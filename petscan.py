"""TEP (PET) injecté : série DICOM -> volume SUV -> statistiques par lésion."""
import argparse
from pathlib import Path

import numpy as np
import pydicom
from scipy import ndimage
from scipy.spatial import ConvexHull


def _seconds(tm):  # DICOM TM "HHMMSS[.ffffff]"
    tm = str(tm).split(".")[0].ljust(6, "0")
    return int(tm[:2]) * 3600 + int(tm[2:4]) * 60 + int(tm[4:6])


def suv_factor(weight_kg, dose_bq, half_life_s, delay_s):
    """Bq/ml -> SUV (poids corporel). La dose est corrigée de la décroissance jusqu'au scan."""
    return weight_kg * 1000 / (dose_bq * 2 ** (-delay_s / half_life_s))


def _scan(folder, modality):
    """{SeriesInstanceUID: [(chemin, en-tête)]} des séries d'une modalité DICOM (fichiers avec ou sans extension)."""
    out = {}
    for p in Path(folder).rglob("*"):
        if not p.is_file():
            continue
        try:
            h = pydicom.dcmread(p, stop_before_pixels=True)
        except pydicom.errors.InvalidDicomError:
            continue
        if h.get("Modality") == modality:
            out.setdefault(h.SeriesInstanceUID, []).append((p, h))
    return out


def load_suv(folder):
    """Retourne (volume SUV [z,y,x], espacement (z,y,x) en mm, origine (z de chaque coupe, y0, x0) en mm)."""
    # Dossier de plusieurs séries : on garde la plus grande série PET (modalité DICOM "PT").
    series = _scan(folder, "PT")
    if not series:
        raise SystemExit(f"Aucune série PET dans {folder}")
    paths = max(series.values(), key=len)
    slices = sorted((pydicom.dcmread(p) for p, _ in paths),
                    key=lambda s: float(s.ImagePositionPatient[2]))
    s0 = slices[0]
    if s0.Units != "BQML":
        raise SystemExit(f"Unités {s0.Units} non gérées (BQML attendu)")
    rp = s0.RadiopharmaceuticalInformationSequence[0]
    start = min(_seconds(s.AcquisitionTime) for s in slices)  # début du scan
    delay = start - _seconds(rp.RadiopharmaceuticalStartTime)
    if delay < 0:  # passage de minuit
        delay += 86400
    f = suv_factor(float(s0.PatientWeight), float(rp.RadionuclideTotalDose),
                   float(rp.RadionuclideHalfLife), delay)
    vol = np.stack([s.pixel_array * float(s.RescaleSlope) + float(s.RescaleIntercept)
                    for s in slices])
    z = np.median(np.diff([float(s.ImagePositionPatient[2]) for s in slices]))
    if [round(v) for v in s0.ImageOrientationPatient] != [1, 0, 0, 0, 1, 0]:
        raise SystemExit("Orientation non axiale non gérée")
    zs = np.array([float(s.ImagePositionPatient[2]) for s in slices])
    origin = (zs, float(s0.ImagePositionPatient[1]), float(s0.ImagePositionPatient[0]))
    return vol * f, (abs(float(z)), float(s0.PixelSpacing[0]), float(s0.PixelSpacing[1])), origin


def ct_series(folder):
    """{uid: [(chemin, en-tête)]} des séries CT avec position (écarte le rapport de dose et les captures)."""
    return {u: hs for u, hs in _scan(folder, "CT").items() if "ImagePositionPatient" in hs[0][1]}


def _ct_overlap(hs, zs):
    z = [float(h.ImagePositionPatient[2]) for _, h in hs]
    return max(0.0, min(max(z), zs[-1]) - max(min(z), zs[0]))


def ct_choices(series, origin):
    """Séries CT proposées : [{uid, desc, n, couverture}], couverture = part de la hauteur du TEP couverte par le CT."""
    zs = origin[0]
    span = (zs[-1] - zs[0]) or 1
    out = [{"uid": u, "desc": str(hs[0][1].get("SeriesDescription", "CT")), "n": len(hs),
            "couverture": round(_ct_overlap(hs, zs) / span, 2)} for u, hs in series.items() if len(hs) > 10]
    return sorted(out, key=lambda c: (-c["couverture"], -c["n"], c["desc"]))


def default_ct_uid(series, origin):
    """Série qui recouvre le mieux le TEP ; à égalité, noyau STANDARD (tissus mous)."""
    return max(series, key=lambda u: (_ct_overlap(series[u], origin[0]), series[u][0][1].get("ConvolutionKernel") == "STANDARD"))


def load_ct(folder, shape, spacing, origin, uid=None, series=None):
    """CT (HU) rééchantillonné sur la grille du TEP [z,y,x] ; -1000 (air) hors du champ du CT.
    uid : série à utiliser (défaut : `default_ct_uid`) ; series : résultat de `ct_series` s'il est déjà scanné.
    Retourne None s'il n'y a pas de CT. Suppose des coupes axiales dans le même repère que le TEP."""
    zs, y0, x0 = origin
    series = series if series is not None else ct_series(folder)
    if not series:
        return None
    uid = uid if uid in series else default_ct_uid(series, origin)
    cts = sorted((pydicom.dcmread(p) for p, _ in series[uid]), key=lambda s: float(s.ImagePositionPatient[2]))
    if [round(v) for v in cts[0].ImageOrientationPatient] != [1, 0, 0, 0, 1, 0]:
        raise SystemExit("CT : orientation non axiale non gérée")
    vol = np.stack([(s.pixel_array * float(s.RescaleSlope) + float(s.RescaleIntercept)).astype(np.int16)
                    for s in cts])
    ct_z = np.array([float(s.ImagePositionPatient[2]) for s in cts])
    dy, dx = (float(v) for v in cts[0].PixelSpacing)
    # indice CT (ligne, colonne) = matrice * indice TEP + décalage
    matrix = [spacing[1] / dy, spacing[2] / dx]
    offset = [(y0 - float(cts[0].ImagePositionPatient[1])) / dy,
              (x0 - float(cts[0].ImagePositionPatient[0])) / dx]
    out = np.full(shape, -1000, np.int16)
    for k, z in enumerate(zs):
        if not ct_z[0] <= z <= ct_z[-1]:
            continue
        f = np.interp(z, ct_z, np.arange(len(ct_z)))  # indice de coupe CT fractionnaire
        i0 = int(f)
        i1 = min(i0 + 1, len(ct_z) - 1)
        sl = (1 - (f - i0)) * vol[i0] + (f - i0) * vol[i1]
        out[k] = ndimage.affine_transform(sl, matrix, offset=offset, output_shape=shape[1:],
                                          order=1, mode="constant", cval=-1000)
    return out


def analyze(suv, voxel_ml, thr=2.5):
    """Lésions = composantes connexes de SUV >= thr. Triées par SUVmax décroissant."""
    labels, n = ndimage.label(suv >= thr)
    out = []
    peaks = ndimage.maximum_position(suv, labels, range(1, n + 1)) if n else []  # position (z,y,x) du pic de chaque lésion
    for i in range(1, n + 1):
        v = suv[labels == i]
        mtv = v.size * voxel_ml
        out.append({"suv_max": float(v.max()), "suv_mean": float(v.mean()),
                    "mtv_ml": mtv, "tlg": float(v.mean() * mtv), "peak": [int(c) for c in peaks[i - 1]]})
    return sorted(out, key=lambda r: -r["suv_max"])


def measure_region(img, sr, sc, r, c, pct=50, radius_mm=8):
    """Mesure automatique sur une coupe 2D autour du pixel (r, c).
    Graine = max local (rayon radius_mm) ; région = pixels connexes >= pct % de ce max.
    sr, sc : taille des pixels en mm (lignes, colonnes).
    length_mm = plus grand diamètre (Feret max) ; width_mm = étendue perpendiculaire à ce diamètre.
    axes/perp = extrémités des deux segments, boundary = pixels de bord, en indices (ligne, colonne).
    None si le max local est <= 0."""
    nr, nc = img.shape
    wr, wc = round(radius_mm / sr), round(radius_mm / sc)
    r0, c0 = max(0, r - wr), max(0, c - wc)
    win = img[r0:r + wr + 1, c0:c + wc + 1]
    dr, dc = np.unravel_index(win.argmax(), win.shape)
    seed = (r0 + dr, c0 + dc)
    peak = float(img[seed])
    if peak <= 0:
        return None
    lab, _ = ndimage.label(img >= min(max(pct, 1), 100) / 100 * peak)
    mask = lab == lab[seed]
    res = _describe(mask, img, sr, sc)
    res["suv_max"], res["suv_mean"] = res.pop("vmax"), res.pop("vmean")
    return res


def measure_tolerance(img, sr, sc, r, c, tol=50):
    """Mesure sur le CT autour du pixel (r, c) : région connexe dont le HU est à ± tol du HU local
    (médiane 3x3). Un seuil de tolérance trop large fait « fuir » la région dans les tissus voisins de même densité."""
    win = img[max(0, r - 1):r + 2, max(0, c - 1):c + 2].astype(np.float32)
    ref = float(np.median(win))
    dr, dc = np.unravel_index(np.abs(win - ref).argmin(), win.shape)  # un pixel de la fenêtre vaut la médiane
    seed = (max(0, r - 1) + dr, max(0, c - 1) + dc)
    lab, _ = ndimage.label(np.abs(img.astype(np.float32) - ref) <= tol)
    res = _describe(lab == lab[seed], img, sr, sc)
    res["hu_ref"], res["hu_mean"], res["hu_max"] = ref, res.pop("vmean"), res.pop("vmax")
    return res


def _describe(mask, img, sr, sc):
    """Géométrie d'une région 2D : plus grand diamètre (Feret max), étendue perpendiculaire, contour, aire, valeurs."""
    edge = np.argwhere(mask & ~ndimage.binary_erosion(mask, border_value=0))
    # coins des pixels de bord : un pixel isolé mesure sr x sc, pas 0
    corners = (edge[:, None, :] + np.array([[-.5, -.5], [-.5, .5], [.5, -.5], [.5, .5]])).reshape(-1, 2)
    mm = corners * [sr, sc]
    hull = mm[ConvexHull(mm).vertices]
    d = np.linalg.norm(hull[:, None] - hull[None], axis=-1)
    a, b = np.unravel_index(d.argmax(), d.shape)
    length = float(d[a, b])
    u = (hull[b] - hull[a]) / length
    n = np.array([-u[1], u[0]])
    proj = hull @ n
    mid = (hull[a] + hull[b]) / 2
    base = mid - n * (mid @ n)  # même position le long du grand axe que le milieu
    vals = img[mask]
    k = np.array([sr, sc])
    return {
        "length_mm": length, "width_mm": float(proj.max() - proj.min()),
        "vmax": float(vals.max()), "vmean": float(vals.mean()), "area_mm2": float(mask.sum() * sr * sc),
        "axes": [(hull[a] / k).tolist(), (hull[b] / k).tolist()],
        "perp": [((base + n * proj.min()) / k).tolist(), ((base + n * proj.max()) / k).tolist()],
        "boundary": edge[:: max(1, len(edge) // 3000)].tolist(),
    }


def save_mip(suv, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.imsave(path, suv.max(axis=1)[::-1], cmap="gray_r", vmin=0, vmax=10)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder", help="dossier de la série PET DICOM")
    ap.add_argument("--thr", type=float, default=2.5, help="seuil SUV (défaut 2.5)")
    ap.add_argument("--mip", help="écrit une MIP coronale PNG")
    a = ap.parse_args()
    suv, spacing, _ = load_suv(a.folder)
    voxel_ml = float(np.prod(spacing)) / 1000
    print(f"SUVmax global : {suv.max():.2f}")
    print(f"{'#':>2} {'SUVmax':>7} {'SUVmean':>8} {'MTV ml':>8} {'TLG':>8}")
    for i, r in enumerate(analyze(suv, voxel_ml, a.thr), 1):
        print(f"{i:>2} {r['suv_max']:>7.2f} {r['suv_mean']:>8.2f} {r['mtv_ml']:>8.1f} {r['tlg']:>8.1f}")
    if a.mip:
        save_mip(suv, a.mip)


if __name__ == "__main__":
    main()
