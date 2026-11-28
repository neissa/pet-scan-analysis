"""TEP (PET) injecté : série DICOM -> volume SUV -> statistiques par lésion."""
import argparse
import os
import sys
from concurrent.futures import BrokenExecutor, ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pydicom
from scipy import ndimage, signal
from scipy.spatial import ConvexHull

from organs import ORGANS, PHYSIO


def _seconds(tm):  # DICOM TM "HHMMSS[.ffffff]"
    tm = str(tm).split(".")[0].ljust(6, "0")
    return int(tm[:2]) * 3600 + int(tm[2:4]) * 60 + int(tm[4:6])


def suv_factor(weight_kg, dose_bq, half_life_s, delay_s):
    """Bq/ml -> SUV (poids corporel). La dose est corrigée de la décroissance jusqu'au scan."""
    return weight_kg * 1000 / (dose_bq * 2 ** (-delay_s / half_life_s))


def _scaled(progress, a, b):
    """Rappel d'avancement ramené de [0, 1] à [a, b] (None reste None)."""
    return None if progress is None else (lambda f: progress(a + (b - a) * f))


def n_workers():
    """Processus de lecture en parallèle, d'après la machine : ses cœurs utilisables moins un (pour le serveur et la page),
    au plus 8 (au-delà, le démarrage des processus coûte plus que le gain) ; 1 = lecture séquentielle."""
    if hasattr(sys, "getandroidapilevel") or "ANDROID_DATA" in os.environ:
        return 1  # Android : pas de processus fils dans l'appli, lecture séquentielle
    try:
        n = len(os.sched_getaffinity(0))  # cœurs réellement alloués au processus (Linux)
    except AttributeError:
        n = os.cpu_count() or 1
    return max(1, min(n - 1, 8))


_CT_TAGS = ("Modality", "SeriesInstanceUID", "ImagePositionPatient", "ConvolutionKernel", "SeriesDescription")


def _header(p):
    try:
        h = pydicom.dcmread(p, stop_before_pixels=True)  # en-tête seul : pas de pixels
    except pydicom.errors.InvalidDicomError:
        return None
    if h.get("Modality") == "PT":  # le patient et le SUV lisent beaucoup de champs du TEP
        return h
    slim = pydicom.Dataset()  # les autres (surtout le CT, des milliers de fichiers) : seulement ce que la sélection de série utilise
    for t in _CT_TAGS:
        if t in h:
            slim[t] = h[t]
    return slim


_pool = None
_headers_cache = {}  # {(dossier, nb de fichiers): [(chemin, en-tête)]} : dernier dossier lu seulement


def _headers(folder, progress=None):
    """[(chemin, en-tête)] de tous les fichiers DICOM du dossier. Une seule lecture sert au TEP, au CT et au patient
    (elle se fait en parallèle sur plusieurs processus) ; le résultat du dernier dossier est gardé tant qu'il ne change pas."""
    global _pool
    files = [p for p in Path(folder).rglob("*") if p.is_file()]
    sig = (str(folder), len(files))  # un dossier dont on remplace des fichiers sans changer leur nombre n'est pas relu : relancer le serveur
    if sig in _headers_cache:
        return _headers_cache[sig]
    heads = None
    if n_workers() > 1 and len(files) >= 200:  # en dessous, le démarrage des processus ne vaut pas le coup
        try:
            _pool = _pool or ProcessPoolExecutor(n_workers())
            heads = []
            for k, h in enumerate(_pool.map(_header, files, chunksize=32)):
                heads.append(h)
                if progress and k % 25 == 0:
                    progress(k / len(files))
        except (BrokenExecutor, OSError):  # processus indisponibles : lecture séquentielle
            _pool, heads = None, None
    if heads is None:
        heads = []
        for k, p in enumerate(files):
            if progress and k % 25 == 0:
                progress(k / len(files))
            heads.append(_header(p))
    _headers_cache.clear()
    _headers_cache[sig] = out = [(p, h) for p, h in zip(files, heads) if h is not None]
    return out


def _read_all(paths, progress=None):
    """Datasets complets (avec pixels) des fichiers, dans l'ordre ; lus en parallèle comme les en-têtes."""
    global _pool
    out = None
    if n_workers() > 1 and len(paths) >= 50:
        try:
            _pool = _pool or ProcessPoolExecutor(n_workers())
            out = []
            for k, d in enumerate(_pool.map(pydicom.dcmread, paths, chunksize=8)):
                out.append(d)
                if progress and k % 10 == 0:
                    progress(k / len(paths))
        except (BrokenExecutor, OSError):
            _pool, out = None, None
    if out is None:
        out = []
        for k, p in enumerate(paths):
            out.append(pydicom.dcmread(p))
            if progress and k % 10 == 0:
                progress(k / len(paths))
    return out


def _scan(folder, modality, progress=None):
    """{SeriesInstanceUID: [(chemin, en-tête)]} des séries d'une modalité DICOM (fichiers avec ou sans extension).
    progress(f) : avancement de 0 à 1 (facultatif, pour l'interface)."""
    out = {}
    for p, h in _headers(folder, progress):
        if h.get("Modality") == modality:
            out.setdefault(h.SeriesInstanceUID, []).append((p, h))
    return out


def load_suv(folder, progress=None):
    """Retourne (volume SUV [z,y,x], espacement (z,y,x) en mm, origine (z de chaque coupe, y0, x0) en mm).
    progress(f) : avancement de 0 à 1 (facultatif)."""
    # Dossier de plusieurs séries : on garde la plus grande série PET (modalité DICOM "PT").
    series = _scan(folder, "PT", _scaled(progress, 0, 0.4))
    if not series:
        raise SystemExit(f"Aucune série PET dans {folder}")
    paths = max(series.values(), key=len)
    read = _read_all([p for p, _ in paths], _scaled(progress, 0.4, 0.9))
    slices = sorted(read, key=lambda s: float(s.ImagePositionPatient[2]))
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
    # float32 : moitié moins de mémoire et de temps de calcul que float64 (500 Mo -> 250 Mo), précision largement suffisante
    vol = np.stack([(s.pixel_array * float(s.RescaleSlope) + float(s.RescaleIntercept)).astype(np.float32)
                    for s in slices])
    z = np.median(np.diff([float(s.ImagePositionPatient[2]) for s in slices]))
    if [round(v) for v in s0.ImageOrientationPatient] != [1, 0, 0, 0, 1, 0]:
        raise SystemExit("Orientation non axiale non gérée")
    zs = np.array([float(s.ImagePositionPatient[2]) for s in slices])
    origin = (zs, float(s0.ImagePositionPatient[1]), float(s0.ImagePositionPatient[0]))
    if progress:
        progress(1.0)
    return vol * f, (abs(float(z)), float(s0.PixelSpacing[0]), float(s0.PixelSpacing[1])), origin


def ct_series(folder, progress=None):
    """{uid: [(chemin, en-tête)]} des séries CT avec position (écarte le rapport de dose et les captures)."""
    return {u: hs for u, hs in _scan(folder, "CT", progress).items() if "ImagePositionPatient" in hs[0][1]}


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


def load_ct(folder, shape, spacing, origin, uid=None, series=None, progress=None):
    """CT (HU) rééchantillonné sur la grille du TEP [z,y,x] ; -1000 (air) hors du champ du CT.
    uid : série à utiliser (défaut : `default_ct_uid`) ; series : résultat de `ct_series` s'il est déjà scanné.
    Retourne None s'il n'y a pas de CT. Suppose des coupes axiales dans le même repère que le TEP."""
    zs, y0, x0 = origin
    series = series if series is not None else ct_series(folder)
    if not series:
        return None
    uid = uid if uid in series else default_ct_uid(series, origin)
    read = _read_all([p for p, _ in series[uid]], _scaled(progress, 0, 0.6))
    cts = sorted(read, key=lambda s: float(s.ImagePositionPatient[2]))
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
        if progress and k % 10 == 0:
            progress(0.6 + 0.4 * k / len(zs))
        if not ct_z[0] <= z <= ct_z[-1]:
            continue
        f = np.interp(z, ct_z, np.arange(len(ct_z)))  # indice de coupe CT fractionnaire
        i0 = int(f)
        i1 = min(i0 + 1, len(ct_z) - 1)
        sl = (1 - (f - i0)) * vol[i0] + (f - i0) * vol[i1]
        out[k] = ndimage.affine_transform(sl, matrix, offset=offset, output_shape=shape[1:],
                                          order=1, mode="constant", cval=-1000)
    return out


def lbm_james(weight_kg, height_m, sex):
    """Masse maigre (kg), formule de James ; None si le sexe n'est ni M ni F ou si la formule sort de son domaine (IMC très élevé)."""
    a, b = {"M": (1.10, 128), "F": (1.07, 148)}.get(sex, (0, 0))
    lbm = a * weight_kg - b * (weight_kg / (height_m * 100)) ** 2
    return lbm if a and lbm > 0 else None


def sul_info(weight_kg, height_m, sex):
    """Facteur SUV -> SUL (= LBM / poids). Poids, taille (m) ou sexe absents ou invalides : facteur 1 (SUV poids), sul=False."""
    try:
        w, h = float(weight_kg), float(height_m)
        lbm = lbm_james(w, h, str(sex).strip().upper()) if w > 0 and 0 < h <= 3 else None
    except (TypeError, ValueError):
        lbm = None
    if lbm is None:
        return {"sul": False, "factor": 1.0, "lbm_kg": None, "unit": "SUV",
                "label": "SUV poids : taille ou sexe absents des en-têtes, SUL indisponible"}
    return {"sul": True, "factor": lbm / w, "lbm_kg": lbm, "unit": "SUL", "label": "SUL (masse maigre, James)"}


PEAK_R_MM = (3000 / (4 * np.pi)) ** (1 / 3)  # sphère de 1 cm³ : rayon ≈ 6,2 mm


def _sphere(radius_mm, spacing, sub=1):
    """Boule de rayon radius_mm centrée sur la grille (z,y,x) d'espacement `spacing` (mm), de taille impaire.
    sub=1 : masque booléen des centres de voxels ; sub>1 : fraction de chaque voxel dans la boule (sous-échantillonnage sub³),
    pour un volume juste même sur une grille grossière (1 cm³ binaire sur 3×2×2 mm ne fait que 0,88 ml)."""
    o = (np.arange(sub) + 0.5) / sub - 0.5 if sub > 1 else np.zeros(1)
    axes = [((np.arange(-n, n + 1)[:, None] + o) * s).ravel() for s in spacing for n in [int(radius_mm // s) + (sub > 1)]]
    g = np.meshgrid(*axes, indexing="ij", sparse=True)
    inside = sum(a ** 2 for a in g) <= radius_mm ** 2
    if sub == 1:
        return inside
    nz, ny, nx = (len(a) // sub for a in axes)
    return inside.reshape(nz, sub, ny, sub, nx, sub).mean(axis=(1, 3, 5))


def liver_reference(suv, labels, spacing, factor=1.0, diam_mm=30, erode_mm=5, liver_id=4):
    """Référence hépatique PERCIST : moyenne et écart-type (SUL) dans une boule de `diam_mm` placée automatiquement
    entièrement dans le foie érodé de `erode_mm`, là où le SUL est le plus homogène (écart-type minimal : évite une lésion).
    Seuil PERCIST = 1,5 × moyenne + 2 × écart-type. None si le foie est absent ou trop petit."""
    liver = labels == liver_id
    if not liver.any():
        return None
    box = ndimage.find_objects(liver.astype(np.uint8))[0]
    # pad : le bord de la boîte compte comme hors du foie
    d = ndimage.distance_transform_edt(np.pad(liver[box], 1), sampling=spacing)[1:-1, 1:-1, 1:-1]
    r = diam_mm / 2
    cand = d > r + erode_mm
    if not cand.any():
        return None
    ker = _sphere(r, spacing)
    img = suv[box].astype(np.float64) * factor
    s1, s2 = (signal.fftconvolve(v, ker.astype(np.float64), mode="same") for v in (img, img * img))  # noyau symétrique
    var = s2 / ker.sum() - (s1 / ker.sum()) ** 2
    c = np.unravel_index(np.argmin(np.where(cand, var, np.inf)), var.shape)
    half = [n // 2 for n in ker.shape]
    vals = img[tuple(slice(ci - h, ci + h + 1) for ci, h in zip(c, half))][ker]  # valeurs exactes de la boule retenue
    mean, sd = float(vals.mean()), float(vals.std())
    return {"mean": mean, "sd": sd, "seuil": 1.5 * mean + 2 * sd, "diam_mm": diam_mm,
            "centre": [int(ci + b.start) for ci, b in zip(c, box)]}


def analyze(suv, voxel_ml, thr=2.5, spacing=None, factor=1.0, organs=None, margin_mm=10):
    """Lésions = composantes connexes de SUV >= thr. Triées par SUVmax décroissant.
    factor : SUV -> SUL (`sul_info`) ; sul_max, sul_mean, tlg_sul en découlent.
    spacing (z,y,x) mm : ajoute sul_peak (moyenne maximale d'une boule de 1 cm³ centrée dans la lésion).
    organs : étiquettes de `organs.segment_organs` ; ajoute physio (organe à fixation physiologique à moins de
    margin_mm, sinon None) et valide (False si physio) : la lésion est marquée, pas retirée."""
    labels, n = ndimage.label(suv >= thr)
    if not n:
        return []
    # Seuls les voxels des lésions sont triés (quelques millions, pas tout le volume) : ndimage.maximum/maximum_position
    # trient les 60 millions de voxels et prennent ~18 s.
    flat = np.flatnonzero(labels)
    lab, val = labels.ravel()[flat], suv.ravel()[flat]
    counts = np.bincount(lab, minlength=n + 1)[1:]
    sums = np.bincount(lab, weights=val, minlength=n + 1)[1:]
    order = np.lexsort((val, lab))                    # par lésion puis par valeur croissante : le dernier de chaque groupe est le pic
    last = order[np.cumsum(counts) - 1]
    peaks = np.unravel_index(flat[last], suv.shape)   # position (z,y,x) du pic de chaque lésion
    out = [{"suv_max": float(val[last[i]]), "suv_mean": float(sums[i] / counts[i]), "mtv_ml": float(counts[i] * voxel_ml),
            "tlg": float(sums[i] * voxel_ml), "peak": [int(p[i]) for p in peaks]} for i in range(n)]
    for r in out:
        r["sul_max"], r["sul_mean"], r["tlg_sul"] = r["suv_max"] * factor, r["suv_mean"] * factor, r["tlg"] * factor
    if spacing is not None:
        ker = _sphere(PEAK_R_MM, spacing, sub=5).astype(np.float32)
        ker /= ker.sum()
        hp = [n // 2 for n in ker.shape]
        hm = [int(np.ceil(margin_mm / s)) for s in spacing]
        grow = lambda sl, h: tuple(slice(max(0, s.start - k), s.stop + k) for s, k in zip(sl, h))
        for i, sl in enumerate(ndimage.find_objects(labels)):
            w = grow(sl, hp)
            m = ndimage.correlate(suv[w].astype(np.float32), ker, mode="nearest")
            out[i]["sul_peak"] = float(m[labels[w] == i + 1].max()) * factor
            if organs is not None:
                w = grow(sl, hm)
                org = organs[w]
                hit = np.isin(org, PHYSIO)
                near = org[hit & (ndimage.distance_transform_edt(labels[w] != i + 1, sampling=spacing) <= margin_mm)] \
                    if hit.any() else org[:0]
                out[i]["physio"] = ORGANS[int(np.bincount(near).argmax())][0] if near.size else None
                out[i]["valide"] = out[i]["physio"] is None
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
