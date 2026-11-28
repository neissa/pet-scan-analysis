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

from organs import ORGANS, PHYSIO, PHYSIO_MARGIN_MM


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
    Seuil PERCIST = 1,5 × moyenne + 2 × écart-type (seuil_suv : le même en SUV, ÷ factor).
    Foie segmenté trop étroit : repli sur une marge nulle, puis sur une boule de 20 mm (sphere_reduite) ; la marge et
    le diamètre réellement utilisés sont rendus (erode_mm, diam_mm) avec un message (msg, vide sans repli).
    None si le foie est absent ou trop petit même pour 20 mm."""
    liver = labels == liver_id
    if not liver.any():
        return None
    box = ndimage.find_objects(liver.astype(np.uint8))[0]
    # pad : le bord de la boîte compte comme hors du foie
    d = ndimage.distance_transform_edt(np.pad(liver[box], 1), sampling=spacing)[1:-1, 1:-1, 1:-1]
    for dm, em in ((diam_mm, erode_mm), (diam_mm, 0), (min(diam_mm, 20), 0)):
        cand = d > dm / 2 + em
        if cand.any():
            break
    else:
        return None
    msg = "" if em == erode_mm else f"Foie segmenté étroit : marge d'érosion {erode_mm} → 0 mm" + (
        "" if dm == diam_mm else f", sphère {diam_mm} → {dm} mm")
    reduite, diam_mm, r = dm < diam_mm, dm, dm / 2
    ker = _sphere(r, spacing)
    img = suv[box].astype(np.float64) * factor
    s1, s2 = (signal.fftconvolve(v, ker.astype(np.float64), mode="same") for v in (img, img * img))  # noyau symétrique
    var = s2 / ker.sum() - (s1 / ker.sum()) ** 2
    c = np.unravel_index(np.argmin(np.where(cand, var, np.inf)), var.shape)
    half = [n // 2 for n in ker.shape]
    vals = img[tuple(slice(ci - h, ci + h + 1) for ci, h in zip(c, half))][ker]  # valeurs exactes de la boule retenue
    mean, sd = float(vals.mean()), float(vals.std())
    seuil = 1.5 * mean + 2 * sd
    return {"mean": mean, "sd": sd, "seuil": seuil, "seuil_suv": seuil / factor, "diam_mm": diam_mm, "erode_mm": em,
            "sphere_reduite": reduite, "msg": msg, "centre": [int(ci + b.start) for ci, b in zip(c, box)]}


def blood_pool_reference(suv, labels, spacing, factor=1.0, erode_mm=4, aorta_id=11):
    """Référence pool sanguin PERCIST (repli quand le foie est atteint) : SUL moyen et écart-type dans l'aorte thoracique
    descendante segmentée, érodée de `erode_mm` (paroi et volume partiel). Seuil PERCIST = 2 × moyenne + 2 × écart-type.
    shortcut: tout le segment d'aorte trouvé au lieu du cylindre de 1 × 2 cm de PERCIST. None sans aorte assez large."""
    ao = labels == aorta_id
    if not ao.any():
        return None
    box = ndimage.find_objects(ao.astype(np.uint8))[0]
    core = ndimage.distance_transform_edt(np.pad(ao[box], 1), sampling=spacing)[1:-1, 1:-1, 1:-1] > erode_mm
    if core.sum() < 10:
        return None
    vals = suv[box][core].astype(np.float64) * factor
    mean, sd = float(vals.mean()), float(vals.std())
    seuil = 2 * mean + 2 * sd
    return {"mean": mean, "sd": sd, "seuil": seuil, "seuil_suv": seuil / factor, "erode_mm": erode_mm}


def analyze(suv, voxel_ml, thr=2.5, spacing=None, factor=1.0, organs=None, margin_mm=10):
    """Lésions = composantes connexes de SUV >= thr. Triées par SUVmax décroissant.
    factor : SUV -> SUL (`sul_info`) ; sul_max, sul_mean, tlg_sul en découlent.
    spacing (z,y,x) mm : ajoute sul_peak (moyenne maximale d'une boule de 1 cm³ centrée dans la lésion).
    organs (avec spacing) : étiquettes de `organs.segment_organs`. Exclusion par voxel : les voxels à moins de margin_mm
    (plafonnée par organe, `PHYSIO_MARGIN_MM` : 5 mm pour le cœur) d'un organe à fixation physiologique (PHYSIO) sont retirés avant l'étiquetage ; une maladie collée à la vessie ou à un
    rein garde ainsi ses autres voxels. physio = organe au contact de la zone retirée (sinon None), simple marque ;
    valide reste True (une lésion entièrement dans la zone retirée n'existe plus)."""
    mask = suv >= thr
    grow = lambda sl, h: tuple(slice(max(0, s.start - k), s.stop + k) for s, k in zip(sl, h))
    if organs is not None and spacing is not None:
        for oid, sl in enumerate(ndimage.find_objects(organs), 1):
            if sl is not None and oid in PHYSIO:  # distance calculée dans la boîte de l'organe seulement (rapide)
                mo = min(margin_mm, PHYSIO_MARGIN_MM.get(oid, margin_mm))
                w = grow(sl, [int(np.ceil(mo / s)) + 1 for s in spacing])
                mask[w] &= ndimage.distance_transform_edt(organs[w] != oid, sampling=spacing) > mo
    labels, n = ndimage.label(mask)
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
        touch = margin_mm + max(spacing)  # au contact de la zone retirée : à moins d'un voxel de plus que la marge
        hm = [int(np.ceil(touch / s)) for s in spacing]
        for i, sl in enumerate(ndimage.find_objects(labels)):
            w = grow(sl, hp)
            m = ndimage.correlate(suv[w].astype(np.float32), ker, mode="nearest")
            out[i]["sul_peak"] = float(m[labels[w] == i + 1].max()) * factor
            if organs is not None:
                w = grow(sl, hm)
                org = organs[w]
                hit = np.isin(org, PHYSIO)
                near = org[hit & (ndimage.distance_transform_edt(labels[w] != i + 1, sampling=spacing) <= touch)] \
                    if hit.any() else org[:0]
                out[i]["physio"] = ORGANS[int(np.bincount(near).argmax())][0] if near.size else None
                out[i]["valide"] = True
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


URINARY = (6, 7, 8)  # reins et vessie : l'urine excrétée y passe


def ureter_paths(labels, spacing):
    """Trajets urétéraux approximatifs : segments (mm) du centre de chaque rein au centre de la vessie.
    shortcut: droite rein -> vessie, pas l'uretère réel ; à remplacer si l'uretère est un jour segmenté."""
    boxes = ndimage.find_objects(labels, max_label=max(URINARY))
    com = {}
    for oid in URINARY:
        sl = boxes[oid - 1]
        if sl is not None:
            c = ndimage.center_of_mass(labels[sl] == oid)
            com[oid] = np.array([(ci + s.start) * sp for ci, s, sp in zip(c, sl, spacing)])
    return [(com[k], com[8]) for k in (6, 7) if k in com and 8 in com]


def _seg_dist(p, a, b):
    ab = b - a
    t = np.clip(np.dot(p - a, ab) / (np.dot(ab, ab) or 1.0), 0, 1)
    return float(np.linalg.norm(p - a - t * ab))


def measure_at(suv, spacing, pos, factor=1.0, labels=None, ureters=(), radius_mm=12, grow_mm=60, pct=50,
               excr_suv=20, excr_mm=20):
    """Mesure d'une lésion mémorisée à la position pos (z,y,x), SANS seuil global : SUVmax dans une boule de radius_mm,
    voxels des organes PHYSIO exclus (urine de la vessie, reins, cœur, cerveau) ; volume = voxels connexes >= pct % de ce max
    (fenêtre ±grow_mm) ; SULpeak (boule de 1 cm³) dans ce volume ; organes à ≤ 5 mm (« au contact »). Drapeau « excrétion possible »
    (la lésion reste mesurée) si le max dépasse excr_suv ou si le pic est à moins de excr_mm d'un rein, de la vessie ou
    d'un trajet urétéral (`ureter_paths`). None si pos est hors du volume."""
    pos = [int(round(p)) for p in pos]
    if not all(0 <= p < n for p, n in zip(pos, suv.shape)):
        return None
    h = [int(np.ceil(grow_mm / s)) for s in spacing]
    lo = [max(0, p - k) for p, k in zip(pos, h)]
    sl = tuple(slice(l, p + k + 1) for l, p, k in zip(lo, pos, h))
    w = suv[sl].astype(np.float32)  # copie
    lab = labels[sl] if labels is not None else None
    if lab is not None:
        w[np.isin(lab, PHYSIO)] = 0
    g = np.ogrid[tuple(slice(0, n) for n in w.shape)]
    c = [p - l for p, l in zip(pos, lo)]
    near = sum(((a - ci) * s) ** 2 for a, ci, s in zip(g, c, spacing)) <= radius_mm ** 2
    vals = np.where(near, w, -np.inf)
    pk = np.unravel_index(int(np.argmax(vals)), w.shape)
    vmax = float(w[pk])
    out = {"suv_max": vmax, "sul_max": vmax * factor, "peak": [int(a + l) for a, l in zip(pk, lo)],
           "dist_mm": float(np.sqrt(sum(((a - ci) * s) ** 2 for a, ci, s in zip(pk, c, spacing)))), "flags": [], "contact": []}
    if vmax <= 0:
        out.update(sul_peak=0.0, vol_ml=0.0)
        return out
    rl, _ = ndimage.label(w >= pct / 100 * vmax)
    region = rl == rl[pk]
    out["vol_ml"] = float(region.sum() * np.prod(spacing) / 1000)
    if any(np.take(region, [0, -1], axis=a).any() for a in range(3)):
        out["flags"].append("volume limité à la fenêtre")
    ker = _sphere(PEAK_R_MM, spacing, sub=5).astype(np.float32)
    ker /= ker.sum()
    box = tuple(slice(max(0, b.start - k // 2), b.stop + k // 2) for b, k in zip(ndimage.find_objects(region.astype(np.uint8))[0], ker.shape))
    out["sul_peak"] = float(ndimage.correlate(w[box], ker, mode="nearest")[region[box]].max()) * factor
    if lab is not None:
        ring = ndimage.binary_dilation(region, iterations=max(1, round(5 / min(spacing)))) & ~region  # « au contact » : ≤ ~5 mm
        out["contact"] = [ORGANS[int(k)][0] for k in np.unique(lab[ring]) if int(k) in ORGANS and k != 10]
        pmm = np.array(out["peak"], float) * spacing
        uri = np.argwhere(np.isin(lab, URINARY))
        close = uri.size and float(np.sqrt((((uri - pk) * spacing) ** 2).sum(axis=1).min())) < excr_mm
        if vmax > excr_suv or close or any(_seg_dist(pmm, a, b) < excr_mm for a, b in ureters):
            out["flags"].append("excrétion possible")
    elif vmax > excr_suv:
        out["flags"].append("excrétion possible")
    return out


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


def ct_diameter(ctv, spacing, pos, vol_ml=None, tol=50):
    """Diamètres RECIST au CT dans le plan axial à la position pos (z,y,x) : plus grand diamètre et petit axe (étendue
    perpendiculaire) de la région par tolérance HU (`measure_tolerance`). None sans CT ; {"why": raison} si la mesure n'est
    pas fiable (densité non tissulaire, région non bornée ou qui fuit bien au-delà de la lésion TEP de volume vol_ml)."""
    if ctv is None:
        return None
    z, y, x = (int(round(p)) for p in pos)
    if not all(0 <= p < n for p, n in zip((z, y, x), ctv.shape)):
        return {"why": "hors du CT"}
    img = ctv[z]
    r = measure_tolerance(img, spacing[1], spacing[2], y, x, tol)
    if not -300 <= r["hu_ref"] <= 300:
        return {"why": f"densité non tissulaire ({r['hu_ref']:.0f} HU)"}
    edge = np.array(r["boundary"])
    if ((edge == 0) | (edge == np.array(img.shape) - 1)).any():
        return {"why": "région non bornée (atteint le bord de la coupe)"}
    d_pet = 20 * (3 * vol_ml / (4 * np.pi)) ** (1 / 3) if vol_ml else None  # diamètre de la sphère de même volume
    if r["length_mm"] > 300 or d_pet and r["length_mm"] > max(2.5 * d_pet, d_pet + 30):
        return {"why": "région fuyante (tissus voisins de même densité)"}
    return {"long_mm": round(r["length_mm"], 1), "short_mm": round(r["width_mm"], 1)}


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
