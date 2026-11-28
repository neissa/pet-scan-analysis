"""Visionneuse web TEP : python viewer.py [dossier_dicom]  (ou --demo)."""
import argparse
import multiprocessing
import importlib.util
import io
import json
import os
import subprocess
import sys
import threading
import webbrowser
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage

from compare import estimate_shift, refine_shift, shift_volume
from export_pdf import build_pdf
from organs import ORGANS, segment_organs
from patient import read_patient_info, read_sul
from petscan import (_scaled, analyze, blood_pool_reference, ct_choices, ct_series, default_ct_uid, liver_reference, load_ct, load_suv,
                     measure_region, measure_tolerance, sul_info)

AXIS = {"axial": 0, "coronal": 1, "sagittal": 2}
CONFIG = Path(os.environ.get("PETSCAN_CONFIG") or Path.home() / ".pet-scan-analysis.json")  # PETSCAN_CONFIG : appli Android ; dossiers récents, dans le dossier utilisateur (jamais dans le repo)
BROWSE = importlib.util.find_spec("tkinter") is not None  # sélecteur de dossier natif disponible ?


def demo():
    """Volume synthétique avec deux « lésions », pour essayer sans DICOM."""
    z, y, x = np.mgrid[:60, :80, :80]
    suv = 1 + 0.5 * np.random.default_rng(0).random((60, 80, 80))
    for (cz, cy, cx), r, v in [((30, 30, 30), 5, 9), ((20, 50, 55), 3, 5)]:
        suv[(z - cz) ** 2 + (y - cy) ** 2 + (x - cx) ** 2 <= r * r] = v
    ctv = np.where(((z - 30) / 28) ** 2 + ((y - 40) / 32) ** 2 + ((x - 40) / 32) ** 2 < 1, 40, -1000)
    return suv, (3.0, 3.0, 3.0), ctv.astype(np.int16)


CMAPS = {"hot", "jet", "gray_r", "bwr"}  # bwr : échelle divergente de la différence entre séances


def bar(cmap):
    """Dégradé 0 -> vmax, pour la légende."""
    buf = io.BytesIO()
    plt.imsave(buf, np.linspace(0, 1, 256)[None, :], cmap=cmap, format="png")
    return buf.getvalue()


UP = 2  # sur-échantillonnage (le navigateur réduit ensuite : plus net qu'un agrandissement)
ROWCOL = {"axial": (1, 2), "coronal": (0, 2), "sagittal": (0, 1), "mip": (0, 2)}  # axes (ligne, colonne)


def _slice(vol, view, i):
    if view == "mip":
        return vol.max(axis=1)[::-1]
    img = np.take(vol, i, axis=AXIS[view], mode="clip")
    return img if view == "axial" else img[::-1]  # tête en haut


def _shift2d(img, dr, dc, fill):
    """out[r, c] = img[r+dr, c+dc] ; `fill` hors de l'image."""
    if dr == 0 and dc == 0:
        return img
    out = np.full_like(img, fill)
    r0, r1 = max(0, -dr), min(img.shape[0], img.shape[0] - dr)
    c0, c1 = max(0, -dc), min(img.shape[1], img.shape[1] - dc)
    if r1 > r0 and c1 > c0:
        out[r0:r1, c0:c1] = img[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
    return out


def _slice_off(vol, view, i, off, fill=0):
    """Comme `_slice` pour un volume décalé de off = (dz, dy, dx) voxels : out[z, y, x] = vol[z+dz, y+dy, x+dx]
    (même convention que `compare.shift_volume`, mais sans copier le volume : on ne décale que la coupe affichée)."""
    if not any(off):
        return _slice(vol, view, i)
    if view == "mip":
        return _shift2d(vol.max(axis=1), off[0], off[2], fill)[::-1]
    ax = AXIS[view]
    k = i + off[ax]
    if 0 <= k < vol.shape[ax]:
        rest = [off[a] for a in range(3) if a != ax]
        img = _shift2d(np.take(vol, k, axis=ax), rest[0], rest[1], fill)
    else:  # coupe hors du volume décalé
        img = np.full(tuple(n for a, n in enumerate(vol.shape) if a != ax), fill, vol.dtype)
    return img if view == "axial" else img[::-1]


def _up(img, view, spacing, box=(0, 0, 1), order=3):
    """Rééchantillonnage cubique à l'espacement réel des voxels (spacing en mm, z,y,x).
    box = (x0, y0, w) : on ne garde que ce carré (fractions de l'image entière) et il remplit
    toute l'image de sortie, donc un zoom augmente vraiment la résolution."""
    r, c = ROWCOL[view]
    m = min(spacing)
    nr, nc = img.shape
    out = (round(UP * nr * spacing[r] / m), round(UP * nc * spacing[c] / m))
    x0, y0, w = box
    sr, sc = w * nr / out[0], w * nc / out[1]  # pas en pixels source par pixel de sortie
    return ndimage.affine_transform(img.astype(np.float32), [sr, sc], output_shape=out, order=order, mode="nearest",
                                    offset=[y0 * nr - 0.5 + sr / 2, x0 * nc - 0.5 + sc / 2])


def measure(suv, spacing, view, i, fx, fy, pct, ctv=None, source="pet", tol=50, mip=None, off=(0, 0, 0)):
    """Mesure auto autour du point (fx, fy), fractions de l'image entière ; résultats en fractions.
    source "pet" : région >= pct % du max local ; source "ct" : région à ± tol HU du HU local.
    mip : image d'un MIP tourné (spacing alors celui de cette image) ; le MIP est toujours mesuré sur le TEP."""
    if view == "mip":
        source = "pet"
    if source == "ct" and ctv is None:
        return {"error": "pas de CT dans ce dossier"}
    img = mip if (view == "mip" and mip is not None) else _slice_off(ctv if source == "ct" else suv, view, i, off,
                                                                      -1000 if source == "ct" else 0)
    nr, nc = img.shape
    rr, cc = ROWCOL[view]
    r, c = min(nr - 1, max(0, int(fy * nr))), min(nc - 1, max(0, int(fx * nc)))
    res = (measure_tolerance(img, spacing[rr], spacing[cc], r, c, tol) if source == "ct"
           else measure_region(img, spacing[rr], spacing[cc], r, c, pct))
    if res is None:
        return {"error": "aucune activité autour du point"}
    frac = lambda pts: [[(c + 0.5) / nc, (r + 0.5) / nr] for r, c in pts]  # (ligne, col) -> [x, y]
    res["axes"], res["perp"], res["boundary"] = frac(res["axes"]), frac(res["perp"]), frac(res["boundary"])
    return res


_POOL = ThreadPoolExecutor(max_workers=8)
_LUTS = {}


def _lut(cmap):
    """Table de 256 couleurs (uint8) d'une palette, calculée une seule fois."""
    if cmap not in _LUTS:
        _LUTS[cmap] = (plt.get_cmap(cmap)(np.linspace(0, 1, 256))[:, :3] * 255).round().astype(np.uint8)
    return _LUTS[cmap]


def _png(rgb):
    """PNG d'une image uint8 (h, w, 3). Compression niveau 1 : 3 à 4 fois plus rapide, fichier ~2 fois plus gros (échange local)."""
    buf = io.BytesIO()
    Image.fromarray(rgb).save(buf, "PNG", compress_level=1)
    return buf.getvalue()


def _organs_overlay(rgb, lab):
    """Teinte translucide + contour de chaque organe (rgb : float32 en 0..255 ; lab : étiquettes de la même taille)."""
    for oid in np.unique(lab):
        if oid in ORGANS:
            col = np.array(matplotlib.colors.to_rgb(ORGANS[oid][1]), np.float32) * 255
            m = lab == oid
            rgb[m] = 0.72 * rgb[m] + 0.28 * col
            rgb[m & ~ndimage.binary_erosion(m, iterations=2)] = col
    return rgb


def _add_b(lesions, bsuv, spacing, radius_mm=12, rel=(0, 0, 0)):
    """Évolution des lésions : pour chacune (séance A), SUVmax de la séance B autour du même point et variation en %.
    Le recalage est approximatif (translation seule) : on cherche le max dans un rayon de `radius_mm`.
    rel = décalage manuel de B moins celui de A (voxels) : le point de B correspondant au pic p de A est p + rel."""
    r = [max(1, round(radius_mm / s)) for s in spacing]
    for lesion in lesions:
        z, y, x = (p + d for p, d in zip(lesion["peak"], rel))
        lo = [max(0, c - n) for c, n in zip((z, y, x), r)]
        hi = [max(0, c + n + 1) for c, n in zip((z, y, x), r)]  # max(0, …) : une fenêtre entièrement hors volume reste vide
        win = bsuv[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
        sb = float(win.max()) if win.size else 0.0
        lesion["suv_b"] = sb
        lesion["delta_pct"] = (lesion["suv_max"] - sb) / sb * 100 if sb >= 0.5 else None
    return lesions


def _organ_at_peak(labels, spacing, peak, radius_mm=15):
    """Organe du pic (« dans X ») ou organe le plus proche à ≤ radius_mm (« près de X », squelette exclu), sinon None.
    Calcul local : seule la fenêtre autour du pic est lue, pas tout le volume."""
    oid = int(labels[tuple(peak)])
    if oid in ORGANS:
        return f"dans {ORGANS[oid][0]}"
    r = [int(np.ceil(radius_mm / s)) for s in spacing]
    lo = [max(0, c - n) for c, n in zip(peak, r)]
    hi = [min(m, c + n + 1) for c, n, m in zip(peak, r, labels.shape)]
    win = labels[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    centre = np.array(peak) - lo
    best = None
    for k in np.unique(win):
        if int(k) not in ORGANS or k == 10:  # le squelette est trop répandu pour dire « près de »
            continue
        pts = np.argwhere(win == k)
        d = float(np.sqrt((((pts - centre) * spacing) ** 2).sum(axis=1).min()))  # distance en mm au voxel le plus proche
        if d <= radius_mm and (best is None or d < best[0]):
            best = (d, int(k))
    return f"près de {ORGANS[best[1]][0]}" if best else None


ANGLES = range(0, 360, 10)  # angles du MIP rotatif


def _angle(i):
    return round(i / 10) * 10 % 360


def render_rgb(suv, view, i, vmax, cmap="jet", spacing=(1, 1, 1), ctv=None, mode="pet", alpha=0.6, box=(0, 0, 1), labels=None,
           mip=None, off=(0, 0, 0)):
    """PNG d'une coupe (axial/coronal/sagittal) ou de la MIP coronale (toujours TEP seul). i est borné.
    mode : "pet", "ct" ou "fusion" (TEP en couleur sur CT ; opacité nulle à SUV 0, pleine dès 25 % de vmax).
    Fenêtre CT fixe -160..240 HU (tissus mous). labels : volume d'organes à superposer (sauf sur la MIP)."""
    # mip : image d'un MIP tourné déjà calculée (spacing = celui de cette image) ; sinon MIP de face au pleine résolution
    show_ct = ctv is not None and mode in ("fusion", "ct") and view != "mip"
    show_organs = labels is not None and view != "mip"
    # les interpolations TEP et CT sont indépendantes : en parallèle (scipy libère le verrou Python)
    # off : décalage manuel (dz, dy, dx) en voxels de CE volume (TEP, CT et organes bougent ensemble)
    pet_f = _POOL.submit(_up, mip if (view == "mip" and mip is not None) else _slice_off(suv, view, i, off), view, spacing, box)
    ct_f = _POOL.submit(_up, _slice_off(ctv, view, i, off, -1000), view, spacing, box) if show_ct else None
    pet = pet_f.result()
    idx = (np.clip(pet / vmax, 0, 1) * 255).astype(np.uint8)  # indice dans la palette
    if not (show_ct or show_organs):
        return _lut(cmap)[idx]
    if show_ct:
        img = np.clip((ct_f.result() + 160) * (255 / 400), 0, 255)[..., None] * np.ones(3, np.float32)
        if mode == "fusion":
            a = (alpha * np.clip(pet / (0.25 * vmax), 0, 1))[..., None]
            img = img * (1 - a) + _lut(cmap)[idx] * a
    else:
        img = _lut(cmap)[idx].astype(np.float32)
    if show_organs:
        img = _organs_overlay(img, np.rint(_up(_slice_off(labels, view, i, off), view, spacing, box, order=0)).astype(np.uint8))
    return np.clip(img, 0, 255).astype(np.uint8)


def render(*args, **kw):
    """PNG d'une coupe ; mêmes paramètres que `render_rgb`."""
    return _png(render_rgb(*args, **kw))


def _diff_rgb(a, b, view, i, vmax, spacing, box, alpha, oa=(0, 0, 0), ob=(0, 0, 0)):
    """Différence de SUV (A − B) : bleu = baisse, rouge = hausse, sur le CT de A ; échelle ±vmax/2, opaque dès 10 % de l'échelle.
    oa, ob : décalages manuels de chaque séance."""
    diff = _up(_slice_off(a[0], view, i, oa).astype(np.float32) - _slice_off(b[0], view, i, ob), view, spacing, box)
    dmax = vmax / 2
    col = _lut("bwr")[((np.clip(diff / dmax, -1, 1) * 0.5 + 0.5) * 255).astype(np.uint8)].astype(np.float32)
    if a[1] is not None and view != "mip":
        base = np.clip((_up(_slice_off(a[1], view, i, oa, -1000), view, spacing, box) + 160) * (255 / 400), 0, 255)[..., None] * np.ones(3, np.float32)
    else:
        base = np.zeros(col.shape, np.float32)
    w = (alpha * np.clip(np.abs(diff) / (0.1 * dmax), 0, 1))[..., None]
    return np.clip(base * (1 - w) + col * w, 0, 255).astype(np.uint8)


def render_study(a, b, cmp, split, view, i, vmax, cmap, spacing, mode, alpha, box, labels=None, mip=None,
                 oa=(0, 0, 0), ob=(0, 0, 0)):
    """Rendu avec comparaison. a, b : (suv, ctv) des deux séances (b déjà recalée sur la grille de a, ou None).
    cmp : "a", "b", "swipe" (A à gauche du curseur `split` de 0 à 1, B à droite) ou "diff" (SUV A − B).
    oa, ob : décalages manuels (dz, dy, dx) en voxels de chaque séance, indépendants.
    Les organes et le MIP tourné ne concernent que la séance A."""
    if b is None or cmp == "a" or mip is not None:
        return _png(render_rgb(a[0], view, i, vmax, cmap, spacing, a[1], mode, alpha, box, labels, mip, oa))
    if cmp == "b":
        return _png(render_rgb(b[0], view, i, vmax, cmap, spacing, b[1], mode, alpha, box, None, None, ob))
    if cmp == "diff":
        return _png(_diff_rgb(a, b, view, i, vmax, spacing, box, alpha, oa, ob))
    ra = render_rgb(a[0], view, i, vmax, cmap, spacing, a[1], mode, alpha, box, labels, None, oa)
    out = render_rgb(b[0], view, i, vmax, cmap, spacing, b[1], mode, alpha, box, None, None, ob)
    x = int(round(split * out.shape[1]))
    out[:, :x] = ra[:, :x]
    if 0 < x < out.shape[1]:
        out[:, max(0, x - 1):x + 1] = 255  # trait de séparation
    return _png(out)


# --- dossiers récents et sélecteur natif ---
def recent():
    try:
        return json.loads(CONFIG.read_text(encoding="utf-8")).get("recent", [])
    except (OSError, ValueError):
        return []


def remember(folder):
    r = [folder] + [f for f in recent() if f != folder]
    try:
        CONFIG.write_text(json.dumps({"recent": r[:8]}), encoding="utf-8")
    except OSError:
        pass


def _ask_folder():
    import tkinter as tk
    from tkinter import filedialog
    r = tk.Tk(); r.withdraw(); r.attributes("-topmost", True)
    print(filedialog.askdirectory(title="Dossier DICOM"))


def browse():
    """Sélecteur de dossier natif (tkinter, dans un sous-processus car il veut le thread principal) ; None si annulé.
    Le sous-processus relance ce programme avec --browse (application packagée : sys.executable est l'exécutable lui-même)."""
    cmd = [sys.executable] + ([] if getattr(sys, "frozen", False) else [__file__]) + ["--browse"]
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=300).stdout.strip() or None
    except (OSError, subprocess.TimeoutExpired):
        return None


class State:
    """Volume courant, remplaçable à chaud : `vol` est un tuple (suv, spacing, ctv, patient, folder) échangé d'un bloc,
    donc lisible sans verrou. Le chargement d'un dossier se fait en arrière-plan ; en cas d'erreur l'ancien volume reste."""

    def __init__(self):
        self.vol = None
        self.organs = None  # (labels, infos) des organes du volume courant, calculés en arrière-plan
        self.mips = None    # {angle: image} du MIP tourné, remplie en arrière-plan
        self.cts = None     # séries CT du dossier : {series, origin, choices, uid}
        self.cmp = None     # séance B recalée sur la grille de A : {suv, ctv, patient, folder, shift, score}
        self.cmp_status, self.cmp_message = "idle", ""
        self._ctcache = {}  # {uid: CT rééchantillonné} : un changement de série déjà chargée est instantané
        self._orgcache = {} # {uid: (labels, infos)} : idem pour les organes
        self.status, self.message = "idle", ""
        self.progress = self.cmp_progress = None  # avancement 0..1 du chargement en cours (None = inconnu)
        self._lock = threading.Lock()

    def set(self, suv, spacing, ctv, patient=None, folder="", sul=None):
        self.sul = sul or sul_info(None, None, None)  # facteur SUV -> SUL ; repli SUV poids sans taille ni sexe
        self._ref = (None, None)  # (clé, référence foie PERCIST) calculée à la demande
        self.organs = self.mips = self.cts = self.cmp = None
        self.cmp_status, self.cmp_message = "idle", ""
        self._ctcache, self._orgcache = {}, {}
        self.vol = (suv, spacing, ctv, patient or {}, folder)
        self.status, self.message = "ready", ""
        threading.Thread(target=self._mips, args=(self.vol,), daemon=True).start()

    def liver_ref(self, suv, labels, spacing):
        """Référence foie PERCIST (`liver_reference`), recalculée seulement si le volume, les organes ou le facteur changent."""
        key = (id(labels), id(suv), self.sul["factor"])
        if self._ref[0] != key:
            self._ref = (key, liver_reference(suv, labels, spacing, self.sul["factor"]))
        return self._ref[1]

    def _organs(self, vol, uid):
        cache = self._orgcache  # celui de ce dossier : si on change de dossier pendant le calcul, on n'y touche plus
        try:
            res = segment_organs(vol[0], vol[2], vol[1])
        except Exception:  # une segmentation ratée ne doit pas gêner la visionneuse
            return
        cache[uid] = res
        if self.vol is vol:  # le volume n'a pas changé pendant le calcul
            self.organs = res

    def _mips(self, vol):
        """MIP de face tourné autour de l'axe vertical, tous les 10° (grille xy réduite ×2 : environ 1 s par image)."""
        try:
            sub = vol[0][:, ::2, ::2]
            frames = {}
            if self.vol is vol:
                self.mips = frames  # remplie au fil du calcul : les premiers angles sont déjà utilisables
            for a in ANGLES:
                if a == 0:
                    continue  # l'angle 0 est servi en pleine résolution, sans cache
                if self.vol is not vol:
                    return
                frames[a] = ndimage.rotate(sub, a, axes=(1, 2), reshape=False, order=1).max(axis=1)[::-1]
        except Exception:
            pass

    def mip_frame(self, i):
        """Image du MIP à l'angle i (arrondi à 10°) ; None pour 0° (pleine résolution) ou si pas encore calculée."""
        a = _angle(i)
        return None if a == 0 else (self.mips or {}).get(a)

    def open_ct(self, uid):
        """Recharge seulement le CT depuis une autre série ; None si la série est inconnue, False si un chargement est en cours."""
        cts = self.cts
        if not cts or uid not in cts["series"]:
            return None
        ctv = self._ctcache.get(uid)
        if ctv is not None:  # série déjà en mémoire : bascule instantanée
            self._swap_ct(uid, ctv)
            return True
        with self._lock:
            if self.status == "loading":
                return False
            self.status, self.message = "loading", "Chargement du CT…"
        threading.Thread(target=self._reload_ct, args=(uid,), daemon=True).start()
        return True

    def _swap_ct(self, uid, ctv):
        """Remplace le CT affiché ; les organes viennent du cache s'ils ont déjà été calculés pour cette série."""
        vol = self.vol
        self.vol = (vol[0], vol[1], ctv, vol[3], vol[4])
        self.cts = {**self.cts, "uid": uid}
        self.organs = self._orgcache.get(uid)
        if self.organs is None:
            threading.Thread(target=self._organs, args=(self.vol, uid), daemon=True).start()

    def _reload_ct(self, uid):
        vol, cts, cache = self.vol, self.cts, self._ctcache
        try:
            ctv = load_ct(vol[4], vol[0].shape, vol[1], cts["origin"], uid, cts["series"], progress=self._set_progress)
        except SystemExit as e:
            self.status, self.message = "error", str(e)
        except Exception as e:
            self.status, self.message = "error", f"{type(e).__name__} : {e}"
        else:
            cache[uid] = ctv
            self._swap_ct(uid, ctv)
            self.status, self.message = "ready", ""
        self.progress = None

    def _preload_cts(self, vol):
        """Charge en arrière-plan les autres séries CT du dossier : les changements de couche sont ensuite instantanés."""
        cache = self._ctcache
        for c in (self.cts or {}).get("choices", []):
            cts = self.cts
            if self.vol is None or self.vol[0] is not vol[0] or not cts:
                return  # autre dossier entre-temps
            if c["uid"] not in cache:
                try:
                    cache[c["uid"]] = load_ct(vol[4], vol[0].shape, vol[1], cts["origin"], c["uid"], cts["series"])
                except Exception:
                    pass

    def open_compare(self, folder):
        """Charge une 2e séance du même patient et la recale sur la grille de la séance courante (en arrière-plan).
        False s'il n'y a pas de volume ou si un chargement de séance B est déjà en cours."""
        with self._lock:
            if self.vol is None or self.cmp_status == "loading":
                return False
            self.cmp_status, self.cmp_message = "loading", "Chargement de la séance B…"
            self.cmp_progress = None
        threading.Thread(target=self._load_compare, args=(folder, self.vol), daemon=True).start()
        return True

    def _load_compare(self, folder, vol):
        try:
            prog = lambda f: setattr(self, "cmp_progress", f)
            suv, spacing, origin = load_suv(folder, _scaled(prog, 0, 0.35))
            if any(abs(a - b) > 0.1 * a for a, b in zip(vol[1], spacing)):
                raise SystemExit("Grilles incompatibles : espacement " + " × ".join(f"{s:.2f}" for s in spacing) +
                                 " mm contre " + " × ".join(f"{s:.2f}" for s in vol[1]) + " mm")
            series = ct_series(folder, _scaled(prog, 0.35, 0.45))
            ctb = load_ct(folder, suv.shape, spacing, origin, default_ct_uid(series, origin) if series else None, series,
                          _scaled(prog, 0.45, 0.8))
            patient = read_patient_info(folder)
            if vol[2] is None or ctb is None:
                shift, score = (0, 0, 0), None  # pas de CT pour recaler : on garde l'alignement des grilles, sans recalage
            else:
                self.cmp_message = "Recalage des deux séances…"
                self.cmp_progress = 0.85
                shift, score = estimate_shift(vol[2], ctb)
            shift = tuple(int(s) for s in shift)
            aligned = (shift_volume(suv, shift, vol[0].shape, 0.0),
                       shift_volume(ctb, shift, vol[0].shape, -1000) if ctb is not None else None)
        except SystemExit as e:
            self.cmp_status, self.cmp_message = "error", str(e)
        except Exception as e:
            self.cmp_status, self.cmp_message = "error", f"{type(e).__name__} : {e}"
        else:
            if self.vol is vol:  # la séance A n'a pas changé pendant le calcul
                self.cmp = {"suv": aligned[0], "ctv": aligned[1], "patient": patient, "folder": folder, "shift": shift, "score": score}
                self.cmp_status, self.cmp_message = "ready", ""
        self.cmp_progress = None

    def open(self, folder):
        """Lance le chargement ; False si un chargement est déjà en cours."""
        with self._lock:
            if self.status == "loading":
                return False
            self.status, self.message = "loading", f"Chargement de {folder}…"
            self.progress = None  # sinon la barre repartirait de la valeur finale du chargement précédent (100 %)
        threading.Thread(target=self._load, args=(folder,), daemon=True).start()
        return True

    def _set_progress(self, f):
        self.progress = f

    def _load(self, folder):
        try:
            suv, spacing, origin = load_suv(folder, _scaled(self._set_progress, 0, 0.45))
            series = ct_series(folder, _scaled(self._set_progress, 0.45, 0.6))
            uid = default_ct_uid(series, origin) if series else None
            ctv = load_ct(folder, suv.shape, spacing, origin, uid, series, _scaled(self._set_progress, 0.6, 1))
            patient, sul = read_patient_info(folder), read_sul(folder)
        except SystemExit as e:  # `load_suv` signale « aucune série PET » par SystemExit : ne pas tuer le serveur
            self.status, self.message = "error", str(e)
        except Exception as e:
            self.status, self.message = "error", f"{type(e).__name__} : {e}"
        else:
            self.set(suv, spacing, ctv, patient, folder, sul)
            if series:
                self.cts = {"series": series, "origin": origin, "choices": ct_choices(series, origin), "uid": uid}
                self._ctcache[uid] = ctv
                threading.Thread(target=self._preload_cts, args=(self.vol,), daemon=True).start()
            remember(folder)
            if ctv is not None:
                threading.Thread(target=self._organs, args=(self.vol, uid), daemon=True).start()
        self.progress = None


def make_handler(state):
    cache, cache_lock = OrderedDict(), threading.Lock()  # images déjà rendues (coupes revues, coupes voisines préchargées)

    def cached(key, make):
        with cache_lock:
            if key in cache:
                cache.move_to_end(key)
                return cache[key]
        png = make()
        with cache_lock:
            cache[key] = png
            while len(cache) > 300:
                cache.popitem(last=False)
        return png

    class H(BaseHTTPRequestHandler):
        def _send(self, body, ct, code=200, extra=None):
            self.send_response(code)
            self.send_header("Content-Type", ct)
            for k, v in (extra or {}).items():
                self.send_header(k, v)
            self.send_header("Cache-Control", "no-store")  # données patient : pas de cache navigateur
            self.end_headers()
            self.wfile.write(body)

        def _mip_source(self, q, spacing):
            """(image du MIP tourné, espacement à utiliser) ; (None, spacing) hors MIP ou à 0° ; (False, …) si pas encore calculée."""
            if q["view"] != "mip":
                return None, spacing
            frame = state.mip_frame(int(q["i"]))
            if frame is not None:
                return frame, (spacing[0], spacing[1] * 2, spacing[2] * 2)  # calculé sur une grille xy réduite ×2
            return (None if _angle(int(q["i"])) == 0 else False), spacing

        @staticmethod
        def _off(q, key, lim=400):
            """Décalage manuel (dz, dy, dx) en voxels, passé sous la forme « z,y,x » ; ValueError si invalide (réponse 400)."""
            v = tuple(int(c) for c in q.get(key, "0,0,0").split(","))
            if len(v) != 3 or any(abs(c) > lim for c in v):
                raise ValueError("décalage invalide")
            return v

        def _json(self, obj, code=200):
            self._send(json.dumps(obj).encode(), "application/json", code)

        def _allowed(self, mutating):
            # Hôte local uniquement (anti DNS-rebinding) ; les actions qui changent l'état refusent les requêtes inter-sites.
            if self.headers.get("Host", "").rsplit(":", 1)[0] not in ("127.0.0.1", "localhost"):
                return False
            return not mutating or self.headers.get("Sec-Fetch-Site", "same-origin") in ("same-origin", "none")

        def do_GET(self):
            u = urlparse(self.path)
            q = {k: v[0] for k, v in parse_qs(u.query).items()}
            if not self._allowed(u.path in ("/open", "/browse", "/ctseries", "/compare", "/autoalign")):
                return self.send_error(403)
            vol = state.vol  # lu une seule fois : cohérent même si un chargement termine entre-temps
            try:
                if u.path == "/":
                    # relu à chaque requête : pas de redémarrage pour modifier la page
                    return self._send(Path(os.environ.get("PETSCAN_INDEX") or Path(__file__).with_name("index.html")).read_bytes(), "text/html; charset=utf-8")
                if u.path == "/status":
                    return self._json({"status": state.status, "message": state.message, "loaded": vol is not None, "progress": state.progress})
                if u.path == "/open":
                    folder = q["folder"]
                    if not Path(folder).is_dir():
                        return self._json({"error": f"Dossier introuvable : {folder}"}, 400)
                    return self._json({"started": state.open(folder)})
                if u.path == "/browse":
                    return self._json({"folder": browse() if BROWSE else None})
                if u.path == "/ctseries":
                    started = state.open_ct(q["uid"])
                    return self._json({"error": "Série CT inconnue"}, 400) if started is None else self._json({"started": started})
                if u.path == "/compare":  # charge une 2e séance du même patient et la recale sur la séance courante
                    folder = q["folder"]
                    if not Path(folder).is_dir():
                        return self._json({"error": f"Dossier introuvable : {folder}"}, 400)
                    return self._json({"started": state.open_compare(folder)})
                if u.path == "/cstatus":
                    c, pa = state.cmp, (vol[3] if vol else {})
                    return self._json({"status": state.cmp_status, "message": state.cmp_message, "ready": c is not None, "progress": state.cmp_progress,
                                       "date_a": pa.get("Date de l'examen"), "date_b": c["patient"].get("Date de l'examen") if c else None,
                                       "weight_a": pa.get("Poids"), "weight_b": c["patient"].get("Poids") if c else None,
                                       "shift_mm": [round(s * v, 1) for s, v in zip(c["shift"], vol[1])] if c and vol else None,
                                       "score": c["score"] if c else None})
                if u.path == "/mips":  # avancement du calcul du MIP rotatif (l'angle 0 est toujours disponible)
                    return self._json({"ready": 0 if state.mips is None else len(state.mips) + 1, "total": len(ANGLES)})
                if u.path == "/bar" and q.get("cmap") in CMAPS:
                    return self._send(bar(q["cmap"]), "image/png")
                if u.path == "/info" and vol is None:
                    return self._json({"loaded": False, "browse": BROWSE, "recent": recent()})
                if vol is None:
                    return self._json({"error": "Aucun volume chargé"}, 503)
                suv, spacing, ctv, patient, folder = vol
                if u.path == "/info":
                    return self._json({"loaded": True, "shape": suv.shape, "spacing": spacing, "suv_max": float(suv.max()),
                                       "ct": ctv is not None, "folder": folder, "browse": BROWSE, "recent": recent(),
                                       "ct_series": (state.cts or {}).get("choices", []), "ct_uid": (state.cts or {}).get("uid")})
                if u.path == "/autoalign":
                    # Aligne B sur A à partir des décalages manuels actuels : renvoie le nouveau décalage de B (convention du serveur)
                    # et le résidu e. scope=local : zone de ±70 mm autour des coupes affichées (c = z,y,x), au voxel près ;
                    # scope=global : volume entier, au pas de 3 voxels.
                    cb = state.cmp
                    if cb is None or ctv is None or cb["ctv"] is None:
                        return self._json({"error": "L'alignement automatique demande le CT des deux séances."})
                    da, db = self._off(q, "oa"), self._off(q, "ob")
                    if q.get("scope") == "local":
                        c = tuple(min(max(v, 0), n - 1) for v, n in zip(self._off(q, "c", 100_000), suv.shape))  # indices de coupes, ramenés dans le volume
                        e, score = refine_shift(ctv, cb["ctv"], da, db, c, tuple(max(8, round(70 / s)) for s in spacing), step=1, margin=14)
                    else:
                        e, score = refine_shift(ctv, cb["ctv"], da, db, step=3, margin=10)
                    return self._json({"ob": [d + x for d, x in zip(db, e)], "e": list(e), "score": score})
                if u.path == "/patient":
                    return self._json(patient)
                if u.path == "/lesions":
                    org = state.organs  # organes pas encore calculés : pas de champs « organe » ni « physio », l'interface relance /lesions
                    if q["thr"] == "percist":  # seuil de la référence hépatique (SUL) ramené en SUV
                        ref = org is not None and state.liver_ref(suv, org[0], spacing)
                        if not ref:
                            return self._json({"error": "Seuil PERCIST indisponible : pas de référence foie"})
                        q["thr"] = ref["seuil_suv"]
                    les = analyze(suv, float(np.prod(spacing)) / 1000, float(q["thr"]), spacing, state.sul["factor"],
                                  org[0] if org is not None else None)
                    if org is not None:
                        for lesion in les:
                            lesion["organe"] = _organ_at_peak(org[0], spacing, lesion["peak"])
                    rel =tuple(b - a for a, b in zip(self._off(q, "oa"), self._off(q, "ob")))  # B relativement à A
                    return self._json(_add_b(les, state.cmp["suv"], spacing, rel=rel) if state.cmp else les)
                if u.path == "/percist":  # unité (SUL ou repli SUV) et référence hépatique ; les lésions restent dans /lesions
                    org, si = state.organs, state.sul
                    ref, msg = None, ("Référence foie indisponible : pas de CT" if ctv is None
                                      else "Référence foie : en attente des organes…" if org is None else "")
                    if org is not None:
                        try:
                            ref = state.liver_ref(suv, org[0], spacing)
                        except Exception as e:
                            return self._json({"error": f"Référence foie : {type(e).__name__} : {e}"})
                        msg = ref["msg"] if ref else "Référence foie indisponible : foie introuvable ou trop petit"
                    pool = blood_pool_reference(suv, org[0], spacing, si["factor"]) if org is not None else None
                    return self._json({**si, "ref": ref, "ref_msg": msg, "pool": pool, "computing": ctv is not None and org is None})
                if u.path == "/img":
                    cmap, mode = q.get("cmap", "jet"), q.get("mode", "pet")
                    z = float(q.get("z", 1))
                    if cmap not in CMAPS or mode not in ("pet", "fusion", "ct") or not 1 <= z <= 32:
                        return self.send_error(400)
                    w = 1 / z
                    box = tuple(min(max(float(q.get(k, 0)), 0), 1 - w) for k in ("x0", "y0")) + (w,)
                    org, cmpd, cmpmode = state.organs, state.cmp, q.get("cmp", "a")
                    if cmpmode not in ("a", "b", "swipe", "diff"):
                        return self.send_error(400)
                    if cmpd is None:
                        cmpmode = "a"
                    if cmpmode != "a" and q["view"] == "mip":
                        frame, sp = None, spacing  # le MIP tourné n'existe que pour la séance A : en comparaison, MIP de face
                    else:
                        frame, sp = self._mip_source(q, spacing)
                    if frame is False:
                        return self._json({"error": "MIP en cours de calcul"}, 503)
                    split = round(min(max(float(q.get("split", 0.5)), 0), 1), 3)
                    b = (cmpd["suv"], cmpd["ctv"]) if cmpd else None
                    labels = org[0] if q.get("organs") == "1" and org else None
                    i = (_angle(int(q["i"])) if cmpmode == "a" else 0) if q["view"] == "mip" else int(q["i"])
                    oa, ob = self._off(q, "oa"), self._off(q, "ob")  # décalages manuels de A et de B
                    key = (q["view"], i, float(q["vmax"]), cmap, mode, float(q.get("alpha", 0.6)), box, id(suv), id(ctv),
                           id(labels), cmpmode, split if cmpmode == "swipe" else 0, id(b[0]) if b else 0, oa, ob)
                    return self._send(cached(key, lambda: render_study(
                        (suv, ctv), b, cmpmode, split, q["view"], int(q["i"]), float(q["vmax"]), cmap, sp, mode,
                        float(q.get("alpha", 0.6)), box, labels, frame, oa, ob)), "image/png")
                if u.path == "/organs":
                    org = state.organs
                    return self._json({"status": "unavailable" if ctv is None else "ready" if org else "computing",
                                       "organes": [dict(id=k, **v) for k, v in (org[1] if org else {}).items()]})
                if u.path == "/organ_at":
                    org = state.organs
                    if not org:
                        return self._json({})
                    lab = _slice_off(org[0], q["view"], int(q["i"]), self._off(q, "oa"))  # les organes suivent la séance A
                    oid = int(lab[min(lab.shape[0] - 1, max(0, int(float(q["fy"]) * lab.shape[0]))),
                                  min(lab.shape[1] - 1, max(0, int(float(q["fx"]) * lab.shape[1])))])
                    return self._json({"nom": ORGANS[oid][0]} if oid in ORGANS else {})
                if u.path == "/measure":
                    cmpd = state.cmp if q.get("cmp") == "b" else None  # mesure sur la séance B (déjà recalée) ou sur A
                    suv_m, ctv_m = (cmpd["suv"], cmpd["ctv"]) if cmpd else (suv, ctv)
                    frame, sp = (None, spacing) if (cmpd and q["view"] == "mip") else self._mip_source(q, spacing)
                    if frame is False:
                        return self._json({"error": "MIP en cours de calcul"})
                    return self._json(measure(suv_m, sp, q["view"], int(q["i"]), float(q["fx"]), float(q["fy"]),
                                              float(q.get("pct", 50)), ctv_m, q.get("source", "pet"), float(q.get("tol", 50)), frame,
                                              self._off(q, "ob" if cmpd else "oa")))
                return self.send_error(404)
            except (KeyError, ValueError):
                return self.send_error(400)

        def do_POST(self):
            if not self._allowed(True):
                return self.send_error(403)
            if urlparse(self.path).path != "/export.pdf":
                return self.send_error(404)
            vol = state.vol
            if vol is None:
                return self._json({"error": "Aucun volume chargé"}, 503)
            suv, spacing, ctv, patient, folder = vol
            try:
                n = int(self.headers.get("Content-Length", 0))
                if not 0 < n <= 2_000_000:
                    return self.send_error(413)
                req = json.loads(self.rfile.read(n))
                st = req["state"]
                if st["cmap"] not in CMAPS or st["mode"] not in ("pet", "fusion", "ct"):
                    raise ValueError("réglage inconnu")
                org = state.organs
                labels = org[0] if st.get("organs") and org else None  # surcouche d'organes comme à l'écran
                cmpd = state.cmp  # comparaison de séances : le PDF reproduit l'affichage (A, B, curseur ou différence)
                cm = st.get("cmp") if cmpd and st.get("cmp") in ("a", "b", "swipe", "diff") else "a"
                split, b = float(st.get("split", 0.5)), ((cmpd["suv"], cmpd["ctv"]) if cmpd else None)
                oa, ob = (tuple(int(c) for c in st.get(k, (0, 0, 0))) for k in ("oa", "ob"))  # décalages manuels affichés
                pdf = build_pdf(lambda view, i, vmax, cmap, mode, alpha, box:
                                render_study((suv, ctv), b, cm, split, view, i, vmax, cmap, spacing, mode, alpha, box, labels,
                                             None, oa, ob),
                                suv.shape, st, req.get("patient") or {})
            except (KeyError, ValueError, TypeError, AttributeError):
                return self.send_error(400)
            self._send(pdf, "application/pdf", extra={"Content-Disposition": 'attachment; filename="tep-ct.pdf"'})

        def log_message(self, fmt, *args):  # pas de journal des requêtes (les URL contiennent des chemins de dossiers)
            pass

    return H


def main():
    if "--browse" in sys.argv:
        return _ask_folder()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder", nargs="?", help="dossier de la série PET DICOM (sinon, à choisir dans la page)")
    ap.add_argument("--demo", action="store_true", help="volume synthétique")
    ap.add_argument("--port", type=int, default=8000)
    a = ap.parse_args()
    state = State()
    if a.demo:
        state.set(*demo())
    elif a.folder:
        state.open(a.folder)
    # 127.0.0.1 : données patient, jamais exposées au réseau
    print(f"http://127.0.0.1:{a.port}", flush=True)
    if getattr(sys, "frozen", False):  # application packagée (double-clic) : on ouvre la page
        webbrowser.open(f"http://127.0.0.1:{a.port}")
    ThreadingHTTPServer(("127.0.0.1", a.port), make_handler(state)).serve_forever()


if __name__ == "__main__":
    multiprocessing.freeze_support()  # application packagée : les processus de lecture DICOM relancent l'exécutable
    main()
