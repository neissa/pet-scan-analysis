"""Visionneuse web TEP : python viewer.py [dossier_dicom]  (ou --demo)."""
import argparse
import importlib.util
import io
import json
import subprocess
import sys
import threading
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

from export_pdf import build_pdf
from organs import ORGANS, segment_organs
from patient import read_patient_info
from petscan import (analyze, ct_choices, ct_series, default_ct_uid, load_ct, load_suv, measure_region,
                     measure_tolerance)

AXIS = {"axial": 0, "coronal": 1, "sagittal": 2}
CONFIG = Path.home() / ".pet-scan-analysis.json"  # dossiers récents, dans le dossier utilisateur (jamais dans le repo)
BROWSE = importlib.util.find_spec("tkinter") is not None  # sélecteur de dossier natif disponible ?


def demo():
    """Volume synthétique avec deux « lésions », pour essayer sans DICOM."""
    z, y, x = np.mgrid[:60, :80, :80]
    suv = 1 + 0.5 * np.random.default_rng(0).random((60, 80, 80))
    for (cz, cy, cx), r, v in [((30, 30, 30), 5, 9), ((20, 50, 55), 3, 5)]:
        suv[(z - cz) ** 2 + (y - cy) ** 2 + (x - cx) ** 2 <= r * r] = v
    ctv = np.where(((z - 30) / 28) ** 2 + ((y - 40) / 32) ** 2 + ((x - 40) / 32) ** 2 < 1, 40, -1000)
    return suv, (3.0, 3.0, 3.0), ctv.astype(np.int16)


CMAPS = {"hot", "jet", "gray_r"}


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


def measure(suv, spacing, view, i, fx, fy, pct, ctv=None, source="pet", tol=50, mip=None):
    """Mesure auto autour du point (fx, fy), fractions de l'image entière ; résultats en fractions.
    source "pet" : région >= pct % du max local ; source "ct" : région à ± tol HU du HU local.
    mip : image d'un MIP tourné (spacing alors celui de cette image) ; le MIP est toujours mesuré sur le TEP."""
    if view == "mip":
        source = "pet"
    if source == "ct" and ctv is None:
        return {"error": "pas de CT dans ce dossier"}
    img = mip if (view == "mip" and mip is not None) else _slice(ctv if source == "ct" else suv, view, i)
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


ANGLES = range(0, 360, 10)  # angles du MIP rotatif


def _angle(i):
    return round(i / 10) * 10 % 360


def render(suv, view, i, vmax, cmap="jet", spacing=(1, 1, 1), ctv=None, mode="pet", alpha=0.6, box=(0, 0, 1), labels=None,
           mip=None):
    """PNG d'une coupe (axial/coronal/sagittal) ou de la MIP coronale (toujours TEP seul). i est borné.
    mode : "pet", "ct" ou "fusion" (TEP en couleur sur CT ; opacité nulle à SUV 0, pleine dès 25 % de vmax).
    Fenêtre CT fixe -160..240 HU (tissus mous). labels : volume d'organes à superposer (sauf sur la MIP)."""
    # mip : image d'un MIP tourné déjà calculée (spacing = celui de cette image) ; sinon MIP de face au pleine résolution
    show_ct = ctv is not None and mode in ("fusion", "ct") and view != "mip"
    show_organs = labels is not None and view != "mip"
    # les interpolations TEP et CT sont indépendantes : en parallèle (scipy libère le verrou Python)
    pet_f = _POOL.submit(_up, mip if (view == "mip" and mip is not None) else _slice(suv, view, i), view, spacing, box)
    ct_f = _POOL.submit(_up, _slice(ctv, view, i), view, spacing, box) if show_ct else None
    pet = pet_f.result()
    idx = (np.clip(pet / vmax, 0, 1) * 255).astype(np.uint8)  # indice dans la palette
    if not (show_ct or show_organs):
        return _png(_lut(cmap)[idx])
    if show_ct:
        img = np.clip((ct_f.result() + 160) * (255 / 400), 0, 255)[..., None] * np.ones(3, np.float32)
        if mode == "fusion":
            a = (alpha * np.clip(pet / (0.25 * vmax), 0, 1))[..., None]
            img = img * (1 - a) + _lut(cmap)[idx] * a
    else:
        img = _lut(cmap)[idx].astype(np.float32)
    if show_organs:
        img = _organs_overlay(img, np.rint(_up(_slice(labels, view, i), view, spacing, box, order=0)).astype(np.uint8))
    return _png(np.clip(img, 0, 255).astype(np.uint8))


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


def browse():
    """Sélecteur de dossier natif (tkinter, dans un sous-processus car il veut le thread principal) ; None si annulé."""
    code = ("import tkinter as tk\nfrom tkinter import filedialog\n"
            "r = tk.Tk(); r.withdraw(); r.attributes('-topmost', True)\n"
            "print(filedialog.askdirectory(title='Dossier DICOM'))")
    try:
        return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300).stdout.strip() or None
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
        self._ctcache = {}  # {uid: CT rééchantillonné} : un changement de série déjà chargée est instantané
        self._orgcache = {} # {uid: (labels, infos)} : idem pour les organes
        self.status, self.message = "idle", ""
        self._lock = threading.Lock()

    def set(self, suv, spacing, ctv, patient=None, folder=""):
        self.organs = self.mips = self.cts = None
        self._ctcache, self._orgcache = {}, {}
        self.vol = (suv, spacing, ctv, patient or {}, folder)
        self.status, self.message = "ready", ""
        threading.Thread(target=self._mips, args=(self.vol,), daemon=True).start()

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
            ctv = load_ct(vol[4], vol[0].shape, vol[1], cts["origin"], uid, cts["series"])
        except SystemExit as e:
            self.status, self.message = "error", str(e)
        except Exception as e:
            self.status, self.message = "error", f"{type(e).__name__} : {e}"
        else:
            cache[uid] = ctv
            self._swap_ct(uid, ctv)
            self.status, self.message = "ready", ""

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

    def open(self, folder):
        """Lance le chargement ; False si un chargement est déjà en cours."""
        with self._lock:
            if self.status == "loading":
                return False
            self.status, self.message = "loading", f"Chargement de {folder}…"
        threading.Thread(target=self._load, args=(folder,), daemon=True).start()
        return True

    def _load(self, folder):
        try:
            suv, spacing, origin = load_suv(folder)
            series = ct_series(folder)
            uid = default_ct_uid(series, origin) if series else None
            ctv = load_ct(folder, suv.shape, spacing, origin, uid, series)
            patient = read_patient_info(folder)
        except SystemExit as e:  # `load_suv` signale « aucune série PET » par SystemExit : ne pas tuer le serveur
            self.status, self.message = "error", str(e)
        except Exception as e:
            self.status, self.message = "error", f"{type(e).__name__} : {e}"
        else:
            self.set(suv, spacing, ctv, patient, folder)
            if series:
                self.cts = {"series": series, "origin": origin, "choices": ct_choices(series, origin), "uid": uid}
                self._ctcache[uid] = ctv
                threading.Thread(target=self._preload_cts, args=(self.vol,), daemon=True).start()
            remember(folder)
            if ctv is not None:
                threading.Thread(target=self._organs, args=(self.vol, uid), daemon=True).start()


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
            if not self._allowed(u.path in ("/open", "/browse", "/ctseries")):
                return self.send_error(403)
            vol = state.vol  # lu une seule fois : cohérent même si un chargement termine entre-temps
            try:
                if u.path == "/":
                    # relu à chaque requête : pas de redémarrage pour modifier la page
                    return self._send(Path(__file__).with_name("index.html").read_bytes(), "text/html; charset=utf-8")
                if u.path == "/status":
                    return self._json({"status": state.status, "message": state.message, "loaded": vol is not None})
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
                if u.path == "/patient":
                    return self._json(patient)
                if u.path == "/lesions":
                    return self._json(analyze(suv, float(np.prod(spacing)) / 1000, float(q["thr"])))
                if u.path == "/img":
                    cmap, mode = q.get("cmap", "jet"), q.get("mode", "pet")
                    z = float(q.get("z", 1))
                    if cmap not in CMAPS or mode not in ("pet", "fusion", "ct") or not 1 <= z <= 32:
                        return self.send_error(400)
                    w = 1 / z
                    box = tuple(min(max(float(q.get(k, 0)), 0), 1 - w) for k in ("x0", "y0")) + (w,)
                    org = state.organs
                    frame, sp = self._mip_source(q, spacing)
                    if frame is False:
                        return self._json({"error": "MIP en cours de calcul"}, 503)
                    labels = org[0] if q.get("organs") == "1" and org else None
                    i = _angle(int(q["i"])) if q["view"] == "mip" else int(q["i"])
                    key = (q["view"], i, float(q["vmax"]), cmap, mode, float(q.get("alpha", 0.6)), box,
                           id(suv), id(ctv), id(labels))  # un autre volume, CT ou jeu d'organes = une autre clé
                    return self._send(cached(key, lambda: render(suv, q["view"], int(q["i"]), float(q["vmax"]), cmap, sp,
                                                                 ctv, mode, float(q.get("alpha", 0.6)), box, labels, frame)),
                                      "image/png")
                if u.path == "/organs":
                    org = state.organs
                    return self._json({"status": "unavailable" if ctv is None else "ready" if org else "computing",
                                       "organes": [dict(id=k, **v) for k, v in (org[1] if org else {}).items()]})
                if u.path == "/organ_at":
                    org = state.organs
                    if not org:
                        return self._json({})
                    lab = _slice(org[0], q["view"], int(q["i"]))
                    oid = int(lab[min(lab.shape[0] - 1, max(0, int(float(q["fy"]) * lab.shape[0]))),
                                  min(lab.shape[1] - 1, max(0, int(float(q["fx"]) * lab.shape[1])))])
                    return self._json({"nom": ORGANS[oid][0]} if oid in ORGANS else {})
                if u.path == "/measure":
                    frame, sp = self._mip_source(q, spacing)
                    if frame is False:
                        return self._json({"error": "MIP en cours de calcul"})
                    return self._json(measure(suv, sp, q["view"], int(q["i"]), float(q["fx"]), float(q["fy"]),
                                              float(q.get("pct", 50)), ctv, q.get("source", "pet"), float(q.get("tol", 50)), frame))
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
                pdf = build_pdf(lambda view, i, vmax, cmap, mode, alpha, box:
                                render(suv, view, i, vmax, cmap, spacing, ctv, mode, alpha, box, labels),
                                suv.shape, st, req.get("patient") or {})
            except (KeyError, ValueError, TypeError, AttributeError):
                return self.send_error(400)
            self._send(pdf, "application/pdf", extra={"Content-Disposition": 'attachment; filename="tep-ct.pdf"'})

        def log_message(self, fmt, *args):  # pas de journal des requêtes (les URL contiennent des chemins de dossiers)
            pass

    return H


def main():
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
    ThreadingHTTPServer(("127.0.0.1", a.port), make_handler(state)).serve_forever()


if __name__ == "__main__":
    main()
