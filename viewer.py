"""Visionneuse web TEP : python viewer.py <dossier_dicom>  (ou --demo)."""
import argparse
import io
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

from petscan import analyze, load_ct, load_suv, measure_region

AXIS = {"axial": 0, "coronal": 1, "sagittal": 2}


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


def _up(img, view, spacing, box=(0, 0, 1)):
    """Rééchantillonnage cubique à l'espacement réel des voxels (spacing en mm, z,y,x).
    box = (x0, y0, w) : on ne garde que ce carré (fractions de l'image entière) et il remplit
    toute l'image de sortie, donc un zoom augmente vraiment la résolution."""
    r, c = ROWCOL[view]
    m = min(spacing)
    nr, nc = img.shape
    out = (round(UP * nr * spacing[r] / m), round(UP * nc * spacing[c] / m))
    x0, y0, w = box
    sr, sc = w * nr / out[0], w * nc / out[1]  # pas en pixels source par pixel de sortie
    return ndimage.affine_transform(img.astype(np.float32), [sr, sc], output_shape=out, order=3, mode="nearest",
                                    offset=[y0 * nr - 0.5 + sr / 2, x0 * nc - 0.5 + sc / 2])


def measure(suv, spacing, view, i, fx, fy, pct):
    """Mesure auto (TEP) autour du point (fx, fy), fractions de l'image entière ; résultats en fractions."""
    img = _slice(suv, view, i)
    nr, nc = img.shape
    rr, cc = ROWCOL[view]
    res = measure_region(img, spacing[rr], spacing[cc], min(nr - 1, max(0, int(fy * nr))),
                         min(nc - 1, max(0, int(fx * nc))), pct)
    if res is None:
        return {"error": "aucune activité autour du point"}
    frac = lambda pts: [[(c + 0.5) / nc, (r + 0.5) / nr] for r, c in pts]  # (ligne, col) -> [x, y]
    res["axes"], res["perp"], res["boundary"] = frac(res["axes"]), frac(res["perp"]), frac(res["boundary"])
    return res


def render(suv, view, i, vmax, cmap="jet", spacing=(1, 1, 1), ctv=None, mode="pet", alpha=0.6, box=(0, 0, 1)):
    """PNG d'une coupe (axial/coronal/sagittal) ou de la MIP coronale (toujours TEP seul). i est borné.
    mode : "pet", "ct" ou "fusion" (TEP en couleur sur CT ; opacité nulle à SUV 0, pleine dès 25 % de vmax).
    Fenêtre CT fixe -160..240 HU (tissus mous)."""
    pet = _up(_slice(suv, view, i), view, spacing, box)
    buf = io.BytesIO()
    if ctv is not None and mode in ("fusion", "ct") and view != "mip":
        img = np.clip((_up(_slice(ctv, view, i), view, spacing, box) + 160) / 400, 0, 1)[..., None] * np.ones(3)
        if mode == "fusion":
            col = plt.get_cmap(cmap)(np.clip(pet / vmax, 0, 1))[..., :3]
            a = (alpha * np.clip(pet / (0.25 * vmax), 0, 1))[..., None]
            img = img * (1 - a) + col * a
        plt.imsave(buf, np.clip(img, 0, 1), format="png")
    else:
        plt.imsave(buf, pet, cmap=cmap, vmin=0, vmax=vmax, format="png")
    return buf.getvalue()


def make_handler(suv, spacing, ctv=None):
    voxel_ml = float(np.prod(spacing)) / 1000

    class H(BaseHTTPRequestHandler):
        def do_GET(self):
            u = urlparse(self.path)
            q = {k: v[0] for k, v in parse_qs(u.query).items()}
            try:
                if u.path == "/":
                    body, ct = Path(__file__).with_name("index.html").read_bytes(), "text/html; charset=utf-8"  # relu à chaque requête : pas de redémarrage pour modifier la page
                elif u.path == "/info":
                    body, ct = json.dumps({"shape": suv.shape, "spacing": spacing, "suv_max": float(suv.max()), "ct": ctv is not None}).encode(), "application/json"
                elif u.path == "/lesions":
                    body, ct = json.dumps(analyze(suv, voxel_ml, float(q["thr"]))).encode(), "application/json"
                elif u.path == "/img":
                    cmap, mode = q.get("cmap", "jet"), q.get("mode", "pet")
                    z = float(q.get("z", 1))
                    if cmap not in CMAPS or mode not in ("pet", "fusion", "ct") or not 1 <= z <= 32:
                        return self.send_error(400)
                    w = 1 / z
                    box = tuple(min(max(float(q.get(k, 0)), 0), 1 - w) for k in ("x0", "y0")) + (w,)
                    body, ct = render(suv, q["view"], int(q["i"]), float(q["vmax"]), cmap, spacing,
                                      ctv, mode, float(q.get("alpha", 0.6)), box), "image/png"
                elif u.path == "/measure":
                    body, ct = json.dumps(measure(suv, spacing, q["view"], int(q["i"]), float(q["fx"]),
                                                  float(q["fy"]), float(q["pct"]))).encode(), "application/json"
                elif u.path == "/bar" and q.get("cmap") in CMAPS:
                    body, ct = bar(q["cmap"]), "image/png"
                else:
                    return self.send_error(404)
            except (KeyError, ValueError):
                return self.send_error(400)
            self.send_response(200)
            self.send_header("Content-Type", ct)
            self.end_headers()
            self.wfile.write(body)

    return H


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder", nargs="?", help="dossier de la série PET DICOM")
    ap.add_argument("--demo", action="store_true", help="volume synthétique")
    ap.add_argument("--port", type=int, default=8000)
    a = ap.parse_args()
    if not (a.folder or a.demo):
        ap.error("donnez un dossier DICOM ou --demo")
    if a.demo:
        suv, spacing, ctv = demo()
    else:
        suv, spacing, origin = load_suv(a.folder)
        ctv = load_ct(a.folder, suv.shape, spacing, origin)
    # 127.0.0.1 : données patient, jamais exposées au réseau
    print(f"http://127.0.0.1:{a.port}")
    ThreadingHTTPServer(("127.0.0.1", a.port), make_handler(suv, spacing, ctv)).serve_forever()


if __name__ == "__main__":
    main()
