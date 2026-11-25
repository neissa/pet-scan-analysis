"""Export PDF (A4) de l'état affiché : les 3 coupes avec leurs mesures, les infos patient, la liste des mesures.
Le PDF est construit en mémoire et renvoyé au navigateur : rien n'est écrit sur disque."""
import io
import math
import textwrap
from datetime import datetime

import numpy as np
from matplotlib import image as mpimg
from matplotlib.figure import Figure

VIEWS = ("axial", "coronal", "sagittal")
W, H, M = 8.27, 11.69, 0.6  # page A4 et marge, en pouces


def _check(state):
    """Valide l'état envoyé par la page ; ValueError si invalide."""
    views = state["views"]
    if [v["view"] for v in views] != list(VIEWS):
        raise ValueError("vues attendues : axial, coronal, sagittal")
    for v in views:
        int(v["i"])
        if not 1 <= float(v["z"]) <= 32:
            raise ValueError("zoom hors limites")
    float(state["vmax"]), float(state["alpha"])
    str(state["cmap"]), str(state["mode"])
    if len(state.get("measures", [])) > 200:
        raise ValueError("trop de mesures")


def _box(v):
    w = 1 / float(v["z"])
    return (min(max(float(v["x0"]), 0), 1 - w), min(max(float(v["y0"]), 0), 1 - w), w)


def _draw_items(ax, items, box, size):
    """Dessine les mesures (coordonnées en fractions de l'image entière) sur un axe de size=(largeur, hauteur) pouces."""
    x0, y0, w = box
    P = lambda p: ((float(p[0]) - x0) / w, 1 - (float(p[1]) - y0) / w)  # -> coordonnées de l'axe (y vers le haut)
    for it in items:
        t = it.get("t")
        if t == "dots" and it.get("pts"):
            xs, ys = zip(*[P(p) for p in it["pts"]])
            ax.plot(xs, ys, ".", ms=1.2, color="yellow", transform=ax.transAxes)
        elif t == "line":
            a, b = P(it["a"]), P(it["b"])
            col = it["color"] if str(it.get("color", "")).startswith("#") else "white"
            ax.plot([a[0], b[0]], [a[1], b[1]], "-", lw=0.9, color=col, transform=ax.transAxes)
            dx, dy = (b[0] - a[0]) * size[0], (b[1] - a[1]) * size[1]  # en pouces, pour des barres vraiment perpendiculaires
            d = math.hypot(dx, dy) or 1
            nx, ny = -dy / d * 0.035 / size[0], dx / d * 0.035 / size[1]
            if it.get("bars"):
                for p in (a, b):
                    ax.plot([p[0] - nx, p[0] + nx], [p[1] - ny, p[1] + ny], "-", lw=0.9, color=col, transform=ax.transAxes)
            right = b[0] < 0.55  # étiquette vers l'intérieur de l'image : jamais coupée par le bord
            ax.text(b[0] + (0.02 if right else -0.02), min(max(b[1], 0.04), 0.96), str(it.get("text", ""))[:60],
                    transform=ax.transAxes, fontsize=5.5, color=col, ha="left" if right else "right", va="center",
                    bbox=dict(facecolor="black", alpha=0.7, pad=1, edgecolor="none"))


def build_pdf(render, shape, state, patient=None):
    """render(view, i, vmax, cmap, mode, alpha, box) -> PNG (octets) du rendu à l'écran ; shape = (nz, ny, nx).
    state : {views: [{view, i, z, x0, y0} x3], vmax, cmap, mode, alpha, measures: [{v, i, text, items}]}.
    patient : {libellé: valeur} (omis s'il est vide). Retourne les octets du PDF."""
    _check(state)
    patient = {str(k)[:60]: str(v)[:120] for k, v in (patient or {}).items()}
    fig = Figure(figsize=(W, H))
    fig.text(M / W, 1 - 0.7 / H, "Compte rendu TEP/CT", fontsize=16, fontweight="bold", va="center")
    fig.text(1 - M / W, 1 - 0.7 / H, datetime.now().strftime("Exporté le %d/%m/%Y à %H:%M"), fontsize=8, ha="right", va="center", color="#555")
    y = 1.05  # position courante depuis le haut de la page, en pouces

    if patient:
        fig.text(M / W, 1 - y / H, "Patient", fontsize=10, fontweight="bold", va="top")
        y += 0.28
        rows = math.ceil(len(patient) / 2)
        for k, (label, value) in enumerate(patient.items()):
            col, row = divmod(k, rows)
            x = M + col * (W - 2 * M) / 2
            fig.text(x / W, 1 - (y + row * 0.19) / H, f"{label} :", fontsize=7.5, color="#555", va="top")
            fig.text((x + 1.35) / W, 1 - (y + row * 0.19) / H, value, fontsize=7.5, va="top")
        y += rows * 0.19 + 0.25

    imgs = []
    for v in state["views"]:
        png = render(v["view"], int(v["i"]), float(state["vmax"]), state["cmap"], state["mode"], float(state["alpha"]), _box(v))
        imgs.append(mpimg.imread(io.BytesIO(png), format="png"))
    aspects = [im.shape[1] / im.shape[0] for im in imgs]
    gap = 0.2
    h = min((W - 2 * M - 2 * gap) / sum(aspects), 3.8)  # même hauteur pour les 3 vues
    x = M + (W - 2 * M - 2 * gap - h * sum(aspects)) / 2
    n = {"axial": shape[0], "coronal": shape[1], "sagittal": shape[2]}
    for v, im, a in zip(state["views"], imgs, aspects):
        ax = fig.add_axes([x / W, 1 - (y + h) / H, h * a / W, h / H])
        ax.imshow(im, aspect="auto", interpolation="bilinear")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_edgecolor("#999")
        for m in state.get("measures", []):  # mesures posées sur cette vue et cette coupe
            if m.get("v") == v["view"] and int(m.get("i", -1)) == int(v["i"]):
                _draw_items(ax, m.get("items", []), _box(v), (h * a, h))
        fig.text((x + h * a / 2) / W, 1 - (y + h + 0.08) / H, f"{v['view']} {int(v['i']) + 1} / {n[v['view']]}",
                 fontsize=8, ha="center", va="top")
        x += h * a + gap
    y += h + 0.4

    cb = fig.add_axes([M / W, 1 - (y + 0.12) / H, 2.6 / W, 0.12 / H])  # barre de couleurs SUV
    cb.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=state["cmap"])
    cb.set_xticks([]); cb.set_yticks([])
    fig.text((M + 2.7) / W, 1 - (y + 0.06) / H, f"SUV 0 → {float(state['vmax']):g}   (affichage : {state['mode']})", fontsize=7.5, va="center")
    y += 0.4

    fig.text(M / W, 1 - y / H, "Mesures", fontsize=10, fontweight="bold", va="top")
    y += 0.26
    lines = [textwrap.fill(f"• {m.get('text', '')}"[:400], 105, subsequent_indent="   ") for m in state.get("measures", [])] or ["Aucune mesure."]
    for ln in lines[:28]:
        fig.text(M / W, 1 - y / H, ln, fontsize=7.5, va="top")
        y += 0.16 * (ln.count("\n") + 1)
    if len(lines) > 28:
        fig.text(M / W, 1 - y / H, f"… {len(lines) - 28} autres mesures non listées", fontsize=7.5, va="top", color="#555")

    fig.text(0.5, 0.5 * M / H, "Recherche — ne constitue pas un dispositif médical.", fontsize=7, ha="center", color="#777")
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf")
    return buf.getvalue()
