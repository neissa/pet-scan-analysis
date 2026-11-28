"""Recalage de deux séances d'un même patient (positions sur la table différentes).
Translation entière estimée par corrélation croisée (FFT) du squelette vu au CT ; pas de rotation ni de déformation.
Grilles d'espacements différents : `resample` ramène un volume sur l'espacement de l'autre avant le recalage."""
import warnings

import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve

from organs import _body

# SciPy < 1.10 avertit à chaque affine_transform à matrice diagonale (1-D), usage voulu ici et dans viewer.py ; filtre global (thread-safe, contrairement à catch_warnings)
warnings.filterwarnings("ignore", message="The behavior of affine_transform with a 1-D array", category=UserWarning)


def _bone(ct, step):
    """Squelette (os dans le corps, sans la table) sur une grille réduite, légèrement flouté pour une corrélation plus robuste."""
    c = ct[::step, ::step, ::step]
    return ndimage.gaussian_filter((_body(c) & (c > 250)).astype(np.float32), 1.0)


def estimate_shift(ct_a, ct_b, step=2):
    """Décalage entier d = (dz, dy, dx) en voxels tel que b[z+dz, y+dy, x+dx] ≈ a[z, y, x], et un score de 0 à 1
    (corrélation normalisée du squelette ; un bon recalage dépasse nettement 0,3). Précision : `step` voxels."""
    fa, fb = _bone(ct_a, step), _bone(ct_b, step)
    corr = fftconvolve(fa, fb[::-1, ::-1, ::-1], mode="full")
    k = np.unravel_index(int(np.argmax(corr)), corr.shape)
    d = tuple((n - 1 - i) * step for n, i in zip(fb.shape, k))
    score = float(corr[k] / np.sqrt((fa ** 2).sum() * (fb ** 2).sum() or 1.0))
    return d, score


def shift_volume(vol, d, shape, fill=0):
    """Volume de forme `shape` tel que out[z, y, x] = vol[z+dz, y+dy, x+dx] ; `fill` hors du volume d'origine."""
    out = np.full(shape, fill, vol.dtype)
    dst, src = [], []
    for n_out, n_src, k in zip(shape, vol.shape, d):
        a, b = max(0, -k), min(n_out, n_src - k)  # indices de sortie dont la source existe
        if b <= a:
            return out
        dst.append(slice(a, b))
        src.append(slice(a + k, b + k))
    out[tuple(dst)] = vol[tuple(src)]
    return out


def same_grid(sp_a, sp_b, tol=0.005):
    """Espacements égaux à tol près (relatif) : pas de rééchantillonnage (1,37 contre 1,35 mm dérive déjà de 4 voxels sur 300)."""
    return all(abs(a - b) <= tol * a for a, b in zip(sp_a, sp_b))


def resample(vol, spacing, new_spacing, order=1, fill=0):
    """Volume sur une grille de même origine (centre du voxel 0, position DICOM inchangée) et d'espacement new_spacing (mm),
    couvrant la même étendue : out[i] = vol(i * new_spacing / spacing). order 1 = linéaire (TEP, CT), 0 = plus proche voisin
    (étiquettes) ; fill hors du volume. Indice de grille d'origine v -> nouvelle grille v * spacing / new_spacing."""
    r = [t / s for s, t in zip(spacing, new_spacing)]
    shape = tuple(int(np.floor((n - 1) / k + 1e-6)) + 1 for n, k in zip(vol.shape, r))
    return ndimage.affine_transform(vol, r, output_shape=shape, order=order, mode="constant", cval=fill)


def _crop(vol, lo, hi, fill):
    """vol[lo:hi] par axe, les indices hors du volume étant remplis avec `fill` (pas d'indices négatifs interprétés)."""
    out = np.full(tuple(h - l for l, h in zip(lo, hi)), fill, vol.dtype)
    src = tuple(slice(max(0, l), min(n, h)) for l, h, n in zip(lo, hi, vol.shape))
    dst = tuple(slice(s.start - l, s.stop - l) for s, l in zip(src, lo))
    if all(s.stop > s.start for s in src):
        out[dst] = vol[src]
    return out


def refine_shift(ct_a, ct_b, da=(0, 0, 0), db=(0, 0, 0), center=None, half=None, step=2, margin=10):
    """Alignement fin automatique de B sur A, en partant des décalages manuels actuels da (A) et db (B), en voxels
    (convention de `shift_volume`). Retourne e = (ez, ey, ex) à ajouter à db, et le score (0 à 1).
    center, half : zone locale (centre et demi-taille en voxels, par axe) ; None = volume entier.
    step : pas de la grille de calcul ; margin : amplitude de recherche en pas (±margin*step voxels autour de db).
    Corrélation normalisée du squelette (les os d'une séance sur l'autre) ; B est cherchée autour de la zone de A."""
    pad = margin * step
    if center is None:  # volume entier
        shape = ct_a.shape
        lo, hi = (0, 0, 0), shape
    else:               # zone locale : on ne calcule le squelette que dans cette zone (et la marge de recherche autour pour B)
        lo = tuple(max(0, c - h) for c, h in zip(center, half))
        hi = tuple(min(n, c + h + 1) for c, h, n in zip(center, half, ct_a.shape))
    # volumes affichés = ct[q + d] ; _crop lit directement la zone voulue, sans copier le volume entier
    ha = _bone(_crop(ct_a, [l + d for l, d in zip(lo, da)], [h + d for h, d in zip(hi, da)], -1000), step)
    hb = _bone(_crop(ct_b, [l - pad + d for l, d in zip(lo, db)], [h + pad + d for h, d in zip(hi, db)], -1000), step)
    corr = fftconvolve(hb, ha[::-1, ::-1, ::-1], mode="valid")
    energy = np.maximum(fftconvolve(hb * hb, np.ones_like(ha), mode="valid"), 0)  # énergie de B sous chaque position (jamais < 0 : arrondi FFT)
    ncc = corr / np.sqrt(energy * float((ha * ha).sum()) + 1e-6)
    ncc[energy < 0.05 * energy.max()] = 0  # positions où B est presque vide : pas de verdict
    k = np.unravel_index(int(np.argmax(ncc)), ncc.shape)
    return tuple(int(i - margin) * step for i in k), float(ncc[k])
