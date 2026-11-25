"""Recalage de deux séances d'un même patient (même grille TEP, positions sur la table différentes).
Translation entière estimée par corrélation croisée (FFT) du squelette vu au CT ; pas de rotation ni de déformation."""
import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve

from organs import _body


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
