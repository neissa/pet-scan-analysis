import numpy as np
import pytest

from organs import ORGANS, segment_organs


def phantom():
    """Fantôme (voxels de 4 mm) : corps cylindrique (40 HU) dans l'air, deux poumons d'air, une sphère chaude de liquide (vessie)."""
    z, y, x = np.mgrid[:60, :70, :70]
    ct = np.full(z.shape, -1000, np.int16)
    ct[((y - 35) ** 2 + (x - 35) ** 2 <= 30 ** 2)] = 40
    suv = np.ones(z.shape, np.float32)
    for cx in (23, 47):  # poumons droit (x faibles) et gauche
        ct[((z - 38) / 14) ** 2 + ((y - 35) / 14) ** 2 + ((x - cx) / 11) ** 2 <= 1] = -800
    bladder = (z - 6) ** 2 + (y - 35) ** 2 + (x - 35) ** 2 <= 5 ** 2
    ct[bladder], suv[bladder] = 10, 15
    return suv, ct


def test_lungs_split_by_side_and_bladder_found():
    suv, ct = phantom()
    lab, info = segment_organs(suv, ct, (4, 4, 4), step=1)
    assert lab.shape == suv.shape
    assert lab[38, 35, 23] == 1 and lab[38, 35, 47] == 2     # poumon droit du patient = x faibles
    assert lab[6, 35, 35] == 8                                # vessie
    assert info[8]["suv_max"] == 15 and info[1]["volume_ml"] > 300
    assert all(info[k]["nom"] == ORGANS[k][0] for k in info)


def test_downsampled_computation_keeps_shape():
    suv, ct = phantom()
    lab, _ = segment_organs(suv, ct, (4, 4, 4), step=2)
    assert lab.shape == suv.shape  # le calcul réduit est remonté à la grille d'origine


def test_no_body_returns_nothing():
    lab, info = segment_organs(np.ones((20, 20, 20), np.float32), np.full((20, 20, 20), -1000, np.int16), (4, 4, 4))
    assert not lab.any() and info == {}


def liver_phantom(hu, seed=0):
    """Fantôme (voxels de 4 mm) : corps graisseux (-80 HU) bruité, poumons, foie ellipsoïde (32×48×48 mm de demi-axes)
    sous la coupole droite, débordant de la ligne médiane, de densité `hu` (injecté ~115, sans contraste ~55)."""
    z, y, x = np.mgrid[:60, :70, :70]
    ct = np.full(z.shape, -1000.0)
    ct[((y - 35) ** 2 + (x - 35) ** 2 <= 30 ** 2)] = -80
    for cx in (23, 47):
        ct[((z - 38) / 14) ** 2 + ((y - 35) / 14) ** 2 + ((x - cx) / 11) ** 2 <= 1] = -800
    liver = ((z - 20) / 8) ** 2 + ((y - 30) / 12) ** 2 + ((x - 26) / 12) ** 2 <= 1
    ct[liver] = hu
    body = ct > -900
    ct[body] += np.random.default_rng(seed).normal(0, 15, body.sum())
    return np.full(z.shape, 2.0, np.float32), ct.astype(np.int16), liver


@pytest.mark.parametrize("hu", [115, 55])
def test_liver_found_injected_and_unenhanced(hu):
    suv, ct, liver = liver_phantom(hu)
    lab, info = segment_organs(suv, ct, (4, 4, 4), step=1)
    found = lab == 4
    assert found[liver].mean() > 0.85 and (found & ~liver).sum() < 0.15 * liver.sum()   # presque tout le foie, peu de fuite
