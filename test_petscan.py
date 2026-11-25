import numpy as np
import pytest

from petscan import analyze, measure_region, suv_factor


def test_suv_factor_decay():
    # 70 kg, 350 MBq, demi-vie F-18 109.77 min, scan 60 min après injection
    f = suv_factor(70, 350e6, 109.77 * 60, 3600)
    assert f == pytest.approx(70000 / (350e6 * 2 ** (-60 / 109.77)))
    assert suv_factor(70, 350e6, 6586.2, 0) == pytest.approx(70000 / 350e6)


def test_analyze_two_lesions():
    suv = np.full((20, 20, 20), 1.0)
    suv[2:4, 2:4, 2:4] = 8.0   # 8 voxels
    suv[10:12, 10:12, 10] = 4.0  # 4 voxels
    r = analyze(suv, voxel_ml=0.5)
    assert [x["suv_max"] for x in r] == [8.0, 4.0]
    assert r[0]["mtv_ml"] == pytest.approx(4.0)
    assert r[0]["tlg"] == pytest.approx(8.0 * 4.0)


def test_measure_single_pixel_uses_pixel_size():
    img = np.zeros((20, 20))
    img[10, 10] = 10
    r = measure_region(img, 2.0, 3.0, 10, 10)  # pixels de 2 mm (lignes) x 3 mm (colonnes)
    assert r["length_mm"] == pytest.approx(np.hypot(2, 3))
    assert r["width_mm"] == pytest.approx(2 * 2 * 3 / np.hypot(2, 3))  # étendue d'un rectangle perpendiculairement à sa diagonale


def test_measure_ellipse_axes_and_click_off_center():
    rr, cc = np.mgrid[:60, :60]
    img = np.where(((rr - 30) / 5) ** 2 + ((cc - 30) / 12) ** 2 <= 1, 8.0, 0.0)  # demi-axes 5 et 12 px
    r = measure_region(img, 1.0, 1.0, 31, 36)  # clic décalé : la graine retrouve le max local
    assert 24 <= r["length_mm"] <= 26 and 10 <= r["width_mm"] <= 12
    assert r["suv_max"] == 8.0
    assert measure_region(np.zeros((9, 9)), 1, 1, 4, 4) is None
