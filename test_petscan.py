import numpy as np
import pytest

from pydicom.dataset import Dataset

from petscan import analyze, ct_choices, default_ct_uid, measure_region, measure_tolerance, suv_factor


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
    z, y, x = r[0]["peak"]  # le pic est dans la lésion de 8, pas dans l'autre
    assert 2 <= z < 4 and 2 <= y < 4 and 2 <= x < 4
    assert r[1]["peak"][2] == 10 and r[1]["peak"][0] in (10, 11)


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


def test_measure_tolerance_disc_and_leak():
    rr, cc = np.mgrid[:80, :80]
    img = np.full((80, 80), -1000.0)
    img[(rr - 40) ** 2 + (cc - 40) ** 2 <= 15 ** 2] = 40.0   # disque de 30 px de diamètre
    r = measure_tolerance(img, 2.0, 2.0, 40, 40, tol=50)    # pixels de 2 mm
    assert 58 <= r["length_mm"] <= 64 and 58 <= r["width_mm"] <= 64 and r["hu_ref"] == 40
    img[(rr - 40) ** 2 + (cc - 40) ** 2 <= 6 ** 2] = 140.0  # zone plus dense au centre
    assert measure_tolerance(img, 2.0, 2.0, 40, 40, tol=30)["length_mm"] < 30   # tolérance étroite : seulement le centre


def _hs(zs, desc, kernel="STANDARD"):
    out = []
    for z in zs:
        d = Dataset()
        d.ImagePositionPatient, d.SeriesDescription, d.ConvolutionKernel = [0, 0, float(z)], desc, kernel
        out.append((f"f{z}", d))
    return out


def test_ct_choices_and_default_series():
    origin = (np.linspace(0, 100, 50), 0.0, 0.0)
    series = {"os": _hs(np.linspace(0, 100, 60), "CT OS", "BONE"), "std": _hs(np.linspace(0, 100, 60), "CT STANDARD"),
              "thorax": _hs(np.linspace(0, 20, 12), "THORAX"), "dose": _hs([50], "Rapport")}
    ch = ct_choices(series, origin)
    assert [c["uid"] for c in ch] == ["os", "std", "thorax"]  # le rapport de dose (1 image) est écarté
    assert ch[0]["couverture"] == 1.0 and ch[2]["couverture"] == 0.2
    assert default_ct_uid(series, origin) == "std"            # à couverture égale : noyau STANDARD
