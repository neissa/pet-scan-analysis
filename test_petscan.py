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


def test_scan_reports_progress(tmp_path):
    from petscan import _scaled, _scan
    for i in range(60):
        (tmp_path / f"f{i}").write_bytes(b"pas du dicom")
    seen = []
    assert _scan(tmp_path, "PT", _scaled(seen.append, 0.2, 0.6)) == {}
    assert seen and seen == sorted(seen) and 0.2 <= seen[0] and seen[-1] < 0.6   # croissant, dans la plage demandée


def test_n_workers_sequential_on_android(monkeypatch):
    from petscan import n_workers
    monkeypatch.setenv("ANDROID_DATA", "/data")  # présent dans l'environnement de toute appli Android
    assert n_workers() == 1


def test_lbm_james_reference_values_and_sul_fallback():
    from petscan import lbm_james, sul_info
    # calculs à la main : 1.10·80 − 128·(80/180)² = 88 − 25.28395 ; 1.07·60 − 148·(60/165)² = 64.2 − 19.57025
    assert lbm_james(80, 1.80, "M") == pytest.approx(62.71605, abs=1e-4)
    assert lbm_james(60, 1.65, "F") == pytest.approx(44.62975, abs=1e-4)
    s = sul_info(80, 1.80, "M")
    assert s["sul"] and s["unit"] == "SUL" and s["factor"] == pytest.approx(62.71605 / 80, abs=1e-6)
    for w, h, sex in [(70, None, "M"), (70, "", "F"), (70, 1.7, "O"), (70, 1.7, None), (0, 1.7, "M"), (70, 170, "M")]:
        f = sul_info(w, h, sex)  # taille absente, en cm, sexe inconnu, poids nul : repli SUV poids signalé, sans planter
        assert f["sul"] is False and f["factor"] == 1.0 and f["unit"] == "SUV" and "indisponible" in f["label"]


def _phantom():
    """Fantôme (z,y,x) anisotrope : foie homogène SUV 2 avec une métastase, lésion sphérique SUV 8 avec un voxel à 20,
    vessie très captante et une lésion proche de la vessie, une lésion isolée loin de tout."""
    sp = (3.0, 2.0, 2.0)
    z, y, x = np.mgrid[:60, :100, :100]
    mm = lambda c: ((z - c[0]) * sp[0]) ** 2 + ((y - c[1]) * sp[1]) ** 2 + ((x - c[2]) * sp[2]) ** 2
    suv = np.full(z.shape, 1.0, np.float32)
    lab = np.zeros(z.shape, np.uint8)
    lab[5:35, 5:55, 5:55] = 4                    # foie : 90 × 100 × 100 mm
    suv[lab == 4] = 2.0
    suv[mm((10, 30, 30)) <= 36] = 6.0            # métastase hépatique (r 6 mm) : la sphère de référence doit l'éviter
    suv[mm((45, 75, 75)) <= 144] = 8.0           # lésion r 12 mm
    suv[45, 75, 75] = 20.0
    lab[mm((45, 30, 80)) <= 100] = 8             # vessie r 10 mm
    suv[lab == 8] = 15.0
    suv[mm((45, 30, 89)) <= 16] = 5.0            # lésion à 4 mm du bord de la vessie, séparée d'elle par un voxel de SUV 1
    suv[mm((50, 90, 20)) <= 16] = 4.0            # lésion isolée
    return suv, sp, lab


def test_sul_peak_physio_and_liver_reference():
    from petscan import PEAK_R_MM, _sphere, liver_reference
    suv, sp, lab = _phantom()
    f = 0.75
    les = analyze(suv, float(np.prod(sp)) / 1000, 2.5, sp, f, lab, margin_mm=10)
    big = next(l for l in les if l["suv_max"] == 20.0)
    w = _sphere(PEAK_R_MM, sp, sub=5)             # poids des voxels (volume partiel)
    assert w.sum() * 12 == pytest.approx(1000, rel=0.03) and w[tuple(n // 2 for n in w.shape)] == 1   # 1 cm³ à 3 % près
    # boule centrée sur le voxel à 20, tous ses autres voxels à 8 : moyenne pondérée = 8 + 12 / somme des poids
    assert big["sul_peak"] == pytest.approx(f * (8 + 12 / w.sum()), rel=1e-5)   # ni le max ni la moyenne
    assert big["sul_max"] == pytest.approx(20 * f) and big["tlg_sul"] == pytest.approx(f * big["tlg"])
    assert big["physio"] is None and big["valide"] is True
    assert not any(l["suv_max"] == 15.0 for l in les)                       # vessie retirée voxel par voxel
    near = next(l for l in les if l["suv_max"] == 5.0)
    assert near["physio"] == "Vessie" and near["valide"] is True             # reste de la lésion au contact : signalé
    assert next(l for l in les if l["suv_max"] == 4.0)["physio"] is None
    assert next(l for l in les if l["suv_max"] == 6.0)["physio"] is None   # foie : pas une fixation physiologique exclue
    assert "physio" not in analyze(suv, 0.012, 2.5, sp, f)[0]               # sans organes : champ absent
    ref = liver_reference(suv, lab, sp, f)
    assert ref["mean"] == pytest.approx(2.0 * f) and ref["sd"] == pytest.approx(0, abs=1e-6)
    assert ref["seuil"] == pytest.approx(1.5 * ref["mean"] + 2 * ref["sd"])
    c = np.array(ref["centre"])
    assert lab[tuple(c)] == 4 and np.sqrt((((c - (10, 30, 30)) * sp) ** 2).sum()) > 15 + 6   # boule hors de la métastase
    noisy = suv + np.where(np.indices(suv.shape).sum(0) % 2, 0.1, -0.1).astype(np.float32)   # damier ±0.1 : écart-type 0.1
    r2 = liver_reference(noisy, lab, sp, 1.0)
    assert r2["sd"] == pytest.approx(0.1, rel=0.05) and r2["seuil"] == pytest.approx(1.5 * r2["mean"] + 0.2, rel=0.01)
    assert liver_reference(suv, np.zeros_like(lab), sp) is None                # foie introuvable
    small = np.zeros_like(lab)
    small[20:26, 40:50, 40:50] = 4
    assert liver_reference(suv, small, sp) is None                             # trop petit pour une boule de 3 cm


def test_liver_reference_fallback_on_thin_liver():
    from petscan import liver_reference
    sp, f = (3.0, 2.0, 2.0), 0.8
    suv = np.full((40, 30, 60), 2.0, np.float32)
    lab = np.zeros(suv.shape, np.uint8)
    lab[5:35, 5:20, 5:55] = 4                    # 30 mm d'épaisseur en y : distance max au bord 16 mm (< 15 + 5)
    ref = liver_reference(suv, lab, sp, f)
    assert ref["erode_mm"] == 0 and ref["diam_mm"] == 30 and not ref["sphere_reduite"] and "0 mm" in ref["msg"]
    assert ref["mean"] == pytest.approx(2.0 * f) and ref["seuil_suv"] == pytest.approx(ref["seuil"] / f)
    lab[:] = 0
    lab[5:35, 5:17, 5:55] = 4                    # 24 mm : seule une sphère de 20 mm tient
    ref = liver_reference(suv, lab, sp, f)
    assert ref["erode_mm"] == 0 and ref["diam_mm"] == 20 and ref["sphere_reduite"] and "20 mm" in ref["msg"]
    full = np.zeros(suv.shape, np.uint8)
    full[5:35, 2:28, 5:55] = 4
    assert liver_reference(suv, full, sp, f)["msg"] == ""                     # foie assez épais : pas de repli


def _fused(organ, organ_suv, lesion_c, bridge=False):
    """Fantôme : organe physiologique (sphère r 10 mm, étiquette `organ`) et lésion SUV 7 (r 10 mm) collée ou reliée."""
    sp = (3.0, 2.0, 2.0)
    z, y, x = np.mgrid[:40, :60, :80]
    mm = lambda c: ((z - c[0]) * sp[0]) ** 2 + ((y - c[1]) * sp[1]) ** 2 + ((x - c[2]) * sp[2]) ** 2
    suv = np.ones(z.shape, np.float32)
    lab = np.zeros(z.shape, np.uint8)
    lab[mm((20, 30, 20)) <= 100] = organ
    suv[lab == organ] = organ_suv
    les = mm(lesion_c) <= 100
    suv[les] = 7.0
    if bridge:
        suv[20, 30, 20:lesion_c[2]] = np.maximum(suv[20, 30, 20:lesion_c[2]], 3.0)
    return suv, sp, lab, les


def test_lesion_fused_with_bladder_measured_without_bladder_voxels():
    suv, sp, lab, les = _fused(8, 15.0, (20, 30, 29))        # centres à 18 mm : la lésion touche la vessie
    vml = float(np.prod(sp)) / 1000
    out = analyze(suv, vml, 2.5, sp, 1.0, lab, margin_mm=10)
    assert len(out) == 1
    r = out[0]
    assert r["suv_max"] == 7.0 and r["physio"] == "Vessie" and r["valide"] is True   # pas la vessie (15), signalée, retenue
    assert 0 < r["mtv_ml"] < les.sum() * vml


def test_retroperitoneal_lesion_near_kidney_not_excluded():
    suv, sp, lab, les = _fused(6, 5.0, (20, 30, 40), bridge=True)   # à 20 mm du bord du rein, reliée par un pont SUV 3
    vml = float(np.prod(sp)) / 1000
    out = analyze(suv, vml, 2.5, sp, 1.0, lab, margin_mm=10)
    r = next(l for l in out if l["suv_max"] == 7.0)
    assert r["valide"] is True and r["mtv_ml"] >= les.sum() * vml                      # lésion entière, retenue
    assert not any(l["suv_max"] == 5.0 for l in out)                                   # le rein lui-même est retiré
