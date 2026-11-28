import numpy as np
import pytest
from scipy import ndimage

from organs import ORGANS, segment_organs, trace_ureter


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


def ball(c, r):
    z, y, x = np.mgrid[:60, :70, :70]
    return sum(((a - ca) / ra) ** 2 for a, ca, ra in zip((z, y, x), c, np.broadcast_to(r, 3))) <= 1


def abdo_phantom():
    """Corps graisseux (-80 HU, SUV 1) avec poumons, sans organe abdominal : on y pose ceux de chaque test."""
    z, y, x = np.mgrid[:60, :70, :70]
    ct = np.full(z.shape, -1000.0)
    ct[((y - 35) ** 2 + (x - 35) ** 2 <= 30 ** 2)] = -80
    for cx in (23, 47):
        ct[ball((38, 35, cx), (14, 14, 11))] = -800
    return np.ones(z.shape, np.float32), ct


def test_liver_posterior_right_lobe_kept():
    suv, ct = abdo_phantom()
    liver = ball((20, 42, 26), (8, 12, 12))    # lobe droit qui descend jusqu'à la paroi postérieure (y élevés)
    ct[liver] = 115
    lab, _ = segment_organs(np.full(ct.shape, 2.0, np.float32), ct, (4, 4, 4), step=1)
    assert (lab == 4)[liver].mean() > 0.85


def test_trachea_is_not_lung():
    suv, ct = phantom()
    ct[36:, 33:38, 33:38] = -900          # trachée : tube d'air médian qui monte des poumons jusqu'au cou
    ct[36:40, 33:38, 30:41] = -900        # rattachée aux deux poumons
    lab, _ = segment_organs(suv, ct, (4, 4, 4), step=1)
    assert lab[56, 35, 35] == 0 and lab[38, 35, 23] == 1


def test_heart_follows_misregistered_hot_myocardium_but_not_a_node():
    suv, ct = phantom()
    myo, node = ball((35, 30, 38), 4), ball((24, 29, 33), 2)
    suv[myo] = 12.0       # myocarde du TEP décalé vers le poumon gauche par rapport au cœur du CT
    suv[node] = 4.0       # ganglion cardio-phrénique, bien moins fixant
    lab, _ = segment_organs(suv, ct, (4, 4, 4), step=1)
    assert (lab == 3)[myo].mean() > 0.9 and lab[35, 30, 41] == 3 and lab[24, 29, 33] != 3


def test_spleen_excludes_adjacent_hot_lesion():
    suv, ct = abdo_phantom()
    lesion = ball((18, 44, 58), 4)
    spleen = ball((18, 44, 50), (6, 8, 8)) & ~lesion
    ct[spleen], suv[spleen] = 120, 2.0
    ct[lesion], suv[lesion] = 110, 6.0     # nodule périsplénique : même densité, nettement plus fixant
    lab, info = segment_organs(suv, ct, (4, 4, 4), step=1)
    found = lab == 5
    assert found[spleen].mean() > 0.8 and not found[lesion].any()


def test_bladder_whole_and_separate_from_pelvic_mass():
    suv, ct = abdo_phantom()
    bladder, mass = ball((8, 35, 35), 5), ball((8, 35, 44), 4) & ~ball((8, 35, 35), 5)
    ct[bladder], suv[bladder] = 45, 40.0   # urine : CT bruité, la densité ne suffit pas
    ct[mass], suv[mass] = 50, 7.0          # masse pelvienne collée, moins fixante que l'urine
    ct += np.random.default_rng(1).normal(0, 20, ct.shape) * (ct > -900)
    lab, info = segment_organs(suv, ct.astype(np.int16), (4, 4, 4), step=1)
    found = lab == 8
    assert found[bladder].mean() > 0.85 and found[mass].sum() < 0.05 * mass.sum()


def test_kidney_needs_enhanced_parenchyma():
    suv, ct = abdo_phantom()
    kidney, pelvis = ball((10, 45, 20), (9, 6, 5)), ball((10, 45, 20), (5, 3, 2))
    ct[kidney], suv[kidney] = 180, 2.5     # parenchyme rehaussé, peu fixant
    ct[pelvis], suv[pelvis] = 10, 10.0     # cavités dilatées, seule zone très fixante (urine)
    mass = ball((8, 45, 50), 6)
    ct[mass], suv[mass] = 30, 6.0          # masse rétropéritonéale hyperfixante, non rehaussée
    lab, _ = segment_organs(suv, ct, (4, 4, 4), step=1)
    assert (lab == 6)[kidney].mean() > 0.8 and not (lab[mass] == 7).any()


def test_free_fluid_found_and_not_on_plain_fat():
    suv, ct = abdo_phantom()
    lab, info = segment_organs(suv, ct, (4, 4, 4), step=1)
    assert 12 not in info
    fluid = ball((10, 35, 25), (4, 7, 7))
    ct[fluid] = 10
    ct += np.random.default_rng(2).normal(0, 15, ct.shape) * (ct > -900)
    lab, info = segment_organs(suv, ct.astype(np.int16), (4, 4, 4), step=1)
    assert (lab == 12)[fluid].mean() > 0.6 and ((lab == 12) & ~fluid).sum() < 0.3 * fluid.sum()


def test_descending_aorta_and_blood_pool_reference():
    from petscan import blood_pool_reference
    suv, ct = phantom()
    tube = np.zeros(ct.shape, bool)
    tube[20:50] = ball((0, 38, 40), (99, 3, 3))[0]   # aorte descendante, à gauche, devant le rachis
    ct[tube], suv[tube] = 160, 1.8
    lab, info = segment_organs(suv, ct, (4, 4, 4), step=1)
    # cherchée à hauteur du cœur seulement (sous la crosse, au-dessus du diaphragme)
    assert 11 in info and (lab == 11)[30:40][tube[30:40]].mean() > 0.5 and not (lab == 11)[~tube].any()
    ref = blood_pool_reference(suv, lab, (4, 4, 4), 0.8)
    assert ref["mean"] == pytest.approx(1.8 * 0.8) and ref["seuil"] == pytest.approx(2 * ref["mean"] + 2 * ref["sd"])
    assert blood_pool_reference(suv, np.zeros_like(lab), (4, 4, 4)) is None


def test_ureter_follows_curved_hot_tube_else_straight_or_none():
    from test_store import curved_ureter
    sp = (4.0, 4.0, 4.0)
    suv, lab, tube, b3 = curved_ureter()
    z, y, x = np.mgrid[:70, :60, :70]
    body = (y - 30) ** 2 + (x - 35) ** 2 <= 28 ** 2
    path = trace_ureter(suv, lab, sp, 6, 35, body)
    on = [p for p in path if tube[p[0]].any() and lab[tuple(p)] == 0]   # hors rein et vessie
    assert len(on) > 40 and all(np.min(np.abs(np.argwhere(tube[p[0]]) - p[1:]).sum(axis=1)) <= 1 for p in on)
    flat = trace_ureter(np.ones_like(suv), lab, sp, 6, 35, body)   # sans fixation : à peu près la droite
    assert np.linalg.norm(np.diff(flat, axis=0) * 4, axis=1).sum() < 1.2 * np.linalg.norm((flat[0] - flat[-1]) * 4)
    assert trace_ureter(suv, lab, sp, 7, 35, body) is None          # pas de rein gauche
    # seul couloir possible : une hélice de 3 tours, plus de 2 fois la ligne droite -> chemin incohérent, repli
    t = np.linspace(0, 1, 2000)
    hx = np.round(np.stack([8 + 52 * t, 35 + 8 * np.sin(6 * np.pi * t), 18 + 8 * np.cos(6 * np.pi * t)], 1)).astype(int)
    helix = np.zeros(body.shape, bool)
    helix[tuple(hx.T)] = True
    lab2 = np.zeros_like(lab)
    lab2[b3(hx[-1], 3)], lab2[b3(hx[0], 3)] = 6, 8
    zone = ndimage.binary_dilation(helix) | (lab2 > 0)
    assert trace_ureter(suv, lab2, sp, 6, 35, zone) is None


def test_segment_organs_adds_estimated_ureter_between_kidney_and_bladder():
    suv, ct = abdo_phantom()
    kidney, pelvis, bladder = ball((12, 45, 20), (9, 6, 5)), ball((12, 45, 20), (5, 3, 2)), ball((3, 35, 38), 3)
    ct[kidney], suv[kidney] = 180, 2.5
    ct[pelvis], suv[pelvis] = 10, 10.0
    ct[bladder], suv[bladder] = 20, 40.0
    lab, info = segment_organs(suv, ct, (4, 4, 4), step=1)
    assert {6, 8, 13} <= set(info) and 14 not in info and info[13]["nom"] == "Uretère droit (estimé)"
    ct[bladder] = -1000                    # sans vessie : pas d'uretère estimé
    assert 13 not in segment_organs(suv, ct, (4, 4, 4), step=1)[1]


@pytest.mark.parametrize("hu", [115, 55])
def test_liver_found_injected_and_unenhanced(hu):
    suv, ct, liver = liver_phantom(hu)
    lab, info = segment_organs(suv, ct, (4, 4, 4), step=1)
    found = lab == 4
    assert found[liver].mean() > 0.85 and (found & ~liver).sum() < 0.15 * liver.sum()   # presque tout le foie, peu de fuite
