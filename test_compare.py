import numpy as np

from compare import estimate_shift, shift_volume


def phantom():
    """CT synthétique : corps cylindrique (40 HU) avec un « squelette » asymétrique (colonne + deux côtes) à 500 HU."""
    z, y, x = np.mgrid[:80, :90, :90]
    ct = np.full(z.shape, -1000, np.int16)
    ct[(y - 45) ** 2 + (x - 45) ** 2 <= 32 ** 2] = 40
    ct[((y - 62) ** 2 + (x - 45) ** 2 <= 6 ** 2)] = 500                        # colonne
    ct[((z - 30) / 5) ** 2 + ((y - 40) / 4) ** 2 + ((x - 22) / 10) ** 2 <= 1] = 500   # côte gauche, plus haute
    ct[((z - 55) / 5) ** 2 + ((y - 40) / 4) ** 2 + ((x - 70) / 10) ** 2 <= 1] = 500   # côte droite
    return ct


def test_shift_volume_copies_and_fills():
    v = np.arange(24, dtype=np.int16).reshape(2, 3, 4) + 1
    out = shift_volume(v, (1, 0, -1), (2, 3, 4), fill=-9)          # out[z, y, x] = v[z+1, y, x-1]
    assert out[0, 1, 1] == v[1, 1, 0] and out[0, 1, 0] == -9 and (out[1] == -9).all()
    assert (shift_volume(v, (5, 0, 0), (2, 3, 4), fill=7) == 7).all()  # décalage hors du volume : tout rempli


def test_estimate_shift_recovers_known_translation():
    a = phantom()
    d0 = (4, -6, 8)
    b = shift_volume(a, d0, a.shape, fill=-1000)       # b[x] = a[x + d0]  =>  b[x - d0] = a[x]  =>  d = -d0
    d, score = estimate_shift(a, b)
    assert all(abs(got - (-want)) <= 2 for got, want in zip(d, d0)), d
    assert score > 0.5
    d_same, score_same = estimate_shift(a, a)
    assert d_same == (0, 0, 0) and score_same > 0.9


def test_aligned_volume_matches_reference():
    a = phantom()
    b = shift_volume(a, (2, 4, -4), a.shape, fill=-1000)
    d, _ = estimate_shift(a, b)
    aligned = shift_volume(b, d, a.shape, fill=-1000)
    inside = (a > 0)
    assert (aligned[inside] == a[inside]).mean() > 0.9   # le squelette redevient superposable


def test_refine_shift_recovers_residual_offset_global_and_local():
    from compare import refine_shift
    a = phantom()
    true_e = (2, -4, 6)                      # B décalée de true_e par rapport à A
    b = shift_volume(a, (-true_e[0], -true_e[1], -true_e[2]), a.shape, fill=-1000)   # b[x] = a[x - e]  =>  b[x + e] = a[x]
    e, score = refine_shift(a, b, step=2, margin=6)
    assert all(abs(g - w) <= 2 for g, w in zip(e, true_e)) and score > 0.5, (e, score)
    e2, score2 = refine_shift(a, b, center=(40, 45, 45), half=(30, 30, 30), step=1, margin=10)   # zone locale, pas de 1 voxel
    assert e2 == true_e and score2 > 0.5, (e2, score2)
    # on part des décalages actuels : si B est déjà décalée de e, le résidu est nul
    e3, _ = refine_shift(a, b, db=true_e, step=1, margin=6, center=(40, 45, 45), half=(30, 30, 30))
    assert e3 == (0, 0, 0)
