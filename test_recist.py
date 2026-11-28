import numpy as np
import pytest

from method_recist import compare
from petscan import ct_diameter

R, D, E = "01/04/2026", "01/07/2026", "01/09/2026"
SP = [3.0, 2.0, 2.0]


def v(long, short=None, flags=(), contact=(), ct=True):
    out = {"suv_max": 8.0, "sul_max": 8.0, "sul_peak": 6.0, "vol_ml": 2.0, "dist_mm": 0.0, "flags": list(flags), "contact": list(contact)}
    if ct is True:
        out["ct"] = None if long is None else {"long_mm": long, "short_mm": short if short is not None else long * 0.8}
    elif ct is not False:
        out["ct"] = ct  # {"why": …}
    return out


def les(i, organe="Foie", since=R, status="valide", ganglion=False):
    return {"id": i, "name": i, "status": status, "since": since, "organe": organe, "ganglion": ganglion}


def run(lesions, *vals, spacing=SP, **extra):
    dates = {d: {"values": x, "score": 0.9, **extra.get(d, {})} for d, x in zip((R, D, E), vals)}
    return compare({"date": R, "lesions": lesions, "unit": "SUL", "spacing": spacing}, dates, {})


# --- mesure au CT : fantômes synthétiques à diamètre connu ---------------------------------------------------------------

def _ct(radius_mm, centre=(10, 50, 50), shape=(20, 100, 100), sp=(3.0, 1.0, 1.0)):
    z, y, x = np.mgrid[:shape[0], :shape[1], :shape[2]].astype(float)
    ct = np.full(shape, -1000, np.int16)
    ct[((y - 50) / 45) ** 2 + ((x - 50) / 45) ** 2 < 1] = -100        # corps : graisse
    ct[((z - centre[0]) * sp[0]) ** 2 + ((y - centre[1]) * sp[1]) ** 2 + ((x - centre[2]) * sp[2]) ** 2 <= radius_mm ** 2] = 40
    return ct, sp


@pytest.mark.parametrize("r", [6, 10, 15])
def test_ct_diametre_sphere_connue(r):
    ct, sp = _ct(r)
    d = ct_diameter(ct, sp, (10, 50, 50), vol_ml=4 / 3 * np.pi * (r / 10) ** 3)
    assert abs(d["long_mm"] - 2 * r) <= 2 and abs(d["short_mm"] - 2 * r) <= 2, d


def test_ct_mesure_non_fiable_jamais_inventee():
    assert ct_diameter(None, SP, (1, 1, 1)) is None
    ct, sp = _ct(8)
    assert "non tissulaire" in ct_diameter(ct, sp, (10, 50, 3))["why"]           # air
    ct[:, 50:53, :] = 40                                                          # bande tissulaire jusqu'au bord de la coupe
    assert "non bornée" in ct_diameter(ct, sp, (10, 50, 50), vol_ml=2.0)["why"]
    ct, sp = _ct(8)
    ct[:, 45:56, 20:80] = 40                                                      # même densité que la lésion, bornée
    assert "fuyante" in ct_diameter(ct, sp, (10, 50, 50), vol_ml=0.5)["why"]


# --- statuts RECIST 1.1 ----------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("before, after, statut", [(40, 27, "RP"), (40, 29, "MS"), (40, 49, "PE"), (40, 47, "MS"), (40, 3, "RC")])
def test_recist_statuts(before, after, statut):
    r = run([les("L1")], {"L1": v(before)}, {"L1": v(after)})
    assert r["statut"] == statut and r["details"][0]["statut"] == statut, r


def test_progression_vs_nadir_et_5_mm():
    r = run([les("L1")], {"L1": v(40)}, {"L1": v(20)}, {"L1": v(25)})                 # nadir 20 → 25 : +25 %, +5 mm
    assert [d["statut"] for d in r["details"]] == ["RP", "PE"] and r["details"][1]["nadir_mm"] == 20
    assert r["details"][1]["var_ref_pct"] == pytest.approx(-37.5)
    r = run([les("L1")], {"L1": v(15)}, {"L1": v(12)}, {"L1": v(15.5)})              # +29 % vs nadir mais +3,5 mm
    assert r["details"][1]["statut"] == "MS"
    t = r["tableau"]
    assert t["colonnes"] == ["Cible", R, D, E] and t["lignes"][-1] == ["Statut", "référence", "MS", "MS"]
    assert t["lignes"][1][1:] == ["15,0", "12,0", "15,5"]


def test_selection_cibles_5_max_2_par_organe_et_seuils():
    lesions = [les("F1"), les("F2"), les("F3"), les("P1", "Poumon droit"), les("P2", "Poumon droit"), les("R1", "Rate"),
               les("G1", None, ganglion=True), les("G2", None, ganglion=True), les("S1", "Rate")]
    vals = {"F1": v(50), "F2": v(30), "F3": v(45), "P1": v(25), "P2": v(20), "R1": v(12), "G1": v(30, 16), "G2": v(30, 12),
            "S1": v(9)}
    r = run(lesions, vals, vals)
    assert r["cibles"] == ["F1", "F3", "P1", "P2", "G1"]                         # F2 : 2 par organe ; R1 : 5 cibles atteintes
    nc = " | ".join(r["non_cibles"])
    for s in ("F2 (2 cibles déjà retenues : Foie)", "R1 (5 cibles déjà retenues)", "G2 (ganglion, petit axe 12,0 mm < 15 mm)",
              "S1 (9,0 mm < 10,0 mm)"):
        assert s in nc, s
    assert r["details"][0]["somme_mm"] == 50 + 45 + 25 + 20 + 16                   # ganglion : petit axe
    assert run([les("L1")], {"L1": v(11)}, {"L1": v(11)}, spacing=[6.0, 2, 2])["cibles"] == []   # coupe de 6 mm : seuil 12 mm


def test_ganglion_rc_sous_10_mm():
    lesions = [les("G1", ganglion=True), les("L1")]
    assert run(lesions, {"G1": v(25, 18), "L1": v(20)}, {"G1": v(12, 9), "L1": v(4)})["statut"] == "RC"
    assert run(lesions, {"G1": v(25, 18), "L1": v(20)}, {"G1": v(14, 11), "L1": v(4)})["statut"] == "RP"


def test_nouvelle_lesion_et_non_cibles_sans_equivoque():
    ref, cur = {"L1": v(40)}, {"L1": v(30), "L9": v(15)}
    assert run([les("L1")], ref, cur)["statut"] == "MS"
    r = run([les("L1"), les("L9", since=D)], ref, cur)
    assert r["statut"] == "PE" and r["details"][0]["nouvelles"] == ["L9"]
    r = run([les("L1"), les("L9", since=D)], ref, {**cur, "L9": v(15, flags=["excrétion possible"])})
    assert r["statut"] == "MS" and "L9 (nouvelle, excrétion possible)" in r["details"][0]["exclues"]
    r = run([les("L1")], ref, cur, **{D: {"pe_non_cible": True}})
    assert r["statut"] == "PE" and r["details"][0]["pe_non_cible"] is True


def test_exclusion_excretion_contact_et_non_mesurable():
    lesions = [les("L1"), les("L2", "Vessie"), les("L3"), les("L4"), les("L5")]
    ref = {"L1": v(40), "L2": v(30, contact=["Vessie"]), "L3": v(30), "L4": v(30, ct={"why": "région fuyante"}), "L5": v(None)}
    cur = {"L1": v(20), "L2": v(30, contact=["Vessie"]), "L3": v(30, flags=["excrétion possible"]), "L4": v(30), "L5": v(None)}
    r = run(lesions, ref, cur)
    assert r["cibles"] == ["L1"] and r["statut"] == "RP"
    assert r["details"][0]["exclues"] == ["L2 (au contact : Vessie)", "L3 (excrétion possible)"]   # signalée à une date : exclue
    w = " | ".join(r["avertissements"])
    for s in ("L4 (région fuyante)", "L5 (pas de CT)", "non diagnostiques", "tolérance HU", "Coupes de 3,0 mm"):
        assert s in w, s


def test_non_evaluable_et_dates_sans_diametre():
    r = run([les("L1")], {"L1": v(40)}, {"L1": v(None, ct=False)})                  # date enregistrée avant la mesure CT
    assert r["statut"] == "NE" and any("diamètre non enregistré" in w for w in r["avertissements"])
    r = run([les("L1")], {"L1": v(40)}, {"L1": v(None)}, **{D: {"score": None}})    # date sans CT, recalage DICOM
    w = " | ".join(r["avertissements"])
    assert r["statut"] == "NE" and "pas de CT" in w and "coordonnées DICOM" in w
    r = run([les("L1")], {"L1": v(40, ct=False)}, {"L1": v(30)})                    # référence sans diamètre : aucune cible
    assert r["statut"] == "NE" and r["cibles"] == []
    assert compare({"date": R, "lesions": [les("L1")]}, {R: {"values": {"L1": v(40)}}}, {})["statut"] == "NE"
