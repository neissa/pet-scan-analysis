import json

import numpy as np
import pytest

import relevance
from relevance import EXCLU_MAX, score

# lésion tissulaire isolée et contrastée (champs de petscan.analyze + organe de viewer)
TISSU = {"suv_max": 8.0, "suv_mean": 4.0, "sul_max": 6.0, "mtv_ml": 4.0, "physio": None, "organe": None}
CTX = {"foie": 1.5, "hu": 40.0}


def _check(r):
    assert 0 <= r["score"] <= 100 and r["classe"] in ("haute", "moyenne", "faible")
    assert sum(c["points"] for c in r["criteres"]) == r["brut"]          # détail cohérent avec le score avant bornage
    assert all({"nom", "valeur", "points"} <= set(c) for c in r["criteres"])
    return r


def test_physiological_and_excretion_score_below_isolated_tissue_lesion():
    iso = _check(score(TISSU, CTX))
    assert iso["classe"] == "haute"
    for les, ctx in (({**TISSU, "physio": "Vessie"}, CTX),                           # au contact de la vessie (analyze)
                     ({**TISSU, "organe": "dans Cœur (approx.)"}, CTX),             # dans le cœur
                     ({**TISSU, "organe": "près de Rein droit"}, {**CTX, "excretion": True}),
                     ({"suv_max": 25.0, "sul_max": 20.0, "vol_ml": 5.0, "flags": ["excrétion possible"], "contact": []}, CTX)):
        r = _check(score(les, ctx))
        assert r["score"] < iso["score"] - 25, r
    assert score({**TISSU, "organe": "près de Rate"}, CTX)["score"] == iso["score"]  # organe non physiologique : neutre


def test_contrast_bonus_dropped_near_physiological_organ():
    hot, warm = {**TISSU, "sul_max": 30.0}, {**TISSU, "sul_max": 3.0}   # contraste 20 → +20 (plafond), 2 → +5
    assert score(hot, CTX)["brut"] - score(warm, CTX)["brut"] == 15
    assert score({**hot, "physio": "Vessie"}, CTX)["brut"] == score({**warm, "physio": "Vessie"}, CTX)["brut"]


def test_small_volume_penalized():
    sc = [score({**TISSU, "mtv_ml": v}, CTX)["score"] for v in (0.1, 0.5, 4.0)]
    assert sc[0] < sc[1] < sc[2]


def test_missing_criteria_are_neutral_and_listed():
    r = _check(score({"suv_max": 8.0, "mtv_ml": 4.0}))   # ni foie, ni CT, ni organes, ni SUVmean
    ne = [c for c in r["criteres"] if c["valeur"].startswith("non évalué")]
    assert {c["nom"] for c in ne} >= {"contraste (SULmax / foie)", "densité CT au pic", "localisation"}
    assert all(c["points"] == 0 for c in ne)
    assert r["brut"] == relevance.BASE + relevance.VOL_OK


def test_no_ct_is_neutral_air_penalized():
    base = score(TISSU, {"foie": 1.5})
    assert base["brut"] == score(TISSU, CTX)["brut"] - relevance.HU_TISSU_PTS
    assert score(TISSU, {"foie": 1.5, "hu": -900})["score"] < base["score"]


def test_follow_status_validated_above_unvalidated_and_excluded_capped():
    plain = score(TISSU, {**CTX, "suivi": {"status": None, "dates": 1}})["score"]
    val = score(TISSU, {**CTX, "suivi": {"status": "valide", "dates": 3, "ref": True}})["score"]
    exc = _check(score(TISSU, {**CTX, "suivi": {"status": "exclu", "dates": 3, "ref": True}}))
    assert val > plain and exc["score"] <= EXCLU_MAX < exc["brut"]


def test_score_bounds_and_sort_order():
    worst = {"suv_max": 2.6, "suv_mean": 2.55, "sul_max": 1.0, "mtv_ml": 0.05, "physio": "Vessie", "organe": "dans Vessie",
             "flags": ["volume limité à la fenêtre", "excrétion possible"]}
    best = {**TISSU, "sul_max": 60.0}
    lo, hi = _check(score(worst, {"foie": 3.0, "hu": -950})), _check(score(best, {**CTX, "suivi": {"status": "valide", "dates": 5, "ref": True}}))
    assert lo["score"] == 0 > lo["brut"] and hi["score"] == 100 < hi["brut"]
    les = [relevance.annotate(dict(l), CTX) for l in (worst, TISSU, {**TISSU, "physio": "Cœur (approx.)"})]
    order = sorted(les, key=lambda l: -l["pertinence"])
    assert order[0] is les[1] and order[-1] is les[0]


def test_lesions_route_adds_relevance():
    from test_petscan import _phantom
    from test_viewer import _get, _server
    from viewer import State
    suv, sp, lab = _phantom()
    ct = np.where(lab > 0, 40, -1000).astype(np.float32)
    ct[50, 90, 20] = 40                       # lésion isolée : tissu au pic
    state = State()
    state.set(suv, sp, ct)
    state.organs = (lab, {})
    srv = _server(state)
    try:
        les = json.loads(_get(srv, "/lesions?thr=2.5")[1])
        assert all({"pertinence", "pertinence_classe", "pertinence_detail"} <= set(l) for l in les)
        assert [l["suv_max"] for l in les] == sorted((l["suv_max"] for l in les), reverse=True)   # ordre serveur inchangé
        vessie, iso = next(l for l in les if l["physio"] == "Vessie"), next(l for l in les if l["suv_max"] == 20.0)
        assert vessie["pertinence"] < iso["pertinence"] - 25 and vessie["pertinence_classe"] == "faible"
        assert any(c["nom"].startswith("contraste") and c["points"] > 0 for c in iso["pertinence_detail"])  # foie trouvé
        assert any("débord" in c["valeur"] for c in vessie["pertinence_detail"])
    finally:
        srv.shutdown()
