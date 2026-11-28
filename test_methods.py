import json

import pytest

import methods
import store
from methods import Method
from method_percist import compare

LIV = {"mean": 1.5, "sd": 0.2}
ACQ = {"Activité injectée": "250 MBq", "Poids": "70 kg", "Incorporation": "60 min", "Machine": "Fab Modèle"}
R, D = "01/04/2026", "01/07/2026"


def v(peak, vol=2.0, flags=(), contact=(), mx=None):
    return {"suv_max": mx or peak * 1.3, "sul_max": mx or peak * 1.3, "sul_peak": peak, "vol_ml": vol, "dist_mm": 0.0,
            "flags": list(flags), "contact": list(contact)}


def les(i, since=R, status="valide"):
    return {"id": i, "name": i, "status": status, "since": since}


def run(lesions, ref_vals, cur_vals, liver=LIV, acq=ACQ, unit="SUL"):
    dates = {R: {"values": ref_vals, "acq": ACQ, "liver": LIV, "unit": "SUL"},
             D: {"values": cur_vals, "acq": acq, "liver": liver, "unit": unit}}
    return compare({"date": R, "lesions": lesions, "unit": "SUL"}, dates, {})


@pytest.mark.parametrize("before, after, statut", [(10, 6, "RMP"), (10, 9, "MMS"), (10, 14, "PMP"), (10, 1.2, "RMC"),
                                                   (3.0, 2.0, "RMP"), (2.6, 1.9, "MMS")])  # -27 % : pas RMP
def test_percist_statuts(before, after, statut):
    r = run([les("L1")], {"L1": v(before)}, {"L1": v(after)})
    assert r["statut"] == statut and r["details"][0]["statut"] == statut, r


def test_seuils_absolus_et_relatifs():
    assert run([les("L1")], {"L1": v(2.0)}, {"L1": v(1.3)}, liver={"mean": 1.0, "sd": 0.1})["statut"] == "MMS"  # -35 % mais -0,7
    assert run([les("L1")], {"L1": v(2.0)}, {"L1": v(2.7)})["statut"] == "MMS"                                  # +35 % mais +0,7


def test_cible_la_plus_active_a_chaque_date():
    r = run([les("L1"), les("L2")], {"L1": v(10), "L2": v(4)}, {"L1": v(3), "L2": v(9)})
    d = r["details"][0]
    assert d["reference"] == ["L1", 10] and d["cible"] == ["L2", 9] and r["statut"] == "MMS"


def test_exclusion_excretion_et_contact_physiologique():
    lesions = [les("L1"), les("L2"), les("L3"), les("L4")]
    ref = {"L1": v(10), "L2": v(20, flags=["excrétion possible"]), "L3": v(12, contact=["Vessie"]), "L4": v(3, contact=["Foie"])}
    cur = {"L1": v(6), "L2": v(60, flags=["excrétion possible"]), "L3": v(30, contact=["Vessie"]), "L4": v(3, contact=["Foie"])}
    d = run(lesions, ref, cur)["details"][0]
    assert d["statut"] == "RMP" and d["reference"][0] == "L1"
    assert d["exclues"] == ["L2 (excrétion possible)", "L3 (au contact : Vessie)"]  # le foie n'est pas un organe physiologique exclu


def test_nouvelle_lesion_retenue_seulement():
    ref, cur = {"L1": v(10)}, {"L1": v(6), "L9": v(8)}
    assert run([les("L1")], ref, cur)["statut"] == "RMP"                                 # proposée, non retenue : ignorée
    assert run([les("L1"), les("L9", D, "exclu")], ref, cur)["statut"] == "RMP"          # retenue puis exclue
    r = run([les("L1"), les("L9", D)], ref, cur)
    assert r["statut"] == "PMP" and r["details"][0]["nouvelles"] == ["L9"]
    r = run([les("L1"), les("L9", D)], ref, {**cur, "L9": v(8, flags=["excrétion possible"])})
    assert r["statut"] == "RMP" and "L9 (nouvelle, excrétion possible)" in r["details"][0]["exclues"]


def test_croissance_nette_et_fenetre():
    assert run([les("L1")], {"L1": v(10, vol=2)}, {"L1": v(10, vol=5)})["statut"] == "PMP"
    r = run([les("L1")], {"L1": v(10, vol=2)}, {"L1": v(10, vol=5, flags=["volume limité à la fenêtre"])})
    assert r["statut"] == "MMS" and any("limité à la fenêtre" in w for w in r["avertissements"])


def test_sulmax_de_repli_et_avertissements():
    r = run([les("L1")], {"L1": v(5, vol=0.3, mx=8)}, {"L1": v(4, vol=0.3, mx=4)}, liver=None,
            acq={**ACQ, "Activité injectée": "400 MBq", "Incorporation": "90 min", "Machine": "Autre"}, unit="SUV")
    d, w = r["details"][0], " | ".join(r["avertissements"])
    assert d["metrique"] == "SULmax" and d["reference"][1] == 8 and r["statut"] == "RMP"
    for s in ("SULpeak non fiable", "SUL de repli", "référence foie absente", "activité injectée", "incorporation", "machine"):
        assert s in w, s
    assert ": poids" not in w
    r = run([les("L1")], {"L1": v(10)}, {"L1": v(6)}, liver={"mean": 2.5, "sd": 0.2}, acq={})
    assert any("foie variable" in x for x in r["avertissements"]) and any("non enregistrées" in x for x in r["avertissements"])


def test_non_evaluable():
    assert compare({"date": R, "lesions": [les("L1")]}, {R: {"values": {"L1": v(10)}}}, {})["statut"] == "NE"
    assert run([les("L1", status="exclu")], {"L1": v(10)}, {"L1": v(20)})["statut"] == "NE"


def test_registre_extensible_recist_declare(monkeypatch):
    cat = methods.catalog()
    assert {c["nom"] for c in cat["cancers"]} >= {"Sarcome desmoplastique à petites cellules rondes", "Autre / général"}
    assert all(c["defaut"] == "percist" for c in cat["cancers"])
    rec = next(m for m in cat["methodes"] if m["id"] == "recist11")
    assert rec["disponible"] is False and rec["donnees"] == ["TEP", "CT"]
    with pytest.raises(ValueError):
        methods.choose("dsrct", "recist11")
    with pytest.raises(ValueError):
        methods.choose("inconnu", "percist")
    fake = Method("factice", "Factice", ("TEP",), ("general",), lambda ref, dates, ctx: {"statut": "OK", "details": [], "avertissements": []})
    monkeypatch.setitem(methods.REGISTRY, "factice", fake)
    assert "factice" in [m["id"] for m in methods.catalog()["methodes"]]
    assert methods.choose("general", "factice") == {"cancer": "general", "methode": "factice"}
    with pytest.raises(ValueError):
        methods.choose("dsrct", "factice")              # hors de ses types de cancer
    m = {"ref_date": R, "lesions": [], "analyse": {"cancer": "general", "methode": "factice"}}
    assert methods.avis(m, {R: {"values": {}}})["statut"] == "OK"


def test_choix_memorise_avec_le_patient_puis_supprime():
    from test_viewer import _get, _post, _server
    from viewer import State, demo
    state = State()
    state.set(*demo())
    key = store.patient_key({"Identifiant": "PH-0002", "Naissance": "01/01/1950"}, create=True)
    store.save_map(key, {"version": 1, "ref_date": R, "unit": "SUL", "lesions": [les("L1")],
                         "dates": {R: {"values": {"L1": v(10)}, "acq": ACQ, "liver": LIV, "unit": "SUL"},
                                   D: {"values": {"L1": v(5)}, "acq": ACQ, "liver": LIV, "unit": "SUL"}}})
    srv = _server(state)
    try:
        assert [m["id"] for m in json.loads(_get(srv, "/methods")[1])["methodes"]][:2] == ["percist", "recist11"]
        a = json.loads(_get(srv, f"/avis?key={key}")[1])
        assert a["statut"] == "RMP" and a["choix"] == {"cancer": "general", "methode": "percist"} and "pas un avis médical" in a["mention"]
        body = {"key": key, "op": "method", "cancer": "dsrct", "methode": "percist"}
        assert _post(srv, "/ref/update", body, site="cross-site") == 403
        assert _post(srv, "/ref/update", {**body, "methode": "recist11"})["error"]
        assert _post(srv, "/ref/update", body) == {"ok": True}
        assert json.loads(_get(srv, f"/ref?key={key}")[1])["analyse"] == {"cancer": "dsrct", "methode": "percist"}
        assert json.loads(_get(srv, f"/avis?key={key}")[1])["choix"]["cancer"] == "dsrct"
        assert "error" in json.loads(_get(srv, "/avis?key=" + "0" * 64)[1])
        st, body_ = _get(srv, "/avis?key=../x")
        assert st == 400 and "error" in json.loads(body_)                        # JSON, jamais la page d'erreur HTML
        assert json.loads(_get(srv, f"/data/delete?what=patient&key={key}&confirm=1")[1]) == {"deleted": 1}
        assert store.load_map(key) is None and "error" in json.loads(_get(srv, f"/avis?key={key}")[1])
    finally:
        srv.shutdown()
