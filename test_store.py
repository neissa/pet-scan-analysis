import json

import numpy as np
import pytest

import store
from petscan import measure_at, ureter_paths


def test_migration_of_recent_then_old_file_removed():
    old = store.LEGACY[0]
    old.write_text(json.dumps({"recent": ["/a", "/b"]}), encoding="utf-8")
    c = store.config()
    assert c["recent"] == ["/a", "/b"] and not old.exists() and (store.DIR / "config.json").is_file()
    assert store.recent() == []                  # récents désactivés par défaut, même migrés
    store.remember("/c")
    assert store.config()["recent"] == ["/a", "/b"]  # rien d'ajouté tant que l'option est coupée
    store.set_config(recent_on=True)
    store.remember("/c")
    assert store.recent() == ["/c", "/a", "/b"]


def test_pseudonym_stable_distinct_and_no_identity_in_files():
    p1 = {"Nom": "Alice TEST", "Identifiant": "ID-123456", "Naissance": "02/03/1961"}
    p2 = {**p1, "Identifiant": "ID-654321"}
    assert store.patient_key(p1) is None and not store.DIR.exists()   # pas de sel : rien écrit avant une référence
    k1 = store.patient_key(p1, create=True)
    assert k1 == store.patient_key(p1) == store.patient_key({**p1, "Nom": "autre"})   # stable (identifiant + naissance)
    assert store.patient_key(p2) != k1 and store.patient_key({**p1, "Naissance": "01/01/1960"}) != k1
    assert store.patient_key({"Nom": "sans id"}) is None
    store.save_map(k1, {"ref_date": "01/04/2026", "lesions": [{"id": "L1"}], "dates": {"01/04/2026": {"values": {}}}})
    blob = b"".join(p.read_bytes() for p in store.DIR.rglob("*") if p.is_file())
    for s in ("Alice", "TEST", "ID-123456", "1961", "02/03"):
        assert s.encode() not in blob
    with pytest.raises(ValueError):
        store.load_map("../../etc/passwd")      # la clé vient de la page : jamais un chemin


def test_inventory_preview_and_delete_patient_recent_all():
    k = store.patient_key({"Identifiant": "X1"}, create=True)
    k2 = store.patient_key({"Identifiant": "X2"})
    for key in (k, k2):
        store.save_map(key, {"ref_date": "01/04/2026", "lesions": [{"id": "L1"}, {"id": "L2"}],
                             "dates": {"01/04/2026": {}, "01/07/2026": {}}})
    store.set_config(recent_on=True)
    store.remember("/x")
    inv = store.inventory()
    assert (inv["patients"], inv["lesions"], inv["dates"], inv["recents"]) == (2, 4, 4, 1) and inv["octets"] > 0
    assert [f["nom"] for f in store.preview("patient", k)] == [f"patients/{k}.json"]
    assert store.delete("patient", k) == 1 and store.load_map(k) is None and store.load_map(k2)
    assert store.delete("recent") == 1 and store.config()["recent"] == []
    assert store.delete("all") >= 2 and not store.DIR.exists()


def test_measure_at_excretion_flag_contact_and_physio_exclusion():
    sp = (2.0, 2.0, 2.0)
    z, y, x = np.mgrid[:50, :60, :60]
    mm = lambda c: ((z - c[0]) * 2) ** 2 + ((y - c[1]) * 2) ** 2 + ((x - c[2]) * 2) ** 2
    suv = np.ones(z.shape, np.float32)
    lab = np.zeros(z.shape, np.uint8)
    lab[mm((25, 30, 20)) <= 100] = 8           # vessie r 10 mm, urine SUV 60
    suv[lab == 8] = 60
    suv[mm((25, 30, 29)) <= 16] = 5.0          # lésion collée à la vessie
    suv[mm((10, 50, 50)) <= 16] = 6.0          # lésion loin de tout
    r = measure_at(suv, sp, (25, 30, 29), labels=lab)
    assert r["suv_max"] == 5.0 and "excrétion possible" in r["flags"] and "Vessie" in r["contact"]   # l'urine est écartée
    assert r["vol_ml"] > 0 and r["sul_peak"] > 0
    far = measure_at(suv, sp, (10, 50, 50), labels=lab)
    assert far["suv_max"] == 6.0 and far["flags"] == [] and far["contact"] == []
    assert "excrétion possible" in measure_at(suv, sp, (25, 30, 20))["flags"]  # sans organes : max très élevé
    assert measure_at(suv, sp, (99, 0, 0)) is None
    lab[10:14, 30:34, 50:54] = 6               # rein : trajet urétéral jusqu'à la vessie
    (a, b), = ureter_paths(lab, sp)
    assert np.allclose(b, np.array([25, 30, 20]) * 2, atol=1)


def curved_ureter():
    """Fantôme (voxels de 4 mm) : rein droit en haut, vessie en bas, uretère courbe fixant (SUV 6) qui bombe en dehors
    (x faibles) à mi-hauteur, à 48 mm de la droite rein -> vessie."""
    z, y, x = np.mgrid[:70, :60, :70]
    ball = lambda c, r: ((z - c[0]) ** 2 + (y - c[1]) ** 2 + (x - c[2]) ** 2) <= r ** 2
    lab = np.zeros(z.shape, np.uint8)
    lab[ball((60, 40, 20), 5)], lab[ball((8, 30, 35), 5)] = 6, 8
    suv = np.ones(z.shape, np.float32)
    suv[lab == 6], suv[lab == 8] = 10, 40
    tube = np.zeros(z.shape, bool)
    for k in range(10, 60):
        t = (k - 8) / 52
        tube[k] = (y[0] - 30 - 10 * t) ** 2 + (x[0] - 35 + 15 * t + 12 * np.sin(np.pi * t)) ** 2 <= 1.5 ** 2
    tube &= lab == 0
    suv[tube] = 6
    return suv, lab, tube, ball


def test_measure_at_excretion_follows_estimated_ureter_not_straight_line():
    sp = (4.0, 4.0, 4.0)
    suv, lab, tube, ball = curved_ureter()
    on, far, line = (34, 35, 16), (34, 25, 16), (34, 35, 27)   # sur l'uretère ; à 40 mm de lui ; sur l'ancienne droite
    suv[ball(on, 2)], suv[ball(far, 2)], suv[ball(line, 2)] = 8, 5, 5
    old = ureter_paths(lab, sp)                                  # pas d'uretère estimé : droite rein -> vessie
    assert len(old) == 1 and "excrétion possible" in measure_at(suv, sp, line, labels=lab, ureters=old)["flags"]
    lab[tube] = 13
    ure = ureter_paths(lab, sp)
    assert len(ure) > 10 and np.allclose(ure[0][0], np.array([60, 40, 20]) * 4) and np.allclose(ure[-1][1], np.array([8, 30, 35]) * 4)
    r = measure_at(suv, sp, on, labels=lab, ureters=ure)
    assert r["suv_max"] == 8 and "excrétion possible" in r["flags"]
    assert measure_at(suv, sp, far, labels=lab, ureters=ure)["flags"] == []
    assert measure_at(suv, sp, line, labels=lab, ureters=ure)["flags"] == []   # près de la droite, loin du trajet courbe
