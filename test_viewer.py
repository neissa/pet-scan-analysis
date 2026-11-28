import pytest

import http.client
import json
import threading
import time
from http.server import ThreadingHTTPServer

from viewer import State, demo, make_handler, measure, render

PNG = b"\x89PNG"


@pytest.mark.parametrize("view", ["axial", "coronal", "sagittal", "mip"])
def test_render_png(view):
    suv, *_ = demo()
    assert render(suv, view, 10, 10).startswith(PNG)


def test_render_anisotropic_spacing():
    suv, *_ = demo()
    assert render(suv, "coronal", 10, 10, spacing=(3, 1, 1)).startswith(PNG)


@pytest.mark.parametrize("mode", ["pet", "ct", "fusion"])
@pytest.mark.parametrize("view", ["axial", "coronal", "mip"])
def test_render_fusion_modes(view, mode):
    suv, _, ctv = demo()
    assert render(suv, view, 30, 10, ctv=ctv, mode=mode).startswith(PNG)


def test_render_zoom_box_keeps_size_and_png():
    suv, *_ = demo()
    full = render(suv, "axial", 30, 10)
    zoomed = render(suv, "axial", 30, 10, box=(0.4, 0.3, 0.25))
    assert zoomed.startswith(PNG) and zoomed != full
    assert zoomed[16:24] == full[16:24]  # même largeur/hauteur de sortie (en-tête IHDR)


def test_measure_demo_lesion():
    suv, spacing, _ = demo()  # sphère de rayon 5 voxels (3 mm) en (30, 30, 30)
    r = measure(suv, spacing, "axial", 30, 30.5 / 80, 30.5 / 80, 50)
    assert 28 <= r["length_mm"] <= 36 and 28 <= r["width_mm"] <= 36 and r["suv_max"] == 9


def test_measure_ct_and_missing_ct():
    suv, spacing, ctv = demo()  # CT synthétique : ellipsoïde de 40 HU dans l'air (-1000)
    r = measure(suv, spacing, "axial", 30, 0.5, 0.5, 50, ctv, "ct", 50)
    assert r["hu_ref"] == 40 and r["length_mm"] > 100 and "error" not in r
    assert measure(suv, spacing, "axial", 30, 0.5, 0.5, 50, None, "ct")["error"]


def test_render_index_clipped():
    suv, *_ = demo()
    assert render(suv, "axial", 10_000, 10).startswith(PNG)


def _server(state):
    srv = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(state))
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv


def _get(srv, path, host=None):
    c = http.client.HTTPConnection("127.0.0.1", srv.server_port)
    c.request("GET", path, headers={"Host": host or f"127.0.0.1:{srv.server_port}"})
    r = c.getresponse()
    return r.status, r.read()


def test_routes_before_and_after_loading_and_host_check():
    state = State()
    srv = _server(state)
    try:
        assert json.loads(_get(srv, "/info")[1])["loaded"] is False
        assert _get(srv, "/img?view=axial&i=0&vmax=10")[0] == 503
        assert _get(srv, "/lesions?thr=2.5")[0] == 503
        state.set(*demo())
        assert json.loads(_get(srv, "/info")[1])["loaded"] is True
        assert _get(srv, "/img?view=axial&i=3&vmax=10")[0] == 200
        assert _get(srv, "/info", host="evil.example")[0] == 403  # anti DNS-rebinding
        assert _get(srv, "/open?folder=/chemin/inexistant")[0] == 400
    finally:
        srv.shutdown()


def test_invalid_folder_keeps_current_volume(tmp_path):
    state = State()
    state.set(*demo())
    before = state.vol
    assert state.open(str(tmp_path))  # dossier sans DICOM
    for _ in range(100):
        if state.status != "loading":
            break
        time.sleep(0.05)
    assert state.status == "error" and "PET" in state.message
    assert state.vol is before  # l'ancien volume reste affiché


def test_render_with_organs_overlay_and_routes():
    import numpy as np
    suv, spacing, ctv = demo()
    lab = np.zeros(suv.shape, np.uint8)
    lab[20:40, 20:40, 20:40] = 3
    plain, tinted = render(suv, "axial", 30, 10, ctv=ctv), render(suv, "axial", 30, 10, ctv=ctv, labels=lab)
    assert tinted.startswith(PNG) and tinted != plain
    assert render(suv, "mip", 0, 10, ctv=ctv, labels=lab) == render(suv, "mip", 0, 10, ctv=ctv)  # pas de surcouche sur la MIP
    state = State()
    state.set(suv, spacing, ctv)
    srv = _server(state)
    try:
        assert json.loads(_get(srv, "/organs")[1])["status"] == "computing"
        state.organs = (lab, {3: {"nom": "Cœur (approx.)", "couleur": "#ff4a4a", "volume_ml": 1.0, "suv_max": 1.0, "suv_mean": 1.0, "hu_mean": 1.0}})
        assert json.loads(_get(srv, "/organs")[1])["organes"][0]["nom"] == "Cœur (approx.)"
        assert json.loads(_get(srv, "/organ_at?view=axial&i=30&fx=0.4&fy=0.4")[1]) == {"nom": "Cœur (approx.)"}
        assert json.loads(_get(srv, "/organ_at?view=axial&i=30&fx=0.02&fy=0.02")[1]) == {}
        assert _get(srv, "/img?view=axial&i=30&vmax=10&organs=1")[0] == 200
    finally:
        srv.shutdown()


def test_rotating_mip_frames_routes_and_measure():
    suv, spacing, ctv = demo()
    state = State()
    state.set(suv, spacing, ctv)
    for _ in range(200):  # calcul en arrière-plan
        if state.mips is not None and len(state.mips) == 35:
            break
        time.sleep(0.05)
    assert len(state.mips) == 35 and state.mip_frame(0) is None and state.mip_frame(94).shape == (60, 40)
    srv = _server(state)
    try:
        assert json.loads(_get(srv, "/mips")[1]) == {"ready": 36, "total": 36}
        assert _get(srv, "/img?view=mip&i=90&vmax=10")[0] == 200 and _get(srv, "/img?view=mip&i=0&vmax=10")[0] == 200
        r = json.loads(_get(srv, "/measure?view=mip&i=90&fx=0.5&fy=0.5&pct=50")[1])
        assert "error" not in r and r["length_mm"] > 0
        assert _get(srv, "/ctseries?uid=inconnue")[0] == 400
    finally:
        srv.shutdown()


def test_compare_modes_lesions_status_and_measure():
    suv, spacing, ctv = demo()
    state = State()
    state.set(suv, spacing, ctv)
    assert json.loads(_get(_server(state), "/cstatus")[1])["ready"] is False
    state.cmp = {"suv": suv * 0.5, "ctv": ctv, "patient": {"Date de l'examen": "04/09/2026"}, "folder": "b",
                 "shift": (1, 0, 0), "score": 0.7}  # séance B : moitié moins de captation, déjà recalée
    srv = _server(state)
    try:
        imgs = {m: _get(srv, f"/img?view=axial&i=30&vmax=10&mode=fusion&cmp={m}&split=0.5")[1] for m in ("a", "b", "swipe", "diff")}
        assert all(v.startswith(PNG) for v in imgs.values()) and len(set(imgs.values())) == 4  # 4 rendus différents
        assert _get(srv, "/img?view=mip&i=90&vmax=10&cmp=diff")[0] == 200   # MIP de face en comparaison
        assert _get(srv, "/img?view=axial&i=30&vmax=10&cmp=nimporte")[0] == 400
        les = json.loads(_get(srv, "/lesions?thr=4")[1])
        assert les[0]["suv_b"] == 4.5 and abs(les[0]["delta_pct"] - 100) < 1e-6     # 9 contre 4,5 : +100 %
        cs = json.loads(_get(srv, "/cstatus")[1])
        assert cs["ready"] and cs["date_b"] == "04/09/2026" and cs["shift_mm"] == [3.0, 0.0, 0.0] and cs["score"] == 0.7
        r = json.loads(_get(srv, "/measure?view=axial&i=30&fx=0.5&fy=0.5&pct=50&cmp=b")[1])
        assert r["suv_max"] <= 4.6                                                  # mesurée sur B, pas sur A
        assert _get(srv, "/compare?folder=/dossier/inexistant")[0] == 400
    finally:
        srv.shutdown()


def test_manual_offsets_slice_semantics_lesions_and_routes():
    import numpy as np
    from compare import shift_volume
    from viewer import _add_b, _slice_off
    rng = np.random.default_rng(1)
    vol = rng.random((12, 20, 24)).astype(np.float32)
    off = (2, -3, 4)
    for view, i in (("axial", 5), ("coronal", 9), ("sagittal", 11)):  # même convention que shift_volume, sans copier le volume
        ref = shift_volume(vol, off, vol.shape, 0.0)
        from viewer import _slice
        assert np.array_equal(_slice_off(vol, view, i, off), _slice(ref, view, i)), view
    off_mip = (2, 0, 4)  # le MIP se projette selon y : le décalage en y ne change que les bords, il est ignoré
    assert np.array_equal(_slice_off(vol, "mip", 0, off_mip), np.max(shift_volume(vol, off_mip, vol.shape, 0.0), axis=1)[::-1])
    assert (_slice_off(vol, "axial", 11, (5, 0, 0)) == 0).all()  # coupe décalée hors du volume : remplie
    # lésions : le point de B correspondant au pic de A est p + rel
    b = np.zeros((30, 30, 30), np.float32)
    b[20, 15, 15] = 8
    les = [{"peak": [10, 15, 15], "suv_max": 4.0}]
    assert _add_b([dict(l) for l in les], b, (2, 2, 2), radius_mm=2)[0]["suv_b"] == 0          # sans décalage : rien en B
    assert _add_b([dict(l) for l in les], b, (2, 2, 2), radius_mm=2, rel=(10, 0, 0))[0]["suv_b"] == 8
    assert _add_b([dict(l) for l in les], b, (2, 2, 2), radius_mm=2, rel=(-40, 0, 0))[0]["suv_b"] == 0  # fenêtre hors volume
    suv, spacing, ctv = demo()
    state = State()
    state.set(suv, spacing, ctv)
    srv = _server(state)
    try:
        a = _get(srv, "/img?view=axial&i=30&vmax=10&mode=pet")[1]
        assert _get(srv, "/img?view=axial&i=30&vmax=10&mode=pet&oa=0,5,0")[1] != a          # A décalée : image différente
        assert _get(srv, "/img?view=axial&i=30&vmax=10&mode=pet&ob=0,5,0")[1] == a          # B seule ne change pas A
        assert _get(srv, "/img?view=axial&i=30&vmax=10&oa=1,2")[0] == 400                    # décalage mal formé
        assert _get(srv, "/img?view=axial&i=30&vmax=10&oa=999,0,0")[0] == 400                # hors bornes
        r = json.loads(_get(srv, "/measure?view=axial&i=30&fx=0.5&fy=0.5&pct=50&oa=0,0,0")[1])
        assert "error" not in r
    finally:
        srv.shutdown()


def test_lesion_organ_location_dans_pres_de_loin_et_absence():
    import numpy as np
    from viewer import _organ_at_peak
    spacing = (3.0, 2.0, 2.0)  # z anisotrope : la distance doit être calculée en mm
    suv = np.zeros((40, 40, 40), np.float32)
    lab = np.zeros(suv.shape, np.uint8)
    peaks = {"dans": (10, 10, 10), "pres": (20, 10, 10), "loin": (30, 30, 30), "bord": (5, 30, 30)}
    for p in peaks.values():
        suv[p] = 5.0
    lab[9:12, 9:12, 9:12] = 4      # Foie : contient le pic « dans »
    lab[20, 10, 16] = 6            # Rein droit : 6 voxels en x = 12 mm du pic « près de » (hors de la marge retirée de 10 mm)
    lab[20, 10, 11] = 10           # squelette à 2 mm : ignoré (trop répandu)
    lab[31, 30, 30] = 10           # squelette à 3 mm du pic « loin » : ignoré
    lab[30, 30, 38] = 9            # Cerveau à 8 voxels en x = 16 mm : trop loin
    lab[10, 30, 30] = 5            # Rate à 5 voxels en z = 15 mm : limite incluse
    assert _organ_at_peak(lab, spacing, peaks["dans"]) == "dans Foie"
    assert _organ_at_peak(lab, spacing, peaks["pres"]) == "près de Rein droit"
    assert _organ_at_peak(lab, spacing, peaks["loin"]) is None
    assert _organ_at_peak(lab, spacing, peaks["bord"]) == "près de Rate"
    state = State()
    state.set(suv, spacing, None)
    srv = _server(state)
    try:
        les = json.loads(_get(srv, "/lesions?thr=2.5")[1])
        assert len(les) == 4 and all("organe" not in l for l in les)      # organes pas encore calculés : champ absent
        state.organs = (lab, {})
        state.cmp = {"suv": suv * 0.5, "ctv": None, "patient": {}, "folder": "b", "shift": (0, 0, 0), "score": None}
        les = {tuple(l["peak"]): l for l in json.loads(_get(srv, "/lesions?thr=2.5")[1])}
        assert les[peaks["dans"]]["organe"] == "dans Foie"
        assert les[peaks["pres"]]["organe"] == "près de Rein droit"
        assert les[peaks["loin"]]["organe"] is None
        assert les[peaks["bord"]]["organe"] == "près de Rate"
        assert les[peaks["pres"]]["suv_max"] == 5.0 and les[peaks["pres"]]["suv_b"] == 2.5  # champs existants intacts
    finally:
        srv.shutdown()


def test_percist_route_sul_fallback_and_lesion_fields():
    from petscan import sul_info
    from test_petscan import _phantom
    suv, sp, lab = _phantom()
    state = State()
    state.set(*demo())
    srv = _server(state)
    try:
        p = json.loads(_get(srv, "/percist")[1])   # démo : pas de taille ni de sexe -> repli SUV poids, signalé
        assert p["sul"] is False and p["factor"] == 1.0 and p["ref"] is None and p["computing"] is True
        state.set(suv, sp, None, sul=sul_info(80, 1.80, "M"))
        p = json.loads(_get(srv, "/percist")[1])
        assert p["sul"] is True and p["ref"] is None and "pas de CT" in p["ref_msg"]
        les = json.loads(_get(srv, "/lesions?thr=2.5")[1])
        assert all("sul_peak" in l and "physio" not in l for l in les)    # organes pas encore calculés
        assert "error" in json.loads(_get(srv, "/lesions?thr=percist")[1])  # pas de référence foie : seuil PERCIST refusé
        state.organs = (lab, {})
        ref = json.loads(_get(srv, "/percist")[1])["ref"]
        assert ref["seuil_suv"] == pytest.approx(ref["seuil"] / state.sul["factor"])   # SUL -> SUV
        lp = json.loads(_get(srv, "/lesions?thr=percist")[1])
        assert lp == json.loads(_get(srv, f"/lesions?thr={ref['seuil_suv']!r}")[1]) and min(l["suv_max"] for l in lp) >= ref["seuil_suv"]
        p = json.loads(_get(srv, "/percist")[1])
        assert p["ref"]["mean"] == pytest.approx(2.0 * p["factor"]) and p["ref_msg"] == ""
        assert p["pool"] is None                                          # pas d'aorte segmentée
        state.organs = ((lab + 7 * (lab == 4)).astype(lab.dtype), {})   # le foie (4) tenu pour l'aorte (11)
        assert json.loads(_get(srv, "/percist")[1])["pool"]["mean"] == pytest.approx(2.0 * p["factor"], rel=0.05)
        les = json.loads(_get(srv, "/lesions?thr=2.5")[1])
        assert isinstance(les, list) and next(l for l in les if l["suv_max"] == 5.0)["physio"] == "Vessie"
        assert les[0]["sul_max"] == pytest.approx(les[0]["suv_max"] * p["factor"])
    finally:
        srv.shutdown()


def test_autoalign_route_global_local_and_errors():
    import numpy as np
    from compare import shift_volume
    from test_compare import phantom
    ct = phantom()
    suv = np.ones(ct.shape, np.float32)
    true_e = (2, -4, 6)  # B décalée de true_e : b[x + e] = a[x]
    b_ct = shift_volume(ct, tuple(-v for v in true_e), ct.shape, -1000)
    state = State()
    state.set(suv, (4.0, 4.0, 4.0), ct)
    srv = _server(state)
    try:
        assert "error" in json.loads(_get(srv, "/autoalign?scope=global")[1])        # pas de séance B
        state.cmp = {"suv": suv, "ctv": b_ct, "patient": {}, "folder": "b", "shift": (0, 0, 0), "score": None}
        g = json.loads(_get(srv, "/autoalign?scope=global")[1])
        assert all(abs(x - w) <= 3 for x, w in zip(g["e"], true_e)) and g["score"] > 0.5, g   # pas de 3 voxels
        loc = json.loads(_get(srv, "/autoalign?scope=local&c=40,45,45")[1])
        assert loc["e"] == list(true_e) and loc["ob"] == list(true_e), loc                  # au voxel près
        # en repartant d'un décalage déjà correct : résidu nul
        again = json.loads(_get(srv, "/autoalign?scope=local&c=40,45,45&ob=2,-4,6")[1])
        assert again["e"] == [0, 0, 0] and again["ob"] == [2, -4, 6]
        assert _get(srv, "/autoalign?scope=local&c=1,2")[0] == 400                          # centre mal formé
        assert _get(srv, "/autoalign?scope=local&c=500,45,45")[0] == 200                    # coupe > 400 : un centre n'est pas un décalage
    finally:
        srv.shutdown()


def test_progress_resets_on_open_and_after_load(tmp_path):
    state = State()
    state.progress = 1.0  # reste du chargement précédent
    assert state.open(str(tmp_path))  # dossier vide : le chargement échoue vite
    assert state.progress is None  # la barre ne doit pas repartir à 100 %
    for _ in range(100):
        if state.status != "loading":
            break
        time.sleep(0.05)
    assert state.status == "error" and state.progress is None


def test_android_env_index_and_config(tmp_path, monkeypatch):
    import os
    import subprocess
    import sys
    page = tmp_path / "index.html"
    page.write_bytes(b"<p>page android</p>")
    monkeypatch.setenv("PETSCAN_INDEX", str(page))
    srv = _server(State())
    try:
        assert _get(srv, "/") == (200, b"<p>page android</p>")  # PETSCAN_INDEX prioritaire
    finally:
        srv.shutdown()
    cfg = tmp_path / "cfg.json"  # PETSCAN_CONFIG : lu à l'import, d'où un processus à part ; ancien fichier migré puis supprimé
    cfg.write_text('{"recent": ["/old"]}', encoding="utf-8")
    code = ("import store, viewer; print(viewer.recent()); store.set_config(recent_on=True); viewer.remember('/d'); "
            "print(viewer.recent(), store.DIR.name)")
    out = subprocess.run([sys.executable, "-c", code], env={**os.environ, "PETSCAN_CONFIG": str(cfg)}, capture_output=True, text=True)
    assert out.stdout.split("\n")[:2] == ["[]", "['/d', '/old'] netscan"], out.stderr
    assert not cfg.exists() and (tmp_path / "netscan" / "config.json").is_file()


def _post(srv, path, obj, site=None):
    c = http.client.HTTPConnection("127.0.0.1", srv.server_port)
    h = {"Host": f"127.0.0.1:{srv.server_port}", "Content-Type": "application/json", **({"Sec-Fetch-Site": site} if site else {})}
    c.request("POST", path, body=json.dumps(obj), headers=h)
    r = c.getresponse()
    return json.loads(r.read()) if r.status != 403 else 403


def _session(d, lesion_suv, extra=None):
    """Séance synthétique décalée de d voxels (b[x + d] = a[x]) : CT, organes (vessie pleine d'urine), lésions."""
    import numpy as np
    from compare import shift_volume
    from test_compare import phantom
    ct = phantom()
    z, y, x = np.mgrid[:80, :90, :90]
    ball = lambda c, r: (z - c[0]) ** 2 + (y - c[1]) ** 2 + (x - c[2]) ** 2 <= r * r
    suv = np.ones(ct.shape, np.float32)
    lab = np.zeros(ct.shape, np.uint8)
    lab[ball((20, 45, 55), 5)] = 8                     # vessie, urine SUV 50
    suv[lab == 8] = 50
    suv[ball((40, 45, 30), 3)] = 0.8 * lesion_suv     # L1 : pic au centre
    suv[40, 45, 30] = lesion_suv
    suv[ball((20, 45, 62), 2)] = 5.0                   # L2 collée à la vessie
    suv[20, 45, 62] = 5.5
    if extra:
        suv[ball(extra, 2)] = 6.0                      # nouvelle fixation
    back = tuple(-v for v in d)
    return (shift_volume(suv, back, ct.shape, 1.0), shift_volume(ct, back, ct.shape, -1000), shift_volume(lab, back, ct.shape, 0))


def test_reference_map_follow_three_dates_excretion_new_fixation_and_deletion():
    import store
    sp = (2.0, 2.0, 2.0)
    pat = lambda date: {"Nom": "Fantôme TEST", "Identifiant": "PH-0001", "Naissance": "01/01/1950", "Date de l'examen": date}
    state = State()
    srv = _server(state)
    try:
        suv, ct, lab = _session((0, 0, 0), 8.0)
        state.set(suv, sp, ct, pat("01/04/2026"), "refA")
        assert "error" in _post(srv, "/ref/define", {"lesions": [{"vox": [40, 45, 30]}]})   # organes pas encore calculés
        state.organs = (lab, {})
        assert json.loads(_get(srv, "/refs")[1]) == {"current": None, "refs": []}
        key = _post(srv, "/ref/define", {"lesions": [{"vox": [40, 45, 30], "name": "rate"}, {"vox": [20, 45, 62]}]})["key"]
        assert _post(srv, "/ref/define", {"lesions": []}) == {"exists": True}      # pas d'écrasement sans confirmation
        assert json.loads(_get(srv, "/refs")[1])["current"] == key
        ct_ref = ct
        for d, date, v, extra in (((3, -4, 5), "01/07/2026", 6.0, (60, 30, 45)), ((-2, 3, -4), "01/09/2026", 4.0, None)):
            suv, ct, lab = _session(d, v, extra)
            state.set(suv, sp, ct, pat(date), "autre")
            state.organs, state._refct = (lab, {}), ("refA", ct_ref)   # CT de référence déjà relu
            f = json.loads(_get(srv, f"/follow?key={key}&thr=3")[1])
            l1, l2 = f["rows"]["L1"], f["rows"]["L2"]
            assert l1["peak"] == [40 + d[0], 45 + d[1], 30 + d[2]] and l1["suv_max"] == v, f        # recalage, lésion retrouvée
            assert l1["flags"] == [] and "excrétion possible" in l2["flags"] and l2["suv_max"] == 5.5  # urine écartée, lésion gardée
            assert f["score"] > 0.5 and f["method"].startswith("squelette")
            news = [n["peak"] for n in f["new"]]   # boule uniforme de rayon 2 : pic n'importe où dedans
            assert len(news) == bool(extra) and all(max(abs(a - b - c) for a, b, c in zip(n, extra, d)) <= 2 for n in news), news
            assert _post(srv, "/ref/record", {"key": key}) == {"ok": True}
        assert _post(srv, "/ref/update", {"key": key, "op": "add", "vox": [40, 45, 30]})["ok"]   # ajout au point courant
        assert _post(srv, "/ref/update", {"key": key, "op": "name", "id": "L2", "value": "vessie ?"})["ok"]
        m = json.loads(_get(srv, f"/ref?key={key}")[1])
        assert m["dates_sorted"] == ["01/04/2026", "01/07/2026", "01/09/2026"] and "folder" not in m and m["folder_name"] == "refA"
        assert [m["dates"][d]["values"]["L1"]["suv_max"] for d in m["dates_sorted"]] == [8.0, 6.0, 4.0]
        assert [l["name"] for l in m["lesions"]] == ["rate", "vessie ?", "L3"] and m["lesions"][2]["vox"] == [42, 42, 34]
        raw = (store.DIR / "patients" / f"{key}.json").read_text(encoding="utf-8")
        assert "PH-0001" not in raw and "Fantôme" not in raw and "1950" not in raw
        assert _get(srv, f"/data/delete?what=patient&key={key}")[0] == 400                      # confirmation exigée
        assert _post(srv, "/ref/record", {"key": key}, site="cross-site") == 403
        c = http.client.HTTPConnection("127.0.0.1", srv.server_port)
        c.request("GET", f"/data/delete?what=all&confirm=1", headers={"Host": f"127.0.0.1:{srv.server_port}", "Sec-Fetch-Site": "cross-site"})
        assert c.getresponse().status == 403                                                      # route qui change l'état
        assert json.loads(_get(srv, "/data")[1])["patients"] == 1
        assert json.loads(_get(srv, f"/data/delete?what=patient&key={key}&confirm=1")[1]) == {"deleted": 1}
        assert json.loads(_get(srv, "/data")[1])["patients"] == 0
    finally:
        srv.shutdown()
