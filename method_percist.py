"""PERCIST 1.0 sur les lésions SUIVIES (carte de référence mesurée à chaque date par `petscan.measure_at`), jamais sur un
MTV/TLG à seuil global (un seuil relatif au foie tronque les lésions peu avides et fait s'effondrer le MTV à tort).
Cible = lésion retenue la plus active à chaque date (SULpeak ; SULmax si la cible fait moins de 1 ml, le SULpeak de 1 cm³
la diluant). Lésions « excrétion possible » ou au contact d'un organe physiologique : exclues du calcul et signalées."""
from organs import ORGANS, PHYSIO
from store import date_key

PHYSIO_NOMS = {ORGANS[k][0] for k in PHYSIO}
LIBELLES = {"RMC": "Réponse métabolique complète", "RMP": "Réponse métabolique partielle", "MMS": "Maladie métabolique stable",
            "PMP": "Progression métabolique", "NE": "Non évaluable"}
EXCR, WIN = "excrétion possible", "volume limité à la fenêtre"
f2 = lambda v: f"{v:.2f}".replace(".", ",")


def _num(s):
    try:
        return float(str(s).split()[0].replace(",", "."))
    except (ValueError, IndexError):
        return None


def _why(v):
    """None : non mesurée ; "" : utilisable ; sinon la raison de l'exclusion."""
    if v is None:
        return None
    if EXCR in v["flags"]:
        return EXCR
    c = PHYSIO_NOMS & set(v["contact"])
    return "au contact : " + ", ".join(sorted(c)) if c else ""


def _acq(a, b, d):
    """Conditions d'acquisition lisibles dans les en-têtes : référence a, date d (b)."""
    if not a or not b:
        return [f"{d} : conditions d'acquisition non enregistrées, comparabilité non vérifiée"]
    out = []
    for k, tol in (("Activité injectée", 0.2), ("Poids", 0.1), ("Incorporation", None)):
        x, y = _num(a.get(k)), _num(b.get(k))
        if x is None or y is None:
            out.append(f"{d} : {k.lower()} illisible, comparabilité non vérifiée")
        elif abs(y - x) > (15 if tol is None else tol * x):  # PERCIST : incorporation à ± 15 min de la référence
            out.append(f"{d} : {k.lower()} {a[k]} → {b[k]} (écart > {'15 min' if tol is None else f'{tol * 100:.0f} %'})")
    if a.get("Machine") != b.get("Machine"):
        out.append(f"{d} : machine différente ({a.get('Machine') or '?'} → {b.get('Machine') or '?'})")
    return out


def _date(d, e, r0, base, valid, warns):
    bv, cv = base["values"], e["values"]
    use, excl = [], []
    for l in valid:
        if l.get("since", r0) != r0:
            continue
        vb, vc = bv.get(l["id"]), cv.get(l["id"])
        w = _why(vb) or _why(vc)  # exclue aux deux dates si signalée à l'une : comparaison homogène
        if w:
            excl.append(f"{l['name']} ({w})")
        elif vb:
            use.append((l["name"], vb, vc))
            if vc is None:
                warns.append(f"{d} : {l['name']} hors champ, non mesurée")
    news = []
    for l in valid:
        s = l.get("since", r0)
        if s != r0 and date_key(s) <= date_key(d):
            w = _why(cv.get(l["id"]))  # retenue par l'utilisateur : comptée même si non mesurée à cette date
            if w:
                excl.append(f"{l['name']} (nouvelle, {w})")
            else:
                news.append(l["name"])
    out = {"date": d, "provisoire": bool(e.get("provisoire")), "exclues": excl, "nouvelles": news}
    cur = [(n, vc) for n, _, vc in use if vc]
    if not use or not cur:
        return {**out, "statut": "NE", "libelle": LIBELLES["NE"],
                "texte": "Aucune lésion de référence utilisable" if not use else "Aucune lésion suivie mesurée à cette date"}
    pick = lambda xs, k: max(xs, key=lambda t: t[1][k])
    tb, tc = pick([(n, vb) for n, vb, _ in use], "sul_peak"), pick(cur, "sul_peak")
    key = "sul_peak" if min(tb[1]["vol_ml"], tc[1]["vol_ml"]) >= 1.0 else "sul_max"
    if key == "sul_max":
        tb, tc = pick([(n, vb) for n, vb, _ in use], key), pick(cur, key)
        warns.append(f"{d} : cible de moins de 1 ml, SULpeak non fiable : SULmax utilisé")
    b, c = tb[1][key], tc[1][key]
    dabs, dpct = c - b, ((c - b) / b * 100 if b > 0 else 0.0)
    lb, lv = base.get("liver"), e.get("liver")
    if lb and lv and abs(lv["mean"] - lb["mean"]) > max(0.2 * lb["mean"], 0.3):
        warns.append(f"{d} : référence foie variable ({f2(lb['mean'])} → {f2(lv['mean'])}), comparaison fragile")
    if not lv:
        warns.append(f"{d} : référence foie absente" + (", celle de la référence est utilisée" if lb else ""))
    if not lb:
        warns.append("Référence foie absente à la séance de référence : mesurabilité de la cible non vérifiée")
    elif b < 1.5 * lb["mean"] + 2 * lb["sd"]:
        warns.append(f"Cible de référence sous le seuil de mesurabilité PERCIST ({f2(1.5 * lb['mean'] + 2 * lb['sd'])})")
    if any(WIN in v["flags"] for _, vb, vc in use for v in (vb, vc) if v):
        warns.append(f"{d} : volume « limité à la fenêtre » pour au moins une lésion (croissance non évaluée pour elle)")
    grow = [n for n, vb, vc in use if vc and WIN not in vb["flags"] + vc["flags"]
            and vc["vol_ml"] >= 1.75 * vb["vol_ml"] and vc["vol_ml"] - vb["vol_ml"] >= 1.0]
    liver = lv or lb
    if news or grow or (dpct >= 30 and dabs >= 0.8):
        st = "PMP"
    elif liver and max(v[key] for _, v in cur) <= liver["mean"]:
        st = "RMC"
    elif dpct <= -30 and dabs <= -0.8:
        st = "RMP"
    else:
        st = "MMS"
        if not liver:
            warns.append(f"{d} : sans référence foie, la réponse complète n'est pas évaluable")
    unit = "SULpeak" if key == "sul_peak" else "SULmax"
    txt = f"Cible {unit} : {tb[0]} {f2(b)} → {tc[0]} {f2(c)} ({dpct:+.0f} %, {'+' if dabs >= 0 else ''}{f2(dabs)})"
    if news:
        txt += " ; nouvelle(s) lésion(s) retenue(s) : " + ", ".join(news)
    if grow:
        txt += " ; croissance nette (volume ≥ +75 % et ≥ +1 ml) : " + ", ".join(grow)
    return {**out, "statut": st, "libelle": LIBELLES[st], "metrique": unit, "reference": [tb[0], b], "cible": [tc[0], c],
            "variation_pct": dpct, "variation_abs": dabs, "croissance": grow, "texte": txt,
            "_t": {k: max(v[k] for _, v in cur) for k in ("sul_peak", "sul_max")}}


def _nadir(x, hist, warns):
    """Nadir de la cible (plus basse valeur depuis la référence incluse, dans la métrique de cette date) et remontée
    ≥ 30 % et ≥ 0,8 depuis ce nadir : signalée, sans changer le statut officiel (calculé vs référence)."""
    k = "sul_peak" if x["metrique"] == "SULpeak" else "sul_max"
    nd, nv = min([(hist[0][0], x["reference"][1])] + [(d, t[k]) for d, t in hist[1:]], key=lambda t: t[1])
    c = x["cible"][1]
    dabs, dpct = c - nv, ((c - nv) / nv * 100 if nv > 0 else 0.0)
    up = dpct >= 30 and dabs >= 0.8
    x.update(nadir=[nd, nv], variation_nadir_pct=dpct, variation_nadir_abs=dabs, remontee_nadir=up)
    x["texte"] += f" ; vs nadir ({nd}, {f2(nv)}) : {dpct:+.0f} %, {'+' if dabs >= 0 else ''}{f2(dabs)}"
    if up:
        x["texte"] += " : remontée depuis le nadir"
        warns.append(f"{x['date']} : remontée depuis le nadir ({nd}, {f2(nv)} → {f2(c)}, {dpct:+.0f} %) ; "
                     "statut PERCIST officiel inchangé (calculé vs référence)")


def compare(ref, dates, ctx):
    """ref : {date, lesions, unit} ; dates : {date: {values, acq, liver, unit, provisoire}} (ordre chronologique)."""
    r0, warns = ref["date"], []
    base = dates.get(r0)
    if base is None:
        return {"statut": "NE", "libelle": LIBELLES["NE"], "details": [], "avertissements": ["Mesures de la référence absentes"]}
    valid = [l for l in ref["lesions"] if l.get("status") == "valide"]
    if any((e.get("unit") or ref.get("unit")) != "SUL" for e in dates.values()):
        warns.append("SUL de repli (SUV poids) pour au moins une date : taille ou sexe absents des en-têtes")
    details, hist = [], [(r0, None)]
    for d, e in dates.items():
        if d != r0:
            warns += _acq(base.get("acq"), e.get("acq"), d)
            x = _date(d, e, r0, base, valid, warns)
            t = x.pop("_t", None)
            if t and date_key(d) > date_key(r0):  # nadir : depuis la référence incluse
                _nadir(x, hist, warns)
                hist.append((d, t))
            details.append(x)
    if not details:
        warns.append("Aucune date de suivi : appliquez le suivi à une autre séance")
    last = details[-1] if details else {"statut": "NE"}
    pc = lambda x, k: f"{x[k]:+.0f} %" if k in x else "—"
    tab = {"colonnes": ["Cible", r0] + [x["date"] for x in details],
           "lignes": [["Valeur (SULpeak / SULmax)", next((f2(x["reference"][1]) for x in details if "reference" in x), "—")] +[f"{x['metrique']} {f2(x['cible'][1])}" if "cible" in x else "—" for x in details],
                      ["% vs référence", "—"] + [pc(x, "variation_pct") for x in details],
                      ["Nadir", "—"] + [f"{f2(x['nadir'][1])} ({x['nadir'][0]})" if "nadir" in x else "—" for x in details],
                      ["% vs nadir", "—"] + [pc(x, "variation_nadir_pct") + (" ⚠ remontée" if x.get("remontee_nadir") else "")
                                             for x in details],
                      ["Statut", "référence"] + [x["statut"] for x in details]]}
    return {"statut": last["statut"], "libelle": LIBELLES[last["statut"]], "details": details, "tableau": tab,
            "avertissements": list(dict.fromkeys(warns))}


METHOD = dict(id="percist", nom="PERCIST 1.0", donnees=("TEP",), cancers=("dsrct", "general"), compare=compare,
              defaut=("dsrct", "general"), note="lésions suivies, SULpeak de la lésion la plus active")
