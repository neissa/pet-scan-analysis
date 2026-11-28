"""RECIST 1.1 sur les lésions TEP SUIVIES : cibles choisies à la référence parmi les lésions retenues, diamètres mesurés au CT
de la TEP-CT à la position de chaque lésion (`petscan.ct_diameter`, enregistrés par date sous "ct"). Aucun diamètre n'est
inventé : CT absent, mesure non fiable ou date enregistrée sans diamètre = lésion non mesurable ou date non évaluable.
Lésions « excrétion possible » ou au contact d'un organe physiologique (à une date au moins) : exclues, comme pour PERCIST."""
from method_percist import _why
from store import date_key

LIBELLES = {"RC": "Réponse complète", "RP": "Réponse partielle", "MS": "Maladie stable", "PE": "Progression", "NE": "Non évaluable"}
MAX_CIBLES, PAR_ORGANE = 5, 2
f1 = lambda v: f"{v:.1f}".replace(".", ",")
pc = lambda v: "—" if v is None else f"{v:+.0f} %"


def _size(v, node):
    """(diamètre RECIST en mm ou None, raison) d'une mesure enregistrée ; ganglion : petit axe, sinon plus grand diamètre."""
    if v is None:
        return None, "hors champ"
    if "ct" not in v:
        return None, "diamètre non enregistré"
    ct = v["ct"]
    if ct is None:
        return None, "pas de CT"
    if "why" in ct:
        return None, ct["why"]
    return (ct["short_mm"] if node else ct["long_mm"]), ""


def _organ(l):
    return "Ganglions" if l.get("ganglion") else (l.get("organe") or "organe non déterminé")


def select(ref, dates, lim=10.0):
    """Cibles à la référence (≤ 5, ≤ 2 par organe, plus grandes d'abord). Mesurable : diamètre ≥ lim (10 mm, ou 2 × l'épaisseur
    de coupe si elle dépasse 5 mm) ; ganglion : petit axe ≥ 15 mm (10 à 15 mm : pathologique non cible, < 10 mm : normal).
    Retourne (cibles, non_cibles, exclues, non_mesurables), les trois dernières en libellés."""
    r0 = ref["date"]
    base = dates[r0]["values"]
    cands, nonc, excl, nm = [], [], [], []
    for l in ref["lesions"]:
        if l.get("status") != "valide" or l.get("since", r0) != r0:
            continue
        w = next((x for x in (_why(e["values"].get(l["id"])) for e in dates.values()) if x), "")
        if w:
            excl.append(f"{l['name']} ({w})")
            continue
        node = bool(l.get("ganglion"))
        s, why = _size(base.get(l["id"]), node)
        if s is None:
            nm.append(f"{l['name']} ({why})")
        elif node and s < 10:
            nonc.append(f"{l['name']} (ganglion non pathologique, petit axe {f1(s)} mm)")
        elif node and s < 15:
            nonc.append(f"{l['name']} (ganglion, petit axe {f1(s)} mm < 15 mm)")
        elif not node and s < lim:
            nonc.append(f"{l['name']} ({f1(s)} mm < {f1(lim)} mm)")
        else:
            cands.append((s, l))
    cib, per = [], {}
    for s, l in sorted(cands, key=lambda t: -t[0]):
        o = _organ(l)
        if len(cib) < MAX_CIBLES and per.get(o, 0) < PAR_ORGANE:
            cib.append(l)
            per[o] = per.get(o, 0) + 1
        else:
            nonc.append(f"{l['name']} ({'5 cibles déjà retenues' if len(cib) >= MAX_CIBLES else f'2 cibles déjà retenues : {o}'})")
    return cib, nonc, excl, nm


def compare(ref, dates, ctx):
    """ref : {date, lesions, unit, spacing} ; dates : {date: {values (dont "ct"), score, pe_non_cible, provisoire}} (ordre chronologique)."""
    r0 = ref["date"]
    if r0 not in dates:
        return {"statut": "NE", "libelle": LIBELLES["NE"], "details": [], "avertissements": ["Mesures de la référence absentes"]}
    warns = ["Diamètres mesurés sur le CT de la TEP-CT (sans injection, basse dose) : non diagnostiques pour RECIST 1.1, "
             "à confirmer sur un scanner diagnostique",
             "Diamètres par tolérance HU (± 50 HU) dans le plan axial : approximatifs, à vérifier sur les images"]
    th = (ref.get("spacing") or [None])[0]
    lim = max(10.0, 2 * th) if th and th > 5 else 10.0
    warns.append(f"Coupes de {f1(th)} mm (CT rééchantillonné sur la grille TEP)" + (f" : lésion mesurable si ≥ {f1(lim)} mm" if th > 5 else "")
                 if th else "Épaisseur de coupe inconnue")
    cib, nonc, excl, nm = select(ref, dates, lim)
    warns += [f"Lésion non mesurable à la référence : {x}" for x in nm]
    if nonc:
        warns.append("Non-cibles (progression « sans équivoque » à signaler par le lecteur) : " + ", ".join(nonc))
    order = list(dates)
    sizes, sums = {}, {}
    for d in order:
        sizes[d] = [_size(dates[d]["values"].get(l["id"]), bool(l.get("ganglion"))) for l in cib]
        miss = [f"{l['name']} ({w})" for l, (s, w) in zip(cib, sizes[d]) if s is None]
        sums[d] = None if miss or not cib else sum(s for s, _ in sizes[d])
        if miss:
            warns.append(f"{d} : diamètre indisponible pour " + ", ".join(miss))
        sc = dates[d].get("score", 1.0)
        if d != r0 and "score" in dates[d] and (sc is None or sc < 0.5):
            warns.append(f"{d} : recalage " + ("sans CT (coordonnées DICOM)" if sc is None else f"peu fiable (corrélation {sc:.2f})".replace(".", ",")))
    details = []
    for i, d in enumerate(order):
        if d == r0 or date_key(d) < date_key(r0):
            continue
        e, s, s0 = dates[d], sums[d], sums[r0]
        prior = [sums[x] for x in order[:i] if sums[x] is not None and date_key(x) >= date_key(r0)]
        nadir = min(prior) if prior else None
        ex, news = list(excl), []
        for l in ref["lesions"]:
            since = l.get("since", r0)
            if l.get("status") == "valide" and since != r0 and date_key(since) <= date_key(d):
                w = _why(e["values"].get(l["id"]))  # retenue par l'utilisateur : comptée même si non mesurée à cette date
                (ex.append(f"{l['name']} (nouvelle, {w})") if w else news.append(l["name"]))
        pe_nc = bool(e.get("pe_non_cible"))
        vr = (s - s0) / s0 * 100 if s is not None and s0 else None
        vn = (s - nadir) / nadir * 100 if s is not None and nadir else None
        grow = s is not None and nadir is not None and s - nadir >= 5 and (nadir == 0 or vn >= 20)
        gone = s is not None and all(v < (10 if l.get("ganglion") else 5) for l, (v, _) in zip(cib, sizes[d]))
        if not cib:
            st, txt = ("PE" if news or pe_nc else "NE"), "Aucune cible mesurable à la référence"
        else:
            st = ("PE" if news or pe_nc or grow else "NE" if s is None else "RC" if gone else "RP" if vr <= -30 else "MS")
            txt = ("Somme des diamètres non calculable (diamètre indisponible pour au moins une cible)" if s is None else
                   f"Somme des diamètres : {f1(s)} mm (référence {f1(s0)} mm, {pc(vr)} ; nadir {f1(nadir)} mm, {pc(vn)}"
                   + (f", {'+' if s - nadir >= 0 else ''}{f1(s - nadir)} mm" if nadir is not None else "") + ")")
        if news:
            txt += " ; nouvelle(s) lésion(s) retenue(s) : " + ", ".join(news)
        if pe_nc:
            txt += " ; progression sans équivoque des non-cibles signalée par le lecteur"
        details.append({"date": d, "provisoire": bool(e.get("provisoire")), "statut": st, "libelle": LIBELLES[st], "texte": txt,
                        "exclues": ex, "nouvelles": news, "somme_mm": s, "nadir_mm": nadir, "var_ref_pct": vr, "var_nadir_pct": vn,
                        "pe_non_cible": pe_nc})
    if not details:
        warns.append("Aucune date de suivi : appliquez le suivi à une autre séance")
    st = {x["date"]: x["statut"] for x in details}
    nad = {x["date"]: x["var_nadir_pct"] for x in details}
    cols = [d for d in order if date_key(d) >= date_key(r0)]
    lignes = [[f"{l['name']} ({_organ(l)}{', petit axe' if l.get('ganglion') else ''})"] +
              [f1(v) if v is not None else "— " + w for v, w in (sizes[d][k] for d in cols)] for k, l in enumerate(cib)]
    s0 = sums[r0]
    lignes += [["Somme (mm)"] + [f1(sums[d]) if sums[d] is not None else "—" for d in cols],
               ["% vs référence"] + [pc((sums[d] - s0) / s0 * 100) if sums[d] is not None and s0 and d != r0 else "—" for d in cols],
               ["% vs nadir"] + [pc(nad.get(d)) for d in cols],
               ["Statut"] + [st.get(d, "référence") for d in cols]]
    last = details[-1]["statut"] if details else "NE"
    return {"statut": last, "libelle": LIBELLES[last], "details": details, "cibles": [l["name"] for l in cib], "non_cibles": nonc,
            "tableau": {"colonnes": ["Cible"] + cols, "lignes": lignes}, "avertissements": list(dict.fromkeys(warns))}


METHOD = dict(id="recist11", nom="RECIST 1.1", donnees=("TEP", "CT"), cancers=("dsrct", "general"), compare=compare,
              note="lésions suivies, diamètres au CT de la TEP-CT (5 cibles, 2 par organe, petit axe des ganglions)")
