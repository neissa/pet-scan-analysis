"""Indice de pertinence d'une lésion (0–100) : trie les fixations de la plus probablement réelle à la moins probable
(fixation physiologique, bruit, excrétion, débord d'un organe voisin). HEURISTIQUE : poids fixés à la main, sans base
annotée ; ce n'est ni une probabilité calibrée ni un avis médical. Score = BASE + somme des points des critères, borné à
[0, 100] (et plafonné à EXCLU_MAX pour une lésion exclue par l'utilisateur). Un critère non évaluable vaut 0 point."""
from organs import ORGANS, PHYSIO

LIBELLE = "indice de pertinence (heuristique, pas une probabilité calibrée ni un avis médical)"
PHYSIO_NOMS = {ORGANS[k][0] for k in PHYSIO}

# --- poids (à recalibrer sur une base annotée) ---
BASE = 50
CONTRASTE_PAR_UNITE, CONTRASTE_ZERO, CONTRASTE_MIN, CONTRASTE_MAX = 10, 1.5, -15, 20  # (SULmax / foie − 1,5) × 10
VOL_BRUIT_ML, VOL_BRUIT = 0.3, -20      # quelques voxels : bruit
VOL_PETIT_ML, VOL_PETIT = 1.0, -8
VOL_GROS_ML, VOL_OK = 50.0, 5           # 1–50 ml : taille plausible ; au-delà, 0 (plateau, organe)
VOL_FENETRE = -10                       # « volume limité à la fenêtre » : région qui fuit (organe, fond)
PIC_NET, PIC_NET_PTS, PIC_MOU, PIC_MOU_PTS, PLATEAU_PTS = 1.6, 8, 1.3, 3, -5  # SUVmax / SUVmean
LOC_DANS, LOC_CONTACT, LOC_PRES = -45, -35, -15  # dans / au contact / à ≤ 15 mm d'un organe à fixation physiologique
EXCRETION = -30
HU_TISSU, HU_TISSU_PTS = (-200, 300), 5
HU_AIR, HU_AIR_PTS = -500, -15          # air : fixation hors tissu (recalage, artefact)
HU_OS, HU_OS_PTS = 700, -5              # os cortical dense
SUIVI_VALIDE, SUIVI_PAR_DATE, SUIVI_DATES_MAX, SUIVI_REF = 25, 5, 10, 5
EXCLU_MAX = 20
CLASSES = ((65, "haute"), (40, "moyenne"), (0, "faible"))

NE = "non évalué"


def _loc(lesion):
    """Pire localisation : (points, texte). organe « dans X » / « près de X » (`viewer._organ_at_peak`), physio (`analyze`),
    contact (`measure_at`). None si aucune information d'organe."""
    if not any(k in lesion for k in ("organe", "physio", "contact")):
        return None
    org, phys, cont = lesion.get("organe"), lesion.get("physio"), lesion.get("contact")
    out = [(0, org or "hors organe physiologique")]
    if org:
        nom = org.split(" ", 1)[1] if org.startswith("dans ") else org[len("près de "):] if org.startswith("près de ") else ""
        if nom in PHYSIO_NOMS:
            out.append((LOC_DANS if org.startswith("dans ") else LOC_PRES, org))
    if phys:
        out.append((LOC_CONTACT, f"au contact : {phys}"))
    c = sorted(PHYSIO_NOMS & set(cont or ()))
    if c:
        out.append((LOC_CONTACT, "au contact : " + ", ".join(c)))
    return min(out, key=lambda t: t[0])


def score(lesion, ctx=None):
    """lesion : champs de `petscan.analyze` (suv_max, suv_mean, mtv_ml, sul_max, physio) et de `viewer` (organe) ou de
    `petscan.measure_at` (vol_ml, flags, contact). ctx (tout facultatif) : foie = SUL moyen de la référence foie (même unité
    que sul_max) ; hu = densité CT au pic (clé absente ou None : pas de CT) ; excretion = bool ; suivi = {status
    ("valide"/"exclu"/None), dates (nb de dates mesurées), ref (présente à la référence)}.
    Retourne {score, brut, classe, criteres: [{nom, valeur, points}]} ; brut = somme des points (BASE comprise)."""
    ctx = ctx or {}
    cr = [{"nom": "base", "valeur": "", "points": BASE}]
    add = lambda nom, valeur, points: cr.append({"nom": nom, "valeur": valeur, "points": int(round(points))})
    flags = lesion.get("flags") or []
    loc = _loc(lesion)
    excr = bool(ctx.get("excretion") or "excrétion possible" in flags)
    debord = excr or (loc is not None and loc[0] < 0)  # près d'un organe physiologique, le SUVmax peut venir de l'organe

    foie, sm = ctx.get("foie"), lesion.get("sul_max", lesion.get("suv_max"))
    if foie and sm is not None:
        q = sm / foie
        p = min(max((q - CONTRASTE_ZERO) * CONTRASTE_PAR_UNITE, CONTRASTE_MIN), CONTRASTE_MAX)
        add("contraste (SULmax / foie)", f"{q:.1f}".replace(".", ",") + (" ; bonus non retenu : débord possible" if debord and p > 0 else ""),
            min(p, 0) if debord else p)
    else:
        add("contraste (SULmax / foie)", NE + " (pas de référence foie)", 0)

    v = lesion.get("mtv_ml", lesion.get("vol_ml"))
    if v is None:
        add("volume", NE, 0)
    else:
        add("volume", f"{v:.2f} ml".replace(".", ","), VOL_BRUIT if v < VOL_BRUIT_ML else VOL_PETIT if v < VOL_PETIT_ML
            else VOL_OK if v < VOL_GROS_ML else 0)
    if "volume limité à la fenêtre" in flags:
        add("volume limité à la fenêtre", "oui", VOL_FENETRE)

    mean = lesion.get("suv_mean")
    if mean:
        q = lesion["suv_max"] / mean
        add("pic isolé (SUVmax / SUVmean)", f"{q:.2f}".replace(".", ","),
            PIC_NET_PTS if q >= PIC_NET else PIC_MOU_PTS if q >= PIC_MOU else PLATEAU_PTS)
    else:
        add("pic isolé (SUVmax / SUVmean)", NE, 0)

    add("localisation", NE + " (organes non calculés)" if loc is None else loc[1], 0 if loc is None else loc[0])

    if excr:
        add("excrétion possible", "oui", EXCRETION)

    hu = ctx.get("hu")
    if hu is None:
        add("densité CT au pic", NE + " (pas de CT)", 0)
    else:
        add("densité CT au pic", f"{hu:.0f} HU", HU_AIR_PTS if hu < HU_AIR else HU_OS_PTS if hu > HU_OS
            else HU_TISSU_PTS if HU_TISSU[0] <= hu <= HU_TISSU[1] else 0)

    s = ctx.get("suivi")
    if s:
        if s.get("status") == "valide":
            add("validée par l'utilisateur", "oui", SUIVI_VALIDE)
        elif s.get("status") == "exclu":
            add("exclue par l'utilisateur", f"score plafonné à {EXCLU_MAX}", 0)
        n = s.get("dates") or 0
        add("persistance", f"{n} date(s)", min(max(n - 1, 0) * SUIVI_PAR_DATE, SUIVI_DATES_MAX))
        if s.get("ref"):
            add("présente à la référence", "oui", SUIVI_REF)

    brut = sum(c["points"] for c in cr)
    sc = min(max(brut, 0), 100)
    if s and s.get("status") == "exclu":
        sc = min(sc, EXCLU_MAX)
    return {"score": sc, "brut": brut, "classe": next(n for t, n in CLASSES if sc >= t), "criteres": cr}


def annotate(lesion, ctx=None):
    """Ajoute à lesion les champs pertinence, pertinence_classe, pertinence_detail (critères) ; renvoie lesion."""
    r = score(lesion, ctx)
    lesion.update(pertinence=r["score"], pertinence_classe=r["classe"], pertinence_detail=r["criteres"])
    return lesion
