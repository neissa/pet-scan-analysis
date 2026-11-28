"""Registre des méthodes d'analyse d'évolution, enfichables par type de cancer.
Une méthode = un module `method_*.py` qui définit METHOD (dict des champs de `Method`) ; l'ajouter = créer le module et
l'importer ci-dessous (import explicite : PyInstaller et la copie Android ne voient que les imports statiques).
viewer.py et index.html ne connaissent aucune méthode : ils lisent `catalog()`, valident par `choose()` et appellent `avis()`."""
from dataclasses import dataclass

CANCERS = {"dsrct": "Sarcome desmoplastique à petites cellules rondes", "general": "Autre / général"}
MENTION = "Aide à la lecture calculée sur les lésions suivies : ce n'est pas un avis médical. Recherche uniquement."


@dataclass(frozen=True)
class Method:
    id: str
    nom: str
    donnees: tuple             # données requises : ("TEP",), ("TEP", "CT")
    cancers: tuple             # identifiants de CANCERS concernés
    compare: object = None     # compare(reference, mesures_par_date, contexte) -> {statut, details, avertissements} ; None = « bientôt »
    defaut: tuple = ()         # types de cancer dont c'est la méthode par défaut
    note: str = ""


REGISTRY = {}


def register(m):
    REGISTRY[m.id] = m
    return m


import method_percist  # noqa: E402  (METHOD = dict des champs de Method : pas d'import circulaire)
from method_percist import _why  # noqa: E402
import method_recist  # noqa: E402

for _mod in (method_percist, method_recist):
    register(Method(**_mod.METHOD))


def default(cancer):
    """Méthode par défaut du type de cancer : celle qui s'y déclare par défaut, sinon la première disponible qui s'y applique."""
    ok = [m for m in REGISTRY.values() if cancer in m.cancers and m.compare]
    return next((m.id for m in ok if cancer in m.defaut), ok[0].id if ok else None)


def catalog():
    return {"cancers": [{"id": k, "nom": v, "defaut": default(k)} for k, v in CANCERS.items()],
            "methodes": [{"id": m.id, "nom": m.nom, "donnees": list(m.donnees), "cancers": list(m.cancers),
                          "disponible": m.compare is not None, "note": m.note} for m in REGISTRY.values()],
            "mention": MENTION}


def choose(cancer, methode):
    """Choix validé {cancer, methode} (ValueError si inconnu, non disponible ou hors du type de cancer)."""
    m = REGISTRY.get(methode)
    if cancer not in CANCERS or m is None or m.compare is None or cancer not in m.cancers:
        raise ValueError("méthode ou type de cancer invalide")
    return {"cancer": cancer, "methode": methode}


def heterogeneite(ref, dates):
    """Réponse hétérogène, commune aux méthodes : à chaque date, chaque lésion de référence retenue (hors excrétion/contact,
    comme PERCIST) classée en baisse / stable / hausse d'après son SULpeak (SULmax sous 1 ml) vs référence, mêmes seuils que
    PERCIST (≥ 30 % et ≥ 0,8). « Mixte » : au moins une en baisse et une en hausse. {date: {...}}"""
    r0 = ref["date"]
    base = (dates.get(r0) or {}).get("values", {})
    out = {}
    for d, e in dates.items():
        if d == r0:
            continue
        cat = {"baisse": [], "stable": [], "hausse": [], "non_evaluable": []}
        for l in ref["lesions"]:
            if l.get("status") != "valide" or l.get("since", r0) != r0:
                continue
            vb, vc = base.get(l["id"]), e["values"].get(l["id"])
            if _why(vb) or _why(vc):
                continue
            if not vb or not vc:
                cat["non_evaluable"].append(l["name"])
                continue
            k = "sul_peak" if min(vb["vol_ml"], vc["vol_ml"]) >= 1.0 else "sul_max"
            dabs = vc[k] - vb[k]
            pct = dabs / vb[k] * 100 if vb[k] > 0 else 0.0
            cat["hausse" if pct >= 30 and dabs >= 0.8 else "baisse" if pct <= -30 and dabs <= -0.8 else "stable"].append(l["name"])
        out[d] = {**cat, "comptes": {k: len(v) for k, v in cat.items()}, "mixte": bool(cat["baisse"] and cat["hausse"])}
    return out


def avis(m, dates, contexte=None):
    """Avis de statut de la carte `m` (store) ; dates : {date: mesures} dans l'ordre chronologique, référence comprise.
    Sans choix mémorisé : type « Autre / général » et sa méthode par défaut."""
    ch = m.get("analyse") or {"cancer": "general", "methode": default("general")}
    meth = REGISTRY.get(ch["methode"])
    if meth is None or meth.compare is None:
        return {"error": f"Méthode indisponible : {ch['methode']}"}
    ref = {"date": m["ref_date"], "lesions": m["lesions"], "unit": m.get("unit"), "spacing": m.get("spacing")}
    r = meth.compare(ref, dates, contexte or {})
    het = heterogeneite(ref, dates)
    for x in r.get("details", []):
        if x["date"] in het:
            x["heterogeneite"] = het[x["date"]]
    return {**r, "choix": ch, "methode": meth.nom, "mention": MENTION}
