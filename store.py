"""Seul stockage sur disque de l'appli : un dossier dédié (réglages, dossiers récents, cartes de lésions de référence).
Patients pseudonymisés : un patient = empreinte SHA-256 (sel aléatoire local) de l'identifiant + date de naissance ;
jamais de nom, d'identifiant ni de date de naissance en clair. Tout se supprime depuis le panneau « Données »."""
import hashlib
import json
import os
import re
import secrets
from pathlib import Path

# PETSCAN_CONFIG (appli Android) : ancien fichier de réglages ; le dossier de données est à côté.
_CFG = os.environ.get("PETSCAN_CONFIG")
DIR = Path(_CFG).parent / "netscan" if _CFG else Path.home() / ".netscan"
LEGACY = [Path(_CFG)] if _CFG else [Path.home() / ".pet-scan-analysis.json"]  # migrés puis supprimés
_KEY = re.compile(r"^[0-9a-f]{64}$")


def _read(path, default):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def _write(path, obj, raw=False):
    """Écriture atomique : fichier temporaire puis remplacement (jamais de fichier à moitié écrit)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(obj if raw else json.dumps(obj, ensure_ascii=False), encoding="utf-8")
    os.replace(tmp, path)


def config():
    """Réglages {recent: [...], recent_on: bool} ; migre l'ancien fichier de dossiers récents (puis le supprime)."""
    cfg = {"recent": [], "recent_on": False, **_read(DIR / "config.json", {})}
    for old in LEGACY:
        if old.is_file():
            rec = _read(old, {}).get("recent", [])
            cfg["recent"] = list(dict.fromkeys(cfg["recent"] + [r for r in rec if isinstance(r, str)]))[:8]
            _write(DIR / "config.json", cfg)
            old.unlink()
    return cfg


def set_config(**kw):
    _write(DIR / "config.json", {**config(), **kw})


def recent():
    """Dossiers récents ; vide tant que l'option n'est pas activée (désactivée par défaut)."""
    c = config()
    return c["recent"] if c["recent_on"] else []


def remember(folder):
    c = config()
    if c["recent_on"]:
        set_config(recent=([folder] + [f for f in c["recent"] if f != folder])[:8])


def _salt(create):
    p = DIR / "salt"
    try:
        return p.read_text(encoding="ascii").strip()
    except OSError:
        if not create:
            return None
        s = secrets.token_hex(32)
        _write(p, s, raw=True)
        return s


def patient_key(patient, create=False):
    """Pseudonyme stable du patient (identifiant + naissance des en-têtes) ; None sans identifiant, ou sans sel local
    si create=False (rien n'est écrit tant qu'aucune référence n'est définie)."""
    pid = str(patient.get("Identifiant") or "").strip()
    s = _salt(create) if pid else None
    return hashlib.sha256(f"{s}\0{pid}\0{patient.get('Naissance', '')}".encode()).hexdigest() if s else None


def _pfile(key):
    if not isinstance(key, str) or not _KEY.match(key):  # la clé vient de la page : jamais de chemin arbitraire
        raise ValueError("patient inconnu")
    return DIR / "patients" / f"{key}.json"


def load_map(key):
    return _read(_pfile(key), None)


def save_map(key, m):
    _write(_pfile(key), m)


def date_key(d):
    """« JJ/MM/AAAA » -> tri chronologique."""
    p = str(d).split("/")
    return tuple(reversed(p)) if len(p) == 3 else (str(d),)


def maps():
    """[(clé, carte)] de tous les patients mémorisés."""
    return [(p.stem, m) for p in sorted((DIR / "patients").glob("*.json")) if _KEY.match(p.stem) and (m := _read(p, None))]


def _files(what="all", key=None):
    if what == "patient":
        p = _pfile(key)
        return [p] if p.is_file() else []
    if what == "recent":
        return [DIR / "config.json"] if (DIR / "config.json").is_file() else []
    return sorted(p for p in DIR.rglob("*") if p.is_file()) if DIR.is_dir() else []


def preview(what, key=None):
    """Fichiers concernés par une suppression : [{nom, octets}] (chemins relatifs au dossier de données, sans identité)."""
    return [{"nom": p.relative_to(DIR).as_posix(), "octets": p.stat().st_size} for p in _files(what, key)]


def inventory():
    ms = maps()
    c = config()
    files = preview("all")
    return {"patients": len(ms), "lesions": sum(len(m.get("lesions", [])) for _, m in ms),
            "dates": sum(len(m.get("dates", {})) for _, m in ms), "fichiers": files,
            "octets": sum(f["octets"] for f in files), "recents": len(c["recent"]), "recent_on": c["recent_on"]}


def delete(what, key=None):
    """Supprime un patient, les dossiers récents ou tout ; renvoie le nombre de fichiers supprimés ou vidés."""
    if what == "recent":
        n = len(_files("recent"))
        if n:
            set_config(recent=[])
        return n
    if what not in ("patient", "all"):
        raise ValueError("suppression inconnue")
    files = _files(what, key)
    for p in files:
        p.unlink()
    if what == "all":
        for d in sorted((p for p in DIR.rglob("*") if p.is_dir()), reverse=True) + [DIR]:
            try:
                d.rmdir()
            except OSError:
                pass
    return len(files)
