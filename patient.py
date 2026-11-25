"""Infos patient/examen lues dans les en-têtes DICOM de la série PET.
Affichage local uniquement : jamais écrites sur disque ni journalisées."""
from petscan import _scan, _seconds

AGE_UNITS = {"Y": "ans", "M": "mois", "W": "sem.", "D": "j"}
SEX = {"M": "Homme", "F": "Femme", "O": "Autre"}


def _date(s):  # YYYYMMDD -> JJ/MM/AAAA
    return f"{s[6:8]}/{s[4:6]}/{s[:4]}" if len(s) >= 8 and s[:8].isdigit() else s


def _time(s):  # HHMMSS[.ffffff] -> HH:MM:SS
    s = s.split(".")[0].ljust(6, "0")
    return f"{s[:2]}:{s[2:4]}:{s[4:6]}"


def _name(s):  # Famille^Prénom -> Prénom FAMILLE
    parts = s.split("^")
    given = parts[1].title() if len(parts) > 1 else ""
    return " ".join(p for p in (given, parts[0].upper()) if p)


def patient_info(headers):
    """headers : en-têtes pydicom de la série PET. Retourne {libellé: texte} ; les champs absents sont omis."""
    if not headers:
        return {}
    h = headers[0]
    g = lambda ds, k: str(ds.get(k, "") or "").strip()
    rp = (h.get("RadiopharmaceuticalInformationSequence") or [None])[0]
    info = {}

    def put(label, value):
        if value:
            info[label] = value

    put("Nom", _name(g(h, "PatientName")))
    put("Identifiant", g(h, "PatientID"))
    put("Naissance", _date(g(h, "PatientBirthDate")))
    put("Sexe", SEX.get(g(h, "PatientSex"), g(h, "PatientSex")))
    age = g(h, "PatientAge")
    put("Âge", f"{int(age[:-1])} {AGE_UNITS.get(age[-1], age[-1])}" if age[:-1].isdigit() else age)
    put("Poids", f"{float(g(h, 'PatientWeight')):g} kg" if g(h, "PatientWeight") else "")
    put("Taille", f"{float(g(h, 'PatientSize')):.2f} m" if g(h, "PatientSize") else "")
    put("Date de l'examen", _date(g(h, "StudyDate")))
    put("Heure de l'examen", _time(g(h, "StudyTime")) if g(h, "StudyTime") else "")
    put("Étude", g(h, "StudyDescription"))
    put("Série", g(h, "SeriesDescription"))
    if rp is not None:
        put("Radiopharmaceutique", g(rp, "Radiopharmaceutical"))
        put("Activité injectée", f"{float(g(rp, 'RadionuclideTotalDose')) / 1e6:.0f} MBq" if g(rp, "RadionuclideTotalDose") else "")
        inj = g(rp, "RadiopharmaceuticalStartTime")
        put("Heure d'injection", _time(inj) if inj else "")
        acq = [g(s, "AcquisitionTime") for s in headers if g(s, "AcquisitionTime")]
        if inj and acq:  # temps d'incorporation = début d'acquisition - injection
            put("Incorporation", f"{((min(_seconds(a) for a in acq) - _seconds(inj)) % 86400) / 60:.0f} min")
    put("Établissement", g(h, "InstitutionName"))
    put("Médecin demandeur", _name(g(h, "ReferringPhysicianName")))
    put("Machine", " ".join(x for x in (g(h, "Manufacturer"), g(h, "ManufacturerModelName")) if x))
    return info


def read_patient_info(folder):
    """Infos de la plus grande série PET du dossier (même règle que `load_suv`) ; {} s'il n'y en a pas."""
    series = _scan(folder, "PT")
    return patient_info([h for _, h in max(series.values(), key=len)]) if series else {}
