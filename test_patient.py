"""Données entièrement fictives : jamais de vraies données patient dans les tests."""
from pydicom.dataset import Dataset

from patient import patient_info


def ds(**kw):
    d = Dataset()
    for k, v in kw.items():
        setattr(d, k, v)
    return d


def test_full_header_is_formatted():
    rp = ds(Radiopharmaceutical="FDG", RadionuclideTotalDose="350000000", RadiopharmaceuticalStartTime="090000.00")
    h = ds(PatientName="DUPONT^JEAN-PIERRE", PatientID="X123", PatientBirthDate="19700131", PatientSex="M",
           PatientAge="054Y", PatientWeight="82.5", PatientSize="1.8", StudyDate="20261125", StudyTime="101500.000",
           AcquisitionTime="100000.00", RadiopharmaceuticalInformationSequence=[rp], Manufacturer="ACME",
           ManufacturerModelName="PET-9", ReferringPhysicianName="MARTIN^ANNE")
    info = patient_info([h])
    assert info["Nom"] == "Jean-Pierre DUPONT"
    assert info["Naissance"] == "31/01/1970" and info["Sexe"] == "Homme" and info["Âge"] == "54 ans"
    assert info["Poids"] == "82.5 kg" and info["Taille"] == "1.80 m"
    assert info["Date de l'examen"] == "25/11/2026" and info["Heure de l'examen"] == "10:15:00"
    assert info["Activité injectée"] == "350 MBq" and info["Heure d'injection"] == "09:00:00"
    assert info["Incorporation"] == "60 min"  # 10:00 - 09:00
    assert info["Médecin demandeur"] == "Anne MARTIN" and info["Machine"] == "ACME PET-9"


def test_absent_fields_are_omitted():
    assert patient_info([ds(PatientID="X1")]) == {"Identifiant": "X1"}
    assert patient_info([]) == {}


def test_uptake_across_midnight():
    rp = ds(RadiopharmaceuticalStartTime="233000")
    h = ds(RadiopharmaceuticalInformationSequence=[rp], AcquisitionTime="003000")
    assert patient_info([h])["Incorporation"] == "60 min"
