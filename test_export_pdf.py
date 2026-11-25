import pytest

from export_pdf import build_pdf
from viewer import demo, render

SUV, SPACING, CT = demo()


def _render(view, i, vmax, cmap, mode, alpha, box):
    return render(SUV, view, i, vmax, cmap, SPACING, CT, mode, alpha, box)


def _state(**kw):
    s = {"views": [{"view": v, "i": 30, "z": 1, "x0": 0, "y0": 0} for v in ("axial", "coronal", "sagittal")],
         "vmax": 10, "cmap": "jet", "mode": "fusion", "alpha": 0.6, "measures": []}
    s.update(kw)
    return s


def test_pdf_with_patient_and_measures():
    meas = [{"v": "axial", "i": 30, "text": "Règle axial 31 : 12.0 mm",
             "items": [{"t": "line", "a": [.2, .3], "b": [.5, .6], "color": "#0f0", "bars": True, "text": "12.0 mm"},
                       {"t": "dots", "pts": [[.3, .3], [.31, .3]], "color": "rgba(255,255,0,.8)"}]}]
    pdf = build_pdf(_render, SUV.shape, _state(measures=meas), {"Nom": "Test PATIENT", "Poids": "70 kg"})
    assert pdf.startswith(b"%PDF") and b"/Count 1" in pdf and len(pdf) > 10_000


def test_pdf_zoomed_without_patient_or_measures():
    views = [{"view": v, "i": 10, "z": 4, "x0": 0.3, "y0": 0.2} for v in ("axial", "coronal", "sagittal")]
    assert build_pdf(_render, SUV.shape, _state(views=views), None).startswith(b"%PDF")


@pytest.mark.parametrize("bad", [
    {"views": [{"view": "axial", "i": 0, "z": 1, "x0": 0, "y0": 0}]},  # il manque 2 vues
    {"views": [{"view": v, "i": 0, "z": 99, "x0": 0, "y0": 0} for v in ("axial", "coronal", "sagittal")]},  # zoom hors limites
])
def test_invalid_state_is_rejected(bad):
    with pytest.raises(ValueError):
        build_pdf(_render, SUV.shape, _state(**bad), None)
