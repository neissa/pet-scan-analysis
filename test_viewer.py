import pytest

from viewer import demo, measure, render

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


def test_render_index_clipped():
    suv, *_ = demo()
    assert render(suv, "axial", 10_000, 10).startswith(PNG)
