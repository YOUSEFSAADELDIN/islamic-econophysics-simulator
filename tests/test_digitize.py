"""The reference values used for comparison must agree with the PDF's own figures."""

import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("PIL")
if shutil.which("pdfimages") is None:
    pytest.skip("pdfimages (poppler) not installed", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
import digitize_kato_figures as dg  # noqa: E402
from islamic_econophysics import kato2022_reference as ref  # noqa: E402


@pytest.fixture(scope="module")
def digitised():
    return dg.digitise(ROOT / "journal.pone.0275113.pdf")


def test_y_axis_calibration_error_is_about_a_pixel(digitised):
    # residual of the label-centre fit = label misalignment, up to ~1 px = ~0.002 Gini
    for fig in ("fig4", "fig5", "fig6", "fig7"):
        assert digitised[fig]["info"]["y_fit_resid_gini"] < 0.003
        assert 0.0015 < digitised[fig]["info"]["gini_per_px"] < 0.0018   # ~0.002 Gini per pixel


def test_fig4_markers_match_reference(digitised):
    s = digitised["fig4"]["series"]
    navy = [p["gini"] for p in s["navy"] if not p["hollow"]]
    blue = [p["gini"] for p in s["blue"] if not p["hollow"]]
    assert np.allclose(navy, ref.FIG4["J dw=0.2"], atol=0.002)
    assert np.allclose(blue, ref.FIG4["J dw=0.1"], atol=0.002)


def test_fig6_and_fig7_markers_match_reference(digitised):
    f6 = [p["gini"] for p in digitised["fig6"]["series"]["purple"] if not p["hollow"]]
    assert np.allclose(f6, ref.FIG6[("L", 10**5)][1:], atol=0.002)
    f7 = [p["gini"] for p in digitised["fig7"]["series"]["blue"] if not p["hollow"]]
    assert np.allclose(f7, ref.FIG7["J0.1"], atol=0.002)
