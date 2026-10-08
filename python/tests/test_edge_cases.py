"""
Tests for edge cases and previously fixed errors.
"""

import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pytest
from scipy.io import loadmat

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from featurecharacterization2d import (  # noqa: E402
    FeatureAttribute,
    Watershed,
    feature_characterization,
    inverse_material_ratio,
    plot_motifs,
)

PROFILES = Path(__file__).resolve().parents[2] / "data" / "profiles"
DX = 0.5e-3  # mm

# simple self-generated profile from the minimal examples
Z_MINIMAL = np.array([
    3.3, 2, 1, 5, 3.8, 4, 1.5, 1.5, 3.5, 2.5, 2, -1, 0, 3, 1.2, 2, -1.2, -5, -4,
    -4.5, -2, -2.3, 1, 3, 3, 3, 4, 4.5, 4.5, 4, 1.5, 1.5, 3.5, 4, 9, 8, -1, -1,
    -1, -1, 7, 7, 7, 0, 0.5, 3, 5, 4, 5, 4.5, 0.5, 1, 2, -1, 0, 3, 5.2, 5, 5.5,
    4, 7,
])
Z_MINIMAL = Z_MINIMAL - np.mean(Z_MINIMAL)


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def run_plot(*args):
    with warnings.catch_warnings():
        # plt.show() with the Agg backend
        warnings.filterwarnings("ignore", "FigureCanvasAgg is non-interactive")
        plot_motifs(*args)


# plot_motifs ---------------------------------------------------------------

@pytest.mark.parametrize(
    "z, TH",
    [(Z_MINIMAL, 2.0), (loadmat(PROFILES / "Bu_1_56_ak.mat")["z"], 3.0)],
    ids=["minimal", "Bu"],
)
def test_plot_motifs_multiple_height_intersections(z, TH):
    # merged motifs (pruning) can have more than one height intersection
    M = Watershed(z, DX, "D", "Wolfprune", TH).motifs()
    assert max(len(ihi) for ihi in M.ihi) > 1
    run_plot(z, DX, M)


@pytest.mark.parametrize("significant", ["Open 0.5", "Closed 0.5", "Closed 50 %"])
def test_plot_motifs_with_threshold(significant):
    _, M, meta = feature_characterization(
        Z_MINIMAL, DX, "D", "Wolfprune 5 %", significant, "HDh", "Mean"
    )
    run_plot(Z_MINIMAL, DX, M, meta["Fsig"], meta["NIsig"])


def test_plot_motifs_without_threshold():
    _, M, meta = feature_characterization(
        Z_MINIMAL, DX, "D", "Wolfprune 5 %", "All", "HDh", "Mean"
    )
    run_plot(Z_MINIMAL, DX, M, meta["Fsig"], meta["NIsig"])


# curvature -----------------------------------------------------------------

@pytest.mark.parametrize("i", [0, 1, 2, 3, 4, 35, 36, 37, 38, 39])
def test_curvature_stencils_exact_for_polynomial(i):
    # finite differences of order 6 are exact for polynomials up to degree 6
    x = np.arange(40.0)
    z = 1e-3 * (x - 20.0) ** 3 + 0.05 * x**2
    dz1 = 3e-3 * (x[i] - 20.0) ** 2 + 0.1 * x[i]
    dz2 = 6e-3 * (x[i] - 20.0) + 0.1
    # dx = 1e-3 mm = 1 µm
    cx = FeatureAttribute.curvature(z, 1e-3, np.array([float(i)]))
    np.testing.assert_allclose(cx, dz2 / (1.0 + dz1**2) ** 1.5, rtol=1e-9)


# empty motif sets ----------------------------------------------------------

def test_monotonic_profile_has_no_motifs():
    M = Watershed(np.arange(10.0), DX, "D", "Wolfprune", 1.0).motifs()
    assert len(M) == 0


def test_pruning_removes_all_motifs():
    M = Watershed(Z_MINIMAL, DX, "D", "Wolfprune", 1000.0).motifs()
    assert len(M) == 0
    with pytest.warns(UserWarning, match="No features detected"):
        xFC, _, meta = feature_characterization(
            Z_MINIMAL, DX, "D", "Wolfprune 1000", "All", "HDh", "Mean"
        )
    assert np.isnan(xFC)
    assert meta["nM"] == 0


# invalid options -----------------------------------------------------------

def test_unknown_statistic_raises():
    with pytest.raises(ValueError, match="Unknown attribute statistic"):
        feature_characterization(Z_MINIMAL, DX, "D", "Wolfprune 1", "All", "HDh", "Median")


def test_unknown_attribute_raises():
    with pytest.raises(ValueError, match="Unknown attribute type"):
        feature_characterization(Z_MINIMAL, DX, "D", "Wolfprune 1", "All", "Foo", "Mean")


# material ratio ------------------------------------------------------------

def test_inverse_material_ratio_bounds():
    assert inverse_material_ratio(Z_MINIMAL, 0.0) == 0.0
    np.testing.assert_allclose(
        inverse_material_ratio(Z_MINIMAL, 100.0), np.min(Z_MINIMAL) - np.max(Z_MINIMAL)
    )
    with pytest.raises(ValueError):
        inverse_material_ratio(Z_MINIMAL, 101.0)
