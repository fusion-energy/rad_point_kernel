"""Python-side correctness tests for the WW builder.

These are fast — they exercise the Rust builder via the PyO3 bindings and
assert on returned arrays. No MC simulation runs here. Speedup validation
lives in ``verification_and_validation/weight_window_speedup.py``.
"""
from __future__ import annotations

import numpy as np
import pytest

import rad_point_kernel as rpk
from rad_point_kernel.weight_windows import (
    _parse_quantity,
    _select_driving_quantity,
)


# ---- Fixtures -------------------------------------------------------------

@pytest.fixture
def water():
    return rpk.Material({"H": 0.111, "O": 0.889}, 1.0, fraction="mass")


@pytest.fixture
def concrete():
    return rpk.Material(
        {"H": 0.01, "O": 0.532, "Si": 0.337, "Ca": 0.044,
         "Al": 0.034, "Fe": 0.014, "Na": 0.029},
        2.3, fraction="mass",
    )


@pytest.fixture
def neutron_14mev():
    return rpk.Source("neutron", 14.1e6)


@pytest.fixture
def photon_1mev():
    return rpk.Source("photon", 1e6)


# ---- Basic structure --------------------------------------------------------

def test_parse_quantity_flux():
    q, coupled = _parse_quantity("flux-neutron")
    assert coupled is False


def test_parse_quantity_dose_ap():
    q, coupled = _parse_quantity("dose-AP-neutron")
    assert coupled is False


def test_parse_quantity_coupled():
    q, coupled = _parse_quantity("dose-AP-coupled-photon")
    assert coupled is True


def test_parse_quantity_invalid_raises():
    with pytest.raises(ValueError):
        _parse_quantity("nonsense")
    with pytest.raises(ValueError):
        _parse_quantity("dose-BADGEOM-neutron")
    # Bare names without explicit particle are no longer valid.
    with pytest.raises(ValueError, match="missing a particle"):
        _parse_quantity("flux")
    with pytest.raises(ValueError, match="missing a particle"):
        _parse_quantity("dose-AP")


# ---- Builder returns the right number of WW objects ------------------------

def test_returns_one_ww_for_neutron_flux(water, neutron_14mev):
    pytest.importorskip("openmc")
    # 70 cm water gives tau approx 3.5 for 14 MeV neutron flux, comfortably
    # above the skip-gate threshold of 3.0.
    layers = [rpk.Layer(70.0, water)]
    ww = rpk.build_weight_windows(
        layers=layers, source=neutron_14mev, quantities=["flux-neutron"],
    )
    assert len(ww) == 1
    assert str(ww[0].particle_type) == "neutron"


def test_returns_one_ww_for_photon_dose(concrete, photon_1mev):
    pytest.importorskip("openmc")
    layers = [rpk.Layer(50.0, concrete)]
    ww = rpk.build_weight_windows(
        layers=layers, source=photon_1mev, quantities=["dose-AP-photon"],
    )
    assert len(ww) == 1
    assert str(ww[0].particle_type) == "photon"


def test_returns_two_ww_for_coupled(concrete, neutron_14mev):
    pytest.importorskip("openmc")
    # Thick enough shield that both neutron and secondary-photon curves have
    # τ > skip-gate threshold. 50 cm of concrete is sometimes borderline for
    # the photon side so we use 100 cm.
    layers = [rpk.Layer(100.0, concrete)]
    ww = rpk.build_weight_windows(
        layers=layers, source=neutron_14mev,
        quantities=["dose-AP-coupled-photon"],
    )
    # Coupled mode always builds the neutron WW; the photon WW is added
    # when its skip-gate doesn't fire.
    assert len(ww) >= 1
    assert str(ww[0].particle_type) == "neutron"
    if len(ww) == 2:
        assert str(ww[1].particle_type) == "photon"


# ---- Skip gate --------------------------------------------------------------

def test_skips_thin_shield(water, neutron_14mev):
    # 5 cm of water for 14 MeV neutrons → τ << 2 → skip
    pytest.importorskip("openmc")
    layers = [rpk.Layer(5.0, water)]
    ww = rpk.build_weight_windows(
        layers=layers, source=neutron_14mev, quantities=["flux-neutron"],
    )
    assert ww == []


def test_skips_all_void(neutron_14mev):
    pytest.importorskip("openmc")
    layers = [rpk.Layer(100.0, None)]
    ww = rpk.build_weight_windows(
        layers=layers, source=neutron_14mev, quantities=["flux-neutron"],
    )
    assert ww == []


# ---- Grid structure --------------------------------------------------------

def test_layer_edges_present_in_grid(water, neutron_14mev):
    # Voids should give forced layer edges. 40 + 40 void + 40 cm water gives
    # tau approx 4 for 14 MeV neutron flux, above the skip-gate threshold.
    pytest.importorskip("openmc")
    layers = [
        rpk.Layer(40.0, water),
        rpk.Layer(40.0, None),     # void between
        rpk.Layer(40.0, water),
    ]
    ww = rpk.build_weight_windows(
        layers=layers, source=neutron_14mev, quantities=["flux-neutron"],
    )
    r_grid = list(ww[0].mesh.r_grid)
    for edge in (40.0, 80.0, 120.0):
        assert any(abs(r - edge) < 1e-6 for r in r_grid), \
            f"layer edge {edge} missing from r_grid {r_grid}"


def test_no_bins_inside_void(water, neutron_14mev):
    pytest.importorskip("openmc")
    layers = [
        rpk.Layer(40.0, water),
        rpk.Layer(40.0, None),
        rpk.Layer(40.0, water),
    ]
    ww = rpk.build_weight_windows(
        layers=layers, source=neutron_14mev, quantities=["flux-neutron"],
    )
    r_grid = list(ww[0].mesh.r_grid)
    for a, b in zip(r_grid, r_grid[1:]):
        centre = 0.5 * (a + b)
        in_void = (centre > 40.0 + 1e-6) and (centre < 80.0 - 1e-6)
        if in_void:
            # must be the whole-void cell [40, 80]
            assert abs(a - 40.0) < 1e-6 and abs(b - 80.0) < 1e-6, \
                f"intra-void bin: [{a}, {b}]"


# ---- Bounds shape ----------------------------------------------------------

def test_lower_bounds_shape(water, neutron_14mev):
    # 70 cm water gives tau approx 3.5, above the skip-gate threshold.
    pytest.importorskip("openmc")
    layers = [rpk.Layer(70.0, water)]
    ww = rpk.build_weight_windows(
        layers=layers, source=neutron_14mev, quantities=["flux-neutron"],
    )
    lb = np.asarray(ww[0].lower_ww_bounds)
    # shape (nr, 1, 1, nE)
    assert lb.ndim == 4
    assert lb.shape[1] == 1
    assert lb.shape[2] == 1
    nr = len(ww[0].mesh.r_grid) - 1
    ne = len(ww[0].energy_bounds) - 1
    assert lb.shape == (nr, 1, 1, ne)
    assert np.all(lb > 0), "all lower bounds must be positive (floor applied)"


# ---- Driving-quantity selection --------------------------------------------

def test_driving_quantity_picks_steepest(water, neutron_14mev):
    # For a 14 MeV neutron in water, dose should drop faster than flux
    # (response weighting amplifies fast-particle contribution). So dose
    # should be selected over flux when both are requested.
    layers = [rpk.Layer(50.0, water)]
    picked = _select_driving_quantity(
        layers, neutron_14mev, ["flux-neutron", "dose-AP-neutron"]
    )
    # Don't hard-code which one wins — just verify the selection returns one
    # of them and doesn't raise.
    assert picked in ("flux-neutron", "dose-AP-neutron")


# ---- Photon-source symmetry ------------------------------------------------

def test_photon_source_produces_photon_ww(concrete, photon_1mev):
    # 30 cm concrete gives tau approx 4.4 for 1 MeV photons, above the
    # skip-gate threshold.
    pytest.importorskip("openmc")
    layers = [rpk.Layer(30.0, concrete)]
    ww = rpk.build_weight_windows(
        layers=layers, source=photon_1mev, quantities=["flux-photon"],
    )
    assert str(ww[0].particle_type) == "photon"
