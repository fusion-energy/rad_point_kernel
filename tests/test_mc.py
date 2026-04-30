"""Monte Carlo integration tests.

These tests require OpenMC and the OPENMC_CROSS_SECTIONS environment variable
pointing at an ENDF/B-VIII.1 (or compatible) cross_sections.xml. They are
skipped automatically when either is missing.
"""

import os

import pytest

import rad_point_kernel as rpk


openmc = pytest.importorskip("openmc")

if not os.environ.get("OPENMC_CROSS_SECTIONS"):
    pytest.skip(
        "OPENMC_CROSS_SECTIONS not set; skipping MC tests",
        allow_module_level=True,
    )


@pytest.fixture(scope="module")
def iron():
    return rpk.Material(composition={"Fe": 1.0}, density=7.874)


def test_compute_buildup_single_geometry(iron):
    layers = [rpk.Layer(thickness=5, material=iron)]
    source = rpk.Source(particle="photon", energy=1e6)

    results = rpk.compute_buildup(
        geometries=[layers],
        source=source,
        quantities=["dose-AP-photon"],
        particles_per_batch=2_000,
        max_batches=20,
        trigger_rel_err=0.1,
    )

    assert len(results) == 1
    r = results[0]
    assert r.buildup["dose-AP-photon"] > 1.0
    assert r.mc["dose-AP-photon"] > 0.0
    assert r.pk["dose-AP-photon"] > 0.0
    assert r.mc_std_dev["dose-AP-photon"] >= 0.0


def test_compute_buildup_multiple_geometries(iron):
    thicknesses = [5, 10]
    geometries = [[rpk.Layer(thickness=t, material=iron)] for t in thicknesses]
    source = rpk.Source(particle="photon", energy=1e6)

    results = rpk.compute_buildup(
        geometries=geometries,
        source=source,
        quantities=["dose-AP-photon"],
        particles_per_batch=2_000,
        max_batches=20,
        trigger_rel_err=0.1,
    )

    assert len(results) == len(thicknesses)
    buildups = [r.buildup["dose-AP-photon"] for r in results]
    # Thicker shield should produce a larger build-up factor for photons.
    assert buildups[1] > buildups[0]


def test_buildup_result_applied_to_calculate_dose(iron):
    layers = [rpk.Layer(thickness=10, material=iron)]
    source = rpk.Source(particle="photon", energy=1e6)

    results = rpk.compute_buildup(
        geometries=[layers],
        source=source,
        quantities=["dose-AP-photon"],
        particles_per_batch=2_000,
        max_batches=20,
        trigger_rel_err=0.1,
    )
    r = results[0]

    corrected = rpk.calculate_dose(
        layers=layers,
        source=source,
        geometry="AP",
        buildup=r,
    ).scale(strength=1e12)
    uncollided = rpk.calculate_dose(
        layers=layers,
        source=source,
        geometry="AP",
    ).scale(strength=1e12)

    assert corrected.dose > uncollided.dose
    assert corrected.buildup_factor == pytest.approx(r.buildup["dose-AP-photon"])


def test_compute_buildup_synthesizes_dose_total(iron):
    """Requesting dose-AP-neutron + dose-AP-coupled-photon for a neutron source
    must auto-add a dose-AP-total quantity that equals the sum of the two doses."""
    layers = [rpk.Layer(thickness=5, material=iron)]
    source = rpk.Source(particle="neutron", energy=14.06e6)

    results = rpk.compute_buildup(
        geometries=[layers],
        source=source,
        quantities=["dose-AP-neutron", "dose-AP-coupled-photon"],
        particles_per_batch=2_000,
        max_batches=20,
        trigger_rel_err=0.1,
    )

    r = results[0]
    assert "dose-AP-total" in r.mc
    assert "dose-AP-total" in r.mc_std_dev
    assert r.mc["dose-AP-total"] == pytest.approx(
        r.mc["dose-AP-neutron"] + r.mc["dose-AP-coupled-photon"], rel=1e-12
    )
    # Variance combines in quadrature.
    expected_std = (
        r.mc_std_dev["dose-AP-neutron"] ** 2
        + r.mc_std_dev["dose-AP-coupled-photon"] ** 2
    ) ** 0.5
    assert r.mc_std_dev["dose-AP-total"] == pytest.approx(expected_std, rel=1e-12)
    # Buildup uses neutron PK as reference.
    assert r.pk["dose-AP-total"] == pytest.approx(r.pk["dose-AP-neutron"])
    assert r.buildup["dose-AP-total"] == pytest.approx(
        r.mc["dose-AP-total"] / r.pk["dose-AP-neutron"], rel=1e-12
    )


def test_compute_buildup_accepts_total_dose_shorthand(iron):
    """quantities=['dose-AP-total'] is shorthand for both halves; the result
    must carry all three keys (neutron, coupled-photon, total)."""
    layers = [rpk.Layer(thickness=5, material=iron)]
    source = rpk.Source(particle="neutron", energy=14.06e6)

    results = rpk.compute_buildup(
        geometries=[layers],
        source=source,
        quantities=["dose-AP-total"],
        particles_per_batch=2_000,
        max_batches=20,
        trigger_rel_err=0.1,
    )

    r = results[0]
    assert "dose-AP-neutron" in r.mc
    assert "dose-AP-coupled-photon" in r.mc
    assert "dose-AP-total" in r.mc
    assert r.mc["dose-AP-total"] == pytest.approx(
        r.mc["dose-AP-neutron"] + r.mc["dose-AP-coupled-photon"], rel=1e-12
    )
