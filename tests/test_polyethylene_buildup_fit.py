"""Regression test: BuildupFit must fit a polyethylene MC dataset cleanly.

Polyethylene at 14.06 MeV is a stress case for B(thickness) in 1D:

  - At small thickness B is dominated by the uncollided 14 MeV neutron
    transmission and approaches 1.
  - From ~0.5 m to ~1.25 m B falls fast as fast neutrons are
    moderated/captured (mostly by H).
  - Past ~1.5 m B turns over and grows because the dose at the outer
    surface becomes dominated by the H(n,γ)D 2.2 MeV secondary-photon
    cascade, while the uncollided neutron dose used as the point-kernel
    reference keeps decaying. The result is a U-shape in B_total over
    [25 cm, 250 cm].

This test asserts:

  1. Each per-quantity fit (`dose-AP-neutron`, `dose-AP-coupled-photon`,
     `dose-AP-total`) reproduces the MC anchor values it was fit to,
     within a few times the per-anchor MC relative error.
  2. The composition property `B_total = B_neutron + B_secondary_photon`,
     which holds exactly for the MC inputs at every anchor, also holds
     in the fit predictions.
  3. The total-dose extrapolation past the deepest anchor stays
     physically bounded.

Polyethylene composition: PNNL Materials Compendium 2021 v2 polyethylene
(C 0.857143 H 0.142857 by mass, density 0.93 g/cc), 25 m radial budget,
spherical geometry, AP irradiation, weight windows on, ENDF/B-VIII.1.
"""

from __future__ import annotations

import math

import pytest
import rad_point_kernel as rpk


# ---------------------------------------------------------------------------
# MC anchors — per source neutron, 14.06 MeV. From OpenMC w/ weight windows,
# converged to <5% relative error on dose-AP-total at every thickness.
# ---------------------------------------------------------------------------
_NEUTRON = "dose-AP-neutron"
_PHOTON = "dose-AP-coupled-photon"
_TOTAL = "dose-AP-total"

MC_ANCHORS = [
    # thickness_cm, mc_n,         mc_p,         mc_t,         pk_n,         std_n,        std_p,        std_t
    (25,  1.634591e-18, 3.127077e-20, 1.665861e-18, 1.690895e-18, 3.688986e-20, 9.953057e-22, 3.690328e-20),
    (50,  2.613481e-19, 1.829294e-20, 2.796411e-19, 4.537608e-19, 1.297116e-20, 3.370850e-22, 1.297554e-20),
    (75,  3.206895e-20, 7.525256e-21, 3.959421e-20, 1.217692e-19, 1.533068e-21, 1.659673e-22, 1.542025e-21),
    (125, 4.069256e-22, 1.196260e-21, 1.603186e-21, 8.769159e-21, 2.026697e-23, 8.944141e-24, 2.215282e-23),
    (175, 3.244167e-24, 2.088413e-22, 2.120854e-22, 6.315076e-22, 6.971206e-25, 2.455867e-24, 2.552893e-24),
    (250, 1.710954e-27, 1.734305e-23, 1.734476e-23, 1.220421e-23, 8.208152e-28, 7.441074e-25, 7.441078e-25),
]


def _build_result(t_cm, mc_n, mc_p, mc_t, pk_n, std_n, std_p, std_t):
    """Construct a BuildupResult dict equivalent to one cache entry."""
    pk_t = pk_n  # point-kernel reference is the uncollided-neutron dose
    return {
        "mc":         {_NEUTRON: mc_n, _PHOTON: mc_p, _TOTAL: mc_t},
        "pk":         {_NEUTRON: pk_n, _TOTAL: pk_t},
        "mc_std_dev": {_NEUTRON: std_n, _PHOTON: std_p, _TOTAL: std_t},
        "buildup":    {_NEUTRON: mc_n / pk_n, _TOTAL: mc_t / pk_t},
        "optical_thickness": 0.0,  # not used by the fit
        "source_strength": 1.0,
    }


@pytest.fixture
def fit():
    points = [{"thickness": t} for t, *_ in MC_ANCHORS]
    results = [rpk.BuildupResult.from_dict(_build_result(*row)) for row in MC_ANCHORS]
    for r in results:
        if hasattr(r, "synthesize_dose_totals"):
            r.synthesize_dose_totals()
    return rpk.BuildupFit(points=points, results=results)


def _b_mc(row, q):
    t, mc_n, mc_p, mc_t, pk_n, *_ = row
    bs = {_NEUTRON: mc_n / pk_n, _PHOTON: mc_p / pk_n, _TOTAL: mc_t / pk_n}
    return bs[q]


def _rel_err(row, q):
    t, mc_n, mc_p, mc_t, pk_n, std_n, std_p, std_t = row
    rel = {_NEUTRON: std_n / mc_n if mc_n > 0 else math.inf,
           _PHOTON: std_p / mc_p if mc_p > 0 else math.inf,
           _TOTAL: std_t / mc_t if mc_t > 0 else math.inf}
    return rel[q]


# ---------------------------------------------------------------------------
# Property 1: each per-quantity fit reproduces its anchors within a few
# times the per-anchor MC relative error.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("quantity", [_NEUTRON, _PHOTON, _TOTAL])
def test_fit_hits_each_anchor(fit, quantity):
    """Per-anchor: |fit - mc| / mc <= max(3 * mc_rel_err, 0.05).

    A 5 % floor handles anchors with sub-1 % MC noise (we don't expect
    the fit to be that perfect on a 6-point dataset).
    """
    misses = []
    for row in MC_ANCHORS:
        t = row[0]
        b_mc = _b_mc(row, quantity)
        if b_mc <= 0:
            continue
        rel_err = _rel_err(row, quantity)
        b_fit = fit.interpolate(quantity=quantity, thickness=t, warn=False).value
        rel_resid = abs(b_fit - b_mc) / b_mc
        tolerance = max(3.0 * rel_err, 0.05)
        if rel_resid > tolerance:
            misses.append(
                f"  t={t:>3} cm: B_mc={b_mc:.4e}, B_fit={b_fit:.4e}, "
                f"rel_resid={rel_resid:.1%} > tol={tolerance:.1%} "
                f"(rel_err={rel_err:.1%})"
            )
    assert not misses, (
        f"\n{quantity} fit misses anchors:\n" + "\n".join(misses)
    )


# ---------------------------------------------------------------------------
# Property 2: total-fit consistency with sum of component fits.
#   B_total(t) = B_neutron(t) + B_secondary_photon(t)
# holds exactly for the MC inputs (B_t = mc_t / pk_n = (mc_n + mc_p) / pk_n
# = B_n + B_p), so the fit predictions should preserve it.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("t", [25, 50, 75, 100, 125, 150, 175, 200, 250])
def test_total_fit_equals_sum_of_components(fit, t):
    n = fit.interpolate(quantity=_NEUTRON, thickness=t, warn=False).value
    p = fit.interpolate(quantity=_PHOTON, thickness=t, warn=False).value
    tot = fit.interpolate(quantity=_TOTAL, thickness=t, warn=False).value
    expected = n + p
    rel_diff = abs(tot - expected) / max(abs(expected), 1e-30)
    assert rel_diff < 0.01, (
        f"At t={t} cm: B_total_fit = {tot:.4f} but B_n_fit + B_p_fit = "
        f"{expected:.4f} ({rel_diff:.1%} apart). MC inputs satisfy "
        "B_t = B_n + B_p exactly so the fits should too."
    )


# ---------------------------------------------------------------------------
# Property 3: extrapolation past the deepest anchor stays physically
# bounded - positive and not blown up to runaway values. Past anchors
# B_total continues to grow because pk_neutron decays exponentially in
# the denominator while mc_photon stays bounded; an independent MC scan
# at t=350 cm gives B_p ~= 82, so the fit can legitimately reach
# triple-digit B at 4 m. The bound here only catches genuine numerical
# blow-up, not physical late-stage growth.
# ---------------------------------------------------------------------------
def test_extrapolation_is_bounded(fit):
    for t in (300, 350, 400):
        b = fit.interpolate(quantity=_TOTAL, thickness=t, warn=False).value
        assert 0.0 < b < 1000.0, (
            f"Extrapolated B_total({t} cm) = {b:.3f} is outside the "
            "physical-plausibility window (0, 1000) for polyethylene."
        )
