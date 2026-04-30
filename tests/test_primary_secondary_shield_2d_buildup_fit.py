"""Regression test: 2D BuildupFit must keep B >= 0 in extrapolation.

Polyethylene + Magnetite + 3% rebar igloo, 14.06 MeV neutrons, 25 m
radial budget. The 2D thin-plate-spline RBF previously extrapolated
linearly past the conc anchor at 150 cm, producing B < 0 around
conc ~= 220 cm because the degree-1 polynomial term in TPS dominates
past the convex hull and the least-squares slope is negative on
attenuating data.

Fix: log-space TPS RBF (`rad_point_kernel_core` 4.1.0+). Fitting
`ln(B)` and exp()'ing on `predict` makes `B = exp(R) >= 0` by
construction; the linear extrapolation in fit-space becomes
exponential decay in B-space - the correct asymptote for shielding
attenuation.

This test pins the post-fix behaviour: the fit must stay non-negative
everywhere within the (poly, conc) plotting range used by the
bulk-shielding study, and reproduce its training anchors within MC
error.

Anchor data and probe grid are from the downstream study cache
(suggestion(1).md section 4.5; ENDF/B-VIII.1, weight windows on,
MAX_BATCHES=12000).
"""

from __future__ import annotations

import pytest
import rad_point_kernel as rpk


_NEUTRON = "dose-AP-neutron"
_PHOTON = "dose-AP-coupled-photon"
_TOTAL = "dose-AP-total"


# (poly_cm, conc_cm) anchors and per-quantity MC values per source
# neutron. PK reference is the uncollided neutron point-kernel at the
# 25 m outer surface for the matching geometry. Real MC tally output
# from the bulk-shielding study cache.
_ANCHORS = [
    # poly,conc,  mc_n,         mc_p,         mc_t,         pk_n,         std_n,        std_p,        std_t
    (  0,  25, 1.366896e-18, 2.549312e-20, 1.392389e-18, 8.764754e-19, 1.969670e-20, 9.752700e-22, 1.972083e-20),
    (  0,  50, 1.723278e-19, 6.343617e-21, 1.786714e-19, 1.219195e-19, 7.267110e-21, 1.760459e-22, 7.269242e-21),
    (  0,  75, 1.632809e-20, 9.696655e-22, 1.729776e-20, 1.695925e-20, 8.138559e-22, 3.712701e-23, 8.147023e-22),
    (  0, 100, 1.598695e-21, 1.430462e-22, 1.741741e-21, 2.359065e-21, 7.852058e-23, 5.094261e-24, 7.868566e-23),
    (  0, 150, 1.117367e-23, 1.897817e-24, 1.307148e-23, 4.564642e-23, 7.554787e-25, 6.264698e-26, 7.580718e-25),
    ( 50,  25, 4.253729e-20, 1.958825e-21, 4.449611e-20, 6.311904e-20, 1.925090e-21, 7.711352e-23, 1.926634e-21),
    ( 50,  50, 4.628619e-21, 2.484210e-22, 4.877040e-21, 8.779985e-21, 2.267100e-22, 1.031996e-23, 2.269448e-22),
    ( 50,  75, 4.970328e-22, 3.577767e-23, 5.328105e-22, 1.221314e-21, 2.477125e-23, 1.147887e-24, 2.479783e-23),
    ( 50, 100, 3.668232e-23, 3.951842e-24, 4.063417e-23, 1.698872e-22, 1.833830e-24, 1.085534e-25, 1.837041e-24),
    ( 50, 150, 2.578733e-25, 3.840407e-26, 2.962774e-25, 3.287209e-24, 7.267985e-26, 5.090343e-27, 7.285789e-26),
    (100,  25, 6.597048e-22, 2.029799e-22, 8.626847e-22, 4.545493e-21, 3.295300e-23, 5.440600e-24, 3.339910e-23),
    (100,  50, 6.774132e-23, 1.338236e-23, 8.112367e-23, 6.322873e-22, 3.174215e-24, 6.670793e-25, 3.243553e-24),
    (100,  75, 5.178051e-24, 1.103079e-24, 6.281130e-24, 8.795243e-23, 3.411050e-25, 8.346684e-26, 3.511685e-25),
    (100, 100, 6.023368e-25, 1.298293e-25, 7.321661e-25, 1.223436e-23, 1.691982e-25, 2.999113e-26, 1.718357e-25),
    (100, 150, 2.023671e-27, 3.967923e-28, 2.420463e-27, 2.367271e-25, 1.039691e-27, 1.655151e-28, 1.052783e-27),
]


def _build_result(poly, conc, mc_n, mc_p, mc_t, pk_n, std_n, std_p, std_t):
    pk_t = pk_n  # PK reference is the uncollided neutron PK for both
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
    points = [{"poly": p, "conc": c} for (p, c, *_) in _ANCHORS]
    results = [rpk.BuildupResult.from_dict(_build_result(*row)) for row in _ANCHORS]
    for r in results:
        if hasattr(r, "synthesize_dose_totals"):
            r.synthesize_dose_totals()
    return rpk.BuildupFit(points=points, results=results)


# Probe points past the deepest concrete anchor (150 cm). Expected to
# all be > 0 after the log-space TPS fix (rad_point_kernel_core 4.1+);
# under earlier linear-space TPS these went negative around conc ~= 220.
@pytest.mark.parametrize("poly, conc", [
    (0, 175), (0, 200), (0, 250), (0, 300), (0, 400),
    (50, 175), (50, 200), (50, 250), (50, 300), (50, 400),
    (100, 175), (100, 200), (100, 250), (100, 300), (100, 400),
])
def test_2d_buildup_stays_nonneg_in_extrapolation(fit, poly, conc):
    b = fit.interpolate(quantity=_TOTAL, poly=poly, conc=conc, warn=False).value
    assert b > 0.0, (
        f"B(poly={poly}, conc={conc}) = {b:.4e}; should stay positive "
        "in extrapolation. Linear-space TPS RBF goes negative past "
        "conc ~= 220 cm because the degree-1 polynomial term "
        "extrapolates linearly. Log-space TPS fixes this."
    )


def test_2d_buildup_within_3sigma_at_anchors(fit):
    """The fit should reproduce its training anchors within MC error."""
    misses = []
    for (poly, conc, _, _, mc_t, pk_n, _std_n, _std_p, std_t) in _ANCHORS:
        b_mc = mc_t / pk_n
        rel = std_t / mc_t if mc_t > 0 else 0.05
        b_fit = fit.interpolate(quantity=_TOTAL, poly=poly, conc=conc, warn=False).value
        rel_resid = abs(b_fit - b_mc) / b_mc
        tol = max(3.0 * rel, 0.05)
        if rel_resid > tol:
            misses.append(
                f"  ({poly}, {conc}): B_mc={b_mc:.4e}, B_fit={b_fit:.4e}, "
                f"rel_resid={rel_resid:.1%} > tol={tol:.1%}"
            )
    assert not misses, "\n" + "\n".join(misses)
