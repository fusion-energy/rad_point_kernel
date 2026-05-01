# Neutron fit form quality across 19 materials

Sibling of the [secondary-photon harness](secondary_photon_fit_quality.md).
Validates the shipped 1D neutron-component fit (Shin-Ishii double
exponential) against the same MC anchor library used for the secondary
photon V&V.

## Setup

Same cache (`verification_and_validation/secondary_photon_fit_cache.json`),
neutron source at 14.06 MeV. The script loads the `dose-AP` field per
material as $B_\mathrm{neutron}(\tau)$, fits the Shin-Ishii form

$$
B(\tau) = A \cdot e^{-\alpha_1 \tau} + (1 - A) \cdot e^{-\alpha_2 \tau}
$$

(`B(0) = 1` by construction), and scores per-material training and
holdout z-scores.

Two materials drop out of the harness because the deepest neutron
anchor underflowed to zero MC counts: `h2o_water` (training) and
`polyethylene_deep` (holdout at 350 cm of polyethylene). Shin-Ishii
isn't meaningful on zero-count anchors.

## Result (19 materials)

| Metric | Value |
|---|---|
| Worst holdout $z$ | 19.77 (magnetite) |
| Mean holdout $z$ | 4.85 |
| Median holdout $z$ | 1.77 |
| Worst training $z$ | 20.22 (magnetite) |

Magnetite is the same outlier as in the secondary-photon harness, for
the same reason: the holdout anchor sits ~50% past the deepest training
anchor in optical thickness, beyond what 6-anchor Shin-Ishii can
extrapolate. Outside that one material, mean holdout $z \approx 3$
across diverse moderators and heavy attenuators.

## Running

```bash
python verification_and_validation/neutron_fit_test_suite.py
```

Prints the per-material residual table sorted best-to-worst.

## Same machinery for flux

The neutron *flux* is also fit with Shin-Ishii in this library; the
harness just hard-codes `dose-AP` as the loaded field. To validate the
flux fit:

1. Regenerate the MC cache including flux:
   `compute_buildup(quantities=["flux-neutron"], ...)`. The `flux-neutron`
   field will land alongside `dose-AP` in each `BuildupResult`.
2. In `_load()`, swap the read of `r["mc"]["dose-AP"]` for
   `r["mc"]["flux-neutron"]` (and the same on `pk`, `mc_std_dev`).

The fit, scoring, and reporting code are identical because both
quantities use the same Shin-Ishii form and have $B(0) = 1$.
