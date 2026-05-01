# Secondary-photon fit form quality across 21 materials

This is a regression-testing harness rather than an analytic verification
script. It exists to answer the question "is my coupled-photon
`B_p(τ)` fit form good enough across the materials this library is
likely to be applied to?", and is the tool that was used to design the
shipped Form K_reg fit.

## Setup

The cache file `verification_and_validation/secondary_photon_fit_cache.json`
holds 6 training anchors plus 1 held-out anchor for each of 21
shielding materials, generated with `compute_buildup` against
ENDF/B-VIII.1 with weight windows on. Materials span:

- **Hydrogenous moderators**: polyethylene, polyethylene_deep,
  borated_poly, lih, d2o, h2o_water (pure water).
- **Mixed light/heavy composites**: WC + 20/50/70 vol% water,
  WB + 20/50/70 vol% water (peak-and-decay $B_p$).
- **Saturating heavy attenuators**: concrete, steel, ss316,
  tungsten, lead, magnetite (heavy concrete), b4c, aluminum,
  concrete_steel10.

Per material, the held-out anchor sits ~30-50% past the deepest
training anchor.

## What the harness reports

For each candidate fit form in `CANDIDATES`, for every material:

- **Training fit z-score** — `max |B_pred - B_mc| / σ` over the 6
  training anchors. Catches in-sample over-/under-fitting.
- **Holdout z-score** — `|B_pred - B_mc| / σ` at the held-out anchor.
  Catches over-fitting that hides in training but shows up in
  extrapolation.

Aggregated across the 21 materials:

| | Current shipped form (Power × sat, 3p) | Form K_reg (5p, L2) |
|---|---|---|
| Worst holdout z | 54.3 | 36.6 |
| Mean holdout z | 14.7 | 6.0 |
| Median holdout z | 7.1 | 2.4 |
| Worst training z | 42.4 | 3.3 |

K_reg dominates baseline on every metric. The wins come from the
late-stage exponential growth tail (`exp(d·x_n + e·x_n²)` with signed
`d, e`) which baseline's pure power-law saturation cannot represent.

## Running

```bash
python verification_and_validation/secondary_photon_fit_test_suite.py
```

The script also writes:

- `secondary_photon_fit_dataset.csv` — flat (material, τ, descriptors, B_p, σ) CSV
  ready for symbolic regression or external analysis (gitignored).
- `/tmp/secondary_photon_fit_suite_plot.png` — per-material comparison panels.

## Adding a new candidate form

Drop a `(fit_fn, eval_fn)` pair into the `CANDIDATES` dict in
`secondary_photon_fit_test_suite.py`. Re-run; the harness scores it across all
21 materials and prints the per-material residuals plus the worst /
mean / median holdout z.

## Same machinery for flux

The harness loads the `dose-AP-coupled-photon` field per material. The
secondary-photon flux is fit with the same Form K_reg in the shipped
library; to validate that path:

1. Regenerate the cache including the flux quantity:
   `compute_buildup(quantities=["dose-AP-coupled-photon", "flux-coupled-photon"], ...)`.
2. In `_load()`, swap the read of `r["mc"]["dose-AP-coupled-photon"]`
   for `r["mc"]["flux-coupled-photon"]` (and the matching `pk`,
   `mc_std_dev`).

Same applies for the [neutron](neutron_fit_quality.md) and
[primary-photon](primary_photon_fit_quality.md) sibling harnesses.
