# Primary-photon fit form quality across 13 materials

Sibling of the [neutron](neutron_fit_quality.md) and
[secondary-photon](secondary_photon_fit_quality.md) harnesses. Validates
the shipped 1D primary-photon fit (Shin-Ishii double exponential) against
photon-source MC anchors.

## Setup

Same cache file as the other two harnesses, but loads only the
photon-source materials (`*_photon` keys: `polyethylene_photon`,
`concrete_photon`, ..., `tungsten_photon`). These come from
`compute_buildup` runs with a 1 MeV photon source — `dose-AP` is then
the *primary-photon* dose (no secondary radiation, no `coupled-photon`
component to worry about).

The fit form is Shin-Ishii:

$$
B(\tau) = A \cdot e^{-\alpha_1 \tau} + (1 - A) \cdot e^{-\alpha_2 \tau}
$$

with $B(0) = 1$ by construction.

## Result (13 photon-source materials)

| Metric | Value |
|---|---|
| Worst holdout $z$ | 5.50 (polyethylene) |
| Mean holdout $z$ | 2.08 |
| Median holdout $z$ | 1.07 |
| Worst training $z$ | 1.69 |

Cleaner than the neutron case because primary-photon physics has no
late-stage thermalisation tail; Shin-Ishii's two decaying exponentials
match the underlying physics well across all materials. Mean holdout
$z = 2.1$ vs $4.9$ for neutron and $6.0$ for secondary photon.

## Running

```bash
python verification_and_validation/primary_photon_fit_test_suite.py
```

## Same machinery for flux

Photon flux uses the same Shin-Ishii form. To validate it:

1. Regenerate the photon-source materials in the cache with
   `compute_buildup(quantities=["flux-photon"], source=photon_src, ...)`.
2. Swap the field read in `_load()` from `r["mc"]["dose-AP"]` to
   `r["mc"]["flux-photon"]` (and the matching `pk`, `mc_std_dev`).

The harness code is identical otherwise.
