# Verification and Validation

This section documents what has been checked to give confidence that the
point-kernel kernel math is correct, and (separately) where its predictions
have been compared against external references.

We split the work using the standard definitions:

- **Verification** ("are we solving the equations right?"): compares code
  output to closed-form analytic solutions or to another code solving the
  same equations (code-to-code benchmarks).
- **Validation** ("are we solving the right equations?"): compares code
  output to **experimental measurements** (e.g. SINBAD, ICSBEP, or ORNL
  shielding benchmarks).

The runnable scripts live in the top-level
[`verification_and_validation/`](https://github.com/fusion-energy/rad_point_kernel/tree/main/verification_and_validation)
folder of the repository.

## Verification status

| Case | Quantity | Reference | Status |
|------|----------|-----------|--------|
| [Uncollided flux in void](uncollided_flux_in_void.md) | $\Phi(r) = S / (4\pi r^2)$ | Analytic | Passing, machine precision |
| [Transmission in void](transmission_in_void.md) | $T = 1$ | Analytic | Passing, exact |
| [Beer-Lambert single slab](beer_lambert_single_slab.md) | $T(t) = e^{-\Sigma t}$ | Analytic (exponential law) | Passing, $< 10^{-15}$ |
| [Neutron fit form quality](neutron_fit_quality.md) | $B_\mathrm{neutron}(\tau)$ across 19 neutron-source materials | MC anchor library (leave-one-out) | Mean holdout $z = 4.9$ |
| [Primary-photon fit form quality](primary_photon_fit_quality.md) | $B_\mathrm{photon}(\tau)$ across 13 photon-source materials | MC anchor library (leave-one-out) | Mean holdout $z = 2.1$ |
| [Secondary-photon fit form quality](secondary_photon_fit_quality.md) | $B_\mathrm{coupled-photon}(\tau)$ across 21 materials | MC anchor library (leave-one-out) | Mean holdout $z = 6.0$ |
| [Harima / ANS-6.4.3 buildup](harima_buildup_comparison.md) | $B_\mathrm{photon}(\mu t)$ across 23 elements + air/water/concrete | Harima 1986 GP fits | Scaffolded; coefficient table not yet populated |

## Validation status

No validation cases have been added yet — validation requires
comparison to **experimental measurements** (reactor / shielding
benchmarks), and those have not been wired into this section.
Published buildup-factor tables (ANS-6.4.3 / Harima) are themselves
*verification* artifacts (derived from moments method, discrete
ordinates, or MC fits) and would be added under the verification
section above if/when included.

## Running the scripts

All verification scripts are self-contained and runnable:

```bash
python verification_and_validation/uncollided_flux_in_void.py
python verification_and_validation/transmission_in_void.py
python verification_and_validation/beer_lambert_single_slab.py
python verification_and_validation/secondary_photon_fit_test_suite.py
python verification_and_validation/neutron_fit_test_suite.py
python verification_and_validation/primary_photon_fit_test_suite.py
python verification_and_validation/harima_buildup_comparison.py
```

Each prints a table of cases and exits with code 0 if all checks pass,
non-zero otherwise. They are suitable for running in CI. The
`secondary_photon_fit_test_suite.py` script doubles as a regression harness when
choosing or tuning a coupled-photon fit form: it loads
`secondary_photon_fit_cache.json` (21-material MC anchor library), fits the
candidate forms in `CANDIDATES`, scores each per-material on the
held-out anchor, and writes a CSV plus a multi-panel comparison plot.
