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

## Validation status

| Case | Reference | Status |
|------|-----------|--------|
| Photon transmission vs shielding experiments | SINBAD measurements | **Planned** |
| Neutron transmission vs shielding experiments | ORNL / SINBAD benchmarks | **Planned** |

Note: comparison to published buildup-factor tables (ANS-6.4.3 / Harima)
is *verification*, not validation, because those tables are themselves
derived from calculation (moments method, discrete ordinates, MC fits).
They belong in the verification table above once added.

## Running the scripts

All verification scripts are self-contained and runnable:

```bash
python verification_and_validation/uncollided_flux_in_void.py
python verification_and_validation/transmission_in_void.py
python verification_and_validation/beer_lambert_single_slab.py
```

Each prints a table of cases and exits with code 0 if all checks pass,
non-zero otherwise. They are suitable for running in CI.
