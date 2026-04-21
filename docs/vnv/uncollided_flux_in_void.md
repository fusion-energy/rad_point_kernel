# Uncollided flux in a void

## What is being verified

For an isotropic point source of strength $S$ (particles/s) in a pure void,
the uncollided scalar flux at distance $r$ is, by geometric spreading alone:

$$
\Phi(r) = \frac{S}{4\pi r^2}
$$

With no attenuating medium the exponential attenuation factor is exactly 1,
so the point-kernel result must reduce to this analytic form to
floating-point precision — independent of particle type, source energy, or
cross-section library.

## What the script does

[`verification_and_validation/uncollided_flux_in_void.py`](https://github.com/shimwell/rad_point_kernel/blob/main/verification_and_validation/uncollided_flux_in_void.py)
sweeps:

- Distances from **1 cm to 5000 cm** (8 values, log-spaced)
- Four sources: **14.06 MeV neutron**, **2.45 MeV neutron**, **1 MeV photon**,
  **662 keV photon**

For each combination it calls `calculate_flux` through a single `Layer` of
void (no `material`) and compares `.uncollided_flux` against
$S / (4\pi r^2)$. It also asserts that `.transmission_fraction == 1.0`
exactly.

## Tolerance

Relative error must be $\leq 10^{-12}$. In practice all 32 cases return
0.00e+00 relative error (the division is bitwise identical to the analytic
formula because both reduce to the same floating-point operations).

## Result

All 32 cases pass. This is a trivial but load-bearing check — if it ever
fails it means the geometric factor has drifted (e.g. a stray density
factor, the wrong $4\pi$ constant, or incorrect handling of zero-thickness
materials).
