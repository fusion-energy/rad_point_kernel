# Transmission in a void

## What is being verified

With no attenuating medium the optical thickness is zero, so the uncollided
transmission fraction is exactly 1:

$$
T = \exp\!\left(-\sum_i \Sigma_{r,i}\, t_i\right) = \exp(0) = 1
$$

for any particle type, energy, or geometry. This is checked via both
`calculate_transmission` (returns the scalar directly) and `calculate_flux`
(which exposes `.transmission_fraction` on its result).

## What the script does

[`verification_and_validation/transmission_in_void.py`](https://github.com/shimwell/rad_point_kernel/blob/main/verification_and_validation/transmission_in_void.py)
covers:

- **Thicknesses**: 0, 1, 10, 100, 1000, 10000 cm
- **Sources**: neutron at 14.06 MeV / 2.45 MeV / thermal (0.0253 eV);
  photon at 1 MeV / 662 keV / 10 keV
- **Multi-layer stacks**: four different void stacks with 1 to 10 layers

The tolerance is **absolute zero** — $T$ must equal $1.0$ exactly, not
"within $10^{-15}$". Any departure from 1.0 would indicate either a stray
cross-section lookup on a void layer or an accumulation of floating-point
error where it shouldn't exist.

## Result

All 40 cases pass. This guards against regressions such as treating a void
as "air" by default, returning $T \approx 1$ instead of $T = 1$ for
zero-thickness layers, or incorrectly summing optical thicknesses across
stacks.
