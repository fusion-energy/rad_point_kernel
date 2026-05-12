# Harima / ANS-6.4.3 buildup verification

## What is being verified

The point-kernel dose buildup factor predicted by `rad_point_kernel`
(Monte Carlo anchor + Shin-Ishii fit) is compared against the
**Geometric Progression (GP)** fits tabulated by Harima et al. for a
**point isotropic photon source in an infinite homogeneous medium**.

Reference: Harima, Sakamoto, Tanaka & Kawai, *Validity of the
Geometric Progression Formula in Approximating Gamma-Ray Buildup
Factors*, **Nuclear Science and Engineering 94** (1986). The same form
is reproduced in ANSI/ANS-6.4.3 (1991, revised 2010).

This is a **verification** case, not validation. The Harima/ANS GP
coefficients are themselves derived from moments-method, discrete
ordinates, or Monte Carlo calculations; agreement with them shows the
two computational chains agree, not that either matches measurement.

## GP form

For a point isotropic source in an infinite homogeneous medium, the
dose buildup factor at optical thickness $\mu t$ is

$$
B(\mu t) = 1 + (b - 1) \frac{K^{\mu t} - 1}{K - 1}, \quad K \neq 1
$$

with

$$
K(\mu t) = c\,(\mu t)^{a} + d\,\frac{\tanh(\mu t / \xi - 2) - \tanh(-2)}{1 - \tanh(-2)}
$$

The five coefficients $(a, b, c, d, \xi)$ are tabulated per element /
compound and per source energy in Harima 1986.

## Geometry match

The GP fit is for an **infinite homogeneous medium**. This script's
MC anchor uses a spherical shell of material with a point isotropic
source at the origin and the dose tally on the outer surface. As long
as the sphere is large enough that the tally radius is well inside the
material (no leakage at the tally point), this geometry is equivalent
to the GP assumption. Slab geometry would be a different quantity and
is **not** used here.

## Coverage

Harima 1986 covers (in their Tables I–VI):

- 23 elements: Be, B, C, N, O, Na, Mg, Al, Si, P, S, Ar, K, Ca, Fe,
  Cu, Mo, Sn, La, Gd, W, Pb, U.
- 3 compounds: air, water, ordinary concrete.
- 25 photon source energies from 15 keV to 15 MeV.
- 0.5 to 40 mean free paths in optical thickness.

This script compares against the full set when the coefficient table
is populated. The `--smoke` flag restricts to the first
`(material, energy)` and a 3-point mfp grid for a quick sanity check.

## Coefficient source file

The GP coefficients are not bundled in source — they live in

```
verification_and_validation/data/harima_1986_gp_coefficients.json
```

with a documented schema. Populating it (from Harima 1986 Tables I–VI
or equivalent reproduction in ANS-6.4.3) is a one-time manual step.
Until at least one entry is filled in, the script reports "table
empty" and exits 0.

## Tolerance

A comparison passes if either

- $|z| \leq 5.0$ (the MC-uncertainty-normalised difference), **or**
- relative error $\leq 5\%$.

The looser of the two is used because Harima's GP fit itself has
~1–5% residuals against the source moments-method buildup factors, so
even a perfect kernel cannot agree more tightly than that.

## Running

```bash
# Smoke test: first material/energy only, 3 mfp values.
python verification_and_validation/harima_buildup_comparison.py --smoke

# Full grid (uses cached MC anchors when present).
python verification_and_validation/harima_buildup_comparison.py

# Re-evaluate the GP <-> cached-MC comparison without re-running MC.
python verification_and_validation/harima_buildup_comparison.py --no-mc
```

Outputs:

- `verification_and_validation/harima_buildup_comparison.csv`
- `verification_and_validation/harima_buildup_comparison.png`
- `verification_and_validation/harima_buildup_comparison_cache.json`
  (cached MC anchors; safe to delete if you want fresh runs).

## Status

Scaffolded; coefficient table not yet populated. Once Harima 1986
Tables I–VI are entered into the JSON file, this case will exercise
all 23 elements + air/water/concrete at 25 energies × 10 mfp values.
