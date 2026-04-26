# Beer-Lambert single slab

## What is being verified

For an uncollided beam traversing a homogeneous slab of thickness $t$, the
transmission obeys Beer-Lambert's law:

$$
T(t) = \exp(-\Sigma\, t)
$$

where $\Sigma$ is the macroscopic (removal) cross section. This
verification does **not** assume any particular value of $\Sigma$ — it
checks the exponential law itself, which holds regardless of which
cross-section library is loaded.

## What the script does

[`verification_and_validation/beer_lambert_single_slab.py`](https://github.com/shimwell/rad_point_kernel/blob/main/verification_and_validation/beer_lambert_single_slab.py)
runs four independent sub-checks over 7 material/source combinations
(iron, water, lead at a range of neutron and photon energies):

### 1. Exponential scaling

$$T(k \cdot t_\text{ref}) = T(t_\text{ref})^k$$

Equivalently: $-\log T$ must be linear in $t$ with a non-negative slope.
Checked for $k \in \{0.5, 1, 2, 3, 5, 10\}$ against a reference thickness
$t_\text{ref} = 5$ cm.

### 2. Multiplicativity of stacked layers

$$T(t_1 + t_2) = T(t_1) \cdot T(t_2)$$

Checked both as a single layer of thickness $t_1 + t_2$ and as two layers
of thicknesses $t_1, t_2$ stacked in one call to `calculate_transmission`.
The stacked form also verifies that multi-layer optical thicknesses are
accumulated correctly.

### 3. Zero-thickness limit

$$T(0) = 1$$

for a material layer (not a void). This protects against a material
lookup being invoked with $t = 0$ and returning something subtly non-unity
(e.g. 0.9999999...).

### 4. Strict monotonicity

$T(t)$ must decrease strictly as $t$ increases. Checked across 9
thicknesses from 0 to 100 cm.

## Tolerance

Relative error must be $\leq 10^{-10}$. In practice all observed errors
are at the floating-point roundoff level ($\sim 10^{-16}$).

## Result

All 84 checks pass. Because sub-check #1 is cross-section-library
independent, this verification catches regressions in the kernel's
exponential math without being sensitive to ENDF/B version changes or
cross-section updates.
