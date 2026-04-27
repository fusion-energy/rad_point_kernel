# API reference

All public classes and functions are accessible from the top-level `rad_point_kernel` module.

```python exec="true" source="material-block" session="api"
import rad_point_kernel as rpk
```

---

## Source

### `Source(particle, energy)`

Defines a radiation source by particle type and energy.

- `particle` - `"neutron"` or `"photon"`.
- `energy` - energy in eV as a single float, or a list of `(energy_eV, weight)` tuples for multi-line spectra.

**Examples:**

```python exec="true" source="material-block" session="api"
# Monoenergetic Cs-137 photon source
source = rpk.Source(particle="photon", energy=662e3)

# Co-60 two-line photon source
source = rpk.Source(particle="photon", energy=[(1173e3, 1.0), (1333e3, 1.0)])

# D-T fusion neutron source
source = rpk.Source(particle="neutron", energy=14.1e6)
```

---

## Materials

### `Material(composition, density, fraction="mass")`

A material defined by composition and density.

- `composition` - dict of component names to fractions. Keys can be element symbols (`"Fe"`), chemical formulas (`"H2O"`), or nuclide names (`"Li6"`).
- `density` - density in g/cm3.
- `fraction` - `"mass"` (default) or `"atom"`.

**Properties:** `density`, `composition`, `fraction`.

**Methods:**

- `nuclide_mass_fractions()` - returns a dict of resolved nuclide-level mass fractions.

### `Material.volume_mix(mat_a, frac_a, mat_b, frac_b)`

Static method. Mix two materials by volume fraction. Returns a new `Material` with volume-weighted density and composition.

---

## Geometry

### `Layer(thickness, material=None)`

A spherical shell layer with thickness in cm and optional material.

- `thickness` - layer thickness in cm.
- `material` - a `Material`, or omit for void (empty space).

**Properties:** `thickness`, `has_material`, `material`.

---

## Calculation functions

All calculation functions take a `Source` object that specifies the particle type and energy. They return a result normalised **per source particle**. Apply an absolute strength via `result.scale(strength=...)` to land in the unit you want: a strength in particles/sec gives Sv/s and particles/cm^2/s; in particles/hr gives Sv/hr; in particles/shot gives Sv/shot. Re-scaling is a numeric post-process - it never re-runs the simulation.

### Transmission fraction

#### `calculate_transmission(layers, source)`

Returns a float in [0, 1]: the pure material attenuation exp(-Sum(Sigma_r * t)) without geometric spreading.

- `layers` - list of `Layer` objects.
- `source` - a `Source` object.

### Uncollided flux

#### `calculate_flux(layers, source, buildup=None)`

Returns a `CalcResult` with flux per source particle. Computes 1 / (4*pi*R^2) * exp(-Sigma*t) * B.

- `layers` - list of `Layer` objects.
- `source` - a `Source` object.
- `buildup` - optional `BuildupModel`, `BuildupResult`, or `InterpolationResult`.

### Dose rate

#### `calculate_dose(layers, source, geometry, buildup=None)`

Returns a `CalcResult` with dose per source particle (Sv/particle). Multiplies the per-particle uncollided flux by the ICRP-116 fluence-to-dose conversion coefficient.

- `layers` - list of `Layer` objects.
- `source` - a `Source` object.
- `geometry` - one of `"AP"`, `"PA"`, `"RLAT"`, `"LLAT"`, `"ROT"`, `"ISO"`.
- `buildup` - optional `BuildupModel`, `BuildupResult`, or `InterpolationResult`.

### Secondary photon dose rate

#### `calculate_secondary_photon_dose_rate(layers, source, geometry, neutron_buildup=None)`

Returns a `SecondaryGammaResult` with dose per source particle. Analytical estimate of the secondary photon dose from neutron capture gammas in each material, plus the neutron dose. For more accurate results, prefer a coupled neutron-photon `compute_buildup` run.

- `layers` - list of `Layer` objects.
- `source` - a monoenergetic neutron `Source` object.
- `geometry` - one of `"AP"`, `"PA"`, `"RLAT"`, `"LLAT"`, `"ROT"`, `"ISO"`.
- `neutron_buildup` - optional build-up applied to the neutron component.

---

## Result types

### `CalcResult`

Returned by flux and dose calculations. Numeric fields are normalised per source particle; multiply through with `.scale(strength=...)` to apply an absolute source strength.

**Properties:**

- `uncollided_flux` - flux at the outer surface (per source particle until scaled).
- `dose_rate` - effective dose (Sv per source particle until scaled). Present for dose calculations only.
- `transmission_fraction` - exp(-Sigma*t).
- `optical_thickness` - Sum(Sigma_r,i * t_i).
- `buildup_factor` - applied build-up factor (1.0 if none).
- `total_distance_cm` - total distance from source to detector (cm).
- `source_strength` - strength baked into `uncollided_flux` and `dose_rate` (1.0 = unscaled).

**Methods:**

- `scale(strength)` - returns a new `CalcResult` with `uncollided_flux` and `dose_rate` multiplied to reflect the given absolute strength. Re-scaling replaces the previous strength rather than compounding. Strength must be > 0.

### `SecondaryGammaResult`

Returned by `calculate_secondary_photon_dose_rate`. Dose fields are per source particle until scaled.

**Properties:**

- `neutron_dose_rate` - neutron dose (per source particle until scaled).
- `secondary_photon_dose_rate` - secondary photon dose (per source particle until scaled).
- `total_dose_rate` - sum of neutron + secondary photon dose.
- `source_strength` - strength baked into the dose fields (1.0 = unscaled).

**Methods:**

- `scale(strength)` - returns a new result with the dose fields multiplied. Same semantics as `CalcResult.scale`.

---

## Build-up

### `BuildupModel`

Analytical build-up factor models. Created via static methods:

- `BuildupModel.none()` - no correction (B = 1).
- `BuildupModel.constant(value)` - fixed scalar B.
- `BuildupModel.linear(a)` - B = 1 + a * mu_t.
- `BuildupModel.polynomial(coeffs)` - B = 1 + a1*(mu_t) + a2*(mu_t)^2 + ...
- `BuildupModel.taylor(a, alpha1, alpha2)` - B = A*exp(-alpha1*mu_t) + (1-A)*exp(-alpha2*mu_t).

### `BuildupResult`

Result of a single Monte Carlo build-up computation from `compute_buildup`. The `mc`, `mc_std_dev`, and `pk` values are stored per source particle until scaled.

**Attributes (dicts keyed by quantity string):**

- `mc` - Monte Carlo tally value (per source particle until scaled).
- `mc_std_dev` - Monte Carlo standard deviation.
- `pk` - point-kernel reference value.
- `buildup` - build-up factor (mc / pk; dimensionless ratio, never affected by scaling).
- `optical_thickness` - total optical thickness (float).
- `source_strength` - strength baked into `mc`, `mc_std_dev`, `pk` (1.0 = unscaled).

Can be passed directly to any calculation function as the `buildup` argument.

**Methods:**

- `scale(strength)` - returns a new `BuildupResult` with `mc`, `mc_std_dev`, and `pk` multiplied to reflect the given absolute strength. Re-scaling replaces the previous strength.
- `to_dict()` / `from_dict(d)` - serialize to/from plain dict (the `source_strength` round-trips).
- `BuildupResult.save(results, path)` - save a list of results to JSON.
- `BuildupResult.load(path)` - load a list of results from JSON.

### `BuildupTable(points, results)`

GP-based interpolation table for build-up factors.

- `points` - list of dicts defining parameter values (e.g. `[{"thickness": 5}, {"thickness": 10}]`).
- `results` - list of `BuildupResult`, one per point.

**Properties:** `axis_names`, `available_quantities`, `axis_ranges`.

**Methods:**

- `interpolate(quantity=None, warn=True, **kwargs)` - returns an `InterpolationResult` at the given parameter values.

### `InterpolationResult`

Result of a `BuildupTable.interpolate()` query.

**Attributes:**

- `value` - predicted build-up factor.
- `sigma` - GP uncertainty (1-sigma standard deviation).
- `is_extrapolated` - True if the query is outside the range of Monte Carlo data.
- `extrapolated_axes` - dict of `{axis: (value, min, max)}` for extrapolated axes.

Can be passed directly to any calculation function as the `buildup` argument.

### `compute_buildup(geometries, source, quantities, ...)`

Run OpenMC Monte Carlo simulations and compute build-up factors. Requires OpenMC; see [Installation](guide/installation.md#openmc-for-monte-carlo-build-up-factors).

**Parameters:**

- `geometries` - list of layer lists, each defining a shield geometry.
- `source` - a `Source` object specifying the particle type and energy.
- `quantities` - list of quantity strings (e.g. `["flux", "dose-AP"]`).
- `particles_per_batch` - particles per batch (default 10,000).
- `max_batches` - maximum batches (default 100).
- `trigger_rel_err` - target relative error (default 0.05).
- `cross_sections` - path to cross_sections.xml or directory (default: uses `OPENMC_CROSS_SECTIONS` env var).

**Returns:** list of `BuildupResult`, one per geometry. MC and PK values are stored per source particle. Apply absolute scaling with `result.scale(strength=...)` to land in your chosen unit.
