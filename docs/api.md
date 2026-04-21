# API reference

All public classes and functions are accessible from the top-level `rad_point_kernel` module.

```python
import rad_point_kernel as rpk
```

---

## Source

### `Source(particle, energy)`

Defines a radiation source by particle type and energy.

- `particle` - `"neutron"` or `"photon"`.
- `energy` - energy in eV as a single float, or a list of `(energy_eV, weight)` tuples for multi-line spectra.

**Examples:**

```python
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

All calculation functions take a `Source` object that specifies the particle type and energy.

### Transmission fraction

#### `calculate_transmission(layers, source)`

Returns a float in [0, 1]: the pure material attenuation exp(-Sum(Sigma_r * t)) without geometric spreading.

- `layers` - list of `Layer` objects.
- `source` - a `Source` object.

### Uncollided flux

#### `calculate_flux(source_strength, layers, source, buildup=None)`

Returns a `CalcResult`. Computes S / (4*pi*R^2) * exp(-Sigma*t) * B.

- `source_strength` - source strength in particles/s.
- `layers` - list of `Layer` objects.
- `source` - a `Source` object.
- `buildup` - optional `BuildupModel`, `BuildupResult`, or `InterpolationResult`.

### Dose rate

#### `calculate_dose(source_strength, layers, source, geometry, buildup=None)`

Returns a `CalcResult`. Computes uncollided flux multiplied by the ICRP-116 fluence-to-dose conversion coefficient.

- `source_strength` - source strength in particles/s.
- `layers` - list of `Layer` objects.
- `source` - a `Source` object.
- `geometry` - one of `"AP"`, `"PA"`, `"RLAT"`, `"LLAT"`, `"ROT"`, `"ISO"`.
- `buildup` - optional `BuildupModel`, `BuildupResult`, or `InterpolationResult`.

### Secondary photon dose rate

#### `calculate_secondary_photon_dose_rate(source_strength, layers, source, geometry, neutron_buildup=None)`

Returns a `SecondaryGammaResult`. Analytical estimate of the secondary photon dose rate from neutron capture gammas in each material, plus the neutron dose. For more accurate results, prefer a coupled neutron-photon `compute_buildup` run.

- `source_strength` - source strength in particles/s.
- `layers` - list of `Layer` objects.
- `source` - a monoenergetic neutron `Source` object.
- `geometry` - one of `"AP"`, `"PA"`, `"RLAT"`, `"LLAT"`, `"ROT"`, `"ISO"`.
- `neutron_buildup` - optional build-up applied to the neutron component.

---

## Result types

### `CalcResult`

Returned by flux and dose calculations.

**Properties:**

- `uncollided_flux` - flux at the outer surface (particles/cm2/s).
- `dose_rate` - effective dose rate in Sv/hr (present for dose calculations).
- `transmission_fraction` - exp(-Sigma*t).
- `optical_thickness` - Sum(Sigma_r,i * t_i).
- `buildup_factor` - applied build-up factor (1.0 if none).
- `total_distance_cm` - total distance from source to detector (cm).

### `SecondaryGammaResult`

Returned by `calculate_secondary_photon_dose_rate`.

**Properties:**

- `neutron_dose_rate` - neutron dose rate in Sv/hr.
- `secondary_photon_dose_rate` - secondary photon dose rate in Sv/hr.
- `total_dose_rate` - sum of neutron + secondary photon dose.

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

Result of a single Monte Carlo build-up computation from `compute_buildup`.

**Attributes (dicts keyed by quantity string):**

- `mc` - Monte Carlo tally value (per source particle).
- `mc_std_dev` - Monte Carlo standard deviation.
- `pk` - point-kernel reference value.
- `buildup` - build-up factor (mc / pk).
- `optical_thickness` - total optical thickness (float).

Can be passed directly to any calculation function as the `buildup` argument.

**Methods:**

- `to_dict()` / `from_dict(d)` - serialize to/from plain dict.
- `BuildupResult.save(results, path)` - save a list of results to JSON.
- `BuildupResult.load(path)` - load a list of results from JSON.

### `BuildupTable(points, results)`

GP-based interpolation table for build-up factors. Requires the `[gp]` optional dependency.

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

Run OpenMC Monte Carlo simulations and compute build-up factors. Requires the `[mc]` optional dependency.

**Parameters:**

- `geometries` - list of layer lists, each defining a shield geometry.
- `source` - a `Source` object specifying the particle type and energy.
- `quantities` - list of quantity strings (e.g. `["flux", "dose-AP"]`).
- `particles_per_batch` - particles per batch (default 10,000).
- `max_batches` - maximum batches (default 100).
- `trigger_rel_err` - target relative error (default 0.05).
- `cross_sections` - path to cross_sections.xml or directory (default: uses `OPENMC_CROSS_SECTIONS` env var).

**Returns:** list of `BuildupResult`, one per geometry.
