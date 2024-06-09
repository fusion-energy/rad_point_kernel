# Build-up factors with Monte Carlo

The point-kernel method computes uncollided flux -- the fraction of particles that travel through the shield without scattering. In reality, scattered particles also contribute to dose behind the shield. The build-up factor B corrects for this:

    Corrected = Point-kernel * B

`compute_buildup` runs OpenMC Monte Carlo simulations on a list of geometries and returns the ratio of Monte Carlo result to point-kernel result for each.

!!! note
    This requires the `[mc]` optional dependency: `pip install rad_point_kernel[mc]`

## Basic example

Compute the photon dose build-up factor for 10 cm of iron:

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [pkc.Layer(thickness=10, material=iron)]

source = pkc.Source("photon", 1e6)
results = pkc.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP"],
)

r = results[0]
print(f"MC dose:     {r.mc['dose-AP']:.4e} (per source particle)")
print(f"PK dose:     {r.pk['dose-AP']:.4e}")
print(f"Build-up B:  {r.buildup['dose-AP']:.3f}")
```

## Multiple geometries at once

Pass a list of layer lists to run Monte Carlo on several shield configurations in one call:

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
thicknesses = [5, 10, 15, 20]
geometries = [[pkc.Layer(thickness=t, material=iron)] for t in thicknesses]

source = pkc.Source("photon", 1e6)
results = pkc.compute_buildup(
    geometries=geometries,
    source=source,
    quantities=["dose-AP"],
)

for t, r in zip(thicknesses, results):
    print(f"{t:>2d} cm: B = {r.buildup['dose-AP']:.3f}")
```

## Quantity strings

The `quantities` argument specifies what to tally. Since the `Source` object already carries the particle type, quantity names don't need to repeat it — `"flux"` means "flux of whatever particle the source emits" and `"dose-AP"` means "dose from the source particle at AP geometry."

The only exception is when you want to tally **secondary photons** produced by neutron interactions. In that case, append `-coupled-photon` to indicate you want the secondary photon contribution rather than the primary particle:

| Quantity string | Description |
|---|---|
| `"flux"` | Flux of the source particle |
| `"dose-AP"` | Dose from the source particle, AP geometry |
| `"dose-ISO"` | Dose from the source particle, ISO geometry |
| `"flux-coupled-photon"` | Secondary photon flux (neutron source only) |
| `"dose-AP-coupled-photon"` | Secondary photon dose, AP geometry (neutron source only) |

You can request multiple quantities in a single call:

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [pkc.Layer(thickness=10, material=iron)]

source = pkc.Source("neutron", 14.1e6)
results = pkc.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP", "flux"],
)

r = results[0]
print(f"Dose B: {r.buildup['dose-AP']:.3f}")
print(f"Flux B: {r.buildup['flux']:.3f}")
```

## BuildupResult

Each element in the returned list is a `BuildupResult` with these dictionaries, keyed by quantity string:

- `r.mc` -- Monte Carlo tally value (per source particle)
- `r.mc_std_dev` -- Monte Carlo standard deviation
- `r.pk` -- point-kernel reference value
- `r.buildup` -- build-up factor (mc / pk)

It also has `r.optical_thickness` (the total optical thickness of the geometry).

## Applying the build-up factor

### Pass BuildupResult directly

The simplest approach. The `calculate_*` functions auto-detect the correct quantity key:

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [pkc.Layer(thickness=10, material=iron)]
SOURCE = 1e12

source = pkc.Source("photon", 1e6)
results = pkc.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP"],
)
r = results[0]

corrected = pkc.calculate_dose(
    SOURCE, layers, source, "AP",
    buildup=r,
)
print(f"Dose with build-up: {corrected.dose_rate:.4e} Sv/hr")
```

### Use BuildupModel.constant() manually

If you have a known build-up factor from another source:

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [pkc.Layer(thickness=10, material=iron)]

# Apply a fixed buildup factor of 2.5 — this multiplies the PK dose
# by 2.5 regardless of thickness. Useful when you have B from an
# external source (literature, previous calculation, etc.)
buildup = pkc.BuildupModel.constant(2.5)

source = pkc.Source("photon", 1e6)
result = pkc.calculate_dose(
    1e12, layers, source, "AP",
    buildup=buildup,
)
print(f"Dose with B=2.5: {result.dose_rate:.4e} Sv/hr")
```

## Cross sections path

If the `OPENMC_CROSS_SECTIONS` environment variable is not set, pass the path directly:

```python
source = pkc.Source("photon", 1e6)
results = pkc.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP"],
    cross_sections="/path/to/cross_sections.xml",
)
```

## Secondary photons (coupled transport)

When using a neutron source, nuclear reactions in the shield produce secondary gamma rays. To include these, add `"dose-AP-coupled-photon"` to your quantities. The Monte Carlo runs with coupled neutron-photon transport automatically:

```python
source = pkc.Source("neutron", 14.1e6)
results = pkc.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP", "dose-AP-coupled-photon"],
)

r = results[0]
S = 1e12
total_dose = (r.mc["dose-AP"] + r.mc["dose-AP-coupled-photon"]) * S
```

Since you're already running Monte Carlo for neutron build-up, coupled transport adds minimal extra cost and gives the secondary photon dose for free.

## Simulation parameters

`compute_buildup` accepts optional parameters to control the Monte Carlo simulation:

| Parameter | Default | Description |
|---|---|---|
| `particles_per_batch` | 10,000 | Particles per batch |
| `max_batches` | 100 | Maximum number of batches |
| `trigger_rel_err` | 0.05 | Target relative error on tallies |
