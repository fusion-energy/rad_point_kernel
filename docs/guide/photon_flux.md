# Photon flux

Three levels of calculation, each building on the previous:

| Function | What it computes | Considers geometry? |
|---|---|---|
| `calculate_transmission` | Material attenuation only (exp(-Sigma*t)) | No |
| **`calculate_flux`** | **1 / (4*pi*R^2) * exp(-Sigma*t) * B per source particle** | **Yes (inverse square law)** |
| `calculate_dose` | Flux * ICRP-116 dose coefficient | Yes (inverse square law) |

`calculate_flux` returns flux per source particle, with inverse-square-law spreading over the total distance and an optional build-up factor. Apply an absolute source strength via `result.scale(strength=...)`; the unit of the result follows the unit of the strength you supply (particles/sec gives /s, particles/shot gives /shot, etc.).

## Example

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),               # 10 m void
    rpk.Layer(thickness=10, material=iron),  # 10 cm iron
]

source = rpk.Source(particle="photon", energy=662e3)
result = rpk.calculate_flux(layers=layers, source=source).scale(strength=1e12)
print(f"Flux: {result.uncollided_flux} photons/cm2/s")
print(f"Transmission: {result.transmission_fraction}")
print(f"Optical thickness: {result.optical_thickness}")
print(f"Build-up factor: {result.buildup_factor}")
print(f"Distance: {result.total_distance_cm} cm")
```

The returned `CalcResult` object has these properties:

- `uncollided_flux` - flux at the outer surface; equals geometry * transmission * B (per source particle until you call `.scale`)
- `transmission_fraction` - exp(-Sigma*t)
- `optical_thickness` - Sum(Sigma_r,i * t_i)
- `buildup_factor` - B (1.0 if no build-up model given)
- `total_distance_cm` - total distance from source to detector
- `source_strength` - the strength the result has been scaled by (1.0 = unscaled)

## With a manual build-up factor

If you already have a build-up factor from tabulated data or a prior run, pass it via `BuildupModel.constant`. See [Calculate build-up with MC](buildup_mc.md) for how to compute B yourself.

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),
    rpk.Layer(thickness=10, material=iron),
]

source = rpk.Source(particle="photon", energy=662e3)

B_flux = 2.0
result = rpk.calculate_flux(
    layers=layers,
    source=source,
    buildup=rpk.BuildupModel.constant(B_flux),
).scale(strength=1e12)
print(f"Flux (B={B_flux}): {result.uncollided_flux} photons/cm2/s")
print(f"Applied build-up:  {result.buildup_factor}")
```
