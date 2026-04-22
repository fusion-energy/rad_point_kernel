# Neutron flux

Three levels of calculation, each building on the previous:

| Function | What it computes | Considers geometry? |
|---|---|---|
| `calculate_transmission` | Material attenuation only (exp(-Sigma*t)) | No |
| **`calculate_flux`** | **S / (4*pi*R^2) * exp(-Sigma*t) * B** | **Yes (inverse square law)** |
| `calculate_dose` | Flux * ICRP-116 dose coefficient | Yes (inverse square law) |

The same `calculate_flux` interface applies to neutrons; just create a neutron `Source`.

## Example

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),               # 10 m void
    rpk.Layer(thickness=10, material=iron),  # 10 cm iron
]

source = rpk.Source(particle="neutron", energy=14.1e6)
result = rpk.calculate_flux(source_strength=1e12, layers=layers, source=source)
print(f"Flux: {result.uncollided_flux} neutrons/cm2/s")
print(f"Transmission: {result.transmission_fraction}")
print(f"Optical thickness: {result.optical_thickness}")
print(f"Distance: {result.total_distance_cm} cm")
```

The returned `CalcResult` has the same properties as in the photon case; see [Flux (photon)](photon_flux.md#example).

## With a manual build-up factor

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),
    rpk.Layer(thickness=10, material=iron),
]

source = rpk.Source(particle="neutron", energy=14.1e6)

B_flux = 3.0
result = rpk.calculate_flux(
    source_strength=1e12,
    layers=layers,
    source=source,
    buildup=rpk.BuildupModel.constant(B_flux),
)
print(f"Flux (B={B_flux}): {result.uncollided_flux} neutrons/cm2/s")
print(f"Applied build-up:  {result.buildup_factor}")
```

See [Calculate build-up with MC](buildup_mc.md) for how to compute B from first principles.
