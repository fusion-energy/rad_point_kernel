# Neutron flux

Three levels of calculation, each building on the previous:

| Function | What it computes | Considers geometry? |
|---|---|---|
| `calculate_transmission` | Material attenuation only (exp(-Sigma*t)) | No |
| **`calculate_flux`** | **1 / (4*pi*R^2) * exp(-Sigma*t) * B per source particle** | **Yes (inverse square law)** |
| `calculate_dose` | Flux * ICRP-116 dose coefficient | Yes (inverse square law) |

The same `calculate_flux` interface applies to neutrons; just create a neutron `Source`. The result is per source particle - apply an absolute strength via `result.scale(strength=...)`. Pulsed neutron devices (e.g. an ICF burn) take a strength in neutrons per shot; steady-state generators take particles per second (or per hour, if you prefer the result in /hr).

## Example

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),               # 10 m void
    rpk.Layer(thickness=10, material=iron),  # 10 cm iron
]

source = rpk.Source(particle="neutron", energy=14.1e6)
# Pulsed DT shot: strength in neutrons per shot, result in /shot
result = rpk.calculate_flux(layers=layers, source=source).scale(strength=1e16)
print(f"Flux: {result.flux} neutrons/cm2/shot")
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
    layers=layers,
    source=source,
    buildup=rpk.BuildupModel.constant(B_flux),
).scale(strength=1e16)
print(f"Flux (B={B_flux}): {result.flux} neutrons/cm2/shot")
print(f"Applied build-up:  {result.buildup_factor}")
```

See [Calculate build-up with MC](buildup_mc.md) for how to compute B from first principles.
