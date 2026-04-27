# Neutron dose

Three levels of calculation, each building on the previous:

| Function | What it computes | Considers geometry? |
|---|---|---|
| `calculate_transmission` | Material attenuation only (exp(-Sigma*t)) | No |
| `calculate_flux` | 1 / (4*pi*R^2) * exp(-Sigma*t) * B per source particle | Yes (inverse square law) |
| **`calculate_dose`** | **Flux * ICRP-116 dose coefficient (per source particle)** | **Yes (inverse square law)** |

The same `calculate_dose` interface applies to neutrons; just create a neutron `Source`. The [irradiation geometries](photon_dose.md#irradiation-geometries) are the same as for photons. As with photons, the dose result is per source particle - apply an absolute strength via `result.scale(strength=...)`. For pulsed neutron devices (e.g. an inertial-confinement burn) supply neutrons per shot to get Sv/shot.

## Example

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),               # 10 m void
    rpk.Layer(thickness=10, material=iron),  # 10 cm iron
]

source = rpk.Source(particle="neutron", energy=14.1e6)
# Pulsed DT shot: 1e16 neutrons per shot, result in Sv/shot
result = rpk.calculate_dose(layers=layers, source=source, geometry="AP").scale(
    strength=1e16,
)
print(f"Neutron dose: {result.dose} Sv/shot")
```

## With a manual build-up factor

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),
    rpk.Layer(thickness=10, material=iron),
]

source = rpk.Source(particle="neutron", energy=14.1e6)

B_dose = 2.5
result = rpk.calculate_dose(
    layers=layers,
    source=source,
    geometry="AP",
    buildup=rpk.BuildupModel.constant(B_dose),
).scale(strength=1e16)
print(f"Neutron dose (B={B_dose}): {result.dose} Sv/shot")
print(f"Applied build-up:        {result.buildup_factor}")
```

See [Calculate build-up with MC](buildup_mc.md) for how to compute B from first principles.

!!! note "Secondary photons"
    Fast neutrons generate secondary photons through inelastic scatter and capture. Those photons can contribute significantly to the dose, especially for high-Z shields. See [Dose (coupled)](coupled_dose.md) for how to include them via Monte Carlo.
