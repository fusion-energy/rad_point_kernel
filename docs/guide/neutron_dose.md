# Neutron dose

Three levels of calculation, each building on the previous:

| Function | What it computes | Considers geometry? |
|---|---|---|
| `calculate_transmission` | Material attenuation only (exp(-Sigma*t)) | No |
| `calculate_flux` | S / (4*pi*R^2) * exp(-Sigma*t) * B | Yes (inverse square law) |
| **`calculate_dose`** | **Flux * ICRP-116 dose coefficient** | **Yes (inverse square law)** |

The same `calculate_dose` interface applies to neutrons; just create a neutron `Source`. The [irradiation geometries](photon_dose.md#irradiation-geometries) are the same as for photons.

## Example

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),               # 10 m void
    rpk.Layer(thickness=10, material=iron),  # 10 cm iron
]

source = rpk.Source(particle="neutron", energy=14.1e6)
result = rpk.calculate_dose(source_strength=1e12, layers=layers, source=source, geometry="AP")
print(f"Neutron dose rate: {result.dose_rate} Sv/hr")
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
    source_strength=1e12,
    layers=layers,
    source=source,
    geometry="AP",
    buildup=rpk.BuildupModel.constant(B_dose),
)
print(f"Neutron dose (B={B_dose}): {result.dose_rate} Sv/hr")
print(f"Applied build-up:        {result.buildup_factor}")
```

See [Calculate build-up with MC](buildup_mc.md) for how to compute B from first principles.

!!! note "Secondary photons"
    Fast neutrons generate secondary photons through inelastic scatter and capture. Those photons can contribute significantly to the dose, especially for high-Z shields. See [Dose (coupled)](coupled_dose.md) for how to include them via Monte Carlo.
