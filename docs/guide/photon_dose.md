# Photon dose

Three levels of calculation, each building on the previous:

| Function | What it computes | Considers geometry? |
|---|---|---|
| `calculate_transmission` | Material attenuation only (exp(-Sigma*t)) | No |
| `calculate_flux` | S / (4*pi*R^2) * exp(-Sigma*t) * B | Yes (inverse square law) |
| **`calculate_dose`** | **Flux * ICRP-116 dose coefficient** | **Yes (inverse square law)** |

`calculate_dose` adds an ICRP-116 fluence-to-effective-dose conversion on top of `calculate_flux`. You must specify the irradiation geometry.

## Example

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),               # 10 m void
    rpk.Layer(thickness=10, material=iron),  # 10 cm iron
]

source = rpk.Source(particle="photon", energy=662e3)
result = rpk.calculate_dose(source_strength=1e12, layers=layers, source=source, geometry="AP")
print(f"Dose rate: {result.dose_rate} Sv/hr")
```

## Irradiation geometries

Six geometries from ICRP-116 are supported:

| Geometry | Description |
|---|---|
| `"AP"` | Anterior-posterior (front exposure) |
| `"PA"` | Posterior-anterior (back exposure) |
| `"RLAT"` | Right lateral |
| `"LLAT"` | Left lateral |
| `"ROT"` | Rotational (averaged over rotation) |
| `"ISO"` | Isotropic (averaged over all directions) |

AP is the most conservative for most scenarios. Compare all six:

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),
    rpk.Layer(thickness=10, material=iron),
]

source = rpk.Source(particle="photon", energy=662e3)
for geo in ["AP", "PA", "RLAT", "LLAT", "ROT", "ISO"]:
    result = rpk.calculate_dose(source_strength=1e12, layers=layers, source=source, geometry=geo)
    print(f"{geo:>4s}: {result.dose_rate} Sv/hr")
```

## With a manual build-up factor

If you already have a dose build-up factor from tabulated data or a prior run, pass it via `BuildupModel.constant`. See [Calculate build-up with MC](buildup_mc.md) for how to compute B yourself.

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),
    rpk.Layer(thickness=10, material=iron),
]

source = rpk.Source(particle="photon", energy=662e3)

B_dose = 1.8
result = rpk.calculate_dose(
    source_strength=1e12,
    layers=layers,
    source=source,
    geometry="AP",
    buildup=rpk.BuildupModel.constant(B_dose),
)
print(f"Dose rate (B={B_dose}): {result.dose_rate} Sv/hr")
print(f"Applied build-up:     {result.buildup_factor}")
```
