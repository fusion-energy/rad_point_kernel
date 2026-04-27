# Photon dose

Three levels of calculation, each building on the previous:

| Function | What it computes | Considers geometry? |
|---|---|---|
| `calculate_transmission` | Material attenuation only (exp(-Sigma*t)) | No |
| `calculate_flux` | 1 / (4*pi*R^2) * exp(-Sigma*t) * B per source particle | Yes (inverse square law) |
| **`calculate_dose`** | **Flux * ICRP-116 dose coefficient (per source particle)** | **Yes (inverse square law)** |

`calculate_dose` adds an ICRP-116 fluence-to-effective-dose conversion on top of `calculate_flux`. You must specify the irradiation geometry. The result is per source particle; apply an absolute strength via `result.scale(strength=...)`. If your strength is in particles/sec the resulting dose is in Sv/s; multiply the strength by 3600 to land in Sv/hr; for a pulsed source pass particles/shot to get Sv/shot.

## Example

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),               # 10 m void
    rpk.Layer(thickness=10, material=iron),  # 10 cm iron
]

source = rpk.Source(particle="photon", energy=662e3)
# Continuous source: 1e12 photons/sec activity, scale by photons per hour
result = rpk.calculate_dose(layers=layers, source=source, geometry="AP").scale(
    strength=1e12 * 3600,
)
print(f"Dose rate: {result.dose} Sv/hr")
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
    result = rpk.calculate_dose(layers=layers, source=source, geometry=geo).scale(
        strength=1e12 * 3600,
    )
    print(f"{geo}: {result.dose} Sv/hr")
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
    layers=layers,
    source=source,
    geometry="AP",
    buildup=rpk.BuildupModel.constant(B_dose),
).scale(strength=1e12 * 3600)
print(f"Dose rate (B={B_dose}): {result.dose} Sv/hr")
print(f"Applied build-up:     {result.buildup_factor}")
```
