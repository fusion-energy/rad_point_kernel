# Flux and dose

There are three levels of calculation, each building on the previous:

| Function | What it computes | Considers geometry? |
|---|---|---|
| `calculate_transmission` | Material attenuation only (exp(-Sigma*t)) | No |
| `calculate_flux` | S / (4*pi*R^2) * exp(-Sigma*t) * B | Yes (inverse square law) |
| `calculate_dose` | Flux * ICRP-116 dose coefficient | Yes (inverse square law) |

- **Transmission fraction**: pure material effect, no source strength or distance.
- **Uncollided flux**: includes source strength (particles/s), inverse-square-law spreading over the total distance, and an optional build-up factor.
- **Dose rate**: converts flux to effective dose rate (Sv/hr) using ICRP Publication 116 conversion coefficients for a specified irradiation geometry.

## Uncollided flux

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),               # 10 m void
    rpk.Layer(thickness=10, material=iron),  # 10 cm iron
]

source = rpk.Source(particle="photon", energy=662e3)
result = rpk.calculate_flux(source_strength=1e12, layers=layers, source=source)
print(f"Flux: {result.uncollided_flux:.4e} photons/cm2/s")
print(f"Transmission: {result.transmission_fraction:.4e}")
print(f"Optical thickness: {result.optical_thickness:.3f}")
print(f"Build-up factor: {result.buildup_factor:.3f}")
print(f"Distance: {result.total_distance_cm:.1f} cm")
```

The returned `CalcResult` object has these properties:

- `uncollided_flux` - flux at the outer surface (particles/cm2/s)
- `transmission_fraction` - exp(-Sigma*t)
- `optical_thickness` - Sum(Sigma_r,i * t_i)
- `buildup_factor` - B (1.0 if no build-up model given)
- `total_distance_cm` - total distance from source to detector

## Dose rate

Dose rate adds an ICRP-116 fluence-to-effective-dose conversion. You must specify the irradiation geometry.

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),               # 10 m void
    rpk.Layer(thickness=10, material=iron),  # 10 cm iron
]

source = rpk.Source(particle="photon", energy=662e3)
result = rpk.calculate_dose(source_strength=1e12, layers=layers, source=source, geometry="AP")
print(f"Dose rate: {result.dose_rate:.4e} Sv/hr")
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
    print(f"{geo:>4s}: {result.dose_rate:.4e} Sv/hr")
```

## Neutron dose rate

The same interface applies to neutrons. Just create a neutron `Source`:

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),
    rpk.Layer(thickness=10, material=iron),
]

source = rpk.Source(particle="neutron", energy=14.1e6)
result = rpk.calculate_dose(source_strength=1e12, layers=layers, source=source, geometry="AP")
print(f"Neutron dose rate: {result.dose_rate:.4e} Sv/hr")
```
