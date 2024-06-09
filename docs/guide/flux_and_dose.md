# Flux and dose

There are three levels of calculation, each building on the previous:

| Function | What it computes | Considers geometry? |
|---|---|---|
| `calculate_transmission` | Material attenuation only (exp(-Sigma*t)) | No |
| `calculate_flux` | S / (4*pi*R^2) * exp(-Sigma*t) * B | Yes (inverse square law) |
| `calculate_dose` | Flux * ICRP-116 dose coefficient | Yes (inverse square law) + dose conversion |

- **Transmission fraction**: pure material effect, no source strength or distance.
- **Uncollided flux**: includes source strength (particles/s), inverse-square-law spreading over the total distance, and an optional build-up factor.
- **Dose rate**: converts flux to effective dose rate (Sv/hr) using ICRP Publication 116 conversion coefficients for a specified irradiation geometry.

## Uncollided flux

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    pkc.Layer(thickness=1000),                    # 10 m void
    pkc.Layer(thickness=10, material=iron),         # 10 cm iron
]

source = pkc.Source("photon", 662e3)
result = pkc.calculate_flux(1e12, layers, source)
print(f"Flux: {result.uncollided_flux:.4e} photons/cm2/s")
print(f"Transmission: {result.transmission_fraction:.4e}")
print(f"Optical thickness: {result.optical_thickness:.3f}")
print(f"Build-up factor: {result.buildup_factor:.3f}")
print(f"Distance: {result.total_distance_cm:.1f} cm")
```

The returned `CalcResult` object has these properties:

- `uncollided_flux` -- flux at the outer surface (particles/cm2/s)
- `transmission_fraction` -- exp(-Sigma*t)
- `optical_thickness` -- Sum(Sigma_r,i * t_i)
- `buildup_factor` -- B (1.0 if no build-up model given)
- `total_distance_cm` -- total distance from source to detector

## Dose rate

Dose rate adds an ICRP-116 fluence-to-effective-dose conversion. You must specify the irradiation geometry.

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    pkc.Layer(thickness=1000),                    # 10 m void
    pkc.Layer(thickness=10, material=iron),         # 10 cm iron
]

source = pkc.Source("photon", 662e3)
result = pkc.calculate_dose(1e12, layers, source, "AP")
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

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    pkc.Layer(thickness=1000),
    pkc.Layer(thickness=10, material=iron),
]

source = pkc.Source("photon", 662e3)
for geo in ["AP", "PA", "RLAT", "LLAT", "ROT", "ISO"]:
    result = pkc.calculate_dose(1e12, layers, source, geo)
    print(f"{geo:>4s}: {result.dose_rate:.4e} Sv/hr")
```

## Neutron dose rate

The same interface applies to neutrons -- just create a neutron `Source`:

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    pkc.Layer(thickness=1000),
    pkc.Layer(thickness=10, material=iron),
]

source = pkc.Source("neutron", 14.1e6)
result = pkc.calculate_dose(1e12, layers, source, "AP")
print(f"Neutron dose rate: {result.dose_rate:.4e} Sv/hr")
```
