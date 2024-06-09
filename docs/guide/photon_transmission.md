# Photon transmission

`calculate_transmission` computes the fraction of uncollided photons that pass through a set of layers. It returns a float between 0 and 1 representing the pure material attenuation exp(-Sum(Sigma_r,i * t_i)), without geometric spreading.

## Single layer example

Calculate photon transmission through 10 cm of iron at 662 keV (Cs-137 gamma line):

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [pkc.Layer(thickness=10, material=iron)]

source = pkc.Source("photon", 662e3)
frac = pkc.calculate_transmission(layers, source)
print(f"Transmission: {frac:.4e}")
```

## Multi-layer example

Transmission through a composite shield:

```python
import rad_point_kernel as pkc

water = pkc.Material(composition={"H2O": 1.0}, density=1.0)
concrete = pkc.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)

layers = [
    pkc.Layer(thickness=30, material=water),
    pkc.Layer(thickness=100, material=concrete),
]

source = pkc.Source("photon", 662e3)
frac = pkc.calculate_transmission(layers, source)
print(f"Transmission: {frac:.4e}")
```

## Notes

- Void layers do not attenuate -- they contribute zero optical thickness.
- Two 5 cm iron layers give the same result as one 10 cm iron layer: exp(-Sigma * 5) * exp(-Sigma * 5) = exp(-Sigma * 10).
- For spectra with multiple energy lines, see [Source spectra](source_spectra.md).
