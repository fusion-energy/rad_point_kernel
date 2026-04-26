# Neutron transmission

`calculate_transmission` computes the fraction of uncollided neutrons that pass through a set of layers. Like the photon version, it returns a float between 0 and 1 representing exp(-Sum(Sigma_r,i * t_i)).

## Single layer example

Calculate neutron transmission through 10 cm of iron at 14.1 MeV (D-T fusion):

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]

source = rpk.Source(particle="neutron", energy=14.1e6)
frac = rpk.calculate_transmission(layers=layers, source=source)
print(f"Transmission: {frac:.4e}")
```

## Comparing shield thicknesses

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
source = rpk.Source(particle="neutron", energy=14.1e6)

for thickness in [5, 10, 20, 50]:
    layers = [rpk.Layer(thickness=thickness, material=iron)]
    frac = rpk.calculate_transmission(layers=layers, source=source)
    print(f"{thickness:>3d} cm iron: {frac:.4e}")
```

## Notes

- Neutron removal cross sections come from ENDF/B-VIII.0.
- The calculation uses the removal cross-section approximation, which is most accurate for fast neutrons and deep penetration problems.
- For spectra with multiple energy lines, see [Source spectra](source_spectra.md).
