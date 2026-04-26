# Source spectra

The `Source` object accepts either a single energy value in eV or a list of `(energy_eV, weight)` tuples for multi-line sources. Weights are normalized internally, so their scale does not matter.

## Single energy line

A monoenergetic source is just a float in eV:

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]

# Cs-137: single 662 keV gamma
source = rpk.Source(particle="photon", energy=662e3)
frac = rpk.calculate_transmission(layers=layers, source=source)
print(f"Cs-137 transmission through 10 cm Fe: {frac}")
```

## Two-line photon source (Co-60)

Co-60 emits two gamma lines of roughly equal intensity at 1.173 MeV and 1.333 MeV. Pass them as a list of `(energy_eV, weight)` tuples:

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]

source = rpk.Source(particle="photon", energy=[(1.173e6, 1.0), (1.333e6, 1.0)])
frac = rpk.calculate_transmission(layers=layers, source=source)
print(f"Co-60 transmission through 10 cm Fe: {frac}")
```

## Mixed neutron source (D-T + D-D)

A fusion source with 95% D-T neutrons (14.06 MeV) and 5% D-D neutrons (2.45 MeV):

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]

source = rpk.Source(particle="neutron", energy=[(14.06e6, 95), (2.45e6, 5)])
frac = rpk.calculate_transmission(layers=layers, source=source)
print(f"D-T + D-D transmission through 10 cm Fe: {frac}")
```

## Weighting

Mixed-energy sources are weighted proportionally internally: weights are normalized, so `(95, 5)` and `(0.95, 0.05)` give identical results. This applies to all calculation functions: transmission, flux, and dose.
