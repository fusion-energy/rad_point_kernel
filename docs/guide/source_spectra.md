# Source spectra

The `Source` object accepts either a single energy value or a list of (energy, weight) tuples for multi-line sources. Weights are normalized internally, so their scale does not matter.

## Single energy line

A monoenergetic source is just a float in eV:

```python
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]

# Cs-137: single 662 keV gamma
source = rpk.Source("photon", 662e3)
frac = rpk.calculate_transmission(layers, source)
```

## Two-line photon source (Co-60)

Co-60 emits two gamma lines of roughly equal intensity at 1.173 MeV and 1.333 MeV. Pass them as a list of `(energy_eV, weight)` tuples:

```python
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]

source = rpk.Source("photon", [(1.173e6, 1.0), (1.333e6, 1.0)])
frac = rpk.calculate_transmission(layers, source)
print(f"Co-60 transmission through 10 cm Fe: {frac:.4e}")
```

## Mixed neutron source (D-T + D-D)

A fusion source with 95% D-T neutrons (14.06 MeV) and 5% D-D neutrons (2.45 MeV):

```python
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]

source = rpk.Source("neutron", [(14.06e6, 95), (2.45e6, 5)])
frac = rpk.calculate_transmission(layers, source)
print(f"D-T + D-D transmission through 10 cm Fe: {frac:.4e}")
```

## How weighting works

Given a spectrum of (E_i, w_i) pairs, the weighted transmission fraction is:

    T = Sum(w_i * T(E_i)) / Sum(w_i)

The weights `(95, 5)` and `(0.95, 0.05)` give identical results. This applies to all calculation functions -- transmission, flux, and dose.
