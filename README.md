# rad_point_kernel

Point-kernel neutron and photon shielding calculator.

Fast engineering estimates for radiation dose and flux behind multi-layer shields, with Monte Carlo build-up factor correction and Gaussian Process extrapolation.

## Installation

```bash
pip install rad_point_kernel
```

For Monte Carlo build-up factors, also install OpenMC:

```bash
python -m pip install --extra-index-url https://shimwell.github.io/wheels openmc
```

## Quick start

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [pkc.Layer(thickness=10, material=iron)]
source = pkc.Source("neutron", 14.1e6)

result = pkc.calculate_dose(1e12, layers, source, "AP")
print(f"Dose rate: {result.dose_rate:.3e} Sv/hr")
```

## Documentation

Full documentation at [https://fusion-energy.github.io/rad_point_kernel](https://fusion-energy.github.io/rad_point_kernel)
