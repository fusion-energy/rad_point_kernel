# point_kernel_calculator

Point-kernel neutron and photon shielding calculator.

Fast engineering estimates for radiation dose and flux behind multi-layer shields, with Monte Carlo build-up factor correction and Gaussian Process extrapolation.

## Installation

```bash
pip install point_kernel_calculator
```

For Monte Carlo build-up factors, also install OpenMC:

```bash
python -m pip install --extra-index-url https://shimwell.github.io/wheels openmc
```

## Quick start

```python
import point_kernel_calculator as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [pkc.Layer(thickness=10, material=iron)]
source = pkc.Source("neutron", 14.1e6)

result = pkc.calculate_dose(1e12, layers, source, "AP")
print(f"Dose rate: {result.dose_rate:.3e} Sv/hr")
```

## Documentation

Full documentation at [https://fusion-energy.github.io/point_kernel_calculator](https://fusion-energy.github.io/point_kernel_calculator)
