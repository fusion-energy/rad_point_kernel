# Point Kernel Calculator

A fast engineering tool for neutron and photon shielding calculations.

## What it does

Point Kernel Calculator estimates radiation dose and flux behind shielding using the point-kernel method. It handles:

- **Photon sources** — monoenergetic (Cs-137) or multi-line spectra (Co-60)
- **Neutron sources** — monoenergetic (D-T) or mixed spectra (D-T + D-D)
- **Secondary photons** — gammas produced by neutron interactions in the shield (coupled neutron-photon transport via OpenMC)
- **Multi-layer shields** — any combination of materials and void
- **Build-up factor correction** — Monte Carlo computed, with Gaussian Process extrapolation and uncertainty
- **ICRP-116 dose coefficients** — 6 irradiation geometries (AP, PA, RLAT, LLAT, ROT, ISO)

## When to use it

The point-kernel method is a fast approximation for shielding design. Use it when:

- You need quick dose estimates for shield thickness scans
- You're doing preliminary design before committing to full Monte Carlo
- You want to screen many material/geometry combinations rapidly
- You need an analytical baseline to compare against transport codes

It is **not** a replacement for Monte Carlo transport codes (OpenMC, MCNP) for final design — but it gets you 80% of the answer in seconds instead of hours.

## How it works

1. **Point kernel** (Rust, instant): Computes uncollided flux through concentric spherical layers using removal cross sections from ENDF/B-VIII.0
2. **Build-up correction** (OpenMC, minutes): Runs Monte Carlo on a few thin shields to determine how much scattered radiation adds to the uncollided estimate
3. **GP extrapolation** ([inference-tools](https://github.com/C-bowman/inference-tools), instant): Fits a Gaussian Process to the Monte Carlo build-up factors and predicts B at any thickness with uncertainty bounds

## Quick start

```python
import rad_point_kernel as pkc

# Define a material and source
concrete = pkc.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3,
    fraction="mass",
)
source = pkc.Source("neutron", 14.1e6)  # D-T fusion

# Define the shield geometry
layers = [
    pkc.Layer(thickness=1000),                        # 10 m void
    pkc.Layer(thickness=100, material=concrete),       # 1 m concrete
]

# Calculate neutron dose rate
result = pkc.calculate_dose(
    source_strength=1e12,
    layers=layers,
    source=source,
    geometry="AP",
)
print(f"Dose rate: {result.dose_rate:.3e} Sv/hr")
```

## Installation

```bash
pip install rad_point_kernel
```

For Monte Carlo build-up factors, also install OpenMC separately:

```bash
python -m pip install --extra-index-url https://shimwell.github.io/wheels openmc
```

See the [installation guide](guide/installation.md) for details on OpenMC and cross section data.
