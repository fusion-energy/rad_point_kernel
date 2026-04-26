# Point Kernel Calculator

A fast engineering tool for neutron and photon shielding calculations.

## What it does

Point Kernel Calculator estimates radiation dose and flux behind shielding using the point-kernel method. It handles:

- **Photon sources** - monoenergetic (e.g. Cs-137) or multi-line spectra (e.g. Co-60)
- **Neutron sources** - monoenergetic (DT) or mixed spectra (DT + DD)
- **Secondary photons** - gammas produced by neutron interactions in the shield (coupled neutron-photon transport)
- **Multi-layer shields** - any combination of materials and voids
- **Build-up factor correction** - Monte Carlo computed (OpenMC), with Gaussian Process extrapolation and uncertainty
- **ICRP-116 dose coefficients** - 6 irradiation directions (AP, PA, RLAT, LLAT, ROT, ISO)

## When to use it

The point-kernel method is a fast approximation for shielding design. Use it when:

- You need quick dose estimates for shield thickness scans
- You're doing preliminary design before committing to full Monte Carlo
- You want to screen many material/geometry combinations rapidly
- You need an analytical baseline to compare against transport codes

It is **not** a replacement for Monte Carlo transport codes for final design, but it gets you closer to an answer in seconds instead of hours.

## How it works

1. **Point kernel** (Rust, instant): Computes uncollided flux through concentric spherical layers using removal cross sections and secondary photon production cross sections generated from ENDF/B-VIII.0
2. **Build-up correction** (OpenMC): Runs Monte Carlo on a few thin shields to determine how much scattered radiation adds to the uncollided estimate
3. **GP extrapolation** ([inference-tools](https://github.com/C-bowman/inference-tools), instant): Fits a Gaussian Process to the Monte Carlo build-up factors and predicts B at any thickness with uncertainty bounds

## Quick start

```python
import rad_point_kernel as rpk

# Define a material and source
concrete = rpk.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3,
    fraction="mass",
)
source = rpk.Source("neutron", 14.1e6)  # D-T fusion

# Define the shield geometry
layers = [
    rpk.Layer(thickness=1000),                    # 10 m void
    rpk.Layer(thickness=100, material=concrete),  # 1 m concrete
]

# Calculate neutron dose rate
result = rpk.calculate_dose(
    source_strength=1e12,
    layers=layers,
    source=source,
    geometry="AP",
)
print(f"Dose rate: {result.dose_rate:.3e} Sv/hr")
```

See the [installation guide](guide/installation.md) for setup, including OpenMC and cross-section data.
