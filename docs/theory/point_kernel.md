# The Point Kernel Method

## Overview

The point-kernel method calculates the uncollided radiation flux from a point source through layers of shielding material. It is the simplest analytical method for shielding calculations and forms the basis for engineering dose estimates.

## The formula

For a point isotropic source of strength $S$ (particles/sec) at the origin, the uncollided flux at distance $R$ through layers of material is:

$$\Phi(R) = \frac{S}{4\pi R^2} \cdot \exp\left(-\sum_i \Sigma_{r,i} \cdot t_i\right) \cdot B$$

Where:

- $\frac{S}{4\pi R^2}$ is the **geometric attenuation** (inverse square law)
- $\exp(-\sum_i \Sigma_{r,i} \cdot t_i)$ is the **material attenuation** (exponential decay through each layer)
- $\Sigma_{r,i}$ is the **macroscopic removal cross section** of layer $i$ (cm$^{-1}$)
- $t_i$ is the **thickness** of layer $i$ (cm)
- $B$ is the **build-up factor** (correction for scattered radiation, $B \geq 1$)

## Transmission fraction

The transmission fraction is the pure material attenuation without geometric spreading:

$$T = \exp\left(-\sum_i \Sigma_{r,i} \cdot t_i\right)$$

This is a dimensionless number between 0 and 1. It answers: "what fraction of particles pass through the shield without interacting?"

## Dose rate

To convert flux to dose rate, multiply by the ICRP-116 dose conversion coefficient $h(E)$:

$$\dot{D} = \Phi \cdot h(E) \cdot 3600 \cdot 10^{-12} \quad \text{[Sv/hr]}$$

Where $h(E)$ is in pSv$\cdot$cm$^2$ and depends on the particle energy and irradiation geometry (AP, PA, etc.).

## Removal cross sections

The removal cross section $\sigma_r$ is not the total cross section. It accounts for the fact that forward-scattered particles effectively continue in the beam direction:

$$\sigma_r(E) = \sigma_{\text{total}}(E) - f_{\text{forward}}(E) \cdot \sigma_{\text{elastic}}(E)$$

Where $f_{\text{forward}}$ is the fraction of elastic scattering into the forward hemisphere, computed from ENDF angular distribution data. This is more accurate than the naive $\sigma_{\text{total}} - \sigma_{\text{elastic}}$ which treats all elastic scattering as forward.

For compound materials, the macroscopic removal cross section is:

$$\Sigma_r = \rho \cdot N_A \cdot \sum_i \frac{w_i \cdot \sigma_{r,i}(E)}{A_i}$$

Where $w_i$ is the mass fraction, $A_i$ the atomic mass, and $\rho$ the density.

## Limitations

- **Uncollided only**: Without build-up factors, the method only counts particles that haven't interacted. For thick shields, this severely underestimates the dose.
- **Spherical geometry**: Layers are concentric spheres. Real geometries with ducts, penetrations, or non-spherical shapes need Monte Carlo.
- **No energy degradation**: Particles that scatter and lose energy are not tracked — the build-up factor corrects for this empirically.
