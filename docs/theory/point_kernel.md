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

The removal cross section $\sigma_r$ is not the total cross section. It accounts for the fact that forward-scattered particles effectively continue in the beam direction and so shouldn't be counted as "removed":

$$\sigma_r(E) = \sigma_t(E) - f_{\text{fwd}}(E) \cdot \sigma_s(E)$$

where $\sigma_s$ is the scattering channel that dominates the forward peak (elastic for neutrons, coherent/Rayleigh for photons), and $f_{\text{fwd}}(E) \in [0, 1]$ is the fraction of that channel scattering into the forward cone $\mu \in [\mu_0, 1]$, with $\mu_0 = \cos\theta_0$. We use $\mu_0 = 0$, i.e. the forward hemisphere ($\theta < 90^\circ$) is treated as "not removed".

### Neutrons

$$\sigma_r(E) = \sigma_t(E) - f_{\text{fwd}}(E) \cdot \sigma_{\text{el}}(E)$$

$f_{\text{fwd}}$ is obtained by integrating the ENDF elastic angular distribution $p(\mu, E)$ from $\mu_0$ to 1:

$$f_{\text{fwd}}(E) = \int_{\mu_0}^{1} p(\mu, E)\, d\mu$$

If no angular distribution is available the scattering is treated as isotropic, giving $f_{\text{fwd}} = (1 - \mu_0)/2 = 0.5$ at $\mu_0 = 0$.

### Photons

$$\sigma_r(E) = \sigma_t(E) - f_{\text{fwd}}(E) \cdot \sigma_{\text{coh}}(E)$$

For photons the forward peak is coherent (Rayleigh) scattering, and $f_{\text{fwd}}$ is computed from the Thomson-weighted atomic form factor $F(x, Z)$:

$$f_{\text{fwd}}(E) = \frac{\displaystyle\int_{\mu_0}^{1} (1 + \mu^2)\,[F(x, Z)]^2\, d\mu}{\displaystyle\int_{-1}^{1} (1 + \mu^2)\,[F(x, Z)]^2\, d\mu}, \qquad x(E, \mu) = \frac{E \sqrt{(1 - \mu)/2}}{hc}$$

Incoherent (Compton) and photoelectric contributions are not subtracted because they genuinely remove the photon from the beam.

### Why not simply $\sigma_t - \sigma_s$?

Subtracting all of the elastic (or coherent) cross section would assume every scattered particle keeps going in the forward direction; this overestimates transmission. Subtracting nothing (just $\sigma_t$) assumes every scatter removes the particle; this underestimates it. Weighting by the angular distribution is between the two and, for anisotropic scatterers at high energy, much closer to reality.

### Compound materials

For a compound the macroscopic removal cross section is:

$$\Sigma_r = \rho \cdot N_A \cdot \sum_i \frac{w_i \cdot \sigma_{r,i}(E)}{A_i}$$

where $w_i$ is the mass fraction, $A_i$ the atomic mass, and $\rho$ the density. The microscopic removal cross sections $\sigma_{r,i}$ are pre-computed per nuclide from ENDF/B-VIII.0 using the [`endf`](https://github.com/shimwell/endf-python) package (see `examples/removal_xs_to_json.py` in that repo for the extraction script).

## Limitations

- **Uncollided only**: Without build-up factors, the method only counts particles that haven't interacted. For thick shields, this severely underestimates the dose.
- **Spherical geometry**: Layers are concentric spheres. Real geometries with ducts, penetrations, or non-spherical shapes need Monte Carlo.
- **No energy degradation**: Particles that scatter and lose energy are not tracked; the build-up factor corrects for this empirically.
