# Installation

## Basic install

```bash
pip install rad_point_kernel
```

This gives you everything for point-kernel calculations, build-up factor interpolation with Gaussian Process extrapolation ([inference-tools](https://github.com/C-bowman/inference-tools)), and plotting.

## OpenMC for Monte Carlo build-up factors

To compute build-up factors from Monte Carlo simulation, you also need OpenMC. OpenMC must be installed separately:

**Pre-built wheels (easiest):**

```bash
python -m pip install --extra-index-url https://shimwell.github.io/wheels openmc
```

**Or follow the official guide:**

See the [OpenMC installation guide](https://docs.openmc.org/en/stable/usersguide/install.html).

## OpenMC cross sections

OpenMC needs nuclear cross section data. If you don't already have it, one easy way to download ENDF/B-VIII.1 is:

```bash
python -m pip install openmc_data
download_endf -r b8.1
```

Other options and pre-packaged libraries are listed on [openmc.org](https://openmc.org/data/#endf-b-viii-1).

Once the data is on disk, point OpenMC at it one of two ways:

**Environment variable:**

Replace "/path/to/cross_sections.xml" with your actual path to where your cross_section.xml is located.

```bash
export OPENMC_CROSS_SECTIONS=/path/to/cross_sections.xml
```

**Pass the path directly:**

```python
results = rpk.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP"],
    cross_sections="/path/to/cross_sections.xml",
)
```
