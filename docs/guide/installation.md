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

OpenMC needs nuclear cross section data. There are two ways to provide it:

**Environment variable:**

```bash
export OPENMC_CROSS_SECTIONS=/path/to/cross_sections.xml
```

**Pass the path directly:**

```python
results = pkc.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP"],
    cross_sections="~/nuclear_data/cross_sections.xml",
)
```

The `~` is expanded automatically. If you pass a directory instead of a file, `/cross_sections.xml` is appended.
