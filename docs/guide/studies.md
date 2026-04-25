# Shielding studies

The `examples/python/` directory contains two study scripts that demonstrate complete workflows: parameter scans with caching, GP extrapolation, and automatic plot generation.

## Single-layer study

**Script:** `examples/python/single_layer_study.py`

Scans shield thickness from thin to thick for different materials, computing total dose (neutron + secondary photon) at each point.

**What it does:**

1. Defines materials (Portland concrete + 3% steel rebar, Magnetite concrete + 3% steel rebar) using `Material.volume_mix()`.
2. Runs coupled neutron-photon Monte Carlo (`dose-AP` + `dose-AP-coupled-photon`) at a few thin thicknesses (10--60 cm).
3. Caches all Monte Carlo results to JSON. On re-run, only simulates missing thicknesses.
4. Computes total build-up as B_total = MC_total / PK_neutron.
5. Uses `BuildupTable` to GP-extrapolate B_total out to 400 cm.
6. Produces a log-scale plot of total dose vs shield thickness, with Monte Carlo points, PK + build-up lines, and uncertainty bands for each material.

**Geometry:** 10 m void + single material layer, 14.1 MeV D-T source, AP dose.

## Multiple-layer study

**Script:** `examples/python/multiple_layer_study.py`

Scans two thickness parameters (water + concrete) simultaneously, demonstrating 2D GP extrapolation.

**What it does:**

1. Defines a geometry: 10 m void + variable water + variable concrete.
2. Runs coupled Monte Carlo at a grid of (water, concrete) thickness pairs.
3. Caches all results incrementally -- only simulates missing grid points.
4. Builds 1D `BuildupTable`s (one per fixed water thickness, one per fixed concrete thickness).
5. GP-extrapolates dose to a dense grid in both dimensions.
6. Produces two plots:
    - Total dose vs concrete thickness (one curve per water thickness)
    - Total dose vs water thickness (one curve per concrete thickness)

**Geometry:** 10 m void + water (0--50 cm) + concrete (10--400 cm), 14.06 MeV, AP dose.

## Running the studies

```bash
cd examples/python

# Single-layer study
python single_layer_study.py

# Multiple-layer study
python multiple_layer_study.py
```

Results and plots are saved to `results/single_layer/` and `results/multiple_layer/` respectively. Cached Monte Carlo data is saved alongside the plots, so subsequent runs are fast.
