# Many thicknesses and extrapolation

Running Monte Carlo at every thickness is expensive. A more practical approach is to run Monte Carlo at a few thin shields, then fit the build-up factor and predict it at all thicknesses.

## Why analytical-form fits?

`BuildupFit` fits the Monte Carlo data points with a closed-form parametric expression chosen to match the physics:

- **1D (single material)**: a Shin-Ishii / Taylor double-exponential
  `B(μt) = A · exp(-α₁·μt) + (1-A) · exp(-α₂·μt)`. Three free parameters,
  with `B(0) = 1` baked in. Captures growth (heavy multiplying materials),
  decay below 1 (hydrogenous materials), peak-and-decay (beryllium), and
  dip-and-recover all with the same form.
- **Multi-D (multi-layer geometries)**: a thin-plate-spline radial-basis-function
  with degree-1 polynomial augmentation. No hyperparameters, works on
  scattered points.

Both forms are weighted by Monte Carlo statistical uncertainty (the fit gives more weight to points with smaller MC error bars).

## Basic example

Run Monte Carlo at 4 thicknesses, then fit and extrapolate to 10:

```python exec="true" source="material-block" result="text" session="extrapolation"
import rad_point_kernel as rpk

concrete = rpk.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)

PARTICLES_PER_HOUR = 1e12 * 3600  # 1e12 photons/sec activity, scale to /hour
source = rpk.Source(particle="photon", energy=1e6)

# Step 1: Monte Carlo at 4 thicknesses
mc_thicknesses = [5, 10, 15, 20]
mc_geometries = [
    [rpk.Layer(thickness=t, material=concrete)]
    for t in mc_thicknesses
]

mc_results = rpk.compute_buildup(
    geometries=mc_geometries,
    source=source,
    quantities=["dose-AP-photon"],
)

for t, r in zip(mc_thicknesses, mc_results):
    print(f"  {t:>2d} cm: B = {r.buildup['dose-AP-photon']}")

# Step 2: Fit
fit = rpk.BuildupFit(
    points=[{"thickness": t} for t in mc_thicknesses],
    results=mc_results,
)

# Step 3: Predict build-up at any thickness
all_thicknesses = mc_thicknesses + [30, 50, 75, 100, 150, 200]
for t in all_thicknesses:
    layers = [rpk.Layer(thickness=t, material=concrete)]
    bi = fit.interpolate(thickness=t)

    result = rpk.calculate_dose(
        layers=layers,
        source=source,
        geometry="AP",
        buildup=bi,
    ).scale(strength=PARTICLES_PER_HOUR)

    status = "EXTRAPOLATED" if bi.is_extrapolated else "interpolated"
    print(f"  {t} cm: dose = {result.dose} Sv/hr, "
          f"B = {bi.value} ({status})")
```

## InterpolationResult

`fit.interpolate()` returns an `InterpolationResult` with:

- `value` - predicted build-up factor
- `is_extrapolated` - True if the query point is outside the range of Monte Carlo data
- `extrapolated_axes` - dict showing which axes are extrapolated and by how much
- `sigma` - **not exposed** (NaN). Predictive uncertainty is on the roadmap; see the project TODO.

## Using InterpolationResult as build-up

You can pass an `InterpolationResult` directly to any calculation function:

```python exec="true" source="material-block" result="text" session="extrapolation"
bi = fit.interpolate(thickness=50)
result = rpk.calculate_dose(
    layers=[rpk.Layer(thickness=50, material=concrete)],
    source=source,
    geometry="AP",
    buildup=bi,
).scale(strength=PARTICLES_PER_HOUR)
print(f"Dose at 50 cm concrete: {result.dose} Sv/hr (B = {bi.value})")
```

## Using your own build-up values

You don't have to use `BuildupFit` - if you have build-up factors from another source (literature, your own fitting, a different interpolation method), you can apply them directly:

```python exec="true" source="material-block" result="text" session="extrapolation"
# Your own B value from any source
my_B = 3.2
layers = [rpk.Layer(thickness=50, material=concrete)]
result = rpk.calculate_dose(
    layers=layers,
    source=source,
    geometry="AP",
    buildup=rpk.BuildupModel.constant(my_B),
).scale(strength=PARTICLES_PER_HOUR)
print(f"Dose with B={my_B}: {result.dose} Sv/hr")
```

This means you can use any interpolation or fitting method you prefer (scipy, scikit-learn, a lookup table, or even a hand-drawn curve) and feed the result into the point-kernel calculation.

## Multi-dimensional fits

`BuildupFit` supports N-dimensional parameter spaces. For multi-layer geometries it switches automatically to a thin-plate-spline RBF interpolator. Example with water and concrete thicknesses as two axes:

```python exec="true" source="material-block" result="text" session="multidim"
import rad_point_kernel as rpk

water = rpk.Material(composition={"H2O": 1.0}, density=1.0)
concrete = rpk.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)

# Monte Carlo at a grid of (water, concrete) thicknesses
mc_water = [0, 10, 20]
mc_conc = [10, 20, 30]
points = []
geometries = []
for w in mc_water:
    for c in mc_conc:
        points.append({"water": w, "conc": c})
        layers = [rpk.Layer(thickness=1000)]
        if w > 0:
            layers.append(rpk.Layer(thickness=w, material=water))
        layers.append(rpk.Layer(thickness=c, material=concrete))
        geometries.append(layers)

source = rpk.Source(particle="neutron", energy=14.06e6)
mc_results = rpk.compute_buildup(
    geometries=geometries,
    source=source,
    quantities=["dose-AP-neutron"],
)

fit = rpk.BuildupFit(points=points, results=mc_results)

# Query at any (water, concrete) combination
bi = fit.interpolate(water=15, conc=25, quantity="dose-AP-neutron")
print(f"B = {bi.value}")
```

## Plotting build-up factors

```python exec="true" source="material-block" html="true" session="extrapolation"
from io import StringIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# After running Monte Carlo and fitting (as above)
thicknesses = np.linspace(5, 200, 100)
b_values = [fit.interpolate(thickness=float(t), warn=False).value for t in thicknesses]

fig, ax = plt.subplots()
ax.plot(thicknesses, b_values, "b-", label="Shin-Ishii fit")

# Monte Carlo points with statistical error bars on B
# sigma_B = mc_std_dev / pk (the same per-point sigma that weights the fit)
mc_bs = [r.buildup["dose-AP-photon"] for r in mc_results]
mc_b_err = [r.mc_std_dev["dose-AP-photon"] / r.pk["dose-AP-photon"] for r in mc_results]
ax.errorbar(mc_thicknesses, mc_bs, yerr=mc_b_err,
            fmt="ko", markersize=7, capsize=3, label="Monte Carlo")

ax.set_xlabel("Thickness (cm)")
ax.set_ylabel("Build-up factor B")
ax.legend()

buf = StringIO()
fig.savefig(buf, format="svg")
plt.close(fig)
print(buf.getvalue())
```

## Plotting dose

```python exec="true" source="material-block" html="true" session="extrapolation"
from io import StringIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

doses = []

for t in all_thicknesses:
    layers = [rpk.Layer(thickness=t, material=concrete)]
    bi = fit.interpolate(thickness=t, warn=False)
    pk = rpk.calculate_dose(
        layers=layers,
        source=source,
        geometry="AP",
    ).scale(strength=PARTICLES_PER_HOUR)
    doses.append(pk.dose * bi.value)

fig, ax = plt.subplots()
ax.plot(all_thicknesses, doses, "b-", label="PK with build-up")

# Monte Carlo reference points
mc_scaled = [r.scale(strength=PARTICLES_PER_HOUR) for r in mc_results]
mc_doses = [r.mc["dose-AP-photon"] for r in mc_scaled]
ax.plot(mc_thicknesses, mc_doses, "ko", markersize=7, label="Monte Carlo")

ax.set_xlabel("Thickness (cm)")
ax.set_ylabel("Dose rate (Sv/hr)")
ax.set_yscale("log")
ax.legend()

buf = StringIO()
fig.savefig(buf, format="svg")
plt.close(fig)
print(buf.getvalue())
```
