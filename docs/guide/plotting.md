# Plotting results

This page shows how to create a dose-vs-thickness plot with Monte Carlo data points and a PK + build-up line using matplotlib.

## Dose vs thickness plot

This example runs Monte Carlo at a few thin concrete shields, fits a `BuildupFit`, and plots the extrapolated dose across a range of thicknesses.

```python exec="true" source="material-block" html="true"
from io import StringIO
import rad_point_kernel as rpk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Setup
concrete = rpk.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)

PARTICLES_PER_HOUR = 1e12 * 3600  # 1e12 photons/sec activity, scale to /hr
source = rpk.Source(particle="photon", energy=1e6)

mc_thicknesses = [5, 10, 15, 20]
all_thicknesses = list(range(5, 205, 5))

# Monte Carlo at thin shields
mc_geometries = [
    [rpk.Layer(thickness=t, material=concrete)]
    for t in mc_thicknesses
]
mc_results = rpk.compute_buildup(
    geometries=mc_geometries,
    source=source,
    quantities=["dose-AP-photon"],
)

# Fit
fit = rpk.BuildupFit(
    points=[{"thickness": t} for t in mc_thicknesses],
    results=mc_results,
)

# Compute dose at all thicknesses
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

# Monte Carlo reference points
mc_scaled = [r.scale(strength=PARTICLES_PER_HOUR) for r in mc_results]
mc_doses = [r.mc["dose-AP-photon"] for r in mc_scaled]
mc_errs = [r.mc_std_dev["dose-AP-photon"] for r in mc_scaled]

# Plot
fig, ax = plt.subplots(figsize=(10, 7))

ax.plot(all_thicknesses, doses, "b-", linewidth=2, label="PK with build-up")
ax.errorbar(
    mc_thicknesses, mc_doses, yerr=mc_errs,
    fmt="ko", markersize=7, capsize=4, zorder=5,
    label="Monte Carlo",
)

ax.set_xlabel("Concrete thickness (cm)", fontsize=13)
ax.set_ylabel("Photon dose rate (Sv/hr)", fontsize=13)
ax.set_title("Photon dose - 1 MeV, AP geometry", fontsize=13)
ax.set_yscale("log")
ax.legend(fontsize=11)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()

buf = StringIO()
fig.savefig(buf, format="svg", bbox_inches="tight")
plt.close(fig)
print(buf.getvalue())
```

## Reading the plot

- **Black points with error bars**: direct Monte Carlo simulation results (ground truth).
- **Blue line**: point-kernel dose multiplied by the fitted build-up factor.
- **Dashed line** (if shown): point-kernel uncollided dose without build-up correction. The gap between this and the solid line shows how much scattered radiation the build-up factor adds.

`BuildupFit` does not currently expose a predictive sigma (the field is NaN). If you need an uncertainty band, propagate the Monte Carlo `mc_std_dev / pk` envelope through your own fit, or watch the project roadmap for a future bootstrap-based predictive sigma.
