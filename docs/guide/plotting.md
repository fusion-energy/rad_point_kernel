# Plotting results

This page shows how to create a dose-vs-thickness plot with Monte Carlo data points, a PK + build-up line, and a GP uncertainty band using matplotlib.

## Dose vs thickness plot

This example runs Monte Carlo at a few thin concrete shields, builds a GP table, and plots the extrapolated dose across a range of thicknesses.

```python
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

SOURCE = 1e12
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
    quantities=["dose-AP"],
)

# Build GP table
table = rpk.BuildupTable(
    points=[{"thickness": t} for t in mc_thicknesses],
    results=mc_results,
)

# Compute dose at all thicknesses
doses = []
doses_lo = []
doses_hi = []

for t in all_thicknesses:
    layers = [rpk.Layer(thickness=t, material=concrete)]
    bi = table.interpolate(thickness=t)
    pk = rpk.calculate_dose(
        source_strength=SOURCE,
        layers=layers,
        source=source,
        geometry="AP",
    )
    doses.append(pk.dose_rate * bi.value)
    doses_lo.append(pk.dose_rate * (bi.value - bi.sigma))
    doses_hi.append(pk.dose_rate * (bi.value + bi.sigma))

# Monte Carlo reference points
mc_doses = [r.mc["dose-AP"] * SOURCE for r in mc_results]
mc_errs = [r.mc_std_dev["dose-AP"] * SOURCE for r in mc_results]

# Plot
fig, ax = plt.subplots(figsize=(10, 7))

ax.plot(all_thicknesses, doses, "b-", linewidth=2, label="PK with build-up")
ax.fill_between(
    all_thicknesses, doses_lo, doses_hi,
    color="blue", alpha=0.15, label="GP uncertainty",
)
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

fig.savefig("dose_vs_thickness.png", dpi=150, bbox_inches="tight")
print("Plot saved to dose_vs_thickness.png")
```

The result looks like this:

![Dose vs thickness](../assets/dose_buildup_comparison.png)

## Reading the plot

- **Black points with error bars**: direct Monte Carlo simulation results (ground truth).
- **Blue line**: point-kernel dose multiplied by the GP-predicted build-up factor.
- **Blue band**: 1-sigma GP uncertainty. The band is tight near Monte Carlo data points and widens as you move further away, especially in the extrapolation region beyond the thickest Monte Carlo shield.
- **Dashed line** (if shown): point-kernel uncollided dose without build-up correction. The gap between this and the solid line shows how much scattered radiation the build-up factor adds.
