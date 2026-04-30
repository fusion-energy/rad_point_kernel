# Multi-layer study

Scans two thickness parameters simultaneously (water then concrete) and shows total dose as a function of each, with the other held at a few representative values. MC is run on a sparse 2D grid; a single 2D `BuildupFit` (thin-plate-spline RBF) interpolates across both axes.

## Setup

```python exec="true" source="material-block" session="multi" result="text"
import json
from pathlib import Path
import rad_point_kernel as rpk

PARTICLES_PER_SHOT = 1e16
GEOMETRY = "AP"
VOID_THICKNESS = 1000  # source-to-shield air gap, cm
source = rpk.Source(particle="neutron", energy=14.06e6)

TOTAL_DOSE = f"dose-{GEOMETRY}-total"

water = rpk.Material(composition={"H2O": 1.0}, density=1.0)
concrete = rpk.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3,
    fraction="mass",
)

mc_water = [0, 5, 10, 15, 20, 25]
mc_conc = [5, 25, 50, 100, 200]
all_water = list(range(0, 30, 5))
all_conc = list(range(5, 205, 5))

CACHE = Path("docs/assets/studies/multi_layer_cache.json")
CACHE.parent.mkdir(parents=True, exist_ok=True)


def make_layers(water_t, conc_t):
    layers = [rpk.Layer(thickness=VOID_THICKNESS)]
    if water_t > 0:
        layers.append(rpk.Layer(thickness=water_t, material=water))
    layers.append(rpk.Layer(thickness=conc_t, material=concrete))
    return layers
```

## Monte Carlo on the 2D grid

```python exec="true" source="material-block" session="multi" result="text"
cached = {}
if CACHE.exists():
    for entry in json.loads(CACHE.read_text()):
        cached[(entry["water"], entry["conc"])] = rpk.BuildupResult.from_dict(entry["result"])

missing = [(w, c) for w in mc_water for c in mc_conc if (w, c) not in cached]
if missing:
    mc_geometries = [make_layers(w, c) for w, c in missing]
    new_results = rpk.compute_buildup(
        geometries=mc_geometries,
        source=source,
        quantities=[TOTAL_DOSE],
        particles_per_batch=10_000,
        max_batches=100,
        trigger_rel_err=0.05,
    )
    for (w, c), r in zip(missing, new_results):
        cached[(w, c)] = r
    cache_data = [
        {"water": w, "conc": c, "result": cached[(w, c)].to_dict()}
        for w, c in sorted(cached)
    ]
    CACHE.write_text(json.dumps(cache_data, indent=2))

print(f"Grid: {len(mc_water)} water x {len(mc_conc)} concrete = {len(cached)} MC points")
```

## Build a 2D buildup table

`compute_buildup` adds `dose-{GEOMETRY}-total` automatically because both the neutron and coupled-photon dose were requested. A single 2D `BuildupFit` keyed on water and concrete thicknesses is enough to interpolate anywhere on the grid:

```python exec="true" source="material-block" session="multi" result="text"
fit = rpk.BuildupFit(
    points=[{"water": wt, "conc": ct} for wt in mc_water for ct in mc_conc],
    results=[cached[(wt, ct)] for wt in mc_water for ct in mc_conc],
)
print(f"Built 2D fit over {fit.axis_ranges}")
```

## Extrapolate and plot

```python exec="true" source="material-block" session="multi" html="true"
import io
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")


def dose_at(water_t, conc_t):
    layers = make_layers(water_t, conc_t)
    bi = fit.interpolate(water=water_t, conc=conc_t, quantity=TOTAL_DOSE, warn=False)
    pk = rpk.calculate_dose(
        layers=layers, source=source, geometry=GEOMETRY,
    ).scale(strength=PARTICLES_PER_SHOT)
    return pk.dose * bi.value


def curve(water_ts, conc_ts):
    return [dose_at(w, c) for w, c in zip(water_ts, conc_ts)]


data_vs_conc = {wt: curve([wt] * len(all_conc), all_conc) for wt in mc_water}
data_vs_water = {ct: curve(all_water, [ct] * len(all_water)) for ct in mc_conc}

colors_water = {0: "#7f7f7f", 5: "#1f77b4", 10: "#ff7f0e",
                15: "#2ca02c", 20: "#d62728", 25: "#9467bd"}
colors_conc = {5: "#1f77b4", 25: "#ff7f0e", 50: "#2ca02c",
               100: "#d62728", 200: "#9467bd"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

for wt in mc_water:
    c = colors_water[wt]
    d = data_vs_conc[wt]
    ax1.plot(all_conc, d, color=c, linewidth=2, label=f"{wt} cm water")
    mc_vals = [cached[(wt, ct)].mc[TOTAL_DOSE] * PARTICLES_PER_SHOT for ct in mc_conc]
    ax1.scatter(mc_conc, mc_vals, color=c, s=30, zorder=5)

ax1.set_xlabel("Concrete thickness (cm)")
ax1.set_ylabel("Total dose (Sv/shot)")
ax1.set_yscale("log")
ax1.set_title("Total dose vs concrete (per water thickness)")
ax1.legend(fontsize=9)
ax1.grid(True, which="both", alpha=0.3)

for ct in mc_conc:
    c = colors_conc[ct]
    d = data_vs_water[ct]
    ax2.plot(all_water, d, color=c, linewidth=2, label=f"{ct} cm concrete")
    mc_vals = [cached[(wt, ct)].mc[TOTAL_DOSE] * PARTICLES_PER_SHOT for wt in mc_water]
    ax2.scatter(mc_water, mc_vals, color=c, s=30, zorder=5)

ax2.set_xlabel("Water thickness (cm)")
ax2.set_ylabel("Total dose (Sv/shot)")
ax2.set_yscale("log")
ax2.set_title("Total dose vs water (per concrete thickness)")
ax2.legend(fontsize=9)
ax2.grid(True, which="both", alpha=0.3)

fig.tight_layout()
buf = io.StringIO()
fig.savefig(buf, format="svg", bbox_inches="tight")
plt.close(fig)
print(buf.getvalue())
```
