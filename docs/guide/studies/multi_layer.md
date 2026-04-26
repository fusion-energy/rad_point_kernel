# Multi-layer study

Scans two thickness parameters simultaneously (water then concrete) and shows total dose as a function of each, with the other held at a few representative values. MC is run on a sparse 2D grid; 1D GP tables extrapolate along each axis.

## Setup

```python exec="true" source="material-block" session="multi" result="text"
import json
from pathlib import Path
import numpy as np
import rad_point_kernel as rpk

SOURCE_STRENGTH = 1e12
GEOMETRY = "AP"
VOID_THICKNESS = 1000
source = rpk.Source(particle="neutron", energy=14.06e6)

N_DOSE = f"dose-{GEOMETRY}"
P_DOSE = f"dose-{GEOMETRY}-coupled-photon"

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
        quantities=[N_DOSE, P_DOSE],
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

## Build per-axis buildup tables

```python exec="true" source="material-block" session="multi" result="text"
def total_buildup_result(r):
    mc_total = r.mc.get(N_DOSE, 0) + r.mc.get(P_DOSE, 0)
    mc_std = np.sqrt(
        r.mc_std_dev.get(N_DOSE, 0) ** 2 + r.mc_std_dev.get(P_DOSE, 0) ** 2
    )
    pk_neutron = r.pk.get(N_DOSE, 0)
    br = rpk.BuildupResult()
    br.mc["total"] = mc_total
    br.mc_std_dev["total"] = mc_std
    br.pk["total"] = pk_neutron
    br.buildup["total"] = mc_total / pk_neutron if pk_neutron > 0 else 1.0
    return br

tables_by_water = {
    wt: rpk.BuildupTable(
        points=[{"conc": ct} for ct in mc_conc],
        results=[total_buildup_result(cached[(wt, ct)]) for ct in mc_conc],
    )
    for wt in mc_water
}
tables_by_conc = {
    ct: rpk.BuildupTable(
        points=[{"water": wt} for wt in mc_water],
        results=[total_buildup_result(cached[(wt, ct)]) for wt in mc_water],
    )
    for ct in mc_conc
}
print(f"Built {len(tables_by_water)} tables by water, {len(tables_by_conc)} tables by concrete")
```

## Extrapolate and plot

```python exec="true" source="material-block" session="multi" html="true"
import io
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")


def extrapolate(table, ts, key, arg_name):
    doses, lo, hi = [], [], []
    for t in ts:
        if arg_name == "conc":
            layers = make_layers(key, t)
            bi = table.interpolate(conc=t, warn=False)
        else:
            layers = make_layers(t, key)
            bi = table.interpolate(water=t, warn=False)
        pk = rpk.calculate_dose(
            source_strength=SOURCE_STRENGTH,
            layers=layers,
            source=source,
            geometry=GEOMETRY,
        )
        doses.append(pk.dose_rate * bi.value)
        lo.append(pk.dose_rate * (bi.value - bi.sigma))
        hi.append(pk.dose_rate * (bi.value + bi.sigma))
    return doses, lo, hi


data_vs_conc = {
    wt: extrapolate(tables_by_water[wt], all_conc, wt, "conc") for wt in mc_water
}
data_vs_water = {
    ct: extrapolate(tables_by_conc[ct], all_water, ct, "water") for ct in mc_conc
}

colors_water = {0: "#7f7f7f", 10: "#1f77b4", 20: "#ff7f0e",
                30: "#2ca02c", 40: "#d62728", 50: "#9467bd"}
colors_conc = {10: "#1f77b4", 50: "#ff7f0e", 100: "#2ca02c",
               200: "#d62728", 400: "#9467bd"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

for wt in mc_water:
    c = colors_water[wt]
    d, lo, hi = data_vs_conc[wt]
    ax1.plot(all_conc, d, color=c, linewidth=2, label=f"{wt} cm water")
    ax1.fill_between(all_conc, lo, hi, color=c, alpha=0.12)
    mc_vals = [
        (cached[(wt, ct)].mc.get(N_DOSE, 0) + cached[(wt, ct)].mc.get(P_DOSE, 0))
        * SOURCE_STRENGTH for ct in mc_conc
    ]
    ax1.scatter(mc_conc, mc_vals, color=c, s=30, zorder=5)

ax1.set_xlabel("Concrete thickness (cm)")
ax1.set_ylabel("Total dose rate (Sv/hr)")
ax1.set_yscale("log")
ax1.set_title("Total dose vs concrete (per water thickness)")
ax1.legend(fontsize=9)
ax1.grid(True, which="both", alpha=0.3)

for ct in mc_conc:
    c = colors_conc[ct]
    d, lo, hi = data_vs_water[ct]
    ax2.plot(all_water, d, color=c, linewidth=2, label=f"{ct} cm concrete")
    ax2.fill_between(all_water, lo, hi, color=c, alpha=0.12)
    mc_vals = [
        (cached[(wt, ct)].mc.get(N_DOSE, 0) + cached[(wt, ct)].mc.get(P_DOSE, 0))
        * SOURCE_STRENGTH for wt in mc_water
    ]
    ax2.scatter(mc_water, mc_vals, color=c, s=30, zorder=5)

ax2.set_xlabel("Water thickness (cm)")
ax2.set_ylabel("Total dose rate (Sv/hr)")
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
