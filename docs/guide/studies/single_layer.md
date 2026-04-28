# Single-layer study

Scans shield thickness from thin to thick for two concrete mixes (Portland + 3 % steel rebar, and Magnetite + 3 % steel rebar), computing total coupled dose (neutron + secondary photon) at each point. Monte Carlo is run only at a sparse set of thin shields and the resulting build-up factor is GP-extrapolated out to 400 cm.

## Setup

```python exec="true" source="material-block" session="single" result="text"
import json
from pathlib import Path
import numpy as np
import rad_point_kernel as rpk

PARTICLES_PER_SHOT = 1e16
GEOMETRY = "AP"
VOID_THICKNESS = 1000  # source-to-shield air gap, cm
source = rpk.Source(particle="neutron", energy=14.1e6)

N_DOSE = f"dose-{GEOMETRY}"
P_DOSE = f"dose-{GEOMETRY}-coupled-photon"
TOTAL_DOSE = f"dose-{GEOMETRY}-total"  # auto-synthesized by compute_buildup

mc_thicknesses = [10, 20, 30, 40, 50, 60]
all_thicknesses = mc_thicknesses + list(range(70, 410, 10))

portland = rpk.Material(
    composition={
        "H": 0.168753, "C": 0.001416, "O": 0.562525, "Na": 0.011838,
        "Mg": 0.0014, "Al": 0.021354, "Si": 0.204119, "K": 0.005656,
        "Ca": 0.018674, "Fe": 0.004264,
    },
    density=2.3,
    fraction="atom",
)
magnetite = rpk.Material(
    composition={
        "H": 0.082377, "O": 0.551, "Mg": 0.010248, "Al": 0.023218,
        "Si": 0.024456, "S": 0.001177, "Ca": 0.047269, "Ti": 0.030274,
        "V": 0.00163, "Cr": 0.000871, "Mn": 0.000962, "Fe": 0.226518,
    },
    density=3.53,
    fraction="atom",
)
steel = rpk.Material(
    composition={"C": 0.003747, "Si": 0.003163, "P": 0.001127,
                 "S": 0.000175, "Mn": 0.000154, "Fe": 0.991634},
    density=7.7,
    fraction="atom",
)

materials = {
    "Portland + 3% steel": rpk.Material.volume_mix(portland, 0.97, steel, 0.03),
    "Magnetite + 3% steel": rpk.Material.volume_mix(magnetite, 0.97, steel, 0.03),
}

CACHE_DIR = Path("docs/assets/studies")
CACHE_DIR.mkdir(parents=True, exist_ok=True)
```

## Monte Carlo on thin shields

One cache file per material; delete the file to re-simulate.

```python exec="true" source="material-block" session="single" result="text"
def cache_path(name):
    safe = name.replace(" ", "_").replace("%", "pct").replace("+", "plus").lower()
    return CACHE_DIR / f"single_layer_{safe}.json"

tables = {}
mc_cached = {}

for name, mat in materials.items():
    cache_file = cache_path(name)

    cached = {}
    if cache_file.exists():
        for entry in json.loads(cache_file.read_text()):
            cached[entry["thickness"]] = rpk.BuildupResult.from_dict(entry["result"])

    missing = [t for t in mc_thicknesses if t not in cached]
    if missing:
        geometries = [
            [rpk.Layer(thickness=VOID_THICKNESS), rpk.Layer(thickness=t, material=mat)]
            for t in missing
        ]
        new_results = rpk.compute_buildup(
            geometries=geometries,
            source=source,
            quantities=[N_DOSE, P_DOSE],
            particles_per_batch=10_000,
            max_batches=100,
            trigger_rel_err=0.05,
        )
        for t, r in zip(missing, new_results):
            cached[t] = r
        cache_data = [
            {"thickness": t, "result": cached[t].to_dict()} for t in sorted(cached)
        ]
        cache_file.write_text(json.dumps(cache_data, indent=2))

    tables[name] = rpk.BuildupTable(
        points=[{"thickness": t} for t in mc_thicknesses],
        results=[cached[t] for t in mc_thicknesses],
    )
    mc_cached[name] = cached

print(f"Loaded MC for {len(materials)} materials")
```

## Extrapolate and plot

```python exec="true" source="material-block" session="single" html="true"
import io
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")

data = {}
for name, mat in materials.items():
    table = tables[name]
    doses, doses_lo, doses_hi = [], [], []
    for t in all_thicknesses:
        layers = [
            rpk.Layer(thickness=VOID_THICKNESS),
            rpk.Layer(thickness=t, material=mat),
        ]
        pk = rpk.calculate_dose(
            layers=layers,
            source=source,
            geometry=GEOMETRY,
        ).scale(strength=PARTICLES_PER_SHOT)
        bi = table.interpolate(thickness=t, quantity=TOTAL_DOSE, warn=False)
        doses.append(pk.dose * bi.value)
        doses_lo.append(pk.dose * (bi.value - bi.sigma))
        doses_hi.append(pk.dose * (bi.value + bi.sigma))
    data[name] = {
        "doses": np.array(doses),
        "doses_lo": np.array(doses_lo),
        "doses_hi": np.array(doses_hi),
    }

colors = {"Portland + 3% steel": "#1f77b4", "Magnetite + 3% steel": "#ff7f0e"}
fig, ax = plt.subplots(figsize=(10, 6))

for name in materials:
    c = colors[name]
    d = data[name]
    valid = d["doses"] > 0
    ax.plot(np.array(all_thicknesses)[valid], d["doses"][valid],
            color=c, linewidth=2, label=name)
    v2 = valid & (d["doses_lo"] > 0) & (d["doses_hi"] > 0)
    ax.fill_between(np.array(all_thicknesses)[v2],
                    d["doses_lo"][v2], d["doses_hi"][v2],
                    color=c, alpha=0.15)

    mc_vals = [mc_cached[name][t].mc[TOTAL_DOSE] * PARTICLES_PER_SHOT for t in mc_thicknesses]
    mc_errs = [mc_cached[name][t].mc_std_dev[TOTAL_DOSE] * PARTICLES_PER_SHOT for t in mc_thicknesses]
    ax.errorbar(mc_thicknesses, mc_vals, yerr=mc_errs,
                fmt="o", color=c, markersize=6, capsize=3, zorder=5)

ax.set_xlabel("Shield thickness (cm)")
ax.set_ylabel("Total dose (Sv/shot)")
ax.set_title(
    f"Total dose (neutron + secondary gamma) - {source.energy/1e6} MeV pulsed DT, {GEOMETRY}"
)
ax.set_yscale("log")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()

buf = io.StringIO()
fig.savefig(buf, format="svg", bbox_inches="tight")
plt.close(fig)
print(buf.getvalue())
```

Monte Carlo points (dots with error bars) are the six thin-shield simulations. Solid lines are PK dose multiplied by the GP-extrapolated build-up factor; shaded bands are the GP +/- one-sigma uncertainty.
