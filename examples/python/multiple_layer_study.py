"""Generate docs/assets/buildup_example.png.

Sweeps a 2D water x concrete grid behind a 14.06 MeV neutron source, runs
coupled neutron-photon Monte Carlo at each grid point, and produces a single
plot of total-dose buildup factor vs concrete thickness with water thickness
as the family axis. The figure illustrates the GP extrapolation story in
docs/theory/buildup_factors.md.

Run with OPENMC_CROSS_SECTIONS pointing at an ENDF/B-VIII.1 cross_sections.xml.
MC results are cached at docs/assets/studies/multi_layer_cache.json - delete
that file to force a re-run.
"""
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import rad_point_kernel as rpk

matplotlib.use("Agg")

GEOMETRY = "AP"
VOID_THICKNESS = 1000  # source-to-shield air gap, cm
TOTAL_DOSE = f"dose-{GEOMETRY}-total"

source = rpk.Source(particle="neutron", energy=14.06e6)

water = rpk.Material(composition={"H2O": 1.0}, density=1.0)
concrete = rpk.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3,
    fraction="mass",
)

mc_water = [0, 10, 20, 30, 40, 50]
mc_conc = [10, 50, 100, 200]
all_conc = list(range(5, 251, 5))

CACHE = Path(__file__).resolve().parent.parent.parent / "docs/assets/studies/multi_layer_cache.json"
CACHE.parent.mkdir(parents=True, exist_ok=True)
OUTPUT = Path(__file__).resolve().parent.parent.parent / "docs/assets/buildup_example.png"


def make_layers(water_t, conc_t):
    layers = [rpk.Layer(thickness=VOID_THICKNESS)]
    if water_t > 0:
        layers.append(rpk.Layer(thickness=water_t, material=water))
    layers.append(rpk.Layer(thickness=conc_t, material=concrete))
    return layers


cached = {}
if CACHE.exists():
    for entry in json.loads(CACHE.read_text()):
        cached[(entry["water"], entry["conc"])] = rpk.BuildupResult.from_dict(entry["result"])

missing = [(w, c) for w in mc_water for c in mc_conc if (w, c) not in cached]
if missing:
    print(f"Running coupled MC for {len(missing)} geometries...")
    new_results = rpk.compute_buildup(
        geometries=[make_layers(w, c) for w, c in missing],
        source=source,
        quantities=[TOTAL_DOSE],
        particles_per_batch=5_000,
        max_batches=200,
        trigger_rel_err=0.05,
    )
    for (w, c), r in zip(missing, new_results):
        cached[(w, c)] = r
    CACHE.write_text(json.dumps(
        [{"water": w, "conc": c, "result": cached[(w, c)].to_dict()}
         for w, c in sorted(cached)],
        indent=2,
    ))

# 2D table over the water x concrete grid
table = rpk.BuildupTable(
    points=[{"water": w, "conc": c} for w in mc_water for c in mc_conc],
    results=[cached[(w, c)] for w in mc_water for c in mc_conc],
)


def b_curve(water_t):
    values, lo, hi = [], [], []
    for c in all_conc:
        ir = table.interpolate(water=water_t, conc=c, quantity=TOTAL_DOSE, warn=False)
        values.append(ir.value)
        lo.append(ir.value - ir.sigma)
        hi.append(ir.value + ir.sigma)
    return np.array(values), np.array(lo), np.array(hi)


fig, ax = plt.subplots(figsize=(10, 6))
colors = {0: "#7f7f7f", 10: "#1f77b4", 20: "#ff7f0e",
          30: "#2ca02c", 40: "#d62728", 50: "#9467bd"}

for w in mc_water:
    c_color = colors[w]
    b, b_lo, b_hi = b_curve(w)
    ax.plot(all_conc, b, color=c_color, linewidth=2, label=f"{w} cm water")
    ax.fill_between(all_conc, b_lo, b_hi, color=c_color, alpha=0.15)

    mc_b = [cached[(w, c)].buildup[TOTAL_DOSE] for c in mc_conc]
    mc_err = [cached[(w, c)].mc_std_dev[TOTAL_DOSE] / cached[(w, c)].pk[TOTAL_DOSE]
              for c in mc_conc]
    ax.errorbar(mc_conc, mc_b, yerr=mc_err, fmt="o",
                color=c_color, markersize=6, capsize=3, zorder=5)

ax.set_xlabel("Concrete thickness (cm)")
ax.set_ylabel("Build-up factor B")
ax.set_yscale("log")
ax.set_title("Total-dose build-up factor (water + concrete, 14.06 MeV neutron, AP)")
ax.legend(fontsize=10)
ax.grid(True, which="both", alpha=0.3)

fig.tight_layout()
fig.savefig(OUTPUT, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Wrote {OUTPUT}")
