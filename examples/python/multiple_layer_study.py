"""Multiple-layer shielding study: total dose (neutron + secondary photon).

Geometry: 10 m void + variable water + variable concrete.
MC (coupled neutron-photon) at thin shields, GP-extrapolated.
Results cached to avoid re-simulation on re-run.

Two plots:
  1. Total dose vs concrete thickness (one line per water thickness)
  2. Total dose vs water thickness (one line per concrete thickness)
"""

import json
import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import rad_point_kernel as rpk

matplotlib.use("Agg")

# --- Parameters ---
SOURCE_STRENGTH = 1e12
GEOMETRY = "AP"
VOID_THICKNESS = 1000
source = rpk.Source("neutron", 14.06e6)

N_DOSE = f"dose-{GEOMETRY}"
P_DOSE = f"dose-{GEOMETRY}-coupled-photon"

water = rpk.Material(composition={"H2O": 1.0}, density=1.0)
concrete = rpk.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3,
    fraction="mass",
)

mc_water = [0, 10, 20, 30, 40, 50]
mc_conc = [10, 50, 100, 200, 400]
all_water = list(range(0, 55, 5))
all_conc = list(range(10, 410, 10))

RESULTS_DIR = Path(os.path.dirname(__file__), "..", "..", "results", "multiple_layer")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
CACHE_FILE = RESULTS_DIR / "mc_cache.json"


def make_layers(water_t, conc_t):
    layers = [rpk.Layer(thickness=VOID_THICKNESS)]
    if water_t > 0:
        layers.append(rpk.Layer(thickness=water_t, material=water))
    layers.append(rpk.Layer(thickness=conc_t, material=concrete))
    return layers


# --- Step 1: MC with cache ---
cached = {}
if CACHE_FILE.exists():
    for entry in json.loads(CACHE_FILE.read_text()):
        key = (entry["water"], entry["conc"])
        cached[key] = rpk.BuildupResult.from_dict(entry["result"])
    print(f"Loaded {len(cached)} cached points")

missing = [(w, c) for w in mc_water for c in mc_conc if (w, c) not in cached]
if missing:
    print(f"Running coupled MC for {len(missing)} geometries...")
    mc_geometries = [make_layers(w, c) for w, c in missing]
    mc_list = rpk.compute_buildup(
        geometries=mc_geometries,
        source=source,
        quantities=[N_DOSE, P_DOSE],
        particles_per_batch=10_000,
        max_batches=100,
        trigger_rel_err=0.05,
    )
    for (w, c), r in zip(missing, mc_list):
        cached[(w, c)] = r

    cache_data = [
        {"water": w, "conc": c, "result": cached[(w, c)].to_dict()}
        for w, c in sorted(cached)
    ]
    CACHE_FILE.write_text(json.dumps(cache_data, indent=2))
    print(f"  Saved {len(cached)} points to {CACHE_FILE}")
else:
    print("All points cached")

# --- Step 2: Build total buildup tables (one per water thickness) ---
print("\nComputing total buildup factors...")

tables_by_water = {}
for wt in mc_water:
    br_list = []
    for ct in mc_conc:
        r = cached[(wt, ct)]
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
        br_list.append(br)
        print(f"  water={wt:>2d}, conc={ct:>3d}: B_total={br.buildup['total']:.3f}")

    tables_by_water[wt] = rpk.BuildupTable(
        points=[{"conc": ct} for ct in mc_conc],
        results=br_list,
    )

# Also build tables by concrete thickness (for the second plot)
tables_by_conc = {}
for ct in mc_conc:
    br_list = []
    for wt in mc_water:
        r = cached[(wt, ct)]
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
        br_list.append(br)

    tables_by_conc[ct] = rpk.BuildupTable(
        points=[{"water": wt} for wt in mc_water],
        results=br_list,
    )

# --- Step 3: Extrapolate ---
print("\nExtrapolating...")

# Data for plot 1: dose vs concrete, one line per water
data_vs_conc = {}
for wt in mc_water:
    table = tables_by_water[wt]
    doses, doses_lo, doses_hi = [], [], []
    for ct in all_conc:
        layers = make_layers(wt, ct)
        pk = rpk.calculate_dose(
            source_strength=SOURCE_STRENGTH,
            layers=layers,
            source=source,
            geometry=GEOMETRY,
        )
        bi = table.interpolate(conc=ct, warn=False)
        doses.append(pk.dose_rate * bi.value)
        doses_lo.append(pk.dose_rate * (bi.value - bi.sigma))
        doses_hi.append(pk.dose_rate * (bi.value + bi.sigma))
    data_vs_conc[wt] = {"doses": doses, "doses_lo": doses_lo, "doses_hi": doses_hi}

# Data for plot 2: dose vs water, one line per concrete
data_vs_water = {}
for ct in mc_conc:
    table = tables_by_conc[ct]
    doses, doses_lo, doses_hi = [], [], []
    for wt in all_water:
        layers = make_layers(wt, ct)
        pk = rpk.calculate_dose(
            source_strength=SOURCE_STRENGTH,
            layers=layers,
            source=source,
            geometry=GEOMETRY,
        )
        bi = table.interpolate(water=wt, warn=False)
        doses.append(pk.dose_rate * bi.value)
        doses_lo.append(pk.dose_rate * (bi.value - bi.sigma))
        doses_hi.append(pk.dose_rate * (bi.value + bi.sigma))
    data_vs_water[ct] = {"doses": doses, "doses_lo": doses_lo, "doses_hi": doses_hi}

# --- Step 4: Plot ---
print("\nGenerating plots...")

colors_water = {
    0: "#7f7f7f", 10: "#1f77b4", 20: "#ff7f0e",
    30: "#2ca02c", 40: "#d62728", 50: "#9467bd",
}
colors_conc = {
    10: "#1f77b4", 50: "#ff7f0e", 100: "#2ca02c",
    200: "#d62728", 400: "#9467bd",
}

# Plot 1: dose vs concrete thickness
fig, ax = plt.subplots(figsize=(10, 7))
for wt in mc_water:
    c = colors_water[wt]
    d = data_vs_conc[wt]
    ax.plot(all_conc, d["doses"], color=c, linewidth=2, label=f"{wt} cm water")
    ax.fill_between(all_conc, d["doses_lo"], d["doses_hi"], color=c, alpha=0.12)

    mc_vals = [
        (cached[(wt, ct)].mc.get(N_DOSE, 0) + cached[(wt, ct)].mc.get(P_DOSE, 0))
        * SOURCE_STRENGTH
        for ct in mc_conc
    ]
    mc_errs = [
        np.sqrt(
            cached[(wt, ct)].mc_std_dev.get(N_DOSE, 0) ** 2
            + cached[(wt, ct)].mc_std_dev.get(P_DOSE, 0) ** 2
        )
        * SOURCE_STRENGTH
        for ct in mc_conc
    ]
    ax.errorbar(
        mc_conc,
        mc_vals,
        yerr=mc_errs,
        fmt="o",
        color=c,
        markersize=6,
        capsize=3,
        zorder=5,
        label="Monte Carlo" if wt == mc_water[0] else None,
    )

ax.set_xlabel("Concrete thickness (cm)", fontsize=13)
ax.set_ylabel("Total dose rate (Sv/hr)", fontsize=13)
ax.set_title(
    f"Total dose vs concrete thickness\n"
    f"({VOID_THICKNESS/100:.0f} m void + water + concrete, {source.energy/1e6:.2f} MeV, {GEOMETRY})",
    fontsize=12,
)
ax.set_yscale("log")
ax.legend(fontsize=10)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(RESULTS_DIR / "dose_vs_concrete.png", dpi=150, bbox_inches="tight")
print(f"  Saved dose_vs_concrete.png")
plt.close(fig)

# Plot 2: dose vs water thickness
fig, ax = plt.subplots(figsize=(10, 7))
for ct in mc_conc:
    c = colors_conc[ct]
    d = data_vs_water[ct]
    ax.plot(all_water, d["doses"], color=c, linewidth=2, label=f"{ct} cm concrete")
    ax.fill_between(all_water, d["doses_lo"], d["doses_hi"], color=c, alpha=0.12)

    mc_vals = [
        (cached[(wt, ct)].mc.get(N_DOSE, 0) + cached[(wt, ct)].mc.get(P_DOSE, 0))
        * SOURCE_STRENGTH
        for wt in mc_water
    ]
    mc_errs = [
        np.sqrt(
            cached[(wt, ct)].mc_std_dev.get(N_DOSE, 0) ** 2
            + cached[(wt, ct)].mc_std_dev.get(P_DOSE, 0) ** 2
        )
        * SOURCE_STRENGTH
        for wt in mc_water
    ]
    ax.errorbar(
        mc_water,
        mc_vals,
        yerr=mc_errs,
        fmt="o",
        color=c,
        markersize=6,
        capsize=3,
        zorder=5,
        label="Monte Carlo" if ct == mc_conc[0] else None,
    )

ax.set_xlabel("Water thickness (cm)", fontsize=13)
ax.set_ylabel("Total dose rate (Sv/hr)", fontsize=13)
ax.set_title(
    f"Total dose vs water thickness\n"
    f"({VOID_THICKNESS/100:.0f} m void + water + concrete, {source.energy/1e6:.2f} MeV, {GEOMETRY})",
    fontsize=12,
)
ax.set_yscale("log")
ax.legend(fontsize=10)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(RESULTS_DIR / "dose_vs_water.png", dpi=150, bbox_inches="tight")
print(f"  Saved dose_vs_water.png")
plt.close(fig)

print(f"\nAll plots saved to {RESULTS_DIR}")
