"""Multiple-layer shielding study: total dose (neutron + secondary photon).

Geometry: 10 m void + variable polyethylene + variable concrete.
MC (coupled neutron-photon) at thin shields, GP-extrapolated.
Results cached to avoid re-simulation on re-run.

Two plots:
  1. Total dose vs concrete thickness (one line per poly thickness)
  2. Total dose vs poly thickness (one line per concrete thickness)
"""

import json
import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import rad_point_kernel as pkc

matplotlib.use("Agg")

# --- Parameters ---
SOURCE_STRENGTH = 1e12
GEOMETRY = "AP"
VOID_THICKNESS = 1000
source = pkc.Source("neutron", 14.06e6)

N_DOSE = f"dose-{GEOMETRY}"
P_DOSE = f"dose-{GEOMETRY}-coupled-photon"

polyethylene = pkc.Material(
    composition={"H": 2, "C": 1}, density=0.94, fraction="atom",
)
concrete = pkc.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3, fraction="mass",
)

mc_poly = [0, 5, 10, 15, 20]
mc_conc = [10, 20, 30, 40]
all_poly = mc_poly + list(range(mc_poly[-1] + 5, 55, 5))
all_conc = mc_conc + list(range(mc_conc[-1] + 10, 210, 10))

RESULTS_DIR = Path(os.path.dirname(__file__), "..", "..", "results", "multiple_layer")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
CACHE_FILE = RESULTS_DIR / "mc_cache.json"


def make_layers(poly_t, conc_t):
    layers = [pkc.Layer(thickness=VOID_THICKNESS)]
    if poly_t > 0:
        layers.append(pkc.Layer(thickness=poly_t, material=polyethylene))
    layers.append(pkc.Layer(thickness=conc_t, material=concrete))
    return layers


# --- Step 1: MC with cache ---
cached = {}
if CACHE_FILE.exists():
    for entry in json.loads(CACHE_FILE.read_text()):
        key = (entry["poly"], entry["conc"])
        cached[key] = pkc.BuildupResult.from_dict(entry["result"])
    print(f"Loaded {len(cached)} cached points")

missing = [(p, c) for p in mc_poly for c in mc_conc if (p, c) not in cached]
if missing:
    print(f"Running coupled MC for {len(missing)} geometries...")
    mc_geometries = [make_layers(p, c) for p, c in missing]
    mc_list = pkc.compute_buildup(
        geometries=mc_geometries,
        source=source,
        quantities=[N_DOSE, P_DOSE],
        particles_per_batch=10_000,
        max_batches=100,
        trigger_rel_err=0.05,
    )
    for (p, c), r in zip(missing, mc_list):
        cached[(p, c)] = r

    cache_data = [
        {"poly": p, "conc": c, "result": cached[(p, c)].to_dict()}
        for p, c in sorted(cached)
    ]
    CACHE_FILE.write_text(json.dumps(cache_data, indent=2))
    print(f"  Saved {len(cached)} points to {CACHE_FILE}")
else:
    print("All points cached")

# --- Step 2: Build total buildup tables (one per poly thickness) ---
print("\nComputing total buildup factors...")

tables_by_poly = {}
for pt in mc_poly:
    br_list = []
    for ct in mc_conc:
        r = cached[(pt, ct)]
        mc_total = r.mc.get(N_DOSE, 0) + r.mc.get(P_DOSE, 0)
        mc_std = np.sqrt(r.mc_std_dev.get(N_DOSE, 0)**2 + r.mc_std_dev.get(P_DOSE, 0)**2)
        pk_neutron = r.pk.get(N_DOSE, 0)

        br = pkc.BuildupResult()
        br.mc["total"] = mc_total
        br.mc_std_dev["total"] = mc_std
        br.pk["total"] = pk_neutron
        br.buildup["total"] = mc_total / pk_neutron if pk_neutron > 0 else 1.0
        br_list.append(br)
        print(f"  poly={pt:>2d}, conc={ct:>3d}: B_total={br.buildup['total']:.3f}")

    tables_by_poly[pt] = pkc.BuildupTable(
        points=[{"conc": ct} for ct in mc_conc], results=br_list,
    )

# Also build tables by concrete thickness (for the second plot)
tables_by_conc = {}
for ct in mc_conc:
    br_list = []
    for pt in mc_poly:
        r = cached[(pt, ct)]
        mc_total = r.mc.get(N_DOSE, 0) + r.mc.get(P_DOSE, 0)
        mc_std = np.sqrt(r.mc_std_dev.get(N_DOSE, 0)**2 + r.mc_std_dev.get(P_DOSE, 0)**2)
        pk_neutron = r.pk.get(N_DOSE, 0)

        br = pkc.BuildupResult()
        br.mc["total"] = mc_total
        br.mc_std_dev["total"] = mc_std
        br.pk["total"] = pk_neutron
        br.buildup["total"] = mc_total / pk_neutron if pk_neutron > 0 else 1.0
        br_list.append(br)

    tables_by_conc[ct] = pkc.BuildupTable(
        points=[{"poly": pt} for pt in mc_poly], results=br_list,
    )

# --- Step 3: Extrapolate ---
print("\nExtrapolating...")

# Data for plot 1: dose vs concrete, one line per poly
data_vs_conc = {}
for pt in mc_poly:
    table = tables_by_poly[pt]
    doses, doses_lo, doses_hi = [], [], []
    for ct in all_conc:
        layers = make_layers(pt, ct)
        pk = pkc.calculate_dose(
            source_strength=SOURCE_STRENGTH, layers=layers,
            source=source, geometry=GEOMETRY,
        )
        bi = table.interpolate(conc=ct, warn=False)
        doses.append(pk.dose_rate * bi.value)
        doses_lo.append(pk.dose_rate * (bi.value - bi.sigma))
        doses_hi.append(pk.dose_rate * (bi.value + bi.sigma))
    data_vs_conc[pt] = {"doses": doses, "doses_lo": doses_lo, "doses_hi": doses_hi}

# Data for plot 2: dose vs poly, one line per concrete
data_vs_poly = {}
for ct in mc_conc:
    table = tables_by_conc[ct]
    doses, doses_lo, doses_hi = [], [], []
    for pt in all_poly:
        layers = make_layers(pt, ct)
        pk = pkc.calculate_dose(
            source_strength=SOURCE_STRENGTH, layers=layers,
            source=source, geometry=GEOMETRY,
        )
        bi = table.interpolate(poly=pt, warn=False)
        doses.append(pk.dose_rate * bi.value)
        doses_lo.append(pk.dose_rate * (bi.value - bi.sigma))
        doses_hi.append(pk.dose_rate * (bi.value + bi.sigma))
    data_vs_poly[ct] = {"doses": doses, "doses_lo": doses_lo, "doses_hi": doses_hi}

# --- Step 4: Plot ---
print("\nGenerating plots...")

colors_poly = {0: "#7f7f7f", 5: "#1f77b4", 10: "#ff7f0e", 15: "#2ca02c", 20: "#d62728"}
colors_conc = {10: "#1f77b4", 20: "#ff7f0e", 30: "#2ca02c", 40: "#d62728"}

# Plot 1: dose vs concrete thickness
fig, ax = plt.subplots(figsize=(10, 7))
for pt in mc_poly:
    c = colors_poly[pt]
    d = data_vs_conc[pt]
    ax.plot(all_conc, d["doses"], color=c, linewidth=2, label=f"{pt} cm poly")
    ax.fill_between(all_conc, d["doses_lo"], d["doses_hi"], color=c, alpha=0.12)

    mc_vals = [
        (cached[(pt, ct)].mc.get(N_DOSE, 0) + cached[(pt, ct)].mc.get(P_DOSE, 0)) * SOURCE_STRENGTH
        for ct in mc_conc
    ]
    mc_errs = [
        np.sqrt(cached[(pt, ct)].mc_std_dev.get(N_DOSE, 0)**2
                + cached[(pt, ct)].mc_std_dev.get(P_DOSE, 0)**2) * SOURCE_STRENGTH
        for ct in mc_conc
    ]
    ax.errorbar(mc_conc, mc_vals, yerr=mc_errs, fmt="o", color=c,
                markersize=6, capsize=3, zorder=5,
                label="Monte Carlo" if pt == mc_poly[0] else None)

ax.set_xlabel("Concrete thickness (cm)", fontsize=13)
ax.set_ylabel("Total dose rate (Sv/hr)", fontsize=13)
ax.set_title(
    f"Total dose vs concrete thickness\n"
    f"({VOID_THICKNESS/100:.0f} m void + poly + concrete, {ENERGY_EV/1e6:.2f} MeV, {GEOMETRY})",
    fontsize=12,
)
ax.set_yscale("log")
ax.legend(fontsize=10)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(RESULTS_DIR / "dose_vs_concrete.png", dpi=150, bbox_inches="tight")
print(f"  Saved dose_vs_concrete.png")
plt.close(fig)

# Plot 2: dose vs poly thickness
fig, ax = plt.subplots(figsize=(10, 7))
for ct in mc_conc:
    c = colors_conc[ct]
    d = data_vs_poly[ct]
    ax.plot(all_poly, d["doses"], color=c, linewidth=2, label=f"{ct} cm concrete")
    ax.fill_between(all_poly, d["doses_lo"], d["doses_hi"], color=c, alpha=0.12)

    mc_vals = [
        (cached[(pt, ct)].mc.get(N_DOSE, 0) + cached[(pt, ct)].mc.get(P_DOSE, 0)) * SOURCE_STRENGTH
        for pt in mc_poly
    ]
    mc_errs = [
        np.sqrt(cached[(pt, ct)].mc_std_dev.get(N_DOSE, 0)**2
                + cached[(pt, ct)].mc_std_dev.get(P_DOSE, 0)**2) * SOURCE_STRENGTH
        for pt in mc_poly
    ]
    ax.errorbar(mc_poly, mc_vals, yerr=mc_errs, fmt="o", color=c,
                markersize=6, capsize=3, zorder=5,
                label="Monte Carlo" if ct == mc_conc[0] else None)

ax.set_xlabel("Polyethylene thickness (cm)", fontsize=13)
ax.set_ylabel("Total dose rate (Sv/hr)", fontsize=13)
ax.set_title(
    f"Total dose vs polyethylene thickness\n"
    f"({VOID_THICKNESS/100:.0f} m void + poly + concrete, {ENERGY_EV/1e6:.2f} MeV, {GEOMETRY})",
    fontsize=12,
)
ax.set_yscale("log")
ax.legend(fontsize=10)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(RESULTS_DIR / "dose_vs_poly.png", dpi=150, bbox_inches="tight")
print(f"  Saved dose_vs_poly.png")
plt.close(fig)

print(f"\nAll plots saved to {RESULTS_DIR}")
