"""Single-layer shielding study: total dose (neutron + secondary photon).

Geometry: 10 m void + single material layer (0-400 cm).
MC (coupled neutron-photon) at thin shields, GP-extrapolated to 400 cm.
Results cached to avoid re-simulation on re-run.

Materials:
  1. Portland concrete + 3% vol steel rebar
  2. Magnetite concrete + 3% vol steel rebar
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
source = rpk.Source("neutron", 14.1e6)

N_DOSE = f"dose-{GEOMETRY}"
P_DOSE = f"dose-{GEOMETRY}-coupled-photon"

mc_thicknesses = [10, 20, 30, 40, 50, 60]
all_thicknesses = mc_thicknesses + list(range(mc_thicknesses[-1] + 10, 410, 10))

RESULTS_DIR = Path(os.path.dirname(__file__), "..", "..", "results", "single_layer")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# --- Materials (PNNL compositions, atom fractions) ---
portland = rpk.Material(
    composition={
        "H": 0.168753,
        "C": 0.001416,
        "O": 0.562525,
        "Na": 0.011838,
        "Mg": 0.0014,
        "Al": 0.021354,
        "Si": 0.204119,
        "K": 0.005656,
        "Ca": 0.018674,
        "Fe": 0.004264,
    },
    density=2.3,
    fraction="atom",
)
magnetite = rpk.Material(
    composition={
        "H": 0.082377,
        "O": 0.551,
        "Mg": 0.010248,
        "Al": 0.023218,
        "Si": 0.024456,
        "S": 0.001177,
        "Ca": 0.047269,
        "Ti": 0.030274,
        "V": 0.00163,
        "Cr": 0.000871,
        "Mn": 0.000962,
        "Fe": 0.226518,
    },
    density=3.53,
    fraction="atom",
)
steel = rpk.Material(
    composition={
        "C": 0.003747,
        "Si": 0.003163,
        "P": 0.001127,
        "S": 0.000175,
        "Mn": 0.000154,
        "Fe": 0.991634,
    },
    density=7.7,
    fraction="atom",
)

materials = {
    "Portland + 3% steel": rpk.Material.volume_mix(portland, 0.97, steel, 0.03),
    "Magnetite + 3% steel": rpk.Material.volume_mix(magnetite, 0.97, steel, 0.03),
}

# --- Step 1: MC with cache ---
tables = {}
mc_cached = {}

for name, mat in materials.items():
    cache_file = RESULTS_DIR / f"mc_{name.replace(' ', '_').lower()}.json"

    cached = {}
    if cache_file.exists():
        for entry in json.loads(cache_file.read_text()):
            cached[entry["thickness"]] = rpk.BuildupResult.from_dict(entry["result"])
        print(f"Loaded {len(cached)} cached thicknesses for {name}")

    missing = [t for t in mc_thicknesses if t not in cached]
    if missing:
        print(f"Running coupled MC for {name} at {missing} cm...")
        layers_list = [
            [rpk.Layer(thickness=VOID_THICKNESS), rpk.Layer(thickness=t, material=mat)]
            for t in missing
        ]
        mc_list = rpk.compute_buildup(
            geometries=layers_list,
            source=source,
            quantities=[N_DOSE, P_DOSE],
            particles_per_batch=10_000,
            max_batches=100,
            trigger_rel_err=0.05,
        )
        for t, r in zip(missing, mc_list):
            cached[t] = r

        cache_data = [
            {"thickness": t, "result": cached[t].to_dict()} for t in sorted(cached)
        ]
        cache_file.write_text(json.dumps(cache_data, indent=2))
        print(f"  Saved {len(cached)} thicknesses to {cache_file}")
    else:
        print(f"All thicknesses cached for {name}")

    # Compute total buildup: B = MC_total / PK_neutron (no analytical secondary photon)
    br_list = []
    for t in mc_thicknesses:
        r = cached[t]
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
        print(f"  {name}: t={t:>3d} cm, B_total={br.buildup['total']:.3f}")

    tables[name] = rpk.BuildupTable(
        points=[{"thickness": t} for t in mc_thicknesses], results=br_list
    )
    mc_cached[name] = cached

# --- Step 2: Extrapolate total dose ---
print("\nExtrapolating...")

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
            source_strength=SOURCE_STRENGTH,
            layers=layers,
            source=source,
            geometry=GEOMETRY,
        )
        bi = table.interpolate(thickness=t, warn=False)
        doses.append(pk.dose_rate * bi.value)
        doses_lo.append(pk.dose_rate * (bi.value - bi.sigma))
        doses_hi.append(pk.dose_rate * (bi.value + bi.sigma))

    data[name] = {
        "doses": np.array(doses),
        "doses_lo": np.array(doses_lo),
        "doses_hi": np.array(doses_hi),
    }

# --- Step 3: Plot ---
print("\nGenerating plot...")

colors = {"Portland + 3% steel": "#1f77b4", "Magnetite + 3% steel": "#ff7f0e"}
fig, ax = plt.subplots(figsize=(10, 7))

for name in materials:
    c = colors[name]
    d = data[name]
    valid = d["doses"] > 0

    ax.plot(
        np.array(all_thicknesses)[valid],
        d["doses"][valid],
        color=c,
        linewidth=2,
        label=name,
    )
    v2 = valid & (d["doses_lo"] > 0) & (d["doses_hi"] > 0)
    ax.fill_between(
        np.array(all_thicknesses)[v2],
        d["doses_lo"][v2],
        d["doses_hi"][v2],
        color=c,
        alpha=0.15,
    )

    mc_ts = mc_thicknesses
    mc_vals = [
        (mc_cached[name][t].mc.get(N_DOSE, 0) + mc_cached[name][t].mc.get(P_DOSE, 0))
        * SOURCE_STRENGTH
        for t in mc_ts
    ]
    mc_errs = [
        np.sqrt(
            mc_cached[name][t].mc_std_dev.get(N_DOSE, 0) ** 2
            + mc_cached[name][t].mc_std_dev.get(P_DOSE, 0) ** 2
        )
        * SOURCE_STRENGTH
        for t in mc_ts
    ]
    ax.errorbar(
        mc_ts,
        mc_vals,
        yerr=mc_errs,
        fmt="o",
        color=c,
        markersize=6,
        capsize=3,
        zorder=5,
    )

ax.set_xlabel("Shield thickness (cm)", fontsize=13)
ax.set_ylabel("Total dose rate (Sv/hr)", fontsize=13)
ax.set_title(
    f"Total dose (neutron + secondary \u03b3) \u2014 {source.energy/1e6:.1f} MeV, {GEOMETRY}\n"
    f"(S={SOURCE_STRENGTH:.0e} n/s, {VOID_THICKNESS/100:.0f} m void + shield)",
    fontsize=12,
)
ax.set_yscale("log")
ax.legend(fontsize=11)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()

output_path = RESULTS_DIR / "total_dose.png"
fig.savefig(output_path, dpi=150, bbox_inches="tight")
print(f"Plot saved to {output_path}")
