"""Photon dose across many thicknesses with GP-extrapolated buildup.

MC at 4 thin concrete thicknesses, BuildupTable GP extrapolation to 4 m.
"""

import os
import matplotlib
import matplotlib.pyplot as plt
import rad_point_kernel as pkc

matplotlib.use("Agg")

concrete = pkc.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3,
    fraction="mass",
)
source = pkc.Source("photon", 662e3)
SOURCE_STRENGTH = 1e12

mc_thicknesses = [5, 10, 15, 20]
all_thicknesses = mc_thicknesses + list(range(25, 410, 5))

# --- MC on thin shields ---
print(f"Running MC for {len(mc_thicknesses)} thicknesses...")
mc_geometries = [[pkc.Layer(thickness=t, material=concrete)] for t in mc_thicknesses]

mc_results = pkc.compute_buildup(
    geometries=mc_geometries,
    source=source,
    quantities=["dose-AP"],
    particles_per_batch=10_000,
    max_batches=50,
    trigger_rel_err=0.05,
)

for t, r in zip(mc_thicknesses, mc_results):
    print(f"  {t:>2d} cm: B = {r.buildup['dose-AP']:.3f}")

# --- BuildupTable ---
table = pkc.BuildupTable(
    points=[{"thickness": t} for t in mc_thicknesses],
    results=mc_results,
)

# --- Dose at all thicknesses ---
all_geometries = [[pkc.Layer(thickness=t, material=concrete)] for t in all_thicknesses]
pk_b_doses, pk_b_lo, pk_b_hi = [], [], []

for t, layers in zip(all_thicknesses, all_geometries):
    bi = table.interpolate(thickness=t)
    pk = pkc.calculate_dose(SOURCE_STRENGTH, layers, source, "AP")
    pk_b_doses.append(pk.dose_rate * bi.value)
    pk_b_lo.append(pk.dose_rate * (bi.value - bi.sigma))
    pk_b_hi.append(pk.dose_rate * (bi.value + bi.sigma))

# --- Plot ---
fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(all_thicknesses, pk_b_doses, "b-", linewidth=2, label="PK with buildup")
ax.fill_between(all_thicknesses, pk_b_lo, pk_b_hi, color="blue", alpha=0.15)

mc_doses = [r.mc["dose-AP"] * SOURCE_STRENGTH for r in mc_results]
mc_errs = [r.mc_std_dev["dose-AP"] * SOURCE_STRENGTH for r in mc_results]
ax.errorbar(
    mc_thicknesses,
    mc_doses,
    yerr=mc_errs,
    fmt="ko",
    markersize=7,
    capsize=4,
    zorder=5,
    label="Monte Carlo simulation",
)

ax.set_xlabel("Concrete thickness (cm)", fontsize=13)
ax.set_ylabel("Photon dose rate (Sv/hr)", fontsize=13)
ax.set_title("Photon dose — 662 keV (Cs-137), AP geometry", fontsize=13)
ax.set_yscale("log")
ax.legend(fontsize=11)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()

output_path = os.path.join(
    os.path.dirname(__file__), "8_single_layer_many_thicknesses_photon_dose.png"
)
fig.savefig(output_path, dpi=150, bbox_inches="tight")
print(f"\nPlot saved to {output_path}")
