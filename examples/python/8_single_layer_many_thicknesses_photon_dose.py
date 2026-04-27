"""Photon dose across many thicknesses with GP-extrapolated buildup.

MC at 4 thin concrete thicknesses, BuildupTable GP extrapolation to 4 m.
Cs-137 is a continuous source; the strength is supplied in photons per hour
so the resulting dose is in Sv/hr.
"""

import os
import matplotlib
import matplotlib.pyplot as plt
import rad_point_kernel as rpk

matplotlib.use("Agg")

concrete = rpk.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3,
    fraction="mass",
)
source = rpk.Source("photon", 662e3)
PARTICLES_PER_HOUR = 1e12 * 3600  # 1e12 photons/sec activity

mc_thicknesses = [5, 10, 15, 20]
all_thicknesses = mc_thicknesses + list(range(25, 410, 5))

# MC on thin shields
print(f"Running MC for {len(mc_thicknesses)} thicknesses...")
mc_geometries = [[rpk.Layer(thickness=t, material=concrete)] for t in mc_thicknesses]

mc_results = rpk.compute_buildup(
    geometries=mc_geometries,
    source=source,
    quantities=["dose-AP"],
    particles_per_batch=10_000,
    max_batches=50,
    trigger_rel_err=0.05,
)

for t, r in zip(mc_thicknesses, mc_results):
    print(f"  {t} cm: B = {r.buildup['dose-AP']}")

# BuildupTable
table = rpk.BuildupTable(
    points=[{"thickness": t} for t in mc_thicknesses],
    results=mc_results,
)

# Dose at all thicknesses
all_geometries = [[rpk.Layer(thickness=t, material=concrete)] for t in all_thicknesses]
pk_b_doses, pk_b_lo, pk_b_hi = [], [], []

for t, layers in zip(all_thicknesses, all_geometries):
    bi = table.interpolate(thickness=t)
    pk = rpk.calculate_dose(layers=layers, source=source, geometry="AP").scale(
        strength=PARTICLES_PER_HOUR
    )
    pk_b_doses.append(pk.dose * bi.value)
    pk_b_lo.append(pk.dose * (bi.value - bi.sigma))
    pk_b_hi.append(pk.dose * (bi.value + bi.sigma))

# Plot
fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(all_thicknesses, pk_b_doses, "b-", linewidth=2, label="PK with buildup")
ax.fill_between(all_thicknesses, pk_b_lo, pk_b_hi, color="blue", alpha=0.15)

mc_scaled = [r.scale(strength=PARTICLES_PER_HOUR) for r in mc_results]
mc_doses = [r.mc["dose-AP"] for r in mc_scaled]
mc_errs = [r.mc_std_dev["dose-AP"] for r in mc_scaled]
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
ax.set_title("Photon dose - 662 keV (Cs-137), AP geometry", fontsize=13)
ax.set_yscale("log")
ax.legend(fontsize=11)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()

output_path = os.path.join(
    os.path.dirname(__file__), "8_single_layer_many_thicknesses_photon_dose.png"
)
fig.savefig(output_path, dpi=150, bbox_inches="tight")
print(f"\nPlot saved to {output_path}")
