"""Generate the two buildup comparison plots for docs/theory/buildup_factors.md.

Produces:
  - docs/assets/dose_buildup_comparison.png
  - docs/assets/flux_buildup_comparison.png

Each plot compares:
  - PK uncollided (dashed, no buildup)  - underestimates
  - PK * B with fit-extrapolated MC buildup (solid)
  - MC reference points

Uses a 662 keV photon source through concrete; the same "canonical" setup
as the dose_buildup_comparison.png shown in the buildup factors theory page.

Nuclear data:
  pip install openmc_data
  download_endf -r b8.1
  export OPENMC_CROSS_SECTIONS=~/nuclear_data/cross_sections.xml
"""

import json
import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import rad_point_kernel as rpk

matplotlib.use("Agg")

# Setup
PARTICLES_PER_HOUR = 1e12 * 3600  # 1e12 photons/sec activity, scale to /hr
PARTICLES_PER_SECOND = 1e12  # for the flux plot, label /s
source = rpk.Source(particle="photon", energy=662e3)
concrete = rpk.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3,
    fraction="mass",
)

mc_thicknesses = [5, 10, 15, 20, 30]
all_thicknesses = mc_thicknesses + list(range(35, 205, 5))

ASSETS_DIR = Path(__file__).resolve().parent
CACHE_FILE = ASSETS_DIR / "buildup_comparison_cache.json"

# Step 1: MC with cache
if CACHE_FILE.exists():
    cached = [rpk.BuildupResult.from_dict(d) for d in json.loads(CACHE_FILE.read_text())]
    mc_results = cached
    print(f"Loaded {len(cached)} cached MC points from {CACHE_FILE.name}")
else:
    print(f"Running MC at {len(mc_thicknesses)} thicknesses (photon, 662 keV, concrete)...")
    mc_geometries = [[rpk.Layer(thickness=t, material=concrete)] for t in mc_thicknesses]
    mc_results = rpk.compute_buildup(
        geometries=mc_geometries,
        source=source,
        quantities=["flux", "dose-AP"],
        particles_per_batch=10_000,
        max_batches=50,
        trigger_rel_err=0.05,
    )
    CACHE_FILE.write_text(json.dumps([r.to_dict() for r in mc_results], indent=2))
    print(f"Saved cache to {CACHE_FILE.name}")

for t, r in zip(mc_thicknesses, mc_results):
    print(f"  {t:>3d} cm: B_flux={r.buildup['flux']}, B_dose={r.buildup['dose-AP']}")


def _build_curve(quantity: str):
    """Return pk_uncoll, pk_b, mc_vals, mc_errs."""
    fit = rpk.BuildupFit(
        points=[{"thickness": t} for t in mc_thicknesses],
        results=mc_results,
    )
    pk_uncoll = []
    pk_b = []
    if quantity == "flux":
        strength = PARTICLES_PER_SECOND
    else:
        strength = PARTICLES_PER_HOUR
    for t in all_thicknesses:
        layers = [rpk.Layer(thickness=t, material=concrete)]
        bi = fit.interpolate(thickness=t, quantity=quantity, warn=False)
        if quantity == "flux":
            pk = rpk.calculate_flux(layers=layers, source=source).scale(
                strength=strength,
            )
            pk_val = pk.flux
        else:
            pk = rpk.calculate_dose(
                layers=layers, source=source, geometry="AP",
            ).scale(strength=strength)
            pk_val = pk.dose
        pk_uncoll.append(pk_val)
        pk_b.append(pk_val * bi.value)
    mc_scaled = [r.scale(strength=strength) for r in mc_results]
    mc_vals = [r.mc[quantity] for r in mc_scaled]
    mc_errs = [r.mc_std_dev[quantity] for r in mc_scaled]
    return pk_uncoll, pk_b, mc_vals, mc_errs


def _plot(quantity: str, ylabel: str, title: str, outfile: str):
    pk_u, pk_b, mc_v, mc_e = _build_curve(quantity)
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(all_thicknesses, pk_u, "k--", linewidth=2, label="PK uncollided (no buildup)")
    ax.plot(all_thicknesses, pk_b, "b-", linewidth=2, label="PK * B (MC-calibrated)")
    ax.errorbar(
        mc_thicknesses, mc_v, yerr=mc_e, fmt="ro", markersize=7, capsize=4,
        ecolor="red", zorder=5, label="Monte Carlo",
    )
    ax.set_xlabel("Concrete thickness (cm)", fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(title, fontsize=13)
    ax.set_yscale("log")
    ax.legend(fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    path = ASSETS_DIR / outfile
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


_plot(
    quantity="dose-AP",
    ylabel="Dose rate (Sv/hr)",
    title="Effect of buildup on dose - 662 keV photon through concrete (AP)",
    outfile="dose_buildup_comparison.png",
)
_plot(
    quantity="flux",
    ylabel="Flux (photons/cm^2/s)",
    title="Effect of buildup on flux - 662 keV photon through concrete",
    outfile="flux_buildup_comparison.png",
)
