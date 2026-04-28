
import json
from pathlib import Path
import numpy as np
import rad_point_kernel as rpk

PARTICLES_PER_HOUR = 1e12 * 3600
GEOMETRY = "AP"
VOID_THICKNESS = 1000
source = rpk.Source(particle="neutron", energy=14.06e6)

N_DOSE = f"dose-{GEOMETRY}"
P_DOSE = f"dose-{GEOMETRY}-coupled-photon"

polyethylene = rpk.Material(composition={"H": 0.143, "C": 0.857}, density=0.94, fraction="mass")
concrete = rpk.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3,
    fraction="mass",
)

mc_poly = [0, 5, 10, 15, 20, 25]
mc_conc = [5, 25, 50, 100, 200]
all_poly = list(range(0, 30, 5))
all_conc = list(range(5, 205, 5))

CACHE = Path("multi_layer_cache.json")
CACHE.parent.mkdir(parents=True, exist_ok=True)


def make_layers(poly_t, conc_t):
    layers = [rpk.Layer(thickness=VOID_THICKNESS)]
    if poly_t > 0:
        layers.append(rpk.Layer(thickness=poly_t, material=polyethylene))
    layers.append(rpk.Layer(thickness=conc_t, material=concrete))
    return layers

cached = {}
if CACHE.exists():
    for entry in json.loads(CACHE.read_text()):
        cached[(entry["poly"], entry["conc"])] = rpk.BuildupResult.from_dict(entry["result"])

missing = [(w, c) for w in mc_poly for c in mc_conc if (w, c) not in cached]
if missing:
    mc_geometries = [make_layers(w, c) for w, c in missing]
    new_results = rpk.compute_buildup(
        geometries=mc_geometries,
        source=source,
        quantities=[N_DOSE, P_DOSE],
        particles_per_batch=100,
        max_batches=100,
        trigger_rel_err=0.05,
        cross_sections="/home/jon/nuclear_data/cross_sections.xml",
    )
    for (w, c), r in zip(missing, new_results):
        cached[(w, c)] = r
    cache_data = [
        {"poly": w, "conc": c, "result": cached[(w, c)].to_dict()}
        for w, c in sorted(cached)
    ]
    CACHE.write_text(json.dumps(cache_data, indent=2))


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

tables_by_poly = {
    wt: rpk.BuildupTable(
        points=[{"conc": ct} for ct in mc_conc],
        results=[total_buildup_result(cached[(wt, ct)]) for ct in mc_conc],
    )
    for wt in mc_poly
}
tables_by_conc = {
    ct: rpk.BuildupTable(
        points=[{"poly": wt} for wt in mc_poly],
        results=[total_buildup_result(cached[(wt, ct)]) for wt in mc_poly],
    )
    for ct in mc_conc
}

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
            bi = table.interpolate(poly=t, warn=False)
        pk = rpk.calculate_dose(
            layers=layers,
            source=source,
            geometry=GEOMETRY,
        ).scale(strength=PARTICLES_PER_HOUR)
        doses.append(pk.dose * bi.value)
        lo.append(pk.dose * (bi.value - bi.sigma))
        hi.append(pk.dose * (bi.value + bi.sigma))
    return doses, lo, hi


data_vs_conc = {
    wt: extrapolate(tables_by_poly[wt], all_conc, wt, "conc") for wt in mc_poly
}
data_vs_poly = {
    ct: extrapolate(tables_by_conc[ct], all_poly, ct, "poly") for ct in mc_conc
}

colors_poly = {0: "#7f7f7f", 5: "#1f77b4", 10: "#ff7f0e",
                15: "#2ca02c", 20: "#d62728", 25: "#9467bd"}
colors_conc = {5: "#1f77b4", 25: "#ff7f0e", 50: "#2ca02c",
               100: "#d62728", 200: "#9467bd"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

for wt in mc_poly:
    c = colors_poly[wt]
    d, lo, hi = data_vs_conc[wt]
    ax1.plot(all_conc, d, color=c, linewidth=2, label=f"{wt} cm polyethylene")
    ax1.fill_between(all_conc, lo, hi, color=c, alpha=0.12)
    mc_vals = [
        (cached[(wt, ct)].mc.get(N_DOSE, 0) + cached[(wt, ct)].mc.get(P_DOSE, 0))
        * PARTICLES_PER_HOUR for ct in mc_conc
    ]
    ax1.scatter(mc_conc, mc_vals, color=c, s=30, zorder=5)

ax1.set_xlabel("Concrete thickness (cm)")
ax1.set_ylabel("Total dose rate (Sv/hr)")
ax1.set_yscale("log")
ax1.set_title("Total dose vs concrete (per polyethylene thickness)")
ax1.legend(fontsize=9)
ax1.grid(True, which="both", alpha=0.3)

for ct in mc_conc:
    c = colors_conc[ct]
    d, lo, hi = data_vs_poly[ct]
    ax2.plot(all_poly, d, color=c, linewidth=2, label=f"{ct} cm concrete")
    ax2.fill_between(all_poly, lo, hi, color=c, alpha=0.12)
    mc_vals = [
        (cached[(wt, ct)].mc.get(N_DOSE, 0) + cached[(wt, ct)].mc.get(P_DOSE, 0))
        * PARTICLES_PER_HOUR for wt in mc_poly
    ]
    ax2.scatter(mc_poly, mc_vals, color=c, s=30, zorder=5)

ax2.set_xlabel("Polyethylene thickness (cm)")
ax2.set_ylabel("Total dose rate (Sv/hr)")
ax2.set_yscale("log")
ax2.set_title("Total dose vs polyethylene (per concrete thickness)")
ax2.legend(fontsize=9)
ax2.grid(True, which="both", alpha=0.3)

fig.tight_layout()
fig.savefig("total_dose.png", bbox_inches="tight")
plt.close(fig)
