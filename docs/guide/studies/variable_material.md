# Variable material study

A 50 cm concrete shield with two composition knobs:

- **Steel volume fraction** in the concrete mix (rebar/steel ball loading)
- **Boron mass fraction** added to that mix (neutron capture additive)

Source is 14.1 MeV D-T neutrons. The study shows three dose maps across the 4x4 composition grid: primary neutron dose, secondary photon dose, and total coupled dose. Steel increases density and attenuates secondary gammas; boron captures thermalised neutrons and suppresses both the neutron dose and the secondary-gamma source.

## Setup

```python exec="true" source="material-block" session="varmat" result="text"
import json
from pathlib import Path
import rad_point_kernel as rpk

PARTICLES_PER_SECOND = 1e12  # steady-state DT generator
GEOMETRY = "AP"
VOID_THICKNESS = 1000
CONCRETE_THICKNESS = 50
source = rpk.Source(particle="neutron", energy=14.1e6)

N_DOSE = f"dose-{GEOMETRY}"
P_DOSE = f"dose-{GEOMETRY}-coupled-photon"

STEEL_VOLS = [0.0, 0.02, 0.05, 0.10]
BORON_WTS = [0.0, 0.01, 0.025, 0.05]

portland = rpk.Material(
    composition={
        "H": 0.168753, "C": 0.001416, "O": 0.562525, "Na": 0.011838,
        "Mg": 0.0014, "Al": 0.021354, "Si": 0.204119, "K": 0.005656,
        "Ca": 0.018674, "Fe": 0.004264,
    },
    density=2.3,
    fraction="atom",
)
steel = rpk.Material(
    composition={"Fe": 0.991634, "C": 0.003747, "Si": 0.003163,
                 "Mn": 0.000154, "P": 0.001127, "S": 0.000175},
    density=7.7,
    fraction="atom",
)
boron = rpk.Material(composition={"B": 1.0}, density=2.34, fraction="mass")

def mix_recipe(steel_vol, boron_wt):
    if steel_vol > 0:
        base = rpk.Material.volume_mix(portland, 1.0 - steel_vol, steel, steel_vol)
    else:
        base = portland
    if boron_wt > 0:
        return rpk.Material.mass_mix(base, 1.0 - boron_wt, boron, boron_wt)
    return base

print(f"Grid: {len(STEEL_VOLS)} x {len(BORON_WTS)} = {len(STEEL_VOLS) * len(BORON_WTS)} MC points")
```

## Monte Carlo on the composition grid

Results are cached to `docs/assets/studies/variable_material_cache.json`; delete the file to re-simulate.

```python exec="true" source="material-block" session="varmat" result="text"
CACHE = Path("docs/assets/studies/variable_material_cache.json")
CACHE.parent.mkdir(parents=True, exist_ok=True)

cached = {}
if CACHE.exists():
    for entry in json.loads(CACHE.read_text()):
        key = (entry["steel_vol"], entry["boron_wt"])
        cached[key] = rpk.BuildupResult.from_dict(entry["result"])
    print(f"Loaded {len(cached)} cached points")

missing = [
    (sv, bw)
    for sv in STEEL_VOLS
    for bw in BORON_WTS
    if (sv, bw) not in cached
]

if missing:
    print(f"Running coupled MC for {len(missing)} compositions...")
    geometries = []
    for sv, bw in missing:
        mat = mix_recipe(sv, bw)
        geometries.append([
            rpk.Layer(thickness=VOID_THICKNESS),
            rpk.Layer(thickness=CONCRETE_THICKNESS, material=mat),
        ])
    new_results = rpk.compute_buildup(
        geometries=geometries,
        source=source,
        quantities=[N_DOSE, P_DOSE],
        particles_per_batch=10_000,
        max_batches=100,
        trigger_rel_err=0.05,
    )
    for (sv, bw), r in zip(missing, new_results):
        cached[(sv, bw)] = r

    cache_data = [
        {"steel_vol": sv, "boron_wt": bw, "result": cached[(sv, bw)].to_dict()}
        for sv, bw in sorted(cached)
    ]
    CACHE.write_text(json.dumps(cache_data, indent=2))
    print(f"Saved {len(cached)} points to {CACHE}")
else:
    print("All grid points cached")
```

## Dose maps

```python exec="true" source="material-block" session="varmat" html="true"
import io
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")

neutron = np.zeros((len(BORON_WTS), len(STEEL_VOLS)))
photon = np.zeros_like(neutron)
for i, bw in enumerate(BORON_WTS):
    for j, sv in enumerate(STEEL_VOLS):
        r = cached[(sv, bw)]
        neutron[i, j] = r.mc[N_DOSE] * PARTICLES_PER_SECOND
        photon[i, j] = r.mc[P_DOSE] * PARTICLES_PER_SECOND
total = neutron + photon

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
panels = [
    (axes[0], neutron, "Neutron dose rate (Sv/s)"),
    (axes[1], photon, "Secondary photon dose rate (Sv/s)"),
    (axes[2], total, "Total dose rate (Sv/s)"),
]

extent = [
    100 * STEEL_VOLS[0] - 1, 100 * STEEL_VOLS[-1] + 1,
    100 * BORON_WTS[0] - 0.5, 100 * BORON_WTS[-1] + 0.5,
]

for ax, data, title in panels:
    im = ax.imshow(
        data, origin="lower", aspect="auto", cmap="viridis",
        norm=matplotlib.colors.LogNorm(vmin=data.min(), vmax=data.max()),
        extent=extent,
    )
    ax.set_xticks([100 * v for v in STEEL_VOLS])
    ax.set_yticks([100 * v for v in BORON_WTS])
    ax.set_xlabel("Steel (vol %)")
    ax.set_ylabel("Boron (wt %)")
    ax.set_title(title)
    for i, bw in enumerate(BORON_WTS):
        for j, sv in enumerate(STEEL_VOLS):
            ax.text(100 * sv, 100 * bw, f"{data[i, j]:.2g}",
                    ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

fig.suptitle(
    f"{CONCRETE_THICKNESS} cm Portland concrete, 14.1 MeV neutron, {GEOMETRY}, "
    f"S = {PARTICLES_PER_SECOND} n/s",
    fontsize=11,
)

buf = io.StringIO()
fig.savefig(buf, format="svg", bbox_inches="tight")
plt.close(fig)
print(buf.getvalue())
```

## Reading the maps

- **Neutron dose** drops steeply with boron; steel alone is nearly neutral (iron has a fast-neutron removal cross section similar to concrete's hydrogen-driven attenuation, just at higher density).
- **Secondary photon dose** drops with boron (fewer thermal captures to produce capture gammas) and drops with steel (iron attenuates the secondary gammas on their way out).
- **Total dose** is dominated by whichever leg wins at that composition - at low boron the secondary-gamma leg is often comparable to or larger than the neutron leg.
