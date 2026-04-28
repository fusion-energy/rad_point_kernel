# Caching Monte Carlo results

Monte Carlo simulations are expensive. `BuildupResult` provides `save` and `load` methods to cache results to JSON files, so you only run Monte Carlo once.

For a flat list of `BuildupResult`, prefer `BuildupResult.save([...])` and `BuildupResult.load(path)`. The lower-level `to_dict()` and `from_dict()` are for when you need to embed a `BuildupResult` inside your own JSON envelope (for example, alongside per-result metadata that the library doesn't track for you).

## Save and load

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

# A list of BuildupResults from a prior compute_buildup call
iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
mc_results = rpk.compute_buildup(
    geometries=[[rpk.Layer(thickness=5, material=iron)]],
    source=rpk.Source(particle="photon", energy=1e6),
    quantities=["dose-AP"],
)

# Save to JSON
rpk.BuildupResult.save(mc_results, "buildup_cache.json")

# Load them back
mc_results = rpk.BuildupResult.load("buildup_cache.json")
print(f"Loaded {len(mc_results)} BuildupResult(s) from cache")
```

## Cache-or-compute pattern

Check if a cache file exists before running MC:

```python exec="true" source="material-block" result="text"
from pathlib import Path
import rad_point_kernel as rpk

concrete = rpk.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)

mc_thicknesses = [5, 10, 15, 20]
CACHE = Path("mc_cache.json")

if CACHE.exists():
    print("Loading cached Monte Carlo results...")
    mc_results = rpk.BuildupResult.load(CACHE)
else:
    print("Running Monte Carlo...")
    geometries = [
        [rpk.Layer(thickness=t, material=concrete)]
        for t in mc_thicknesses
    ]
    source = rpk.Source(particle="photon", energy=1e6)
    mc_results = rpk.compute_buildup(
        geometries=geometries,
        source=source,
        quantities=["dose-AP"],
    )
    rpk.BuildupResult.save(mc_results, CACHE)
    print(f"Saved to {CACHE}")

for t, r in zip(mc_thicknesses, mc_results):
    print(f"  {t:>2d} cm: B = {r.buildup['dose-AP']}")
```

## Incremental caching

When you need to add more thicknesses to an existing study, load the cache, run Monte Carlo only for the missing points, and save the combined results. `BuildupResult.save/load` is for the result payload; a small sidecar JSON keeps the thickness-keyed metadata aligned by position.

```python exec="true" source="material-block" result="text"
import json
from pathlib import Path
import rad_point_kernel as rpk

concrete = rpk.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)

mc_thicknesses = [5, 10, 15, 20, 30, 40]
CACHE = Path("mc_incremental.json")
KEYS = CACHE.with_suffix(".thicknesses.json")

# Load existing cache (results + parallel list of thicknesses)
if CACHE.exists() and KEYS.exists():
    cached_thicknesses = json.loads(KEYS.read_text())
    cached_results = rpk.BuildupResult.load(CACHE)
    print(f"Loaded {len(cached_thicknesses)} cached thicknesses")
else:
    cached_thicknesses, cached_results = [], []

cached = dict(zip(cached_thicknesses, cached_results))
missing = [t for t in mc_thicknesses if t not in cached]

if missing:
    print(f"Running Monte Carlo for {missing}...")
    new_results = rpk.compute_buildup(
        geometries=[[rpk.Layer(thickness=t, material=concrete)] for t in missing],
        source=rpk.Source(particle="photon", energy=1e6),
        quantities=["dose-AP"],
    )
    cached.update(zip(missing, new_results))

    # Re-save: results + parallel sidecar of keys
    keys_sorted = sorted(cached)
    rpk.BuildupResult.save([cached[t] for t in keys_sorted], CACHE)
    KEYS.write_text(json.dumps(keys_sorted))
    print(f"Saved {len(cached)} thicknesses")
else:
    print("All thicknesses already cached")

for t in sorted(cached):
    print(f"  {t:>2d} cm: B = {cached[t].buildup['dose-AP']}")
```

## Caches and version upgrades

`BuildupResult` round-trips through JSON faithfully, but a cache only contains the quantities that existed when it was written. If a later release adds new quantities (for example `dose-{geo}-total` was added in 1.4.0 alongside `dose-{geo}` and `dose-{geo}-coupled-photon`), loading an older cache gives you a `BuildupResult` *without* the new keys. Two consequences to watch for:

- Reading a missing key by name raises `KeyError` (e.g. `r.mc["dose-AP-total"]` on a 1.3.x cache).
- Building a `BuildupTable` from an older cache produces a table whose `available_quantities` reflects only what was tallied at the time, so `interpolate()` may pick a different default than you expect.

The simplest cure after a minor version bump is to **delete the cache JSON and re-run Monte Carlo**. The cache filename is yours to manage; if your workflow can tolerate it, namespacing by version (e.g. `mc_cache_v1.4.json`) lets old and new caches coexist while you transition.
