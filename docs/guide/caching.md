# Caching Monte Carlo results

Monte Carlo simulations are expensive. `BuildupResult` provides `save` and `load` methods to cache results to JSON files, so you only run Monte Carlo once.

## Save and load

```python
import rad_point_kernel as rpk

# Save a list of BuildupResults
rpk.BuildupResult.save(mc_results, "buildup_cache.json")

# Load them back
mc_results = rpk.BuildupResult.load("buildup_cache.json")
```

## Cache-or-compute pattern

Check if a cache file exists before running MC:

```python
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
    print(f"  {t:>2d} cm: B = {r.buildup['dose-AP']:.3f}")
```

## Incremental caching

When you need to add more thicknesses to an existing study, load the cache, run Monte Carlo only for the missing points, and save the combined results.

```python
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

# Load existing cache as {thickness: BuildupResult}
cached = {}
if CACHE.exists():
    for entry in json.loads(CACHE.read_text()):
        cached[entry["thickness"]] = rpk.BuildupResult.from_dict(entry["result"])
    print(f"Loaded {len(cached)} cached thicknesses")

# Find missing thicknesses
missing = [t for t in mc_thicknesses if t not in cached]

if missing:
    print(f"Running Monte Carlo for {missing}...")
    geometries = [
        [rpk.Layer(thickness=t, material=concrete)]
        for t in missing
    ]
    source = rpk.Source(particle="photon", energy=1e6)
    new_results = rpk.compute_buildup(
        geometries=geometries,
        source=source,
        quantities=["dose-AP"],
    )
    for t, r in zip(missing, new_results):
        cached[t] = r

    # Save all results
    cache_data = [
        {"thickness": t, "result": cached[t].to_dict()}
        for t in sorted(cached)
    ]
    CACHE.write_text(json.dumps(cache_data, indent=2))
    print(f"Saved {len(cached)} thicknesses")
else:
    print("All thicknesses already cached")

# Use all cached results
for t in sorted(cached):
    print(f"  {t:>2d} cm: B = {cached[t].buildup['dose-AP']:.3f}")
```
