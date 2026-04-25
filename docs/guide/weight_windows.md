# Weight window generation for Monte Carlo

When `compute_buildup` runs OpenMC, it often has to push particles through many
mean free paths of shielding. Without variance reduction, most particles are
absorbed before reaching the detector surface and the tally sits at zero. Even
with long run times, the relative error stays unusable.

`rad_point_kernel` generates **weight windows (WW)** automatically for every
MC run. The WW come from the same point-kernel physics the rest of the
package uses, so no extra MC pre-iteration is needed — they are built
analytically in milliseconds before the MC starts.

## What weight windows do

A weight window is an importance map over space and energy. When a particle
enters a region of higher importance (deeper into the shield), the MC code
splits it into multiple daughters with smaller weights — many particles now
sample the region that would otherwise see almost no statistics. When a
particle enters a region of lower importance, it is Russian-rouletted,
probabilistically killed (with weight compensated on the survivors) so that
CPU is not spent on histories that won't contribute to the tally.

Splitting and rouletting are *unbiased* — they change variance, not the
expected value of the tally.

## How `rad_point_kernel` chooses the map

Three things need to be decided:

1. **The radial mesh** — where to put the boundaries of the weight-window
   cells in radius.
2. **The per-cell importance values** — the lower weight-window bound in each
   cell and energy group.
3. **The energy bin structure** — so that fast and slow particles can be
   biased differently.

All three are driven by the **point-kernel importance curve** for the
quantity being tallied (flux, dose, or coupled secondary-photon dose). The
point-kernel is evaluated along the radius, and the mesh is placed so that
each cell represents roughly a factor of $e$ drop in expected tally
contribution. Cell values are the point-kernel importance evaluated at each
bin centre, normalised so the innermost bin has lower bound 1.

### Layer interfaces are always mesh edges

Voids and material changes get a mesh boundary automatically. A void cell
stays as a single cell (no attenuation inside it, so there's nothing to
bias). A thin material layer — thinner than one mean free path — also gets
exactly one cell. Thick layers get subdivided into $\lceil \tau / \ln(e)
\rceil$ cells of roughly equal optical thickness.

### Energy bins come from the source

The builder chooses log-spaced energy bins automatically from the peak
source energy, with a thermal floor for neutrons and a keV floor for
photons. Per-energy bounds are evaluated as if the source were
monoenergetic at each bin's geometric mean.

For coupled neutron → photon simulations, a second `WeightWindows` object
is built for the photon population, using the secondary-gamma importance
curve.

## Limitations — why this isn't "perfect"

The weight windows are built from the **bare attenuation** physics — the
same $\exp(-\Sigma_r \, t)$ attenuation the point-kernel computes — **with
the build-up factor set to 1**. The build-up factor isn't used during WW
construction because:

- If the WW map already incorporated build-up, it would flatten the
  importance curve and split too little at depth.
- The MC run itself is what measures the true build-up.

This means the analytic WW **underestimate the flux** at depth when
scattering (build-up) is significant — e.g. dose behind a hydrogenous
shield. In practice the WW still produce enormous speedups (often 10× to
1000×) because the $\exp(-\Sigma_r t)$ shape dominates, and the under-bias
only costs a modest amount of variance efficiency on the scattered tail.

Similar caveats apply to fissile materials (point-kernel doesn't model
neutron multiplication) and to strongly-downscattered spectra (the
per-energy bound uses a mono-at-bin-centre approximation). For those cases,
MAGIC-style iterative refinement seeded from the analytic WW is a future
upgrade path — but for the fast-neutron / photon shielding problems
`rad_point_kernel` targets, the analytic WW is "good enough" and costs
nothing.

## Controlling WW generation

Weight windows are **on by default** in `compute_buildup`. Turn them off
with:

```python
results = rpk.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP"],
    use_weight_windows=False,   # default True
)
```

The builder skips itself when it wouldn't help — if the geometry has
total $\tau < 2$ (trivial attenuation) or produces a mesh with $\leq 2$
cells, it returns an empty WW list and `compute_buildup` just runs plain
MC. An INFO log is emitted so you know this happened.

For advanced control, the builder is callable directly:

```python
from rad_point_kernel.weight_windows import build_weight_windows

ww_list = build_weight_windows(
    layers=layers,
    source=source,
    quantities=["dose-AP"],
    log_ratio_per_bin=1.0,         # factor of e per bin (default)
    min_bin_width_cm=0.5,
    upper_bound_ratio=5.0,
)
# Attach to an OpenMC settings object directly:
my_openmc_settings.weight_windows = ww_list
```

## When to tally both flux and dose

If `quantities` contains more than one tally (e.g. `["flux", "dose-AP"]`),
the builder uses the importance curve that drops the fastest across the
shield to drive the WW. That curve needs the most variance help; WW sized
for it is automatically adequate for the gentler quantity. The cost is a
modest efficiency reduction on the gentler tally, never variance
starvation on the steeper one.
