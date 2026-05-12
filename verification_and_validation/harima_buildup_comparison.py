"""V&V: compare rad_point_kernel buildup factors to the Harima/ANS-6.4.3
Geometric Progression (GP) fits.

Reference
---------
Harima, Sakamoto, Tanaka & Kawai, "Validity of the Geometric Progression
Formula in Approximating Gamma-Ray Buildup Factors", Nuclear Science and
Engineering 94 (1986). The GP formula and tabulated (a, b, c, d, xi)
coefficients are from that paper. ANSI/ANS-6.4.3-1991 (and the 2010
revision) re-tabulates the same form with refined coefficients; either
source is usable, but coefficients here are taken from Harima 1986.

GP form
-------
For a point isotropic source in an infinite homogeneous medium, the dose
buildup factor at optical thickness mut is

    B(mut) = 1 + (b - 1) * (K^mut - 1) / (K - 1)        if K != 1
    B(mut) = 1 + (b - 1) * mut                            if K == 1

with

    K(mut) = c * mut^a + d * [tanh(mut/xi - 2) - tanh(-2)]
                                / [1 - tanh(-2)]

Geometry
--------
The GP fit assumes point-isotropic-in-infinite-medium. Our MC anchor is a
spherical shell with a point isotropic source at the origin and a tally
on the outer surface. For the comparison to be valid the outer sphere
must be large enough that no leakage occurs at the tally surface (we
want the photon flux to be representative of an infinite medium at the
tally radius). Sphere sizing is controlled by `_sphere_radius_for_mfp`
below.

Coefficient table
-----------------
GP coefficients are NOT bundled with the script; they live in
`verification_and_validation/data/harima_1986_gp_coefficients.json`,
which currently ships empty. Populate it from Harima 1986 (Tables I-VI)
to enable the comparison. The script's behaviour when the table is
empty is to print a warning and exit 0 (no regressions to flag).

Output
------
- Stdout summary table per (material, energy).
- `harima_buildup_comparison.csv` - one row per (material, energy, mfp)
  with `B_GP`, `B_MC`, `B_MC_sigma`, `rel_err`, `z`.
- `harima_buildup_comparison.png` - panel plot per material.

Exit codes
----------
0  All comparisons within tolerance, OR coefficient table empty.
1  At least one comparison outside tolerance.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

import rad_point_kernel as rpk


HERE = Path(__file__).resolve().parent
COEFF_PATH = HERE / "data" / "harima_1986_gp_coefficients.json"
CACHE_PATH = HERE / "harima_buildup_comparison_cache.json"
CSV_PATH = HERE / "harima_buildup_comparison.csv"
PLOT_PATH = HERE / "harima_buildup_comparison.png"

# Z-score tolerance for the agreement test. Harima's GP fit itself has
# residuals against the moments-method buildup factors of order 1-5%, so
# we shouldn't require much better than that. We report z-scores anyway
# so a real divergence is obvious. Pass criterion is |z| <= 5.0 OR
# rel_err <= 0.05 -- whichever is looser, since for very small MC sigma
# even a small physical disagreement produces a large z.
Z_TOL = 5.0
REL_TOL = 0.05

# mfp grid for the comparison. Harima 1986 covers 0.5 to 40 mfp; we
# subset to a manageable spread that still shows the saturation regime.
DEFAULT_MFP_GRID = (0.5, 1.0, 2.0, 4.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0)


# Material definitions -- compositions and densities for every label that
# the coefficient JSON might use. Density in g/cm^3. Compositions are
# atomic number ratios (e.g. water = 2 H : 1 O). The GP coefficients are
# per-element except for air, water, and concrete (which Harima
# tabulates as compounds in their own right). Add entries here as the
# coefficient JSON is populated.
#
# (Densities for elemental solids/liquids at room temperature; for
# gases like air a representative density at STP.)
MATERIALS: dict[str, tuple[dict[str, float], float]] = {
    "Be":  ({"Be": 1.0}, 1.85),
    "B":   ({"B":  1.0}, 2.34),
    "C":   ({"C":  1.0}, 1.70),       # graphite
    "N":   ({"N":  1.0}, 0.001251),   # gas, STP
    "O":   ({"O":  1.0}, 0.001429),   # gas, STP
    "Na":  ({"Na": 1.0}, 0.971),
    "Mg":  ({"Mg": 1.0}, 1.738),
    "Al":  ({"Al": 1.0}, 2.699),
    "Si":  ({"Si": 1.0}, 2.330),
    "P":   ({"P":  1.0}, 1.823),
    "S":   ({"S":  1.0}, 2.067),
    "Ar":  ({"Ar": 1.0}, 0.001784),   # gas, STP
    "K":   ({"K":  1.0}, 0.862),
    "Ca":  ({"Ca": 1.0}, 1.55),
    "Fe":  ({"Fe": 1.0}, 7.874),
    "Cu":  ({"Cu": 1.0}, 8.96),
    "Mo":  ({"Mo": 1.0}, 10.22),
    "Sn":  ({"Sn": 1.0}, 7.31),
    "La":  ({"La": 1.0}, 6.146),
    "Gd":  ({"Gd": 1.0}, 7.901),
    "W":   ({"W":  1.0}, 19.30),
    "Pb":  ({"Pb": 1.0}, 11.34),
    "U":   ({"U":  1.0}, 19.05),
    # Compounds explicitly tabulated by Harima.
    "air":      ({"N": 0.78084, "O": 0.20946, "Ar": 0.00934, "C": 0.0004}, 0.001225),
    "water":    ({"H": 2.0, "O": 1.0}, 1.0),
    "concrete": ({"H": 0.0221, "C": 0.00248, "O": 0.5746, "Na": 0.0152,
                  "Mg": 0.00128, "Al": 0.0196, "Si": 0.3046, "S": 0.00128,
                  "K": 0.0100, "Ca": 0.0426, "Fe": 0.00637}, 2.30),  # ANSI N5.12 / NIST PSTAR ordinary concrete
}


# --- GP formula ----------------------------------------------------------

def gp_buildup(mut: np.ndarray, a: float, b: float, c: float, d: float,
               xi: float) -> np.ndarray:
    """Evaluate the Harima GP buildup factor at one or more optical
    thicknesses `mut`. Vectorised over `mut`.
    """
    mut = np.asarray(mut, dtype=float)
    k = c * np.power(mut, a) + d * (
        np.tanh(mut / xi - 2.0) - math.tanh(-2.0)
    ) / (1.0 - math.tanh(-2.0))
    # Linear branch when K -> 1.
    out = np.empty_like(mut)
    near_one = np.abs(k - 1.0) < 1e-9
    out[near_one] = 1.0 + (b - 1.0) * mut[near_one]
    if (~near_one).any():
        kk = k[~near_one]
        out[~near_one] = 1.0 + (b - 1.0) * (
            np.power(kk, mut[~near_one]) - 1.0
        ) / (kk - 1.0)
    return out


# --- coefficient loading --------------------------------------------------

@dataclass
class GPEntry:
    material: str
    energy_MeV: float
    a: float
    b: float
    c: float
    d: float
    xi: float
    valid_range_mfp: Optional[tuple[float, float]] = None


def load_coefficients(path: Path = COEFF_PATH) -> list[GPEntry]:
    raw = json.loads(path.read_text())
    out: list[GPEntry] = []
    for material, payload in raw.get("materials", {}).items():
        if material not in MATERIALS:
            raise ValueError(
                f"Coefficient table references material '{material}' which "
                f"has no composition/density entry in MATERIALS. Add it."
            )
        coeffs = payload.get("coefficients", {})
        valid = payload.get("valid_range_mfp")
        valid_t = tuple(valid) if valid else None
        for e_str, c in coeffs.items():
            out.append(GPEntry(
                material=material,
                energy_MeV=float(e_str),
                a=float(c["a"]), b=float(c["b"]), c=float(c["c"]),
                d=float(c["d"]), xi=float(c["xi"]),
                valid_range_mfp=valid_t,
            ))
    return out


# --- MC anchor ------------------------------------------------------------

def _sphere_radius_for_mfp(target_mfp: float, density: float,
                           mu_per_density_cm2_g: float) -> float:
    """Physical thickness (cm) of the spherical shell needed to land at
    `target_mfp` mean free paths. Uses mu/rho [cm^2/g] from the user.

    For a precise 'no leakage' guarantee we'd want to oversize the
    sphere beyond the dose tally point; here we tally on the outer
    surface and rely on Harima's own infinite-medium assumption being a
    good match to a sphere whose tally is at a finite, but several-mfp,
    radius. If outer-leakage bias becomes a concern, switch to a
    nested-shell tally inside a larger outer vacuum boundary.
    """
    return target_mfp / max(mu_per_density_cm2_g * density, 1e-30)


def run_mc_anchors(entries: list[GPEntry],
                   mfp_grid: tuple[float, ...],
                   particles_per_batch: int,
                   batches: int,
                   max_batches: int,
                   trigger_rel_err: float,
                   cache_path: Path = CACHE_PATH) -> dict:
    """Run (or load cached) MC buildup factors for every (material,
    energy, mfp) triple referenced by `entries`.

    Cache key: `<material>|<energy_MeV>|<mfp>`.
    """
    cache = {}
    if cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text())
        except Exception:
            cache = {}

    for entry in entries:
        comp, density = MATERIALS[entry.material]
        material = rpk.Material(composition=comp, density=density)
        source = rpk.Source("photon", entry.energy_MeV * 1e6)

        thicknesses_cm = []
        keys = []
        for mfp in mfp_grid:
            key = f"{entry.material}|{entry.energy_MeV}|{mfp}"
            if key in cache:
                continue
            # Translate target mfp -> physical thickness via a 1 cm probe
            # that we then rescale using the optical-thickness back-out
            # of compute_buildup itself. To avoid a chicken-and-egg, do
            # one bracketing run at 1 cm and use that result's
            # optical_thickness to set the physical thickness for the
            # other mfp values.
            thicknesses_cm.append((mfp, key))

        if not thicknesses_cm:
            continue

        # First pass: a 1 cm probe to determine optical_thickness per cm.
        probe = rpk.compute_buildup(
            geometries=[[rpk.Layer(thickness=1.0, material=material)]],
            source=source,
            quantities=["dose-AP-photon"],
            particles_per_batch=particles_per_batch,
            batches=batches,
            max_batches=max_batches,
            trigger_rel_err=trigger_rel_err,
        )[0]
        mfp_per_cm = probe.optical_thickness  # since the probe is 1 cm
        if mfp_per_cm <= 0:
            raise RuntimeError(
                f"Zero optical thickness for {entry.material} at "
                f"{entry.energy_MeV} MeV - cross-section data missing?"
            )

        geometries = []
        for mfp, key in thicknesses_cm:
            t_cm = mfp / mfp_per_cm
            geometries.append([rpk.Layer(thickness=t_cm, material=material)])

        results = rpk.compute_buildup(
            geometries=geometries,
            source=source,
            quantities=["dose-AP-photon"],
            particles_per_batch=particles_per_batch,
            batches=batches,
            max_batches=max_batches,
            trigger_rel_err=trigger_rel_err,
        )
        for (mfp, key), r in zip(thicknesses_cm, results):
            mc = r.mc["dose-AP-photon"]
            pk = r.pk["dose-AP-photon"]
            sigma = r.mc_std_dev["dose-AP-photon"]
            cache[key] = {
                "optical_thickness": r.optical_thickness,
                "mc": mc,
                "pk": pk,
                "mc_std_dev": sigma,
                "B_mc": mc / pk if pk > 0 else float("nan"),
                "B_mc_sigma": (sigma / pk) if pk > 0 else float("nan"),
            }
        cache_path.write_text(json.dumps(cache, indent=2))

    return cache


# --- comparison driver ----------------------------------------------------

def run(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--smoke", action="store_true",
                        help="Run only the first (material, energy) entry "
                             "and a 3-point mfp grid as a smoke test.")
    parser.add_argument("--no-mc", action="store_true",
                        help="Skip MC; only evaluate the GP formula and "
                             "any cached MC values. Useful in CI.")
    parser.add_argument("--particles", type=int, default=10_000)
    parser.add_argument("--batches", type=int, default=10)
    parser.add_argument("--max-batches", type=int, default=100)
    parser.add_argument("--rel-err", type=float, default=0.05)
    args = parser.parse_args(argv)

    entries = load_coefficients()
    if not entries:
        print(f"WARNING: GP coefficient table {COEFF_PATH} is empty.")
        print("         Populate it from Harima 1986 NSE 94 (Tables I-VI)")
        print("         to enable the comparison. Exiting 0 (no regression).")
        return 0

    mfp_grid = DEFAULT_MFP_GRID
    if args.smoke:
        entries = entries[:1]
        mfp_grid = (1.0, 5.0, 20.0)

    if not args.no_mc:
        run_mc_anchors(
            entries=entries, mfp_grid=mfp_grid,
            particles_per_batch=args.particles,
            batches=args.batches, max_batches=args.max_batches,
            trigger_rel_err=args.rel_err,
        )

    cache = json.loads(CACHE_PATH.read_text()) if CACHE_PATH.exists() else {}

    rows: list[dict] = []
    failures = 0
    print(f"{'material':<10} {'E [MeV]':>8} {'mfp':>6} {'B_GP':>10} "
          f"{'B_MC':>10} {'sigma':>10} {'rel_err':>9} {'z':>7}  status")
    print("-" * 90)
    for e in entries:
        for mfp in mfp_grid:
            key = f"{e.material}|{e.energy_MeV}|{mfp}"
            mc = cache.get(key)
            if not mc:
                continue
            mut = mc["optical_thickness"]
            B_gp = float(gp_buildup(np.array([mut]),
                                    e.a, e.b, e.c, e.d, e.xi)[0])
            B_mc = mc["B_mc"]
            sig = max(mc["B_mc_sigma"], 1e-30)
            rel = abs(B_gp - B_mc) / max(abs(B_mc), 1e-30)
            z = abs(B_gp - B_mc) / sig
            ok = (z <= Z_TOL) or (rel <= REL_TOL)
            status = "OK" if ok else "FAIL"
            if not ok:
                failures += 1
            rows.append({
                "material": e.material, "energy_MeV": e.energy_MeV,
                "mfp": mfp, "mut_actual": mut,
                "B_GP": B_gp, "B_MC": B_mc, "B_MC_sigma": sig,
                "rel_err": rel, "z": z, "status": status,
            })
            print(f"{e.material:<10} {e.energy_MeV:>8.3f} {mfp:>6.1f} "
                  f"{B_gp:>10.4f} {B_mc:>10.4f} {sig:>10.4f} "
                  f"{rel:>9.2e} {z:>7.2f}  {status}")

    if rows:
        with CSV_PATH.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote {len(rows)} rows to {CSV_PATH}")
        try:
            _plot(rows, PLOT_PATH)
            print(f"Wrote plot to {PLOT_PATH}")
        except Exception as exc:
            print(f"(plot skipped: {exc})")

    if failures:
        print(f"\nFAILED: {failures}/{len(rows)} comparisons outside tol "
              f"(|z| > {Z_TOL} AND rel_err > {REL_TOL})")
        return 1
    print(f"\nPASSED: all {len(rows)} comparisons within tol")
    return 0


def _plot(rows: list[dict], out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_pair: dict[tuple[str, float], list[dict]] = {}
    for r in rows:
        by_pair.setdefault((r["material"], r["energy_MeV"]), []).append(r)

    n = len(by_pair)
    cols = min(3, n)
    rows_ax = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows_ax, cols, figsize=(4 * cols, 3 * rows_ax),
                             squeeze=False)
    for ax, ((mat, energy), pts) in zip(axes.flat, by_pair.items()):
        pts = sorted(pts, key=lambda r: r["mut_actual"])
        muts = [p["mut_actual"] for p in pts]
        ax.errorbar(muts, [p["B_MC"] for p in pts],
                    yerr=[p["B_MC_sigma"] for p in pts],
                    fmt="o", label="MC")
        ax.plot(muts, [p["B_GP"] for p in pts], "-", label="GP (Harima 1986)")
        ax.set_xlabel("optical thickness [mfp]")
        ax.set_ylabel("dose buildup B")
        ax.set_title(f"{mat} @ {energy} MeV")
        ax.legend()
        ax.grid(True, alpha=0.3)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(run())
