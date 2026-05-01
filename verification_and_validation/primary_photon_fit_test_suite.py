"""V&V harness for the primary-photon dose fit (Shin-Ishii form).

Sibling of `secondary_photon_fit_test_suite.py` and
`neutron_fit_test_suite.py`. Validates the Shin-Ishii fit shipped in
`rad_point_kernel_core` for the primary-photon dose response when the
source is a photon (not neutron):

    B(tau) = A * exp(-alpha1 * tau) + (1 - A) * exp(-alpha2 * tau)

The `_photon` cache entries (e.g. `polyethylene_photon`,
`concrete_photon`) come from `compute_buildup` runs with a 1 MeV photon
source; the `dose-AP` field is then the primary-photon dose, not the
secondary-photon dose.

Same machinery applies to flux: regenerate the cache with
`quantities=["flux-photon"]` and load the `flux-photon` field instead
of `dose-AP`.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


CACHE = Path(__file__).resolve().parent / "secondary_photon_fit_cache.json"

PHOTON_MATERIALS = [
    "polyethylene_photon", "borated_poly_photon", "lih_photon", "d2o_photon",
    "b4c_photon", "concrete_photon", "aluminum_photon",
    "concrete_steel10_photon", "steel_photon", "ss316_photon",
    "tungsten_photon", "lead_photon", "beryllium_photon",
]


@dataclass
class MaterialData:
    name: str
    tau_train: np.ndarray
    b_train: np.ndarray
    sigma_train: np.ndarray
    tau_holdout: float
    b_holdout: float
    sigma_holdout: float


def _load() -> dict[str, MaterialData]:
    raw = json.loads(CACHE.read_text())
    out = {}
    for name in PHOTON_MATERIALS:
        if name not in raw or f"{name}_holdout" not in raw:
            continue
        rs = raw[name]
        h = raw[f"{name}_holdout"][0]
        tau = np.array([r["optical_thickness"] for r in rs], dtype=float)
        pk = np.array([r["pk"]["dose-AP"] for r in rs], dtype=float)
        mc = np.array([r["mc"]["dose-AP"] for r in rs], dtype=float)
        sd = np.array([r["mc_std_dev"]["dose-AP"] for r in rs], dtype=float)
        h_mc = float(h["mc"]["dose-AP"])
        h_pk = float(h["pk"]["dose-AP"])
        if not (mc > 0).all() or h_mc <= 0:
            continue
        b = mc / pk
        sigma = np.maximum(sd / pk, 1e-30)
        out[name] = MaterialData(
            name=name,
            tau_train=tau,
            b_train=b,
            sigma_train=sigma,
            tau_holdout=float(h["optical_thickness"]),
            b_holdout=h_mc / h_pk,
            sigma_holdout=max(float(h["mc_std_dev"]["dose-AP"]) / h_pk, 1e-30),
        )
    return out


PRIMARY_PHOTON_DATA = _load()


def shin_ishii_eval(tau, A, alpha1, alpha2):
    t = np.asarray(tau, dtype=float)
    return A * np.exp(-alpha1 * t) + (1 - A) * np.exp(-alpha2 * t)


def fit_shin_ishii(tau, b, sigma):
    def loss(p):
        A, a1, a2 = p
        if not np.isfinite(p).all():
            return 1e30
        pred = shin_ishii_eval(tau, A, a1, a2)
        return float(np.sum(((pred - b) / sigma) ** 2))

    starts = [
        [0.5, 0.1, 0.5], [0.3, 0.05, 1.0], [0.7, 0.5, 2.0],
        [0.5, 1.0, 0.1], [0.5, 2.0, 0.5], [0.9, 0.05, 1.5],
    ]
    best, best_l = None, np.inf
    for s in starts:
        r = minimize(loss, s, method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 30000})
        if r.fun < best_l:
            best_l, best = r.fun, r.x
    return best


def run_suite(name: str = "Shin-Ishii (3p)") -> dict:
    per_mat = {}
    for material, m in PRIMARY_PHOTON_DATA.items():
        params = fit_shin_ishii(m.tau_train, m.b_train, m.sigma_train)
        pred_train = shin_ishii_eval(m.tau_train, *params)
        z_fit = float(np.max(np.abs((pred_train - m.b_train) / m.sigma_train)))
        pred_h = float(shin_ishii_eval(np.array([m.tau_holdout]), *params)[0])
        z_h = abs(pred_h - m.b_holdout) / max(m.sigma_holdout, 1e-30)
        per_mat[material] = {
            "z_holdout": z_h, "z_fit_max": z_fit,
            "b_pred_holdout": pred_h, "b_mc_holdout": m.b_holdout,
        }
    z_h = np.array([v["z_holdout"] for v in per_mat.values()])
    z_f = np.array([v["z_fit_max"] for v in per_mat.values()])
    return {
        "name": name, "per_material": per_mat,
        "n_materials": len(per_mat),
        "worst_z_holdout": float(z_h.max()),
        "mean_z_holdout": float(z_h.mean()),
        "median_z_holdout": float(np.median(z_h)),
        "worst_z_fit": float(z_f.max()),
    }


if __name__ == "__main__":
    print(f"Primary-photon fit V&V over {len(PRIMARY_PHOTON_DATA)} photon-source materials")
    print()
    report = run_suite()
    print(f"{'material':>26} {'tau_h':>6} {'B_mc':>9} {'B_pred':>9} {'z_holdout':>10} {'z_fit_max':>10}")
    print("-" * 80)
    for mat in sorted(report["per_material"], key=lambda k: report["per_material"][k]["z_holdout"]):
        v = report["per_material"][mat]
        m = PRIMARY_PHOTON_DATA[mat]
        short = mat[:-len("_photon")] if mat.endswith("_photon") else mat
        print(f"{short:>26} {m.tau_holdout:>6.2f} {v['b_mc_holdout']:>9.4f} "
              f"{v['b_pred_holdout']:>9.4f} {v['z_holdout']:>10.2f} {v['z_fit_max']:>10.2f}")
    print("-" * 80)
    print(f"Worst holdout z: {report['worst_z_holdout']:.2f}  "
          f"Mean: {report['mean_z_holdout']:.2f}  "
          f"Median: {report['median_z_holdout']:.2f}")
    print(f"Worst training z: {report['worst_z_fit']:.2f}")
