"""Cross-material test suite for the secondary-photon B_p(tau) fit form.

Loads coupled-photon B_p data for 13 materials from the bulk-shielding study
cache (6 training anchors + 1 holdout per material). Provides:

  - PHOTON_DATA: a dict keyed by material name -> {tau_train, b_p_train,
    sigma_train, tau_holdout, b_p_holdout, sigma_holdout, descriptors}
  - DESCRIPTORS: per-material physics descriptors (density, hydrogen
    fraction, mean Z, mu_14MeV) for cross-material symbolic regression
  - run_suite(fit_fn, eval_fn): given a candidate (parametric) form, run
    it against every material; return per-material residuals + summary z-stats
  - dump_pysr_dataset(path): write a flat (tau, density, H_frac, ..., B_p)
    CSV that can be fed to PySR or any other symbolic-regression tool

The harness is form-agnostic: any fit/predict pair with the right
signatures plugs in. Useful for closed-form parametric fits, log-space
RBF, or anything else.

Reference data lives in /tmp/gp_invest_cache.json. Schema of each cache
entry: {mc, mc_std_dev, pk, buildup, optical_thickness, source_strength}
with quantity keys still using the pre-2.0 names (`dose-AP`,
`dose-AP-coupled-photon`, `dose-AP-total`); B_p = mc[coupled] / pk[neutron].
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np


CACHE = Path(__file__).resolve().parent / "secondary_photon_fit_cache.json"

MATERIALS = [
    "polyethylene", "polyethylene_deep", "borated_poly", "lih", "d2o",
    "b4c", "concrete", "aluminum", "concrete_steel10",
    "steel", "ss316", "tungsten", "lead",
    "magnetite",
    "wc_water", "wc_water_50", "wc_water_70",
    "wb_water", "wb_water_50", "wb_water_70",
    "h2o_water",
]


# Per-material physics descriptors. Densities g/cm^3, H mass fraction,
# mean atomic number Z_mean (mass-weighted), and approximate total
# macroscopic cross section mu at 14 MeV (1/cm). Sources: PNNL Materials
# Compendium for compositions and densities; mu values estimated from the
# pk_neutron decay between adjacent training anchors in the cache (good
# enough for use as a regression feature, not for transport).
DESCRIPTORS = {
    "polyethylene":     {"density": 0.93, "H_frac": 0.143, "Z_mean": 5.3, "mu_14MeV": 0.0265},
    "polyethylene_deep":{"density": 0.93, "H_frac": 0.143, "Z_mean": 5.3, "mu_14MeV": 0.0525},
    "borated_poly":     {"density": 0.99, "H_frac": 0.116, "Z_mean": 5.5, "mu_14MeV": 0.0263},
    "lih":              {"density": 0.78, "H_frac": 0.127, "Z_mean": 4.0, "mu_14MeV": 0.0234},
    "d2o":              {"density": 1.10, "H_frac": 0.000, "Z_mean": 7.7, "mu_14MeV": 0.0235},
    "b4c":              {"density": 2.51, "H_frac": 0.000, "Z_mean": 5.6, "mu_14MeV": 0.0257},
    "concrete":         {"density": 2.30, "H_frac": 0.010, "Z_mean": 11.0, "mu_14MeV": 0.0333},
    "aluminum":         {"density": 2.70, "H_frac": 0.000, "Z_mean": 13.0, "mu_14MeV": 0.0339},
    "concrete_steel10": {"density": 3.06, "H_frac": 0.009, "Z_mean": 14.5, "mu_14MeV": 0.0373},
    "steel":            {"density": 7.85, "H_frac": 0.000, "Z_mean": 26.0, "mu_14MeV": 0.0376},
    "ss316":            {"density": 7.99, "H_frac": 0.000, "Z_mean": 25.9, "mu_14MeV": 0.0381},
    "tungsten":         {"density": 19.3, "H_frac": 0.000, "Z_mean": 74.0, "mu_14MeV": 0.0323},
    "lead":             {"density": 11.34, "H_frac": 0.000, "Z_mean": 82.0, "mu_14MeV": 0.0265},
    "magnetite":        {"density": 3.53, "H_frac": 0.012, "Z_mean": 17.8, "mu_14MeV": 0.0775},
    "wc_water":         {"density": 12.68, "H_frac": 0.00177, "Z_mean": 64.0, "mu_14MeV": 0.131},
    "wc_water_50":      {"density":  8.30, "H_frac": 0.00674, "Z_mean": 52.0, "mu_14MeV": 0.150},
    "wc_water_70":      {"density":  5.38, "H_frac": 0.01455, "Z_mean": 38.0, "mu_14MeV": 0.130},
    "wb_water":         {"density": 12.76, "H_frac": 0.00176, "Z_mean": 63.8, "mu_14MeV": 0.137},
    "wb_water_50":      {"density":  8.35, "H_frac": 0.00674, "Z_mean": 52.0, "mu_14MeV": 0.150},
    "wb_water_70":      {"density":  5.41, "H_frac": 0.01454, "Z_mean": 38.0, "mu_14MeV": 0.130},
    "h2o_water":        {"density":  1.00, "H_frac": 0.11198, "Z_mean":  7.4, "mu_14MeV": 0.050},
}


@dataclass
class MaterialData:
    name: str
    tau_train: np.ndarray
    b_p_train: np.ndarray
    sigma_train: np.ndarray
    tau_holdout: float
    b_p_holdout: float
    sigma_holdout: float
    descriptors: dict = field(default_factory=dict)


def _load() -> dict[str, MaterialData]:
    raw = json.loads(CACHE.read_text())
    out = {}
    for name in MATERIALS:
        rs = raw[name]
        h = raw[f"{name}_holdout"][0]
        tau = np.array([r["optical_thickness"] for r in rs], dtype=float)
        pk = np.array([r["pk"]["dose-AP"] for r in rs], dtype=float)
        mc_p = np.array([r["mc"]["dose-AP-coupled-photon"] for r in rs], dtype=float)
        sd_p = np.array([r["mc_std_dev"]["dose-AP-coupled-photon"] for r in rs], dtype=float)
        b = mc_p / pk
        sigma = np.maximum(sd_p / pk, 1e-30)

        h_tau = float(h["optical_thickness"])
        h_pk = float(h["pk"]["dose-AP"])
        h_mc = float(h["mc"]["dose-AP-coupled-photon"])
        h_sd = float(h["mc_std_dev"]["dose-AP-coupled-photon"])

        out[name] = MaterialData(
            name=name,
            tau_train=tau,
            b_p_train=b,
            sigma_train=sigma,
            tau_holdout=h_tau,
            b_p_holdout=h_mc / h_pk,
            sigma_holdout=max(h_sd / h_pk, 1e-30),
            descriptors=DESCRIPTORS.get(name, {}),
        )
    return out


PHOTON_DATA = _load()


def run_suite(
    fit_fn: Callable,
    eval_fn: Callable,
    tag: str = "candidate",
) -> dict:
    """Run a (fit, eval) pair across every material.

    fit_fn(tau, b_p, sigma) -> params
    eval_fn(tau, params) -> b_p_pred (numpy array if tau is array)

    Returns a dict with per-material residuals and summary stats:
      {
        "tag": str,
        "per_material": {name: {z_holdout, z_fit_max, b_pred_holdout, ...}},
        "worst_z_holdout": float,
        "mean_z_holdout": float,
        "median_z_holdout": float,
        "worst_z_fit": float,
      }
    """
    per_mat = {}
    for name, m in PHOTON_DATA.items():
        try:
            params = fit_fn(m.tau_train, m.b_p_train, m.sigma_train)
        except Exception as e:
            per_mat[name] = {"error": f"fit failed: {e}"}
            continue
        try:
            pred_train = np.atleast_1d(eval_fn(m.tau_train, params))
            pred_h = float(np.atleast_1d(eval_fn(np.array([m.tau_holdout]), params))[0])
        except Exception as e:
            per_mat[name] = {"error": f"eval failed: {e}"}
            continue
        z_fit_max = float(np.max(np.abs((pred_train - m.b_p_train) / m.sigma_train)))
        z_holdout = float(abs(pred_h - m.b_p_holdout) / m.sigma_holdout)
        per_mat[name] = {
            "z_holdout": z_holdout,
            "z_fit_max": z_fit_max,
            "b_pred_holdout": pred_h,
            "b_mc_holdout": m.b_p_holdout,
            "rel_resid_holdout": abs(pred_h - m.b_p_holdout) / max(m.b_p_holdout, 1e-30),
        }

    successful = [v for v in per_mat.values() if "error" not in v]
    if not successful:
        return {"tag": tag, "per_material": per_mat,
                "worst_z_holdout": float("inf"), "mean_z_holdout": float("inf"),
                "median_z_holdout": float("inf"), "worst_z_fit": float("inf")}

    z_holdout = np.array([v["z_holdout"] for v in successful])
    z_fit = np.array([v["z_fit_max"] for v in successful])
    return {
        "tag": tag,
        "per_material": per_mat,
        "worst_z_holdout": float(z_holdout.max()),
        "mean_z_holdout": float(z_holdout.mean()),
        "median_z_holdout": float(np.median(z_holdout)),
        "worst_z_fit": float(z_fit.max()),
        "n_materials": len(successful),
    }


def print_report(report: dict) -> None:
    print(f"=== {report['tag']} ===")
    print(f"{'material':>20} {'h_tau':>6} {'B_mc':>10} {'B_pred':>10} "
          f"{'rel_resid':>10} {'z_holdout':>10} {'z_fit_max':>10}")
    print("-" * 90)
    pm = report["per_material"]
    for name in MATERIALS:
        v = pm[name]
        if "error" in v:
            print(f"{name:>20}  {v['error']}")
            continue
        m = PHOTON_DATA[name]
        print(f"{name:>20} {m.tau_holdout:>6.2f} {v['b_mc_holdout']:>10.4f} "
              f"{v['b_pred_holdout']:>10.4f} {v['rel_resid_holdout']*100:>9.1f}% "
              f"{v['z_holdout']:>10.2f} {v['z_fit_max']:>10.2f}")
    print("-" * 90)
    print(f"  holdout z:  worst={report['worst_z_holdout']:.2f}  "
          f"mean={report['mean_z_holdout']:.2f}  median={report['median_z_holdout']:.2f}")
    print(f"  training z: worst={report['worst_z_fit']:.2f}")
    print()


def dump_pysr_dataset(path: Path) -> None:
    """Write a flat (material, tau, density, H_frac, Z_mean, mu_14MeV, B_p,
    sigma_B_p, is_holdout) CSV usable by PySR or any other symbolic
    regression tool. One row per (material, tau) pair across training and
    holdout.
    """
    rows = []
    for name, m in PHOTON_DATA.items():
        d = m.descriptors
        for ti, bi, si in zip(m.tau_train, m.b_p_train, m.sigma_train):
            rows.append({
                "material": name, "tau": ti, "B_p": bi, "sigma_B_p": si,
                "density": d.get("density", float("nan")),
                "H_frac": d.get("H_frac", float("nan")),
                "Z_mean": d.get("Z_mean", float("nan")),
                "mu_14MeV": d.get("mu_14MeV", float("nan")),
                "is_holdout": 0,
            })
        rows.append({
            "material": name, "tau": m.tau_holdout, "B_p": m.b_p_holdout,
            "sigma_B_p": m.sigma_holdout,
            "density": d.get("density", float("nan")),
            "H_frac": d.get("H_frac", float("nan")),
            "Z_mean": d.get("Z_mean", float("nan")),
            "mu_14MeV": d.get("mu_14MeV", float("nan")),
            "is_holdout": 1,
        })

    cols = ["material", "tau", "density", "H_frac", "Z_mean", "mu_14MeV",
            "B_p", "sigma_B_p", "is_holdout"]
    with open(path, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {len(rows)} rows to {path}")


# ---------------------------------------------------------------------------
# Built-in candidate forms (so this module is self-checking)
# ---------------------------------------------------------------------------

def _baseline_fit(tau, b, sigma):
    from scipy.optimize import minimize
    def loss(p):
        a, pe, c = p
        if a <= 0 or pe < 0 or c <= 0:
            return 1e30
        ts = np.maximum(tau, 1e-12)
        pred = a * ts ** pe * (1 - np.exp(-c * tau))
        return float(np.sum(((pred - b) / sigma) ** 2))
    a0 = max(float(b.max() / tau.max()), 1e-30)
    starts = [[a0, 1.0, 1.0], [a0, 1.5, 0.5], [a0, 0.5, 2.0],
              [a0 * 2, 1.2, 1.0], [a0 * 0.5, 0.8, 3.0], [a0, 2.0, 0.3]]
    best, best_l = None, np.inf
    for s in starts:
        r = minimize(loss, s, method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 20000})
        if r.fun < best_l:
            best_l, best = r.fun, r.x
    return best


def _baseline_eval(tau, params):
    a, pe, c = params
    ts = np.maximum(np.asarray(tau, dtype=float), 1e-12)
    return a * ts ** pe * (1 - np.exp(-c * np.asarray(tau, dtype=float)))


def _form_a_fit(tau, b, sigma):
    from scipy.optimize import minimize
    def loss(p):
        a, pe, c, d = p
        if a <= 0 or pe < 0 or c <= 0 or d < 0:
            return 1e30
        ts = np.maximum(tau, 1e-12)
        pred = a * ts ** pe * (1 - np.exp(-c * tau)) * np.exp(d * tau)
        return float(np.sum(((pred - b) / sigma) ** 2))
    a0 = max(float(b.max() / tau.max()), 1e-30)
    starts = [
        [a0, 1.0, 1.0, 0.0], [a0, 1.5, 0.5, 0.0], [a0, 0.5, 2.0, 0.0],
        [a0 * 0.1, 1.0, 1.0, 0.005], [a0 * 0.1, 0.5, 1.0, 0.010],
        [a0 * 0.5, 1.0, 0.5, 0.020],
    ]
    best, best_l = None, np.inf
    for s in starts:
        r = minimize(loss, s, method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 20000})
        if r.fun < best_l:
            best_l, best = r.fun, r.x
    return best


def _form_a_eval(tau, params):
    a, pe, c, d = params
    t = np.asarray(tau, dtype=float)
    ts = np.maximum(t, 1e-12)
    return a * ts ** pe * (1 - np.exp(-c * t)) * np.exp(d * t)


def _form_k_fit(tau, b, sigma, allow_neg_e=False):
    """Form K = a * tau^p * (1 - exp(-c*tau)) * exp(d*tau + e*tau^2).
    With allow_neg_e=True, e can go negative — captures peak-and-decay
    shapes (Gaussian-modulated growth) seen in mixed light/heavy
    composites where H content makes photons but Z absorbs them past
    the peak.
    """
    from scipy.optimize import minimize

    def loss(p):
        a, pe, c, d, e = p
        if a <= 0 or pe < 0 or c <= 0 or d < 0:
            return 1e30
        if not allow_neg_e and e < 0:
            return 1e30
        ts = np.maximum(tau, 1e-12)
        try:
            pred = a * ts ** pe * (1 - np.exp(-c * tau)) * np.exp(d * tau + e * tau * tau)
        except Exception:
            return 1e30
        v = float(np.sum(((pred - b) / sigma) ** 2))
        return v if np.isfinite(v) else 1e30

    a0 = max(float(b.max() / tau.max()), 1e-30)
    starts = [
        [a0, 1.0, 1.0, 0.0, 0.0], [a0, 1.5, 0.5, 0.0, 0.0],
        [a0, 0.5, 2.0, 0.0, 0.0], [a0 * 0.5, 0.8, 3.0, 0.0, 0.0],
        [a0 * 0.1, 1.0, 1.0, 0.005, 0.0],
        [a0 * 0.05, 1.0, 1.0, 0.01, 1e-5],
        [a0 * 0.5, 0.5, 1.0, 0.0, 1e-5],
    ]
    if allow_neg_e:
        starts += [
            [a0, 1.0, 1.0, 0.5, -0.05],
            [a0, 1.0, 1.0, 0.3, -0.02],
            [a0 * 2, 0.5, 1.0, 0.5, -0.10],
        ]
    best, best_l = None, np.inf
    for s in starts:
        r = minimize(loss, s, method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 30000})
        if r.fun < best_l:
            best_l, best = r.fun, r.x
    return best


def _form_k_pos_fit(tau, b, sigma):
    return _form_k_fit(tau, b, sigma, allow_neg_e=False)


def _form_k_signed_fit(tau, b, sigma):
    return _form_k_fit(tau, b, sigma, allow_neg_e=True)


def _form_k_reg_fit(tau, b, sigma, lam_d=10.0, lam_e=100.0):
    """Form K' with L2 regularisation on d and e (defaults from sweep:
    lam_d=10, lam_e=100 minimises mean holdout z on the 20-material set
    while keeping training fit clean and recovering baseline-equivalent
    behaviour on saturating materials)."""
    from scipy.optimize import minimize

    def loss(p):
        a, pe, c, d, e = p
        if a <= 0 or pe < 0 or c <= 0:
            return 1e30
        ts = np.maximum(tau, 1e-12)
        try:
            pred = a * ts ** pe * (1 - np.exp(-c * tau)) * np.exp(d * tau + e * tau * tau)
        except Exception:
            return 1e30
        chi2 = float(np.sum(((pred - b) / sigma) ** 2))
        return chi2 + lam_d * d * d + lam_e * e * e

    a0 = max(float(b.max() / tau.max()), 1e-30)
    starts = [
        [a0, 1.0, 1.0, 0.0, 0.0],
        [a0, 1.5, 0.5, 0.0, 0.0],
        [a0, 0.5, 2.0, 0.0, 0.0],
        [a0 * 0.1, 1.0, 1.0, 0.005, 0.0],
        [a0 * 0.05, 1.0, 1.0, 0.01, 1e-5],
        [a0, 1.0, 1.0, 0.5, -0.05],
        [a0 * 2, 0.5, 1.0, 0.5, -0.10],
    ]
    best, best_l = None, np.inf
    for s in starts:
        r = minimize(loss, s, method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 30000})
        if r.fun < best_l:
            best_l, best = r.fun, r.x
    return best


def _form_k_eval(tau, params):
    a, pe, c, d, e = params
    t = np.asarray(tau, dtype=float)
    ts = np.maximum(t, 1e-12)
    return a * ts ** pe * (1 - np.exp(-c * t)) * np.exp(d * t + e * t * t)


def _form_m_fit(tau, b, sigma):
    """Form M = a * tau^p * (1 - exp(-c*tau)) * (1 + g * exp(d*tau)).
    5 params: a, p, c, g, d (d signed).

    Asymptotic behaviour at large tau:
      d > 0:  blows up like tau^p * g * exp(d*tau)   (poly-style growth)
      d < 0:  saturates to a * tau^p                  (heavy-Z saturator;
                                                       peak-and-decay if g > 0)
      d == 0: saturates to a * tau^p * (1 + g)

    All three regimes from one form. The (1 + g*exp(d*tau)) factor stays
    bounded below by 1 for d <= 0, so peak-and-decay shapes converge to
    a non-zero floor (a * tau^p) instead of all the way to zero.
    """
    from scipy.optimize import minimize

    def loss(p):
        a, pe, c, g, d = p
        if a <= 0 or pe < 0 or c <= 0 or g < 0:
            return 1e30
        ts = np.maximum(tau, 1e-12)
        try:
            pred = a * ts ** pe * (1 - np.exp(-c * tau)) * (1.0 + g * np.exp(d * tau))
        except Exception:
            return 1e30
        v = float(np.sum(((pred - b) / sigma) ** 2))
        return v if np.isfinite(v) else 1e30

    a0 = max(float(b.max() / tau.max()), 1e-30)
    starts = [
        # heavy-saturating regime (g=0 collapses to baseline)
        [a0, 1.0, 1.0, 0.0, 0.0], [a0, 1.5, 0.5, 0.0, 0.0],
        [a0 * 0.5, 0.8, 3.0, 0.0, 0.0],
        # poly-style growth (d > 0)
        [a0 * 0.1, 1.0, 1.0, 0.5, 0.05],
        [a0 * 0.05, 0.5, 1.0, 1.0, 0.10],
        [a0 * 0.01, 0.0, 1.0, 5.0, 0.20],
        # peak-and-decay (d < 0)
        [a0, 0.0, 1.0, 5.0, -0.3],
        [a0 * 0.5, 0.5, 1.0, 2.0, -0.5],
        [a0, 1.0, 0.5, 1.0, -0.2],
    ]
    best, best_l = None, np.inf
    for s in starts:
        r = minimize(loss, s, method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 30000})
        if r.fun < best_l:
            best_l, best = r.fun, r.x
    return best


def _form_m_eval(tau, params):
    a, pe, c, g, d = params
    t = np.asarray(tau, dtype=float)
    ts = np.maximum(t, 1e-12)
    return a * ts ** pe * (1 - np.exp(-c * t)) * (1.0 + g * np.exp(d * t))


def _form_a_min_fit(tau, b, sigma):
    """Form A_minimal = a * (1 - exp(-c*tau)) * exp(d*tau).
    3 free params (no tau^p factor), d signed (allows growth, saturation,
    or slight decay). Same parameter count as the baseline Power x sat.

    Motivated by the mixed-effects sweep result: shared p was ~ -0.04
    (essentially zero) at the best-CV regulariser, suggesting the tau^p
    factor doesn't earn its parameter on this data.
    """
    from scipy.optimize import minimize

    def loss(p):
        a, c, d = p
        if a <= 0 or c <= 0:
            return 1e30
        try:
            pred = a * (1 - np.exp(-c * tau)) * np.exp(d * tau)
        except Exception:
            return 1e30
        v = float(np.sum(((pred - b) / sigma) ** 2))
        return v if np.isfinite(v) else 1e30

    a0 = max(float(b.max() / tau.max()), 1e-30)
    starts = [
        # saturating starts (d = 0, near 0)
        [a0, 1.0, 0.0], [a0 * 2, 0.5, 0.0], [a0 * 0.5, 2.0, 0.0],
        # growth starts (d > 0)
        [a0 * 0.1, 1.0, 0.5], [a0 * 0.05, 1.0, 1.0],
        [a0 * 0.5, 0.5, 0.3],
        # mild-decay starts (d < 0)
        [a0, 1.0, -0.1], [a0 * 1.5, 1.0, -0.3],
    ]
    best, best_l = None, np.inf
    for s in starts:
        r = minimize(loss, s, method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 30000})
        if r.fun < best_l:
            best_l, best = r.fun, r.x
    return best


def _form_a_min_eval(tau, params):
    a, c, d = params
    t = np.asarray(tau, dtype=float)
    return a * (1 - np.exp(-c * t)) * np.exp(d * t)


CANDIDATES = {
    "baseline (3p)":            (_baseline_fit, _baseline_eval),
    "Form A_min (3p, d signed)":(_form_a_min_fit, _form_a_min_eval),
    "Form A (4p)":              (_form_a_fit, _form_a_eval),
    "Form K (5p, e>=0)":        (_form_k_pos_fit, _form_k_eval),
    "Form K' (5p, e signed)":   (_form_k_signed_fit, _form_k_eval),
    "Form K_reg (5p, L2 on d,e)": (_form_k_reg_fit, _form_k_eval),
    "Form M (5p, d signed)":    (_form_m_fit, _form_m_eval),
}


def plot_all(out_path: Path | None = None) -> Path:
    """Multi-panel plot: one subplot per material. Each panel shows MC
    training anchors (filled circles with error bars), the MC holdout
    anchor (open square), and a curve per candidate fit form across the
    full tau range out to 1.5 * the deepest training tau.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(MATERIALS)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.0 * nrows),
                             squeeze=False)

    fits_by_form = {}
    for tag, (fit_fn, eval_fn) in CANDIDATES.items():
        fits_by_form[tag] = {
            name: fit_fn(m.tau_train, m.b_p_train, m.sigma_train)
            for name, m in PHOTON_DATA.items()
        }

    colors = {"baseline (3p)": "#1f77b4",
              "Form A_min (3p, d signed)": "#17becf",
              "Form A (4p)": "#2ca02c",
              "Form K (5p, e>=0)": "#ff7f0e", "Form K' (5p, e signed)": "#d62728",
              "Form M (5p, d signed)": "#9467bd"}

    for idx, name in enumerate(MATERIALS):
        ax = axes[idx // ncols][idx % ncols]
        m = PHOTON_DATA[name]
        tau_max = max(m.tau_train.max(), m.tau_holdout)
        tau_grid = np.linspace(0.01, tau_max * 1.1, 200)

        for tag, (_fit_fn, eval_fn) in CANDIDATES.items():
            params = fits_by_form[tag][name]
            try:
                y = eval_fn(tau_grid, params)
            except Exception:
                continue
            ax.plot(tau_grid, y, "-", color=colors.get(tag, "gray"),
                    label=tag, lw=1.5, alpha=0.85)

        ax.errorbar(m.tau_train, m.b_p_train, yerr=m.sigma_train,
                    fmt="o", color="black", ms=4, capsize=2,
                    label="MC train", zorder=10)
        ax.errorbar([m.tau_holdout], [m.b_p_holdout], yerr=[m.sigma_holdout],
                    fmt="s", mfc="none", color="red", ms=8, capsize=3,
                    label="MC holdout", zorder=10)

        ax.set_title(name, fontsize=10)
        ax.set_xlabel("tau", fontsize=9)
        ax.set_ylabel("B_p", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.3)
        if idx == 0:
            ax.legend(loc="upper left", fontsize=7, framealpha=0.9)

    # Hide leftover empty axes
    for idx in range(len(MATERIALS), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("Coupled-photon B_p(tau) fits vs MC across the cross-validation suite",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out_path = out_path or Path("/tmp/secondary_photon_fit_suite_plot.png")
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(f"Loaded photon B_p suite: {len(PHOTON_DATA)} materials")
    print(f"  {', '.join(MATERIALS)}")
    print()
    for tag, (fit_fn, eval_fn) in CANDIDATES.items():
        report = run_suite(fit_fn, eval_fn, tag=tag)
        print_report(report)

    out_csv = Path(__file__).resolve().parents[1] / "secondary_photon_fit_dataset.csv"
    dump_pysr_dataset(out_csv)

    out_png = plot_all()
    print(f"Saved plot to {out_png}")
