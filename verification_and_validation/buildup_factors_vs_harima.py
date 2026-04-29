"""Verification: MC-derived photon dose buildup factors vs ANS-6.4.3 / Harima G-P fits.

ANSI/ANS-6.4.3-1991 publishes Geometric-Progression-form coefficients
(a, b, c, d, X_k) for ~25 elements plus standard concretes and water at
photon energies 0.015 - 15 MeV. The full B(mu*t) curve is

    K(mu*t) = c*(mu*t)^a + d*[tanh(mu*t/X_k - 2) - tanh(-2)] / [1 - tanh(-2)]
    B(mu*t) = 1 + (b - 1) * (K^mu*t - 1) / (K - 1)        if K != 1
    B(mu*t) = 1 + (b - 1) * mu*t                          if K == 1

This script runs `compute_buildup` for a small set of (material, energy)
cases, extracts B(mu*t) from the Monte Carlo dose, and checks that it
agrees with the published Harima curve within 3 sigma of the MC error.

Note: per the TODO this is *verification*, not validation - the Harima
tables are themselves computational, so we are checking code-to-code
consistency, not agreement with measurement.

Coefficients are placeholders read from public secondary sources; before
treating this script as production verification, cross-check each row
against the primary tables in ANS-6.4.3-1991 or JAERI-M 86-058
(Harima, Sakamoto, Tanaka, Kawai, 1986).

Requires OpenMC and a cross_sections.xml available via the
OPENMC_CROSS_SECTIONS environment variable.
"""
import math

import rad_point_kernel as rpk


def harima_gp_form(mu_t: float, a: float, b: float, c: float, d: float, X_k: float) -> float:
    """ANS-6.4.3 Geometric-Progression form."""
    K = (
        c * mu_t ** a
        + d * (math.tanh(mu_t / X_k - 2.0) - math.tanh(-2.0))
        / (1.0 - math.tanh(-2.0))
    )
    if abs(K - 1.0) < 1e-6:
        return 1.0 + (b - 1.0) * mu_t
    return 1.0 + (b - 1.0) * (K ** mu_t - 1.0) / (K - 1.0)


# Materials matching the compositions used in the Harima tables.
WATER = rpk.Material(composition={"H2O": 1.0}, density=1.0)
IRON = rpk.Material(composition={"Fe": 1.0}, density=7.874)
LEAD = rpk.Material(composition={"Pb": 1.0}, density=11.34)
CONCRETE = rpk.Material(  # ANS-6.4.3 ordinary concrete composition
    composition={
        "H": 0.0056, "C": 0.0149, "O": 0.4983, "Na": 0.0171,
        "Mg": 0.0048, "Al": 0.0456, "Si": 0.3158, "K": 0.0192,
        "Ca": 0.0834, "Fe": 0.0124,
    },
    density=2.30, fraction="mass",
)


# (material, label, energy_eV, harima coefficients, thicknesses_cm to test)
# Coefficients below are placeholders - cross-check against ANS-6.4.3-1991.
CASES = [
    (WATER, "water", 1.0e6,
     {"a": 0.046, "b": 1.881, "c": 0.518, "d": -0.045, "X_k": 14.43},
     [10.0, 25.0, 50.0]),
    (IRON, "iron", 1.0e6,
     {"a": 0.115, "b": 1.483, "c": 0.506, "d": -0.051, "X_k": 14.13},
     [3.0, 6.0, 10.0]),
    (LEAD, "lead", 1.0e6,
     {"a": 0.116, "b": 1.379, "c": 0.570, "d": -0.054, "X_k": 14.06},
     [1.0, 3.0, 6.0]),
    (CONCRETE, "concrete", 1.0e6,
     {"a": 0.082, "b": 1.610, "c": 0.515, "d": -0.054, "X_k": 13.21},
     [5.0, 15.0, 30.0]),
]

# Pass/fail tolerance: how many sigma of MC noise we'll accept between
# B_MC and B_Harima before flagging a row.
TOL_SIGMA = 3.0


def run() -> int:
    header = (
        f"{'material':>10}  {'E(MeV)':>7}  {'t(cm)':>7}  {'mu*t':>6}  "
        f"{'B_MC':>9}  {'sigma_B':>9}  {'B_Harima':>9}  {'rel_err':>9}  {'n_sig':>6}"
    )
    print(header)
    print("-" * len(header))

    failures = 0
    total = 0
    for material, label, energy, coeffs, thicks in CASES:
        source = rpk.Source(particle="photon", energy=energy)
        geometries = [
            [rpk.Layer(thickness=1000.0), rpk.Layer(thickness=t, material=material)]
            for t in thicks
        ]
        try:
            results = rpk.compute_buildup(
                geometries=geometries,
                source=source,
                quantities=["dose-AP"],
                particles_per_batch=10_000,
                max_batches=200,
                trigger_rel_err=0.02,
            )
        except Exception as exc:
            print(f"{label:>10}: MC failed - {exc}")
            failures += len(thicks)
            total += len(thicks)
            continue

        for t, r in zip(thicks, results):
            mu_t = r.optical_thickness
            pk_dose = r.pk["dose-AP"]
            mc_dose = r.mc["dose-AP"]
            sd_dose = r.mc_std_dev["dose-AP"]
            b_mc = mc_dose / pk_dose if pk_dose > 0 else float("nan")
            sigma_b = sd_dose / pk_dose if pk_dose > 0 else float("nan")
            b_ref = harima_gp_form(mu_t, **coeffs)
            rel_err = (b_mc - b_ref) / b_ref
            n_sig = abs(b_mc - b_ref) / max(sigma_b, 1e-12)
            status = "OK" if n_sig <= TOL_SIGMA else "FAIL"
            if status == "FAIL":
                failures += 1
            total += 1
            print(
                f"{label:>10}  {energy / 1e6:>7.3f}  {t:>7.1f}  {mu_t:>6.2f}  "
                f"{b_mc:>9.4f}  {sigma_b:>9.4f}  {b_ref:>9.4f}  "
                f"{rel_err * 100:>+7.2f}%  {n_sig:>6.2f}  {status}"
            )

    print()
    if failures:
        print(f"FAILED: {failures}/{total} case(s) exceeded {TOL_SIGMA:.0f}-sigma tolerance")
        return 1
    print(f"PASSED: all {total} case(s) agree with Harima G-P fits within {TOL_SIGMA:.0f} sigma")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
