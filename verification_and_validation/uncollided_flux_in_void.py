"""Verification: uncollided flux from an isotropic point source in a pure void.

Reference: for an isotropic point source of strength S (particles/s) at the
origin and no attenuating medium, the uncollided scalar flux at distance r is

    phi(r) = S / (4 * pi * r^2)

This script computes phi(r) with `calculate_transmission` is not needed here -
we use `calculate_flux` through a single void Layer of thickness r, for a
logarithmic sweep of r. The point-kernel result must match the analytic
formula to machine precision for any source energy and particle type (void
implies zero cross section, so the transmission factor is exactly 1).
"""

import math

import rad_point_kernel as rpk


SOURCE_STRENGTH = 1e12  # particles/s
TOL_REL = 1e-12         # void case should be exact up to floating-point noise
DISTANCES_CM = [1.0, 10.0, 50.0, 100.0, 250.0, 500.0, 1000.0, 5000.0]

CASES = [
    ("neutron", 14.06e6),
    ("neutron", 2.45e6),
    ("photon", 1.0e6),
    ("photon", 6.62e5),
]


def analytic(s: float, r: float) -> float:
    return s / (4.0 * math.pi * r * r)


def run() -> int:
    header = f"{'particle':>8}  {'E (eV)':>10}  {'r (cm)':>8}  {'PK':>14}  {'analytic':>14}  {'rel err':>10}"
    print(header)
    print("-" * len(header))

    failures = 0
    for particle, energy in CASES:
        source = rpk.Source(particle, energy)
        for r in DISTANCES_CM:
            result = rpk.calculate_flux(SOURCE_STRENGTH, [rpk.Layer(thickness=r)], source)
            pk = result.flux
            ref = analytic(SOURCE_STRENGTH, r)
            rel = abs(pk - ref) / ref
            status = "OK" if rel <= TOL_REL else "FAIL"
            if rel > TOL_REL:
                failures += 1
            print(f"{particle:>8}  {energy:>10.3e}  {r:>8.1f}  {pk:>14.6e}  {ref:>14.6e}  {rel:>10.2e}  {status}")

            assert result.transmission_fraction == 1.0, (
                f"void transmission must be exactly 1.0, got {result.transmission_fraction}"
            )

    print()
    if failures:
        print(f"FAILED: {failures} case(s) exceeded tol {TOL_REL:.0e}")
        return 1
    print(f"PASSED: all {len(CASES) * len(DISTANCES_CM)} cases within tol {TOL_REL:.0e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
