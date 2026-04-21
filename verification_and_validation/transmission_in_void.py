"""Verification: transmission fraction through a pure void.

With no attenuating medium the optical thickness is zero, so the uncollided
transmission fraction must be exactly 1.0 for any particle type, energy, or
distance:

    T = exp(-sum_i(Sigma_i * t_i)) = exp(0) = 1

This is checked via both `calculate_transmission` (returns the scalar
directly) and `calculate_flux` (which exposes `.transmission_fraction` on its
result). Also verifies that a multi-layer void stack still gives T = 1.
"""

import rad_point_kernel as rpk


TOL_ABS = 0.0  # void transmission must be exactly 1.0
THICKNESSES_CM = [0.0, 1.0, 10.0, 100.0, 1000.0, 10_000.0]

CASES = [
    ("neutron", 14.06e6),
    ("neutron", 2.45e6),
    ("neutron", 0.0253),
    ("photon", 1.0e6),
    ("photon", 6.62e5),
    ("photon", 1.0e4),
]

MULTILAYER_STACKS = [
    [50.0],
    [10.0, 20.0, 30.0],
    [1.0] * 10,
    [100.0, 200.0, 50.0, 75.0, 25.0],
]


def run() -> int:
    header = f"{'particle':>8}  {'E (eV)':>10}  {'t (cm)':>10}  {'T (trans)':>12}  {'T (flux)':>12}"
    print(header)
    print("-" * len(header))

    failures = 0
    total = 0

    for particle, energy in CASES:
        source = rpk.Source(particle, energy)
        for t in THICKNESSES_CM:
            layers = [rpk.Layer(thickness=t)]
            t_trans = rpk.calculate_transmission(layers, source)
            # calculate_flux needs r > 0 for inverse-square; use t+1 so we can still read transmission_fraction
            r = t if t > 0 else 1.0
            t_flux = rpk.calculate_flux(1e12, [rpk.Layer(thickness=r)], source).transmission_fraction

            ok_trans = abs(t_trans - 1.0) <= TOL_ABS
            ok_flux = abs(t_flux - 1.0) <= TOL_ABS
            status = "OK" if (ok_trans and ok_flux) else "FAIL"
            if not (ok_trans and ok_flux):
                failures += 1
            total += 1
            print(f"{particle:>8}  {energy:>10.3e}  {t:>10.1f}  {t_trans:>12.10f}  {t_flux:>12.10f}  {status}")

    print()
    print("multi-layer void stacks (neutron 14.06 MeV):")
    source = rpk.Source("neutron", 14.06e6)
    for stack in MULTILAYER_STACKS:
        layers = [rpk.Layer(thickness=t) for t in stack]
        t_trans = rpk.calculate_transmission(layers, source)
        ok = abs(t_trans - 1.0) <= TOL_ABS
        status = "OK" if ok else "FAIL"
        if not ok:
            failures += 1
        total += 1
        print(f"  thicknesses={stack}  T={t_trans:.15f}  {status}")

    print()
    if failures:
        print(f"FAILED: {failures}/{total} case(s)")
        return 1
    print(f"PASSED: all {total} cases have T == 1.0 exactly")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
