"""Verification: Beer-Lambert single-slab attenuation.

For an uncollided beam traversing a homogeneous slab, the transmission is

    T(t) = exp(-Sigma * t)

where Sigma is the macroscopic total cross section (monoenergetic). This
verification does NOT assume any particular Sigma value — it verifies the
*exponential law* itself, which is independent of cross-section data:

  1. Self-similarity across thickness:
        T(k * t_ref) == T(t_ref) ** k    for any k > 0
     Equivalently: -log(T) must be linear in t (slope = Sigma >= 0).

  2. Multiplicativity of stacked layers of the same material:
        T(t1 + t2) == T(t1) * T(t2)
     This is checked both as (a) a single layer of thickness t1+t2 vs.
     (b) two stacked layers of thicknesses t1 and t2.

  3. Zero-thickness limit: T(0) == 1 for a material layer (not a void).

  4. Strict monotonicity: T decreases as t increases.

If any of these fail, the kernel's attenuation math is wrong regardless of
which cross-section library is loaded.
"""

import math

import rad_point_kernel as pkc


TOL_REL = 1e-10  # double-precision exp + one cross-section lookup round-trip

# (label, Material, Source)
def _mat(comp, density):
    return pkc.Material(composition=comp, density=density)

CASES = [
    ("iron / n 14.06 MeV",   _mat({"Fe": 1.0}, 7.874),              pkc.Source("neutron", 14.06e6)),
    ("iron / n 2.45 MeV",    _mat({"Fe": 1.0}, 7.874),              pkc.Source("neutron", 2.45e6)),
    ("iron / photon 1 MeV",  _mat({"Fe": 1.0}, 7.874),              pkc.Source("photon", 1.0e6)),
    ("water / n 14.06 MeV",  _mat({"H": 2.0, "O": 1.0}, 1.0),       pkc.Source("neutron", 14.06e6)),
    ("water / photon 1 MeV", _mat({"H": 2.0, "O": 1.0}, 1.0),       pkc.Source("photon", 1.0e6)),
    ("lead / photon 662 keV",_mat({"Pb": 1.0}, 11.34),              pkc.Source("photon", 6.62e5)),
    ("lead / photon 1 MeV",  _mat({"Pb": 1.0}, 11.34),              pkc.Source("photon", 1.0e6)),
]

T_REF = 5.0  # cm — reference thickness for the exponential-scaling check
K_VALUES = [0.5, 1.0, 2.0, 3.0, 5.0, 10.0]  # multiples of t_ref
SPLIT_PAIRS = [(2.0, 8.0), (4.0, 4.0), (1.0, 9.0), (3.5, 6.5)]


def _T(material, source, thickness):
    return pkc.calculate_transmission([pkc.Layer(thickness=thickness, material=material)], source)


def _T_stacked(material, source, t1, t2):
    return pkc.calculate_transmission(
        [pkc.Layer(thickness=t1, material=material), pkc.Layer(thickness=t2, material=material)],
        source,
    )


def _relerr(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def run() -> int:
    failures = 0
    total = 0

    print("1. Exponential scaling: T(k*t_ref) vs T(t_ref)**k")
    print(f"   t_ref = {T_REF} cm")
    print(f"   {'case':<24}  {'k':>5}  {'T(k*t_ref)':>14}  {'T(t_ref)^k':>14}  {'rel err':>10}")
    print("   " + "-" * 78)
    for label, material, source in CASES:
        T_ref = _T(material, source, T_REF)
        for k in K_VALUES:
            T_k = _T(material, source, k * T_REF)
            T_expected = T_ref ** k
            err = _relerr(T_k, T_expected)
            ok = err <= TOL_REL
            status = "OK" if ok else "FAIL"
            if not ok:
                failures += 1
            total += 1
            print(f"   {label:<24}  {k:>5.1f}  {T_k:>14.6e}  {T_expected:>14.6e}  {err:>10.2e}  {status}")

    print()
    print("2. Multiplicativity: T(t1+t2) vs T(t1)*T(t2), single-layer vs stacked")
    print(f"   {'case':<24}  {'t1':>5}  {'t2':>5}  {'T_sum':>14}  {'T1*T2':>14}  {'T_stack':>14}  {'max err':>10}")
    print("   " + "-" * 98)
    for label, material, source in CASES:
        for t1, t2 in SPLIT_PAIRS:
            T_sum = _T(material, source, t1 + t2)
            T1 = _T(material, source, t1)
            T2 = _T(material, source, t2)
            T_stack = _T_stacked(material, source, t1, t2)
            e1 = _relerr(T_sum, T1 * T2)
            e2 = _relerr(T_stack, T1 * T2)
            err = max(e1, e2)
            ok = err <= TOL_REL
            status = "OK" if ok else "FAIL"
            if not ok:
                failures += 1
            total += 1
            print(f"   {label:<24}  {t1:>5.1f}  {t2:>5.1f}  {T_sum:>14.6e}  {T1*T2:>14.6e}  {T_stack:>14.6e}  {err:>10.2e}  {status}")

    print()
    print("3. Zero-thickness limit: T(0) == 1 exactly for a material layer")
    print(f"   {'case':<24}  {'T(0)':>14}  {'status':>6}")
    print("   " + "-" * 50)
    for label, material, source in CASES:
        T0 = _T(material, source, 0.0)
        ok = T0 == 1.0
        status = "OK" if ok else "FAIL"
        if not ok:
            failures += 1
        total += 1
        print(f"   {label:<24}  {T0:>14.12f}  {status:>6}")

    print()
    print("4. Strict monotonicity: T(t) decreases with t")
    print(f"   {'case':<24}  {'status':>6}")
    print("   " + "-" * 36)
    thicknesses = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0]
    for label, material, source in CASES:
        Ts = [_T(material, source, t) for t in thicknesses]
        ok = all(Ts[i] > Ts[i + 1] for i in range(len(Ts) - 1))
        status = "OK" if ok else "FAIL"
        if not ok:
            failures += 1
        total += 1
        print(f"   {label:<24}  {status:>6}")

    print()
    if failures:
        print(f"FAILED: {failures}/{total} case(s)")
        return 1
    print(f"PASSED: all {total} Beer-Lambert checks within tol {TOL_REL:.0e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
