"""Weight-window generation for OpenMC using the point-kernel importance.

Users typically don't call this module directly — `compute_buildup` invokes it
under the hood when ``use_weight_windows=True`` (the default). See
``docs/guide/weight_windows.md`` for a user-facing overview.

This module is pure glue over ``rad_point_kernel_core``'s Rust mesh builder.
No physics lives here.
"""
from __future__ import annotations

import logging
from typing import Iterable

import numpy as np

from rad_point_kernel_core import (
    Quantity,
    build_weight_windows_maybe,
    importance_at,
)

from rad_point_kernel.buildup import _parse_quantity as _parse_quantity_str

_log = logging.getLogger("rad_point_kernel.weight_windows")


def _parse_quantity(q: str) -> tuple[Quantity, bool]:
    """Parse a quantity string into a Rust ``Quantity`` enum + coupled flag.

    Shares vocabulary with ``rad_point_kernel.buildup._parse_quantity`` so
    that string validation and the dose-geometry whitelist live in one place.
    """
    measure, geo, coupled = _parse_quantity_str(q)
    if measure == "flux":
        return Quantity.flux(), coupled
    return Quantity.dose(geo), coupled


def _geom_from(q: str) -> str:
    """Pull the dose geometry code out of a quantity string."""
    base = q.removesuffix("-coupled-photon")
    if not base.startswith("dose-"):
        raise ValueError(f"quantity '{q}' has no dose geometry")
    return base[5:]


def _select_driving_quantity(layers, source, quantities: Iterable[str]) -> str:
    """Pick the quantity with the steepest importance drop across the shield.

    See ``ww-plan.md §6a`` for why this works. Essentially a free call — the
    point-kernel evaluations are microseconds each.
    """
    qs = list(quantities)
    if len(qs) == 1:
        return qs[0]

    r_outer = sum(layer.thickness for layer in layers)
    if r_outer <= 0.0:
        return qs[0]

    def tau_total(q: str) -> float:
        q_enum, _ = _parse_quantity(q)
        try:
            i_in = importance_at(layers, source, q_enum, r_outer * 1e-6)
            i_out = importance_at(layers, source, q_enum, r_outer)
        except Exception as e:
            _log.warning("importance_at failed for quantity %r: %s", q, e)
            return 0.0
        return float(np.log(max(i_in, 1e-300) / max(i_out, 1e-300)))

    return max(qs, key=tau_total)


def _ww_from_plan(
    plan,
    *,
    particle_type: str,
    upper_bound_ratio: float,
):
    """Convert a ``WeightWindowPlan`` into an ``openmc.WeightWindows`` object."""
    import openmc  # lazy import — openmc is optional at package install time

    nr = len(plan.r_grid) - 1
    ne = len(plan.energy_bounds) - 1

    # plan.lower_bounds is shape (nE, nr) — transpose to (nr, nE) then reshape
    # to (nr, 1, 1, nE) for OpenMC's SphericalMesh storage order.
    lower = np.asarray(plan.lower_bounds, dtype=float).T  # (nr, nE)
    lower = lower.reshape(nr, 1, 1, ne)

    mesh = openmc.SphericalMesh(r_grid=list(plan.r_grid))
    return openmc.WeightWindows(
        mesh=mesh,
        energy_bounds=list(plan.energy_bounds),
        lower_ww_bounds=lower,
        upper_bound_ratio=upper_bound_ratio,
        particle_type=particle_type,
    )


def build_weight_windows(
    layers,
    source,
    *,
    quantities: Iterable[str],
    log_ratio_per_bin: float = 1.0,
    min_bin_width_cm: float = 0.5,
    upper_bound_ratio: float = 5.0,
) -> list:
    """Build openmc.WeightWindows for a spherical-shell simulation.

    Parameters
    ----------
    layers
        List of ``rad_point_kernel.Layer`` describing the spherical shells.
    source
        ``rad_point_kernel.Source`` — particle type and energy.
    quantities
        One or more tally-quantity strings, same vocabulary as
        ``compute_buildup``: ``"flux"``, ``"dose-AP"``, etc., optionally with
        ``"-coupled-photon"`` suffix.
    log_ratio_per_bin
        Target log drop in importance between adjacent radial bins. Default
        1.0 (factor of e per bin, adjacent-bin ratio ≈ 2.7). Smaller values
        (0.5 → ratio 1.6) give finer splitting; larger values (2.0 → ratio 7)
        are more aggressive and risk correlated-daughter CPU waste.
    min_bin_width_cm
        Drop mesh boundaries that would produce bins narrower than this.
    upper_bound_ratio
        Ratio of upper to lower WW bound. OpenMC's default is 5.0.

    Returns
    -------
    list[openmc.WeightWindows]
        Length 0 if the skip-gate fired (τ < 3 or n_bins ≤ 2).
        Length 1 for a non-coupled simulation (single primary particle).
        Length 2 for coupled n→γ: first element is the neutron WW, second is
        the photon WW driven by secondary-photon dose importance.
    """
    qs = list(quantities)
    if not qs:
        return []

    # Pick the steepest quantity for the primary WW (§6a)
    driving = _select_driving_quantity(layers, source, qs)
    q_primary, _ = _parse_quantity(driving)
    primary_particle = source.particle  # "neutron" or "photon"

    plan = build_weight_windows_maybe(
        layers, source, q_primary, primary_particle,
        log_ratio_per_bin, min_bin_width_cm,
    )
    if plan is None:
        _log.info(
            "Weight windows skipped: geometry too thin (τ < 3 or n_bins ≤ 2) "
            "for driving quantity %r.",
            driving,
        )
        return []

    _log.info(
        "Weight windows built for %r: %d radial bins, %d energy bins, τ ≈ %.2f.",
        driving, len(plan.r_grid) - 1, len(plan.energy_bounds) - 1,
        plan.tau_total,
    )

    ww_list = [_ww_from_plan(
        plan,
        particle_type=primary_particle,
        upper_bound_ratio=upper_bound_ratio,
    )]

    # Coupled n→γ: always add photon WW if any quantity is coupled,
    # regardless of which won the driving selection
    coupled_q = next((q for q in qs if q.endswith("-coupled-photon")), None)
    if coupled_q is not None and primary_particle == "neutron":
        q_secondary = Quantity.secondary_photon_dose(_geom_from(coupled_q))
        plan_g = build_weight_windows_maybe(
            layers, source, q_secondary, "photon",
            log_ratio_per_bin, min_bin_width_cm,
        )
        if plan_g is None:
            _log.info(
                "Photon WW skipped: τ < 3 or n_bins ≤ 2 for the secondary-gamma "
                "importance curve.",
            )
        else:
            _log.info(
                "Photon WW built: %d radial bins, %d energy bins, τ ≈ %.2f.",
                len(plan_g.r_grid) - 1, len(plan_g.energy_bounds) - 1,
                plan_g.tau_total,
            )
            ww_list.append(_ww_from_plan(
                plan_g,
                particle_type="photon",
                upper_bound_ratio=upper_bound_ratio,
            ))

    return ww_list
