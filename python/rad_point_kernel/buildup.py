"""Build-up factor computation and interpolation.

Quantity naming convention for compute_buildup:
    flux                        Flux of the source particle
    dose-AP                     Dose of the source particle at AP geometry
    flux-coupled-photon         Secondary photon flux (neutron source only)
    dose-AP-coupled-photon      Secondary photon dose at AP (neutron source only)
    dose-AP-total               Auto-synthesized when both 'dose-AP' and
                                'dose-AP-coupled-photon' are requested. Equals
                                the neutron + secondary-photon dose; the
                                buildup factor uses the neutron PK as
                                reference. Same rule applies to PA, RLAT,
                                LLAT, ROT, ISO.
"""

import math
import sys
import tempfile

import numpy as np

from rad_point_kernel_core import (
    BuildupResult,
    InterpolationResult,
    calculate_dose,
    calculate_flux,
)

VALID_DOSE_GEOMETRIES = {"AP", "PA", "RLAT", "LLAT", "ROT", "ISO"}


def _parse_quantity(q):
    """Parse a quantity string. Source carries the particle type.

    Returns (measure, geometry, coupled).
    """
    coupled = q.endswith("-coupled-photon")
    base = q[: -len("-coupled-photon")] if coupled else q

    if base == "flux":
        return ("flux", None, coupled)

    if base.startswith("dose-"):
        geo = base[5:]
        if geo not in VALID_DOSE_GEOMETRIES:
            raise ValueError(
                f"Invalid quantity '{q}': geometry must be one of {VALID_DOSE_GEOMETRIES}"
            )
        return ("dose", geo, coupled)

    raise ValueError(f"Invalid quantity '{q}': must be 'flux' or 'dose-{{geometry}}'")


class BuildupTable:
    """GP-based interpolation table for build-up factors.

    Axis values are layer thicknesses in cm (matching Layer.thickness).

    Example::

        table = BuildupTable(
            points=[{"poly": 5, "conc": 40}, ...],
            results=mc_results,
        )
        result = table.interpolate(poly=7, conc=60)
    """

    def __init__(self, points, results):
        if len(points) != len(results):
            raise ValueError(
                f"points ({len(points)}) and results ({len(results)}) must have same length"
            )
        if len(points) < 2:
            raise ValueError("Need at least 2 data points for interpolation")

        self._axis_names = sorted(points[0].keys())
        for p in points:
            if sorted(p.keys()) != self._axis_names:
                raise ValueError(
                    f"All points must have the same axes. "
                    f"Expected {self._axis_names}, got {sorted(p.keys())}"
                )

        self._points = points
        self._results = results
        self._n_points = len(points)
        self._n_dims = len(self._axis_names)

        self._x = np.array(
            [[p[ax] for ax in self._axis_names] for p in points], dtype=float
        )

        self._axis_ranges = {
            ax: (float(self._x[:, i].min()), float(self._x[:, i].max()))
            for i, ax in enumerate(self._axis_names)
        }

        available = set()
        r0 = results[0]
        for q in r0.buildup:
            if all(q in r.buildup and r.buildup[q] is not None for r in results):
                available.add(q)

        self._kept_indices = {}
        for q in sorted(available):
            kept = [i for i, r in enumerate(results) if r.buildup[q] != 0.0]
            dropped = [points[i] for i, r in enumerate(results) if r.buildup[q] == 0.0]
            if dropped:
                print(
                    f"WARNING: BuildupTable dropping {len(dropped)} point(s) "
                    f"with build-up = 0 for quantity '{q}' at {dropped}. "
                    f"Usually the MC tally returned 0 (insufficient statistics); "
                    f"the GP would be pulled toward 0 and extrapolation unreliable. "
                    f"Re-run MC with more particles or use variance reduction.",
                    file=sys.stderr,
                )
            if len(kept) < 2:
                print(
                    f"WARNING: Quantity '{q}' has only {len(kept)} non-zero "
                    f"point(s); need >= 2 to fit a GP. Removing from table.",
                    file=sys.stderr,
                )
                continue
            self._kept_indices[q] = kept

        self._available = set(self._kept_indices.keys())
        if not self._available:
            raise ValueError(
                "BuildupTable has no quantity with >= 2 non-zero build-up "
                "points. All MC tallies may have returned 0 - re-run with "
                "more particles or use variance reduction."
            )

        sorted_qs = sorted(self._available)
        self._default_quantity = sorted_qs[0]
        self._gps = {}

    @property
    def axis_names(self):
        return list(self._axis_names)

    @property
    def available_quantities(self):
        return set(self._available)

    @property
    def axis_ranges(self):
        return dict(self._axis_ranges)

    def _get_gp(self, quantity):
        if quantity in self._gps:
            return self._gps[quantity]

        from inference.gp import GpRegressor
        from inference.gp.covariance import SquaredExponential

        keep = self._kept_indices[quantity]
        x_kept = self._x[keep]
        kept_results = [self._results[i] for i in keep]

        y = np.array([r.buildup[quantity] for r in kept_results], dtype=float)
        y_err = np.array(
            [
                r.mc_std_dev.get(quantity, 0.01) / r.pk[quantity]
                if r.pk.get(quantity, 0) > 0
                else 0.01
                for r in kept_results
            ],
            dtype=float,
        )
        y_err = np.maximum(y_err, 1e-10)

        hyperpar_bounds = [(-10, 10)]
        for i in range(self._n_dims):
            coords = np.sort(np.unique(x_kept[:, i]))
            if len(coords) > 1:
                min_spacing = np.min(np.diff(coords))
            else:
                min_spacing = 1.0
            min_log_ls = np.log(max(min_spacing * 1.5, 1e-6))
            hyperpar_bounds.append((min_log_ls, 15))

        kernel = SquaredExponential(hyperpar_bounds=hyperpar_bounds)
        x_input = x_kept if self._n_dims > 1 else x_kept.ravel()
        gp = GpRegressor(x_input, y, y_err=y_err, kernel=kernel)
        self._gps[quantity] = gp
        return gp

    def interpolate(self, quantity=None, warn=True, **kwargs):
        if quantity is None:
            quantity = self._default_quantity
        if quantity is None:
            raise ValueError("No quantities available in this table")
        if quantity not in self._available:
            raise ValueError(
                f"Quantity '{quantity}' not available. Available: {sorted(self._available)}"
            )

        provided = sorted(kwargs.keys())
        if provided != self._axis_names:
            raise ValueError(f"Expected axes {self._axis_names}, got {provided}")

        x_query = np.array([[kwargs[ax] for ax in self._axis_names]], dtype=float)
        if self._n_dims == 1:
            x_query = x_query.ravel()

        is_extrapolated = False
        extrapolated_axes = {}
        for i, ax in enumerate(self._axis_names):
            val = kwargs[ax]
            lo, hi = self._axis_ranges[ax]
            if val < lo or val > hi:
                is_extrapolated = True
                extrapolated_axes[ax] = (val, lo, hi)

        if is_extrapolated and warn:
            parts = [
                f"{ax}={val} is outside simulated range [{lo}, {hi}]"
                for ax, (val, lo, hi) in extrapolated_axes.items()
            ]
            print(
                f"WARNING: Extrapolating: {'; '.join(parts)}. "
                "Build-up factor may be inaccurate.",
                file=sys.stderr,
            )

        gp = self._get_gp(quantity)
        mu, sigma = gp(x_query)

        return InterpolationResult(
            value=float(mu[0]),
            sigma=float(sigma[0]),
            is_extrapolated=is_extrapolated,
            extrapolated_axes=extrapolated_axes,
        )

    def __repr__(self):
        return (
            f"BuildupTable(axes={self._axis_names}, "
            f"n_points={self._n_points}, "
            f"quantities={sorted(self._available)})"
        )


def compute_buildup(
    geometries,
    source,
    quantities=("flux",),
    particles_per_batch=10_000,
    batches=10,
    max_batches=100,
    trigger_rel_err=0.05,
    cross_sections=None,
    use_weight_windows=True,
    max_history_splits=100,
):
    """Run OpenMC MC on a list of geometries and compute build-up factors.

    Args:
        geometries: List of layer lists, each defining a shield geometry
            (Layer thicknesses in cm).
        source: Source object with particle type and energy (eV).
        quantities: List of quantity strings. Examples: "flux", "dose-AP",
            "flux-coupled-photon", "dose-AP-coupled-photon". When both
            "dose-{geo}" and "dose-{geo}-coupled-photon" are requested for
            the same geometry, a synthetic "dose-{geo}-total" quantity is
            added to each result (sum of the two doses; buildup factor uses
            the neutron PK as reference).
        particles_per_batch: Particles per batch (default 10,000).
        batches: Minimum number of batches before the trigger can stop the
            run early (default 10). Gives the tally enough statistics to
            evaluate the relative-error trigger reliably.
        max_batches: Safety cap on number of batches (default 100).
        trigger_rel_err: Target relative error on tallies (default 0.05).
        cross_sections: Path to cross_sections.xml or directory containing it.

    Returns:
        List of BuildupResult, one per geometry. MC flux and dose are returned
        per source particle. Apply an absolute strength via
        `result.scale(strength)` to land in the unit you want (Sv/hr for
        particles/sec, Sv/shot for particles/shot, etc.).
    """
    if isinstance(quantities, str):
        quantities = [quantities]
    quantities = list(quantities)

    if cross_sections is not None:
        import openmc
        from pathlib import Path
        path = Path(cross_sections).expanduser()
        if path.is_dir():
            path = path / "cross_sections.xml"
        openmc.config["cross_sections"] = str(path)

    # Parse quantities
    parsed = [(q, _parse_quantity(q)) for q in quantities]

    # Determine transport mode
    needs_coupled = any(p[2] for _, p in parsed)

    results = []
    for layers in geometries:
        result = BuildupResult()

        # Optical thickness from PK
        pk_flux = calculate_flux(layers=layers, source=source)
        result.optical_thickness = pk_flux.optical_thickness

        # Run MC
        mc_data = _run_mc(
            layers, source, needs_coupled,
            parsed, particles_per_batch, batches, max_batches, trigger_rel_err,
            use_weight_windows=use_weight_windows,
            max_history_splits=max_history_splits,
        )
        _populate_result(result, layers, source, parsed, mc_data)
        _synthesize_dose_totals(result, parsed)
        results.append(result)

    return results


def _synthesize_dose_totals(result, parsed):
    """Add 'dose-{geo}-total' for each geometry where both the neutron and
    coupled-photon dose were tallied. Total dose = neutron + secondary photon;
    the buildup factor uses the neutron PK as reference (secondary photons
    have no analytic PK)."""
    quantities = {q_name for q_name, _ in parsed}

    for q_name, (measure, geo, coupled) in parsed:
        if measure != "dose" or coupled:
            continue
        coupled_name = f"dose-{geo}-coupled-photon"
        if coupled_name not in quantities:
            continue

        total_name = f"dose-{geo}-total"
        n_mc = result.mc.get(q_name, 0.0)
        p_mc = result.mc.get(coupled_name, 0.0)
        n_std = result.mc_std_dev.get(q_name, 0.0)
        p_std = result.mc_std_dev.get(coupled_name, 0.0)

        result.mc[total_name] = n_mc + p_mc
        result.mc_std_dev[total_name] = math.sqrt(n_std ** 2 + p_std ** 2)

        pk_n = result.pk.get(q_name)
        if pk_n is not None:
            result.pk[total_name] = pk_n
            if pk_n > 0:
                result.buildup[total_name] = (n_mc + p_mc) / pk_n


def _run_mc(layers, source, coupled, quantities, particles_per_batch, batches, max_batches, trigger_rel_err,
            *, use_weight_windows=True, max_history_splits=100):
    """Run a single MC simulation. Returns dict: quantity_name -> (mean, std_dev)."""
    import openmc

    source_particle = source.particle
    source_energy = source.energy  # float or list of (energy, weight) tuples

    # Materials
    omc_materials = []
    layer_mats = []
    for layer in layers:
        mat = layer.material
        if mat is not None:
            omc_mat = openmc.Material(
                components=mat.nuclide_mass_fractions(),
                density=mat.density,
                density_units="g/cm3",
                percent_type="wo",
            )
            omc_materials.append(omc_mat)
            layer_mats.append(omc_mat)
        else:
            layer_mats.append(None)

    # Geometry
    r = 0.0
    radii = []
    for layer in layers:
        r += layer.thickness
        radii.append(r)
    outer_radius = radii[-1]
    surface_area = 4.0 * math.pi * outer_radius ** 2

    spheres = []
    for i, radius in enumerate(radii):
        btype = "vacuum" if i == len(radii) - 1 else "transmission"
        spheres.append(openmc.Sphere(r=radius, boundary_type=btype))

    cells = []
    for i in range(len(layers)):
        if i == 0:
            region = -spheres[0]
        else:
            region = +spheres[i - 1] & -spheres[i]
        cells.append(openmc.Cell(region=region, fill=layer_mats[i]))

    omc_geometry = openmc.Geometry(openmc.Universe(cells=cells))

    # Source
    omc_source = openmc.IndependentSource()
    omc_source.space = openmc.stats.Point((0, 0, 0))
    omc_source.angle = openmc.stats.Isotropic()
    omc_source.particle = source_particle

    if isinstance(source_energy, (list, tuple)):
        energies = [e for e, _ in source_energy]
        weights = [w for _, w in source_energy]
        omc_source.energy = openmc.stats.Discrete(energies, weights)
    else:
        omc_source.energy = openmc.stats.Discrete([source_energy], [1.0])

    # Settings
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = particles_per_batch
    settings.batches = batches
    settings.source = [omc_source]
    if coupled:
        settings.photon_transport = True
    settings.trigger_active = True
    settings.trigger_max_batches = max_batches
    settings.trigger_batch_interval = 5

    # Weight-window generation — driven by the user's tally quantities.
    if use_weight_windows:
        from .weight_windows import build_weight_windows
        ww_list = build_weight_windows(
            layers, source,
            quantities=[q for q, _ in quantities],
        )
        if ww_list:
            settings.weight_windows = ww_list
            settings.max_history_splits = max_history_splits

    # Tallies
    surface_filter = openmc.SurfaceFilter([spheres[-1]])
    tallies_obj = openmc.Tallies()
    tally_map = {}

    for q_name, (measure, geo, is_coupled) in quantities:
        tally_particle = "photon" if is_coupled else source_particle

        particle_filter = openmc.ParticleFilter([tally_particle])

        if measure == "flux":
            tally = openmc.Tally(name=q_name)
            tally.filters = [surface_filter, particle_filter]
            tally.scores = ["current"]
            tally.triggers.append(openmc.Trigger("rel_err", trigger_rel_err))
        else:
            energy, coeffs = openmc.data.dose_coefficients(tally_particle, geo)
            tally = openmc.Tally(name=q_name)
            tally.filters = [
                surface_filter, particle_filter,
                openmc.EnergyFunctionFilter(energy, coeffs),
            ]
            tally.scores = ["current"]
            tally.triggers.append(openmc.Trigger("rel_err", trigger_rel_err))

        tallies_obj.append(tally)
        tally_map[q_name] = tally

    # Run
    model = openmc.Model(
        geometry=omc_geometry,
        materials=openmc.Materials(omc_materials),
        settings=settings,
        tallies=tallies_obj,
    )

    mc_data = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        model.run(cwd=tmpdir, output=False, apply_tally_results=True)

        for q_name, tally in tally_map.items():
            mean = tally.mean.flatten()[0]
            std = tally.std_dev.flatten()[0]

            _, (measure, geo, _coupled) = next(
                (qn, p) for qn, p in quantities if qn == q_name
            )

            if measure == "flux":
                mc_data[q_name] = (mean / surface_area, std / surface_area)
            else:
                # Dose tally with EnergyFunctionFilter is in pSv*cm^2 per
                # source particle integrated over the surface; convert to Sv
                # per source particle and divide out the surface area to get
                # the per-cm^2 dose. No time factor: scale(strength) sets the
                # time unit at the call site.
                mc_data[q_name] = (
                    mean / surface_area * 1e-12,
                    std / surface_area * 1e-12,
                )

    return mc_data


def _populate_result(result, layers, source, quantities, mc_data):
    """Fill a BuildupResult with MC data, PK references, and buildup factors."""
    for q_name, (measure, geo, coupled) in quantities:
        mc_val, mc_std = mc_data[q_name]
        result.mc[q_name] = mc_val
        result.mc_std_dev[q_name] = mc_std

        if mc_val == 0.0:
            thicknesses = [layer.thickness for layer in layers]
            print(
                f"WARNING: MC tally '{q_name}' returned 0 for geometry with "
                f"layer thicknesses {thicknesses} cm - not enough particles "
                f"penetrated the shield. Increase particles_per_batch or "
                f"max_batches, or use variance reduction. Build-up factor "
                f"from this result will be unreliable.",
                file=sys.stderr,
            )

        # PK reference - only for primary particle quantities (not coupled secondary)
        if not coupled:
            if measure == "flux":
                pk = calculate_flux(layers=layers, source=source)
                pk_val = pk.flux
            else:
                pk = calculate_dose(layers=layers, source=source, geometry=geo)
                pk_val = pk.dose

            result.pk[q_name] = pk_val
            if pk_val and pk_val > 0:
                result.buildup[q_name] = mc_val / pk_val
