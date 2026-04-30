"""Build-up factor computation and interpolation.

Quantity naming convention for compute_buildup:
    flux                        Flux of the source particle
    dose-AP                     Dose of the source particle at AP geometry
    flux-coupled-photon         Secondary photon flux (neutron source only)
    dose-AP-coupled-photon      Secondary photon dose at AP (neutron source only)
    dose-AP-total               Neutron + secondary-photon dose. Accepted as
                                input shorthand (expands to both halves) and
                                auto-added to the output. The buildup factor
                                uses the neutron PK as reference. Same rule
                                applies to PA, RLAT, LLAT, ROT, ISO.
"""

import math
import tempfile
import warnings

import numpy as np

from rad_point_kernel_core import (
    BuildupFit,
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


def _expand_total_requests(quantities, source):
    """Expand 'dose-{geo}-total' shorthand into the two halves it implies:
    'dose-{geo}' and 'dose-{geo}-coupled-photon'. The auto-synth on the
    output side fills the total key back in, so the user can request and
    read the same string. Order-preserving and dedup'd.
    """
    out = []
    seen = set()

    def add(q):
        if q not in seen:
            seen.add(q)
            out.append(q)

    for q in quantities:
        if q.startswith("dose-") and q.endswith("-total"):
            geo = q[len("dose-") : -len("-total")]
            if geo not in VALID_DOSE_GEOMETRIES:
                raise ValueError(
                    f"Invalid quantity '{q}': geometry must be one of "
                    f"{sorted(VALID_DOSE_GEOMETRIES)}"
                )
            if source.particle != "neutron":
                raise ValueError(
                    f"'{q}' requires a neutron source; got {source.particle!r}. "
                    f"Photon sources don't produce secondary photons in this tool."
                )
            add(f"dose-{geo}")
            add(f"dose-{geo}-coupled-photon")
        else:
            add(q)
    return out


def compute_buildup(
    geometries,
    source,
    quantities,
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
        quantities: Required. List of quantity strings. Examples: "flux",
            "dose-AP", "flux-coupled-photon", "dose-AP-coupled-photon".
            "dose-{geo}-total" is accepted as a shorthand for both halves
            (requires a neutron source); the result still carries all
            three keys. When both halves are present in the result a
            "dose-{geo}-total" entry is auto-synthesized (sum of the two
            doses; buildup factor uses the neutron PK as reference).
        particles_per_batch: Particles per batch (default 10,000).
        batches: Minimum number of batches before the trigger can stop the
            run early (default 10). Gives the tally enough statistics to
            evaluate the relative-error trigger reliably.
        max_batches: Safety cap on number of batches (default 100).
        trigger_rel_err: Target relative error on tallies (default 0.05).
        cross_sections: Path to cross_sections.xml or directory containing it.
            Restored to its previous value when this call returns.

    Returns:
        List of BuildupResult, one per geometry. MC flux and dose are returned
        per source particle. Apply an absolute strength via
        `result.scale(strength)` to land in the unit you want (Sv/hr for
        particles/sec, Sv/shot for particles/shot, etc.).
    """
    if isinstance(quantities, str):
        quantities = [quantities]
    quantities = list(quantities)
    if not quantities:
        raise ValueError("compute_buildup: 'quantities' must be non-empty")
    quantities = _expand_total_requests(quantities, source)

    geometries = list(geometries)
    if not geometries:
        raise ValueError("compute_buildup: 'geometries' must be non-empty")
    for i, layers in enumerate(geometries):
        if not layers:
            raise ValueError(
                f"compute_buildup: geometries[{i}] has no layers - need at "
                f"least one Layer to define an outer surface."
            )
        total_thickness = sum(layer.thickness for layer in layers)
        if total_thickness <= 0:
            raise ValueError(
                f"compute_buildup: geometries[{i}] has zero or negative total "
                f"thickness ({total_thickness} cm); detector surface would be "
                f"at the origin."
            )

    sentinel = object()
    prev_xs = sentinel
    if cross_sections is not None:
        import openmc
        from pathlib import Path
        path = Path(cross_sections).expanduser()
        if path.is_dir():
            path = path / "cross_sections.xml"
        prev_xs = openmc.config.get("cross_sections", sentinel)
        openmc.config["cross_sections"] = str(path)

    # Parse quantities
    parsed = [(q, _parse_quantity(q)) for q in quantities]

    # Determine transport mode
    needs_coupled = any(p[2] for _, p in parsed)

    results = []
    try:
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
            result.synthesize_dose_totals()
            results.append(result)
    finally:
        if cross_sections is not None:
            import openmc
            if prev_xs is sentinel:
                openmc.config.pop("cross_sections", None)
            else:
                openmc.config["cross_sections"] = prev_xs

    return results


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
            warnings.warn(
                f"MC tally '{q_name}' returned 0 for geometry with "
                f"layer thicknesses {thicknesses} cm - not enough particles "
                f"penetrated the shield. Increase particles_per_batch or "
                f"max_batches, or use variance reduction. Build-up factor "
                f"from this result will be unreliable.",
                UserWarning,
                stacklevel=2,
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
