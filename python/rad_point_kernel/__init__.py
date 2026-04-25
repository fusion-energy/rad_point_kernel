"""Point-kernel neutron and photon transmission/dose calculator."""

from rad_point_kernel_core import (
    BuildupModel,
    CalcResult,
    Layer,
    Material,
    SecondaryGammaResult,
    Source,
)
from rad_point_kernel_core import (
    calculate_transmission as _calculate_transmission,
    calculate_flux as _calculate_flux,
    calculate_dose as _calculate_dose,
    calculate_secondary_photon_dose_rate as _calculate_secondary_photon_dose_rate,
)
from rad_point_kernel.buildup import (
    BuildupResult,
    BuildupTable,
    InterpolationResult,
    compute_buildup,
)


def _resolve_buildup(buildup, source, quantity_type, geometry=None):
    """Convert BuildupResult or InterpolationResult to BuildupModel.

    Builds the quantity key from quantity type and geometry. The particle
    type is carried by the source, so it is not part of the key — this
    matches the convention in buildup.py (keys are "flux", "dose-AP", ...).
    """
    if buildup is None:
        return None
    if isinstance(buildup, BuildupModel):
        return buildup
    if isinstance(buildup, InterpolationResult):
        return BuildupModel.constant(buildup.value)
    if isinstance(buildup, BuildupResult):
        key = f"dose-{geometry}" if geometry else "flux"
        b = buildup.buildup.get(key)
        if b is None:
            available = sorted(buildup.buildup.keys())
            raise ValueError(
                f"BuildupResult doesn't contain '{key}'. "
                f"Available: {available}"
            )
        return BuildupModel.constant(b)
    raise TypeError(
        f"buildup must be BuildupModel, BuildupResult, or InterpolationResult, "
        f"got {type(buildup).__name__}"
    )


def calculate_transmission(layers, source):
    """Calculate transmission fraction through layers.

    Args:
        layers: List of Layer objects.
        source: Source object with particle type and energy.

    Returns:
        Transmission fraction (float in [0, 1]).
    """
    return _calculate_transmission(layers, source)


def calculate_flux(source_strength, layers, source, buildup=None):
    """Calculate uncollided flux at the outer surface.

    Returns S / (4piR^2) * exp(-sum(Sigma_r * t)) * B, in particles/cm^2/s.

    Args:
        source_strength: Source strength in particles/sec.
        layers: List of Layer objects.
        source: Source object with particle type and energy.
        buildup: BuildupModel, BuildupResult, or InterpolationResult (optional).

    Returns:
        CalcResult with uncollided_flux, transmission_fraction, etc.
    """
    buildup = _resolve_buildup(buildup, source, "flux")
    return _calculate_flux(source_strength, layers, source, buildup)


def calculate_dose(source_strength, layers, source, geometry, buildup=None):
    """Calculate dose rate at the outer surface.

    Args:
        source_strength: Source strength in particles/sec.
        layers: List of Layer objects.
        source: Source object with particle type and energy.
        geometry: Dose geometry ("AP", "PA", "RLAT", "LLAT", "ROT", "ISO").
        buildup: BuildupModel, BuildupResult, or InterpolationResult (optional).

    Returns:
        CalcResult with dose_rate in Sv/hr.
    """
    buildup = _resolve_buildup(buildup, source, "dose", geometry)
    return _calculate_dose(source_strength, layers, source, geometry, buildup)


# TODO: Consider removing from public API — coupled MC is more accurate.
def calculate_secondary_photon_dose_rate(
    source_strength, layers, source, geometry, neutron_buildup=None,
):
    """Calculate secondary photon dose rate (analytical, no MC).

    Args:
        source_strength: Source strength in particles/sec.
        layers: List of Layer objects.
        source: Source object (must be a monoenergetic neutron source).
        geometry: Dose geometry ("AP", "PA", "RLAT", "LLAT", "ROT", "ISO").
        neutron_buildup: Optional BuildupModel, BuildupResult, or InterpolationResult.
    """
    if source.particle != "neutron":
        raise ValueError(
            f"calculate_secondary_photon_dose_rate requires a neutron source, "
            f"got {source.particle!r}"
        )
    neutron_buildup = _resolve_buildup(neutron_buildup, source, "dose", geometry)
    return _calculate_secondary_photon_dose_rate(
        source_strength, layers, source.energy, geometry, neutron_buildup,
    )


__all__ = [
    "BuildupModel",
    "BuildupResult",
    "BuildupTable",
    "CalcResult",
    "InterpolationResult",
    "Layer",
    "Material",
    "SecondaryGammaResult",
    "Source",
    "calculate_dose",
    "calculate_flux",
    "calculate_secondary_photon_dose_rate",
    "calculate_transmission",
    "compute_buildup",
]
