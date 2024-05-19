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

    Builds the quantity key from source particle, quantity type, and geometry.
    """
    if buildup is None:
        return None
    if isinstance(buildup, BuildupModel):
        return buildup
    if isinstance(buildup, InterpolationResult):
        return BuildupModel.constant(buildup.value)
    if isinstance(buildup, BuildupResult):
        # Build key: "flux-neutron", "dose-AP-photon", etc.
        particle = source.particle
        if geometry:
            key = f"dose-{geometry}-{particle}"
        else:
            key = f"flux-{particle}"
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

