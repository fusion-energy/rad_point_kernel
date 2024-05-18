"""Build-up factor computation and interpolation.

Quantity naming convention for compute_buildup:
    flux                        Flux of the source particle
    dose-AP                     Dose of the source particle at AP geometry
    flux-coupled-photon         Secondary photon flux (neutron source only)
    dose-AP-coupled-photon      Secondary photon dose at AP (neutron source only)
"""

from dataclasses import dataclass, field
import math
import sys
import tempfile

import numpy as np

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


@dataclass
class BuildupResult:
    """Result of a single MC build-up factor computation."""

    mc: dict = field(default_factory=dict)
    mc_std_dev: dict = field(default_factory=dict)
    pk: dict = field(default_factory=dict)
    buildup: dict = field(default_factory=dict)
    optical_thickness: float = 0.0

    def to_dict(self):
        return {
            "mc": self.mc,
            "mc_std_dev": self.mc_std_dev,
            "pk": self.pk,
            "buildup": self.buildup,
            "optical_thickness": self.optical_thickness,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            mc=d.get("mc", {}),
            mc_std_dev=d.get("mc_std_dev", {}),
            pk=d.get("pk", {}),
            buildup=d.get("buildup", {}),
            optical_thickness=d.get("optical_thickness", 0.0),
        )

    @staticmethod
    def save(results, path):
        import json
        from pathlib import Path
        Path(path).write_text(json.dumps([r.to_dict() for r in results], indent=2))

    @staticmethod
    def load(path):
        import json
        from pathlib import Path
        return [BuildupResult.from_dict(d) for d in json.loads(Path(path).read_text())]


@dataclass
class InterpolationResult:
    """Result of a BuildupTable interpolation query."""

    value: float
    sigma: float
    is_extrapolated: bool = False
