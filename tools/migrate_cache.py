"""Migrate a BuildupResult cache JSON from <2.0 quantity names to >=2.0.

Background
----------
Versions <2.0 of `rad_point_kernel` used `flux` and `dose-{geo}` as quantity
keys and inferred the particle from the source. That made cache files
ambiguous in isolation -- you could not tell whether a `flux` value was a
neutron flux or a photon flux without remembering which Source produced it.

>=2.0 names the particle in the key:

    flux                   ->  flux-neutron     (neutron source)
                          or  flux-photon      (photon source)
    dose-{geo}             ->  dose-{geo}-neutron  (neutron source)
                          or  dose-{geo}-photon   (photon source)
    dose-{geo}-coupled-photon, dose-{geo}-total   unchanged

This script rewrites a cache JSON in-place (or to a new file) using the new
names. It auto-detects the source particle when the cache contains any
`*-coupled-photon` keys (only neutron sources produce secondary photons in
this tool). Otherwise pass `--particle neutron` or `--particle photon`.

Usage
-----
    python tools/migrate_cache.py CACHE.json
    python tools/migrate_cache.py CACHE.json --output NEW.json
    python tools/migrate_cache.py CACHE.json --particle photon
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable


_DOSE_GEOMETRIES = ("AP", "PA", "RLAT", "LLAT", "ROT", "ISO")
_RENAMABLE_BARE_DOSE_NAMES = {f"dose-{g}" for g in _DOSE_GEOMETRIES}
_RESULT_SECTIONS = ("mc", "mc_std_dev", "pk", "buildup")


def _walk_result_sections(obj):
    """Yield every nested ``{quantity_key: value}`` dict that lives under
    one of the BuildupResult section names. Handles both the canonical save
    format (list of result dicts at the top level) and study formats that
    wrap each result in extra metadata (e.g. ``{"water": ..., "result": {...}}``).
    """
    if isinstance(obj, list):
        for item in obj:
            yield from _walk_result_sections(item)
    elif isinstance(obj, dict):
        for k, v in obj.items():
            if k in _RESULT_SECTIONS and isinstance(v, dict):
                yield obj, k, v
            else:
                yield from _walk_result_sections(v)


def _rename_dict(d: dict, particle: str) -> dict:
    """Return a new dict with old keys renamed to particle-explicit form."""
    out: dict = {}
    for k, v in d.items():
        if k == "flux":
            out[f"flux-{particle}"] = v
        elif k in _RENAMABLE_BARE_DOSE_NAMES:
            out[f"{k}-{particle}"] = v
        else:
            # Already explicit, or one of: -coupled-photon, -total. Leave alone.
            out[k] = v
    return out


def _detect_particle(data) -> str | None:
    """A cache with any `*-coupled-photon` key must come from a neutron source;
    only neutron sources produce secondary photons. Otherwise we can't tell.
    """
    for _parent, _name, section in _walk_result_sections(data):
        if any(k.endswith("-coupled-photon") for k in section.keys()):
            return "neutron"
    return None


def _has_old_keys(data) -> bool:
    for _parent, _name, section in _walk_result_sections(data):
        for k in section.keys():
            if k == "flux" or k in _RENAMABLE_BARE_DOSE_NAMES:
                return True
    return False


def migrate(input_path: Path, output_path: Path, particle: str | None) -> None:
    text = input_path.read_text()
    data = json.loads(text)

    if not _has_old_keys(data):
        print(f"{input_path}: already in new format; nothing to do.")
        if input_path != output_path:
            output_path.write_text(text)
        return

    if particle is None:
        particle = _detect_particle(data)
        if particle is None:
            raise SystemExit(
                f"{input_path}: cannot auto-detect the source particle "
                f"(no *-coupled-photon keys present). Re-run with "
                f"--particle neutron or --particle photon."
            )
        print(
            f"{input_path}: detected source particle = {particle!r} "
            f"(secondary photons present)"
        )

    if particle not in ("neutron", "photon"):
        raise SystemExit(
            f"--particle must be 'neutron' or 'photon'; got {particle!r}"
        )

    n_sections = 0
    for parent, section_name, section_dict in _walk_result_sections(data):
        parent[section_name] = _rename_dict(section_dict, particle)
        n_sections += 1
    output_path.write_text(json.dumps(data, indent=2))
    print(
        f"{input_path}: migrated {n_sections} sections; wrote {output_path}"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Migrate a BuildupResult cache JSON from <2.0 quantity names "
            "(flux, dose-{geo}) to >=2.0 names (flux-{particle}, "
            "dose-{geo}-{particle})."
        )
    )
    parser.add_argument("input", type=Path, help="Cache JSON to migrate.")
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Write to this path instead of overwriting INPUT.",
    )
    parser.add_argument(
        "--particle",
        choices=("neutron", "photon"),
        default=None,
        help=(
            "Source particle the cache was produced from. Auto-detected if "
            "the cache contains any *-coupled-photon keys (must be a "
            "neutron source); required otherwise."
        ),
    )
    args = parser.parse_args(argv)

    output = args.output if args.output is not None else args.input
    migrate(args.input, output, args.particle)
    return 0


if __name__ == "__main__":
    sys.exit(main())
