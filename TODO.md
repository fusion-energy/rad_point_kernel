# TODO

## Docs

- Regenerate `docs/assets/buildup_example.png` using water + concrete (currently still shows the old polyethylene + concrete grid). Driven by a 2D MC run: `mc_water = [0, 10, 20, 30, 40, 50]` x `mc_conc = [10, 50, 100, 200, 400]` = 30 coupled neutron-photon geometries. Uses `examples/python/multiple_layer_study.py` as the MC driver; needs a new plotting pass that outputs B (not dose) vs concrete thickness with water as the family axis, matching the caption in `docs/theory/buildup_factors.md`. Nuclear data setup: `pip install openmc_data && download_endf -r b8.1` (see `docs/index.md`), then `export OPENMC_CROSS_SECTIONS=~/nuclear_data/cross_sections.xml`.

## Verification and validation

- Add a verification case comparing point-kernel buildup factors against published ANS-6.4.3 / Harima polynomial fits (this is verification, not validation - those tables are themselves computational).
- Add validation against measured shielding benchmarks (SINBAD, ICSBEP, ORNL) once source data is pulled in.
