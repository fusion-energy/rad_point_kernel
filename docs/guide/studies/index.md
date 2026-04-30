# Shielding studies

Three worked examples that exercise the full pipeline - material mixing, coupled-transport Monte Carlo, analytical-form build-up fitting, and plotting. Each page is self-contained: the code runs as part of the docs build and the plot below it is generated from a committed MC cache so rebuilds are fast.

- **[Single-layer study](single_layer.md)** - dose vs shield thickness for two concrete mixes, extrapolated from thin-shield MC.
- **[Multi-layer study](multi_layer.md)** - water + concrete, 2D thickness sweep with `BuildupFit` (thin-plate-spline RBF) extrapolation.
- **[Variable material study](variable_material.md)** - fixed 50 cm concrete with varying steel volume fraction and boron weight fraction, shown as three dose maps.

To re-run the Monte Carlo from scratch, delete the relevant `docs/assets/studies/*_cache.json` and rebuild.
