"""Photon transmission through a single material layer (Cs-137 source)."""

import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]
source = rpk.Source("photon", 662e3)

frac = rpk.calculate_transmission(layers, source)
print(f"Photon transmission through 10 cm Fe at 662 keV (Cs-137): {frac:.4e}")
