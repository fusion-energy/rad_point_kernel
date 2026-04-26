"""Neutron transmission through a single material layer."""

import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]
source = rpk.Source("neutron", 14.1e6)

frac = rpk.calculate_transmission(layers, source)
print(f"Neutron transmission through 10 cm Fe at 14.1 MeV: {frac:.4e}")
