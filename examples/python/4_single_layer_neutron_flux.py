"""Neutron flux with MC-computed buildup factor (D-T + D-D source)."""

import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [pkc.Layer(thickness=10, material=iron)]
source = pkc.Source("neutron", [(14.06e6, 95), (2.45e6, 5)])
SOURCE_STRENGTH = 1e12

# PK flux (no buildup)
pk = pkc.calculate_flux(SOURCE_STRENGTH, layers, source)
print(f"PK flux (no buildup): {pk.uncollided_flux:.4e} n/cm2/s")

# MC buildup
print("Running MC with D-T + D-D spectrum...")
results = pkc.compute_buildup(
    geometries=[layers], source=source, quantities=["flux"],
    particles_per_batch=10_000, max_batches=50, trigger_rel_err=0.05,
)

r = results[0]
print(f"Buildup factor: {r.buildup['flux']:.4f}")

# Apply buildup
corrected = pkc.calculate_flux(SOURCE_STRENGTH, layers, source, buildup=r)
print(f"PK flux with buildup: {corrected.uncollided_flux:.4e} n/cm2/s")
print(f"MC flux (reference):  {r.mc['flux'] * SOURCE_STRENGTH:.4e} n/cm2/s")
