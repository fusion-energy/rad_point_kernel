"""Photon flux with MC-computed buildup factor (Co-60 source).

Co-60 is a continuous source (a radioisotope), so the strength is in
photons per second and the resulting flux lands in photons/cm^2/s.
"""

import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]
source = rpk.Source("photon", [(1173e3, 1.0), (1333e3, 1.0)])
PARTICLES_PER_SECOND = 1e12  # photons/s

# PK flux (no buildup)
pk = rpk.calculate_flux(layers=layers, source=source).scale(strength=PARTICLES_PER_SECOND)
print(f"PK flux (no buildup): {pk.uncollided_flux} photons/cm2/s")

# MC buildup
print("Running MC with Co-60 spectrum...")
results = rpk.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["flux"],
    particles_per_batch=10_000,
    max_batches=50,
    trigger_rel_err=0.05,
)

r = results[0].scale(strength=PARTICLES_PER_SECOND)
print(f"Buildup factor: {r.buildup['flux']}")

# Apply buildup
corrected = rpk.calculate_flux(layers=layers, source=source, buildup=r).scale(
    strength=PARTICLES_PER_SECOND
)
print(f"PK flux with buildup: {corrected.uncollided_flux} photons/cm2/s")
print(f"MC flux (reference):  {r.mc['flux']} photons/cm2/s")
