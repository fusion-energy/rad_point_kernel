"""Photon dose with MC-computed buildup factor (Cs-137 source).

Cs-137 is a continuous source. The activity is in photons per second, so we
multiply by 3600 to scale by photons per hour and land in the conventional
Sv/hr dose-rate unit.
"""

import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]
source = rpk.Source("photon", 662e3)
PARTICLES_PER_HOUR = 1e12 * 3600  # 1e12 photons/sec activity

pk = rpk.calculate_dose(layers=layers, source=source, geometry="AP").scale(
    strength=PARTICLES_PER_HOUR
)
print(f"PK dose (no buildup): {pk.dose_rate} Sv/hr")

print("Running MC...")
results = rpk.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP"],
    particles_per_batch=10_000,
    max_batches=50,
    trigger_rel_err=0.05,
)

r = results[0].scale(strength=PARTICLES_PER_HOUR)
print(f"Buildup factor: {r.buildup['dose-AP']}")

corrected = rpk.calculate_dose(
    layers=layers, source=source, geometry="AP", buildup=r
).scale(strength=PARTICLES_PER_HOUR)
print(f"PK dose with buildup: {corrected.dose_rate} Sv/hr")
print(f"MC dose (reference):  {r.mc['dose-AP']} Sv/hr")
