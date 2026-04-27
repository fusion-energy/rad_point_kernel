"""Neutron dose with MC-computed buildup factor (pulsed DT source).

A 14.1 MeV neutron pulse from an inertial-confinement shot. The strength
is in neutrons per shot, so the resulting dose is in Sv per shot.
"""

import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]
source = rpk.Source("neutron", 14.1e6)
PARTICLES_PER_SHOT = 1e16

pk = rpk.calculate_dose(layers=layers, source=source, geometry="AP").scale(
    strength=PARTICLES_PER_SHOT
)
print(f"PK dose (no buildup): {pk.dose_rate} Sv/shot")

print("Running MC...")
results = rpk.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP"],
    particles_per_batch=10_000,
    max_batches=50,
    trigger_rel_err=0.05,
)

r = results[0].scale(strength=PARTICLES_PER_SHOT)
print(f"Buildup factor: {r.buildup['dose-AP']}")

corrected = rpk.calculate_dose(
    layers=layers, source=source, geometry="AP", buildup=r
).scale(strength=PARTICLES_PER_SHOT)
print(f"PK dose with buildup: {corrected.dose_rate} Sv/shot")
print(f"MC dose (reference):  {r.mc['dose-AP']} Sv/shot")
