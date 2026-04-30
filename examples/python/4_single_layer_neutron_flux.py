"""Neutron flux with MC-computed buildup factor (D-T + D-D source).

This represents a pulsed-fusion shot, so the strength is in neutrons per shot
and the resulting flux is in neutrons/cm^2 per shot. Swap the strength unit
(without re-running MC) to model a steady-state generator instead.
"""

import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]
source = rpk.Source("neutron", [(14.06e6, 95), (2.45e6, 5)])
PARTICLES_PER_SHOT = 1e16

# PK flux (no buildup)
pk = rpk.calculate_flux(layers=layers, source=source).scale(strength=PARTICLES_PER_SHOT)
print(f"PK flux (no buildup): {pk.flux} n/cm2/shot")

# MC buildup
print("Running MC with D-T + D-D spectrum...")
results = rpk.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["flux-neutron"],
    particles_per_batch=10_000,
    max_batches=50,
    trigger_rel_err=0.05,
)

r = results[0].scale(strength=PARTICLES_PER_SHOT)
print(f"Buildup factor: {r.buildup['flux-neutron']}")

# Apply buildup
corrected = rpk.calculate_flux(layers=layers, source=source, buildup=r).scale(
    strength=PARTICLES_PER_SHOT
)
print(f"PK flux with buildup: {corrected.flux} n/cm2/shot")
print(f"MC flux (reference):  {r.mc['flux-neutron']} n/cm2/shot")
