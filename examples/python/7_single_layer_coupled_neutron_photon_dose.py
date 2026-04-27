"""Secondary photon dose from coupled neutron-photon MC (pulsed DT shot).

The neutron source is a single ICF-style burn, so we report dose per shot.
"""

import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]
source = rpk.Source("neutron", 14.1e6)
PARTICLES_PER_SHOT = 1e16

print("Running coupled MC...")
results = rpk.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP-coupled-photon"],
    particles_per_batch=10_000,
    max_batches=50,
    trigger_rel_err=0.05,
)

r = results[0].scale(strength=PARTICLES_PER_SHOT)
print(f"Secondary photon dose (MC): {r.mc['dose-AP-coupled-photon']} Sv/shot")
