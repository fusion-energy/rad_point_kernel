"""Secondary photon dose from coupled neutron-photon MC."""

import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]
source = rpk.Source("neutron", 14.1e6)
SOURCE_STRENGTH = 1e12

print("Running coupled MC...")
results = rpk.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP-coupled-photon"],
    particles_per_batch=10_000,
    max_batches=50,
    trigger_rel_err=0.05,
)

r = results[0]
mc_photon = r.mc["dose-AP-coupled-photon"] * SOURCE_STRENGTH
print(f"Secondary photon dose (MC): {mc_photon:.4e} Sv/hr")
