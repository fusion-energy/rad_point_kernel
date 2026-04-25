"""Secondary photon dose from coupled neutron-photon MC."""

import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [pkc.Layer(thickness=10, material=iron)]
source = pkc.Source("neutron", 14.1e6)
SOURCE_STRENGTH = 1e12

print("Running coupled MC...")
results = pkc.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP-coupled-photon"],
    particles_per_batch=10_000,
    max_batches=50,
    trigger_rel_err=0.05,
    cross_sections="/home/jon/nuclear_data/cross_sections.xml",
)

r = results[0]
mc_photon = r.mc["dose-AP-coupled-photon"] * SOURCE_STRENGTH
print(f"Secondary photon dose (MC): {mc_photon:.4e} Sv/hr")
