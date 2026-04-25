"""Photon dose with MC-computed buildup factor."""

import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
layers = [pkc.Layer(thickness=10, material=iron)]
source = pkc.Source("photon", 662e3)
SOURCE_STRENGTH = 1e12

pk = pkc.calculate_dose(SOURCE_STRENGTH, layers, source, "AP")
print(f"PK dose (no buildup): {pk.dose_rate:.4e} Sv/hr")

print("Running MC...")
results = pkc.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP"],
    particles_per_batch=10_000,
    max_batches=50,
    trigger_rel_err=0.05,
)

r = results[0]
print(f"Buildup factor: {r.buildup['dose-AP']:.4f}")

corrected = pkc.calculate_dose(SOURCE_STRENGTH, layers, source, "AP", buildup=r)
print(f"PK dose with buildup: {corrected.dose_rate:.4e} Sv/hr")
print(f"MC dose (reference):  {r.mc['dose-AP'] * SOURCE_STRENGTH:.4e} Sv/hr")
