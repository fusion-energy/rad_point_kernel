# Coupled neutron and photon dose

A fast neutron source inside a shield generates secondary photons through inelastic scatter and neutron capture. For a complete dose estimate you need both the primary neutron dose and the secondary photon dose. There are two ways to get this:

- **Coupled Monte Carlo** via `compute_buildup` with a `-coupled-photon` quantity - the most accurate option.
- **Analytical PK estimate** via `calculate_secondary_photon_dose_rate` - derives the secondary-gamma source per layer from neutron fluence and capture cross-sections, then transports the gammas to the detector. Faster and no MC required, but less accurate.

## Coupled MC

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [rpk.Layer(thickness=10, material=iron)]

source = rpk.Source(particle="neutron", energy=14.1e6)
results = rpk.compute_buildup(
    geometries=[layers],
    source=source,
    quantities=["dose-AP", "dose-AP-coupled-photon"],
)

# Pulsed DT shot: 1e16 neutrons per shot, dose lands in Sv/shot
r = results[0].scale(strength=1e16)
neutron = r.mc["dose-AP"]
gamma = r.mc["dose-AP-coupled-photon"]
print(f"Neutron dose:    {neutron} Sv/shot")
print(f"Secondary gamma: {gamma} Sv/shot")
print(f"Total:           {neutron + gamma} Sv/shot")
```

## With a manual build-up factor

`calculate_secondary_photon_dose_rate` accepts `neutron_buildup` as a `BuildupModel`. If you already have a neutron build-up factor from tabulated data or a prior run, pass it in:

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
layers = [
    rpk.Layer(thickness=1000),
    rpk.Layer(thickness=10, material=iron),
]

source = rpk.Source(particle="neutron", energy=14.1e6)

B_neutron = 2.5
result = rpk.calculate_secondary_photon_dose_rate(
    layers=layers,
    source=source,
    geometry="AP",
    neutron_buildup=rpk.BuildupModel.constant(B_neutron),
).scale(strength=1e16)
print(f"Neutron dose (B={B_neutron}): {result.neutron_dose_rate} Sv/shot")
print(f"Secondary gamma dose:       {result.secondary_photon_dose_rate} Sv/shot")
print(f"Total dose:                 {result.total_dose_rate} Sv/shot")
```

### Why the manual B is neutron-only

The secondary-gamma pathway is a derived quantity: neutron fluence reaches each material layer, capture and inelastic-scatter cross-sections determine how many secondary photons are born there, and the gammas are transported to the detector analytically. A user-supplied build-up factor only makes sense on the neutron leg - that is where "some of the scattered neutrons also reach this layer" is a first-class thing you tabulate. The secondary-photon leg has no analogous user-facing B input; its attenuation and geometry are built into the capture-gamma model itself. If you need to correct the secondary-photon side independently, prefer the [Coupled MC](#coupled-mc) path instead.

## Notes

- Coupled transport only makes sense for a neutron source; photon sources do not produce secondary neutrons in this tool.
- `dose-AP-coupled-photon` is available alongside any irradiation geometry suffix: `dose-PA-coupled-photon`, `dose-ISO-coupled-photon`, etc.

## Next: extracting build-up factors from MC

The `compute_buildup` call above tallies neutron and secondary-photon doses for one specific geometry. The same machinery can produce reusable build-up factors across a range of thicknesses, energies, or layer stacks - numbers you can plug back into `calculate_flux` or `calculate_dose` as `BuildupModel.constant(B)` on later runs without re-simulating. [Calculate build-up with MC](buildup_mc.md) walks through that workflow.
