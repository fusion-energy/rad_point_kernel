"""Making materials, layers, and sources."""

import rad_point_kernel as rpk

# --- Materials from elements (mass fractions) ---
iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)

concrete = rpk.Material(
    composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
    density=2.3,
    fraction="mass",
)

# --- Materials from elements (atom fractions) ---
polyethylene = rpk.Material(composition={"H": 2, "C": 1}, density=0.94, fraction="atom")

# --- Materials from chemical formulas ---
water = rpk.Material(composition={"H2O": 1.0}, density=1.0)

# --- Materials from specific nuclides ---
enriched_lithium = rpk.Material(
    composition={"Li6": 0.6, "Li7": 0.4}, density=0.534, fraction="atom"
)

# --- Mixing materials by volume ---
steel = rpk.Material(
    composition={"C": 0.004, "Fe": 0.996}, density=7.87, fraction="atom"
)
rebar_concrete = rpk.Material.volume_mix(concrete, 0.97, steel, 0.03)
print(f"Rebar concrete density: {rebar_concrete.density:.3f} g/cm3")

# --- Sources ---
cs137 = rpk.Source("photon", 662e3)
co60 = rpk.Source("photon", [(1173e3, 1.0), (1333e3, 1.0)])
dt_neutron = rpk.Source("neutron", 14.06e6)
dt_dd = rpk.Source("neutron", [(14.06e6, 95), (2.45e6, 5)])

print(f"Sources: {cs137}, {co60}, {dt_neutron}, {dt_dd}")

# --- Layers ---
single = [rpk.Layer(thickness=100, material=concrete)]
multi = [
    rpk.Layer(thickness=500),
    rpk.Layer(thickness=20, material=polyethylene),
    rpk.Layer(thickness=100, material=concrete),
]

for name, layers in [("Single", single), ("Multi", multi)]:
    total = sum(l.thickness for l in layers)
    print(f"{name}: {len(layers)} layers, {total:.0f} cm total")
