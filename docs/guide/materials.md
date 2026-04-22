# Materials

## Creating materials

A `Material` needs a composition dictionary and a density in g/cm3. The composition keys can be element symbols, chemical formulas, or specific nuclide names.

### From elements (mass fractions)

Mass fractions are the default. The values are relative; they are normalized internally.

```python exec="true" source="material-block"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)

concrete = rpk.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)
```

### From elements (atom fractions)

Set `fraction="atom"` to specify atom (number) fractions instead.

```python exec="true" source="material-block"
import rad_point_kernel as rpk

polyethylene = rpk.Material(
    composition={"H": 2, "C": 1},
    density=0.94,
    fraction="atom",
)
```

### From chemical formulas

Formulas like `H2O` are expanded automatically.

```python exec="true" source="material-block"
import rad_point_kernel as rpk

water = rpk.Material(composition={"H2O": 1.0}, density=1.0)
```

### From specific nuclides

Use nuclide names (element symbol + mass number) when isotopic composition matters.

```python exec="true" source="material-block"
import rad_point_kernel as rpk

enriched_lithium = rpk.Material(
    composition={"Li6": 0.6, "Li7": 0.4},
    density=0.534,
    fraction="atom",
)
```

## Mixing materials by volume

`Material.volume_mix()` combines two materials by volume fraction. The resulting density and composition are computed from the volume-weighted combination.

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

concrete = rpk.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)
steel = rpk.Material(
    composition={"C": 0.004, "Fe": 0.996},
    density=7.87,
    fraction="atom",
)

# 97% concrete + 3% steel rebar by volume
rebar_concrete = rpk.Material.volume_mix(concrete, 0.97, steel, 0.03)
print(f"Density: {rebar_concrete.density} g/cm3")
```

## Mixing materials by mass

`Material.mass_mix()` combines two materials by mass fraction. Use this for additives that are specified by weight percent - boron in borated concrete, tungsten in heavy-loaded polyethylene, etc.

The mixed density is volume-additive (each component retains its own bulk density):

    1 / rho_mix = w_a / rho_a + w_b / rho_b

```python exec="true" source="material-block" result="text"
import rad_point_kernel as rpk

concrete = rpk.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)
boron = rpk.Material(composition={"B": 1.0}, density=2.34, fraction="mass")

# Concrete with 3 wt% natural boron
borated = rpk.Material.mass_mix(concrete, 0.97, boron, 0.03)
print(f"Density: {borated.density} g/cm3")
```
