# Materials and layers

## Creating materials

A `Material` needs a composition dictionary and a density in g/cm3. The composition keys can be element symbols, chemical formulas, or specific nuclide names.

### From elements (mass fractions)

Mass fractions are the default. The values are relative -- they are normalized internally.

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)

concrete = pkc.Material(
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

```python
import rad_point_kernel as pkc

polyethylene = pkc.Material(
    composition={"H": 2, "C": 1},
    density=0.94,
    fraction="atom",
)
```

### From chemical formulas

Formulas like `H2O` are expanded automatically.

```python
import rad_point_kernel as pkc

water = pkc.Material(composition={"H2O": 1.0}, density=1.0)
```

### From specific nuclides

Use nuclide names (element symbol + mass number) when isotopic composition matters.

```python
import rad_point_kernel as pkc

enriched_lithium = pkc.Material(
    composition={"Li6": 0.6, "Li7": 0.4},
    density=0.534,
    fraction="atom",
)
```

## Mixing materials by volume

`Material.volume_mix()` combines two materials by volume fraction. The resulting density and composition are computed from the volume-weighted combination.

```python
import rad_point_kernel as pkc

concrete = pkc.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)
steel = pkc.Material(
    composition={"C": 0.004, "Fe": 0.996},
    density=7.87,
    fraction="atom",
)

# 97% concrete + 3% steel rebar by volume
rebar_concrete = pkc.Material.volume_mix(concrete, 0.97, steel, 0.03)
print(f"Density: {rebar_concrete.density:.3f} g/cm3")
```

## Creating layers

A `Layer` represents a spherical shell with a thickness (in cm) and an optional material. Layers stack outward from the source.

```python
import rad_point_kernel as pkc

iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)

# A 10 cm iron layer
shield = pkc.Layer(thickness=10, material=iron)
```

### Void layers

Omit the material to create a void (empty space). Void layers contribute distance for the inverse-square-law calculation but no attenuation.

```python
import rad_point_kernel as pkc

void = pkc.Layer(thickness=1000)  # 10 m of empty space
```

## Multi-layer geometry

A geometry is simply a list of layers. They are evaluated outward from the source point.

```python
import rad_point_kernel as pkc

water = pkc.Material(composition={"H2O": 1.0}, density=1.0)
concrete = pkc.Material(
    composition={
        "H": 0.01, "O": 0.53, "Si": 0.34,
        "Ca": 0.04, "Al": 0.03, "Fe": 0.01,
    },
    density=2.3,
    fraction="mass",
)
iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)

layers = [
    pkc.Layer(thickness=1000),                    # 10 m void
    pkc.Layer(thickness=5, material=iron),         # 5 cm iron
    pkc.Layer(thickness=30, material=water),        # 30 cm water
    pkc.Layer(thickness=100, material=concrete),    # 100 cm concrete
]
```
