# Layers

A `Layer` represents a spherical shell with a thickness (in cm) and an optional material. Layers stack outward from the source. For how to build materials, see [Materials](materials.md).

## Creating layers

```python exec="true" source="material-block"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)

# A 10 cm iron layer
shield = rpk.Layer(thickness=10, material=iron)
```

## Void layers

Omit the material to create a void (empty space). Void layers contribute distance for the inverse-square-law calculation but no attenuation.

```python exec="true" source="material-block"
import rad_point_kernel as rpk

void = rpk.Layer(thickness=1000)  # 10 m of empty space
```

## Multi-layer geometry

A geometry is simply a list of layers. They are evaluated outward from the source point.

```python exec="true" source="material-block"
import rad_point_kernel as rpk

iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
lead = rpk.Material(composition={"Pb": 1.0}, density=11.34)

layers = [
    rpk.Layer(thickness=1000),                # 10 m void
    rpk.Layer(thickness=5, material=iron),    # 5 cm iron
    rpk.Layer(thickness=10, material=lead),   # 10 cm lead
]
```
