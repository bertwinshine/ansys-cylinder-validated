# ANSYS validation of the pressure vessel case (Lamé check)

I replicated the cylindrical pressure vessel case in ANSYS and checked the
result against the exact Lamé solution for a thick-walled cylinder. The idea
is to have an independent reference for this case that's verified against the
analytical answer, not just against another solver.

The ANSYS stresses match Lamé to about 0.05% mean error across the field.

## Case

- Inner radius 0.9 m, outer radius 1.0 m, length 3.0 m
- Internal pressure 3e5 Pa (3 bar)
- E = 3 GPa, nu = 0.3 (placeholder values from the case I replicated)
- Linear static structural, pressure on the inner wall

Note: E and nu are placeholder values, not a real material. Stress doesn't
depend on E anyway (only displacement does), so the check holds.

## Results vs Lamé

| | ANSYS | Lamé |
|---|---|---|
| Hoop, inner wall | 2.857 MPa | 2.858 MPa |
| Hoop, outer wall | 2.557 MPa | 2.558 MPa |
| Radial, inner wall | -0.299 MPa | -0.300 MPa (= -p) |
| Radial, outer wall | ~0 | 0 |
| Von Mises peak | 3.021 MPa | 3.02 MPa |

Mean error over the field: 0.05%.

## Files

- `merge_and_validate.py`: merges the ANSYS exports, runs the Lamé check, writes the CSV
- `pressure_vessel_reference.csv`: the final field (node, x, y, z, radius, radial, hoop, von_mises, deformation)
- `ansys_exports/`: the raw ANSYS text exports (coordinates + stresses + deformation)

## How to run

Needs Python with numpy and pandas.

```
python merge_and_validate.py
```

It reads the files in `ansys_exports/`, prints the validation result, and
writes `pressure_vessel_reference.csv`. Geometry (inner/outer radius) is read
from the data, so it works if the size changes.

## Limitations

- Mesh was kept under the ANSYS Student node limit, so it's not super fine, but
  the error is already 0.05%. A finer mesh on a full license would tighten it.
- This is the open-end case (zero axial stress), which is what Lamé assumes.
- It's one operating point (one E, p). For surrogate training you'd want a set
  of these across the parameter range, done the same way.

## Credit

The ANSYS model and validation are my work. I wrote the Python with help from an
AI assistant and have checked and understand all of it.

Bertwin Kurisinkal Shine, TU Chemnitz
