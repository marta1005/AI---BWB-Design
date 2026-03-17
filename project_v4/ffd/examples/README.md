# FFD Examples

Use this folder for isolated executable examples related to:

- building an FFD box around the wing
- embedding a surface or mesh
- testing deformation variables
- checking mesh quality after deformation

## Available examples

### 1. Build and load the reference semispan FFD box with DVGeometry

```bash
python project_v4/ffd/examples/build_reference_ffd_box.py
```

Outputs:

- `project_v4/ffd/example_outputs/reference_ffd_box/reference_semispan_ffd.xyz`
- `project_v4/ffd/example_outputs/reference_ffd_box/reference_semispan_ffd.json`
- `project_v4/ffd/example_outputs/reference_ffd_box/reference_semispan_baseline_ffd.dat`
- `project_v4/ffd/example_outputs/reference_ffd_box/reference_semispan_baseline_surface_wing_surface.dat`
- `project_v4/ffd/example_outputs/reference_ffd_box/reference_semispan_ffd.png`

The script does four things:

- generates the Plot3D FFD box
- loads it with `DVGeometry`
- embeds a reference wing surface pointset
- writes Tecplot outputs and a PNG preview

### 2. Build and load the reference full-wing FFD box with DVGeometry

```bash
python project_v4/ffd/examples/build_fullwing_ffd_box.py
```

Outputs:

- `project_v4/ffd/example_outputs/reference_ffd_box/reference_fullwing_ffd.xyz`
- `project_v4/ffd/example_outputs/reference_ffd_box/reference_fullwing_ffd.json`
- `project_v4/ffd/example_outputs/reference_ffd_box/reference_fullwing_baseline_ffd.dat`
- `project_v4/ffd/example_outputs/reference_ffd_box/reference_fullwing_baseline_surface_wing_surface.dat`
- `project_v4/ffd/example_outputs/reference_ffd_box/reference_fullwing_ffd.png`

### 3. Apply a small local deformation with DVGeometry

This script requires the local `pySpline` + `pyGeo` stack to be working.

```bash
python project_v4/ffd/examples/deform_reference_ffd_dvgeo.py
python project_v4/ffd/examples/deform_reference_ffd_dvgeo.py --full-wing
```

Outputs:

- baseline and deformed Plot3D FFD files
- `DVGeometry` Tecplot files for the FFD lattice
- baseline and deformed embedded surface point sets
- a combined PNG comparing the undeformed and deformed wing/FFD
- a dedicated 3D PNG with control points and deformed points
- an animated GIF showing the deformation progression

## What you need to visualize this

### Baseline and deformed FFD preview

You need:

- the local `project/pyspline` build to work with your Python
- the local `project/pygeo` package available
- `mpi4py`
- `matplotlib`

Optional but useful:

- Tecplot, if you want to inspect the `.dat` files written by `DVGeometry`

If you do not have Tecplot, the scripts still generate PNG previews from the `DVGeometry`-updated pointsets.
