# AI BWB Design

This repository contains the current codebase for geometric parametrization of blended-wing-body (BWB) aircraft.

## Repository Structure

- `parametrization/`: active aircraft and BWB parametrization packages
- `project_v4/`: previous active BWB workflow kept for reference
- `project/`: local `pyGeo`/`pySpline`/`baseclasses` dependencies used by the geometry builders

## Current Focus

`parametrization` is the main working area. It includes:

- planform parametrization based on span partitions, chord law, and segment sweeps
- CST-based airfoil sections with independent coefficients at each control section
- spanwise interpolation with `pySpline`
- loft and IGES export with `pyGeo`
- design-space tools and GEMSEO integration

## Documentation

Start here:

- `parametrization/README.md`
- `parametrization/CTA/README.md`

Additional notes:

- `parametrization/CTA/docs/README.md`
- `project_v4/docs/geometric_parametrization_tables.md`
- `project_v4/docs/runtime_requirements.md`
- `project_v4/docs/gemseo_install.md`

## Note

The old standalone `cst/` utilities have been folded into `parametrization/shared`.
