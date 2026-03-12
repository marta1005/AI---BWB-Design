# AI BWB Design

This repository contains the current codebase for geometric parametrization of blended-wing-body (BWB) aircraft.

## Repository Structure

- `project_v4/`: active BWB geometry framework
- `cst/`: CST-related utilities and airfoil support code

## Current Focus

`project_v4` is the main working version. It includes:

- planform parametrization based on span partitions, chord law, and segment sweeps
- CST-based airfoil sections with independent coefficients at each control section
- spanwise interpolation with `pySpline`
- loft and IGES export with `pyGeo`
- design-space tools and GEMSEO integration

## Documentation

Start here:

- `project_v4/README.md`

Additional notes:

- `project_v4/docs/geometric_parametrization_tables.md`
- `project_v4/docs/runtime_requirements.md`
- `project_v4/docs/gemseo_install.md`

## Note

Only the active `project_v4` workflow and the `cst` utilities are included in this repository snapshot.
