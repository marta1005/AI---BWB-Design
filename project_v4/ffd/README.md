# FFD Workspace

This folder is intentionally separated from the current `project_v4` geometric parametrization and design-space workflow.

Its purpose is to host future work related to:

- FFD box generation
- `pyGeo` / `DVGeo`-based mesh deformation
- mesh-morphing experiments
- local deformation studies on top of an already defined BWB geometry

## Scope

This folder is **not** part of the active `project_v4` design-space pipeline.

The current main workflow remains outside this folder:

- geometric parametrization
- CST-based section definition
- GEMSEO design space
- LHS sampling
- post-evaluated geometric constraints

The FFD work should remain isolated here until the workflow is mature enough to connect back to the main project.

## Proposed internal structure

- `docs/`: notes, assumptions, FFD architecture decisions
- `examples/`: small standalone FFD or mesh-deformation experiments
- `example_outputs/`: generated files from those experiments

## Current status

This workspace now contains an isolated FFD workflow based on the local
`pyGeo`/`DVGeometry` stack:

- Plot3D FFD box generation for the current BWB reference geometry
- `DVGeometry` loading of the generated FFD
- embedding of a reference wing pointset
- local deformation tests written with `DVGeometry`
- Tecplot outputs and PNG previews for inspection

It is still separated from the main `project_v4` design-space workflow.
