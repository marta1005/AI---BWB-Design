# MADS CTA Migration Plan

## Purpose

This document defines the target architecture and migration plan to absorb the
current `BWB` / `CTA` geometric parameterization workflow into `MADS`.

The end goal is:

- `MADS` becomes the single framework for geometry parameterization,
  multidisciplinary execution, and optimization.
- `BWB` and `CTA` disappear as active runtime dependencies.
- the current `BWB` / `CTA` codebase is treated as a temporary reference source
  for logic, data, and validation only.

This document is intentionally rigorous and prescriptive. It is meant to be the
baseline for implementation decisions.

## Scope

This migration covers:

- CST-based wing parameterization
- spanwise interpolation laws
- section interpolation and local profile transforms
- geometry build and resolved stations
- geometric metrics
- export to `pyGeo` / IGES
- internal volume constraint parsing and evaluation
- the CTA case definition as a reusable case preset

This migration does not initially cover:

- aerodynamic solver changes
- thermal / propulsion / mission solver changes
- deep changes to `GEMSEO` orchestration
- CPACS-driven workflow beyond what is already structurally available in `MADS`

## Architectural Principles

1. `MADS` is the final home of the geometry workflow.
2. `BWB` and `CTA` are reference implementations during migration only.
3. The geometry kernel must be generic.
4. CTA-specific logic must not pollute the generic kernel.
5. Solvers remain the unit of execution.
6. Disciplines remain thin wrappers over solvers.
7. Components represent physical entities; parameters are variables owned by
   components.
8. Complex geometry state is cached in components and not exposed as raw
   `GEMSEO` outputs.
9. Volume constraints are evaluated in a defined CAD reference frame.
10. No early modifications to `MADS` core classes unless proven necessary.

## Decisions Already Fixed

The following decisions are considered agreed and should be preserved:

- The wing parameterization must become native to `MADS`.
- The geometry build must be implemented as one or more `MADS` solvers.
- `pyGeo` export must be handled by a solver.
- Internal volume constraints must be handled by a solver.
- CTA must become a case preset, not a parallel framework.
- The initial implementation must avoid modifying the `MADS` scenario,
  discipline, and base solver core unless integration truly requires it.

## No-Touch Areas In Phase 1

The following files should not be modified in the first migration phase:

- `MADS/src/multiads/scenario/__init__.py`
- `MADS/src/multiads/disciplines/__init__.py`
- `MADS/src/multiads/solvers/__init__.py`

The following areas should also remain untouched initially:

- existing aerodynamic solvers
- existing mission solvers
- existing thermal / propulsion / power supply solvers

## Target Directory Tree

The target additions inside `MADS` are:

```text
MADS/
├── assets/
│   └── cta/
│       ├── bwb_glider.geo
│       ├── B359_internal_volume_constraints_set1.csv
│       └── cta_defaults.json
├── examples/
│   ├── cta_geometry_only.py
│   ├── cta_geometry_constraints.py
│   └── cta_geometry_export.py
└── src/
    └── multiads/
        ├── assembly/
        │   └── parametric_wing.py
        ├── geometry/
        │   ├── __init__.py
        │   ├── model.py
        │   ├── cst_profiles.py
        │   ├── spanwise_laws.py
        │   ├── sections.py
        │   ├── builder.py
        │   ├── state.py
        │   ├── metrics.py
        │   ├── pygeo_export.py
        │   └── internal_volume_constraints.py
        ├── solvers/
        │   └── geometry/
        │       ├── __init__.py
        │       ├── cst_geometry.py
        │       ├── pygeo_export.py
        │       └── internal_volume_constraints.py
        ├── disciplines/
        │   ├── geometry.py
        │   └── packaging.py
        └── cases/
            ├── __init__.py
            ├── cta.py
            ├── cta_laws.py
            └── cta_constraints.py
```

## Responsibilities By File

### `src/multiads/assembly/parametric_wing.py`

Primary responsibility:

- declare the parametric wing component model

Contents:

- `ParametricWing`
- `WingAnchorSection`
- `WingSpanSegment`
- internal cache for resolved geometry state

Owned data:

- planform variables
- CST coefficient variables
- twist and leading-edge vertical offset variables
- interpolation choices
- export options
- cached prepared geometry
- cached export results
- cached packaging results

### `src/multiads/geometry/model.py`

Primary responsibility:

- define reusable geometry configuration dataclasses

Contents:

- `WingGeometryConfig`
- `AnchorSectionDefinition`
- `SamplingConfig`
- `InterpolationConfig`
- optional export and constraint evaluation config containers

### `src/multiads/geometry/cst_profiles.py`

Primary responsibility:

- generate and evaluate CST profiles

Contents:

- CST basis evaluation
- upper / lower profile generation
- leading-edge and trailing-edge handling
- thickness / camber helper functions

### `src/multiads/geometry/spanwise_laws.py`

Primary responsibility:

- host reusable spanwise interpolation laws

Contents:

- linear interpolation
- `PCHIP`
- custom interpolation hooks
- law composition utilities

### `src/multiads/geometry/sections.py`

Primary responsibility:

- resolve section-to-section interpolation and local profile transforms

Contents:

- interpolation between anchor sections
- local upper / lower profile transforms
- optional TE / LE local corrections
- support for case-specific custom laws without embedding case data here

### `src/multiads/geometry/state.py`

Primary responsibility:

- define resolved geometry state containers

Contents:

- `PreparedGeometry`
- `ResolvedStation`
- `GeometryEnvelope`
- resolved section coordinates and caches

### `src/multiads/geometry/builder.py`

Primary responsibility:

- build the full wing geometry from the configuration

Contents:

- `prepare_geometry(...)`
- resolved station generation
- section interpolation orchestration
- envelope extraction
- preparation for export and constraint evaluation

### `src/multiads/geometry/metrics.py`

Primary responsibility:

- compute geometry-derived scalar outputs

Contents:

- span
- area
- enclosed volume
- root chord
- tip chord
- section-based metrics used as solver outputs

### `src/multiads/geometry/pygeo_export.py`

Primary responsibility:

- export resolved geometry to `pyGeo` and IGES

Contents:

- station selection policy
- `pyGeo` adapter logic
- optional meshing frame export support
- IGES writing helpers

### `src/multiads/geometry/internal_volume_constraints.py`

Primary responsibility:

- parse and evaluate internal packaging constraints

Contents:

- `CadReferenceFrame`
- `IndicatorSurfaceSpec`
- `InternalVolumeConstraintSet`
- CSV parsing
- triangulation and sampling
- pointwise constraint evaluation

### `src/multiads/solvers/geometry/cst_geometry.py`

Primary responsibility:

- run the geometry generation workflow as a `MADS` solver

Main class:

- `CSTGeometrySolver`

Expected behavior:

- `parse_variables(...)` collects all geometry-driving variables from the wing
  component
- `_run()` builds the resolved geometry and stores it in component state
- `compute_output()` exposes numeric geometry outputs

Expected outputs:

- wing span
- wing area
- enclosed volume
- root chord
- tip chord
- optional additional geometry indicators

### `src/multiads/solvers/geometry/pygeo_export.py`

Primary responsibility:

- run export as a `MADS` solver

Main class:

- `PyGeoExportSolver`

Expected behavior:

- reuses the prepared geometry stored on the wing component
- performs export to `pyGeo` / IGES
- stores export artifacts and status

Expected outputs:

- export success flag
- number of stations exported
- number of surfaces exported

### `src/multiads/solvers/geometry/internal_volume_constraints.py`

Primary responsibility:

- evaluate internal packaging margins as a `MADS` solver

Main class:

- `InternalVolumeConstraintSolver`

Expected behavior:

- reuses prepared geometry stored on the wing component
- loads the active constraint set
- evaluates all margins in the CAD frame

Expected outputs:

- minimum global margin
- per-surface margins
- optional binary feasibility flag

### `src/multiads/disciplines/geometry.py`

Primary responsibility:

- semantic wrapper over the geometry solver

Main class:

- `Geometry`

Implementation style:

- thin `MADSDiscipline` wrapper, similar to existing discipline files

### `src/multiads/disciplines/packaging.py`

Primary responsibility:

- semantic wrapper over packaging constraint evaluation

Main class:

- `Packaging`

### `src/multiads/cases/cta.py`

Primary responsibility:

- define the CTA case preset

Contents:

- default geometry values
- default active variables
- bounds for optimization
- case-level convenience builders

### `src/multiads/cases/cta_laws.py`

Primary responsibility:

- store CTA-specific geometric laws

Contents:

- any special upper / lower / TE / LE law retained from CTA
- no generic geometry kernel code

### `src/multiads/cases/cta_constraints.py`

Primary responsibility:

- define CTA-specific packaging references

Contents:

- path to constraint CSV
- CTA CAD reference frame
- helper loader for CTA constraint sets

### `assets/cta/*`

Primary responsibility:

- case data only

Contents:

- geometry reference files
- CSV constraints
- optional defaults JSON

### `examples/cta_*.py`

Primary responsibility:

- executable demonstrations and regression entrypoints

Examples:

- geometry only
- geometry plus constraints
- geometry plus export

## Migration Mapping From Current Code

### Generic logic currently in `BWB`

| Current source | Target location |
| --- | --- |
| `parametrization/bwb/specs.py` | `multiads/geometry/model.py` |
| `parametrization/bwb/spanwise_laws.py` | `multiads/geometry/spanwise_laws.py` |
| `parametrization/bwb/sections.py` | `multiads/geometry/sections.py` |
| `parametrization/bwb/builder.py` | `multiads/geometry/builder.py` |
| `parametrization/bwb/exporters.py` | `multiads/geometry/pygeo_export.py` |
| `parametrization/bwb/internal_volume_constraints.py` | `multiads/geometry/internal_volume_constraints.py` |

### CTA-specific logic currently in `CTA`

| Current source | Target location |
| --- | --- |
| `parametrization/CTA/case.py` generic geometry parts | `multiads/geometry/*` |
| `parametrization/CTA/case.py` CTA case values and laws | `multiads/cases/cta.py` and `multiads/cases/cta_laws.py` |
| `parametrization/CTA/design_space.py` | `multiads/cases/cta.py` |
| `parametrization/CTA/internal_volume_constraints.py` case-specific pieces | `multiads/cases/cta_constraints.py` |

## Implementation Phases

### Phase 0: Freeze and reference

Goal:

- treat the current `BWB/CTA` implementation as the reference baseline

Deliverables:

- this migration plan
- preserved reference outputs and plots
- explicit list of case data files

Acceptance criteria:

- no code migrated yet
- target architecture agreed

### Phase 1: Create generic geometry kernel in `MADS`

Goal:

- create `multiads/geometry/*`

Deliverables:

- configuration dataclasses
- CST profile generation
- spanwise laws
- section interpolation
- prepared geometry state

Acceptance criteria:

- geometry kernel runs independently inside `MADS`
- no CTA-specific constants embedded in the generic kernel

### Phase 2: Create the parametric wing component

Goal:

- create `ParametricWing`

Deliverables:

- component class
- variable ownership model
- internal geometry cache

Acceptance criteria:

- all geometry-driving parameters are represented as component-owned variables
- component can hold resolved geometry state without leaking nonnumeric data to
  the optimization layer

### Phase 3: Create the geometry solver

Goal:

- create `CSTGeometrySolver`

Deliverables:

- solver implementation
- numeric outputs
- state caching on component

Acceptance criteria:

- solver runs through the standard `BaseSolver` / `MADSDiscipline` path
- outputs are numerically stable and usable by `MADSScenario`

### Phase 4: Create packaging constraints workflow

Goal:

- create generic internal volume constraint evaluation in `MADS`

Deliverables:

- parser
- CAD reference frame handling
- margin evaluation solver

Acceptance criteria:

- the same CTA CSV can be evaluated fully inside `MADS`
- per-surface margins are reproduced consistently

### Phase 5: Create export workflow

Goal:

- create `PyGeoExportSolver`

Deliverables:

- export adapter
- IGES generation path

Acceptance criteria:

- export can run directly from the resolved geometry cache
- no dependency on the old standalone export scripts

### Phase 6: Recreate CTA as a case preset

Goal:

- add CTA as a case definition, not a framework

Deliverables:

- `cases/cta.py`
- `cases/cta_laws.py`
- `cases/cta_constraints.py`
- `assets/cta/*`

Acceptance criteria:

- CTA can be constructed entirely from within `MADS`
- no runtime dependency on `parametrization/CTA`

### Phase 7: Examples and regression checks

Goal:

- add reproducible execution entrypoints

Deliverables:

- examples for geometry, constraints, and export

Acceptance criteria:

- example runs generate expected outputs
- geometry and margins are traceable against the old reference

## Validation Strategy

The migration must be validated progressively.

Validation levels:

1. local profile reproduction
2. spanwise envelope reproduction
3. station-by-station geometry comparison
4. volume constraint margin comparison
5. `pyGeo` export comparison

Recommended checks:

- profile overlays against current reference stations
- front envelope comparison
- side envelope comparison
- per-surface packaging margins
- exported IGES visual comparison

## Rigor Rules

The following rules must be enforced during implementation:

1. No direct dependency from the final `MADS` implementation back to
   `parametrization/bwb` or `parametrization/CTA`.
2. Generic code must not import CTA-specific case files.
3. CTA-specific laws must stay in `cases/cta_laws.py`.
4. Internal volume constraints must remain case-configurable through
   case-level helpers.
5. Geometry state must be cached in the component, not pushed as arbitrary
   Python objects through the optimization outputs.
6. Every migrated block must have at least one executable example.
7. Every phase must be validated against the current reference before moving on.

## Immediate Next Step

The next implementation step should be:

1. create `multiads/geometry/`
2. create `multiads/assembly/parametric_wing.py`
3. define the minimal dataclasses and component model
4. do not touch solvers yet until those interfaces are stable

That is the safest starting point because it fixes the data model before any
execution logic is migrated.
