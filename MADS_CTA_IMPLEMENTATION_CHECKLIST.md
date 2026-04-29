# MADS CTA Implementation Checklist

## Purpose

This document translates the architectural migration plan into an execution
checklist.

It is intentionally strict:

- exact order of implementation
- file-by-file creation sequence
- dependency rules
- validation gates
- explicit stop conditions

This checklist should be followed before any broad refactor starts.

## Governing Rules

1. Do not modify `MADS` core orchestration first.
2. Do not create runtime dependencies from the final `MADS` code back to
   `parametrization/bwb` or `parametrization/CTA`.
3. Do not mix CTA-specific laws into generic geometry files.
4. Do not expose complex Python geometry objects through optimizer outputs.
5. Do not proceed to a new phase before the previous phase is validated.

## No-Touch Files In Initial Phases

These files must remain unchanged until explicitly justified:

- `MADS/src/multiads/scenario/__init__.py`
- `MADS/src/multiads/disciplines/__init__.py`
- `MADS/src/multiads/solvers/__init__.py`

## Phase Overview

The implementation must follow this exact phase order:

1. Skeleton and namespaces
2. Generic geometry data model
3. Parametric wing component
4. Generic geometry kernel
5. Geometry solver
6. Volume-constraint kernel
7. Volume-constraint solver
8. `pyGeo` export kernel
9. `pyGeo` export solver
10. CTA case preset
11. Example scripts and regression checks
12. Decommission dependency on `BWB/CTA`

---

## Phase 1: Skeleton And Namespaces

### Goal

Create the empty target structure inside `MADS` without moving logic yet.

### Files To Create In Order

1. `MADS/src/multiads/geometry/__init__.py`
2. `MADS/src/multiads/solvers/geometry/__init__.py`
3. `MADS/src/multiads/cases/__init__.py`
4. `MADS/src/multiads/disciplines/geometry.py`
5. `MADS/src/multiads/disciplines/packaging.py`
6. `MADS/src/multiads/assembly/parametric_wing.py`
7. `MADS/src/multiads/geometry/model.py`
8. `MADS/src/multiads/geometry/state.py`
9. `MADS/src/multiads/geometry/metrics.py`
10. `MADS/src/multiads/geometry/cst_profiles.py`
11. `MADS/src/multiads/geometry/spanwise_laws.py`
12. `MADS/src/multiads/geometry/sections.py`
13. `MADS/src/multiads/geometry/builder.py`
14. `MADS/src/multiads/geometry/pygeo_export.py`
15. `MADS/src/multiads/geometry/internal_volume_constraints.py`
16. `MADS/src/multiads/solvers/geometry/cst_geometry.py`
17. `MADS/src/multiads/solvers/geometry/pygeo_export.py`
18. `MADS/src/multiads/solvers/geometry/internal_volume_constraints.py`
19. `MADS/src/multiads/cases/cta.py`
20. `MADS/src/multiads/cases/cta_laws.py`
21. `MADS/src/multiads/cases/cta_constraints.py`
22. `MADS/examples/cta_geometry_only.py`
23. `MADS/examples/cta_geometry_constraints.py`
24. `MADS/examples/cta_geometry_export.py`
25. `MADS/assets/cta/cta_defaults.json`

### Minimal Content Allowed

At this phase, these files may contain:

- docstrings
- empty dataclasses
- empty class shells
- import-safe placeholders

They must not yet contain:

- copied CTA logic
- direct imports from `parametrization/bwb`
- direct imports from `parametrization/CTA`

### Exit Criteria

- all files exist
- all imports are syntactically valid
- no circular imports introduced
- no logic migrated yet

### Stop Condition

Do not start Phase 2 if any created file has unclear ownership or overlapping
responsibility.

---

## Phase 2: Generic Geometry Data Model

### Goal

Define the reusable geometry configuration and resolved-state model.

### Files To Fill

1. `MADS/src/multiads/geometry/model.py`
2. `MADS/src/multiads/geometry/state.py`
3. `MADS/src/multiads/geometry/metrics.py`

### Exact Work Order

1. Define configuration dataclasses in `model.py`
2. Define resolved-state dataclasses in `state.py`
3. Define metric result containers in `metrics.py`
4. Review names and ownership before any solver logic is added

### Classes To Introduce

In `model.py`:

- `WingGeometryConfig`
- `AnchorSectionDefinition`
- `SamplingConfig`
- `InterpolationConfig`
- `ExportConfig`
- `ConstraintEvaluationConfig`

In `state.py`:

- `ResolvedStation`
- `GeometryEnvelope`
- `PreparedGeometry`

In `metrics.py`:

- `GeometryMetricSet`

### Validation

Check that the model is capable of representing:

- anchor stations
- CST coefficients
- spanwise interpolation choices
- local profile transforms
- export settings
- constraint evaluation settings

### Exit Criteria

- data model covers all current CTA geometry inputs
- no CTA-specific constants are embedded
- no solver code is required to instantiate the objects

### Stop Condition

Do not proceed if any current CTA input still lacks a clear home in the data
model.

---

## Phase 3: Parametric Wing Component

### Goal

Create the component that will own all wing parameters and all geometry state.

### File To Fill

1. `MADS/src/multiads/assembly/parametric_wing.py`

### Exact Work Order

1. Define `WingAnchorSection`
2. Define `WingSpanSegment`
3. Define `ParametricWing`
4. Add state/cache fields to `ParametricWing`
5. Add helper methods to build configuration objects from component variables

### Required Responsibilities

`ParametricWing` must own:

- planform variables
- CST variables
- twist variables
- LE vertical offset variables
- interpolation settings
- cached prepared geometry
- cached export results
- cached packaging results

### Validation

Confirm that a CTA wing can be represented by `ParametricWing` without needing
external mutable global state.

### Exit Criteria

- all geometry-driving parameters belong to the wing component
- all nonnumeric state is stored as component cache, not as optimization output

### Stop Condition

Do not proceed if any geometry parameter still lives only in ad hoc case logic.

---

## Phase 4: Generic Geometry Kernel

### Goal

Implement the reusable geometry logic, independent of CTA.

### Files To Fill

1. `MADS/src/multiads/geometry/cst_profiles.py`
2. `MADS/src/multiads/geometry/spanwise_laws.py`
3. `MADS/src/multiads/geometry/sections.py`
4. `MADS/src/multiads/geometry/builder.py`

### Exact Work Order

1. Port CST profile evaluation into `cst_profiles.py`
2. Port generic spanwise laws into `spanwise_laws.py`
3. Port section interpolation logic into `sections.py`
4. Port prepared-geometry build orchestration into `builder.py`

### Source Mapping

Use as reference only:

- `parametrization/bwb/specs.py`
- `parametrization/bwb/spanwise_laws.py`
- `parametrization/bwb/sections.py`
- `parametrization/bwb/builder.py`

### Rules

- strip CTA case constants out of the generic kernel
- keep law hooks generic
- keep profile transform interfaces generic

### Validation

At the end of this phase, build a minimal synthetic wing case entirely from the
new `MADS` geometry kernel.

### Exit Criteria

- `prepare_geometry(...)`-equivalent behavior exists in `MADS`
- the generic kernel can build a resolved wing without any CTA-specific code

### Stop Condition

Do not proceed if the new geometry kernel still imports `parametrization/CTA`.

---

## Phase 5: Geometry Solver

### Goal

Wrap the generic geometry kernel in a `MADS` solver.

### File To Fill

1. `MADS/src/multiads/solvers/geometry/cst_geometry.py`

### Exact Work Order

1. Define `CSTGeometrySolver`
2. Implement `parse_variables(...)`
3. Implement `set_state(...)` only if needed
4. Implement `_run()`
5. Implement `compute_output()`
6. Stub `compute_sensitivities()`

### Required Inputs

The solver must read from `ParametricWing`:

- planform variables
- CST variables
- law settings
- export / constraint settings only if needed downstream

### Required Outputs

The solver must expose numeric outputs only, such as:

- `wing.span`
- `wing.area`
- `wing.volume`
- `wing.root_chord`
- `wing.tip_chord`

### Validation

Create and run a minimal example:

- instantiate wing component
- instantiate geometry solver
- run through `MADSDiscipline`

### Exit Criteria

- geometry builds correctly through `BaseSolver` + `MADSDiscipline`
- component cache contains resolved geometry after execution

### Stop Condition

Do not proceed if solver outputs still depend on ad hoc script-level state.

---

## Phase 6: Volume-Constraint Kernel

### Goal

Port the generic packaging-constraint engine into `MADS`.

### File To Fill

1. `MADS/src/multiads/geometry/internal_volume_constraints.py`

### Exact Work Order

1. Port `CadReferenceFrame`
2. Port `IndicatorSurfaceSpec`
3. Port `InternalVolumeConstraintSet`
4. Port CSV parser
5. Port triangulation and sampling
6. Port bound evaluation
7. Port result dataclasses

### Source Mapping

Use as reference only:

- `parametrization/bwb/internal_volume_constraints.py`

### Validation

Load the current CTA CSV and verify:

- number of surfaces
- number of points per surface
- reference frame is represented correctly

### Exit Criteria

- generic packaging constraints can be loaded without any CTA runtime import

### Stop Condition

Do not proceed if the generic kernel contains CTA-specific file paths.

---

## Phase 7: Volume-Constraint Solver

### Goal

Expose packaging evaluation as a solver.

### File To Fill

1. `MADS/src/multiads/solvers/geometry/internal_volume_constraints.py`

### Exact Work Order

1. Define `InternalVolumeConstraintSolver`
2. Implement `parse_variables(...)`
3. Implement `_run()` using cached geometry
4. Implement `compute_output()`
5. Stub `compute_sensitivities()`

### Required Outputs

At minimum:

- `minimum_margin`
- `feasibility_flag`

Recommended per-surface outputs:

- `cabin_floor_margin`
- `door_ceiling_margin`
- `cabin_aisle_ceiling_margin`
- `cabin_seat_ceiling_margin`
- `cargo1_floor_margin`
- `cargo2_floor_margin`
- `cargo2_ceiling_margin`
- `nlg_margin`
- `mlg1_margin`
- `mlg2_margin`

### Validation

Run the solver on the CTA preset and compare margins against the current
reference workflow.

### Exit Criteria

- margins are reproducible inside `MADS`
- feasibility status is available as numeric solver output

### Stop Condition

Do not proceed if packaging evaluation requires old plotting or export scripts.

---

## Phase 8: `pyGeo` Export Kernel

### Goal

Port export logic into a generic reusable geometry module.

### File To Fill

1. `MADS/src/multiads/geometry/pygeo_export.py`

### Exact Work Order

1. port generic station-to-`pyGeo` adapter logic
2. port IGES writing helpers
3. port optional meshing-frame export helpers

### Source Mapping

Use as reference only:

- `parametrization/bwb/exporters.py`

### Validation

Confirm that the export kernel can:

- use resolved stations
- export all stations or selected anchor stations
- write IGES consistently

### Exit Criteria

- export logic is reusable and independent from CTA scripts

### Stop Condition

Do not proceed if export logic still depends on old CTA plotting modules.

---

## Phase 9: `pyGeo` Export Solver

### Goal

Expose geometry export as a solver.

### File To Fill

1. `MADS/src/multiads/solvers/geometry/pygeo_export.py`

### Exact Work Order

1. define `PyGeoExportSolver`
2. implement `parse_variables(...)`
3. implement `_run()` using cached geometry
4. implement `compute_output()`
5. stub `compute_sensitivities()`

### Required Outputs

- `export_ok`
- `n_sections_exported`
- `n_surfaces_exported`

### Validation

Run the solver using the CTA preset and verify that the exported IGES matches
the reference process.

### Exit Criteria

- export works fully from inside `MADS`

### Stop Condition

Do not proceed if export still relies on standalone CTA scripts.

---

## Phase 10: CTA Case Preset

### Goal

Recreate CTA as a case definition, not as a framework.

### Files To Fill

1. `MADS/src/multiads/cases/cta.py`
2. `MADS/src/multiads/cases/cta_laws.py`
3. `MADS/src/multiads/cases/cta_constraints.py`
4. `MADS/assets/cta/cta_defaults.json`
5. `MADS/assets/cta/bwb_glider.geo`
6. `MADS/assets/cta/B359_internal_volume_constraints_set1.csv`

### Exact Work Order

1. define default data in `cta_defaults.json`
2. define CTA default component builder in `cta.py`
3. move CTA-specific laws into `cta_laws.py`
4. define CTA constraint loader and reference frame in `cta_constraints.py`
5. copy static assets into `assets/cta/`

### Rules

- `cta.py` may assemble the pieces
- `cta_laws.py` may contain special-case laws
- generic geometry files must remain CTA-free

### Validation

Instantiate CTA entirely from `MADS` assets and case helpers.

### Exit Criteria

- CTA can be built without importing `parametrization/CTA`

### Stop Condition

Do not proceed if CTA setup still requires direct access to the old repo.

---

## Phase 11: Examples And Regression Checks

### Goal

Create reproducible entrypoints and validation scripts.

### Files To Fill

1. `MADS/examples/cta_geometry_only.py`
2. `MADS/examples/cta_geometry_constraints.py`
3. `MADS/examples/cta_geometry_export.py`

### Exact Work Order

1. geometry-only example
2. geometry-plus-constraints example
3. geometry-plus-export example

### Validation Outputs

The examples must be able to reproduce:

- geometry build
- profile comparison
- envelope comparison
- packaging margins
- IGES export

### Exit Criteria

- every major feature has one runnable example

### Stop Condition

Do not proceed to decommissioning until examples are stable.

---

## Phase 12: Remove Runtime Dependence On `BWB/CTA`

### Goal

Ensure the new `MADS` implementation is self-contained.

### Checklist

- remove temporary imports from `parametrization/bwb`
- remove temporary imports from `parametrization/CTA`
- ensure all assets are local to `MADS`
- ensure examples run without the old package

### Exit Criteria

- `MADS` is the sole runtime owner of geometry, packaging, and export

### Stop Condition

Do not declare the migration complete while any runtime path still points back
to `BWB` or `CTA`.

---

## File Creation Order Summary

If we need the shortest strict sequence, it is this:

1. `src/multiads/geometry/__init__.py`
2. `src/multiads/solvers/geometry/__init__.py`
3. `src/multiads/cases/__init__.py`
4. `src/multiads/assembly/parametric_wing.py`
5. `src/multiads/geometry/model.py`
6. `src/multiads/geometry/state.py`
7. `src/multiads/geometry/metrics.py`
8. `src/multiads/geometry/cst_profiles.py`
9. `src/multiads/geometry/spanwise_laws.py`
10. `src/multiads/geometry/sections.py`
11. `src/multiads/geometry/builder.py`
12. `src/multiads/solvers/geometry/cst_geometry.py`
13. `src/multiads/geometry/internal_volume_constraints.py`
14. `src/multiads/solvers/geometry/internal_volume_constraints.py`
15. `src/multiads/geometry/pygeo_export.py`
16. `src/multiads/solvers/geometry/pygeo_export.py`
17. `src/multiads/disciplines/geometry.py`
18. `src/multiads/disciplines/packaging.py`
19. `src/multiads/cases/cta.py`
20. `src/multiads/cases/cta_laws.py`
21. `src/multiads/cases/cta_constraints.py`
22. `assets/cta/*`
23. `examples/cta_geometry_only.py`
24. `examples/cta_geometry_constraints.py`
25. `examples/cta_geometry_export.py`

## First Real Coding Step

When implementation starts, the first real coding step should be:

1. fill `src/multiads/geometry/model.py`
2. fill `src/multiads/assembly/parametric_wing.py`
3. fill `src/multiads/geometry/state.py`

Only after those three are stable should solver implementation begin.
