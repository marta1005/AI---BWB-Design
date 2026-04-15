# BWB Core Guide

`parametrization/bwb` is the generic geometry core for sectioned blended-wing-
body wings.

Its job is not to define one specific aircraft.
Its job is to provide the machinery needed to declare a BWB wing, reconstruct
its planform and profiles, and prepare a sampled geometry that can then be
plotted, validated, or exported.

This package is meant to be reused by case-specific layers such as CTA.

## What The Core Does

At a high level, the core takes:

- a topology
- a planform specification
- a section/profile family
- spanwise laws
- sampling/export settings

and turns them into a `SectionedBWBModelConfig`, then into a prepared geometry.

The prepared geometry contains:

- continuous LE and TE definitions
- reconstructed profiles at any spanwise station
- spanwise twist/camber/vertical laws
- sampled stations and loft-ready data

## Main Modules

- `topology.py`
  Spanwise topology and station placement.

- `specs.py`
  Main dataclasses used to define a BWB case.

- `case_definition.py`
  Utilities to declare cases, solve relations, and build configs from either
  compact variables or explicit sections.

- `design_variables.py`
  Compact legacy design vector for four main sections.

- `design_space.py`
  Generic design-space helpers: bounds, metadata, flattening and sampling.

- `planform.py`
  Reconstruction of LE/TE curves from sections, helper points, and blend rules.

- `sections.py`
  Reconstruction of the airfoil family along the span.

- `spanwise_laws.py`
  Twist, camber, and vertical-offset interpolation.

- `builder.py`
  End-to-end geometry preparation.

- `validation.py`
  Geometry checks.

- `volume.py`
  Volume-related evaluation helpers.

- `exporters.py`
  Export helpers such as IGES/pyGeo-facing sampling.

## Mental Model

The core is easiest to understand if you split it into five blocks:

1. Topology
   Where are the main sections and profile anchors along the span?

2. Planform
   What are the section chords, LE positions, TE positions, and local shaping
   rules?

3. Profiles
   What airfoils exist at the anchor stations, and how are they interpolated?

4. Spanwise laws
   How do twist, camber delta, and vertical offset vary with span?

5. Sampling/export
   How densely do we sample the geometry for plots, lofts, or export?

## Two Main Ways To Define A Geometry

The core supports two complementary declaration styles.

### 1. Compact design-vector mode

This is the older, smaller, more parameterized mode.

It uses:

- `SectionedBWBDesignVariables`
- semispan ratios `b1_span_ratio`, `b2_span_ratio`, `b3_span_ratio`
- root chord + chord ratios
- sweep angles per main segment

This is useful when:

- you want a compact optimization vector
- your geometry still fits the “4 main sections / 3 main segments” pattern

### 2. Explicit section mode

This is the more general mode.

It uses:

- `SectionedBWBTopologySpec(section_y=...)`
- `PlanformSpec(section_le_x=..., section_chords_override=...)`
- optional explicit `leading_edge_control_points`
- optional explicit `trailing_edge_control_points`

This is useful when:

- you want to choose the number of sections explicitly
- you want to add one more section without changing the core
- you want to decide exactly how LE and TE are joined
- you want a case layer to be more declarative

In practice, this is the recommended mode for case-specific packages.

## Core Dataclasses

## `SectionedBWBTopologySpec`

Defined in [topology.py](/Users/martaarnabatmartin/Desktop/BWB/parametrization/bwb/topology.py).

It declares the spanwise structure.

Main fields:

- `span`
- `b1_span_ratio`
- `b2_span_ratio`
- `b3_span_ratio`
- `section_y`
- `anchor_y`
- `topology_name`

Two patterns are supported:

- compact mode:
  - set `span`
  - set `b1_span_ratio`, `b2_span_ratio`, `b3_span_ratio`
  - leave `section_y=None`

- explicit mode:
  - set `span`
  - set `section_y=(...)`
  - leave the ratios unused

Rules:

- section stations must be strictly increasing
- first station must be `0.0`
- last station must equal `span`
- anchor stations must also be strictly increasing
- anchors must lie inside `[0, span]`
- anchors must start at `0.0` and end at `span`

## `PlanformSpec`

Defined in [specs.py](/Users/martaarnabatmartin/Desktop/BWB/parametrization/bwb/specs.py).

This is the heart of the planform declaration.

In compact mode, the main fields are:

- `le_root_x`
- `c1_root_chord`
- `c2_c1_ratio`
- `c3_c1_ratio`
- `c4_c1_ratio`
- `s1_deg`
- `s2_deg`
- `s3_deg`

In explicit mode, the important fields are:

- `section_le_x`
- `section_chords_override`
- `leading_edge_control_points`
- `trailing_edge_control_points`

The shaping controls available in all modes include:

- `body_le_fixed_points`
- `med_3_te_sweep_deg`
- `med_3_te_helper_fraction`
- `te_exact_segments`
- `le_exact_segments`
- `te_spline_bridge`
- `le_spline_bridge`
- `blend_fraction`
- `te_blend_fraction`
- `min_linear_core_fraction`
- `te_min_linear_core_fraction`
- `le_linear_start_index`
- `symmetry_blend_y`

This means the core can express:

- exact linear segments
- blended transitions
- spline bridges
- TE helper points
- symmetry-side smoothing

## `SectionCSTSpec`

Defined in [specs.py](/Users/martaarnabatmartin/Desktop/BWB/parametrization/bwb/specs.py).

This declares one airfoil anchor section.

Fields:

- `upper_coeffs`
- `lower_coeffs`
- `tc_max`
- `x_tmax`
- `te_thickness`

These are bundled into a `SectionFamilySpec`.

## `SectionFamilySpec`

Defined in [specs.py](/Users/martaarnabatmartin/Desktop/BWB/parametrization/bwb/specs.py).

This declares the profile family used by the wing.

Main fields:

- `cst_degree`
- `n1`
- `n2`
- `shared_leading_edge`
- `profile_generation_mode`
- `section_specs_override`
- `profile_relations`

So this block answers:

- what the anchor airfoils are
- how CST is interpreted
- how one profile may reuse information from an earlier one

## `SpanwiseLawSpec`

Defined in [specs.py](/Users/martaarnabatmartin/Desktop/BWB/parametrization/bwb/specs.py).

This declares spanwise laws for:

- `twist_deg`
- `camber_delta`
- `vertical_offset_z`

In most case layers these are built from `AnchoredSpanwiseLaw`.

## `SamplingSpec`

Defined in [specs.py](/Users/martaarnabatmartin/Desktop/BWB/parametrization/bwb/specs.py).

This controls numerical resolution.

Important fields:

- `num_airfoil_points`
- `num_base_stations`
- `section_curve_n_ctl`
- `section_interpolation`
- `airfoil_distribution_mode`
- `k_span`

## `SectionedBWBModelConfig`

Defined in [specs.py](/Users/martaarnabatmartin/Desktop/BWB/parametrization/bwb/specs.py).

This is the top-level object that the builder expects.

It groups:

- `topology`
- `planform`
- `sections`
- `spanwise`
- `sampling`
- `export`
- `volume`

## Typical Build Flow

The standard pipeline is:

1. declare the case
2. build a `SectionedBWBModelConfig`
3. call `prepare_geometry(config)`

That happens in [builder.py](/Users/martaarnabatmartin/Desktop/BWB/parametrization/bwb/builder.py).

The result is a prepared geometry object with:

- `planform`
- `section_model`
- `spanwise_laws`
- `loft`

Those are what plotting/export code should usually consume.

## Case Utilities

The generic helpers in [case_definition.py](/Users/martaarnabatmartin/Desktop/BWB/parametrization/bwb/case_definition.py)
exist to make case-specific packages easier to write.

The key pieces are:

- `Relation`
  One derived relation between declared values.

- `CanonicalSectionDeclaration`
  A compact declaration for four main sections.

- `ExplicitSectionSpec`
  One explicitly declared section with:
  - `name`
  - `y_m`
  - `le_x_m`
  - `chord_m`

- `CaseTemplate`
  A reusable template that contains profile anchors, CST family, laws, and
  sampling settings.

- `solve_relations(...)`
  Resolves a graph of declared and derived values.

- `build_case_config_from_resolved(...)`
  Builds a config from resolved compact values.

- `build_case_config_from_explicit_sections(...)`
  Builds a config from explicit sections and optional explicit LE/TE controls.

## Example A: Compact Legacy Configuration

This is the smaller, ratio/sweep-driven path:

```python
from parametrization.bwb.design_variables import SectionedBWBDesignVariables
from parametrization.bwb.builder import prepare_geometry

design = SectionedBWBDesignVariables.reference_seed()
config = design.to_model_config()
prepared = prepare_geometry(config)
```

Use this when:

- the compact design vector is enough
- your case really fits the four-main-section pattern

## Example B: Explicit Sections

This is the more general path:

```python
from parametrization.bwb.case_definition import (
    ExplicitSectionSpec,
    build_case_config_from_explicit_sections,
)
from parametrization.bwb.builder import prepare_geometry

sections = (
    ExplicitSectionSpec("S0", 0.0, 0.0, 40.0),
    ExplicitSectionSpec("S1", 8.0, 16.0, 14.0),
    ExplicitSectionSpec("S2", 13.0, 22.0, 8.0),
    ExplicitSectionSpec("S3", 39.5, 36.0, 0.8),
)

config = build_case_config_from_explicit_sections(
    sections_definition=sections,
    template=template,
    topology_name="my_bwb_case",
)

prepared = prepare_geometry(config)
```

Use this when:

- you want to insert another section
- you want section positions to be declared directly
- you want LE/TE control to be explicit

## Example C: Explicit Sections With Control Points

You can go one step further and drive LE/TE through control points:

```python
config = build_case_config_from_explicit_sections(
    sections_definition=sections,
    template=template,
    topology_name="my_bwb_case",
    leading_edge_control_points=(
        (0.0, 0.0),
        (2.0, 2.0),
        (10.0, 6.0),
        (16.0, 8.0),
        (22.0, 13.0),
        (36.0, 39.5),
    ),
    trailing_edge_control_points=(
        (40.0, 0.0),
        (32.0, 6.0),
        (30.0, 8.0),
        (29.5, 13.0),
        (36.8, 39.5),
    ),
    te_exact_segments=(0, 3),
    te_spline_bridge=(1, 3),
)
```

This is the path to use when you want the case layer to control:

- where the main sections are
- how the LE is joined
- how the TE is joined
- which segments must remain exact
- where spline bridges must be applied

## What You Can Declare

At the generic BWB-core level, you can declare:

### Topology

- total semispan
- section stations
- anchor stations
- or legacy span ratios

### Planform

- section LE positions
- section chords
- root-chord + ratios
- LE sweeps
- TE helper rules
- LE/TE control points
- exact segments
- spline bridges
- blend fractions

### Profiles

- anchor CST sections
- CST family defaults (`n1`, `n2`, `shared_leading_edge`, generation mode)
- profile inheritance rules

### Spanwise laws

- twist by anchor
- camber delta by anchor
- vertical offset by anchor

### Sampling/export

- airfoil resolution
- section interpolation mode
- station count
- spline control resolution

## What The Core Does Not Assume

The core does not need to know:

- project-specific names like `C0`, `S1`, `Bw`
- whether a section is called “transition wing” or “outer wing”
- whether a given public parameter is absolute or ratio-based

Those are case-layer decisions.

So the right split is:

- `bwb` = geometry machinery
- case package = naming, relations, project-specific parameterization

## When To Use Each Layer

Use `bwb` directly when:

- you are building a new generic BWB case
- you want a new section layout
- you want explicit section/control-point declaration

Use a case package on top of `bwb` when:

- you need project-facing names
- you need public driver relations
- you need specific bounds/plots/exports for one aircraft case

## Recommended Pattern For New Cases

If you are building a new case on top of `bwb`, a good structure is:

1. one JSON or Python declaration file containing:
   - canonical sections
   - public drivers if needed
   - template

2. one `case.py` that:
   - loads the declaration
   - resolves relations
   - builds explicit sections and control points
   - returns a `SectionedBWBModelConfig`

3. plotting/export scripts that consume only the built config

That keeps the core generic and the case layer clean.
