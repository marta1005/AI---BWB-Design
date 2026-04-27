# CTA Guide

`parametrization/CTA` defines one concrete wing case on top of the generic
`parametrization/bwb` machinery.

This package is useful when we want two things at the same time:

- keep the generic BWB core reusable
- declare one specific case with its own geometry, naming, and outputs

The source of truth for the case is:

- `parametrization/CTA/data/cta_case_inputs.json`
- `parametrization/CTA/case.py`

The JSON declares the case.
The Python code explains how those declared values are turned into a
`SectionedBWBModelConfig`.

## Mental Model

The CTA case is split into three layers:

1. `canonical_sections`
   This is the direct geometric declaration of the wing planform by section.

2. `public_drivers`
   This is an alternative declaration using project-facing parameters such as
   `Bw`, `C0`, `C3`, `C4/C3`, `S1`, `S2`, etc.

3. `template`
   This contains everything that is not part of the main planform driver set:
   profile anchors, CST coefficients, twist anchors, LE z anchors, interpolation
   choices, helper points, blend rules, and sampling options.

The idea is:

- `canonical_sections` says what the geometry is
- `public_drivers` says how the same geometry may be driven
- `template` says how to loft, interpolate, and sample it

## Package Layout

- `parametrization/CTA/case.py`
  Main case-definition logic.

- `parametrization/CTA/data/cta_case_inputs.json`
  Declarative input file for the CTA case.

- `parametrization/CTA/design_space.py`
  Optional public-variable layer used for bounds, GIFs, variation plots, and
  public sampling.

- `parametrization/CTA/codes/plotting`
  Plotting and animation scripts for the active CTA wing case.

- `parametrization/CTA/codes/exports`
  IGES and bounds exports.

## Step 1: Declare The Canonical Geometry

The first block in the JSON is `canonical_sections`.

It declares the wing through four main planform stations:

```json
"canonical_sections": {
  "y_0_m": 0.0,
  "y_1_m": 8.041,
  "y_2_m": 12.5081007083,
  "y_3_m": 39.4995,
  "le_x_0_m": 0.0,
  "le_x_1_m": 15.785947611133833,
  "le_x_2_m": 21.946845426354557,
  "le_x_3_m": 36.103616972154164,
  "chord_0_m": 41.17952274,
  "chord_1_m": 13.9269627,
  "chord_2_m": 7.76845979406,
  "chord_3_m": 0.8
}
```

These are the fields that can be declared here:

- `y_0_m`, `y_1_m`, `y_2_m`, `y_3_m`
- `le_x_0_m`, `le_x_1_m`, `le_x_2_m`, `le_x_3_m`
- `chord_0_m`, `chord_1_m`, `chord_2_m`, `chord_3_m`

Interpretation:

- `y_i_m` is the spanwise station position
- `le_x_i_m` is the leading-edge x position of that station
- `chord_i_m` is the local chord at that station

So this block is the direct declaration of:

- section positions
- leading-edge positions
- section chords

From these values, the code derives:

- semispan
- B1/B2/B3 span ratios
- LE sweeps
- internal chord ratios

The code path is:

- `load_cta_canonical_declaration(...)`
- `resolve_cta_from_sections(...)`

## Step 2: Declare The Public Driver View

The second block is `public_drivers`.

This is not the direct geometry.
This is the project-facing parameterization used when we want to talk in terms
of `Bw`, `C0`, `C3`, `C4/C3`, `S1`, `S2`, etc.

Example:

```json
"public_drivers": {
  "wing_span_m": 31.458499999999997,
  "b2_over_wing_span": 0.1419998,
  "c0_m": 41.17952274,
  "c3_m": 13.9269627,
  "wing_med_3_tr": 0.5578,
  "c5_m": 0.8,
  "s1_50_deg": 34.6,
  "s2_25_deg": 24.7,
  "med_3_te_sweep_deg": 0.0
}
```

These are the public fields that can be declared here:

- `wing_span_m`
- `b2_over_wing_span`
- `c0_m`
- `c3_m`
- `wing_med_3_tr`
- `c5_m`
- `s1_50_deg`
- `s2_25_deg`
- `med_3_te_sweep_deg`

Interpretation:

- `wing_span_m` = `B2 + B3`
- `b2_over_wing_span` = `B2 / (B2 + B3)`
- `c0_m` = `C0`
- `c3_m` = `C3`
- `wing_med_3_tr` = `C4 / C3`
- `c5_m` = `C5`
- `s1_50_deg` = transition-wing sweep measured at 50% chord
- `s2_25_deg` = outer-wing sweep measured at 25% chord
- `med_3_te_sweep_deg` = explicit TE shaping on the `C3 -> C4` segment

From these values, the code derives:

- `y_1_m`, `y_2_m`, `y_3_m`
- `span_m`
- `b1_span_ratio`, `b2_span_ratio`, `b3_span_ratio`
- section chords
- internal LE sweeps `s1_deg`, `s2_deg`, `s3_deg`
- `le_x_1_m`, `le_x_2_m`, `le_x_3_m`

The code path is:

- `load_cta_public_declaration(...)`
- `resolve_cta_from_public(...)`

## Step 3: Declare The Template

The `template` block contains everything that tells the core how to build the
final 3D wing once the main section geometry is known.

It includes four groups of data.

### 3.1 Profile anchor stations

```json
"anchor_y_m": [0.0, 1.9, 5.694, 8.041, 12.5081007083, 39.4995],
"anchor_le_z_m": [0.25865, 1.03474, 0.48145, 0.67693, 0.72992, 1.97538],
"anchor_twist_deg": [0.778, -0.342, 0.371, -0.249, 0.483, 3.177],
"anchor_camber_delta": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
```

These define:

- where profile anchors sit along the span
- the vertical LE offset at each anchor
- the twist at each anchor
- the camber delta at each anchor

### 3.2 CST section definitions

The `anchor_sections` list defines the actual fitted profiles.

Each entry declares:

- `upper_coeffs`
- `lower_coeffs`
- `tc_max`
- `x_tmax`
- `te_thickness`

So this is the aerodynamic profile family of the case.

### 3.3 Profile generation defaults

These fields define how the CST family is interpreted:

- `cst_n1`
- `cst_n2`
- `shared_leading_edge`
- `profile_generation_mode`

For the current CTA case, the defaults are:

- `cst_n1 = 0.35`
- `cst_n2 = 1.1`
- `shared_leading_edge = false`
- `profile_generation_mode = "cst_only"`

### 3.4 Planform helper and blending rules

These fields control how LE and TE are reconstructed between the declared
stations:

- `body_le_fixed_points`
- `te_exact_segments`
- `te_c1_span_fraction`
- `te_inboard_blend_fraction`
- `te_inboard_blend_dx`
- `te_inboard_radius_factor`
- `med_3_te_helper_fraction`
- `te_outer_blend_fraction`
- `blend_fraction`
- `min_linear_core_fraction`
- `te_blend_fraction`
- `te_min_linear_core_fraction`
- `te_spline_bridge`
- `symmetry_blend_y`

This is the part that lets the case describe:

- exact segments
- local TE shaping
- helper points
- symmetry-side blending
- spline bridge zones

### 3.5 Sampling and interpolation defaults

These fields control the generated geometry resolution:

- `section_interpolation`
- `airfoil_distribution_mode`
- `num_airfoil_points`
- `num_base_stations`
- `section_curve_n_ctl`
- `k_span`

For the current case, the main defaults are:

- `section_interpolation = "pchip"`
- `airfoil_distribution_mode = "anchors"`
- `num_airfoil_points = 241`
- `num_base_stations = 41`
- `section_curve_n_ctl = 41`
- `k_span = 4`

### 3.6 Current spanwise laws

The current CTA case uses these laws between anchors:

| Quantity | Current law | Where it is declared | Affects IGES? |
| --- | --- | --- | --- |
| Profile interpolation between anchor sections | `pchip` | `template.section_interpolation` in `cta_case_inputs.json` | Yes |
| Twist law | `pchip` | `config.spanwise.twist_deg` in `case.py` | Yes |
| Vertical `LE_z` / `vertical_offset_z` law | `pchip` | `config.spanwise.vertical_offset_z` in `case.py` | Yes |
| Camber-delta law | `pyspline` | `config.spanwise.camber_delta` in `case.py` | Yes |
| Plot-only front-view smoothing | whatever is coded in `plot_cta_views.py` | plotting only | No |

This distinction is important:

- anything declared in `case.py` / `cta_case_inputs.json` feeds the real geometry
- anything added only in `plot_cta_views.py` changes the figure, not the IGES

## What Is Fixed, What Is Derived

Within the current CTA setup:

- the case JSON fixes the current CTA instance
- `canonical_sections` is the direct geometric declaration
- `public_drivers` is an alternative driver declaration of the same case
- the template fixes the profile family and the interpolation/loft rules

Derived internally from the canonical/public data:

- `span_m`
- `b1_span_ratio`
- `b2_span_ratio`
- `b3_span_ratio`
- `c2_c1_ratio`
- `c3_c1_ratio`
- `c4_c1_ratio`
- `s1_deg`
- `s2_deg`
- `s3_deg`
- helper span fractions tied to fixed absolute positions

## Main Build Paths

There are two main ways to build the CTA case.

### Path A: build from canonical sections

Use this when the section geometry is the source of truth.

```python
from parametrization.CTA.case import (
    load_cta_case_payload,
    load_cta_canonical_declaration,
    load_cta_case_template,
    resolve_cta_from_sections,
    build_cta_case_config_from_resolved,
)

payload = load_cta_case_payload()
canonical = load_cta_canonical_declaration(payload)
template = load_cta_case_template(payload, canonical)
resolved = resolve_cta_from_sections(canonical, template)
config = build_cta_case_config_from_resolved(resolved, template)
```

### Path B: build from public drivers

Use this when `Bw`, `C0`, `C3`, `C4/C3`, `S1`, `S2`, etc. are the source of
truth.

```python
from parametrization.CTA.case import (
    load_cta_case_payload,
    load_cta_canonical_declaration,
    load_cta_public_declaration,
    load_cta_case_template,
    resolve_cta_from_sections,
    resolve_cta_from_public,
    build_cta_case_config_from_resolved,
)

payload = load_cta_case_payload()
canonical = load_cta_canonical_declaration(payload)
public = load_cta_public_declaration(payload)
template = load_cta_case_template(payload, canonical)

fixed_internal = resolve_cta_from_sections(canonical, template)
resolved = resolve_cta_from_public(
    public,
    canonical,
    template,
    fixed_internal_s1_deg=fixed_internal["s1_deg"],
)
config = build_cta_case_config_from_resolved(resolved, template)
```

### Path C: build the current declared target directly

If you just want the current declared CTA case:

```python
from parametrization.CTA.case import build_cta_case_target_config

config = build_cta_case_target_config()
```

This is the shortest and cleanest entrypoint when you want the current active
case exactly as declared.

## How The Final Config Is Built

The final object is a `SectionedBWBModelConfig`.

The chain is:

1. load case payload
2. resolve either canonical or public declaration
3. create explicit sections and control points
4. build `SectionedBWBModelConfig`
5. pass it to `prepare_geometry(...)`

The explicit-section helpers in `case.py` are:

- `build_cta_explicit_sections_from_resolved(...)`
- `build_cta_leading_edge_control_points(...)`
- `build_cta_trailing_edge_control_points(...)`

These are important because they are the bridge between:

- high-level CTA declarations
- and the generic `bwb` geometry core

## If You Want To Use pyGeo

When we export CTA as IGES, `pyGeo` does **not** start from the plotting views.
It starts from the resolved BWB geometry.

The actual export path is:

1. `build_cta_design()`
2. `to_cta_model_config(...)`
3. `prepare_geometry(config)`
4. `build_surface(config)`
5. `pyGeo("liftingSurface", ...)`

So the flow is:

```text
CTA JSON/case.py
    -> resolved CTA design
    -> SectionedBWBModelConfig
    -> prepared loft definition
    -> pyGeo lifting-surface loft
    -> IGES
```

The object passed to `pyGeo` is not a front-view sketch.
It is a spanwise set of sections plus their placement laws.

Concretely, `pyGeo` receives:

- `xsections`
  Airfoil files written at anchor/sample stations.
- `scale`
  Local chord at each station.
- `x`
  Leading-edge x position at each station.
- `y`
  Vertical position at each station.
- `z`
  Spanwise position at each station.
- `rotZ`
  Local twist at each station.
- `offset`
  Local section offset.
- `teHeightScaled`
  Local trailing-edge thickness.
- `nCtl`, `kSpan`
  Loft-control settings.

For CTA, the export script is:

- `parametrization/CTA/codes/exports/run_cta_iges.py`

and the core handoff to `pyGeo` happens in:

- `parametrization/bwb/exporters.py`

## What This Means For Transition Control

This is the key point:

- `pyGeo` gives us a continuous loft
- but `pyGeo` does **not** decide what the transition should be

The transition shape is already largely determined **before** the export, by:

- the chosen profile anchors
- the section interpolation mode
- the LE/TE control points and blend settings
- the twist law
- the vertical offset law
- the spanwise station sampling

So, yes, we do have control over the transitions, but it is **indirect**.
We do **not** usually control them as one final explicit hand-drawn curve.
Instead, we control the ingredients that define the loft.

In practice:

- if you change something only in `plot_cta_views.py`, that affects the figure
  only
- if you want the IGES to change, the change must happen in `case.py` / the CTA
  config / the BWB build inputs

So, for the current CTA case:

- the loft uses `pchip` for profile interpolation
- the twist law uses `pchip`
- the vertical law uses `pchip`
- the camber delta uses `pyspline`
- and `airfoil_distribution_mode = "anchors"` means the exported `.dat` files
  are the anchor sections only, while `pyGeo` interpolates the missing spanwise
  stations during loft construction

That is why a front-view tweak can look right in a plot and still not belong to
the exported geometry.

## What You Can Change Safely

If you want to edit the current case without changing the structure:

- change `canonical_sections` to redefine the planform directly
- change `public_drivers` to redefine the case via public parameters
- change `anchor_*` arrays to modify the spanwise profile/twist/z laws
- change `anchor_sections` to update the airfoil family
- change template blend/helper values to alter the LE/TE behavior
- change sampling values to alter discretization only

## What You Should Keep Consistent

Some arrays are expected to stay aligned:

- `anchor_y_m`
- `anchor_le_z_m`
- `anchor_twist_deg`
- `anchor_camber_delta`
- `anchor_sections`

They all describe the same set of profile anchors, so their lengths and ordering
must match.

Also:

- `y_0_m < y_1_m < y_2_m < y_3_m`
- all `chord_i_m` should stay positive
- `wing_med_3_tr` should stay positive
- spanwise fractions should remain physically valid

## Role Of `design_space.py`

`case.py` is enough to define and build the CTA case.

`design_space.py` is optional.
It exists only when we want a public-variable layer for:

- bounds tables
- animation sampling
- parameter-variation plots
- frame traces

So the split is:

- `case.py` = case declaration and geometry build
- `design_space.py` = public exploration layer

## Current Useful Outputs

The main wing outputs are written to:

- `parametrization/CTA/outputs/wing`

Useful files there include:

- `cta_layout.png`
- `cta_front.png`
- `cta_front_half_le_te.png`
- `cta_profiles.png`
- `cta_definition_scheme.png`
- `cta_3d.png`
- `cta_planform_parameter_variations_grid.png`
- `cta_planform_parameter_sweep.gif`
- `cta_planform_parameter_frame_traces.png`
- `cta_planform_parameter_frame_traces.csv`

## Typical Commands

Build wing plots:

```bash
./.venv/bin/python parametrization/CTA/codes/plotting/plot_cta_views.py
./.venv/bin/python parametrization/CTA/codes/plotting/plot_cta_planform_parameter_variations.py
./.venv/bin/python parametrization/CTA/codes/plotting/animate_cta_planform_parameter_sweep.py
./.venv/bin/python parametrization/CTA/codes/plotting/plot_cta_vs_glider_vertical_overlay.py
```

Export:

```bash
./.venv/bin/python parametrization/CTA/codes/exports/export_cta_bounds_table.py
./.venv/bin/python parametrization/CTA/codes/exports/run_cta_iges.py
```
