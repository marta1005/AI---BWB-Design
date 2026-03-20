# Aircraft Parametrization

`parametrization/aircraft` is the generic conventional-aircraft geometry
package. It depends on `parametrization/shared` for reusable CST and dependency
helpers, but it does not depend on `parametrization/bwb` or
`parametrization/CTA`.

This package is for section-based wings, tails, fuselages, and simple aircraft
assembly. CTA itself should be modeled in `parametrization/bwb`, because CTA is
a BWB and needs centerbody and nose logic that does not belong here.

## Current Scope

The package currently contains:

- `profiles.py`
  Airfoil profile definitions based on CST, iCST, and an intuitive iCST wrapper.
- `laws.py`
  Spanwise interpolation and scalar-law definitions.
- `sections.py`
  Generic section placement data structures.
- `lifting_surface.py`
  Generic lifting-surface preparation, section sampling, and IGES export.
- `wing.py`
  High-level wing API built from section definitions.
- `wing_designer.py`
  Higher-level wing designer API built on top of `WingSpec`.
- `vertical_tail.py`
  High-level vertical-tail wrapper built on the same lifting-surface engine.
- `fuselage.py`
  Section-based fuselage definition and export.
- `aircraft.py`
  Lightweight whole-aircraft assembly of wings, tails, and fuselages.
- `plotting.py`
  2D and 3D plotting utilities for the prepared geometries.

## What You Can Define Right Now

### Wings

The wing workflow is the most mature part of the package.

There are two wing APIs:

- `WingDesignerSpec`
  Recommended user-facing API for aircraft design work.
- `WingSpec`
  Lower-level API when you want direct control of the resolved section geometry.

Right now you can define:

- an arbitrary number of physical wing sections
- a different airfoil profile at each section
- local chord, twist, pitch, roll, and thickness scale per section
- leading-edge position section by section
- wing vertical position section by section
- either absolute section placement or incremental sweep/dihedral placement
- hard transitions, local `C1` transitions, or local `C2` transitions between sections
- a globally curved spanwise spline only when you explicitly request it

The intended workflow is section-based, not point-by-point reconstruction.

Use:

- `WingDesignerSpec`
  Define a wing from global quantities such as `span`, `root_chord`, root
  placement, section chord ratios, quarter-chord sweep, and dihedral.
- `WingStationDesignerSpec`
  Define designer-facing stations such as `root`, `kink`, and `tip`.
- `WingTransitionSpec`
  Define the transition behavior between sections.
- `WingSpanStationSpec`
  Define sections with physical `span_position`, `x_le`, `vertical_y`, `chord`,
  `twist_deg`, `pitch_deg`, `roll_deg`, and `profile_id`.
- `WingSpec.from_span_sections(...)`
  Build a low-level wing from those physical sections.
- `WingStationSpec`
  Use this lower-level form if you prefer normalized `eta` stations and direct
  absolute geometry control.

### Airfoil Profiles

For wing and tail sections you can currently define:

- direct CST profiles with `CSTAirfoilProfileSpec`
- iCST-constrained profiles with `ICSTAirfoilProfileSpec`
- intuitive iCST profiles with `IntuitiveAirfoilProfileSpec`

The recommended user-facing option is `IntuitiveAirfoilProfileSpec`.

Typical profile inputs are:

- `leading_edge_radius`
- `max_thickness`
- `x_tmax`
- `max_camber`
- `x_cmax`
- `trailing_edge_wedge_angle_deg`
- `trailing_edge_camber_angle_deg`
- `aft_control_x`
- `te_thickness`

### Spanwise Transitions

Spanwise behavior is controlled with `InterpolationSpec`.

Supported methods:

- `InterpolationMethod.LINEAR`
  Hard piecewise-linear transitions.
- `InterpolationMethod.SEGMENTED`
  Straight spanwise segments with only a local smooth transition around the
  interior stations.
- `InterpolationMethod.PYSPLINE`
  Globally curved interpolation.

For conventional wings, the recommended default is:

```python
InterpolationSpec(
    method=InterpolationMethod.SEGMENTED,
    continuity=ContinuityOrder.C2,
    blend_fraction=0.18,
)
```

That means:

- the wing stays straight between sections
- only the neighborhood of the section junction is blended
- the user chooses whether the local blend is `C1` or `C2`

### Scalar Laws

You can override spanwise laws explicitly with `ScalarLawSpec`.

Supported scalar-law names for lifting surfaces are:

- `x`
- `y`
- `z`
- `chord`
- `roll_deg`
- `pitch_deg`
- `twist_deg`
- `thickness_scale`
- `te_thickness`

### Vertical Tails

`VerticalTailSpec` is available and uses the same lifting-surface machinery.

At the moment it is best understood as a wrapper around the generic lifting
surface builder, not a fully separate vertical-tail design language.

### Fuselages

`FuselageSpec` currently supports section-based fuselage definition using
closed superellipse-like sections.

Right now you can define:

- section position along the body
- width and height of each section
- section center offsets
- top, bottom, and side shape exponents

This is good for generic transport-style fuselages, but it is still a simple
body model.

### Aircraft Assembly

`AircraftAssemblySpec` can assemble:

- wings
- vertical tails
- fuselages

This is currently a geometry-organization and plotting/export layer. It is not
yet a watertight multi-body CAD assembly with boolean trimming or true fairings.

### Outputs

The package can currently generate:

- section `.dat` files
- lifting-surface `.igs` exports
- fuselage `.igs` exports
- 2D overview figures
- 3D figures

## Recommended Wing Workflow

1. Define a profile catalog.
2. Define a compact set of designer stations such as `root`, `kink`, and `tip`.
3. Choose a transition mode.
4. Convert the designer wing into a `WingSpec`.
5. Export `.igs` or generate plots.

### Minimal Wing Example

```python
from parametrization.aircraft import (
    ContinuityOrder,
    IntuitiveAirfoilProfileSpec,
    ProfileCatalog,
    WingDesignerSpec,
    WingStationDesignerSpec,
    WingTransitionSpec,
)

profiles = ProfileCatalog(
    (
        IntuitiveAirfoilProfileSpec(
            profile_id="root_airfoil",
            leading_edge_radius=0.018,
            max_thickness=0.145,
            x_tmax=0.34,
            max_camber=0.024,
            x_cmax=0.44,
        ),
        IntuitiveAirfoilProfileSpec(
            profile_id="tip_airfoil",
            leading_edge_radius=0.008,
            max_thickness=0.090,
            x_tmax=0.29,
            max_camber=0.008,
            x_cmax=0.36,
        ),
    )
)

segmented_c2 = WingTransitionSpec(
    continuity=ContinuityOrder.C2,
    blend_fraction=0.18,
)

wing = WingDesignerSpec(
    wing_id="main_wing",
    span=29.0,
    root_chord=6.2,
    symmetric=True,
    root_le_x=8.35,
    root_vertical_y=0.10,
    stations=(
        WingStationDesignerSpec(
            station_id="root",
            profile_id="root_airfoil",
            eta=0.0,
            chord_ratio=1.0,
            twist_deg=2.5,
        ),
        WingStationDesignerSpec(
            station_id="tip",
            profile_id="tip_airfoil",
            eta=1.0,
            chord_ratio=1.55 / 6.2,
            twist_deg=-2.4,
            sweep_qc_deg=25.9,
            dihedral_deg=5.5,
        ),
    ),
    section_transition=segmented_c2,
    spine_transition=segmented_c2,
)

wing_spec = wing.to_wing_spec()
```

If `symmetric=True`, `span` is interpreted as the full wingspan. Internally,
`to_wing_spec()` still produces the half-wing representation used by the
current wing builder. Use `to_full_span_component_spec()` when you want a full
mirrored wing component directly. If `symmetric=False`, `span` is interpreted
as the physical span of the single lifting surface you are defining.

## Wing Parameters You Can Define Right Now

### Per Designer Wing

- `wing_id`
- `span`
- `symmetric`
- `root_chord`
- `root_le_x`
- `root_vertical_y`
- `root_z`
- `spine_transition`
- `section_transition`
- `scalar_laws`
- `metadata`

### Per Designer Station

- `station_id`
- `eta`
- `profile_id`
- `chord_ratio`
- `twist_deg`
- `pitch_deg`
- `roll_deg`
- `thickness_scale`
- `sweep_qc_deg`
- `dihedral_deg`

### Per Profile

- `profile_id`
- `leading_edge_radius`
- `max_thickness`
- `x_tmax`
- `max_camber`
- `x_cmax`
- `trailing_edge_wedge_angle_deg`
- `trailing_edge_camber_angle_deg`
- `aft_control_x`
- `te_thickness`

### Per Wing Section

- `station_id`
- `span_position` or `eta`
- `profile_id`
- `chord`
- `x_le`
- `vertical_y`
- `twist_deg`
- `pitch_deg`
- `roll_deg`
- `thickness_scale`

If you use `WingStationSpec`, you may also define:

- `sweep_le_deg`
- `dihedral_deg`

### Per Wing

- `wing_id`
- `semispan`
- `root_x`
- `root_y`
- `root_z`
- `section_interpolation`
- `spine_interpolation`
- `scalar_laws`
- `metadata`

## What Still Needs Improvement

The package is useful, but it is not finished.

### Wing-Level Improvements

- a sizing-oriented wing API based on area, aspect ratio, taper ratio,
  and quarter-chord sweep targets
- better handling of mirrored full-wing construction as a first-class concept
- more explicit support for multiple breaks, flaps, and control-surface regions
- stronger validation of geometric feasibility between neighboring sections
- more direct control of trailing-edge closure and manufacturability constraints

### Fuselage Improvements

- richer nose and tail shaping parameters
- better cockpit and radome control
- section families beyond the current superellipse-like body model

### Aircraft Assembly Improvements

- true fairings
- component intersection handling
- watertight export of multi-component assemblies
- better placement rules for horizontal tails and pylons/nacelles

## Out Of Scope For This Package

The following should not be implemented in `parametrization/aircraft`:

- CTA-specific geometry reconstruction
- BWB centerbody and nose logic
- CTA-specific reference matching

Those belong in `parametrization/bwb` and `parametrization/CTA`.

## Examples

The generic examples that remain relevant are:

- `python parametrization/aircraft/examples/run_demo_wing_iges.py`
- `python parametrization/aircraft/examples/plot_demo_wing_views.py`
- `python parametrization/aircraft/examples/plot_demo_wing_transition_comparison.py`
- `python parametrization/aircraft/examples/show_demo_wing_3d.py`

Wing examples are kept intentionally lean for now. Fuselage and full-aircraft
examples have been removed until those APIs are refined further.

## Summary

Today, `aircraft` is strongest as a section-based wing and lifting-surface
package, plus a simple fuselage and assembly layer.

If the task is:

- generic wing or tail parametrization: use `aircraft`
- generic fuselage sections: use `aircraft`
- BWB centerbody, nose, or CTA-specific geometry: use `bwb` / `CTA`
