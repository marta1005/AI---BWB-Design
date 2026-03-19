# Wing Parametrization

This folder contains the generic aircraft geometry package. At this stage, the
wing workflow is the most mature part of the package, and this document focuses
on how a wing is defined, parameterized, built, and exported.

The core design goal is:

- keep the wing definition intuitive for aircraft design work
- allow different wing layouts, not only one CTA-like geometry
- keep spanwise segments constructively straight when needed
- make the transition around interior stations smooth with either `C1` or `C2`
- avoid a global spline bending the whole wing unless the user explicitly wants
  that behavior

## Geometry Convention

Global axes:

- `x`: streamwise / longitudinal
- `y`: vertical
- `z`: spanwise / lateral

Wing span parameter:

- `eta = 0.0` at the wing root
- `eta = 1.0` at the wing tip

The current high-level wing API defines one physical semi-span. A mirrored full
wing can be assembled later at aircraft level.

## Wing Definition Stack

A wing is built in five layers:

1. Define airfoil profiles.
2. Define wing stations.
3. Define spanwise interpolation behavior.
4. Optionally define scalar spanwise laws.
5. Build the lifting surface and export geometry or plots.

The main classes involved are:

- `IntuitiveAirfoilProfileSpec`
- `ProfileCatalog`
- `WingStationSpec`
- `WingSpec`
- `InterpolationSpec`
- `ScalarLawSpec`
- `LiftingSurfaceBuildOptions`

## 1. Airfoil Profiles

Wing sections are identified by `profile_id`. The recommended profile input is
`IntuitiveAirfoilProfileSpec`, which internally resolves to the iCST-based CST
representation used by the geometry engine.

### Recommended Profile Parameters

Required:

- `profile_id`
- `leading_edge_radius`
- `max_thickness`
- `x_tmax`

Common optional parameters:

- `max_camber`
- `x_cmax`
- `trailing_edge_wedge_angle_deg`
- `trailing_edge_camber_angle_deg`
- `aft_control_x`
- `te_thickness`

### Example

```python
IntuitiveAirfoilProfileSpec(
    profile_id="root_airfoil",
    leading_edge_radius=0.018,
    max_thickness=0.145,
    x_tmax=0.34,
    max_camber=0.024,
    x_cmax=0.44,
    trailing_edge_wedge_angle_deg=13.5,
    trailing_edge_camber_angle_deg=-0.8,
    aft_control_x=0.74,
    te_thickness=0.0018,
)
```

All profiles used by the wing are stored in a `ProfileCatalog`.

## 2. Wing Stations

`WingStationSpec` defines each anchor section of the wing.

Each station describes:

- where the section is along the span
- which profile it uses
- the local chord
- the local orientation
- and optionally the local placement rule

### Required `WingStationSpec` Parameters

| Parameter | Meaning |
| --- | --- |
| `station_id` | Unique name of the station |
| `eta` | Normalized span location in `[0, 1]` |
| `profile_id` | Airfoil profile used at the station |
| `chord` | Local chord length |

### Common Optional `WingStationSpec` Parameters

| Parameter | Meaning |
| --- | --- |
| `twist_deg` | Local geometric twist |
| `pitch_deg` | Additional pitch rotation |
| `roll_deg` | Additional roll rotation |
| `thickness_scale` | Multiplicative thickness scale |
| `x_le` | Absolute leading-edge `x` position |
| `vertical_y` | Absolute vertical `y` position |
| `sweep_le_deg` | Segment sweep from previous station |
| `dihedral_deg` | Segment dihedral from previous station |

### Placement Rule

For each non-root station, use one of these approaches:

- absolute placement:
  - `x_le`
  - `vertical_y`
- segment-based placement:
  - `sweep_le_deg`
  - `dihedral_deg`

Do not mix the absolute and segment-based versions of the same quantity on the
same station.

### Station Rules

- the first station must be at `eta = 0.0`
- the last station must be at `eta = 1.0`
- the root station position is controlled by `WingSpec.root_x`,
  `WingSpec.root_y`, and `WingSpec.root_z`
- wing stations must be strictly increasing in `eta`

### Example

```python
stations = (
    WingStationSpec(
        station_id="root",
        eta=0.0,
        profile_id="root_airfoil",
        chord=6.2,
        twist_deg=2.5,
    ),
    WingStationSpec(
        station_id="kink",
        eta=0.38,
        profile_id="kink_airfoil",
        chord=4.0,
        twist_deg=0.8,
        sweep_le_deg=20.0,
        dihedral_deg=3.5,
    ),
    WingStationSpec(
        station_id="tip",
        eta=1.0,
        profile_id="tip_airfoil",
        chord=1.55,
        twist_deg=-2.4,
        sweep_le_deg=29.0,
        dihedral_deg=5.5,
    ),
)
```

## 3. `WingSpec`: The High-Level Wing Object

`WingSpec` gathers the whole wing definition.

### `WingSpec` Parameters

| Parameter | Meaning |
| --- | --- |
| `wing_id` | Wing identifier |
| `semispan` | Semi-span length |
| `stations` | Tuple of `WingStationSpec` anchor sections |
| `root_x` | Root leading-edge `x` reference |
| `root_y` | Root vertical position |
| `root_z` | Root spanwise position |
| `section_interpolation` | How the airfoil family changes between stations |
| `spine_interpolation` | How planform-like quantities evolve along the span |
| `scalar_laws` | Optional user-defined spanwise laws |
| `mirrored` | Reserved field for future symmetry handling |
| `metadata` | Free-form descriptive metadata |

### What `spine_interpolation` Controls

`spine_interpolation` drives the spanwise evolution of:

- `x`
- `y`
- `z`
- `chord`
- `twist_deg`
- `pitch_deg`
- `roll_deg`
- `thickness_scale`

### What `section_interpolation` Controls

`section_interpolation` drives the transition of the airfoil family itself,
which means the resolved iCST/CST profile coefficients and trailing-edge
thickness values.

## 4. Interpolation Modes

Interpolation is defined with `InterpolationSpec`.

### `InterpolationSpec` Parameters

| Parameter | Meaning |
| --- | --- |
| `method` | Interpolation method |
| `continuity` | Requested continuity order |
| `blend_fraction` | Local transition width for segmented interpolation |

### Available Methods

#### `InterpolationMethod.LINEAR`

Use this when you want a true hard break between stations.

Behavior:

- piecewise linear
- no local smoothing
- good for intentionally sharp kinks

#### `InterpolationMethod.SEGMENTED`

This is the recommended mode for conventional aircraft wings.

Behavior:

- each spanwise segment stays straight away from the station transitions
- only a local window around each interior station is smoothed
- the user chooses whether that local blend is `C1` or `C2`
- the whole wing does not become a global spline

This is the key mode for constructible wings with root-kink-tip style layouts.

#### `InterpolationMethod.PYSPLINE`

Use this only when you explicitly want a globally curved spanwise evolution.

Behavior:

- global spline through all anchor stations
- useful for intentionally smooth and curved planforms
- not recommended for standard panel-like wings unless the geometry truly
  requires it

### Continuity Choice for `SEGMENTED`

- `ContinuityOrder.C1`
  - tangent-continuous local transition
- `ContinuityOrder.C2`
  - curvature-continuous local transition

### `blend_fraction`

`blend_fraction` defines the local transition half-width around each interior
station as a fraction of the neighboring span segment.

Practical interpretation:

- smaller values keep the straight segment dominant and make the transition more
  localized
- larger values spread the transition over a wider neighborhood

Typical starting values:

- `0.10` to `0.16` for tighter local transitions
- `0.16` to `0.22` for smoother and broader transitions

### Recommended Default for Aircraft Wings

```python
InterpolationSpec(
    method=InterpolationMethod.SEGMENTED,
    continuity=ContinuityOrder.C2,
    blend_fraction=0.18,
)
```

## 5. Optional Scalar Laws

`ScalarLawSpec` lets the user override specific spanwise quantities explicitly.

This is useful when:

- station values are not enough
- you want direct control of twist or thickness evolution
- you want a law with a different interpolation behavior than the main spine

Supported wing scalar laws are:

- `x`
- `y`
- `z`
- `chord`
- `roll_deg`
- `pitch_deg`
- `twist_deg`
- `thickness_scale`
- `te_thickness`

### Example: Twist Law

```python
ScalarLawSpec(
    name="twist_deg",
    anchors=(
        ScalarLawAnchor(0.0, 2.5),
        ScalarLawAnchor(0.38, 0.8),
        ScalarLawAnchor(1.0, -2.4),
    ),
    interpolation=InterpolationSpec(
        method=InterpolationMethod.SEGMENTED,
        continuity=ContinuityOrder.C2,
        blend_fraction=0.16,
    ),
)
```

## 6. Complete Wing Example

```python
from parametrization.aircraft import (
    ContinuityOrder,
    InterpolationMethod,
    InterpolationSpec,
    IntuitiveAirfoilProfileSpec,
    LiftingSurfaceBuildOptions,
    ProfileCatalog,
    ScalarLawAnchor,
    ScalarLawSpec,
    WingSpec,
    WingStationSpec,
    export_lifting_surface_iges,
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
            profile_id="kink_airfoil",
            leading_edge_radius=0.013,
            max_thickness=0.118,
            x_tmax=0.32,
            max_camber=0.016,
            x_cmax=0.41,
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

segmented_c2 = InterpolationSpec(
    method=InterpolationMethod.SEGMENTED,
    continuity=ContinuityOrder.C2,
    blend_fraction=0.16,
)

wing = WingSpec(
    wing_id="main_wing",
    semispan=14.5,
    stations=(
        WingStationSpec("root", 0.0, "root_airfoil", chord=6.2, twist_deg=2.5),
        WingStationSpec(
            "kink",
            0.38,
            "kink_airfoil",
            chord=4.0,
            twist_deg=0.8,
            sweep_le_deg=20.0,
            dihedral_deg=3.5,
        ),
        WingStationSpec(
            "tip",
            1.0,
            "tip_airfoil",
            chord=1.55,
            twist_deg=-2.4,
            sweep_le_deg=29.0,
            dihedral_deg=5.5,
        ),
    ),
    root_x=8.35,
    root_y=0.10,
    root_z=0.0,
    section_interpolation=segmented_c2,
    spine_interpolation=segmented_c2,
    scalar_laws=(
        ScalarLawSpec(
            name="twist_deg",
            anchors=(
                ScalarLawAnchor(0.0, 2.5),
                ScalarLawAnchor(0.38, 0.8),
                ScalarLawAnchor(1.0, -2.4),
            ),
            interpolation=segmented_c2,
        ),
    ),
)

result = export_lifting_surface_iges(
    wing.to_component_spec(),
    profiles,
    out_dir="parametrization/aircraft/example_outputs/my_wing",
    options=LiftingSurfaceBuildOptions(
        station_count=13,
        include_anchor_sections=True,
        airfoil_sample_count=301,
        blunt_te=True,
        tip_style="rounded",
    ),
)
```

## 7. Build and Export Options

The lifting-surface builder uses `LiftingSurfaceBuildOptions`.

### Main Parameters

| Parameter | Meaning |
| --- | --- |
| `station_count` | Number of generated spanwise stations |
| `station_etas` | Optional explicit station grid |
| `include_anchor_sections` | Force anchor stations into the build grid |
| `airfoil_sample_count` | Number of points used to sample each airfoil |
| `fit_n_ctl` | Optional pyGeo fitting control count |
| `k_span` | Optional spanwise spline order for pyGeo fitting |
| `tip_style` | `"rounded"` or `"pinched"` |
| `blunt_te` | Export blunt trailing edge |
| `rounded_te` | Use rounded trailing edge when blunt trailing edge is enabled |

### Notes

- if you use `SEGMENTED`, the builder automatically inserts extra support
  stations around interior wing breaks so the local `C1` or `C2` transition is
  actually represented in the loft
- you do not need to manually add those extra stations yourself

## 8. Typical Design Patterns

### Conventional Transport Wing

Use:

- 3 to 5 anchor stations
- `SEGMENTED`
- `C2`
- moderate `blend_fraction`

Good starting point:

```python
InterpolationSpec(
    method=InterpolationMethod.SEGMENTED,
    continuity=ContinuityOrder.C2,
    blend_fraction=0.16,
)
```

### Wing with a Hard Structural Kink

Use:

- `LINEAR`

This keeps the break explicit and does not smooth it.

### Intentionally Curved Wing

Use:

- `PYSPLINE`

Only use this if a global curved planform is actually desired.

## 9. Example Scripts

The current example scripts for the wing workflow are:

- `python parametrization/aircraft/examples/run_demo_wing_iges.py`
  - exports the demo wing to IGES
- `python parametrization/aircraft/examples/plot_demo_wing_views.py`
  - generates top, side, front, and 3D wing plots
- `python parametrization/aircraft/examples/plot_demo_wing_transition_comparison.py`
  - compares the same wing with `C1` and `C2` segmented transitions

Generated outputs are written to:

- `parametrization/aircraft/example_outputs/demo_main_wing/`

## 10. Current Scope and Limitations

Current strengths:

- generic root-kink-tip style wing definition
- iCST-driven airfoil profiles
- straight-by-segment wing layouts
- local `C1` or `C2` transitions between sections
- IGES export
- real-scale plotting

Current limitations:

- no dedicated higher-level planform API yet for area, aspect ratio, taper, or
  quarter-chord sweep targets
- no automatic structural reference surfaces yet
- no dedicated flap, slat, or control-surface parameterization yet

## 11. Implementation Files

If you want to inspect the implementation directly, the key files are:

- `parametrization/aircraft/profiles.py`
- `parametrization/aircraft/wing.py`
- `parametrization/aircraft/laws.py`
- `parametrization/aircraft/lifting_surface.py`
- `parametrization/aircraft/examples/common_demo.py`

## 12. Quick Summary

The current wing parameterization is based on:

- profile definitions
- anchor sections
- straight spanwise segments by default
- optional local `C1` or `C2` transitions at interior stations
- explicit user control over twist, thickness, chord, sweep, and dihedral

That makes the workflow suitable for conventional aircraft wings while still
leaving the door open for more curved geometries when needed.
