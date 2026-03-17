# project_v4 Geometric Parametrization Tables

This document separates the current `project_v4` geometry controls into:

1. **Active design-space parameters** to be used for DOE, surrogate modeling, generative AI, and optimization.
2. **Fixed internal modeling parameters** that are not part of the current design vector, but can still be retuned if the modeling strategy changes.

The intent of the current design-space definition is to keep the exploration focused on:

- spanwise topology
- chord law
- leading-edge sweep law
- sectional twist
- full CST coefficients for each control section

The current active design-space definition keeps `dihedral = 0` fixed and leaves thickness targets, CST class exponents, root-blend settings, and other modeling controls constant.

## 1. Active Design-Space Parameters

| Group | Parameter(s) | Symbol | Units | Role in Geometry |
|---|---|---|---|---|
| Semi-span topology | `span` | `b/2` | m | Half-span of the aircraft. Sets the global geometric scale in the spanwise direction. |
| Semi-span topology | `b1_span_ratio` | `B1 / (b/2)` | `-` | Inboard span partition, normalized by semi-span. |
| Semi-span topology | `b2_span_ratio` | `B2 / (b/2)` | `-` | Transition-wing span partition, normalized by semi-span. |
| Semi-span topology | `b3_span_ratio` | `B3 / (b/2)` | `-` | Outboard span partition, normalized by semi-span. |
| Chord law | `c1_root_chord` | `C1` | m | Root chord. Main absolute chord reference. |
| Chord law | `c2_c1_ratio` | `C2 / C1` | `-` | Chord ratio at control section `C2`. |
| Chord law | `c3_c1_ratio` | `C3 / C1` | `-` | Chord ratio at control section `C3`. |
| Chord law | `c4_c1_ratio` | `C4 / C1` | `-` | Chord ratio at control section `C4`. |
| LE sweep law | `s1_deg` | `S1` | deg | Leading-edge sweep angle of the first segment (`C0 -> C1`). |
| LE sweep law | `s2_deg` | `S2` | deg | Leading-edge sweep angle of the second segment (`C1 -> C2`). |
| LE sweep law | `s3_deg` | `S3` | deg | Leading-edge sweep angle of the third segment (`C2 -> C3`). |
| Twist law | `twist_c1_deg` | `twist_C1` | deg | Twist angle at control section `C1`. |
| Twist law | `twist_c2_deg` | `twist_C2` | deg | Twist angle at control section `C2`. |
| Twist law | `twist_c3_deg` | `twist_C3` | deg | Twist angle at control section `C3`. |
| Twist law | `twist_c4_deg` | `twist_C4` | deg | Twist angle at control section `C4`. |
| CST shape at section C1 | `c1_upper_cst_0..5` | `A^u_{C1,0..5}` | `-` | Upper-surface Kulfan/CST coefficients at `C1`. |
| CST shape at section C1 | `c1_lower_cst_0..5` | `A^l_{C1,0..5}` | `-` | Lower-surface Kulfan/CST coefficients at `C1`. |
| CST shape at section C2 | `c2_upper_cst_0..5` | `A^u_{C2,0..5}` | `-` | Upper-surface Kulfan/CST coefficients at `C2`. |
| CST shape at section C2 | `c2_lower_cst_0..5` | `A^l_{C2,0..5}` | `-` | Lower-surface Kulfan/CST coefficients at `C2`. |
| CST shape at section C3 | `c3_upper_cst_0..5` | `A^u_{C3,0..5}` | `-` | Upper-surface Kulfan/CST coefficients at `C3`. |
| CST shape at section C3 | `c3_lower_cst_0..5` | `A^l_{C3,0..5}` | `-` | Lower-surface Kulfan/CST coefficients at `C3`. |
| CST shape at section C4 | `c4_upper_cst_0..5` | `A^u_{C4,0..5}` | `-` | Upper-surface Kulfan/CST coefficients at `C4`. |
| CST shape at section C4 | `c4_lower_cst_0..5` | `A^l_{C4,0..5}` | `-` | Lower-surface Kulfan/CST coefficients at `C4`. |

### Current design-space scope

The current recommended **primary design space** is:

- `span`
- `B1, B2, B3` ratios
- `C1, C2/C1, C3/C1, C4/C1`
- `S1, S2, S3`
- `twist_C1, twist_C2, twist_C3, twist_C4`
- full `6 + 6` CST coefficients at each control section

This gives direct control over:

- global wing size
- spanwise partitioning
- planform chord distribution
- segment-by-segment sweep
- spanwise twist law
- local airfoil shape at each control section

## 2. Fixed Internal Modeling Parameters That Can Be Retuned

These parameters are currently **fixed in the implementation**. They are not part of the main design vector, but they still affect the model behavior and can be changed if the geometry strategy needs to be retuned.

| Group | Parameter | Current Value | Role in Geometry / Modeling |
|---|---|---:|---|
| Planform continuity | `continuity_order` | `2` | Uses `C2` continuity at planform junctions. If changed to `1`, the model becomes `C1` continuous. |
| Planform junction smoothing | `blend_fraction` | `0.10` effective | Sets the spanwise width of the spline transition zone around each planform junction. |
| Planform constructibility | `min_linear_core_fraction` | `0.75` | Guarantees that each planform segment keeps a linear core. Curvature is concentrated near the junctions only. |
| Trailing-edge exact segment rule | `te_exact_segments` | `(0,)` | Forces the first TE segment (`C0 -> C1`) to remain exactly straight. |
| Root nose blending | `symmetry_blend_y` | driven by `nose_blend_y` | Controls how the root region is rounded close to the symmetry plane. |
| CST class function | `cst_n1` / `cst_n2` | `0.50 / 1.00` | Kulfan class exponents. The current default is NACA-like behavior. |
| Thickness evaluation window | `x_tc_window` | `(0.15, 0.65)` | Chordwise window used to evaluate and enforce maximum thickness. |
| Geometry validation window | `x_valid_window` | `(0.02, 0.98)` | Chordwise window used for robust profile validity checks away from LE/TE numerical singularities. |
| Camber shaping mode | `camber_mode_center` | `3.5` | Center of the modal CST-based camber perturbation. |
| Camber shaping mode | `camber_mode_width` | `2.0` | Width of the modal CST-based camber perturbation. |
| Twist interpolation | `twist_deg.interpolation` | `pyspline` | Section twist is interpolated spanwise with `pySpline`. |
| Camber interpolation | `camber_delta.interpolation` | `pyspline` | Camber deltas are interpolated spanwise with `pySpline`. |
| Section interpolation | `section_interpolation` | `pyspline` | CST coefficients, `tc_max`, `x_tmax`, and `te_thickness` are interpolated spanwise with `pySpline`. |
| Airfoil sampling | `num_airfoil_points` | `241` | Number of chordwise points used to discretize each 2D profile. |
| Loft sampling | `num_base_stations` | `41` | Number of spanwise base stations used before lofting/export. |
| Section curve resolution | `section_curve_n_ctl` | `41` | Resolution used for section curves in the loft preparation. |
| Span spline order | `k_span` | `4` | Order of the spanwise `pySpline` interpolation. |
| Airfoil distribution mode | `airfoil_distribution_mode` | `all` | All sampled airfoil points are used in the loft process, not only sparse anchors. |
| pyGeo tip treatment | `tip_style` | `rounded` | How the tip is closed during `pyGeo` surface generation. |
| pyGeo TE treatment | `blunt_te` | `True` | Preserves a finite trailing-edge thickness in the exported loft instead of collapsing it to a sharp edge. |
| Topology identifier | `topology_name` | `cbs_bwb_v4` | Internal label of the current topology convention. |

## Recommended message for presentations

> The current `project_v4` design space is intentionally restricted to span, span partition, chord law, leading-edge sweep by segment, twist at the control sections, and full CST coefficients at each control section, while dihedral is fixed to zero. Continuity, spline blending, trailing-edge exactness, sampling, and export behavior are controlled by fixed internal modeling parameters that can still be retuned if needed.
