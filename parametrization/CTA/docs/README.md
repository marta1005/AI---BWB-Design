# CTA Parametrization (`parametrization/CTA`)

This package isolates the CTA-specific design-space logic on top of the generic
`parametrization/bwb` geometry builder.

The public CTA layer uses the legacy naming agreed for the project:

- `S`  -> internal `s1_deg` and fixed
- `S1` -> internal `s2_deg` and active
- `S2` -> internal `s3_deg` and active
- `C0` -> active body chord
- `C3` -> active transition-wing chord
- `C4` -> active outer-wing chord
- `C5` -> active wing-tip chord
- `B1` -> fixed in meters
- `B2` -> active as a fraction of wing span, i.e. `B2/(B2+B3)`
- `wing span` -> active and equal to `B2 + B3`

## CTA Public Design Space

Fixed or constrained:

- `S`
- `B1`
- `C1`, derived from `C0` with fixed `S` and straight `TE(C0->C1)`
- `C2`, kept only as an internal spline helper on the inboard trailing edge
- `C3` spanwise location fixed at `y = 8.041 m`
- twist constant from root to `C3`

Active public variables:

- `wing span = B2 + B3`
- `B2/(B2+B3)`
- `C0`
- `C3`
- `C4`
- `C5`
- `S1`
- `S2`
- twist at `C0/C3`, `C4` and `C5`
- full CST coefficients for `C0`, `C3`, `C4` and `C5`

Geometry note:

- `TE(C3->C4)` is not constrained to stay straight anymore because `C3` and `C4`
  are both active design variables.

## Useful Commands

Public CTA bounds:

```bash
python parametrization/CTA/codes/exports/export_cta_bounds_table.py
```

Geometry export:

```bash
python parametrization/CTA/codes/exports/run_cta_iges.py
python parametrization/CTA/codes/plotting/plot_cta_views.py
python parametrization/CTA/codes/plotting/plot_cta_vs_glider_vertical_overlay.py
```

## How CTA Parameters Reach pyGeo

The CTA layer does not pass semantic parameters such as sweep names, helper
locations, or `te_inboard_radius_factor` directly into `pyGeo`.

Instead, the pipeline is:

1. `parametrization/CTA/case.py`
   Builds the CTA-specific `SectionedBWBModelConfig`.
2. `parametrization/bwb/specs.py`
   Converts CTA public sections and helper points into explicit LE/TE control points.
3. `parametrization/bwb/planform.py`
   Resolves those control points into continuous spanwise functions `le_x(y)` and `te_x(y)`.
4. `parametrization/bwb/builder.py`
   Samples the geometry along the span to compute:
   - `leading_edge_x`
   - `trailing_edge_x`
   - `chord`
   - `twist_deg`
   - `vertical offsets`
5. `parametrization/bwb/exporters.py`
   Sends those sampled arrays to `pyGeo("liftingSurface", ...)`.

So a CTA planform parameter affects `pyGeo` only after it has been converted
into explicit sampled geometry.
