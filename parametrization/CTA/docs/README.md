# CTA Parametrization (`parametrization/CTA`)

This package isolates the CTA-focused BWB parametrization on top of
`parametrization/bwb`, with reusable low-level helpers living in
`parametrization/shared`, so the next step (global aircraft parametrization with
fuselage, VTP, additional sections, etc.) can be developed in parallel without
mixing concerns.

## Sweep Naming

The requested naming remap is implemented here:

- `S`  -> legacy `s1_deg` (fixed)
- `S1` -> legacy `s2_deg` (variable)
- `S2` -> legacy `s3_deg` (variable)

## Fixed Parameters in CTA AI Design Space

- `S` (`s1_deg`)
- `C1` (`c1_root_chord`)
- `C2` (`c2_c1_ratio`)
- straight TE segment `C0->C1`
- straight TE segment `C3->C4`

## Variable Parameters in CTA AI Design Space

- semi-span and B ratios:
  - `span`
  - `b1_span_ratio`
  - `b2_span_ratio`
  - `b3_span_ratio`
- sweeps:
  - `S1` (`s2_deg`)
  - `S2` (`s3_deg`)
- chord family:
  - `c3_c1_ratio`
  - `c4_c1_ratio` (used for outboard/tip family, i.e., C4/C5 in CTA figure notation)
- twists:
  - `twist_c1_deg`
  - `twist_c2_deg`
  - `twist_c3_deg`
  - `twist_c4_deg`
- CST:
  - full upper/lower 6+6 coefficients for C1, C2, C3 and C4

## Usage

- Export bounds table and metadata:
  - `python parametrization/CTA/examples/export_cta_bounds_table.py`
- Export CTA reference IGES from this package:
  - `python parametrization/CTA/examples/run_reference_iges.py`

Generated tables in `docs/` and generated plots or IGES files in `example_outputs/`
are intentionally not treated as source files; regenerate them from the scripts above.
