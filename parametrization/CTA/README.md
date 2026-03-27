# CTA Parametrization Package

This folder contains the CTA-only parametrization and design-space setup on top of
`parametrization/bwb`.

- Shared utilities: `parametrization/shared/`
- Generic BWB layer used by CTA: `parametrization/bwb/`
- Core reference adapters: `parametrization/CTA/reference.py`
- CTA AI design-space (fixed + variable split): `parametrization/CTA/design_space.py`
- CTA GEMSEO adapter: `parametrization/CTA/gemseo_space.py`
- CTA GEMSEO requirements: `parametrization/CTA/requirements-gemseo.txt`
- Docs: `parametrization/CTA/docs/README.md`
- Generated bounds tables and plots are rebuilt on demand from the helper scripts in `parametrization/CTA/codes/`.
- Generated files live under `parametrization/CTA/outputs/`.
- Codes:
  - `parametrization/CTA/codes/export_cta_bounds_table.py`
  - `parametrization/CTA/codes/export_cta_gemseo_bounds_table.py`
  - `parametrization/CTA/codes/inspect_cta_gemseo_design_space.py`
  - `parametrization/CTA/codes/sample_cta_gemseo_doe.py`
  - `parametrization/CTA/codes/profile_parametrization_lab.py`
  - `parametrization/CTA/codes/plot_cta_reference_views.py`
  - `parametrization/CTA/codes/show_cta_reference_3d.py`
  - `parametrization/CTA/codes/run_reference_iges.py`

## How CTA Geometry Reaches pyGeo

CTA-specific parameters are not passed directly to `pyGeo` as semantic knobs.
They are first converted into explicit geometry:

1. `parametrization/CTA/reference.py`
   Defines the public CTA reference choices and maps them into a
   `SectionedBWBModelConfig`.
2. `parametrization/bwb/specs.py`
   Builds the public control points for the leading and trailing edges,
   including hidden helper points such as the inboard `TE(C1->C3)` spline helper.
3. `parametrization/bwb/planform.py`
   Turns those points into continuous spanwise axes (`le_x(y)` and `te_x(y)`),
   using segmented C1/C2 transitions or the dedicated spline bridge.
4. `parametrization/bwb/builder.py`
   Samples the resolved geometry on span stations and computes:
   - `leading_edge_x`
   - `trailing_edge_x`
   - `chord = TE - LE`
   - `twist_deg`
   - `vertical offsets`
5. `parametrization/bwb/exporters.py`
   Passes those sampled arrays to `pyGeo("liftingSurface", ...)` as:
   - `x = leading_edge_x`
   - `scale = chord`
   - `z = span`
   - `y = vertical offsets`
   - `rotZ = twist`

So, for example, the intensity of the inboard trailing-edge curve leaving `C1`
changes the final `te_x(y)` distribution first, and only then reaches `pyGeo`
through the sampled `x` and `scale` arrays.
