# CTA Parametrization Package

This folder contains the CTA-only parametrization and design-space setup on top of `parametrization/core`.

- Shared geometry core: `parametrization/core/`
- Core reference adapters: `parametrization/CTA/reference.py`
- CTA AI design-space (fixed + variable split): `parametrization/CTA/design_space.py`
- Docs: `parametrization/CTA/docs/README.md`
- Generated bounds tables and plots are rebuilt on demand from the example scripts.
- Examples:
  - `parametrization/CTA/examples/export_cta_bounds_table.py`
  - `parametrization/CTA/examples/profile_parametrization_lab.py`
  - `parametrization/CTA/examples/plot_cta_reference_views.py`
  - `parametrization/CTA/examples/show_cta_reference_3d.py`
  - `parametrization/CTA/examples/run_reference_iges.py`
