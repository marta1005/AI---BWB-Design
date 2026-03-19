# Parametrization Layout

This repository now separates reusable geometry utilities from architecture-specific
parametrizations.

## Packages

- `parametrization/shared`
  Shared utilities reused by multiple parametrization families.
  This package contains dependency loading for `pyspline`/`pyGeo`, CST helpers,
  and airfoil I/O helpers.

- `parametrization/bwb`
  Generic sectioned blended-wing-body tooling.
  This is the home of the planform, design variables, topology, validation,
  exporters, and geometry builder used by CTA and any future BWB variants.

- `parametrization/CTA`
  CTA-specific presets, examples, docs, and reference outputs.
  CTA depends on `parametrization.bwb`, not on `aircraft`.

- `parametrization/aircraft`
  Generic whole-aircraft building blocks.
  This package uses only `parametrization.shared` helpers plus its own local
  wing/fuselage/component logic.

- `parametrization/core`
  Backward-compatible import shim.
  Existing scripts can still import `parametrization.core.*`, but new code should
  use `parametrization.shared.*` or `parametrization.bwb.*` directly.

## Rule Of Thumb

- If it is reusable across aircraft families, put it in `shared`.
- If it is specific to sectioned BWB logic, put it in `bwb`.
- If it is CTA-specific, keep it in `CTA`.
- If it is generic aircraft assembly logic, keep it in `aircraft`.
