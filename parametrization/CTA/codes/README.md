# CTA Codes

Helper scripts are grouped by topic:

- [airfoils](/Users/martaarnabatmartin/Desktop/BWB/parametrization/CTA/codes/airfoils)
  Raw glider inspection and CST fitting.

- [plotting](/Users/martaarnabatmartin/Desktop/BWB/parametrization/CTA/codes/plotting)
  Shared wing plotting helpers, vertical comparison plots, and animations.

- [exports](/Users/martaarnabatmartin/Desktop/BWB/parametrization/CTA/codes/exports)
  Bounds tables and IGES export.

Removed during cleanup:

- the old interactive 3D viewer script, because the main plot generator
  already produces the saved 3D figure
- the exploratory profile lab script, because it was not part of the active CTA
  workflow
- the old validation package name, because the CTA package itself is already a
  concrete parametrization case
- the explicit reproduction/consistency script layer, because the CTA case is
  already defined directly in `case.py`
- the GEMSEO-specific layer, because this CTA package is no longer using that
  workflow
