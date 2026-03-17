# Section Profile Relations

`project_v4` now separates two use cases:

1. `AI design space`: all control sections are treated as independent unless the user deliberately ties them.
2. `Manual parametrization`: sections can share part or all of their profile definition.

The relation layer is applied on top of the base CST definitions of `C1`, `C2`, `C3` and `C4`.

## Supported relations

Each section can define any of these links to an earlier section:

- `shape_source_index`: copy the full profile definition from another section
  - upper CST coefficients
  - lower CST coefficients
  - `tc_max`
  - `x_tmax`
  - `te_thickness`
- `tie_cst(...)`: copy only upper and lower CST coefficients from another section
- `upper_source_index`: copy only the upper-surface CST coefficients
- `lower_source_index`: copy only the lower-surface CST coefficients
- `te_source_index`: copy only the trailing-edge thickness

These links are resolved from inboard to outboard. A section can only reference an earlier section.

## What makes geometric sense

The following combinations are supported and coherent:

- fully independent profiles
- same full shape as another section, but with a different local chord
- same upper surface as another section
- same lower surface as another section
- same CST coefficients, with a different local chord

The following idea is deliberately **not** treated as a native relation:

- "same CST coefficients but different thickness"

That is not a pure CST relation. If thickness must change independently, the section should either:

- use its own CST coefficients, or
- share the CST and change only the physical scale through the local chord

## API

The main dataclass is:

- `SectionProfileRelationSpec`

Utility helpers are provided in:

- `independent_profile_relations()`
- `tie_cst(...)`
- `tie_shape(...)`
- `tie_upper(...)`
- `tie_lower(...)`
- `tie_te_thickness(...)`

## Example

```python
from project_v4 import (
    SectionedBWBDesignVariables,
    independent_profile_relations,
    tie_cst,
    tie_upper,
)

design = SectionedBWBDesignVariables.reference_seed()

relations = independent_profile_relations()
relations = tie_cst(relations, target_index=1, source_index=0)     # C2 = same CST coeffs as C1
relations = tie_upper(relations, target_index=3, source_index=2)   # C4 upper = C3 upper

config = design.to_model_config(
    profile_generation_mode="cst_only",
    profile_relations=relations,
)
```

In that example:

- `C2` uses the same upper/lower CST coefficients as `C1`
- `C2` can still have a different chord because chord belongs to the planform, not to the profile relation
- `C4` only shares the upper surface with `C3`
