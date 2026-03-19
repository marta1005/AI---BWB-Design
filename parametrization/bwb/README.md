# Sectioned BWB Package

`parametrization/bwb` contains the generic sectioned blended-wing-body
parametrization logic.

This is the package used by CTA and by any future BWB-family geometry.

Main modules:

- `topology.py`: semispan segmentation and anchor placement
- `specs.py`: configuration dataclasses
- `design_variables.py`: editable BWB parameter vector
- `design_space.py`: bounds, flattening, metadata, and sampling helpers
- `planform.py`: smooth sectioned planform reconstruction
- `spanwise_laws.py`: twist/camber/offset interpolation
- `sections.py`: spanwise section profile reconstruction
- `validation.py`: section and loft checks
- `volume.py`: volume constraint evaluation
- `exporters.py`: airfoil and IGES export helpers
- `builder.py`: end-to-end geometry preparation/export
