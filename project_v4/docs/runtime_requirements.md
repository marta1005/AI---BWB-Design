# project_v4 runtime requirements

## Recommended Python

- Recommended: `Python 3.11.x`
- Acceptable for GEMSEO 6.x: `Python >=3.10,<3.14`
- Not recommended here: `Python 3.9`

Reason:

- `project_v4` itself can be kept mostly Python 3.9-friendly, but the intended stack
  now includes `GEMSEO 6.x`, which requires Python `>=3.10`.
- The local `pySpline + pyGeo + GEMSEO` stack was designed around a Python 3.11 virtual environment.

## Python packages

### Minimum runtime used by project_v4

These are the direct Python packages needed to run the geometry, plotting and IGES export workflow:

```text
numpy==2.4.2
scipy==1.17.1
matplotlib==3.10.8
mpi4py==4.1.1
mdolab-baseclasses==1.8.4
packaging==26.0
```

File:

- `project_v4/requirements-runtime.txt`

### GEMSEO extension

For design-space work with GEMSEO:

```text
gemseo>=6,<7
```

Files:

- `project_v4/requirements-gemseo.txt`
- `project_v4/requirements-full.txt`

## Local MDOLab geometry stack

`project_v4` does not depend on `pySpline` and `pyGeo` from PyPI. It uses the local copies in:

- `project/pyspline`
- `project/pygeo`

`pySpline` is compiled locally and `pyGeo` is then imported from the local source tree.

Relevant package metadata in this repository:

- `project/pyspline/setup.py` requires:
  - `numpy>=1.21`
  - `scipy>=1.7`
- `project/pygeo/setup.py` requires:
  - `numpy>=1.21`
  - `pyspline>=1.1`
  - `scipy>=1.7`
  - `mpi4py>=3.1.5`
  - `mdolab-baseclasses`
  - `packaging`

## System requirements

To compile and use the local `pySpline` library, you also need:

- `make`
- `gfortran`
- a C compiler (`gcc` or `cc`)
- `python3-config` matching the active interpreter
- `f2py` from the active NumPy installation
- an MPI runtime compatible with `mpi4py`

## Recommended installation sequence

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r project_v4/requirements-full.txt
```

Then build the local spline library:

```bash
cd project/pyspline
cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk
make clean
make -B module \
  PYTHON=$(which python) \
  PYTHON-CONFIG=$(which python3-config) \
  F2PY=$(which f2py)
```

## Validation

After installation:

```bash
python - <<'PY'
from project_v4.dependency_setup import load_pyspline_curve, load_pygeo_class
print(load_pyspline_curve(rebuild_if_needed=False))
print(load_pygeo_class(rebuild_if_needed=False))
PY
```

and:

```bash
python project_v4/examples/run_reference_example.py
```

