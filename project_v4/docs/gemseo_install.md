# GEMSEO installation for project_v4

## Recommended environment

Use the existing project virtual environment based on Python 3.11:

```bash
source .venv/bin/activate
python --version
```

Expected result:

```text
Python 3.11.x
```

GEMSEO 6 requires Python `>=3.10,<3.14`. Do not try to install it in an older
Python 3.9 environment.

## Install GEMSEO

```bash
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r project_v4/requirements-gemseo.txt
```

## Common installation failure

If `pip` reports a conflict such as:

```text
The user requested (constraint) numpy==1.21.6
```

then the problem does not come from `project_v4`. It means the current shell,
CI job or `pip` configuration is forcing an external NumPy constraint.

Check it with:

```bash
python --version
echo "$PIP_CONSTRAINT"
python -m pip config list -v
```

If you see Python `3.9` or any old NumPy constraint, fix that first:

```bash
unset PIP_CONSTRAINT
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r project_v4/requirements-gemseo.txt
```

If your current virtual environment is still Python `3.9`, create a Python `3.11`
virtual environment and reinstall the project stack there.

## Quick validation

```bash
python - <<'PY'
import gemseo
print(gemseo.__version__)
PY
```

## Notes for this repository

- `project_v4` already supports GEMSEO through `project_v4/gemseo_space.py`.
- If the Python interpreter changes, rebuild or revalidate the local `pyspline/pyGeo` stack before exporting IGES.
- The initial recommended GEMSEO preset is `presentation_core`.
