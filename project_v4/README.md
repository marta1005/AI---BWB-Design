# project_v4

Base nueva para geometrías BWB con esta parametrización de planta:

- `C1`: cuerda raíz
- `C2/C1`, `C3/C1`, `C4/C1`: cuerdas relativas en las 3 secciones exteriores
- `B1`, `B2`, `B3`: tramos spanwise normalizados con la semienvergadura
- `S1`, `S2`, `S3`: sweep del `LE` por tramo

La geometría mantiene continuidad `C2` en planta y usa:

- `CST` para los perfiles
- `pySpline` para las leyes spanwise, twist por sección e interpolación de secciones
- `pyGeo` para el loft y export `IGES`

Los perfiles usan ahora `CST` completos e independientes en `C1`, `C2`, `C3` y `C4`.
Ahora mismo cada sección usa:

- `6` coeficientes CST en extradós
- `6` coeficientes CST en intradós
- `tc_max` por sección
- `x_tmax` por sección
- `te_thickness` por sección

La clase CST está fijada por defecto a `N1=0.5`, `N2=1.0` para una familia tipo `NACA`.
El modo geométrico por defecto es ahora `cst_only`:

- la forma 2D se genera directamente desde los coeficientes `CST`
- `tc_max` y `x_tmax` se conservan como referencias/constraints geométricas posteriores
- `te_thickness` sigue disponible como parámetro geométrico separado del borde de salida

Sigue existiendo un modo legacy `enforce_targets` para reconstrucciones donde se quiera imponer
`tc_max` y `x_tmax` durante la generación del perfil.

La cabina no forma parte del modelo geométrico de `v4`.

## Módulos principales

- `topology.py`: semienvergadura y segmentos `B1/B2/B3`
- `specs.py`: specs de planta, familias de perfil, leyes spanwise, sampling y export
- `planform.py`: construcción `LE/TE` suave con continuidad `C1/C2`
- `sections.py`: construcción CST e interpolación spanwise
- `spanwise_laws.py`: leyes `thickness`, `twist`, `camber` y diedro
- `validation.py`: validación geométrica
- `volume.py`: evaluación de la constraint de volumen encerrado
- `exporters.py`: escritura de perfiles y exportación `pyGeo`
- `dependency_setup.py`: bootstrap y recompilación automática de `pyspline` si hace falta
- `builder.py`: ensamblado completo de la geometría
- `design_variables.py`: vector de variables para IA/optimización
- `design_space.py`: metadata, grupos, presets y muestreo local del espacio de diseño
- `gemseo_space.py`: adaptador de `project_v4` a `GEMSEO DesignSpace`

## Requisitos

Python recomendado:

- `Python 3.11.x`

Ficheros de requisitos:

- `project_v4/requirements-runtime.txt`: geometría, plotting y export base
- `project_v4/requirements-gemseo.txt`: capa GEMSEO
- `project_v4/requirements-full.txt`: stack completo recomendado

Guía resumida:

- `project_v4/docs/runtime_requirements.md`

## Estructura de ejemplos

Los scripts reales están ahora organizados así:

- `examples/reference/`: caso base, plots de referencia y comparación blunt-TE
- `examples/design_space/`: exploración de variantes, export de bounds e inspección GEMSEO
- `examples/legacy/`: nombres antiguos mantenidos como compatibilidad

Los scripts de primer nivel en `examples/` son wrappers para no romper tus comandos actuales.

## Ejemplo base

- `.venv/bin/python project_v4/examples/run_reference_example.py`
- `.venv/bin/python project_v4/examples/plot_reference_planform.py`

## Exploración del espacio de diseño

Hay una capa explícita de espacio de diseño en `design_space.py` para empezar a probar
variantes sin abrir todavía todo el vector geométrico. La idea es trabajar por grupos:

- `topology`: `B1/B2/B3`
- `planform`: posición raíz, cuerdas y sweeps
- `twist`: twist por sección
- `thickness_targets`: `tc_max`, `x_tmax`, `te_thickness`
- `camber_mode`: control adicional de camber
- `cst_c1` ... `cst_c4`: coeficientes CST completos por sección

Presets disponibles:

- `presentation_core`: planform + twist + thickness targets
- `planform_only`
- `aero_sections`
- `ai_geometry_core`: spans + chords + sweeps + nose blend + twists + full CST coefficients, while dihedral and TE thickness remain fixed by the reference design
- `gemseo_geometry_core`: spans + chords + sweeps + twists + full CST coefficients
- `section_shapes`
- `full_geometry`

Ejemplo para generar variantes válidas y compararlas:

- `.venv/bin/python project_v4/examples/explore_design_space.py`
- `.venv/bin/python project_v4/examples/explore_design_space.py --preset full_geometry --count 8 --variation-scale 0.15`
- `.venv/bin/python project_v4/examples/design_space/export_bounds_table.py`
- `.venv/bin/python project_v4/examples/design_space/inspect_gemseo_design_space.py --preset ai_geometry_core`
- `.venv/bin/python project_v4/examples/design_space/export_gemseo_bounds_table.py --preset ai_geometry_core`
- `.venv/bin/python project_v4/examples/design_space/sample_gemseo_lhs.py --preset ai_geometry_core --n-samples 32 --algo LHS --seed 7`

La salida se escribe en `project_v4/example_outputs/design_space_<preset>/` e incluye:

- overlay de plantas
- comparación de perfiles `C1` y `C4`
- CSV con las variables activas de cada variante
- JSON resumen

Para empezar a validar la parametrización, el preset recomendado es `presentation_core`,
porque mueve geometría global y leyes básicas sin abrir todavía todos los coeficientes CST.

La tabla completa de bounds se genera en:

- `project_v4/docs/design_variable_bounds.csv`
- `project_v4/docs/design_variable_bounds.md`
- `project_v4/docs/geometric_parametrization_tables.md`
- `project_v4/docs/profile_relations.md`

## GEMSEO

El adaptador GEMSEO está en `gemseo_space.py` y construye un `DesignSpace` agrupado por bloques
geométricos, no como 84 escalares sueltos. La topología `B1/B2/B3` se pasa como pesos positivos
que el adaptador re-proyecta a ratios que suman `1.0`.

Instalación recomendada:

```bash
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r project_v4/requirements-gemseo.txt
```

Importante:

- `gemseo 6.x` requiere Python `>=3.10,<3.14`
- si `pip` muestra un conflicto del tipo `(constraint) numpy==1.21.6`, ese pin viene
  de un constraint externo del entorno, no de `project_v4`
- antes de instalar GEMSEO, comprueba `python --version`, `echo $PIP_CONSTRAINT`
  y `python -m pip config list -v`

Inspección rápida:

- `.venv/bin/python project_v4/examples/design_space/inspect_gemseo_design_space.py`
- `.venv/bin/python project_v4/examples/design_space/inspect_gemseo_design_space.py --preset full_geometry`
- `.venv/bin/python project_v4/examples/design_space/sample_gemseo_lhs.py --preset ai_geometry_core --n-samples 32`

Guía corta de instalación:

- `project_v4/docs/gemseo_install.md`

## Compatibilidad

Se mantienen wrappers en `examples/` con los nombres antiguos para no romper flujos previos,
pero la referencia limpia de `v4` es `run_reference_example.py`.
