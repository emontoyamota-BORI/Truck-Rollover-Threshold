# Model C: multibody simulation in PyChrono

Multibody model (Chrono::Vehicle, Python bindings) of the 4x2 tractor + two-axle semitrailer of Yoo et al. (2025)
on the accident ramp of Zhang et al. (2024), with two suspensions: **multi-leaf steel springs** and
**trailing-arm air suspension**. It is used to verify the static roll-plane model (Model B) of the article.

Code author: Luis Montoya (UANL-FIME). License: MIT (see the repository root).
The file and variable names and the inline comments of this folder are in Spanish; this README documents the
workflow in English. The original Spanish notes are kept in `LEEME.md`.

## Requirements

- A conda environment with PyChrono (version 10 was used) and its Irrlicht module for the 3D view, e.g.
  `conda create -n chrono -c projectchrono -c conda-forge pychrono` (see the PyChrono installation guide).
  See `environment.yml`.
- numpy and matplotlib; optional for MP4 videos and Excel tables: `pip install imageio imageio-ffmpeg openpyxl`.
- Chrono cannot open paths that contain accented characters: the scripts copy `datos/` to an ASCII folder in the
  user's home directory and copy the results back to `resultados/`.

## How to run

1. Activate the environment (`conda activate chrono`) and run `paso0_verificar.py`: it checks that all modules
   are available.
2. Run `correr_todo.py`. It executes the whole workflow: input data → tilt table → ramp runs → comparison
   videos → plots and tables. When the 3D view opens, do not close it; it closes automatically at the end of
   each run.

Each step can also be run on its own:

| File | Needs PyChrono | Purpose |
|---|---|---|
| `config.py` | No | All simulation parameters: masses, geometry, suspensions (K_phi per axle, bounce frequency, spring spacing), tires, run cases (`CASOS_RAMPA`) |
| `paso0_verificar.py` | No | Environment check |
| `paso1_preparar.py` | No | Ramp mesh and vehicle JSON files for both suspensions; parameter tables |
| `geometria_rampa.py` | No | Ramp geometry (Zhang et al., 2024, Table 3) and terrain mesh |
| `vehiculo_json.py` | No | Chrono::Vehicle JSON files of the tractor and semitrailer for one suspension: king-pin load for the axle-load split of Yoo et al. (2025), CG of the sprung masses, spring stiffness from the bounce frequency, damping, TMeasy tires |
| `paso2_mesa_basculante.py [susp] [--ver]` | Yes | Static loads and virtual tilt table: first lift-off, lift-off of each axle, roll gradient |
| `paso3_rampa.py susp mode v_kmh mu video` | Yes | One run on the ramp, e.g. `python paso3_rampa.py aire constante 40 0.8 1` |
| `chrono_modelo.py` | Yes | Builds the model for one suspension: tractor (`WheeledVehicle`) + semitrailer (`WheeledTrailer`) from the JSON files, TMeasy tires, powertrain, torsional springs with the auxiliary roll stiffness, terrain (ramp mesh or flat ground) and telemetry functions |
| `lanzar.py` | No | Runs each simulation in a separate Python process, so that Chrono and Irrlicht start clean in every run |
| `video.py [run]` | No | Side-by-side GIF/MP4 (multi-leaf / air) |
| `paso4_resultados.py` | No | All comparison plots and tables |
| `tablas.py` | No | Helpers to read and write CSV (and Excel) tables |

Suspension keys: `multihoja` = multi-leaf steel springs, `aire` = trailing-arm air suspension.

## Outputs (`resultados/`)

```
resultados/
  comparacion/tablas/   comparison tables used by the article and by tractor_semitrailer_srt.py
  multihoja/ and aire/
    mesa/               tilt table: mesa_basculante.csv, resumen_mesa.json
    telemetria/         <run>.csv (every 0.05 s) and <run>_resumen.json
    tablas/             <run>_cada5m.csv (telemetry every 5 m of station)
```

Videos, frames and plots (about 150 MB) are not included in the repository; `correr_todo.py` regenerates them.

Main tables in `resultados/comparacion/tablas/`:

| Table | Content | Used in |
|---|---|---|
| `tabla_mesa_curvas.csv` | LTR of each axle, whole-vehicle LTR, tractor and semitrailer roll, Fz per side, every 0.01 g | Part B.3 (`CHRONO_TILT`), article Sec. 3.2 |
| `tabla_mesa_basculante.csv` | First lift-off, lift-off acceleration of each axle, lift-off order, roll gradient [deg/g] | Part B.3, article Sec. 3.2 |
| `tabla_rampa_estacionario.csv` | Mean values in the constant-radius curves (R = 125 m and 60 m) of each run: speed, a_y, LTR, roll, steering, articulation, slip angle, Fz and Fy per axle | Part B.3 (`CHRONO_R60`), Part D.4, article Sec. 3.4 |
| `tabla_eventos.csv` | Lift-off and rollover events of each run | article Sec. 3.4 |
| `tabla_validacion_ref11.csv` | LTR in the 60 m curve at 30/40/50 km/h versus TruckSim (Zhang et al., 2024) | article Sec. 3.5 |
| `tabla_parametros_suspension.csv`, `tabla_masas_y_cg.csv` | Model parameters | article Secs. 2.3–2.4 |

## Sign convention (the ramp turns to the right)

Side L = outer side of the curve, side R = inner side. a_y, F_y and steering are positive towards the centre of
the curve; roll is positive towards the outside; LTR = (outer − inner)/(outer + inner). On the ramp, the tractor
and semitrailer roll angles are relative to the road surface (superelevation removed), and the suspension roll
of each axle is the chassis rotation relative to that axle. Fz is the load normal to the road and Fy the lateral
force in the road plane.

## Known limitation

The roll-centre height is the same for both suspensions, because it is fixed by the Chrono axle template
(`SolidThreeLinkAxle`). This is stated in the article.
