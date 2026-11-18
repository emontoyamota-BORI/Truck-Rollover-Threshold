# Physics-based static rollover threshold of a tractor–semitrailer: multi-leaf versus air suspension

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24329653.svg)](https://doi.org/10.5281/zenodo.24329653)

Companion code and data of the article

> L. Montoya, J. Ibarra, R. A. Ramirez-Mendoza, *A Physics-Based Static Rollover Threshold with Measured
> Suspension Roll Stiffness for Safe Curve-Speed Estimation of a Tractor–Semitrailer: Multi-Leaf versus Air
> Suspension* (manuscript submitted for publication).

The repository contains everything needed to reproduce the numerical results, tables and Figures 3, 4, 7, 8
and 9 of the article.

## Models

| Model | Description | Where |
|---|---|---|
| **A** | Reproduction of the quasi-static safe-speed model of Yoo, Kim & Lee (2025): imposed SRT, fixed roll-stiffness split per axle. | `tractor_semitrailer_srt.py`, Part A |
| **B** | Proposed model: axle-by-axle static roll-plane model. The SRT is computed from the measured roll stiffness of each suspension, the tire vertical stiffness and the lateral CG shift. | `tractor_semitrailer_srt.py`, Part B |
| **C** | Multibody simulation in PyChrono (Chrono::Vehicle), used to verify Model B: virtual tilt table and constant-speed runs on the accident ramp. | `multibody_pychrono/` |

## Repository structure

```
tractor_semitrailer_srt.py   Models A and B, comparison with Model C, sensitivity studies (single file)
make_article_figures.py      Regenerates Figures 3, 4, 7, 8 and 9 of the article into figures/
figures/                     Figures produced by make_article_figures.py
requirements.txt             Python dependencies of the two scripts above
multibody_pychrono/          Model C: PyChrono scripts, input data and output tables (see its README)
CITATION.cff                 Citation metadata
LICENSE                      MIT license
```

## Quick start

**Google Colab.** Upload `tractor_semitrailer_srt.py` to a new notebook (or paste its content into one cell) and
run it. All tables are printed and all diagnostic plots are shown in the notebook.

**Local Python (≥ 3.9).**

```bash
pip install -r requirements.txt
python tractor_semitrailer_srt.py      # prints all tables (about 10 s)
python make_article_figures.py         # writes the article figures to figures/
```

Set `SHOW_MODEL_A_FIGURES = True` at the top of the main script to also plot the reproduction of Figures 4–8 of
Yoo et al. (2025).

## Correspondence between the code and the article

| Code section | Content | Article |
|---|---|---|
| A.1–A.2 | Vehicle, friction and ramp parameters | Secs. 2.2–2.3 |
| A.3–A.6 | Model A: loads, LTR, SRT calibration, safe speed, forward–backward speed profile | Sec. 2.5 |
| B.1–B.2 | Tire and suspension parameters of Model B | Secs. 2.3–2.4 |
| B.3 | Results of Model C (hard-coded from the PyChrono output tables) | Sec. 2.8 |
| B.4 | Model B: roll-plane equilibrium, axle lift-off, SRT as the peak of the stability curve | Sec. 2.6 |
| C.1–C.2 | Parameter tables and consistency checks | Secs. 2.3–2.4 |
| C.3–C.4 | Tilt table: axle LTR, first lift-off, SRT, roll angle and CG shift | Sec. 3.2, Figs. 3–4 |
| C.5 | LTR and rollover speed in the R = 60 m curve | Sec. 3.5, Fig. 7 |
| C.6 | Safe-speed profile along the ramp | Sec. 3.6, Fig. 8 |
| C.7, C.7b | Sensitivity of the SRT to roll centre, payload height, off-centre load and roll-stiffness ranges | Sec. 3.3 |
| C.8 | Comparison with the ORNL track tests | Discussion |
| D.2 | Friction law: power law versus 1987 Magic Formula peak, effect on the safe speed | Secs. 2.7, 3.7, Fig. 9 |
| D.3 (second half), D.4 | Supplementary slip-angle analysis, **not reported in the article** | — |

## Data sources

All input data come from published sources cited in the article: vehicle and friction parameters from
Yoo et al. (2025); ramp geometry and TruckSim results from Zhang et al. (2024); measured suspension roll
stiffness and roll-centre heights from Winkler, Karamihas & Bogard (1992) and Fu & Cebon (2002); track-test
lift-off accelerations from Oak Ridge National Laboratory (2009). The truck-tire friction points were digitized
from Fig. 5 of Yoo et al. (2025) (original data of Segel & Ervin). No proprietary data are used.

## Model C (PyChrono)

The multibody model is not needed to run the main script: its results are stored in
`multibody_pychrono/resultados/` and copied into Part B.3. To rerun it, follow
[`multibody_pychrono/README.md`](multibody_pychrono/README.md) (requires a conda environment with PyChrono).

## How to cite

Please cite the article above. To cite this code specifically, use the metadata in `CITATION.cff`
(GitHub shows a "Cite this repository" button) or the Zenodo DOI.

## Acknowledgments

Parts of the code and of its documentation were developed with the assistance of Claude (Anthropic), as declared
in the article. The authors reviewed and verified all code and results.

## License

MIT, see [`LICENSE`](LICENSE).
