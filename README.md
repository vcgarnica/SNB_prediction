# Hierarchical Bayesian generalized linear models for predicting Stagonospora nodorum blotch in winter wheat
## Part I: Model development and evaluation

## Citation
The following source code accompanies our publication, currently under review in *Phytopathology* Journal:

```
@article{garnica2025bayesian,
  author = {Vinicius C. Garnica and Peter S. Ojiambo},
  title = {Hierarchical Bayesian generalized linear models for predicting Stagonospora nodorum blotch in winter wheat. Part I: Model development and evaluation},
  year = {2025},
  doi = {},
  journal = {Phytopathology}
}
```

## Introduction

This repository contains data and code for Bayesian hierarchical modeling of Stagonospora nodorum blotch (SNB) disease severity in winter wheat. 
The study develops predictive models that integrate agronomic management factors, cultivar resistance, and pre-anthesis weather variables while 
accounting for hierarchical data structures and environmental variability.

Nine Bayesian beta regression models with increasing complexity were fitted to disease severity data from commercial soft red winter 
wheat cultivars planted across 14 environments in North Carolina from 2021 to 2024. Models were implemented using `brms` in 
R with Stan backend, employing Hamiltonian Monte Carlo sampling with literature-informed priors.

The repository includes the complete analysis pipeline, detailing the R scripts and data required to reproduce the study. 
The folder structure is as follows:

```
SNB_bayesian/
├── code/
│   ├── 9_relationship_matrix.R
│   ├── 10_simulation_plots.R
│   └── Window pane/
│       ├── code/
│       ├── data/
│       ├── figures/
│       └── results/
├── data/
│   ├── anthesis.csv
│   ├── euclimat.RData
│   ├── fa1_library.Rdata
│   ├── fa2_library.Rdata
│   ├── fa3_library.Rdata
│   ├── locations.csv
│   ├── matrix_predictors.csv
│   ├── pre_anthesis_variables.Rdata
│   ├── SNB.Rdata
│   └── weather_hour.Rdata
├── results/
│   ├── summary_metrics.csv
│   ├── weather_coefficients.csv
│   ├── models/
│   │   ├── M1.rds
│   │   ├── M2.rds
│   │   ├── M3.rds
│   │   ├── M4.rds
│   │   ├── M5.rds
│   │   ├── M6.rds
│   │   ├── M7.rds
│   │   ├── M8.rds
│   │   └── M9.rds
│   └── figures/
│       ├── fig1.tiff
│       ├── fig2.tiff
│       ├── fig3.tiff
│       ├── fig4.tiff
│       ├── fig5.tiff
│       ├── figS1.tiff
│       ├── figS2.tiff
│       ├── figS3.tiff
│       ├── figS4.tiff
│       ├── figS5.tiff
│       ├── figS6.tiff
│       └── figS7.tiff
├── simulation/
│   ├── data_final_sev.RData
│   ├── K_obs.RData
│   ├── sev_c.RData
│   ├── sev_u.RData
│   ├── Simulation.Rmd
│   └── Simulation.html
├── Bayesian_SNB_model.Rmd
└── Bayesian_SNB_model.html
```

## Analysis Pipeline

The complete analysis workflow is fully described in `Bayesian_SNB_model.Rmd` and rendered in `Bayesian_SNB_model.html`. The pipeline integrates data preparation, model fitting, diagnostics, and visualization to develop hierarchical Bayesian models for SNB severity prediction. 

In the `code/Window pane/` folder, all pre-processing steps were conducted to generate the factor-analytic weather libraries (`fa1_library.Rdata`, `fa2_library.Rdata`, `fa3_library.Rdata`) and the final disease dataset (`SNB.Rdata`).

Nine Bayesian beta regression models (M1-M9) were fitted using the `brms` package with literature-informed priors and Hamiltonian Monte Carlo sampling. Fitted models are stored in `results/models/` as .rds files. Model comparison was performed using LOO-CV as well as data reserved from entire new environments (simulating deployment scenario), with results summarized in `results/summary_metrics.csv`. Weather variable coefficients and their importance rankings are provided in `results/weather_coefficients.csv`.

Publication-ready figures are available in `results/figures/`, including main figures (`fig1.tiff` through `fig5.tiff`) and supplementary figures (`figS1.tiff` through `figS7.tiff`).

### Additional Scripts

The `9_relationship_matrix.R` script generates correlation matrices and visualizations of relationships between weather variables and disease severity across environments.

The `10_simulation_plots.R` script creates visualizations for simulation studies validating model performance under various scenarios.

### Simulation Study

The `simulation/` directory contains a complete simulation study evaluating model performance:
- `Simulation.Rmd`: R Markdown document with simulation code and analysis
- `Simulation.html`: Rendered report of simulation results
- Supporting data files: `data_final_sev.RData`, `K_obs.RData`, `sev_c.RData`, `sev_u.RData`

### Window Pane Analysis

The `code/Window pane/` directory contains a self-contained sub-analysis extracting disease-associated weather variables 
through windowed temporal analysis, with its own code, data, figures, and results subdirectories. You can find further information in
[Garnica and Ojiambo 2025)](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2025.1637130/full) and the project's [repository](https://github.com/vcgarnica/SNB_window_pane).

## References

* Garnica, V. C., Shalizi, M. N., and Ojiambo, P. S. (2025). Performance and stability of winter wheat cultivars to Stagonospora nodorum blotch epidemics in multi-environment trials. Phytopathology.

* Garnica, V. C., and Ojiambo, P. S. (2025). Leveraging window-pane analysis with environmental factor loadings of genotype-by-environment interaction to identify high-resolution weather-based variables associated with plant disease. Frontiers in Plant Science, 16, 1637130.

* Hobbs, N. T., and Hooten, M. B. (2015). Bayesian models: a statistical primer for ecologists. Princeton University Press.

* Bürkner, P. C. (2017). brms: An R package for Bayesian multilevel models using Stan. Journal of Statistical Software, 80(1), 1-28.

* Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing, 27(5), 1413-1432.