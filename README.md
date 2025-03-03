
# Human impacts and conservation legacy on the fish community contributions

<!-- badges: start -->

[![License: GPL (\>=
2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)
<!-- badges: end -->

Research Compendium of the project **Modelling the reef fish contributions, human impacts, and the legacy of conservation**

### How to cite

Please cite this compendium as:

> **Flandrin et al. 20XX**

### Content

This repository is structured as follow:

- [`data/raw_data`](https://github.com/FlandrinU/Contributions_prediction/tree/master/data/raw_data):
  contains all raw data required to perform analyses, do NOT modify this folder.
  
- [`data/derived_data`](https://github.com/FlandrinU/Contributions_prediction/tree/master/data/derived_data):
  contains all intermediate results for further analysis.

- [`analyses/`](https://github.com/FlandrinU/Contributions_prediction/tree/master/analyses/):
  contains R scripts to run each step of the workflow, structured in 3 parts.

- [`outputs/`](https://github.com/FlandrinU/Contributions_prediction/tree/master/outputs):
  contains all the results created during the workflow

- [`figures/`](https://github.com/FlandrinU/Contributions_prediction/tree/master/figures):
  contains all the figures created during the workflow, structured following the analysis parts

- [`R/`](https://github.com/FlandrinU/Contributions_prediction/tree/master/R): contains
  R functions developed especially for this project, and called in the analysis scripts.

- [`DESCRIPTION`](https://github.com/FlandrinU/Contributions_prediction/tree/master/DESCRIPTION):
  contains project metadata (author, date, dependencies, etc.), lists all the needed 
  packages to during the workflow

- [`make.R`](https://github.com/FlandrinU/Contributions_prediction/tree/master/make.R):
  main R script to run the entire project by calling each R script
  stored in the `analyses/` folder. Be careful, this script would take several
  days to run and need important computational ressources



### Workflow

#### 1) Species traits and contributions

In `analyses/1_species_traits_and_contributions` folder, the scripts `1a_` to `1c_` extract the needed species traits and species-level contributions to assess fish contributions, and infer of the missing data with Random Forest.

#### 2) Fish community contributions

In `analyses/2_contributions_of_RLS_surveys` folder, the scripts `2a_` to `2h_` assess the fish contributions at the survey scale, and aggregate all contributions at the site level, following the workflow detailed in Flandrin et al. 2024 (10.1016/j.oneear.2024.09.011).

#### 3) Modelling the contributions

In `analyses/3_predict_contributions` folder, the scripts `3a_` to `3e_` use covariates to model the reef fish contributions at the site scale, using the Bayesian framework 'HMSC', to investigate the drivers of community contributions. `3d_` runs countercatual scenarios to deepen our understanding of the human footprint on ecosystems, and explore the potential of conservation efforts. `3d_` and `3e` reproduce all figures presented in the paper Flandrin et al. associated to this project.



### Usage

Clone the repository, open R/RStudio and run:

``` r
source("make.R")
```

### Notes

- All required packages, listed in the `DESCRIPTION` file, will be
  installed (if necessary)
- All required packages and R functions will be loaded
- Some analyses listed in the `make.R` might take time and resources. You should 
run lines 13-19 of make.R and then go to the other script. Note that the lines 35-44 let to assess each contributions in each surveys; you can skip this part and run the lines 52-65 to reproduce all the figures of the paper Flandrin et al.
