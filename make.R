#' Run the Entire Project
#'
#' 
#' @description 
#' This script runs the entire project and produces all figures present in the
#' Flandrin U. _et al._ 2025.
#' 
#' @author Ulysse Flandrin \email{ulysse.flandrin@gmail.com}
#' 
#' @date 25-02-2025


#-----------------Install required dependencies---------------------
if (!("remotes" %in% utils::installed.packages()))
  install.packages("remotes")
if (!("renv" %in% utils::installed.packages()))
  install.packages("renv")
renv::install("r-lib/devtools") #to install packages from github
renv::install() #install all packages noted in file DESCRIPTION



#-----------------Run the project---------------------
### 1) Extract and complete species traits ###

source(here::here('analyses', '1_species_traits_and_contributions', '1a_extract_species_contributions.R'))
source(here::here('analyses', '1_species_traits_and_contributions', '1b_filter_tropical_data.R'))
source(here::here('analyses', '1_species_traits_and_contributions', '1c_infer_species_traits.R'))




### 2) Assess all fish contributions at the community level ###

source(here::here('analyses', '2_contributions_of_RLS_surveys', '2a_PA_and_abundance_matrix.R'))
source(here::here('analyses', '2_contributions_of_RLS_surveys', '2b_biodiversity_indices.R')) # /!\ long time
source(here::here('analyses', '2_contributions_of_RLS_surveys', '2c_biomass_distribution.R')) # /!\ long time
source(here::here('analyses', '2_contributions_of_RLS_surveys', '2d_biochemical_indices.R')) # /!\ long time
source(here::here('analyses', '2_contributions_of_RLS_surveys', '2e_Food_Web_indices.R'))
source(here::here('analyses', '2_contributions_of_RLS_surveys', '2f_food_intake.R'))
source(here::here('analyses', '2_contributions_of_RLS_surveys', '2g_cultural_contributions.R'))

## Merge all contributions in a unique table
source(here::here('analyses', '2_contributions_of_RLS_surveys', '2h_merge_contributions.R'))




### 3) Predict contributions and run couterfactual scenarios ###

## Clean the contributions and covariates data
source(here::here('analyses', '3_predict_contributions', '3a_preping_data_to_predict.R'))

## Fit Bayesian models
# -> need to set-up the HMSC-HPC package. See the file HMSC_package/Set_up_HMSC_package.Rmd
# /!\ run all models will take several days -> open the file and select the wanted model
source(here::here('analyses', '3_predict_contributions', '3b_fit_model_hmsc.R')) 

## Plot model outputs (Fig. 1 and Fig. 2)
# Plot all figures of the Full model. Run line 'file_name <- gsub("output_", "", list_files[5])' to see another model.
source(here::here('analyses', '3_predict_contributions', '3c_results_hmsc.R'))

## Run and plot counterfactual scenarios (Fig. 3 and Fig. 4)
source(here::here('analyses', '3_predict_contributions', '3d_conterfactuals.R'))
source(here::here('analyses', '3_predict_contributions', '3e_sensibility_analysis.R'))

