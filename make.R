#' template: A Research Compendium
#' 
#' @description 
#' A paragraph providing a full description of the project and describing each 
#' step of the workflow.
#' 
#' @author Ulysse Flandrin \email{ulysse.flandrin@gmail.com}
#' 
#' @date 2024

#-----------------Install required dependencies---------------------
# # renv::restore()
# # renv::init()
# # renv::install()
# # renv::status()
# # renv::snapshot()
# 
# if (!("remotes" %in% utils::installed.packages())) 
#   install.packages("remotes")
# if (!("renv" %in% utils::installed.packages())) 
#   install.packages("renv")
# renv::install("r-lib/devtools") #to install packages from github
# renv::install() #install all packages noted in file DESCRIPTION
# renv::snapshot()


## Load packages & functions + Setup project ----

devtools::load_all(here::here())



#-----------------Run the project---------------------
## 1) Complete species traits (from github project 'RLS_data')
source(here::here("analyses/1_species_traits_and_contributions/1a_extract_species_contributions.R"))
source(here::here("analyses/1_species_traits_and_contributions/1b_filter_tropical_data.R"))
source(here::here("analyses/1_species_traits_and_contributions/1c_infer_species_traits.R"))

## 2) Assess contributions at the survey level
source(here::here("analyses/2_contributions_of_RLS_surveys/2a_PA_and_abundance_matrix.R")) #OK
source(here::here("analyses/2_contributions_of_RLS_surveys/2b_biodiversity_indices.R")) #OK /!\ long time
source(here::here("analyses/2_contributions_of_RLS_surveys/2c_biomass_distribution.R")) # OK /!\ long time
source(here::here("analyses/2_contributions_of_RLS_surveys/2d_biochemical_indices.R")) # OK /!\ long time
source(here::here("analyses/2_contributions_of_RLS_surveys/2e_Food_Web_indices.R")) # OK
source(here::here("analyses/2_contributions_of_RLS_surveys/2f_food_intake.R")) # OK
source(here::here("analyses/2_contributions_of_RLS_surveys/2g_cultural_contributions.R"))
source(here::here("analyses/2_contributions_of_RLS_surveys/2h_merge_contributions.R"))

## 3) Predict contributions and run couterfactual scenarios
source(here::here("analyses/3_predict_contributions/3a_preping_data_to_predict.R"))
#-> need to set-up the HMSC-HPC package. See the file HMSC_package/Set_up_HMSC_package.Rmd
source(here::here("analyses/3_predict_contributions/3b_fit_model_hmsc.R")) #/!\ run all models will take several days..

