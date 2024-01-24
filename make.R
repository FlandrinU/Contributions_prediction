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
source(here::here("analyses/1_species_traits_and_contributions/1b_complete_tropical_species_traits.R"))

## 2) Assess contributions at the survey level

