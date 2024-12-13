################################################################################
##
## Script to determine the list of species in each surveys of RLS data  
## (occurrence matrix) and their relative biomass
##
## 2a_PA_and_abundance_matrix.R
##
## 23/01/2024
##
## Ulysse Flandrin
##
################################################################################

## cleaning memory
rm(list=ls())

##------------------- loading datasets-------------------
#RLS observations
load(file = here::here("data/derived_data/rls_actino_trop.Rdata"))


# ------------------- total biomass and abundance in each survey  -------------------
rls_actino_trop <- rls_actino_trop|>
  dplyr::mutate(biomass = raw_biomass) |> #IN THIS PROJECT, WE CONSIDER THE BIOMASS
  # OF FISHES ESTIMATES AS a*Size^b , WITH THE SIZE REPORTED BY DIVERS, AND NOT CORRECTED.
  dplyr::group_by(survey_id) |>
  dplyr::mutate(
    abundance_tot_survey = sum(total), #total abundance in the survey
    biomass_tot_survey = sum(biomass))  #total biomass in the survey

# ------------------- biomass of species in each survey  -------------------
surveys_sp_biom <- rls_actino_trop |> 
  dplyr::select(survey_id, rls_species_name, biomass) |>
  dplyr::group_by(survey_id, rls_species_name) |>
  dplyr::summarize( sp_biom = sum(biomass) ) |>
  tidyr::pivot_wider(names_from = rls_species_name, values_from = sp_biom, values_fill = 0) |>
  tibble::column_to_rownames(var="survey_id") 


dim(surveys_sp_biom) # OK: 1655 actino according RLS -> 1637 unique species


# ------------------- occurrence matrix of species in surveys  -------------------
surveys_sp_occ <- surveys_sp_biom
surveys_sp_occ[surveys_sp_occ!=0] <- 1

# ------------------- relative biomass of species in surveys  -------------------
surveys_sp_pbiom <- surveys_sp_biom / apply(surveys_sp_biom,1,sum)

# ------------------- save matrices  -------------------
save(rls_actino_trop, file=here::here("data", "derived_data", "2_rls_actino_trop.Rdata"))
save(surveys_sp_occ, file=here::here("data", "derived_data", "2_occurrence_matrix_sp_survey.Rdata"))
save(surveys_sp_pbiom, file=here::here("data", "derived_data", "2_relative_biom_matrix_sp_survey.Rdata"))
save(surveys_sp_biom, file=here::here("data", "derived_data", "2_biom_matrix_sp_survey.Rdata"))


