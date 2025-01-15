rm(list=ls())

#### Check data aesthetic and public interest ####
# data
load(file = here::here("data", "raw_data", "RLS_actinopterygii_data.Rdata"))
load(file = here::here("data", "raw_data", "environmental_covariates", "TEMP_list_survey_sst.Rdata"))
load(file = here::here("outputs","2g_cultural_all_species.Rdata"))


# load( file = here::here("data", "derived_data", "RLS_species_traits_contrib.Rdata"))
# load(file = here::here("data", "raw_data", "RLS_species_traits.Rdata"))


#list of fish in tropical reefs
rls_trop <- RLS_actinopterygii_data |> 
  dplyr::filter(survey_id %in% list_survey_sst$survey_id) |> 
  dplyr::select(survey_id, rls_species_name, spec_code) |> 
  unique()

species_trop <- rls_trop |>
  dplyr::select(rls_species_name, spec_code) |> 
  unique()

#traits of theses species / Remove pleuro and syngnatiformes
traits_sp_trop <- dplyr::select(species_trop, rls_species_name) |> 
  dplyr::left_join(cultural_all_species) |> 
  dplyr::filter(!order %in% c("Pleuronectiformes", "Syngnathiformes"))

# Missing data
missing_interest_data <- traits_sp_trop |> 
  dplyr::select(rls_species_name:order, family:fishbase_name, public_attention, academic_knowledge) |> 
  dplyr::filter(is.na(public_attention))

missing_aesthetic_data <- traits_sp_trop |> 
  dplyr::select(rls_species_name:order, family:fishbase_name, esthe_score)|> 
  dplyr::filter(is.na(esthe_score))

# Remove species according to Mouquet 2024
to_remove <- c("Elagatis bipinnulata", "Euthynnus affinis", "Euthynnus lineatus",
              "Gadus morhua","Mola mola", "Salmo salar", "Thunnus albacares")

missing_interest_data <- missing_interest_data |> 
  dplyr::filter(!fishbase_name %in% to_remove)

#save sp list
write.csv(missing_aesthetic_data, here::here("data","list_species_missing_aesthetic_data.csv"),
          row.names = F)
write.csv(missing_interest_data, here::here("data","list_species_missing_interest_data.csv"),
          row.names = F)





#### Photoquadrats ####
metadata_survey <- read.csv(here::here("data", "raw_data","ep_survey_list.csv"))
all_survey <- metadata_survey 
# |> 
#   dplyr::filter(has_pq_scores_in_db == "true" | has_pqs_catalogued_in_db == "true")

#liste des anciens surveys avec PQ: 
load("/home/u_flandrin/Bureau/RLS_data/data/derived_data/covariates/raw_covariates_all_surveys.Rdata")
old_pq  <- raw_covariates_all |> 
  dplyr::filter(!is.na(coral)) |> 
  dplyr::select(survey_id)

#liste des surveys actuels
studied_survey <- list_survey_sst$survey_id

# PQ data ?
pq_data <- all_survey |> 
  dplyr::filter(survey_id %in% studied_survey) |>  # all survey with PQ data
  dplyr::filter(!survey_id %in% old_pq$survey_id)

# => PAS DE NOUVEAU PQ DATA DEPUIS DEBUT 2023


## Proportion de surveys finalement sélectionné sans PQ:
load("/home/u_flandrin/Bureau/RLS_data/data/derived_data/covariates/raw_covariates_all_surveys.Rdata")
load(file = here::here("data", "derived_data", "3_all_covariates_to_predict.Rdata"))

raw_cov <- raw_covariates_all |> 
  dplyr::filter(survey_id %in% rownames(covariates_final))
summary(raw_cov$coral)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.3027  0.5348  0.6424  0.8254  1.5708    1779 
1779/4829
