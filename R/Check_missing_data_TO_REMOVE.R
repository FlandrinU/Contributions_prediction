rm(list=ls())

### Check data aesthetic and public interest ###
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



