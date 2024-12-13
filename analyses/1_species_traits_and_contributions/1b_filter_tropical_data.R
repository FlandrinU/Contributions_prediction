################################################################################
##
##  
##
## 1b_filter_tropical_data.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("here", "dplyr", "ggplot2", "funbiogeo", "patchwork","rnaturalearth",
#           "tibble", "stringr")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

library(funbiogeo)
library(ggplot2)

##-------------loading data and functions-------------
#Species traits
load( file = here::here("data", "derived_data", "RLS_species_traits_contrib.Rdata"))

#RLS observations
load( file = here::here("data", "raw_data", "RLS_actinopterygii_data.Rdata"))
load( file = here::here("data", "raw_data", "RLS_elasmobranchii_data.Rdata"))

#survey metadata
load(here::here("data", "raw_data", "environmental_covariates", 
                "all_covariates_benthos_inferred_tropical_surveys.Rdata"))

#coastline shapefile
worldmap <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')

# functions #
source(here::here("R","evaluation_prediction_model.R"))

##------------- Filtering tropical fishes -------------
### filtering to keep only RLS data for sub-tropical and tropical sites ######
##  i.e with min(SST)>=17°C over the 5 years before the survey

sites_upper_17SST <- all_covariates_benthos_inferred |>
   dplyr::filter(min_5year_analysed_sst >= 17) |> 
   dplyr::filter(latitude > -30) # Approximatelly the limit of the Allen Atlas ([-31.56, 29.5])

# check filter
NA_on_map(data=sites_upper_17SST, variable = "mean_1year_analysed_sst",
                      xlim = c(-180,180),ylim = c(-50,50),
                      jitter = 0, lat_line = 30, 
                      priority= "no") #which points are display on the others


 ### filter observations ###
 rls_actino_trop <- RLS_actinopterygii_data |>
   dplyr::filter(survey_id %in% sites_upper_17SST$survey_id) 
 
 # diversity remaining
 dplyr::n_distinct(rls_actino_trop$survey_id) # 6936 surveys
 dplyr::n_distinct(rls_actino_trop$site_code) # 2194 sites
 dplyr::n_distinct(rls_actino_trop$rls_species_name) # 1655 taxa
 dplyr::n_distinct(rls_actino_trop$family) # 84 families


 ### Elasmobranch observations ###
 rls_elasmo_trop <- RLS_elamsobranchii_data |>
   dplyr::filter(survey_id %in% unique(rls_actino_trop$survey_id))
 
 dplyr::n_distinct(rls_elasmo_trop$rls_species_name) # 60 taxa
 dplyr::n_distinct(rls_elasmo_trop$family) # 16 families
 


##-------------filter observed species in tropical reefs (Actino + Elasmo) -------------
observed_species <- dplyr::distinct(
 rbind(dplyr::select(rls_actino_trop, -raw_biomass),
                     rls_elasmo_trop), 
 rls_species_name, .keep_all = TRUE)

tropical_species_traits <- species_traits_contrib |> 
  dplyr::filter(rls_species_name %in% observed_species$rls_species_name)

# diversity remaining
dplyr::n_distinct(tropical_species_traits$family) # 99 families (84 actino)
dplyr::n_distinct(tropical_species_traits$rls_species_name) # 1715 taxa (1655 actino)
dplyr::n_distinct(tropical_species_traits$fishbase_name) # 1695 taxa -> 20 duplicates
dplyr::n_distinct(tropical_species_traits$fishbase_name[
   tropical_species_traits$class != "Elasmobranchii"]) # 1637 taxa
dplyr::n_distinct(tropical_species_traits$fishbase_name[
   tropical_species_traits$class == "Elasmobranchii"]) # 58 taxa

##-------------Observe data with missing values-------------
# Explore data
species_traits <- dplyr::rename(tropical_species_traits,
                                species = rls_species_name) |> 
  dplyr::select(-phylum, -class, -order, -family, -spec_code, -worms_id)

fb_plot_species_traits_completeness(species_traits)
ggsave(plot= last_plot(), file= here::here("figures", "1_traits_completedness_tropical.png"), width = 15, height = 7)
fb_plot_number_species_by_trait(species_traits, threshold_species_proportion = 1)
ggsave(plot= last_plot(), file= here::here("figures", "1_percent_species_per_traits_tropical.png"), width = 8, height = 8)

# ### espèces manquantes pour l'esthétique
# sp_missing_aest <- tropical_species_traits[is.na(tropical_species_traits$aesthetic),]
# added_sp_aest <- readr::read_csv(here::here("outputs/2g_aesthetic_inference_new_sp.csv")) |> 
#    dplyr::mutate(name = gsub("([A-Za-z]+)\\_([A-Za-z]+).*", "\\1 \\2", image_name))
# a_faire <- dplyr::filter(sp_missing_aest, !fishbase_name %in% added_sp_aest$name) |> 
#    dplyr::filter(class != "Elasmobranchii",
#                  order != "Pleuronectiformes",
#                  order != "Syngnathiformes")
# readr::write_csv(a_faire, file="~/Téléchargements/missing_species_aesthetic.csv")
# ###
# ### espèces manquantes pour le public attention
# sp_missing_hum_int <- tropical_species_traits[is.na(tropical_species_traits$public_interest),]
# readr::write_csv(sp_missing_hum_int, file="~/Téléchargements/missing_species_human_interst.csv")
# ###


##-------------Save tropical data-------------
# RLS data #
save(rls_actino_trop, file = here::here("data/derived_data/rls_actino_trop.Rdata"))
save(rls_elasmo_trop, file = here::here("data/derived_data/rls_elasmo_trop.Rdata"))

# load(file = here::here("data/derived_data/rls_actino_trop.Rdata"))
# load(file = here::here("data/derived_data/rls_elasmo_trop.Rdata"))

# species traits #
save(tropical_species_traits, file = here::here("data/derived_data/tropical_species_traits.Rdata"))

