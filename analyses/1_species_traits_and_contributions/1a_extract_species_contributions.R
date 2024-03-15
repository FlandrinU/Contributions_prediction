################################################################################
##
##  
##
## 1a_extract_species_contributions.R
##
## 10/01/2024
##
## Ulysse Flandrin
##
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("here", "dplyr")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

##-------------loading data -------------
#All RLS species
load(file = here::here("data", "raw_data", "RLS_species_traits.Rdata"))

#Nutrients data
load(file = here::here("data", "raw_data", "data_rfishbase.RData")) #/!\ large file

#Aesthethic and interests data
aesth_interest_table <- read.csv(here::here("data", "raw_data", "Mouquet2023_Human_Interest_final_table.csv"),
                                 sep=";", dec=",")

#Recycling data -> see project "Tropical reef contributions" for assessement details
load(file = here::here("data", "raw_data", "flux_final_data_species.Rdata"))


##----------------------------Nutrients supply----------------------------
variables <- colnames(data_rfishbase)
variables[order(variables)]

nutrient_list <- c("SpecCode", "Species",
                   "Calcium", "Iron", "Omega3", "Selenium", "VitaminA", "Zinc")


nutrients <- data_rfishbase |> 
  dplyr::select(dplyr::all_of(nutrient_list)) |> 
  dplyr::filter(SpecCode %in% species_traits_final$spec_code) |> 
  dplyr::distinct() |> 
  dplyr::rename(spec_code = SpecCode,
                fishbase_name = Species)

# Clear environnement
rm("data_rfishbase")

# Merge data
species_nutrients <- species_traits_final |> 
  dplyr::left_join(nutrients)

# # Save data 
# save(species_nutrients, file = here::here("data", "derived_data", "2a_species_nutrients.Rdata"))
# # load(file = here::here("data", "derived_data", "2a_species_nutrients.Rdata"))


##---------------------------Aesthetic and interests----------------------------
aesth_interest <- dplyr::select(aesth_interest_table,
                                fb_sci_name, public, acad, Aesthetic) |> 
  dplyr::mutate(fb_sci_name = gsub("_", " ", fb_sci_name)) |> 
  dplyr::rename(fishbase_name = fb_sci_name,
                public_interest = public,
                academic_knowledge = acad,
                aesthetic = Aesthetic)


# Merge data
species_cultural_contrib <- species_nutrients |> 
  dplyr::left_join(aesth_interest) 



# # Save data
# save(species_cultural_contrib, file = here::here("data", "derived_data", "2b_species_cultural_contrib.Rdata"))
# # load( file = here::here("data", "derived_data", "2b_species_cultural_contrib.Rdata"))




##---------------------------Recycling data----------------------------
recycling <- flux_sp |> dplyr::rename(rls_species_name = species)
recycling$rls_species_name <- gsub("_", " ", recycling$rls_species_name)

# Merge data
species_traits_contrib <- species_cultural_contrib |> 
  dplyr::left_join(recycling) 



##---------------------------Explore data----------------------------
library(funbiogeo)
species_traits <- dplyr::rename(species_traits_contrib,
                                species = rls_species_name) |> 
  dplyr::select(-phylum, -class, -order, -family, -spec_code, -worms_id)

fb_plot_species_traits_completeness(species_traits)
fb_plot_number_species_by_trait(species_traits, threshold_species_proportion = 0.75)

species_traits_numeric <- species_traits |> 
  dplyr::select(-fishbase_name, -IUCN_category, -Importance, -PriceCateg, 
                -UsedforAquaculture, -Aquarium, -Schooling, -DemersPelag,
                -trophic_guild, -IUCN_inferred_Loiseau23 )
fb_plot_trait_correlation(species_traits_numeric)


##---------------------------Save data----------------------------
save(species_traits_contrib, file = here::here("data", "derived_data", 
                                               "RLS_species_traits_contrib.Rdata"))
# load( file = here::here("data", "derived_data", "RLS_species_traits_contrib.Rdata"))
