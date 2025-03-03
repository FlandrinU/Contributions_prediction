################################################################################
##
##  This script gathers all information on species traits, from fishbase and 
##   other published sources, needed to assess species contributions
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

#aesthetic values of species from Langlois et al. 2022
aesth_2022 <- read.csv(
  "data/raw_data/aesthetic_deep_model_Langlois_2022/Langlois_el_al_2022.csv")

#Recycling data -> see project "Tropical reef contributions" for assessement details
load(file = here::here("data", "raw_data", "flux_final_data_species.Rdata"))

## Load functions
source("R/check_scientific_names.R")

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

### AESTHETIC  NEW INFERENCE ###

# # Install missing librairies for python
# reticulate::py_install("pytz")
# reticulate::py_install("sympy")
# reticulate::py_install("Pillow")
# reticulate::py_install(c("torch", "torchvision"))


#run aesthetic for new species
# reticulate::source_python("R/inference_aesthetic_score_Langlois2022.py")
#IF IMPOSSIBLE TO RUN RETICULATE: open a terminal in the R/ folder, and run the
# python code via bash: > python3 04_inference.py

## Observe new inference:
aesth_new <- read.csv("outputs/2g_aesthetic_inference_new_sp.csv") 

# Check inference with old scores of Langlois 2022
check <- aesth_new |> 
  dplyr::mutate(old = ifelse(grepl("Langlois22", image_name), 1, 0)) |> 
  dplyr::mutate(file_name = gsub(".png", "", image_name)) |> 
  dplyr::mutate(file_name = gsub("_Langlois22", "", file_name)) |> 
  dplyr::mutate(sp_name = gsub("^(.*?)_[A-Z]_[0-9]+$", "\\1", file_name)) |> 
  dplyr::left_join(aesth_2022)

plot(check$predicted_score~check$esthe_score)
abline(a=0, b=1)
cor.test(check$predicted_score,check$esthe_score)
# INFERED SCORES CONSISTENT FOR THE 55 PICURES IN COMMON -> WE CAN INFER NEW 
# SCORES FROM PICTURES AND ADD THEM TO THE FILE OF LANGLOIS ET AL. 2022

asthetic_inferred <- check |> 
  dplyr::filter(old == 0) |> 
  dplyr::select(file_name = image_name, sp_name, predicted_score) |> 
  dplyr::group_by(sp_name) |> 
  dplyr::mutate(esthe_score = max(predicted_score)) |> 
  dplyr::filter(esthe_score == predicted_score) |> 
  dplyr::select(-predicted_score)

## Merge all scores
aesthetic_score <- aesth_2022 |> 
  dplyr::select(file_name, sp_name, esthe_score) |> 
  dplyr::bind_rows(asthetic_inferred)


### PUBLIC ATTENTION ###

human_interest <- aesth_interest_table |>
  dplyr::select(public_attention = public,
                academic_knowledge = acad,
                sp_name = fb_sci_name)

### CHECK SCIENTIFIC NAMES ###

cultural_sp_contrib <- dplyr::full_join(aesthetic_score, human_interest) |> 
  dplyr::mutate(sp_name = gsub("_", " ", sp_name)) 

cultural_sp_contrib_checked <- code_sp_check(cultural_sp_contrib, original_name = 'sp_name',
                                     mc_cores = 15)

cultural_sp_contrib <- cultural_sp_contrib_checked |> 
  dplyr::select(-worms_id, -sp_name, -file_name, -check, -spec_code)|>
  dplyr::group_by(fishbase_name) |> 
  tidyr::fill(tidyr::everything(), .direction = 'updown') |> 
  dplyr::mutate(esthe_score = max(esthe_score),
                public_attention = max(public_attention),
                academic_knowledge = max(academic_knowledge)) |>  
  dplyr::rename("aesthetic" = esthe_score) |> 
  unique()



cultural_all_species <- species_nutrients |> 
  dplyr::left_join(cultural_sp_contrib) 

save(cultural_all_species, file = here::here("outputs","1a_cultural_all_species.Rdata"))



##---------------------------Recycling data----------------------------
recycling <- flux_sp |> dplyr::rename(rls_species_name = species)
recycling$rls_species_name <- gsub("_", " ", recycling$rls_species_name)

# Merge data
species_traits_contrib <- cultural_all_species |> 
  dplyr::left_join(recycling) 



##---------------------------Explore data----------------------------
## Check merge data ##
colnames(species_traits_contrib) # OK
species_traits_contrib$rls_species_name[
  duplicated(species_traits_contrib$rls_species_name)] #OK
  
length(unique(species_traits_contrib$rls_species_name)) # OK
length(unique(species_traits_contrib$fishbase_name)) # 33 synonyms


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
