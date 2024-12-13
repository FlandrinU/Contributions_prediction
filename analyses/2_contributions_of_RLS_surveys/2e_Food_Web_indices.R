#################################################################################
#'
#'This script infers the trophic interaction metaweb among all the RLS species
#' and extract different trophic indicators of the local webs, defined at the 
#' survey scale
#' 
#'
#'@author Ulysse Flandrin, \email{ulysse.flandrin@@gmail.com}
#'
#' @date 2024/02/16
#'
###############################################################################"
##-----------------cleaning memory-------------------
rm(list=ls())


##-----------------Loading packages-------------------
# pkgs <- c("here", "tidyverse", "rfishbase", "dplyr", "DHARMa", "gtools",
#           "parallel", "car", "PresenceAbsence",  "igraph", "NetIndices",
#           "ggplot2", "ggsignif", "ggraph" )
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))


## ----------------- Loading data -------------------
#Species traits
load(file= here::here("outputs", "RLS_species_traits_inferred.Rdata"))

#Occurence and biomass matrices
load(file=here::here("data", "derived_data", "2_relative_biom_matrix_sp_survey.Rdata"))


##----------------- Metaweb construction ---------------------####
# Construct the metaweb from species size (and Trophic level) from an allometric
#  niche model (Albouy et al. 2019), calibrated with observed interactions from 
# Barnes et al. 2008 data
source(here::here("R", "metaweb_construction.R")) 


##----------------- cross validation niche model ---------------------####
# Test the confidence of the allometric niche model with a Boyce index
# source(here::here("R", "cross_validation_niche_model.R")) #OK but long time to run


##----------------- Metaweb analysis ---------------------####
# Figures of trophic web analysis -> not necessary, here to evaluate the metaweb
source(here::here("R", "metaweb_figures.R")) 


##----------------- extraction of local web indices ---------------------####
# Extract local trophic web and trophic indicators of each sites
source(here::here("R", "extract_local_web_trophic_indicators.R")) 



##----------------- Assess the weighted mTL from fishbase's Troph ---------------------####
#reorder the species 
traits <- inferred_species_traits |> 
  dplyr::filter(!class == "Elasmobranchii") |> 
  dplyr::select(Troph, Length) |> 
  tidyr::drop_na() |> 
  as.matrix()

# traits <- as.matrix(traits[colnames(surveys_sp_pbiom),])
surveys_sp_pbiom <- as.matrix(surveys_sp_pbiom[,rownames(traits)]) #reorder

which(rownames(traits) != colnames(surveys_sp_pbiom)) #OK

#Mean trophic level: Trophic level weihted by the relative biomass of each species
troph_mTL <- as.matrix(surveys_sp_pbiom) %*% traits |> #MATRICIAL PRODUCT
  as.data.frame() |> 
  dplyr::select(troph_mTL = Troph) |> 
  tibble::rownames_to_column("survey_id")



##----------------- observation ---------------------####
load(file= here::here("outputs", "2e_trophic_indicators_survey.Rdata"))
food_web_indices <- trophic_indicators_survey |> 
  dplyr::select(survey_id, Species, Connectance, b_power_law, weighted_mTL) |> 
  dplyr::mutate(survey_id = as.character(survey_id)) |> 
  dplyr::full_join(troph_mTL)


# Check distributions
library(ggplot2)
ggplot(data=tidyr::pivot_longer(food_web_indices,
                                cols = -survey_id,
                                names_to = "index", values_to = "values"), 
       aes(x=values, group=index, fill=index)) +
  geom_histogram(aes(y = ..density..), bins = 20, color = "grey40", fill ="white") +
  geom_density(aes(fill = index), alpha = 0.2) +
  hrbrthemes::theme_ipsum() +
  facet_wrap(~index, scales = "free") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )
#

plot(food_web_indices$b_power_law~food_web_indices$Species)
plot(food_web_indices$weighted_mTL~food_web_indices$Species)
plot(food_web_indices$troph_mTL~food_web_indices$Species)
plot(food_web_indices$weighted_mTL~food_web_indices$troph_mTL)



##----------------- save data ---------------------####

save(food_web_indices, file = here::here("outputs", "2e_food_web_indices_surveys.Rdata" ))
# load( file = here::here("outputs", "2e_food_web_indices_surveys.Rdata" ))
