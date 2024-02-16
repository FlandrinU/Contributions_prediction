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
#           "ggplot2", "ggsignif" )
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

## ----------------- Loading data -------------------
#Species traits
load(file= here::here("outputs", "RLS_species_traits_inferred.Rdata"))


##----------------- Metaweb construction ---------------------####
# Construct the metaweb from species size (and Trophic level) from an allometric
#  niche model (Albouy et al. 2019), calibrated with observed interactions from 
# Barnes et al. 2008 data
source(here::here("R", "metaweb_construction.R")) 


##----------------- Species data ---------------------####
# Test the confidence of the size niche model with a Boyce index
source(here::here("trophic_web", "R", "cross_validation_niche_model.R")) 

##----------------- Species data ---------------------####
# Extract local trophic web and trophic indicators of each sites
source(here::here("trophic_web", "R", "extract_local_web_trophic_indicators.R")) 

##----------------- Species data ---------------------####
# Figures of trophic web analysis -> not necessary, here to evaluate the metaweb
source(here::here("trophic_web", "R", "metaweb_figures.R")) 
