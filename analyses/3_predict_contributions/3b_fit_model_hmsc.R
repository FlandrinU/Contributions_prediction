################################################################################
##
##  
##
## 3c_fit_model_hmsc.R
##
## 07/05/2024
##
## Ulysse Flandrin
##
################################################################################
##---------------------------- cleaning memory ---------------------------------
rm(list=ls())

##----------------------Loading packages and functions--------------------------
pkgs <- c("here", "Hmsc", "coda", "ggmcmc")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

source("R/HMSC_function.R")

##------------------------------- load data ------------------------------------
## Survey scale, all covariates
load(here::here("data/derived_data/3_all_contributions_to_predict.Rdata"))
load(here::here("data/derived_data/3_all_covariates_to_predict.Rdata"))
load( here::here("data", "derived_data", "3_datasets_for_predict_CV_80_20.Rdata"))
Y_data =  observations_final
X_data = covariates_final[rownames(Y_data),]

## Without Allen
load(here::here("data/derived_data/3_contributions_without_Allen_to_predict.Rdata"))
load(here::here("data/derived_data/3_covariates_without_Allen_to_predict.Rdata"))
Y_data_wo_allen =  observations_final_without_Allen
X_data_wo_allen = covariates_final_without_Allen[rownames(Y_data_wo_allen),]

## Site scale
load(here::here("data/derived_data/3_sites_contributions_to_predict.Rdata"))
load(here::here("data/derived_data/3_sites_covariates_to_predict.Rdata"))
Y_data_site =  observations_site_final
X_data_site = covariates_site_final[rownames(Y_data_site),]


## Only Australia
X_data_aust = covariates_final |> dplyr::filter(country == "Australia")
Y_data_aust =  observations_final[rownames(X_data_aust),]


##----------------------------- Set-up parameters ------------------------------
nSamples = 100
thin = 2000
nChains = 2 
verbose = 1000 
nb_neighbours = 10
transient = nSamples * thin

set_shrink = NULL
response_distribution <- rep("normal", ncol(Y_data))
#response_distribution[colnames(Y_data) == "iucn_species_richness"] <- "poisson"

quadratic_effects = NULL
# quadratic_effects = colnames(
#   dplyr::select(X_data,-longitude,-latitude,-year,-country, -ecoregion,-realm,
#                 -hdi,-marine_ecosystem_dependency,-ngo,-natural_ressource_rent))


save_path = here::here("outputs/models/hmsc")

##----------------------------- Fit HMSC models ------------------------------

#### FULL MODEL SURVEYS ####
name = "FULL_survey_site_country_in_rL"
random_factors = c("sample_unit", "site", "country")


hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data,
              X_data,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data,
                         X_data = X_data,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, name,
                         run_python = T, save_path)


#### FULL MODEL SITES ####
name = "FULL_SITE_SCALE_site_country_in_rL"
random_factors = c("sample_unit", "country")

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data_site,
                         X_data = X_data_site,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, name,
                         run_python = T, save_path)



#### FULL MODEL SITES, QUADRATIC EFFECTS ####
quadratic_effects_all = colnames(
  dplyr::select(X_data,-longitude,-latitude,-year,-country, -ecoregion,-realm,
                -hdi,-marine_ecosystem_dependency,-ngo,-natural_ressource_rent))

name = "quadratic_SITE_SCALE_site_country_in_rL"

hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, 
              quadratic_effects = quadratic_effects_all,
              random_factors, nb_neighbours, set_shrink, name,
              run_python = T, save_path)







# ##--------------------------- Sensitivity analyses -----------------------------
#  source("R/evaluation_prediction_model.R")
# ##### Quadratic effects ####
# full_data <- cbind(X_data, Y_data) |> 
#   tidyr::pivot_longer(cols = colnames(Y_data), names_to = "contributions")
# 
# plot_interaction(full_data,
#                  var_facet_wrap = "contributions",
#                  X_values = "median_5year_analysed_sst",
#                  Y_values = "value",
#                  add_curve = T) # (mgcv::gam() is used with formula = y ~ s(x, bs = "cs") with method = "REML".)
# 
# ggsave(filename = here::here("figures/models/hmsc/sensitivity_analyses/contributions_VS_SST.jpg"),
#        width = 15, height = 8)
# 
# 
# 
# 
# ##### fit on several realms ####
# 
# nSamples = 300 #1000 
# thin = 3000 #100
# nChains = 3 
# verbose = 100 
# nb_neighbours = 10
# 
# random_factors = c("sample_unit","country")
# response_distribution <- rep("normal", ncol(Y_data))
# 
# 
# # Run on the realms separately
# realms <- unique(X_data$realm)
# 
# for(rlm in realms){
#   cat("Fit model for realm ", rlm, "\n")
#   
#   X_data_realm <- dplyr::filter(X_data, realm == rlm)
#   Y_data_realm <- Y_data[rownames(X_data_realm),]
# 
#   name = paste0("test_REALM_", gsub(" ", "_", rlm), "_asso_and_country_in_rL")
#   
#   hmsc_function(nSamples,
#                 thin,
#                 nChains,
#                 verbose,
#                 transient = nSamples * thin,
#                 Y_data_realm,
#                 X_data_realm,
#                 response_distribution,
#                 random_factors,
#                 nb_neighbours,
#                 name,
#                 run_python = TRUE,
#                 save_path = here::here("outputs/models/hmsc/sensitivity_analysis"))
#   }




