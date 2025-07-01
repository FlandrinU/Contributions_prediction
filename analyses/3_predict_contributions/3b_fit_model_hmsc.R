################################################################################
##
##  This script uses the HMSC framework to fit a hierarchical Bayesian linear 
##   model on the reef fish contributions. It runs the full model and all 
##   alternative models used for sensitivity analyses.
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
# pkgs <- c("here", "Hmsc", "coda", "ggmcmc")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

source("R/HMSC_function.R")

##------------------------------- load data ------------------------------------

## Site scale
load(here::here("data/derived_data/3_sites_contributions_to_predict.Rdata"))
load(here::here("data/derived_data/3_sites_covariates_to_predict.Rdata"))
Y_data_site =  observations_site_final
X_data_site = covariates_site_final[rownames(Y_data_site),
                -which(colnames(covariates_site_final) == "protection_status_detailed")]


##----------------------------- time frame ------------------------------
ids <- rownames(X_data_site)
date <- stringr::word(ids, 2, sep = "_")
date <- date[order(date)]
head(date,1) ; tail(date,1) #"2006-09-02" -> "2024-08-22"

##----------------------------- Set-up parameters ------------------------------
nSamples = 200
thin = 1000
nChains = 4
verbose = 500 
nb_neighbours = 10
transient = nSamples * thin

set_shrink = NULL
response_distribution <- rep("normal", ncol(Y_data_site))

quadratic_effects = NULL
test_null_model = NULL

save_path = here::here("outputs/models/hmsc")

##----------------------------- Fit HMSC FULL model ------------------------------
#### FULL MODEL SITES ####
name = "FULL_model_SITE_SCALE"
random_factors = c("sample_unit", "country")

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data_site,
                         X_data = X_data_site,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model, name,
                         run_python = T, save_path)



##--------------------------- Sensitivity analyses -----------------------------
## Survey scale, all covariates
load(here::here("data/derived_data/3_all_contributions_to_predict.Rdata"))
load(here::here("data/derived_data/3_all_covariates_to_predict.Rdata"))
Y_data =  observations_final
X_data = covariates_final[rownames(Y_data),
                          -which(colnames(covariates_final) == "protection_status_detailed")]


## Without Allen
load(here::here("data/derived_data/3_contributions_without_Allen_to_predict.Rdata"))
load(here::here("data/derived_data/3_covariates_without_Allen_to_predict.Rdata"))
Y_data_wo_allen =  observations_final_without_Allen
X_data_wo_allen = covariates_final_without_Allen[rownames(Y_data_wo_allen),
                       -which(colnames(covariates_final_without_Allen) ==
                                "protection_status_detailed")]

## Only/Without Australia
X_data_aust = covariates_site_final |> 
  dplyr::filter(country == "Australia") |> 
  dplyr::select(-protection_status_detailed,
                -Natural_ressource_rent, -Marine_ecosystem_dependency,
                -HDI) # No variability inside Australia
Y_data_aust =  observations_site_final[rownames(X_data_aust),]
summary(X_data_aust)


X_data_no_aust = covariates_site_final |> 
  dplyr::filter(country != "Australia") |> 
  dplyr::select(-protection_status_detailed,
                -Patch_reefs)
Y_data_no_aust =  observations_site_final[rownames(X_data_no_aust),]

summary(X_data_no_aust)

## Covariates before benthic composition inference
load(file = here::here("data", "derived_data", 
                       "3_sites_without_NA_in_PQ_covariates_to_predict.Rdata"))
load(file = here::here("data", "derived_data",
                       "3_sites_without_NA_in_PQ_contributions_to_predict.Rdata"))

## Habitat covariates summarized into PCA dimensions
load(file = here::here("data", "derived_data",
                       "3_sites_summarized_PCA_covariates_to_predict.Rdata"))




#### FULL MODEL SURVEYS ####
name = "FULL_survey_site_country_in_rL"
random_factors = c("sample_unit", "site", "country")


hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data,
              X_data,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data,
                         X_data = X_data,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model, name,
                         run_python = T, save_path)




#### FULL MODEL SITES, QUADRATIC EFFECTS ####
quadratic_effects_all = colnames(
  dplyr::select(X_data_site,- site_code, -longitude,-latitude,-year,-country, 
                -ecoregion,-realm))

random_factors = c("sample_unit", "country")

name = "Quadratic_SITE_SCALE_site_country_in_rL"

hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, 
              quadratic_effects = quadratic_effects_all,
              random_factors, nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)



#### FULL MODEL SITES HIGH SHRINK ####
name = "test_1000shrink_var_SITE_SCALE_site_country_in_rL"
random_factors = c("sample_unit", "country")
set_shrink = 1000
quadratic_effects = NULL

hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = F, save_path)






#### NULL MODEL, SITES ####
name = "Null_model_SITE_SCALE_site_country_in_rL"
random_factors = c("sample_unit", "country")
set_shrink = NULL
test_null_model = TRUE

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink,
              test_null_model = test_null_model,
              name, run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data_site,
                         X_data = X_data_site,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model = test_null_model, 
                         name, run_python = T, save_path)






#### PCA DIMENSIONS FOR HABITAT COV, SITES ####
name = "Habitat_cov_in_PCA_dimensions_model_SITE_SCALE"
random_factors = c("sample_unit", "country")

X_data_site_pca <- covariates_site_PCA_hab

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site_pca,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)



#### WITHOUT PQ COVARIATES MODEL SITES ####
name = "Without_PQ_model_SITE_SCALE"
random_factors = c("sample_unit", "country")
X_data_pq <- X_data_site |> 
  dplyr::select(-Coral_RLS, -Sand_RLS, -Other_sessile_invert_RLS, -Rock_RLS, 
                -Coralline_algae_RLS, -Coral_rubble_RLS)

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_pq,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data_site,
                         X_data = X_data_pq,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model, name,
                         run_python = T, save_path)



#### WITHOUT REEFS WITH NA IN PQ DATA - MODEL SITES ####
name = "Without_NA_PQ_SITE_SCALE"
random_factors = c("sample_unit", "country")
X_data_pq <- covariates_site_without_NA_in_PQ
Y_data_pq <- observations_site_without_NA_in_PQ[rownames(covariates_site_without_NA_in_PQ),]
  

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_pq,
              X_data = X_data_pq,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)







#### WITHOUT AUSTRALIA - MODEL SITES ####
name = "Sensitivity_without_Australia_SITE_SCALE"
random_factors = c("sample_unit", "country")

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_no_aust,
              X_data = X_data_no_aust,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)



##--------------------------- Sensitivity analyses - REVIEW NEE -----------------------------

#### ONLY COUNTRY IN RANDOM ####
name = "Sensitivity_only_country_in_RL_SITE_SCALE"
random_factors = c("country")

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data_site,
                         X_data = X_data_site,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model, name,
                         run_python = T, save_path)



#### ONLY ECOREGION IN RANDOM ####
name = "Sensitivity_only_ecoregion_in_RL_SITE_SCALE"
random_factors = c("ecoregion")

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data_site,
                         X_data = X_data_site,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model, name,
                         run_python = T, save_path)


#### ECOREGION AND SAMPLE UNIT IN RANDOM ####
name = "Sensitivity_ecoregion&SU_in_RL_SITE_SCALE"
random_factors = c("sample_unit","ecoregion")

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data_site,
                         X_data = X_data_site,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model, name,
                         run_python = T, save_path)


#### NO RANDOM EFFECTS ####
name = "Sensitivity_NO_random_factors_SITE_SCALE"
random_factors = c()

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data_site,
                         X_data = X_data_site,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model, name,
                         run_python = T, save_path)

#### DETAILS MPA - FULL MODEL SITES ####
Y_data_site =  observations_site_final
X_data_site = covariates_site_final[rownames(Y_data_site),
        -which(colnames(covariates_site_final) == "protection_status")] |> 
  dplyr::rename(protection_status = protection_status_detailed)


name = "Sensitivity_MPA_features_detailed_10y_test2_SITE_SCALE"
random_factors = c("sample_unit", "country")

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)


#### MODEL SITES WITHOUT UNDERSAMPLED COUNTRIES ####
# name = "Sensitivity_countries_with_full_and_out_samples_model_SITE_SCALE"
# random_factors = c("sample_unit", "country")
# 
# selected_countries <- X_data_site |> 
#   dplyr::group_by(country, protection_status) |> 
#   dplyr::summarise(n_sample = dplyr::n()) |> 
#   dplyr::group_by(country) |> 
#   dplyr::mutate(n_conditions = length(unique(protection_status)),
#                 both_conditions = sum(c("out", "full") %in% unique(protection_status))) |> 
#   dplyr::filter(n_conditions > 1) |> 
#   dplyr::filter(both_conditions == 2) |> 
#   dplyr::pull(country)
# 
# # selected_countries <- X_data_site |> 
# #   dplyr::filter(country == "Australia") |> 
# #   dplyr::group_by(ecoregion, protection_status) |> 
# #   dplyr::summarise(n_sample = dplyr::n()) |> 
# #   dplyr::group_by(ecoregion) |> 
# #   dplyr::mutate(n_conditions = length(unique(protection_status)),
# #                 both_conditions = sum(c("out", "full") %in% unique(protection_status))) |> 
# #   dplyr::filter(n_conditions > 1) |> 
# #   dplyr::filter(both_conditions == 2)
# 
# X_data_site_subset <- X_data_site |> 
#   dplyr::filter(country %in% selected_countries)
# 
# Y_data_site_subset <- Y_data_site[rownames(X_data_site_subset),]
# 
# #Fit full model
# hmsc_function(nSamples, thin, nChains, verbose, transient,
#               Y_data = Y_data_site_subset,
#               X_data = X_data_site_subset,
#               response_distribution, quadratic_effects,random_factors,
#               nb_neighbours, set_shrink, test_null_model, name,
#               run_python = T, save_path)

#### MODEL SITES WITHOUT UNDERSAMPLED COUNTRIES AND OVERSAMPLED ECOREGION ####
name = "Sensitivity_without_undersampled_countries_and_oversampled_ecoregion_samples_model_SITE_SCALE"
random_factors = c("sample_unit", "country")

# Identify undersampled countries
selected_countries <- X_data_site |> 
  dplyr::group_by(country, protection_status) |> 
  dplyr::summarise(n_sample = dplyr::n()) |> 
  dplyr::group_by(country) |> 
  dplyr::mutate(n_conditions = length(unique(protection_status)),
                both_conditions = sum(c("out", "full") %in% unique(protection_status))) |> 
  dplyr::filter(n_conditions > 1) |> 
  dplyr::filter(both_conditions == 2) |> 
  dplyr::pull(country)

# Identify oversampled ecoregion
oversampled_ecoregion <- X_data_site |>
  dplyr::group_by(ecoregion, protection_status) |>
  dplyr::summarise(n_sample = dplyr::n()) |>
  dplyr::group_by(ecoregion) |>
  dplyr::mutate(n_tot_ecoregion = sum(n_sample),
                n_conditions = length(unique(protection_status)),
                both_conditions = sum(c("out", "full") %in% unique(protection_status))) |> 
  dplyr::ungroup() |> 
  dplyr::distinct(ecoregion, .keep_all = TRUE) |>
  dplyr::slice_max(order_by = n_tot_ecoregion, n = 3)|> 
  dplyr::pull(ecoregion)


X_data_site_subset <- X_data_site |> 
  dplyr::filter(country %in% selected_countries, #this remove 28 countries, and 444 samples
                !ecoregion %in% oversampled_ecoregion) #by removing the 3 most sampled ecoregion: remove 1194 samples

Y_data_site_subset <- Y_data_site[rownames(X_data_site_subset),]

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site_subset,
              X_data = X_data_site_subset,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)
