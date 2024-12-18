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


## Only/Without Australia
X_data_aust = covariates_site_final |> dplyr::filter(country == "Australia")
Y_data_aust =  observations_site_final[rownames(X_data_aust),]

X_data_no_aust = covariates_site_final |> dplyr::filter(country != "Australia")
Y_data_no_aust =  observations_site_final[rownames(X_data_no_aust),]

##----------------------------- Set-up parameters ------------------------------
nSamples = 200
thin = 1000
nChains = 2 
verbose = 1000 
nb_neighbours = 10
transient = nSamples * thin

set_shrink = NULL
response_distribution <- rep("normal", ncol(Y_data))
#response_distribution[colnames(Y_data) == "iucn_species_richness"] <- "poisson"

quadratic_effects = NULL
test_null_model = NULL
# quadratic_effects = colnames(
#   dplyr::select(X_data,-longitude,-latitude,-year,-country, -ecoregion,-realm,
#                 -hdi,-marine_ecosystem_dependency,-ngo,-natural_ressource_rent))


save_path = here::here("outputs/models/hmsc")

##----------------------------- Fit HMSC models ------------------------------
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



# #### FULL MODEL SITES test mpa 2 ####
# name = "test_mpa2_FULL_model_SITE_SCALE"
# random_factors = c("sample_unit", "country")
# 
# table(X_data_site$protection_status2)
# X_data_site$protection_status <- NULL
# 
# #Fit full model
# hmsc_function(nSamples, thin, nChains, verbose, transient,
#               Y_data = Y_data_site,
#               X_data = X_data_site,
#               response_distribution, quadratic_effects,random_factors,
#               nb_neighbours, set_shrink, test_null_model, name,
#               run_python = T, save_path)
# 
# #Fit crossvalidation
# fit_hmsc_crossvalidation(k_fold = 5, 
#                          nSamples, thin, nChains, verbose, transient,
#                          Y_data = Y_data_site,
#                          X_data = X_data_site,
#                          response_distribution, quadratic_effects,random_factors,
#                          nb_neighbours, set_shrink, test_null_model, name,
#                          run_python = T, save_path)
# 

##--------------------------- Sensitivity analyses -----------------------------

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
              run_python = F, save_path)



#### FULL MODEL SITES, SPATIAL IN RRANDOM ####
name = "Spatial_effect_SITE_SCALE"
random_factors = c("sample_unit", "country", "spatial")

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_site,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
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
              name, run_python = F, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data_site,
                         X_data = X_data_site,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model = test_null_model, 
                         name, run_python = F, save_path)





#### UNIVARIATES MODEL SITES ####
# # Y <- Y_data_site |> dplyr::select(invertivores_biomass, public_interest )
# Y <- Y_data_site |> dplyr::select(actino_richness)
# response_distribution <- rep("normal", ncol(Y))
# name = "test_univariate_SITE_SCALE_site_country_in_rL"
# random_factors = c("sample_unit", "country")
# 
# #Fit model -> around 1h40 for each response
# hmsc_function(nSamples, thin, nChains, verbose, transient,
#               Y_data = Y,
#               X_data = X_data_site,
#               response_distribution, quadratic_effects,random_factors,
#               nb_neighbours, set_shrink, test_null_model, name,
#               run_python = T, save_path)
# 
# #Fit crossvalidation
# fit_hmsc_crossvalidation(k_fold = 5, 
#                          nSamples, thin, nChains, verbose, transient,
#                          Y_data = Y,
#                          X_data = X_data_site,
#                          response_distribution, quadratic_effects,random_factors,
#                          nb_neighbours, set_shrink, test_null_model, name,
#                          run_python = T, save_path)
# 
# 
# 
#### RANDOM ASSOCIATIONS BETWEEN CONTRIBUTIONS ####
name = "test_RANDOM_contrib_except_actino_SITE_SCALE"
random_factors = c("sample_unit", "country")

Y_random <- Y_data_site |> 
  dplyr::mutate(across(-actino_richness, .fns = sample, .names = "{.col}"))
cor.test(Y_random$actino_richness, Y_random$mean_endemism)
cor.test(Y_data_site$actino_richness, Y_data_site$mean_endemism)
cor.test(Y_random$actino_richness, Y_data_site$actino_richness)


#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_random,
              X_data = X_data_site,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_random,
                         X_data = X_data_site,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model, name,
                         run_python = T, save_path)






#### AUSTRALIA MODEL SITES ####
name = "Australia_only_SITE_SCALE"
random_factors = c("sample_unit")
X_data_aust <- X_data_aust |> dplyr::select(-hdi, -marine_ecosystem_dependency,
                                            -natural_ressource_rent)
#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_aust,
              X_data = X_data_aust,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)

#### WITHOUT AUSTRALIA MODEL SITES ####
name = "Without_Australia_SITE_SCALE"
random_factors = c("sample_unit", "country")
X_data_no_aust <- dplyr::select(X_data_no_aust, -Patch_Reefs_500m)

#Fit full model
hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data = Y_data_no_aust,
              X_data = X_data_no_aust,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)



#### WITHOUT ALLEN MODEL SURVEYS ####
name = "Without_Allen_survey_site_country_in_rL"
random_factors = c("sample_unit", "site", "country")


hmsc_function(nSamples, thin, nChains, verbose, transient,
              Y_data_wo_allen,
              X_data_wo_allen,
              response_distribution, quadratic_effects,random_factors,
              nb_neighbours, set_shrink, test_null_model, name,
              run_python = T, save_path)

#Fit crossvalidation
fit_hmsc_crossvalidation(k_fold = 5, 
                         nSamples, thin, nChains, verbose, transient,
                         Y_data = Y_data_wo_allen,
                         X_data = X_data_wo_allen,
                         response_distribution, quadratic_effects,random_factors,
                         nb_neighbours, set_shrink, test_null_model, name,
                         run_python = T, save_path)




# ##### fit on several realms ####
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

#### Test simple lm ####
data <- cbind(X_data_site, Y_data_site)

cov <- colnames(X_data_site)[8:46]
formula <- as.formula(paste("available_biomass ~", paste(cov, collapse = " + ")))
model <- lm(formula, data = data)

summary(model)

ggplot(data)+
  geom_violin(aes(x=effectiveness, y=available_biomass))+
  # geom_boxplot(aes(x=effectiveness, y=available_biomass))+
  geom_hline(yintercept = median(data[data$effectiveness == "out", "available_biomass"]),
             color = "red")

ggplot(data)+
  geom_boxplot(aes(x=effectiveness, y=n_fishing_vessels))+
  geom_violin(aes(x=effectiveness, y=n_fishing_vessels))+
  geom_hline(yintercept = median(data[data$effectiveness == "out", "n_fishing_vessels"]),
             color = "red")

ggplot(data)+
  geom_boxplot(aes(x=effectiveness, y=gravtot2))+
  geom_violin(aes(x=effectiveness, y=gravtot2))+
  geom_hline(yintercept = median(data[data$effectiveness == "out", "gravtot2"]),
             color = "red")

# ##--------------------------- Sensitivity analyses -----------------------------
#  source("R/evaluation_prediction_model.R")
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