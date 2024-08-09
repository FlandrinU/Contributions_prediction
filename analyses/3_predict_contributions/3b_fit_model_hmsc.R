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

##-----------------------------Loading packages---------------------------------
pkgs <- c("here", "Hmsc", "coda", "ggmcmc")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

source("R/HMSC_function.R")

##------------------------------- load data ------------------------------------
load(here::here("data/derived_data/3_all_contributions_to_predict.Rdata"))
load(here::here("data/derived_data/3_all_covariates_to_predict.Rdata"))
load( here::here("data", "derived_data", "3_datasets_for_predict_CV_80_20.Rdata"))
Y_data =  observations_final#[sample(1:nrow(observations_final),1000),] ####################### reduce data
X_data = covariates_final[rownames(Y_data),] 
# |> 
#   dplyr::filter(country == "Australia")

Y_data =  observations_final[rownames(X_data),]
rownames(X_data) <- rownames(Y_data)

##----------------------------- Set-up parameters ------------------------------
nSamples = 50 #1000 
thin = 10 #100
nChains = 2 
verbose = 100 
nb_neighbours = 10

random_factors = c("associations","country")
response_distribution <- rep("normal", ncol(Y_data))
#response_distribution[colnames(Y_data) == "iucn_species_richness"] <- "poisson"

name = "test_spatial_full_and_country_in_rL"


##----------------------------- Run HMSC function ------------------------------

hmsc_function(nSamples,
              thin,
              nChains,
              verbose,
              transient = nSamples * thin,
              Y_data,
              X_data,
              response_distribution,
              random_factors,
              nb_neighbours,
              name,
              run_python = TRUE,
              save_path = here::here("outputs/models/hmsc"))


##----------------------------- Crossvalidation: predictive power ------------------------------
nSamples = 100 #1000 
thin = 2000 #100
nChains = 2 
verbose = 100 
nb_neighbours = 10
random_factors = c("country", "associations")
response_distribution <- rep("normal", ncol(Y_data))
name = "model_train"



predictions_cv <- parallel::mclapply(1:length(datasets), 
                                     mc.cores = 5,
                                     function(i){
  
  # Set parameters
  fold <- datasets[[i]]
  train <- fold[[1]]
  test <- fold[[2]]
  name_cv <- paste0(name, "_CV",i)
  
  Y_data_cv <- train
  X_data_cv <- covariates_final[rownames(Y_data_cv),]
  
  # Train model
  hmsc_function(nSamples,
                thin,
                nChains,
                verbose,
                transient = nSamples * thin,
                Y_data_cv,
                X_data_cv,
                response_distribution,
                random_factors,
                nb_neighbours,
                name_cv,
                run_python = T,
                save_path = here::here("outputs/models/hmsc/cross_validation"))
  
  # Import model output
  save_init <- paste0(save_path,"/init_multi/")
  save_out <- paste0(save_path, "out_multi/")
  localDir <- paste0(save_path, "multivariate/")
  file_name <- paste0(name_cv, "_", paste(nChains, "chains",
                                               thin, "thin",
                                               nSamples, "samples",
                                               sep = "_"),".rds")
  
  load( file.path(localDir, paste0("model_fit_", file_name, ".Rdata")))
  file_rds <- readRDS(file = paste0(save_out, "output_", file_name))[[1]]
  importFromHPC <- jsonify::from_json(file_rds)
  postList <- importFromHPC[1:nChains]
  model_fit_mcmc <- Hmsc::importPosteriorFromHPC(model_fit,
                                                 postList, 
                                                 nSamples, 
                                                 thin, 
                                                 transient = nSamples * thin)
  
  
  # Perform crossvalidation
  Y_data_test <- test
  X_data_test <- covariates_final[rownames(Y_data_test),]
  
  studyDesign <- data.frame(associations = as.factor(rownames(X_data_test)),
                            spatial = as.factor(rownames(X_data_test)),
                            year = as.factor(X_data_test$year),
                            country = as.factor(X_data_test$country),
                            ecoregion = as.factor(X_data_test$ecoregion))
  
  rL_asso = Hmsc::HmscRandomLevel(units = studyDesign$associations)
  rL_year = Hmsc::HmscRandomLevel(units = studyDesign$year)
  rL_ecoregion = Hmsc::HmscRandomLevel(units = studyDesign$ecoregion)
  rL_country =  Hmsc::HmscRandomLevel(units = studyDesign$country)
  
  ranLevels = list(associations = rL_asso,
                   year = rL_year,
                   ecoregion = rL_ecoregion,
                   country = rL_country)
  
  studyDesign <- studyDesign |> dplyr::select(all_of(random_factors))
  ranLevels <- ranLevels[random_factors]
  
  predY = predict(model_fit_mcmc, 
                  XData = X_data_test, 
                  studyDesign = studyDesign,
                  ranLevels = ranLevels, 
                  expected = T)
  
  PredY_mean <- as.data.frame(Reduce("+", predY) / length(predY))
  
  
  list(test, PredY_mean)
})

save(predictions_cv, file = here::here("outputs/models/hmsc/cross_validation/predictions_crossval_5folds.Rdata"))


##----------------------------- Sensitivity analyses: fit on several realms ------------------------------
nSamples = 300 #1000 
thin = 3000 #100
nChains = 3 
verbose = 100 
nb_neighbours = 10

random_factors = c("associations","country")
response_distribution <- rep("normal", ncol(Y_data))


# Run on the realms separately
realms <- unique(X_data$realm)

for(rlm in realms){
  cat("Fit model for realm ", rlm, "\n")
  
  X_data_realm <- dplyr::filter(X_data, realm == rlm)
  Y_data_realm <- Y_data[rownames(X_data_realm),]

  name = paste0("test_REALM_", gsub(" ", "_", rlm), "_asso_and_country_in_rL")
  
  hmsc_function(nSamples,
                thin,
                nChains,
                verbose,
                transient = nSamples * thin,
                Y_data_realm,
                X_data_realm,
                response_distribution,
                random_factors,
                nb_neighbours,
                name,
                run_python = TRUE,
                save_path = here::here("outputs/models/hmsc/sensitivity_analysis"))
  }




