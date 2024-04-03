################################################################################
##
##  
##
## 3b_prediction_sjSDM.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
################################################################################
##----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("sjSDM")


##-------------loading data and functions-------------
# #full data
# load(file = here::here("data", "derived_data", "3_all_contributions_to_predict.Rdata"))
# 
# #datasets to predict
# load( file = here::here("data", "derived_data", "3_datasets_for_predict_CV_80_20.Rdata"))
# 
# #covariates
# load(file = here::here("data", "derived_data", "3_all_covariates_to_predict.Rdata"))

## SITE SCALE
load(file = here::here("data", "derived_data", "3_sites_contributions_to_predict.Rdata"))
load( file = here::here("data", "derived_data", "3_sites_datasets_for_predict_CV_80_20.Rdata"))
load(file = here::here("data", "derived_data", "3_sites_covariates_to_predict.Rdata"))

#load functions
source("R/evaluation_prediction_model.R")


##------------- sjSDM with DNN-------------
# Sys.unsetenv("GITHUB_PAT")
# devtools::install_github("https://github.com/TheoreticalEcology/s-jSDM", subdir = "sjSDM", ref = "master")
# 
# install.packages("sjSDM")
# #sjSDM::install_sjSDM(method = "gpu")
# sjSDM::install_sjSDM(method = "cpu")

# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:2,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  #Moran's eigenvectors map predictors
  XYcoords <- covariates_final |>
    dplyr::select(longitude, latitude) |> 
    as.matrix()
  SPeigen <- sjSDM::generateSpatialEV(XYcoords)[, 1:30] #reduce the dimension to always have 30 columns
  rownames(SPeigen) <- rownames(covariates_final)
  
  # TRAIN
  Y_train <- as.matrix(data[[1]])
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    dplyr::select(-longitude, -latitude) |> 
    as.matrix()
  
  SPeigen_train <- SPeigen[rownames(Y_train),]
  
  
  #Fit model
  model <- sjSDM::sjSDM(Y_train, 
                        env = sjSDM::DNN(X_train,
                                         hidden = c(50L, 50L, 50L),
                                         activation = "selu"),
                        spatial = sjSDM::DNN(SPeigen_train, 
                                             hidden = c(5L, 5L),
                                             ~ 0+.),
                        # env = sjSDM::linear(X_train, ~.),
                        # spatial = sjSDM::linear(SPeigen_train, ~0+.),
                        
                        control = sjSDM::sjSDMControl(optimizer = sjSDM::RMSprop(),
                                                      early_stopping_training = 10,
                                                      scheduler = 5,
                                                      lr_reduce_factor = 0.5),
                        learning_rate = 0.01, #good with 0.01
                        step_size = round(0.1*nrow(X_train)),
                        # sampling = 100L,
                        iter = 200L, 
                        
                        se = F,
                        family = gaussian(),
                        device = "cpu")
  # plot(model)
  
  
  #Predict on new data
  Y_test <- as.matrix(data[[2]])
  
  X_test <- covariates_final[rownames(Y_test),] |> 
    dplyr::select(-longitude, -latitude) |> 
    as.matrix()
  
  SPeigen_test <- SPeigen[rownames(Y_test),]
  
  preds <- predict(model, newdata = X_test, SP = SPeigen_test)
  
  #Compare prediction
  colnames(preds) <- colnames(Y_test)
  rownames(preds) <- rownames(Y_test)
  
  preds_long <- as.data.frame(preds) |> 
    tibble::rownames_to_column("survey_id") |> 
    tidyr::pivot_longer(cols = -survey_id ,
                        names_to = "variable",
                        values_to = "imputed" )
  
  eval <- as.data.frame(Y_test) |> 
    tibble::rownames_to_column("survey_id") |> 
    tidyr::pivot_longer(cols = -survey_id ,
                        names_to = "variable",
                        values_to = "observed" ) |> 
    dplyr::full_join(preds_long) |> 
    dplyr::mutate(model = i)
  
  ##obs
  eval_act <- dplyr::filter(eval, variable == "herbivores_biomass")
  plot(eval_act$imputed ~ eval_act$observed)
  cor.test(eval_act$imputed, eval_act$observed)
  
  eval
}) #END OF LAPPLY ON CROSSVALIDATION


# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models", "Pred_sjSMD_DNN_3x50hidden_lr001.jpg"),
       width = 12, height = 7)

#density plot
density_prediction(all_res)
ggsave(filename = here::here("figures", "models", "Pred_density_sjSMD_DNN_3x50hidden_lr001.jpg"),
       width = 20, height = 10)


