################################################################################
##
##  
##
## 3b_prediction_gllvm.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
################################################################################
##----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("gllvm")
install.packages('gllvm')

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

##------------- gllvm multi-contributions-------------
library(gllvm)

# data(antTraits)
# y <- as.matrix(antTraits$abund)
# X <- scale(as.matrix(antTraits$env))
# TR <- antTraits$traits


# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:1,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- as.matrix(data[[1]])
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    dplyr::select(-country, -ecoregion, -longitude, -latitude) |> 
    dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
                                                "out" = 0,
                                                "Low" = 1,
                                                "Medium" = 2,
                                                "High" = 3)) |>  #effectiveness into quantitative values
    as.matrix()
  
  
  #Fit model
  model <- gllvm::gllvm(Y_train,X_train,
                        family = "gaussian")
  
  #Predict on new data
  Y_test <- as.matrix(data[[2]][,colnames(Y_train)])
  
  X_test <- covariates_final[rownames(Y_test), colnames(X_train)] |> 
    dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
                                                "out" = 0,
                                                "Low" = 1,
                                                "Medium" = 2,
                                                "High" = 3)) |> 
    as.matrix()
  
  preds <- predict(model, newX = X_test, level =0)
  
  
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
  
  # ##obs
  eval_act <- dplyr::filter(eval, variable == "herbivores_biomass")
  plot(eval_act$imputed ~ eval_act$observed)
  cor.test(eval_act$imputed, eval_act$observed)
  
  eval
}) #END OF LAPPLY ON CROSSVALIDATION


# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models", "Pred_gllvm_multivariate_without_long_lat.jpg"),
       width = 12, height = 7)


density_prediction(all_res)
ggsave(filename = here::here("figures", "models", "Pred_density_gllvm_multivariate.jpg"),
       width = 20, height = 10)



