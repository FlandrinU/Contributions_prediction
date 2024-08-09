################################################################################
##
##  
##
## 3b_prediction_mgpr.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
################################################################################
##----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("mgpr")
library(mgpr)

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

##------------- mgpr multi-contributions-------------

# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:2,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- as.matrix(data[[1]])[,1:20]
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    dplyr::select(-country, -ecoregion) |> 
    dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
                                                "out" = 0,
                                                "Low" = 1,
                                                "Medium" = 2,
                                                "High" = 3)) |>  #effectiveness into quantitative values
    as.matrix()
  
  
  #Fit model
  model <- mgpr::mgpr(datay = Y_train, datax = X_train,
                      kernel = "matern32",
                      # kernel ="rbf",
                      # kernpar = list(sigma = 0.37, corlen = 5.91, errorvar = 0.1),
                      meanf = "avg",
                      verbose=T
  )

  #Predict on new data
  Y_test <- as.matrix(data[[2]][,colnames(Y_train)])
  
  X_test <- covariates_final[rownames(Y_test), colnames(X_train)] |> 
    dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
                                                "out" = 0,
                                                "Low" = 1,
                                                "Medium" = 2,
                                                "High" = 3)) |> 
    as.matrix()
  
  preds <- predict(model, newdata = X_test,
                   fixneg = F)
  
  
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
  # eval_act <- dplyr::filter(eval, variable == "herbivores_biomass")
  # plot(eval_act$imputed ~ eval_act$observed)
  # cor.test(eval_act$imputed, eval_act$observed)
  
  eval
}) #END OF LAPPLY ON CROSSVALIDATION


# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models", "Pred_mgpr_multivariate_mattern32_20var.jpg"),
       width = 12, height = 7)


density_prediction(all_res)
ggsave(filename = here::here("figures", "models", "Pred_density_mgpr_multivariate_mattern32_20var.jpg"),
       width = 20, height = 10)



