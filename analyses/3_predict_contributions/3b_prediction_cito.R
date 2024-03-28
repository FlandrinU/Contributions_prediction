################################################################################
##
##  
##
## 3b_prediction_cito.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
################################################################################
##----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("cito")


##-------------loading data and functions-------------
#full data
load(file = here::here("data", "derived_data", "3_all_contributions_to_predict.Rdata"))

#datasets to predict
load( file = here::here("data", "derived_data", "3_datasets_for_predict_CV_80_20.Rdata"))

#covariates
load(file = here::here("data", "derived_data", "3_all_covariates_to_predict.Rdata"))

#load functions
source("R/evaluation_prediction_model.R")

##------------- cito multi-contributions-------------
# ## installation
# remove.packages("torch")
# remove.packages("cito")
# 
# # check package 
# if(!require('torch',quietly = TRUE)) install.packages('torch')
# # library('torch') 
# 
# #install torch
# torch::install_torch()
# 
# Sys.unsetenv("GITHUB_PAT")
# remotes::install_github("mlverse/torch")
# 
# #install cito
# install.packages("cito")
# 
# if(!require('devtools', quietly = TRUE)) install.packages('devtools')
# devtools::install_github('citoverse/cito')
# Sys.unsetenv("GITHUB_PAT")


###  IF PROBLEM: ###
# -> run cito::dnn(...)
# -> Restart R
# -> run torch::install_torch()
# => cito should work, but is in conflict with sjSDM
###


# cross_val <- lapply(1:length(datasets),FUN = function(i){
   cross_val <- lapply(1:2,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- as.matrix(data[[1]])
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    # dplyr::select(-longitude, -latitude) |> 
    as.matrix()
  
  
  #Fit model
  model <- cito::dnn(X = X_train, Y = Y_train,
                     hidden = c(50L,50L,50L),
                     #hidden = round(mean(ncol(X_train), ncol(Y_train))),
                     epochs =100,
                     lr = 0.01, #learning rate
                     # loss = "gaussian",
                     activation = "selu",
                     lr_scheduler = cito::config_lr_scheduler("reduce_on_plateau", 
                                                              patience = 5,
                                                              factor = 0.5), #divide by 2 lr if the loss is constant over 5 epochs
                     validation = 0.2, #take 20% of the data to test the model
                     early_stopping = 10,
                     lambda = 0.001,
                     alpha = 0.2
                     # ,
                     # batchsize = round(nrow(X_train) * 0.1)
                     # ,
                     # bootstrap = 10
                     
  )
  
  # #Evaluate training
  # plot(model)
  # cito::analyze_training(model)
  # summary(model)
  
  #Predict on new data
  Y_test <- as.matrix(data[[2]][,colnames(Y_train)])
  
  X_test <- covariates_final[rownames(Y_test), colnames(X_train)]|> 
    as.matrix()
  
  preds <- predict(model, newdata = X_test)
  
  
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
  
  eval
}) #END OF LAPPLY ON CROSSVALIDATION


# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models", "Pred_cito_multivariate_3x50hiddens.jpg"),
       width = 12, height = 7)


density_prediction(all_res)
ggsave(filename = here::here("figures", "models", ".jpg"),
       width = 20, height = 10)




##------------- cito NN and NP-------------

# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:5,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- as.matrix(data[[1]][,c("NN_score", "NP_score")])
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    dplyr::select(-longitude, -latitude) |> 
    as.matrix()
  
  
  #Fit model
  model <- cito::dnn(X = X_train, Y = Y_train,
                     hidden = c(25L),
                     epochs = 100,
                     lr = 0.0001,
                     lr_scheduler = cito::config_lr_scheduler("reduce_on_plateau", 
                                                              patience = 5,
                                                              factor = 0.5),
                     validation = 0.1,
                     early_stopping = 10, lambda = 0.001, alpha = 0.2
                     
  )
  
  # #Evaluate training
  # plot(model)
  # cito::analyze_training(model)
  # summary(model)
  
  #Predict on new data
  Y_test <- as.matrix(data[[2]][,colnames(Y_train)])
  
  X_test <- covariates_final[rownames(Y_test), colnames(X_train)]|> 
    as.matrix()
  
  preds <- predict(model, newdata = X_test)
  
  
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
  
  eval
}) #END OF LAPPLY ON CROSSVALIDATION


# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models", "Pred_cito_NN_NP.jpg"),
       width = 12, height = 7)


density_prediction(all_res)
ggsave(filename = here::here("figures", "models", ".jpg"),
       width = 20, height = 10)



##------------- cito univariate-------------

# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:2,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y <- as.matrix(data[[1]])
  
  X_train <- covariates_final[rownames(Y),] |> 
    dplyr::select(-longitude, -latitude) |> 
    as.matrix()
  
  preds_raw <- lapply(colnames(Y), FUN = function(contrib){
    
    Y_train <- Y[,contrib]
    
    #Fit model
    model <- cito::dnn(X = X_train, Y = Y_train,
                       hidden = c(25L),
                       epochs = 100,
                       lr = 0.01,
                       lr_scheduler = cito::config_lr_scheduler("reduce_on_plateau", 
                                                                patience = 10,
                                                                factor = 0.5),
                       validation = 0.2,
                       early_stopping = 20, lambda = 0.001, alpha = 0.2
                       
                        )
    
    #Predict on new data
    Y_test <- as.matrix(data[[2]])[,contrib]
    
    X_test <- covariates_final[names(Y_test), colnames(X_train)]|> 
      as.matrix()
    
     preds_contrib <- predict(model, newdata = X_test)
     # plot(preds_contrib~Y_test)
     preds_contrib
  })
  
  #extract result
  preds <- do.call(cbind, preds_raw)
  colnames(preds) <- colnames(Y)
  rownames(preds) <- rownames(data[[2]])
  
  #Compare prediction
  preds_long <- as.data.frame(preds) |> 
    tibble::rownames_to_column("survey_id") |> 
    tidyr::pivot_longer(cols = -survey_id ,
                        names_to = "variable",
                        values_to = "imputed" )
  
  eval <- data[[2]] |> 
    tibble::rownames_to_column("survey_id") |> 
    tidyr::pivot_longer(cols = -survey_id ,
                        names_to = "variable",
                        values_to = "observed" ) |> 
    dplyr::full_join(preds_long) |> 
    dplyr::mutate(model = i)
  
  eval
  # ##obs
  # eval_act <- dplyr::filter(eval, variable == "mean_endemism")
  # eval_act <- dplyr::filter(eval, variable == "actino_richness")
  # plot(eval_act$imputed ~ eval_act$observed)
  # cor.test(eval_act$imputed, eval_act$observed)
  
}) #END OF LAPPLY ON CROSSVALIDATION


# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models", "Pred_cito_univariates_25hidden.jpg"),
       width = 12, height = 7)


density_prediction(all_res)
ggsave(filename = here::here("figures", "models", ".jpg"),
       width = 20, height = 10)

