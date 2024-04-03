################################################################################
##
##  
##
## 3b_prediction_Random_Forest.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
################################################################################
##----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("SpatialML", "randomforest)


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
source("R/grf_function_package_SpatialML.R")


##------------- RF univariate-------------

# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:4,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y <- data[[1]]
  
  X_train <- covariates_final[rownames(Y),] |> 
    dplyr::select(-longitude, -latitude) 
  
  # PREDICT EACH CONTRIBUTION INDIVIDUALLY
  preds_raw <- pbmcapply::pbmclapply(colnames(Y),
                                     # mc.cores = parallel::detectCores()-10,
                                     mc.cores = 15,
                                     FUN = function(contrib){
                                       # contrib <- "actino_richness"
                                       Y_train <- Y[rownames(X_train),contrib]
                                       
                                       # train <- cbind(X_train, Y_train)
                                       # fmla <- as.formula(paste("Y_train", " ~ ", paste(colnames(X_train), collapse= "+")))
                                       
                                       Y_test <- as.matrix(data[[2]])[,contrib]
                                       X_test <- covariates_final[names(Y_test), colnames(X_train)]
                                       
                                       #Fit model
                                       model <- randomForest::randomForest(x = X_train,
                                                                           y = Y_train,
                                                                           # formula = fmla,
                                                                           # data = train,
                                                                           ntree = 100,
                                                                           xtest = X_test,
                                                                           ytest = Y_test
                                       )
                                       
                                       
                                       # #Predict on new data
                                       # preds_contrib <- predict(model, new.data = X_test, type = 'response')
                                       # preds_contrib
                                       model[["test"]][["predicted"]]
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
ggsave(filename = here::here("figures", "models", "Pred_RF_100tree.jpg"),
       width = 12, height = 7)


density_prediction(all_res)
ggsave(filename = here::here("figures", "models", "Pred_density_RF_100tree.jpg"),
       width = 20, height = 10)



##------------- spRF univariate-------------

# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:2,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y <- data[[1]]
  
  X_train <- covariates_final[rownames(Y),] |> 
    dplyr::select(-longitude, -latitude) 
  
  # PREDICT EACH CONTRIBUTION INDIVIDUALLY
  preds_raw <- pbmcapply::pbmclapply(colnames(Y),
                                     mc.cores = parallel::detectCores()-10,
                                     FUN = function(contrib){
    # contrib <- "NN_score"
    Y_train <- Y[rownames(X_train),contrib]
    
    train <<- cbind(X_train, Y_train)
    fmla <<- as.formula(paste("Y_train", " ~ ", paste(colnames(X_train), collapse= "+")))
    coords <<- covariates_final[rownames(train),] |> 
      dplyr::select(longitude, latitude) 
    
    #Fit model
    model <- GRF(formula = fmla,
                            dframe = train,
                            bw = 10, #15
                            kernel = "adaptive",
                            coords = coords,
                            ntree = 100, #1000 or 100
                            geo.weighted = FALSE)
    
    
    #Predict on new data
    Y_test <- as.matrix(data[[2]])[,contrib]
    
    X_test <- covariates_final[names(Y_test), c("longitude", "latitude",colnames(X_train))] 
    
    preds_contrib <- SpatialML::predict.grf(model, new.data = X_test, 
                                            x.var.name = "longitude" ,
                                            y.var.name = "latitude",
                                            local.w=0.25, global.w=.75)
    # plot(preds_contrib~Y_test)
    # summary(lm(preds_contrib~Y_test))
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
ggsave(filename = here::here("figures", "models", "Pred_spRF_10tree_0.5_0.5.jpg"),
       width = 12, height = 7)


density_prediction(all_res)
ggsave(filename = here::here("figures", "models", "Pred_density_spRF_10tree_0.5_0.5.jpg"),
       width = 20, height = 10)



