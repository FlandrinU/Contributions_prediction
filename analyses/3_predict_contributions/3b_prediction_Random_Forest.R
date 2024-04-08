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
ggsave(filename = here::here("figures", "models", "Pred_spRF_100tree_0.25_0.75.jpg"),
       width = 12, height = 7)


density_prediction(all_res)
ggsave(filename = here::here("figures", "models", "Pred_density_spRF_100tree_0.25_0.75.jpg"),
       width = 20, height = 10)





##------------- Multivariate RF -> impossible to run on a large matrix -------------

# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:2, FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- as.matrix(data[[1]][1:300, 1:2])
  
  # #Scale Y
  # Ymeans <- colMeans(data[[1]])
  # Ysd <- apply(data[[1]], 2, sd)
  # Y_train <- as.matrix(data[[1]][1:100, 1:5])
  # for(i in 1:ncol(Y_train)) Y_train[,i] <- (Y_train[,i] - Ymeans[i]) / Ysd[i]
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    dplyr::select(-longitude, -latitude, -country, - ecoregion) |> 
    dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
                                                "out" = 0,
                                                "Low" = 1,
                                                "Medium" = 2,
                                                "High" = 3)) |>  #effectiveness into quantitative values
    as.matrix()
  
  # XYcoords_train <- covariates_final[rownames(Y_train),] |>
  #   dplyr::select(longitude, latitude) |> 
  #   as.matrix()
  # 
  # #Moran's eigenvectors map predictors
  # SPeigen_train <- sjSDM::generateSpatialEV(XYcoords_train)[, 1:30] #reduce the dimension to always have 30 columns
  
  
  #Fit model and predict
  
  Y_test <- as.matrix(data[[2]][1:300, 1:2])
  # for(i in 1:ncol(Y_test)) Y_test[,i] <- (Y_test[,i] - Ymeans[i]) / Ysd[i]
  
  X_test <- covariates_final[rownames(Y_test),] |> 
    dplyr::select(-longitude, -latitude, -country, - ecoregion) |> 
    dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
                                                "out" = 0,
                                                "Low" = 1,
                                                "Medium" = 2,
                                                "High" = 3)) |>  #effectiveness into quantitative values
    as.matrix()
  
  preds <- MultivariateRandomForest::build_forest_predict(trainX = X_train,
                                                          trainY = Y_train,
                                                          n_tree = 100,
                                                          m_feature = round(ncol(X_train)/3,0),
                                                          min_leaf = 10,
                                                          testX = X_test)
  
  
  
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







##------------- Generalised Least Square RF (RF-GLS) -> problem of package updating?? -------------
install.packages(RandomForestsGLS)
library(RandomForestsGLS)


cross_val <- lapply(1:1,FUN = function(i){
  
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  
  XYcoords <- covariates_final |>
    dplyr::select(longitude, latitude) 
  
  Y_train <- data[[1]]
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    dplyr::select(-longitude, -latitude, -country, -ecoregion, -effectiveness) |> 
    as.matrix()
  
  coords_train <- XYcoords[rownames(Y_train),] |> as.matrix()

    
  #run model on each variable
  preds_raw <- lapply(colnames(Y_train),
                      # pbmcapply::pbmclapply(colnames(Y_train), mc.cores = parallel::detectCores()-3,
                      FUN = function(contrib){
                        
      # contrib = "herbivores_biomass"
      cat("contribution", which(contrib == colnames(Y_train)), "/", 
          length(colnames(Y_train)), ":", contrib, "\n")
      
      
      model <- RFGLS_estimate_spatial(coords = coords_train,
                                      y = Y_train[, contrib],
                                      X = X_train,
                                      Xtest = X_test,
                                      ntree = 50, cov.model = "exponential",
                                      nthsize = 20, param_estimate = TRUE)
      
      
      #New data
      Y_test <- data[[2]][,contrib]
      
      X_test <- covariates_final[rownames(data[[2]]), colnames(X_train)]
      
      coords_test <- XYcoords[rownames(X_test),]
        
      #Predict on new data
       preds <- RFGLS_predict_spatial(model, coords_test, X_test,
                                      h = 1, verbose = FALSE)
       preds

      # #Obs
      # plot(preds ~ Y_test[,contrib])
      # cor.test(preds, Y_test[,contrib])
      
      
    }) #END OF LAPPLY ON EACH VARIABLE
  
  #extract result
  preds <- do.call(cbind, preds_raw)
  colnames(preds) <- colnames(Y_train)
  rownames(preds) <- rownames(data[[2]])
  
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
  
  # ##obs
  # eval_act <- dplyr::filter(eval, variable == "herbivores_biomass")
  # plot(eval_act$imputed ~ eval_act$observed)
  # cor.test(eval_act$imputed, eval_act$observed)
  
  #Return dataframe
  eval
  
}) #END OF LAPPLY ON CROSSVALIDATION

# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models",
                             "Pred_spaMM_univariate_methodML_Matern_spatial_ONLY.jpg"),
       width = 12, height = 7)


#density plot
density_prediction(all_res)
ggsave(filename = here::here("figures", "models", 
                             "Pred_density_spaMM_univariate_methodML_Matern_spatial_ONLY.jpg"),
       width = 20, height = 10)





##------------- drf-------------
# install.packages('drf')
library(drf)

# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:2,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]]
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    dplyr::select(longitude, latitude)
  
  # dplyr::select(-longitude, -latitude, -country, - ecoregion)
  # |> 
  #   dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
  #                                               "out" = 0,
  #                                               "Low" = 1,
  #                                               "Medium" = 2,
  #                                               "High" = 3)) 
  
  
  #Fit model
  model <- drf::drf(X = X_train, Y = Y_train,
                    num.trees = 1000,
                    mtry = 15)
  
  # variableImportance(
  #   model,
  #   h = NULL,
  #   response.scaling = TRUE,
  #   type = "difference"
  # ) # /!\ long to run...

  #Predict on new data
  Y_test <- data[[2]][,colnames(Y_train)]
  
  X_test <- covariates_final[rownames(Y_test), colnames(X_train)]
  # |> 
  #   dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
  #                                               "out" = 0,
  #                                               "Low" = 1,
  #                                               "Medium" = 2,
  #                                               "High" = 3)) |> 
  #   as.matrix()
  
  predictions <- predict(model, newdata = X_test, functional = "mean")
  
  #Compare prediction
  preds <- predictions[["mean"]]
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
ggsave(filename = here::here("figures", "models", 
                             paste0("Pred_drf_multivariate_1000tree_mtry15.jpg")),
       width = 12, height = 7)


density_prediction(all_res)



##------------- CovRegRF -> useless to predict ?-------------
# install.packages('CovRegRF')
library(CovRegRF)

options(rf.cores=1, mc.cores=1)
## load generated example data
data(data, package = "CovRegRF")
xvar.names <- colnames(data$X)
yvar.names <- colnames(data$Y)
data1 <- data.frame(data$X, data$Y)
## define train/test split
set.seed(2345)
smp <- sample(1:nrow(data1), size = round(nrow(data1)*0.6), replace = FALSE)
traindata <- data1[smp,,drop=FALSE]
testdata <- data1[-smp, xvar.names, drop=FALSE]
## formula object
formula <- as.formula(paste(paste(yvar.names, collapse="+"), ".", sep=" ~ "))
## train covregrf
covregrf.obj <- covregrf(formula, traindata, params.rfsrc = list(ntree = 50),
                         importance = TRUE)





# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:2,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]]
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    dplyr::select(-longitude, -latitude, -country, - ecoregion)
  # |> 
  #   dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
  #                                               "out" = 0,
  #                                               "Low" = 1,
  #                                               "Medium" = 2,
  #                                               "High" = 3)) 
  
  train <- cbind(X_train, Y_train)
  
  formula <- as.formula(paste(paste(colnames(Y_train), collapse="+"), ".", sep=" ~ "))
  
  ## train covregrf
  covregrf.obj <- covregrf(formula, train, 
                           params.rfsrc = list(ntree = 50),
                           importance = F)
  
  
  #Predict on new data
  Y_test <- data[[2]][,colnames(Y_train)]
  
  X_test <- covariates_final[rownames(Y_test), colnames(X_train)]
  # |> 
  #   dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
  #                                               "out" = 0,
  #                                               "Low" = 1,
  #                                               "Medium" = 2,
  #                                               "High" = 3)) |> 
  #   as.matrix()
  
  test <- cbind(X_test, Y_test)
  
  predictions <- predict(covregrf.obj, newdata = test)
  
  #Compare prediction
  preds <- predictions$predicted
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
ggsave(filename = here::here("figures", "models", 
                             paste0("Pred_drf_multivariate_1000tree_mtry15.jpg")),
       width = 12, height = 7)


density_prediction(all_res)
