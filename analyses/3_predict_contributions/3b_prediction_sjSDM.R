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
#full data
load(file = here::here("data", "derived_data", "3_all_contributions_to_predict.Rdata"))

#datasets to predict
load( file = here::here("data", "derived_data", "3_datasets_for_predict_CV_80_20.Rdata"))

#covariates
load(file = here::here("data", "derived_data", "3_all_covariates_to_predict.Rdata"))

#load functions
source("R/evaluation_prediction_model.R")


##------------- sjSDM with DNN-------------


cross_val <- lapply(1:length(datasets),FUN = function(i){
# cross_val <- lapply(1:4,FUN = function(i){
  
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  # i=1
  data = datasets[[i]] 
  
  # Y_train <- as.matrix(data[[1]])
  
  #Scale Y
  Ymeans <- colMeans(data[[1]])
  Ysd <- apply(data[[1]], 2, sd)
  Y_train <- as.matrix(data[[1]])
  for(i in 1:ncol(data[[1]])) Y_train[,i] <- (Y_train[,i] - Ymeans[i]) / Ysd[i]
  
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    dplyr::select(-longitude, -latitude) |> 
    as.matrix()
  
  XYcoords_train <- covariates_final[rownames(Y_train),] |>
    dplyr::select(longitude, latitude) |> 
    as.matrix()

  #Moran's eigenvectors map predictors
  SPeigen_train <- sjSDM::generateSpatialEV(XYcoords_train)[, 1:30] #reduce the dimension to always have 30 columns
  
  #Fit model
  model <- sjSDM::sjSDM(Y_train, 
                       env = sjSDM::DNN(X_train,
                                        hidden = c(50L,50L),
                                        activation = "relu"),
                       spatial = sjSDM::DNN(SPeigen_train, ######################################  coordinates or eigenvectors ?
                                            hidden = c(5L, 5L),
                                            ~ 0+.),
                       # env = sjSDM::linear(X_train, ~.),
                       # spatial = sjSDM::linear(SPeigen_train, ~0+.),
                       
                       control = sjSDM::sjSDMControl(optimizer = sjSDM::RMSprop(),
                                                     early_stopping_training = 10,
                                                     scheduler = 10,
                                                     lr_reduce_factor = 0.5),
                       learning_rate = 0.01,
                       iter = 500L, 
                       
                       se = F,
                       family = gaussian())
  
  
  #Predict on new data
  Y_test <- as.matrix(data[[2]])
  
  #scale Y
  for(i in 1:ncol(data[[2]])) Y_test[,i] <- (Y_test[,i] - Ymeans[i]) / Ysd[i]
  
  X_test <- covariates_final[rownames(Y_test),]|> 
    dplyr::select(-longitude, -latitude) |> 
    as.matrix()
  XYcoords_test <- covariates_final[rownames(Y_test),] |>
    dplyr::select(longitude, latitude) |> 
    as.matrix()
  SPeigen_test <- sjSDM::generateSpatialEV(XYcoords_test)[, 1:30]
  
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
  
  eval
}) #END OF LAPPLY ON CROSSVALIDATION


# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)



# create basic plot
ggplot(all_res) +
  geom_density_2d_filled(aes(x = observed, y = imputed),
                         contour_var = 'count',
                         contour = F, n = 100, bins= 10, colour = 'transparent') +
  scale_fill_viridis_d(option = 'viridis', begin = 0.2, end = 0.9,
                       name = 'Count') +
  # geom_point(aes(x = observed, y = predicted), alpha = 0.2)
  theme_bw()+ 
  theme(panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'grey90', colour = 'grey90'), 
        aspect.ratio = 1) + 
  facet_wrap(~variable, scales = "free") +
  geom_abline() +
  labs(x = "Observed", y = "Predicted") +
  # geom_text(label = paste0("r² = ", round(r2,3), "\nestimate = ", round(a,2)),
  #           x= 0.3, y=0.9, size = 5)+
  theme(
    axis.text=element_text(size=10),
    axis.title=element_text(size=15),
    legend.text=element_text(size=8), 
    legend.title=element_text(size=8),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())




##------------- Multivariate RF -> impossible to run on a large matrix -------------

# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:4,FUN = function(i){
  
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  # i=1
  data = datasets[[i]] 
  
  #Scale Y
  Ymeans <- colMeans(data[[1]])
  Ysd <- apply(data[[1]], 2, sd)
  Y_train <- as.matrix(data[[1]][1:100, 1:5])
  for(i in 1:ncol(Y_train)) Y_train[,i] <- (Y_train[,i] - Ymeans[i]) / Ysd[i]
  
  X_train <- covariates_final[rownames(Y_train),1:10] |> 
    dplyr::select(-longitude, -latitude) |> 
    as.matrix()
  
  # XYcoords_train <- covariates_final[rownames(Y_train),] |>
  #   dplyr::select(longitude, latitude) |> 
  #   as.matrix()
  # 
  # #Moran's eigenvectors map predictors
  # SPeigen_train <- sjSDM::generateSpatialEV(XYcoords_train)[, 1:30] #reduce the dimension to always have 30 columns
  
  
  #Fit model and predict
  
  Y_test <- as.matrix(data[[2]][1:100, 1:5])
  for(i in 1:ncol(Y_test)) Y_test[,i] <- (Y_test[,i] - Ymeans[i]) / Ysd[i]
  
  X_test <- covariates_final[rownames(Y_test),1:10]|> 
    dplyr::select(-longitude, -latitude) |> 
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


##------------- ..? -------------
