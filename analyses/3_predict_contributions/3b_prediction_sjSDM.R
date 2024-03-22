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
# Sys.unsetenv("GITHUB_PAT")
# devtools::install_github("https://github.com/TheoreticalEcology/s-jSDM", subdir = "sjSDM", ref = "master")
# 
# install.packages("sjSDM")
# #sjSDM::install_sjSDM(method = "gpu")
# sjSDM::install_sjSDM(method = "cpu")

cross_val <- lapply(1:length(datasets),FUN = function(i){
# cross_val <- lapply(1:4,FUN = function(i){

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
                                        hidden = c(50L, 50L),
                                        activation = "relu"),
                       spatial = sjSDM::DNN(SPeigen_train, 
                                            hidden = c(5L, 5L),
                                            ~ 0+.),
                       # env = sjSDM::linear(X_train, ~.),
                       # spatial = sjSDM::linear(SPeigen_train, ~0+.),
                       
                       control = sjSDM::sjSDMControl(optimizer = sjSDM::RMSprop(),
                                                     early_stopping_training = 15,
                                                     scheduler = 5,
                                                     lr_reduce_factor = 0.5),
                       learning_rate = 0.01,
                       iter = 500L, 
                       
                       se = F,
                       family = gaussian(),
                       device = "cpu")
  
  
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
  
  eval
}) #END OF LAPPLY ON CROSSVALIDATION


# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models", "Pred_sjSMD_DNN_50hidden.jpg"),
       width = 12, height = 7)



# DENSITY PLOT
res_scaled <- all_res |> 
  dplyr::group_by(variable) |> 
  dplyr::mutate(min_obs = min(observed), max_obs = max(observed),
                min_imp = min(imputed), max_imp = max(imputed)) |> 
  dplyr::mutate(observed = (observed - min_obs)/(max_obs - min_obs),
                imputed = (imputed - min_imp)/(max_imp- min_imp))

ggplot(res_scaled, aes(x = observed, y = imputed)) +
  geom_density_2d_filled(aes(x = observed, y = imputed),
                         contour_var = 'ndensity',
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
ggsave(filename = here::here("figures", "models", "Pred_density_sjSMD_DNN_50hidden.jpg"),
       width = 20, height = 10)


##------------- spaMM -------------
# install.packages("spaMM")

#IMPOSSIBLE TO LOAD PKG AFTER VERSION 4.0.0
# install.packages("data/raw_data/spaMM_3.13.0.tar.gz", repos = NULL, type = "source")

library(spaMM)

cross_val <- lapply(1:length(datasets),FUN = function(i){
  # cross_val <- lapply(1:4,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]]
  
  X_train <- covariates_final[rownames(Y_train),] 
  
  train <- cbind(X_train, Y_train) |> 
    tibble::rownames_to_column("survey_id")
  
  formulas = list(
    
      # formula=cbind(actino_richness, elasmobranch_richness, functional_distinctiveness,
      #               iucn_species_richness, mean_endemism, evolutionary_distinctiveness,
      #               functional_entropy, phylogenetic_entropy, herbivores_biomass,
      #               invertivores_biomass, piscivores_biomass, trophic_web_robustness,
      #               mean_trophic_level, available_biomass, selenium,
      #               iron, vitamin_A, available_biomass_turnover,
      #               NN_score, NP_score) ~ 1 +
      #               depth + year + median_7days_chl +median_5year_chl+
      #               median_7days_degree_heating_week+median_5year_degree_heating_week+
      #               (1|survey_id)+
      #               Matern(1|longitude+latitude),
      # family=gaussian())

    actino_richness=list(
      formula = actino_richness ~ 1 + year + Matern(1|longitude+latitude),
      family=gaussian()),

    functional_distinctiveness=list(
      formula = functional_distinctiveness ~ 1 + year + Matern(1|longitude+latitude),
      # + corrMatrix(0+mv(1,2)|survey_id) +
      #   (0+mv(1,2)|survey_id)+ Matern(1|longitude+latitude),
      family=gaussian())
    )

  #Fit model
  model <- spaMM::fitmv(submodels = formulas,
                        data = train,
                        method = "ML")
  
  
  #Predict on new data
  Y_test <- data[[2]]
  
  X_test <- covariates_final[rownames(Y_test),]
  
  test <- cbind(X_test, Y_test) |> 
    tibble::rownames_to_column("survey_id")
  
  preds <- predict(model, newdata = test)
  
  #Compare prediction
  
  # #make a table
  # cuts <- c(0, seq(from = nrow(Y_test), to = nrow(preds), by = nrow(Y_test)))
  # preds_table <- sapply(2:length(cuts), 
  #                       function(i) preds[seq(cuts[i-1]+1, cuts[i]),1])
  # 
  # colnames(preds_table) <- names(formulas)
  # rownames(preds_table) <- rownames(Y_test)
  # 
  # preds_long <- as.data.frame(preds) |> 
  #   tibble::rownames_to_column("survey_id") |> 
  #   tidyr::pivot_longer(cols = -survey_id ,
  #                       names_to = "variable",
  #                       values_to = "imputed" )
  
  preds_long <- as.data.frame(preds) |> 
    dplyr::mutate(survey_id = rep(rownames(Y_test), times = nrow(preds)/nrow(Y_test))) |> 
    dplyr::mutate(variable = rep(names(formulas), each = nrow(Y_test))) |> 
    dplyr::rename(imputed = V1)
    
  eval <- as.data.frame(Y_test) |> 
    tibble::rownames_to_column("survey_id") |> 
    tidyr::pivot_longer(cols = names(formulas) ,
                        names_to = "variable",
                        values_to = "observed" ) |> 
    dplyr::full_join(preds_long) |> 
    dplyr::mutate(model = i)
  
  eval
}) #END OF LAPPLY ON CROSSVALIDATION



##------------- cito -------------
## installation
remove.packages("torch")
remove.packages("cito")

# check package 
if(!require('torch',quietly = TRUE)) install.packages('torch')
library('torch') 

#install torch
torch::install_torch()
torch::install_torch(torch_lantern = TRUE)

Sys.unsetenv("GITHUB_PAT")
remotes::install_github("mlverse/torch")


#install cito
install.packages("cito")

if(!require('devtools', quietly = TRUE)) install.packages('devtools')
devtools::install_github('citoverse/cito')
Sys.unsetenv("GITHUB_PAT")




cross_val <- lapply(1:length(datasets),FUN = function(i){
  # cross_val <- lapply(1:4,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- as.matrix(data[[1]])
  
  X_train <- covariates_final[rownames(Y_train),] |> 
    # dplyr::select(-longitude, -latitude) |> 
    as.matrix()
  

  #Fit model
  model <- cito::dnn(X = Xtrain, Y = Y_train)
  
  
  #Predict on new data
  Y_test <- as.matrix(data[[2]])
  
  X_test <- covariates_final[rownames(Y_test),]|> 
    # dplyr::select(-longitude, -latitude) |> 
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





##------------- Multivariate RF -> impossible to run on a large matrix -------------

# cross_val <- lapply(1:length(datasets),FUN = function(i){
cross_val <- lapply(1:4, FUN = function(i){
  
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  # i=1
  data = datasets[[i]] 
  
  Y_train <- as.matrix(data[[1]][1:100, 1:5])
  
  # #Scale Y
  # Ymeans <- colMeans(data[[1]])
  # Ysd <- apply(data[[1]], 2, sd)
  # Y_train <- as.matrix(data[[1]][1:100, 1:5])
  # for(i in 1:ncol(Y_train)) Y_train[,i] <- (Y_train[,i] - Ymeans[i]) / Ysd[i]
  
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
  # for(i in 1:ncol(Y_test)) Y_test[,i] <- (Y_test[,i] - Ymeans[i]) / Ysd[i]
  
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




