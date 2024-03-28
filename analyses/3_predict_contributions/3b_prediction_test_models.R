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
                       learning_rate = 0.1,
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

#density plot
density_prediction(all_res)
ggsave(filename = here::here("figures", "models", ".jpg"),
       width = 20, height = 10)



##------------- spaMM -------------
# install.packages("spaMM")

#IF IMPOSSIBLE TO LOAD PKG AFTER VERSION 4.0.0
# install.packages("data/raw_data/spaMM_3.13.0.tar.gz", repos = NULL, type = "source")

# library(spaMM)

# cross_val <- pbmcapply::pbmclapply(1:length(datasets), mc.cores = 5,
#                                    FUN = function(i){
# cross_val <- lapply(1:2,FUN = function(i){
cross_val <- pbmcapply::pbmclapply(1:5, mc.cores = 5,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]][1:300,] ################################" reduce data
  
  X_train <- covariates_final[rownames(Y_train),]
  
  train <- cbind(X_train, Y_train) |> 
    tibble::rownames_to_column("survey_id")
  
  covariates <- colnames(dplyr::select(X_train, -longitude, -latitude))
  
  
  # List of submodels
  
    # submodels = list(
      # list(formula=cbind(actino_richness, elasmobranch_richness, functional_distinctiveness,
      #               iucn_species_richness, mean_endemism, evolutionary_distinctiveness,
      #               functional_entropy, phylogenetic_entropy, herbivores_biomass,
      #               invertivores_biomass, piscivores_biomass, trophic_web_robustness,
      #               mean_trophic_level, available_biomass, selenium,
      #               iron, vitamin_A, available_biomass_turnover,
      #               NN_score, NP_score) ~ 1 +
      #               depth + year + median_7days_chl +median_5year_chl+
      #               median_7days_degree_heating_week+median_5year_degree_heating_week+
      #               # (1|survey_id)+
      #               Matern(1|longitude+latitude),
      # family=gaussian())
    # )
  
  responses <- colnames(Y_train)[5:7]
  # responses <- c("actino_richness", "mean_endemism", "available_biomass_turnover",
  #                "NN_score", "iron", "available_biomass", "functional_entropy")##############################"
  submodels <- list()
  
  for (response in responses) {
    form <- as.formula(
      paste(response, "~ ",
            paste(c("1", covariates, "Matern(1|longitude+latitude)"
                    # , paste0("spaMM::corrMatrix(0+spaMM::mv(",
                    #          paste(seq(1,length(responses)),collapse = ","),
                    #          ")|survey_id)")
                    # , paste0("(0+spaMM::mv(",
                    #          paste(seq(1,length(responses)),collapse = ","),
                    #          ")|survey_id)")
            ) , 
            collapse = "+")))
    submodels[[response]] <- list(formula = form, family = gaussian())
  }
  
  #New data
  Y_test <- data[[2]][,responses]
  
  X_test <- covariates_final[rownames(Y_test),]
  
  test <- cbind(X_test, Y_test) |> 
    tibble::rownames_to_column("survey_id")
  
  
tryCatch({
    #Fit the model
    model <- spaMM::fitmv(submodels,
                          data = train,
                          method = "ML",
                          control.HLfit = list(
                            NbThreads = parallel::detectCores()-5))
    
  
    #Predict on new data
    preds <- predict(model, newdata = test)

    }, error = function(e){
      
      preds <<- data.frame(V1 = rep(NA, nrow(Y_test) * ncol(Y_test)))

  })#End of tryCatch
    
  
  #Compare prediction
  var_to_predict <-  names(
    unlist( 
      lapply(submodels, function(mod) {
        form <- mod[[1]]
        return(all.vars(form$formula)[1])
      })))
  
  preds_long <- as.data.frame(preds) |> 
    dplyr::mutate(survey_id = rep(rownames(Y_test), times = nrow(preds)/nrow(Y_test))) |> 
    dplyr::mutate(variable = rep(var_to_predict, each = nrow(Y_test))) |> 
    dplyr::rename(imputed = V1)
    
  eval <- as.data.frame(Y_test) |> 
    tibble::rownames_to_column("survey_id") |> 
    tidyr::pivot_longer(cols = all_of(var_to_predict) ,
                        names_to = "variable",
                        values_to = "observed" ) |> 
    dplyr::full_join(preds_long) |> 
    dplyr::mutate(model = i)
    
  
  # ##obs
  # eval_act <- dplyr::filter(eval, variable == "mean_endemism")
  # eval_act <- dplyr::filter(eval, variable == "actino_richness")
  # plot(eval_act$imputed ~ eval_act$observed)
  # cor.test(eval_act$imputed, eval_act$observed)

  #Return dataframe
  eval
  
}) #END OF LAPPLY ON CROSSVALIDATION

# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models", "Pred_spaMM.jpg"),
       width = 12, height = 7)


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




