################################################################################
##
##  
##
## 3b_prediction_spaMM.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
################################################################################
##----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("spaMM")


##-------------loading data and functions-------------

## SITE SCALE
load( file = here::here("data", "derived_data", "3_sites_datasets_for_predict_CV_80_20.Rdata"))
load(file = here::here("data", "derived_data", "3_sites_covariates_to_predict.Rdata"))

#load functions
source("R/evaluation_prediction_model.R")


##------------- spaMM multivariate -------------


# cross_val <- lapply(1:length(datasets), FUN = function(i){
# cross_val <- lapply(1:2,FUN = function(i){

  i=1 # test on one iteration
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]][1:300,] ################################" reduce data if needed
  
  X_train <- covariates_final[rownames(Y_train),]
  
  train <- cbind(X_train, Y_train) |> 
    tibble::rownames_to_column("survey_id")
  
  covariates <- colnames(dplyr::select(X_train, -longitude, -latitude,
                                       -country, 
                                       -year,
                                       -voice,
                                       -no_violence))
  
  # List of submodels
  responses <- colnames(Y_train)[1:2] ###### reduce the number of responses if needed

  submodels <- list()
  
  for (response in responses) {
    form <- as.formula(
      paste(response, "~ ",
            paste(c("1", covariates 
                    , "(1|year)"
                    # , "(1|country)"
                    
                    # ,"Matern(1|longitude+latitude)"
                    , paste0("Matern(0+mv(",
                             paste(seq(1,length(responses)),collapse = ","),
                             ")|longitude+latitude)")

            ) , 
            collapse = "+")))
    
    submodels[[response]] <- list(formula = form, family = gaussian())
  }
  
  
  #Fit the model
  library(spaMM)
  model <- spaMM::fitmv(submodels,
                        data = train,
                        method = "ML",
                        verbose=c(TRACE=TRUE),
                        control.HLfit = list(
                          NbThreads = parallel::detectCores()-5))
    
  model
    
  # Predict on new data
  Y_test <- data[[2]][,responses]

  X_test <- covariates_final[rownames(Y_test),]

  test <- cbind(X_test, Y_test) |>
    tibble::rownames_to_column("survey_id")

  
  preds <- predict(model, newdata = test)
    
  

  
  #Compare prediction and observed values
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
  # eval_act <- dplyr::filter(eval, variable == "herbivores_biomass")
  # eval_act <- dplyr::filter(eval, variable == "actino_richness")
  # plot(eval_act$imputed ~ eval_act$observed)
  # cor.test(eval_act$imputed, eval_act$observed)
  
  #Return dataframe
  eval
  
# }) #END OF LAPPLY ON CROSSVALIDATION
