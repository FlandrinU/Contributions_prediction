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
# pkgs <- c("spaMM", "DHARMa)


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


#Quick and good explanation on linear mixed effect models:
# https://www.youtube.com/watch?v=QCqF-2E86r0


##------------- spaMM multivariate -------------
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
  
  Y_train <- data[[1]]#[1:300,] ################################" reduce data
  
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
  
  responses <- colnames(Y_train)[5:12]
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





##------------- spaMM univariate -------------
distribution_plot(covariates_final, longer = T,
                  cols_not_plot = c("longitude", "latitude", "country",
                                    "ecoregion", "effectiveness"))
distribution_plot(observations_final, longer = T,
                  cols_plot = colnames(observations_final))

# cross_val <- pbmcapply::pbmclapply(1:length(datasets), mc.cores = 5,
#                                    FUN = function(i){
cross_val <- lapply(1:1,FUN = function(i){
  
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]]
  
  X_train <- covariates_final[rownames(Y_train),]
  
  train <- cbind(X_train, Y_train) |> 
    tibble::rownames_to_column("survey_id")
  
  covariates <- colnames(dplyr::select(X_train, -longitude, -latitude, -ecoregion, -year))
  
  #run model on each variable
  preds_raw <- lapply(colnames(Y_train),
              # pbmcapply::pbmclapply(colnames(Y_train), mc.cores = parallel::detectCores()-3,
                      FUN = function(contrib){
                        
    # contrib = "herbivores_biomass"
    cat("contribution", which(contrib == colnames(Y_train)), "/", 
        length(colnames(Y_train)), ":", contrib, "\n")
    
    form <- as.formula(
        paste(
          contrib, "~ ", paste(c("1", 
                                 covariates,
                                 "Matern(1|longitude+latitude)"
                                 , "(1|country)" #intercept random effect on country
                                 # # , "(effectiveness|country)"
                                 , "(1|year)"
                                 ),
                               collapse = "+")
              )
        )
    
    
    #New data
    Y_test <- data[[2]]
    
    X_test <- covariates_final[rownames(Y_test),]
    
    test <- cbind(X_test, Y_test) |> 
      tibble::rownames_to_column("survey_id")
    
    
    tryCatch({
      #Fit the model
      model <- spaMM::fitme(form,
                            data = train,
                            family = gaussian(),
                            control.HLfit = list(NbThreads = parallel::detectCores()-2),
                            method = "ML"
                            )


      #Predict on new data
      preds <- predict(model, newdata = test)
      preds
      }, error = function(e){
        data.frame(V1 = rep(NA, nrow(Y_test)))
        # error <<- data.frame(V1 = rep(NA, nrow(Y_test)))
        # error
      })#End of tryCatch
    
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
  "Pred_spaMM_univariate_methodML_random_effect_year_country.jpg"),
       width = 12, height = 7)


#density plot
density_prediction(all_res)
ggsave(filename = here::here("figures", "models", 
 "Pred_density_spaMM_univariate_methodML_random_effect_year_country.jpg"),
       width = 20, height = 10)



# #model quality
# library(DHARMa)
# summary(model)
# testDispersion(model)
# simulationOutput <- simulateResiduals(fittedModel = model, plot = F)
# residuals(simulationOutput)
# residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
# plot(simulationOutput)
# plotResiduals(simulationOutput, X_train$longitude)
# plotResiduals(simulationOutput, testData$Environment2)




##------------- spaMM variable importance -------------

# cross_val <- pbmcapply::pbmclapply(1:length(datasets), mc.cores = 5,
#                                    FUN = function(i){
cross_val <- lapply(1:2,FUN = function(i){
  
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]]
  
  X_train <- covariates_final[rownames(Y_train),]
  
  train <- cbind(X_train, Y_train) |> 
    tibble::rownames_to_column("survey_id")
  
  covariates <- colnames(dplyr::select(X_train, -longitude, -latitude))
  
  #run model on each variable
  preds_raw <- lapply(colnames(Y_train),
                      # pbmcapply::pbmclapply(colnames(Y_train), mc.cores = parallel::detectCores()-3,
                      FUN = function(contrib){
                        
    cat("contribution", which(contrib == colnames(Y_train)), "/", 
        length(colnames(Y_train)), ":", contrib, "\n")
    
    form <- as.formula(
      paste(
        contrib, "~ ", paste(c("1", covariates, 
                               "Matern(1|longitude+latitude)",
                               "(1|country)", #intercept random effect on country
                               "(1|ecoregion)",
                               "(1|year)"),
                             collapse = "+")
      )
    )
    
    
    #New data
    Y_test <- data[[2]]
    
    X_test <- covariates_final[rownames(Y_test),]
    
    test <- cbind(X_test, Y_test) |> 
      tibble::rownames_to_column("survey_id")
    
    
    #Fit the model
    model <- spaMM::fitme(form,
                          data = train,
                          family = gaussian(),
                          control.HLfit = list(NbThreads = parallel::detectCores()-5),
                          method = "ML"
    )
      
      
      #OBSERVE IMPORTANCE OF VARIABLES
      #Coefficient plot
      coefficients = summary(model, details=list(p_value=TRUE))$beta_table  |> 
        as.data.frame() |>
        tibble::rownames_to_column("term") |>
        janitor::clean_names() |>
        dplyr::mutate(significance = ifelse(p_value >= 0.001, "Not significant p-value","Significant p-value")) |>
        dplyr::mutate(t_value_signif = ifelse(t_value > 1.96 | t_value < - 1.96, "Significant t-value","Not significant t-value")) |>
        dplyr::filter(term != "(Intercept)") |> 
        dplyr::filter(!grepl("country", term)) |> 
        dplyr::filter(!grepl("ecoregion", term)) 
        
      
      #Plot
      library(ggplot2)
      (coefficient_plot = coefficients |>
          ggplot(aes(x=reorder(term, - estimate), y=estimate)) +
          geom_errorbar(aes(ymin=estimate-cond_se,ymax=estimate+cond_se),width=0.2)+
          geom_point(size = 3,aes(color = significance, shape = t_value_signif)) +
          harrypotter::scale_color_hp_d(option = "Ravenclaw") +
          labs(title = "(A) Coefficients of regression model fixed effects with 95% confidence intervals",
               color = "Significance of p-value (P < 0.001)",
               shape = "Significance of t-value (t < -1.96 or > 1.96)",
               x = " ",
               y = " ") +
          coord_flip()+
          theme_minimal(base_size = 15) +
          theme(legend.position = "top",legend.box="vertical", legend.margin=margin()) +
          theme(plot.caption = element_text(hjust = 0)))
      ggsave(filename = here::here("figures", "models","variable_importance",
                                   paste0("spaMM_univariate_", contrib, "_crossval_",
                                   i,".jpg")),
             width = 15, height = 10)
    
      
      #Predict on new data
      preds <- predict(model, newdata = test)
      preds
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
  
  ##obs
  eval_act <- dplyr::filter(eval, variable == "herbivores_biomass")
  plot(eval_act$imputed ~ eval_act$observed)
  cor.test(eval_act$imputed, eval_act$observed)
  
  #Return dataframe
  eval
  
}) #END OF LAPPLY ON CROSSVALIDATION

# Observe result
all_res <- do.call(rbind, cross_val)
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models",
                             "Pred_spaMM_univariate_methodML_random_effect_year_country_ecoregion.jpg"),
       width = 12, height = 7)


#density plot
density_prediction(all_res)
ggsave(filename = here::here("figures", "models", 
                             "Pred_density_spaMM_univariate_methodML_random_effect_year_country_ecoregion.jpg"),
       width = 20, height = 10)




