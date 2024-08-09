################################################################################
##
##  
##
## 3b_prediction_spaMM.R
##
## 23/10/2023
##
## Ulysse Flandrin
##
################################################################################
##----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("spaMM", "DHARMa")
# install.packages(here::here("data", "raw_data","spaMM_4.4.40.tar.gz") )

##-------------loading data and functions-------------
#full data
load(file = here::here("data", "derived_data", "3_all_contributions_to_predict.Rdata"))

#datasets to predict
load( file = here::here("data", "derived_data", "3_datasets_for_predict_CV_80_20.Rdata"))

#covariates
load(file = here::here("data", "derived_data", "3_all_covariates_to_predict.Rdata"))

# ## SITE SCALE
# load(file = here::here("data", "derived_data", "3_sites_contributions_to_predict.Rdata"))
# load( file = here::here("data", "derived_data", "3_sites_datasets_for_predict_CV_80_20.Rdata"))
# load(file = here::here("data", "derived_data", "3_sites_covariates_to_predict.Rdata"))

#load functions
source("R/evaluation_prediction_model.R")


#Quick and good explanation on linear mixed effect models:
# https://www.youtube.com/watch?v=QCqF-2E86r0


##------------- spaMM univariate + spatial-------------
distribution_plot(covariates_final, longer = T,
                  cols_not_plot = c("longitude", "latitude", 
                                    "country", "ecoregion", "realm",
                                    "effectiveness"))
distribution_plot(observations_final, longer = T,
                  cols_plot = colnames(observations_final))

# cross_val <- pbmcapply::pbmclapply(1:length(datasets), mc.cores = 5,
#                                    FUN = function(i){
cross_val <- lapply(1:1,FUN = function(i){
  
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]]#[1:2000,]
  
  X_train <- covariates_final[rownames(Y_train),]
  # |>
  #   dplyr::filter(realm == "Eastern Indo-Pacific")
  # Y_train <- Y_train[rownames(X_train),]
  
  train <- cbind(X_train, Y_train) |> 
    tibble::rownames_to_column("survey_id")
  
  covariates <- colnames(dplyr::select(X_train, -longitude, -latitude, 
                                       -year,
                                       -country,
                                       -ecoregion,
                                       -realm))
  
  #run model on each variable
  preds_raw <- #lapply(colnames(Y_train),
    pbmcapply::pbmclapply(colnames(Y_train), mc.cores = parallel::detectCores()-15,
                          FUN = function(contrib){
                            
     # contrib = "herbivores_biomass"
     # contrib = "piscivores_biomass"
     cat("contribution", which(contrib == colnames(Y_train)), "/", 
         length(colnames(Y_train)), ":", contrib, "\n")
     
     form <- as.formula(
       paste(
         contrib, "~ ", paste(c("1", 
                                covariates,
                                "Matern(1|longitude+latitude)"
                                # , "(1|country)" #intercept random effect on country
                                # # , "(effectiveness|country)"
                                , "(1|year)"
                                # , "(1|depth)"
         ),
         collapse = "+")
       )
     )
     
     
     #New data
     Y_test <- data[[2]]
     
     X_test <- covariates_final[rownames(Y_test),]
     # |>
     #   dplyr::filter(realm == "Eastern Indo-Pacific")
     # Y_test <- Y_test[rownames(X_test),]
     
     test <- cbind(X_test, Y_test) |> 
       tibble::rownames_to_column("survey_id")
     
     
     tryCatch({
       #Fit the model
       library(spaMM)
       model <- spaMM::fitme(form,
                             data = train,
                             family = gaussian(),
                             control.HLfit = list(NbThreads = 1),
                             method = "ML"
       )
       
       # save(model, file = here::here("outputs", "models",
       #                               paste0("spaMM_UNIvariate_2000_rows_", contrib, ".Rdata")))
       
       
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
  
  Y_test <- data[[2]]
  
  X_test <- covariates_final[rownames(Y_test),]
  # |>
  #   dplyr::filter(realm == "Eastern Indo-Pacific")
  # Y_test <- Y_test[rownames(X_test),]
  
  #extract result
  preds <- do.call(cbind, preds_raw)
  colnames(preds) <- colnames(Y_test)
  rownames(preds) <- rownames(Y_test)
  
  preds_long <- as.data.frame(preds) |> 
    tibble::rownames_to_column("survey_id") |> 
    tidyr::pivot_longer(cols = -survey_id ,
                        names_to = "variable",
                        values_to = "imputed" )
  
  eval <- data[[2]] |>
    # eval <- data[[1]] |>  # try on train
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
save(all_res, file = here::here("outputs", "models",
                                "CROSS_VAL_spaMM_univariate_year_spatial_Rdeffect.Rdata"))
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models",
                             "Pred_spaMM_univariate_year_spatial_Rdeffect.jpg"),
       width = 12, height = 7)


#density plot
plot_interaction(all_res, var_facet_wrap = "variable", X_values = "observed",
                 Y_values = "imputed", xlabel = "Observed contributions",
                 ylabel = "Predicted contributions")
density_prediction(all_res)
ggsave(filename = here::here("figures", "models", 
                             "Pred_density_spaMM_univariate_random_effect_year_spatial.jpg"),
       width = 20, height = 10)






##------------- spaMM univariate without spatial-------------
distribution_plot(covariates_final, longer = T,
                  cols_not_plot = c("longitude", "latitude", 
                                    "country", "ecoregion", "realm",
                                    "effectiveness"))
distribution_plot(observations_final, longer = T,
                  cols_plot = colnames(observations_final))

cross_val <- pbmcapply::pbmclapply(1:length(datasets), mc.cores = 5,
                                   FUN = function(i){
# cross_val <- lapply(1:10,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]]
  
  X_train <- covariates_final[rownames(Y_train),]
  
  
  train <- cbind(X_train, Y_train) |> 
    tibble::rownames_to_column("survey_id")
  
  covariates <- colnames(dplyr::select(X_train, -longitude, -latitude, 
                                       -year,
                                       -country,
                                       -ecoregion,
                                       -realm))
  
  #run model on each variable
  preds_raw <- #lapply(colnames(Y_train),
    pbmcapply::pbmclapply(colnames(Y_train), mc.cores = parallel::detectCores()-15,
                          FUN = function(contrib){
                            
  # contrib = "herbivores_biomass"
  # contrib = "piscivores_biomass"
  cat("contribution", which(contrib == colnames(Y_train)), "/", 
      length(colnames(Y_train)), ":", contrib, "\n")
  
  form <- as.formula(
    paste(
      contrib, "~ ", paste(c("1", 
                             covariates
                             # ,"Matern(1|longitude+latitude)"
                             # , "(1|country)" #intercept random effect on country
                             , "(1|year)"
                            , "(1|ecoregion)"
                            # ,"(1|realm)"
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
    library(spaMM)
    model <- spaMM::fitme(form,
                          data = train,
                          family = gaussian(),
                          control.HLfit = list(NbThreads = 1),
                          method = "ML"
    )
    
    # save(model, file = here::here("outputs", "models",
    #                               paste0("spaMM_UNIvariate_2000_rows_", contrib, ".Rdata")))
    
    
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
  
  Y_test <- data[[2]]
  
  X_test <- covariates_final[rownames(Y_test),]

  
  #extract result
  preds <- do.call(cbind, preds_raw)
  colnames(preds) <- colnames(Y_test)
  rownames(preds) <- rownames(Y_test)
  
  preds_long <- as.data.frame(preds) |> 
    tibble::rownames_to_column("survey_id") |> 
    tidyr::pivot_longer(cols = -survey_id ,
                        names_to = "variable",
                        values_to = "imputed" )
  
  eval <- data[[2]] |>
    # eval <- data[[1]] |>  # try on train
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
ggsave(filename = here::here("figures", "models", "spaMM_univar_without_spatial",
                             "Pred_spaMM_univar_Rdeffect(year-ecoregion).jpg"),
       width = 12, height = 7)


#density plot
plot_interaction(all_res, var_facet_wrap = "variable", X_values = "observed",
                 Y_values = "imputed", xlabel = "Observed contributions",
                 ylabel = "Predicted contributions")
density_prediction(all_res)
ggsave(filename = here::here("figures", "models", "spaMM_univar_without_spatial",
                             "Pred_density_spaMM_univar_Rdeffect(year-ecoregion).jpg"),
       width = 20, height = 10)






##------------- spaMM multivariate -------------
# install.packages("spaMM")

#IF IMPOSSIBLE TO LOAD PKG AFTER VERSION 4.0.0
# install.packages("data/raw_data/spaMM_3.13.0.tar.gz", repos = NULL, type = "source")

# library(spaMM)


## Correlation matrix
cor_matrix <- cor(observations_final, method = "pearson")
pairwise_corr <- cor_matrix[upper.tri(cor_matrix)]
summary(pairwise_corr) 
high_correlation_indices <- which(abs(cor_matrix) > 0.5 & cor_matrix < 1, arr.ind = TRUE)
correlated_pairs <- cbind(rownames(cor_matrix)[high_correlation_indices[, "row"]],
                          colnames(cor_matrix)[high_correlation_indices[, "col"]])
unique(correlated_pairs[,1])

#Observe the highly correlated pairs:
network_graph <- igraph::graph_from_edgelist(correlated_pairs, directed = FALSE)
plot(network_graph, layout = igraph::layout_with_fr(network_graph),
     edge.arrow.size = 1, vertex.label.cex = 0.8)


groups <- list(
  c("functional_entropy", "phylogenetic_entropy"),
  
  c("trophic_web_robustness", "vitamin_A", "mean_trophic_level"),
  
  c("calcium", "NP_score", "available_biomass_turnover", "iron", "zinc"),
  
  c("actino_richness", "mean_endemism", "omega_3", "functional_distinctiveness",
    "evolutionary_distinctiveness", "selenium"),
  
  c("herbivores_biomass", "available_biomass", "invertivores_biomass",
    "NN_score", "piscivores_biomass", "iucn_species_richness")
)

#test group
group <- groups[[5]]

# cross_val <- pbmcapply::pbmclapply(1:length(datasets), mc.cores = 5,
#                                    FUN = function(i){
# cross_val <- lapply(1:1,FUN = function(i){
cross_val <- pbmcapply::pbmclapply(1:2, mc.cores = 1,FUN = function(i){
  
  # i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]][1:2000,] ################################" reduce data
  
  X_train <- covariates_final[rownames(Y_train),]
  
  train <- cbind(X_train, Y_train) |> 
    tibble::rownames_to_column("survey_id")
  
  covariates <- colnames(dplyr::select(X_train, -longitude, -latitude,
                                       -country, 
                                       -year))
  
  
  # List of submodels

  # responses <- colnames(Y_train)
  responses <- group
  # responses <- c("actino_richness", "mean_endemism", "available_biomass_turnover",
  #                "NN_score", "iron", "available_biomass", "functional_entropy")##############################"
  submodels <- list()
  
  for (response in responses) {
    form <- as.formula(
      paste(response, "~ ",
            paste(c("1", covariates
                    # , "Matern(1|longitude+latitude)"
                    , paste0("Matern( 0 + mv(",
                             paste(seq(1,length(responses)),collapse = ","),
                             ")|longitude+latitude)")
                    , "(1|year)"
                    # , "(1|country)"
                    # , paste0("(0+spaMM::mv(",
                    #          paste(seq(1,length(responses)),collapse = ","),
                    #          ")|survey_id)")
            ) , 
            collapse = "+")))
    
    submodels[[response]] <- list(formula = form, family = gaussian())
  }
  
  
tryCatch({
    #Fit the model
    library(spaMM) #spaMM need to be attached to hahe the verbose.
    model <- spaMM::fitmv(submodels,
                          data = train,
                          method = "ML",
                          verbose=c(TRACE=TRUE),
                          control.HLfit = list(
                            NbThreads = 10))
                            # NbThreads = parallel::detectCores()-5))
    
    save(model, file = here::here("outputs", "models",
                                  paste0("spaMM_multivariate_2000_rows_", paste(group, collapse ="-"), ".Rdata")))
    

    #New data
    Y_test <- data[[2]][,responses]
    
    X_test <- covariates_final[rownames(Y_test), ] |> 
      # dplyr::select(-all_of(useless))
      dplyr::filter(year %in% unique(X_train$year) &
                      country %in% unique(X_train$country))
    
    Y_test <- Y_test[rownames(X_test),]
    
    test <- cbind(X_test, Y_test) |> 
      tibble::rownames_to_column("survey_id")
    
    #Predict on new data
    preds <- predict(model, newdata = test) 

    }, error = function(e){
      
      preds <<- data.frame(V1 = rep(NA, nrow(Y_test) * ncol(Y_test)))

  })#End of tryCatch
    
  # Distriburtion of predictions
  # distibution_plot(preds, )
  
  #Compare prediction
  var_to_predict <-  names(
    unlist( 
      lapply(submodels, function(mod) {
        form <- mod[[1]]
        return(all.vars(form$formula)[1])
      })))
  
  preds_long <- as.data.frame(preds) |> 
    dplyr::mutate(survey_id = rep(rownames(Y_train), times = nrow(preds)/nrow(Y_train))) |> 
    dplyr::mutate(variable = rep(var_to_predict, each = nrow(Y_train))) |> 
    dplyr::rename(imputed = V1)
    
  eval <- as.data.frame(Y_train) |> 
    tibble::rownames_to_column("survey_id") |> 
    tidyr::pivot_longer(cols = all_of(var_to_predict) ,
                        names_to = "variable",
                        values_to = "observed" ) |> 
    dplyr::full_join(preds_long) |> 
    dplyr::mutate(model = i)
    
  
  # ##obs
  # eval_act <- dplyr::filter(eval, variable == "functional_entropy")
  # eval_act <- dplyr::filter(eval, variable == "actino_richness")
  # plot(eval_act$imputed ~ eval_act$observed)
  # cor.test(eval_act$imputed, eval_act$observed)

  #Return dataframe
  eval
  
}) #END OF LAPPLY ON CROSSVALIDATION

# Observe result
all_res <- do.call(rbind, cross_val)
save(all_res, file = here::here("outputs", "models",
     aste0("CROSS_VAL_spaMM_multivariate_3000_rows_", paste(group, collapse ="-"), ".Rdata")))
result <- extract_result_contrib(cross_val)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models",
                             paste0("Pred_spaMM_multivariate_3000_rows_", paste(group, collapse ="-"),
                             ".jpg")),
       width = 12, height = 7)

#observe distribution of predicted values
all_res_poisson <- eval
res_test <- all_res_poisson[all_res_poisson$variable == "iucn_species_richness",]
plot(res_test$imputed~res_test$observed)
cor.test(res_test$imputed, res_test$observed)

all_res_gaussian<- eval
res_test <- all_res_gaussian[all_res_gaussian$variable == "iucn_species_richness",]
plot(res_test$imputed~res_test$observed)
cor.test(res_test$imputed, res_test$observed)

all_res <- all_res_poisson
all_res <- all_res_gaussian

ggplot(all_res)+
  geom_histogram(aes(x=imputed, group=variable, y = after_stat(density)), 
                 bins = 20, color = "grey40", fill ="white") +
  geom_density(aes(x=observed, group=variable, fill = variable), alpha = 0.2) +
  hrbrthemes::theme_ipsum(axis_title_size = 11, strip_text_size = 14,
                          strip_text_face = "bold",) +
  xlab("Contribution values") + 
  ylab("Histogram of imputed values and density of observed values")+
  facet_wrap(~variable, scales = "free") +
  theme(legend.position="none",panel.spacing = unit(0.1, "lines"),
        axis.ticks.x=element_blank())


#observe residuals
residuals <- all_res |> 
  dplyr::mutate(residuals = observed - imputed)

distribution_plot(residuals, longer = F, index_values = c("variable", "residuals"))





##------------- spaMM multivariate without spatial -> no interest =  univariate spaMM-------------

groups <- list(
  c("functional_entropy", "phylogenetic_entropy"),
  
  c("NN_score", "NP_score"),
  
  c("trophic_web_robustness", "vitamin_A", "mean_trophic_level"),
  
  c("calcium", "available_biomass_turnover", "iron", "zinc"),
  
  c("herbivores_biomass", "available_biomass", "invertivores_biomass",
    "piscivores_biomass", "iucn_species_richness"),
  
  c("actino_richness", "mean_endemism", "omega_3", "functional_distinctiveness",
    "evolutionary_distinctiveness", "selenium")
)


models_multi <- lapply(groups, FUN = function(group){
  # pbmcapply::pbmclapply(groups, mc.cores = 5, FUN = function(group){
  
  cat( "\n", "Multivariate fit with:", paste(group, collapse = " / "), "\n")
  
  # List of submodels
  #group = c("functional_entropy", "phylogenetic_entropy")
  submodels <- list()
  
# # cross_val <- lapply(1:1,FUN = function(i){
# cross_val <- pbmcapply::pbmclapply(1:2, mc.cores = 1,FUN = function(i){
  
  i=1
  cat("Crossvalidation ", i,"/", length(datasets), "\n")
  data = datasets[[i]] 
  
  Y_train <- data[[1]]#[1:2000,] ################################" reduce data
  
  X_train <- covariates_final[rownames(Y_train),]
  
  train <- cbind(X_train, Y_train) |> 
    tibble::rownames_to_column("survey_id")
  
  covariates <- colnames(dplyr::select(X_train, -longitude, -latitude,
                                       -country,
                                       -ecoregion,
                                       -realm,
                                       -year))
  
  
  # List of submodels
  
  # responses <- colnames(Y_train)
  responses <- group
  # responses <- c("actino_richness", "mean_endemism", "available_biomass_turnover",
  #                "NN_score", "iron", "available_biomass", "functional_entropy")##############################"
  submodels <- list()
  
  for (response in responses) {
    form <- as.formula(
      paste(response, "~ ",
            paste(c("1", covariates
                    # , paste0("Matern( 0 + mv(",
                    #          paste(seq(1,length(responses)),collapse = ","),
                    #          ")|longitude+latitude)")
                    , "(1|year)"
                    , "(1|ecoregion)"
            ) , 
            collapse = "+")))
    
    submodels[[response]] <- list(formula = form, family = gaussian())
  }
  
  
  tryCatch({
    #Fit the model
    library(spaMM) #spaMM need to be attached to hahe the verbose.
    model <- spaMM::fitmv(submodels,
                          data = train,
                          method = "ML",
                          # verbose=c(TRACE=TRUE),
                          control.HLfit = list(NbThreads = 5))
    # NbThreads = parallel::detectCores()-5))
    
    # save(model, file = here::here("outputs", "models",
    #                               paste0("spaMM_multivariate_2000_rows_", paste(group, collapse ="-"), ".Rdata")))
    
    
    #New data
    Y_test <- data[[2]][,responses]
    
    X_test <- covariates_final[rownames(Y_test), ] |> 
      # dplyr::select(-all_of(useless))
      dplyr::filter(year %in% unique(X_train$year) &
                      ecoregion %in% unique(X_train$ecoregion))
    
    Y_test <- Y_test[rownames(X_test),]
    
    test <- cbind(X_test, Y_test) |> 
      tibble::rownames_to_column("survey_id")
    
    #Predict on new data
    preds <- predict(model, newdata = test) 
    
  }, error = function(e){
    
    preds <<- data.frame(V1 = rep(NA, nrow(Y_test) * ncol(Y_test)))
    
  })#End of tryCatch
  
  # Distriburtion of predictions
  # distibution_plot(preds, )
  
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
  # eval_act <- dplyr::filter(eval, variable == "herbivores_biomass")
  # eval_act <- dplyr::filter(eval, variable == "actino_richness")
  # plot(eval_act$imputed ~ eval_act$observed)
  # cor.test(eval_act$imputed, eval_act$observed)
  
  #Return dataframe
  eval
  
}) #END OF LAPPLY ON groups of contributions

# Observe result
all_res <- do.call(rbind, models_multi)
result <- extract_result_contrib(models_multi)
estimates_boxplot(result)
ggsave(filename = here::here("figures", "models", "spaMM_multivariate",
                             "Pred_spaMM_multivariate_without_spatial.jpg"),
       width = 12, height = 7)