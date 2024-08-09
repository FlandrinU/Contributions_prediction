################################################################################
##
##  
##
## 3c_fit_model_spaMM.R
##
## 17/04/2024
##
## Ulysse Flandrin
##
################################################################################
##----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("spaMM", "DHARMa")
# install.packages(here::here("data", "raw_data","spaMM_4.4.40.tar.gz") )

library(spaMM)
##-------------loading data and functions-------------
#full data
load(file = here::here("data", "derived_data", "3_all_contributions_to_predict.Rdata"))

#covariates
load(file = here::here("data", "derived_data", "3_all_covariates_to_predict.Rdata"))


##------------- spaMM UNIVARIATE -------------

contributions <- colnames(observations_final)

fixed_effects <- colnames(dplyr::select(covariates_final,
                                           -longitude, -latitude, 
                                           -year,
                                           -country))

dataset <- cbind(observations_final, covariates_final)


#fit the model on each variable
models_uni <- #lapply(contributions,
              pbmcapply::pbmclapply(contributions, 
                                mc.cores = 10,
                                FUN = function(contrib){
                            
   # contrib = "herbivores_biomass"

  cat("contribution", which(contrib == contributions), "/", 
       length(contributions), ":", contrib, "\n")
   
   form <- as.formula(
     paste(
       contrib, "~ ", paste(c("1", 
                              fixed_effects
                              , "Matern(1|longitude+latitude)"
                              , "(1|year)"
       ),collapse = "+") ) )
   
   
   tryCatch({
     #Fit the model
     model <- spaMM::fitme(form,
                           data = dataset,
                           family = gaussian(),
                           control.HLfit = list(NbThreads = 1),
                           method = "ML"
     )

     model
     
   }, error = function(e){ model = NA }) #End of tryCatch
   
 }) #END OF MCLAPPLY ON EACH CONTRIBUTION
  

save(models_uni, file = here::here("outputs", "models", "Model_fit_spaMM_univariate.Rdata"))


##------------- spaMM UNIVARIATE witout spatial -------------

contributions <- colnames(observations_final)

fixed_effects <- colnames(dplyr::select(covariates_final,
                                        -longitude, -latitude, 
                                        -year,
                                        -country,
                                        -ecoregion,
                                        -realm))

dataset <- cbind(observations_final, covariates_final)


#fit the model on each variable
models_uni_without_spatial <- pbmcapply::pbmclapply(contributions, 
                                                    mc.cores = 5,
                                                    FUN = function(contrib){
   
   # contrib = "herbivores_biomass"
   
   cat("contribution", which(contrib == contributions), "/", 
       length(contributions), ":", contrib, "\n")
   
   form <- as.formula(
     paste(
       contrib, "~ ", paste(c("1", 
                              fixed_effects
                              , "(1|year)"
                              # , "(1|country)"
                              , "(1|ecoregion)"
       ),collapse = "+") ) )
   
   
   tryCatch({
     #Fit the model
     model <- spaMM::fitme(form,
                           data = dataset,
                           family = gaussian(),
                           control.HLfit = list(NbThreads = 1),
                           method = "ML"
     )
     
     model
     
   }, error = function(e){ model = NA }) #End of tryCatch
   
 }) #END OF MCLAPPLY ON EACH CONTRIBUTION


save(models_uni_without_spatial, 
     file = here::here("outputs", "models", "Model_fit_spaMM_univariate_without_spatial.Rdata"))




##------------- spaMM UNIVARIATE Central Indo-Pacific -------------

contributions <- colnames(observations_final)

fixed_effects <- colnames(dplyr::select(covariates_final,
                                        -longitude, -latitude, 
                                        -year,
                                        -country,
                                        -realm))

dataset <- cbind(observations_final, covariates_final) |> 
  dplyr::filter(realm == "Central Indo-Pacific")


#fit the model on each variable
models_uni_central_indo_pacific <- #lapply(contributions,
  pbmcapply::pbmclapply(contributions, 
                        mc.cores = 5,
                        FUN = function(contrib){
  
  # contrib = "herbivores_biomass"
  
  cat("contribution", which(contrib == contributions), "/", 
      length(contributions), ":", contrib, "\n")
  
  form <- as.formula(
    paste(
      contrib, "~ ", paste(c("1", 
                             fixed_effects
                             , "Matern(1|longitude+latitude)"
                             , "(1|year)"
      ),collapse = "+") ) )
  
  
  tryCatch({
    #Fit the model
    model <- spaMM::fitme(form,
                          data = dataset,
                          family = gaussian(),
                          control.HLfit = list(NbThreads = 1),
                          method = "ML"
    )
    
    model
    
  }, error = function(e){ model = NA }) #End of tryCatch
  
}) #END OF MCLAPPLY ON EACH CONTRIBUTION


save(models_uni_central_indo_pacific,
     file = here::here("outputs", "models", 
                       "Model_fit_spaMM_univariate_central_indo_pacific.Rdata"))



##------------- spaMM UNIVARIATE SUB-MODELS -------------

contributions <- colnames(observations_final)

fixed_effects <- colnames(dplyr::select(covariates_final,
                                        -longitude, -latitude, 
                                        -year,
                                        -country,
                                        -realm))
random_effect <- c( "Matern(1|longitude+latitude)" , "(1|year)")

predictors <- c(fixed_effects, random_effect)

dataset <- cbind(observations_final, covariates_final)

# models_without_one_predictor <- 
lapply(predictors, function(pred){
  # pred <- "Matern(1|longitude+latitude)"
  
  cat("Submodel", which(pred == predictors), "/", 
      length(predictors), ": remove", pred, "\n")
  
  predict <- predictors[predictors != pred]
  
  #fit the model on each response variable
  submodels_on_contributions <- pbmcapply::pbmclapply(contributions, 
                                                      mc.cores = 5,
                                                      FUN = function(contrib){
                            
      # contrib = "herbivores_biomass"
      
      form <- as.formula( paste( contrib, "~ ", paste(predict,collapse = "+") ) )
      
      tryCatch({
        #Fit the model
        model <- spaMM::fitme(form,
                              data = dataset,
                              family = gaussian(),
                              control.HLfit = list(NbThreads = 1),
                              method = "ML"
        )
        
        model
        
      }, error = function(e){ model = NA }) #End of tryCatch
      
    }) #END OF MCLAPPLY ON EACH CONTRIBUTION
  
  names(submodels_on_contributions) <- contributions
  
  save(submodels_on_contributions, 
       file = here::here("outputs", "models", 
       paste0("Model_fit_spaMM_univariate_WITHOUT_", pred,".Rdata")))
  
  rm(list = "submodels_on_contributions" )
  
})

# names(models_without_one_predictor) <- predictors
# 
# save(models_without_one_predictor, 
#      file = here::here("outputs", "models", "All_models_fit_spaMM_univariate_WITHOUT_ONE_predictor.Rdata"))



##------------- spaMM MULTIVARIATE -------------
# ## Correlation matrix
# cor_matrix <- cor(observations_final, method = "pearson")
# pairwise_corr <- cor_matrix[upper.tri(cor_matrix)]
# summary(pairwise_corr) 
# high_correlation_indices <- which(abs(cor_matrix) > 0.3 & cor_matrix < 1, arr.ind = TRUE)
# correlated_pairs <- cbind(rownames(cor_matrix)[high_correlation_indices[, "row"]],
#                           colnames(cor_matrix)[high_correlation_indices[, "col"]])
# unique(correlated_pairs[,1])
# 
# #Observe the highly correlated pairs:
# network_graph <- igraph::graph_from_edgelist(correlated_pairs, directed = FALSE)
# plot(network_graph, layout = igraph::layout_with_fr(network_graph),
#      edge.arrow.size = 1, vertex.label.cex = 0.8)


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

# Data
contributions <- colnames(observations_final)

fixed_effects <- colnames(dplyr::select(covariates_final,
                                        -longitude, -latitude, 
                                        -year,
                                        -country))

dataset <- cbind(observations_final, covariates_final)

# Check all contributions are considered
length(unlist(groups)) == length(contributions)


models_multi <- lapply(groups, FUN = function(group){
               # pbmcapply::pbmclapply(groups, mc.cores = 5, FUN = function(group){
  
  cat( "\n", "Multivariate fit with:", paste(group, collapse = " / "), "\n")
  
  # List of submodels
  #group = c("functional_entropy", "phylogenetic_entropy")
  submodels <- list()
  
  for (response in group) {
    form <- as.formula(
      paste(response, "~ ",
            paste(c("1", fixed_effects
                    , paste0("Matern( 0 + mv(",
                             paste(seq(1,length(group)),collapse = ","),
                             ")|longitude+latitude)")
                    , "(1|year)"
            ) , collapse = "+")))
  
    submodels[[response]] <- list(formula = form, family = gaussian())
  }
  
  
  tryCatch({
    #Fit the model
    model <- spaMM::fitmv(submodels,
                          data = dataset,
                          method = "ML",
                          verbose=c(TRACE=TRUE),
                          control.HLfit = list(NbThreads = 5))
    
    save(model, file = here::here("outputs", "models",
         paste0("Model_fit_spaMM_multivariate_", paste(group, collapse = "-"), ".Rdata")))
    
    model
  }, error = function(e){ model = NA }) 
  
}) #END OF LAPPLY ON groups of contributions

save(models_multi, file = here::here("outputs", "models", "Model_fit_spaMM_multivariate.Rdata"))


