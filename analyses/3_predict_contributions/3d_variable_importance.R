################################################################################
##
##  
##
## 3d_variable_importance.R
##
## 18/04/2024
##
## Ulysse Flandrin
##
################################################################################
##----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("DHARMa")

##-------------loading data and functions-------------
#full model spaMM univariate
load(file = here::here("outputs", "models", "Model_fit_spaMM_univariate.Rdata"))

# spaMM univariate without spatial in Matern(...)
load(file = here::here("outputs", "models", "Model_fit_spaMM_univariate_without_spatial.Rdata"))

# spaMM univariate without Rd effect
load(file = here::here("outputs", "models", "Model_fit_spaMM_univariate_WITHOUT_rd_effect.Rdata"))

# spaMM univariate in central indo pacific only
load(file = here::here("outputs", "models", "Model_fit_spaMM_univariate_central_indo_pacific.Rdata"))

#full model spaMM multivariate
load( file = here::here("outputs", "models", "Model_fit_spaMM_multivariate.Rdata"))
# load( file = here::here("outputs", "models", "Model_fit_spaMM_multivariate_calcium-available_biomass_turnover-iron-zinc.Rdata"))

#covariates
load(file = here::here("data", "derived_data", "3_all_covariates_to_predict.Rdata"))

#load functions
source("R/evaluation_prediction_model.R")


###############################################################################"
##
##                         #### CHECK MODELS ####
##
###############################################################################"

##---------------------- spaMM univariate ----------------------


# Analyse each models

results <- do.call(rbind, pbmcapply::pbmclapply(models_uni, 
                                                mc.cores = 10,
                                                function(model){
  # model <- models_uni[[1]]
   preds <- predict(model, newdata = covariates_final)
  
   #observe residuals
   residuals <- data.frame(variable =  as.character(model[["predictor"]][[2]]),
                           observed = model[["data"]][[model[["predictor"]][[2]]]],
                           imputed = preds,
                           res = spaMM::residuals.HLfit(model) ) |>   #residuals from spaMM package
      tibble::rownames_to_column("survey_id")
   
   residuals
}))

## Predictions distributions -> same patterns as observations ?
distribution_plot(results, longer = F, index_values = c("variable", "imputed"), 
                  xlabel = "Predictions", ylabel = "frequency" )
ggsave(filename = here::here("figures", "models", "spaMM_univariate",
                             "Distributions_of_predicted_contributions.jpg"),
       width = 20, height = 10)


## Predictions vs observations -> model quality ?

# trace(ggpubr:::.stat_lm, edit = TRUE) #change signif() in line 9 to have more digits
ggplot(results)+
  geom_point(aes(x = observed, y = imputed, fill = variable),
             color = "grey40", alpha = 0.2, shape = 21) +
  hrbrthemes::theme_ipsum() +
  xlab("Observed contributions") + ylab("imputed")+
  geom_abline(slope = 1) + 
  # annotate("text", x = 2.2, y = 1.8, label = "1:1")+
  ggpubr::stat_regline_equation(data = results,
                                aes(x = observed, y = imputed, label = ..rr.label..))   +
  facet_wrap(~variable, scales = "free") +
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"))

# plot_interaction(results, var_facet_wrap = "variable", X_values = "observed",
#                 Y_values = "imputed", xlabel = "Observed contributions",
#                 ylabel = "imputed") + geom_abline(slope = 1)

ggsave(filename = here::here("figures", "models", "spaMM_univariate",
                             "Predictions_VS_observations.jpg"),
       width = 20, height = 10)


## Residuals distributions -> Normal ? zero-centered ?
distribution_plot(results, longer = F, index_values = c("variable", "res"), 
                 xlabel = "Residuals", ylabel = "frequency" )
ggsave(filename = here::here("figures", "models", "spaMM_univariate",
                             "Residuals_distributions.jpg"),
       width = 20, height = 10)


## QQ-plot of residuals -> Normality of residuals ?
jpeg(filename = here::here("figures", "models", "spaMM_univariate", "residuals_QQplots.jpg"),
     units = "cm", width = 40, height = 30, res = 400 )
par(mfrow = c(5,5))
for(model in models_uni){
  # model = models_uni[[1]]
  simulationOutput <- DHARMa::simulateResiduals(fittedModel = model, plot = F)
  DHARMa::plotQQunif(simulationOutput) 
  title(main = model[["predictor"]][[2]] , line = -0.5, cex.main = 1.5)
}
dev.off()


## Residuals according to Y values -> Homoscedasticity of residuals ?
plot_interaction(results, var_facet_wrap = "variable", X_values = "imputed",
                 Y_values = "res", xlabel = "Predicted contributions",
                 ylabel = "Residuals") + geom_hline(yintercept = 0)
ggsave(filename = here::here("figures", "models", "spaMM_univariate",
                             "Residuals_VS_predicted_contributions.jpg"),
       width = 20, height = 10)


## Residuals according to covariates -> Patterns ?
colnames(covariates_final)
cov = "latitude"
res_cov <- results |> 
  dplyr::left_join( tibble::rownames_to_column(covariates_final, "survey_id"))
plot_interaction(res_cov, var_facet_wrap = "variable", X_values = cov,
                 Y_values = "res", xlabel = cov, ylabel = "Residuals")






## Simulated residuals vs predicted
jpeg(filename = here::here("figures", "models", "spaMM_univariate", "residualsVSpredicted.jpg"),
     units = "cm", width = 40, height = 25, res = 400 )
  par(mfrow = c(5,5))
  for(model in models_uni){
    # model = models_uni[[1]]
    simulationOutput <- DHARMa::simulateResiduals(fittedModel = model, plot = F)
    DHARMa::plotResiduals(simulationOutput) 
    title(main = model[["predictor"]][[2]] , line = 0.5, cex.main = 1.5)
  }
dev.off()

   
   
# DHARMa::plotResiduals(simulationOutput, covariates_final$latitude)
   
   


##---------------------- spaMM multivariates ----------------------


# Analyse each models

results <- do.call(rbind,lapply(models_multi, function(model){
  
  # model <- models_multi[[5]]
  preds <- predict(model, newdata = covariates_final) 
  
  #observe residuals
  responses <- c()
  for(i in 1:length(model[["predictor"]])){
    responses <- c(responses, as.character(model[["predictor"]][[i]][[2]]))
  }
  
  observed <- lapply(responses, function(i) model[["data"]][[i]] ) 
  
  residuals <- data.frame(variable =  rep(responses, each = length(observed[[1]])),
                          observed = unlist(observed),
                          imputed = preds,
                          res = spaMM::residuals.HLfit(model) ) |>   #residuals from spaMM package
    tibble::rownames_to_column("survey_id") |> 
    dplyr::mutate(survey_id = gsub("[^0-9.]", "", survey_id)) |> 
    dplyr::mutate(survey_id = gsub("\\..*", "", survey_id))
    
  residuals
}))
# results = residuals

## Predictions distributions -> same patterns as observations ?
distribution_plot(results, longer = F, index_values = c("variable", "imputed"), 
                  xlabel = "Predictions", ylabel = "frequency" )
ggsave(filename = here::here("figures", "models", "spaMM_multivariate",
                             "Distributions_of_predicted_contributions.jpg"),
       width = 20, height = 10)


## Predictions vs observations -> model quality ?

# trace(ggpubr:::.stat_lm, edit = TRUE) #change signif() in line 9 to have more digits
ggplot(results)+
  geom_point(aes(x = observed, y = imputed, fill = variable),
             color = "grey40", alpha = 0.2, shape = 21) +
  hrbrthemes::theme_ipsum() +
  xlab("Observed contributions") + ylab("imputed")+
  geom_abline(slope = 1) + 
  # annotate("text", x = 2.2, y = 2, label = "1:1")+
  ggpubr::stat_regline_equation(data = results,
              aes(x = observed, y = imputed, label = ..rr.label..))   +
  facet_wrap(~variable, scales = "free") +
  theme(legend.position="none",panel.spacing = unit(0.1, "lines"))

# plot_interaction(results, var_facet_wrap = "variable", X_values = "observed",
#                  Y_values = "imputed", xlabel = "Observed contributions",
#                  ylabel = "imputed") 

ggsave(filename = here::here("figures", "models", "spaMM_multivariate",
                             "Predictions_VS_observations.jpg"),
       width = 20, height = 10)


## Residuals distributions -> Normal ? zero-centered ?
distribution_plot(results, longer = F, index_values = c("variable", "res"), 
                  xlabel = "Residuals", ylabel = "frequency" )
ggsave(filename = here::here("figures", "models", "spaMM_multivariate",
                             "Residuals_distributions.jpg"),
       width = 20, height = 10)


# ## QQ-plot of residuals -> Normality of residuals ?
# jpeg(filename = here::here("figures", "models", "spaMM_multivariate", "residuals_QQplots.jpg"),
#      units = "cm", width = 40, height = 30, res = 400 )
# par(mfrow = c(5,5))
# for(model in models_multi){
#   # model = models_multi[[1]]
#   simulationOutput <- DHARMa::simulateResiduals(fittedModel = model, plot = F)
#   DHARMa::plotQQunif(simulationOutput) 
#   title(main = model[["predictor"]][[2]] , line = -0.5, cex.main = 1.5)
# }
# dev.off()


## Residuals according to Y values -> Homoscedasticity of residuals ?
plot_interaction(results, var_facet_wrap = "variable", X_values = "imputed",
                 Y_values = "res", xlabel = "Predicted contributions",
                 ylabel = "Residuals") + geom_hline(yintercept = 0)
ggsave(filename = here::here("figures", "models", "spaMM_multivariate",
                             "Residuals_VS_predicted_contributions.jpg"),
       width = 20, height = 10)



## Residuals according to covariates -> Patterns ?
colnames(covariates_final)
cov = "latitude"
res_cov <- results |> 
  dplyr::left_join( tibble::rownames_to_column(covariates_final, "survey_id"))
plot_interaction(res_cov, var_facet_wrap = "variable", X_values = cov,
                 Y_values = "res", xlabel = cov, ylabel = "Residuals")



# ## Simulated residuals vs predicted
# jpeg(filename = here::here("figures", "models", "spaMM_multivariate", "residualsVSpredicted.jpg"),
#      units = "cm", width = 40, height = 25, res = 400 )
# par(mfrow = c(5,5))
# for(model in models_uni){
#   # model = models_multi[[1]]
#   simulationOutput <- DHARMa::simulateResiduals(fittedModel = model, plot = F)
#   DHARMa::plotResiduals(simulationOutput) 
#   title(main = model[["predictor"]][[2]] , line = 0.5, cex.main = 1.5)
# }
# dev.off()



# DHARMa::plotResiduals(simulationOutput, covariates_final$latitude)



##---------------------- spaMM univar WITHOUT spatial ----------------------


# Analyse each models

results <- do.call(rbind, pbmcapply::pbmclapply(models_uni_without_spatial, 
                                                mc.cores = 10,
                                                function(model){
  # model <- models_uni_without_spatial[[1]]
  preds <- predict(model, newdata = covariates_final)
  
  #observe residuals
  residuals <- data.frame(variable =  as.character(model[["predictor"]][[2]]),
                          observed = model[["data"]][[model[["predictor"]][[2]]]],
                          imputed = preds,
                          res = spaMM::residuals.HLfit(model) ) |>   #residuals from spaMM package
    tibble::rownames_to_column("survey_id")
  
  residuals
}))

#Mean R-squared
results |> 
  dplyr::group_by(variable) |> 
  dplyr::summarise(r2 = summary(lm(imputed ~ observed))$r.squared) |> 
  dplyr::pull(r2) |> 
  summary()
# Min.   1st Qu.  Median    Mean  3rd Qu.   Max. 
# 0.1804  0.2410  0.2950  0.3657  0.3342  0.9400


## Predictions distributions -> same patterns as observations ?
distribution_plot(results, longer = F, index_values = c("variable", "imputed"), 
                  xlabel = "Predictions", ylabel = "frequency" )
ggsave(filename = here::here("figures", "models", "spaMM_univar_without_spatial",
                             "Distributions_of_predicted_contributions.jpg"),
       width = 20, height = 10)



## Predictions vs observations -> model quality ?

# trace(ggpubr:::.stat_lm, edit = TRUE) #change signif() in line 9 to have more digits
ggplot(results)+
  geom_point(aes(x = observed, y = imputed, fill = variable),
             color = "grey40", alpha = 0.2, shape = 21) +
  hrbrthemes::theme_ipsum() +
  xlab("Observed contributions") + ylab("imputed")+
  geom_abline(slope = 1) + 
  # annotate("text", x = 2.2, y = 1.8, label = "1:1")+
  ggpubr::stat_regline_equation(data = results,
                                aes(x = observed, y = imputed, label = ..rr.label..))   +
  facet_wrap(~variable, scales = "free") +
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"))

ggsave(filename = here::here("figures", "models", "spaMM_univar_without_spatial",
                             "Predictions_VS_observations.jpg"),
       width = 20, height = 10)



## Residuals distributions -> Normal ? zero-centered ?
distribution_plot(results, longer = F, index_values = c("variable", "res"), 
                  xlabel = "Residuals", ylabel = "frequency" )
ggsave(filename = here::here("figures", "models", "spaMM_univar_without_spatial",
                             "Residuals_distributions.jpg"),
       width = 20, height = 10)



## QQ-plot of residuals -> Normality of residuals ?
jpeg(filename = here::here("figures", "models", "spaMM_univar_without_spatial", 
    "residuals_QQplots.jpg"), units = "cm", width = 40, height = 30, res = 400 )
par(mfrow = c(5,5))
for(model in models_uni_without_spatial){
  # model = models_uni[[1]]
  simulationOutput <- DHARMa::simulateResiduals(fittedModel = model, plot = F)
  DHARMa::plotQQunif(simulationOutput) 
  title(main = model[["predictor"]][[2]] , line = -0.5, cex.main = 1.5)
}
dev.off()



## Residuals according to Y values -> Homoscedasticity of residuals ?
plot_interaction(results, var_facet_wrap = "variable", X_values = "imputed",
                 Y_values = "res", xlabel = "Predicted contributions",
                 ylabel = "Residuals") + geom_hline(yintercept = 0)
ggsave(filename = here::here("figures", "models", "spaMM_univar_without_spatial",
                             "Residuals_VS_predicted_contributions.jpg"),
       width = 20, height = 10)


## Residuals according to covariates -> Patterns ?
colnames(covariates_final)
cov = "latitude"
res_cov <- results |> 
  dplyr::left_join( tibble::rownames_to_column(covariates_final, "survey_id"))
plot_interaction(res_cov, var_facet_wrap = "variable", X_values = cov,
                 Y_values = "res", xlabel = cov, ylabel = "Residuals")




## Simulated residuals vs predicted
jpeg(filename = here::here("figures", "models", "spaMM_univar_without_spatial", 
    "residualsVSpredicted.jpg"), units = "cm", width = 40, height = 25, res = 400 )
par(mfrow = c(5,5))
for(model in models_uni_without_spatial){
  # model = models_uni[[1]]
  simulationOutput <- DHARMa::simulateResiduals(fittedModel = model, plot = F)
  DHARMa::plotResiduals(simulationOutput) 
  title(main = model[["predictor"]][[2]] , line = 0.5, cex.main = 1.5)
}
dev.off()





##---------------------- spaMM univar WITHOUT Rd effect ----------------------


# Analyse each models

results <- do.call(rbind, pbmcapply::pbmclapply(models_uni_without_rd_effect, 
                                                mc.cores = 10,
                                                function(model){
  # model <- models_uni[[1]]
  preds <- predict(model, newdata = covariates_final)
  
  #observe residuals
  residuals <- data.frame(variable =  as.character(model[["predictor"]][[2]]),
                          observed = model[["data"]][[model[["predictor"]][[2]]]],
                          imputed = preds,
                          res = spaMM::residuals.HLfit(model) ) |>   #residuals from spaMM package
    tibble::rownames_to_column("survey_id")
  
  residuals
}))


## Predictions distributions -> same patterns as observations ?
distribution_plot(results, longer = F, index_values = c("variable", "imputed"), 
                  xlabel = "Predictions", ylabel = "frequency" )
ggsave(filename = here::here("figures", "models", "spaMM_univariate_without_Random",
                             "Distributions_of_predicted_contributions.jpg"),
       width = 20, height = 10)



## Predictions vs observations -> model quality ?

# trace(ggpubr:::.stat_lm, edit = TRUE) #change signif() in line 9 to have more digits
ggplot(results)+
  geom_point(aes(x = observed, y = imputed, fill = variable),
             color = "grey40", alpha = 0.2, shape = 21) +
  hrbrthemes::theme_ipsum() +
  xlab("Observed contributions") + ylab("imputed")+
  geom_abline(slope = 1) + 
  # annotate("text", x = 2.2, y = 1.8, label = "1:1")+
  ggpubr::stat_regline_equation(data = results,
                                aes(x = observed, y = imputed, label = ..rr.label..))   +
  facet_wrap(~variable, scales = "free") +
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"))

ggsave(filename = here::here("figures", "models", "spaMM_univariate_without_Random",
                             "Predictions_VS_observations.jpg"),
       width = 20, height = 10)



## Residuals distributions -> Normal ? zero-centered ?
distribution_plot(results, longer = F, index_values = c("variable", "res"), 
                  xlabel = "Residuals", ylabel = "frequency" )
ggsave(filename = here::here("figures", "models", "spaMM_univariate_without_Random",
                             "Residuals_distributions.jpg"),
       width = 20, height = 10)



## QQ-plot of residuals -> Normality of residuals ?
jpeg(filename = here::here("figures", "models", "spaMM_univariate_without_Random", "residuals_QQplots.jpg"),
     units = "cm", width = 40, height = 30, res = 400 )
par(mfrow = c(5,5))
for(model in models_uni){
  # model = models_uni[[1]]
  simulationOutput <- DHARMa::simulateResiduals(fittedModel = model, plot = F)
  DHARMa::plotQQunif(simulationOutput) 
  title(main = model[["predictor"]][[2]] , line = -0.5, cex.main = 1.5)
}
dev.off()



## Residuals according to Y values -> Homoscedasticity of residuals ?
plot_interaction(results, var_facet_wrap = "variable", X_values = "imputed",
                 Y_values = "res", xlabel = "Predicted contributions",
                 ylabel = "Residuals") + geom_hline(yintercept = 0)
ggsave(filename = here::here("figures", "models", "spaMM_univariate_without_Random",
                             "Residuals_VS_predicted_contributions.jpg"),
       width = 20, height = 10)


## Residuals according to covariates -> Patterns ?
colnames(covariates_final)
cov = "latitude"
res_cov <- results |> 
  dplyr::left_join( tibble::rownames_to_column(covariates_final, "survey_id"))
plot_interaction(res_cov, var_facet_wrap = "variable", X_values = cov,
                 Y_values = "res", xlabel = cov, ylabel = "Residuals")




## Simulated residuals vs predicted
jpeg(filename = here::here("figures", "models", "spaMM_univariate_without_Random", "residualsVSpredicted.jpg"),
     units = "cm", width = 40, height = 25, res = 400 )
par(mfrow = c(5,5))
for(model in models_uni){
  # model = models_uni[[1]]
  simulationOutput <- DHARMa::simulateResiduals(fittedModel = model, plot = F)
  DHARMa::plotResiduals(simulationOutput) 
  title(main = model[["predictor"]][[2]] , line = 0.5, cex.main = 1.5)
}
dev.off()




##---------------------- spaMM univar Central Indo-pacific only ----------------------

# Analyse each models

results <- do.call(rbind, pbmcapply::pbmclapply(models_uni_central_indo_pacific, 
                                                mc.cores = 10,
                                                function(model){
  # model <- models_uni_central_indo_pacific[[1]]
  preds <- predict(model, newdata = covariates_final|> 
                     dplyr::filter(realm == "Central Indo-Pacific"))
  
  #observe residuals
  residuals <- data.frame(variable =  as.character(model[["predictor"]][[2]]),
                          observed = model[["data"]][[model[["predictor"]][[2]]]],
                          imputed = preds,
                          res = spaMM::residuals.HLfit(model) ) |>   #residuals from spaMM package
    tibble::rownames_to_column("survey_id")
  
  residuals
}))


## Predictions distributions -> same patterns as observations ?
distribution_plot(results, longer = F, index_values = c("variable", "imputed"), 
                  xlabel = "Predictions", ylabel = "frequency" )
ggsave(filename = here::here("figures", "models", "spaMM_univariate_central_indo_pacific",
                             "Distributions_of_predicted_contributions.jpg"),
       width = 20, height = 10)



## Predictions vs observations -> model quality ?

# trace(ggpubr:::.stat_lm, edit = TRUE) #change signif() in line 9 to have more digits
ggplot(results)+
  geom_point(aes(x = observed, y = imputed, fill = variable),
             color = "grey40", alpha = 0.2, shape = 21) +
  hrbrthemes::theme_ipsum() +
  xlab("Observed contributions") + ylab("imputed")+
  geom_abline(slope = 1) + 
  # annotate("text", x = 2.2, y = 1.8, label = "1:1")+
  ggpubr::stat_regline_equation(data = results,
                                aes(x = observed, y = imputed, label = after_stat(rr.label)))   +
  facet_wrap(~variable, scales = "free") +
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"))

ggsave(filename = here::here("figures", "models", "spaMM_univariate_central_indo_pacific",
                             "Predictions_VS_observations.jpg"),
       width = 20, height = 10)



## Residuals distributions -> Normal ? zero-centered ?
distribution_plot(results, longer = F, index_values = c("variable", "res"), 
                  xlabel = "Residuals", ylabel = "frequency" )
ggsave(filename = here::here("figures", "models", "spaMM_univariate_central_indo_pacific",
                             "Residuals_distributions.jpg"),
       width = 20, height = 10)



## QQ-plot of residuals -> Normality of residuals ?
jpeg(filename = here::here("figures", "models", "spaMM_univariate_central_indo_pacific",
                           "residuals_QQplots.jpg"),
     units = "cm", width = 40, height = 30, res = 400 )
par(mfrow = c(5,5))
for(model in models_uni){
  # model = models_uni[[1]]
  simulationOutput <- DHARMa::simulateResiduals(fittedModel = model, plot = F)
  DHARMa::plotQQunif(simulationOutput) 
  title(main = model[["predictor"]][[2]] , line = -0.5, cex.main = 1.5)
}
dev.off()



## Residuals according to Y values -> Homoscedasticity of residuals ?
plot_interaction(results, var_facet_wrap = "variable", X_values = "imputed",
                 Y_values = "res", xlabel = "Predicted contributions",
                 ylabel = "Residuals") + geom_hline(yintercept = 0)
ggsave(filename = here::here("figures", "models", "spaMM_univariate_central_indo_pacific",
                             "Residuals_VS_predicted_contributions.jpg"),
       width = 20, height = 10)


## Residuals according to covariates -> Patterns ?
colnames(covariates_final)
cov = "latitude"
res_cov <- results |> 
  dplyr::left_join( tibble::rownames_to_column(covariates_final, "survey_id"))
plot_interaction(res_cov, var_facet_wrap = "variable", X_values = cov,
                 Y_values = "res", xlabel = cov, ylabel = "Residuals")




## Simulated residuals vs predicted
jpeg(filename = here::here("figures", "models", "spaMM_univariate_central_indo_pacific", 
                           "residualsVSpredicted.jpg"),
     units = "cm", width = 40, height = 25, res = 400 )
par(mfrow = c(5,5))
for(model in models_uni){
  # model = models_uni[[1]]
  simulationOutput <- DHARMa::simulateResiduals(fittedModel = model, plot = F)
  DHARMa::plotResiduals(simulationOutput) 
  title(main = model[["predictor"]][[2]] , line = 0.5, cex.main = 1.5)
}
dev.off()





###############################################################################
##
##               #### RELATIVE IMPORTANCE OF RANDOM EFFECTS ####
##
###############################################################################"

# model <- models_uni[[1]]
random_effect <- spaMM::ranef(model)
fixed_effect <- spaMM::fixef(model)
spaMM::pseudoR2(model)

install.packages("MuMIn")
model <- models_uni[[1]]
MuMIn::r.squaredGLMM(model)
install.packages("glmm.hp")
glmm.hp::glmm.hp(model)

###############################################################################"
##
##                       #### VARIABLE IMPORTANCE ####
##
###############################################################################"

##------------- AIC... -------------

test <- spaMM::get_any_IC(model, short.names = F)
spaMM::extractAIC.HLfit(model)


##------------- spaMM univariate - variables importance -------------

all_coefficients <- do.call(rbind, 
                            parallel::mclapply(models_uni,
                                               mc.cores = 10,
                                               function(model){
             
    contrib <- model[["call"]][["formula"]][[2]]
      
    #Coefficient plot
    coefficients = summary(model, details=list(p_value=TRUE))$beta_table  |> 
      as.data.frame() |>
      tibble::rownames_to_column("term") |>
      janitor::clean_names() |>
      dplyr::mutate(significance = ifelse(p_value >= 0.05, "Not significant p-value","Significant p-value")) |>
      # dplyr::mutate(t_value_signif = ifelse(t_value > 1.96 | t_value < - 1.96, "Significant t-value","Not significant t-value")) |>
      dplyr::filter(term != "(Intercept)") |> 
      dplyr::mutate(contribution = as.character(contrib))
    

    
    #Plot
    library(ggplot2)
    coefficients |>
        ggplot(aes(x=reorder(term, - estimate), y=estimate)) +
        geom_errorbar(aes(ymin=estimate-cond_se,ymax=estimate+cond_se),width=0.2)+
        geom_point(size = 3,aes(color = significance
                                # , shape = t_value_signif
                                )) +
        harrypotter::scale_color_hp_d(option = "Ravenclaw") +
        labs(title = paste0( contrib, 
                           " : Regression estimates of fixed effects (and 95% CI)"),
             color = "Significance of p-value (P < 0.05)",
             # shape = "Significance of t-value (t < -1.96 or > 1.96)",
             x = "", y = "") +
        coord_flip()+
        theme_minimal(base_size = 15) +
        theme(legend.position = "top",legend.box="vertical", legend.margin=margin()) +
        theme(plot.caption = element_text(hjust = 0))
    
    ggsave(filename = here::here("figures", "models","variable_importance",
                                 paste0("spaMM_univariate_", contrib, ".jpg")),
           width = 15, height = 10)
    
    coefficients
}))


# Classification of covariates
covariates <- unique(all_coefficients$term)
envir <- c(grep("median", covariates, value = T))
habitat <- c(grep("500m", covariates, value = T),
             "depth", "algae", "coral", "Sand", "seagrass", "microalgal_mats",
             "other_sessile_invert", "Rock", "coralline_algae", "coral_rubble")
human <- setdiff(covariates, c(envir, habitat))

# Classification of contributions
cont_list <- unique(all_coefficients$contribution)
NP <- c("available_biomass", "selenium", "zinc", "omega_3" , "calcium",  "iron",                       
        "vitamin_A", "available_biomass_turnover", "NP_score")
NN <- setdiff(cont_list, NP)

#Resume data
coeff_plot <- all_coefficients |> 
  dplyr::mutate(cov_class = ifelse(term %in% envir, "environmental",
                                   ifelse(term %in% habitat, "habitat", "human"))) |> 
  dplyr::mutate(contrib_class = ifelse(contribution %in% NP, "Nature-for-People",
                                    "Nature-for-Nature"))

# order covariates
df_ordered <- dplyr::filter(coeff_plot, contrib_class ==  "Nature-for-Nature") |> 
  dplyr::group_by(term) |> 
  dplyr::summarise(mean_estimate = median(estimate)) |> 
  dplyr::arrange(mean_estimate) |> 
  dplyr::pull(term)


## Boxplot of estimates
ggplot(coeff_plot) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_boxplot(aes(y=factor(term, levels = df_ordered), 
                   x = estimate, fill = cov_class),
               alpha = 0.7) +
  hrbrthemes::theme_ipsum() +
  # scale_fill_manual(values = c("#9467bd", "#ff7f0e", "#2ca02c"))+
  harrypotter::scale_fill_hp_d(option = "Ravenclaw") +
  facet_wrap(~ contrib_class) +
  xlab("Regression coefficient estimates") +
  ylab("Predictors") +
  theme( axis.title.x = element_text(size = 10),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         legend.position = "bottom" 
  )

ggsave(filename = here::here("figures", "models","variable_importance",
                             paste0("spaMM_univariate_ALL_contrib_boxplot.jpg")),
       width = 15, height = 10)




ggplot(coeff_plot) +
  stat_summary(aes(y = factor(term, levels = df_ordered), x = estimate, color = cov_class),
               fun.data = function(x) {
                 mean_val <- mean(x)
                 q1 <- quantile(x, 0.25)
                 q3 <- quantile(x, 0.75)
                 data.frame(y = mean_val, ymin = q1, ymax = q3)
               }, geom = "errorbar", width = 0.2) +
  stat_summary(aes(y = factor(term, levels = df_ordered), x = estimate, color = cov_class),
               fun = mean, geom = "point", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  hrbrthemes::theme_ipsum() +
  harrypotter::scale_color_hp_d(option = "Ravenclaw") +
  facet_wrap(~ contrib_class) +
  xlab("Regression coefficient estimates") +  ylab("Predictors") +
  theme( axis.title.x = element_text(size = 10),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         legend.position = "bottom" 
  )

ggsave(filename = here::here("figures", "models","variable_importance",
                             paste0("spaMM_univariate_ALL_contrib_geom_point.jpg")),
       width = 15, height = 10)


## Heat map of estimates
palette <- rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788"))
coef_matrix <- dplyr::select(coeff_plot, term, estimate, contribution) |> 
  tidyr::pivot_wider(names_from = "contribution", values_from = "estimate") |> 
  tibble::column_to_rownames("term") |> 
  as.matrix()

jpeg(filename = here::here("figures", "models", "variable_importance", "spaMM_univariate_heatmap_estimates.jpg"),
     units = "cm", width = 25, height = 15, res = 400 )
stats::heatmap(coef_matrix,  Colv=T,
               hclustfun=function(x) hclust(x, method="ward.D2"), 
               scale='none', col=palette, cexCol=0.6)
dev.off()


ggplot(coeff_plot)+
  aes(x = contribution, y = term, fill = estimate) +
  geom_tile() +
  scale_fill_gradient( low = "#FEE090", high = "#542788") +  
  theme_minimal() +
  labs(x = "", y = "Term", fill = "Contribution")




##------------- spaMM univar - WITHOUT spatial -------------

all_coefficients <- do.call(rbind, 
                            parallel::mclapply(models_uni_without_spatial,
                                               mc.cores = 10,
                                               function(model){
   # model = models_uni_without_spatial[[1]]
   contrib <- model[["call"]][["formula"]][[2]]
   
   #Coefficient plot
   coefficients = summary(model, details=list(p_value=TRUE))$beta_table  |> 
     as.data.frame() |>
     tibble::rownames_to_column("term") |>
     janitor::clean_names() |>
     dplyr::mutate(significance = ifelse(p_value >= 0.05, "Not significant p-value","Significant p-value")) |>
     # dplyr::mutate(t_value_signif = ifelse(t_value > 1.96 | t_value < - 1.96, "Significant t-value","Not significant t-value")) |>
     dplyr::filter(term != "(Intercept)") |> 
     dplyr::mutate(contribution = as.character(contrib))
   
   
   
   #Plot
   library(ggplot2)
   coefficients |>
     ggplot(aes(x=reorder(term, - estimate), y=estimate)) +
     geom_errorbar(aes(ymin=estimate-cond_se,ymax=estimate+cond_se),width=0.2)+
     geom_point(size = 3,aes(color = significance
                             # , shape = t_value_signif
     )) +
     harrypotter::scale_color_hp_d(option = "Ravenclaw") +
     labs(title = paste0( contrib, 
                          " : Regression estimates of fixed effects (and 95% CI)"),
          color = "Significance of p-value (P < 0.05)",
          # shape = "Significance of t-value (t < -1.96 or > 1.96)",
          x = "", y = "") +
     coord_flip()+
     theme_minimal(base_size = 15) +
     theme(legend.position = "top",legend.box="vertical", legend.margin=margin()) +
     theme(plot.caption = element_text(hjust = 0))
   
    ggsave(filename = here::here("figures", "models","spaMM_univar_without_spatial",
                                "variable_importance",
                                paste0("spaMM_univariate_", contrib, ".jpg")),
          width = 15, height = 10)
   
   coefficients
 }))


# Classification of covariates
covariates <- unique(all_coefficients$term)
envir <- c(grep("median", covariates, value = T))
habitat <- c(grep("500m", covariates, value = T),
             "depth", "algae", "coral", "Sand", "seagrass", "microalgal_mats",
             "other_sessile_invert", "Rock", "coralline_algae", "coral_rubble")
human <- setdiff(covariates, c(envir, habitat))

# Classification of contributions
cont_list <- unique(all_coefficients$contribution)
NP <- c("available_biomass", "selenium", "zinc", "omega_3" , "calcium",  "iron",                       
        "vitamin_A", "available_biomass_turnover", "NP_score")
NN <- setdiff(cont_list, NP)

#Resume data
coeff_plot <- all_coefficients |> 
  dplyr::mutate(cov_class = ifelse(term %in% envir, "environmental",
                                   ifelse(term %in% habitat, "habitat", "human"))) |> 
  dplyr::mutate(contrib_class = ifelse(contribution %in% NP, "Nature-for-People",
                                       "Nature-for-Nature"))

# order covariates
df_ordered <- dplyr::filter(coeff_plot, contrib_class ==  "Nature-for-Nature") |> 
  dplyr::group_by(term) |> 
  dplyr::summarise(mean_estimate = median(estimate)) |> 
  dplyr::arrange(mean_estimate) |> 
  dplyr::pull(term)


## Boxplot of estimates
ggplot(coeff_plot) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_boxplot(aes(y=factor(term, levels = df_ordered), 
                   x = estimate, fill = cov_class),
               alpha = 0.7) +
  hrbrthemes::theme_ipsum() +
  # scale_fill_manual(values = c("#9467bd", "#ff7f0e", "#2ca02c"))+
  harrypotter::scale_fill_hp_d(option = "Ravenclaw") +
  facet_wrap(~ contrib_class) +
  xlab("Regression coefficient estimates") +
  ylab("Predictors") +
  theme( axis.title.x = element_text(size = 10),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         legend.position = "bottom" 
  )

ggsave(filename = here::here("figures", "models","spaMM_univar_without_spatial",
                             "variable_importance",
                             paste0("spaMM_univariate_ALL_contrib_boxplot.jpg")),
       width = 15, height = 10)




ggplot(coeff_plot) +
  stat_summary(aes(y = factor(term, levels = df_ordered), x = estimate, color = cov_class),
               fun.data = function(x) {
                 mean_val <- mean(x)
                 q1 <- quantile(x, 0.25)
                 q3 <- quantile(x, 0.75)
                 data.frame(y = mean_val, ymin = q1, ymax = q3)
               }, geom = "errorbar", width = 0.2) +
  stat_summary(aes(y = factor(term, levels = df_ordered), x = estimate, color = cov_class),
               fun = mean, geom = "point", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  hrbrthemes::theme_ipsum() +
  harrypotter::scale_color_hp_d(option = "Ravenclaw") +
  facet_wrap(~ contrib_class) +
  xlab("Regression coefficient estimates") +  ylab("Predictors") +
  theme( axis.title.x = element_text(size = 10),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         legend.position = "bottom" 
  )

ggsave(filename = here::here("figures", "models","spaMM_univar_without_spatial",
                             "variable_importance",
                             paste0("spaMM_univariate_ALL_contrib_geom_point.jpg")),
       width = 15, height = 10)


## Heat map of estimates
palette <- rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788"))
coef_matrix <- dplyr::select(coeff_plot, term, estimate, contribution) |> 
  tidyr::pivot_wider(names_from = "contribution", values_from = "estimate") |> 
  tibble::column_to_rownames("term") |> 
  as.matrix()

jpeg(filename = here::here("figures", "models","spaMM_univar_without_spatial",
                           "variable_importance", "spaMM_univariate_heatmap_estimates.jpg"),
     units = "cm", width = 25, height = 15, res = 400 )
stats::heatmap(coef_matrix,  Colv=T,
               hclustfun=function(x) hclust(x, method="ward.D2"), 
               scale='none', col=palette, cexCol=0.6)
dev.off()


ggplot(coeff_plot)+
  aes(x = contribution, y = term, fill = estimate) +
  geom_tile() +
  scale_fill_gradient( low = "#FEE090", high = "#542788") +  
  theme_minimal() +
  labs(x = "", y = "Term", fill = "Contribution")





##------------- spaMM univariate CENTRAL INDO-PACIFIC - variables importance -------------

all_coefficients <- do.call(rbind, 
                            parallel::mclapply(models_uni_central_indo_pacific,
                                               mc.cores = 10,
                                               function(model){
  
  contrib <- model[["call"]][["formula"]][[2]]
  
  #Coefficient plot
  coefficients = summary(model, details=list(p_value=TRUE))$beta_table  |> 
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    janitor::clean_names() |>
    dplyr::mutate(significance = ifelse(p_value >= 0.05, "Not significant p-value","Significant p-value")) |>
    # dplyr::mutate(t_value_signif = ifelse(t_value > 1.96 | t_value < - 1.96, "Significant t-value","Not significant t-value")) |>
    dplyr::filter(term != "(Intercept)") |> 
    dplyr::mutate(contribution = as.character(contrib))
  
  
  
  #Plot
  library(ggplot2)
  coefficients |>
    ggplot(aes(x=reorder(term, - estimate), y=estimate)) +
    geom_errorbar(aes(ymin=estimate-cond_se,ymax=estimate+cond_se),width=0.2)+
    geom_point(size = 3,aes(color = significance
                            # , shape = t_value_signif
    )) +
    harrypotter::scale_color_hp_d(option = "Ravenclaw") +
    labs(title = paste0( contrib, 
                         " : Regression estimates of fixed effects (and 95% CI)"),
         color = "Significance of p-value (P < 0.05)",
         # shape = "Significance of t-value (t < -1.96 or > 1.96)",
         x = "", y = "") +
    coord_flip()+
    theme_minimal(base_size = 15) +
    theme(legend.position = "top",legend.box="vertical", legend.margin=margin()) +
    theme(plot.caption = element_text(hjust = 0))
  
  ggsave(filename = here::here("figures", "models","variable_importance_central_indo_pacific",
                               paste0("spaMM_univariate_", contrib, ".jpg")),
         width = 15, height = 10)
  
  coefficients
}))


# Classification of covariates
covariates <- unique(all_coefficients$term)
envir <- c(grep("median", covariates, value = T))
habitat <- c(grep("500m", covariates, value = T),
             "depth", "algae", "coral", "Sand", "seagrass", "microalgal_mats",
             "other_sessile_invert", "Rock", "coralline_algae", "coral_rubble")
human <- setdiff(covariates, c(envir, habitat))

# Classification of contributions
cont_list <- unique(all_coefficients$contribution)
NP <- c("available_biomass", "selenium", "zinc", "omega_3" , "calcium",  "iron",                       
        "vitamin_A", "available_biomass_turnover", "NP_score")
NN <- setdiff(cont_list, NP)

#Resume data
coeff_plot <- all_coefficients |> 
  dplyr::mutate(cov_class = ifelse(term %in% envir, "environmental",
                                   ifelse(term %in% habitat, "habitat", "human"))) |> 
  dplyr::mutate(contrib_class = ifelse(contribution %in% NP, "Nature-for-People",
                                       "Nature-for-Nature"))

# order covariates
df_ordered <- dplyr::filter(coeff_plot, contrib_class ==  "Nature-for-Nature") |> 
  dplyr::group_by(term) |> 
  dplyr::summarise(mean_estimate = median(estimate)) |> 
  dplyr::arrange(mean_estimate) |> 
  dplyr::pull(term)


## Boxplot of estimates
ggplot(coeff_plot) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_boxplot(aes(y=factor(term, levels = df_ordered), 
                   x = estimate, fill = cov_class),
               alpha = 0.7) +
  xlim(-1,1)+
  hrbrthemes::theme_ipsum() +
  # scale_fill_manual(values = c("#9467bd", "#ff7f0e", "#2ca02c"))+
  harrypotter::scale_fill_hp_d(option = "Ravenclaw") +
  facet_wrap(~ contrib_class) +
  xlab("Regression coefficient estimates") +
  ylab("Predictors") +
  theme( axis.title.x = element_text(size = 10),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         legend.position = "bottom" 
  )

ggsave(filename = here::here("figures", "models","variable_importance_central_indo_pacific",
                             paste0("spaMM_univariate_ALL_contrib_boxplot.jpg")),
       width = 15, height = 10)




ggplot(coeff_plot) +
  stat_summary(aes(y = factor(term, levels = df_ordered), x = estimate, color = cov_class),
               fun.data = function(x) {
                 mean_val <- mean(x)
                 q1 <- quantile(x, 0.25)
                 q3 <- quantile(x, 0.75)
                 data.frame(y = mean_val, ymin = q1, ymax = q3)
               }, geom = "errorbar", width = 0.2) +
  stat_summary(aes(y = factor(term, levels = df_ordered), x = estimate, color = cov_class),
               fun = mean, geom = "point", size = 3) +
  xlim(-0.7,0.7)+
  geom_vline(xintercept = 0, linetype = "dashed")+
  hrbrthemes::theme_ipsum() +
  harrypotter::scale_color_hp_d(option = "Ravenclaw") +
  facet_wrap(~ contrib_class) +
  xlab("Regression coefficient estimates") +  ylab("Predictors") +
  theme( axis.title.x = element_text(size = 10),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         legend.position = "bottom" 
  )

ggsave(filename = here::here("figures", "models","variable_importance_central_indo_pacific",
                             paste0("spaMM_univariate_ALL_contrib_geom_point.jpg")),
       width = 15, height = 10)


## Heat map of estimates
palette <- rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788"))
coef_matrix <- dplyr::select(coeff_plot, term, estimate, contribution) |> 
  tidyr::pivot_wider(names_from = "contribution", values_from = "estimate") |> 
  tibble::column_to_rownames("term") |> 
  as.matrix()

jpeg(filename = here::here("figures", "models", "variable_importance_central_indo_pacific",
                           "spaMM_univariate_heatmap_estimates.jpg"),
     units = "cm", width = 25, height = 15, res = 400 )
stats::heatmap(coef_matrix,  Colv=T,
               hclustfun=function(x) hclust(x, method="ward.D2"), 
               scale='none', col=palette, cexCol=0.6)
dev.off()


ggplot(coeff_plot)+
  aes(x = contribution, y = term, fill = estimate) +
  geom_tile() +
  scale_fill_gradient( low = "#FEE090", high = "#542788") +  
  theme_minimal() +
  labs(x = "", y = "Term", fill = "Contribution")


