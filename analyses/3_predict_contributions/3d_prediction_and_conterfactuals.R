###############################################################################"
##
##
##
## 3d_prediction_and_conterfactuals.R
##
## 28/06/2024
##
## Ulysse Flandrin
##
###############################################################################"

##---------------------------- cleaning memory ---------------------------------
rm(list=ls())

##-----------------------------Loading packages---------------------------------
pkgs <- c("here", "Hmsc", "coda", "ggmcmc", "jsonify", "FactoMineR", "factoextra",
          "ggplot2")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))
library(ggplot2)
library(patchwork)

##------------------------------- load data ------------------------------------
load(here::here("data/derived_data/3_all_contributions_to_predict.Rdata"))
load(here::here("data/derived_data/3_all_covariates_to_predict.Rdata"))
# response =  observations_final#[sample(1:nrow(observations_final),1000),] ####################### reduce data
# covariates = covariates_final[rownames(response),]

# Crossvalidation result
load( here::here("outputs/models/hmsc/cross_validation/predictions_crossval_5folds.Rdata"))


# PATHS
save_init <- here::here("outputs/models/hmsc/init_multi/")
save_out <- here::here("outputs/models/hmsc/out_multi/")
localDir <- here::here("outputs/models/hmsc/multivariate")

# plot functions
source(here::here("R","evaluation_prediction_model.R"))
source(here::here("R/HMSC_function.R"))


##----------------------------- Predictive power in crossvalidation ------------------------------

### Predictive power ###
obs <- NULL
preds <- NULL

for(i in 1:length(predictions_cv)){
  temp <- predictions_cv[[i]]
  obs <- rbind(obs, temp[[1]])
  preds <- rbind(preds, temp[[2]])
}

obs_long <- obs |> 
  tibble::rownames_to_column("survey_id") |> 
  tidyr::pivot_longer(cols = -survey_id, names_to = "contribution", values_to = "observation")

predictive_power <- preds |> 
  tibble::rownames_to_column("survey_id") |> 
  tidyr::pivot_longer(cols = -survey_id, 
                      names_to = "contribution", 
                      values_to = "prediction") |> 
  dplyr::left_join(obs_long)


ggplot(predictive_power)+
  geom_point(aes(x = observation, y = prediction, fill = contribution),
             color = "grey40", alpha = 0.2, shape = 21) +
  hrbrthemes::theme_ipsum() +
  xlab("Observed contributions") + ylab("imputed")+
  geom_abline(slope = 1) + 
  ggpubr::stat_regline_equation(data = predictive_power,
                                aes(x = observation, y = prediction,
                                    label = after_stat(rr.label)))   +
  facet_wrap(~contribution, scales = "free") +
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"))

ggsave(filename = paste0("figures/models/hmsc/Crossvalidation_predictive_power.jpg"),
       width = 15, height = 8)


predictive_power_summary <- predictive_power |> 
  dplyr::group_by(contribution) |> 
  dplyr::summarise(r_squared = summary(lm(prediction ~ observation))[["r.squared"]])

png("figures/models/hmsc/predictive_power_CV_5_folds.png",
    width = 10, height = 10, units = "cm", res = 300)
hist(predictive_power_summary$r_squared, xlim = c(0,1), 
     main=paste0("Mean = ", round(mean(predictive_power_summary$r_squared),2)))
dev.off()

##----------------------------- import hmsc output ------------------------------

## List all files in the directory and choose the model
list_files <- list.files(save_out) 
list_files
file_name <- gsub("output_", "", list_files[9]) #choose the wanted file


## Import initial object

#Model design
# load(file = file.path(localDir, paste0( file_name))) #old models
load(file = file.path(localDir, paste0("model_fit_", file_name, ".Rdata")))

#Initial hmsc object
init_obj_rds <- readRDS(paste0(save_init, paste0("init_", file_name)))
init_obj <- jsonify::from_json(init_obj_rds)

nSamples = init_obj[["samples"]]
thin = init_obj[["thin"]]
nChains = init_obj[["nChains"]]
transient = init_obj[["transient"]]


## Import posterior probability
file_rds <- readRDS(file = paste0(save_out, paste0("output_", file_name)))[[1]]
importFromHPC <- jsonify::from_json(file_rds)
postList <- importFromHPC[1:nChains]

##Export and merge chains
model_fit_mcmc <- Hmsc::importPosteriorFromHPC(model_fit,
                                               postList, 
                                               nSamples, 
                                               thin, 
                                               transient)



## Estimates for each chains
mpost <- Hmsc::convertToCodaObject(model_fit_mcmc)

postBeta <- Hmsc::getPostEstimate(model_fit_mcmc, parName = "Beta")


##----------------------------- Predict contributions -----------------------------------
predictions <- run_hmsc_prediction(X_data = covariates_final,
                                   model_fit_mcmc = model_fit_mcmc)[[1]]


distribution_plot(predictions, longer = T, cols_plot = colnames(predictions))
#Quite good distributions compared to observations, except for IUCN species.
ggsave(width=15, height= 10, here::here("figures/models/hmsc/conterfactuals", 
                             "predicted_contributions_distribution.jpg"))

residuals <- predictions - observations_final
distribution_plot(residuals, longer = T,cols_plot = colnames(residuals))
# ~zero centered.


##### TO DO: test spatial autocorrelations of residuals. ######



##----------------------------- Counterfactual scenarios -----------------------------------

#### Change initial conditions ###
summary(covariates_final)

#(1) Change effectiveness only: from "out" to "high protection"
X_new_mpa <- covariates_final
new_surveys_mpa <- rownames(X_new_mpa |> dplyr::filter(effectiveness == "out"))
X_new_mpa[new_surveys_mpa, "effectiveness"] <- as.factor("High")

#(2) Change fishing pressure only
X_new_vessels <- covariates_final
new_surveys_vessels <- rownames(
  X_new_vessels[X_new_vessels$n_fishing_vessels != min(X_new_vessels$n_fishing_vessels),]
)
X_new_vessels[new_surveys_vessels, "n_fishing_vessels"] <- min(X_new_vessels$n_fishing_vessels)

#(3) Change fishing pressure and MPA = real protection
X_new_protected <- covariates_final
new_surveys_protected <- unique(c(new_surveys_mpa, new_surveys_vessels))
X_new_protected[new_surveys_protected, "effectiveness"] <- as.factor("High")
X_new_protected[new_surveys_protected, "n_fishing_vessels"] <- min(X_new_protected$n_fishing_vessels)

#(4) Change human pressure: no gravity and high neartt
X_new_far <- covariates_final
new_surveys_far <- rownames(
  X_new_far[X_new_far$gravtot2 != min(X_new_far$gravtot2) |
              X_new_far$neartt != max(X_new_far$neartt),])
X_new_far[new_surveys_far, "gravtot2"] <- min(X_new_far$gravtot2)
X_new_far[new_surveys_far, "neartt"] <- max(X_new_far$neartt)

#(5) Pristine sites
X_new_pristine <- X_new_protected
new_surveys_pristine <- rownames(X_new_pristine)
X_new_pristine[new_surveys_pristine, "gravtot2"] <- min(X_new_pristine$gravtot2)
X_new_pristine[new_surveys_pristine, "neartt"] <- max(X_new_pristine$neartt)



#### Run model with new conditions ###
X_new <- X_new_mpa ; new_survey <- new_surveys_mpa
X_new <- X_new_vessels ; new_survey <- new_surveys_vessels
X_new <- X_new_protected ; new_survey <- new_surveys_protected
X_new <- X_new_far ; new_survey <- new_surveys_far
X_new <- X_new_pristine ; new_survey <- new_surveys_pristine


new_predictions <- run_hmsc_prediction(X_data = covariates_final,
                                       model_fit_mcmc = model_fit_mcmc,
                                       X_new_data = X_new,
                                       new_surveys = new_survey)

preds <- new_predictions[["predictions"]]
conterfactual <- new_predictions[["new_scenario"]]
effective_change <- new_predictions[["effective_change"]]


distribution_plot(conterfactual, longer = T,
                  cols_plot = colnames(conterfactual))


## Check differences ##
distrib_boxplot <- function(data, x, y, fill, hline = 0, title = NULL){
  ggplot(data) +
    aes_string(x= x, y= y, fill = fill)+
    geom_boxplot() +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "text", 
                 aes(label = paste(round(after_stat(y), 2))), 
                 position = position_dodge(width = 0.75), vjust = 0.5, size = 3,
                 color = "grey50") + 
    geom_hline(yintercept = hline, linetype = "dashed", color = "coral3") +
    xlab("") + ylab("Contributions change in counterfactual scenarios") +
    labs(title=title)+
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 10))+
    coord_flip()
}

distribution_plot(effective_change, longer = F, cols_plot = colnames(effective_change))

effective_change <- effective_change  |> 
  dplyr::mutate(index = reorder(index, values, FUN = median))

distrib_boxplot(effective_change, x = "index", y = "values", fill = "index",
                title = "No vessels")
ggsave(width=8, height= 8, here::here("figures/models/hmsc/conterfactuals", 
                                        "changes_in_new_vessels.jpg"))


#obs heterogeneity by country
data_1_contrib <- effective_change |>  
  dplyr::filter(index == "iucn_species_richness") |> 
  dplyr::mutate(country = reorder(country, values, FUN = median))

distrib_boxplot(data_1_contrib, x = "country", y = "values", fill = "country",
                hline = mean(data_1_contrib$values))




##----------------------------- Plot changes -----------------------------------
selected_country <- covariates_final |> 
  dplyr::count(country) |> 
  dplyr::filter(n > 20) |> 
  dplyr::pull(country)
  
X_new <- X_new_mpa ; new_survey <- new_surveys_mpa
X_new <- X_new_vessels ; new_survey <- new_surveys_vessels
X_new <- X_new_protected ; new_survey <- new_surveys_protected
X_new <- X_new_far ; new_survey <- new_surveys_far
X_new <- X_new_pristine ; new_survey <- new_surveys_pristine


new_predictions <- run_hmsc_prediction(X_data = covariates_final,
                                       model_fit_mcmc = model_fit_mcmc,
                                       X_new_data = X_new,
                                       new_surveys = new_survey)


#Original PCA
pca <- FactoMineR::PCA(new_predictions[[1]], scale.unit = T, graph=F, ncp=15)

factoextra::fviz_screeplot(pca, ncp=15)

pca_plot <- factoextra::fviz_pca_biplot(pca,
                            title="",
                            label = "var",
                            labelsize = 4, 
                            geom=c("point"), 
                            pointshape=21,
                            stroke=0, pointsize=2,
                            alpha.ind = 0.7,
                            alpha.var = 0.7,
                            fill.ind = "grey",    
                            repel = TRUE)

# Calculate barycenters (centroids) for each country
barycenters <- aggregate(pca$ind$coord[, 1:2], 
                         by = list(X_new$country), FUN = mean) |> 
  dplyr::filter(Group.1 %in% selected_country)
colnames(barycenters) <- c("country", "PC1", "PC2")

# Add barycenters to the plot
pca_plot_barycenter <- pca_plot + 
  geom_point(data = barycenters, aes(x = PC1, y = PC2, fill = country), 
             shape = 23, size = 4, color = "black")

### Obs new scenario ###
# (1) Project the new predictions in the same space
projected_points <- predict(pca, newdata = new_predictions[[2]])

# (2) calculate new barycenters
barycenters_new <- aggregate(projected_points$coord[, 1:2], 
                         by = list(X_new$country), FUN = mean)
colnames(barycenters_new) <- c("country", "PC1", "PC2")

barycenters_movement <- merge(barycenters, barycenters_new, 
                              by = "country", suffixes = c("_old", "_new"))

pca_plot_barycenter + 
  # geom_point(data = barycenters_new, aes(x = PC1, y = PC2, fill = country), 
  #                                  shape = 23, size = 3, color = "black")
  geom_segment(data = barycenters_movement, 
               aes(x = PC1_old, y = PC2_old, 
                   xend = PC1_new, yend = PC2_new,
                   color = country), 
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1)+
  labs(title = "No human and total protection")

ggsave(width=15, height= 8, here::here("figures/models/hmsc/conterfactuals", 
                                      "barycenters_movment_pristine.jpg"))



## General barycenter:
barycenter_tot <- apply(pca$ind$coord[, 1:2], 2, mean) # 0, 0 by construction
barycenters_new <- t(as.data.frame(apply(projected_points$coord[, 1:2], 2, mean)))
pca_plot + geom_segment(data = barycenters_new, 
                        aes(x = 0, y = 0, 
                            xend = Dim.1, yend = Dim.2), 
                        arrow = arrow(length = unit(0.3, "cm")),
                        size = 1)+
  labs(title = "No human and total protection")

