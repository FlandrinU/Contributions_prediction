################################################################################
##
## 
##
## 3a_preping_data_to_predict.R
##
## 11/03/2024
##
## Ulysse Flandrin
##
###############################################################################"

## cleaning memory
rm(list=ls())

##------------------- Loading datasets-------------------
# Surveys metadata
load(file = here::here("data", "raw_data", "environmental_covariates",
                       "all_covariates_benthos_inferred_tropical_surveys.Rdata"))

# Contributions matrix
load(file = here::here("outputs", "2_all_contributions.Rdata"))


##------------------- Clean observations and covariates-------------------
# FILTER OBSERVATIONS
colnames(contributions)

observations <- contributions |> 
  dplyr::select(-N_recycling, -P_recycling) |> 
  tidyr::drop_na()


# FILTER COVARIATES
colnames(all_covariates_benthos_inferred)

df <- all_covariates_benthos_inferred[,-c(1:8,11,13,14,15,16,176)]
pca <- FactoMineR::PCA(df, scale = T, graph=F, ncp=30)
factoextra::fviz_screeplot(pca, ncp = 20)
factoextra::fviz_pca_var(pca, col.var = "contrib",repel = TRUE)
factoextra::fviz_pca_var(pca, col.var = "contrib", axes=c(3,4), repel = TRUE)

covariates <- all_covariates_benthos_inferred |> 
  dplyr::select(survey_id, latitude, longitude, depth, year,
                
                median_7days_chl, median_5year_chl, median_7days_degree_heating_week,
                median_5year_degree_heating_week, median_7days_nppv, median_5year_nppv,
                median_7days_o2, median_5year_o2, median_5year_ph, median_7days_so_mean,
                median_5year_so_mean, median_7days_analysed_sst, median_5year_analysed_sst,
                
                coral_algae_500m, Rock_500m, Sand_500m, Seagrass_500m, Plateau_500m,
                algae, coral, Sand, seagrass, microalgal_mats, "other sessile invert",
                Rock, "coralline algae", "coral rubble",
                
                control_of_corruption, gdp, gravtot2, hdi, marine_ecosystem_dependency,
                effectiveness, natural_ressource_rent, neartt, ngo, no_violence,
                voice, n_fishing_vessels, 
  ) 

funbiogeo::fb_plot_species_traits_completeness(dplyr::rename(covariates, species = survey_id))
           
#PREPING COVARIATES                                       
covariates_final <- covariates |> 
  dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
                                              "out" = 0,
                                              "Low" = 1,
                                              "Medium" = 2,
                                              "High" = 3)) |>  #effectiveness into quantitative values
  tidyr::drop_na() |> #30% of loss notably due to the Allen Atlas
  dplyr::filter(survey_id %in% rownames(observations)) |> 
  tibble::column_to_rownames("survey_id") |> 
  dplyr::mutate(across(-c(longitude, latitude), scale)) |> #SCALE ALL COVARIATES
  dplyr::mutate(across(everything(), as.numeric))
  # dplyr::mutate(across(-c(longitude, latitude, effectiveness), scale)) |> #SCALE ALL COVARIATES
  # dplyr::mutate(across(-effectiveness, as.numeric))

#Distribution of covariates
library(ggplot2)
ggplot(data=tidyr::pivot_longer(covariates_final,
                                cols = -all_of(c("longitude", "latitude", "effectiveness")),
                                names_to = "index", values_to = "values"))+
  aes(x=values, group=index, fill=index) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, color = "grey40", fill ="white") +
  geom_density(aes(fill = index), alpha = 0.2) +
  hrbrthemes::theme_ipsum() +
  facet_wrap(~index, scales = "free") +
  theme(legend.position="none",panel.spacing = unit(0.1, "lines"),
        axis.ticks.x=element_blank())


#COMMON SURVEYS FOR COVARIATES AND OBSERVATIONS
observations_final <- observations[rownames(covariates_final),]
dim(observations_final) #4427 SURVEYS, 21 CONTRIBUTIONS



##------------------- Rough cross-validation -------------------
n_CV = 10

datasets <- lapply(c(1:n_CV), FUN = function(i){
  sample <- sample.int(nrow(observations_final), round(0.2*nrow(observations_final),0))
  
  train <- observations_final[-sample,]
  test <- observations_final[sample,]
  
  list(train, test)
})



##------------------- save datasets -------------------
save(covariates_final, file = here::here("data", "derived_data", "3_all_covariates_to_predict.Rdata"))
save(observations_final, file = here::here("data", "derived_data", "3_all_contributions_to_predict.Rdata"))

save(datasets, file = here::here("data", "derived_data", "3_datasets_for_predict_CV_80_20.Rdata"))
