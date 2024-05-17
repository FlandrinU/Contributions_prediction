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
load(file = here::here("outputs", "2_all_contributions_with_synthetic_score.Rdata"))

# Contribution matrix at the site scale
load(file = here::here("outputs", "2_contributions_site&date.Rdata"))

## Load functions ##
source(here::here("R","evaluation_prediction_model.R"))

###############################################################################"
##
##                            #### SURVEY SCALE ####
##
###############################################################################"

##------------------- Clean observations and covariates-------------------

# FILTER OBSERVATIONS

# colnames(contributions)
colnames(contributions_with_synthetic_score)

# observations <- contributions |> 
observations <- contributions_with_synthetic_score |> 
  dplyr::select(-N_recycling, -P_recycling, -elasmobranch_richness) |> 
  tidyr::drop_na() |> 
  dplyr::mutate(across(-c(iucn_species_richness), scale)) |> 
  dplyr::mutate(across(everything(), as.numeric))




# FILTER COVARIATES
colnames(all_covariates_benthos_inferred)
colnames(all_covariates_benthos_inferred) <- gsub(" ", "_", colnames(all_covariates_benthos_inferred))

## Too many variables -> look at correlations:
sapply(all_covariates_benthos_inferred, class) 
metadata_col <- c("survey_id", "country", "area", "ecoregion", "realm", "location",
                  "site_code", "site_name", "latitude", "longitude", "survey_date",
                  "program", "hour", "effectiveness")

data_to_filter <- all_covariates_benthos_inferred[,-which(names(
  all_covariates_benthos_inferred) %in% metadata_col)] |> 
  tidyr::drop_na() |> 
  scale()

## Correlation matrix
cor_matrix <- cor(data_to_filter, method = "pearson")
pairwise_corr <- cor_matrix[upper.tri(cor_matrix)]
summary(pairwise_corr) 
high_correlation_indices <- which(cor_matrix > 0.7 & cor_matrix < 1, arr.ind = TRUE)
correlated_pairs <- cbind(rownames(cor_matrix)[high_correlation_indices[, "row"]],
                          colnames(cor_matrix)[high_correlation_indices[, "col"]])
unique(correlated_pairs[,1])

#Observe the highly correlated pairs:
network_graph <- igraph::graph_from_edgelist(correlated_pairs, directed = FALSE)
plot(network_graph, layout = igraph::layout_with_fr(network_graph),
     edge.arrow.size = 1, vertex.label.cex = 0.8)

### REMOVE OBVIOUS CORRELATIONS AND RE-RUN LINES 61 to 73
data_to_filter <- as.data.frame(data_to_filter) |> 
  dplyr::select(-grep("q95", colnames(data_to_filter)),
                -grep("q05", colnames(data_to_filter)),
                -grep("mean", colnames(data_to_filter)),
                -grep("max", colnames(data_to_filter)),
                -grep("min", colnames(data_to_filter))
  ) |> as.matrix()


#SELECT VARIABLES
colnames(data_to_filter)
pca <- FactoMineR::PCA(data_to_filter, scale = T, graph=F, ncp=30)
factoextra::fviz_pca_var(pca, col.var = "contrib",repel = TRUE)
factoextra::fviz_pca_var(pca, col.var = "contrib", axes=c(3,4), repel = TRUE)

covariates <- all_covariates_benthos_inferred |> 
  dplyr::select(survey_id, latitude, longitude, depth, year,
                country, realm, ecoregion, 
                
                median_7days_chl, median_5year_chl,
                median_7days_degree_heating_week, median_5year_degree_heating_week,
                median_7days_nppv, median_5year_nppv,
                median_7days_o2, median_5year_o2, 
                median_5year_ph,
                median_7days_so_mean, median_5year_so_mean, 
                median_7days_analysed_sst, median_5year_analysed_sst,
                
                coral_algae_500m, Rock_500m, Sand_500m, Seagrass_500m, Plateau_500m,
                algae, coral, Sand, seagrass, microalgal_mats, other_sessile_invert,
                Rock, coralline_algae, coral_rubble,
                
                control_of_corruption, gdp, gravtot2, hdi, marine_ecosystem_dependency,
                effectiveness, natural_ressource_rent, neartt, ngo, no_violence,
                voice, n_fishing_vessels, 
  ) 
# |> 
#   dplyr::select(-grep("7days", colnames(data_to_filter)), -grep("500m", colnames(data_to_filter))) # REDUCE DATA TO TEST MODELS######################################################"
  
  
  
funbiogeo::fb_plot_species_traits_completeness(dplyr::rename(covariates, species = survey_id))
           
#PREPING FINAL COVARIATES                                       
covariates_final <- covariates |> 
  # dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
  #                                             "out" = 0,
  #                                             "Low" = 1,
  #                                             "Medium" = 2,
  #                                             "High" = 3)) |>  
  
  #Change the order of levels of effectiveness for the GLM
  dplyr::mutate(effectiveness = factor(effectiveness, 
                                       levels = c("out", "Low", "Medium", "High"))) |> 
  tidyr::drop_na() |> #30% of loss notably due to the Allen Atlas
  dplyr::filter(survey_id %in% rownames(observations)) |> 
  tibble::column_to_rownames("survey_id") |> 
  # dplyr::mutate(across(-c(longitude, latitude), scale)) |> #SCALE ALL COVARIATES
  # dplyr::mutate(across(everything(), as.numeric))
  dplyr::mutate(across(-c(longitude, latitude,
                          effectiveness,
                          country,
                          realm,
                          ecoregion), scale)) |> #SCALE ALL COVARIATES
  dplyr::mutate(across(-c(effectiveness, country, realm, ecoregion), as.numeric))



#Distribution of covariates
distribution_plot(covariates_final, longer = T,
                  cols_not_plot = c("longitude", "latitude", "effectiveness",
                                    "country", "realm", "ecoregion") )
ggsave(plot = last_plot(), width=15, height= 10,
       filename = here::here("figures", "3_covariates_distribution_scaled.jpg"))

#COMMON SURVEYS FOR COVARIATES AND OBSERVATIONS
observations_final <- observations[rownames(covariates_final),]
dim(observations_final) #4420 SURVEYS, 21 CONTRIBUTIONS + 2 SCORES

#Distribution of observations
distribution_plot(observations_final, longer = T,
                  cols_plot = colnames(observations_final)) #OK: log transformed and scaled values
pca <- FactoMineR::PCA(observations_final, scale.unit = T, graph=F, ncp=15,
                       quanti.sup = c("NN_score", "NP_score"))
factoextra::fviz_pca_biplot(pca, repel = TRUE, geom="point", pointshape=21,
                            stroke=0, pointsize=3, alpha.ind = 0.7, 
                            fill.ind = "grey", col.quanti.sup = "firebrick")


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




##------------------- Chose statistical distribution -------------------
cov <- covariates_final |> dplyr::select(-longitude, -latitude, -country, -year)

distrib <- list()
distri_df <- data.frame(row.names = colnames(observations_final))

for( contrib in colnames(observations_final)){
  fmla <- as.formula( paste(contrib, "~ ", paste(c("1", colnames(cov)) , 
          collapse = "+")))
  fit <- glm(formula = fmla, data = cbind(cov, observations_final))
  
  distribution <- performance::check_distribution(fit)
  which_max <- distribution$p_Response == max(distribution$p_Response)
  
  distrib[[contrib]] <- distribution
  
  distri_df[contrib, "distribution"] <- paste(distribution$Distribution[which_max],
                                              collapse = "_OR_")
  distri_df[ contrib, "proba"] <- unique(distribution$p_Response[which_max])
}


###############################################################################"
##
##                            #### SITE SCALE ####
##   (WE MEAN CONTRIBUTIONS AND COVARIATES IN THE SAME SITE AT THE SAME DATE)
##
###############################################################################"
##------------------- Clean observations -------------------
# FILTER OBSERVATIONS

colnames(contributions_sites_date)

observations <- contributions_sites_date |> 
  dplyr::select(-N_recycling, -P_recycling) |> 
  tidyr::drop_na() |> 
  dplyr::mutate(across(everything(), scale)) |> 
  dplyr::mutate(across(everything(), as.numeric))


##------------------- Mean covariates -------------------
# Mean the covariates at the site scale, for a given date: all surveys in the
# same place, observed at the same date are merged
covariates_sites <- all_covariates_benthos_inferred |> 
  dplyr::select(-survey_id, -area, -location, -site_name, -program, -visibility, -hour) |> 
  dplyr::group_by(site_code, latitude, longitude, country, ecoregion, realm, 
                  survey_date, year, effectiveness) |> 
  dplyr::summarise(across(.cols = everything(),
                          .fns = ~mean(., na.rm = TRUE), .names = "{.col}")) |> 
  dplyr::mutate(survey_id = paste0(site_code, "_", survey_date)) |> 
  dplyr::ungroup()

nrow(covariates_sites) # 3405 mean_sites


# FILTER COVARIATES
colnames(covariates_sites)
colnames(covariates_sites) <- gsub(" ", "_", colnames(covariates_sites))

## Too many variables -> look at correlations:
sapply(covariates_sites, class) 
metadata_col <- c("survey_id", "country", "area", "ecoregion", "realm", "location",
                  "site_code", "site_name", "latitude", "longitude", "survey_date",
                  "program", "hour", "effectiveness")

data_to_filter <- covariates_sites[,-which(names(covariates_sites)
                                           %in% metadata_col)] |> 
  tidyr::drop_na() |> 
  scale()

## Correlation matrix
cor_matrix <- cor(data_to_filter, method = "pearson")
pairwise_corr <- cor_matrix[upper.tri(cor_matrix)]
summary(pairwise_corr) 
high_correlation_indices <- which(cor_matrix > 0.7 & cor_matrix < 1, arr.ind = TRUE)
correlated_pairs <- cbind(rownames(cor_matrix)[high_correlation_indices[, "row"]],
                          colnames(cor_matrix)[high_correlation_indices[, "col"]])
unique(correlated_pairs[,1])

#Observe the highly correlated pairs:
network_graph <- igraph::graph_from_edgelist(correlated_pairs, directed = FALSE)
plot(network_graph, layout = igraph::layout_with_fr(network_graph),
     edge.arrow.size = 1, vertex.label.cex = 0.8)

### REMOVE OBVIOUS CORRELATIONS AND RE-RUN LINES 228 to 241
data_to_filter <- as.data.frame(data_to_filter) |> 
  dplyr::select(-grep("q95", colnames(data_to_filter)),
                -grep("q05", colnames(data_to_filter)),
                -grep("mean", colnames(data_to_filter)),
                -grep("max", colnames(data_to_filter)),
                -grep("min", colnames(data_to_filter))
  ) |> as.matrix()


#SELECT VARIABLES
colnames(data_to_filter)
pca <- FactoMineR::PCA(data_to_filter, scale = T, graph=F, ncp=30)
factoextra::fviz_pca_var(pca, col.var = "contrib",repel = TRUE)
factoextra::fviz_pca_var(pca, col.var = "contrib", axes=c(3,4), repel = TRUE)

covariates <- covariates_sites |> 
  dplyr::select(survey_id, latitude, longitude, depth, year,
                country, ecoregion, realm,
                
                median_7days_chl, median_5year_chl,
                median_7days_degree_heating_week, median_5year_degree_heating_week,
                median_7days_nppv, median_5year_nppv,
                median_7days_o2, median_5year_o2, 
                median_5year_ph,
                median_7days_so_mean, median_5year_so_mean, 
                median_7days_analysed_sst, median_5year_analysed_sst,
                
                coral_algae_500m, Rock_500m, Sand_500m, Seagrass_500m, Plateau_500m,
                algae, coral, Sand, seagrass, microalgal_mats, other_sessile_invert,
                Rock, coralline_algae, coral_rubble,
                
                control_of_corruption, gdp, gravtot2, hdi, marine_ecosystem_dependency,
                effectiveness, natural_ressource_rent, neartt, ngo, no_violence,
                voice, n_fishing_vessels
  ) 


funbiogeo::fb_plot_species_traits_completeness(dplyr::rename(covariates, species = survey_id))


## SEE DISTRIBUTION
distribution_plot(covariates, longer = T,
                  cols_not_plot = c("longitude", "latitude", 
                                     "effectiveness", "country", "ecoregion",
                                     "realm", "survey_id") )

# to_log <- c("coral_algae_500m", "coral_rubble","coralline_algae" ,"gdp" ,"gravtot2",
#             "hdi", "median_5year_chl","median_5year_nppv", "median_7days_chl", 
#             "median_7days_degree_heating_week", "median_7days_nppv", "microalgal_mats",
#             "n_fishing_vessels", "neartt", "other_sessile_invert", "Plateau_500m",
#             "Seagrass_500m", "seagrass")



#PREPING FINAL COVARIATES                                       
covariates_final <- covariates |> 
  # dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
  #                                             "out" = 0,
  #                                             "Low" = 1,
  #                                             "Medium" = 2,
  #                                             "High" = 3)) |>  #effectiveness into quantitative values
 
  #Change the order of levels of effectiveness for the GLM
  dplyr::mutate(effectiveness = factor(effectiveness, 
                                       levels = c("out", "Low", "Medium", "High"))) |> 
  tidyr::drop_na() |> #30% of loss notably due to the Allen Atlas
  dplyr::filter(survey_id %in% rownames(observations)) |> 
  tibble::column_to_rownames("survey_id") |> 
  # # Log transformation
  # dplyr::mutate(across(.cols = all_of(to_log),
  #                      .fns = ~ .x +1 , .names = "{.col}")) |>
  # dplyr::mutate(across(.cols = all_of(to_log),
  #                      .fns = log10 , .names = "{.col}")) |> 
  # dplyr::mutate(across(-c(longitude, latitude), scale)) |> #SCALE ALL COVARIATES
  # dplyr::mutate(across(everything(), as.numeric))
  dplyr::mutate(across(-c(longitude, latitude, effectiveness, ecoregion,
                          country, realm), scale)) |> #SCALE ALL COVARIATES
  dplyr::mutate(across(-c(effectiveness, country, ecoregion, realm), as.numeric))


#Distribution of covariates
distribution_plot(covariates_final, longer = T,
                  cols_not_plot = c("longitude", "latitude", 
                                    "effectiveness", "country",
                                    "ecoregion", "realm") )



#COMMON SURVEYS FOR COVARIATES AND OBSERVATIONS
observations_final <- observations[rownames(covariates_final),]


## TRY FILTERS ON OBSERVATIONS
# -> FAKE COLUMNS TO CHECK MODELS
# observations_final$test_sum <- covariates_final$median_5year_analysed_sst + covariates_final$coral
# observations_final$test_prod <- covariates_final$median_5year_chl * covariates_final$gravtot2
# remove elasmobranchs (impossible to predict)
observations_final$elasmobranch_richness <- NULL

# # 1) remove well predict observations
# observations_final$mean_endemism <- NULL
# observations_final$functional_distinctiveness <- NULL
# observations_final$omega_3 <- NULL
# observations_final$actino_richness <- NULL
# observations_final$evolutionary_distinctiveness <- NULL
# 2) take only correlated contributions togethers:
# observations_final <- dplyr::select(observations_final, )


#Distribution of covariates
dim(observations_final) #2467 SITES, 20 CONTRIBUTIONS + 2 SCORES
distribution_plot(observations_final, longer = T,
                  cols_plot = colnames(observations_final) ) #OK: log transformed and scales values
pca <- FactoMineR::PCA(observations_final, scale.unit = T, graph=F, ncp=15,
                       quanti.sup = c("NN_score", "NP_score"))
factoextra::fviz_pca_biplot(pca, repel = TRUE, geom="point", pointshape=21,
                            stroke=0, pointsize=3, alpha.ind = 0.7, 
                            fill.ind = "grey", col.quanti.sup = "firebrick")





##------------------- Rough cross-validation -------------------
n_CV = 10

datasets <- lapply(c(1:n_CV), FUN = function(i){
  sample <- sample.int(nrow(observations_final), round(0.2*nrow(observations_final),0))
  
  train <- observations_final[-sample,]
  test <- observations_final[sample,]
  
  list(train, test)
})



##------------------- save datasets -------------------
save(covariates_final, file = here::here("data", "derived_data", "3_sites_covariates_to_predict.Rdata"))
save(observations_final, file = here::here("data", "derived_data", "3_sites_contributions_to_predict.Rdata"))

save(datasets, file = here::here("data", "derived_data", "3_sites_datasets_for_predict_CV_80_20.Rdata"))

