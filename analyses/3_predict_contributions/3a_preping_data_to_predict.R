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

##-----------------Loading packages-------------------
# Sys.unsetenv("GITHUB_PAT")
# remotes::install_github("FRBCesab/funbiogeo")
# pkgs <- c("here", "dplyr", "funbiogeo", "performance", "rsample", "car")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)

##------------------- Loading datasets-------------------
# Surveys metadata
load(file = here::here("data", "raw_data", "environmental_covariates",
                       "all_covariates_benthos_inferred_tropical_surveys.Rdata"))

# Contributions matrix
load(file = here::here("outputs", "2_all_contributions.Rdata"))
load(file = here::here("outputs", "2_all_contributions_with_synthetic_score.Rdata"))
load(file = here::here("outputs", "2_metadata_backtransformation_contrib.Rdata") )


# Contribution matrix at the site scale
load(file = here::here("outputs", "2_contributions_site&date.Rdata"))

#Coastline
coast <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')

## Load functions ##
source(here::here("R","evaluation_prediction_model.R"))

###############################################################################"
##
##                            #### SURVEY SCALE ####
##
###############################################################################"

##------------------- Clean observations and covariates-------------------

#### FILTER OBSERVATIONS ####

# colnames(contributions)
colnames(contributions_with_synthetic_score)
funbiogeo::fb_plot_species_traits_completeness(
  contributions_with_synthetic_score |> tibble::rownames_to_column("species"))
  


# observations <- contributions |> 
observations <- contributions_with_synthetic_score |> 
  dplyr::select(-N_recycling, -P_recycling, -elasmobranch_richness,) |> 
  tidyr::drop_na() |> 
  dplyr::mutate(across(-c(iucn_species_richness), scale)) |> 
  dplyr::mutate(across(everything(), as.numeric))

means <- contributions_with_synthetic_score[rownames(observations),] |> 
  dplyr::summarise(across(-iucn_species_richness, mean)) |> 
  tidyr::pivot_longer(cols = everything(), names_to = "contribution", values_to = "mean_survey")

sd <- contributions_with_synthetic_score[rownames(observations),] |> 
  dplyr::summarise(across(-iucn_species_richness, ~sd(.)))|> 
  tidyr::pivot_longer(cols = everything(), names_to = "contribution", values_to = "sd_survey")

#Metadata back-scaling
contributions_transformation <- contributions_transformation |> 
  dplyr::mutate(scaled = ifelse(contribution != "iucn_species_richness", T, F)) |> 
  dplyr::left_join(means)|> 
  dplyr::left_join(sd)
  

#### FILTER COVARIATES ####
colnames(all_covariates_benthos_inferred)
colnames(all_covariates_benthos_inferred) <- gsub(" ", "_", colnames(all_covariates_benthos_inferred))

## Too many variables -> look at correlations:
sapply(all_covariates_benthos_inferred, class) 
metadata_col <- c("survey_id", "country", "area", "ecoregion", "realm", "location",
                  "site_code", "site_name", "latitude", "longitude", "survey_date",
                  "program", "hour", "effectiveness")

data_to_filter <- all_covariates_benthos_inferred[,-which(names(
  all_covariates_benthos_inferred) %in% metadata_col)] |> 
  tidyr::drop_na()

## Correlation matrix
env <- as.data.frame(data_to_filter) |> 
  dplyr::select(
    grep("q95", colnames(data_to_filter)),
    grep("q05", colnames(data_to_filter)),
    # grep("mean", colnames(data_to_filter)), #SAME AS MEDIAN
    grep("median", colnames(data_to_filter)),
    # grep("max", colnames(data_to_filter)), #SAME AS 95% QUANTILE
    # grep("min", colnames(data_to_filter)),
    -q05_5year_degree_heating_week
    )
M <- cor(env)

png(filename = here::here("figures/models/covariates","corr_matrix_env_covariates.png"), 
    width= 40, height = 30, units = "cm", res = 1000)
  corrplot::corrplot(M, order = 'AOE', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r', tl.cex = 0.5)
dev.off() 

# More observations of DHW:
dhw <- as.data.frame(data_to_filter) |> 
  dplyr::select(grep("degree_heating_week", colnames(data_to_filter)))
dhw <- dhw[, sapply(dhw, function(x) var(x, na.rm = TRUE) != 0)]
corrplot::corrplot(cor(dhw), order = 'AOE', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r', tl.cex = 0.5)

distribution_plot(dhw, cols_plot = colnames(dhw)) #NEED TO BE LOG-TRANSFORMED

# More observation of 7 days variables
days <- as.data.frame(data_to_filter) |> 
  dplyr::select(
    grep("q95_7", colnames(data_to_filter)),
    grep("q05_7", colnames(data_to_filter)),
    grep("median_7", colnames(data_to_filter)),
  )
corrplot::corrplot(cor(days), order = 'AOE', tl.pos = 'tp',
                   tl.srt = 60, cl.pos = 'r', tl.cex = 0.5) #MEDIAN AND QUANTILES ARE THE SAME

#Other variables
others <- as.data.frame(data_to_filter) |> 
  dplyr::select(-grep("q95", colnames(data_to_filter)),
                -grep("q05", colnames(data_to_filter)),
                -grep("mean", colnames(data_to_filter)),
                -grep("median", colnames(data_to_filter)),
                -grep("max", colnames(data_to_filter)),
                -grep("min", colnames(data_to_filter)))
M <- cor(others)
png(filename = here::here("figures/models/covariates","corr_matrix_others_covariates.png"), 
    width= 40, height = 30, units = "cm", res = 1000)
corrplot::corrplot(M, order = 'AOE', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r', tl.cex = 0.7)
dev.off() 

#Importance of "country"
ggplot(all_covariates_benthos_inferred) +
geom_boxplot(aes(x = reorder(country, median_5year_o2), y = median_5year_o2)) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.text = element_text(size = 6))


ggplot(all_covariates_benthos_inferred) +
  geom_boxplot(aes(x = reorder(country, median_5year_analysed_sst), y = median_5year_analysed_sst)) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.text = element_text(size = 6))
ggsave(filename = paste0("figures/models/covariates/country_VS_SST.jpg"),
       width = 8, height = 6)


### REMOVE OBVIOUS CORRELATIONS 
data_to_filter <- as.data.frame(data_to_filter) |> 
  dplyr::select(-grep("mean_", colnames(data_to_filter)),
                -grep("max_", colnames(data_to_filter)),
                -grep("min_", colnames(data_to_filter)),
                -q05_5year_degree_heating_week
                ) |> 
  as.matrix()


### FIRST SELECTION OF VARIABLES
cov <- colnames(data_to_filter)
cov <- cov[order(cov)]

metadata_to_select <- c("survey_id", "site_code", "latitude", "longitude",
                        "country", "ecoregion", "realm",
                        "depth", "year",
                        "effectiveness")

cov_to_select <- c(#Environment
                     cov[grep("median_5year", cov)],
                     cov[grep("degree_heating_week", cov)],
                     cov[grep("q05_5year", cov)],
                     cov[grep("q95_5year", cov)],
                     
                     cov[grep("median_7days", cov)],
                     
                   #Habitat
                     cov[grep("500m", cov)],
                     "algae", "coral", "Sand", 
                     "seagrass", "microalgal_mats", "other_sessile_invert",
                     "Rock", "coralline_algae", "coral_rubble",
                     
                   #Human
                     "control_of_corruption", "gdp", "gravtot2", "hdi", 
                     "marine_ecosystem_dependency", 
                     "natural_ressource_rent", "neartt", "ngo",
                     "no_violence", "voice", "n_fishing_vessels"
  )

cov_to_select

### OBSERVE COVARIATES THAT ARE STILL CORRELATED
cor_matrix <- cor(as.data.frame(data_to_filter) |> 
                    dplyr::select(all_of(cov_to_select)))
png(filename = here::here("figures/models/covariates","corr_matrix_first_selection_covariates.png"),
    width= 40, height = 30, units = "cm", res = 1000)
corrplot::corrplot(cor_matrix, order = 'AOE', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r', tl.cex = 0.7)
dev.off()

# Covariates network -> THRESHOLD r = 0.7 of correlation.
pairwise_corr <- cor_matrix[upper.tri(cor_matrix)]
high_correlation_indices <- which((cor_matrix > 0.7 & cor_matrix < 1) |
                                    cor_matrix < -0.7 , arr.ind = TRUE)
correlated_pairs <- cbind(rownames(cor_matrix)[high_correlation_indices[, "row"]],
                          colnames(cor_matrix)[high_correlation_indices[, "col"]])
network_graph <- igraph::graph_from_edgelist(correlated_pairs, directed = FALSE)

png("figures/models/covariates/Selection of covariates_higly_correlated_cov_(0.7).png",
    width = 30, height = 30, units = "cm", res = 300)
plot(network_graph, 
     layout = igraph::layout_with_kk(network_graph),
     edge.arrow.size = 1, 
     vertex.label.cex = 0.8)
dev.off()


### REFINE THE SELECTION
cov_to_select2 <- 
  cov_to_select[-which(cov_to_select %in% 
                         c("median_5year_nppv",
                           "median_7days_nppv", # correlated to chlorophyll
                           "median_7days_o2",
                           "median_5year_o2",
                           "median_7days_so_mean",
                           # "median_5year_ph", #correlated with sst, low ecological meaning
                           "q05_5year_nppv",
                           "q05_5year_chl",
                           "q95_5year_chl",
                           "q95_5year_o2",
                           "q05_5year_o2",
                           "q95_5year_ph",
                           "q05_5year_ph", # correlated at the median_5year_pH
                           "q95_5year_nppv",
                           "q95_5year_so_mean",
                           "q05_5year_so_mean",
                           "q95_5year_analysed_sst",
                           "q05_5year_analysed_sst",
                           "q05_7days_degree_heating_week",
                           "q95_7days_degree_heating_week",
                           "q05_1year_degree_heating_week", #no influence in models
                           "q95_1year_degree_heating_week", # low contribution to variance explained, correlated at 0.5 with q95_5year_DHW
                           "median_7days_analysed_sst", #correlated with median_5year_sst, low ecological relevance
                           
                           "Outer_Reef_Flat_500m", # correlated to Rock_500m
                           "Shallow_Lagoon_500m", #correlated to Sand_500m
                           "seagrass", #no influence in models
                           "microalgal_mats", #no influence in models
                           "algae", #highly negatively correlated with coral, and less influence in variance explained
                           
                           "no_violence", #correlated to hdi
                           "control_of_corruption",
                           "voice",
                           # "hdi",  #country level covariates
                           # "marine_ecosystem_dependency",  #country level covariates
                           "ngo" #country level covariates
                           # "natural_ressource_rent" #country level covariates
                         ))]


selection2 <- as.data.frame(data_to_filter) |> dplyr::select(all_of(cov_to_select2))
png("figures/models/covariates/Selected_covariates_correlation.png",
    width = 35, height = 25, units = "cm", res = 300)
corrplot::corrplot(cor(selection2), order = 'AOE', tl.pos = 'tp', tl.srt = 60, 
                   cl.pos = 'r', tl.cex = 0.7, tl.offset = 0.5)
corrplot::corrplot(cor(selection2), add = TRUE, type = 'upper', method = 'number',
                   order = 'AOE', insig = 'p-value', diag = FALSE, tl.pos = 'n', 
                   cl.pos = 'n', number.digits = 1, number.cex = 0.7)
dev.off()

pca <- FactoMineR::PCA(selection2, scale = T, graph=F, ncp=30)
factoextra::fviz_pca_var(pca, col.var = "contrib",repel = TRUE)

### VIF ANALYSIS
cov_selected <- all_covariates_benthos_inferred |> 
  dplyr::select(all_of(c("survey_id",cov_to_select2))) |> 
  tibble::column_to_rownames("survey_id") 

model <- lm(observations[rownames(cov_selected), "available_biomass"] ~ ., data = cov_selected)

vif_values <- car::vif(model)
print(vif_values[order(vif_values)]) #No covariates with VIF > 5.

### DISTRIBUTION OF COVARIATES
distribution_plot(selection2, longer = T,  cols_plot = cov_to_select2 )
ggsave( width=15, height= 10,
        filename = here::here("figures/models/covariates", "3_raw_covariates_distribution.jpg"))

cov_to_log_transformed <- 
  c("Back_Reef_Slope_500m", "Deep_Lagoon_500m", "gdp", "gravtot2",
    "Inner_Reef_Flat_500m", 
    "marine_ecosystem_dependency", "natural_ressource_rent",
    "median_1year_degree_heating_week", "median_5year_chl",
    "median_5year_degree_heating_week", "median_7days_chl",
    "median_7days_degree_heating_week", "median_7days_degree_heating_week",
    "Microalgal_Mats_500m", "n_fishing_vessels", 
    "neartt", "Patch_Reefs_500m", "Plateau_500m" ,
    "Reef_Crest_500m","Reef_Slope_500m",  
    "Rock_500m", "Rubble_500m", "Sand_500m", "Seagrass_500m",
    "Sheltered_Reef_Slope_500m", "Terrestrial_Reef_Flat_500m"
  )


#### FINAL SELECTION OF COVARIATES ####
metadata_to_select
cov_to_select2

covariates <- all_covariates_benthos_inferred |> 
  dplyr::select(all_of(c(metadata_to_select,cov_to_select2))) |> 
  dplyr::mutate(across(.cols = all_of(cov_to_log_transformed),
                       .fns = ~ .x +1 , .names = "{.col}")) |>
  dplyr::mutate(across(.cols = all_of(cov_to_log_transformed),
                       .fns = log10 , .names = "{.col}"))          # log(x+1) to avoid -Inf values


funbiogeo::fb_plot_number_species_by_trait(dplyr::rename(covariates, species = survey_id))
funbiogeo::fb_plot_species_traits_completeness(dplyr::rename(covariates, species = survey_id))
ggsave(plot = last_plot(), width=15, height= 10,
       filename = here::here("figures/models/covariates", "covariates_completedness.jpg"))

#FINAL COVARIATES                                       
covariates_final <- covariates |> 
  # dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
  #                                             "out" = 0,
  #                                             "Low" = 1,
  #                                             "Medium" = 2,
  #                                             "High" = 3)) |>  
  
  #Change the order of levels of effectiveness for the GLM
  dplyr::mutate(effectiveness = factor(effectiveness, 
                                       levels = c("out", "Low", "Medium", "High"))) |> 
  # tidyr::drop_na() |> #30% of loss notably due to the Allen Atlas
  dplyr::filter(survey_id %in% rownames(observations)) |> 
  tibble::column_to_rownames("survey_id") |>
  dplyr::mutate(across(-c(longitude, latitude,
                          site_code,
                          effectiveness,
                          country,
                          realm,
                          ecoregion), scale)) |> #SCALE ALL COVARIATES
  dplyr::mutate(across(-c(site_code, effectiveness, country, realm, ecoregion), as.numeric))



#Distribution of covariates
distribution_plot(covariates_final, longer = T, cols_plot = cov_to_select2 )
ggsave(plot = last_plot(), width=15, height= 10,
       filename = here::here("figures/models/covariates", 
                             "covariates_distribution_log_transformed_and_scaled.jpg"))

# ## Check covariates with low heterogeneity
# low_hetero <- c("median_7days_degree_heating_week", "Patch_Reefs_500m",
#                 "Seagrass_500m", "Terrestrial_Reef_Flat_500m")
# summary(covariates_final[,low_hetero])
# distribution_plot(covariates_final, longer = T, cols_plot = low_hetero )+
#   ylim(c(0,5))


#### COMMON SURVEYS FOR COVARIATES AND OBSERVATIONS ####
covariates_final_without_Allen <- covariates_final[rownames(observations),] |> 
  dplyr::select(-grep("_500m", colnames(covariates_final))) |> 
  tidyr::drop_na() 

covariates_final <- covariates_final[rownames(observations),] |> 
  tidyr::drop_na() #30% of loss notably due to the Allen Atlas


observations_final <- observations[rownames(covariates_final),] |> 
  dplyr::select(-NN_score, -NP_score)
dim(observations_final) #4423 SURVEYS, 22 CONTRIBUTIONS

observations_final_aggregated_score <- observations[rownames(covariates_final),] |> 
  dplyr::select(NN_score, NP_score)

observations_final_without_Allen <- observations[rownames(covariates_final_without_Allen),] |> 
  dplyr::select(-NN_score, -NP_score)
dim(observations_final_without_Allen) #5170 SURVEYS, 22 CONTRIBUTIONS


#Distribution of observations
distribution_plot(observations_final, longer = T,
                  cols_plot = colnames(observations_final)) #OK: log transformed and scaled values
ggsave(plot = last_plot(), width=15, height= 10,
       filename = here::here("figures/models/covariates", 
                             "contributions_distribution_log_transformed_and_scaled.jpg"))


pca <- FactoMineR::PCA(cbind(observations_final,observations_final_aggregated_score),
                       scale.unit = T, graph=F, ncp=15,
                       quanti.sup = c("NN_score", "NP_score"))
factoextra::fviz_pca_biplot(pca, repel = TRUE, geom="point", pointshape=21,
                            stroke=0, pointsize=3, alpha.ind = 0.7, 
                            fill.ind = "grey", col.quanti.sup = "firebrick")


##------------------- save datasets -------------------
save(covariates_final, file = here::here("data", "derived_data", "3_all_covariates_to_predict.Rdata"))
save(observations_final, file = here::here("data", "derived_data", "3_all_contributions_to_predict.Rdata"))
save(observations_final_aggregated_score, 
     file = here::here("data", "derived_data", "3_NN_NP_scores_to_predict.Rdata"))

# save(datasets, file = here::here("data", "derived_data", "3_datasets_for_predict_CV_66_33.Rdata"))

save(covariates_final_without_Allen, file = here::here("data", "derived_data", "3_covariates_without_Allen_to_predict.Rdata"))
save(observations_final_without_Allen, file = here::here("data", "derived_data", "3_contributions_without_Allen_to_predict.Rdata"))



# ##------------------- Chose statistical distribution -------------------
# cov <- covariates_final |> dplyr::select(-longitude, -latitude, -country, -year)
# 
# distrib <- list()
# distri_df <- data.frame(row.names = colnames(observations_final))
# 
# for( contrib in colnames(observations_final)){
#   fmla <- as.formula( paste(contrib, "~ ", paste(c("1", colnames(cov)) , 
#           collapse = "+")))
#   fit <- glm(formula = fmla, data = cbind(cov, observations_final))
#   
#   distribution <- performance::check_distribution(fit)
#   which_max <- distribution$p_Response == max(distribution$p_Response)
#   
#   distrib[[contrib]] <- distribution
#   
#   distri_df[contrib, "distribution"] <- paste(distribution$Distribution[which_max],
#                                               collapse = "_OR_")
#   distri_df[ contrib, "proba"] <- unique(distribution$p_Response[which_max])
# }


###############################################################################"
##
##                            #### SITE SCALE ####
##   (WE MEAN CONTRIBUTIONS AND COVARIATES IN THE SAME SITE AT THE SAME DATE)
##
###############################################################################"
##------------------- Clean observations -------------------
colnames(contributions_sites_date)

observations_site <- contributions_sites_date |>
  dplyr::select(-N_recycling, -P_recycling,) |> 
  tidyr::drop_na() |> 
  dplyr::mutate(across(-c(iucn_species_richness), scale)) |> 
  dplyr::mutate(across(everything(), as.numeric))

means <- contributions_sites_date[rownames(observations_site),] |> 
  dplyr::summarise(across(-iucn_species_richness, mean)) |> 
  tidyr::pivot_longer(cols = everything(), names_to = "contribution", values_to = "mean_site")

sd <- contributions_sites_date[rownames(observations_site),] |> 
  dplyr::summarise(across(-iucn_species_richness, ~sd(.)))|> 
  tidyr::pivot_longer(cols = everything(), names_to = "contribution", values_to = "sd_site")

#Metadata back-scaling
contributions_transformation <- contributions_transformation |> 
  dplyr::left_join(means)|> 
  dplyr::left_join(sd)


##------------------- Mean covariates -------------------
# Mean the covariates at the site scale, for a given date: all surveys in the
# same place, observed at the same date are merged
#### FINAL SELECTION OF COVARIATES ####
metadata_to_select
cov_to_select2

covariates_site <- all_covariates_benthos_inferred |> 
  dplyr::select(all_of(c(metadata_to_select, "survey_date",cov_to_select2))) |> 
  #Agregate at the site scale
  dplyr::select(-survey_id) |> 
  dplyr::group_by(site_code, latitude, longitude, country, ecoregion, realm, 
                  survey_date, year, effectiveness) |> 
  dplyr::summarise(across(.cols = everything(),
                          .fns = ~mean(., na.rm = TRUE), .names = "{.col}")) |> 
  dplyr::mutate(across(.cols = all_of(cov_to_select2),
                       .fns = ~ifelse(is.nan(.), NA, .), .names = "{.col}")) |> 
  dplyr::mutate(id = paste0(site_code, "_", survey_date)) |> 
  dplyr::ungroup() |> 
  dplyr::select(id, everything(), -survey_date) |> 
  #log-transform
  dplyr::mutate(across(.cols = all_of(cov_to_log_transformed),
                       .fns = ~ .x +1 , .names = "{.col}")) |>
  dplyr::mutate(across(.cols = all_of(cov_to_log_transformed),
                       .fns = log10 , .names = "{.col}"))      # log(x+1) to avoid -Inf values


funbiogeo::fb_plot_number_species_by_trait(dplyr::rename(covariates_site, species = id))
funbiogeo::fb_plot_species_traits_completeness(dplyr::rename(covariates_site, species = id))
# ggsave(plot = last_plot(), width=15, height= 10,
#        filename = here::here("figures/models/covariates", "covariates_completedness.jpg"))


#FINAL COVARIATES                                       
covariates_site_final <- covariates_site |> 
  # dplyr::mutate(effectiveness = dplyr::recode(effectiveness,
  #                                             "out" = 0,
  #                                             "Low" = 1,
  #                                             "Medium" = 2,
  #                                             "High" = 3)) |>  
  
  #Change the order of levels of effectiveness for the GLM
  dplyr::mutate(effectiveness = factor(effectiveness, 
                                       levels = c("out", "Low", "Medium", "High"))) |> 
  # tidyr::drop_na() |> #30% of loss notably due to the Allen Atlas
  dplyr::filter(id %in% rownames(observations_site)) |> 
  tibble::column_to_rownames("id") |>
  dplyr::mutate(across(-c(longitude, latitude,
                          site_code,
                          effectiveness,
                          country,
                          realm,
                          ecoregion), scale)) |> #SCALE ALL COVARIATES
  dplyr::mutate(across(-c(site_code, effectiveness, country, realm, ecoregion), as.numeric))



#Distribution of covariates
distribution_plot(covariates_site_final, longer = T, cols_plot = cov_to_select2 )
# ggsave(plot = last_plot(), width=15, height= 10,
#        filename = here::here("figures/models/covariates", 
#                              "covariates_distribution_log_transformed_and_scaled.jpg"))





#### COMMON SITES FOR COVARIATES AND OBSERVATIONS ####
covariates_site_final <- covariates_site_final[rownames(observations_site),] |> 
  tidyr::drop_na() # loss notably due to the Allen Atlas


observations_site_final <- observations_site[rownames(covariates_site_final),] |> 
  dplyr::select(-NN_score, -NP_score)
dim(observations_site_final) #2488 SITE/DATE, 22 CONTRIBUTIONS


#Distribution of observations
distribution_plot(observations_site_final, longer = T,
                  cols_plot = colnames(observations_site_final)) #OK: log transformed and scaled values
# ggsave(plot = last_plot(), width=15, height= 10,
#        filename = here::here("figures/models/covariates", 
#                              "contributions_distribution_log_transformed_and_scaled.jpg"))


pca <- FactoMineR::PCA(observations_site_final,
                       scale.unit = T, graph=F, ncp=15
                       )
factoextra::fviz_pca_biplot(pca, repel = TRUE, geom="point", pointshape=21,
                            stroke=0, pointsize=3, alpha.ind = 0.7, 
                            fill.ind = "grey")



##------------------- Map and caracterize RLS sites -------------------

table(covariates_site_final$effectiveness)
length(unique(covariates_site_final$country))

plot_mpa <-function(covariates_site_final, xlim=c(-180,180), ylim = c(-36, 31),
                    legend_pos = "none", jitter = 0.2){
  ggplot(covariates_site_final) +
    geom_sf(data = coast, color = "grey30", fill = "lightgrey",
            aes(size=0.1)) +
    
    geom_point(size = 2, na.rm = T,
               position = position_jitter(width =jitter, height =jitter),
               alpha = 0.8,
               colour = "black",
               stroke=0.1,
               shape = 21,
               aes(x = longitude, y = latitude,
                   fill=effectiveness)) +
    
    coord_sf(xlim, ylim , expand = FALSE) +
    guides(alpha = "none", size = "none", colour = "none") +
    # scale_shape_manual(values=c(21,24,23))+
    
    theme_minimal()+
    labs(title = "",
         x="", y= "") +
    theme(legend.position = legend_pos,
          plot.title = element_text(size=15, face="bold"),
          legend.text = element_text(size=13),
          legend.title = element_text(size=15),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
}

mpa <-  plot_mpa(covariates_site_final , xlim=c(-180,180), ylim = c(-36, 31),
                 legend_pos = "none", jitter = 1)+
  geom_rect(aes(xmin = 110, xmax = 160, ymin = -32, ymax = -7), color = "black", fill= "transparent")+
  geom_rect(aes(xmin = -95, xmax = -67, ymin = -3, ymax = 18), color = "black", fill= "transparent")

gold_coast_mpa<- plot_mpa( covariates_site_final,ylim = c(-32, -7),
                           xlim= c(110,160), legend_pos = "right", jitter = 0.5)

caraib_mpa <- plot_mpa( covariates_site_final, ylim = c(-3, 18),
                        xlim= c(-95,-67), legend_pos = "none", jitter = 0.5)

ggpubr::ggarrange(mpa, # First row with world map
                  ggpubr::ggarrange(caraib_mpa, gold_coast_mpa,  
                                    ncol = 2, labels = c("B", "C"), widths = c(1, 1.3)), # Second row with zooms
                  nrow = 2, labels = "A") 
ggsave(plot = last_plot(), width = 14, height = 7,
       filename = here::here("figures","RLS sites with protection and zoom.jpg"))
       


##------------------- save datasets -------------------
save(covariates_site_final, file = here::here("data", "derived_data", "3_sites_covariates_to_predict.Rdata"))
save(observations_site_final, file = here::here("data", "derived_data", "3_sites_contributions_to_predict.Rdata"))

# load(file = here::here("data", "derived_data", "3_sites_covariates_to_predict.Rdata"))

save(contributions_transformation,
     file = here::here("outputs", "3_metadata_backtransformation_contrib.Rdata") )
