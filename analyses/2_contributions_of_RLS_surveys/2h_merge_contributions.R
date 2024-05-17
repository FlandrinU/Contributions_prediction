################################################################################
##
## 
##
## 2h_merge_contributions.R
##
## 26/02/2024
##
## Ulysse Flandrin
##
###############################################################################"

## cleaning memory
rm(list=ls())

##------------------- loading datasets-------------------
# Surveys metadata
load(file = here::here("data", "raw_data", "environmental_covariates",
                       "all_covariates_benthos_inferred_tropical_surveys.Rdata"))

# Biodiversity indices
load(file = here::here("outputs", "2b_biodiv_indices_surveys.Rdata" ))

# Biomass distribution
load(file = here::here("outputs", "2b_biomass_indices_surveys.Rdata" ))

# Biochemichal flows
load(file=here::here( "outputs", "2d_surveys_fluxes.Rdata") )

# Food web stability
load( file = here::here("outputs", "2e_food_web_indices_surveys.Rdata" ))

# Food intake
load(file = here::here("outputs", "2f_food_indices_surveys.Rdata" ))

#Cultural values


#Coastline
coast <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')

## Load functions ##
source(here::here("R","evaluation_prediction_model.R"))


##------------------- Merge all estimated contributions -------------------####
metadata <- all_covariates_benthos_inferred |> 
  dplyr::select(survey_id, country, ecoregion, realm, site_code, latitude,
                longitude, survey_date, depth, year, sst  = mean_5year_analysed_sst,
                effectiveness, coral
                )
var_metadata <- colnames(metadata)

contributions_surveys <- metadata |> 
  dplyr::left_join(biodiv_indices_surveys) |> 
  dplyr::left_join(biomass_indices_surveys) |> 
  dplyr::left_join(surveys_fluxes_final) |> 
  dplyr::left_join(food_web_indices) |> 
  dplyr::left_join(food_indices_surveys) |> 
  # dplyr::left_join(cultural_indices) |> 
  dplyr::select(all_of(var_metadata),
                actino_richness,
                elasmobranch_richness, 
                functional_distinctiveness,
                iucn_species_richness = iucn_species,
                mean_endemism, 
                evolutionary_distinctiveness = ED_mean,
                functional_entropy = funct_entropy,
                phylogenetic_entropy = phylo_entropy_Mean,
                herbivores_biomass = biom_lowTL,
                invertivores_biomass = biom_mediumTL, 
                piscivores_biomass = biom_highTL,
                total_biomass,
                N_recycling = recycling_N, 
                P_recycling = recycling_P,
                trophic_web_robustness = b_power_law,
                mean_trophic_level = troph_mTL,
                available_biomass, 
                selenium = Selenium_C, 
                zinc = Zinc_C, 
                omega_3 = Omega3_C, 
                calcium = Calcium_C, 
                iron = Iron_C,
                vitamin_A = VitaminA_C,
                available_biomass_turnover = productivity)

summary(contributions_surveys)
#299 RLS SITES WITH NO FISHES KEPT FOR THE ANALYSIS -> TO REMOVE
contributions_surveys <- dplyr::filter(contributions_surveys, !is.na(actino_richness))

# diversity remaining
dplyr::n_distinct(contributions_surveys$survey_id) # 6002 surveys
dplyr::n_distinct(contributions_surveys$site_code) # 1977 sites
dplyr::n_distinct(contributions_surveys$country) # 41 countries





###############################################################################"
##
##                            #### SURVEY SCALE ####
##
###############################################################################"

##------------------- Check distributions and transform values -------------------####
summary(contributions_surveys)

# distributions
distribution_plot(contributions_surveys, longer = T,
                  cols_not_plot = var_metadata)

ggsave(plot = last_plot(), width=15, height= 10,
       filename = here::here("figures", "2_contributions_distribution_raw.jpg"))


# high magnitude distribution
median <-  robustbase::colMedians(as.matrix(
  dplyr::select(contributions_surveys, actino_richness:available_biomass_turnover)),
  na.rm = T)
max <- apply(dplyr::select(contributions_surveys,
                           actino_richness:available_biomass_turnover),2,
             function(x) max(x, na.rm=T))
max/median
which(max/median > 5)

# IDENTIFICATION OF RIGHT-SKEWED DISTRIBUTIONS:
to_log <- c("herbivores_biomass", "invertivores_biomass", "piscivores_biomass",
            "total_biomass", "N_recycling" , "P_recycling", "available_biomass",
            "available_biomass_turnover",
            
            "vitamin_A", "omega_3", "iron", "calcium")

zero_values <- c("herbivores_biomass", "piscivores_biomass")


# Log transformation
contributions_surveys_log <- contributions_surveys |>
  dplyr::mutate(across(.cols = all_of(zero_values),
                       .fns = ~ .x +1 , .names = "{.col}")) |>
  dplyr::mutate(across(.cols = all_of(to_log),
                       .fns = log10 , .names = "{.col}")) # log(x+1) to avoid -Inf values

# CHECK DISTRIBUTIONS
summary(contributions_surveys_log)

distribution_plot(contributions_surveys_log, longer = T,
                  cols_not_plot = var_metadata)

ggsave(plot = last_plot(), width=15, height= 10,
       filename = here::here("figures", "2_contributions_distribution_log_transformed.jpg"))



##------------------- Save log-transformed contributions -------------------####
#remove total biomass (redundant information)
contributions <- contributions_surveys_log |> 
  tibble::column_to_rownames("survey_id") |> 
  dplyr::select(-total_biomass, -any_of(var_metadata))

# Save #
save(contributions, file = here::here("outputs", "2_all_contributions.Rdata"))



##------------------- Check NA's distribution -------------------####
library(funbiogeo)
df <- dplyr::rename(contributions_surveys_log, species = survey_id) |> 
  dplyr::select(-any_of(var_metadata))
fb_plot_species_traits_completeness(df)
ggsave(plot = last_plot(), width=15, height= 10,
       filename = here::here("figures", "2_contributions_completedness.jpg"))
fb_plot_number_species_by_trait(df)

#NA position on map
NA_on_map(data=contributions_surveys_log, variable = "iucn_species_richness",
          xlim = c(-180,180),ylim = c(-40,40),
          jitter = 2, lat_line = 30, 
          priority= "no") #which points are display on the others


##------------------- Correlations between contributions -------------------####
## Corr-matrix for all Contributions (log transformed)
M <- cor(dplyr::select(contributions_surveys_log, 
                       -total_biomass, -any_of(var_metadata)) |> 
           tidyr::drop_na()) #Many loss due to recycling

# M <- cor(dplyr::select(contributions_surveys_log,-total_biomass, -any_of(var_metadata),
#                        -N_recycling, -P_recycling) |> 
#            tidyr::drop_na()) 

png(filename = here::here("figures","2h_corr_matrix_contributions.png"), 
    width= 40, height = 30, units = "cm", res = 1000)
print({
  corrplot::corrplot(M, order = 'AOE', type = 'lower', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r')
  corrplot::corrplot(M, add = TRUE, type = 'upper', method = 'number', order = 'AOE', insig = 'p-value',
                     diag = FALSE, tl.pos = 'n', cl.pos = 'n', number.digits = 1)
})
dev.off() 

## Study correlogram
pairwise_corr <- M[upper.tri(M)]
summary(pairwise_corr) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.65994 -0.08566  0.07111  0.08510  0.20663  0.99300 
summary(abs(pairwise_corr))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001281 0.072355 0.144830 0.203633 0.285199 0.993001  



##Classify contributions into 2 categories: Nature for Nature (NN) and Nature 
## for People (NP)
colnames(contributions)
grp_NN_NP <- as.factor(c(actino_richness = "NN",
                         elasmobranch_richness = "NN",
                         functional_distinctiveness = "NN",
                         iucn_species_richness = "NN", 
                         mean_endemism = "NN",
                         evolutionary_distinctiveness = "NN",
                         functional_entropy = "NN",
                         phylogenetic_entropy = "NN",
                         herbivores_biomass = "NN",
                         invertivores_biomass = "NN",
                         piscivores_biomass = "NN",
                         N_recycling = "NN",
                         P_recycling = "NN",
                         trophic_web_robustness = "NN",
                         mean_trophic_level = "NN",
                    
                         available_biomass = "NP",
                         selenium = "NP",
                         zinc = "NP",
                         omega_3 = "NP",
                         calcium = "NP",
                         iron = "NP",
                         vitamin_A = "NP",
                         available_biomass_turnover = "NP")) # /!\ the order matter

### PCA with weigths according categories and items (see table 1, Flandrin et al. 2023)
w <- c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6,
       1/5, 1/5, 1/5, 1/5, 1/5,
       1/2, 1/2, 
       1/2, 1/2,
       1/2, 
       1/6, 1/6, 1/6, 1/6, 1/6, 1/6,
       1/2)
names(w) <- colnames(contributions)
w

pca <- FactoMineR::PCA(contributions, col.w = w, scale.unit = T, graph=F, ncp=15)

factoextra::fviz_screeplot(pca, ncp=15)

color_grad = c("#A50026", "#D73027", "#F46D43", "#FDAE61","#FEE090", "#D8DAEB",
               "#B2ABD2", "#8073AC", "#6D469C", "#603692", "#542788",
               "#473C8B")

png(filename = here::here("figures","2h_PCA_contributions_axes1-2.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
factoextra::fviz_pca_biplot(pca, col.var = grp_NN_NP, 
                            title="",
                            palette = c("forestgreen", "dodgerblue3"),
                            labelsize = 4, 
                            repel = TRUE,
                            geom="point", pointshape=21,
                            stroke=0, pointsize=3,
                            alpha.ind = 0.7,
                            fill.ind = contributions_surveys_log$total_biomass,
                            gradient.cols = rev(color_grad))+
  labs(col = "Reef Contributions", fill = "log(total biomass)")+
  scale_color_discrete(type=c("forestgreen", "dodgerblue3"),labels = c("NN", "NP"))+
  theme(legend.position = "bottom")

dev.off()



##-------------plot Contributions on map-------------
#plot function
plot_Contrib_on_world_map <- function(NCP = "Taxonomic_Richness", xlim=c(-180,180), ylim = c(-36, 31),
                                      title="world map with ", jitter=1.5, pt_size=2){
  data <- contributions_surveys_log[order(contributions_surveys_log[,NCP]),]
  library(ggplot2)
  map <- ggplot(data) +
    geom_sf(data = coast, color = "grey30", fill = "lightgrey",
            size=0.1) +
    
    geom_point(data=data,
               size = pt_size, shape = 20,
               position=position_jitter(width=jitter, height = jitter),
               aes(x = longitude, y = latitude,
                   colour= data[,NCP],
                   alpha = 0.7)) +
    scale_colour_gradientn(name  = NCP,
                           colours = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu")))+
    
    
    coord_sf(xlim, ylim, expand = FALSE) +
    guides(alpha= "none" ) +
    # scale_size_continuous(range = c(0.5, 4), guide = "none") +
    theme_minimal()+
    labs(#title = paste0(NCP, " geographic distribution"),
      x="", y= "") +
    theme(legend.position = "bottom",
          plot.title = element_text(size=10, face="bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
  
  ggsave( here::here("figures", "map_contributions", paste0( title , NCP, ".jpg")), plot = map, width=15, height = 7 )
  #map
}

# save world maps
parallel::mclapply(colnames(contributions), function(NCP){
  plot_Contrib_on_world_map(NCP, xlim=c(-180,180), ylim = c(-36, 31), 
                            title="world_map_with_", jitter=1.5, pt_size=2)
},mc.cores=parallel::detectCores()-5)



#-------------Compute weighted mean of Contributions: NN and NP scores-------------
scaled_contrib <- scale(contributions)


grp_NN_NP <- as.factor(c(actino_richness = "NN",
                         elasmobranch_richness = "removed from mean",
                         functional_distinctiveness = "NN",
                         iucn_species_richness = "NN", 
                         mean_endemism = "NN",
                         evolutionary_distinctiveness = "NN",
                         functional_entropy = "NN",
                         phylogenetic_entropy = "NN",
                         herbivores_biomass = "NN",
                         invertivores_biomass = "NN",
                         piscivores_biomass = "NN",
                         N_recycling = "NN",
                         P_recycling = "NN",
                         trophic_web_robustness = "NN",
                         mean_trophic_level = "NN",
                         
                         available_biomass = "NP",
                         selenium = "NP",
                         zinc = "NP",
                         omega_3 = "NP",
                         calcium = "NP",
                         iron = "NP",
                         vitamin_A = "NP",
                         available_biomass_turnover = "NP")) # /!\ the order matter


## Nature to Nature (NN) score ##
NN_names <- names(grp_NN_NP)[ grp_NN_NP=="NN" ]
contrib_NN <- scaled_contrib[,NN_names]

colnames(contrib_NN)
weighting_par <- c(1/5, 1/5, 1/5, 1/5, 1/5, #Biodiversity
                   1/5, 1/5, 1/5, 1/5, 1/5, #Biomass distribution
                   1/2, 1/2, #recycling
                   1/2,1/2) # trophic web
names(weighting_par) <- colnames(contrib_NN)
weighting_par

mean_NN <- c()
for( site in 1:nrow(contrib_NN)){
  EDS_site <- sum(weighting_par * contrib_NN[site,], na.rm=T) / sum(weighting_par)
  mean_NN <- c(mean_NN, EDS_site) }

summary(mean_NN)


## Nature to People (NP) score ##
NP_names <- names(grp_NN_NP)[ grp_NN_NP=="NP" ]
contrib_NP <- scaled_contrib[,NP_names]

colnames(contrib_NP)
weighting_par <- c(1/2, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/2)
names(weighting_par) <- colnames(contrib_NP)
weighting_par

mean_NP <- c()
for( site in 1:nrow(contrib_NP)){
  EDS_site <- sum(weighting_par * contrib_NP[site,], na.rm=T) / sum(weighting_par)
  mean_NP <- c(mean_NP, EDS_site) }

summary(mean_NP)


contributions_with_synthetic_score <- as.data.frame( contributions ) |> 
  dplyr::mutate(NN_score = mean_NN,
                NP_score = mean_NP)

# plot(mean_NP ~ mean_NN)
# cor.test(mean_NP, mean_NN)

# Save #
save(contributions_with_synthetic_score, 
     file = here::here("outputs", "2_all_contributions_with_synthetic_score.Rdata"))

# load(file = here::here("outputs", "2_all_contributions_with_synthetic_score.Rdata"))




###############################################################################"
##
##                            #### SITE SCALE ####
##
###############################################################################"

#------------------- Mean contribution values at the site X year scale -------------------####

#TO REDUCE NOISES IN DATA, WE MEANS THE CONTRIBUTIONS VALUES OF SURVEYS, 
# IN A GIVEN SITE, OBSERVED IN THE SAME DATE

# Mean the contribution at the site scale, for a given date: all surveys in the
# same place, observed at the same date are merged



contributions_sites <- contributions_surveys |> 
  dplyr::select(-survey_id) |> 
  dplyr::group_by(site_code, latitude, longitude, country, ecoregion, realm, 
                  survey_date, year, effectiveness) |> 
  dplyr::summarise(across(.cols = everything(),
                          .fns = ~mean(., na.rm = TRUE), .names = "{.col}")) |> 
  dplyr::mutate(across(.cols = everything(),
                       .fns = ~ifelse(is.nan(.), NA, .), .names = "{.col}")) |> 
  dplyr::mutate(id = paste0(site_code, "_", survey_date)) |> 
  dplyr::ungroup() 

nrow(contributions_sites) # 3325 mean_surveys

distribution_plot(contributions_sites, longer = T,cols_plot = colnames(contributions_sites)[13:36])


# Log transformation
contributions_sites_log <- contributions_sites |>
  dplyr::mutate(across(.cols = all_of(zero_values),
                       .fns = ~ .x +1 , .names = "{.col}")) |>
  dplyr::mutate(across(.cols = all_of(to_log),
                       .fns = ~ifelse(is.na(.), NA, log10(.)) , .names = "{.col}")) # log(x+1) to avoid -Inf values

#remove total biomass (redundant information)
contributions_sites_date <- contributions_sites_log |> 
  tibble::column_to_rownames("id") |> 
  dplyr::select(-total_biomass, -any_of(var_metadata))

#Check NA and distributions
df <- tibble::rownames_to_column(contributions_sites_date, "id") |> 
  dplyr::rename(species = id)
funbiogeo::fb_plot_species_traits_completeness(df)


distribution_plot(contributions_sites_date, longer = T,
                  cols_plot = colnames(contributions_sites_date))


#------------------- Compute weighted mean of contributions -------------------####
to_mean <- scale(contributions_sites_date)

## Nature to Nature (NN) score ##
NN_names <- names(grp_NN_NP)[ grp_NN_NP=="NN" ]
contrib_NN <- to_mean[,NN_names]

colnames(contrib_NN)
weighting_par <- c(1/5, 1/5, 1/5, 1/5, 1/5, #Biodiversity
                   1/5, 1/5, 1/5, 1/5, 1/5, #Biomass distribution
                   1/2, 1/2, #recycling
                   1/2,1/2) # trophic web
names(weighting_par) <- colnames(contrib_NN)
weighting_par

mean_NN <- c()
for( site in 1:nrow(contrib_NN)){
  EDS_site <- sum(weighting_par * contrib_NN[site,], na.rm=T) / sum(weighting_par)
  mean_NN <- c(mean_NN, EDS_site) }

summary(mean_NN)


## Nature to People (NP) score ##
NP_names <- names(grp_NN_NP)[ grp_NN_NP=="NP" ]
contrib_NP <- to_mean[,NP_names]

colnames(contrib_NP)
weighting_par <- c(1/2, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/2)
names(weighting_par) <- colnames(contrib_NP)
weighting_par

mean_NP <- c()
for( site in 1:nrow(contrib_NP)){
  EDS_site <- sum(weighting_par * contrib_NP[site,], na.rm=T) / sum(weighting_par)
  mean_NP <- c(mean_NP, EDS_site) }

summary(mean_NP)

# plot(mean_NP ~ mean_NN)
# cor.test(mean_NP, mean_NN)


contributions_sites_date <- contributions_sites_date |> 
  dplyr::mutate(NN_score = mean_NN,
                NP_score = mean_NP)


# Save #
save(contributions_sites_date, file = here::here("outputs", "2_contributions_site&date.Rdata"))
# load(file = here::here("outputs", "2_contributions_site&date.Rdata"))


## Check values
distribution_plot(contributions_sites_date, cols_plot = colnames(contributions_sites_date))
pca <- FactoMineR::PCA(contributions_sites_date, scale.unit = T, graph=F, ncp=15,
                       quanti.sup = c("NN_score", "NP_score"))
factoextra::fviz_pca_biplot(pca, repel = TRUE, geom="point", pointshape=21,
                            stroke=0, pointsize=3, alpha.ind = 0.7, 
                            fill.ind = "grey", col.quanti.sup = "firebrick")



#------------------- Mean contribution values at the site scale -------------------####
# Mean the contribution at the site scale, whatever the date:
#  we offset the temporal information = look only on the spatial patterns.

synthetic_scores <- contributions_with_synthetic_score[,c( "NN_score", "NP_score")] |> 
  tibble::rownames_to_column("survey_id")

contributions_sites_no_date <- contributions_surveys |> 
  dplyr::left_join(synthetic_scores) |> 
  dplyr::select(-survey_id, -survey_date) |> 
  dplyr::group_by(site_code, latitude, longitude, country, ecoregion,
                  realm, effectiveness) |> 
  dplyr::summarise(across(.cols = everything(),
                          .fns = ~mean(., na.rm = TRUE), .names = "{.col}")) |> 
  dplyr::mutate(across(.cols = everything(),
                       .fns = ~ifelse(is.nan(.), NA, .), .names = "{.col}")) |> 
  dplyr::mutate(id = paste0(site_code, "_", effectiveness)) |>
  dplyr::ungroup()

# /!\ SOME SITES HAVE CHANGED PROTECTION STATUS BETWEEN DIFFERENT SURVEYS

nrow(contributions_sites_no_date) # 2128 sites

contributions_sites_no_date <- contributions_sites_no_date |>
  dplyr::mutate(across(.cols = all_of(zero_values),
                       .fns = ~ .x +1 , .names = "{.col}")) |>
  dplyr::mutate(across(.cols = all_of(to_log),
                       .fns = ~ifelse(is.na(.), NA, log10(.)) , .names = "{.col}")) |>  # log(x+1) to avoid -Inf values
  tibble::column_to_rownames("id") |> 
  dplyr::select(-total_biomass, -any_of(var_metadata))

# Save #
save(contributions_sites_no_date, file = here::here("outputs", "2_contributions_site_NO_date.Rdata"))
# load( file = here::here("outputs", "2_contributions_site_NO_date.Rdata"))

