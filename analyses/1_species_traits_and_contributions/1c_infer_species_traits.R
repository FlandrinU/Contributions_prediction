################################################################################
##
##  
##
## 1c_infer_species_traits.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
###############################################################################"
#----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("here", "dplyr", "ggplot2", "funbiogeo", "missForest", "slam",
#           "pbmcapply", "patchwork", "ggplot2", "maditr","rnaturalearth",
#           "tibble", "stringr")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

library(funbiogeo)
library(ggplot2)

##-------------loading data and functions-------------
# tropical species traits #
load(file = here::here("data","derived_data","tropical_species_traits.Rdata"))

# functions #
source(here::here("R","MissForest_functions.R"))
source(here::here("R","evaluation_prediction_model.R"))


##-------------Check for traits distribution-------------
numeric_cols <- colnames(tropical_species_traits)[sapply(tropical_species_traits, class) == "numeric"]
df <- tidyr::pivot_longer(tropical_species_traits, cols = all_of(numeric_cols) ,
                          names_to = "trait", values_to = "value")

ggplot(data = df)+
  aes(x=value, group=trait, fill=trait) +
  geom_histogram(aes(y = ..density.., fill = trait), bins = 20,
                 color = "grey40", alpha = 0.2) +
  hrbrthemes::theme_ipsum() +
  facet_wrap(~trait, scales = "free") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )

#DEPTH MIN AND DEPTH MAX + GEOGRAPHIC RANGE ARE HIGLY RIGHT-SKEWED -> LOG TRANSFORMATION
cols <- c("depth_min", "depth_max", "geographic_range_Albouy19", "range_n_cells_01")

tropical_species_traits[,cols] <- log10(tropical_species_traits[,cols] +1) # log10(X+1)
colnames(tropical_species_traits)[colnames(tropical_species_traits) %in% cols] <- 
  c("log(depth_min)", "log(depth_max)", "log(geographic_range_Albouy19)", "log10(range_n_cells_01)")


##-------------Complete IUCN data from known categories-------------
## IUCN data: IUCN categories used to complete inference by Loiseau et al. (2023)
## mainly for Elasmobranch (that are not taken in Loiseau's study)
iucn <- tropical_species_traits |> 
  tibble::column_to_rownames("rls_species_name") |> 
  dplyr::select(class, fishbase_name,
                iucn_inferred = IUCN_inferred_Loiseau23, 
                iucn_redlist = IUCN_category) 

missing_values <- dplyr::filter(iucn, is.na(iucn$iucn_inferred)) |>
  dplyr::mutate(iucn_redlist = dplyr::recode(iucn_redlist,  
                                             "NE" = "No Status",
                                             "DD" = "No Status",
                                             "LC" = "Non Threatened",
                                             "NT" = "Non Threatened",
                                             "VU" = "Threatened",
                                             "EN" = "Threatened",
                                             "CR" = "Threatened")) 

iucn[rownames(missing_values), "iucn_inferred"] <- missing_values$iucn_redlist

tropical_species_traits[, "IUCN_inferred_Loiseau23"] <- 
  iucn[tropical_species_traits$rls_species_name, "iucn_inferred"]




trait_completedness <- dplyr::rename(tropical_species_traits,
                                species = rls_species_name) |> 
  dplyr::select(-phylum, -class, -order, -family, -spec_code, -worms_id)

fb_plot_species_traits_completeness(trait_completedness)
fb_plot_number_species_by_trait(species_traits, threshold_species_proportion = 1)


################################################################################
##
##         Integrate phylogeny with Principal Component Analyses
##
###############################################################################"

##-------------Choose the number of dimensions-------------
phylogeny <- tropical_species_traits |> 
  dplyr::ungroup() |> 
  dplyr::select(rls_species_name, phylum, class, order, family) |> 
  dplyr::mutate(genus = stringr::word(rls_species_name,1)) |> 
  tibble::column_to_rownames(var="rls_species_name") |> 
  as.matrix()

#Integrate phylogeny with MCA
dimensionality <- FactoMineR::MCA(phylogeny,ncp=10) 
factoextra::fviz_eig(dimensionality, choice = "variance",
                     addlabels = F, ylim = c(0, 0.75), barcolor = "darkred",
                     barfill = "darkred", ncp=550)


# Choose the number of dimensions
dimensions <- c(0,2,5,10,15,20,25,30,35,40,45,50,55,60,75,100,200,300)


## Run missforest for each dimension of phylogeny
test_dimensionality <- lapply(dimensions, FUN = function(dim){
  #dim=2
  cat(dim, "dimensions:", "\n")
  
  if( dim == 0 ){
    
    #No phylogeny
    traits_and_phylo <- tropical_species_traits
    #prep data
    data_for_missforest <- preping_data(species_traits_df =traits_and_phylo)
    
    data_to_infer <- data_for_missforest[[1]]
    traits_data_factors <- data_for_missforest[[2]]
    traits_data_num <- data_for_missforest[[3]]
    factor_length <- data_for_missforest[[4]]
    
  }else{
    
    dimensionality <- FactoMineR::MCA(phylogeny,ncp=dim) 
    
    phylo_space <- dimensionality[["ind"]][["coord"]]
    colnames(phylo_space)<- gsub(" ", "", colnames(phylo_space))
    
    traits_and_phylo <- cbind(phylo_space, tropical_species_traits)
    rownames(traits_and_phylo) <- NULL
    
    data_for_missforest <- preping_data(species_traits_df =traits_and_phylo)
    
    data_to_infer <- data_for_missforest[[1]]
    traits_data_factors <- data_for_missforest[[2]]
    traits_data_num <- dplyr::select(data_for_missforest[[3]],-colnames(phylo_space))
    factor_length <- data_for_missforest[[4]]
    
  } # data prepped
  
  
  # Run Missforest
  model_eval_missforest_MCA <- fct_missforest_evaluation_MCA(
    data_to_infer, traits_data_factors = traits_data_factors,
    traits_data_num = traits_data_num, factor_length = factor_length,
    model_iteration = 10, #100
    maxiter_mf=10, #10
    ntree_mf=100,  #100
    prop_NA = 0.2)

  # Extract results
  results <- extract_model_perf(raw_result = model_eval_missforest_MCA)
  
  traits_performance <- results[[1]] |> 
    dplyr::mutate(dimensions = dim)
  
  traits_performance
}) #END of lapply on dimension number

# #time of computation:
# time <- c("0 dim 04:46", "2 dim 04:31", "5 dim 03:50", "10 dim 05:44", "25 dim 04:52",
#           "50 dim 06:16", "75 dim 07:35", "100 dim 10:05", "200 dim 13:05")

dimensionality_res <- do.call(rbind, test_dimensionality) |> 
  dplyr::mutate(dimensions = as.factor(dimensions))  
# dplyr::filter(variable == "K")

save(dimensionality_res, file = here::here("outputs","1c_choose_dimensionality_phylogeny_mF.Rdata"))
# load(file = here::here("outputs","1c_choose_dimensionality_phylogeny_mF.Rdata"))

resumed_data <- dimensionality_res |>
  dplyr::group_by(variable, dimensions) |>
  dplyr::summarise(median_estimate = median(estimate))

## Plot performance against the number of dimensions
ggplot(resumed_data) +
  aes(x= dimensions, y= median_estimate, fill= dimensions)+
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.15, alpha=0.2)+
  xlab("") + ylab("Assessement quality (median of R-squared and Accuracy among traits)") +
  theme_minimal() +
  scale_fill_manual(values = grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(9, "OrRd"))(length(unique(resumed_data$dimensions)))) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10))
  
ggsave(filename = here::here("figures","1_species_traits", "1_phylogeny_importance_in_MF_perf.jpg"),
       plot = last_plot(), width = 10, height =8 )


## Plot performance of prediction for each dimensionality ##
dim = 35
res <- dimensionality_res |> dplyr::filter(dimensions == dim)
# Estimate distributions (boxplot)
boxplot_missforest <- estimates_boxplot(df_estimates = res)
boxplot_missforest
ggsave(filename = here::here("figures","1_species_traits", paste0("1_Missforest_performance_boxplot_traits_phylogeny_",dim,".jpg")),
       boxplot_missforest, width = 12, height =8 )

# Estimate distribution (Histograms)
hist_missforest <- estimates_histogramm(data = res)
hist_missforest




##-------------Assess error in MissForest-------------

## number of dimensions chosen 
dim = 35
##

## Preping data
phylogeny <- tropical_species_traits |> 
  dplyr::ungroup() |> 
  dplyr::select(rls_species_name, phylum, class, order, family) |> 
  dplyr::mutate(genus = stringr::word(rls_species_name,1)) |> 
  tibble::column_to_rownames(var="rls_species_name") |> 
  as.matrix()

dimensionality <- FactoMineR::MCA(phylogeny,ncp=dim) 

phylo_space <- dimensionality[["ind"]][["coord"]]
colnames(phylo_space)<- gsub(" ", "", colnames(phylo_space))

traits_and_phylo <- cbind(phylo_space, tropical_species_traits)
rownames(traits_and_phylo) <- NULL

data_for_missforest <- preping_data(species_traits_df = traits_and_phylo)

data_to_infer <- data_for_missforest[[1]]
traits_data_factors <- data_for_missforest[[2]]
traits_data_num <- dplyr::select(data_for_missforest[[3]],-colnames(phylo_space))
factor_length <- data_for_missforest[[4]]


## Run Missforest evaluation
model_eval_missforest_MCA <- fct_missforest_evaluation_MCA(
  data_to_infer, traits_data_factors , traits_data_num, factor_length,
  model_iteration = 50, #100
  maxiter_mf=10, #10
  ntree_mf=100,  #100
  prop_NA = 0.2)


save(model_eval_missforest_MCA, file = here::here("outputs", "predictive_model_eval_species_traits.Rdata"))
# load(file = here::here("outputs", "predictive_model_eval_species_traits.Rdata"))


## Extract results
results <- extract_model_perf(raw_result = model_eval_missforest_MCA)
traits_performance <- results[[1]]
order_performance <- results[[2]]
raw_factors_perf <- results[[3]]
raw_numeric_perf <- results[[4]]


## Plot performance of prediction
# Estimate distributions (boxplot)
boxplot_missforest <- estimates_boxplot(df_estimates = traits_performance)
boxplot_missforest
ggsave(filename = here::here("figures","1_species_traits", "1_Missforest_final_performance_boxplot_traits.jpg"),
       boxplot_missforest, width = 12, height =8 )

#Choose variable to plot
var_to_plot <- c("Length", "K", "trophic_guild",
                  "Calcium", "Iron", "Omega3", "Selenium", "VitaminA", "Zinc",
                  "log(geographic_range_Albouy19)",
                  "public_attention", "aesthetic")

#Plot performance for these traits
boxplot_missforest <- estimates_boxplot(df_estimates = dplyr::filter(
  traits_performance, variable %in% var_to_plot),
  add_mean = F)
boxplot_missforest
ggsave(filename = here::here("figures","1_species_traits", "1_Missforest_performance_distribution_selected_traits.jpg"),
       boxplot_missforest, width = 10, height =7 )

# Estimate distributions (Histograms)
hist_missforest <- estimates_histogramm(data = traits_performance)
hist_missforest
ggsave(filename = here::here("figures","1_species_traits", "1_Missforest_performance_distribution.png"),
       hist_missforest, width = 22, height =14 )

boxplot_missforest_order <- estimates_boxplot_per_order(df_estimates = order_performance)
boxplot_missforest_order #very large plot
ggpubr::ggexport(boxplot_missforest_order, 
                 filename = here::here("figures","1_species_traits", "1_Missforest_performance_boxplot_order.pdf"),
                 width = 16, height =11)


## Sum up missforest performance
max_proportion <- apply(data_to_infer, 2, FUN = function(x){
  table(x)[which.max(table(x))]
}) / nrow(data_to_infer)

resumed_data <- traits_performance |> 
  dplyr::group_by(variable) |> 
  dplyr::summarise(median_estimate = median(estimate))

prop_NA <- c()
for( i in resumed_data$variable){
  p_NA <- sum(is.na(tropical_species_traits[,i]))/nrow(tropical_species_traits)
  prop_NA <- c(prop_NA, p_NA)
}

#assess the risk: proportion of NA * (1- Rsquared or Accuracy) ~= proportion of error in the final matrix
resumed_data <- resumed_data |> 
  dplyr::mutate(prop_NA = prop_NA,
                prop_error = prop_NA * (1-median_estimate))

useful <- c("a", "b", "DemersPelag", "Length", "Calcium", "Iron", "Omega3",
            "Selenium", "VitaminA", "Zinc", "K", "ClimVuln_SSP585", "public_interest",
            "aesthetic", "range_n_cells_01", "IUCN_category", "trophic_guild", 
            "recycling_P", "recycling_N")
hist(resumed_data$prop_error, breaks = 20)
hist(resumed_data[resumed_data$variable %in% useful,]$prop_error, breaks = 20)
resumed_data[resumed_data$prop_error < 0.05,]$variable



##------------- Infer data on chosen variables-------------

#Choose variable to infer
var_to_infer <- c("Length", "K", "trophic_guild",
                  "Calcium", "Iron", "Omega3", "Selenium", "VitaminA", "Zinc",
                  "ClimVuln_SSP585", "geographic_range_Albouy19",
                  "public_attention", "aesthetic")

#Number of dimensions chosen 
dim = 35

#Preping data
phylogeny <- tropical_species_traits |> 
  dplyr::ungroup() |> 
  dplyr::select(rls_species_name, phylum, class, order, family) |> 
  dplyr::mutate(genus = stringr::word(rls_species_name,1)) |> 
  tibble::column_to_rownames(var="rls_species_name") |> 
  as.matrix()

dimensionality <- FactoMineR::MCA(phylogeny,ncp=dim) 

phylo_space <- dimensionality[["ind"]][["coord"]]
colnames(phylo_space)<- gsub(" ", "", colnames(phylo_space))

traits_and_phylo <- cbind(phylo_space, tropical_species_traits)
rownames(traits_and_phylo) <- NULL

data_for_missforest <- preping_data(species_traits_df = traits_and_phylo)

data_to_infer <- data_for_missforest[[1]]
traits_data_factors <- data_for_missforest[[2]]
traits_data_num <- dplyr::select(data_for_missforest[[3]],-colnames(phylo_space))
factor_length <- data_for_missforest[[4]]


#Run several missforest, and keep the most frequent inference if most models 
# converge to the same result
inferred_data <- 
  missforest_applied(
    data_to_infer, 
    factor_length,
    traits_data_factors,
    traits_data_num,
    var_to_infer,
    confidence_threshold = 0.8, #factors: proportion of consistency with the norm, 
                                # or numeric: 1-Coefficient of variation > 0.8
    model_iteration = 50,
    maxiter_mf=10, #missforest parameter
    ntree_mf=100) #missforest parameter


# Explore data
species_traits <- inferred_data |> 
  tibble::rownames_to_column(var = "species") |> 
  dplyr::select(-phylum, -class, -order)

fb_plot_species_traits_completeness(species_traits)
fb_plot_number_species_by_trait(species_traits, threshold_species_proportion = 1)
ggsave(plot= last_plot(), width = 8, height = 8, 
       file= here::here("figures","1_species_traits","1_percent_species_INFERRED_TRAITS_tropical.png"))

#check trophic guild inference
df <- data.frame(troph_fishbase  = data_to_infer$Troph[which(is.na(data_to_infer$trophic_guild))],
                 trophic_guild_inferred = inferred_data$trophic_guild[which(is.na(data_to_infer$trophic_guild))])

ggplot(data=df, aes(x=troph_fishbase, group=trophic_guild_inferred, fill=trophic_guild_inferred)) +
  geom_histogram(aes(y = ..density..), bins = 20, color = "grey40", fill ="white") +
  geom_density(aes(fill = trophic_guild_inferred), alpha = 0.2) +
  hrbrthemes::theme_ipsum() +
  facet_wrap(~trophic_guild_inferred, scales = "free_y") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )
ggsave(filename = here::here("figures","1_species_traits", "1_trophic_guild_inference_VS_troph.jpg"),
       width = 12, height = 8)


##-------------save data-------------
inferred_species_traits <- inferred_data |> 
  tibble::rownames_to_column("rls_species_name") |> 
  dplyr::left_join(
    dplyr::select(tropical_species_traits, 
                  rls_species_name, family, spec_code, worms_id, fishbase_name)) |> 
  tibble::column_to_rownames("rls_species_name")
  

save(inferred_species_traits, file= here::here("outputs", "RLS_species_traits_inferred.Rdata"))
# load(file= here::here("outputs", "RLS_species_traits_inferred.Rdata"))



#Save used traits completedness
used_traits <- c("fishbase_name", "IUCN_category", "Length", "Importance",
                 "Troph", "a", "b", "K", "DemersPelag", "trophic_guild",
                 "IUCN_inferred_Loiseau23", 
                 "log(depth_min)", "log(depth_max)", "log(geographic_range_Albouy19)", "Iron",
                 "VitaminA" , "Zinc", "Calcium", "Omega3", "Selenium" ,
                 "aesthetic", "public_attention")
species_traits <- dplyr::select(inferred_species_traits, all_of(used_traits)) |> 
  tibble::rownames_to_column("species")
fb_plot_number_species_by_trait(species_traits, threshold_species_proportion = 1)
ggsave(plot= last_plot(), width = 8, height = 8, 
       file= here::here("figures","1_species_traits","1_percent_species_INFERRED_used_TRAITS_tropical.png"))
