###############################################################################"
##
##
##
## 3d_conterfactuals.R
##
## 28/06/2024
##
## Ulysse Flandrin
##
###############################################################################"

##---------------------------- cleaning memory ---------------------------------
rm(list=ls())

##-----------------------------Loading packages---------------------------------
# pkgs <- c("here", "Hmsc", "ggplot2")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# # ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))
library(ggplot2)
library(patchwork)

source(here::here("R/HMSC_function.R"))
source(here::here("R","evaluation_prediction_model.R"))


##------------------------------- load data ------------------------------------
## Survey scale, all covariates
load(here::here("data/derived_data/3_all_covariates_to_predict.Rdata"))
X_data = covariates_final

# ## Without Allen
# load(here::here("data/derived_data/3_covariates_without_Allen_to_predict.Rdata"))
# X_data_wo_allen = covariates_final_without_Allen

## Site scale
load(here::here("data/derived_data/3_sites_covariates_to_predict.Rdata"))
X_data_site = covariates_site_final


# ## Only Australia
# X_data_aust = covariates_final |> dplyr::filter(country == "Australia")
# Y_data_aust =  observations_final[rownames(X_data_aust),]


##load metadata 
load(file = here::here("data", "raw_data", "environmental_covariates",
                       "all_covariates_benthos_inferred_tropical_surveys.Rdata"))

metadata <- all_covariates_benthos_inferred |> 
  dplyr::select(survey_id:year) |> 
  tibble::column_to_rownames("survey_id")

## If site scale
metadata_sites <- all_covariates_benthos_inferred |>
  dplyr::select(country:year, -depth, -visibility, -hour) |>
  dplyr::mutate(id = paste0(site_code, "_", survey_date)) |> 
  unique()
rownames(metadata_sites) <- NULL
metadata_sites <- tibble::column_to_rownames(metadata_sites, "id")


# Path to model folder
path = here::here("outputs/models/hmsc")



##----------------------------- Choose model ------------------------------

## List all files in the directory and choose the model
list_files <- list.files(file.path(path, "out_multi")) 
list_files
model_name <- gsub("output_", "", list_files[14]) #choose the wanted file
concatenate_chains = F


##------------- New conditions in counterfactual scenarios ---------------------

### Initial conditions ###
X <- X_data_site
metadata = metadata_sites
# X <- X_data


#### Change initial conditions ###
summary(X)


# (1) HUMAN FOOTPRINT - Change unprotected sites to: full protection, 
#  min fishing vessels, min gravity and max neartt
X_pristine_conditions <- X
new_pristine <- rownames(X_pristine_conditions |> dplyr::filter(protection_status2 == "out"))

X_pristine_conditions[new_pristine, "protection_status2"] <- as.factor("full")
X_pristine_conditions[new_pristine, "n_fishing_vessels"] <- min(X_pristine_conditions$n_fishing_vessels)
X_pristine_conditions[new_pristine, "gravtot2"] <- min(X_pristine_conditions$gravtot2)
X_pristine_conditions[new_pristine, "neartt"] <- min(X_pristine_conditions$neartt)


# (2) CONSERVATION LEGACY - Change protected sites to: out MPA, mean fishing vessels of "out" sites
X_conservation_legacy <- X
# new_conserv_legacy <- rownames(X_conservation_legacy |> dplyr::filter(protection_status == "full"))
new_conserv_legacy <- rownames(X_conservation_legacy |> dplyr::filter(protection_status2 != "out"))

#protection status
X_conservation_legacy[new_conserv_legacy, "protection_status2"] <- as.factor("out")

#fishing vessels
mean_fishing_out <- mean(X_conservation_legacy[new_pristine, "n_fishing_vessels"])
mean_country_out <- X_conservation_legacy[new_pristine, ] |> 
  dplyr::group_by(country) |> 
  dplyr::summarise(fishing_out = mean(n_fishing_vessels))
missing_rows <- data.frame(
  country = unique(X$country)[!unique(X$country) %in% mean_country_out$country],
  fishing_out = mean_fishing_out)
mean_country_out <- rbind(mean_country_out, missing_rows)


X_conservation_legacy <- X_conservation_legacy |> 
  tibble::rownames_to_column("id") |> 
  dplyr::left_join(mean_country_out) |> 
  dplyr::mutate(n_fishing_vessels = dplyr::case_when(
    id %in% new_conserv_legacy ~ fishing_out,
    TRUE ~ n_fishing_vessels)) |> 
  tibble::column_to_rownames("id") |> 
  dplyr::select(-fishing_out)
# X_conservation_legacy[new_conserv_legacy, "n_fishing_vessels"] <- mean_fishing_out





#### Explore other changes ###

#(3) Change MPA only: from "out" to "full protection"
X_new_mpa <- X
new_mpa <- rownames(X_new_mpa |> dplyr::filter(protection_status2 == "out"))
X_new_mpa[new_mpa, "protection_status2"] <- as.factor("full")

#(4) Change fishing pressure only: set the number of fishing vessel to the minimum known.
X_no_vessels <- X
new_vessels <- rownames(
  X_no_vessels[X_no_vessels$n_fishing_vessels != min(X_no_vessels$n_fishing_vessels),]
)
X_no_vessels[new_vessels, "n_fishing_vessels"] <- min(X_no_vessels$n_fishing_vessels)

#(5) Conservation potential: unprotected sites are placed in reserves, without fishing pressure.
X_new_mpa_no_vessels <- X
new_mpa_no_vessels <- unique(c(new_mpa, new_vessels))
X_new_mpa_no_vessels[new_mpa_no_vessels, "protection_status2"] <- as.factor("full")
X_new_mpa_no_vessels[new_mpa_no_vessels, "n_fishing_vessels"] <- min(X_new_mpa_no_vessels$n_fishing_vessels)

#(6) Low gravity
X_low_gravity <- X
new_low_gravity <- rownames(
  X_low_gravity[X_low_gravity$gravtot2 != min(X_low_gravity$gravtot2) ,])
X_low_gravity[new_low_gravity, "gravtot2"] <- min(X_low_gravity$gravtot2)

#(7) No human population
X_no_human <- X_low_gravity
new_no_human <-  unique(rownames(X_no_human[X_no_human$neartt != max(X_no_human$neartt),]),
                   new_low_gravity)
X_no_human[new_no_human, "neartt"] <- max(X_no_human$neartt)


#(8) Change MPA only: from "out" to medium protection ("restricted")
X_new_medium_mpa <- X
new_medium_mpa <- rownames(X_new_medium_mpa |> dplyr::filter(protection_status2 == "out"))
X_new_medium_mpa[new_medium_mpa, "protection_status2"] <- as.factor("restricted")


##----------------------------- Plot changes -----------------------------------
# set_ids = new_mpa #look only at currently unprotected sites
# set_ids = new_mpa_no_vessels
# set_ids = new_pristine
# set_ids = new_no_human

# Select countries we want to plot
threshold_nb_sites = 20
selected_countries <- metadata |> 
  dplyr::count(country) |> 
  dplyr::filter(n > threshold_nb_sites) |> 
  dplyr::pull(country)



# (1) HUMAN FOOTPRINT 
human_footprint <- 
  plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                               X_new_data = X_pristine_conditions,
                               metadata,
                               save_name = "Human_footprint",
                               selected_countries,
                               plot_responders_on_map = F,
                               is_counterfactual = TRUE)


# (2) CONSERVATION LEGACY 
conservation_legacy <- 
  plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                               X_new_data = X_conservation_legacy,
                               metadata,
                               save_name = "Conservation_legacy",
                               selected_countries,
                               plot_responders_on_map = F,
                               is_counterfactual = TRUE)





#### Explore other changes ###

#(3) Change effectiveness only: from "out" to "full protection"
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_new_mpa,
                             metadata,
                             save_name = "full_protection_potential_gains_of_out_sites",
                             selected_countries,
                             plot_responders_on_map = T,
                             is_counterfactual = F)

#(4) Change fishing pressure only: set the number of fishing vessel to the minimum known.
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_no_vessels,
                             metadata,
                             save_name = "fishing_pressure_footprint",
                             selected_countries,
                             is_counterfactual = TRUE)

#(5) Conservation potential: unprotected sites are placed in reserves, without fishing pressure.
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_new_mpa_no_vessels,
                             metadata,
                             save_name = "Conservation_potential_gains-no-vessels-full-MPA",
                             selected_countries,
                             plot_responders_on_map = T,
                             is_counterfactual = F)

#(6) Change human "pollution": minimal gravity 
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_low_gravity,
                             metadata,
                             save_name = "gravity_footprint",
                             selected_countries,
                             is_counterfactual = TRUE)

#(7) Change human pressure: minimal gravity and maximal neartt (travel time to the nearest market)
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_no_human,
                             metadata,
                             save_name = "Population_footprint-gravity&neartt",
                             selected_countries,
                             is_counterfactual = TRUE)



#(8) Change effectiveness only: from "out" to "restricted protection"
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_new_medium_mpa,
                             metadata,
                             save_name = "restricted_MPA_potential_gains",
                             selected_countries,
                             is_counterfactual = F)




# ##----------------------------- Plot panel -----------------------------------
# set_ids = new_mpa_no_vessels #look only at currently unprotected sites
# # set_ids = new_pristine
# # set_ids = new_no_human
# 
# 
# protected <- plot_conterfactual_scenarios(path, model_name, concatenate_chains,
#                                           X_new_data = X_new_mpa_no_vessels,
#                                           metadata, save_name = "",
#                                           selected_countries,
#                                           set_ids = set_ids)
# 
# no_human <- plot_conterfactual_scenarios(path, model_name, concatenate_chains,
#                                          X_new_data = X_no_human,
#                                          metadata, save_name = "",
#                                          selected_countries,
#                                          set_ids = set_ids)
# 
# pristine <- plot_conterfactual_scenarios(path, model_name, concatenate_chains,
#                                          X_new_data = X_new_pristine,
#                                          metadata, save_name = "",
#                                          selected_countries,
#                                          set_ids = set_ids)
# 
# 
# A <- pristine[[1]] + labs(title = "")+
#   theme(legend.position = "none")+
#   xlim(-8,8)+
#   annotate("text", x = -7.3, y = Inf, 
#            label = expression(bold("A)")~"Pristine conditions"), 
#            size = 6, vjust = 1, fontface = "bold")
# 
# 
# B <- protected[[2]] + labs(title = "") + 
#   theme(legend.position = "none")+
#   xlim(-8,8)+
#   annotate("text", x = -6.7, y = Inf, 
#            label = expression(bold("B)")~"Protected sites"), 
#            size = 6, vjust = 1, fontface = "bold")
# 
# C <- no_human[[2]] + labs(title = "") +
#   theme(legend.position = "none")+
#   xlim(-8,8)+
#   annotate("text", x = -5.7, y = Inf, 
#            label = expression(bold("C)")~"No human degradations"), 
#            size = 6, vjust = 1, fontface = "bold")
# 
# 
# # Get legend
# pristine_with_larger_legend <- pristine[[1]] +
#   theme(legend.text = element_text(size = 15),
#         legend.title =  element_text(size = 20),
#         legend.key.height = unit(1, 'cm'))  
# legend <- ggpubr::get_legend(pristine_with_larger_legend)
# 
# 
# 
# 
# panel <- (A + legend + plot_layout(widths = c(4, 1))) / 
#   (B | C) + 
#   # plot_annotation(tag_levels = "A") +
#   plot_layout(heights = c(1.6, 1)) &
#   theme(
#     legend.text = element_text(size = 13),
#     plot.title = element_text(size = 16),       
#     plot.tag = element_text(face = 'bold', size = 14)
#   )
# 
# # Save panel
# path_file <- here::here("figures","models","hmsc", "conterfactuals", gsub(".rds", "", model_name))    
# 
# ggsave(filename = file.path( path_file, paste0("Panel_conterfactuals_",model_name,".jpg")),
#        width = 20, height =15 )

folder_name <- gsub(".rds", "", model_name)
path_file <- here::here("figures","models","hmsc", "conterfactuals", folder_name)    


conserv_legacy <- conservation_legacy[[3]] |> 
  dplyr::select(id, contribution, conserv_legacy_change = raw_change_percent) |> 
  dplyr::group_by(contribution) |> 
  dplyr::summarise(mean_change_percent_conserv_legacy = mean(conserv_legacy_change),
                   median_change_percent_conserv_legacy = median(conserv_legacy_change),
                   sd_change_percent_conserv_legacy = sd(conserv_legacy_change))


hum_footprint <- human_footprint[[3]]|> 
  dplyr::select(id, contribution, hum_footprint_change = raw_change_percent) |> 
  dplyr::group_by(contribution) |> 
  dplyr::summarise(mean_change_percent_hum_footprint = mean(hum_footprint_change),
                   median_change_percent_hum_footprint = median(hum_footprint_change),
                   sd_change_percent_hum_footprint = sd(hum_footprint_change))


all_changes <- conservation_legacy[[3]] |> 
  dplyr::select(id, contribution, change = raw_change_percent) |> 
  dplyr::mutate(counterfactual = "conservation_legacy") |> 
  dplyr::bind_rows(
    human_footprint[[3]]|> 
      dplyr::select(id, contribution, change = raw_change_percent) |> 
      dplyr::mutate(counterfactual = "human_footprint")
  )

test_zero <- all_changes |> 
  dplyr::group_by(counterfactual, contribution) |> 
  dplyr::summarise(
    t_test = list(t.test(change, mu = 0)), .groups = "drop") |> 
  dplyr::mutate(
    t_results = purrr::map(t_test, broom::tidy))  |> 
  tidyr::unnest(t_results)

test_between_couterfactuals <- all_changes |> 
  dplyr::group_by(contribution) |> 
  dplyr::summarise(
    t_test = list(t.test(
      change ~ counterfactual,  # Test t entre les deux countrefactuels
      data = dplyr::cur_data()       # Applique à chaque contribution
    )), .groups = "drop")|> 
  dplyr::mutate(t_results = purrr::map(t_test, broom::tidy))  |> 
   tidyr::unnest(t_results)
 

## Loliplot
change_percent <- conserv_legacy |> 
  dplyr::left_join(hum_footprint) |>
  dplyr::mutate(contribution = forcats::fct_reorder(contribution,
                                                    median_change_percent_hum_footprint, .fun = max, .desc = TRUE)) |>
  tidyr::pivot_longer(
    cols = starts_with("median") | starts_with("sd"),
    names_to = c("stat", "type"),
    names_pattern = "(median|sd)_change_percent_(.*)",
    values_to = "value" ) |> 
  tidyr::pivot_wider(
    names_from = "stat",
    values_from = "value") 


ggplot(change_percent, aes(x = median, y = contribution, color = type)) +
  geom_segment(aes(x = 0, xend = median, y = contribution, yend = contribution), 
               size = 1) + 
  geom_point(size = 4) + 
  scale_x_continuous(
    breaks = c(-250,-200,-150,-100,-50,0,50),      
    labels = c("-250%","-200%","-150%","-100%","-50%","0","+50%"))+
  scale_color_manual(values = c("conserv_legacy" = "darkseagreen3", 
                                "hum_footprint" = "firebrick1"),
                     labels = c("conserv_legacy" = "Conservation legacy", 
                                "hum_footprint" = "Human footprint")) + 
  labs(
    x = "Contribution changes in counterfactual scenarios", 
    y = "Contributions",
    color = "Counterfactual") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

ggsave(width = 10, height = 8, filename = file.path(
  path_file,paste0("loliplot_human_footprint_conserv_legacy","_", folder_name,".jpg")))




## Boxplot with log scale
change_percent_log <- conserv_legacy |> 
  dplyr::left_join(hum_footprint) |>
  dplyr::mutate(contribution = forcats::fct_reorder(contribution,
                                                    median_change_percent_hum_footprint, .fun = max, .desc = TRUE)) |>
  tidyr::pivot_longer(
    cols = starts_with("median") | starts_with("sd"),
    names_to = c("stat", "type"),
    names_pattern = "(median|sd)_change_percent_(.*)",
    values_to = "value" ) |> 
  tidyr::pivot_wider(
    names_from = "stat",
    values_from = "value") |> 
  dplyr::mutate(median = dplyr::case_when( median > 1 ~ log10(median),
                                           median < -1 ~ -log10(-(median)),
                                           T ~ median),
                sd = dplyr::case_when( sd > 1 ~ log10(sd),
                                       sd < -1 ~ -log10(-(sd)),
                                       T ~ sd)) 


ggplot(change_percent_log, aes(x = median, y = contribution, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +  
  geom_errorbar(aes(xmin = median, xmax = median + sd), 
                position = position_dodge(width = 0.9), 
                width = 0.2, color = "black") +
  geom_errorbar(aes(xmin = median - sd, xmax = median), 
                position = position_dodge(width = 0.9), 
                width = 0.2, color = "black")+
  scale_fill_manual(values = c("conserv_legacy" = "darkseagreen3", 
                                "hum_footprint" = "firebrick"),
                     labels = c("conserv_legacy" = "Conservation legacy", 
                                "hum_footprint" = "Human footprint")) + 
  scale_x_continuous(
    breaks = c(-3, -2, -1, 0, 1, 2, 3),      
    labels = c("-1000%", "-100%", "-10%", "0", "+10%", "+100%", "+1000%"))+  
  # xlim(-3,3)+
  labs(
    x = "Contribution changes in counterfactual scenarios", 
    y = "Contribution",
    color = "Counterfactual") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

ggsave(width = 10, height = 8, filename = file.path(
  path_file,paste0("barplot_human_footprint_conserv_legacy","_", folder_name,".jpg")))

