###############################################################################"
##
##  This script uses the HMSC output of the full model to build counterfactual
##   scenarios. By altering chosen covariates, it assesses the human footprint 
##   and the conservation legacy on reef fish contributions. This script produces
##   Fig. 3 and 4 of the paper flandrin et al.
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

## Site scale
load(here::here("data/derived_data/3_sites_covariates_to_predict.Rdata"))
X_data_site = covariates_site_final


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

model_name <-"FULL_model_SITE_SCALE_4_chains_1000_thin_200_samples.rds"
# model_name <- gsub("output_", "", list_files[5]) #choose the wanted file

concatenate_chains = F
##------------- New conditions in counterfactual scenarios ---------------------

### Initial conditions ###
X <- X_data_site
metadata = metadata_sites
# X <- X_data

## Identify "pristine" sites
low_grav_sites <- rownames(X[ X$Gravity == min(X$Gravity), ])
max_neartt_sites <- rownames(X[ X$Travel_time == min(X$Travel_time), ])

#### Change initial conditions ###
summary(X)


# (1) HUMAN FOOTPRINT - Change unprotected sites to: full protection, 
#  min fishing vessels, min Gravity and max Travel_time
X_pristine_conditions <- X
new_pristine <- rownames(X_pristine_conditions |> dplyr::filter(protection_status == "out"))
# new_pristine <- rownames(X_pristine_conditions |> dplyr::filter(protection_status == "full"))

X_pristine_conditions[new_pristine, "protection_status"] <- as.factor("full")
X_pristine_conditions[new_pristine, "Fishing_vessel_density"] <- min(X_pristine_conditions$Fishing_vessel_density)
X_pristine_conditions[new_pristine, "Gravity"] <- min(X_pristine_conditions$Gravity)
X_pristine_conditions[new_pristine, "Travel_time"] <- min(X_pristine_conditions$Travel_time)


# (2) CONSERVATION LEGACY - Change protected sites to: out MPA, mean fishing vessels of "out" sites

#fishing vessels
mean_fishing_out <- mean(X[new_pristine, "Fishing_vessel_density"])
mean_country_out <- X[new_pristine, ] |> 
  dplyr::group_by(country) |> 
  dplyr::summarise(fishing_out = mean(Fishing_vessel_density))
missing_rows <- data.frame(
  country = unique(X$country)[!unique(X$country) %in% mean_country_out$country],
  fishing_out = mean_fishing_out)
mean_country_out <- rbind(mean_country_out, missing_rows)


## Convervation legacy of full MPA only
X_conservation_legacy_full <- X
new_conserv_legacy_full_mpa <- rownames(X_conservation_legacy_full |> 
                                          dplyr::filter(protection_status == "full"))

X_conservation_legacy_full[new_conserv_legacy_full_mpa, "protection_status"] <- as.factor("out")

X_conservation_legacy_full <- X_conservation_legacy_full |> 
  tibble::rownames_to_column("id") |> 
  dplyr::left_join(mean_country_out) |> 
  dplyr::mutate(Fishing_vessel_density = dplyr::case_when(
    id %in% new_conserv_legacy_full_mpa ~ fishing_out,
    TRUE ~ Fishing_vessel_density)) |> 
  tibble::column_to_rownames("id") |> 
  dplyr::select(-fishing_out)



## Convervation legacy of all MPAs
X_conservation_legacy_all <- X
new_conserv_legacy_all_mpa <- rownames(X_conservation_legacy_all |> 
                                         dplyr::filter(protection_status != "out"))

X_conservation_legacy_all[new_conserv_legacy_all_mpa, "protection_status"] <- as.factor("out")

X_conservation_legacy_all <- X_conservation_legacy_all |> 
  tibble::rownames_to_column("id") |> 
  dplyr::left_join(mean_country_out) |> 
  dplyr::mutate(Fishing_vessel_density = dplyr::case_when(
    id %in% new_conserv_legacy_all_mpa ~ fishing_out,
    TRUE ~ Fishing_vessel_density)) |> 
  tibble::column_to_rownames("id") |> 
  dplyr::select(-fishing_out)



#### Explore other changes ####

#(3) Change MPA only: from "out" to "full protection"
X_new_mpa <- X
new_mpa <- rownames(X_new_mpa |> dplyr::filter(protection_status == "out"))
X_new_mpa[new_mpa, "protection_status"] <- as.factor("full")

#(4) Change fishing pressure only: set the number of fishing vessel to the minimum known.
X_no_vessels <- X
new_vessels <- rownames(
  X_no_vessels[X_no_vessels$Fishing_vessel_density != min(X_no_vessels$Fishing_vessel_density),]
)
X_no_vessels[new_vessels, "Fishing_vessel_density"] <- min(X_no_vessels$Fishing_vessel_density)

#(5) Conservation potential: unprotected sites are placed in reserves, without fishing pressure.
X_new_mpa_no_vessels <- X
new_mpa_no_vessels <- unique(c(new_mpa, new_vessels))
X_new_mpa_no_vessels[new_mpa_no_vessels, "protection_status"] <- as.factor("full")
X_new_mpa_no_vessels[new_mpa_no_vessels, "Fishing_vessel_density"] <- 
  min(X_new_mpa_no_vessels$Fishing_vessel_density)

#(6) Low Gravity
X_low_Gravity <- X
new_low_Gravity <- rownames(
  X_low_Gravity[X_low_Gravity$Gravity != min(X_low_Gravity$Gravity) ,])
X_low_Gravity[new_low_Gravity, "Gravity"] <- min(X_low_Gravity$Gravity)

#(7) No human population
X_no_human <- X_low_Gravity
new_no_human <-  unique(rownames(X_no_human[X_no_human$Travel_time != max(X_no_human$Travel_time),]),
                   new_low_Gravity)
X_no_human[new_no_human, "Travel_time"] <- max(X_no_human$Travel_time)


#(8) Change MPA only: from "out" to medium protection ("restricted")
X_new_medium_mpa <- X
new_medium_mpa <- rownames(X_new_medium_mpa |> dplyr::filter(protection_status == "out"))
X_new_medium_mpa[new_medium_mpa, "protection_status"] <- as.factor("restricted")


#(9) No gravity and Full protection
X_no_grav_protected <- X_low_Gravity
new_no_grav_protected <-  unique(rownames(X_no_grav_protected |> 
                                            dplyr::filter(protection_status == "out"))) #,
                        # new_low_Gravity)
X_no_grav_protected[new_no_grav_protected,  "protection_status"] <- as.factor("full")


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
order <- human_footprint[[3]]  |> 
  dplyr::mutate(contribution = reorder(contribution, raw_change_percent,
                                       FUN = median, decreasing = T )) 
set_order_boxplot <- levels(order$contribution)



conservation_legacy_full_mpa <- 
  plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                               X_new_data = X_conservation_legacy_full,
                               metadata,
                               save_name = "Conservation_legacy_full_MPA",
                               selected_countries,
                               plot_responders_on_map = F,
                               is_counterfactual = TRUE,
                               set_order_boxplot = set_order_boxplot)


conservation_legacy_all_mpa <- 
  plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                               X_new_data = X_conservation_legacy_all,
                               metadata,
                               save_name = "Conservation_legacy_all_MPA",
                               selected_countries,
                               plot_responders_on_map = F,
                               is_counterfactual = TRUE,
                               set_order_boxplot = set_order_boxplot)



## Save conterfactuals ##
conterfactuals <- list(human_footprint, conservation_legacy_all_mpa, conservation_legacy_full_mpa)
save(conterfactuals, file = paste0(here::here("figures","3_models","hmsc", "conterfactuals"),
                                   "/", gsub(".rds", "", model_name), 
                                   "/counterfactuals_to_plot.Rdata"))


#### Plot other changes ####

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

#(6) Change human "pollution": minimal Gravity 
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_low_Gravity,
                             metadata,
                             save_name = "Gravity_footprint",
                             selected_countries,
                             is_counterfactual = TRUE,
                             set_ids = new_mpa)

#(7) Change human pressure: minimal Gravity and maximal Travel_time (travel time to the nearest market)
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_no_human,
                             metadata,
                             save_name = "Population_footprint-Gravity&Travel_time",
                             selected_countries,
                             is_counterfactual = TRUE,
                             set_ids = new_mpa)



#(8) Change effectiveness only: from "out" to "restricted protection"
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_new_medium_mpa,
                             metadata,
                             save_name = "restricted_MPA_potential_gains",
                             selected_countries,
                             is_counterfactual = F)


#(9) Change human pressure: minimal Gravity and full protection
sensitivity_HF <- plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_no_grav_protected,
                             metadata,
                             save_name = "Sensitivity_analysis_HF_low_Gravity&full_protection",
                             selected_countries,
                             is_counterfactual = T,
                             set_ids = new_mpa)

##----------------------------- Plot Loliplot -----------------------------------
#Counterfactuals colors:
Cl_color = "darkseagreen3"
HF_color = "firebrick3"


folder_name <- gsub(".rds", "", model_name)
path_file <- here::here("figures","3_models","hmsc", "conterfactuals", folder_name)    

# ## load file
# load( file = paste0(here::here("figures","3_models","hmsc", "conterfactuals"),
#                     "/", gsub(".rds", "", model_name),
#                     "/counterfactuals_to_plot.Rdata"))
# human_footprint <- conterfactuals[[1]]
# conservation_legacy_all_mpa <- conterfactuals[[2]]
# conservation_legacy_full_mpa <- conterfactuals[[3]]

## choose the type of conservation legacy 
conservation_legacy <- conservation_legacy_full_mpa


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

## Statistical test
check_normality <- all_changes |> 
  dplyr::group_by(counterfactual, contribution) |> 
  dplyr::summarise(
    shapiro =  list(shapiro.test(change)))|> 
  dplyr::mutate(
    shapiro_results = purrr::map(shapiro, broom::tidy))  |> 
  tidyr::unnest(shapiro_results) |> 
  dplyr::mutate(normal_distrib = ifelse(p.value < 0.05, "no", "yes"))
      
check_homoscedasticity <- all_changes |> 
  dplyr::mutate(counterfactual = as.factor(counterfactual)) |>
  dplyr::group_by(contribution) |> 
  dplyr::summarise(
    levene = list(car::leveneTest(change ~ counterfactual))) |> 
  dplyr::mutate(
    levene_results = purrr::map(levene, broom::tidy)) |> 
  tidyr::unnest(levene_results) |> 
  dplyr::mutate(homoscedastic = ifelse(p.value < 0.05, "no", "yes"))


# Non-Parametric tests      
test_zero <- all_changes |> 
  dplyr::group_by(counterfactual, contribution) |> 
  dplyr::summarise(
    wilcox_test = list(wilcox.test(change, mu = 0)), .groups = "drop") |> 
  dplyr::mutate(
    wilcox_results = purrr::map(wilcox_test, broom::tidy))  |> 
  tidyr::unnest(wilcox_results) |> 
  dplyr::mutate(different_from_zero = ifelse(p.value < 0.05, "yes", "no")) |> 
  dplyr::select(contribution, counterfactual, different_from_zero)|> 
  # Arrange names
  dplyr::mutate(contribution = gsub("_", " ", contribution),
                contribution = dplyr::case_when(
                  contribution == "Iucn species richness" ~ "IUCN species richness",
                  TRUE ~ contribution
                ))

test_between_couterfactuals <- all_changes |> 
  dplyr::mutate(change = dplyr::case_when(
    counterfactual == "conservation_legacy" ~ -change, # take the opposite of conservation gains to compare with human footprint
    TRUE ~ change)) |> 
  dplyr::group_by(contribution) |> 
  dplyr::summarise(
    wilcox_test = list(wilcox.test(change ~ counterfactual, exact = FALSE)), .groups = "drop") |> 
  dplyr::mutate(
    wilcox_results = purrr::map(wilcox_test, broom::tidy)) |> 
  tidyr::unnest(wilcox_results) |> 
  dplyr::mutate(footprint_different_from_gain = ifelse(p.value < 0.05, "yes", "no")) |> 
  dplyr::select(contribution, footprint_different_from_gain)|> 
  # Arrange names
  dplyr::mutate(contribution = gsub("_", " ", contribution),
                contribution = dplyr::case_when(
                  contribution == "Iucn species richness" ~ "IUCN species richness",
                  TRUE ~ contribution
                ))




### PLOT LOLLIPLOT (FIGURE 3)  ####
grp_NN_NP <- data.frame(
  contribution = c("Actinopterygian richness","Functional distinctiveness",
                   "IUCN species richness" ,"Endemism",
                   "Evolutionary distinctiveness","Functional entropy",
                   "Phylogenetic entropy","Herbivores biomass",
                   "Invertivores biomass",  "Piscivores biomass",
                   "Trophic web robustness", "Mean trophic level",
                   
                   "Public attention", "Aesthetic",
                   "Available biomass", "Selenium",
                   "Zinc",   "Omega 3", "Calcium",
                   "Iron","Vitamin A", "Available biomass turnover"),
  group = c(rep("NN", 12), 
            rep("NC", 2),
            rep("NS", 8)))

## Arrange data
change_percent <- conserv_legacy |> 
  dplyr::left_join(hum_footprint) |> 
  # Arrange names
  dplyr::mutate(contribution = gsub("_", " ", contribution),
                contribution = dplyr::case_when(
                  contribution == "Iucn species richness" ~ "IUCN species richness",
                  TRUE ~ contribution
                )) |> 
  # Compare medians
  dplyr::mutate(contrib_order = median_change_percent_hum_footprint,
                bigger_effect = 
                  ifelse(abs(median_change_percent_conserv_legacy) > abs(median_change_percent_hum_footprint),
                         "conservation_legacy", "human_footprint")) |> 
  # Arrange table
  tidyr::pivot_longer(
    cols = starts_with("median") | starts_with("sd"),
    names_to = c("stat", "counterfactual"),
    names_pattern = "(median|sd)_change_percent_(.*)",
    values_to = "value" ) |> 
  tidyr::pivot_wider(
    names_from = "stat",
    values_from = "value") |> 
  dplyr::mutate(counterfactual = dplyr::recode(counterfactual, 
                                        "conserv_legacy" = "conservation_legacy", 
                                        "hum_footprint" = "human_footprint")) |> 
  dplyr::left_join(test_zero) |> 
  dplyr::left_join(test_between_couterfactuals) |> 
  # Keep significant differences between counterfactuals
  dplyr::mutate(bigger_effect = ifelse(footprint_different_from_gain == "yes" &
                                         bigger_effect == counterfactual,
                                       "yes", "no"))|> 
  # Reapply order to `contribution` after reshaping
  dplyr::left_join(grp_NN_NP) |> 
  dplyr::mutate(contribution = forcats::fct_reorder(contribution,
             contrib_order, .fun = max, .desc = TRUE))

 

## Broke  X-axis
x_break = -42
xticks <- c(signif(min(change_percent$median),2),-40,-20, 0, 15)

# Transform the data
trans_x <- function(x) {
  pmax(x, x_break) + 0.12 * pmin(x + abs(x_break), 0)
}
trans_x(xticks)
change_percent$median_t <- trans_x(change_percent$median)

##Plot
lolliplot <- ggplot(change_percent)+
  aes(x = median_t, y = contribution, color = counterfactual) +
  geom_segment(aes(x = 0, xend = median_t, y = contribution, yend = contribution,
                   alpha = different_from_zero), 
               size = 1) + 
  geom_point(aes(alpha = different_from_zero, size = bigger_effect)) + 
  scale_alpha_manual(values = c("yes" = 1, "no" = 0.3)) +
  scale_size_manual(values = c("yes" = 5, "no" = 2)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
  
  #broken axis
  geom_rect(aes(xmin = trans_x(x_break-16), xmax = x_break-4, ymin = -Inf, ymax = Inf), 
            fill = "white", color = NA) +  # Mask broken area
  annotate("text", x = x_break-3, y = max(as.numeric(factor(change_percent$contribution))), 
           label = "//", size = 7) + # Add a visual cue for the break
  
  scale_x_continuous(
    breaks = trans_x(xticks),
    labels = paste0(xticks, "%") )+
  scale_color_manual(values = c("conservation_legacy" = Cl_color, 
                                "human_footprint" = HF_color)) + 
  labs(
    x = "Contribution changes in counterfactual scenarios", 
    y = "", #"Contributions",
    color = "Counterfactuals") +
  theme_minimal() +
  
  # Custom labels with colored points
  scale_y_discrete(labels = function(labels) {
    sapply(seq_along(labels), function(i) {
      group_color <- unique(change_percent$group[change_percent$contribution == labels[i]])
      color <- ifelse(group_color == "NN", "forestgreen",
                      ifelse(group_color == "NC", "darkgoldenrod2", "dodgerblue3")) # Adjust colors
      paste0("<span style='color:", color, "; font-size: 13px;'>&#9679;</span> ", labels[i])
    })
  }) +
  
  theme(axis.text = ggtext::element_markdown(size = 12), #axis.text.y = element_text(size = 10), #, color = change_percent$group),
        axis.title = element_text(size = 13),
        legend.position = "bottom",
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0,1,0,0, unit = "cm"),
        legend.title = element_text(face="bold", size = 11)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE),
         alpha = "none",#guide_legend(nrow = 2, byrow = TRUE),
         size = "none")
lolliplot

# ggsave(width = 9, height = 8, plot = lolliplot, filename = file.path(
#   path_file,paste0("loliplot_human_footprint_conserv_legacy_full_MPA","_", folder_name,".jpg")))



## Test column proportion mitigated ####
change_prop <- change_percent |> 
  dplyr::group_by(contribution) |> 
  dplyr::mutate(mitigation = 
                  (median[counterfactual == "conservation_legacy"] / 
                  -median[counterfactual == "human_footprint"])*100 )


x_text_pos <- max(change_percent$median_t) * 1.3
existing_x_breaks <- trans_x(xticks)
existing_x_labels <- paste0(xticks, "%")

lolliplot_prop_column <- lolliplot +
  geom_tile(data = change_prop, aes(x = x_text_pos, y = contribution, fill = mitigation), 
            color = "black", width = 6, height = 1, alpha = 0.5) +  
    geom_text(data = change_prop, aes(x = x_text_pos, y = contribution, 
                                    label = paste0(round(mitigation, 0), " %")), 
            color = "black", size = 3) + 
  # scale_x_continuous(breaks = c( trans_x(xticks), x_text_pos),  
  #                    labels = c(paste0(xticks, "%"), 
  #                               "<b style='font-size:14px;'>Mitigation (%)</b>")) +
  scale_fill_gradient2(
    low = "grey40", mid = "grey50", high = "grey90", 
    midpoint = 0, 
    limits = c(min(change_prop$mitigation), 100),
    na.value = "grey95"
  ) +
  theme( #axis.text.x = ggtext::element_markdown(size = 10, angle = 45,  color = "black")
   axis.title.x = element_text(vjust = -5, hjust = 0.2))+
  guides(fill = "none")
lolliplot_prop_column

# ## Test Loliplot comparative plot ####
# 
# arrange_change <- change_percent |> 
#   dplyr::group_by(contribution) |> 
#   dplyr::mutate(x_start = ifelse(counterfactual == "conservation_legacy", median[counterfactual == "human_footprint"], 0),
#                 x_end = ifelse(counterfactual == "conservation_legacy", median[counterfactual == "human_footprint"]+
#                                  median[counterfactual == "conservation_legacy"], median[counterfactual == "human_footprint"])) |> 
#   dplyr::ungroup()|> 
#   dplyr::mutate(counterfactual = factor(counterfactual, levels = c("human_footprint", "conservation_legacy"))) |> 
#   dplyr::arrange(counterfactual)
# 
# x_break = -42
# xticks <- c(-150, -140,-40,-20, 0, 15)
# trans_x_linear <- function(x) {
#   x + (x < x_break)*90
# }
# trans_x_linear(xticks)
# arrange_change$x_start <- trans_x_linear(arrange_change$x_start)
# arrange_change$x_end <- trans_x_linear(arrange_change$x_end)
# 
# ggplot(arrange_change) +
#   aes(x = x_end, y = contribution, color = counterfactual, linetype = counterfactual) +
#   geom_segment(aes(x = x_start, xend = x_end, 
#                    y = as.numeric(contribution) + ifelse(counterfactual == "conservation_legacy", 0, -0), 
#                    yend = as.numeric(contribution) + ifelse(counterfactual == "conservation_legacy", 0, -0),
#                    alpha = different_from_zero), 
#                size = 1) + 
#   geom_point(aes(y = as.numeric(contribution) + ifelse(counterfactual == "conservation_legacy", 0, -0),
#                  alpha = different_from_zero, size = bigger_effect, shape = counterfactual)) + 
#   scale_alpha_manual(values = c("yes" = 1, "no" = 1)) +
#   scale_size_manual(values = c("yes" = 5, "no" = 5)) +
#   scale_linetype_manual(values= c(1,1))+
#   scale_y_continuous(breaks = seq_along(unique(arrange_change$contribution)), 
#                      labels = unique(arrange_change$contribution))+
#   
#   geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
#   
#   #broken axis
#   geom_rect(aes(xmin =x_break-2, xmax = x_break-4, ymin = -Inf, ymax = Inf),
#             fill = "white", color = NA) +  # Mask broken area
#   annotate("text", x = x_break-3, y = max(as.numeric(factor(change_percent$contribution))),
#            label = "//", size = 7) + # Add a visual cue for the break
#   
#   scale_x_continuous(
#     breaks = trans_x_linear(xticks),
#     labels = paste0(xticks, "%") )+
#   scale_color_manual(values = c("human_footprint" = HF_color, 
#                                 "conservation_legacy" = Cl_color )) + 
#   labs(
#     x = "Contribution changes in counterfactual scenarios", 
#     y = "", #"Contributions",
#     color = "Counterfactuals") +
#   theme_minimal() +
#   
#   # Custom labels with colored points
#   scale_y_discrete(labels = function(labels) {
#     sapply(seq_along(labels), function(i) {
#       group_color <- unique(change_percent$group[change_percent$contribution == labels[i]])
#       color <- ifelse(group_color == "NN", "forestgreen",
#                       ifelse(group_color == "NC", "darkgoldenrod2", "dodgerblue3")) # Adjust colors
#       paste0("<span style='color:", color, "; font-size: 13px;'>&#9679;</span> ", labels[i])
#     })
#   }) +
#   
#   theme(axis.text.y = ggtext::element_markdown(size = 10), #axis.text.y = element_text(size = 10), #, color = change_percent$group),
#         legend.position = "bottom",
#         legend.key.spacing.y = unit(0.1, "cm"),
#         legend.margin = margin(0,1,0,0, unit = "cm"),
#         legend.title = element_text(face="bold", size = 11)) +
#   guides(color = guide_legend(nrow = 1, byrow = TRUE),
#          alpha = "none",#guide_legend(nrow = 2, byrow = TRUE),
#          size = "none")
# 
# 
# ggsave(width = 9, height = 8, filename = file.path(
#   path_file,paste0("loliplot_test","CL_according_to_HF",".jpg")))

##----------------------------- Plot  Mean of loliplot -----------------------------------
grp_NN_NP <- data.frame(
  contribution = c("Actinopterygian richness","Functional distinctiveness",
                   "IUCN species richness" ,"Endemism",
                   "Evolutionary distinctiveness","Functional entropy",
                   "Phylogenetic entropy","Herbivores biomass",
                   "Invertivores biomass",  "Piscivores biomass",
                   "Trophic web robustness", "Mean trophic level",
                   
                   "Public attention", "Aesthetic",
                   "Available biomass", "Selenium",
                   "Zinc",   "Omega 3", "Calcium",
                   "Iron","Vitamin A", "Available biomass turnover"),
  group = c(rep("NN", 12), 
            rep("NC", 2),
            rep("NS", 8))) |> 
  dplyr::right_join(dplyr::select(change_percent, -group))

mean_NN_NP <- grp_NN_NP |> 
  dplyr::group_by(group, counterfactual) |> 
  dplyr::summarise(mean = mean(median)) |> 
  dplyr::mutate(group = factor(group, levels = c("NC", "NS", "NN")))

get_val <- function(grp, type){
  dplyr::filter(mean_NN_NP, group == grp & counterfactual == type)$mean
}

mean_plot <- ggplot(mean_NN_NP)+
  aes(x = mean , y = group, fill = group) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.5) +
  scale_fill_manual(values = c("NS" = "dodgerblue3", 
                               "NC" = "darkgoldenrod2", 
                               "NN" = "forestgreen")) +
  theme_test() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size=22), 
    axis.text.x = element_text(size = 20),
    legend.position = "none",
    plot.margin = unit(c(0,1,0,0), "cm")
  ) +
  labs(x = "Mean changes (%)")+

  #Add labels
  geom_text(aes(x = 0, y = "NN"), label = "Nature-for-Nature", color = "white", 
            hjust = 1.65, vjust = 0.5, size = 8) + 
  geom_text(aes(x = 0, y = "NS"), label = "Nature-for-Society", color = "black", 
            hjust = 1.6, vjust = 0.5, size = 8)+
  geom_text(aes(x = 0, y = "NC"), label = "Nature-as-Culture", color = "black",
            hjust = 1.64, vjust = 0.5, size = 8)+
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.7)

mean_plot


mean_plot_with_title <-  mean_plot +
  geom_text(aes(x = min(mean_NN_NP$mean), y = 3, 
                label = "Human footprint"), 
            color = HF_color, size = 7,
            hjust = -0.4, vjust = -5.2, fontface = "bold") +
  geom_text(aes(x = max(mean_NN_NP$mean) * 0.6, y =3, 
                label = "Conservation\n legacy"),
            color = Cl_color, size = 7, 
            hjust = 0.45, vjust = -1.4, fontface = "bold") +
  # make some space for the title
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 2, r = 1, b = 0, l = 0, unit = "cm"))


ggsave(filename =  paste0(path_file,"/Mean_changes_in_NFF.jpg"),
       plot = mean_plot_with_title,  width = 7, height = 5)
# ggsave(filename =  paste0(path_file,"/Mean_changes_in_NFF_sensitivity_HF_only_grav&mpafull.jpg"),
#        plot = mean_plot_with_title,  width = 7, height = 5)



## Inser Mean inside lolliplot
inser <- cowplot::ggdraw() + 
  cowplot::draw_image(paste0(path_file,"/Mean_changes_in_NFF.jpg"))+
  theme(plot.margin =unit(c(0,0,0,0), 'cm'))

plot_merged <- lolliplot +
  theme(legend.position = "none")+
  annotation_custom(ggplotGrob(inser),
                     xmin = -50,
                     xmax = -10,
                     ymin = 1,
                     ymax = 11)
plot_merged
# ggsave(filename =  paste0(path_file,"/Lolliplot_with_insert.jpg"),
#        plot = plot_merged,  width = 9, height = 8)

## With proportion collumn
title_column <- ggplot() + theme_void() +
  theme(legend.position = "none") +
  geom_text(aes(x = 0, y = 0.1, label = "Restoration \n by MPAs"), size = 3.5,
           color = "black", hjust = 0.5, vjust = 0.6, fontface = "bold", angle = 50)

plot_merged_prop <- lolliplot_prop_column +
  theme(plot.margin =unit(c(0,0,1.3,0), 'cm'))+
  theme(legend.position = "none")+
  annotation_custom(ggplotGrob(inser),
                    xmin = -50,
                    xmax = -10,
                    ymin = 1,
                    ymax = 11)+
  annotation_custom(ggplotGrob(title_column),
                    xmin = 13,
                    xmax = 27,
                    ymin = -2.5,
                    ymax = 0)
plot_merged_prop
ggsave(filename =  paste0(path_file,"/Lolliplot_prop_mitigated.jpg"),
       plot = plot_merged_prop,  width = 9, height = 8)
# ggsave(filename =  paste0(path_file,"/Lolliplot_prop_mitigated_sensitivity_HF_only_grav&mpafull.jpg"),
#        plot = plot_merged_prop,  width = 9, height = 8)


##----------------------------- Plot panel PCA -----------------------------------
PCA_HF <- human_footprint[[1]] +
  labs(title = "Human footprint")+
  coord_cartesian(xlim= c(-7.9,8))+
  theme(legend.position = "none",
        plot.title = element_text(size = 20, colour = HF_color, face = "bold", hjust = 0.5),
        axis.title.x = element_text(vjust = 123, hjust = 0, size = 12),
        axis.title.y = element_text(vjust = -143, hjust = 0, angle = 90, size = 12))


PCA_CL <- conservation_legacy[[1]]+
  labs(title = "Conservation legacy")+
  coord_cartesian(xlim= c(-7.9,8))+
  theme(legend.position = "none",
        plot.title = element_text(size = 20, colour = Cl_color, face = "bold", hjust = 0.5),
        axis.title.x = element_text(vjust = 123, hjust = 0, size = 12),
        axis.title.y = element_text(vjust = -143, hjust = 0, angle = 90, size = 12))

legend <- ggpubr::get_legend(human_footprint[[1]] + 
   theme(legend.text = element_text(size = 15),
         legend.title = element_text(size = 16, face = "bold", 
                                     margin = margin(r=1, unit = "cm")),
         legend.spacing.x = unit(2, "cm"),
         legend.key.spacing.x = unit(.8, "cm") )+
   guides(
     fill = guide_legend(nrow =3, override.aes = list(shape = 21), order = 2),  
     shape = guide_legend(nrow = 2, override.aes = list(fill = NA), order = 1),
     color = "none"
   ))

  

panel_PCA <- (PCA_HF + PCA_CL) / ggpubr::as_ggplot(legend) +
  plot_layout(heights = c(8,2))

ggsave(filename =  paste0(path_file,"/panel_PCA_HF_CL.jpg"),
       plot = panel_PCA,  width = 20, height = 11)


##----------------------------- FINAL PANEL FIG 3 -----------------------------------
method <- cowplot::ggdraw() + 
  cowplot::draw_image(here::here("report/Current reefs.jpg"))+
  theme(plot.margin =unit(c(0,0,0,0), 'cm'))

# loliplot_inser <- cowplot::ggdraw() +
#   cowplot::draw_image(paste0(path_file,"/Lolliplot_with_insert.jpg"))+
#   theme(plot.margin =unit(c(0,0,0,0), 'cm'))
loliplot_inser <- cowplot::ggdraw() +
  cowplot::draw_image(paste0(path_file,"/Lolliplot_prop_mitigated.jpg"))+
  theme(plot.margin =unit(c(0,0,0,0), 'cm'))

PCA <- cowplot::ggdraw() + 
  cowplot::draw_image(paste0(path_file,"/panel_PCA_HF_CL.jpg"))+
  theme(plot.margin =unit(c(0,0,0,0), 'cm'))



final_panel <- cowplot::plot_grid(
  cowplot::plot_grid(
    method, 
    loliplot_inser, 
    rel_widths = c(1.1, 2), 
    labels = c("A)", "B)"), 
    label_size = 16,      
    label_fontface = "bold"),
  PCA,
  rel_heights = c(1.1, 1),
  ncol = 1,
  labels = c("", "C)"),  
  label_size = 16,
  label_fontface = "bold")
# final_panel

# ggsave(filename =  paste0(path_file,"/Final_panel_Figure_3.jpg"),
#        plot = final_panel,  width = 13, height = 15)
ggsave(filename =  paste0(path_file,"/Final_panel_Figure_3_mitigation_prop.jpg"),
       plot = final_panel,  width = 13, height = 15)



##----------------------------- DENSITY PLOTS -----------------------------------

grp <- data.frame(
  contribution = c("Actinopterygian_richness","Functional_distinctiveness",
                   "Iucn_species_richness" ,"Endemism",
                   "Evolutionary_distinctiveness","Functional_entropy",
                   "Phylogenetic_entropy","Herbivores_biomass",
                   "Invertivores_biomass",  "Piscivores_biomass",
                   "Trophic_web_robustness", "Mean_trophic_level",
                   
                   "Public_attention", "Aesthetic",
                   "Available_biomass", "Selenium",
                   "Zinc",   "Omega_3", "Calcium",
                   "Iron","Vitamin_A", "Available_biomass_turnover"),
  group = c(rep("NN", 12), 
            rep("NC", 2),
            rep("NS", 8)))


conserv_legacy <- conservation_legacy[[3]] |> 
  dplyr::select(id, contribution, country, changes = raw_change_values, raw_original_prediction) |> 
  dplyr::mutate(counterfactual = "CL")

hum_footprint <- human_footprint[[3]] |> 
  dplyr::select(id, contribution, country, changes = raw_change_values, raw_original_prediction) |> 
  dplyr::mutate(counterfactual = "HF")

changes_df <- hum_footprint |> 
  dplyr::bind_rows(conserv_legacy) |> 
  dplyr::mutate(
    changes = ifelse(contribution %in% c("Available_biomass", "Herbivores_biomass",
                                         "Invertivores_biomass", "Piscivores_biomass"),
                     changes*20/1000, changes),
    raw_original_prediction = ifelse(contribution %in% c("Available_biomass", "Herbivores_biomass",
                                                         "Invertivores_biomass", "Piscivores_biomass"),
                                     raw_original_prediction*20/1000, raw_original_prediction)) |> 
  dplyr::mutate(name = dplyr::recode(
    contribution,
    "Actinopterygian_richness" = "Actinopterygian richness (species/500m²)",
    "Functional_distinctiveness" = "Functional distinctiveness",
    "Iucn_species_richness" = "IUCN species richness (species/500m²)",
    "Endemism" = "Endemism",
    "Evolutionary_distinctiveness" = "Evolutionary distinctiveness",
    "Functional_entropy" = "Functional entropy",
    "Phylogenetic_entropy" = "Phylogenetic entropy",
    "Herbivores_biomass" = "Herbivores biomass (kg/ha)",
    "Invertivores_biomass" = "Invertivores biomass (kg/ha)",
    "Piscivores_biomass" = "Piscivores biomass (kg/ha)",
    "Trophic_web_robustness" = "Trophic web robustness",
    "Mean_trophic_level" = "Mean trophic level",
    "Public_attention" = "Public attention (abritrary unit)",
    "Aesthetic" = "Aesthetic (abritrary unit)",
    "Available_biomass" = "Available biomass (kg/ha)",
    "Selenium" = "Selenium (µg/100g of fish)",
    "Zinc" = "Zinc (mg/100g of fish)",
    "Omega_3" = "Omega 3 (g/100g of fish)",
    "Calcium" = "Calcium (mg/100g of fish)",
    "Iron" = "Iron (mg/100g of fish)",
    "Vitamin_A" = "Vitamin A (µg/100g of fish)",
    "Available_biomass_turnover" = "Available biomass turnover (proportion of biomass renewed per day)"
  )
  )

summary_changes <- changes_df |> 
  dplyr::group_by(contribution, counterfactual) |> 
  dplyr::summarise(median = quantile(changes, 0.5),
                quantile_5 = quantile(changes, 0.05),
                quantile_95 = quantile(changes, 0.95))
density_plot_fct <- 
  function(contribution_to_plot =  c("Iucn_species_richness","Actinopterygian_richness"),
           add_quantile = 0.5,
           log_transform = T,
           facet_ncol = 1,
           scales = "free",
           x_label = "",
           x_text_CL = 0,
           x_text_HF = 0,
           y_text = 0.02,
           n_round = 1,
           legend.pos = "none"
           ){
    
    data <- changes_df|>
      dplyr::filter(contribution  %in% contribution_to_plot)  |> 
      dplyr::left_join(grp) |> 
      dplyr::mutate(color = ifelse(group == "NN", "forestgreen",
                                    ifelse(group == "NC", "darkgoldenrod2",
                                           "dodgerblue3")),
                    name = paste0("<span style='color:",
                                  color,
                                  "; font-size: 35px;'>&#9679;</span> ",
                                  name))

    
  
    
    plot <- ggplot(data)+
      aes(x=changes, group=counterfactual, fill=counterfactual) +
      geom_density(alpha = 0.4)+
      # geom_density(aes(x=raw_original_prediction, group=counterfactual, fill=counterfactual),alpha = 0.1)+
      
      #Add medians
      geom_point(data = data |> 
                   dplyr::group_by(name, counterfactual) |> 
                   dplyr::summarise(median_value = median(changes, na.rm = TRUE),
                                    .groups = 'drop'), 
                 aes(x = median_value, y = 0, fill = counterfactual), 
                 size = 6, shape = 21, color = "black") +
      
      geom_text(data = data |> 
                  dplyr::group_by(name, counterfactual) |>
                  dplyr::summarise(median_value = median(changes, na.rm = TRUE),
                                   sd = sd(changes, na.rm = TRUE),
                                   .groups = 'drop')|>
                  dplyr::mutate(x_text = dplyr::case_when(
                    counterfactual == "CL" ~ x_text_CL, 
                    counterfactual == "HF" ~ x_text_HF   
                  )),
                aes(x = median_value + sign(median_value)*sd*x_text,
                    y = y_text, 
                    label = round(median_value, n_round), color = counterfactual),
                size = 8, vjust = -1,
                show.legend = FALSE) +
      
      scale_fill_manual(
        values = c("firebrick3","darkseagreen3"),
        labels = c("HF" = "Human footprint", "CL" = "Conservation legacy"),
        limits = c("HF", "CL"),
        name = "Contribution changes in counterfactuals:") +
      scale_color_manual( values = c(colorspace::darken("darkseagreen3", 0.3), 
                                     colorspace::darken("firebrick3", 0.1)) ) +

            geom_vline(xintercept = 0)+ # linetype = "dashed", alpha = 0.5,
      geom_hline(yintercept = 0)+
      
      hrbrthemes::theme_ipsum(grid = F, axis = "x", ticks = T,
                              base_size = 17,
                              strip_text_size = 20,
                              axis_title_size = 20) +
      xlab(x_label) + ylab("")+
      theme(legend.position=legend.pos, panel.spacing = unit(0.1, "lines"),
            legend.text = element_text(size=20),
            legend.title = element_text(size=20),
            # legend.key.spacing = unit(0.7, "cm"),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_text(color = "black"),
            plot.margin = margin(0, 0, 0, 0),
            strip.text = ggtext::element_markdown(size = 20))+
      facet_wrap(~name, ncol = facet_ncol, scales = scales)
    
    
    # Add observed value
    if(is.numeric(add_quantile )){
      plot <- plot +
        geom_vline(data = data |>
                     dplyr::group_by(name, counterfactual) |>
                     dplyr::summarise(quantile = quantile(raw_original_prediction, 
                                                          probs = add_quantile),
                                      .groups = 'drop') |>
                     dplyr::group_by(name) |>
                     dplyr::mutate(quantile = max(quantile)) ,
                   aes(xintercept = quantile),
                   linewidth = 3,color = "black", linetype = "81", alpha = 0.4)
    }
    
    # Log transform
    if(log_transform){
      breaks = c(seq(-10000, -2000, 1000), seq(-1000, -200, 100), 
                 seq(-100, 100, 10), 
                 seq(200, 1000, 100), seq(2000, 10000, 1000))
                       
      plot <- plot +
        scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                           breaks = breaks,
                           # minor_breaks = c(seq(-10000, -2000, 1000), seq(-1000, -200, 100),
                           #                  seq(-100, 100, 10),
                           #                  seq(200, 1000, 100), seq(2000, 10000, 1000)),
                           labels = ifelse((grepl("^1", as.character(abs(breaks))) |
                                              breaks == 0), as.character(breaks), "")
        )
    }
    
    plot
  }

unique(changes_df$contribution)


#Biomass
density_plot_fct(contribution_to_plot =  c( "Available_biomass", "Herbivores_biomass",
                                            "Invertivores_biomass", "Piscivores_biomass"),
                 add_quantile = 0.5,
                 log_transform = T,
                 facet_ncol = 1,
                 scales = "free_y",
                 y_text = 0.07,
                 x_label = "Biomass changes (kg/ha)")
ggsave(filename =  paste0(path_file,"/density_changes_biomass.jpg"), width = 8, height = 12)


# Aesthetic
density_plot_fct(contribution_to_plot =  c("Aesthetic"),
                 add_quantile = 0.5,
                 x_text_CL = 11,
                 y_text = 0.03,
                 log_transform = T)
ggsave(filename =  paste0(path_file,"/density_changes_aesth.jpg"), width = 8, height = 3)


# IUCN
density_plot_fct(contribution_to_plot =  c("Iucn_species_richness"),
                 add_quantile = 0.5,
                 log_transform = F,
                 n_round = 2,
                 x_text_CL = 0.1,
                 y_text = 0.5)
ggsave(filename =  paste0(path_file,"/density_changes_IUCN.jpg"), width = 8, height = 3)


# Richness
density_plot_fct(contribution_to_plot =  c("Actinopterygian_richness"),
                 add_quantile = 0.5,
                 log_transform = F,
                 x_text_CL = 5.5,
                 x_text_HF = 2.7,
                 n_round = 2)
ggsave(filename =  paste0(path_file,"/density_changes_richness.jpg"), width = 8, height = 3)


# Biomass turnover
density_plot_fct(contribution_to_plot =  c("Available_biomass_turnover"),
                 add_quantile = 0.5,
                 log_transform = F,
                 x_text_CL = 0.5,
                 x_text_HF = 0.4,
                 y_text = 3,
                 n_round = 2)
ggsave(filename =  paste0(path_file,"/density_changes_turnover.jpg"), width = 8, height = 3)


# Nutrients
density_plot_fct(contribution_to_plot =  c("Selenium",     
                                           "Zinc", "Omega_3", "Calcium", "Iron",                   
                                           "Vitamin_A"),
                 add_quantile = 0.5,
                 log_transform = F,
                 facet_ncol = 1,
                 scales = "free",
                 x_text_CL = 2,
                 x_text_HF = 2,
                 n_round = 2)
ggsave(filename =  paste0(path_file,"/density_changes_nutrients_with_median.jpg"), width = 8, height = 10)


density_plot_fct(contribution_to_plot =  c("Selenium",     
                                           "Zinc", "Omega_3", "Calcium", "Iron",                   
                                           "Vitamin_A"),
                 add_quantile = F,
                 log_transform = F,
                 facet_ncol = 2,
                 scales = "free",
                 x_text_CL = 0.5,
                 x_text_HF = 0.5,
                 n_round = 2)
ggsave(filename =  paste0(path_file,"/density_changes_nutrients_without_max.jpg"), width = 13, height = 10)



#Others
density_plot_fct(contribution_to_plot =  c( "Functional_distinctiveness", "Evolutionary_distinctiveness",
                                            "Endemism", "Functional_entropy", "Phylogenetic_entropy",
                                            "Trophic_web_robustness", "Mean_trophic_level",
                                            "Public_attention"),
                 add_quantile = F,
                 log_transform = F,
                 facet_ncol = 2,
                 scales = "free",
                 n_round = 2)

ggsave(filename =  paste0(path_file,"/density_changes_others.jpg"), width = 17, height = 13)


#Panel SUPP.
density_plot_fct(contribution_to_plot =  c( "Selenium",     
                                            "Zinc", "Omega_3", "Calcium", "Iron",                   
                                            "Vitamin_A",
                                            "Functional_distinctiveness", "Evolutionary_distinctiveness",
                                            "Endemism", "Functional_entropy", "Phylogenetic_entropy",
                                            "Trophic_web_robustness", "Mean_trophic_level",
                                            "Public_attention"),
                 add_quantile = F,
                 log_transform = F,
                 facet_ncol = 2,
                 scales = "free",
                 x_text_CL = 0.5,
                 x_text_HF = 0.5,
                 n_round = 2,
                 legend.pos = "bottom")

ggsave(filename =  paste0(path_file,"/Supp_Panel_density_plot.jpg"), width = 15, height = 15)


### PANEL DENSITY PLOT ####

biomass <- cowplot::ggdraw() + 
  cowplot::draw_image(paste0(path_file,"/density_changes_biomass.jpg"))+
  theme(plot.margin =unit(c(0,0,0,0.5), 'cm'))

iucn <- cowplot::ggdraw() + 
  cowplot::draw_image(paste0(path_file,"/density_changes_IUCN.jpg"))+
  theme(plot.margin =unit(c(0,0,0,0.5), 'cm'))

richness <- cowplot::ggdraw() + 
  cowplot::draw_image(paste0(path_file,"/density_changes_richness.jpg"))+
  theme(plot.margin =unit(c(0,0,0,0.5), 'cm'))

aesth <- cowplot::ggdraw() + 
  cowplot::draw_image(paste0(path_file,"/density_changes_aesth.jpg"))+
  theme(plot.margin =unit(c(0,0,0,0.5), 'cm'))

turnover <- cowplot::ggdraw() + 
  cowplot::draw_image(paste0(path_file,"/density_changes_turnover.jpg"))+
  theme(plot.margin =unit(c(0,0,0,0.5), 'cm'))

legend <- ggpubr::as_ggplot(
  ggpubr::get_legend(
    density_plot_fct(contribution_to_plot =  c("Aesthetic"),legend.pos = "bottom"))
)

final_panel <- cowplot::plot_grid(
  biomass,
  
  cowplot::plot_grid(
    iucn, 
    richness,
    aesth,
    turnover,
    ncol = 1, 
    labels = c("B)", "C)", "D)", "E)"), 
    label_size = 13,      
    label_fontface = "bold"),
  
  ncol = 2,
  labels = c("A)", ""),  
  label_size = 13,
  label_fontface = "bold") / legend +
  plot_layout(heights = c(10,1))

# final_panel

ggsave(filename =  paste0(path_file,"/Figure_4_density_plot.jpg"),
       plot = final_panel,  width = 16, height = 12.5)