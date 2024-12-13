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
model_name <- gsub("output_", "", list_files[15]) #choose the wanted file
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
new_pristine <- rownames(X_pristine_conditions |> dplyr::filter(protection_status == "out"))

X_pristine_conditions[new_pristine, "protection_status"] <- as.factor("full")
X_pristine_conditions[new_pristine, "n_fishing_vessels"] <- min(X_pristine_conditions$n_fishing_vessels)
X_pristine_conditions[new_pristine, "gravity"] <- min(X_pristine_conditions$gravity)
X_pristine_conditions[new_pristine, "neartt"] <- min(X_pristine_conditions$neartt)


# (2) CONSERVATION LEGACY - Change protected sites to: out MPA, mean fishing vessels of "out" sites

#fishing vessels
mean_fishing_out <- mean(X[new_pristine, "n_fishing_vessels"])
mean_country_out <- X[new_pristine, ] |> 
  dplyr::group_by(country) |> 
  dplyr::summarise(fishing_out = mean(n_fishing_vessels))
missing_rows <- data.frame(
  country = unique(X$country)[!unique(X$country) %in% mean_country_out$country],
  fishing_out = mean_fishing_out)
mean_country_out <- rbind(mean_country_out, missing_rows)


## Convervation legacy of full MPA only
X_conservation_legacy_full <- X
new_conserv_legacy_full_mpa <- rownames(X_conservation_legacy_full |> dplyr::filter(protection_status == "full"))

X_conservation_legacy_full[new_conserv_legacy_full_mpa, "protection_status"] <- as.factor("out")

X_conservation_legacy_full <- X_conservation_legacy_full |> 
  tibble::rownames_to_column("id") |> 
  dplyr::left_join(mean_country_out) |> 
  dplyr::mutate(n_fishing_vessels = dplyr::case_when(
    id %in% new_conserv_legacy_full_mpa ~ fishing_out,
    TRUE ~ n_fishing_vessels)) |> 
  tibble::column_to_rownames("id") |> 
  dplyr::select(-fishing_out)



## Convervation legacy of all MPAs
X_conservation_legacy_all <- X
new_conserv_legacy_all_mpa <- rownames(X_conservation_legacy_all |> dplyr::filter(protection_status != "out"))

X_conservation_legacy_all[new_conserv_legacy_all_mpa, "protection_status"] <- as.factor("out")

X_conservation_legacy_all <- X_conservation_legacy_all |> 
  tibble::rownames_to_column("id") |> 
  dplyr::left_join(mean_country_out) |> 
  dplyr::mutate(n_fishing_vessels = dplyr::case_when(
    id %in% new_conserv_legacy_all_mpa ~ fishing_out,
    TRUE ~ n_fishing_vessels)) |> 
  tibble::column_to_rownames("id") |> 
  dplyr::select(-fishing_out)

# X_conservation_legacy[new_conserv_legacy, "n_fishing_vessels"] <- mean_fishing_out





#### Explore other changes ####

#(3) Change MPA only: from "out" to "full protection"
X_new_mpa <- X
new_mpa <- rownames(X_new_mpa |> dplyr::filter(protection_status == "out"))
X_new_mpa[new_mpa, "protection_status"] <- as.factor("full")

#(4) Change fishing pressure only: set the number of fishing vessel to the minimum known.
X_no_vessels <- X
new_vessels <- rownames(
  X_no_vessels[X_no_vessels$n_fishing_vessels != min(X_no_vessels$n_fishing_vessels),]
)
X_no_vessels[new_vessels, "n_fishing_vessels"] <- min(X_no_vessels$n_fishing_vessels)

#(5) Conservation potential: unprotected sites are placed in reserves, without fishing pressure.
X_new_mpa_no_vessels <- X
new_mpa_no_vessels <- unique(c(new_mpa, new_vessels))
X_new_mpa_no_vessels[new_mpa_no_vessels, "protection_status"] <- as.factor("full")
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
new_medium_mpa <- rownames(X_new_medium_mpa |> dplyr::filter(protection_status == "out"))
X_new_medium_mpa[new_medium_mpa, "protection_status"] <- as.factor("restricted")


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


conservation_legacy_all_mpa <- 
  plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                               X_new_data = X_conservation_legacy_all,
                               metadata,
                               save_name = "Conservation_legacy_all_MPA",
                               selected_countries,
                               plot_responders_on_map = F,
                               is_counterfactual = TRUE,
                               set_order_boxplot = set_order_boxplot)


conservation_legacy_full_mpa <- 
  plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                               X_new_data = X_conservation_legacy_full,
                               metadata,
                               save_name = "Conservation_legacy_full_MPA",
                               selected_countries,
                               plot_responders_on_map = F,
                               is_counterfactual = TRUE,
                               set_order_boxplot = set_order_boxplot)




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




##----------------------------- Plot panel -----------------------------------

folder_name <- gsub(".rds", "", model_name)
path_file <- here::here("figures","models","hmsc", "conterfactuals", folder_name)    

## choose the type of conservation legacy ###########################################"
conservation_legacy <- conservation_legacy_full_mpa
# conservation_legacy <- conservation_legacy_all_mpa


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

# # Parametric tests      
# test_zero <- all_changes |> 
#   dplyr::group_by(counterfactual, contribution) |> 
#   dplyr::summarise(
#     t_test = list(t.test(change, mu = 0)), .groups = "drop") |> 
#   dplyr::mutate(
#     t_results = purrr::map(t_test, broom::tidy))  |> 
#   tidyr::unnest(t_results) |> 
#   dplyr::mutate(different_from_zero = ifelse(p.value < 0.05, "yes", "no")) |> 
#   dplyr::select(contribution, counterfactual, different_from_zero)
# 
# test_between_couterfactuals <- all_changes |> 
#   dplyr::mutate(change = dplyr::case_when(
#     counterfactual == "conservation_legacy" ~ - change, # take the opposite of conservation gains to compare with human footprint
#     TRUE ~ change)) |> 
#   dplyr::group_by(contribution) |> 
#   dplyr::summarise(
#     welch_t = list(t.test(change ~ counterfactual, var.equal = FALSE))) |> 
#   dplyr::mutate(
#     welch_t_results = purrr::map(welch_t, broom::tidy)) |> 
#   tidyr::unnest(welch_t_results) |> 
#   dplyr::mutate(footprint_different_from_gain = ifelse(p.value < 0.05, "yes", "no"))|> 
#   dplyr::select(contribution, footprint_different_from_gain)

# Non-Parametric tests      
test_zero <- all_changes |> 
  dplyr::group_by(counterfactual, contribution) |> 
  dplyr::summarise(
    wilcox_test = list(wilcox.test(change, mu = 0)), .groups = "drop") |> 
  dplyr::mutate(
    wilcox_results = purrr::map(wilcox_test, broom::tidy))  |> 
  tidyr::unnest(wilcox_results) |> 
  dplyr::mutate(different_from_zero = ifelse(p.value < 0.05, "yes", "no")) |> 
  dplyr::select(contribution, counterfactual, different_from_zero)

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
  dplyr::select(contribution, footprint_different_from_gain)

## Loliplot
change_percent <- conserv_legacy |> 
  dplyr::left_join(hum_footprint) |>
  dplyr::mutate(contribution =forcats::fct_reorder(contribution,
     median_change_percent_hum_footprint, .fun = max, .desc = TRUE)) |>
  # Compare medians
  dplyr::mutate(bigger_effect = 
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
                                       "yes", "no"))


## Broke  X-axis
x_break = -60
xticks <- c(-120, -60,-40,-20, 0, 20)

# Transform the data
trans_x <- function(x) {
  pmax(x, x_break) + 0.12 * pmin(x + abs(x_break), 0)
}
trans_x(xticks)
change_percent$median_t <- trans_x(change_percent$median)


ggplot(change_percent, aes(x = median_t, y = contribution, color = counterfactual)) +
  geom_segment(aes(x = 0, xend = median_t, y = contribution, yend = contribution,
                   alpha = different_from_zero), 
               size = 1) + 
  geom_point(aes(alpha = different_from_zero, size = bigger_effect)) + 
  scale_alpha_manual(values = c("yes" = 1, "no" = 0.3)) +
  scale_size_manual(values = c("yes" = 5, "no" = 2)) +
  geom_rect(aes(xmin = trans_x(x_break-15), xmax = x_break-4, ymin = -Inf, ymax = Inf), 
            fill = "white", color = NA) +  # Mask broken area
  annotate("text", x = x_break-3, y = max(as.numeric(factor(change_percent$contribution))), 
           label = "//", size = 7) + # Add a visual cue for the break
  scale_x_continuous(
    breaks = trans_x(xticks),
    labels = paste0(xticks, "%") )+
  scale_color_manual(values = c("conservation_legacy" = "darkseagreen3", 
                                "human_footprint" = "firebrick3")) + 
  labs(
    x = "Contribution changes in counterfactual scenarios", 
    y = "Contributions",
    color = "Counterfactuals") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0,1,0,0, unit = "cm"),
        legend.title = element_text(face="bold", size = 11)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE),
         alpha = "none",#guide_legend(nrow = 2, byrow = TRUE),
         size = "none")

ggsave(width = 9, height = 8, filename = file.path(
  path_file,paste0("loliplot_human_footprint_conserv_legacy_full_MPA","_", folder_name,".jpg")))


# ggplot(change_percent, aes(x = median, y = contribution, color = counterfactual)) +
#   geom_segment(aes(x = 0, xend = median, y = contribution, yend = contribution,
#                    alpha = different_from_zero),
#                size = 1) +
#   geom_point(aes(alpha = different_from_zero, size = bigger_effect)) +
#   scale_alpha_manual(values = c("yes" = 1, "no" = 0.3)) +
#   scale_size_manual(values = c("yes" = 5, "no" = 3)) +
#   # Stars for significant differences
#   # geom_text(
#   #   data = change_percent |> dplyr::filter(footprint_different_from_gain == "yes") |>
#   #     dplyr::distinct(contribution, .keep_all = TRUE),
#   #   aes(x = 25, y = contribution, label = "*"),
#   #   color = "black", hjust = 0, size = 5
#   # ) +
#   scale_x_continuous(
#     breaks = c(-250,-200,-150,-100,-50,0,50),
#     labels = c("-250%","-200%","-150%","-100%","-50%","0","+50%"))+
#   scale_color_manual(values = c("conservation_legacy" = "darkseagreen3",
#                                 "human_footprint" = "firebrick3")) +
#   labs(
#     x = "Contribution changes in counterfactual scenarios",
#     y = "Contributions",
#     color = "Counterfactual") +
#   theme_minimal() +
#   theme(axis.text.y = element_text(size = 10),
#         legend.position = "bottom",
#         legend.key.spacing.y = unit(0.1, "cm"),
#         legend.margin = margin(0,1,0,0, unit = "cm"),
#         legend.title = element_text(face="bold", size = 11)) +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE),
#          alpha = guide_legend(nrow = 2, byrow = TRUE),
#          size = "none")



## Barplot with log scale
change_percent_log <- conserv_legacy |> 
  dplyr::left_join(hum_footprint) |>
  dplyr::mutate(contribution = forcats::fct_reorder(contribution,
                                                    median_change_percent_hum_footprint, .fun = max, .desc = TRUE)) |>
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
  dplyr::mutate(median = dplyr::case_when( median > 1 ~ log10(median),
                                           median < -1 ~ -log10(-(median)),
                                           T ~ median),
                sd = dplyr::case_when( sd > 1 ~ log10(sd),
                                       sd < -1 ~ -log10(-(sd)),
                                       T ~ sd)) 


ggplot(change_percent_log, aes(x = median, y = contribution, fill = counterfactual)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +  
  geom_errorbar(aes(xmin = median, xmax = median + sd), 
                position = position_dodge(width = 0.9), 
                width = 0.2, color = "black") +
  geom_errorbar(aes(xmin = median - sd, xmax = median), 
                position = position_dodge(width = 0.9), 
                width = 0.2, color = "black")+
  scale_fill_manual(values = c("conservation_legacy" = "darkseagreen3", 
                                "human_footprint" = "firebrick3")) +
  scale_x_continuous(
    breaks = c(-3, -2, -1, 0, 1, 2, 3),      
    labels = c("-1000%", "-100%", "-10%", "0", "+10%", "+100%", "+1000%"))+  
  labs(
    x = "Contribution changes in counterfactual scenarios", 
    y = "Contribution"
    ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

ggsave(width = 10, height = 8, filename = file.path(
  path_file,paste0("barplot_human_footprint_conserv_legacy_full_MPA","_", folder_name,".jpg")))



## Boxplot 
all_changes_log_transformed <- all_changes|> 
  dplyr::mutate(change = dplyr::case_when(
    change > 1 ~ log10(change),
    change < -1 ~ -log10(-(change)),
    T ~ change
  )) |> 
  dplyr::mutate(contribution = reorder(contribution, change,
                                       FUN = median))
ggplot(all_changes) +
  geom_hline(yintercept = 0, color = "grey", size = 1) +
  aes(x= change, y= contribution, fill = counterfactual)+
  scale_fill_manual(values = c("conservation_legacy" = "darkseagreen3", 
                                "human_footprint" = "firebrick3")) +
  # geom_violin(trim = FALSE, position = position_dodge(width =1), alpha = 0.7) +
  geom_boxplot(alpha = 0.7, outliers = F) +

  ylab("") + xlab("Contributions change in counterfactual scenarios") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 15, margin = margin(r = 20)),
        panel.spacing = unit(0.3, "lines"),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 0, hjust = 0.5,size = 13))

ggsave(width = 10, height = 8, filename = file.path(
  path_file,paste0("boxplot_human_footprint_conserv_legacy_full_MPA","_", folder_name,".jpg")))






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
