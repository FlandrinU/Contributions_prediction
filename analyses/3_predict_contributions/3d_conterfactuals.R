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
model_name <- gsub("output_", "", list_files[3]) #choose the wanted file
concatenate_chains = F


##------------- New conditions in counterfactual scenarios ---------------------

### Initial conditions ###
X <- X_data_site
metadata = metadata_sites
# X <- X_data


#### Change initial conditions ###
summary(X)

#(1) Change effectiveness only: from "out" to "high protection"
X_new_mpa <- X
new_mpa <- rownames(X_new_mpa |> dplyr::filter(effectiveness == "out"))
X_new_mpa[new_mpa, "effectiveness"] <- as.factor("High")

#(2) Change fishing pressure only: set the number of fishing vessel to the minimum known.
X_no_vessels <- X
new_vessels <- rownames(
  X_no_vessels[X_no_vessels$n_fishing_vessels != min(X_no_vessels$n_fishing_vessels),]
)
X_no_vessels[new_vessels, "n_fishing_vessels"] <- min(X_no_vessels$n_fishing_vessels)

#(3) Change fishing pressure (minimual fishing vessels) and high effectiveness MPA => real protection
X_new_mpa_no_vessels <- X
new_mpa_no_vessels <- unique(c(new_mpa, new_vessels))
X_new_mpa_no_vessels[new_mpa_no_vessels, "effectiveness"] <- as.factor("High")
X_new_mpa_no_vessels[new_mpa_no_vessels, "n_fishing_vessels"] <- min(X_new_mpa_no_vessels$n_fishing_vessels)

#(4) Change human "pollution": minimal gravity 
X_low_gravity <- X
new_low_gravity <- rownames(
  X_low_gravity[X_low_gravity$gravtot2 != min(X_low_gravity$gravtot2) ,])
X_low_gravity[new_low_gravity, "gravtot2"] <- min(X_low_gravity$gravtot2)

#(5) Change human pressure: minimal gravity and maximal neartt (travel time to the nearest market)
X_no_human <- X_low_gravity
new_no_human <-  unique(rownames(X_no_human[X_no_human$neartt != max(X_no_human$neartt),]),
                   new_low_gravity)
X_no_human[new_no_human, "neartt"] <- max(X_no_human$neartt)

#(6) Pristine sites: total protection and no human around
X_new_pristine <- X_new_mpa_no_vessels
new_pristine <- rownames(X_new_pristine)
X_new_pristine[new_pristine, "gravtot2"] <- min(X_new_pristine$gravtot2)
X_new_pristine[new_pristine, "neartt"] <- max(X_new_pristine$neartt)

#(7) Change effectiveness only: from "out" to "low protection"
X_new_low_mpa <- X
new_low_mpa <- rownames(X_new_low_mpa |> dplyr::filter(effectiveness == "out"))
X_new_low_mpa[new_low_mpa, "effectiveness"] <- as.factor("Low")

#(8) Change effectiveness only: from "out" to "medium protection"
X_new_medium_mpa <- X
new_medium_mpa <- rownames(X_new_medium_mpa |> dplyr::filter(effectiveness == "out"))
X_new_medium_mpa[new_medium_mpa, "effectiveness"] <- as.factor("Medium")

#(9) Effect of current MPAs: from protected to "out"
X_remove_mpa <- X
new_rm_mpa <- rownames(X_remove_mpa |> dplyr::filter(effectiveness == "High"))
X_remove_mpa[new_rm_mpa, "effectiveness"] <- as.factor("out")

##----------------------------- Plot changes -----------------------------------

# Select countries we want to plot
threshold_nb_sites = 10
selected_countries <- metadata |> 
  dplyr::count(country) |> 
  dplyr::filter(n > threshold_nb_sites) |> 
  dplyr::pull(country)


#(1) Change effectiveness only: from "out" to "high protection"
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_new_mpa,
                             metadata,
                             save_name = "High_effectiveness",
                             selected_countries,
                             plot_responders_on_map = T)

#(2) Change fishing pressure only: set the number of fishing vessel to the minimum known.
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_no_vessels,
                             metadata,
                             save_name = "Low_fishing_pressure",
                             selected_countries)

#(3) Change fishing pressure (minimal fishing vessels) and high effectiveness MPA => real protection
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_new_mpa_no_vessels,
                             metadata,
                             save_name = "Total_protection-MPA&vessels",
                             selected_countries,
                             plot_responders_on_map = T)

#(4) Change human "pollution": minimal gravity 
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_low_gravity,
                             metadata,
                             save_name = "Low_gravity",
                             selected_countries)

#(5) Change human pressure: minimal gravity and maximal neartt (travel time to the nearest market)
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_no_human,
                             metadata,
                             save_name = "No_human-gravity&neartt",
                             selected_countries)

#(6) Pristine sites: total protection and no human around
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_new_pristine,
                             metadata,
                             save_name = "Pristine-lowgravity&neartt&MPA&vessels",
                             selected_countries,
                             plot_responders_on_map = T)

# #(7) Change effectiveness only: from "out" to "low protection" -> only 79 sites to fit the model
# plot_conterfactual_scenarios(path, model_name, concatenate_chains,
#                              X_new_data = X_new_low_mpa,
#                              metadata,
#                              save_name = "Low_effectiveness",
#                              selected_countries)

#(8) Change effectiveness only: from "out" to "Medium protection"
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_new_medium_mpa,
                             metadata,
                             save_name = "Medium_effectiveness",
                             selected_countries)


#(9) Effect of current MPAs: from protected to "out"
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_remove_mpa,
                             metadata,
                             save_name = "Remove_high_effectiveness_MPAs",
                             selected_countries)

##----------------------------- Plot panel -----------------------------------
set_ids = new_mpa_no_vessels #look only at currently unprotected sites
# set_ids = new_pristine
# set_ids = new_no_human

protected <- plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                                          X_new_data = X_new_mpa_no_vessels,
                                          metadata, save_name = "",
                                          selected_countries,
                                          set_ids = set_ids)

no_human <- plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                                         X_new_data = X_no_human,
                                         metadata, save_name = "",
                                         selected_countries,
                                         set_ids = set_ids)

pristine <- plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                                         X_new_data = X_new_pristine,
                                         metadata, save_name = "",
                                         selected_countries,
                                         set_ids = set_ids)


A <- pristine[[1]] + labs(title = "")+
  theme(legend.position = "none")+
  xlim(-8,8)+
  annotate("text", x = -7.3, y = Inf, 
           label = expression(bold("A)")~"Pristine conditions"), 
           size = 6, vjust = 1, fontface = "bold")


B <- protected[[2]] + labs(title = "") + 
  theme(legend.position = "none")+
  xlim(-8,8)+
  annotate("text", x = -6.7, y = Inf, 
           label = expression(bold("B)")~"Protected sites"), 
           size = 6, vjust = 1, fontface = "bold")

C <- no_human[[2]] + labs(title = "") +
  theme(legend.position = "none")+
  xlim(-8,8)+
  annotate("text", x = -5.7, y = Inf, 
           label = expression(bold("C)")~"No human degradations"), 
           size = 6, vjust = 1, fontface = "bold")


# Get legend
pristine_with_larger_legend <- pristine[[1]] +
  theme(legend.text = element_text(size = 15),
        legend.title =  element_text(size = 20),
        legend.key.height = unit(1, 'cm'))  
legend <- ggpubr::get_legend(pristine_with_larger_legend)




panel <- (A + legend + plot_layout(widths = c(4, 1))) / 
  (B | C) + 
  # plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1.6, 1)) &
  theme(
    legend.text = element_text(size = 13),
    plot.title = element_text(size = 16),       
    plot.tag = element_text(face = 'bold', size = 14)
  )

# Save panel
path_file <- here::here("figures","models","hmsc", "conterfactuals", gsub(".rds", "", model_name))    

ggsave(filename = file.path( path_file, paste0("Panel_conterfactuals_",model_name,".jpg")),
       width = 20, height =15 )

