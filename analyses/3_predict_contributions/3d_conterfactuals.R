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
model_name <- gsub("output_", "", list_files[2]) #choose the wanted file
concatenate_chains = F


##----------------------- Counterfactual scenarios -----------------------------

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
                             selected_countries)

#(2) Change fishing pressure only: set the number of fishing vessel to the minimum known.
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_no_vessels,
                             metadata,
                             save_name = "Low_fishing_pressure",
                             selected_countries)

#(3) Change fishing pressure (minimal fishing vessels) and high effectiveness MPA => real protection
protected <- plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                                          X_new_data = X_new_mpa_no_vessels,
                                          metadata,
                                          save_name = "Total_protection-MPA&vessels",
                                          selected_countries)

#(4) Change human "pollution": minimal gravity 
plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                             X_new_data = X_low_gravity,
                             metadata,
                             save_name = "Low_gravity",
                             selected_countries)

#(5) Change human pressure: minimal gravity and maximal neartt (travel time to the nearest market)
no_human <- plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                                         X_new_data = X_no_human,
                                         metadata,
                                         save_name = "No_human-gravity&neartt",
                                         selected_countries)

#(6) Pristine sites: total protection and no human around
pristine <- plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                                         X_new_data = X_new_pristine,
                                         metadata,
                                         save_name = "Pristine-lowgravity&neartt&MPA&vessels",
                                         selected_countries)


##----------------------------- Plot panel -----------------------------------
A <- pristine[[1]] + labs(title = "Pristine conditions")+
  xlim(-6,8)

B <- protected[[2]] + labs(title = "Protected sites") + 
  theme(legend.position = "none")+
  xlim(-6,8)

C <- no_human[[2]] + labs(title = "No human pollutions") +
  theme(legend.position = "none")+
  xlim(-6,8)

panel <- A / (B | C) + 
  plot_annotation(tag_levels = "A")+
  plot_layout(heights = c(3, 2), guides = "collect") &
  theme(
    legend.text =  element_text(size = 13),
    plot.title = element_text(size = 16),       
    plot.tag = element_text(face = 'bold', size = 14)  
  )


# Save panel
path_file <- here::here("figures","models","hmsc", "conterfactuals", gsub(".rds", "", model_name))    

ggsave(filename = file.path( path_file, paste0("Panel_conterfactuals_",model_name,".jpg")),
       width = 20, height =15 )

