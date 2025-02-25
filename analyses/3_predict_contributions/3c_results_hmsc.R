################################################################################'
##
##  
##
## 3d_results_hmsc.R
##
## 13/06/2024
##
## Ulysse Flandrin
##
################################################################################'

##---------------------------- cleaning memory ---------------------------------
rm(list=ls())

##-----------------------------Loading packages---------------------------------
pkgs <- c("here", "Hmsc", "coda", "ggmcmc", "jsonify", "sf")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))


source("R/HMSC_function.R")

##------------------------------- load data ------------------------------------
## Contributions to predict
load(here::here("data/derived_data/3_all_contributions_to_predict.Rdata"))

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

##----------------------------- import hmsc output ------------------------------
## List all files in the directory and choose the model

list_files <- list.files(file.path(path, "out_multi")) 
list_files
file_name <- gsub("output_", "", list_files[5]) #choose the wanted file
concatenate_chains = F

##----------------------------- Plot hmsc results ------------------------------

plot_hmsc_result(metadata = metadata_sites,
                 file_name = file_name,
                 path = path,
                 concatenate_chains = concatenate_chains,
                 plot_convergence = T,
                 plot_explanatory_power = T,
                 plot_variance_partitioning = T,
                 plot_residual_associations = T,
                 plot_estimates = T,
                 plot_partial_graph = T,
                 check_residuals = T,
                 check_spatial_autocorrelation = F,
                 latent_factors = T,
                 drivers_to_plot =  list(
                   c("protection_statusfull", "Fishing_vessel_density","Gravity","GDP"),
                   c("protection_statusfull", "protection_statusrestricted",
                     "Fishing_vessel_density"),
                   # c("protection_statushigh","Gravity", "Fishing_vessel_density"),
                   # c("protection_statushigh", "protection_statusmedium",
                   #   "protection_statuslow", "Fishing_vessel_density"),
                   # 
                   c("protection_statusfull", "protection_statusrestricted",
                     "Marine_ecosystem_dependency", "Natural_ressource_rent",
                     "HDI","Travel_time"),
                   c("SST_5_years", "Coral_RLS", "Chlorophyll_5_years",
                     "DHW_quantile95_5_years")
                 )
)


### Run prediction ###
folder_name <- gsub(".rds", "", file_name)

make_crossval_prediction_hmsc(path = path,
                              folder_name = folder_name,
                              concatenate_chains = concatenate_chains,
                              conditional_prediction = T,
                              mcmcStep_conditional = 100, 
                              marginal_responses = colnames(observations_final)
)


### Plot predictive power ###

plot_predictive_power(path = path,
                      folder_name = folder_name)




