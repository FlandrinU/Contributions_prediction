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
pkgs <- c("here", "Hmsc", "coda", "ggmcmc", "jsonify")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))


source("R/HMSC_function.R")

##------------------------------- load data ------------------------------------
load(here::here("data/derived_data/3_all_contributions_to_predict.Rdata"))
load(here::here("data/derived_data/3_all_covariates_to_predict.Rdata"))
response =  observations_final#[sample(1:nrow(observations_final),1000),] ####################### reduce data
covariates = covariates_final[rownames(response),]


# PATHS
save_init <- here::here("outputs/models/hmsc/init_multi/")
save_out <- here::here("outputs/models/hmsc/out_multi/")
localDir <- here::here("outputs/models/hmsc/multivariate")

##----------------------------- import hmsc output ------------------------------
## List all files in the directory and choose the model

list_files <- list.files(save_out) 
list_files
file_name <- gsub("output_", "", list_files[9]) #choose the wanted file


##----------------------------- Plot hmsc results ------------------------------

plot_hmsc_result(covariates = covariates,
                 response = response,
                 file_name = file_name,
                 save_init = save_init,
                 save_out = save_out,
                 localDir = localDir,
                 plot_convergence = T,
                 plot_explanatory_power = T,
                 plot_variance_partitioning =T,
                 plot_residual_associations = T,
                 plot_estimates = T,
                 check_residuals = T,
                 drivers_to_plot =  list(
                   c("gravtot2", "gdp", "neartt", "n_fishing_vessels"),
                   c("effectivenessHigh", "effectivenessMedium",
                     "effectivenessLow", "n_fishing_vessels"),
                   c("q95_5year_degree_heating_week", "median_5year_analysed_sst",
                     "median_5year_ph"),
                   c( "coral", "algae", "Terrestrial_Reef_Flat_500m", "depth")
                 )
)
