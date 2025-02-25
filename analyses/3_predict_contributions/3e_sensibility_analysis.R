################################################################################'
##
##  
##
## 3e_sensibility_analysis.R
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
library(ggplot2)
library(patchwork)

##----------------------------- data and hmsc outputs ------------------------------

# Path to model folder
path = here::here("outputs/models/hmsc")
## List all files in the directory and choose the model
list_files <- list.files(file.path(path, "out_multi")) 
list_files

##----------------------------- functions ------------------------------

extract_estimates <- function(file_name = full_model_file){
  # PATHS
  # file_name = full_model_file
  
  save_init <- file.path(path, "init_multi/")
  save_out <- file.path(path, "out_multi/")
  localDir <- file.path(path, "multivariate/")
  save_name <- gsub(".rds", "", file_name)
  
  #Model design
  load(file = paste0(localDir, paste0("model_fit_", save_name, ".rds.Rdata")))
  
  #Initial hmsc object
  init_obj_rds <- readRDS(paste0(save_init, paste0("init_", save_name, ".rds")))
  init_obj <- jsonify::from_json(init_obj_rds)
  
  nSamples = init_obj[["samples"]]
  thin = init_obj[["thin"]]
  nChains = init_obj[["nChains"]]
  transient = init_obj[["transient"]]
  
  
  ## Import posterior probability
  all_chains <- readRDS(file = paste0(save_out, paste0("output_", file_name)))[[1]]
  importFromHPC <- jsonify::from_json(all_chains)
  chainList <- importFromHPC[1:nChains]
  
  ##Export and merge chains
  model_fit_mcmc <- Hmsc::importPosteriorFromHPC(model_fit,
                                                 chainList, 
                                                 nSamples, 
                                                 thin, 
                                                 transient)
  
  ## Estimates for each chains
  mpost <- Hmsc::convertToCodaObject(model_fit_mcmc)
  
  postBeta <- Hmsc::getPostEstimate(model_fit_mcmc, parName = "Beta")
  
  
  estimates <-  as.data.frame(postBeta[["mean"]]) 
  rownames(estimates) <- model_fit_mcmc[["covNames"]]
  
  support_estimates <- postBeta[["support"]]
  rownames(support_estimates) <- model_fit_mcmc[["covNames"]]
  support_estimates <- as.data.frame(support_estimates) |> 
    tibble::rownames_to_column("covariate") |> 
    tidyr::pivot_longer(cols = -covariate,
                        names_to = "response", values_to = "support") |> 
    dplyr::mutate(support = dplyr::case_when(support > 0.95 ~ 1, 
                                             support < 0.05 ~ 1,
                                             TRUE ~ 0 ))|> 
    tidyr::pivot_wider(names_from = response, values_from = support) |> 
    tibble::column_to_rownames("covariate") 
  
  toPlot <- support_estimates * estimates 
  
  toPlot <- toPlot[-1,]
  
  
  list(estimates = estimates, support_estimates = support_estimates, 
       toPlot = toPlot, n_sampled = nrow(model_fit_mcmc[["Y"]]))
}  # END of extract_estimates function

##------------------ Photoquadrat sensibility analysis ---------------------

#### TEST WITH MODEL WITHOUT PQ COVARIATES ####
full_model_file <- "FULL_model_SITE_SCALE_2_chains_1000_thin_200_samples.rds"
PQ_model_file <- "Without_PQ_model_SITE_SCALE_2_chains_1000_thin_200_samples.rds"



full_model_estimates <- extract_estimates(full_model_file)
PQ_model_estimates <- extract_estimates(PQ_model_file)



PQ_model_estimates_matrix <- PQ_model_estimates[["estimates"]] |> 
  tibble::rownames_to_column("covariate") |> 
  tidyr::pivot_longer(-covariate, names_to = "contribution", 
                      values_to = "effect_sizes_without_PQ")

full_model_estimates_matrix <- full_model_estimates[["estimates"]]|> 
  tibble::rownames_to_column("covariate") |> 
  tidyr::pivot_longer(-covariate, names_to = "contribution",
                      values_to = "effect_sizes_full_model")


compare_estimate <- PQ_model_estimates_matrix |> 
  dplyr::left_join(full_model_estimates_matrix) |> 
  dplyr::filter(covariate != "(Intercept)")


corr <- cor(compare_estimate$effect_sizes_without_PQ, compare_estimate$effect_sizes_full_model)

ggplot(compare_estimate, aes(x=effect_sizes_full_model, 
                             y=effect_sizes_without_PQ)) +
  geom_point(shape = 21, alpha = 0.5,
             fill = "grey50", color= "black")+
  theme_bw()+
  
  ylab("Effect sizes in model without Photo-quadrat covariates") +
  
  xlab(paste0("Effect sizes in full model (N_reef = ",
              full_model_estimates[["n_sampled"]], ")"))+
  ylab(paste0("Effect sizes in model without Photo-quadrat covariates (N_reef = ",
              PQ_model_estimates[["n_sampled"]], ")")) +
  
  guides(colour = guide_colourbar(title.position="top")) +
  theme(legend.position = "none",
        legend.key.size = unit(0.5, 'cm'),
        legend.direction = "vertical",
        legend.title = element_text( size = 11),
        legend.background = element_rect(fill='transparent'),
        axis.title=element_text(size=12),
        axis.text=element_text(size=12)) +
  geom_abline(intercept = 0, slope = 1,color="coral3",linetype = "dashed",
              linewidth = 1)+ 
  annotate("text", x = -Inf , y = -Inf, label = "1:1", color = "coral3", 
           size = 4, hjust = -1, vjust = -0.5, angle = 45)+
  
  geom_label(aes(label = paste0("Correlation between effect sizes:\n r = ",
                                round(corr, 2)),
                 y = 0.4, x = -0.3), size = 5,
             color = "black", hjust = 0, fill = "grey95")+
  ggrepel::geom_label_repel(
    data=dplyr::filter(compare_estimate, abs(effect_sizes_full_model - effect_sizes_without_PQ) >
                         quantile(abs(effect_sizes_full_model - effect_sizes_without_PQ), 0.99) ),
    aes(label= paste(covariate, "-", contribution)),
    size=4, fill = "white",
    min.segment.length = 0.2,
    color = "black", alpha = 0.8,
    direction = "both",
    seed = 1968)+
  scale_x_continuous(limits = c(-0.3, 0.4)) +  
  scale_y_continuous(limits = c(-0.3, 0.4)) +
  coord_equal()


ggsave(filename = here::here("figures", "3_species_traits",
                             "correlation_effect_sizes_model_without_PQ_covariates.jpg"),
       width = 8, height = 7)




#### TEST WITH MODEL WITHOUT INFERRED PQ ####
full_model_file <- "FULL_model_SITE_SCALE_2_chains_1000_thin_200_samples.rds"
PQ_model_file <- "Without_NA_PQ_SITE_SCALE_2_chains_1000_thin_200_samples.rds"



full_model_estimates <- extract_estimates(full_model_file)
PQ_model_estimates <- extract_estimates(PQ_model_file)



PQ_model_estimates_matrix <- PQ_model_estimates[["estimates"]] |> 
  tibble::rownames_to_column("covariate") |> 
  tidyr::pivot_longer(-covariate, names_to = "contribution", 
                      values_to = "effect_sizes_without_inferred_PQ")

full_model_estimates_matrix <- full_model_estimates[["estimates"]]|> 
  tibble::rownames_to_column("covariate") |> 
  tidyr::pivot_longer(-covariate, names_to = "contribution",
                      values_to = "effect_sizes_full_model")


compare_estimate <- PQ_model_estimates_matrix |> 
  dplyr::left_join(full_model_estimates_matrix) |> 
  dplyr::filter(covariate != "(Intercept)")


corr <- cor(compare_estimate$effect_sizes_without_inferred_PQ, compare_estimate$effect_sizes_full_model)

ggplot(compare_estimate, aes(x=effect_sizes_full_model, 
                             y=effect_sizes_without_inferred_PQ)) +
  geom_point(shape = 21, alpha = 0.5,
             fill = "grey50", color= "black")+
  theme_bw()+
  
  xlab(paste0("Effect sizes in full model (N_reef = ",
              full_model_estimates[["n_sampled"]], ")"))+
  ylab(paste0("Effect sizes in model without Photo-quadrat inferred data (N_reef = ",
              PQ_model_estimates[["n_sampled"]], ")")) +
  
  guides(colour = guide_colourbar(title.position="top")) +
  theme(legend.position = "none",
        legend.key.size = unit(0.5, 'cm'),
        legend.direction = "vertical",
        legend.title = element_text( size = 11),
        legend.background = element_rect(fill='transparent'),
        axis.title=element_text(size=12),
        axis.text=element_text(size=12)) +
  geom_abline(intercept = 0, slope = 1,color="coral3",linetype = "dashed",
              linewidth = 1)+ 
  annotate("text", x = -Inf , y = -Inf, label = "1:1", color = "coral3", 
           size = 4, hjust = -1, vjust = -0.5, angle = 45)+
  
  geom_label(aes(label = paste0("Correlation between effect sizes:\n r = ",
                                round(corr, 2)),
                 y = 0.4, x = -0.3), size = 5,
             color = "black", hjust = 0, fill = "grey95")+
  ggrepel::geom_label_repel(
    data=dplyr::filter(compare_estimate, abs(effect_sizes_full_model - effect_sizes_without_inferred_PQ) >
                         quantile(abs(effect_sizes_full_model - effect_sizes_without_inferred_PQ), 0.99) ),
    aes(label= paste(covariate, "-", contribution)),
    size=4, fill = "white",
    min.segment.length = 0.2,
    color = "black", alpha = 0.8,
    direction = "both",
    seed = 1968)+
  scale_x_continuous(limits = c(-0.3, 0.4)) +  
  scale_y_continuous(limits = c(-0.3, 0.4)) +
  coord_equal()


ggsave(filename = here::here("figures", "3_species_traits",
                             "correlation_effect_sizes_model_without_inferred_PQ.jpg"),
       width = 8, height = 7)

##------------------ Surveys / Sites sensibility analysis ---------------------
site_model_file <- "FULL_model_SITE_SCALE_2_chains_1000_thin_200_samples.rds"
survey_model_file <- "FULL_survey_site_country_in_rL_2_chains_1000_thin_200_samples.rds"

site_model_estimates <- extract_estimates(file_name = site_model_file)
survey_model_estimates <- extract_estimates(file_name = survey_model_file)



survey_model_estimates_matrix <- survey_model_estimates[["estimates"]] |> 
  tibble::rownames_to_column("covariate") |> 
  tidyr::pivot_longer(-covariate, names_to = "contribution", 
                      values_to = "effect_sizes_survey_scale")

site_model_estimates_matrix <- site_model_estimates[["estimates"]]|> 
  tibble::rownames_to_column("covariate") |> 
  tidyr::pivot_longer(-covariate, names_to = "contribution",
                      values_to = "effect_sizes_reef_scale")


compare_estimate <- survey_model_estimates_matrix |> 
  dplyr::left_join(site_model_estimates_matrix) |> 
  dplyr::filter(covariate != "(Intercept)")


corr <- cor(compare_estimate$effect_sizes_survey_scale, compare_estimate$effect_sizes_reef_scale)

ggplot(compare_estimate, aes(x=effect_sizes_reef_scale, 
                             y=effect_sizes_survey_scale)) +
  geom_point(shape = 21, alpha = 0.5,
             fill = "grey50", color= "black")+
  theme_bw()+
  xlab("Effect sizes in reef scale full model")+
  ylab("Effect sizes in survey scale full model") +
  
  xlab(paste0("Effect sizes in reef scale full model (N = ",
              site_model_estimates[["n_sampled"]], ")"))+
  ylab(paste0("Effect sizes in survey scale full model (N = ",
              survey_model_estimates[["n_sampled"]], ")")) +
  
  guides(colour = guide_colourbar(title.position="top")) +
  theme(legend.position = "none",
        legend.key.size = unit(0.5, 'cm'),
        legend.direction = "vertical",
        legend.title = element_text( size = 11),
        legend.background = element_rect(fill='transparent'),
        axis.title=element_text(size=12),
        axis.text=element_text(size=12)) +
  geom_abline(intercept = 0, slope = 1,color="coral3",linetype = "dashed",
              linewidth = 1)+ 
  annotate("text", x = -Inf , y = -Inf, label = "1:1", color = "coral3", 
           size = 4, hjust = -1, vjust = -0.5, angle = 45)+
  
  geom_label(aes(label = paste0("Correlation between effect sizes:\n r = ",
                                round(corr, 2)),
                 y = 0.4, x = -0.3), size = 5,
             color = "black", hjust = 0, fill = "grey95")+
  ggrepel::geom_label_repel(
    data=dplyr::filter(compare_estimate, abs(effect_sizes_reef_scale - effect_sizes_survey_scale) >
                         quantile(abs(effect_sizes_reef_scale - effect_sizes_survey_scale), 0.99) ),
    aes(label= paste(covariate, "-", contribution)),
    size=4, fill = "white",
    min.segment.length = 0.2,
    color = "black", alpha = 0.8,
    direction = "both",
    seed = 1968)+
  scale_x_continuous(limits = c(-0.3, 0.4)) +  
  scale_y_continuous(limits = c(-0.3, 0.4)) + 
  coord_equal()


ggsave(filename = here::here("figures", "3_species_traits",
                             "correlation_effect_sizes_surveys_site_models.jpg"),
       width = 8, height = 7)

