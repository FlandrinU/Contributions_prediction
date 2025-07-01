################################################################################'
##
##  This script uses different HMSC outputs tu run sensitivity analysis, by
##   comparing the the covariate effect sizes on contributions of the full model
##   and alternative models.
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

source("R/HMSC_function.R")
source(here::here("R","evaluation_prediction_model.R"))

##----------------------------- data and hmsc outputs ------------------------------
#Coastline
coast <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')

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




plot_boxplot <- function(df = conterfact_CL,
                         val = "raw_change_percent",
                         prop_outliers = 0.015){
  
  df <- df |> 
    dplyr::mutate(contribution = reorder(contribution, .data[[val]],
                                         FUN = median, decreasing = F ))
  ggplot(df) +
    geom_hline(yintercept = 0, color = "grey", size = 1) +
    aes(x= contribution, y= .data[[val]] , fill = group, color = method)+
    scale_fill_manual(values = c( "darkgoldenrod2", "forestgreen", "dodgerblue3")
                      # ,
                      # labels = c("NN" = "Nature-for-Nature",
                      #            "NS" = "Nature-for-Society",
                      #            "NC" = "Nature-as-Culture")
    )+
    scale_color_manual(values = c( "grey70",  "grey10"))+
    
    geom_boxplot(alpha = 0.7, outliers = T) +
    
    xlab("") + ylab("Contributions change in counterfactual scenarios") +
    labs(title="", fill = "")+
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 15, margin = margin(r = 20)),
          panel.spacing = unit(0.3, "lines"),
          axis.text.y = element_text(size = 13),
          axis.title = element_text(size = 13),
          axis.text.x = element_text(angle = 0, hjust = 0.5,size = 13))+
    
    #Flip coordinates and cut outliers
    coord_flip(ylim=c(quantile(df[, val], prop_outliers, na.rm = T), 
                      quantile(df[, val], 1-prop_outliers, na.rm = T)))
  
}



density_plot_fct <- 
  function(changes_df = observed_conservation_legacy ,
           contribution_to_plot =  c("Iucn_species_richness","Actinopterygian_richness"),
           log_transform = T,
           facet_ncol = 1,
           scales = "free",
           x_label = "",
           x_text_CL = 0,
           x_text_HF = 0,
           y_text = 0.02,
           n_round = 1,
           legend.pos = "bottom",
           prop_outliers = 0,
           add_random_distri = F
  ){
    
    data <- changes_df|>
      dplyr::filter(contribution  %in% contribution_to_plot)  |> 
      dplyr::mutate(color = ifelse(group == "NN", "forestgreen",
                                   ifelse(group == "NC", "darkgoldenrod2",
                                          "dodgerblue3")),
                    name = paste0("<span style='color:",
                                  color,
                                  "; font-size: 35px;'>&#9679;</span> ",
                                  contribution))
    
    
    
    
    plot <- ggplot(data)+
      aes(x=raw_change, group=method, fill=method) +
      geom_density(alpha = 0.6)+
      # geom_density(aes(x=raw_original_prediction, group=method, fill=method),alpha = 0.1)+
      
      #Add medians
      geom_point(data = data |> 
                   dplyr::group_by(name, method) |> 
                   dplyr::summarise(median_value = median(raw_change, na.rm = TRUE),
                                    .groups = 'drop'), 
                 aes(x = median_value, y = 0, fill = method), 
                 size = 6, shape = 21, color = "black") +
      
      geom_text(data = data |>
                  dplyr::group_by(name, method) |>
                  dplyr::summarise(median_value = median(raw_change, na.rm = TRUE),
                                   sd = sd(raw_change, na.rm = TRUE),
                                   .groups = 'drop')|>
                  dplyr::mutate(x_text = dplyr::case_when(
                    method == "counterfactual" ~ x_text_CL,
                    method == "observed_IN_OUT" ~ x_text_HF
                  )),
                aes(x = median_value + sign(median_value)*abs(median_value)*x_text,#*sd*x_text,
                    y = y_text,
                    label = round(median_value, n_round), color = method),
                size = 8, vjust = -1,
                show.legend = FALSE) +
      
      scale_fill_manual(
        values = c("observed_IN_OUT" = "chartreuse3",
                   "counterfactual" = "#0077B6"),
        labels = c("Observed conservation legacy \n (in/out design)" , 
                   "Model-based conservation legacy \n (counterfactuals)"),
        name = "Changes:") +
      
      scale_color_manual(
        values = c("observed_IN_OUT" = colorspace::darken("chartreuse3", 0.4),
                   "counterfactual" = colorspace::darken("#0077B6", 0.3)))  +
      
      geom_vline(xintercept = 0)+ # linetype = "dashed", alpha = 0.5,
      geom_hline(yintercept = 0)+
      
      hrbrthemes::theme_ipsum(grid = F, axis = "x", ticks = T,
                              base_size = 17,
                              strip_text_size = 20,
                              axis_title_size = 20) +
      xlab(x_label) + ylab("")+
      # xlim(c(c(quantile(data[, "raw_change"], prop_outliers, na.rm = T), 
      #          quantile(data[, "raw_change"], 1-prop_outliers, na.rm = T)))) + 
      theme(legend.position=legend.pos, panel.spacing = unit(0.1, "lines"),
            legend.text = element_text(size=20),
            legend.title = element_text(size=20),
            # legend.key.spacing = unit(0.7, "cm"),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_text(color = "black"),
            plot.margin = margin(0, 0, 0, 0),
            strip.text = ggtext::element_markdown(size = 20))+
      facet_wrap(~name, ncol = facet_ncol, scales = "free")

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
    
    if(add_random_distri){
      
      df <- data |> 
        dplyr::ungroup() |> 
        dplyr::mutate(raw_change = runif(nrow(data), 
                                         quantile(data$raw_change, .01),
                                         quantile(data$raw_change, .99)))
        
        
      plot <- plot +
        geom_density(data = df, aes(x = raw_change, fill = "Uniform distribution"), group =1,
                     color = "grey50", alpha = 0.3)+
        scale_fill_manual(
          values = c("Uniform distribution" = "grey",
                     "observed_IN_OUT" = "chartreuse3",
                     "counterfactual" = "#0077B6"),
          labels = c("Observed conservation legacy \n (in/out design)" , 
                     "Model-based conservation legacy \n (counterfactuals)",
                     "Uniform distribution"),
          name = "Changes:") 
    }
    plot
  }


plot_mpa <-function(covariates_site_final, xlim=c(-180,180), ylim = c(-36, 31),
                    legend_pos = "none", jitter = 0.2){
  ggplot(covariates_site_final) +
    geom_sf(data = coast, color = "grey30", fill = "lightgrey",
            aes(size=0.1)) +
    
    geom_point(size = 3, na.rm = T,
               position = position_jitter(width =jitter, height =jitter),
               alpha = 0.8,
               colour = "black",
               stroke=0.1,
               shape = 21,
               aes(x = longitude, y = latitude,
                   fill=protection_status)) +
    
    coord_sf(xlim, ylim , expand = FALSE) +
    guides(alpha = "none", size = "none", colour = "none") +
    # scale_shape_manual(values=c(21,24,23))+
    
    theme_minimal()+
    theme(legend.position = legend_pos,
          plot.title = element_text(size=15, face="bold"),
          legend.text = element_text(size=13),
          legend.title = element_text(size=15),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
}

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



##------------------------- Subset IN/OUT design - Use full model ----------------------------

## Matching conditions
time_buffer = 60 #days
spatial_buffer = 50 #km
protection = c("full")#, "restricted")
long_min = -180
long_max = 180

path_file = here::here("figures","3_sensibility_analysis", "counterfactuals")
parameters = paste0(time_buffer, "days_", spatial_buffer, "km_", 
                    paste0(protection, collapse = "_"), 
                    "_MPA_", long_min, "_", long_max, "_longitudes")


#### 1) Site selection ----

## Site data
load(here::here("data/derived_data/3_sites_covariates_to_predict.Rdata"))

## MPA data
mpa_csv <- read.csv("data/raw_data/rls_mpa_MASTER_09012025.csv", header = TRUE)


## Match inside/outside mpa sites 
metadata <- covariates_site_final |> 
  dplyr::mutate(date = stringr::word(rownames(covariates_site_final),
                                     2, sep = "_")) |> 
  dplyr::select(site_code, date, latitude, longitude, country, ecoregion,
                realm, year, protection_status, depth) |> 
  dplyr::mutate(year = stringr::word(date, sep = "-"),
                date = as.Date(date)) |> 
  tibble::rownames_to_column("id")


list_mpa <- metadata |> 
  dplyr::filter(protection_status %in% protection) |>
  dplyr::filter(longitude > long_min & longitude < long_max) |> 
  dplyr::pull(id)

match_sites <- parallel::mclapply(list_mpa, mc.cores = 15,
                                  FUN = function(mpa){
 # mpa = list_mpa[1]
 # mpa = "ACEH1_2009-02-22"
 mpa_data <- metadata |> 
   dplyr::filter(id == mpa)
 coord_mpa <- mpa_data |> 
   dplyr::select(longitude, latitude)
 
 time_lag <- metadata |> 
   dplyr::filter(id != mpa,
                 protection_status == "out") |>
                 # protection_status != "full") |> 
   dplyr::mutate(lag_days = abs(difftime(date, mpa_data$date, units = "days"))) |> 
   dplyr::filter(lag_days <= time_buffer )
 
 
 if(nrow(time_lag) > 0){
   
   distances <- geodist::geodist(
     x = time_lag |> dplyr:: select(longitude, latitude),
     y = coord_mpa,
     measure = "haversine")
   
   matching_sites <- time_lag |> 
     dplyr::mutate(distance_km = as.vector(distances) / 1000) |> 
     dplyr::filter(distance_km <= spatial_buffer)
   
   if(nrow(matching_sites) > 0){
     mpa_df <- mpa_data |> dplyr::mutate(lag_days = 0, distance_km = 0)
     match_df <- rbind(mpa_df, matching_sites) |> 
       dplyr::mutate(match = paste0(mpa_data$id, "_match"))
     
     match_df
   }
 }
})

match_sites_df <- do.call(rbind, match_sites)
#Number of matched sites in MPA
in_mpa <- length(unique(match_sites_df$match)) ; in_mpa
#Number of match outside MPAs
out_mpa <- sum(!unique(match_sites_df$id) %in% list_mpa) ; out_mpa
#Number of MPAs considered
n_mpa <- dplyr::filter(match_sites_df, protection_status != "out") |> 
  dplyr::select(site_code) |> 
  unique() |> 
  dplyr::left_join(mpa_csv) |> 
  dplyr::pull(rls_mpa) |> 
  unique() ; n_mpa

n_mpa_total <- dplyr::filter(metadata, protection_status == "full") |> 
  dplyr::select(site_code) |> 
  unique() |> 
  dplyr::left_join(mpa_csv) |> 
  dplyr::pull(rls_mpa) |> 
  unique() ; n_mpa_total


## Selected sites on map
plot_mpa(match_sites_df,  xlim=c(-180,180), ylim = c(-36, 31),
         legend_pos = "bottom", jitter = 0.3)+
  labs(title = paste0("In/Out-design full MPA (time buffer: ", time_buffer,
                      " days, spatial buffer: ", spatial_buffer, " km)"),
       subtitle = paste0("N = ", length(n_mpa), " MPAs surveyed, ", in_mpa, " samples in MPA, ",
                         out_mpa, " samples in MPA"),
       x="", y= "") 

# ylim = c(5.4, 6),xlim= c(95,95.8)
ggsave(filename = paste0(path_file,"/Selected_sites_IN_OUT_design_",
                         parameters,".jpg"),
       width = 15, height = 8)


#### 2) Check covariates matching ----
load(here::here("data/derived_data/3_sites_covariates_to_predict.Rdata"))

cov_to_check <- c("SST_5_years", "Chlorophyll_5_years", "pH_5_years", "Sand", 
                  "Gravity", "depth", "Terrestrial_reef_flat", "Rubble", "Coral_RLS")
  
  
X <- covariates_site_final |> 
  dplyr::select(all_of(cov_to_check)) |> 
  tibble::rownames_to_column("id")

cov_match <- match_sites_df |> 
  dplyr::left_join(X) |> 
  tidyr::pivot_longer(cols = all_of(cov_to_check),
                      names_to = "covariate",
                      values_to = "val")

# Standardized Mean Difference total group
SDM <- cov_match |> 
  dplyr::group_by(covariate) |> 
  dplyr::summarise(
    mean_out = mean(val[protection_status == "out"], na.rm = TRUE),
    mean_full = mean(val[protection_status == "full"], na.rm = TRUE),
    var_out = var(val[protection_status == "out"], na.rm = TRUE),
    var_full = var(val[protection_status == "full"], na.rm = TRUE),
    pooled_sd = sqrt((var_out + var_full) / 2),
    SMD = (mean_out - mean_full) / pooled_sd
  )
ggplot(SDM, aes(x = reorder(covariate, abs(SMD)), y = SMD)) +
  geom_col(fill = "darkslateblue") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    x = "Covariates",
    y = "Standardized Mean Difference (SMD)",
    title = "Love Plot : Déséquilibre des covariables"
  )


cov <- cov_match |>
  # dplyr::mutate(id = paste0(id, match)) |> 
  dplyr::select(-c(site_code:realm), -year, -lag_days, -distance_km) |> 
  dplyr::group_by(covariate) |> 
  dplyr::mutate(val = as.numeric(scale(val)))

cov_full <- cov |> 
  dplyr::filter(protection_status == "full") |> 
  dplyr::rename(val_full = val) |> 
  dplyr::select(-id, - protection_status)

cov_out <- cov |> 
  dplyr::filter(protection_status == "out") |> 
  dplyr::rename(val_out = val) |> 
  dplyr::select(- protection_status) |> 
  dplyr::group_by(covariate, match) |> 
  dplyr::summarise(val_out = mean(val_out))

cov_diff <- cov_out |> 
  dplyr::left_join(cov_full) |> 
  dplyr::mutate(diff = val_full - val_out)

ggplot(cov_diff, aes(x = diff, y  = covariate)) +
  geom_boxplot(aes(fill = covariate), color = "grey20", alpha = 0.7, outliers = F) +
  geom_jitter(aes(color = covariate),  color = "grey50",alpha = 0.3, height = 0.2) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "firebrick", size = 0.5) +
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "firebrick", size = 0.5) +
  
  labs(
    y = "Covariates",
    x = "Standardized differences of intra-pair covariates"
  ) +
  theme_bw()


#### 3) In/Out comparison in original prediction ----

## Full model prediction
model_name <-"FULL_model_SITE_SCALE_4_chains_1000_thin_200_samples.rds"
load( file = paste0(here::here("figures","3_models","hmsc", "conterfactuals"),
                    "/", gsub(".rds", "", model_name),
                    "/counterfactuals_to_plot.Rdata"))

out_reefs <- conterfactuals[[1]][[3]]

if(all(protection == "full")){
  in_reefs <- conterfactuals[[3]][[3]]
}else{in_reefs <-  conterfactuals[[2]][[3]]}

all_prediction <- rbind(out_reefs, in_reefs) |> 
  dplyr::select(id, contribution, original_prediction, conterfactual,
                raw_original_prediction, raw_conterfactual, raw_change_values,
                raw_change_percent, change_percent_log, group) |> 
  dplyr::filter(id %in% unique(match_sites_df$id)) |> 
  dplyr::left_join(
    dplyr::select(match_sites_df, id, protection_status) |> unique()
  ) |> 
  dplyr::group_by(contribution) |> 
  dplyr::mutate(mean = mean(raw_original_prediction),
                sd = sd(raw_original_prediction))

### In/Out comparison
preds <-  match_sites_df |> 
  dplyr::select(id, match, distance_km) |> 
  dplyr::left_join(all_prediction, relationship = "many-to-many") |>
  dplyr::mutate(site_code = stringr::word( id, sep = "_",1),
                date = stringr::word( id, sep = "_",2),
                year = stringr::word( date, sep = "-",1))|> 
  dplyr::left_join(
    dplyr::select(mpa_csv, site_code, rls_mpa, ps_name, year_of_protection))

in_mpa_preds <- preds |> 
  dplyr::filter(protection_status != "out") |> 
  dplyr::select(match,contribution,
                original_prediction_in_mpa = original_prediction, 
                raw_original_prediction_in_mpa = raw_original_prediction,
                group, protection_status,
                mean, sd, rls_mpa, year)

out_mpa_preds <- preds |> 
  dplyr::filter(protection_status == "out") |> 
  dplyr::select(id, match, contribution,
                original_prediction_out_mpa = original_prediction, 
                raw_original_prediction_out_mpa = raw_original_prediction,
                mean, sd)
                 

observed_conservation_legacy <- in_mpa_preds |> 
  dplyr::left_join(out_mpa_preds) |> 
  dplyr::mutate(change = original_prediction_in_mpa - original_prediction_out_mpa,
                raw_change = raw_original_prediction_in_mpa - raw_original_prediction_out_mpa,
                raw_change_percent = (raw_original_prediction_in_mpa - 
                                        raw_original_prediction_out_mpa)/
                  raw_original_prediction_out_mpa *100,
                method = "observed_IN_OUT",
                
                scaled_pred_in_mpa = (raw_original_prediction_in_mpa - mean)/sd,
                scaled_pred_out_mpa = (raw_original_prediction_out_mpa - mean)/sd,
                scaled_change = scaled_pred_in_mpa - scaled_pred_out_mpa
                )
# |>                                           ########## MEAN BY MATCHS????
#   dplyr::group_by(match, contribution)|> 
#   dplyr::mutate(dplyr::across(where(is.numeric), mean, na.rm = TRUE)) |> 
#   unique()
  

# ## Mean per MPA -> QUITE THE SAME
# averaged_observed_conservation_legacy <- in_mpa_preds |> 
#   dplyr::left_join(out_mpa_preds) |> 
#   dplyr::group_by(rls_mpa, year, contribution, group, mean, sd) |> 
#   dplyr::summarise(original_prediction_in_mpa = mean(original_prediction_in_mpa),
#                    raw_original_prediction_in_mpa = mean(raw_original_prediction_in_mpa),
#                    
#                    original_prediction_out_mpa = mean(original_prediction_out_mpa),
#                    raw_original_prediction_out_mpa = mean(raw_original_prediction_out_mpa)) |> 
#   dplyr::ungroup() |> 
#   dplyr::mutate(change = original_prediction_in_mpa - original_prediction_out_mpa,
#                 raw_change = raw_original_prediction_in_mpa - raw_original_prediction_out_mpa,
#                 raw_change_percent = (raw_original_prediction_in_mpa - 
#                                         raw_original_prediction_out_mpa)/
#                   raw_original_prediction_out_mpa *100,
#                 method = "observed_IN_OUT",
#                 
#                 scaled_pred_in_mpa = (raw_original_prediction_in_mpa - mean)/sd,
#                 scaled_pred_out_mpa = (raw_original_prediction_out_mpa - mean)/sd,
#                 scaled_change = scaled_pred_in_mpa - scaled_pred_out_mpa
#   )


# plot_boxplot(observed_conservation_legacy, val = "raw_change_percent", prop_outliers = 0.015 )
# 
# 
# ggsave(filename = paste0(path_file,
#                          "/Observed_conservation_legacy_IN_OUT_design_", 
#                          parameters,".jpg"),
#        width = 10, height = 8)
# 
# plot_boxplot(observed_conservation_legacy, val = "scaled_change", prop_outliers = 0.015 )
# 




unique(observed_conservation_legacy$contribution)

#Biomass
density_plot_fct(changes_df = observed_conservation_legacy ,
                 contribution_to_plot =  c( "Available_biomass", "Herbivores_biomass",
                                            "Invertivores_biomass", "Piscivores_biomass"),
                 log_transform = T,
                 facet_ncol = 1,
                 scales = "free_y",
                 y_text = 0.07,
                 x_label = "Biomass changes (kg/ha)",
                 add_random_distri = T)


#### 3') Effect sizes log ratio responses ----

## Control/impzct effect sizes
observed_effect_size <- in_mpa_preds |> 
  dplyr::left_join(out_mpa_preds)|> 
  dplyr::group_by(contribution, rls_mpa) |> 
  dplyr::summarise(value_in_mpa = mean(raw_original_prediction_in_mpa),
                   value_out_mpa = mean(raw_original_prediction_out_mpa)) |> 
  dplyr::mutate(log_ratio_response_observed = log10(value_in_mpa/value_out_mpa)) |> 
  dplyr::select(-value_in_mpa, -value_out_mpa)

conterfact_effect_size <- all_prediction |> 
  dplyr::filter(protection_status != "out") |> 
  dplyr::mutate(site_code = stringr::word( id, sep = "_",1)) |> 
  dplyr::left_join(
    dplyr::select(mpa_csv, site_code, rls_mpa, ps_name, year_of_protection)) |> 
  dplyr::group_by(contribution, rls_mpa) |> 
  dplyr::summarise(value_in_mpa = mean(raw_original_prediction),
                   value_out_mpa = mean(raw_conterfactual)) |> 
  dplyr::mutate(log_ratio_response_counterfactual = log10(value_in_mpa/value_out_mpa)) |> 
  dplyr::select(-value_in_mpa, -value_out_mpa)


effect_sizes <- observed_effect_size |> 
  dplyr::left_join(conterfact_effect_size)

ggplot(data = effect_sizes, 
       aes(x = log_ratio_response_observed,
           y = log_ratio_response_counterfactual)) +
  geom_point(alpha = 0.7, size = 3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ contribution, scales = "free") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Observed effect size (log10 ratio)",
    y = "Counterfactual effect size (log10 ratio)",
    title = "Comparison of observed and counterfactual effect sizes"
  )

#### 4) Protected/Unprotected comparison in counterfactuals ----

conterfact_CL <- all_prediction |> 
  dplyr::filter(protection_status != "out") |> 
  dplyr::mutate(method = "counterfactual",
                
                scaled_pred_in_mpa = (raw_original_prediction - mean)/sd,
                scaled_counterfact= (raw_conterfactual - mean)/sd,
                scaled_change = scaled_pred_in_mpa - scaled_counterfact,
                
                raw_change = raw_change_values
  )
  
# plot_boxplot(conterfact_CL, val = "raw_change_percent", prop_outliers = 0.015 )
# 
# 
# ggsave(filename = paste0(path_file,
#                          "/Counterfactual_conservation_legacy_IN_OUT_design_", 
#                          parameters,".jpg"),
#        width = 10, height = 8)


plot_boxplot(conterfact_CL, val = "scaled_change", prop_outliers = 0.015 )


density_plot_fct(changes_df = conterfact_CL ,
                 contribution_to_plot =  c( "Available_biomass", "Herbivores_biomass",
                                            "Invertivores_biomass", "Piscivores_biomass"),
                 log_transform = T,
                 facet_ncol = 1,
                 scales = "free_y",
                 y_text = 0.07,
                 x_label = "Biomass changes (kg/ha)",
                 add_random_distri = T)

#### 5) Compare methods ----

observed_conservation_legacy
conterfact_CL


both_methods <- observed_conservation_legacy |> 
  dplyr::select(id, contribution, group, protection_status, method, raw_change_percent,
                scaled_change, raw_change) |> 
  dplyr::bind_rows(
    dplyr::select(conterfact_CL,
                  id, contribution, group, protection_status, method, raw_change_percent,
                  scaled_change, raw_change)
  ) |> 
  dplyr::mutate(
    raw_change = ifelse(contribution %in% c("Available_biomass", "Herbivores_biomass",
                                         "Invertivores_biomass", "Piscivores_biomass"),
                        raw_change*20/1000, raw_change)) |> 
  dplyr::mutate(site_code = stringr::word(id, sep="_", 1)) |> 
  dplyr::group_by(id, contribution, group, protection_status, method) |> 
  dplyr::summarise(raw_change_percent = median(raw_change_percent),
                   scaled_change = median(scaled_change),
                   raw_change = median(raw_change),
                   n_out = dplyr::n())


# number of comparison in In/Out design
summary(dplyr::filter(both_methods, method == "observed_IN_OUT")$n_out)


# Statistical tests
results <- both_methods |> 
  dplyr::group_by(contribution) |> 
  dplyr::summarise(
    p_value = tryCatch({
      wilcox.test(raw_change ~ method)$p.value
    }, error = function(e) NA),
    n_methods = dplyr::n_distinct(method)) |> 
  dplyr::ungroup()

test <- both_methods |> 
       dplyr::group_by(contribution)  |> 
       dplyr::summarise(
             delta = effsize::cliff.delta(raw_change ~ method)$estimate,
             p_value = wilcox.test(raw_change ~ method)$p.value )


## BOXPLOT OF CHANGES
plot_boxplot(both_methods, val = "raw_change_percent", prop_outliers = 0.025 )


ggsave(filename = paste0(path_file,
                         "/Comparison_methods_conservation_legacy_IN_OUT_design_", 
                         parameters,".jpg"),
       width = 10, height = 8)


plot_boxplot(both_methods, val = "scaled_change", prop_outliers = 0.015 )

# ggsave(filename = paste0(path_file,
#                          "/Comparison_SCALED_methods_conservation_legacy_IN_OUT_design_", 
#                          parameters,".jpg"),
#        width = 10, height = 8)




## DENSITY PLOT OF CHANGES
unique(observed_conservation_legacy$contribution)

biom <- density_plot_fct(changes_df = both_methods ,
                 contribution_to_plot =  c( "Available_biomass", "Herbivores_biomass",
                                            "Invertivores_biomass", "Piscivores_biomass"),
                 log_transform = T,
                 facet_ncol = 1,  scales = "free",
                 y_text = 0.07, x_text_CL = -0.5, x_text_HF = 1,
                 x_label = "Biomass changes (kg/ha)", add_random_distri = T,
                 legend.pos = "none")
biom

sp <- density_plot_fct(changes_df = both_methods ,
                 contribution_to_plot =  c("Actinopterygian_richness") ,
                 log_transform = F, facet_ncol = 1,
                 scales = "free", y_text = 0.07,
                 x_label = "",  x_text_CL = 0.5, x_text_HF = -0.5,
                 legend.pos = "none")
sp


iucn <- density_plot_fct(changes_df = both_methods ,
                       contribution_to_plot =  c("Iucn_species_richness") ,
                       log_transform = F, facet_ncol = 1, scales = "free",
                       y_text = 0.07, x_label = "", x_text_CL = -0.1, x_text_HF = 0.1,
                       legend.pos = "none")
iucn

aest <- density_plot_fct(changes_df = both_methods ,
                       contribution_to_plot =  c("Aesthetic") ,
                       log_transform = F, facet_ncol = 1, scales = "free",
                       y_text = 0.001, x_label = "", n_round = 1,
                       x_text_CL = 0.4, x_text_HF = -0.4,
                       legend.pos = "none")
aest

turn <- density_plot_fct(changes_df = both_methods ,
                         contribution_to_plot =  c("Available_biomass_turnover") ,
                         log_transform = F, facet_ncol = 1, scales = "free",
                         y_text = 0.001, x_label = "", n_round = 2,
                         x_text_CL = -2, x_text_HF = 0.5,
                         legend.pos = "none")
turn


legend <- ggpubr::as_ggplot(
  ggpubr::get_legend(
    density_plot_fct(changes_df = both_methods, contribution_to_plot =  c("Aesthetic"),
                     legend.pos = "bottom", add_random_distri = T)))

(biom + ( iucn /sp / aest / turn)) / legend + plot_layout(heights = c(10,1))
ggsave(filename = paste0(path_file,
                         "/Density_plot_comparison_conservation_legacy_biom_richnes_", 
                         parameters,".jpg"),
       width = 16, height = 12)


others <- density_plot_fct(changes_df = both_methods ,
                 contribution_to_plot =  c( "Selenium",     
                                            "Zinc", "Omega_3", "Calcium", "Iron",                   
                                            "Vitamin_A", "Functional_distinctiveness", "Evolutionary_distinctiveness",
                                            "Endemism", "Functional_entropy", "Phylogenetic_entropy",
                                            "Trophic_web_robustness", "Mean_trophic_level",
                                            "Public_attention"),
                 log_transform = F,
                 facet_ncol = 3,
                 scales = "free",
                 y_text = 0,
                 x_label = "",
                 n_round = 2,
                 x_text_CL = 2, x_text_HF = -2)
others
ggsave(filename = paste0(path_file,
                         "/Density_plot_comparison_conservation_legacy_others_", 
                         parameters,".jpg"),
       width = 16, height = 12)




##------------------------- MPA features sensitivity ----------------------------
library(ggplot2)
library(patchwork)

source(here::here("R/HMSC_function.R"))
source(here::here("R","evaluation_prediction_model.R"))

load(here::here("data/derived_data/3_sites_covariates_to_predict.Rdata"))

color_grad = c("#A50026", "#D73027", "#F46D43", "#FDAE61","#FEE090", "#D8DAEB",
               "#B2ABD2", "#8073AC", "#6D469C", "#603692", "#542788",
               "#473C8B")

## Covariates
X_data_site = covariates_site_final[, -which(colnames(covariates_site_final) ==
                                               "protection_status")] |>
  dplyr::rename(protection_status = protection_status_detailed)

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


## List all files in the directory and choose the model
list_files <- list.files(file.path(path, "out_multi")) 
list_files

model_name <-"Sensitivity_MPA_features_detailed_10y_SITE_SCALE_4_chains_1000_thin_200_samples.rds"
# model_name <- gsub("output_", "", list_files[13]) #choose the wanted file

concatenate_chains = F



##------- New conditions in counterfactual scenarios ---------"

### Initial conditions ###
X <- X_data_site
metadata = metadata_sites


# (1) HUMAN FOOTPRINT - Change unprotected sites to: full protection, 
#  min fishing vessels, min Gravity and max Travel_time
X_pristine_conditions <- X
new_pristine <- rownames(X_pristine_conditions |> dplyr::filter(protection_status == "out"))

X_pristine_conditions[new_pristine, "protection_status"] <- as.factor("full_large_old")
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


## Convervation legacy of full large old MPA only
X_conservation_legacy_full <- X
new_conserv_legacy_full_mpa <- rownames(X_conservation_legacy_full |> 
                                          dplyr::filter(protection_status == "full_large_old"))

X_conservation_legacy_full[new_conserv_legacy_full_mpa, "protection_status"] <- as.factor("out")

X_conservation_legacy_full <- X_conservation_legacy_full |> 
  tibble::rownames_to_column("id") |> 
  dplyr::left_join(mean_country_out) |> 
  dplyr::mutate(Fishing_vessel_density = dplyr::case_when(
    id %in% new_conserv_legacy_full_mpa ~ fishing_out,
    TRUE ~ Fishing_vessel_density)) |> 
  tibble::column_to_rownames("id") |> 
  dplyr::select(-fishing_out)

## Convervation legacy of others full MPA only
X_conservation_legacy_others <- X
new_conserv_legacy_others_mpa <- rownames(X_conservation_legacy_others |> 
                                          dplyr::filter(protection_status == "full_others"))
# | protection_status == "full_large_old"))

X_conservation_legacy_others[new_conserv_legacy_others_mpa, "protection_status"] <- as.factor("out")

X_conservation_legacy_others <- X_conservation_legacy_others |> 
  tibble::rownames_to_column("id") |> 
  dplyr::left_join(mean_country_out) |> 
  dplyr::mutate(Fishing_vessel_density = dplyr::case_when(
    id %in% new_conserv_legacy_others_mpa ~ fishing_out,
    TRUE ~ Fishing_vessel_density)) |> 
  tibble::column_to_rownames("id") |> 
  dplyr::select(-fishing_out)


## Convervation legacy of ALL full MPA
X_conservation_legacy_all <- X
new_conserv_legacy_all_mpa <- rownames(X_conservation_legacy_all |> 
                                            dplyr::filter(protection_status == "full_others" |
                                                            protection_status == "full_large_old"))

X_conservation_legacy_all[new_conserv_legacy_all_mpa, "protection_status"] <- as.factor("out")

X_conservation_legacy_all <- X_conservation_legacy_all |> 
  tibble::rownames_to_column("id") |> 
  dplyr::left_join(mean_country_out) |> 
  dplyr::mutate(Fishing_vessel_density = dplyr::case_when(
    id %in% new_conserv_legacy_all_mpa ~ fishing_out,
    TRUE ~ Fishing_vessel_density)) |> 
  tibble::column_to_rownames("id") |> 
  dplyr::select(-fishing_out)


##---------- Plot changes ----------"
# Select countries we want to plot
threshold_nb_sites = 20
selected_countries <- metadata |> 
  dplyr::count(country) |> 
  dplyr::filter(n > threshold_nb_sites) |> 
  dplyr::pull(country)


# (2) CONSERVATION LEGACY
conservation_legacy_full_mpa <- 
  plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                               X_new_data = X_conservation_legacy_full,
                               metadata,
                               save_name = "Conservation_legacy_full&large&old_MPA",
                               selected_countries,
                               plot_responders_on_map = F,
                               is_counterfactual = TRUE)



order <- conservation_legacy_full_mpa[[3]]  |> 
  dplyr::mutate(contribution = reorder(contribution, raw_change_percent,
                                       FUN = median, decreasing = T )) 
set_order_boxplot <- levels(order$contribution)

conservation_legacy_other_full_mpa <- 
  plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                               X_new_data = X_conservation_legacy_others,
                               metadata,
                               save_name = "Conservation_legacy_others_full_MPA",
                               selected_countries,
                               plot_responders_on_map = F,
                               is_counterfactual = TRUE,
                               set_order_boxplot = set_order_boxplot)



conservation_legacy_all_full_mpa <- 
  plot_conterfactual_scenarios(path, model_name, concatenate_chains,
                               X_new_data = X_conservation_legacy_all,
                               metadata,
                               save_name = "Conservation_legacy_all_full_MPA",
                               selected_countries,
                               plot_responders_on_map = F,
                               is_counterfactual = TRUE,
                               set_order_boxplot = set_order_boxplot)



## Save conterfactuals ##
conterfactuals <- list(conservation_legacy_other_full_mpa, conservation_legacy_full_mpa,
                       conservation_legacy_all_full_mpa)
save(conterfactuals, file = paste0(here::here("figures","3_models","hmsc", "conterfactuals"),
                                   "/", gsub(".rds", "", model_name), 
                                   "/counterfactuals_to_plot.Rdata"))



##--------------- DENSITY PLOTS ----------------"
#Counterfactuals colors:
Cl_color = "#0077B6"#"darkslategrey"
HF_color = "chartreuse3"# "#52B788"#
all_full_mpa_color = "#6C757D"

folder_name <- gsub(".rds", "", model_name)
path_file <- here::here("figures","3_models","hmsc", "conterfactuals", folder_name)    

# ## load file
# load(file = paste0(here::here("figures","3_models","hmsc", "conterfactuals"),
#                                                        "/", gsub(".rds", "", model_name),
#                                                       "/counterfactuals_to_plot.Rdata"))
# conservation_legacy_other_full_mpa <- conterfactuals[[1]]
# conservation_legacy_full_mpa <- conterfactuals[[2]] 
# conservation_legacy_all_full_mpa <- conterfactuals[[3]]


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


conserv_legacy <- conservation_legacy_full_mpa[[3]] |> 
  dplyr::select(id, contribution, country, changes = raw_change_values, raw_original_prediction) |> 
  dplyr::mutate(counterfactual = "CL")

hum_footprint <- conservation_legacy_other_full_mpa[[3]] |> 
  dplyr::select(id, contribution, country, changes = raw_change_values, raw_original_prediction) |> 
  dplyr::mutate(counterfactual = "HF")

all_mpa_legacy <- conservation_legacy_all_full_mpa[[3]] |> 
  dplyr::select(id, contribution, country, changes = raw_change_values, raw_original_prediction) |> 
  dplyr::mutate(counterfactual = "all_full_mpa")

changes_df <- hum_footprint |> 
  dplyr::bind_rows(conserv_legacy) |> 
  # dplyr::bind_rows(all_mpa_legacy) |>
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
           add_quantile = 0.5, log_transform = T, facet_ncol = 1, scales = "free",
           x_label = "",x_text_CL = 0, x_text_HF = 0,  y_text_CL=0, y_text_HF = 0, 
           n_round = 1, legend.pos = "none" ){
    
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
                    counterfactual == "HF" ~ x_text_HF , 
                    counterfactual == "all_full_mpa" ~ x_text_HF),
                    y_text = dplyr::case_when(
                      counterfactual == "CL" ~ y_text_CL, 
                      counterfactual == "HF" ~ y_text_HF )),
                aes(x = median_value + sign(median_value)*sd*x_text,
                    y = y_text, 
                    label = round(median_value, n_round), color = counterfactual),
                size = 8, vjust = -1,
                show.legend = FALSE) +
      
      scale_fill_manual(
        values = c(HF_color,Cl_color), #, all_full_mpa_color),
        labels = c("HF" = "Others full MPA", "CL" = "Full MPA > 10 years"),
                   # "all_full_mpa" = "All full MPA"),
        limits = c("HF", "CL"),#, "all_full_mpa"),
        name = "Contribution changes in counterfactuals:") +
      scale_color_manual( values = c(#colorspace::darken(all_full_mpa_color, 0.5),
                                     colorspace::darken(Cl_color, 0.2), 
                                     colorspace::darken(HF_color, 0.4)) ) +
      
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
                 y_text_HF = -0.05, y_text_CL = 0.1,
                 x_label = "Biomass changes (kg/ha)")
ggsave(filename =  paste0(path_file,"/density_changes_biomass.jpg"), width = 8, height = 12)


# Aesthetic
density_plot_fct(contribution_to_plot =  c("Aesthetic"),
                 add_quantile = 0.5,
                 x_text_CL = 5,x_text_HF = -2.5,
                 y_text_HF = 0.0, y_text_CL = 0,
                 log_transform = T)
ggsave(filename =  paste0(path_file,"/density_changes_aesth.jpg"), width = 8, height = 3)


# IUCN
density_plot_fct(contribution_to_plot =  c("Iucn_species_richness"),
                 add_quantile = 0.5,
                 log_transform = F,
                 n_round = 2,
                 x_text_CL = 0,
                 y_text_HF = 0.0, y_text_CL = 0)
ggsave(filename =  paste0(path_file,"/density_changes_IUCN.jpg"), width = 8, height = 3)


# Richness
density_plot_fct(contribution_to_plot =  c("Actinopterygian_richness"),
                 add_quantile = 0.5,
                 log_transform = F,
                 x_text_CL = -3, x_text_HF = +2.7,
                 y_text_HF = 0.0, y_text_CL = 0,
                 n_round = 2)
ggsave(filename =  paste0(path_file,"/density_changes_richness.jpg"), width = 8, height = 3)


# Biomass turnover
density_plot_fct(contribution_to_plot =  c("Available_biomass_turnover"),
                 add_quantile = 0.5,
                 log_transform = F,
                 x_text_CL = -1,  x_text_HF = 2,
                 y_text_HF = 0.0, y_text_CL = 0,
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


### PANEL DENSITY PLOT ###

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

ggsave(filename =  paste0(path_file,"/Figure_4_density_plot", "_", model_name, ".jpg"),
       plot = final_panel,  width = 16, height = 12.5)









##------------------------- Spatial autocorrelation - random factors ----------------------------

## Spatial autocorrelation in full model
path_file <- here::here("figures", "3_models", "hmsc",
                        "FULL_model_SITE_SCALE_4_chains_1000_thin_200_samples")
load(file = paste0(path_file,"/Residuals_spatial_correlation_Moran_indices.Rdata"))

moranI_df_full <- moranI_df |> dplyr::rename(Moran_index_full_model = moran_index)

## Spatial autocorrelation in model without reef in random
path_file <- here::here("figures", "3_models", "hmsc",
                        "Sensitivity_only_country_in_RL_SITE_SCALE_4_chains_1000_thin_200_samples")
load(file = paste0(path_file,"/Residuals_spatial_correlation_Moran_indices.Rdata"))

moranI_df_only_country <- moranI_df |> dplyr::rename(Moran_index_only_country_in_RL = moran_index)


## Spatial autocorrelation in model without any random effect
path_file <- here::here("figures", "3_models", "hmsc",
                        "Sensitivity_NO_random_factors_SITE_SCALE_4_chains_1000_thin_200_samples")
load(file = paste0(path_file,"/Residuals_spatial_correlation_Moran_indices.Rdata"))

moranI_df_NO_random <- moranI_df |> dplyr::rename(Moran_index_NO_RL = moran_index)



## Compare spatial autocorrelation
moranI_df <- dplyr::left_join(moranI_df_full, moranI_df_only_country)|> 
  dplyr::left_join(moranI_df_NO_random) |> 
  dplyr::arrange(Moran_index_NO_RL) |> 
  dplyr::mutate(response = factor(response, levels = response))

ggplot(moranI_df, aes(y = response)) +
  geom_point(aes(x = Moran_index_full_model, color = "full_model"), size = 4) + 
  geom_point(aes(x = Moran_index_only_country_in_RL, color = "only_country"), size = 4) + 
  geom_point(aes(x = Moran_index_NO_RL, color = "no_rd"), size = 4) + 
  
  geom_vline(xintercept = mean(moranI_df$Moran_index_full_model),
             linetype = "dashed", color = "skyblue") + 
  geom_vline(xintercept = mean(moranI_df$Moran_index_only_country_in_RL),
             linetype = "dashed", color = "orange") +
  geom_vline(xintercept = mean(moranI_df$Moran_index_NO_RL),
             linetype = "dashed", color = "grey50") +
  
  annotate("text", x = mean(moranI_df$Moran_index_full_model) - 0.02, y = -Inf, 
           label = round(mean(moranI_df$Moran_index_full_model), 2), 
           vjust = -1.5, color = "skyblue", size = 6) +
  annotate("text", x = mean(moranI_df$Moran_index_only_country_in_RL) + 0.02, y = -Inf, 
           label = round(mean(moranI_df$Moran_index_only_country_in_RL), 2), 
           vjust = -1.5, color = "orange", size = 6) +
  annotate("text", x = mean(moranI_df$Moran_index_NO_RL) + 0.02, y = -Inf, 
           label = round(mean(moranI_df$Moran_index_NO_RL), 2), 
           vjust = -1.5, color = "grey50", size = 6) +
  
  labs(x = "Moran index", y = "Contributions", 
       title = "Spatial autocorrelation in model residuals",
       color = "Moran indices") +
  scale_color_manual(values = c("full_model" = "skyblue",
                                "only_country" = "orange",
                                "no_rd" = "grey50"),
                     labels = c("\nFull model \n(country and sampled \nreef in random)\n",
                                 "\n \n No random effect \n \n", "Only country in random")) +
  # xlim(c(0,1))+
  theme_bw()+
  theme(title = element_text(size = 17),
        legend.text =  element_text(size = 13),
        axis.text = element_text(size = 15),     
        axis.title = element_text(size = 15),   
        strip.text = element_text(size = 13))


ggsave(filename = paste0(here::here("figures/3_sensibility_analysis"), 
                         "/Spatial_autocorrelation_in_residuals_according_models.jpg"),
       width = 15, height = 8)
