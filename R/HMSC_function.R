################################################################################
##
##  This function contains all code necessary to run Hmsc models and plot results
##
## HMSC_function.R
##
## 07/05/2024
##
## Ulysse Flandrin
##
################################################################################

# load(here::here("data/derived_data/3_all_contributions_to_predict.Rdata"))
# load(here::here("data/derived_data/3_all_covariates_to_predict.Rdata"))
# Y_data =  observations_final#[sample(1:nrow(observations_final),1000),] ####################### reduce data
# X_data = covariates_final[rownames(Y_data),]
# rownames(X_data) <- rownames(Y_data)
# 
# nSamples = 10 #1000
# thin = 20 #100
# nChains = 2 #2
# verbose = 100 #100
# transient = nSamples * thin # 1500 #100 * thin
# nb_neighbours = 10
# random_factors = c("spatial","country")
# name = "test"
# response_distribution <- rep("normal", ncol(Y_data))
# #response_distribution[colnames(Y_data) == "iucn_species_richness"] <- "poisson"
# run_python = TRUE
# save_path = here::here("outputs/models/hmsc")


#' HMSC function.
#'
#' @param nSamples Number of samples in the markov chain
#' @param thin Number of iterations between each sample
#' @param nChains Number of independant chains
#' @param verbose Frequence of messages diplayed
#' @param transient Number of iteration before sample begining
#' @param Y_data dataframe of responses
#' @param X_data dataframe of covariates
#' @param response_distribution statistical distribution of responses ("normal", "poisson", or...)
#' @param random_factors list of random factors as character
#' @param nb_neighbours number of neighbours taken in spatial model
#' @param name file name saved
#' @param run_python if FALSE, the function creates the model and init object but doesn't run the markov chains
#'
#' @return save hmsc outputs in save_out directory
#' @export
#'
#' @examples
#' 



hmsc_function <- function(nSamples,
                          thin,
                          nChains,
                          verbose,
                          transient,
                          Y_data,
                          X_data,
                          response_distribution,
                          random_factors = NULL,
                          nb_neighbours = NULL,
                          name,
                          run_python = TRUE,
                          save_path = here::here("outputs/models/hmsc")
                          ){
  
  set.seed(0612)
  
  # PATHS
  save_init <- paste0(save_path,"/init_multi")
  save_out <- paste0(save_path, "/out_multi")
  localDir <- paste0(save_path, "/multivariate")
  
  
  cat("-------------- Initialisation of Hmsc model --------------\n")
  
  ## Set fixed effects ##
  fixed_effects <- colnames(dplyr::select(X_data,
                                          -longitude, -latitude, 
                                          -year,
                                          -country,
                                          -ecoregion,
                                          -realm,
                                          
                                          # country level covariates
                                          -hdi,
                                          -marine_ecosystem_dependency,
                                          -ngo,
                                          -natural_ressource_rent
  ))
  
  formula <- as.formula(paste("~ ", paste(fixed_effects, collapse = "+") ) )
  
  response_distribution <- rep("normal", ncol(Y_data))
  # response_distribution[colnames(Y_data) == "iucn_species_richness"] <- "poisson"
  
  
  ## Set random effects ##
  
  # setting model structure with spatial structure ‘Nearest Neighbour Gaussian Process (NNGP)’
  # (1) add noise to coordinates for them to differ sltitly
  noise_magnitude <- 0.00001 #noise in coordinates
  
  X <- X_data |>
    dplyr::rowwise() |> 
    dplyr::mutate(latitude = latitude + runif(1, -noise_magnitude, noise_magnitude),
                  longitude = longitude + runif(1, -noise_magnitude, noise_magnitude)) |>
    dplyr::ungroup()
  rownames(X) <- rownames(X_data)
  
  # (2) create a matrix with coordinates : xycoords is a matrix with 2 columns "x-coordinate","y-coordinate" and row names with spygen_code
  xycoords <- X |>
    dplyr::select(longitude, latitude)
  rownames(xycoords) <- rownames(X)
  
  
  
  # (3) create spatial random effect and study design (rows in Y)
  rL.nngp = Hmsc::HmscRandomLevel(sData = xycoords,
                                  sMethod = 'NNGP',
                                  nNeighbours = nb_neighbours,
                                  # sMethod = 'Full', #too much memory needed...
                                  longlat = F) #Should be True but doesn't work after
  
  rL.nngp = Hmsc::setPriors(rL.nngp, nfMin=1,nfMax=1)
  
  studyDesign <- data.frame(associations = as.factor(rownames(X)),
                            spatial = as.factor(rownames(X)),
                            year = as.factor(X$year),
                            country = as.factor(X$country),
                            ecoregion = as.factor(X$ecoregion))
  
  # (4) Set random levels
  rL_asso = Hmsc::HmscRandomLevel(units = studyDesign$associations)
  rL_year = Hmsc::HmscRandomLevel(units = studyDesign$year)
  rL_ecoregion = Hmsc::HmscRandomLevel(units = studyDesign$ecoregion)
  rL_country =  Hmsc::HmscRandomLevel(units = studyDesign$country)
  
  ranLevels = list(associations = rL_asso,
                   year = rL_year,
                   ecoregion = rL_ecoregion,
                   country = rL_country, 
                   spatial = rL.nngp)
  
  
  ## Construct HMSC model strucure ##
  if(is.null(random_factors)){
    model_fit = Hmsc::Hmsc(Y = Y_data,
                           XData = X[, fixed_effects],
                           XFormula = formula,
                           distr = response_distribution) 
  }else{
    studyDesign <- studyDesign |> dplyr::select(all_of(random_factors))
    ranLevels <- ranLevels[random_factors]
    
    model_fit = Hmsc::Hmsc(Y = Y_data,
                           XData = X[, fixed_effects],
                           XFormula = formula,
                           studyDesign = studyDesign,
                           ranLevels = ranLevels,
                           distr = response_distribution) # see help for other distributions
  }
  
  
  
  # create object for computation on HMSC-HPC
  init_obj = Hmsc::sampleMcmc(model_fit,
                              samples=nSamples,
                              thin=thin,
                              transient=transient, 
                              nChains=nChains,
                              verbose=verbose,
                              nParallel = nChains,
                              engine="HPC") # HPC : try to use the GPU under python command: on github version only (v3.1-2)
  
  
  # save it locally
  file_name <- sprintf(paste0(name, "_", paste(nChains, "chains",
                                          thin, "thin",
                                          nSamples, "samples",
                                          sep = "_"),
                              ".rds"))
  
    
    
  init_file_path = file.path(save_init, paste0("init_", file_name))
  saveRDS(jsonify::to_json(init_obj), file=init_file_path)
  
  save(model_fit, file = file.path(localDir, paste0("model_fit_", file_name, ".Rdata")))
  
  

  
  cat("-------------- Initialisation of python for Hmsc --------------\n")
  
  python = file.path(getwd(),"HMSC_package","hmsc-venv", "bin", "python")  # hmsc-venv for Linux and macOS
  
  # Define the output file path
  post_file_path <- file.path(save_out, paste0("output_", file_name))
  
  # Construct the Python command
  python_cmd_args <- paste("-m hmsc.run_gibbs_sampler",
                                      "--input", shQuote(init_file_path),
                                      "--output", shQuote(post_file_path),
                                      "--samples", nSamples,
                                      "--transient", format(transient, scientific = FALSE),
                                      "--thin", thin,
                                      "--verbose", verbose)
  
  cat(paste(shQuote(python), python_cmd_args), "\n")
  
  
  if(run_python){
    cat("-------------- Run Hmsc with python --------------\n")
    system2(python, python_cmd_args)
  }
  
}







#' Plot results Hmsc.
#'
#' @param covariates dataframe of covariates
#' @param response dataframe of responses
#' @param file_name name of the model fitted
#' @param save_init directory for models
#' @param save_out directory for models
#' @param localDir directory for models
#' @param plot_convergence TRUE or FALSE, if you want to plot figures related to convergence of the model
#' @param plot_explanatory_power TRUE or FALSE, if you want to plot figures related to explanatory power
#'
#' @return save hmsc plots in figures/hmsc/file_name directory
#' @export
#'
#' @examples
#' 



plot_hmsc_result <- function(covariates = covariates,
                             response = response,
                             file_name = file_name,
                             save_init = save_init,
                             save_out = save_out,
                             localDir = localDir,
                             plot_convergence = T,
                             plot_explanatory_power = T,
                             plot_variance_partitioning = T,
                             plot_residual_associations = T,
                             plot_estimates = T,
                             check_residuals = T,
                             drivers_to_plot =  list(
                               c("n_fishing_vessels", "gravtot2", "gdp", "neartt"))
                            ){
  
  
  ##----Source packages and functions----
  
  library(ggplot2)
  library(patchwork)
  
  source("R/evaluation_prediction_model.R")
  
  # Hmsc::computeVariancePartitioning creates bug in Rstudio because it can fill the console of warning.
  # to clear it you can run 'cat("\014")' or make Ctrl + L , to continue.
  # Otherwise, use the modified following function of computeVariancePartitioning
  source("R/HMSC_computeVariancePartitioning.R")
  
  
  ##---- Import initial object ----
  cat("Import HMSC fitted model... \n")
  
  save_name <- gsub(".rds", "", file_name)
  
  #Model design
  load(file = file.path(localDir, paste0("model_fit_", file_name, ".Rdata")))
  
  #Initial hmsc object
  init_obj_rds <- readRDS(paste0(save_init, paste0("init_", file_name)))
  init_obj <- jsonify::from_json(init_obj_rds)
  
  nSamples = init_obj[["samples"]]
  thin = init_obj[["thin"]]
  nChains = init_obj[["nChains"]]
  transient = init_obj[["transient"]]
  
  
  ## Import posterior probability
  file_rds <- readRDS(file = paste0(save_out, paste0("output_", file_name)))[[1]]
  importFromHPC <- jsonify::from_json(file_rds)
  postList <- importFromHPC[1:nChains]
  ### /!\ CHECK THE STRUCTURE OF THE RESULT: MATRICES IN PYTHON VS LIST IN R ###
  
  ##Result model
  cat(sprintf("fitting time %.1f h\n", importFromHPC[[nChains+1]] / 3600))
  
  ##Export and merge chains
  model_fit_mcmc <- Hmsc::importPosteriorFromHPC(model_fit,
                                                 postList, 
                                                 nSamples, 
                                                 thin, 
                                                 transient)
  
  
  ## Estimates for each chains
  mpost <- Hmsc::convertToCodaObject(model_fit_mcmc)
  
  postBeta <- Hmsc::getPostEstimate(model_fit_mcmc, parName = "Beta")
  
  
  ##---- create directory to save figures ----
  path_file <- here::here("figures","models","hmsc", save_name)    
  dir.exists(path_file)
  
  if(!dir.exists(path_file)) dir.create(path_file)

  
  ##----- Check convergence ----
  
  if(plot_convergence){
    cat("Plot convergence of the model... \n")
  
  ### Effective size ###
  
  # Estimate convergence -> effective size should be around nSamples*nChains, and psrf around 1
  # ESS = the number of posterior row effectively independent. low ESS => poor prediction
  png(paste0(path_file,"/convergence_estimates_effective_size_", save_name,".png"),
      width = 20, height = 20, units = "cm", res = 200)
  par(mfrow=c(2,2), mar = c(4,4,4,4))
  hist(coda::effectiveSize(mpost$Beta), main="ess(beta)", breaks = 30)
  abline(v= nSamples*nChains, col = "red4", lty = 2, lwd = 3)
  hist(coda::gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)", breaks = 30)
  abline(v= 1, col = "red4", lty = 2, lwd = 3)
  hist(coda::effectiveSize(mpost$Omega[[1]]), main="ess(omega)", breaks = 30) #also check associations between y
  abline(v= nSamples*nChains, col = "red4", lty = 2, lwd = 3)
  hist(coda::gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf, main="psrf(omega)", breaks = 30)
  abline(v= 1, col = "red4", lty = 2, lwd = 3)
  dev.off()
  
  
  
  ### Chains convergence ###
  # install.packages("ggmcmc")
  # library("ggmcmc")
  # see help at : http://xavier-fim.net/post/using_ggmcmc/
  
  S <- ggmcmc::ggs(mpost$Beta)
  cov <- "n_fishing_vessels"
  S <- S[grep(paste(cov, collapse = "|"), S$Parameter),]
  
  # # quick look on the distribution of the values and the shape of the posterior distribution.
  # ggmcmc::ggs_histogram(S)+facet_wrap(~ Parameter, ncol = 5)
  # 
  # # ggs_density() allows comparing the target distribution by chains and whether each
  # # chain has converged in a similar space.
  # ggmcmc::ggs_density(S)+facet_wrap(~ Parameter, ncol = 5)
  
  # assess convergence and diagnose chain problems. the expected outcome is to produce 
  # “white noise”. + a good tool to assess within-chain convergence
  ggmcmc::ggs_traceplot(S)+facet_wrap(~ Parameter, ncol = 5)
  ggsave(filename = paste0(path_file,"/traceplot_", save_name,".jpg"),
         width = 15, height = 10)
  
  # The expected output is a line that quickly approaches the overall mean, in addition
  # to the fact that all chains are expected to have the same mean
  ggmcmc::ggs_running(S)+facet_wrap(~ Parameter, ncol = 5)
  ggsave(filename = paste0(path_file,"/Chain_convergence_", save_name,".jpg"),
         width = 15, height = 10)
  
  
  # # Ideally, the initial and final parts of the chain have to be sampling in the same 
  # # target distribution
  # ggmcmc::ggs_compare_partial(S)+facet_wrap(~ Parameter, ncol = 5)
  
  # The autocorrelation plot expectes a bar at one in the first lag, but no autocorrelation
  # beyond it. While autocorrelation is not per se a signal of lack of convergence,
  # it may indicate some misbehaviour of several chains or parameters, or indicate that
  # a chain needs more time to converge. 
  ggmcmc::ggs_autocorrelation(S)+facet_wrap(~ Parameter, ncol = 5)
  ggsave(filename = paste0(path_file,"/autocorrelation_", save_name,".jpg"),
         width = 15, height = 10)
  
  # # diagnose potential problems of convergence due to highly correlated parameters
  # ggmcmc::ggs_crosscorrelation(S)
  # 
  # # compare the between-chain variation with the within-chain variation. It is expected to be close to 1.
  # ggmcmc::ggs_Rhat(S) + xlab("R_hat")
  # ggmcmc::ggs_caterpillar(S)
  
  } #END OF PLOT CONVERGENCE
  
  
  ##----------------------------- Explanatory power -----------------------------------
  
  if(plot_explanatory_power){
    cat("Plot explanatory power of the model... \n")
    
    ### Explanatory  power ###
    preds <- Hmsc::computePredictedValues(model_fit_mcmc)
    MF <- Hmsc::evaluateModelFit(hM=model_fit_mcmc, predY=preds)
    
    png(paste0(path_file,"/explanatory_power_", save_name,".png"),
        width = 20, height = 10, units = "cm", res = 300)
    par(mfrow=c(1,2), mar = c(4,4,4,4))
    hist(MF$R2, xlim = c(0,1), main=paste0("Mean R2 = ", round(mean(MF$R2),2)))
    hist(MF$RMSE, xlim = c(0,1), main=paste0("Mean RMSE = ", round(mean(MF$RMSE),2)))
    dev.off()
  } #END OF PLOT EXPLANATORY POWER
  
  ##----------------------------- Covariates estimates -----------------------------------
  
  if(plot_variance_partitioning){
    cat("Complute variance partitioning of the model... \n")
    
    #### Variance partitioning ####
    VP <- computeVariancePartitioning(model_fit_mcmc)
    VP$R2T$Beta ############################################################## why 0????
    
    # svg("figures/models/hmsc/variance_partitioning.svg")
    # Hmsc::plotVariancePartitioning(model_fit_mcmc, VP = VP)
    # dev.off()
    
    ##### Preping VP table #####
    VP_long <- as.data.frame(VP[["vals"]]) |> 
      tibble::rownames_to_column(var = "Covariate") |> 
      tidyr::pivot_longer(
        cols = - Covariate,
        names_to = "Response",
        values_to = "Value"
      )
    
    #classify covariates
    human <- c("gdp", "gravtot2", "effectiveness", "natural_ressource_rent", 
               "neartt","n_fishing_vessels")
    habitat <- c("depth", "algae", "coral", "Sand", "seagrass", "microalgal_mats",
                 "other_sessile_invert", "Rock", "coralline_algae", "coral_rubble",
                 as.character(unique(VP_long$Covariate)[grepl("500m", unique(VP_long$Covariate))]))
    envir <-  c(as.character(unique(VP_long$Covariate)[grepl("median", unique(VP_long$Covariate))]),
                "q05_1year_degree_heating_week", "q95_1year_degree_heating_week", 
                "q95_5year_degree_heating_week")
    random <- unique(VP_long$Covariate)[grepl("Random", unique(VP_long$Covariate))]
    
    # cat_colors <- list(
    #   human = colorRampPalette(c("#D0E1F9", "#1E429F"))(length(human)),
    #   habitat = colorRampPalette(c("#F8D7DA", "#EC7063"))(length(habitat)),
    #   envir = colorRampPalette(c("#D4EDDA", "#58D68D"))(length(envir)),
    #   random = colorRampPalette(c("grey90", "grey50"))(length(random))
    # )
    
    # harrypotter::harrypotter(8, option = "ronweasley2")
    # c("#CC6229", "#FFC199")
    cat_colors <- list(
      human = colorRampPalette(c("#452C52FF", "#D5B3D1"))(length(human)),
      habitat = colorRampPalette(c("#FFA00AFF", "#FFDD99"))(length(habitat)),
      envir = colorRampPalette(c("#FF7A33FF", "#FFC199"))(length(envir)),
      random = colorRampPalette(c("grey95", "grey50"))(length(random))
    )
    
    VP_long <- VP_long |> 
      dplyr::mutate(category = dplyr::case_when(
        Covariate %in% human ~ "human",
        Covariate %in% habitat ~ "habitat",
        Covariate %in% envir ~ "envir",
        Covariate %in% random ~ "random",
        TRUE ~ NA_character_
      )) |> 
      dplyr::rowwise() |> 
      dplyr::mutate(color = if (!is.na(category)) cat_colors[[category]][match(Covariate, get(category))] else NA_character_) |> 
      dplyr::ungroup()
    
    # Order covariates
    VP_long$Covariate <- forcats::fct_relevel(VP_long$Covariate, 
                                              c(human, habitat, envir, random))
    
    # measure where to place labels in relative variance partitioning
    VP_long <- VP_long |> dplyr::arrange(Covariate)
    VP_long$Symbol <- as.character(as.numeric(factor(VP_long$Covariate)))
    VP_long$labels <- paste0(VP_long$Symbol, ": ", VP_long$Covariate)
    VP_long$mid_y <- NA
    
    for(resp in unique(VP_long$Response)){
      df_resp <- VP_long |> dplyr::filter(Response == resp)
      labs <- rev(df_resp$labels)
      
      mid <- 0
      sum <- 0
      for(cov in labs){
        val <- dplyr::filter(df_resp, labels == cov)$Value
        mid <- sum + val/2
        sum <- sum + val
        
        VP_long$mid_y[VP_long$Response == resp & VP_long$labels == cov] <- mid
      }
    }
    
    
    # Aggregated
    VP_aggregated <- VP_long |> 
      dplyr::group_by(category, Response) |> 
      dplyr::summarise(prop_variance = sum(Value)) |> 
      dplyr::group_by(category) |> 
      dplyr::summarise(sd = sd(prop_variance),
                       prop_variance = mean(prop_variance))
    
    # VP absolute
    variance_explained <- data.frame(Response = model_fit_mcmc[["spNames"]],
                                     R2 = MF$R2)
    VP_long_absolute <- VP_long |> 
      dplyr::left_join(variance_explained) |> 
      dplyr::mutate(Value = Value * R2)
    
    for(resp in unique(VP_long_absolute$Response)){
      df_resp <- VP_long_absolute |> dplyr::filter(Response == resp)
      labs <- rev(df_resp$labels)
      
      mid <- 0
      sum <- 0
      for(cov in labs){
        val <- dplyr::filter(df_resp, labels == cov)$Value
        mid <- sum + val/2
        sum <- sum + val
        
        VP_long_absolute$mid_y[VP_long_absolute$Response == resp &
                                 VP_long_absolute$labels == cov] <- mid
      }
    }
    
    ##### Plot VP #####
    
    ## Sum of VP per category
    
    ggplot(VP_aggregated) +
      geom_col(aes(x = reorder(category, prop_variance),
                   y = prop_variance, fill = category)) +
      geom_errorbar(aes(x = category, y = prop_variance,
                        ymin=prop_variance-sd, ymax=prop_variance+sd),
                    width=.1, position=position_dodge(.9)) +
      scale_fill_manual(values = c("random" =  "#CBCBCB",
                                   "envir" = "#FFA976",
                                   "habitat" = "#FFCF7A",
                                   "human" = "#9B7D9E")) +
      theme_minimal() +
      coord_flip() +
      labs(y = "Proportion in the variance explained", x = "") +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12)) 
    ggsave(filename = paste0(path_file,"/variance_explained_per_category_",
                             save_name,".jpg"), width = 15, height = 10)
    
    
    
    ## Relative variance partitioning
    ggplot(VP_long, aes(x = Response, y = Value, fill = Covariate)) +
      geom_bar(stat = "identity", position = "stack", color = "black", size = 0.1) +
      scale_fill_manual(values = setNames(VP_long$color, VP_long$Covariate),
                        labels = unique(VP_long$labels)) +
      labs(title = "", x = "", y = "") +
      theme_classic(base_size = 11,
                    base_line_size = 0.1) +
      theme(
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 6))+
      geom_text(data = dplyr::filter(VP_long, Value > 0.05),
                aes(y = mid_y, label = Symbol), size = 2, color = "black") 
    
    
    ggsave(filename =  paste0(path_file,"/variance_partitioning_", save_name,".jpg"),
           width = 15, height = 10)
    
    
    
    ## VP with absolute value of variance explained
    variance_explained 
    
    ggplot(VP_long_absolute, aes(x = reorder(Response,-R2),
                                 y = Value, fill = Covariate)) +
      geom_bar(stat = "identity", position = "stack", color = "black", size = 0.1) +
      scale_fill_manual(values = setNames(VP_long_absolute$color, VP_long_absolute$Covariate),
                        labels = unique(VP_long_absolute$labels)) +
      labs( title = "",x = "",y = "") +
      theme_classic(base_size = 11,
                    base_line_size = 0.1) +
      theme(
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 6))+
      geom_text(data = dplyr::filter(VP_long_absolute, Value > 0.03),
                aes(y = mid_y, label = Symbol), size = 2, color = "black") 
    
    ggsave(filename =  paste0(path_file,"/variance_partitioning_absolute_values_", save_name,".jpg"),
           width = 15, height = 10)
    
    
    
    ## Mean contributions of covariates
    covariate_contrib <- VP_long |> 
      dplyr::group_by(Covariate, color, category) |> 
      dplyr::summarise(sd = sd(Value),
                       contribution = mean(Value)) |> 
      dplyr::filter(!grepl("Random",Covariate))
    
    ggplot(covariate_contrib) +
      geom_col(aes(x = reorder(Covariate, contribution),
                   y = contribution, fill = category)) +
      # geom_errorbar(aes(x = Covariate, y = contribution,
      #                   ymin=contribution-sd,ymax=contribution+sd),
      #               width=.1, position=position_dodge(.9)) +
      scale_fill_manual(values = c("random" =  "#CBCBCB",
                                   "envir" = "#FFA976",
                                   "habitat" = "#FFCF7A",
                                   "human" = "#9B7D9E")) +
      theme_minimal() +
      coord_flip() +
      labs(y = "Proportion in the variance explained", x = "") +
      theme(legend.position = "right") +
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12)) 
    ggsave(filename =  paste0(path_file,"/contributions_of_covariate_in_variance_",
                              save_name,".jpg"), width = 15, height = 10)
    
  } #END OF PLOT VARIANCE PARTITION
  
  
  #### Residual associations ####
  
  if(plot_residual_associations){
    cat("Check residual associations between responses... \n")
    
    # Check associations
    OmegaCor <- Hmsc::computeAssociations(model_fit_mcmc)
    supportLevel <- 0.95
    
    toPlot <- ((OmegaCor[[1]]$support>supportLevel) +
                 (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
    
    png(paste0(path_file,"/residual_associations_", save_name,".png"),
        width = 25, height = 15, units = "cm", res = 300)
    # par(mar = c(10,10,2,2))
    corrplot::corrplot(toPlot, method = "color",
                       col=colorRampPalette(c("blue","white","red"))(200),
                       tl.cex=.6, tl.col="black",
                       title=paste("random effect level:", model_fit_mcmc$rLNames[1]), 
                       mar=c(0,0,1,0))
    #residual associations among species can be generated by correlated responses to missing covariates,
    # or by ecological interactions.
    dev.off()
    
  } #END OF PLOT RESIDUAL ASSOCIATIONS
  
  
  #### Beta estimates ####
  
  if(plot_estimates){
    cat("Plot driver's estimates of the model... \n")
    
    ##### Support level  #####
    png(paste0(path_file,"/estimate_significance_", save_name,".png"),
        width = 25, height = 15, units = "cm", res = 300)
    par(mar = c(10,10,2,2))
    # Hmsc::plotBeta(model_fit_mcmc, post = postBeta, param = "Mean", supportLevel = 0.95)
    Hmsc::plotBeta(model_fit_mcmc, post = postBeta, param = "Sign", supportLevel = 0.95)
    dev.off()
    
    
    ##### ridges plot  #####
    S <- ggmcmc::ggs(mpost$Beta)
    
    # Extract one or some covariate(s)
    unique(S$Parameter)
    # drivers_to_plot <-  list(
    #   c("gravtot2", "n_fishing_vessels", "q95_5year_degree_heating_week"),
    #   c("effectiveness", "n_fishing_vessels")
    #   )
    

    for(drivers in drivers_to_plot){
      
      df <- S[grep(paste(drivers, collapse = "|"), S$Parameter),] |> 
        dplyr::mutate(Parameter = as.character(Parameter),
                      resp_name = sapply(strsplit(sapply(strsplit(Parameter, ","),
                                                         function(x) x[2]), " "), function(x) x[2]),
                      driver = sapply(strsplit(sapply(strsplit(sapply(strsplit(Parameter, ","),
                                                                      function(x) x[1]), "\\["), function(x) x[2]), " "),
                                      function(x) x[1])) |> 
        dplyr::filter(driver %in% drivers) #remove composed names as "coral_rubble"
      
      
      medians <- df  |> 
        dplyr::filter(driver == drivers[1]) |> 
        dplyr::group_by(resp_name)  |> 
        dplyr::summarise(median_value = median(value)) |> 
        dplyr::arrange(median_value)
      
      df <- df |> 
        dplyr::mutate(resp_name = factor(resp_name, levels = medians$resp_name))
      
      ggplot(df) +
        aes(y = resp_name, x = value,  fill = resp_name) +
        ggridges::geom_density_ridges(alpha=0.3)+ #, bandwidth = 0.005) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6)+
        hrbrthemes::theme_ipsum( axis_title_size = 0 ) +
        theme(
          legend.position="none",
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size = 15)
        ) +
        xlab(drivers) + ylab("Nature Contributions to People and Nature")+
        facet_wrap(~driver, ncol = length(df$driver)
                   , scales = "free_x"
        )
      
      ggsave(filename = paste0(path_file,"/posterior_distribution_of_estimates_", save_name,
                               paste(drivers, collapse = "-"), ".jpg"),
             width = 12, height = 8)
      
    }
    
    
    ##### Estimates Heatmap  #####
    S_df <- S |> 
      dplyr::mutate(Parameter = as.character(Parameter),
                    resp_name = sapply(strsplit(sapply(strsplit(Parameter, ","),
                                                       function(x) x[2]), " "), function(x) x[2]),
                    cov = sapply(strsplit(sapply(strsplit(sapply(strsplit(Parameter, ","),
                                                                 function(x) x[1]), "\\["), function(x) x[2]), " "),
                                 function(x) x[1]))
    
    mean_estimate <- S_df |> 
      dplyr::group_by(cov, resp_name) |> 
      dplyr::summarise(mean_posterior_distrib  = mean(value)) |> 
      dplyr::filter(cov != "(Intercept)")
    
    summary(mean_estimate$mean_posterior_distrib)
    
    ## Heat map of estimates
    palette <- rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788"))
    coef_matrix <- mean_estimate |> 
      tidyr::pivot_wider(names_from = "resp_name", values_from = "mean_posterior_distrib") |> 
      tibble::column_to_rownames("cov") |> 
      as.matrix()
    
    jpeg(filename = paste0(path_file,"/heatmap_estimates_", save_name,".jpg"),
         units = "cm", width = 20, height = 20, res = 400 )
    stats::heatmap(coef_matrix,  Colv=T,
                   hclustfun=function(x) hclust(x, method="ward.D2"), 
                   scale='none', col=palette, cexCol=0.6, margins = c(7,10))
    dev.off()
    
    
    
    ##### Estimate distribution #####
    # Classification of covariates
    covariates <- unique(mean_estimate$cov)
    envir <- c(grep("median", covariates, value = T))
    habitat <- c(grep("500m", covariates, value = T),
                 "depth", "algae", "coral", "Sand", "seagrass", "microalgal_mats",
                 "other_sessile_invert", "Rock", "coralline_algae", "coral_rubble")
    human <- setdiff(covariates, c(envir, habitat))
    
    # Classification of contributions
    cont_list <- unique(mean_estimate$resp_name)
    NP <- c("available_biomass", "selenium", "zinc", "omega_3" , "calcium",  "iron",                       
            "vitamin_A", "available_biomass_turnover", "NP_score")
    NN <- setdiff(cont_list, NP)
    
    #Resume data
    coeff_plot <- mean_estimate |> 
      dplyr::mutate(cov_class = ifelse(cov %in% envir, "environmental",
                                       ifelse(cov %in% habitat, "habitat", "human"))) |> 
      dplyr::mutate(contrib_class = ifelse(resp_name %in% NP, "Nature-for-People",
                                           "Nature-for-Nature"))
    
    # Plot estimate distribution
    
    plot_distri <- function(category = "Nature-for-Nature"){
      
      data <- dplyr::filter(coeff_plot, contrib_class ==  category)
      
      # order covariates
      df_ordered <- data |> 
        dplyr::group_by(cov) |> 
        dplyr::summarise(mean_estimate = median(mean_posterior_distrib)) |> 
        dplyr::arrange(mean_estimate) |> 
        dplyr::pull(cov)
      
      
      # Boxplot of estimates
      ggplot(data) +
        geom_vline(xintercept = 0, linetype = "dashed")+
        geom_boxplot(aes(y=factor(cov, levels = df_ordered), 
                         x = mean_posterior_distrib, fill = cov_class),
                     alpha = 0.7) +
        hrbrthemes::theme_ipsum() +
        # scale_fill_manual(values = c("#9467bd", "#ff7f0e", "#2ca02c"))+
        harrypotter::scale_fill_hp_d(option = "Ravenclaw") +
        # facet_wrap(~ contrib_class, scales = "free_x") +
        xlab("Regression coefficient estimates") +
        ylab("Predictors") +
        labs(title = paste("Estimate distribution on", category, " contributions"))+
        theme( axis.title.x = element_text(size = 10),
               axis.title.y = element_text(size = 10),
               axis.text.y = element_text(size = 10),
               plot.title = element_text(size = 12),
               legend.position = "bottom" 
        )
    }
    
    plot_distri(category = "Nature-for-Nature") + plot_distri(category = "Nature-for-People")
    ggsave(filename =paste0(path_file,"/mean_estimate_distribution_", save_name,".jpg"),
           width = 15, height = 8)
    
  
  } # END OF PLOT ESTIMATES
  
  
  
  # #### Partial graph ####
  # Gradient = Hmsc::constructGradient(model_fit_mcmc,
  #                                    focalVariable = "n_fishing_vessels")
  # predY = predict(model_fit_mcmc, 
  #                 XData = Gradient$XDataNew, 
  #                 studyDesign = Gradient$studyDesignNew,
  #                 ranLevels = Gradient$rLNew, 
  #                 expected=TRUE)
  # plotGradient(model_fit_mcmc, Gradient, pred=predY, measure="Y", showData = TRUE,
  #              index=13)
  # plotGradient(model_fit_mcmc, Gradient, pred=predY, measure="S", showData = TRUE)
  # 
  # 
  # png("figures/models/hmsc/partial_plot.png",
  #     width = 15, height = 15, units = "cm", res = 300)
  # plotGradient(model_fit_mcmc, Gradient, pred=predY, measure="Y", showData = TRUE,
  #              index=13)
  # dev.off()
  
  
  
  
  ##----------------------------- Residuals -----------------------------------
  
  if(check_residuals){
    cat("Plot predictive resituals of the model... \n")
    
    #Run predictions on all data
    X <- covariates_final
    studyDesign <- model_fit_mcmc$studyDesign
    ranLevels <- model_fit_mcmc$ranLevels
    response_long <- response |> 
      tibble::rownames_to_column("survey_id") |> 
      tidyr::pivot_longer(cols = -survey_id, names_to = "response",
                          values_to = "observation") 
    
    predY = predict(model_fit_mcmc, 
                    XData = X, 
                    studyDesign = studyDesign,
                    ranLevels = ranLevels, 
                    expected = T)
    
    PredY_mean <- as.data.frame(Reduce("+", predY) / length(predY))
    
    residuals <- response - PredY_mean
    
    compare_pred <- PredY_mean |> 
      tibble::rownames_to_column("survey_id") |> 
      tidyr::pivot_longer(cols = -survey_id, names_to = "response",
                          values_to = "prediction") |> 
      dplyr::left_join(response_long)
    
    compare_residuals <- residuals |> 
      tibble::rownames_to_column("survey_id") |> 
      tidyr::pivot_longer(cols = -survey_id, names_to = "response",
                          values_to = "residuals") |> 
      dplyr::left_join(response_long)
    
    
    #histograms of residuals
    distribution_plot(residuals, cols_plot = colnames(residuals))
    
    #Predictions vs observations
    ggplot(compare_pred)+
      geom_point(aes(x = observation, y = prediction, fill = response),
                 color = "grey40", alpha = 0.2, shape = 21) +
      hrbrthemes::theme_ipsum() +
      xlab("Observed contributions") + ylab("imputed")+
      geom_abline(slope = 1) + 
      ggpubr::stat_regline_equation(data = compare_pred,
                                    aes(x = observation, y = prediction, label = after_stat(rr.label)))   +
      facet_wrap(~response, scales = "free") +
      theme(legend.position="none", panel.spacing = unit(0.1, "lines"))
    ggsave(filename = paste0(path_file,"/Predictions_VS_Observations_", save_name,".jpg"),
           width = 15, height = 8)
    
    #residuals vs observations
    plot_interaction(compare_residuals, var_facet_wrap = "response", 
                     X_values = "observation", Y_values = "residuals")+
      geom_hline(yintercept = 0)
    ggsave(filename = paste0(path_file,"/Residual_VS_Observations_", save_name,".jpg"),
           width = 15, height = 8)
    
  } #END OF CHECK RESIDUALS
                
  
} ## END OF FUNCTION PLOT_HMSC_RESULT
                          

#' Run hmsc prediction.
#'
#' @param X_data dataframe of covariates
#' @param model_fit_mcmc hmsc fitted model
#' @param X_new_data same dataframe as X_data, but with different values to run conterfactual scenarios
#' @param new_surveys list of survey_id that have been changed
#'
#' @return save hmsc plots in figures/hmsc/file_name directory
#' @export
#'
#' @examples
#' 



run_hmsc_prediction <- function(X_data = X,
                                model_fit_mcmc = model_fit_mcmc,
                                new_surveys = NULL,
                                X_new_data = NULL
){
  ## Run prediction on original data
  studyDesign <- model_fit_mcmc$studyDesign
  ranLevels <- model_fit_mcmc$ranLevels
  
  
  predY <- predict(model_fit_mcmc, 
                   XData = X_data, 
                   studyDesign = studyDesign,
                   ranLevels = ranLevels, 
                   expected = T)
  
  PredY_mean <- as.data.frame(Reduce("+", predY) / length(predY))
  
  
  ## Run prediction in conterfactual scenarios
  if(!is.null(X_new_data)){
    
    predY_new_scenario <- predict(model_fit_mcmc, 
                                  XData = X_new_data, 
                                  studyDesign = studyDesign,
                                  ranLevels = ranLevels, 
                                  expected = F)
    
    PredY_new_scenario_mean <- as.data.frame(
      Reduce("+", predY_new_scenario) / length(predY_new_scenario))
    
    
    ## Observe changes
    
    changes_matrix <- PredY_new_scenario_mean - PredY_mean
    
    effective_change <- changes_matrix[new_surveys,] |> 
      tibble::rownames_to_column("survey_id") |> 
      tidyr::pivot_longer(cols = -survey_id, names_to = "index", values_to = "values")|> 
      dplyr::left_join( data.frame(survey_id = rownames(X_data),
                                   country = X_data$country)
      )
    
    }else{
      PredY_new_scenario_mean <- NULL
      effective_change <- NULL
    } # END of conterfactual scenarios
  
  list(predictions = PredY_mean, 
       new_scenario = PredY_new_scenario_mean,
       effective_change =effective_change)
} ## END OF FUNCTION PLOT_HMSC_RESULT
