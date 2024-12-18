###############################################################################'
##
##  This function contains all code necessary to run Hmsc models and plot results
##
## HMSC_function.R
##
## 07/05/2024
##
## Ulysse Flandrin
##
###############################################################################'


##--------------------------FIT HMSC MODELS-----------------------------------##


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
                          quadratic_effects = NULL,
                          random_factors = NULL,
                          nb_neighbours = NULL,
                          set_shrink = NULL,
                          test_null_model = NULL,
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
  if(!is.null(quadratic_effects)) cat("quadratic effects:", quadratic_effects, "\n")
  
  fixed_effects <- colnames(dplyr::select(X_data,
                                          -site_code,
                                          -longitude, -latitude, 
                                          -year,
                                          -country,
                                          -ecoregion,
                                          -realm
                                          # ,                     #country level covariates
                                          # -hdi,
                                          # -marine_ecosystem_dependency
                                          # -ngo,
                                          # -natural_ressource_rent
  ))
  
  linear_effects <- fixed_effects[!fixed_effects %in% quadratic_effects]
    
  if(!is.null(quadratic_effects)){
    formula <- as.formula(
      paste("~ ", 
            paste(
              c(linear_effects,
                paste0("poly(", quadratic_effects, ",degree = 2,raw = TRUE)")),
              collapse = "+")      
      ))
  }else{ 
    formula <- as.formula(paste("~ ", paste(linear_effects, collapse = "+")))
  }
  
  ## Test null model: only random levels
  if(!is.null(test_null_model)){ formula <- as.formula("~1")}
  
  # response_distribution <- rep("normal", ncol(Y_data))
  # response_distribution[colnames(Y_data) == "iucn_species_richness"] <- "poisson"
  
  
  ## Set random effects ##
  
  # setting model structure with spatial structure ‘Nearest Neighbour Gaussian Process (NNGP)’
  # (1) add noise to coordinates for them to differ sltitly
  noise_magnitude <- 0.00001 #noise in coordinates
  
  X <- X_data |>
    tibble::rownames_to_column("survey_id") |> 
    dplyr::rowwise() |> 
    dplyr::mutate(latitude = latitude + runif(1, -noise_magnitude, noise_magnitude),
                  longitude = longitude + runif(1, -noise_magnitude, noise_magnitude)) |>
    dplyr::ungroup() |> 
    tibble::column_to_rownames("survey_id")

  # (2) create a matrix with coordinates : xycoords is a matrix with 2 columns "x-coordinate","y-coordinate" and row names with spygen_code
  xycoords <- X |>
    dplyr::select(longitude, latitude)

  
  # (3) create spatial random effect and study design (rows in Y)
  rL.spatial = Hmsc::HmscRandomLevel(sData = xycoords,
                                  sMethod = 'NNGP',
                                  nNeighbours = nb_neighbours,
                                  # sMethod = 'Full', #too much memory needed...
                                  longlat = F) #Should be True but doesn't work after
  
  # Knots = Hmsc::constructKnots(xycoords, nKnots = 60) #regular nodes in areas with surveys points
  # rL.spatial = Hmsc::HmscRandomLevel(sData = xycoords,
  #                                 sMethod = "GPP",
  #                                 sKnot = Knots,
  #                                 longlat = F
  #                                 )
  
  
  rL.spatial = Hmsc::setPriors(rL.spatial, nfMin=1,nfMax=1) #max latent factor = 1
  
  studyDesign <- data.frame(sample_unit = as.factor(rownames(X)),
                            site = as.factor(X$site_code),
                            spatial = as.factor(rownames(X)),
                            year = as.factor(X$year),
                            country = as.factor(X$country),
                            ecoregion = as.factor(X$ecoregion))
  
  # (4) Set random levels
  rL_sample = Hmsc::HmscRandomLevel(units = unique(studyDesign$sample_unit))
  rL_site = Hmsc::HmscRandomLevel(units = unique(studyDesign$site))
  rL_year = Hmsc::HmscRandomLevel(units = unique(studyDesign$year))
  rL_ecoregion = Hmsc::HmscRandomLevel(units = unique(studyDesign$ecoregion))
  rL_country =  Hmsc::HmscRandomLevel(units = unique(studyDesign$country))
  
  ranLevels = list(sample_unit = rL_sample,
                   site = rL_site,
                   year = rL_year,
                   ecoregion = rL_ecoregion,
                   country = rL_country, 
                   spatial = rL.spatial)
  
  # (Set shrinkage
  if(!is.null(set_shrink)){
    rL_country=Hmsc::setPriors(rL_country, a1=set_shrink, a2=set_shrink)
    rL_sample=Hmsc::setPriors(rL_sample, a1=set_shrink, a2=set_shrink)
    rL_site=Hmsc::setPriors(rL_site, a1=set_shrink, a2=set_shrink)
    rL_ecoregion=Hmsc::setPriors(rL_ecoregion, a1=set_shrink, a2=set_shrink)
    rL_year=Hmsc::setPriors(rL_year, a1=set_shrink, a2=set_shrink)
  }
  #)
  
  
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
  if(ncol(Y_data) == 1){ engine="R"}else{engine="HPC"}
  init_obj = Hmsc::sampleMcmc(model_fit,
                              samples=nSamples,
                              thin=thin,
                              transient=transient, 
                              nChains=nChains,
                              verbose=verbose,
                              nParallel = nChains,
                              engine=engine) # HPC : try to use the GPU under python command: on github version only (v3.1-2)
  
  
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
    if(ncol(Y_data) == 1){
      saveRDS(jsonify::to_json(init_obj), file=post_file_path)
      }else{
      cat("-------------- Run Hmsc with python --------------\n")
      system2(python, python_cmd_args)
    }
  }
  
} # END OF hmsc_function



##--------------------------FIT HMSC CROSSVALIDATION-----------------------------------##


#' HMSC function.
#'
#' @param k_fold Number of fold in crossvalidations
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



fit_hmsc_crossvalidation <- function(k_fold = 5, 
                                 nSamples,
                                 thin,
                                 nChains,
                                 verbose,
                                 transient,
                                 Y_data,
                                 X_data,
                                 response_distribution,
                                 quadratic_effects = NULL,
                                 random_factors = NULL,
                                 nb_neighbours = NULL,
                                 set_shrink = NULL,
                                 test_null_model = NULL,
                                 name,
                                 run_python = TRUE,
                                 save_path = here::here("outputs/models/hmsc")
){
  
  save_path_cv <- file.path(save_path, "cross_validation/")
  
  ### Data sets ###
  folds <- caret::groupKFold(X_data$site_code, k = k_fold)
  
  datasets <- lapply(c(1:k_fold), FUN = function(i){
    
    train_indices <- folds[[i]]

    train <- Y_data[train_indices,] 
    test <- Y_data[-train_indices,]

    list(train, test)
  })
  
  
  ### Train models ###
  train_model <- parallel::mclapply(1:k_fold, 
                                    mc.cores = 5,
                                    function(i){
      # Data
      fold <- datasets[[i]]
      Y_train <- fold[[1]]
      Y_test <- fold[[2]]
      name_cv <- paste0( "CV",i, "_", name)
      
      X_train <- X_data[rownames(Y_train),]
      X_test <- X_data[rownames(Y_test),]
      
      # Train model
      hmsc_function(nSamples, thin, nChains, verbose, transient,
                    Y_data = Y_train,
                    X_data = X_train,
                    response_distribution, quadratic_effects,random_factors,
                    nb_neighbours, set_shrink, test_null_model, name_cv,
                    run_python = run_python, save_path_cv)

      list(X_train, Y_train, X_test, Y_test)
      
      }) # END OF TRAIN CROSS-VAL MODEL
      
  
  
  ## CLEAN FOLDER ###
  save_init <- paste0(save_path_cv,"init_multi/")
  save_out <- paste0(save_path_cv, "out_multi/")
  localDir <- paste0(save_path_cv, "multivariate/")
  
  model_name <- paste0(name, "_", paste(nChains, "chains",
                                        thin, "thin",
                                        nSamples, "samples",
                                        sep = "_"))
  
  paths <- c(save_init, localDir, save_out)
  for(path in paths){

    dir.create(paste0(path, model_name), showWarnings = F)
    
    list_files <- list.files(path)
    crossval <- list_files[grep(paste0("CV.*", model_name), list_files)]
    
    #deplace all files and rename them
    for (file in crossval) {
      destination_file <- file.path(paste0(path, model_name),
                                    paste0(stringr::str_extract(file, "CV[0-9]+"), ".rds"))
      
      file.rename(file.path(path, file),
                  destination_file
                  )
      } # end of for each cross val
    } # end for each path
  
  
  #Save X and Y data associated with crossval
  save(train_model, file = file.path(paste0(save_out, model_name),
                                     "dataset_crossvalidation.Rdata"))
  
} # END OF fit_hmsc_crossvalidation





##-----------------------------PLOT HMSC RESULTS-----------------------------------##

#' Plot results Hmsc.
#'
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



plot_hmsc_result <- function(metadata = metadata,
                             #response = response,
                             file_name = file_name,
                             path = path,
                             concatenate_chains = F,
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
                               c("n_fishing_vessels", "gravity", "gdp", "neartt"))
                            ){
  
  color_grad = rev(c("#A50026", "#D73027", "#FDAE61", "#FEE090","white",
                 "#D8DAEB", "#B2ABD2", "#8073AC", "#542788"))
  color_grad_soft <- colorspace::desaturate(color_grad, amount = 0.3)
  RdBu = rev(c("#67001F", "#B2182B", 
           "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", 
           "#92C5DE", "#4393C3", "#2166AC", "#053061"))
  
  #Coastline
  coast <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')
  
  # Contributions classifications
  grp_NN_NP <- data.frame(
    contribution = c("actino_richness","functional_distinctiveness",
                      "iucn_species_richness" ,"mean_endemism",
                      "evolutionary_distinctiveness","functional_entropy",
                      "phylogenetic_entropy","herbivores_biomass",
                      "invertivores_biomass",  "piscivores_biomass",
                      "trophic_web_robustness", "mean_trophic_level",
      
                      "public_attention", "aesthetic",
                      "available_biomass", "selenium",
                      "zinc",   "omega_3", "calcium",
                      "iron","vitamin_A", "available_biomass_turnover"),
    group = c(rep("NN", 12), 
              rep("NC", 2),
              rep("NS", 8)))
  
  ##----Source packages and functions----
  
  library(ggplot2)
  library(patchwork)
  library(stringr)
  
  source("R/evaluation_prediction_model.R")
  
  # Hmsc::computeVariancePartitioning creates bug in Rstudio because it can fill the console of warning.
  # to clear it you can run 'cat("\014")' or make Ctrl + L , to continue.
  # Otherwise, use the modified following function of computeVariancePartitioning
  source("R/HMSC_computeVariancePartitioning.R")
  
  # PATHS
  save_init <- file.path(path, "init_multi/")
  save_out <- file.path(path, "out_multi/")
  localDir <- file.path(path, "multivariate/")
  
  ##---- Import initial object ----
  cat("Import HMSC fitted model... \n")
  
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
  if(concatenate_chains){
    chainList = vector("list", nChains)
    importFromHPC <- vector("list", nChains+1)
    for(cInd in 1:nChains){
      chain_file_path = file.path(
        paste0(save_out,"output_", file_name, "/chain_", cInd, ".rds"))
      chainList[[cInd]] = jsonify::from_json(readRDS(file = chain_file_path)[[1]])[[1]]
    }
  }else{
    all_chains <- readRDS(file = paste0(save_out, paste0("output_", file_name)))[[1]]
    importFromHPC <- jsonify::from_json(all_chains)
    chainList <- importFromHPC[1:nChains]
    ##Result model
    cat(sprintf("fitting time %.1f h\n", importFromHPC[[nChains+1]] / 3600))
    
  }
  ### /!\ CHECK THE STRUCTURE OF THE RESULT: MATRICES IN PYTHON VS LIST IN R ###
  
  ##Export and merge chains
  model_fit_mcmc <- Hmsc::importPosteriorFromHPC(model_fit,
                                                 chainList, 
                                                 nSamples, 
                                                 thin, 
                                                 transient)
  
  ## Initial data
  Y_data <- as.data.frame(model_fit_mcmc$Y)
  X_data <- model_fit_mcmc$XData
  metadata_model <- metadata[rownames(X_data),]
  
  
  ## Estimates for each chains
  mpost <- Hmsc::convertToCodaObject(model_fit_mcmc)
  
  postBeta <- Hmsc::getPostEstimate(model_fit_mcmc, parName = "Beta")
  
  S_raw <- ggmcmc::ggs(mpost$Beta)
  
  #Arrange S table
  S_arranged <- S_raw |>
    dplyr::mutate(
      # Extract the covariates
      covariate = dplyr::case_when(
        # Handle polynomial terms like poly(reef_500m, degree = 2)
        grepl("poly\\(", Parameter) ~ str_replace(
          str_replace(str_extract(Parameter, "poly\\(([^,]+),"), "poly\\(", ""), 
          "\\,", ""),
        # Handle with Intercept
        grepl("Intercept", Parameter) ~ str_extract(Parameter, "(?<=\\()[a-zA-Z]+(?=\\))"),
        # Handle other cases without poly
        TRUE ~ str_extract(Parameter, "[a-zA-Z0-9_]+(?=\\s*\\()")
      ),
      # Extract the responses
      response = dplyr::case_when(
        # Handle polynomial terms like poly(reef_500m, omega_3)
        grepl("poly\\(", Parameter) ~ str_replace(
          str_replace(
            str_extract(Parameter, ",\\s*([a-zA-Z0-9_]+)\\s*\\("), ",\\s*", ""), 
          " \\(", ""),
        # Handle other cases without poly
        TRUE ~ str_extract(Parameter, "(?<=,\\s)[a-zA-Z0-9_]+")
      ),
      #Extract the polynomial degree
      degree = dplyr::case_when(
        grepl("poly\\(", Parameter) ~ str_replace(str_extract(Parameter, "\\)\\d+"), "\\)", ""),
        TRUE ~ ""
      ))|>
    dplyr::mutate(covariate = dplyr::case_when(
      grepl("poly\\(", Parameter) ~ paste0(covariate,"_Deg",degree),
      TRUE ~ covariate),
      Parameter = paste0(covariate, "-", response))
  
  
  
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
  # Beta are the estimates of covariables on the responses = fixed effects
  # Omega is the matrix of species-to-species residual covariances.
  png(paste0(path_file,"/convergence_estimates_effective_size_", save_name,".png"),
      width = 20, height = 20, units = "cm", res = 200)
  par(mfrow=c(2,2), mar = c(4,4,4,4))
  hist(coda::effectiveSize(mpost$Beta), main="fixed_effects: ess(beta)", breaks = 30)
  abline(v= nSamples*nChains, col = "red4", lty = 2, lwd = 3)
  hist(coda::gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)", breaks = 30)
  abline(v= 1, col = "red4", lty = 2, lwd = 3)
  hist(coda::effectiveSize(mpost$Omega[[1]]), main="random_effects: ess(omega)", breaks = 30) #also check associations between y
  abline(v= nSamples*nChains, col = "red4", lty = 2, lwd = 3)
  hist(coda::gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf, main="psrf(omega)", breaks = 30)
  abline(v= 1, col = "red4", lty = 2, lwd = 3)
  dev.off()
  
  
  
  ### Chains convergence ###
  # install.packages("ggmcmc")
  # library("ggmcmc")
  # see help at : http://xavier-fim.net/post/using_ggmcmc/
  
  cov <- "median_5year_analysed_sst"
  S <- S_arranged |> dplyr::filter(grepl(cov, covariate))
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
    
    ### AIC and computing time ###
    comput_time <- round(importFromHPC[[nChains+1]] / 3600, 1)
    AIC <- round(Hmsc::computeWAIC(model_fit_mcmc), 2)
    
    ### Explanatory  power ###
    preds <- Hmsc::computePredictedValues(model_fit_mcmc)
    MF <- Hmsc::evaluateModelFit(hM=model_fit_mcmc, predY=preds)
    
    MF_table <- as.data.frame(MF)
    MF_table$responses <- model_fit_mcmc[["spNames"]]
    save(MF_table, file = paste0(path_file,"/explanatory_power_data.Rdata"))
    
    png(paste0(path_file,"/explanatory_power_", save_name,".png"),
        width = 20, height = 13, units = "cm", res = 300)
    par(mfrow=c(1,2), mar = c(4,4,10,4))
    hist(MF$R2, xlim = c(0,1), main=paste0("Mean R2 = ", round(mean(MF$R2),2)))
    hist(MF$RMSE, xlim = c(0,1), main=paste0("Mean RMSE = ", round(mean(MF$RMSE),2)))
    mtext(paste0("WAIC = ", AIC, "     Computing Time: ", comput_time, "h"),
          side = 3, line = -2.5, outer = TRUE, cex = 1.1, font = 2)
    dev.off()
    
    
    
    ### R² per ecoregions ###
    predY <- predict(model_fit_mcmc, 
                     XData = X_data, 
                     studyDesign = model_fit_mcmc$studyDesign,
                     ranLevels = model_fit_mcmc$ranLevels)
    
    PredY_mean <- as.data.frame(Reduce("+", predY) / length(predY))
    
    predictions <- PredY_mean |> 
      tibble::rownames_to_column("id") |> 
      tidyr::pivot_longer(cols = -id, names_to = "contributions", values_to = "predicted") |> 
      dplyr::left_join(
        tidyr::pivot_longer(
          tibble::rownames_to_column(Y_data, "id"),
          cols = -id, names_to = "contributions", values_to = "observed")
      ) |> 
      dplyr::left_join(
        tibble::rownames_to_column(metadata_model, "id")
      ) 
    
    
    ### R² per ecoregions ###
    sf::sf_use_s2(FALSE)
    ecoregions <- sf::st_read(here::here("data/raw_data/MEOW/meow_ecos.shp"))
   
    points_sf <- sf::st_as_sf(predictions, 
                              coords = c("longitude", "latitude"), crs = 4326) |> 
      dplyr::select(id, geometry) |> 
      unique()
    
    points_in_ecoreg <- sf::st_join(points_sf, ecoregions, join = sf::st_within)
    
    predictions_ecoreg <- predictions |> 
      dplyr::left_join(points_in_ecoreg) |> 
      dplyr::group_by(ECOREGION, contributions) |> 
      dplyr::summarise(R_squared = cor(observed, predicted)^2)|> 
      dplyr::group_by(ECOREGION) |> 
      dplyr::summarise(R_squared_ecoreg = mean(R_squared, na.rm =T))
    
    ecoregions_Rsqu <- ecoregions |> 
      dplyr::left_join(predictions_ecoreg)
    
    ggplot() +
      geom_sf(data = ecoregions_Rsqu, aes(fill = R_squared_ecoreg), color = "grey50")+
      geom_sf(data = coast, fill = "grey60", color = "grey40")+
      coord_sf(c(-180,180), c(-40, 40) , expand = FALSE) +
      scico::scale_fill_scico(palette = "davos", direction = -1, 
                              na.value = "lightgrey", end = 0.7 )+
      theme_minimal()+
      labs(title = "",x="", y= "", fill = "Mean R²") +
      theme(plot.title = element_text(size=15, face="bold"),
            legend.text = element_text(size=13),
            legend.title = element_text(size=15),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0,0,0,0), units = , "cm")
      )
    
    ggsave(filename = paste0(path_file,"/Explanatory_power_per_ecoregion_",
                             save_name,".jpg"), width = 15, height = 5)
    
    
    
    # ### R² per contry ###
    # predictions_country <- predictions |> 
    #   dplyr::group_by(country, contributions) |> 
    #   dplyr::summarise(R_squared = cor(observed, predicted)^2) |> 
    #   dplyr::group_by(country) |> 
    #   dplyr::summarise(R_squared_country = mean(R_squared, na.rm =T)) |> 
    #   dplyr::mutate(country = dplyr::recode(country,
    #                                         "United States" = "United States of America",
    #                                         "Caribbean Sea" ="Colombia")) |> 
    #   dplyr::rename(name_en = "country") 
    # 
    # coast_rsquared <- coast|> 
    #   dplyr::left_join(predictions_country)
    # 
    # ggplot() +
    #   geom_sf(data = ecoregions_Rsqu, aes(fill = R_squared_ecoreg), color = "grey30")+
    #   geom_sf(data = coast_rsquared, aes(fill = R_squared_country), color = "grey30")+
    #   coord_sf(c(-180,180), c(-36, 33) , expand = FALSE) +
    #   scico::scale_fill_scico(palette = "davos", direction = -1, 
    #                           na.value = "lightgrey", end = 0.7 )+
    #   theme_minimal()+
    #   labs(title = "",x="", y= "", fill = "Mean R²") +
    #   theme(plot.title = element_text(size=15, face="bold"),
    #         legend.text = element_text(size=13),
    #         legend.title = element_text(size=15),
    #         axis.text.x = element_blank(),
    #         axis.ticks.x = element_blank(),
    #         plot.margin = unit(c(0,0,0,0), units = , "cm")
    #   )
    # 
    # ggsave(filename = paste0(path_file,"/Explanatory_power_per_country_",
    #                          save_name,".jpg"), width = 15, height = 5)
   
    
  } #END OF PLOT EXPLANATORY POWER
  
  ##----------------------------- Covariates importance -----------------------------------
  
  if(plot_variance_partitioning){
    cat("Complute variance partitioning of the model... \n")
    
    preds <- Hmsc::computePredictedValues(model_fit_mcmc)
    MF <- Hmsc::evaluateModelFit(hM=model_fit_mcmc, predY=preds)
    
    #### Variance partitioning ####
    VP <- computeVariancePartitioning(model_fit_mcmc)

    ##### Preping VP table #####
    VP_long <- as.data.frame(VP[["vals"]]) |> 
      tibble::rownames_to_column(var = "Covariate") |> 
      tidyr::pivot_longer(
        cols = - Covariate,
        names_to = "Response",
        values_to = "Value"
      )
    
    #classify covariates
    human <- c("gdp", "gravity", "protection_status", "natural_ressource_rent",
               # "protection_status2",
               "neartt","n_fishing_vessels", "hdi", "marine_ecosystem_dependency")
    
    habitat <- c("depth", #"algae",
                 "coral", "Sand", #"seagrass", "microalgal_mats",
                 "other_sessile_invert", "Rock", "coralline_algae", "coral_rubble",
                 "Back_Reef_Slope_500m", "coral_algae_500m",  "Deep_Lagoon_500m", 
                 "Inner_Reef_Flat_500m", "Microalgal_Mats_500m", "Patch_Reefs_500m",   
                 "Plateau_500m", "Reef_Crest_500m", "Reef_Slope_500m", 
                 "Rock_500m", "Rubble_500m", "Sand_500m", "Seagrass_500m", 
                 "Sheltered_Reef_Slope_500m", "Terrestrial_Reef_Flat_500m")
    
    envir <-  c("median_5year_analysed_sst", "median_5year_chl", 
                "median_5year_degree_heating_week", #"q05_5year_ph",
                "median_5year_so_glor", "median_1year_degree_heating_week",
                "median_7days_degree_heating_week", #"median_7days_analysed_sst",
                "median_7days_chl", #"q95_1year_degree_heating_week", 
                "median_5year_ph",
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
      dplyr::mutate(
        Covariate = dplyr::case_when(
          grepl("poly\\(", Covariate) ~ str_replace(
            str_replace(str_extract(Covariate, "poly\\(([^,]+),"), "poly\\(", ""), 
            "\\,", ""),
          TRUE ~ Covariate)
        # 
        # degree = dplyr::case_when(
        #   grepl("poly\\(", Covariate) ~ str_replace(str_extract(Covariate, "\\)\\d+"), "\\)", ""),
        #   TRUE ~ "")
        )|>
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
    
    
    # VP_long <- VP_long |> 
    #   dplyr::mutate(category = dplyr::case_when(
    #     Covariate %in% human ~ "human",
    #     Covariate %in% habitat ~ "habitat",
    #     Covariate %in% envir ~ "envir",
    #     Covariate %in% random ~ "random",
    #     TRUE ~ NA_character_
    #   )) |> 
    #   dplyr::rowwise() |> 
    #   dplyr::mutate(color = if (!is.na(category)) cat_colors[[category]][match(Covariate, get(category))] else NA_character_) |> 
    #   dplyr::ungroup()
    

    # Order covariates
    VP_long$Covariate <- forcats::fct_relevel(VP_long$Covariate, 
                                              c(human, habitat, envir, random))
    
    #Rename covariates
    VP_long <- VP_long |> 
      dplyr::mutate(Covariate = dplyr::recode(Covariate,
                                              "median_5year_analysed_sst" = "SST_5years",
                                              "median_5year_degree_heating_week" = "DHW_5years",
                                              "median_5year_so_glor" = "Salinity_5_years",
                                              "median_1year_degree_heating_week" = "DHW_1year",
                                              "median_7days_degree_heating_week" = "DHW_7days",
                                              "median_7days_chl" = "Chlorophyll_7days",
                                              "median_5year_ph" = "pH_5years",
                                              "q95_5year_degree_heating_week" = "DHW_quartile95_5years",
                                              "median_5year_chl" = "Chlorophyll_5years"))
    
    # Mean contribution of covariates
    mean_contrib <-  VP_long |> 
      dplyr::group_by(Covariate, color, category) |> 
      dplyr::summarise(Value = mean(Value)) |> 
      dplyr::mutate(Response = "Mean contribution")
    
    add_gap <- mean_contrib|> 
      dplyr::mutate(Response = "", Value = 0)
    
    VP_long <- VP_long |>
      dplyr::bind_rows(mean_contrib) |>
      dplyr::bind_rows(add_gap)
    
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
      dplyr::filter(!Response %in% c("","Mean contribution")) |> 
      dplyr::group_by(category, Response) |> 
      dplyr::summarise(prop_variance = sum(Value)) |> 
      dplyr::group_by(category) |> 
      dplyr::summarise(sd = sd(prop_variance),
                       prop_variance = mean(prop_variance))
    
    
    # VP absolute
    variance_explained <- data.frame(Response = c(model_fit_mcmc[["spNames"]],"", "Mean contribution"),
                                     R2 = c(MF$R2,0,1))
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
    VP_random <- VP_long |> 
      dplyr::filter(!Response %in% c("","Mean contribution")) |> 
      dplyr::filter(category == "random") |> 
      dplyr::group_by(Covariate, color) |> 
      dplyr::summarise(mean_imp = mean(Value)) 
    
    ggplot(VP_aggregated) +
      geom_col(aes(x = reorder(category, prop_variance),
                   y = prop_variance, fill = category),
               color = "grey20") +
      
      #Add precision on the random bar
      annotate("rect", 
               xmin = which(levels(reorder(VP_aggregated$category,
                            VP_aggregated$prop_variance)) == "random") - 0.45,
               xmax = which(levels(reorder(VP_aggregated$category,
                            VP_aggregated$prop_variance)) == "random") + 0.45, 
               ymin = 0, 
               ymax = VP_random$mean_imp[2], 
               fill = VP_random$color[2])+
      annotate("text", label = stringr::word(VP_random$Covariate[2], sep = " ", 2),
               x = which(levels(reorder(VP_aggregated$category,
                                    VP_aggregated$prop_variance)) == "random"),
               y = VP_random$mean_imp[2]/2)+
      
      annotate("rect", 
               xmin = which(levels(reorder(VP_aggregated$category, 
                          VP_aggregated$prop_variance)) == "random") - 0.45,
               xmax = which(levels(reorder(VP_aggregated$category, 
                          VP_aggregated$prop_variance)) == "random") + 0.45, 
               ymin = VP_random$mean_imp[2], 
               ymax = VP_random$mean_imp[2] + VP_random$mean_imp[1], 
               fill = VP_random$color[1], alpha = 0.5)+
      annotate("text", label = stringr::word(VP_random$Covariate[1], sep = " ", 2),
               x = which(levels(reorder(VP_aggregated$category,
                                        VP_aggregated$prop_variance)) == "random"),
               y = VP_random$mean_imp[2] + VP_random$mean_imp[1]/2.5)+
      
      #plot parameters
      geom_errorbar(aes(x = category, y = prop_variance,
                        ymin=prop_variance-sd, ymax=prop_variance+sd),
                    width=.1, position=position_dodge(.9),
                    color = "grey20") +
      scale_fill_manual(values = c("random" =  "#CBCBCB",
                                   "envir" = "#FFA976",
                                   "habitat" = "#FFCF7A",
                                   "human" = "#9B7D9E")) +
    
      scale_x_discrete(labels = c("random" = "Random Effects", 
                                  "envir" = "Environment", 
                                  "habitat" = "Habitat", 
                                  "human" = "Human")) +
      theme_minimal() +
      coord_flip() +
      labs(y = "Proportion in the variance explained", x = "") +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12)) 
    ggsave(filename = paste0(path_file,"/variance_explained_per_category_",
                             save_name,".jpg"), width = 15, height = 10)
    
    
    
    ## Relative variance partitioning
    ggplot(dplyr::filter(VP_long,!Response %in% c("","Mean contribution"))) +    
      aes(x = Response, y = Value, fill = Covariate)+
      geom_bar(stat = "identity", position = "stack", 
               color = "black",linewidth = 0.1) +
      scale_fill_manual(values = setNames(VP_long$color, VP_long$Covariate),
                        labels = unique(VP_long$labels)) +
      labs(title = "", x = "", y = "", fill ="Covariates") +
      theme_classic(base_size = 11,
                    base_line_size = 0.1) +
      theme(
        axis.text.x = element_text(size=12, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=12),
        legend.position = "right",
        legend.title = element_text(size = 15, hjust = 0),
        legend.text = element_text( size = 10))+
      guides(fill = guide_legend(ncol = 2)) + 
      geom_text(data = dplyr::filter(VP_long, Response != "Mean contribution" & Value > 0.04),
                aes(y = mid_y, label = Symbol), size = 3, color = "black") 

    ggsave(filename =  paste0(path_file,"/variance_partitioning_", save_name,".jpg"),
           width = 15, height = 10)
    
    
    
    ## VP with absolute value of variance explained
    variance_explained 
    
    VP_absolute_aggregated <- VP_long_absolute |> 
      dplyr::filter(Response != "Mean contribution") |> 
      dplyr::group_by(category, Response) |> 
      dplyr::summarise(prop_variance = sum(Value)) |> 
      dplyr::group_by(category) |> 
      dplyr::summarise(sd = sd(prop_variance),
                       prop_variance = mean(prop_variance))
    
    VP_long_absolute[VP_long_absolute$Response == "Mean contribution", "R2"] <- -1
    
    
    
    ## Plot Figure 1
    VP_plot_absolute <- 
      ggplot(dplyr::filter(VP_long_absolute, !Response %in% c("", "Mean contribution")))+
      aes(x = reorder(Response,-R2), y = Value, fill = Covariate) +
      geom_bar(stat = "identity", position = "stack",
               color = "black", linewidth = 0.1) +
      scale_fill_manual(values = setNames(VP_long_absolute$color, VP_long_absolute$Covariate),
                        labels = unique(VP_long_absolute$labels)) +
      labs(title = "", x = "", y = "", fill ="Covariates") +
      theme_classic(base_size = 11,
                    base_line_size = 0.1) +
      theme(
        axis.text.x = element_text(size=15, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=15),
        legend.position = "bottom",
        legend.title = element_text(size = 20, hjust = 0.5, vjust= 5, angle = 90),
        legend.text = element_text( size = 14),
        legend.key.spacing.y = unit(0, units = "cm"))+
      guides(fill = guide_legend(nrow = 9)) +
      geom_text(data = dplyr::filter(VP_long_absolute, Value > 0.02 &
                                       !Response %in% c("", "Mean contribution")),
                aes(y = mid_y, label = Symbol), size = 4, color = "black") 
    
    # #Custom mean bar
    # geom_rect(aes(xmin = length(unique(Response))-1 - 0.5, # Hide gap bar
    #               xmax = length(unique(Response))-1 + 0.5,
    #               ymin = -Inf, ymax = Inf),
    #           inherit.aes = FALSE, fill = "white", color = NA)+ 
    # geom_rect(aes(xmin = length(unique(Response))- 0.5, #higlight mean bar
    #               xmax = length(unique(Response))+ 0.5,
    #               ymin = 0, ymax = 1.001),
    #           inherit.aes = FALSE, fill = NA, color = "black", linewidth = 2)+
    # # Increase font size of "Average effects"
    # scale_x_discrete(labels = c("Mean contribution" = expression(bold("Average effects"))))+

    
    ## Mean importance column
    VP_plot_mean_imp <- 
      ggplot(dplyr::filter(VP_long_absolute, Response=="Mean contribution"))+
      aes(x = Response, y = Value, fill = Covariate) +
      geom_bar(stat = "identity", position = "stack",
               color = "black", linewidth = 0.1) +
      scale_fill_manual(values = setNames(VP_long_absolute$color, VP_long_absolute$Covariate),
                        labels = unique(VP_long_absolute$labels)) +
      theme_classic(base_size = 11, base_line_size = 0.1) +
      scale_x_discrete(labels = c("Mean contribution" = "Average importance"))+
      theme(
        axis.text.x = element_text(size=18, angle = 50, face= "bold", 
                                   hjust = 1, vjust = 1),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),)+
      geom_text(data = dplyr::filter(VP_long_absolute, Value > 0.02 &
                                       Response == "Mean contribution"),
                aes(y = mid_y, label = Symbol), size = 4, color = "black")+
      geom_rect(aes(xmin = 1- 0.5, 
                    xmax = 1+ 0.5,
                    ymin = 0, ymax = 1.003),
                inherit.aes = FALSE, fill = NA, color = "black", linewidth = 1.8)
  
    
    #Merge the plot
    library(patchwork)
    VP_plot_absolute + plot_spacer() + VP_plot_mean_imp + plot_layout(widths = c(20, 0.2, 1))
    ggsave(filename =  paste0(path_file,"/variance_partitioning_absolute_values_", save_name,".jpg"),
           width = 18, height = 16)
    
    
    
    ## Mean contributions of covariates
    covariate_contrib <- VP_long |> 
      dplyr::filter(Response != "Mean contribution") |> 
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
                                   "human" = "#9B7D9E"),
                        labels = c("envir" = "Environment", 
                                   "habitat" = "Habitat", 
                                   "human" = "Human")) +
      theme_minimal() +
      coord_flip() +
      labs(y = "Proportion in the variance explained", x = "") +
      theme(legend.position = "right") +
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.title = element_text(size = 15, hjust = 0),
            legend.text = element_text( size = 12)) 
    ggsave(filename =  paste0(path_file,"/contributions_of_covariate_in_variance_",
                              save_name,".jpg"), width = 15, height = 10)
    
  } #END OF PLOT VARIANCE PARTITION
  
  
  #### Residual associations ####
  
  if(plot_residual_associations){
    cat("Check residual associations between responses... \n")
    
    # Check associations
    OmegaCor <- Hmsc::computeAssociations(model_fit_mcmc)
    supportLevel <- 0.95
    
    list_corr <- list()
    for(i in 1:length(OmegaCor)){
      #i=1
      resid <- OmegaCor[[i]]$mean
      cor_contrib <- cor(Y_data)[colnames(resid), colnames(resid)]
      
      mantel_test <- vegan::mantel(resid, cor_contrib)
      
      paste0("Mantel statistic:\n r = ", round(mantel_test[["statistic"]], 3),
             ", p-val = ", mantel_test[["signif"]])
      
      toPlot <- ((OmegaCor[[i]]$support>supportLevel) +
                   (OmegaCor[[i]]$support<(1-supportLevel))>0)*OmegaCor[[i]]$mean
      
      png(paste0(path_file,"/residual_associations_", model_fit_mcmc$rLNames[i], save_name,".png"),
          width = 25, height = 15, units = "cm", res = 300)
      # par(mar = c(10,10,2,2))
      corrplot::corrplot(toPlot, method = "color",
                         col=colorRampPalette(RdBu)(200),
                         tl.cex=.6, tl.col="black",
                         order = "hclust",
                         title=paste("random effect level:", model_fit_mcmc$rLNames[i]), 
                         mar=c(0,0,1,0))
      #residual associations among species can be generated by correlated responses to missing covariates,
      # or by ecological interactions/similar life story and constraints.
      mtext( paste0("Mantel statistic with initial correlation\n of contributions: r = ",
                    round(mantel_test[["statistic"]], 2),
                    ", p-val = ", mantel_test[["signif"]]),
             side = 3, adj=0, line = 1, cex = 0.8)
      dev.off()
      
      list_corr[[i]] <- toPlot
      names(list_corr)[i] <- model_fit_mcmc$rLNames[i]
    }
    
    # etaPost=getPostEstimate(model_fit_mcmc, "Eta")
    # lambdaPost=getPostEstimate(model_fit_mcmc, "Lambda")
    # Hmsc::biPlot(model_fit_mcmc, etaPost = etaPost, lambdaPost = lambdaPost, 
    #              factors = c(1,2))
    
    
    ## Plot 2 random levels together
    corr_1 <- list_corr[[1]]
    corplot_1 <- corrplot::corrplot(corr_1,order = "hclust", mar=c(0,0,0,0))
    
    corr_2 <- list_corr[[2]]
    order_name <- colnames(corplot_1[["corr"]]) ; n <- nrow(corr_1) 
    corr_2_ordered <- corr_2[match(order_name, rownames(corr_2)), 
                             match(order_name, colnames(corr_2))]
    
    png(paste0(path_file,"/residual_associations_both_levels", save_name,".png"),
        width = 25, height = 25, units = "cm", res = 300)
      corrplot::corrplot(corr_1, method = "color", diag = T,
                         col=colorRampPalette(color_grad_soft)(200),
                         tl.cex=.6, tl.col="black", tl.srt = 60,
                         order = "hclust",
                         #title=paste("random effect level:", model_fit_mcmc$rLNames[1]), 
                         mar=c(4,5,5,0))   
      corrplot::corrplot(corr_2_ordered, method = "color", type = "lower", diag = F,
                         col=colorRampPalette(color_grad_soft)(200),
                         tl.pos='n', tl.col="black",order = "original",
                         cl.pos = 'n',
                         #title=paste("random effect level:", model_fit_mcmc$rLNames[2]), 
                         add= T)
      polygon(x = c(0.5, n + 0.5), y = c(n + 0.5, 0.5), border = "white",lwd = 4)
      polygon(x = c(0.6, n + 0.5, n + 0.5), y = c(n + 0.5, n + 0.5, 0.6), 
              border = "grey10", #border = color_grad[1], 
              lwd = 3)
      polygon(x = c(0.5, 0.5, n + 0.4), y = c(0.5, n + 0.4, 0.5), 
              border = "grey50", #border = color_grad[length(color_grad)], 
              lwd = 3)
      
      mtext(paste("random effect level:", model_fit_mcmc$rLNames[1]), side = 3, 
            line = 0, cex = 1.5, font = 2, 
            col ="grey10" )# color_grad[1]) 
      mtext(paste("random effect level:", model_fit_mcmc$rLNames[2]), side = 2,
            line = 0, cex = 1.5, font = 2, las = 3,
            col = "grey50")# color_grad[length(color_grad)]) 
 
    dev.off()

    
  } #END OF PLOT RESIDUAL ASSOCIATIONS
  
  
  #### Beta estimates ####
  
  if(plot_estimates){
    cat("Plot driver's estimates of the model... \n")
    
    ##### Support level  #####
    png(paste0(path_file,"/estimate_significance_", save_name,".png"),
        width = 30, height = 20, units = "cm", res = 300)
    par(mar = c(15,10,2,2))
    # Hmsc::plotBeta(model_fit_mcmc, post = postBeta, param = "Mean", supportLevel = 0.5)
    Hmsc::plotBeta(model_fit_mcmc, post = postBeta, param = "Sign",
                   supportLevel = 0.95, colors = colorRampPalette(c("#2166AC", "white", "#A50026")),
                   mar = c(0,0,0,0) )
    dev.off()
    
    ## Keep support levels
    support_estimates <- postBeta[["support"]]
    rownames(support_estimates) <- model_fit_mcmc[["covNames"]]
    support_estimates <- as.data.frame(support_estimates) |> 
      tibble::rownames_to_column("covariate") |> 
      tidyr::pivot_longer(cols = -covariate,
                          names_to = "response", values_to = "support") |> 
      dplyr::mutate(support = dplyr::case_when(support > 0.95 ~ 1, 
                                               support < 0.05 ~ 1,
                                               TRUE ~ 0 ))|> 
      dplyr::mutate(support = factor(support, levels = c(0, 1)))
    
    
    ##### ridges plot  #####
    for(drivers in drivers_to_plot){
      
      #drivers = drivers_to_plot[[1]]
      drivers_name <- drivers
      all_drivers <- c(drivers, paste0(drivers, "_Deg1"),  paste0(drivers, "_Deg2"))
      
      #Filter estimates table
      df <- S_arranged |>
        dplyr::filter(covariate %in% all_drivers) |> 
        dplyr::left_join(support_estimates)
      
      medians <- df  |> 
        dplyr::filter(covariate %in% c(all_drivers[1], paste0(all_drivers[1], "_Deg1"))) |> 
        dplyr::group_by(response)  |> 
        dplyr::summarise(median_value = median(value)) |> 
        dplyr::arrange(median_value)
      
      df <- df |> 
        dplyr::left_join(dplyr::rename(grp_NN_NP, response = "contribution"))|> 
        dplyr::mutate(covariate = factor(covariate, levels = all_drivers)) |> 
        dplyr::mutate(response = factor(response, levels = medians$response)) 
      
      new_titles <- c(
        "protection_statusfull" = "Full protection",
        "n_fishing_vessels" = "Fishing pressure",
        "gravity" = "Gravity",
        "gdp" = "GDP",
        
        "protection_statusrestricted" = "Restricted MPA",
        "median_5year_analysed_sst" = "SST_5years",
        "median_5year_chl" = "Chlorophyll_5years",
        "q95_5year_degree_heating_week" = "DHW_5years"
      )
      
      ridges_plot <- ggplot(df) +
        aes(y = response, x = value,  fill = group) +
        ggridges::geom_density_ridges(aes(alpha = support), linewidth = 0.3)+#alpha=0.5, bandwidth = 0.005) +
        scale_fill_manual(values = c( "darkgoldenrod2", "forestgreen", "dodgerblue3"),
                          labels = c("NN" = "Nature-for-Nature",
                                     "NS" = "Nature-for-Society",
                                     "NC" = "Nature-as-Culture"))+
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
        scale_alpha_manual(values = c("0" = 0.25, "1" = 0.7), guide = "none") +
        hrbrthemes::theme_ipsum( axis_title_size = 0 ) +
        theme(
          legend.position="bottom",
          legend.text = element_text(size = 15, margin = margin(r = 20)),
          panel.spacing = unit(0.3, "lines"),
          strip.text.x = element_text(size = 16, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13)
          ) +
        labs(fill = "")+
        xlab(all_drivers) + ylab("Nature Contributions to People and Nature")+
        facet_wrap(~covariate, ncol = length(df$covariate),
                   scales = "free_x",
             labeller = labeller(covariate = new_titles)
        )
      
      # if(drivers == c("protection_statusfull", "n_fishing_vessels","gravity","gdp"))
      
      ggsave(filename = paste0(path_file,"/posterior_distribution_of_estimates_", save_name,
                               paste(drivers_name, collapse = "-"), ".jpg"),
             plot = ridges_plot, width = 12, height = 8)
      
    }
    
    
    ##### Estimates Heatmap  #####
    # S_df <- S |> 
    #   dplyr::mutate(Parameter = as.character(Parameter),
    #                 response = sapply(strsplit(sapply(strsplit(Parameter, ","),
    #                                                    function(x) x[2]), " "), function(x) x[2]),
    #                 cov = sapply(strsplit(sapply(strsplit(sapply(strsplit(Parameter, ","),
    #                                                              function(x) x[1]), "\\["), function(x) x[2]), " "),
    #                              function(x) x[1]))
    
    # support_estimates <- postBeta[["support"]]
    # rownames(support_estimates) <- model_fit_mcmc[["covNames"]]
    # support_estimates[support_estimates > 0.95] <- 1
    # support_estimates[support_estimates < 0.95] <- NA
    
    mean_estimate <- S_arranged |> 
      dplyr::group_by(covariate, response) |> 
      dplyr::summarise(mean_posterior_distrib  = mean(value)) |> 
      dplyr::filter(covariate != "Intercept")
    
    summary(mean_estimate$mean_posterior_distrib)
    
    ## Heat map of estimates
    palette <- rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788"))
    coef_matrix <- mean_estimate |> 
      tidyr::pivot_wider(names_from = "response", values_from = "mean_posterior_distrib") |> 
      tibble::column_to_rownames("covariate") |> 
      as.matrix()
    
    # for(i in rownames(coef_matrix)){
    #   for(j in colnames(coef_matrix)){
    #     coef_matrix[i,j] <- coef_matrix[i,j] * support_estimates[i,j]
    #   }
    # }
    
    jpeg(filename = paste0(path_file,"/heatmap_estimates_", save_name,".jpg"),
         units = "cm", width = 20, height = 20, res = 400 )
    stats::heatmap(coef_matrix,  Colv=T,
                   hclustfun=function(x) hclust(x, method="ward.D2"), 
                   scale='none', col=palette, cexCol=0.6, margins = c(7,10))
    dev.off()
    
    
    
    # ##### Estimate distribution #####
    # # Classification of covariates
    # covariates <- unique(mean_estimate$covariate)
    # envir <- c(grep("median", covariates, value = T))
    # habitat <- c(grep("500m", covariates, value = T),
    #              "depth", "algae", "coral", "Sand", "seagrass", "microalgal_mats",
    #              "other_sessile_invert", "Rock", "coralline_algae", "coral_rubble")
    # human <- setdiff(covariates, c(envir, habitat))
    # 
    # # Classification of contributions
    # cont_list <- unique(mean_estimate$response)
    # NP <- c("available_biomass", "selenium", "zinc", "omega_3" , "calcium",  "iron",                       
    #         "vitamin_A", "available_biomass_turnover", "NP_score")
    # NN <- setdiff(cont_list, NP)
    # 
    # #Resume data
    # coeff_plot <- mean_estimate |> 
    #   dplyr::mutate(cov_class = ifelse(covariate %in% envir, "environmental",
    #                                    ifelse(covariate %in% habitat, "habitat", "human"))) |> 
    #   dplyr::mutate(contrib_class = ifelse(response %in% NP, "Nature-for-People",
    #                                        "Nature-for-Nature"))
    # 
    # # Plot estimate distribution
    # 
    # plot_distri <- function(category = "Nature-for-Nature"){
    #   
    #   data <- dplyr::filter(coeff_plot, contrib_class ==  category)
    #   
    #   # order covariates
    #   df_ordered <- data |> 
    #     dplyr::group_by(covariate) |> 
    #     dplyr::summarise(mean_estimate = median(mean_posterior_distrib)) |> 
    #     dplyr::arrange(mean_estimate) |> 
    #     dplyr::pull(covariate)
    #   
    #   
    #   # Boxplot of estimates
    #   ggplot(data) +
    #     geom_vline(xintercept = 0, linetype = "dashed")+
    #     geom_boxplot(aes(y=factor(covariate, levels = df_ordered), 
    #                      x = mean_posterior_distrib, fill = cov_class),
    #                  alpha = 0.7) +
    #     hrbrthemes::theme_ipsum() +
    #     # scale_fill_manual(values = c("#9467bd", "#ff7f0e", "#2ca02c"))+
    #     harrypotter::scale_fill_hp_d(option = "Ravenclaw") +
    #     # facet_wrap(~ contrib_class, scales = "free_x") +
    #     xlab("Regression coefficient estimates") +
    #     ylab("Predictors") +
    #     labs(title = paste("Estimate distribution on", category, " contributions"))+
    #     theme( axis.title.x = element_text(size = 10),
    #            axis.title.y = element_text(size = 10),
    #            axis.text.y = element_text(size = 10),
    #            plot.title = element_text(size = 12),
    #            legend.position = "bottom" 
    #     )
    # }
    # 
    # plot_distri(category = "Nature-for-Nature") + plot_distri(category = "Nature-for-People")
    # ggsave(filename =paste0(path_file,"/mean_estimate_distribution_", save_name,".jpg"),
    #        width = 15, height = 8)
    
  
  } # END OF PLOT ESTIMATES
  
  #### Partial graph ####
  if(plot_partial_graph){
    
  Gradient = Hmsc::constructGradient(model_fit_mcmc,
                                     focalVariable = "median_5year_analysed_sst")
  predY = predict(model_fit_mcmc,
                  XData = Gradient$XDataNew,
                  studyDesign = Gradient$studyDesignNew,
                  ranLevels = Gradient$rLNew,
                  expected=TRUE)
                  
  # Hmsc::plotGradient(model_fit_mcmc, Gradient, pred = predY, measure = "Y",
  #                    showData = TRUE, index = 13)
  
  png(paste0(path_file,"/Partial_plot_SST_", save_name,".jpg"),
      width = 35, height = 25, units = "cm", res = 300)
  par(mfrow = c(6, 4), mar = c(4, 4, 2, 2)) 
  for (i in 1:22) {
    Hmsc::plotGradient(model_fit_mcmc, Gradient, pred = predY, measure = "Y", 
                       showPosteriorSupport =F, xlabel="",
                       showData = TRUE, index = i)
  } 
  dev.off()
  
  }
  
  ##----------------------------- Residuals -----------------------------------
  
  if(check_residuals){
    cat("Plot predictive residuals of the model... \n")
    
    #Run predictions on all data
    studyDesign <- model_fit_mcmc$studyDesign
    ranLevels <- model_fit_mcmc$ranLevels
    response_long <- Y_data |> 
      tibble::rownames_to_column("survey_id") |> 
      tidyr::pivot_longer(cols = -survey_id, names_to = "response",
                          values_to = "observation") 
    
    predY = predict(model_fit_mcmc, 
                    XData = X_data, 
                    studyDesign = studyDesign,
                    ranLevels = ranLevels, 
                    expected = T)
    
    PredY_mean <- as.data.frame(Reduce("+", predY) / length(predY))
    
    residuals <- Y_data - PredY_mean
    rownames(PredY_mean) <- rownames(Y_data)
      
    compare_pred <- PredY_mean |> 
      tibble::rownames_to_column("survey_id") |> 
      tidyr::pivot_longer(cols = -survey_id, names_to = "response",
                          values_to = "prediction") |> 
      dplyr::left_join(response_long) |> 
      dplyr::mutate(residuals = observation-prediction)
    
    
    # #histograms of residuals
    # distribution_plot(residuals, cols_plot = colnames(residuals))
    
    #Predictions vs observations
    ggplot(compare_pred)+
      geom_point(aes(x = observation, y = prediction, fill = response),
                 color = "grey40", alpha = 0.2, shape = 21) +
      hrbrthemes::theme_ipsum() +
      xlab("Observed contributions") + ylab("imputed")+
      geom_abline(slope = 1) + 
      ggpubr::stat_regline_equation(data = compare_pred,
                                    aes(x = observation, y = prediction,
                                        label = after_stat(rr.label)))   +
      facet_wrap(~response, scales = "free") +
      theme(legend.position="none", panel.spacing = unit(0.1, "lines"))
    ggsave(filename = paste0(path_file,"/Predictions_VS_Observations_Inference_", 
                             save_name,".jpg"),
           width = 15, height = 8)
    
    #residuals vs observations
    plot_interaction(compare_pred, var_facet_wrap = "response", 
                     X_values = "observation", Y_values = "residuals")+
      geom_hline(yintercept = 0)
    ggsave(filename = paste0(path_file,"/Residual_VS_Observations_", save_name,".jpg"),
           width = 15, height = 8)
    
    
    
    
    ## Check for spatial structure:
    if(check_spatial_autocorrelation){
      cat("... Measure Moran Indices... \n") 
      noise_magnitude <- 0.00001 #noise in coordinates
      
      geo_points <- metadata_model |> 
        dplyr::select(longitude, latitude)|>
        dplyr::rowwise() |> 
        dplyr::mutate(latitude = latitude + runif(1, -noise_magnitude, noise_magnitude),
                      longitude = longitude + runif(1, -noise_magnitude, noise_magnitude)) |>
        dplyr::ungroup()
      
      site_dist <- as.matrix(
        geodist::geodist(geo_points[,c("longitude", "latitude")], measure = "geodesic"))
      site_dist_inv <- 1/site_dist
      diag(site_dist_inv) <- 0
      
      moranI <- parallel::mclapply(1:ncol(residuals), mc.cores = 22, FUN = function(i){
        ape::Moran.I(residuals[,i], site_dist_inv)$observed
      })
      names(moranI) <- colnames(residuals)
      
      moranI_df <- as.data.frame(t(as.data.frame(moranI))) |> 
        dplyr::rename(moran_index = V1) |> 
        tibble::rownames_to_column("response")
      
      
      residual_coord <- compare_pred |> 
        dplyr::left_join(
          tibble::rownames_to_column(
            dplyr::select(metadata_model, longitude, latitude), "survey_id")) |> 
        dplyr::left_join(moranI_df)
      
      
      #Plot residuals
      plot_interaction(residual_coord, var_facet_wrap = "response", 
                       X_values = "latitude", Y_values = "residuals")+
        geom_hline(yintercept = 0)+
        geom_text( aes(x = Inf, y = Inf, 
                  label = paste0("Moran's I = ", round(moran_index, 2))), 
                  hjust = 1.1, vjust = 1.1, size = 3, color = "black")
      
      ggsave(filename = paste0(path_file,"/residuals_VS_latitude_MoranI_", save_name,".jpg"),
             width = 15, height = 8)
      
      # Plot all variogramms
      data <- metadata_model |> 
        dplyr::select(longitude, latitude) |> 
        dplyr::bind_cols(residuals)
      
      increment = 300
      resamp= 30
      contributions <- colnames(residuals)[order(colnames(residuals))]
      
      correlog_plot <- function(contrib = "calcium"){
        
        spatial_cor <- ncf::correlog(x= data$longitude, 
                                     y= data$latitude, 
                                     z = data[,contrib],
                                     increment = increment, 
                                     resamp= resamp, 
                                     latlon = TRUE) #distance and increment are in km
        ggplot() +
          geom_point(aes(y = spatial_cor$correlation[ which(spatial_cor$mean.of.class < 15000) ], 
                         x = spatial_cor$mean.of.class[ which(spatial_cor$mean.of.class < 15000) ]),
                     alpha = 0.6, size = 1) +
          geom_smooth(aes(y = spatial_cor$correlation[ which(spatial_cor$mean.of.class < 15000) ], 
                          x = spatial_cor$mean.of.class[ which(spatial_cor$mean.of.class < 15000) ]))+
          
          geom_point(aes(y = spatial_cor$correlation[ which(spatial_cor$p <= 0.05 &
                                                              spatial_cor$mean.of.class < 15000) ] ,
                         x = spatial_cor$mean.of.class[ which(spatial_cor$p <= 0.05 &
                                                                spatial_cor$mean.of.class < 15000) ]),
                     alpha = 0.8, size = 1, col ="red") +
          geom_vline(xintercept = spatial_cor[["x.intercept"]], linetype="dashed", 
                     color = "black", linewidth=1)+
          
          xlab("Distance class (in km)") + ylab("Moran Index")+
          labs(title =  contrib,
               subtitle = paste0("X intercept = ", round(spatial_cor[["x.intercept"]],1), "km"))+
          # ylim(-0.7,0.9)+
          theme_bw()
      }
      
      cat("... Plot Variogramm ... \n") 
      plots <- pbmcapply::pbmclapply(contributions, FUN = correlog_plot, 
                                     mc.cores = 10)
      
      all_plot <- Reduce(`+`, plots) +
        plot_layout(axis_titles = "collect")+
        plot_annotation(tag_levels = "a",
                        title = paste0("Increment = ", increment, "km, ", "permutations = ", resamp)) &
        theme(plot.tag = element_text(face = 'bold'),
              axis.title.x = element_text(size = 14),   
              axis.title.y = element_text(size = 14),)
      
      ggsave(filename = paste0(path_file,"/Variogramms_residuals_", save_name,".jpg"),
             all_plot, width = 22, height =14 )
      
    }#END OF CHECK SPATIAL AUTOCORRELATION
    
  } #END OF CHECK RESIDUALS
  
  ##----------------------------- Latent factors -----------------------------------
  
  if(latent_factors){
    cat("Plot latent factors of the model... \n")
    #site loadings (η) and species loadings (λ) on the latent factors
    postEta <- Hmsc::getPostEstimate(model_fit_mcmc, parName = "Eta")
    
    postEtamean1 <- postEta$mean[,1] # Latent factor 1
    postEtamean2 <- postEta$mean[,2] # Latent factor 2
    postEtamean3 <- postEta$mean[,3] # Latent factor 3
    
    id <- c()
    id_name <- NULL
    for(rd_level in 1:length(model_fit_mcmc[["rL"]])){
      id_temp <- model_fit_mcmc[["rL"]][[rd_level]][["pi"]]
      if(length(id_temp) > length(id)){ 
        id <- id_temp
        id_name <- names(model_fit_mcmc[["rL"]])[rd_level] 
        }
    }
    
    if(id_name == "site") id_name <- "site_code"
    if(id_name == "sample_unit") id_name <- "survey_id"
    
    
    data_lf <- data.frame(
      # survey_id = model_fit_mcmc[["ranLevels"]][["sample_unit"]][["pi"]],
      id_latent_factor = id,
      latent_factor_1 = postEtamean1,
      latent_factor_2 = postEtamean2,
      latent_factor_3 = postEtamean3)
    colnames(data_lf)[1] <- id_name
    
    data_lf <- data_lf|> 
      dplyr::left_join(
        dplyr::select(
          tibble::rownames_to_column(metadata_model, "survey_id"),
          all_of(id_name), longitude, latitude),
        keep = F) |> 
      unique()
    
    source(here::here("R/evaluation_prediction_model.R"))
    
    lf1 <- plot_Contrib_on_world_map(data = data_lf,
                              "latent_factor_1",
                              xlim=c(-180,180), ylim = c(-36, 31),
                              title="",jitter=1.5, pt_size=2, save=F)
    lf2 <- plot_Contrib_on_world_map(data = data_lf,
                                     "latent_factor_2",
                                     xlim=c(-180,180), ylim = c(-36, 31),
                                     title="",jitter=1.5, pt_size=2, save=F)
    lf3 <- plot_Contrib_on_world_map(data = data_lf,
                                     "latent_factor_3",
                                     xlim=c(-180,180), ylim = c(-36, 31),
                                     title="",jitter=1.5, pt_size=2, save=F)
    library(patchwork)
    latent_factors <- lf1 / lf2 / lf3
    latent_factors
    
    ggsave(latent_factors, 
           filename = paste0(path_file,"/Latent_factors_", save_name,".jpg"),
           width = 12, height = 10)
    
    # check_lf <- metadata_model |> 
    #   tibble::rownames_to_column("survey_id") |> 
    #   dplyr::left_join(data_lf) |> 
    #   dplyr::left_join(tibble::rownames_to_column(X_data, "survey_id"))
    # 
    # plot(check_lf$latent_factor_1 ~ check_lf$latitude)
    
    }#END OF CHECK LATENT FACTORS
  
                
  
} ## END OF FUNCTION PLOT_HMSC_RESULT
       


##--------------------------RUN HMSC PREDICTION------------------------------##


#' Run Hmsc prediction from the fitted model and new X data. Can run conditional 
#' prediction based on chosen responses.
#'
#' @param path file path to the model folder
#' @param folder_name name of the folder ontaining crossvalidations
#' @param conditional_prediction Boolean, if conditional prediction should be computed
#' @param mcmcStep_conditional Number of additionnal Mcmc step in case of conditional prediction
#' @param marginal_responses name of the responses predicted as marginal and used for conditional prediction
#'
#' @return a list of dataframe with observed, marginal predicted, and conditional predicted responses, for each crossvalidation
#' @export
#'
#' @examples
#' 
#' 
make_crossval_prediction_hmsc <- function(path = here::here("outputs/models/hmsc"),
                                 folder_name = "test",
                                 concatenate_chains = F,
                                 conditional_prediction = T,
                                 mcmcStep_conditional = 1, #"should be set high enough to obtain appropriate conditional predictions."
                                 marginal_responses = c("actino_richness",
                                                        "functional_distinctiveness")
                                 ){
  # PATHS
  save_init <- file.path(path, "cross_validation/init_multi")
  save_out <- file.path(path, "cross_validation/out_multi")
  localDir <- file.path(path, "cross_validation/multivariate")
  

  ## Extract data and CV names
  cross_val <- list.files(file.path(save_out,folder_name))
  load(file.path(save_out,folder_name, "dataset_crossvalidation.Rdata"))
  cv_files <- cross_val[grep("CV", cross_val)]
  
  if(concatenate_chains) cv_files <- unique(sub("_.*", "", cv_files))
  
  ## For each fold of crossvalidation, predict the Y responses on test dataset
  predictions_cv <- pbmcapply::pbmclapply(cv_files,
                                          mc.cores = 1, #parallelize the predict instead
                                          FUN = function(cv){
    # cv <-cv_files[4]
                                            
    if(concatenate_chains) cv <- paste0(cv, ".rds")
      
    # Initial data
    data <- train_model[[as.numeric(stringr::str_extract(cv, "(?<=CV)[0-9]+"))]]
    X_train <- data[[1]]
    X_test <- data[[3]]
    Y_test <- data[[4]]
    
    
    # Import initial model
    #Model design
    load(file = file.path(localDir, folder_name, cv))
    
    #Initial hmsc object
    init_obj_rds <- readRDS(file.path(save_init, folder_name, cv))
    init_obj <- jsonify::from_json(init_obj_rds)
    
    nSamples = init_obj[["samples"]]
    thin = init_obj[["thin"]]
    nChains = init_obj[["nChains"]]
    transient = init_obj[["transient"]]
    
    
    # Import posterior probability
    if(concatenate_chains){
      chainList = vector("list", nChains)
      importFromHPC <- vector("list", nChains+1)
      for(cInd in 1:nChains){
        chain_file_path = file.path(
          save_out, folder_name, paste0(gsub(".rds","",cv), "_chain_", cInd, ".rds")) ### To check ###
        chainList[[cInd]] = jsonify::from_json(readRDS(file = chain_file_path)[[1]])[[1]]
      }
    }else{
      all_chains <- readRDS(file.path(save_out, folder_name, cv))[[1]]
      importFromHPC <- jsonify::from_json(all_chains)
      chainList <- importFromHPC[1:nChains]
      #Result model
      cat(sprintf("fitting time %.1f h\n", importFromHPC[[nChains+1]] / 3600))
    }
    
    #Export and merge chains
    model_fit_mcmc <- Hmsc::importPosteriorFromHPC(model_fit,
                                                   chainList, 
                                                   nSamples, 
                                                   thin, 
                                                   transient)
    
    # Check data
    X_model <- model_fit_mcmc$XData
    if( sum(rownames(X_model) != rownames(X_train)) > 0 ){
      cat("Model and data doesn't match, check data.")
    } 
    
    ## Set study design of new data
    random_factors <- colnames(model_fit_mcmc[["studyDesign"]])
    studyDesign <- data.frame(sample_unit = as.factor(rownames(X_test)),
                              site = as.factor(X_test$site_code),
                              spatial = as.factor(rownames(X_test)),
                              year = as.factor(X_test$year),
                              country = as.factor(X_test$country),
                              ecoregion = as.factor(X_test$ecoregion))
    
    rL_asso = Hmsc::HmscRandomLevel(units = studyDesign$sample_unit)
    rL_site = Hmsc::HmscRandomLevel(units = unique(studyDesign$site))
    rL_year = Hmsc::HmscRandomLevel(units = studyDesign$year)
    rL_ecoregion = Hmsc::HmscRandomLevel(units = studyDesign$ecoregion)
    rL_country =  Hmsc::HmscRandomLevel(units = studyDesign$country)
    
    ranLevels = list(sample_unit = rL_asso,
                     site = rL_site,
                     year = rL_year,
                     ecoregion = rL_ecoregion,
                     country = rL_country)
  
    studyDesign <- studyDesign |> dplyr::select(all_of(random_factors))
    ranLevels <- ranLevels[random_factors]
    
    ## MAKE PREDICTIONS
    PredY_marginal <- PredY_conditional <- NULL
    
    #Marginal prediction
    predY_test <- predict(model_fit_mcmc, 
                          XData = X_test[, colnames(X_model)], 
                          studyDesign = studyDesign,
                          ranLevels = ranLevels, 
                          expected = F)
    
    PredY_marginal <- as.data.frame(Reduce("+", predY_test) / length(predY_test))
    
    if(conditional_prediction){
      
      #Predict some responses in marginal
      pred_marg <- predict(model_fit_mcmc, 
                           XData = X_test[, colnames(X_model)], 
                           studyDesign = studyDesign,
                           ranLevels = ranLevels, 
                           expected = F)
      
      PredY_marg <- PredY_marginal |> as.matrix()
      
      #use only selected reponses as backbone for other conditional predictions
      partition.sp <- as.numeric(colnames(PredY_marg) %in% marginal_responses)
      
      ############ TEST ##################
      # add_noise <- function(col, target_r2 = 0.2) {
      #   noise <- scale(rnorm(length(col)))
      #   new_col <- sqrt(target_r2) * col + sqrt(1 - target_r2) * noise
      #   return(new_col)
      # }
      # Y_test_noise <- apply(Y_test, 2, add_noise)
      ####################################
      
      conditional_prediction <- predict(model_fit_mcmc, 
                                        XData = X_test[, colnames(X_model)], 
                                        studyDesign = studyDesign,
                                        ranLevels = ranLevels, 
                                        expected = F,
                                        #use conditional prediction
                                        # Yc = PredY_marg,
                                        Yc = as.matrix(Y_test),
                                        # Yc = Y_test_noise,
                                        partition.sp = partition.sp,
                                        mcmcStep = mcmcStep_conditional,
                                        nParallel = parallel::detectCores()-4)
      
      PredY_conditional <- as.data.frame(
        Reduce("+", conditional_prediction) / length(conditional_prediction))
      
      # #Replace conditional prediction for the initial predicted responses
      # PredY_conditional[,marginal_responses] <- PredY_marg[,marginal_responses]
    }
    
    #save data
    list(observed = Y_test, 
         marginal_prediction = PredY_marginal, 
         conditional_prediction = PredY_conditional)
  }) #END OF PREDICTION ON EACH CROSSVAL
  
  
  save(predictions_cv, file = file.path(save_out,folder_name, "prediction_results.Rdata"))
  # load(file = file.path(save_out,folder_name, "prediction_results.Rdata"))
  
} # END OF FUNCION make_crossval_prediction_hmsc
                   



##--------------------------PLOT PREDICTIVE POWER------------------------------##


#' 
#'
#' @param path file path to the model folder
#' @param folder_name name of the folder ontaining crossvalidations
#' @param conditional_prediction Boolean, if conditional prediction should be computed
#' @param mcmcStep_conditional Number of additionnal Mcmc step in case of conditional prediction
#' @param marginal_responses name of the responses predicted as marginal and used for conditional prediction
#'
#' @return a list of dataframe with observed, marginal predicted, and conditional predicted responses, for each crossvalidation
#' @export
#'
#' @examples
#' 
#' 
plot_predictive_power <- function(path = here::here("outputs/models/hmsc"),
                                 folder_name = "test"){
  
  
  source("R/evaluation_prediction_model.R")
  
  path_file <- here::here("figures","models","hmsc", folder_name)    
  dir.exists(path_file)
  
  if(!dir.exists(path_file)) dir.create(path_file)
  
  
  # load explanatory power
  load(file.path(path_file, "explanatory_power_data.Rdata"))    
  
  # load predictions
  save_out <- file.path(path, "cross_validation/out_multi")
  load(file = file.path(save_out,folder_name, "prediction_results.Rdata"))
  
  ## Merge crossvalidations
  observed_contributions <- do.call(rbind, lapply(predictions_cv, `[[`, "observed"))
  marginal_prediction <- do.call(rbind, lapply(predictions_cv, `[[`, "marginal_prediction"))
  conditional_prediction <- do.call(rbind, lapply(predictions_cv, `[[`, "conditional_prediction"))
  
  ## Merge datasets
  
  observed_long <- observed_contributions |> 
    tibble::rownames_to_column("id") |> 
    tidyr::pivot_longer(-id, names_to = "responses", values_to = "observed")
  
  marginal_long <- marginal_prediction |> 
    tibble::rownames_to_column("id") |> 
    tidyr::pivot_longer(-id, names_to = "responses", values_to = "marginal_prediction")
  
  conditional_long <- conditional_prediction |> 
    tibble::rownames_to_column("id") |> 
    tidyr::pivot_longer(-id, names_to = "responses", values_to = "conditional_prediction")
  
  prediction <- observed_long |> 
    dplyr::left_join(marginal_long) |> 
    dplyr::left_join(conditional_long) |> 
    dplyr::mutate(residuals_mar = marginal_prediction - observed,
                  residuals_cond = conditional_prediction- observed)
  
  predictive_power_summary <- prediction |> 
    dplyr::group_by(responses) |> 
    dplyr::summarise(r_squared_marginal = summary(lm(marginal_prediction ~ observed))[["r.squared"]],
                     r_squared_conditional = summary(lm(conditional_prediction ~ observed))[["r.squared"]]) |> 
    dplyr::left_join(MF_table)
  
  
  ### Predictive power ###
  
  ## 1) Marginal predictions
  ggplot(prediction)+
    geom_point(aes(x = observed, y = marginal_prediction, fill = responses),
               color = "grey40", alpha = 0.2, shape = 21) +
    hrbrthemes::theme_ipsum() +
    xlab("Observed contributions") + ylab("Joint predictions")+
    geom_abline(slope = 1) + 
    ggpubr::stat_regline_equation(data = prediction,
                                  aes(x = observed, y = marginal_prediction,
                                      label = after_stat(rr.label)))   +
    facet_wrap(~responses, scales = "free") +
    theme(legend.position="none", panel.spacing = unit(0.1, "lines"))
  
  ggsave(filename = paste0(path_file, "/Joint_predictions_", folder_name, ".jpg"),
         width = 15, height = 8)
  
  
  
  
  # ## 2) Conditional predictions
  # ggplot(prediction)+
  #   geom_point(aes(x = observed, y = conditional_prediction, fill = responses),
  #              color = "grey40", alpha = 0.2, shape = 21) +
  #   hrbrthemes::theme_ipsum() +
  #   xlab("Observed contributions") + ylab("Conditional conditions")+
  #   geom_abline(slope = 1) + 
  #   ggpubr::stat_regline_equation(data = prediction,
  #                                 aes(x = observed, y = conditional_prediction,
  #                                     label = after_stat(rr.label)))   +
  #   facet_wrap(~responses, scales = "free") +
  #   theme(legend.position="none", panel.spacing = unit(0.1, "lines"))
  # 
  # ggsave(filename = paste0(path_file, "/Conditional_predictions_", folder_name, ".jpg"),
  #        width = 15, height = 8)
  
  
  ## 3) Summary
  png(paste0(path_file, "/Predictive_power_", folder_name, ".jpg"),
      width = 30, height = 15, units = "cm", res = 300)
  par(mfrow = c(1, 2))
  hist(predictive_power_summary$r_squared_marginal, xlim = c(0,1), 
       main=paste0("R_squared marginal Mean = ", 
                   round(mean(predictive_power_summary$r_squared_marginal),2)))
  hist(predictive_power_summary$r_squared_conditional, xlim = c(0,1), 
       main=paste0("R_squared conditional Mean = ", 
                   round(mean(predictive_power_summary$r_squared_conditional),2)))
  
  dev.off()
  
  
  ## 4) residuals
  plot_interaction(prediction, var_facet_wrap = "responses", 
                   X_values = "observed", Y_values = "residuals_mar")+
    geom_hline(yintercept = 0)
  
  plot_interaction(prediction, var_facet_wrap = "responses", 
                   X_values = "observed", Y_values = "residuals_cond")+
    geom_hline(yintercept = 0)
  # ggsave(filename = paste0(path_file, "/Residuals_conditional_predictions_", folder_name, ".jpg"),
  #        width = 15, height = 8)
  
  
  ## 5) R_squared summary
  
  predictive_power_summary <- predictive_power_summary |> 
    dplyr::arrange(R2) |> 
    dplyr::mutate(responses = factor(responses, levels = responses))
  
  ggplot(predictive_power_summary, aes(y = responses)) +
    geom_point(aes(x = r_squared_marginal, color = "Marginal prediction"), size = 3) + 
    geom_point(aes(x = r_squared_conditional, color = "Conditional prediction"), size = 3) + 
    geom_point(aes(x = R2, color = "Explanatory power"), size = 3) + 
    
    geom_vline(xintercept = mean(predictive_power_summary$r_squared_marginal),
               linetype = "dashed", color = "skyblue") + 
    geom_vline(xintercept = mean(predictive_power_summary$r_squared_conditional),
               linetype = "dashed", color = "orange") +
    geom_vline(xintercept = mean(predictive_power_summary$R2),
               linetype = "dashed", color = "grey50") +
    
    annotate("text", x = mean(predictive_power_summary$r_squared_marginal) + 0.02, y = -Inf, 
             label = round(mean(predictive_power_summary$r_squared_marginal), 2), 
             vjust = -1.5, color = "skyblue", size = 4) +
    annotate("text", x = mean(predictive_power_summary$r_squared_conditional) - 0.02, y = -Inf, 
             label = round(mean(predictive_power_summary$r_squared_conditional), 2), 
             vjust = -1.5, color = "orange", size = 4) +
    annotate("text", x = mean(predictive_power_summary$R2)+ 0.02, y = -Inf, 
             # label = paste0("Explanatory power: R² = ",round(mean(predictive_power_summary$R2), 2)), 
             label = round(mean(predictive_power_summary$R2), 2), 
             vjust =-1.5, color = "grey50", size = 4) +
    
    labs(x = "R_squared", y = "Responses", 
         title = "Explanatory and predictive power",
         color = "R_squared Type") +
    scale_color_manual(values = c("Marginal prediction" = "skyblue",
                                  "Conditional prediction" = "orange",
                                  "Explanatory power" = "grey50")) +
    theme_minimal()
  
  ggsave(filename = paste0(path_file, "/Joint_vs_Conditional_pred_power_", folder_name, ".jpg"),
         width = 15, height = 8)
  
  ##### TO DO: test spatial autocorrmapping = ##### TO DO: test spatial autocorrelations of residuals ######
  
  
} # END OF FUNCION plot_predictive_power


##--------------------------HMSC CONTERFACTUAL PREDICTION-----------------------------------##

#' Run hmsc prediction.
#'
#' @param path 
#' @param model_name 
#' @param concatenate_chains 
#' @param X_new_data
#' @param metadata 
#' @param is_counterfactual Precise if it is a couterfactual (change the past) or a
#'  scenario (change for the future). If TRUE, change = observed_value - new_value,
#'  i.e. current value - the one you should have in the "past". If FALSE, 
#'  change = new_value - observed_value, i.e. future value - current value.
#' 
#' 
#' @return save hmsc plots in figures/hmsc/file_name directory
#' @export
#'
#' @examples
#' 



run_hmsc_prediction <- function(path = path,
                                model_name = file_name,
                                concatenate_chains =F,
                                X_new_data = NULL,
                                metadata = NULL,
                                is_counterfactual = TRUE
){
  ## Paths
  save_init <- file.path(path, "init_multi")
  save_out <- file.path(path, "out_multi")
  localDir <- file.path(path, "multivariate")
  
  
  ## Import initial object
  
  #Model design
  load(file = file.path(localDir, paste0("model_fit_", model_name, ".Rdata")))
  
  #Initial hmsc object
  init_obj_rds <- readRDS(file.path(save_init, paste0("init_", model_name)))
  init_obj <- jsonify::from_json(init_obj_rds)
  
  nSamples = init_obj[["samples"]]
  thin = init_obj[["thin"]]
  nChains = init_obj[["nChains"]]
  transient = init_obj[["transient"]]
  
  
  ## Import posterior probability
  if(concatenate_chains){
    chainList = vector("list", nChains)
    importFromHPC <- vector("list", nChains+1)
    for(cInd in 1:nChains){
      chain_file_path = file.path(
        paste0(save_out,"output_", model_name, "/chain_", cInd, ".rds"))
      chainList[[cInd]] = jsonify::from_json(readRDS(file = chain_file_path)[[1]])[[1]]
    }
  }else{
    all_chains <- readRDS(file = file.path(save_out, paste0("output_", model_name)))[[1]]
    importFromHPC <- jsonify::from_json(all_chains)
    chainList <- importFromHPC[1:nChains]
    ##Result model
    cat(sprintf("fitting time %.1f h\n", importFromHPC[[nChains+1]] / 3600))
    
  }
  
  ##Export and merge chains
  model_fit_mcmc <- Hmsc::importPosteriorFromHPC(model_fit,
                                                 chainList, 
                                                 nSamples, 
                                                 thin, 
                                                 transient)
  
  ## Initial data
  Y_train <- as.data.frame(model_fit_mcmc$Y)
  X_train <- model_fit_mcmc$XData
  # metadata_model <- metadata[rownames(X_train),]
  
  
  ## Run prediction on original data
  studyDesign <- model_fit_mcmc$studyDesign
  ranLevels <- model_fit_mcmc$ranLevels
  
  
  predY <- predict(model_fit_mcmc, 
                   XData = X_train, 
                   studyDesign = studyDesign,
                   ranLevels = ranLevels)
  
  PredY_mean <- as.data.frame(Reduce("+", predY) / length(predY))
  
  
  ## Run prediction in conterfactual scenarios
  if(!is.null(X_new_data)){
    # X_new_data = X_train
    # metadata = metadata_sites
    
    ## Identify new conditions
    X_new_data <- X_new_data[, colnames(X_train)]
    rows_with_changes <- names(which(rowSums(X_new_data != X_train)>0))
    cat("Conditions have been changed in", length(rows_with_changes), "locations \n")
    
    ## Run new predictions
    predY_new_scenario <- predict(model_fit_mcmc, 
                                  XData = X_new_data, 
                                  studyDesign = studyDesign,
                                  ranLevels = ranLevels)
    
    PredY_new_scenario_mean <- as.data.frame(
      Reduce("+", predY_new_scenario) / length(predY_new_scenario))
    
    
    ## Observe changes
    effective_change <- PredY_mean[rows_with_changes,] |> 
      tibble::rownames_to_column("id") |> 
      tidyr::pivot_longer(cols = -id, names_to = "contribution", values_to = "original_prediction")|> 
      dplyr::left_join( 
        PredY_new_scenario_mean[rows_with_changes,] |> 
          tibble::rownames_to_column("id") |> 
          tidyr::pivot_longer(cols = -id, names_to = "contribution", values_to = "conterfactual")
        ) |> 
      dplyr::mutate(change = (-1)^is_counterfactual * (conterfactual - original_prediction))|> #if counterfactual: change = original_prediction-counterfactual, if it's ascenario, change = conterfactual - original_prediction.
      dplyr::left_join( data.frame(id = rownames(metadata),
                                   country = metadata$country)
      )
    
    }else{
      PredY_new_scenario_mean <- NULL
      effective_change <- NULL
    } # END of conterfactual scenarios
  
  list(predictions = PredY_mean, 
       new_scenario = PredY_new_scenario_mean,
       effective_change =effective_change)
} ## END OF FUNCTION RUN_HMSC_PREDICTION









#' Plot conterfactual scenarios
#'
#' @param 
#'
#' @return save hmsc plots in figures/hmsc/file_name directory
#' @export
#'
#' @examples
#' 



plot_conterfactual_scenarios <- function(path = path,
                                         model_name = model_name,
                                         concatenate_chains =F,
                                         X_new_data = X_new_mpa,
                                         metadata = metadata_sites,
                                         save_name = "effectiveness_high",
                                         selected_countries = selected_countries,
                                         plot_responders_on_map = FALSE,
                                         set_ids = NULL,
                                         is_counterfactual = TRUE,
                                         set_order_boxplot = NULL
){
  #Contributions group
  grp_NN_NP <- data.frame(
    contribution = c("actino_richness","functional_distinctiveness",
                      "iucn_species_richness" ,"mean_endemism",
                      "evolutionary_distinctiveness","functional_entropy",
                      "phylogenetic_entropy","herbivores_biomass",
                      "invertivores_biomass",  "piscivores_biomass",
                      "trophic_web_robustness", "mean_trophic_level",
      
                      "public_attention", "aesthetic",
                      "available_biomass", "selenium",
                      "zinc",   "omega_3", "calcium",
                      "iron","vitamin_A", "available_biomass_turnover"),
    group = c(rep("NN", 12), 
              rep("NC", 2),
              rep("NS", 8)))
  
  #Coastline
  coast <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')
  
  RdBu = c("#67001F", "#B2182B", 
           "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", 
           "#92C5DE", "#4393C3", "#2166AC", "#053061")
  
  # Create folder
  folder_name <- gsub(".rds", "", model_name)
  path_file <- here::here("figures","models","hmsc", "conterfactuals", folder_name)    

  if(!dir.exists(path_file)) dir.create(path_file)
  
  # Predict all contributions
  new_predictions <- run_hmsc_prediction(path,
                                         model_name,
                                         concatenate_chains,
                                         X_new_data,
                                         metadata,
                                         is_counterfactual)
  
  preds <- new_predictions[["predictions"]]
  conterfactual <- new_predictions[["new_scenario"]]
  effective_change <- new_predictions[["effective_change"]] |> 
    dplyr::left_join(grp_NN_NP)
  
  if(!is.null(set_ids)){
    effective_change <- effective_change |> dplyr::filter(id %in% set_ids)
  }
  
  #Back-transformation to original data to plot percent change
  load(file = here::here("outputs", "3_metadata_backtransformation_contrib.Rdata") )

  raw_change <- effective_change |> 
    dplyr::left_join(contributions_transformation) |> 
    dplyr::mutate(raw_original_prediction = ifelse(scaled, 
                                                   original_prediction*sd_site+mean_site,
                                                   original_prediction)) |> 
    dplyr::mutate(raw_original_prediction = ifelse(log_transformed, 
                                                   10^raw_original_prediction,
                                                   raw_original_prediction)) |> 
    dplyr::mutate(raw_original_prediction = ifelse(add_1, 
                                                   raw_original_prediction+1,
                                                   raw_original_prediction)) |> 
    dplyr::mutate(raw_conterfactual = ifelse(scaled, 
                                             conterfactual*sd_site+mean_site,
                                             conterfactual)) |> 
    dplyr::mutate(raw_conterfactual = ifelse(log_transformed, 
                                                   10^raw_conterfactual,
                                             raw_conterfactual)) |> 
    dplyr::mutate(raw_conterfactual = ifelse(add_1, 
                                             raw_conterfactual+1,
                                             raw_conterfactual)) |> 
    dplyr::mutate(raw_change_values = 
                    (-1)^is_counterfactual * (raw_conterfactual - raw_original_prediction),
                  raw_change_percent = 
                  ((-1)^is_counterfactual * (raw_conterfactual - raw_original_prediction))
                  /raw_original_prediction *100) 
  
  # distribution_plot(raw_change, longer = F, index_values = c("contribution","raw_change_percent"))
  
  
      
  # #---- Plot conterfactual VS original predictions ---
  # plot_interaction(effective_change, var_facet_wrap = "contribution",
  #                  X_values = "original_prediction", Y_values = "conterfactual",
  #                  xlabel = "Original predictions", 
  #                  ylabel = paste("Conterfactual predictions in", save_name))
  # 
  # ggsave( width = 15, height = 8, filename = file.path(
  #   path_file,paste0("Conterfactual_VS_predictions_", save_name, folder_name, ".jpg"))
  #       )
  
  ###---- Plot histogramms of changes ---####
  distribution_plot(raw_change, longer = F,
                    index_values =  c("contribution", "raw_change_values"))

  ggsave( width = 15, height = 8, filename = file.path(
    path_file,paste0("Raw_changes_histogramms_", save_name, folder_name, ".jpg"))
        )
  
  ###---- Plot distributions of contribution changes ---####
  distrib_boxplot <- function(data, x, y, fill, hline = 0, title = NULL,
                              prop_outliers = 0,
                              include_stat = TRUE){
    p <- ggplot(data) +
          geom_hline(yintercept = 0, color = "grey", size = 1) +
          aes_string(x= x, y= y, fill = fill)+
          # scale_fill_manual(values = c("forestgreen", "dodgerblue3"),
          #                   labels = c("NN" = "Nature-for-Nature contributions",
          #                              "NP" = "Nature-for-People contributions"))+
          scale_fill_manual(values = c( "darkgoldenrod2", "forestgreen", "dodgerblue3"),
                            labels = c("NN" = "Nature-for-Nature",
                                       "NS" = "Nature-for-Society",
                                       "NC" = "Nature-as-Culture"))+
          # geom_violin(trim = FALSE, position = position_dodge(width =1), alpha = 0.7) +
          geom_boxplot(alpha = 0.7, outliers = T) +
          
          xlab("") + ylab("Contributions change in counterfactual scenarios") +
          labs(title=title, fill = "")+
          theme_minimal() +
          theme(legend.position = "bottom",
                panel.grid.minor = element_blank(),
                legend.text = element_text(size = 15, margin = margin(r = 20)),
                panel.spacing = unit(0.3, "lines"),
                axis.text.y = element_text(size = 13),
                axis.title = element_text(size = 13),
                axis.text.x = element_text(angle = 0, hjust = 0.5,size = 13))+
          
          #Flip coordinates and cut outliers
          coord_flip(ylim=c(quantile(data[,y][[1]], prop_outliers), 
                            quantile(data[,y][[1]], 1-prop_outliers)))
    
    if(include_stat){ 
      p <- p + 
        stat_summary(fun = "median", geom = "text",
                     aes(label = round(after_stat(y), 1)),
                     position = position_dodge(width = 0.75), vjust = 0.5,
                     hjust = -0.2, size = 4, color = "grey50") +
        
        geom_hline(yintercept = hline,
                   linetype = "dashed", color = "coral3") +
        annotate("text", y = hline + 0.02, x = -Inf, 
                 label = round(hline, 2), vjust = 0.05, hjust = -0.1,
                 vjust = -1.5, color = "coral3", size = 4)
    }
  
    p
  }
  
  effective_change <- effective_change  |> 
    dplyr::mutate(contribution = reorder(contribution, change, FUN = median))
  
  mean <- dplyr::group_by(effective_change, contribution) |> 
    dplyr::summarise(change = median(change)) |> 
    dplyr::summarise(mean = mean(change)) |> 
    as.numeric()
  
  distrib_boxplot(effective_change, x = "contribution", y = "change",
                  fill = "group", title = save_name, hline = mean,
                  prop_outliers = 0)
  
  ggsave(width = 8, height = 8, filename = file.path(
    path_file,paste0("Changes_distrib_", save_name, "_", folder_name, ".jpg")))
  
  
  
  ## Change in percent on raw values
  raw_change <- raw_change  |> 
    dplyr::mutate(contribution = reorder(contribution, raw_change_percent,
                                         FUN = median, decreasing = T ))
  # set_order_boxplot <- levels(raw_change$contribution)
  if(!is.null(set_order_boxplot)){
    raw_change$contribution <- factor(raw_change$contribution,
                                         levels = set_order_boxplot)
  }
  
  mean <- dplyr::group_by(raw_change, contribution) |> 
    dplyr::summarise(change = median(raw_change_percent)) |> 
    dplyr::summarise(mean = mean(change)) |> 
    as.numeric()
  
  distrib_boxplot(raw_change, x = "contribution", y = "raw_change_percent",
                  fill = "group", title = paste(save_name, "(change in percent)"),
                  hline = mean, prop_outliers = 0.015) +
    ylab("Percent change in contributions in counterfactual scenarios (%)")
  
  ggsave(width = 12, height = 10, filename = file.path(
    path_file,paste0("Percent_changes_distrib_", save_name, "_", folder_name, ".jpg")))
    
  
  
  ## Log scale ##
  raw_change_log_transformed <- raw_change |> 
    dplyr::mutate(change_percent_log = dplyr::case_when(
      raw_change_percent > 1 ~ log10(raw_change_percent),
      raw_change_percent < -1 ~ -log10(-(raw_change_percent)),
      T ~ raw_change_percent
    )) |> 
    dplyr::mutate(contribution = reorder(contribution, change_percent_log,
                                         FUN = median, decreasing=T))
  
  if(!is.null(set_order_boxplot)){
    raw_change_log_transformed$contribution <- 
      factor(raw_change_log_transformed$contribution, levels = set_order_boxplot)
  }
  
  median_labels <- raw_change |> 
    dplyr::group_by(contribution, group) |> 
    dplyr::summarize(median_raw = median(raw_change_percent)) |> 
    dplyr::mutate(median_position = dplyr::case_when(
      median_raw > 0 ~ log10(median_raw+1),
      median_raw < 0 ~ -log10(-(median_raw-1)),
      T ~ 0))
  
  hline = sign(mean(median_labels$median_raw)) * 
    abs(log10(abs(mean(median_labels$median_raw))))

  distrib_boxplot(raw_change_log_transformed, x = "contribution", y = "change_percent_log",
                 fill = "group", title = paste(save_name, "(change in percent, LOG-scale)"),
                  hline = 0, prop_outliers = 0.01,
                 include_stat = F) + 
    geom_text(data = median_labels,
              aes(x = contribution, y = median_position, 
                  label = paste(signif(median_raw, 2), "")),
              vjust = 0.5, hjust = -0.3, size = 4, color = "grey30") +
    
    geom_hline(yintercept = hline,
               linetype = "dashed", color = "coral3") +
    annotate("text", y = hline , x = -Inf, 
             label = paste(round(mean(median_labels$median_raw), 2), "%"),
             vjust = 0, hjust = -0.2, color = "coral3", size = 4)+
    scale_y_continuous(
      breaks = c(-3, -2, -1, 0, 1, 2, 3),      
      labels = c("-1000%", "-100%", "-10%", "0", "+10%", "+100%", "+1000%")
    )+
    ylab("Percent change in contributions in counterfactuals (%, log-scale)")
  
  ggsave(width = 12, height = 10, filename = file.path(
    path_file,paste0("Percent_changes_LOG_scale_", save_name, "_", folder_name, ".jpg")))
  
  ###---- Plot changes in each countries ---####
  
  # #obs heterogeneity by country
  # data_1_contrib <- effective_change |>
  #   dplyr::filter(contribution == "available_biomass") |>
  #   dplyr::mutate(country = reorder(country, change, FUN = median))
  #
  # distrib_boxplot(data_1_contrib, x = "country", y = "change", fill = "country",
  #                 hline = mean(data_1_contrib$change))


  
  # Original PCA
  pca <- FactoMineR::PCA(preds, scale.unit = T, graph=F, ncp=15)
  
  # factoextra::fviz_screeplot(pca, ncp=15)
  
  pca_plot_no_labels <- factoextra::fviz_pca_biplot(pca,
                                          axes = c(1,2),title="",
                                          label = "none",
                                          geom=c("point"), 
                                          select.var = list(cos2 = 0.3),
                                          pointshape=21,
                                          stroke=0, pointsize=2,
                                          alpha.ind = 0.7,
                                          alpha.var = 0.7,
                                          fill.ind = "grey",  
                                          col.var = "grey20",
                                          repel = TRUE)
  
  #label position
  data <- data.frame(obsnames=row.names(pca$ind$coord), pca$ind$coord)
  datapc <- data.frame(varnames=rownames(pca$var$coord), pca$var$coord)
  mult <- min(
    (max(data[,"Dim.2"]) - min(data[,"Dim.2"])/(max(datapc[,"Dim.2"])-min(datapc[,"Dim.2"]))),
    (max(data[,"Dim.1"]) - min(data[,"Dim.1"])/(max(datapc[,"Dim.1"])-min(datapc[,"Dim.1"])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (datapc$Dim.1),
                      v2 = .7 * mult * (datapc$Dim.2) )
  
  #Add labels
  pca_plot <- pca_plot_no_labels +
    ggrepel::geom_text_repel(
      data= datapc[which(factoextra::facto_summarize(
                             pca, element= "var",axes = c(1,2))$cos2 > 0.3)
                   ,], #extract cos2 of variables
      aes(x= v1*1.05, y= v2*1.03, label= gsub("_"," ", varnames)),
      size=5,
      color =  "grey20", 
      force_pull = 5,
      direction = "both")
    
  pca_plot
  
  # Calculate barycenters (centroids) for each country in the 6 first dimensions
  if(!is.null(set_ids)){
    new_conditions <- set_ids
  }else{ new_conditions <- unique(effective_change$id) }
  
  coord_old_points <- as.data.frame(pca$ind$coord[, 1:6]) |> 
    tibble::rownames_to_column("id") |> 
    dplyr::filter(id %in% new_conditions) |> #Considering only points that will change
    tibble::column_to_rownames("id")
  
  countries <- metadata[rownames(coord_old_points),]$country
  
  barycenters <- aggregate(coord_old_points, 
                           by = list(countries), FUN = mean) |> 
    dplyr::filter(Group.1 %in% selected_countries)
  colnames(barycenters) <- c("country", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
  
  if(is.null(set_ids)){
    absent_countries <- selected_countries[!selected_countries %in% countries]
    missing_rows <- data.frame(country = absent_countries,
                               PC1 = NA, PC2 = NA, PC3 = NA, PC4 = NA, PC5 = NA, PC6 = NA)
    barycenters <- rbind(barycenters, missing_rows)
  }
  
  # Add barycenters to the plot
  pca_plot_barycenter <- pca_plot + 
    geom_point(data = barycenters, aes(x = PC1, y = PC2, fill = country), 
               shape = 21, size = 5, color = "black", stroke = 1)
  
  
  ### Observe new scenario ###
  # (1) Project the new predictions in the same space
  projected_points <- predict(pca, newdata = conterfactual)
  
  # (2) calculate new barycenters
  coord_new_points <- as.data.frame(projected_points$coord[, 1:6]) |> 
    tibble::rownames_to_column("id") |> 
    dplyr::filter(id %in% new_conditions) |> 
    tibble::column_to_rownames("id")
  
  barycenters_new <- aggregate(coord_new_points, 
                               by = list(countries), FUN = mean) |> 
    dplyr::filter(Group.1 %in% selected_countries)
  colnames(barycenters_new) <- c("country", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
  if(is.null(set_ids)){
    barycenters_new <- rbind(barycenters_new, missing_rows)
  }
  barycenters_movement <- merge(barycenters, barycenters_new, 
                                by = "country", suffixes = c("_old", "_new"))
  
  
  plot_mvt <- pca_plot_barycenter + 
    geom_segment(data = barycenters_movement, 
                 aes(xend = PC1_old, yend = PC2_old, 
                     x = PC1_new, y = PC2_new,
                     color = country), 
                 # arrow = arrow(length = unit(0.3, "cm")),
                 linewidth = 1, linetype = "longdash")+ 
    geom_point(data = barycenters, 
               aes(x = PC1, y = PC2, fill = country, shape = "Current conditions"), 
               size = 5, color = "black", stroke = 1) +
    geom_point(data = barycenters_new, 
               aes(x = PC1, y = PC2, fill = country, shape = "Counterfactual conditions"), 
               size = 4, color = "black", alpha = 1) +
    coord_cartesian(xlim= c(-7.7,7))+
    labs(title = save_name, fill = "Country", shape = "Conditions")+#, y = "", x="") +
    scale_shape_manual(values = c("Current conditions" = 21, "Counterfactual conditions" = 24)) +
    guides(
      fill = guide_legend(nrow = 3, override.aes = list(shape = 21)),  
      shape = guide_legend(nrow = 2, override.aes = list(fill = NA)),
      color = "none"
    )+
    theme(legend.position = "bottom",
          axis.text = element_text(size = 0),
          axis.ticks = element_line(linewidth = 0),
          axis.title.x = element_text(vjust = 127, hjust = 0, size = 12),
          axis.title.y = element_text(vjust = -197, hjust = 0, angle = 90, size = 12),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9, face = "bold", 
                                      margin = margin(r=.5, unit = "cm")),
          legend.spacing.x = unit(1, "cm"),
          legend.key.spacing.x = unit(.1, "cm") )
  plot_mvt
  
  # if(is.null(set_ids)){ # if the ID are fixed to plot the barycenters, don't save theses
    #barycenters alone because they don't necessarily represent the reefs with changes conditions.
  ggsave(plot = plot_mvt, width=13, height= 10, filename = file.path(
    path_file,paste0("Barycenters_movement_", save_name,"_", folder_name, ".jpg"))
  )
  # }
  
  
  # ## General barycenter:
  # barycenter_tot <- apply(pca$ind$coord[, 1:2], 2, mean) # 0, 0 by construction
  # barycenters_new <- t(as.data.frame(apply(projected_points$coord[, 1:2], 2, mean)))
  # pca_plot + geom_segment(data = barycenters_new, 
  #                         aes(x = 0, y = 0, 
  #                             xend = Dim.1, yend = Dim.2), 
  #                         arrow = arrow(length = unit(0.3, "cm")),
  #                         linewidth = 1)+
  #   labs(title = save_name)
  no_labels <- pca_plot_no_labels +
    geom_segment(data = barycenters_movement, 
                 aes(xend = PC1_old, yend = PC2_old, 
                     x = PC1_new, y = PC2_new,
                     color = country), 
                 # arrow = arrow(length = unit(0.3, "cm")),
                 linewidth = 1, linetype = "longdash")+ 
    geom_point(data = barycenters, 
               aes(x = PC1, y = PC2, fill = country, shape = "Current conditions"), 
               size = 5, color = "black", stroke = 1) +
    geom_point(data = barycenters_new, 
               aes(x = PC1, y = PC2, fill = country, shape = "Counterfactual conditions"), 
               size = 4, color = "black", alpha = 1) +
    coord_cartesian(xlim= c(-7.5,7))+
    labs(title = save_name, fill = "Country", shape = "Conditions", y = "", x="") +
    scale_shape_manual(values = c("Current conditions" = 21, "Counterfactual conditions" = 24)) +
    guides(
      fill = guide_legend(nrow = 2, override.aes = list(shape = 21)),  
      shape = guide_legend(nrow = 2, override.aes = list(fill = NA)),
      color = "none"
    )+
    theme(legend.position = "bottom",
          axis.text = element_text(size = 0),
          axis.ticks = element_line(linewidth = 0))
  
  
  ###---- Responders ---####
  
  ## Plot the distances in the contribution multidimensional space
  if(plot_responders_on_map){
      
    distances <- data.frame(
      id = rownames(coord_old_points),
      distance = apply(coord_old_points - coord_new_points, 1, function(row) sqrt(sum(row^2)))
    ) |> 
      dplyr::left_join(tibble::rownames_to_column(metadata, "id") ) |> 
      dplyr::arrange(distance)
    
    
    plot_mpa <-function(distances, xlim=c(-180,180), ylim = c(-36, 31),
                        legend_pos = "right", jitter = 2){
      ggplot(distances) +
        geom_sf(data = coast, color = "grey30", fill = "lightgrey",
                aes(size=0.1)) +
        
        geom_point( #size = 2,
                    na.rm = T,
                   position = position_jitter(width =jitter, height =jitter),
                   alpha = 0.8,
                   colour = "black",
                   stroke=0.1,
                   shape = 21,
                   aes(x = longitude, y = latitude,
                       fill=distance, size = distance)) +
        
        coord_sf(xlim, ylim , expand = FALSE) +
        guides(alpha = "none", size = "none", colour = "none") +
        # scale_shape_manual(values=c(21,24,23))+
        scale_fill_gradientn(colors = RdBu)+
        scale_size(range = c(1, 3), guide = "legend") +
        
        theme_minimal()+
        labs(title = "",
             x="", y= "") +
        theme(legend.position = legend_pos,
              plot.title = element_text(size=15, face="bold"),
              legend.text = element_text(size=13),
              legend.title = element_text(size=15),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
        )
    }
    
    map <-  plot_mpa(distances , xlim=c(-180,180), ylim = c(-36, 31),
                     legend_pos = "none", jitter = 1)+
      geom_rect(aes(xmin = 110, xmax = 160, ymin = -32, ymax = -7), color = "black", fill= "transparent")+
      geom_rect(aes(xmin = -95, xmax = -67, ymin = -3, ymax = 18), color = "black", fill= "transparent")
    
    gold_coast_map <- plot_mpa( distances,ylim = c(-32, -7),
                               xlim= c(110,160), legend_pos = "right", jitter = 0.5)
    
    caraib_map <- plot_mpa( distances, ylim = c(-3, 18),
                            xlim= c(-95,-67), legend_pos = "none", jitter = 0.5)
    
    ggpubr::ggarrange(map, # First row with world map
                      ggpubr::ggarrange(caraib_map, gold_coast_map,  
                                        ncol = 2, labels = c("B", "C"), widths = c(1, 1.3)), # Second row with zooms
                      nrow = 2, labels = "A") 
  
    ggsave(plot = last_plot(), width=15, height= 8, filename = file.path(
        path_file,paste0("High_responders_", save_name, folder_name, ".jpg"))
      )
  } # End of if(plot_responders_on_map)
  
  
  #keep the plots in memory
  return(list(plot_mvt, no_labels, raw_change_log_transformed))
} #END OF plot_conterfactual_scenarios



