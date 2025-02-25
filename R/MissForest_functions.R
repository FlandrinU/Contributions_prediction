################################################################################
##
##  
##
## MissForest_functions.R
##
## 29/11/2022
##
## Ulysse Flandrin
##
################################################################################
 # ##-----------------Loading packages-------------------
 # pkgs <- c("here", "dplyr", "missForest", "pbmcapply", "patchwork", "ggplot2",
 #           "maditr", "slam")
 # nip <- pkgs[!(pkgs %in% installed.packages())]
 # nip <- lapply(nip, install.packages, dependencies = TRUE)
 # ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))


##-------------preping dataset-------------
preping_data <- function(species_traits_df = tropical_species_traits){
  
  traits_data <- species_traits_df |>  
    dplyr::select(-fishbase_name,-spec_code, -worms_id) |> 
    tibble::column_to_rownames("rls_species_name") |> 
    dplyr::mutate(across(where(is.character), as.factor)) 
  
  ## Check data
  factor_length = data.frame()
  for (i in 1:length(traits_data[,])){
    #If the trait is a factor the get the length
    if (is.factor(traits_data[,i])==TRUE){
      data_frame = data.frame(id = colnames(traits_data[i]),
                              length = length(unique(traits_data[,i])))
      factor_length = rbind(factor_length,data_frame)
    }
  } #OK
  
  #IF any of the factors has over 53 categories filter it out of the data
  if(any(factor_length>53)){
    
    over53 <- factor_length |> 
      dplyr::filter(length>=53)
    
    data_to_infer <- traits_data |> 
      dplyr::select(-(all_of(over53$id)))
    # warning("Some of your traits had more than 53 categories. These traits were
    #         filtered out during the missForest and then re-added to your data later") 
    
    #Else keep it as it is
  }else{
    data_to_infer <- traits_data 
  }
  
  
  ## long format table ##
  
  traits_data_factors <- traits_data  |>
    #extract only categorial variables 
    dplyr::select(all_of(factor_length$id)) |>
    tibble::rownames_to_column("rls_species_name")  |> 
    #long format table
    tidyr::pivot_longer(cols = IUCN_category:last_col(),
                        names_to = "variable",
                        values_to = "observed")
  
  traits_data_num <- traits_data  |>
    #extract only categorial variables 
    dplyr::select(-all_of(unique(traits_data_factors$variable))) |>
    tibble::rownames_to_column("rls_species_name")  |> 
    #long format table
    tidyr::pivot_longer(cols = c(Length:last_col()),
                        names_to = "variable",
                        values_to = "observed")
  
  return(list(data_to_infer, traits_data_factors, traits_data_num, factor_length))
         
} # END of function preping_data


##-------------MCA: Evaluate inference with Missforest - Phylogeny with MCA-------------

fct_missforest_evaluation_MCA <- function(data_to_infer,
                                          traits_data_factors = traits_data_factors,
                                          traits_data_num = traits_data_num,
                                          factor_length = factor_length,
                                          model_iteration = 5,
                                          maxiter_mf=3, #missforest parameter #10
                                          ntree_mf=10,  #missforest parameter #100
                                          prop_NA = 0.2){ # Proportion of NA created to evaluate the model
  
  
  model_eval_missforest <- pbmcapply::pbmclapply(c(1:model_iteration), 
                                                 mc.cores = parallel::detectCores() - 3, 
                                                 FUN = function(N){
   
                                                   
   categorials <- dplyr::filter(factor_length,length<53)$id
   categorials <- categorials[! categorials %in% c("phylum","class","order")]
   Dims <- colnames(data_to_infer)[grep("Dim",colnames(data_to_infer))]
                                   
   phylo <- data_to_infer |> dplyr::select(all_of(Dims))
   traits <- data_to_infer |> dplyr::select(-phylum, -class, -order, -all_of(Dims))
   
   
   #Produce 20% of NA in the matrix
   data_withNA <- dplyr::bind_cols(
     phylo,
     missForest::prodNA(traits, noNA = prop_NA)
   )
   
   #long format Data frame with NA
   data_withNA_factors <- data_withNA  |>
     dplyr::select(all_of(categorials)) |>
     tibble::rownames_to_column("rls_species_name")  |> 
     tidyr::pivot_longer(cols = c(IUCN_category:last_col()),
                         names_to = "variable",
                         values_to = "obs_NA")
   
   data_withNA_num <- data_withNA  |>
     dplyr::select(-all_of(c(unique(data_withNA_factors$variable), Dims))) |>
     tibble::rownames_to_column("rls_species_name")  |> 
     tidyr::pivot_longer(cols = c(Length:last_col()),
                         names_to = "variable",
                         values_to = "obs_NA") 
   
   
   ## Imputing data ##
   impute = missForest::missForest(data_withNA, 
                                   maxiter = maxiter_mf, ntree = ntree_mf,
                                   variablewise = T, verbose = T)
   
   
   #long format Data frame with imputed values 
   imputed_factors <- impute$ximp  |>
     dplyr::select(all_of(categorials)) |>
     tibble::rownames_to_column("rls_species_name")  |> 
     tidyr::pivot_longer(cols = c(IUCN_category:last_col()),
                         names_to = "variable",
                         values_to = "imputed")
   
   imputed_num <- impute$ximp  |>
     dplyr::select(-all_of(c(unique(data_withNA_factors$variable), Dims))) |>
     tibble::rownames_to_column("rls_species_name")  |> 
     tidyr::pivot_longer(cols = c(Length:last_col()),
                         names_to = "variable",
                         values_to = "imputed") 
   
   
   
   ### Model evaluation ###
   ##Factors imputation
   eval_factors <- traits_data_factors |> 
     dplyr::left_join(data_withNA_factors) |> 
     dplyr::left_join(imputed_factors) |> 
     dplyr::mutate(missForest = N) |> #identify the iteration N of the mclapply
     dplyr::filter(is.na(obs_NA) & !is.na(observed)) # keep only rows where we can assess the imputation accuracy
   
   
   ##Numeric imputation
   eval_num <- traits_data_num |> 
     dplyr::left_join(data_withNA_num) |> 
     dplyr::left_join(imputed_num) |> 
     dplyr::mutate(missForest = N) |> #identify the iteration N of the mclapply
     dplyr::filter(is.na(obs_NA) & !is.na(observed)) # keep only rows where we can assess the imputation precision
   
   
   result <- list(eval_factors, eval_num)
   return(result)
   
 }) ## END OF MCLAPPLY ON MISSFOREST
  
  model_eval_missforest
  
} ## END OF fct_missforest_evaluation_MCA


## ----------------------- Missforest application -----------------

#' missforest_applied
#'
#' @param data_to_infer the species X traits dataframe with some missing values
#' @param factor_length a two column df with the names of factoral variables of 'data_to_infer'
#'                      and the number of categories inside each.
#' @param traits_data_factors same as 'data_to_infer' but in long format with only factor variables
#' @param traits_data_num same as 'data_to_infer' but in long format with only numerical variables
#' @param var_to_infer a list of the variable we want to infer
#' @param confidence_threshold a threshold of consistency between the different missforest: 
#' proportion of missforest with the same result for factors, and 1-sd/median for numericals.
#' @param model_iteration number of missforest to do, to test the consistency of the model
#' @param maxiter_mf max number of iteration in each missforest 
#' @param ntree_mf number of tree in each missforest 
#'
#' @return the species X traits dataframe with inferred values by missforest for
#'  the selected variables, and values for which most models (>confidence_threshold)
#'  converged towards the same result
#'
#' @examples


missforest_applied <- function(data_to_infer, 
                               factor_length,
                               traits_data_factors,
                               traits_data_num,
                               var_to_infer= c("Length", "K", "IUCN_category", "trophic_guild"),
                               confidence_threshold = 0.8,
                               model_iteration = 2,
                               maxiter_mf=1, #missforest parameter #10
                               ntree_mf=10){ #missforest parameter #10
                                      
  
  res_missforest <- pbmcapply::pbmclapply(c(1:model_iteration), 
                                          mc.cores = parallel::detectCores() - 5, 
                                          FUN = function(N){
                                            
    Dims <- colnames(data_to_infer)[grep("Dim",colnames(data_to_infer))]
                                            
    ## Imputing data ##
    impute = missForest::missForest(data_to_infer, 
                                    maxiter = maxiter_mf, ntree = ntree_mf,
                                    # maxiter = 10, ntree = 100, ######################## maxiter = 10, ntree = 100 when it's OK
                                    variablewise = T, verbose = T)
    
    
    #long format Data frame with imputed values 
    imputed_factors <- impute$ximp  |>
      dplyr::select(all_of(dplyr::filter(factor_length,length<53)$id)) |>
      tibble::rownames_to_column("rls_species_name")  |> 
      tidyr::pivot_longer(cols = c(IUCN_category:last_col()),
                          names_to = "variable",
                          values_to = "imputed")
    
    imputed_num <- impute$ximp  |>
      dplyr::select(-all_of(c(unique(imputed_factors$variable),Dims))) |>
      tibble::rownames_to_column("rls_species_name")  |> 
      tidyr::pivot_longer(cols = c(Length:last_col()),
                          names_to = "variable",
                          values_to = "imputed") 

    
    ### Result ###
    ##Factors imputation
    eval_factors <- traits_data_factors |> 
      dplyr::left_join(imputed_factors) |> 
      dplyr::filter(is.na(observed) &         
                   variable %in% var_to_infer) |>  # keep only rows where we want to infer the data 
      dplyr::select(-observed)
    
    
    ##Numeric imputation
    eval_num <- traits_data_num |> 
      dplyr::left_join(imputed_num) |> 
      dplyr::filter(is.na(observed) &         
                      variable %in% var_to_infer) |>  # keep only rows where we want to infer the data 
      dplyr::select(-observed)
      
    
    result <- list(eval_factors, eval_num)
    return(result)
    
  }) ## END OF MCLAPPLY ON MISSFOREST
  save(res_missforest, file = here::here("outputs", "Missforest_application_raw_resulst.Rdata"))
  # load( file = here::here("outputs", "Missforest_application_raw_resulst.Rdata"))
  
  ### Extract infered data ###
  flat_list <- unlist(res_missforest, recursive = F)
  
  #gather all factors estimations
  estimates_factors <- flat_list[[1]]
  for(i in seq(3,length(flat_list),2)){
    estimates_factors <- estimates_factors |> 
      dplyr::left_join(flat_list[[i]],
                       by=c("rls_species_name","phylum","class","order","family","variable"),
                       suffix = c("", paste(".",i)))
  }
  
  #gather all numerical estimations
  estimates_num <- flat_list[[2]]
  for(i in seq(4,length(flat_list),2)){
    estimates_num <- estimates_num |> 
      dplyr::left_join(flat_list[[i]],
                       by=c("rls_species_name","phylum","class","order","family","variable"),
                       suffix = c("", paste(".",i)))
  }
  
  
  ##Assess Missforest confidence
  # factors
  imputed <- estimates_factors[, grep("imputed", colnames(estimates_factors))]
  norm <- apply(imputed, 1, FUN = function(x){names(which.max(table(x)))})
  agreement <- apply(imputed, 1, FUN = function(x){ 
    mean(x == names(which.max(table(x))))
  })
  
  final_estimation_factors <- estimates_factors |> 
    dplyr::mutate(norm = norm, confidence = agreement) |> 
    dplyr::select(-grep("imputed", colnames(estimates_factors)))
  
  
  # numerical and all imputation
  imputed <- estimates_num[, grep("imputed", colnames(estimates_num))]
  median <- apply(imputed, 1, median)
  deviation <- apply(imputed, 1, FUN = function(x){ 1- sd(x) / mean(x) }) #Coefficient of variation = sd/mean => we want 1-CV > confidence_threshold
  
  final_estimation_num <- estimates_num |> 
    dplyr::mutate(norm = median, confidence = deviation) |> 
    dplyr::select(-grep("imputed", colnames(estimates_factors)))    
  
  
  ### Insert imputed data into initial data frame ###
  Dims <- colnames(data_to_infer)[grep("Dim",colnames(data_to_infer))]
  
  final_imputation <- rbind(final_estimation_factors, final_estimation_num)
  infered_data <- data_to_infer[, !colnames(data_to_infer) %in% Dims]
  
  for( i in 1:nrow(final_imputation)){
    if(final_imputation$confidence[i] > confidence_threshold){ # if most of missforest converged into the same result
      sp <- final_imputation$rls_species_name[i]
      var <- final_imputation$variable[i]
      infered_data[sp, var] <- final_imputation$norm[i]
    }
  }
  
  infered_data <- infered_data |>  
    dplyr::mutate(across(where(is.character), as.numeric)) 
  
  return(infered_data)
} ## END OF missforest_applied
