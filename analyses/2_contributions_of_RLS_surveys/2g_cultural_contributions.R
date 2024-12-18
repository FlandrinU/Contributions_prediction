################################################################################
##
## 
##
## 2g_cultural_contributions.R
##
## 26/02/2024
##
## Ulysse Flandrin
##
###############################################################################"

## cleaning memory
rm(list=ls())

##-----------------Loading packages-------------------
# pkgs <- c("reticulate")

##------------------- loading datasets-------------------
#Species traits
load(file= here::here("outputs", "RLS_species_traits_inferred.Rdata"))


#Occurrence matrix of communities
load("data/derived_data/2_occurrence_matrix_sp_survey.Rdata")

#RLS observations
load(file = here::here("data/derived_data/2_rls_actino_trop.Rdata"))

## Load functions
source("R/check_scientific_names.R")

# ##------------------- 1) Aesthetic and public attention at the species level -------------------####
# 
# ### AESTHETIC ###
# 
# # # Install missing librairies for python
# # reticulate::py_install("pytz")
# # reticulate::py_install("sympy")
# # reticulate::py_install("Pillow")
# # reticulate::py_install(c("torch", "torchvision"))
# 
# 
# #run aesthetic for new species
# # reticulate::source_python("R/inference_aesthetic_score_Langlois2022.py")
# #IF IMPOSSIBLE TO RUN RETICULATE: open a terminal in the R/ folder, and run the
# # python code via bash: > python3 04_inference.py
# 
# ## Observe new inference:
# aesth_new <- read.csv("outputs/2g_aesthetic_inference_new_sp.csv") 
# 
# # Check inference with old scores of Langlois 2022
# check <- aesth_new |> 
#   dplyr::mutate(old = ifelse(grepl("Langlois22", image_name), 1, 0)) |> 
#   dplyr::mutate(file_name = gsub(".png", "", image_name)) |> 
#   dplyr::mutate(file_name = gsub("_Langlois22", "", file_name)) |> 
#   dplyr::mutate(sp_name = gsub("^(.*?)_[A-Z]_[0-9]+$", "\\1", file_name)) |> 
#   dplyr::left_join(aesth_2022)
# 
# plot(check$predicted_score~check$esthe_score)
# abline(a=0, b=1)
# cor.test(check$predicted_score,check$esthe_score)
# # INFERED SCORES CONSISTENT FOR THE 55 PICURES IN COMMON -> WE CAN INFER NEW 
# # SCORES FROM PICTURES AND ADD THEM TO THE FILE OF LANGLOIS ET AL. 2022
# 
# asthetic_inferred <- check |> 
#   dplyr::filter(old == 0) |> 
#   dplyr::select(file_name = image_name, sp_name, predicted_score) |> 
#   dplyr::group_by(sp_name) |> 
#   dplyr::mutate(esthe_score = max(predicted_score)) |> 
#   dplyr::filter(esthe_score == predicted_score) |> 
#   dplyr::select(-predicted_score)
# 
# ## Merge all scores
# aesthetic_score <- aesth_2022 |> 
#   dplyr::select(file_name, sp_name, esthe_score) |> 
#   dplyr::bind_rows(asthetic_inferred)
# 
# 
# ### PUBLIC ATTENTION ###
# 
# human_interest <- hum_int |>
#   dplyr::select(public_attention = public,
#                 academic_knowledge = acad,
#                 sp_name = fb_sci_name)
# 
# ### CHECK SCIENTIFIC NAMES ###
# 
# cultural_sp_contrib <- dplyr::full_join(aesthetic_score, human_interest) |> 
#   dplyr::mutate(sp_name = gsub("_", " ", sp_name)) 
# 
# cultural_sp_contrib <- code_sp_check(cultural_sp_contrib, original_name = 'sp_name',
#                                      mc_cores = 15)|> 
#   dplyr::select(-worms_id, -sp_name, -file_name, -check)
# 
# cultural_all_species <- inferred_species_traits |> 
#   dplyr::select(-public_interest, -academic_knowledge ) |> 
#   tibble::rownames_to_column("rls_species_name") |> 
#   dplyr::left_join(cultural_sp_contrib) |> 
#   dplyr::group_by(fishbase_name) |> 
#   tidyr::fill(tidyr::everything(), .direction = 'updown') |> 
#   dplyr::mutate(esthe_score = max(esthe_score)) |> 
#   unique()
# 
# save(cultural_all_species, file = here::here("outputs","2g_cultural_all_species.Rdata"))
#   

##------------------- 2) Aesthetic scores (survey level) -------------------####

#' The following code computes community aesthe values at the survey level. It 
#' comes from Mclean et al. 2025 "Conserving the beauty of fish communities"
#' @author Matthew McLean, \email {mcleamj@@gmail.com}
#' 
# INIT ----

# coefficient a and b from Tribot, A.S, Deter, J., Claverie, T., Guillhaumon, F., Villeger, S., & Mouquet, N. (2019). Species diversity and composition drive the aesthetic value of coral reef fish assemblages. Biology letters, 15, 20190703, doi:10.1098/rsbl.2019.0703
# will be used to compute the aesthetic scores of assemblages
intercept_sr <- 7.0772149
slope_sr     <- 0.20439752

#data
aesthe_species <- inferred_species_traits |> 
  tibble::rownames_to_column("rls_species_name")
sp_pres_matrix <- surveys_sp_occ

# Compute the aesthe contribution of each species
# with parameters from Tribot, A.S, Deter, J., Claverie, T., Guillhaumon, F., Villeger, S., & Mouquet, N. (2019). Species diversity and composition drive the aesthetic value of coral reef fish assemblages. Biology letters, 15, 20190703, doi:10.1098/rsbl.2019.0703
# positive and negative effect are relative to the espected effect computed with the species 

## Computing the aesthe_effect
aesthe_species$aesthe_effect <- (log(aesthe_species$aesthetic) - 7.3468679)/7.937672

# ----

# COMPUTE AESTHETICS ----
#aesthe_survey : predicted aesthe 
#aesthe_SR_survey : predicted aesthe based only on species richness

surveyid_vect <- row.names(sp_pres_matrix)  

survey_aesth <- do.call(rbind, pbmcapply::pbmclapply((1:length(surveyid_vect)), function(i){
  
  #i=1
  # surveyID <- "1003646"
  surveyID <-surveyid_vect[i]
  # presence_absence of the species of the survey
  
  vector_abs_pres <- sp_pres_matrix[surveyID,]
  
  # species present
  sp_survey <- colnames(vector_abs_pres)[vector_abs_pres[1,]>0]
  
  # number of species of the survey
  nb_species <- length(sp_survey)
  
  # species component in community aesthetic
  species_effect <- sum(
    aesthe_species$aesthe_effect[aesthe_species$rls_species_name %in%  gsub("_"," ",sp_survey)],
    na.rm = T
    )
  
  #Prevalence of NA in the survey
  missing_species <- aesthe_species[aesthe_species$rls_species_name %in%  gsub("_"," ",sp_survey),] |> 
    dplyr::filter(is.na(aesthetic))
  
  if(nrow(missing_species)>0){
    missing_rls_obs <- rls_actino_trop |> 
      dplyr::filter(survey_id == surveyID,
                    rls_species_name %in% missing_species$rls_species_name)
    
    prop_biom_aesth <- 1 - sum(missing_rls_obs$raw_biomass) / unique(missing_rls_obs$biomass_tot_survey)
    prop_abund_aesth <- 1 - sum(missing_rls_obs$total ) / unique(missing_rls_obs$abundance_tot_survey)
  }else{
    prop_biom_aesth <- 1 
    prop_abund_aesth <- 1
  }
  
  # aesthe of the survey
  E <- intercept_sr + slope_sr * log(nb_species) + species_effect
  score <-  exp(E)
  
  E <-  intercept_sr + slope_sr * log(nb_species)
  score_SR  <-  exp(E) 
  
  # compute the number of species with positive and negative effect in each survey
  
  vect  <-  (aesthe_species$aesthe_effect * vector_abs_pres)
  nb_sp_pos_survey   <-  length(which(vect>0))
  nb_sp_neg_survey   <-  length(which(vect<0))
  
  cbind.data.frame(survey_id=surveyid_vect[i],
                   nb_species=nb_species,
                   aesthe_survey=score,
                   aesthe_SR_survey=score_SR,
                   nb_sp_pos_survey=nb_sp_pos_survey,
                   nb_sp_neg_survey=nb_sp_neg_survey,
                   prop_biom_aesth = prop_biom_aesth, 
                   prop_abund_aesth = prop_abund_aesth)
  
}, mc.cores = parallel::detectCores()-1))

plot(survey_aesth$nb_species,survey_aesth$aesthe_survey)

write.csv(survey_aesth, here::here("outputs", "survey_aesth.csv"), row.names = FALSE)

# ---- 

#Check outputs
old_surveys_aesth <- read.csv(here::here("data/derived_data/survey_aesth_McLean2024.csv")) |> 
  dplyr::mutate(SurveyID = as.character(SurveyID)) |> 
  dplyr::select(survey_id = SurveyID,
                aesthe_old = aesthe_survey, 
                aesthe_SR_old = aesthe_SR_survey) |> 
  dplyr::full_join(survey_aesth)

plot(old_surveys_aesth$aesthe_survey ~ old_surveys_aesth$aesthe_old) ; abline(a=0, b=1)
plot(old_surveys_aesth$aesthe_SR_survey ~ old_surveys_aesth$aesthe_SR_old) ; abline(a=0, b=1)


## REMOVE NON REPRESENTATIVE ESTIMATIONS: (less than 80% abundance or biomass with known species)
survey_aesth_filtered <- survey_aesth |> 
  dplyr::filter(prop_biom_aesth > 0.8 & prop_abund_aesth > 0.8 )


##------------------- 2) Public attention (survey level) -------------------####
cultural <- inferred_species_traits
cultural_survey <-  pbmcapply::pbmclapply(rownames(surveys_sp_occ), function(id){
  cat(id, "\n")
  sp <- names(surveys_sp_occ[id, which(surveys_sp_occ[id,] > 0)])
  if( length(sp) > 0){
    n <- length(sp)
    
    mean_academic_knowledge <- mean(cultural[sp, "academic_knowledge"], na.rm = T)
    mean_public_attention <- mean(cultural[sp, "public_attention"], na.rm = T)
    mean_esthe <- mean(cultural[sp, "aesthetic"], na.rm = T)
    
    academic_knowledge <- quantile(cultural[sp, "academic_knowledge"], 0.75, na.rm = T)
    public_attention <- quantile(cultural[sp, "public_attention"], 0.75, na.rm = T)
    quantile_esthe <- quantile(cultural[sp, "aesthetic"], 0.75, na.rm = T)
    
    cult <- as.numeric(c(id, academic_knowledge, public_attention, quantile_esthe,
                         mean_academic_knowledge, mean_public_attention, mean_esthe,
                         n))
  }else{cult <- c(id, NA, NA, NA, NA, NA, NA, NA)}
  
  names(cult) <- c("survey_id", "academic_knowledge", "public_attention", "quantile_esthe",
                   "mean_academic_knowledge", "mean_public_attention", "mean_esthe",
                   "diversity")
  cult
}, mc.cores = parallel::detectCores()-1)

public_contrib_survey <- data.frame(do.call(rbind, cultural_survey)) |>
  dplyr::mutate( across(everything(), ~as.numeric(.))) |>
  dplyr::mutate(survey_id = as.character(survey_id))


## Check the importance of NA 
rows_with_na <- which(!complete.cases(inferred_species_traits[, "public_attention"]))
sp_na <- rownames(inferred_species_traits[rows_with_na,])

NA_prop <- rls_actino_trop |>
  dplyr::group_by(survey_id) |>
  dplyr::filter(rls_species_name %in% sp_na) |> 
  dplyr::group_by(survey_id, abundance_tot_survey, biomass_tot_survey) |> 
  dplyr::summarise(total = sum(total),
                   biomass = sum(biomass)) |> 
  dplyr::mutate(prop_abund = total / abundance_tot_survey,
                prop_biom = biomass / biomass_tot_survey)|> #proportion of species with NA
  dplyr::filter(prop_abund > 0.2 | prop_biom  > 0.2)

#SELECT WHEN "NAs" SPECIES REPRESENT MORE THAN 20% OF THE SURVEY, IN ABUNDANCE OR BIOMASS:
na_survey <- unique(NA_prop$survey_id)

## REMOVE NON REPRESENTATIVE ESTIMATIONS:
if(length(na_survey)>0){
  public_contrib_survey[
  public_contrib_survey$survey_id == na_survey, "public_attention"] <- NA
}



##------------------- 3) Save cultural contributions -------------------####

cultural_contributions <- public_contrib_survey |>
  dplyr::full_join(survey_aesth_filtered) |> 
  dplyr::select(survey_id, public_attention, aesthe_survey)

fb_plot_species_traits_completeness(dplyr::rename(cultural_contributions,
                                                  species = survey_id))

save(cultural_contributions, file = here::here("outputs", "2g_cultural_contributions.Rdata"))
