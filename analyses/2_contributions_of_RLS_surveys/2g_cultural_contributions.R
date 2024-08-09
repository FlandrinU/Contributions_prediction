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

#cultural values of species
hum_int <- read.csv2("data/raw_data/Mouquet2023_Human_Interest_final_table.csv", dec = ",")

#aesthetic values of species from Langlois et al. 2022
aesth_2022 <- read.csv(
  "data/raw_data/aesthetic_deep_model_Langlois_2022/Langlois_el_al_2022.csv")

#Occurrence matrix of communities
load("data/derived_data/2_occurrence_matrix_sp_survey.Rdata")

#RLS observations
load(file = here::here("data/derived_data/2_rls_actino_trop.Rdata"))

## Load functions
source("R/check_scientific_names.R")

##------------------- 1) Aesthetic and public attention at the species level -------------------####

### AESTHETIC ###

# # Install missing librairies for python
# reticulate::py_install("pytz")
# reticulate::py_install("sympy")
# reticulate::py_install("Pillow")
# reticulate::py_install(c("torch", "torchvision"))


#run aesthetic for new species
reticulate::source_python("R/inference_aesthetic_score_Langlois2022.py")
#IF IMPOSSIBLE TO RUN RETICULATE: open a terminal in the R/ folder, and run the
# python code via bash: > python3 04_inference.py

## Observe new inference:
aesth_new <- read.csv("outputs/2g_aesthetic_inference_new_sp.csv") 

# Check inference with old scores of Langlois 2022
check <- aesth_new |> 
  dplyr::mutate(old = ifelse(grepl("Langlois22", image_name), 1, 0)) |> 
  dplyr::mutate(file_name = gsub(".png", "", image_name)) |> 
  dplyr::mutate(file_name = gsub("_Langlois22", "", file_name)) |> 
  dplyr::mutate(sp_name = gsub("^(.*?)_[A-Z]_[0-9]+$", "\\1", file_name)) |> 
  dplyr::left_join(aesth_2022)

plot(check$predicted_score~check$esthe_score)
abline(a=0, b=1)
cor.test(check$predicted_score,check$esthe_score)
# INFERED SCORES CONSISTENT FOR THE 55 PICURES IN COMMON -> WE CAN INFER NEW 
# SCORES FROM PICTURES AND ADD THEM TO THE FILE OF LANGLOIS ET AL. 2022

asthetic_inferred <- check |> 
  dplyr::filter(old == 0) |> 
  dplyr::select(file_name = image_name, sp_name, predicted_score) |> 
  dplyr::group_by(sp_name) |> 
  dplyr::mutate(esthe_score = max(predicted_score)) |> 
  dplyr::filter(esthe_score == predicted_score) |> 
  dplyr::select(-predicted_score)

## Merge all scores
aesthetic_score <- aesth_2022 |> 
  dplyr::select(file_name, sp_name, esthe_score) |> 
  dplyr::bind_rows(asthetic_inferred)


### PUBLIC ATTENTION ###

human_interest <- hum_int |>
  dplyr::select(public_interest = public,
                academic_knowledge = acad,
                sp_name = fb_sci_name)

### CHECK SCIENTIFIC NAMES ###

cultural_sp_contrib <- dplyr::full_join(aesthetic_score, human_interest) |> 
  dplyr::mutate(sp_name = gsub("_", " ", sp_name))

cultural_sp_contrib <- code_sp_check(cultural_sp_contrib, original_name = 'sp_name', mc_cores = 15)

cultural_all_species <- inferred_species_traits |> 
  dplyr::select(-public_interest, -academic_knowledge ) |> 
  tibble::rownames_to_column("rls_species_name") |> 
  dplyr::left_join(cultural_sp_contrib) |> 
  dplyr::select(-sp_name, -file_name) |> 
  dplyr::group_by(fishbase_name) |> 
  tidyr::fill(tidyr::everything(), .direction = 'updown') |> 
  dplyr::mutate(esthe_score = max(esthe_score)) |> 
  unique()
  

##------------------- 2) Aesthetic scores (survey level) -------------------####
survey_aesth_all <- read.csv(here::here("data", "derived_data", "survey_aesth.csv"))


survey_aesth <- survey_aesth_all |> 
  dplyr::select(survey_id = SurveyID, aesthe_survey) |> 
  dplyr::mutate(survey_id = as.character(survey_id))


## TO DO: Check NAs in aesthetic #######################"

##------------------- 3) Public attention (survey level) -------------------####
cultural <- tibble::column_to_rownames(cultural_all_species, var ="rls_species_name")
cultural_survey <- lapply( rownames(surveys_sp_occ), function(id){
  cat(id, "\n")
  sp <- names(surveys_sp_occ[id, which(surveys_sp_occ[id,] >0)])
  if( length(sp) > 0){
    mean_academic_knowledge <- mean(cultural[sp, "academic_knowledge"], na.rm = T)
    mean_public_interest <- mean(cultural[sp, "public_interest"], na.rm = T)
    
    academic_knowledge <- quantile(cultural[sp, "academic_knowledge"], 0.75, na.rm = T)
    public_interest <- quantile(cultural[sp, "public_interest"], 0.75, na.rm = T)
    
    cult <- as.numeric(c(id, academic_knowledge, public_interest,
                         mean_academic_knowledge, mean_public_interest))
  }else{cult <- c(id, NA, NA, NA, NA)}
  
  names(cult) <- c("SurveyID", "academic_knowledge", "public_interest", 
                   "mean_academic_knowledge", "mean_public_interest")
  cult
})

public_contrib_survey <- data.frame(do.call(rbind, cultural_survey)) |>
  dplyr::mutate( across(everything(), ~as.numeric(.))) |>
  dplyr::mutate(survey_id = as.character(SurveyID))


## Check the importance of NA 
rows_with_na <- which(!complete.cases(cultural_all_species[, "public_interest"]))
sp_na <- cultural_all_species$rls_species_name[rows_with_na]

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

# ## REMOVE NON REPRESENTATIVE ESTIMATIONS:
# public_contrib_survey[
#   public_contrib_survey$survey_id == na_survey, "public_interst"] <- NA

####### TO DO WHEN ALL PUBLIC INTEREST WILL BE CALCULATED 


##------------------- 4) Save cultural contributions -------------------####

cultural_contributions <- public_contrib_survey |>
  dplyr::full_join(survey_aesth) |> 
  dplyr::select(survey_id, public_interest, aesthe_survey)

save(cultural_contributions, file = here::here("outputs", "2g_cultural_contributions.Rdata"))
