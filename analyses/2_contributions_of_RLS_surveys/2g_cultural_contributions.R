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


##------------------- 1) Aesthetic scores -------------------####
# # Install missing librairies for python
# reticulate::py_install("pytz")
# reticulate::py_install("sympy")
# reticulate::py_install("Pillow")
# reticulate::py_install(c("torch", "torchvision"))
# 
# 
# #run aesthetic for new species
# reticulate::source_python("R/04_inference.py")
# #IF IMPOSSIBLE TO RUN RETICULATE: open a terminal in the R/ folder, and run the 
# # python code via bash: > python3 04_inference.py 

#Observe new inference:
aesth_new <- read.csv("outputs/2g_aesthetic_inference_new_sp.csv") |> 
  dplyr::mutate(old = ifelse(grepl("Langlois22", image_name), 1, 0)) |> 
  dplyr::mutate(file_name = gsub(".png", "", image_name)) |> 
  dplyr::mutate(file_name = gsub("_Langlois22", "", file_name)) |> 
  dplyr::mutate(sp_name = gsub("^(.*?)_[A-Z]_[0-9]+$", "\\1", file_name)) 
  

#compare inference
aesth <- aesth_new |> 
  dplyr::left_join(aesth_2022)
plot(aesth$predicted_score~aesth$esthe_score)
abline(a=0, b=1)
cor.test(aesth$predicted_score,aesth$esthe_score)
df <- aesth[aesth$old == 1,]
hist(rank(df$predicted_score)-rank(df$esthe_score))
