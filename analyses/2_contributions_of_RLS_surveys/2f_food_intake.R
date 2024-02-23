################################################################################
##
## 
##
## 2f_food_intake.R
##
## 23/01/2024
##
## Ulysse Flandrin
##
###############################################################################"

## cleaning memory
rm(list=ls())

##------------------- loading datasets-------------------
#Species traits
load(file= here::here("outputs", "RLS_species_traits_inferred.Rdata"))
species_traits <- tibble::rownames_to_column(inferred_species_traits, "species_name")
  
#RLS observations
load(file = here::here("data/derived_data/2_rls_actino_trop.Rdata"))

#Survey metadata
load(file = here::here("data", "raw_data", "environmental_covariates", "RLS_surveys_covariates.Rdata"))


##------------------- loading functions-------------------
source(here::here("R", "calc_prod_rfishprod.R"))
source(here::here("R", "calc_prod_transect.R"))



## -------------------1) Available biomass -------------------####

# This script determines the total biomass of fishery species in each transect,
#  using the list of families which are fished in the world 
#  (expert opinion, Cinner et al. 2020)

max_size <- species_traits |> 
  dplyr::select(species_name, max_length=Length)

data_surveys_fishery <- rls_actino_trop |>
  dplyr::left_join(max_size)


#families in which all species are targeted in the family: (Cinner et al. 2020 + other expert opinion )
all_sp <- c("Acanthuridae", "Caesionidae", "Carangidae", "Ephippidae", "Haemulidae", "Kyphosidae",
            "Labridae", "Lethrinidae", "Lutjanidae", "Mullidae", "Nemipteridae", "Scaridae",
            "Scombridae", "Serranidae", "Siganidae", "Sparidae", "Sphyraenidae",
            "Mugilidae", "Bothidae", "Scorpaenidae", "Sciaenidae")

#families with species targeted if larger than 20cm: (Cinner et al. 2020 + other expert opinion )
target_larger_20cm <- c("Balistidae", "Holocentridae", "Pomacanthidae", "Priacanthidae")

#non-targeted families: (Cinner et al. 2020 + other expert opinion )
untargeted_fam <- c("Chaetodontidae","Cirrhitidae", "Diodontidae", "Grammatidae", 
                    "Monacanthidae", "Pempheridae", "Pinguipedidae", "Pseudochromidae", 
                    "Synodontidae", "Tetraodontidae", "Zanclidae","Fistulariidae",
                    "Ostraciidae", "Pomacentridae")

unclassified_fam <- setdiff(unique(data_surveys_fishery$family), 
                            c(all_sp, target_larger_20cm, untargeted_fam))
unclassified_fam #OK


#keep only targeted fishes
data_surveys_fishery <- data_surveys_fishery |> 
  dplyr::filter(family %in% all_sp | (family %in% target_larger_20cm & max_length > 20))




### biomass of fishery species in each survey ###
surveys_fishery_sp_biom <- data_surveys_fishery |> 
  dplyr::group_by(survey_id, species_name) |>
  dplyr:: summarize( available_biomass = sum(biomass) ) |>
  tidyr::pivot_wider(names_from = species_name, values_from = available_biomass, 
                     values_fill = 0) |>
  tibble::column_to_rownames(var="survey_id") 



surveys_fishery_biom <- data.frame(
  survey_id = rownames(surveys_fishery_sp_biom),
  available_biomass = rowSums(surveys_fishery_sp_biom))




## -------------------2) Available nutrients -------------------####

#aggregate by species in each transect
# (keep only species targeted by fisheries)
data_surveys_nutrient <- rls_actino_trop |>
  dplyr::group_by(survey_id, species_name) |>
  dplyr::summarise(biomass_sp = sum(biomass)) |> 
  dplyr::filter(species_name %in% unique(data_surveys_fishery$species_name)) 


#nutrients per species per transect
RLS_nut_sp_surv <- data_surveys_nutrient |>
  dplyr::left_join(
    dplyr::select(species_traits,
                  species_name, Calcium, Iron, Omega3, Selenium, VitaminA, Zinc)) |>
  dplyr::mutate(Selenium_tot  = Selenium  * biomass_sp,
                Zinc_tot      = Zinc      * biomass_sp,
                Omega3_tot   = Omega3     * biomass_sp,
                Calcium_tot   = Calcium   * biomass_sp,
                Iron_tot      = Iron      * biomass_sp,
                VitaminA_tot = VitaminA   * biomass_sp)


#sum of nutrients in each transect, and concentration of it  
RLS_nut_surv <- RLS_nut_sp_surv |>
  dplyr::group_by(survey_id) |>
  dplyr::summarise(biom_tot    = sum(biomass_sp),
                   Selenium_C  = sum(Selenium_tot) / biom_tot,
                   Zinc_C      = sum(Zinc_tot)     / biom_tot,
                   Omega3_C   = sum(Omega3_tot)  / biom_tot,
                   Calcium_C   = sum(Calcium_tot)  / biom_tot,
                   Iron_C      = sum(Iron_tot)     / biom_tot,
                   VitaminA_C = sum(VitaminA_tot)/ biom_tot ) |> 
  dplyr::select(-biom_tot)




## Check the importance of NA 
cols <- c("Calcium", "Iron", "Omega3", "Selenium", "VitaminA", "Zinc")
rows_with_na <- which(!complete.cases(species_traits[, cols]))
cols_with_na <- apply(is.na(species_traits[rows_with_na, cols]), 1, function(x) cols[x])
sp_na <- species_traits$species_name[rows_with_na]

NA_prop <- rls_actino_trop |>
  dplyr::group_by(survey_id) |>
  dplyr::mutate(abundance_tot_survey = sum(total),
                   biomass_tot_survey = sum(biomass)) |> 
  dplyr::filter(species_name %in% sp_na) |> 
  dplyr::group_by(survey_id, abundance_tot_survey, biomass_tot_survey) |> 
  dplyr::summarise(total = sum(total),
                   biomass = sum(biomass)) |> 
  dplyr::mutate(prop_abund = total / abundance_tot_survey,
                prop_biom = biomass / biomass_tot_survey)|> #proportion of species with NA
  dplyr::filter(prop_abund > 0.2 | prop_biom  > 0.2)

#SELECT WHEN "NAs" SPECIES REPRESENT MORE THAN 20% OF THE SURVEY, IN ABUNDANCE OR BIOMASS:
na_survey <- unique(NA_prop$survey_id)


## REMOVE NON REPRESENTATIVE ESTIMATIONS:
RLS_nut_surv[RLS_nut_surv$survey_id == na_survey, -1] <- NA


## -------------------3) Turnover of biomass -------------------####

##--- Prepping RLS data ---
#rename variables
surveys <- dplyr::select(rls_actino_trop,
                         survey_id, Species = species_name, Num = total, 
                         Sizeclass = size_class, Biomass = biomass)

sp <- dplyr::select(species_traits,
                    Species = species_name, Family = family,
                    MaxLength = Length, lwa = a, lwb = b)

metadata <- dplyr::select(all_covariates,
                          survey_id, Temperature = mean_5year_analysed_sst) |> 
  dplyr::mutate(survey_id= as.character(survey_id))

#Merge variables
RLS_clean <- surveys |> 
  dplyr::left_join(sp) |> 
  dplyr::left_join(metadata) |> 
  dplyr::rename(SurveyID = survey_id)

#add new variables
RLS_clean$Mmax = RLS_clean$lwa*RLS_clean$MaxLength^RLS_clean$lwb
RLS_clean$logMmax  = log10(RLS_clean$Mmax)
RLS_clean$logLmax  = log10(RLS_clean$MaxLength)

data_final <- RLS_clean  |> dplyr::mutate(Area=50*10)  |> dplyr::arrange(SurveyID)


#Filtering only fishable fish 
data_final_fishery <- data_final |> 
  dplyr::filter(Family %in% all_sp | (Family %in% target_larger_20cm & MaxLength > 20))

unique(data_final_fishery$Family)


##--- Run functions to predict productivity ---

#individual level#
RLS_prod_indiv = calc_prod_rfishprod(data_final_fishery)

#transect level#
RLS_prod_transect = calc_prod_transect(RLS_prod_indiv, all_covariates)
#columns: 'Prod' = biomass producted in one day on the surface of the survey
#         'Productivity' = percent of biomass producted in a day

productivuty_surveys <- RLS_prod_transect |> 
  dplyr::select(survey_id = SurveyID, Prod, productivity = Productivity)

## ------------------- Merge and save food intake indices -------------------####

food_indices_surveys <- surveys_fishery_biom |> 
  dplyr::full_join(RLS_nut_surv) |>
  dplyr::full_join(productivuty_surveys) 


## Check indices measures
summary(food_indices_surveys)

library(funbiogeo)
fb_plot_species_traits_completeness(dplyr::rename(food_indices_surveys, species = survey_id))

pca <- FactoMineR::PCA(food_indices_surveys[,-1], scale = T, graph=F, ncp=30)
factoextra::fviz_screeplot(pca, ncp = 20)
factoextra::fviz_pca_var(pca, col.var = "contrib", repel = TRUE)
factoextra::fviz_pca_var(pca, col.var = "contrib", axes=c(3,4), repel = TRUE)

# Check distributions
ggplot(data=tidyr::pivot_longer(food_indices_surveys,
                                cols = -survey_id,
                                names_to = "index", values_to = "values"), 
       aes(x=values, group=index, fill=index)) +
  geom_histogram(aes(y = ..density..), bins = 20, color = "grey40", fill ="white") +
  geom_density(aes(fill = index), alpha = 0.2) +
  hrbrthemes::theme_ipsum() +
  facet_wrap(~index, scales = "free") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )
#



save(food_indices_surveys, file = here::here("outputs", "2f_food_indices_surveys.Rdata" ))
# load(file = here::here("outputs", "2f_food_indices_surveys.Rdata" ))




