################################################################################
##
## 
##
## 2b_biodiversity_indices.R
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

#RLS observations
load(file = here::here("data/derived_data/2_rls_actino_trop.Rdata"))

#Elasmobranch presence (RLS observations)
load(file = here::here("data/derived_data/rls_elasmo_trop.Rdata"))

#SurveysXspecies matrices
load( file=here::here("data", "derived_data", "2_occurrence_matrix_sp_survey.Rdata"))
load( file=here::here("data", "derived_data", "2_relative_biom_matrix_sp_survey.Rdata"))

#Phylogenetic tree
  #chronogram of ray-finned fishes from Rabosky 2018
  # distribution of 100 all-taxon assembled (ATA) time-calibrated trees with 31,526 taxa
trees <- lapply(list.files(
  here::here("data/raw_data/phylogeny/dryad_Rabosky_2018/trees/full_trees/"), full.names = T),
  FUN = function(i){ape::read.tree(i)})
                

#Coastline map
coast <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')

#Survey metadata
load(file = here::here("data", "raw_data", "environmental_covariates",
                       "all_covariates_benthos_inferred_tropical_surveys.Rdata"))


##------------------- loading functions-------------------
source(here::here("R", "check_scientific_names.R"))


## -------------------1) taxonomic diversity in surveys -------------------####

# taxo richness of actinopterygii
taxo_richness = apply(surveys_sp_occ, 1,  sum)

# taxo richness of elasmobranchii
elasmo_by_surveys <- rls_elasmo_trop |>
  dplyr::mutate( number = 1) |>
  dplyr::group_by(survey_id) |>
  dplyr::distinct(rls_species_name, .keep_all=TRUE) |> #keep only one size class per species
  dplyr::summarise(elasmobranch_richness = sum(number))

table(elasmo_by_surveys$elasmobranch_richness)
# nb of elasmobranch species: 1   2   3   4 
# nb of surveys             946 122  21   3 


#merge richness data
surveys_richness <- data.frame(survey_id = names(taxo_richness), 
                               actino_richness = taxo_richness) |> 
  dplyr::left_join(elasmo_by_surveys)
surveys_richness$elasmobranch_richness[is.na(surveys_richness$elasmobranch_richness)] <- 0

summary(surveys_richness) #Actino:2-123, Elasmo: 0-4





##------------------- 2) functional distinctiveness and richness in surveys -------------------####

# dataframe with species as row names and traits as variables: select traits
sp_traits <- inferred_species_traits|> 
  dplyr::filter(class != "Elasmobranchii") |> 
  dplyr::select(Length, trophic_guild, Troph, DemersPelag, K, Schooling)
#Choice for traits: cf McLean 2024, Mouillot 2014, Auber 2022


# type of traits (for mFD package)
traits_cat<-data.frame(trait_name=c("Length", "trophic_guild", "Troph", "DemersPelag", "K", "Schooling"),
                       trait_type=c("Q", "N", "Q", "N", "Q", "N"))


# Gower distance between species ---
sp_gower <- mFD::funct.dist(sp_tr = sp_traits, tr_cat = traits_cat, metric = "gower",
                            stop_if_NA = F)
summary(as.numeric(sp_gower)) # most distances are small (Q3=0.40) 
# warning message means some species have same trait values
# => not an issue for computation of Chao et al 2019 indices
save(sp_gower, file = here::here("data/derived_data/2_funct_gower_distance_species.Rdata"))
# load(file = here::here("data/derived_data/2_funct_gower_distance_species.Rdata"))


## COMPUTING FUNCTIONNAL DISTINCTIVENESS ---
funct_distinctiveness_sp <- apply( as.matrix(sp_gower), 1, sum) / (ncol(as.matrix(sp_gower))-1) #Grenié et al. 2018

funct_distinct_surveys_raw <- lapply(1:nrow(surveys_sp_occ),function(i){
  Sp <- colnames(surveys_sp_occ) [which(surveys_sp_occ[i,]==1)]
  mean(funct_distinctiveness_sp[Sp])
})

surveys_funct_distinctiveness <-  do.call(rbind,funct_distinct_surveys_raw)[,1]


## COMPUTING FUNCTIONAL RICHNESS ---

# applying Chao 2019 framework with package mFD
# richness on species occurrences with q=0
surveys_funct_richness<-mFD::alpha.fd.hill(asb_sp_w = as.matrix(surveys_sp_occ), 
                                           sp_dist = sp_gower,
                                           q=0, 
                                           tau="mean",
                                           details_returned =FALSE)[,1]



## MERGING
surveys_functional <- data.frame(
  survey_id = names(surveys_funct_richness), 
  functional_richness = surveys_funct_richness) |> 
  dplyr::bind_cols(functional_distinctiveness =surveys_funct_distinctiveness)

summary(surveys_functional)






##------------------- 3) IUCN index (Actino + Elasmo) -------------------####
## Merge all rls observations (only the presence of species matter, not their abundance)
rls_trop <- rls_actino_trop |> dplyr::select(-raw_biomass) |> 
  dplyr::bind_rows(rls_elasmo_trop)

## IUCN data: IUCN categories, completed by inference by Loiseau et al. (2023)
iucn <- inferred_species_traits |> 
  dplyr::select(class, fishbase_name,
                iucn_inferred = IUCN_inferred_Loiseau23, 
                iucn_redlist = IUCN_category) 

# missing_values <- dplyr::filter(iucn, is.na(iucn$iucn_inferred)) |>
#   dplyr::mutate(iucn_redlist = dplyr::recode(iucn_redlist,  
#                                              "NE" = "No Status",
#                                              "DD" = "No Status",
#                                              "LC" = "Non Threatened",
#                                              "NT" = "Non Threatened",
#                                              "VU" = "Threatened",
#                                              "EN" = "Threatened",
#                                              "CR" = "Threatened")) 
# 
# iucn[rownames(missing_values), "iucn_inferred"] <- missing_values$iucn_redlist

## remaining NAs -> recode them into "No Status" to fit with Loiseau 2023.
iucn[is.na(iucn$iucn_inferred),] #27 species including 8 elasmobranchs

iucn$iucn_inferred[is.na(iucn$iucn_inferred)] <- "No Status"


## Merge survey data and iucn category
rls_trop_iucn <- tibble::rownames_to_column(iucn, "rls_species_name") |> 
  dplyr::right_join(rls_trop) 

summary(rls_trop_iucn$iucn_inferred)
summary(rls_trop_iucn$iucn_inferred[rls_trop_iucn$class == "Elasmobranchii"])
#Most of Elasmobranches are Threathned


## Number of IUCN species per surveys
iucn_by_surveys <- rls_trop_iucn|>
  dplyr::group_by(survey_id) |>
  dplyr::distinct(fishbase_name, .keep_all=TRUE) |> #keep only one size class per species
  dplyr::summarise(iucn_species = sum(iucn_inferred == "Threatened"),
                   elasmobranch = sum(class == "Elasmobranchii"))

table(iucn_by_surveys$iucn_species)
#nb of iucn species:   0    1    2    3    4    5    6    7 
# nb of surveys     3295 1814  624  202   49   14    3    1 



## Check the importance of NA (i.e. No Status)
NA_prop <- rls_trop_iucn|>
  dplyr::group_by(survey_id) |>
  dplyr::mutate(
    abund_tot = sum(total), #total abundance in the survey, including elasmobranchii
    biom_tot = sum(biomass)) |>  #total biomass in the survey, including elasmobranchii
  dplyr::group_by(survey_id, iucn_inferred) |> 
  dplyr::mutate(
    abund_per_iucn = sum(total),  #total abundance of TH, NT, and NS 
    biom_per_iucn = sum(biomass), #total biomass of TH, NT, and NS 
    prop_abund_per_iucn = abund_per_iucn / abund_tot, #proportion of TH, NT, and NS 
    prop_biom_per_iucn = biom_per_iucn / biom_tot) |> 
  dplyr::filter(iucn_inferred == "No Status",
                prop_abund_per_iucn > 0.2 | prop_biom_per_iucn  > 0.2) 
#SELECT WHEN "NO STATUS" SPECIES REPRESENT MORE THAN 20% OF THE SURVEY, IN ABUNDANCE OR BIOMASS:
na_survey <- unique(NA_prop$survey_id)
  


## REMOVE NON REPRESENTATIVE ESTIMATIONS:
iucn_by_surveys[iucn_by_surveys$survey_id %in% na_survey,c("iucn_species")] <- NA

table(iucn_by_surveys$iucn_species)
#nb of iucn species:   0    1    2    3    4    5    6    7 
# nb of surveys     3196 1726  602  196   47   13    3    1 
plot(c(iucn_by_surveys$iucn_species + rnorm(nrow(iucn_by_surveys), sd=0.1))~
       c(iucn_by_surveys$elasmobranch+rnorm(nrow(iucn_by_surveys), sd=0.1)))





##------------------- 4) Endemism -------------------####

## Species endemism

endemism_sp <- inferred_species_traits |>
  tibble::rownames_to_column("rls_species_name") |> 
  dplyr::select(rls_species_name, log_range = `log(geographic_range_Albouy19)`) |> 
  dplyr::mutate(range = 10^(log_range)) |> 
  dplyr::mutate(endemism_log = (max(log_range, na.rm=T) - log_range)/
                  (max(log_range, na.rm=T) - min(log_range, na.rm=T)),
                
                endemism_raw = (max(range, na.rm=T) - range)/
                  (max(range, na.rm=T) - min(range, na.rm=T))) #observe without logging the ranges

hist(endemism_sp$log_range, breaks = 20) 
hist(endemism_sp$range, breaks = 20) #highly right-skewed

hist(endemism_sp$endemism_log, breaks = 20) #good
hist(endemism_sp$endemism_raw, breaks = 20) #highly influenced by 'max(range)'

plot(endemism_sp$endemism_log~endemism_sp$endemism_raw)
#KEEP ENDEMISM ASSESSED FROM LOG-TRANSFORMED RANGE SIZE



## Survey endemic scores
endemism_survey <- as.data.frame(apply(surveys_sp_occ, 1,  sum))
colnames(endemism_survey) <- "sp_richness"

occ_endem <- surveys_sp_occ |>
  tibble::rownames_to_column("survey_id") |> 
  tidyr::pivot_longer( cols=2:last_col(),
                                names_to = "rls_species_name", values_to = "occ") |> 
  dplyr::filter(occ == 1) |> #remove useless rows: when the species is absent
  dplyr::left_join(endemism_sp) |> 
  dplyr::select(survey_id, rls_species_name, endemism_log)
  
endemism_survey <- occ_endem |> 
  dplyr::group_by(survey_id) |> 
  dplyr::summarise(sp_richness = dplyr::n(),
                   mean_endemism = mean(endemism_log, na.rm = T),
                   sum_endemism = sum(endemism_log, na.rm = T),
                   Q3_endemism = quantile(endemism_log, probs = 0.75, na.rm = T)) |> 
  dplyr::mutate(residuals_endemism_richness = lm(sum_endemism~sp_richness)[["residuals"]])


#check and plot survey endemism
hist(endemism_survey$mean_endemism, breaks = 20)
hist(endemism_survey$sum_endemism, breaks = 20)
hist(endemism_survey$Q3_endemism, breaks = 20)
hist(endemism_survey$residuals_endemism_richness, breaks = 20)

endem_metadata <- endemism_survey |> 
  dplyr::left_join(all_covariates_benthos_inferred)

library(ggplot2)
ggplot(endem_metadata) + geom_point(aes(sp_richness,sum_endemism, color = realm ))
ggplot(endem_metadata) + geom_point(aes(sp_richness,mean_endemism, color = realm ))
ggplot(endem_metadata) + geom_point(aes(sp_richness,residuals_endemism_richness, color = realm ))
ggplot(endem_metadata) + geom_point(aes(sp_richness,Q3_endemism, color = realm ))
ggplot(endem_metadata) + geom_point(aes(mean_endemism,residuals_endemism_richness, color = realm ))
#KEEP THE MEAN ENDEMISM: LESS CORRELATED WITH SP RICHNESS, AND MORE CONSISTENT
# WITH OTHER VARIABLES


## Check the importance of NA 
NA_prop <- occ_endem |> 
  dplyr::filter(is.na(endemism_log)) |> 
  dplyr::left_join(rls_actino_trop) |> 
  dplyr::group_by(survey_id, abundance_tot_survey, biomass_tot_survey) |> #keep totals at the survey scale
  dplyr::summarise(total = sum(total),
                biomass = sum(biomass)) |> 
  dplyr::mutate(prop_abund = total / abundance_tot_survey,
                prop_biom = biomass / biomass_tot_survey)|> #proportion of species with NA
  dplyr::filter(prop_abund > 0.2 | prop_biom  > 0.2) 
#SELECT WHEN "NAs" SPECIES REPRESENT MORE THAN 20% OF THE SURVEY, IN ABUNDANCE OR BIOMASS:
na_survey <- unique(NA_prop$survey_id)


## REMOVE NON REPRESENTATIVE ESTIMATIONS:
endemism <- endemism_survey |> 
  dplyr::select(survey_id, mean_endemism) |> 
  dplyr::mutate(mean_endemism = as.numeric(ifelse(survey_id %in% na_survey, NA, mean_endemism)))




##------------------- 5) Evolutionnary distinctiveness -------------------####

## extract trees with RLS reef species

#Check names in phylogeny
treesp <- trees[[1]][["tip.label"]]
names_phylogeny <- data.frame(rls_species_name =gsub("_", " ", treesp),
                              tree = 1)

# /!\ very long time to run
# subset1 <- code_sp_check(names_phylogeny[c(1:10000),], original_name = "rls_species_name", mc_cores = 12) 
# subset2 <- code_sp_check(names_phylogeny[c(10001:20000),], original_name = "rls_species_name", mc_cores = 12) 
# subset3 <- code_sp_check(names_phylogeny[c(20001:nrow(names_phylogeny)),], original_name = "rls_species_name", mc_cores = 15) 
# 
# names_phylogeny_fishbase <- rbind(subset1, subset2, subset3)
# summary(names_phylogeny_fishbase$check) #should be only 1, few mistakes
# names_phylogeny_fishbase[is.na(names_phylogeny_fishbase$fishbase_name),] #freshwater fishes, wrong names, ...
# save(names_phylogeny_fishbase, file = here::here("outputs", "2b_names_phylogeny_fishbase.Rdata"))
load(file = here::here("outputs", "2b_names_phylogeny_fishbase.Rdata"))

#Check NA
NA_fishbase <- names_phylogeny_fishbase[is.na(names_phylogeny_fishbase$fishbase_name)|
                names_phylogeny_fishbase$check == 2,] #mainly sub-species


#Extract names in the trees
names <- rep(NA, ncol(surveys_sp_occ))
sp_wrong_name <- c()

for( i in 1:ncol(surveys_sp_occ)){
  
  rls_name <- colnames(surveys_sp_occ)[i]
  fb_name <- inferred_species_traits[rls_name, "fishbase_name"]
  tree_name <- names_phylogeny_fishbase[
    which(names_phylogeny_fishbase$fishbase_name == fb_name), "rls_species_name"]
  
  if (length(tree_name) == 1) {
    names[i] <- tree_name
  } else if (length(tree_name) == 2) {
    names[i] <- tree_name[1] #choose one of the accepted name, we suppose they are very close in the phylogeny
    cat("Replace ", rls_name, " by ", tree_name[1], "\n")
    } else {
    sp_wrong_name <- c(sp_wrong_name, rls_name)
      }
}
  
 
length(sp_wrong_name) #11 species are still missing
#none of these species are in NA_fishbase


#SOME SPECIES ARE ABSENT FROM THE PHYLOGENY, AND SOME NAMES ARE DUPLICATED:
sp_occ_phylogeny <- as.data.frame(t(surveys_sp_occ)) |> 
  dplyr::mutate(species_tree = gsub(" ", "_", names)) |> 
  dplyr::filter(!is.na(species_tree)) |> #remove the 11 species not in the phylogeny
  dplyr::group_by(species_tree) |>
  dplyr::summarise(across(everything(), max)) |>  #merge identical species by keeping all occurrences
  tibble::column_to_rownames("species_tree") |> 
  t() |> 
  as.data.frame()

sp_pbiom_phylogeny <- as.data.frame(t(surveys_sp_pbiom)) |> 
  dplyr::mutate(species_tree = gsub(" ", "_", names)) |> 
  dplyr::filter(!is.na(species_tree)) |> 
  dplyr::group_by(species_tree) |>
  dplyr::summarise(across(everything(), max)) |> 
  tibble::column_to_rownames("species_tree") |> 
  t() |> 
  as.data.frame()

save(sp_occ_phylogeny, file = here::here("outputs", "2b_sp_occ_phylogeny_matrix.Rdata"))
save(sp_pbiom_phylogeny, file = here::here("outputs", "2b_sp_pbiom_phylogeny_matrix.Rdata"))

#Extract trees
phylo_100 <- list()
for (i in 1:100) {
  phylo_100[[i]] <- picante::match.phylo.comm(trees[[i]], sp_occ_phylogeny)
}

#Occurrence matrix in RLS surveys
occ_matrix <- as.matrix(phylo_100[[1]][["comm"]])


## Compute Evolutionary distinctivness: ED  Isaac et al. method

# by species -> long time to run
ED_species_raw <- parallel::mclapply(phylo_100, mc.cores=10, function(x) {
  picante::evol.distinct(x$phy, type = c("fair.proportion"), scale = FALSE, use.branch.lengths = TRUE)})

save(ED_species_raw, file = here::here("outputs", "2b_evolutionary_distinctivness_species.Rdata"))
# load(file = here::here("outputs", "2b_evolutionary_distinctivness_species.Rdata"))

ED_species <- as.data.frame(do.call(cbind, lapply(ED_species_raw,function(y){y[,2]})))
rownames(ED_species) <- phylo_100[[1]]$phy$tip.label

ED_species_summary <- t(apply(ED_species, 1, summary))
ED_species_summary <- cbind( ED_species_summary, sd = apply(ED_species, 1, sd))


#by surveys
median_ED_sp <- ED_species_summary[,"Median"]
ED_surveys_raw <- parallel::mclapply(1:nrow(occ_matrix), mc.cores=parallel::detectCores()-5 ,function(i){
  if(sum(occ_matrix[i,])==0){
    rep(0,6)
  }else{
    Sp <- colnames(occ_matrix) [which(occ_matrix[i,]==1)]
    c(rownames(occ_matrix)[i], summary(median_ED_sp[Sp]))
  }
})

ED_surveys <- do.call(rbind,ED_surveys_raw)
colnames(ED_surveys) <- c("survey_id", "ED_min", "ED_Q1", "ED_median", "ED_mean",
                          "ED_Q3", "ED_max")



## Check the importance of NA 
NA_prop <- rls_actino_trop|> 
  dplyr::filter(rls_species_name %in% sp_wrong_name) |> 
  dplyr::group_by(survey_id, abundance_tot_survey, biomass_tot_survey) |> #keep totals at the survey scale
  dplyr::summarise(total = sum(total),
                   biomass = sum(biomass)) |> 
  dplyr::mutate(prop_abund = total / abundance_tot_survey,
                prop_biom = biomass / biomass_tot_survey)|> #proportion of species with NA
  dplyr::filter(prop_abund > 0.2 | prop_biom  > 0.2) 
#SELECT WHEN "NAs" SPECIES REPRESENT MORE THAN 20% OF THE SURVEY, IN ABUNDANCE OR BIOMASS:
na_survey <- unique(NA_prop$survey_id)


## REMOVE NON REPRESENTATIVE ESTIMATIONS:
evol_distinct <- as.data.frame(ED_surveys) |> 
  dplyr::select(survey_id, ED_mean) |> 
  dplyr::mutate(ED_mean = as.numeric(ifelse(survey_id %in% na_survey, NA, ED_mean)))

hist(evol_distinct$ED_mean)





##------------- Merge and save all biodiversity indices -------------####

biodiv_indices_surveys <- surveys_richness |> 
  dplyr::full_join(surveys_functional) |> 
  dplyr::full_join(iucn_by_surveys) |> 
  dplyr::full_join(endemism) |> 
  dplyr::full_join(evol_distinct) |> 
  dplyr::select(-elasmobranch, -functional_richness)

## Check indices measures
summary(biodiv_indices_surveys)

library(funbiogeo)
fb_plot_species_traits_completeness(dplyr::rename(biodiv_indices_surveys, species = survey_id))

pca <- FactoMineR::PCA(biodiv_indices_surveys[,-1], scale = T, graph=F, ncp=30)
factoextra::fviz_screeplot(pca, ncp = 20)
factoextra::fviz_pca_var(pca, col.var = "contrib", repel = TRUE)
factoextra::fviz_pca_var(pca, col.var = "contrib", axes=c(3,4), repel = TRUE)

# Check distributions
ggplot(data=tidyr::pivot_longer(biodiv_indices_surveys,
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



save(biodiv_indices_surveys, file = here::here("outputs", "2b_biodiv_indices_surveys.Rdata" ))
# load(file = here::here("outputs", "2b_biodiv_indices_surveys.Rdata" ))

  
  
  
  