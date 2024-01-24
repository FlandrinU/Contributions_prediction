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
################################################################################

## cleaning memory
rm(list=ls())

##------------------- loading datasets-------------------
#Species traits
load(file= here::here("outputs", "RLS_species_traits_inferred.Rdata"))

#RLS observations
load(file = here::here("data/derived_data/rls_actino_trop.Rdata"))

#Elasmobranch presence (RLS observations)
load(file = here::here("data/derived_data/rls_elasmo_trop.Rdata"))

#SurveysXspecies matrices
load( file=here::here("data", "derived_data", "2_occurrence_matrix_sp_survey.Rdata"))
load( file=here::here("data", "derived_data", "2_relative_biom_matrix_sp_survey.Rdata"))

#Phylogenetic tree
trees <- ape::read.tree(file = here::here("biodiversity", "data","Code&Data_Siquiera2020",
                                          "TACT", "Reef_fish_all_combined.trees")) #chronogram of ray-finned fishes, Siquiera 2020 from Rabosky 2018
#Coastline map
coast <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')



### -------------------1) taxonomic diversity in surveys -------------------####

# taxo richness of actinopterygii
taxo_richness = apply(surveys_sp_occ, 1,  sum)

# taxo richness of elasmobranchii
elasmo_by_surveys <- rls_elasmo_trop |>
  dplyr::mutate( number = 1) |>
  dplyr::group_by(survey_id) |>
  dplyr::distinct(species_name, .keep_all=TRUE) |> #keep only one size class per species
  dplyr::summarise(elasmobranch_diversity = sum(number))

table(elasmo_by_surveys$elasmobranch_diversity)
# nb of elasmobranch species: 1   2   3   4 
# nb of surveys             576  92  17   1 


#merge richness data
surveys_richness <- data.frame(survey_id = names(taxo_richness), 
                               taxo_richness = taxo_richness) |> 
  dplyr::left_join(elasmo_by_surveys)
surveys_richness$elasmobranch_diversity[is.na(surveys_richness$elasmobranch_diversity)] <- 0



######################################################################################
##------------------- 2) functional diversity in surveys -------------------####

# dataframe with species as row names and traits as variables
sp_traits <- infered_data
sp_traits <- sp_traits |> 
  dplyr::select(Length, trophic_guild, Position, Activity)

# type of traits
traits_cat<-data.frame(trait_name=c("Size", "Diet", "DemersPelag", "Activity"),
                       trait_type=c("Q", "N", "O", "N")
)

# Gower distance between species ----
sp_gower <- mFD::funct.dist(sp_tr = sp_traits, tr_cat = traits_cat, metric = "gower")
summary(as.numeric(sp_gower)) # most distances are small (Q3=0.32) 
# warning message means some species have same trait values
# => not an issue for computation of Chao et al 2019 indices


# computing functional distinctiveness
funct_distinctiveness_sp <- apply( as.matrix(sp_gower), 1, sum) / (ncol(as.matrix(sp_gower))-1) #Grenié et al. 2018

funct_distinct_surveys_raw <- lapply(1:nrow(surveys_sp_occ),function(i){
  Sp <- colnames(surveys_sp_occ) [which(surveys_sp_occ[i,]==1)]
  mean(funct_distinctiveness_sp[Sp])
})

surveys_biodiversity$funct_distinctiveness <-  do.call(rbind,funct_distinct_surveys_raw)[,1]

# computing functional richness and functional entropy ----
# applying Chao 2019 framework with package mFD

# richness on species occurrences with q=0
surveys_biodiversity$funct_richness<-mFD::alpha.fd.hill(asb_sp_w = surveys_sp_occ, 
                                                        sp_dist = sp_gower,
                                                        q=0, 
                                                        tau="mean",
                                                        details_returned =FALSE)[,1]

# richness on species relative biomass with q=1
surveys_biodiversity$funct_entropy<-mFD::alpha.fd.hill(asb_sp_w = surveys_sp_pbiom, 
                                                       sp_dist = sp_gower,
                                                       q=1, 
                                                       tau="mean",
                                                       details_returned =FALSE)[,1]

summary(surveys_biodiversity)








# merging 
surveys_biodiversity <- surveys_biodiversity |>
  tibble::rownames_to_column("SurveyID")  |>
  dplyr::left_join(surveys_size) |>
  dplyr::left_join(surveys_biomTL)
summary(surveys_biodiversity)






##------------------- 3) IUCN index -------------------####


## merge data
all_sp <- rbind( dplyr::select(data_species, species, species_corrected),
                 dplyr::select(data_species_elasmo, species, species_corrected))

names <- questionr::na.rm(all_sp$species_corrected) # 3 elasmobranch species identifies only at the genus level

## IUCN redlist data
iucn_name <- gsub("_", " ", names)
iucn_data_raw <- parallel::mclapply(iucn_name, mc.cores = parallel::detectCores()-5,
                                    function(i){rredlist ::rl_search(name = i,  key= IUCN_KEY)}) #/!\ long time to run
iucn_data <- do.call(rbind, lapply(iucn_data_raw, "[[", "result"))

table(iucn_data$category) # 3 CR, 35 DD, 10 EN, 819 LC, 18 NT, 34 VU

save(iucn_data, file = here::here("biodiversity", "outputs", "iucn_data_all_species.Rdata"))
# load(file = here::here("biodiversity", "outputs", "iucn_data_all_species.Rdata"))

#keep the iucn category
iucn_category <- iucn_data |>
  dplyr::mutate(scientific_name = gsub(" ", "_", iucn_data$scientific_name)) |>
  dplyr::select( species_corrected = scientific_name, category) |>
  dplyr::right_join(all_sp)

table(iucn_category$category, useNA = "always") # 3 CR, 35 DD, 10 EN, 818 LC, 18 NT, 34 VU, 166 <NA>


#Complete category by random forrest and deep learning (cf N.Loiseau)
dat_network$IUCN_final <- as.character(dat_network$IUCN_final)

for(i in which(is.na(iucn_category$category))){
  name <- iucn_category$species[i]
  name_corrected <- iucn_category$species_corrected[i]
  
  if(length(which(dat_network$species == name))==1){
    iucn_category$category[i] <- dat_network$IUCN_final[which(dat_network$species == name)] 
  }else{ if(length(which(dat_network$species == name_corrected))==1){
    iucn_category$category[i] <- dat_network$IUCN_final[which(dat_network$species == name_corrected)]
  }}
}

table(iucn_category$category, useNA = "always") 
# 3 CR, 34 DD,  10 EN, 755 LC,  6 No Status, 202 Non Threatened, 18 NT, 17 Threatened, 32 VU, 7 NA



#merge survey data and iucn category
all_surveys_iucn <- rbind( dplyr::select(data_surveys, SurveyID, species, size_class, number, biomass),
                           dplyr::select(data_surveys_elasmo, SurveyID, species, size_class, number, biomass)) |>
  dplyr::left_join(iucn_category) 


#remove species without iucn category
all_surveys_iucn <- all_surveys_iucn |>
  dplyr::filter(is.na(all_surveys_iucn$category) == F,
                category != "No Status")

## Number of IUCN species per surveys
iucn_category_surveys <- all_surveys_iucn |>
  dplyr::mutate(category = forcats::fct_recode(category, 
                                               "0" = "LC",
                                               "0" = "NE",
                                               "0" = "DD",
                                               "0" = "NT",
                                               "0" = "Non Threatened",
                                               "1" = "VU",
                                               "1" = "EN",
                                               "1" = "CR",
                                               "1" = "Threatened")) 


iucn_by_surveys <- iucn_category_surveys|>
  dplyr::mutate(category = as.numeric(as.character(category))) |>
  dplyr::group_by(SurveyID) |>
  dplyr::distinct(species, .keep_all=TRUE) |> #keep only one size class per species
  dplyr::summarise(iucn_species = sum(category))

table(iucn_by_surveys$iucn_species)
#nb of iucn species:   0    1    2    3    4    5    6 
# nb of surveys     2244  979  304   78   18    3    1  







##------------------- 4) Endemism -------------------####



##-------------compute species endemism-------------
# Deal with species names
names <- gsub(" ", "_", colnames(mat_PA_teleo[, -c(1,2)]))

for( i in 1:length(names)){
  if(names[i] %in% data_species$species){
    names[i] <- data_species$species_corrected[which(data_species$species == names[i])]
  }
}

colnames(mat_PA_teleo) <- c("Longitude", "Latitude", names)
mat_PA_rls <- mat_PA_teleo[, which(colnames(mat_PA_teleo) %in% data_species$species_corrected)]


## species range
range_sp <- as.data.frame(colSums(mat_PA_rls)) 
colnames(range_sp) <- "range"


## Endemism
endemism_sp <- range_sp |>
  dplyr::mutate(endemism = (max(range) - range)/(max(range) - min(range)))
hist(endemism_sp$endemism, breaks = 20)


##-------------survey endemic score------------- = mean of species endemism in a given survey
endemism_survey <- rep(NA, nrow(surveys_sp_occ))

for(i in 1:nrow(surveys_sp_occ)){
  Names <- names(surveys_sp_occ[i,which(surveys_sp_occ[i,] == 1)])
  corrected_names <- dplyr::filter(data_species, species %in% Names)$species_corrected
  endemism_survey[i] <- mean(endemism_sp[corrected_names, "endemism"], na.rm=T)
}

endemism_survey_rls <- data.frame(SurveyID = rownames(surveys_sp_occ), Endemism = endemism_survey)
save(endemism_survey_rls, file = here::here("biodiversity", "outputs", "survey_endemism_score.Rdata"))


# #check and plot survey endemism
# hist(endemism_survey_rls$Endemism, breaks = 20)
# endemism_survey_map <- endemism_survey_rls |>
#   dplyr::left_join( dplyr::select(metadata_surveys, SurveyID, SiteLongitude, SiteLatitude))
# 
# ggplot(endemism_survey_map) +
#   geom_sf(data = coast, color = "grey30", fill = "lightgrey",
#           aes(size=0.1)) +
#   geom_point(data=endemism_survey_map,
#              size = 4, shape = 20,
#              aes(x = SiteLongitude, y = SiteLatitude,
#                  colour= Endemism)) +
#   scale_colour_gradient(low = "dodgerblue", high="darkred",
#                         na.value="black") +
#   guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
#   theme_minimal()+
#   labs(title = paste0("Endemism", " geographic distribution"),
#        x="", y= "") +
#   theme(legend.position = "right",
#         plot.title = element_text(size=10, face="bold"),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
#   )
# 
# ggsave(plot = last_plot(), filename = here::here("biodiversity", "figures", "endemism_on_world_map.jpg"),
#        width = 15, height = 10)




##------------------- 5) Evolutionnary distinctiveness -------------------####


## -------------extract trees with RLS reef species-------------
## Change names in matrix to fit with names in phylogeny
treesp <- trees[[1]][["tip.label"]]
names <- rep(NA, ncol(surveys_sp_occ))
sp_wrong_name <- c()
for( i in 1:ncol(surveys_sp_occ)){
  if(colnames(surveys_sp_occ)[i] %in% treesp){
    names[i] <- colnames(surveys_sp_occ)[i]
    
  }else if(data_species$species_corrected[
    which(data_species$species == colnames(surveys_sp_occ)[i])] %in% treesp){
    names[i] <- data_species$species_corrected[which(data_species$species == colnames(surveys_sp_occ)[i])]
    
  }else{
    sp_wrong_name <- c(sp_wrong_name, colnames(surveys_sp_occ)[i])
    cat(colnames(surveys_sp_occ)[i], " not in phylogeny \n")
  }
}

#Is there a synonym names in the tree?
for( sp in sp_wrong_name){
  sp <- stringr::str_replace(sp, "_", " ")
  tab_syn <- taxize::synonyms(sp , db="worms")
  if( nrow(tab_syn[[sp]]) > 0){
    sp_syn <- stringr::str_replace(tab_syn[[sp]][["scientificname"]], " ", "_")
    cat(length(sp_syn), "synonyms", "\n" )
    cat( "new name : ", sp_syn[which(is.element(sp_syn, treesp)==T)], "\n" )
  }
} #No synonyms found in the tree

colnames(surveys_sp_occ) <- colnames(surveys_sp_pbiom) <- names

surveys_sp_occ <- as.data.frame(surveys_sp_occ) |>
  dplyr::select( dplyr::contains( "_")) #Remove species absent form phylogeny

surveys_sp_pbiom <- as.data.frame(surveys_sp_pbiom) |>
  dplyr::select( dplyr::contains( "_")) |> #Remove species absent form phylogeny
  as.matrix()


#Extract trees
phylo_100<-list()
for (i in 1:100) {
  phylo_100[[i]]<-picante::match.phylo.comm(trees[[i]], surveys_sp_occ)
}

#Occurrence matrix in RLS surveys
occ_matrix <- as.matrix(phylo_100[[1]][["comm"]])


##-------------Compute phylogenetic indices-------------

## Evolutionary distinctivness: ED  Isaac et al. method
#by species
ED_species_raw<-parallel::mclapply(phylo_100, mc.cores=parallel::detectCores()-5, function(x) {
  picante::evol.distinct(x$phy, type = c("fair.proportion"), scale = FALSE, use.branch.lengths = TRUE)})

ED_species <- as.data.frame(do.call(cbind, lapply(ED_species_raw,function(y){y[,2]})))
rownames(ED_species) <- phylo_100[[1]]$phy$tip.label

ED_species_summary <- t(apply(ED_species, 1, summary))
ED_species_summary <- cbind( ED_species_summary, sd = apply(ED_species, 1, sd))
colnames(ED_species_summary) <- paste0("ED_", colnames(ED_species_summary))

save(ED_species_summary, file = here::here("biodiversity", "outputs", "evolutionary_distinctivness_species.Rdata"))

#by surveys
mean_ED_sp <- apply(ED_species, 1, mean)
ED_surveys_raw <- parallel::mclapply(1:nrow(occ_matrix), mc.cores=parallel::detectCores()-5 ,function(i){
  if(sum(occ_matrix[i,])==0){
    rep(0,6)
  }else{
    Sp <- colnames(occ_matrix) [which(occ_matrix[i,]==1)]
    summary(mean_ED_sp[Sp])
  }
})

ED_surveys <- do.call(rbind,ED_surveys_raw)
row.names(ED_surveys) <- row.names(occ_matrix)
colnames(ED_surveys) <- paste0("ED_", colnames(ED_surveys))

save(ED_surveys, file = here::here("biodiversity", "outputs", "evolutionary_distinctivness_surveys.Rdata"))



## PD and SES.PD quantification with Phyloregion
#PD
sparse_occ_matrix <- methods::as( as.matrix(occ_matrix), "sparseMatrix")
PD_surveys_raw <- lapply(phylo_100, function (x){ phyloregion::PD(sparse_occ_matrix, x$phy) })
PD_surveys_100 <- t( do.call(rbind, PD_surveys_raw ) )

PD_surveys_summary <- t(apply(PD_surveys_100, 1, summary))
colnames(PD_surveys_summary) <- paste0("PD_", colnames(PD_surveys_summary))


#Residuals of PD ~ taxonomic richness
Mean_PD <- apply(PD_surveys_100,1,mean)
taxo_richness <- apply(occ_matrix, 1, sum)
residuals_PD_richness <- residuals(lm(Mean_PD ~ taxo_richness)) # Approach using residuals

# #SES.PD
# SES_PD_raw <- parallel::mclapply(phylo_100, mc.cores= parallel::detectCores()-5, function (x){
#   phyloregion::PD_ses(sparse_occ_matrix, x$phy, model = c("tipshuffle"), reps= 1000)$zscore #Change reps for a quickier analyse
# }) # Approach using SES.PD based on a null model shuffling tip labels                        
# 
# SES_PD_surveys_100 <- do.call(cbind, SES_PD_raw)
# SES_PD_surveys_summary <- t(apply(SES_PD_surveys_100, 1, summary))
# colnames(SES_PD_surveys_summary) <- paste0("SES_PD_", colnames(SES_PD_surveys_summary))
# 

PD_surveys <- cbind( PD_surveys_summary, residuals_PD_richness = residuals_PD_richness)#, SES_PD_surveys_summary)
save(PD_surveys, file = here::here("biodiversity", "outputs", "phylogenetic_diversity_surveys.Rdata"))


## PE: Phylogenetic endemism (pkg Pyloregion)
PE_surveys_raw <- parallel::mclapply(phylo_100, mc.cores=parallel::detectCores()-5, function(x) {
  phyloregion::phylo_endemism(sparse_occ_matrix, x$phy, weighted = TRUE)
})

PE_surveys_100 <- do.call(cbind, PE_surveys_raw)
PE_surveys_summary <- t(apply(PE_surveys_100, 1, summary))
colnames(PE_surveys_summary) <- paste0("PE_", colnames(PE_surveys_summary))

save(PE_surveys_summary, file = here::here("biodiversity", "outputs", "phylogenetic_endemism_surveys.Rdata"))






##-------------Save phylogenetic indices-------------
phylo_indices_surveys_all <- cbind(ED_surveys, PD_surveys, PE_surveys_summary, phylo_entropy_summary )
phylo_indices_surveys_all <- as.data.frame(phylo_indices_surveys_all)
phylo_indices_surveys_all <- tibble::rownames_to_column(phylo_indices_surveys_all,var="SurveyID")





##-------------merge and save all biodiversity indices-------------
surveys_richness