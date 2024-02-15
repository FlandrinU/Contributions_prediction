################################################################################
##
## 
##
## 2c_biomass_distribution.R
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


##-------------------biomass belonging to low and high trophic levels -------------------####
summary(data_species$Diet)

# species with low TL = herbivores_microvores_detritivores diet
species_lowTL<- data_species |>
  dplyr::filter(Diet=="herbivores_microvores_detritivores") |>
  dplyr::pull(species) 
length(species_lowTL) # 200 species

biom_lowTL <- rowSums(surveys_sp_biom[,species_lowTL]) 
summary(biom_lowTL) # from 0 to 1123316 , Q1=3643, median=21833, Q3=24280

# species with medium TL = all type of invertivorous diet (including planktivores and corrallivores)
species_mediumTL<- data_species |>
  dplyr::filter(Diet %in% c("corallivores", "sessile_invertivores", 
                            "microinvertivores", "macroinvertivores", "crustacivores", 
                            "planktivores")  ) |>
  dplyr::pull(species) 
length(species_mediumTL) # 746 species

biom_mediumTL <- rowSums(surveys_sp_biom[,species_mediumTL]) 
summary(biom_mediumTL) # from 0 to 5629533, Q1=5653, median=13216, Q3=28245


# species with high TL = piscivores diet
species_highTL<- data_species |>
  dplyr::filter(Diet=="piscivores") |>
  dplyr::pull(species) 
length(species_highTL) # 78 species

biom_highTL <- rowSums(surveys_sp_biom[,species_highTL])
summary(biom_highTL) # from 0 to 474937.1, Q1=7.7, median=685.6, Q3=2511.0

# merging
surveys_biomTL<-data.frame(biom_lowTL, biom_mediumTL, biom_highTL) |>
  tibble::rownames_to_column("SurveyID")
summary(surveys_biomTL)


# size_structure of fishes (weighted by number of individuals)
surveys_size<- data_surveys |>
  dplyr::select(SurveyID, species, size_class, number) |>
  tidyr::uncount( number ) |>
  dplyr::select( - species) |>
  dplyr::group_by(SurveyID) |>
  dplyr::summarize( size_Q1=quantile(size_class,0.25) ,
                    size_median=quantile(size_class,0.5) ,
                    size_Q3=quantile(size_class, 0.75),
                    size_p90=quantile(size_class, 0.90)
  )
summary(surveys_size)  



## phylogenetic entropy: Marcon 2015, from Allen 2009  -> very long time to run: run only on 10 trees due to negligible variations
library('entropart')
phylo_entropy_raw <- lapply( phylo_100[1:10], function(x) {
  cat("Start computation for one tree... \n")
  list <- parallel::mclapply(rownames(surveys_sp_pbiom), mc.cores = parallel::detectCores()-5, function(survey){
    community_biom <- surveys_sp_pbiom[survey,]
    phylo_entropy <- entropart::PhyloEntropy(Ps=community_biom,
                                             q = 1,
                                             Tree = x[["phy"]], 
                                             Normalize = T) #phylogenetic entropy of order 1 of relative biomass vector.
    phylo_entropy[["Total"]]
  })
  cat("one more tree done.  \n")
  do.call(rbind, list)
})


phylo_entropy_surveys_10 <- do.call(cbind, phylo_entropy_raw)
# test variability
sd <- apply(phylo_entropy_surveys_10, 1, sd)
summary(sd) # Median 0.0000918      Mean 0.0011747  3rd Qu. 0.0010619      Max. 0.0405119
#

phylo_entropy_summary <- t(apply(phylo_entropy_surveys_10, 1, summary))
colnames(phylo_entropy_summary) <- paste0("phylo_entropy_", colnames(phylo_entropy_summary))
rownames(phylo_entropy_summary) <- rownames(surveys_sp_pbiom)

save(phylo_entropy_surveys_10, file = here::here("biodiversity", "outputs", "phylogenetic_entropy_surveys_10trees.Rdata"))
save(phylo_entropy_summary, file = here::here("biodiversity", "outputs", "phylogenetic_entropy_surveys.Rdata"))

##-------------Save phylogenetic indices-------------
phylo_indices_surveys_all <- cbind(ED_surveys, PD_surveys, PE_surveys_summary, phylo_entropy_summary )
phylo_indices_surveys_all <- as.data.frame(phylo_indices_surveys_all)
phylo_indices_surveys_all <- tibble::rownames_to_column(phylo_indices_surveys_all,var="SurveyID")


save(phylo_indices_surveys_all, file = here::here("biodiversity", "outputs", "phylogenetic_indices_surveys.Rdata"))


###------------------- 2) functional entropy in surveys -------------------####
# functional entropy / richness on species relative biomass with q=1
surveys_biodiversity$funct_entropy<-mFD::alpha.fd.hill(asb_sp_w = surveys_sp_pbiom, 
                                                       sp_dist = sp_gower,
                                                       q=1, 
                                                       tau="mean",
                                                       details_returned =FALSE)[,1]