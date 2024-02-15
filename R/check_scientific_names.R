#' check_scientific_names.R
#'
#' @param data A data frame with species in rows and a column "species_name", with names in this format: "Mola mola"
#'
#' @return
#' @export
#'
#' @examples


#data = species_list
#var <- sample(unique(data$species_name), 30)
# check_ind <- synonyms(str_replace("Acanthobrama_terraesanctae", "_", " "),version = "19.04")
#i = 5908
#

code_sp_check <- function(data, mc_cores = 1){

  var <- unique(data$species_name)
  
  check <-  do.call(rbind, pbmcapply::pbmclapply(1:length(var),function(i){ 
    
    if( i %% 500 == 0 ){ Sys.sleep(sample(1:60, 1))} # Avoid problem with API request
    
    # For sub species keep the first and the last names
    if(is.na(stringr::word(var[i],3,sep=" "))){
      # sp <- "Raja africana"
      sp <- var[i]
      
    }else{
      
      sp <- paste0(stringr::word(var[i], 1,sep=" "),
                   " ", stringr::word(var[i], 3,sep=" "))}
    
    print(paste0(i,"/",length(var),", ",round(i/length(var),3)*100,"%")) #need to use only 1 core to have this message
    
    check_ind <- rfishbase::synonyms(sp ) #, version = "19.04")
    
    if(sum(is.na(check_ind[1,])) == 9){
      check_ind <-  data.frame(species_name = var[i],
                               spec_code = NA, 
                               worms_id = NA,
                               fishbase_name = NA,
                               check = NA)
      
    }else{
      
      check_ind  <- check_ind[check_ind$Status == "accepted name" |
                                check_ind$Status == "synonym" |
                                check_ind$Status == "provisionally accept"|
                                check_ind$Status == "provisionally accepted name"|
                                check_ind$Status == "ambiguous synonym",]
      
      if(nrow(check_ind) == 1 && check_ind$Species == "Genus Species") {
        check_ind$Species  <- check_ind$synonym }
      if(sum(check_ind$Species %in%  "Genus Species") >0) {
        check_ind  <- check_ind[check_ind$Species != "Genus Species",] }
      if(nrow(check_ind)>1 && sum(check_ind$Status  %in% "accepted name") == 1) {
        check_ind  <- check_ind[check_ind$Status == "accepted name",] }
      if(nrow(check_ind)>1 && sum(!check_ind$Status  %in% "accepted name") >0 && sum(check_ind$Status  %in% "provisionally accept")>0) {
        check_ind  <- check_ind[check_ind$Status == "provisionally accept",] }
      if(nrow(check_ind)>1 && sum(!check_ind$Status  %in% "accepted name")>0 && sum(check_ind$Status  %in% "provisionally accepted name")>0) {
        check_ind  <- check_ind[check_ind$Status == "provisionally accepted name",] }
      if(nrow(check_ind)>1 && sum(!check_ind$Status  %in% "accepted name")>0 && sum(!check_ind$Status  %in% "provisionally accept")>0) {
        check_ind  <- check_ind[check_ind$Status == "synonym",] }
      if(nrow(check_ind)>1 && sum(!check_ind$Status  %in% "accepted name")>0 && sum(!check_ind$Status  %in% "provisionally accepted name")>0) {
        check_ind  <- check_ind[check_ind$Status == "synonym",] }
     
      if(nrow(check_ind) == 0){
        check_ind <-  data.frame(species_name = var[i],
                                 SpecCode = NA, 
                                 WoRMS_ID = NA,
                                 Species = NA) }
      
      check_ind <-  data.frame(species_name = var[i],
                               spec_code = check_ind$SpecCode,
                               worms_id = check_ind$WoRMS_ID,
                               fishbase_name = check_ind$Species,
                               check = nrow(check_ind)) # check if very not a problem in the filter accepted, synonym...should be 1
    } # end of else
    
  }, mc.cores = mc_cores)) #end of lapply
  
  res <- merge(data,check, by = "species_name")
  return(res)
} # end of function code_sp_check
