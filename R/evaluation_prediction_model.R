################################################################################
##
##  
##
## evaluation_prediction_model.R
##
## 29/11/2022
##
## Ulysse Flandrin
##
################################################################################
# ##-----------------Loading packages-------------------
# pkgs <- c("here", "dplyr", "missForest", "pbmcapply", "patchwork", "ggplot2")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

library(patchwork)
library(ggplot2)

##-------------1) Extract model result -------------

extract_model_perf <- function(raw_result = model_eval_missforest){
  flat_list <- unlist(raw_result, recursive = F)
  
  ### Raw data ###
  raw_estimates_factor <- do.call(rbind,flat_list[seq(1,length(flat_list),2)])
  raw_estimates_num <- do.call(rbind,flat_list[seq(2,length(flat_list),2)])
  
  ### Extract evaluation measure per traits ###
  traits_perf = list()
  
  ##Factors
  for(i in seq(1,length(flat_list),2)){ 
    model_i <- list()
    temp_fact <- flat_list[[i]]
    
    for( var in unique(temp_fact$variable)){
      temp_fact_trait <- temp_fact |> dplyr::filter(variable == var)
      model_number <- unique(temp_fact_trait$missForest)
      
      m_test <- mean(temp_fact_trait$observed==temp_fact_trait$imputed) #Measure the 'Accuracy' of the model
      
      #number of NA created in the trait
      na_created <- nrow(temp_fact_trait)
      
      #Saving accuracy measure into list
      model_i[[var]] <- c(m_test, na_created, model_number)
    }
    
    res <- as.data.frame(do.call(rbind,model_i)) |> 
      dplyr::rename(estimate = V1, NA_created = V2, model = V3) |> 
      dplyr::mutate(method = "accuracy") |> 
      tibble::rownames_to_column(var = "variable") 
    
    traits_perf[[i]] <- res
  } #End factor estimates per trait
  
  
  #Numeric
  for(i in seq(2,length(flat_list),2)){ 
    model_i <- list()
    temp_num <- flat_list[[i]]
    
    for( var in unique(temp_num$variable)){
      temp_num_trait <- temp_num |> dplyr::filter(variable == var)
      model_number <- unique(temp_num_trait$missForest)
      
      # lm_test = cor(temp_num_trait$observed, temp_num_trait$imputed) # Pearson correlation between observed and inferred
      lm_test = summary(lm(temp_num_trait$observed ~ temp_num_trait$imputed))[["r.squared"]] #R-squared evaluation
      
      #number of NA created in the trait
      na_created = nrow(temp_num_trait)
      
      #Saving accuracy measure into list
      model_i[[var]] <- c(lm_test, na_created, model_number)
    }
    
    res <- as.data.frame(do.call(rbind,model_i)) |> 
      dplyr::rename(estimate = V1, NA_created = V2, model = V3) |> 
      dplyr::mutate(method = "r_squared") |> 
      tibble::rownames_to_column(var = "variable")
    
    traits_perf[[i]] <- res
  } #End num estimates per trait
  
  traits_performance = as.data.frame(do.call(rbind,traits_perf))
   
   
  
  ### Extract evaluation measure per order ###
  order_perf = list()
  
  ##Factors
  for(i in seq(1,length(flat_list),2)){ 
    model_i <- list()
    temp_fact <- flat_list[[i]]
    
    for( taxa in unique(temp_fact$order)){
      temp_fact_taxa <- temp_fact |> dplyr::filter(order == taxa)
      model_number <- unique(temp_fact_taxa$missForest)
      
      m_test <- mean(temp_fact_taxa$observed==temp_fact_taxa$imputed) #Measure the 'Accuracy' of the model
      
      #number of NA created in the trait
      na_created <- nrow(temp_fact_taxa)
      
      #Saving accuracy measure into list
      model_i[[taxa]] <- c(m_test, na_created,model_number)
    }
    
    res <- as.data.frame(do.call(rbind,model_i)) |> 
      dplyr::rename(estimate = V1, NA_created = V2, model = V3) |> 
      dplyr::mutate(method = "accuracy") |> 
      tibble::rownames_to_column(var = "order") 
    
    order_perf[[i]] <- res
  } #End factor estimates per order
  
  
  #Numeric
  for(i in seq(2,length(flat_list),2)){ 
    model_i <- list()
    temp_num <- flat_list[[i]] |> 
      dplyr::group_by(variable) |> 
      dplyr::mutate(mean=mean(observed),
                    sd = sd(observed)) |> 
      dplyr::mutate(observed= (observed - mean) / sd ,
                    imputed = (imputed - mean) / sd ) #scale independently all variables to be comparable
    
    for( taxa in unique(temp_num$order)){
      temp_num_taxa <- temp_num |> dplyr::filter(order == taxa)
      model_number <- unique(temp_num_taxa$missForest)
      
      lm_test = NA
      if(nrow(temp_num_taxa)>2){
        # lm_test = cor(temp_num_taxa$observed, temp_num_taxa$imputed) # Pearson correlation between observed and inferred
        lm_test = summary(lm(temp_num_taxa$observed ~ temp_num_taxa$imputed))[["r.squared"]] #R-squared evaluation
        }
      
      #number of NA created in the trait
      na_created = nrow(temp_num_taxa)
      
      #Saving accuracy measure into list
      model_i[[taxa]] <- c( lm_test, na_created, model_number)
    }
    
    res <- as.data.frame(do.call(rbind,model_i)) |> 
      dplyr::rename(estimate = V1, NA_created = V2, model = V3) |> 
      dplyr::mutate(method = "r_squared") |> 
      tibble::rownames_to_column(var = "order")
    
    order_perf[[i]] <- res
  } #End num estimates per order
  
  order_performance = as.data.frame(do.call(rbind,order_perf)) 
   
  
  ### All results ###
  all_res <- list(traits_performance, order_performance, 
                raw_estimates_factor, raw_estimates_num)
  return(all_res)
  
} # END of function extract_model_perf


##-------------2) Plot estimates boxplot-------------

## BOXPLOT BY TRAITS
estimates_boxplot <- function(df_estimates = traits_performance,
                              add_mean = T){
  library(ggplot2)
  #Colors
  col <- rev(fishualize::fish(n = length(unique(df_estimates$variable)), 
                          option = "Ostracion_whitleyi", begin = 0.2, end = 0.9))
  
  #Mean estimate
  mean_estimate <- mean(df_estimates$estimate)
  
  #order boxplot
  order_trait <- df_estimates |>
    dplyr::group_by(variable) |>
    dplyr::summarise(median_estimate = median(estimate))
  
  order_boxplot <- order_trait$variable[order(order_trait$median_estimate)]
  
  #merge data
  data <- merge(df_estimates, order_trait)
  data$variable <- factor(data$variable, levels = order_boxplot)
  
  plot <- ggplot(data) +
    aes(x= variable, y= estimate, fill = variable, col = method)+
    geom_boxplot() +
    ylim(0,1)+
    
    # Add the mean above each boxplot
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "text", 
                 aes(label = paste(round(after_stat(y), 2))), 
                 position = position_dodge(width = 0.75), vjust = -2, size = 4,
                 color = "grey50") + 

    
    #Set ggplot theme
    scale_color_manual(values = c("grey", "black"))+
    scale_fill_manual(values = col)+
    xlab("") + ylab("Predictive power") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 12))
  
  #Add the mean prediction
  if(add_mean){
    plot <- plot +
    geom_hline(yintercept = mean_estimate, linetype = "dashed", 
               color = "coral3") +
      annotate("text", x = length(unique(data$variable))-1, y = mean_estimate, 
               label = paste("Mean:", round(mean_estimate, 2)), vjust = -1,
               color = "coral3") 
  }
  
  plot
  
} # END of estimates_boxplot (2a)

## BOXPLOT BY ORDER
estimates_boxplot_per_order <- function(df_estimates = order_performance){
  #Colors
  col <- rev(fishualize::fish(n = length(unique(df_estimates$order)), 
                              option = "Ostracion_whitleyi", begin = 0.2, end = 0.9))
  
  ### Continuous traits
  df_estimates_num <- dplyr::filter(df_estimates, method == "r_squared")
  #order boxplot
  order_taxa <- df_estimates_num |>
    dplyr::group_by(order) |>
    dplyr::summarise(median_estimate = median(estimate, na.rm = T))
  
  order_boxplot <- order_taxa$order[order(order_taxa$median_estimate)]
  
  #merge data
  data <- merge(df_estimates_num, order_taxa)
  data$order <- factor(data$order, levels = order_boxplot)
  
  plot_r_squared <-ggplot(data) +
    aes(x= order, y= estimate, fill = order)+
    geom_boxplot(col = "black") +
    scale_fill_manual(values = col)+
    xlab("") + ylab("Assessement quality (R-squared)") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 0))
  
  ### Categorial traits
  df_estimates_trait <- dplyr::filter(df_estimates, method == "accuracy")

  #merge data
  data <- merge(df_estimates_trait, order_taxa)
  data$order <- factor(data$order, levels = order_boxplot)
  
  plot_accuracy <- ggplot(data) +
    aes(x= order, y= estimate, fill = order)+
    geom_boxplot(col = "grey") +
    scale_fill_manual(values = col)+
    xlab("") + ylab("Assessement quality (Accuracy)") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 10))
  
  plot <- plot_r_squared / plot_accuracy
  plot
  
} # END of estimates_boxplot (2a)

##-------------3) Plot estimates histogram-------------
estimates_histogramm <- function(data = traits_performance){

  plot_distribution <- function(trait, data){
    col <- fishualize::fish(n = length(unique(data$variable)), option = "Ostracion_whitleyi", begin = 0.2, end = 0.9)
    names(col) <- unique(data$variable)[order(unique(data$variable))]
    
    data <- dplyr::filter(data, variable == trait)
    
    ggplot(data) +
      aes(x = estimate) +
      geom_histogram(bins = 30,
                     fill = col[trait][[1]],
                     col = "black") +
      xlim(0,1)+
      labs(title = trait) +
      xlab("") + ylab("") +
      theme_minimal() +
      theme( legend.position = "none")+
      geom_vline(xintercept = median( data$estimate),
                 linetype="dashed", 
                 color = "black", linewidth=1) +
      annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1,
               label=paste("median = ", round(median( data$estimate), 2) ),
               col="firebrick4")

  }
  
  
  plots <- lapply(unique(data$variable)[order(unique(data$variable))], FUN = plot_distribution, data )

  vect_plots <- c()
  for(i in c(1:length(plots))){ vect_plots <- c(vect_plots, paste0("plots[[", i, "]]")) }
  code_expr <- paste(vect_plots, collapse = "+")
  
  all_plot <- eval(parse(text=code_expr)) +
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  
  all_plot
} # END of estimates_histogramm (1b)



estimates_initial_NA <- function(NA_proportion = NA_proportion){
  data <- NA_proportion[[1]]
  
  #Colors
  col <- rev(fishualize::fish(n = 2, option = "Ostracion_whitleyi", begin = 0.2, end = 0.9))
  
  P1 <- ggplot(data) +
    geom_point(aes(x= prop_NA, y= precision, col = type)) +
    geom_smooth(aes(x= prop_NA, y= precision, col = type))+
    scale_color_manual(values = col)+
    xlab("Proportion of NAs") + ylab("Assessement quality (R-squared, or Accuracy)") +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(size = 10))
  
  data_traits <- NA_proportion[[2]]
  
  P2 <-ggplot(data_traits) +
    geom_point(aes(x= prop_NA, y= precision, col = type)) +
    geom_smooth(aes(x= prop_NA, y= precision, col = type))+
    scale_color_manual(values = col)+
    xlab("Proportion of NAs") + ylab("Assessement quality (R-squared, or Accuracy)") +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(size = 10))+
    facet_wrap(~variable)
  
  list(P1, P2)
}  #END OF FUNCTION estimates_initial_NA
  




estimates_according_NA_sp <- function(order_performance = order_performance,
                                   original_data = species_traits_final){
  # Number of initial NAs per order
  order <- unique(original_data$order)[order(unique(original_data$order))]
  ratio_NA <- list()
  
  for(i in order){
    temp <- dplyr::filter(original_data, order == i) |> 
      dplyr::select(-c(rls_species_name:fishbase_name))
    prop_NA <- mean(is.na(temp))
    nb_species <- nrow(temp)
    
    ratio_NA[[i]] <- c(prop_NA, nb_species)
  }
  
  NA_order <- do.call(rbind, ratio_NA) 
  colnames(NA_order) <- c("NA_proportion", "species_nb")
  NA_order <- tibble::rownames_to_column(as.data.frame(NA_order), var="order")

  
  # mean performance per order
  perf_order <- list()
  for(i in order){
    temp <- dplyr::filter(order_performance, order == i) |> 
      dplyr::group_by(method) |> 
      dplyr::summarise(mean = mean(estimate, na.rm=T),
                       sd = sd(estimate, na.rm =T)) |> 
      # tidyr::pivot_wider(names_from = "method", values_from = c("mean", "sd")) |> 
      dplyr::mutate(order = i)
    
    perf_order[[i]] <- temp
  }
  
  res <- do.call(rbind, perf_order) |> 
    dplyr::left_join(NA_order)
  
  
  # Plot
  P1 <- ggplot(res)+
    geom_point(aes(x=NA_proportion, y = mean, col = method))+
    # geom_smooth(aes(x=NA_proportion, y = mean, col = method))+
    xlab("Na's proportion per order") + ylab(paste("Mean per order (on",
                                                   max(order_performance$model),
                                                   "models)")) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10))
  
  
  P2 <- ggplot(res)+
    geom_point(aes(x=log10(species_nb), y = mean, col = method))+
    xlab("log10(Number of species per order)") + ylab(paste("Mean per order (on",
                                                            max(order_performance$model),
                                                            "models)")) +
    theme_minimal() +
    theme(legend.position = "right",
          legend.background = element_rect(color ="black"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10))
  
  P1+P2
  } # END of estimates_according_NA_sp (6)






##-------------4) Check NA on map-------------

NA_on_map <- function(data=covariates,
                      variable = "coral_algae_10km",
                      xlim = c(-180,180),
                      ylim = c(-60,60),
                      jitter = 1,
                      lat_line = 30,
                      priority= "no"){
  
  worldmap <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')
  not_priority <- ifelse(priority=="no", "yes", "no") #set the order of the points in ggplot
  
  df <- data |> 
    dplyr::select(longitude, latitude, all_of(variable)) |> 
    dplyr::mutate(presence_of_data = ifelse(!is.na(data[,variable]), "yes", "no"))
  
  ggplot() +
    geom_sf(data = worldmap, color = NA, fill = "grey70") +
    geom_jitter(data=df[df$presence_of_data == not_priority,], 
                aes(x=longitude, y=latitude, color= presence_of_data),
                width = jitter, height = jitter)+
    geom_jitter(data=df[df$presence_of_data == priority,], 
                aes(x=longitude, y=latitude, color= presence_of_data),
                width = jitter, height = jitter)+
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "PuOr")[c(1,3)]) +  # Change colors here
    
    coord_sf(ylim= ylim, xlim =xlim , expand = FALSE) +
    
    geom_hline(aes(yintercept=c(-lat_line,lat_line)), linetype = "dashed", linewidth=0.3)+
    theme_bw()+
    labs(x="Longitude", y= "Latitude", title = variable) +
    theme(axis.title.x = element_text(face = "bold",
                                      size = 15),
          axis.title.y = element_text(face = "bold",
                                      size = 15),
          axis.text = element_text(size=13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          plot.title = element_text(size=10, face="bold"),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
  
}


##-------------5) Density plot predictive models -------------

density_prediction <- function(raw_result = all_res){
  
  #rescale the data between 0 and 1
  res_scaled <- raw_result |> 
  dplyr::group_by(variable) |> 
  dplyr::mutate(min_obs = min(observed), max_obs = max(observed),
                min_imp = min(imputed), max_imp = max(imputed)) |> 
  dplyr::mutate(observed = (observed - min_obs)/(max_obs - min_obs),
                imputed = (imputed - min_imp)/(max_imp- min_imp))
  
  #density plot
  ggplot(res_scaled, aes(x = observed, y = imputed)) +
    geom_density_2d_filled(aes(x = observed, y = imputed),
                           contour_var = 'ndensity',
                           contour = F, n = 100, bins= 10, colour = 'transparent') +
    scale_fill_viridis_d(option = 'viridis', begin = 0.2, end = 0.9,
                         name = 'Count') +
    # geom_point(aes(x = observed, y = predicted), alpha = 0.2)
    theme_bw()+ 
    theme(panel.grid = element_blank(), 
          strip.background = element_rect(fill = 'grey90', colour = 'grey90'), 
          aspect.ratio = 1) + 
    facet_wrap(~variable, scales = "free") +
    geom_abline() +
    labs(x = "Observed", y = "Predicted") +
    # geom_text(label = paste0("r² = ", round(r2,3), "\nestimate = ", round(a,2)),
    #           x= 0.3, y=0.9, size = 5)+
    theme(
      axis.text=element_text(size=10),
      axis.title=element_text(size=15),
      legend.text=element_text(size=8), 
      legend.title=element_text(size=8),
      strip.text.x = element_text(size = 10),
      strip.text.y = element_text(size = 10),
      strip.background = element_blank(),
      panel.background = element_rect(fill = "white", colour = "grey50",
                                      linewidth = 1, linetype = "solid"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank())

} # END of density_prediction




##-------------6) plot distributions of variables -------------

distribution_plot <- function(dataframe = covariates_final,
                              longer = TRUE,
                              cols_plot = c("depth", "gdp"),
                              cols_not_plot = NULL,
                              index_values = NULL,  #if longer = F, give the names of the 'index' and 'values' columns in a vector
                              xlabel = NULL,
                              ylabel = NULL,
                              axis_title_size = NULL,
                              strip_txt_size = NULL,
                              strip_txt_face = NULL
                              )
  { 
  
  library(ggplot2)
  
  if(longer){
    
    if(!is.null(cols_not_plot)){
      
      data <- tidyr::pivot_longer(dataframe,
                               cols = -all_of(cols_not_plot),
                               names_to = "index",
                               values_to = "values")
    }else{
      
      data <- tidyr::pivot_longer(dataframe,
                                  cols = all_of(cols_plot),
                                  names_to = "index",
                                  values_to = "values")
    }
    
  }else{
    data <- dataframe |> 
      dplyr::rename( index = index_values[1],
                     values = index_values[2])
    
  }
  
  
  #plot distributions
  ggplot(data)+
    aes(x=values, group=index, fill=index) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, color = "grey40", fill ="white") +
    geom_density(aes(fill = index), alpha = 0.2) +
    hrbrthemes::theme_ipsum(axis_title_size = axis_title_size,
                            strip_text_size = strip_txt_size,
                            strip_text_face = strip_txt_size) +
    xlab(xlabel) + ylab(ylabel)+
    facet_wrap(~index, scales = "free", ncol =4) +
    theme(legend.position="none",panel.spacing = unit(0.1, "lines"),
          axis.ticks.x=element_blank())
  
} #END OF FUNCTION 'distribution_plot'



##-------------7) plot relationship of X and Y -------------

plot_interaction <- function(dataframe = residuals,
                              var_facet_wrap = NULL, #give the name of the column giving the different elements of the facet wrap
                              X_values = NULL,
                              Y_values = NULL,
                              xlabel = NULL,
                              ylabel = NULL,
                              axis_title_size = NULL,
                              strip_txt_size = NULL,
                              strip_txt_face = NULL,
                             add_curve = NULL
)
{ 
  
  library(ggplot2)
  
    data <- dataframe |> 
      dplyr::rename( index = var_facet_wrap,
                     X = X_values,
                     Y = Y_values)
    
  
  #plot relationship
  plot <- ggplot(data)+
    geom_point(aes(x = X, y = Y, fill = index),
               color = "grey40", alpha = 0.2, shape = 21) +
    hrbrthemes::theme_ipsum(
      axis_title_size = axis_title_size,
                            strip_text_size = strip_txt_size,
                            strip_text_face = strip_txt_size
      ) +
    xlab(xlabel) + ylab(ylabel)+
    facet_wrap(~index, scales = "free") +
    theme(legend.position="none",panel.spacing = unit(0.1, "lines"),
          axis.ticks.x=element_blank())
  
  if(!is.null(add_curve)){
    plot <- plot +
      geom_smooth(aes(X,Y),method = "gam", se = FALSE,color = "grey50")
  }
  
  plot
} #END OF FUNCTION 'interaction_plot'



##-------------8) plot Contributions on map-------------

plot_Contrib_on_world_map <-function(data=contributions_surveys_log[order(contributions_surveys_log[,NCP]),],
                                     NCP = "Taxonomic_Richness",
                                     xlim=c(-180,180), ylim = c(-36, 31),
                                     title="world map with ",
                                     jitter=1.5, 
                                     pt_size=2,
                                     save=T){
  
  library(ggplot2)
  
  #Coastline
  coast <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')
  
  map <- ggplot(data) +
    geom_sf(data = coast, color = "grey30", fill = "lightgrey",
            size=0.1) +
    
    geom_point(data=data,
               size = pt_size, shape = 20,
               position=position_jitter(width=jitter, height = jitter),
               aes(x = longitude, y = latitude,
                   colour= data[,NCP],
                   alpha = 0.7)) +
    scale_colour_gradientn(name  = NCP,
                           colours = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu")))+
    
    
    coord_sf(xlim, ylim, expand = FALSE) +
    guides(alpha= "none" ) +
    # scale_size_continuous(range = c(0.5, 4), guide = "none") +
    theme_minimal()+
    labs(#title = paste0(NCP, " geographic distribution"),
      x="", y= "") +
    theme(legend.position = "bottom",
          plot.title = element_text(size=10, face="bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
  
  if(save){
    ggsave( here::here("figures", "2_map_contributions", paste0( title , NCP, ".jpg")), plot = map, width=15, height = 7 )
  }else{map}
}
