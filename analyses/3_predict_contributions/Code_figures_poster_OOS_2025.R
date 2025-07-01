#Counterfactuals colors:
Cl_color = "darkseagreen3"
HF_color = "firebrick3"

Cl_color= "#86b386"
HF_color = "#951d1d"
# # install.packages("colorspace")
# library(colorspace)
# 
# # Couleurs assombries
# Cl_color <- darken(Cl_color, amount = 0.1)  # 20% plus sombre
# HF_color <- darken(HF_color, amount = 0.1)


##Plot
lolliplot <- ggplot(change_percent)+
  aes(x = median_t, y = contribution, color = counterfactual) +
  geom_segment(aes(x = 0, xend = median_t, y = contribution, yend = contribution,
                   alpha = different_from_zero), 
               size = 1) + 
  geom_point(aes(alpha = different_from_zero, size = bigger_effect)) + 
  scale_alpha_manual(values = c("yes" = 1, "no" = 0.3)) +
  scale_size_manual(values = c("yes" = 5, "no" = 2)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
  
  #broken axis
  geom_rect(aes(xmin = trans_x(x_break-16), xmax = x_break-4, ymin = -Inf, ymax = Inf), 
            fill = "white", color = NA) +  # Mask broken area
  annotate("text", x = x_break-3, y = max(as.numeric(factor(change_percent$contribution))), 
           label = "//", size = 7) + # Add a visual cue for the break
  
  scale_x_continuous(
    breaks = trans_x(xticks),
    labels = paste0(xticks, "%") )+
  scale_color_manual(values = c("conservation_legacy" = Cl_color, 
                                "human_footprint" = HF_color)) + 
  labs(
    x = "Contribution changes in counterfactual scenarios", 
    y = "", #"Contributions",
    color = "Counterfactuals") +
  theme_bw() +
  
  # Custom labels with colored points
  scale_y_discrete(labels = function(labels) {
    sapply(seq_along(labels), function(i) {
      group_color <- unique(change_percent$group[change_percent$contribution == labels[i]])
      color <- ifelse(group_color == "NN", "forestgreen",
                      ifelse(group_color == "NC", "darkgoldenrod2", "dodgerblue3")) # Adjust colors
      paste0("<span style='color:", color, "; font-size: 13px;'>&#9679;</span> ", labels[i])
    })
  }) +
  
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text = ggtext::element_markdown(size = 12, color = "black"), #axis.text.y = element_text(size = 10), #, color = change_percent$group),
        axis.title = element_text(size = 13),
        legend.position = "bottom",
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0,1,0,0, unit = "cm"),
        legend.title = element_text(face="bold", size = 11),
        panel.grid.major = element_line(color = "grey20", size = 0.1),
        panel.grid.minor = element_line(color = "grey20", size = 0.1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE),
         alpha = "none",#guide_legend(nrow = 2, byrow = TRUE),
         size = "none")
lolliplot


lolliplot_prop_column <- lolliplot +
  theme( #axis.text.x = ggtext::element_markdown(size = 10, angle = 45,  color = "black")
    axis.title.x = element_text(vjust = -5, hjust = 0.2))+
  guides(fill = "none")
lolliplot_prop_column

plot_merged_prop <- lolliplot_prop_column +
  theme(plot.margin =unit(c(0.1,0.1,1.3,0.1), 'cm'))+
  theme(legend.position = "none")
plot_merged_prop

ggsave(filename =  paste0(path_file,"/Lolliplot_OOS_poster.png"),
       plot = plot_merged_prop,  width = 9, height = 8,
       bg = "transparent")






## Without grid

##Plot
lolliplot <- ggplot(change_percent)+
  aes(x = median_t, y = contribution, color = counterfactual) +
  geom_segment(aes(x = 0, xend = median_t, y = contribution, yend = contribution,
                   alpha = different_from_zero), 
               size = 1) + 
  geom_point(aes(alpha = different_from_zero, size = bigger_effect)) + 
  scale_alpha_manual(values = c("yes" = 1, "no" = 0.3)) +
  scale_size_manual(values = c("yes" = 5, "no" = 2)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
  
  #broken axis
  geom_rect(aes(xmin = trans_x(x_break-16), xmax = x_break-4, ymin = -Inf, ymax = Inf), 
            fill = "white", color = NA) +  # Mask broken area
  annotate("text", x = x_break-3, y = max(as.numeric(factor(change_percent$contribution))), 
           label = "//", size = 7) + # Add a visual cue for the break
  
  scale_x_continuous(
    breaks = trans_x(xticks),
    labels = paste0(xticks, "%") )+
  scale_color_manual(values = c("conservation_legacy" = Cl_color, 
                                "human_footprint" = HF_color)) + 
  labs(
    x = "Contribution changes in counterfactual scenarios", 
    y = "", #"Contributions",
    color = "Counterfactuals") +
  theme_bw() +
  
  # Custom labels with colored points
  scale_y_discrete(labels = function(labels) {
    sapply(seq_along(labels), function(i) {
      group_color <- unique(change_percent$group[change_percent$contribution == labels[i]])
      color <- ifelse(group_color == "NN", "forestgreen",
                      ifelse(group_color == "NC", "darkgoldenrod2", "dodgerblue3")) # Adjust colors
      paste0("<span style='color:", color, "; font-size: 13px;'>&#9679;</span> ", labels[i])
    })
  }) +
  
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text = ggtext::element_markdown(size = 12, color = "black"), #axis.text.y = element_text(size = 10), #, color = change_percent$group),
        axis.title = element_text(size = 13),
        legend.position = "bottom",
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0,1,0,0, unit = "cm"),
        legend.title = element_text(face="bold", size = 11),
        panel.grid.major = element_line(color = "grey20", size = 0),
        panel.grid.minor = element_line(color = "grey20", size = 0)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE),
         alpha = "none",#guide_legend(nrow = 2, byrow = TRUE),
         size = "none")
lolliplot


lolliplot_prop_column <- lolliplot +
  theme( #axis.text.x = ggtext::element_markdown(size = 10, angle = 45,  color = "black")
    axis.title.x = element_text(vjust = -5, hjust = 0.2))+
  guides(fill = "none")
lolliplot_prop_column

plot_merged_prop <- lolliplot_prop_column +
  theme(plot.margin =unit(c(0.1,0.1,1.3,0.1), 'cm'))+
  theme(legend.position = "none")
plot_merged_prop

ggsave(filename =  paste0(path_file,"/Lolliplot_OOS_poster_without_grid.png"),
       plot = plot_merged_prop,  width = 9, height = 8,
       bg = "transparent")





### fig 2 -----
drivers_to_plot =  list(
  c("protection_statusfull","Gravity","GDP")
)

drivers = drivers_to_plot[[1]]
drivers_name <- drivers
all_drivers <- c(drivers, paste0(drivers, "_Deg1"),  paste0(drivers, "_Deg2"))

support_estimates <- postBeta[["support"]]
rownames(support_estimates) <- model_fit_mcmc[["covNames"]]
support_estimates <- as.data.frame(support_estimates) |> 
  tibble::rownames_to_column("covariate") |> 
  tidyr::pivot_longer(cols = -covariate,
                      names_to = "response", values_to = "support") |> 
  dplyr::mutate(support = dplyr::case_when(support > 0.95 ~ 1, 
                                           support < 0.05 ~ 1,
                                           TRUE ~ 0 ))|> 
  dplyr::mutate(support = factor(support, levels = c(0, 1)))


#Filter estimates table
df <- S_arranged |>
  dplyr::filter(covariate %in% all_drivers) |> 
  dplyr::left_join(support_estimates) |> 
  dplyr::left_join(dplyr::rename(grp_NN_NP, response = "contribution")) |>
  dplyr::mutate(response = gsub("_", " ", response),
                response = dplyr::case_when(
                  response == "Iucn species richness" ~ "IUCN species richness",
                  TRUE ~ response
                )) 

medians <- df  |> 
  dplyr::filter(covariate %in% c(all_drivers[1], paste0(all_drivers[1], "_Deg1"))) |> 
  dplyr::group_by(response)  |> 
  dplyr::summarise(median_value = median(value),
                   group = unique(group)) |>
  dplyr::arrange(factor(group, levels = c("NS", "NN", "NC")), median_value)

df <- df |> 
  dplyr::mutate(covariate = factor(covariate, levels = all_drivers)) |> 
  dplyr::mutate(response = factor(response, levels = medians$response))

new_titles <- c(
  "protection_statusfull" = "Full protection",
  "Fishing_vessel_density" = "Fishing pressure",
  "Gravity" = "Gravity",
  "GDP" = "GDP",
  
  "protection_statusrestricted" = "Restricted MPA",
  "SST_5_years" = "SST (5 years)",
  "Chlorophyll_5_years" = "Chlorophyll (5 years)",
  "DHW_quantile95_5_years" = "DHW (5 years)",
  "Marine_ecosystem_dependency" = "MED",
  "Travel_time" = "Travel time",
  "HDI" = "HDI",
  "Natural_ressource_rent" = "Ressource rent"
)


ridges_plot <- ggplot(df) +
  aes(y = response, x = value,  fill = group) +
  ggridges::geom_density_ridges(aes(alpha = support), linewidth = 0.3,
                                scale = 3,  # >1 pour du chevauchement
                                rel_min_height = 0  # permet aux courbes de "déborder"
  )+#alpha=0.5, bandwidth = 0.005) +
  scale_fill_manual(values = c( "darkgoldenrod2", "forestgreen", "dodgerblue3"),
                    labels = c("NN" = "Nature-for-Nature",
                               "NS" = "Nature-for-Society",
                               "NC" = "Nature-as-Culture"))+
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
  scale_alpha_manual(values = c("0" = 0.25, "1" = 0.7), guide = "none") +
  hrbrthemes::theme_ipsum( axis_title_size = 0 ) +
  theme(
    legend.position="bottom",
    legend.text = element_text(size = 12, margin = margin(r = 20)),
    panel.spacing = unit(0.3, "lines"),
    strip.text.x = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, vjust = -0.2, color = "black"),
    panel.grid.major = element_line(color = "grey20", size = 0.1),
    panel.grid.minor = element_line(color = "grey20", size = 0.1)) +
  labs(fill = "")+
  xlab(all_drivers) + ylab("Nature Contributions to People and Nature")+
  facet_wrap(~covariate, ncol = length(df$covariate),
             scales = "free_x",
             labeller = labeller(covariate = new_titles)
  )

# if(drivers == c("protection_statusfull", "Fishing_vessel_density","Gravity","GDP"))
ridges_plot
ggsave(filename = paste0(path_file,"/OOS_poster_posterior_distribution_of_estimates_", save_name,
                         paste(drivers_name, collapse = "-"), ".png"),
       plot = ridges_plot, width = 9, height = 7)



## Fig 3 -----
grp_NN_NP <- data.frame(
  contribution = c("Actinopterygian richness","Functional distinctiveness",
                   "IUCN species richness" ,"Endemism",
                   "Evolutionary distinctiveness","Functional entropy",
                   "Phylogenetic entropy","Herbivores biomass",
                   "Invertivores biomass",  "Piscivores biomass",
                   "Trophic web robustness", "Mean trophic level",
                   
                   "Public attention", "Aesthetic",
                   "Available biomass", "Selenium",
                   "Zinc",   "Omega 3", "Calcium",
                   "Iron","Vitamin A", "Available biomass turnover"),
  group = c(rep("NN", 12), 
            rep("NC", 2),
            rep("NS", 8))) |> 
  dplyr::right_join(dplyr::select(change_percent, -group))

mean_NN_NP <- grp_NN_NP |> 
  dplyr::group_by(group, counterfactual) |> 
  dplyr::summarise(mean = mean(median)) |> 
  dplyr::mutate(group = factor(group, levels = c("NC", "NS", "NN")))

get_val <- function(grp, type){
  dplyr::filter(mean_NN_NP, group == grp & counterfactual == type)$mean
}

mean_plot <- ggplot(mean_NN_NP)+
  aes(x = mean , y = group, fill = group) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.5) +
  scale_fill_manual(values = c("NS" = "dodgerblue3", 
                               "NC" = "darkgoldenrod2", 
                               "NN" = "forestgreen")) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size=22), 
    axis.text.x = element_text(size = 20),
    legend.position = "none",
    plot.margin = unit(c(0,1,0,0), "cm"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_line( size = 0),
    panel.grid.minor = element_line( size = 0)) +
  labs(x = "Mean changes (%)")+
  
  #Add labels
  geom_text(aes(x = 0, y = "NN"), label = "Nature-for-Nature", color = "white", 
            hjust = 1.65, vjust = 0.5, size = 8) + 
  geom_text(aes(x = 0, y = "NS"), label = "Nature-for-Society", color = "black", 
            hjust = 1.6, vjust = 0.5, size = 8)+
  geom_text(aes(x = 0, y = "NC"), label = "Nature-as-Culture", color = "black",
            hjust = 1.64, vjust = 0.5, size = 8)+
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.7)

mean_plot


mean_plot_with_title <-  mean_plot +
  geom_text(aes(x = min(mean_NN_NP$mean), y = 3, 
                label = "Human footprint"), 
            color = HF_color, size = 7,
            hjust = -0.4, vjust = -5.2, fontface = "bold") +
  geom_text(aes(x = max(mean_NN_NP$mean) * 0.6, y =3, 
                label = "Conservation\n legacy"),
            color = Cl_color, size = 7, 
            hjust = 0.45, vjust = -1.4, fontface = "bold") +
  # make some space for the title
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 2, r = 1, b = 0, l = 0, unit = "cm"))


ggsave(filename =  paste0(path_file,"/OOS_poster_Mean_changes_in_NFF.png"),
       plot = mean_plot_with_title,  width = 7, height = 5)



## Map ----
load(file = here::here("data", "derived_data", "3_sites_covariates_to_predict.Rdata"))


table(covariates_site_final$protection_status)
# out restricted       full 
# 1132        802        869 
length(unique(covariates_site_final$country))

plot_mpa <-function(covariates_site_final, xlim=c(-180,180), ylim = c(-36, 31),
                    legend_pos = "none", jitter = 0.2){
  ggplot(covariates_site_final) +
    geom_sf(data = coast, color = "grey30", fill = "lightgrey",
            aes(size=0.1)) +
    
    geom_point(size = 2.5, na.rm = T,
               position = position_jitter(width =jitter, height =jitter),
               alpha = 0.6,
               colour = "#03045e",
               fill="#03045e",
               stroke=0.1,
               shape = 21,
               aes(x = longitude, y = latitude)) +
    
    coord_sf(xlim, ylim , expand = FALSE) +
    guides(alpha = "none", size = "none", colour = "none") +
    # scale_shape_manual(values=c(21,24,23))+
    
    theme_bw()+
    labs(title = "",
         x="", y= "") +
    theme(legend.position = legend_pos,
          plot.title = element_text(size=15, face="bold"),
          legend.text = element_text(size=13),
          legend.title = element_text(size=15),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.major = element_line( size = 0),
          panel.grid.minor = element_line( size = 0),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
}

plot_mpa(covariates_site_final , xlim=c(-180,180), ylim = c(-45, 45),
                 legend_pos = "none", jitter = 1.5)
ggsave(plot = last_plot(), width = 14, height = 7,
       filename = here::here("figures","OOS_poster_2025_RLS sites.png"))
