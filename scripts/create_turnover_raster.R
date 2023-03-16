create_turnover_raster <- function(turnover_dat = species_turnover_allcells,
                                   turnover_metric = "beta_tot",
                                   create_plot = TRUE,
                                   display_plot = FALSE,
                                   plot_layers = NULL,
                                   basemap = world,
                                   template_raster = r){
  
  turnover_r <- turnover_dat %>% 
    select(lon, lat, t_comp, t_comp_lab, all_of(turnover_metric)) %>% 
    arrange(t_comp) %>% 
    select(-t_comp) %>% 
    pivot_wider(names_from = t_comp_lab,
                values_from = all_of(turnover_metric)) %>% 
    rast(type = "xyz", crs = crs(template_raster)) %>% 
    project(template_raster)
  
  # create the plot if required
  if(create_plot == TRUE){
    if(is_null(plot_layers)){
      r_plot <- turnover_r
    } else {
      r_plot <- select(turnover_r, all_of(plot_layers))
    }
    turnover_plot <- ggplot() +
      geom_spatraster(data = r_plot) +
      scale_fill_viridis_c(name = turnover_metric) +
      geom_sf(data = basemap, colour = "grey95", fill = "grey85", lwd = 0.1, alpha = 0.5) +
      xlim(as.vector(ext(template_raster)[1:2])) +
      ylim(as.vector(ext(template_raster)[3:4])) +
      coord_sf(expand = FALSE) +
      facet_wrap(~lyr) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  } else {
    turnover_plot <- NULL
  }
  
  if(display_plot == TRUE & !is.null(turnover_plot)){
    print(turnover_plot)
  }
  
  if(is.null(turnover_plot)){
    ret <- turnover_r
  } else {
    ret <- list(turnover_rast = turnover_r, turnover_plot = turnover_plot)
  }
  
  ret
  
}