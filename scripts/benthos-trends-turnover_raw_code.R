## ----setup, include=FALSE------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- load_packages, message = FALSE-------------------------------------------------------
# basic data manipulation and visualisation
library(tidyverse)
library(here)
library(janitor)
library(viridis)
library(biscale)
# spatial data processing and mapping
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(RNetCDF)
# diversity and turnover analysis
library(vegan)
library(BAT)


## ---- read_data, message = FALSE-----------------------------------------------------------
sample_events <- read_csv(here("data",
                               "derived_data/sample_events.csv"))
pres_df <- read_csv(here("data",
                         "derived_data/pres_df.csv"))


## ---- exclude_nas--------------------------------------------------------------------------
sample_events_includingNA <- sample_events
sample_events <- sample_events %>% filter(includes_na == FALSE)


## ---- filter_samples_to_species_data-------------------------------------------------------
sample_events <- sample_events %>% filter(sample_id %in% unique(pres_df$sample_id))


## ---- five-yr-blocs------------------------------------------------------------------------
sample_events <- sample_events %>% 
  mutate(year_5 = round(year/5)*5)
ggplot(sample_events) + geom_bar(aes(x = year_5))



## ---- samples_x_time_x_space---------------------------------------------------------------
ggplot(st_as_sf(sample_events, coords = c("lon", "lat"), crs = 4326)) +
  geom_sf(size = 0.1) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 6)) +
  facet_wrap(~ year_5)


## ---- define_time_periods------------------------------------------------------------------
sample_events <- sample_events %>%
  mutate(time_slice = case_when(
    year < 1990 ~ 1,
    year >= 1990 & year <2000 ~ 2,
    year >= 2000 & year <2005 ~ 3,
    year >= 2005 & year <2010 ~ 4,
    year >= 2010 & year <2015 ~ 5,
    year >= 2015 ~ 6
  ))

sample_events %>% count(time_slice)


## ---- show_spatial_extent------------------------------------------------------------------
sample_events %>% select(lon, lat) %>% summary()


## ---- create_grid--------------------------------------------------------------------------
extent_tb <- tibble(lon = c(-34, 59), lat = c(28, 82)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
r <- rast(ext(extent_tb), resolution = 1, crs = crs(extent_tb))


## ---- create_samples_by_time_rast----------------------------------------------------------
# create 6-layer raster
samp_event_by_time_r <- rep(r, 6)

# vector of names for time slices
time_slices <- c("pre1990", "1990s", "2000-2004", "2005-2009", "2010-2014", "post2015")

# fill the raster from sample_events
for(i in 1:nlyr(samp_event_by_time_r)){
  samp_event_by_time_r[[i]] <- sample_events %>%
    filter(time_slice == i) %>% 
    dplyr::select(lon, lat) %>% 
    as.matrix() %>% 
    rasterize(r, fun = "length") %>%
    setNames(time_slices[i])
}



## ---- get_world----------------------------------------------------------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_transform(crs = crs(r))


## ---- map_samples_through_time-------------------------------------------------------------
ggplot() +
  geom_spatraster(data = log10(samp_event_by_time_r)) +
  scale_fill_viridis_c(name = "log N Samples") +
  geom_sf(data = world, colour = "grey95", fill = "grey85", lwd = 0.1, alpha = 0.5) +
  xlim(as.vector(ext(r)[1:2])) +
  ylim(as.vector(ext(r)[3:4])) +
  coord_sf(expand = FALSE) +
  facet_wrap(~lyr)



## ---- add_cellid_sample_events-------------------------------------------------------------
sample_events <- sample_events %>% 
  mutate(sample_cell = cellFromXY(r,
                                  as.matrix(dplyr::select(sample_events, lon, lat))))



## ---- get_samples_per_cell_per_time, message = FALSE---------------------------------------
sample_events %>%
  dplyr::select(sample_cell, time_slice) %>%
  distinct() %>%
  count(sample_cell) %>%
  count(n) %>% 
  rename(n_time_periods = n, n_grid_cells = nn) %>% 
  mutate(p_cells = round(n_grid_cells / sum(n_grid_cells), 2))


## ------------------------------------------------------------------------------------------
cells_single_time <- 
  sample_events %>% select(sample_cell, time_slice) %>% 
  distinct() %>% 
  count(sample_cell) %>% filter(n == 1) %>%
  pull(sample_cell)

cells_multi_time <- 
  sample_events %>% select(sample_cell, time_slice) %>% 
  distinct() %>%  
  count(sample_cell) %>% filter(n > 1) %>%
  pull(sample_cell)

single_multi_time <- c(r,r) %>% setNames(c("single_t", "multi_t"))

values(single_multi_time[[1]]) <- NA
values(single_multi_time[[1]])[cells_single_time] <- 1

values(single_multi_time[[2]]) <- NA
values(single_multi_time[[2]])[cells_multi_time] <- 1

(single_v_multi_t_cell_map <- ggplot() +
  geom_spatraster(data = single_multi_time) +
  scale_fill_viridis_c(option = "turbo") +
  geom_sf(data = world, colour = "grey95", fill = "grey85", lwd = 0.1, alpha = 0.5) +
  xlim(as.vector(ext(r)[1:2])) +
  ylim(as.vector(ext(r)[3:4])) +
  coord_sf(expand = FALSE) +
  theme(legend.position = "none") +
  facet_wrap(~lyr)
)


## ---- create_gridded_occs, message = FALSE-------------------------------------------------
gridded_occs <- pres_df %>%
  left_join(sample_events, join_by(sample_id)) %>%
  group_by(aphia_id, time_slice, sample_cell) %>% 
  summarise(sp_occs = n()) %>% 
  arrange(sample_cell, time_slice) %>% ungroup()



## ---- filter_occs_multi_time---------------------------------------------------------------
gridded_occs <- gridded_occs %>%
  filter(sample_cell %in% cells_multi_time)


## ---- get_occs_by_time_1sp-----------------------------------------------------------------
gridded_occs %>%
  filter(aphia_id == 101160) %>% arrange(sample_cell, time_slice)



## ---- message = FALSE----------------------------------------------------------------------
gridded_occs %>%
  filter(aphia_id == 101160) %>%
  count(sample_cell) %>% 
  count(n) %>% 
  rename(n_time_periods = n, n_grid_cells = nn)


## ---- get_example_grid_cell----------------------------------------------------------------
gridded_occs %>% dplyr::select(time_slice, sample_cell) %>%
  distinct() %>%
  count(sample_cell) %>%
  arrange(desc(n)) %>% 
  head()



## ---- example_community_mat----------------------------------------------------------------
(comm_mat_eg <- gridded_occs %>% filter(sample_cell == 1544) %>% 
  mutate(pa = ifelse(sp_occs > 0, 1, 0)) %>% 
  dplyr::select(time_slice, aphia_id, pa) %>% 
  pivot_wider(names_from = aphia_id, values_from = pa, values_fill = 0) %>% 
  column_to_rownames(var = "time_slice")
)


## ---- get_beta_example---------------------------------------------------------------------
b_eg <- comm_mat_eg %>% 
  as.matrix() %>% beta()
abc_eg <- comm_mat_eg %>% betadiver(method = NULL)



## ---- get_species_turnover-----------------------------------------------------------------
source(here("scripts", "get_species_turnover.R"))


## ------------------------------------------------------------------------------------------
gridded_occs %>% filter(sample_cell == 1544) %>% get_species_turnover()



## ---- get_species_turnover_all_cells, message = FALSE--------------------------------------
species_turnover_allcells <- gridded_occs %>%
  split(.$sample_cell) %>%
  purrr::map(get_species_turnover, .progress = TRUE) %>% 
  bind_rows()


## ---- tidy_species_turnover_all_cells------------------------------------------------------
(species_turnover_allcells  <- species_turnover_allcells  %>% 
  mutate(
    t_start_long = case_match(
      time_start,
      1 ~ "pre1990",
      2 ~ "1990s",
      3 ~ "2000-2004",
      4 ~ "2005-2009",
      5 ~ "2012-2014",
      6 ~ "post2015"
    ),
    t_comp_long = case_match(
      time_comp,
      1 ~ "pre1990",
      2 ~ "1990s",
      3 ~ "2000-2004",
      4 ~ "2005-2009",
      5 ~ "2012-2014",
      6 ~ "post2015"),
    t_comp = paste(time_start, time_comp, sep = "_"),
    t_comp_lab = paste(t_start_long, t_comp_long, sep = "_v_"),
    t_steps = time_comp - time_start,
    b_prop = b / (a + b + c), c_prop = c / (a + b + c)
  ) %>% 
  select(sample_cell, t_comp, t_steps, beta_tot:c, b_prop, c_prop, t_comp_lab)
)


## ---- add_sample_meta_to_turnover, message = FALSE-----------------------------------------
# Get sample cell metadata
sample_cell_meta <- sample_events %>%
  group_by(sample_cell) %>%
  summarise(n_samps = n_distinct(sample_id)) %>% 
  bind_cols(xyFromCell(r, .$sample_cell)) %>% 
  rename(lon = x, lat = y)

# Then join to species turnover data
species_turnover_allcells <- species_turnover_allcells %>% 
  left_join(sample_cell_meta, join_by(sample_cell)) %>% 
  select(lon, lat, everything())



## ---- create_beta_tot_raster---------------------------------------------------------------
beta_tot_r <- species_turnover_allcells %>% 
  select(lon, lat, t_comp, t_comp_lab, beta_tot) %>% 
  arrange(t_comp) %>% 
  select(-t_comp) %>% 
  pivot_wider(names_from = t_comp_lab,
              values_from = beta_tot) %>% 
  rast(type = "xyz", crs = crs(r)) %>% 
  project(r)



## ---- plot_beta_tot_raster-----------------------------------------------------------------

(beta_tot_map <- ggplot() +
  geom_spatraster(data = beta_tot_r) +
  scale_fill_viridis_c(name = "beta (total)") +
  geom_sf(data = world, colour = "grey95", fill = "grey85", lwd = 0.1, alpha = 0.5) +
  xlim(as.vector(ext(r)[1:2])) +
  ylim(as.vector(ext(r)[3:4])) +
  coord_sf(expand = FALSE) +
  facet_wrap(~lyr) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
)


## ---- plot_beta_tot_raster_layer-----------------------------------------------------------
ggplot() +
  geom_spatraster(data = select(beta_tot_r, 1, 5)) +
  scale_fill_viridis_c(name = "beta (total)") +
  geom_sf(data = world, colour = "grey95", fill = "grey85", lwd = 0.1, alpha = 0.5) +
  xlim(as.vector(ext(r)[1:2])) +
  ylim(as.vector(ext(r)[3:4])) +
  coord_sf(expand = FALSE) +
  facet_wrap(~lyr) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


## ---- create_turnover_raster---------------------------------------------------------------
source(here("scripts", "create_turnover_raster.R"))


## ---- get_turnover_raster_cprop------------------------------------------------------------
turnover_r_cprop <- create_turnover_raster(turnover_metric = "c_prop",
                       create_plot = TRUE,
                       display_plot = FALSE,
                       plot_layers = c("2000-2004_v_2005-2009",
                                       "2000-2004_v_2012-2014",
                                       "2000-2004_v_post2015")
                       )


## ---- display_cprop_raster-----------------------------------------------------------------
turnover_r_cprop$turnover_rast


## ---- display_cprop_plot-------------------------------------------------------------------
turnover_r_cprop$turnover_plot


## ---- create_all_turnover_rasters----------------------------------------------------------
beta_tot_r <- create_turnover_raster(turnover_metric = "beta_tot", create_plot = FALSE)
beta_repl_r <- create_turnover_raster(turnover_metric = "beta_repl", create_plot = FALSE)
beta_rich_r <- create_turnover_raster(turnover_metric = "beta_rich", create_plot = FALSE)

sp_shared_r <- create_turnover_raster(turnover_metric = "a", create_plot = FALSE)
sp_lost_r <- create_turnover_raster(turnover_metric = "b", create_plot = FALSE)
sp_gained_r <- create_turnover_raster(turnover_metric = "c", create_plot = FALSE)

p_sp_lost_r <- create_turnover_raster(turnover_metric = "b_prop", create_plot = FALSE)
p_sp_gained_r <- create_turnover_raster(turnover_metric = "c_prop", create_plot = FALSE)



## ---- write_turnover_tifs------------------------------------------------------------------
writeRaster(beta_tot_r, filename = here("product", "beta_tot_r.tif"), overwrite = TRUE)
writeRaster(beta_repl_r, filename = here("product", "beta_repl_r.tif"), overwrite = TRUE)
writeRaster(beta_rich_r, filename = here("product", "beta_rich_r.tif"), overwrite = TRUE)
writeRaster(sp_shared_r, filename = here("product", "sp_shared_r.tif"), overwrite = TRUE)
writeRaster(sp_lost_r, filename = here("product", "sp_lost_r.tif"), overwrite = TRUE)
writeRaster(sp_gained_r, filename = here("product", "sp_gained_r.tif"), overwrite = TRUE)
writeRaster(p_sp_lost_r, filename = here("product", "p_sp_lost_r.tif"), overwrite = TRUE)
writeRaster(p_sp_gained_r, filename = here("product", "p_sp_gained_r.tif"), overwrite = TRUE)


## ---- get_netcdfify------------------------------------------------------------------------
source(here("scripts", "netcdfify.R"))


## ---- read_diveristy_rasters---------------------------------------------------------------
beta_total <- rast(here("product", "beta_tot_r.tif"))
beta_repl <- rast(here("product", "beta_repl_r.tif"))
beta_rich <- rast(here("product", "beta_rich_r.tif"))
sp_shared <- rast(here("product", "sp_shared_r.tif"))
sp_lost <- rast(here("product", "sp_lost_r.tif"))
sp_gained <- rast(here("product", "sp_gained_r.tif"))
p_sp_lost <- rast(here("product", "p_sp_lost_r.tif"))
p_sp_gained <- rast(here("product", "p_sp_gained_r.tif"))



## ---- add_layername_prefix-----------------------------------------------------------------
names(beta_total) <- paste0("beta_total_", names(beta_total))
names(beta_repl) <- paste0("beta_repl_", names(beta_repl))
names(beta_rich) <- paste0("beta_rich_", names(beta_rich))
names(sp_shared) <- paste0("sp_shared_", names(sp_shared))
names(sp_lost) <- paste0("sp_lost_", names(sp_lost))
names(sp_gained) <- paste0("sp_gained_", names(sp_gained))
names(p_sp_lost) <- paste0("p_sp_lost_", names(p_sp_lost))
names(p_sp_gained) <- paste0("p_sp_gained_", names(p_sp_gained))



## ---- combine_rasters----------------------------------------------------------------------
all_measures_r <- c(
  beta_total, beta_repl, beta_rich,
  sp_shared, sp_lost, sp_gained,
  p_sp_lost, p_sp_gained
)



## ---- write_all_measures_tiff--------------------------------------------------------------
writeRaster(all_measures_r,
            filename = here("product", "all_diversity_measures.tif"),
            overwrite = TRUE)



## ------------------------------------------------------------------------------------------
all_measures_minmax <- all_measures_r %>% 
  minmax() %>% 
  as_tibble() %>% 
  mutate(value = c("min", "max")) %>% 
  select(value, everything()) %>% 
  pivot_longer(cols = -value, names_to = "diversity_measure", values_to = "val") %>% 
  pivot_wider(names_from = value, values_from = val)

write_csv(all_measures_minmax,
          here::here("data", "derived_data/all_measures_minmax.csv"))



## ---- set_global_attrs---------------------------------------------------------------------
global_attr <- list(
  title = "Diversity_measures",
  summary = "Eight measures of changes in species compositions of European macrobenthic communities between different time periods based on presence / absence in  on a 1 degree grid",                       
  Conventions = "CF-1.8",
  # id = "",
  naming_authority = "emodnet-biology.eu",
  history = "https://github.com/EMODnet/benthos-trends",
  source = "https://github.com/EMODnet/benthos-trends",
  # processing_level = "",
  # comment = "", 
  # acknowledgment = "",
  license = "CC-BY",
  standard_name_vocabulary = "CF Standard Name Table v1.8",
  date_created = as.character(Sys.Date()),
  creator_name = "Tom Webb",
  creator_email = "t.j.webb@sheffield.ac.uk",
  creator_url = "https://www.sheffield.ac.uk/biosciences/people/academic-staff/tom-webb",
  institution = "The University of Sheffield",
  project = "EMODnet-Biology",
  publisher_name = "EMODnet-Biology",                 
  publisher_email = "bio@emodnet.eu",                
  publisher_url = "www.emodnet-biology.eu",                  
  # geospatial_bounds = "",              
  # geospatial_bounds_crs = "",          
  # geospatial_bounds_vertical_crs = "", 
  geospatial_lat_min = ext(all_measures_r)[3],
  geospatial_lat_max = ext(all_measures_r)[4],
  geospatial_lon_min = ext(all_measures_r)[1],
  geospatial_lon_max = ext(all_measures_r)[2],
  # geospatial_vertical_min = "",        
  # geospatial_vertical_max = "",        
  # geospatial_vertical_positive = "",  
  # time_coverage_start = "1911",            
  # time_coverage_end = "2016",              
  # time_coverage_duration = "",         
  # time_coverage_resolution = "",       
  # uuid = "",                           
  # sea_name = "",                       
  # creator_type = "",                   
  creator_institution = "The University of Sheffield",            
  # publisher_type = "",                 
  publisher_institution = "Flanders Marine Institute (VLIZ)",        
  # program = "",                        
  # contributor_name = "",               
  # contributor_role  = "",              
  geospatial_lat_units = "degrees_north",           
  geospatial_lon_units = "degrees_east",           
  # geospatial_vertical_units   = "",    
  # date_modified = "",               
  # date_issued = "",                    
  # date_metadata_modified   = "",       
  # product_version = "",            
  # keywords_vocabulary = "",          
  # platform  = "",              
  # platform_vocabulary = "",          
  # instrument = "",          
  # instrument_vocabulary  = "",        
  # featureType = "Point",                  
  # metadata_link = "",                  
  # references = "",
  comment = "Uses attributes recommended by http://cfconventions.org",
  license = "CC-BY", 
  publisher_name = "EMODnet Biology Data Management Team",
  citation = "Webb, T.J. (2023). Temporal turnover of macrobenthos in European seas.",
  acknowledgement = "European Marine Observation Data Network (EMODnet) Biology project (EMFF/2019/1.3.1.9/Lot 6/SI2.837974), funded by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund"
)



## ---- netcdfify_all_measures---------------------------------------------------------------
netcdfify(focal_rast = all_measures_r,
          output_fname = "all_diversity_measures",
          add_measure_name = FALSE,
          global_atts = global_attr)


## ---- netcdfify_single_measure-------------------------------------------------------------
global_attr$title <- "Beta diversity"
global_attr$summary <- "Total beta diversity between different time periods based on presence / absence in European macrobenthic communities on a 1 degree grid"
beta_total <- rast(here("product", "beta_tot_r.tif"))
netcdfify(focal_rast = beta_total,
          output_fname = "beta_total",
          add_measure_name = TRUE,
          global_atts = global_attr)


## ---- reproducibility_date, echo = FALSE---------------------------------------------------
lubridate::now()


## ---- session_info, echo = FALSE-----------------------------------------------------------
sessionInfo()

