
# get the data from the netcdf version for the full European database
library(tidyverse)
library(ncdf4)
library(here)
library(worrms)

# Open the connection to the nc file
nc_fil <- nc_open(here("data", "EMODnet-Biology-Benthos-European-Seas/Macrobenthos_Eur_Seas_Pres_Abs_v0-4.nc"))

presabs_crs <- ncatt_get(nc_fil, "crs")



# get sampling events into a data frame

# assmeble sampling events
sample_events <- tibble(
  lat = ncvar_get(nc_fil, "lat"),
  lon = ncvar_get(nc_fil,"lon"),
  date = lubridate::ymd(as.Date(ncvar_get(nc_fil, "Date"), origin = "1970-01-01 00:00:00")),
  year = lubridate::year(as.Date(ncvar_get(nc_fil, "Date"), origin = "1970-01-01 00:00:00"))
)

# NB:
distinct(sample_events)
# is slightly shorter than sample_events
# Also:
sample_events %>% filter(date > lubridate::now())


# get taxon names and aphia ids
taxon_name <- ncvar_get(nc_fil, "Taxon_Name")
aphia_id <- ncvar_get(nc_fil,"AphiaID")
# get aphia_id from the full lsid:
# NB parse_number(aphia_id[1]) does not work, so:
aphia_id <- as.integer(word(aphia_id, start = -1, sep = ":"))

# DECISION: use species rank only. For this we need WoRMS classification:

get_worms_rank <- function(aphiaid){
  classif <- wm_classification(aphiaid) %>% slice(n())
  classif
}
# one species
get_worms_rank(aphia_id[1000])
# over species
aphia_id[1:10] %>% map_df(get_worms_rank, .progress = TRUE)

# all species:
taxo <- aphia_id %>% map_df(get_worms_rank, .progress = TRUE)

# this takes ~40 mins - so write to file:
write_csv(taxo, file = here("data", "EMODnet-Biology-Benthos-European-Seas/benthos_taxo_rank.csv"))

taxo <- read_csv(here("data", "EMODnet-Biology-Benthos-European-Seas/benthos_taxo_rank.csv"))

taxo %>% count(rank) %>% print(n = 100)

# Subset to species only
spp <- taxo %>% filter(rank == "Species")
# need to get a vector indicating their position in the bigger aphiaid vector
spp <- spp %>% mutate(aphia_index = match(.$AphiaID, taxo$AphiaID))


# Get Individual species occurrence by Event
# get presence-absence for a single species
presabs <- ncvar_get(nc_fil,"Pres_abs", start = c(1, spp$aphia_index[1]), count = c(-1,1))

# in a function:
get_presabs_sp <- function(sp_index, sp_df = spp, samps = sample_events){
  
  # get pres-abs data
  presabs <- ncvar_get(nc_fil, "Pres_abs",
                       start = c(1, sp_index), count = c(-1,1))
  # assemble data frame
  ret <- tibble(aphia_id = spp$AphiaID[spp$aphia_index == sp_index],
                presabs = presabs)
  ret <- samps %>% bind_cols(ret)
  
  ret
}

test_pa <- get_presabs_sp(sp_index = spp$aphia_index[1])
test_pa %>% count(presabs)
# NB: NA = not looked for, 0 = looked for but absent
# Maybe exclude events with NAs?
# How to do this - flag them in samp events as NA for a species?
# Use mapping to assemble big long dataset (might be really long - too long - 2.4bn rows…)
exp(log(nrow(sample_events)) + log(nrow(spp)))
# could using a sparse matrix help? But unsure how to assemble that easily… https://www.geeksforgeeks.org/working-with-sparse-matrices-in-r-programming/

# Try just getting a summary in a function:
get_presabs_sp_summ <- function(sp_index, sp_df = spp){
  
  # sp_index = spp$aphia_index[1]
  # sp_df = spp
  # samps = sample_events
  
  # get pres-abs data
  presabs <- ncvar_get(nc_fil, "Pres_abs",
                       start = c(1, sp_index), count = c(-1,1))
  # assemble data frame
    ret <- tibble(presabs = presabs) %>%
      count(presabs) %>%
      mutate(aphia_id = spp$AphiaID[spp$aphia_index == sp_index]) %>% 
      dplyr::select(aphia_id, everything())
  
  ret
}

get_presabs_sp_summ(sp_index = spp$aphia_index[1])
get_presabs_sp_summ(sp_index = 7)

pres_abs_summ <- spp %>% pull(aphia_index) %>%
  purrr::map(get_presabs_sp_summ, .progress = TRUE) %>%
  bind_rows()

pres_abs_summ %>% group_by(presabs) %>% summarise(n_cases = sum(n))

# manageable number of presences; too many NAs and 0s. 
# Interesting headline numbers though: ~2bn absences, 3.3M presences, 388M NAs.
# How many sample events would we bin if we canned all those with NAs? Check this next.
# Think about how to make this work…

sample_events <- sample_events %>% mutate(includes_na = FALSE)

get_na_events <- function(sp_index, samps = sample_events){
  
  # get pres-abs data
  presabs <- ncvar_get(nc_fil, "Pres_abs",
                       start = c(1, sp_index), count = c(-1,1))
  # flag NAs in sample_events
  samps <- samps %>% mutate(includes_na = ifelse(is.na(presabs), TRUE, includes_na))
  samps
}

na_events <- sample_events
for(i in spp$aphia_index){
  na_events <- get_na_events(sp_index = i, samps = na_events)
}

na_events %>% count(includes_na)
# So getting rid of these would remove ~16.5% of sample events
# quick look at dates:
range(na_events$year[na_events$includes_na == TRUE])
range(na_events$year[na_events$includes_na == FALSE])
ggplot(na_events) + geom_histogram(aes(x = year, fill = includes_na), alpha = 2/3)

sample_events <- na_events

# add a sample id variable:
sample_events <- sample_events %>% mutate(sample_id = row_number()) %>% 
  dplyr::select(sample_id, everything())

# we also need to exclude sites with sample dates in the future:
sample_events %>% filter(date > lubridate::now())

# simple way to do this is to set includes_na to TRUE for these sites:
sample_events <- sample_events %>%
  mutate(includes_na = ifelse(date > lubridate::now(), TRUE, includes_na))

# write this to file
write_csv(sample_events, file = here("data", "EMODnet-Biology-Benthos-European-Seas/sample_events.csv"))

sample_events <- read_csv(here("data", "EMODnet-Biology-Benthos-European-Seas/sample_events.csv"))


# Assemble the long species presence x site df, for non-na sites only. Need a 'sites to include' vector.

# in a function:
get_presabs_sp <- function(sp_index, sp_df = spp, samps = sample_events){
  
  # get pres-abs data
  presabs <- ncvar_get(nc_fil, "Pres_abs",
                       start = c(1, sp_index), count = c(-1,1))
  # assemble data frame
  ret <- tibble(aphia_id = spp$AphiaID[spp$aphia_index == sp_index],
                presabs = presabs)
  ret <- samps %>% bind_cols(ret)
  
  ret <- ret %>% filter(presabs == 1 & includes_na == FALSE) %>% 
    dplyr::select(sample_id, aphia_id)
  
  ret
}

get_presabs_sp(sp_index = spp$aphia_index[1])

# generate full dataframe of species presence at each site (excluding sites with NA values)
pres_df <- spp %>% pull(aphia_index) %>%
  purrr::map(get_presabs_sp, .progress = TRUE) %>%
  bind_rows()

sp_occ_summary <- pres_df %>% count(aphia_id) %>% arrange(n)
ggplot(sp_occ_summary) + geom_histogram(aes(x = n))
sum(sp_occ_summary$n == 1)
sum(sp_occ_summary$n > 1)

# get some info on CRS:
presabs_crs <- ncatt_get(nc_fil, "crs")

nc_close(nc_fil)

# Then proceed as in previous trends code
# First save the full list of sample events (including NAs) to a new object:
sample_events_includingNA <- sample_events
# Then filter sample_events to exclude samples which have NAs
sample_events <- sample_events %>% filter(includes_na == FALSE)

# simple frequency plot:
ggplot(sample_events, aes(x = year)) + geom_histogram()

# five year blocks:
sample_events <- sample_events %>% 
  mutate(year_5 = round(year/5)*5)

sample_events %>% count(year_5) %>% print(n = 100)
ggplot(sample_events) + geom_bar(aes(x = year_5))

# samples per grid square per unit time?
library(sf)
sample_events_s <- st_as_sf(sample_events,
                              coords = c("lon", "lat"),
                              crs = 4326)
st_bbox(sample_events_s)

# rough plot:
ggplot(sample_events_s) + geom_sf()

# define some time periods
sample_events <- sample_events %>%
  mutate(time_slice = case_when(
    year < 1990 ~ 1,
    year >= 1990 & year <2000 ~ 2,
    year >= 2000 & year <2005 ~ 3,
    year >= 2005 & year <2010 ~ 4,
    year >= 2010 & year <2015 ~ 5,
    year >= 2015 ~ 6
  ))

#################

# create raster - this is 0.5 degree resolution, but product uses 1 degree - see below. Product also uses terra and tidyterra rather than raster as here


extent_tb <- tibble(lon = c(-33.5, 59), lat = c(28, 81.5)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
r <- raster(extent(extent_tb), resolution = 0.5, crs = st_crs(extent_tb)$proj4string)

# all years
sample_events_r <- rasterize(dplyr::select(sample_events, lon, lat), r,
                               fun = function(x,...) length(x))
plot(sample_events_r, col = viridis(64))
plot(log10(sample_events_r), col = viridis(64))

# by year
samp_ev_1 <- rasterize(dplyr::select(filter(sample_events, time_slice == 1),
                                         lon, lat), r,
                           fun = function(x,...) length(x))
samp_ev_2 <- rasterize(dplyr::select(filter(sample_events, time_slice == 2),
                                     lon, lat), r,
                       fun = function(x,...) length(x))
samp_ev_3 <- rasterize(dplyr::select(filter(sample_events, time_slice == 3),
                                     lon, lat), r,
                       fun = function(x,...) length(x))
samp_ev_4 <- rasterize(dplyr::select(filter(sample_events, time_slice == 4),
                                     lon, lat), r,
                       fun = function(x,...) length(x))
samp_ev_5 <- rasterize(dplyr::select(filter(sample_events, time_slice == 5),
                                     lon, lat), r,
                       fun = function(x,...) length(x))
samp_ev_6 <- rasterize(dplyr::select(filter(sample_events, time_slice == 6),
                                     lon, lat), r,
                       fun = function(x,...) length(x))

# combine into a stack, and then remove the individual slice rasters
samp_ev_by_time_r <- raster::stack(samp_ev_1, samp_ev_2,
                                   samp_ev_3, samp_ev_4,
                                   samp_ev_5, samp_ev_6)

rm(samp_ev_1, samp_ev_2, samp_ev_3, samp_ev_4, samp_ev_5, samp_ev_6)

# simple plot
names(samp_ev_by_time_r) <- c("pre 1990", "1990s", "2000-04", "2005-09",
                                  "2010-14", "post 2015")
plot(log10(samp_ev_by_time_r), col = viridis(64))

# get grid cell + time slice for each sample
# group by this
# for each group:
## create species x time matrix
## calculate beta / turnover stats
## store somewhere with grid ID

#sample_cell <- cellFromXY(r, as.matrix(dplyr::select(sample_events, lon, lat)))


sample_events <- sample_events %>% 
  mutate(sample_cell = cellFromXY(r,
                                  as.matrix(dplyr::select(sample_events, lon, lat))))

# join this to the big species presence dataset and summarise by grid cell
gridded_occs <- pres_df %>% left_join(sample_events, join_by(sample_id)) %>% 
  group_by(aphia_id, time_slice, sample_cell) %>% 
  summarise(sp_occs = n()) %>% 
  arrange(sample_cell, time_slice) %>% ungroup()

# look at number of time slices per sampling cell:
gridded_occs %>% dplyr::select(sample_cell, time_slice) %>% distinct() %>% count(sample_cell) %>% count(n)

cells_single_time <- gridded_occs %>% dplyr::select(sample_cell, time_slice) %>% distinct() %>%
  count(sample_cell) %>% filter(n == 1) %>% 
  left_join(sample_events, join_by(sample_cell)) %>% 
  distinct(sample_cell, .keep_all = TRUE) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

cells_multiple_time <- gridded_occs %>% dplyr::select(sample_cell, time_slice) %>% distinct() %>%
  count(sample_cell) %>% filter(n > 1) %>% 
  left_join(sample_events, join_by(sample_cell)) %>% 
  distinct(sample_cell, .keep_all = TRUE) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

single_time_plot <- ggplot(cells_single_time) + geom_sf()
multi_time_plot <- ggplot(cells_multiple_time) + geom_sf()

single_time_plot + multi_time_plot

# nearly half of grid cells (730) have samples from only one time slice
# Increase time slices, increase grid size, or both?

######################
# Re-do at 1 degree resolution to increase number of repeat-sampled cells
# Also use terra

library(terra)
library(tidyterra)
library(viridis)

extent_tb <- tibble(lon = c(-34, 59), lat = c(28, 82)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

r <- rast(ext(extent_tb), resolution = 1, crs = crs(extent_tb))

# all years
sample_events_r <- sample_events %>%
  dplyr::select(lon, lat) %>% 
  as.matrix() %>% 
  rasterize(r, fun = "length")

plot(sample_events_r, col = viridis(64))
plot(log10(sample_events_r), col = viridis(64))

world <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") %>% 
  st_transform(crs = crs(sample_events_r))



ggplot() +
  geom_spatraster(data = log10(sample_events_r)) +
  scale_fill_viridis_c(name = "log N Samples") +
  geom_sf(data = world, colour = "grey95", fill = NA, size = 0.1) +
  xlim(as.vector(ext(sample_events_r)[1:2])) +
  ylim(as.vector(ext(sample_events_r)[3:4])) +
  coord_sf(expand = FALSE)


# by year
samp_ev_1 <- sample_events %>%
  filter(time_slice == 1) %>% 
  dplyr::select(lon, lat) %>% 
  as.matrix() %>% 
  rasterize(r, fun = "length") %>% 
  rename(time1 = lyr.1)

samp_ev_2 <- sample_events %>%
  filter(time_slice == 2) %>% 
  dplyr::select(lon, lat) %>% 
  as.matrix() %>% 
  rasterize(r, fun = "length") %>% 
  rename(time2 = lyr.1)

samp_ev_3 <- sample_events %>%
  filter(time_slice == 3) %>% 
  dplyr::select(lon, lat) %>% 
  as.matrix() %>% 
  rasterize(r, fun = "length") %>% 
  rename(time3 = lyr.1)

samp_ev_4 <- sample_events %>%
  filter(time_slice == 4) %>% 
  dplyr::select(lon, lat) %>% 
  as.matrix() %>% 
  rasterize(r, fun = "length") %>% 
  rename(time4 = lyr.1)

samp_ev_5 <- sample_events %>%
  filter(time_slice == 5) %>% 
  dplyr::select(lon, lat) %>% 
  as.matrix() %>% 
  rasterize(r, fun = "length") %>% 
  rename(time5 = lyr.1)

samp_ev_6 <- sample_events %>%
  filter(time_slice == 6) %>% 
  dplyr::select(lon, lat) %>% 
  as.matrix() %>% 
  rasterize(r, fun = "length") %>% 
  rename(time6 = lyr.1)


# combine into a stack, and then remove the individual slice rasters
samp_ev_by_time_r <- c(samp_ev_1, samp_ev_2, samp_ev_3, samp_ev_4, samp_ev_5, samp_ev_6)

rm(samp_ev_1, samp_ev_2, samp_ev_3, samp_ev_4, samp_ev_5, samp_ev_6)
names(samp_ev_by_time_r) <- c("pre 1990", "1990s", "2000-04", "2005-09",
                              "2010-14", "post 2015")

# simple plot
plot(log10(samp_ev_by_time_r), col = viridis(64))

ggplot() +
  geom_spatraster(data = log10(samp_ev_by_time_r)) +
  scale_fill_viridis_c(name = "log N Samples") +
  geom_sf(data = world, colour = "grey95", fill = NA, size = 0.1) +
  xlim(as.vector(ext(sample_events_r)[1:2])) +
  ylim(as.vector(ext(sample_events_r)[3:4])) +
  coord_sf(expand = FALSE) +
  facet_wrap(~lyr)


# get grid cell + time slice for each sample
# group by this
# for each group:
## create species x time matrix
## calculate beta / turnover stats
## store somewhere with grid ID

sample_events <- sample_events %>% 
  mutate(sample_cell = cellFromXY(r,
                                  as.matrix(dplyr::select(sample_events, lon, lat))))

# join this to the big species presence dataset and summarise by grid cell
gridded_occs <- pres_df %>% left_join(sample_events, join_by(sample_id)) %>% 
  group_by(aphia_id, time_slice, sample_cell) %>% 
  summarise(sp_occs = n()) %>% 
  arrange(sample_cell, time_slice) %>% ungroup()

# kind of interesting to look at individual species turnover within and between cells, e.g.
# gridded_occs %>% filter(aphia_id == 101160) %>% arrange(sample_cell, time_slice) %>% print(n = 1000)

# look at number of time slices per sampling cell:
gridded_occs %>% dplyr::select(sample_cell, time_slice) %>% distinct() %>% count(sample_cell) %>% count(n)

cells_single_time <- gridded_occs %>% dplyr::select(sample_cell, time_slice) %>% distinct() %>%
  count(sample_cell) %>% filter(n == 1) %>% 
  left_join(sample_events, join_by(sample_cell)) %>% 
  distinct(sample_cell, .keep_all = TRUE) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

cells_multiple_time <- gridded_occs %>% dplyr::select(sample_cell, time_slice) %>% distinct() %>%
  count(sample_cell) %>% filter(n > 1) %>% 
  left_join(sample_events, join_by(sample_cell)) %>% 
  distinct(sample_cell, .keep_all = TRUE) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

single_time_plot <- ggplot(cells_single_time) + geom_sf()
multi_time_plot <- ggplot(cells_multiple_time) + geom_sf()

single_time_plot + multi_time_plot

# still >40% cells have only a single time slice, but more have multiple (e.g.6)
gridded_occs %>% dplyr::select(sample_cell) %>% distinct()

  
# get rid of cells with only one sampling period
multi_time_cells <- gridded_occs %>%
  # get distinct combinations of sample_cell and time_slice from gridded_occs
  dplyr::select(sample_cell, time_slice) %>%
  distinct() %>%
  # identify sample cells with samples from >1 time slice
  count(sample_cell) %>%
  filter(n > 1)

# filter gridded_occs to these cells only
gridded_occs <- gridded_occs %>% filter(sample_cell %in% multi_time_cells$sample_cell)
grid_coords <- xyFromCell(r, gridded_occs$sample_cell) %>% as_tibble() %>% rename(lon = x, lat = y)
gridded_occs <- gridded_occs %>% bind_cols(grid_coords)

# try a sample cell:
gridded_occs %>% dplyr::select(time_slice, sample_cell) %>%
  distinct() %>%
  count(sample_cell) %>%
  arrange(desc(n))

gridded_occs %>% filter(sample_cell == 1729) %>% count(time_slice)

# 5 time slice example:
# get a matrix of species presence/absence for a specific grid cell:
test_mat <- gridded_occs %>% filter(sample_cell == 1729) %>% 
  mutate(pa = ifelse(sp_occs > 0, 1, 0)) %>% 
  dplyr::select(time_slice, aphia_id, pa) %>% 
  pivot_wider(names_from = aphia_id, values_from = pa, values_fill = 0) %>% 
  column_to_rownames(var = "time_slice")
#
test_b <- test_mat %>% 
  as.matrix() %>% beta()

# To get a, b and c for each comparison:
test_abc <- test_mat %>% betadiver(method = NULL)

sp_shared <- data.frame(as.matrix(test_abc$a)) %>% 
  rownames_to_column("time_start") %>% 
  pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "a") %>% 
  mutate(time_comp = str_sub(time_comp, start = 2),
         time_start = as.numeric(time_start),
         time_comp = as.numeric(time_comp)) %>% 
  filter(time_start < time_comp)

sp_lost <- data.frame(as.matrix(test_abc$b)) %>% 
  rownames_to_column("time_start") %>% 
  pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "b") %>% 
  mutate(time_comp = str_sub(time_comp, start = 2),
         time_start = as.numeric(time_start),
         time_comp = as.numeric(time_comp)) %>% 
  filter(time_start < time_comp)

sp_gained <- data.frame(as.matrix(test_abc$c)) %>% 
  rownames_to_column("time_start") %>% 
  pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "c") %>% 
  mutate(time_comp = str_sub(time_comp, start = 2),
         time_start = as.numeric(time_start),
         time_comp = as.numeric(time_comp)) %>% 
  filter(time_start < time_comp)

beta_tot <- data.frame(as.matrix(test_b$Btotal)) %>% rownames_to_column("time_start") %>% 
  pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "beta_tot") %>% 
  mutate(time_comp = str_sub(time_comp, start = 2),
         time_start = as.numeric(time_start),
         time_comp = as.numeric(time_comp)) %>% 
  filter(time_start < time_comp)

beta_repl <- data.frame(as.matrix(test_b$Brepl)) %>% rownames_to_column("time_start") %>% 
  pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "beta_repl") %>% 
  mutate(time_comp = str_sub(time_comp, start = 2),
         time_start = as.numeric(time_start),
         time_comp = as.numeric(time_comp)) %>% 
  filter(time_start < time_comp)

beta_rich <- data.frame(as.matrix(test_b$Brich)) %>% rownames_to_column("time_start") %>% 
  pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "beta_rich") %>% 
  mutate(time_comp = str_sub(time_comp, start = 2),
         time_start = as.numeric(time_start),
         time_comp = as.numeric(time_comp)) %>% 
  filter(time_start < time_comp)

# function this all up…
get_species_turnover <- function(comm_dat){
  
  # This gets some measures of species turnover for a community dataset
  # representing a grid cell through time.
  # Community data comm_dat should have columns time_slice, aphia_id, sp_occs
  # (i.e. in the format of gridded_occs) but only have data for a single sample_cell
  
  # first record sample cell
  sample_cell <- unique(comm_dat$sample_cell)
  
  comm_dat <- comm_dat %>% 
    mutate(pa = ifelse(sp_occs > 0, 1, 0)) %>% 
    dplyr::select(time_slice, aphia_id, pa) %>% 
    pivot_wider(names_from = aphia_id, values_from = pa, values_fill = 0) %>% 
    column_to_rownames(var = "time_slice")
  
  
  # Get incidence-weighted jaccard beta diversity
  # partitioned into total, replacement, richness components
  beta_div <- comm_dat %>% as.matrix() %>% BAT::beta()
  
  # Get each time-slice comparison in a tidy format - here for total beta
  beta_tot <- data.frame(as.matrix(beta_div$Btotal)) %>% rownames_to_column("time_start") %>% 
    pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "beta_tot") %>% 
    mutate(time_comp = str_sub(time_comp, start = 2),
           time_start = as.numeric(time_start),
           time_comp = as.numeric(time_comp)) %>% 
    filter(time_start < time_comp)
  
  # Then for replacement beta
  beta_repl <- data.frame(as.matrix(beta_div$Brepl)) %>% rownames_to_column("time_start") %>% 
    pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "beta_repl") %>% 
    mutate(time_comp = str_sub(time_comp, start = 2),
           time_start = as.numeric(time_start),
           time_comp = as.numeric(time_comp)) %>% 
    filter(time_start < time_comp)
  
  # Finally for richness beta:
  beta_rich <- data.frame(as.matrix(beta_div$Brich)) %>% rownames_to_column("time_start") %>% 
    pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "beta_rich") %>% 
    mutate(time_comp = str_sub(time_comp, start = 2),
           time_start = as.numeric(time_start),
           time_comp = as.numeric(time_comp)) %>% 
    filter(time_start < time_comp)
  
  # We also want a, b, and c (species in common, species lost, species gained) for each comparison
  beta_abc <- comm_dat %>% vegan::betadiver(method = NULL)
  
  # Get shared species (a) in tidy format
  sp_shared <- data.frame(as.matrix(beta_abc$a)) %>% 
    rownames_to_column("time_start") %>% 
    pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "a") %>% 
    mutate(time_comp = str_sub(time_comp, start = 2),
           time_start = as.numeric(time_start),
           time_comp = as.numeric(time_comp)) %>% 
    filter(time_start < time_comp)
  
  # Get species lost between time slices (b)
  sp_lost <- data.frame(as.matrix(beta_abc$b)) %>% 
    rownames_to_column("time_start") %>% 
    pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "b") %>% 
    mutate(time_comp = str_sub(time_comp, start = 2),
           time_start = as.numeric(time_start),
           time_comp = as.numeric(time_comp)) %>% 
    filter(time_start < time_comp)
  
  # Get species gained between time slices (c)
  sp_gained <- data.frame(as.matrix(beta_abc$c)) %>% 
    rownames_to_column("time_start") %>% 
    pivot_longer(cols = -time_start, names_to = "time_comp", values_to = "c") %>% 
    mutate(time_comp = str_sub(time_comp, start = 2),
           time_start = as.numeric(time_start),
           time_comp = as.numeric(time_comp)) %>% 
    filter(time_start < time_comp)
  
  
  # join all of these together:
  
  joined <- list(beta_tot, beta_repl, beta_rich,
                 sp_shared, sp_lost, sp_gained) %>% 
    purrr::reduce(left_join, by = c("time_start", "time_comp"))
  
  
  joined <- joined %>% mutate(sample_cell = sample_cell) %>% 
    dplyr::select(sample_cell, everything())
  
  joined
}

test <- gridded_occs %>% filter(sample_cell == 1729) 
# run for one cell
get_species_turnover(test)
rm(test)

# run over all cells
species_turnover_allcells <- gridded_occs %>%
  split(.$sample_cell) %>%
  purrr::map(get_species_turnover, .progress = TRUE) %>% 
  bind_rows()

# get a time comparison variable
species_turnover_allcells <- species_turnover_allcells %>% 
  mutate(t_comp = paste(time_start, time_comp, sep = "_"),
         t_steps = time_comp - time_start) %>% 
  dplyr::select(sample_cell, t_comp, t_steps, beta_tot:c)

# add lon and lat for the cells
cell_latlon <- gridded_occs %>% dplyr::select(sample_cell, lon, lat) %>% distinct()
# get total distinct samples per cell
sample_events_per_cell <- sample_events %>% group_by(sample_cell) %>% summarise(n_samps = n_distinct(sample_id))
# join to cell_latlon
cell_latlon <- cell_latlon %>% left_join(sample_events_per_cell, join_by(sample_cell))

species_turnover_allcells <- species_turnover_allcells %>% 
  left_join(cell_latlon, join_by(sample_cell))
  

# restrict to single time step changes, create proportional species gain / loss, and make spatial
species_turnover_single_t <- species_turnover_allcells %>% filter(t_steps == 1) %>% 
  mutate(b_prop = b / (a + b + c), c_prop = c / (a + b + c)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

ggplot(filter(species_turnover_single_t, t_comp == "4_5")) +
  coord_sf() +
  geom_raster(aes(x = lon, y = lat, fill = b_prop)) +
  scale_fill_viridis()


# try biscale map
sp_dat <- bi_class(species_turnover_single_t, x = b_prop, y = c_prop, style = "quantile", dim = 3)

# get basemap
europe <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf') %>% 
  st_transform(crs = st_crs(sample_events_s))
europe <- europe %>% st_make_valid() %>% 
  st_crop(xmin = -33, xmax = 44, ymin = 20, ymax = 85)
ggplot(europe) + geom_sf() +
  ylim(28, 79)

# lookup table for labels of time comparisons
time_comp <- c(
  "1_2" = "1990s v pre 1990",
  "2_3" = "2000-2004 v 1990s",
  "3_4" = "2005-2009 v 2000-2004",
  "4_5" = "2010-2014 V 2005-2009",
  "5_6" = "post 2015 v 2010-2014"
)

# plot b_prop, transparency scaled by total number of sampling events per cell
ggplot() +
  geom_sf(data = europe) +
  xlim(-33, 44) +
  ylim(28, 79) +
  geom_raster(data = species_turnover_single_t, aes(x = lon, y = lat, fill = b_prop, alpha = log10(n_samps))) +
  xlab("") + ylab("") +
  scale_fill_viridis() +
  facet_wrap(~t_comp, labeller = labeller(t_comp = time_comp)) +
  theme_bw()


ggplot() +
  geom_sf(data = europe) +
  xlim(-33, 44) +
  ylim(28, 79) +
  geom_raster(data = species_turnover_single_t, aes(x = lon, y = lat, fill = beta_tot, alpha = log10(n_samps))) +
  xlab("") + ylab("") +
  scale_fill_viridis() +
  facet_wrap(~t_comp, labeller = labeller(t_comp = time_comp)) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.2))


# this uses species gained and lost as a proportion of all species recorded in the cell across the two time periods
loss_gain <- ggplot() +
  geom_sf(data = europe) +
  xlim(-33, 44) +
  ylim(28, 79) +
  geom_raster(data = sp_dat, aes(x = lon, y = lat, fill = bi_class), show.legend = FALSE) +
  xlab("") + ylab("") +
  bi_scale_fill(pal = "DkBlue2", dim = 4) +
  facet_wrap(~t_comp, labeller = labeller(t_comp = time_comp)) +
  theme_bw()

lg_legend <-  bi_legend(pal = "DkBlue2", dim = 4, xlab = "species lost", ylab = "species gained")

loss_gain + inset_element(lg_legend, left = 0.66, bottom = 0.1, right = 1, top = 0.43)

sp_change_maps <- loss_gain + inset_element(lg_legend, left = 0.66, bottom = 0.1, right = 1, top = 0.43)
ggsave(filename = "sp_change_maps_1deg.png", plot = sp_change_maps, width = 24, height = 16, units = "cm")

## Might be useful to provide species x site x time too as a simpler 4D array


