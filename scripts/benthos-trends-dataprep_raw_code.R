## ----setup, include=FALSE------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- load_packages, message = FALSE-------------------------------------------------------
library(tidyverse)
library(ncdf4)
library(here)
library(worrms)


## ---- open_nc_connection-------------------------------------------------------------------
nc_fil <- nc_open(here("data", "raw_data/Macrobenthos_Eur_Seas_Pres_Abs_v0-4.nc"))


## ---- get_sampling_events------------------------------------------------------------------
sample_events <- tibble(
  lat = as.vector(ncvar_get(nc_fil, "lat")),
  lon = as.vector(ncvar_get(nc_fil,"lon")),
  date = lubridate::ymd(as.Date(ncvar_get(nc_fil, "Date"), origin = "1970-01-01 00:00:00")),
  year = lubridate::year(as.Date(ncvar_get(nc_fil, "Date"), origin = "1970-01-01 00:00:00"))
)


## ---- future_samples, echo = FALSE---------------------------------------------------------
sample_events %>% filter(date > lubridate::now())


## ---- get_taxa-----------------------------------------------------------------------------
taxon_name <- ncvar_get(nc_fil, "Taxon_Name")
aphia_id <- ncvar_get(nc_fil,"AphiaID") %>% word(start = -1, sep = ":") %>% as.integer()


## ---- get_taxonrank------------------------------------------------------------------------
get_worms_rank <- function(aphiaid){
  classif <- wm_classification(aphiaid) %>% slice(n())
  classif
}


## ---- get_taxonrank_example----------------------------------------------------------------
get_worms_rank(aphia_id[1000])


## ---- get_taxonrank_alltaxa, eval = FALSE--------------------------------------------------
## taxo <- aphia_id %>% map_df(get_worms_rank, .progress = TRUE)


## ---- read_taxo, message = FALSE-----------------------------------------------------------
taxo <- read_csv(here("data", "derived_data/benthos_taxo_rank.csv"))


## ---- print_taxonranks---------------------------------------------------------------------
taxo %>% count(rank) %>% arrange(desc(n)) %>% print(n = 100)


## ---- get_speciesrank----------------------------------------------------------------------
spp <- taxo %>%
  filter(rank == "Species") %>% 
  mutate(aphia_index = match(.$AphiaID, taxo$AphiaID))


## ---- get_presabs_singlespecies------------------------------------------------------------
presabs <- ncvar_get(nc_fil,"Pres_abs",
                     start = c(1, spp$aphia_index[1]), count = c(-1,1))


## ---- get_presabs_species------------------------------------------------------------------
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



## ---- get_presabs_eg-----------------------------------------------------------------------

presabs <- get_presabs_sp(sp_index = spp$aphia_index[1])
presabs %>% count(presabs)



## ---- get_presabs_summary------------------------------------------------------------------
get_presabs_sp_summ <- function(sp_index, sp_df = spp){
  
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



## ---- get_presabs_summary_eg---------------------------------------------------------------
get_presabs_sp_summ(sp_index = spp$aphia_index[1])


## ---- get_presabs_sunmary_allspp, message = FALSE------------------------------------------
pres_abs_summ <- spp %>% pull(aphia_index) %>%
  purrr::map(get_presabs_sp_summ, .progress = TRUE) %>%
  bind_rows()

pres_abs_summ %>% group_by(presabs) %>% summarise(n_cases = sum(n))


## ---- create_includes_na-------------------------------------------------------------------
sample_events <- sample_events %>% mutate(includes_na = FALSE)


## ---- get_na_events------------------------------------------------------------------------
get_na_events <- function(sp_index, samps = sample_events){
  
  # get pres-abs data
  presabs <- ncvar_get(nc_fil, "Pres_abs",
                       start = c(1, sp_index), count = c(-1,1))
  # flag NAs in sample_events
  samps <- samps %>% mutate(includes_na = ifelse(is.na(presabs), TRUE, includes_na))
  samps
}



## ---- get_na_events_allspp-----------------------------------------------------------------
na_events <- sample_events
for(i in spp$aphia_index){
  na_events <- get_na_events(sp_index = i, samps = na_events)
}



## ---- summarise_includes_na----------------------------------------------------------------
na_events %>% count(includes_na)


## ---- echo = FALSE, message = FALSE--------------------------------------------------------
ggplot(na_events) + geom_histogram(aes(x = year, fill = includes_na), alpha = 2/3)


## ---- finalise_sample_events---------------------------------------------------------------

sample_events <- na_events %>% 
  mutate(sample_id = row_number()) %>% 
  dplyr::select(sample_id, everything()) %>% 
  mutate(includes_na = ifelse(date > lubridate::now(), TRUE, includes_na))



## ---- read_sample_events-------------------------------------------------------------------
sample_events <- read_csv(here("data",
                               "derived_data/sample_events.csv"))



## ---- get_presabs_sp_na_rm-----------------------------------------------------------------

get_presabs_sp_na_rm <- function(sp_index, sp_df = spp, samps = sample_events){
  
  # get pres-abs data
  presabs <- ncvar_get(nc_fil, "Pres_abs",
                       start = c(1, sp_index), count = c(-1,1))
  # assemble data frame
  ret <- tibble(aphia_id = spp$AphiaID[spp$aphia_index == sp_index],
                presabs = presabs)
  ret <- samps %>% bind_cols(ret)
  
  # filter to presences and non-NA sites
  ret <- ret %>% filter(presabs == 1 & includes_na == FALSE) %>% 
    dplyr::select(sample_id, aphia_id)
  
  ret
}



## ---- read_species_list, message = FALSE---------------------------------------------------
taxo <- read_csv(here("data", "derived_data/benthos_taxo_rank.csv"))
spp <- taxo %>% 
  filter(rank == "Species") %>% 
  mutate(aphia_index = match(.$AphiaID, taxo$AphiaID))


## ---- get_presabs_sp_na_rm_eg--------------------------------------------------------------
get_presabs_sp_na_rm(sp_index = spp$aphia_index[1])


## ---- get_presabs_sp_na_rm_allspp, message = FALSE-----------------------------------------

pres_df <- spp %>% pull(aphia_index) %>%
  purrr::map(get_presabs_sp_na_rm, .progress = TRUE) %>%
  bind_rows()



## ---- get_sp_occ_summary, message = FALSE--------------------------------------------------
sp_occ_summary <- pres_df %>% count(aphia_id) %>% arrange(n)
ggplot(sp_occ_summary) + geom_histogram(aes(x = n))


## ---- get_presabs_crs----------------------------------------------------------------------
presabs_crs <- ncatt_get(nc_fil, "crs")


## ---- read_pres_df, eval = FALSE-----------------------------------------------------------
## pres_df <- read_csv(here("data",
##                          "derived_data/pres_df.csv"))
## 


## ---- close_nc-----------------------------------------------------------------------------
nc_close(nc_fil)


## ---- reproducibility_date, echo = FALSE---------------------------------------------------
lubridate::now()


## ---- session_info, echo = FALSE-----------------------------------------------------------
sessionInfo()

