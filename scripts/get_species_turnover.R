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
