netcdfify <- function(focal_rast,
                      output_fname,
                      output_dir = "product",
                      global_atts){
  
  # focal_rast is the focal raster, in SpatRaster format
  # output_fname is the name of the output file to produce (with no file extension)
  # output_dir is the directory (relative to the RProject source) in which to save the output file
  # global_atts is a list of the global attributes for the output NetCDF
  
  # packages required
  require(tidyverse)
  require(here)
  require(janitor)
  require(terra)
  require(RNetCDF)
  
  # add diversity variable to names of raster layers and clean names
  names(focal_rast) <- paste(output_fname, names(focal_rast))
  focal_rast <- focal_rast %>% janitor::clean_names()
  
  # create full filename of output file
  output_fname <- paste0(output_fname, ".nc")
  
  # Create nc file
  nc <- create.nc(here(output_dir, output_fname))
  
  # turn raster to dataframe
  focal_df <- focal_rast %>%
    as.data.frame(xy = TRUE, na.rm = FALSE) %>%
    as_tibble() %>% 
    rename(lon = x, lat = y) %>% 
    arrange(lon, lat)
  
  # get lon and lat variables
  focal_lon <- focal_df %>% select(lon) %>% distinct() %>% arrange(lon) %>% pull(lon)
  focal_lat <- focal_df %>% select(lat) %>% distinct() %>% arrange(lat) %>% pull(lat)
  
  # Define lon dimension
  dim.def.nc(nc, dimname = "lon", dimlength = length(focal_lon)) 
  # Define lon variable
  var.def.nc(nc, varname = "lon", vartype = "NC_DOUBLE", dimensions = "lon")
  # Add attributes
  att.put.nc(nc, variable = "lon", name = "units", type = "NC_CHAR", value = "degrees_east")
  att.put.nc(nc, variable = "lon", name = "standard_name", type = "NC_CHAR", value = "longitude")
  att.put.nc(nc, variable = "lon", name = "long_name", type = "NC_CHAR", value = "Longitude")
  # Put data
  var.put.nc(nc, variable = "lon", data = focal_lon) 
  
  # Define lat dimension
  dim.def.nc(nc, dimname = "lat", dimlength = length(focal_lat)) 
  # Define lat variable
  var.def.nc(nc, varname = "lat", vartype = "NC_DOUBLE", dimensions = "lat")
  # Add attributes
  att.put.nc(nc, variable = "lat", name = "units", type = "NC_CHAR", value = "degrees_north")
  att.put.nc(nc, variable = "lat", name = "standard_name", type = "NC_CHAR", value = "latitude")
  att.put.nc(nc, variable = "lat", name = "long_name", type = "NC_CHAR", value = "Latitude")
  # Put data
  var.put.nc(nc, variable = "lat", data = focal_lat) 
  
  # Define non-dimensional crs variable 
  var.def.nc(nc, varname = "crs", vartype = "NC_CHAR", dimensions = NA)
  
  # Add attributes
  att.put.nc(nc, variable = "crs", name = "long_name", type = "NC_CHAR", value = "Coordinate Reference System")
  att.put.nc(nc, variable = "crs", name = "geographic_crs_name", type = "NC_CHAR", value = "WGS 84")
  att.put.nc(nc, variable = "crs", name = "grid_mapping_name", type = "NC_CHAR", value = "latitude_longitude")
  att.put.nc(nc, variable = "crs", name = "reference_ellipsoid_name", type = "NC_CHAR", value = "WGS 84")
  att.put.nc(nc, variable = "crs", name = "horizontal_datum_name", type = "NC_CHAR", value = "WGS 84")
  att.put.nc(nc, variable = "crs", name = "prime_meridian_name", type = "NC_CHAR", value = "Greenwich")
  att.put.nc(nc, variable = "crs", name = "longitude_of_prime_meridian", type = "NC_DOUBLE", value = 0.)
  att.put.nc(nc, variable = "crs", name = "semi_major_axis", type = "NC_DOUBLE", value = 6378137.)
  att.put.nc(nc, variable = "crs", name = "semi_minor_axis", type = "NC_DOUBLE", value = 6356752.314245179)
  att.put.nc(nc, variable = "crs", name = "inverse_flattening", type = "NC_DOUBLE", value = 298.257223563)
  att.put.nc(nc, variable = "crs",name = "spatial_ref", type = "NC_CHAR",
             value = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]')
  att.put.nc(nc, variable = "crs", name = "GeoTransform", type = "NC_CHAR",
             value = '-180 0.08333333333333333 0 90 0 -0.08333333333333333 ')
  
  # get variable names
  vnames <- names(focal_rast)
  
  # create and fill variables for each comparison
  for(i in vnames){
    # Create the diversity metric variable defined by the four dimensions
    var.def.nc(nc, varname = i, vartype = "NC_DOUBLE", dimensions = c("lon", "lat"))
    # Add attributes
    att.put.nc(nc, variable = i, name = "_FillValue", type = "NC_DOUBLE", value = -99999)
    vnames_val <- paste(output_fname, i)
    att.put.nc(nc, variable = i, name = "long_name", type = "NC_CHAR", value = vnames_val)
    
    # add data
    focal_array <- array(
      pull(focal_df, all_of(i)),
      dim = c(length(focal_lat), length(focal_lon))
    )
    var.put.nc(nc, variable = i, data = t(focal_array))
  }
  
  # Define function that detects if the data type should be character of 
  # integer and add to global attributes
  add_global_attributes <- function(nc, attributes){
    
    stopifnot(is.list(attributes))
    
    for(i in 1:length(attributes)){
      if(is.character(attributes[[i]])){
        type <- "NC_CHAR"
      }else if(is.numeric(attributes[[i]])){
        type <- "NC_DOUBLE"
      }
      att.put.nc(nc, variable = "NC_GLOBAL", name = names(attributes[i]), type = type, value = attributes[[i]])
    }
    sync.nc(nc)
  }
  
  # Add attributes
  add_global_attributes(nc, attributes = global_atts)
  
  sync.nc(nc)
  close.nc(nc)
  
}