# Temporal Turnover in European Macrobenthos Communities

## Introduction

This product builds on the EMODnet Biology data product [Presence/absence data of macrozoobenthos in the European Seas](https://github.com/EMODnet/EMODnet-Biology-Benthos-European-Seas) to derive estimates of temporal turnover in benthic communities on a spatial grid across European seas. This product only uses species-level records, and only uses sampling events where the full macrobenthic community was surveyed (i.e. where there are no 'NA' values in the presence/absence dataset for any species). Six time periods are considered, based on data availability: before 1990, 1990-1999, 2000-2004, 2005-2009, 2010-2014, and 2015 and after. A 1 degree grid is used to obtain reasonable numbers of repeat samples per grid cell. The code below could be adapted to set different time periods and/or a different grid resolution. This readme describes the product structure, including the workflow to generate the required derived datasets and the process for turning them into gridded maps of community turnover.

## Directory structure

```
EMODnet_benthos_trends/
├── analysis
├── data/
│   ├── derived_data/
│   └── raw_data/
├── docs/
├── product/
└── scripts/
```

* **analysis** - Markdown or Jupyter notebooks
* **data** - Raw and derived data
* **docs** - Rendered reports
* **product** - Output product files
* **scripts** - Reusable code

## Data series

This product uses the EMODnet Biology data product [Presence/absence data of macrozoobenthos in the European Seas](https://github.com/EMODnet/EMODnet-Biology-Benthos-European-Seas). This is available as a single NetCDF file - because that file is large (~3.2GB) it is not included here. Rather, this product uses two datasets that are derived from the NetCDF file, and which are available in `data/derived_data`: `sample_events.csv` is a table of unique sampling events, including latitude, longitude, and sampling date; and `pres_df.csv` is a table of species presences, referenced by sample id (linked to `sample_events`) and species identity given as [WoRMS AphiaID](https://marinespecies.org/about.php). The full workflow for deriving these datasets from the raw NetCDF file is described in `docs/benthos-trends-dataprep`.

## Data product

The product (available in the `product` directory in this repository) includes a series of gridded estimates of macrobenthic community turnover in European seas. The six time periods are as defined above (before 1990, 1990-1999, 2000-2004, 2005-2009, 2010-2014, and 2015 and after), and turnover metrics are calcuated at the cell level for all possible pairwise time comparisions - from a minimum of 1 (where a cell contains observations from only two time periods) to a maxmimum of 15 (where a cell contains observations in all six time periods). 

The diversity measures used are: 

*Beta diversity* - these metrics are derived using the `beta` function in the `BAT` package ([Cardoso et al. 2022](https://CRAN.R-project.org/package=BAT)). There are three measures: Btotal (total beta diversity), which is then decomposed into Brepl (turnover due to species replacement) and Brich (turnover due to changes in species richness), such that Btotal = Brepl + Brich. Filenames for these three measures begin with "beta_tot", "beta_repl", and "beta_rich" for Btotal, Brepl and Brich respectively. Each of these measures can vary between 0 and 1.

We also calculate three measures of species turnover using the `betadiver` function in the `vegan` package ([Oksanen et al. 2022](https://CRAN.R-project.org/package=vegan)) with `method = NULL`, to give *species in common (a)* in a grid cell between two time periods, as well as *species lost (b)* and *species gained (c)*. Filenames for these measures begin with "sp_shared", "sp_lost", and "sp_gained" for a, b and c respectively. Values are in numbers of species, with a minimum of 0 and a maximum of the total number of species within a grid cell across both time periods under comparison.

Finally we also provide b and c as proportions of all species observed in the grid cell across both time periods (i.e. *proportion of species lost, b / (a + b + c)* and *proportion of species gained, c / (a + b + c)*) as measures of relative species turnover. Filnames for these measures begin with "p_sp_lost" and "p_sp_gained" for the proportion of species lost and gained respectively. Values of both can vary between 0 and 1.

Gridded estimates of each diversity measure are produced in R ([R Core Team 2022](https://www.R-project.org/)) as `SpatRaster` objects using the `terra` package ([Hijmans 2023](https://CRAN.R-project.org/package=terra)). They are provided here as A NetCDF file (a single file with each diversity measure x time comparison combination included as a separate variable, giving 120 variables in total) and as GeoTiff files (one file per diversity measure, with one layer per time comparision within each file, as well as a single file of all combinations matching the structure of the NetCDF file).

## More information:

The full workflow for deriving the products from the datasets derived from the European presence/absence data product is set out in `docs/benthos-trends-turnover`.

### References

Cardoso P, Mammola S, Rigal F, Carvalho J (2022). _BAT: Biodiversity Assessment Tools_. R package version 2.9.2, https://CRAN.R-project.org/package=BAT

Herman, P M J (2022) _Summary presence/absence maps of macro-endobenthos in European Seas, based on the EMODNET Biology database_. Integrated data products created under the European Marine Observation Data Network (EMODnet) Biology project Phase IV (EMFF/2019/1.3.1.9/Lot 6/SI2.837974), funded by the by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund. https://www.vliz.be/imis?dasid=8216 

Hijmans R (2023). _terra: Spatial Data Analysis_. R package version 1.7-3, https://CRAN.R-project.org/package=terra

Oksanen J, et al. (2022). _vegan: Community Ecology Package_. R package version 2.6-4, https://CRAN.R-project.org/package=vegan

R Core Team (2022). _R: A language and environment for statistical computing_. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

WoRMS Editorial Board (2023). _World Register of Marine Species_. Available from https://www.marinespecies.org at VLIZ. doi:10.14284/170

### Code and methodology

See `docs/benthos-trends-dataprep` and `docs/benthos-trends-turnover` in this repository for a full description of methodology including relevant code.

### Citation and download link

This product should be cited as:

Webb, T.J. (2023) _Temporal turnover of macrobenthos in European seas_. Integrated data products created under the European Marine Observation Data Network (EMODnet) Biology project Phase IV (EMFF/2019/1.3.1.9/Lot 6/SI2.837974), funded by the by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund.

Available to download in:

https://www.vliz.be/imis?dasid=8229

### Authors

Tom Webb
