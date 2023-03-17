# To render the RMarkdown docs from the analysis folder into the docs folder:


# Data preparation
rmarkdown::render(here::here("analysis", "benthos-trends-dataprep.Rmd"),
                  output_file = here::here("docs", "benthos-trends-dataprep.html"))

# Turnover calculations
rmarkdown::render(here::here("analysis", "benthos-trends-turnover.Rmd"),
                  output_file = here::here("docs", "benthos-trends-turnover.html"))

# Extract code
knitr::purl(input = here::here("analysis", "benthos-trends-dataprep.Rmd"),
            output = here::here("scripts", "benthos-trends-dataprep_raw_code.R"))

knitr::purl(input = here::here("analysis", "benthos-trends-turnover.Rmd"),
            output = here::here("scripts", "benthos-trends-turnover_raw_code.R"))