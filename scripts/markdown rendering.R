# To render the RMarkdown docs from the analysis folder into the docs folder:


# Data preparation
rmarkdown::render(here::here("analysis", "benthos-trends-dataprep.Rmd"),
                  output_file = here::here("docs", "benthos-trends-dataprep.html"))

# Turnover calculations
rmarkdown::render(here::here("analysis", "benthos-trends-turnover.Rmd"),
                  output_file = here::here("docs", "benthos-trends-turnover.html"))