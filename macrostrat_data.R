# import Macrostrat database files
library(jsonlite)
library(dplyr)

macrostrat_data <- function(filename, data_url){
  if (!file.exists(filename)) {
    download.file(data_url, filename, method="auto")
  }
  data <- fromJSON(filename)
  return(data[["success"]][["data"]])
}

sedimentary_rocks <- function() {
  rocks <- macrostrat_data("liths.json", "https://macrostrat.org/api/defs/lithologies?all")
  rocks <- rocks %>% filter(class == 'sedimentary')
  rocks <- rocks$name
  return(rocks)
}


