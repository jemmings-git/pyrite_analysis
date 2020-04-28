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

sedimentary_rocks <- c(sedimentary_rocks(), "sediment")

meta_sedimentary_rocks <- function() {
  rocks <- macrostrat_data("liths.json", "https://macrostrat.org/api/defs/lithologies?all")
  rocks <- rocks %>% filter(class == 'sedimentary' | type == 'metasedimentary')
  rocks <- rocks$name
  return(rocks)
}

meta_sedimentary_rocks <- c(meta_sedimentary_rocks(), "sediment")

mud_rocks <- function() {
  rocks <- macrostrat_data("liths.json", "https://macrostrat.org/api/defs/lithologies?all")
  rocks <- rocks %>% filter(class == 'sedimentary' | type == 'metasedimentary')
  rocks <- rocks$name
  return(rocks)
}

mud_rocks <- mud_rocks()
mud_rocks <- c(mud_rocks[4:9], "mudrock")


environments <- function() {
  rocks <- macrostrat_data("envs.json", "https://macrostrat.org/api/defs/environments?all")
  return(rocks)
}

environments <- environments()

non.mar <- subset(environments, class == "non-marine")
non.mar <- c("eolian", "fluvial", "glacial", "terrestrial", "non-marine", "nonmarine") # simple approach, avoids false hits

marine <- subset(environments, class == "marine")
marine <- marine$name

deep_water <- environments[c(11:13,15:20,28,32:37),]
deep_water <- c(deep_water$name, "deep water","deep marine", "deep-water", "deep-marine")

