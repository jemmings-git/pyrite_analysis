
library(dplyr) # plyr is not needed (note loading plyr after dplyr will prevent execution of PART 2)
library(tidyr)
library(readr)
library(jsonlite)
library(reshape2)
library(docstring)
library(devtools) # run once
#devtools::install_github("moodymudskipper/safejoin") # needed for coalesce join # run once
library(safejoin)

collect_strat_units <- function() {
  #' Collect reference data for stratigraphic units (age and location)
  #' Uses data from the Macrostrat API at http://macrostrat.org/api/

  # import TWO Macrostrat database files - download and keep local copies
  macrostrat_data <- function(filename, data_url){
    #' Download data from API, checking for local cache first
    file_location <- paste(getwd(), "data", filename, sep="/")
    if (!file.exists(file_location)) {
      download.file(data_url, filename, method="auto")
    }
    data <- fromJSON(file_location)
    return(data[["success"]][["data"]])
  }

  # first file is "strat names" including ALL strat_names_ID
  strat <- macrostrat_data("strat.json", "https://macrostrat.org/api/defs/strat_names?all&format=json&response=long")

  # second file is "units" which includes environment of deposition info, palaeolatitude, etc.
  # but note many strat_names_ID listed in 'strat' are not listed in 'units'
  units <- macrostrat_data("units.json", "https://macrostrat.org/api/units?age_top=0&age_bottom=4540000000&response=long&format=json")

  # merge the two Macrostrat databases into one - match to strat_names_ID
  strat$strat_name_id  = as.character(strat$strat_name_id)
  units$strat_name_id  = as.character(units$strat_name_id)

  # merge the GDD extracts and Macrostrat files, using 'strat_name_id' as the unique identifier
  # there are presently two options -

  # OPTION 1 - merge strat and units, and propagate all 'units'
  # normalisation to total 'marine units' should correct for spatial bias
  # but it might also introduce uncertainty - i.e., how can we be sure a given strat ID
  # can be joined to all units associated with the ID?

  # that said, OPTION 1 is presently favoured (and appears to be adopted in Peters et al. 2017)
  strat_units <- safe_right_join(units, strat, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)

  return(strat_units)
}

import_and_strat_align <- function() {
  #' Load our collection of results from GeoDeepDive
  #' Merge with stratigraphic data

  extracts <- read_csv(paste(getwd(), "data", "results.csv", sep="/"))

  # remove unresolved strat_name_id hits

  extracts <- extracts[!grepl(pattern = "\\~", extracts$strat_name_id),]

  # OR OPTION 2 - merged datasets without propagation of 'units'
  #strat_units <- safe_right_join(units, strat, by = c("strat_name_id", "strat_name_id"), conflict = ~.y)
  #strat_units <- strat_units[!duplicated(strat_units[c("strat_name_id", "t_age", "b_age")]),] # remove duplicate IDs but only for those with matching age ranges
  # produces a more 'blocky' and apparently lower resolution output

  extracts$strat_name_id  = as.character(extracts$strat_name_id)

  strat_units <- collect_strat_units()

  data <- safe_right_join(extracts, strat_units, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)
  return(data)

}

  # Options could now be params to function above

  # optional export of the merged data file
  # write.table(data, "data5.txt", sep="\t")

  # optional remove extracts and strat from working environment
  # rm(extracts, strat)

subset_and_bin <- function(data) {
  #' Prepare dataset for normalisation to Macrostrat units #####
  #' - Addition of Phanerozoic time bins between top and base ages ####
  #' subset dataframe into phanerozoic & precambrian (ultimately bins = 1 and 10, respectively)

  cutoff <- 541 # Ma
  phanerozoic_increment <- 1 # in Ma
  precambrian_increment <- 10 # in Ma


  add_all_bins <- function(df, increment) {
    ncols <- ceiling((max(df$b_age - df$t_age))/increment)

    df <- df %>% mutate(bin1 = ifelse(b_age-increment > t_age+increment, t_age+increment, NA))
    for(n in 2:(ncols-1)){
      next_bin <- paste('bin', n, sep='')
      bin <- paste0('bin', (n-1), sep='')
      # See https://stackoverflow.com/a/49311813/4319767 for why this syntax
      df <- df %>% mutate(!!next_bin := ifelse(UQ(rlang::sym(bin)) < b_age-increment, UQ(rlang::sym(bin))+increment, NA))
    }
    return(df)
  }

  data_precambrian <- subset(data, t_age > cutoff)
  data_precambrian <- add_all_bins(data_precambrian, precambrian_increment)

  # including units spanning Phanerozoic-Precambrian boundary

  data_phanerozoic1 <- subset(data, t_age < cutoff & b_age < cutoff)
  data_phanerozoic1 <- add_all_bins(data_phanerozoic1, phanerozoic_increment)

  data_phanerozoic2 <- subset(data, t_age < cutoff & b_age > cutoff)
  data_phanerozoic2$b_age <- cutoff
  data_phanerozoic2 <- add_all_bins(data_phanerozoic2, phanerozoic_increment)

  data_phanerozoic3 <- subset(data, t_age < cutoff & b_age > cutoff)
  data_phanerozoic3$t_age <- cutoff
  data_phanerozoic3 <- add_all_bins(data_phanerozoic3, precambrian_increment)

  # Bind the Phanerozoic and Precambrian datasets

  data <- bind_rows(data_phanerozoic1, data_phanerozoic2, data_phanerozoic3, data_precambrian)

  # export the merged file in order to avoid re-running the above code during analysis

  # first convert the nested lists to character strings (otherwise the .txt export fails)
  return(data)
}

export_data <- function(data) {
  #' Export the data after segment_and_bin phase into a flat file
  data$lith  = as.character(data$lith)
  data$environ  = as.character(data$environ)
  data$econ  = as.character(data$econ)
  data$measure = as.character(data$measure)
  data$units_above = as.character(data$units_above)
  data$units_below = as.character(data$units_below)
  data$refs = as.character(data$refs)
  # export the data
  write.table(data, "data_comp2.txt", sep="\t")
}

simplify_lith_env <- function(data) {
  #' Simplify datafile, simplify lith and environ variables ####

  # convert the dataframe from wide to long format (with respect to the age increments)

  data <- melt(data, id.vars = c(1:20, 23:72), na.rm=TRUE)

  data$value <- as.numeric(data$value)

  # optional compilation of target words
  phrases <- as.data.frame(as.array(data$target_word))

  # simplification of environment of deposition - in this case into Marine and Non-Marine units
  data$envCODE <- gsub("TRUE", "Non-Marine", grepl("non-marine", data$environ))
  data$envCODE <- gsub("FALSE", "Marine", data$envCODE)
  data$envCODE[is.na(data$unit_name)] <- NA

  # simplification of lithology - in this case into sedimentary & non-sedimentary rocks
  data$lithCODE <- gsub("TRUE", "Sedimentary", grepl("sedimentary", data$lith))
  data$lithCODE <- gsub("FALSE", "Non-Sed", data$lithCODE)
  data$lithCODE[is.na(data$unit_name)] <- NA
  return(data)
}
