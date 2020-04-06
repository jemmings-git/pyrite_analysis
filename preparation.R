##### This R script was written by Joe Emmings and Jo Walsh (British Geological Survey) ####
##### This script is designed to manipulate and process GeoDeepDive extractions
##### interfaced with the Macrostrat database

##### PART 1 - Load and merge the GDD extractions and Macrostrat databases ####

# install the following libraries (and dependencies)

# run once

library(devtools)


# attach

library(dplyr) # plyr is not needed (note loading plyr after dplyr will prevent execution of PART 2)
library(tidyr)
library(readr)
library(jsonlite)
library(reshape2)
library(ggplot2)
library(safejoin)
library(pammtools)
library(gridExtra)

# clear working environment if re-running this script from Part 1 in the same R session
rm(list=ls())

# set working directory
project_home <- 'N:/Data/xGDD/analysis'



source('macrostrat_data.R') # use choices below to find strat flag hits not in concepts database

tryCatch({
  setwd(project_home)
}, error = function(err) {
  if (dir.exists('./data')) {
    setwd('./data') }
}
)
 
rocks <- sedimentary_rocks # normalisation to all sedimentary rocks
rocks <- meta_sedimentary_rocks # normalisation to all sedimentary and metamorphosed sedimentary rocks
rocks <- mud_rocks # normalisation to mudstones - perhaps this is the most robust approach


# import GDD extractions # (note 'results2' is with new pyrite terms 24/10/19

if (!file.exists('results2.csv')) {
  source_data <- 'https://geodeepdive.org/app_output/jemmings_with_pyrite_24Oct2019.zip'
  download.file(source_data, 'jemmings_etal.zip', method='auto')
  unzip('jemmings_etal.zip')
}

extracts <- read_csv("results.csv")

# remove unresolved strat_name_id hits

extracts <- extracts[!grepl(pattern = "\\~", extracts$strat_name_id),]

length(unique(extracts$docid)) # metrics - no. of documents
length(unique(extracts$strat_name_id)) # metrics - no. of strat packages
length(unique(extracts$result_id))

# import strat package Macrostrat database

strat <- macrostrat_data("strat.json", "https://macrostrat.org/api/defs/strat_names?all&format=json&response=long")

# import strat concepts Macrostrat database

concepts <- macrostrat_data("concepts.json", "https://macrostrat.org/api/defs/strat_name_concepts?all&format=json&response=long")

## output 1 - merge of strat and concepts
strat_concepts <- right_join(strat, concepts, by = c("concept_id", "concept_id")) # right join - on the basis strat packages with concept = 0 are captured below

# find entries without concept but can be included on the basis of strat_flag def

strat_flags <- subset(strat, grepl(paste(rocks, collapse = "|"), strat$strat_name_long, ignore.case=TRUE))

strat_flags <- subset(strat_flags, concept_id == 0)

# join for complete lith output

strat_concepts <- safe_full_join(strat_flags, strat_concepts, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)

# import units Macrostrat database

units <- macrostrat_data("units.json", "https://macrostrat.org/api/units?age_top=0&age_bottom=4540000000&response=long&format=json")

# merge the two Macrostrat databases into one - match using strat_names_ID

strat_concepts$strat_name_id  = as.character(strat_concepts$strat_name_id)
strat$strat_name_id  = as.character(strat$strat_name_id)
units$strat_name_id  = as.character(units$strat_name_id)

## output 2 - composite database where units are propagated and prioritised (where available)

# coalesce units and strat databases (units prioritised)

units_strat <- safe_right_join(units, strat, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)

units_strat_concepts <- right_join(concepts, units_strat, by = c("concept_id", "concept_id"))

# extract units

units2 <- units_strat_concepts[!is.na(units_strat_concepts$unit_id),] # carry forward all units (lith info is subsetted during later analysis)

# extract strat units

strat2 <- units_strat_concepts[is.na(units_strat_concepts$unit_id),]

strat_flags2 <- subset(strat2, grepl(paste(rocks, collapse = "|"), strat2$strat_name_long, ignore.case=TRUE))

strat_flags2 <- subset(strat_flags2, concept_id == 0) # this brings forward all strat packages with lith mentioned in the strat name (e.g., flag) but without a concept

strat2 <- subset(strat2, concept_id != 0) #| concept_id != NA) this part not necessary. This brings forward all strat packages with a concept

strat_units_concepts <- rbind(strat2, strat_flags2, units2) # combine

## next step is to join both dataframes (options 1 and 2) to the xDD extracts 

extracts$strat_name_id  = as.character(extracts$strat_name_id)
strat_concepts$strat_name_id  = as.character(strat_concepts$strat_name_id)

# generate output files, part1 excluding units, part 2 including units

data_part1 <- safe_right_join(extracts, strat_concepts, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)
data_part2 <- safe_right_join(extracts, strat_units_concepts, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)

##### PART 2 - Prepare dataset for normalisation to Macrostrat units #####
##### This section interpolates time bins between top and base ages ####
##### seperately defining phanerozoic & precambrian bin size ####

# for output 1

cutoff <- 541 # Ma
phanerozoic_increment <- 1 # in Ma
precambrian_increment <- 10 # in Ma

# add placeholder numbered bins to a dataframe, by increment, for later population
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

data_part1_precambrian <- subset(data_part1, t_age > cutoff)
data_part1_precambrian <- add_all_bins(data_part1_precambrian, precambrian_increment)

# including units spanning Phanerozoic-Precambrian boundary

data_part1_phanerozoic1 <- subset(data_part1, t_age < cutoff & b_age < cutoff)
data_part1_phanerozoic1 <- add_all_bins(data_part1_phanerozoic1, phanerozoic_increment)

data_part1_phanerozoic2 <- subset(data_part1, t_age < cutoff & b_age > cutoff)
data_part1_phanerozoic2$b_age <- cutoff
data_part1_phanerozoic2 <- add_all_bins(data_part1_phanerozoic2, phanerozoic_increment)

data_part1_phanerozoic3 <- subset(data_part1, t_age < cutoff & b_age > cutoff)
data_part1_phanerozoic3$t_age <- cutoff
data_part1_phanerozoic3 <- add_all_bins(data_part1_phanerozoic3, precambrian_increment)

# Bind the Phanerozoic and Precambrian datasets

data_part1 <- bind_rows(data_part1_phanerozoic1, data_part1_phanerozoic2, data_part1_phanerozoic3, data_part1_precambrian)

# remove any nested carriage returns, which can cause problems during export

data_part1 <- as.data.frame(sapply(data_part1, function(x) gsub("\n", "", x)))
data_part1 <- as.data.frame(sapply(data_part1, function(x) gsub("\r", "", x)))

# export output1 (csv works well here)

write.csv(data_part1, "data_part1_comp.csv")

# repeat again for output 2 (composite dataset)

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

data_part2_precambrian <- subset(data_part2, t_age > cutoff)
data_part2_precambrian <- add_all_bins(data_part2_precambrian, precambrian_increment)

# including units spanning Phanerozoic-Precambrian boundary

data_part2_phanerozoic1 <- subset(data_part2, t_age < cutoff & b_age < cutoff)
data_part2_phanerozoic1 <- add_all_bins(data_part2_phanerozoic1, phanerozoic_increment)

data_part2_phanerozoic2 <- subset(data_part2, t_age < cutoff & b_age > cutoff)
data_part2_phanerozoic2$b_age <- cutoff
data_part2_phanerozoic2 <- add_all_bins(data_part2_phanerozoic2, phanerozoic_increment)

data_part2_phanerozoic3 <- subset(data_part2, t_age < cutoff & b_age > cutoff)
data_part2_phanerozoic3$t_age <- cutoff
data_part2_phanerozoic3 <- add_all_bins(data_part2_phanerozoic3, precambrian_increment)

# Bind the Phanerozoic and Precambrian datasets

data_part2 <- bind_rows(data_part2_phanerozoic1, data_part2_phanerozoic2, data_part2_phanerozoic3, data_part2_precambrian)

#  convert the nested lists to character strings (otherwise can cause problems during export)

data_part2$lith  = as.character(data_part2$lith)
data_part2$environ  = as.character(data_part2$environ)
data_part2$econ  = as.character(data_part2$econ)
data_part2$measure = as.character(data_part2$measure)
data_part2$units_above = as.character(data_part2$units_above)
data_part2$units_below = as.character(data_part2$units_below)
data_part2$refs.x = as.character(data_part2$refs.x)
data_part2$refs.y = as.character(data_part2$refs.y)

# remove any nested carriage returns, which can cause problems during export
data_part2 <- as.data.frame(sapply(data_part2, function(x) gsub("\n", "", x)))
data_part2 <- as.data.frame(sapply(data_part2, function(x) gsub("\r", "", x)))

# export output 2 (again csv works well here)

write.csv(data_part2, "data_part2_comp.csv")

# collect additional data sources

zaffos_et_al = 'https://raw.githubusercontent.com/UW-Macrostrat/PNAS_201702297/master/FinalData/ContinuousTimeSeries.csv'
download.file(zaffos_et_al, 'Zaffos_et_al.txt', method="auto")


##### END ####
