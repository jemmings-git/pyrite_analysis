### xDD snippets extraction file written by J Emmings ##
## purpose: to manipulate Wilkin-framboid xDD snippets ##

setwd("D:/Pyrite_backup/pyrite_analysis-master/pyrite_analysis-master") ### set your working directory

library(jsonlite)
library(ggplot2)
library(tidyr)
library(readr)
library(reshape2)
library(safejoin)
library(pammtools)
library(gridExtra)
library(dplyr)

# download Wilkin-framboid snippets (all pages)
# https://xdd.wisc.edu/api/snippets?term=Wilkin,framboid&full_results=true&inclusive=true&clean&known_terms=stratigraphic_names

# bind snippets pages 1-6

#data <- fromJSON("snippets.json")
#data <- data[["success"]][["data"]]

#data2 <- fromJSON("snippets2.json")
#data2 <- data2[["success"]][["data"]]

#data3 <- fromJSON("snippets3.json")
#data3 <- data3[["success"]][["data"]]

#data4 <- fromJSON("snippets4.json")
#data4 <- data4[["success"]][["data"]]

#data5 <- fromJSON("snippets5.json")
#data5 <- data5[["success"]][["data"]]

#data6 <- fromJSON("snippets6.json")
#data6 <- data6[["success"]][["data"]]

#wilkin <- rbind(data, data2, data3, data4, data5, data6)

#write_json(wilkin, "wilkin_snippets.json")

# read file

wilkin <- fromJSON("wilkin_snippets.json")

out <- hoist(wilkin, known_terms, "stratigraphic_names")

out <- subset(out, stratigraphic_names != "NULL")

out <- out[,-11]

out$stratigraphic_names <- gsub("list(), ", "", out$stratigraphic_names, fixed = TRUE)
out$stratigraphic_names <- gsub("list(", "", out$stratigraphic_names, fixed = TRUE)
out$stratigraphic_names <- gsub("c(", "", out$stratigraphic_names, fixed = TRUE)
out$stratigraphic_names <- gsub(")", "", out$stratigraphic_names, fixed = TRUE)
out$stratigraphic_names <- gsub("\"", "", out$stratigraphic_names)

test <- strsplit(out$stratigraphic_names, ",")
test <- do.call("rbind", strsplit(out$stratigraphic_names, ","))
test <- t(apply(test, 1, function(x) replace(x, duplicated(x), NA)))

out <- cbind(out, test)

stratigraphic_names <- melt(out, id.vars = c(1:10))

stratigraphic_names <- stratigraphic_names[!is.na(stratigraphic_names$value),]
stratigraphic_names <- subset(stratigraphic_names, value != " ")
stratigraphic_names$value <- trimws(stratigraphic_names$value, which = "left")

stratigraphic_names <- stratigraphic_names %>%
  group_by(`_gddid`) %>%
  mutate(count = n())

# remove duplicate strat names

stratigraphic_names <- stratigraphic_names[!duplicated(stratigraphic_names$value),]

# import macrostrat 

source('macrostrat_data.R')

rocks <- meta_sedimentary_rocks

strat <- macrostrat_data("strat.json", "https://macrostrat.org/api/defs/strat_names?all&format=json&response=long")

# remove any non-unique strat names

strat <- strat %>%
  group_by(strat_name_long) %>%
  filter(n() == 1)

concepts <- macrostrat_data("concepts.json", "https://macrostrat.org/api/defs/strat_name_concepts?all&format=json&response=long")

## output 1 - merge of strat and concepts 
strat_concepts <- right_join(strat, concepts, by = "concept_id") # right join - on the basis strat packages with concept = 0 are captured below

# find entries without concept but can be included on the basis of strat_flag def

strat_flags <- subset(strat, grepl(paste(rocks, collapse = "|"), strat$strat_name_long, ignore.case=TRUE))

strat_flags <- subset(strat_flags, concept_id == 0)

# join for complete lith output

strat_concepts <- safe_full_join(strat_flags, strat_concepts, by = "strat_name_id", conflict = coalesce)

# import units Macrostrat database

units <- macrostrat_data("units.json", "https://macrostrat.org/api/units?age_top=0&age_bottom=4540000000&response=long&format=json")

# merge the two Macrostrat databases into one - match using strat_names_ID

strat_concepts$strat_name_id  = as.character(strat_concepts$strat_name_id)
strat$strat_name_id  = as.character(strat$strat_name_id)
units$strat_name_id  = as.character(units$strat_name_id)

## output 2 - composite database where units are propagated and prioritised (where available)

# coalesce units and strat databases (units prioritised)

units_strat <- safe_right_join(units, strat, by = "strat_name_id", conflict = coalesce)

units_strat_concepts <- right_join(concepts, units_strat, by = "concept_id")

# extract units

units2 <- units_strat_concepts[!is.na(units_strat_concepts$unit_id),] # carry forward all units (lith info is subsetted during later analysis)

# extract strat units

strat2 <- units_strat_concepts[is.na(units_strat_concepts$unit_id),]

strat_flags2 <- subset(strat2, grepl(paste(rocks, collapse = "|"), strat2$strat_name_long, ignore.case=TRUE))

strat_flags2 <- subset(strat_flags2, concept_id == 0) # this brings forward all strat packages with lith mentioned in the strat name (e.g., flag) but without a concept

strat2 <- subset(strat2, concept_id != 0) #| concept_id != NA) this part not necessary. This brings forward all strat packages with a concept

strat_units_concepts <- rbind(strat2, strat_flags2, units2) # combine

## next step is to join both dataframes (options 1 and 2) to the xDD extracts 

# generate output files, part1 excluding units, part 2 including units

stratigraphic_names <- as.data.frame(stratigraphic_names)
names(stratigraphic_names)[12] <- "strat_name_long"

data_part1 <- left_join(stratigraphic_names, strat_concepts, by = "strat_name_long", keep = TRUE)
data_part1 <- data_part1[!duplicated(data_part1$strat_name_long.x),]
data_part1 <- data_part1[!is.na(data_part1$t_age),]

data_part2 <- left_join(stratigraphic_names, strat_units_concepts, by = "strat_name_long", keep = TRUE)
data_part2 <- data_part2[!is.na(data_part2$t_age),]

#### inflate ###

cutoff <- 541 # Ma 
phanerozoic_increment <- 1 # in Ma
precambrian_increment <- 10 # in Ma

# precambrian

data_part1_precambrian <- subset(data_part1, t_age > cutoff)
data_part1_precambrian$t_age2 <- round(data_part1_precambrian$t_age, digits = -1)+(cutoff-round(cutoff, digits = -1))
data_part1_precambrian$b_age2 <- round(data_part1_precambrian$b_age, digits = -1)+(cutoff-round(cutoff, digits = -1))

# phanerozoic & overlapping units

data_part1_phanerozoic1 <- subset(data_part1, t_age < cutoff & b_age < cutoff)
data_part1_phanerozoic1$t_age2 <- round(data_part1_phanerozoic1$t_age, digits = 0)
data_part1_phanerozoic1$b_age2 <- round(data_part1_phanerozoic1$b_age, digits = 0)


data_part1_phanerozoic2 <- subset(data_part1, t_age < cutoff & b_age > cutoff)
data_part1_phanerozoic2$b_age <- cutoff
data_part1_phanerozoic2$t_age2 <- round(data_part1_phanerozoic2$t_age, digits = 0)
data_part1_phanerozoic2$b_age2 <- round(data_part1_phanerozoic2$b_age, digits = 0)


data_part1_phanerozoic3 <- subset(data_part1, t_age < cutoff & b_age > cutoff)
data_part1_phanerozoic3$t_age <- cutoff
data_part1_phanerozoic3$t_age2 <- round(data_part1_phanerozoic3$t_age, digits = 0)
data_part1_phanerozoic3$b_age2 <- round(data_part1_phanerozoic3$b_age, digits = -1)+(cutoff-round(cutoff, digits = -1))

# Bind the Phanerozoic and Precambrian datasets

data_part1 <- bind_rows(data_part1_phanerozoic1, data_part1_phanerozoic2, data_part1_phanerozoic3, data_part1_precambrian)

# then repeat for part2

# precambrian

data_part2_precambrian <- subset(data_part2, t_age > cutoff)
data_part2_precambrian$t_age2 <- round(data_part2_precambrian$t_age, digits = -1)+(cutoff-round(cutoff, digits = -1))
data_part2_precambrian$b_age2 <- round(data_part2_precambrian$b_age, digits = -1)+(cutoff-round(cutoff, digits = -1))

# phanerozoic & overlapping units

data_part2_phanerozoic1 <- subset(data_part2, t_age < cutoff & b_age < cutoff)
data_part2_phanerozoic1$t_age2 <- round(data_part2_phanerozoic1$t_age, digits = 0)
data_part2_phanerozoic1$b_age2 <- round(data_part2_phanerozoic1$b_age, digits = 0)


data_part2_phanerozoic2 <- subset(data_part2, t_age < cutoff & b_age > cutoff)
data_part2_phanerozoic2$b_age <- cutoff
data_part2_phanerozoic2$t_age2 <- round(data_part2_phanerozoic2$t_age, digits = 0)
data_part2_phanerozoic2$b_age2 <- round(data_part2_phanerozoic2$b_age, digits = 0)


data_part2_phanerozoic3 <- subset(data_part2, t_age < cutoff & b_age > cutoff)
data_part2_phanerozoic3$t_age <- cutoff
data_part2_phanerozoic3$t_age2 <- round(data_part2_phanerozoic3$t_age, digits = 0)
data_part2_phanerozoic3$b_age2 <- round(data_part2_phanerozoic3$b_age, digits = -1)+(cutoff-round(cutoff, digits = -1))

# Bind

data_part2 <- bind_rows(data_part2_phanerozoic1, data_part2_phanerozoic2, data_part2_phanerozoic3, data_part2_precambrian)


# for output 1


add_all_bins <- function(df, increment) {
  ncols <- ceiling((max(df$b_age2 - df$t_age2))/increment)
  
  df <- df %>% mutate(bin1 = ifelse(b_age2-increment > t_age2+increment, t_age2+increment, NA))
  for(n in 2:(ncols-1)){
    next_bin <- paste('bin', n, sep='')
    bin <- paste0('bin', (n-1), sep='')
    # See https://stackoverflow.com/a/49311813/4319767 for why this syntax 
    df <- df %>% mutate(!!next_bin := ifelse(UQ(rlang::sym(bin)) < b_age2-increment, UQ(rlang::sym(bin))+increment, NA))
  }  
  return(df)
}

data_part1_precambrian <- subset(data_part1, t_age2 > cutoff)
data_part1_precambrian <- add_all_bins(data_part1_precambrian, precambrian_increment)

# including units spanning Phanerozoic-Precambrian boundary

data_part1_phanerozoic1 <- subset(data_part1, t_age2 < cutoff & b_age2 < cutoff)
data_part1_phanerozoic1 <- add_all_bins(data_part1_phanerozoic1, phanerozoic_increment)

data_part1_phanerozoic2 <- subset(data_part1, t_age2 < cutoff & b_age2 == cutoff)
data_part1_phanerozoic2 <- add_all_bins(data_part1_phanerozoic2, phanerozoic_increment)

data_part1_phanerozoic3 <- subset(data_part1, t_age2 == cutoff & b_age2 > cutoff)
data_part1_phanerozoic3 <- add_all_bins(data_part1_phanerozoic3, precambrian_increment)

# Bind the Phanerozoic and Precambrian datasets

data_part1 <- bind_rows(data_part1_phanerozoic1, data_part1_phanerozoic2, data_part1_phanerozoic3, data_part1_precambrian)

# remove any nested carriage returns, which can cause problems during export

data_part1 <- as.data.frame(sapply(data_part1, function(x) gsub("\n", "", x)))
data_part1 <- as.data.frame(sapply(data_part1, function(x) gsub("\r", "", x)))

# export output1 (csv works well here)

write.csv(data_part1, "wilkin_framboids_strat.csv")

# repeat again for output 2 (composite dataset)

cutoff <- 541 # Ma 
phanerozoic_increment <- 1 # in Ma
precambrian_increment <- 10 # in Ma


add_all_bins <- function(df, increment) {
  ncols <- ceiling((max(df$b_age2 - df$t_age2))/increment)
  
  df <- df %>% mutate(bin1 = ifelse(b_age2-increment > t_age2+increment, t_age2+increment, NA))
  for(n in 2:(ncols-1)){
    next_bin <- paste('bin', n, sep='')
    bin <- paste0('bin', (n-1), sep='')
    # See https://stackoverflow.com/a/49311813/4319767 for why this syntax 
    df <- df %>% mutate(!!next_bin := ifelse(UQ(rlang::sym(bin)) < b_age2-increment, UQ(rlang::sym(bin))+increment, NA))
  }  
  return(df)
}

data_part2_precambrian <- subset(data_part2, t_age2 > cutoff)
data_part2_precambrian <- add_all_bins(data_part2_precambrian, precambrian_increment)

# including units spanning Phanerozoic-Precambrian boundary

data_part2_phanerozoic1 <- subset(data_part2, t_age2 < cutoff & b_age2 < cutoff)
data_part2_phanerozoic1 <- add_all_bins(data_part2_phanerozoic1, phanerozoic_increment)

data_part2_phanerozoic2 <- subset(data_part2, t_age2 < cutoff & b_age2 == cutoff)
data_part2_phanerozoic2 <- add_all_bins(data_part2_phanerozoic2, phanerozoic_increment)

data_part2_phanerozoic3 <- subset(data_part2, t_age2 == cutoff & b_age2 > cutoff)
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

write.csv(data_part2, "wilkin_framboids_comp.csv")
