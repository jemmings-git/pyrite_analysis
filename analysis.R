##### This R script was written by Joe Emmings, Jo Walsh and Kathryn Leeming (British Geological Survey) ####
##### This script is designed to plot GeoDeepDive extractions
##### interfaced with the Macrostrat database


rm(list=ls())

library(dplyr) # plyr is not needed (note loading plyr after dplyr will prevent execution of PART 2)
library(tidyr)
library(readr)
library(jsonlite)
library(reshape2)
library(ggplot2)
library(safejoin)
library(pammtools)
library(gridExtra)
library(tidyverse)
library(deeptime)
library(ggfittext)
library(grid)

project_home <- '...' # insert directory
tryCatch({
  setwd(project_home)
}, error = function(err) { 
  if (dir.exists('./data')) {
    setwd('./data') }
}
)


# choices for analysis:

source('macrostrat_data.R')

# choose one of the following rocks and rocks2 defs

#rocks <- sedimentary_rocks # normalisation to all sedimentary rocks
#rocks2 <- sedimentary_rocks

rocks <- meta_sedimentary_rocks # normalisation to all sedimentary and metamorphosed sedimentary rocks - choose this option to replicate manuscript output
rocks2 <- meta_sedimentary_rocks

#rocks <- mud_rocks # normalisation to mudstones
#rocks2 <- mud_rocks

###### PART 1 - using 'output 2' ('p2') - composite analysis of units and strat packages #####

data_p2 <-  read.csv(file = "data_part2_comp.csv", row.names = 1)

length(unique(data_p2$result_id))


# generate complete list of sedimentary units for part 2 - 
# this is done by propagating unit lith info to strat IDs
# so, for any given strat package, if the majority of units within that package
# are 'sedimentary', then the overall package is considered sedimentary
# this info is then joined to the other sediment definitions derived from concepts descriptions

# run this block below - this is required for normalisation to all sediments
# Further subdivision, such as 'marine', etc. can be implemented by via 'rocks' above

data_p2$lith2 <- (grepl(paste(rocks, collapse = "|"), data_p2$lith, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), data_p2$other, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), data_p2$strat_phrase_root, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), data_p2$strat_flag, ignore.case=TRUE) |
  grepl(paste(rocks, collapse = "|"), data_p2$phrase, ignore.case=TRUE) |
  grepl(paste(rocks, collapse = "|"), data_p2$environ, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), data_p2$lith, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), data_p2$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), data_p2$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), data_p2$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), data_p2$phrase, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), data_p2$environ, ignore.case=TRUE))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# keep this file for part 2

sed.list <- aggregate(lith2 ~ strat_name_id, data_p2, Mode)

# block end

# optional compilation of target words and summary
phrases <- as.data.frame(as.array(data_p2$target_word))
list.phrases <- as.data.frame(summary(data_p2$target_word))

# compile framboid mentions

# block start

# includes search in phrase since 'framboid' is sufficiently unique

framboids <- grepl("framboid", data_p2$target_word, ignore.case=TRUE) |
  grepl("framboid", data_p2$phrase, ignore.case=TRUE)

framboids <- subset(data_p2, framboids)

length(unique(framboids$result_id)) # no of framboid results
length(unique(framboids$strat_name_id)) # no of strat packages
length(unique(framboids$unit_id)) # no of units

# include any framboid mentions in concept description (not likely)

framboids1 <- grepl("framboid", data_p2$other, ignore.case=TRUE)
framboids1 <- subset(data_p2, framboids1)
framboids <- rbind(framboids, framboids1)

# note some duplicate results_id should not be removed - these are present for strat packages spanning the Precambrian-Phanerozoic, e.g., Tindir Gp
# therefore must subdivide the dataset before removing duplicates associated with multiple phrase hits

# the block below cleans the framboid dataframe, removing genuine duplicate hits
# and limiting unit matches to sedimentary units (as defined in lith)

cutoff <- 541 # Ma 
phanerozoic_increment <- 1 # in Ma
precambrian_increment <- 10 # in Ma

framboids1 <- subset(framboids, t_age2 > cutoff & b_age2 > cutoff)
framboids2 <- framboids1[!duplicated(framboids1$unit_id), ]
framboids3 <- framboids1[is.na(framboids1$unit_id),]
framboids2 <- framboids2[!is.na(framboids2$unit_id),] 
framboids3 <- framboids3[!duplicated(framboids3$strat_name_id),]
framboids1 <- rbind(framboids2, framboids3)

framboids2 <- subset(framboids, t_age2 < cutoff & b_age2 < cutoff)
framboids3 <- framboids2[!duplicated(framboids2$unit_id), ]
framboids4 <- framboids2[is.na(framboids2$unit_id),]
framboids3 <- framboids3[!is.na(framboids3$unit_id),]
framboids4 <- framboids4[!duplicated(framboids4$strat_name_id),]
framboids2 <- rbind(framboids3, framboids4)

framboids3 <- subset(framboids, t_age2 < cutoff & b_age2 == cutoff)
framboids4 <- framboids3[!duplicated(framboids3$unit_id), ]
framboids5 <- framboids3[is.na(framboids3$unit_id),]
framboids4 <- framboids4[!is.na(framboids4$unit_id),]
framboids5 <- framboids5[!duplicated(framboids5$strat_name_id),]
framboids3 <- rbind(framboids4, framboids5)

framboids4 <- subset(framboids, t_age2 == cutoff & b_age2 > cutoff)
framboids5 <- framboids4[!duplicated(framboids4$unit_id), ]
framboids6 <- framboids4[is.na(framboids4$unit_id),]
framboids5 <- framboids5[!is.na(framboids5$unit_id),]
framboids6 <- framboids6[!duplicated(framboids6$strat_name_id),]
framboids4 <- rbind(framboids5, framboids6)

framboids <- rbind(framboids1, framboids2, framboids3, framboids4)

sediments <- (grepl(paste(rocks, collapse = "|"), framboids$lith, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), framboids$other, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), framboids$environ, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), framboids$strat_phrase_root, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), framboids$phrase, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), framboids$strat_flag, ignore.case=TRUE) |
  grepl(paste(rocks, collapse = "|"), framboids$environ, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), framboids$lith, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), framboids$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), framboids$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), framboids$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), framboids$phrase, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), framboids$environ, ignore.case=TRUE))

framboids <- subset(framboids, sediments) # subset to 'rocks' of interest

# block end

# repeat the same process for any other 'similar' phrases - in this case 'pyrite nodules'

# block start

nodules <- grepl("nodul", data_p2$target_word, ignore.case=TRUE) | grepl("concretion", data_p2$target_word, ignore.case=TRUE) 

nodules <- subset(data_p2, nodules)

length(unique(nodules$result_id)) # no. of nodule+concretion results
length(unique(nodules$strat_name_id)) # no of strat packages
length(unique(nodules$unit_id)) # no of units

# include nodule mentions pulled from lexicons (concepts)

nodules1 <- grepl("pyrite nodul", data_p2$other, ignore.case=TRUE) | grepl("pyrite concretion", data_p2$other, ignore.case=TRUE) |
  grepl("pyritic nodul", data_p2$other, ignore.case=TRUE) | grepl("pyritic concretion", data_p2$other, ignore.case=TRUE) |
  grepl("nodules of pyrite", data_p2$other, ignore.case=TRUE) | grepl("concretions of pyrite", data_p2$other, ignore.case=TRUE) | 
  grepl("nodular pyrite", data_p2$other, ignore.case=TRUE) | grepl("concretionary pyrite", data_p2$other, ignore.case=TRUE)

nodules1 <- subset(data_p2, nodules1)
nodules <- rbind(nodules, nodules1)

nodules1 <- subset(nodules, t_age2 > cutoff & b_age2 > cutoff)
nodules2 <- nodules1[!duplicated(nodules1$unit_id), ]
nodules3 <- nodules1[is.na(nodules1$unit_id),]
nodules2 <- nodules2[!is.na(nodules2$unit_id),]
nodules3 <- nodules3[!duplicated(nodules3$strat_name_id),]
nodules1 <- rbind(nodules2, nodules3)

nodules2 <- subset(nodules, t_age2 < cutoff & b_age2 < cutoff)
nodules3 <- nodules2[!duplicated(nodules2$unit_id), ]
nodules4 <- nodules2[is.na(nodules2$unit_id),]
nodules3 <- nodules3[!is.na(nodules3$unit_id),]
nodules4 <- nodules4[!duplicated(nodules4$strat_name_id),]
nodules2 <- rbind(nodules3, nodules4)

nodules3 <- subset(nodules, t_age2 < cutoff & b_age2 == cutoff)
nodules4 <- nodules3[!duplicated(nodules3$unit_id), ]
nodules5 <- nodules3[is.na(nodules3$unit_id),]
nodules4 <- nodules4[!is.na(nodules4$unit_id),]
nodules5 <- nodules5[!duplicated(nodules5$strat_name_id),]
nodules3 <- rbind(nodules4, nodules5)

nodules4 <- subset(nodules, t_age2 == cutoff & b_age2 > cutoff)
nodules5 <- nodules4[!duplicated(nodules4$unit_id), ]
nodules6 <- nodules4[is.na(nodules4$unit_id),]
nodules5 <- nodules5[!is.na(nodules5$unit_id),]
nodules6 <- nodules6[!duplicated(nodules6$strat_name_id),]
nodules4 <- rbind(nodules5, nodules6)

nodules <- rbind(nodules1, nodules2, nodules3, nodules4)

sediments <- (grepl(paste(rocks, collapse = "|"), nodules$lith, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), nodules$other, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), nodules$environ, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), nodules$strat_phrase_root, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), nodules$phrase, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), nodules$strat_flag, ignore.case=TRUE) |
  grepl(paste(rocks, collapse = "|"), nodules$environ, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), nodules$lith, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), nodules$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), nodules$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), nodules$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), nodules$phrase, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), nodules$environ, ignore.case=TRUE))

nodules <- subset(nodules, sediments) # subset to 'rocks' of interest

# block end

# the next section processes undifferentiated pyrite mentions - this is handled slightly differently

# block start

# start by including pyritic strat concepts pulled from lexicons (undifferentiated)

pyritic_strat <- grepl("pyrit", data_p2$other, ignore.case=TRUE)  |  grepl("pyrit", data_p2$lith, ignore.case=TRUE)

pyritic_strat <- subset(data_p2, pyritic_strat)
false <- !grepl("non-pyrit", pyritic_strat$other, ignore.case=TRUE) & !grepl("non-pyrit", pyritic_strat$lith, ignore.case=TRUE)
pyritic_strat <- subset(pyritic_strat, false)

# now need to carry forward only units not already listed in framboids or nodules record

# therefore first bind framboid and nodules datasets

all0 <- rbind(framboids, nodules)

all1 <- all0[!duplicated(all0$unit_id), ]
all2 <- all0[is.na(all0$unit_id),]
all1 <- all1[!is.na(all1$unit_id),]
all2 <- all2[!duplicated(all2$strat_name_id),]
all0 <- rbind(all1, all2)

# right join framboids+nodules to the pyritic strat record

pyritic_strat <- safe_right_join(all0, pyritic_strat, by = "strat_name_id", conflict = coalesce)
pyritic_strat <- pyritic_strat[is.na(pyritic_strat$target_word),]

pyrite_undif <- grepl("\\<pyrite\\>", data_p2$target_word, ignore.case=TRUE) | grepl("\\<pyritic\\>", data_p2$target_word, ignore.case=TRUE)

pyrite_undif <- subset(data_p2, pyrite_undif)

pyrite_undif <- rbind(pyritic_strat, pyrite_undif)

# pyrite mineralisation and veins searches - including framboid- and nodule-bearing rocks

pyrite_all <- pyrite_undif

pyrite_all1 <- subset(pyrite_all, t_age2 > cutoff & b_age2 > cutoff)
pyrite_all2 <- pyrite_all1[!duplicated(pyrite_all1$unit_id), ]
pyrite_all3 <- pyrite_all1[is.na(pyrite_all1$unit_id),]
pyrite_all2 <- pyrite_all2[!is.na(pyrite_all2$unit_id),]
pyrite_all3 <- pyrite_all3[!duplicated(pyrite_all3$strat_name_id),]
pyrite_all1 <- rbind(pyrite_all2, pyrite_all3)

pyrite_all2 <- subset(pyrite_all, t_age2 < cutoff & b_age2 < cutoff)
pyrite_all3 <- pyrite_all2[!duplicated(pyrite_all2$unit_id), ]
pyrite_all4 <- pyrite_all2[is.na(pyrite_all2$unit_id),]
pyrite_all3 <- pyrite_all3[!is.na(pyrite_all3$unit_id),]
pyrite_all4 <- pyrite_all4[!duplicated(pyrite_all4$strat_name_id),]
pyrite_all2 <- rbind(pyrite_all3, pyrite_all4)

pyrite_all3 <- subset(pyrite_all, t_age2 < cutoff & b_age2 == cutoff)
pyrite_all4 <- pyrite_all3[!duplicated(pyrite_all3$unit_id), ]
pyrite_all5 <- pyrite_all3[is.na(pyrite_all3$unit_id),]
pyrite_all4 <- pyrite_all4[!is.na(pyrite_all4$unit_id),]
pyrite_all5 <- pyrite_all5[!duplicated(pyrite_all5$strat_name_id),]
pyrite_all3 <- rbind(pyrite_all4, pyrite_all5)

pyrite_all4 <- subset(pyrite_all, t_age2 == cutoff & b_age2 > cutoff)
pyrite_all5 <- pyrite_all4[!duplicated(pyrite_all4$unit_id), ]
pyrite_all6 <- pyrite_all4[is.na(pyrite_all4$unit_id),]
pyrite_all5 <- pyrite_all5[!is.na(pyrite_all5$unit_id),]
pyrite_all6 <- pyrite_all6[!duplicated(pyrite_all6$strat_name_id),]
pyrite_all4 <- rbind(pyrite_all5, pyrite_all6)

pyrite_all <- rbind(pyrite_all1, pyrite_all2, pyrite_all3, pyrite_all4)

sediments <- (grepl(paste(rocks, collapse = "|"), pyrite_all$lith, ignore.case=TRUE) | 
                grepl(paste(rocks, collapse = "|"), pyrite_all$other, ignore.case=TRUE) | 
                grepl(paste(rocks, collapse = "|"), pyrite_all$environ, ignore.case=TRUE) | 
                grepl(paste(rocks, collapse = "|"), pyrite_all$strat_phrase_root, ignore.case=TRUE) | 
                grepl(paste(rocks, collapse = "|"), pyrite_all$phrase, ignore.case=TRUE) | 
                grepl(paste(rocks, collapse = "|"), pyrite_all$strat_flag, ignore.case=TRUE) |
                grepl(paste(rocks, collapse = "|"), pyrite_all$environ, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), pyrite_all$lith, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), pyrite_all$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), pyrite_all$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), pyrite_all$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), pyrite_all$phrase, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), pyrite_all$environ, ignore.case=TRUE))

pyrite_all <- subset(pyrite_all, sediments)

# extract nearby mentions of veins or mineralisation

veins2 <- grepl("vein", pyrite_all$phrase, ignore.case=TRUE)  |  grepl("mineralisation", pyrite_all$phrase, ignore.case=TRUE)

veins2 <- subset(pyrite_all, veins2)

# use anti_join to again carry forward only packages not already identified as containing framboids or nodules

pyrite_undif <- anti_join(pyrite_undif, all0, by = "strat_name_id") 

pyrite_undif1 <- subset(pyrite_undif, t_age2 > cutoff & b_age2 > cutoff)
pyrite_undif2 <- pyrite_undif1[!duplicated(pyrite_undif1$unit_id), ]
pyrite_undif3 <- pyrite_undif1[is.na(pyrite_undif1$unit_id),]
pyrite_undif2 <- pyrite_undif2[!is.na(pyrite_undif2$unit_id),]
pyrite_undif3 <- pyrite_undif3[!duplicated(pyrite_undif3$strat_name_id),]
pyrite_undif1 <- rbind(pyrite_undif2, pyrite_undif3)

pyrite_undif2 <- subset(pyrite_undif, t_age2 < cutoff & b_age2 < cutoff)
pyrite_undif3 <- pyrite_undif2[!duplicated(pyrite_undif2$unit_id), ]
pyrite_undif4 <- pyrite_undif2[is.na(pyrite_undif2$unit_id),]
pyrite_undif3 <- pyrite_undif3[!is.na(pyrite_undif3$unit_id),]
pyrite_undif4 <- pyrite_undif4[!duplicated(pyrite_undif4$strat_name_id),]
pyrite_undif2 <- rbind(pyrite_undif3, pyrite_undif4)

pyrite_undif3 <- subset(pyrite_undif, t_age2 < cutoff & b_age2 == cutoff)
pyrite_undif4 <- pyrite_undif3[!duplicated(pyrite_undif3$unit_id), ]
pyrite_undif5 <- pyrite_undif3[is.na(pyrite_undif3$unit_id),]
pyrite_undif4 <- pyrite_undif4[!is.na(pyrite_undif4$unit_id),]
pyrite_undif5 <- pyrite_undif5[!duplicated(pyrite_undif5$strat_name_id),]
pyrite_undif3 <- rbind(pyrite_undif4, pyrite_undif5)

pyrite_undif4 <- subset(pyrite_undif, t_age2 == cutoff & b_age2 > cutoff)
pyrite_undif5 <- pyrite_undif4[!duplicated(pyrite_undif4$unit_id), ]
pyrite_undif6 <- pyrite_undif4[is.na(pyrite_undif4$unit_id),]
pyrite_undif5 <- pyrite_undif5[!is.na(pyrite_undif5$unit_id),]
pyrite_undif6 <- pyrite_undif6[!duplicated(pyrite_undif6$strat_name_id),]
pyrite_undif4 <- rbind(pyrite_undif5, pyrite_undif6)

pyrite_undif <- rbind(pyrite_undif1, pyrite_undif2, pyrite_undif3, pyrite_undif4)

sediments <- (grepl(paste(rocks, collapse = "|"), pyrite_undif$lith, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), pyrite_undif$other, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), pyrite_undif$environ, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), pyrite_undif$strat_phrase_root, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), pyrite_undif$phrase, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), pyrite_undif$strat_flag, ignore.case=TRUE) |
  grepl(paste(rocks, collapse = "|"), pyrite_undif$environ, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), pyrite_undif$lith, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), pyrite_undif$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), pyrite_undif$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), pyrite_undif$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), pyrite_undif$phrase, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), pyrite_undif$environ, ignore.case=TRUE))

pyrite_undif <- subset(pyrite_undif, sediments)

# extract nearby mentions of veins or mineralisation - not including framboid- or nodule- bearing rocks

veins <- grepl("vein", pyrite_undif$phrase, ignore.case=TRUE)  |  grepl("mineralisation", pyrite_undif$phrase, ignore.case=TRUE)

veins <- subset(pyrite_undif, veins)

# block end

# the following section melts all datasets (i.e., convert wide to long format)

# block start

# for example framboids

framboids <- melt(framboids, id.vars = c(1:83,626), na.rm=TRUE)
framboids$value <- as.numeric(framboids$value)
framboids <- subset(framboids, variable != "b_age2")

nodules <- melt(nodules, id.vars = c(1:83,626), na.rm=TRUE)
nodules$value <- as.numeric(nodules$value)
nodules <- subset(nodules, variable != "b_age2")

pyrite_undif <- melt(pyrite_undif, id.vars = c(1:83,626), na.rm=TRUE)
pyrite_undif$value <- as.numeric(pyrite_undif$value)
pyrite_undif <- subset(pyrite_undif, variable != "b_age2")

veins <- melt(veins, id.vars = c(1:83,626), na.rm=TRUE)
veins$value <- as.numeric(veins$value)
veins <- subset(veins, variable != "b_age2")

veins2 <- melt(veins2, id.vars = c(1:83,626), na.rm=TRUE)
veins2$value <- as.numeric(veins2$value)
veins2 <- subset(veins2, variable != "b_age2")

# block end

# bin the datasets

# block start

phanerozoic_increment <- 1 # in Ma
precambrian_increment <- 10 # in MA

# for example framboids

# Phanerozoic bins

framboids_bins1 <- hist(framboids$value[framboids$value >= 0 & framboids$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
framboids_bins1 <- as.data.frame(cbind(framboids_bins1$counts, framboids_bins1$breaks))
framboids_bins1 <- framboids_bins1[1:542,]

# Precambrian bins
framboids_bins2 <- hist(framboids$value[framboids$value >= cutoff & framboids$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
framboids_bins2 <- as.data.frame(cbind(framboids_bins2$counts, framboids_bins2$breaks))

# Bind
framboids_bins <- rbind(framboids_bins1,framboids_bins2)

framboids_bins <- framboids_bins[c(1:541,543:889),] # remove duplicate 541

# repeat for the others

# nodules

nodules_bins1 <- hist(nodules$value[nodules$value >= 0 & nodules$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
nodules_bins1 <- as.data.frame(cbind(nodules_bins1$counts, nodules_bins1$breaks))
nodules_bins1 <- nodules_bins1[1:542,]

nodules_bins2 <- hist(nodules$value[nodules$value >= cutoff & nodules$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
nodules_bins2 <- as.data.frame(cbind(nodules_bins2$counts, nodules_bins2$breaks))

nodules_bins <- rbind(nodules_bins1,nodules_bins2)
nodules_bins <- nodules_bins[c(1:541,543:889),] # remove duplicate 541

# undifferentiated pyrite mentions

pyrite_undif_bins1 <- hist(pyrite_undif$value[pyrite_undif$value >= 0 & pyrite_undif$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
pyrite_undif_bins1 <- as.data.frame(cbind(pyrite_undif_bins1$counts, pyrite_undif_bins1$breaks))
pyrite_undif_bins1 <- pyrite_undif_bins1[1:542,]

pyrite_undif_bins2 <- hist(pyrite_undif$value[pyrite_undif$value >= cutoff & pyrite_undif$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
pyrite_undif_bins2 <- as.data.frame(cbind(pyrite_undif_bins2$counts, pyrite_undif_bins2$breaks))

pyrite_undif_bins <- rbind(pyrite_undif_bins1,pyrite_undif_bins2)
pyrite_undif_bins <- pyrite_undif_bins[c(1:541,543:889),] # remove duplicate 541

# undifferentiated pyrite mentions - associated with veins and mineralisation mentions

veins_bins1 <- hist(veins$value[veins$value >= 0 & veins$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
veins_bins1 <- as.data.frame(cbind(veins_bins1$counts, veins_bins1$breaks))
veins_bins1 <- veins_bins1[1:542,]

veins_bins2 <- hist(veins$value[veins$value >= cutoff & veins$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
veins_bins2 <- as.data.frame(cbind(veins_bins2$counts, veins_bins2$breaks))

veins_bins <- rbind(veins_bins1,veins_bins2)
veins_bins <- veins_bins[c(1:541,543:889),] # remove duplicate 541

veins2_bins1 <- hist(veins2$value[veins2$value >= 0 & veins2$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
veins2_bins1 <- as.data.frame(cbind(veins2_bins1$counts, veins2_bins1$breaks))
veins2_bins1 <- veins2_bins1[1:542,]

veins2_bins2 <- hist(veins2$value[veins2$value >= cutoff & veins2$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
veins2_bins2 <- as.data.frame(cbind(veins2_bins2$counts, veins2_bins2$breaks))

veins2_bins <- rbind(veins2_bins1,veins2_bins2)
veins2_bins <- veins2_bins[c(1:541,543:889),] # remove duplicate 541

# block end

# now to normalise to sedimentary packages - including nodules and framboids if not identified as 'sedimentary'

# block start

# run this block for sediment normalisation

sediments <- (grepl(paste(rocks, collapse = "|"), data_p2$lith, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), data_p2$other, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), data_p2$environ, ignore.case=TRUE) |
  grepl(paste(rocks, collapse = "|"), data_p2$strat_phrase_root, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), data_p2$phrase, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), data_p2$strat_flag, ignore.case=TRUE) |
  grepl(paste(rocks, collapse = "|"), data_p2$environ, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), data_p2$lith, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), data_p2$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), data_p2$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), data_p2$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), data_p2$phrase, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), data_p2$environ, ignore.case=TRUE))

sediments <- subset(data_p2, sediments)

sediments1 <- subset(sediments, t_age2 > cutoff & b_age2 > cutoff)
sediments2 <- sediments1[!duplicated(sediments1$unit_id), ]
sediments3 <- sediments1[is.na(sediments1$unit_id),]
sediments2 <- sediments2[!is.na(sediments2$unit_id),]
sediments3 <- sediments3[!duplicated(sediments3$strat_name_id),]
sediments1 <- rbind(sediments2, sediments3)

sediments2 <- subset(sediments, t_age2 < cutoff & b_age2 < cutoff)
sediments3 <- sediments2[!duplicated(sediments2$unit_id), ]
sediments4 <- sediments2[is.na(sediments2$unit_id),]
sediments3 <- sediments3[!is.na(sediments3$unit_id),]
sediments4 <- sediments4[!duplicated(sediments4$strat_name_id),]
sediments2 <- rbind(sediments3, sediments4)

sediments3 <- subset(sediments, t_age2 < cutoff & b_age2 == cutoff)
sediments4 <- sediments3[!duplicated(sediments3$unit_id), ]
sediments5 <- sediments3[is.na(sediments3$unit_id),]
sediments4 <- sediments4[!is.na(sediments4$unit_id),]
sediments5 <- sediments5[!duplicated(sediments5$strat_name_id),]
sediments3 <- rbind(sediments4, sediments5)

sediments4 <- subset(sediments, t_age2 == cutoff & b_age2 > cutoff)
sediments5 <- sediments4[!duplicated(sediments4$unit_id), ]
sediments6 <- sediments4[is.na(sediments4$unit_id),]
sediments5 <- sediments5[!is.na(sediments5$unit_id),]
sediments6 <- sediments6[!duplicated(sediments6$strat_name_id),]
sediments4 <- rbind(sediments5, sediments6)

sediments <- rbind(sediments1, sediments2, sediments3, sediments4)

sediments <- melt(sediments, id.vars = c(1:83,626), na.rm=TRUE)
sediments$value <- as.numeric(sediments$value)
sediments <- subset(sediments, variable != "b_age2")

# Phanerozoic bins

sediments_bins1 <- hist(sediments$value[sediments$value >= 0 & sediments$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
sediments_bins1 <- as.data.frame(cbind(sediments_bins1$counts, sediments_bins1$breaks))
sediments_bins1 <- sediments_bins1[1:542,]

# Precambrian bins
sediments_bins2 <- hist(sediments$value[sediments$value >= cutoff & sediments$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
sediments_bins2 <- as.data.frame(cbind(sediments_bins2$counts, sediments_bins2$breaks))

# Bind
sediments_bins <- rbind(sediments_bins1,sediments_bins2)
sediments_bins <- sediments_bins[c(1:541,543:889),] # remove duplicate 541

# block end

### import Wilkin-framboids ##

wilkin_comp <-  read.csv(file = "wilkin_framboids_comp.csv", row.names = 1)

sediments <- (grepl(paste(rocks, collapse = "|"), wilkin_comp$lith, ignore.case=TRUE) | 
                grepl(paste(rocks, collapse = "|"), wilkin_comp$other, ignore.case=TRUE) | 
                grepl(paste(rocks, collapse = "|"), wilkin_comp$environ, ignore.case=TRUE) | 
                grepl(paste(rocks, collapse = "|"), wilkin_comp$highlight, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), wilkin_comp$strat_name_long.x, ignore.case=TRUE)) &
  (grepl(paste(rocks2, collapse = "|"), wilkin_comp$lith, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), wilkin_comp$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), wilkin_comp$environ, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), wilkin_comp$highlight, ignore.case=TRUE) | 
  grepl(paste(rocks2, collapse = "|"), wilkin_comp$strat_name_long.x, ignore.case=TRUE))

wilkin_comp <- subset(wilkin_comp, sediments) # subset to 'rocks' of interest

wilkin_comp <- melt(wilkin_comp, id.vars = c(1:87), na.rm=TRUE)
wilkin_comp$value <- as.numeric(wilkin_comp$value)
wilkin_comp <- subset(wilkin_comp, variable != "b_age2")

# all strat names mentioned in Wilkin-pyrite papers

# Phanerozoic bins

wilkin_comp_bins1 <- hist(wilkin_comp$value[wilkin_comp$value >= 0 & wilkin_comp$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
wilkin_comp_bins1 <- as.data.frame(cbind(wilkin_comp_bins1$counts, wilkin_comp_bins1$breaks))
wilkin_comp_bins1 <- wilkin_comp_bins1[1:542,]

# Precambrian bins
wilkin_comp_bins2 <- hist(wilkin_comp$value[wilkin_comp$value >= cutoff & wilkin_comp$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
wilkin_comp_bins2 <- as.data.frame(cbind(wilkin_comp_bins2$counts, wilkin_comp_bins2$breaks))

# Bind
wilkin_comp_bins <- rbind(wilkin_comp_bins1,wilkin_comp_bins2)

wilkin_comp_bins <- wilkin_comp_bins[c(1:541,543:889),] # remove duplicate 541

# cull to 1 strat names mentioned in Wilkin-pyrite papers

wilkin_comp.one <- subset(wilkin_comp, count <= 1)

# Phanerozoic bins

wilkin_comp.one_bins1 <- hist(wilkin_comp.one$value[wilkin_comp.one$value >= 0 & wilkin_comp.one$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
wilkin_comp.one_bins1 <- as.data.frame(cbind(wilkin_comp.one_bins1$counts, wilkin_comp.one_bins1$breaks))
wilkin_comp.one_bins1 <- wilkin_comp.one_bins1[1:542,]

# Precambrian bins
wilkin_comp.one_bins2 <- hist(wilkin_comp.one$value[wilkin_comp.one$value >= cutoff & wilkin_comp.one$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
wilkin_comp.one_bins2 <- as.data.frame(cbind(wilkin_comp.one_bins2$counts, wilkin_comp.one_bins2$breaks))

# Bind
wilkin_comp.one_bins <- rbind(wilkin_comp.one_bins1,wilkin_comp.one_bins2)

wilkin_comp.one_bins <- wilkin_comp.one_bins[c(1:541,543:889),] # remove duplicate 541

# cull to 3 strat names mentioned in Wilkin-pyrite papers

wilkin_comp.three <- subset(wilkin_comp, count <= 3)

# Phanerozoic bins

wilkin_comp.three_bins1 <- hist(wilkin_comp.three$value[wilkin_comp.three$value >= 0 & wilkin_comp.three$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
wilkin_comp.three_bins1 <- as.data.frame(cbind(wilkin_comp.three_bins1$counts, wilkin_comp.three_bins1$breaks))
wilkin_comp.three_bins1 <- wilkin_comp.three_bins1[1:542,]

# Precambrian bins
wilkin_comp.three_bins2 <- hist(wilkin_comp.three$value[wilkin_comp.three$value >= cutoff & wilkin_comp.three$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
wilkin_comp.three_bins2 <- as.data.frame(cbind(wilkin_comp.three_bins2$counts, wilkin_comp.three_bins2$breaks))

# Bind
wilkin_comp.three_bins <- rbind(wilkin_comp.three_bins1,wilkin_comp.three_bins2)

wilkin_comp.three_bins <- wilkin_comp.three_bins[c(1:541,543:889),] # remove duplicate 541

# cull to 5 strat names mentioned in Wilkin-pyrite papers

wilkin_comp.five <- subset(wilkin_comp, count <= 5)

# Phanerozoic bins

wilkin_comp.five_bins1 <- hist(wilkin_comp.five$value[wilkin_comp.five$value >= 0 & wilkin_comp.five$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
wilkin_comp.five_bins1 <- as.data.frame(cbind(wilkin_comp.five_bins1$counts, wilkin_comp.five_bins1$breaks))
wilkin_comp.five_bins1 <- wilkin_comp.five_bins1[1:542,]

# Precambrian bins
wilkin_comp.five_bins2 <- hist(wilkin_comp.five$value[wilkin_comp.five$value >= cutoff & wilkin_comp.five$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
wilkin_comp.five_bins2 <- as.data.frame(cbind(wilkin_comp.five_bins2$counts, wilkin_comp.five_bins2$breaks))

# Bind
wilkin_comp.five_bins <- rbind(wilkin_comp.five_bins1,wilkin_comp.five_bins2)

wilkin_comp.five_bins <- wilkin_comp.five_bins[c(1:541,543:889),] # remove duplicate 541

# now for normalisation to sedimentary bins, as stacked output

# block start

# optional additional time points for plotting

HEATT <- c(56, 66, 93, 116, 183, 200, 251, 359, 372, 383, 444, 514, 541) 
ME <- c(445, 372, 252, 201, 65)
OAEs <- c(183, 120, 111, 93) # from Jenkyns 2010
PETM <- 55.8 # from Jenkyns 2010
Sturt <- c(716, 663) # Sturtian glaciation from ???
Es <- c(517, 502, 405, 393, 388, 382, 359, 330, 249, 240, 230, 220, 188, 145)
Extras <- c(2500, 1600, 1000, 850, 635, 541, 485.4, 443.8, 419.2, 358.9, 298.9, 
            251.9, 201.3, 145, 66, 23.03, 2.58) # Chronostrat divisions
supercontinents <- c(320, 170, 900, 700, 1600, 1400) # from Li et al. 2019 Precambrian Research
carb.intervals <- c(541, 465, 372, 323, 265, 227, 164, 133)
lows <- c(323, 299, 201, 170) # carb.intervals and lows from Riding et al. 2019
glacials <- cbind(c(2360, 2400,2450,717,650,467,330,315,307,305,300,34),
                  c(2380,2440,2460,660,635,444,317,311,308,306,288,0))
glacials <- as.data.frame(glacials)
BIF <- cbind(c(3000,1000, 541), c(1800, 542, 527))
BIF <- as.data.frame(BIF)
HEBS <- cbind(c(130,310, 335,330,362,388,490,505,2050), c(89, 305, 336, 323,359,383,491,506,2060))
HEBS <- as.data.frame(HEBS)

SEDEX <- read.csv("SEDEX.csv")

MVT <- subset(SEDEX, prevtype == "MVT")
SEDEX <- subset(SEDEX, prevtype == "SEDEX")

# plots

Top <- -0 # top age for plot (in Ma), set at -0.5 in order to centre bins
Bottom <- 3000 # bottom age for plot (in Ma)
phanerozoic_increment <- 1
precambrian_increment <- 10

# use below for plot 1 - S indicates stacked

framboidsS <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1)/sediments_bins$V1,nodules_bins$V2))
nodulesS <- as.data.frame(cbind((nodules_bins$V1)/sediments_bins$V1,nodules_bins$V2))
veinsS <- as.data.frame(cbind((veins_bins$V1+framboids_bins$V1+nodules_bins$V1)/sediments_bins$V1,veins_bins$V2))
veins2S <- as.data.frame(cbind((veins2_bins$V1+framboids_bins$V1+nodules_bins$V1)/sediments_bins$V1,veins2_bins$V2))
undifS <- as.data.frame(cbind((pyrite_undif_bins$V1+framboids_bins$V1+nodules_bins$V1)/sediments_bins$V1,nodules_bins$V2))

#### output - Figure 1A ###

# also add thick precambrian-cambrian boundary..

a <- ggplot() + theme_bw() +
  geom_linerange(data = glacials, aes(xmin= V1, xmax = V2, y = 0.38), colour = "deepskyblue", size = 2) +
   annotate(x=c(3000,2850,2350,2100, 1900),y=0.4,label= c("BIF", "Major glacials", "SEDEX", "MVT", "HEBS"),
           vjust=0,hjust = 0, geom="text", colour = c("orange", "deepskyblue", "black", "grey50", "red"), fontface = 2, size = 3.5) +
  geom_linerange(data = BIF, aes(xmin= V1, xmax = V2, y = 0.39), colour = "orange", size = 2) +
  geom_point(data = SEDEX, aes(agemy, y = 0.36)) +
  geom_point(data = MVT, aes(agemy, y = 0.37), colour = "grey50") +
  geom_point(data = HEBS, aes(V1, y = 0.37), colour = "red") +
  geom_point(data = HEBS, aes(V2, y = 0.37), colour = "red") +
  geom_linerange(data = HEBS, aes(xmax= V1, xmin = V2, y = 0.37), colour = "red", size = 2) +
  scale_x_reverse(limits = c(Bottom, 542), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0, 0.35), clip="off") +
  geom_vline(xintercept = c(Extras)) +
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  geom_stepribbon(data = undifS, aes(V2-0.5, ymin = 0, ymax = V1, fill = "Undifferentiated pyrite"), colour = "black") +
  geom_stepribbon(data = framboidsS, aes(V2-0.5, ymin = 0, ymax = V1, fill = "Framboids"), colour = "black") +
  geom_stepribbon(data = nodulesS, aes(V2-0.5, ymin = 0, ymax = V1, fill = "Nodules/concretions"), colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.95,0,0.11,0), "cm"),
        axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        plot.tag.position = c(0.01,0.99),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = 'white'),
        axis.text.y=element_text(size = 12),
        legend.position = c(0.44, 0.9),
        legend.background = element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.key.size = unit(0.7, 'cm')) +
  ylab("Proportion of\npyrite-bearing\n(meta)sedimentary\nrocks") +
  labs(tag = "A")  +
  annotate(x=c(HEATT)-5,y=0.35,label=c(letters[length(HEATT):1]),
           vjust=1,geom="label", colour = "black", fontface = 2, size = 3.5,
           label.padding = unit(0.18, "lines"),
           label.r = unit(0.11, "lines"),
           label.size = NA) +
  scale_fill_manual(values = c("Undifferentiated pyrite" = "grey90",
                               "Framboids" = "red",
                               "Nodules/concretions" = "blue"),
                     name = NULL, guide = 'legend', 
                     labels = c('Undifferentiated pyrite','Framboids', 'Nodules/concretions'))

b <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(542, Top), expand = c(0,0)) +
  geom_linerange(data = BIF, aes(xmin= V1, xmax = V2, y = 0.39), colour = "orange", size = 2) +
  geom_linerange(data = glacials, aes(xmin= V1, xmax = V2, y = 0.38), colour = "deepskyblue", size = 2) +
  geom_point(data = SEDEX, aes(agemy, y = 0.36)) +
  geom_point(data = MVT, aes(agemy, y = 0.37), colour = "grey50") +
  geom_point(data = HEBS, aes(V1, y = 0.37), colour = "red") +
  geom_point(data = HEBS, aes(V2, y = 0.37), colour = "red") +
  geom_linerange(data = HEBS, aes(xmax= V1, xmin = V2, y = 0.37), colour = "red", size = 2) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0, 0.35), clip="off") +
  geom_vline(xintercept = c(Extras)) +
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  geom_stepribbon(data = undifS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "grey90", colour = "black") +
  geom_stepribbon(data = framboidsS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "red", colour = "black") +
  geom_stepribbon(data = nodulesS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "blue", colour = "black") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.95,0.15,0.11,-0.16), "cm"),
        panel.background = element_rect(fill = 'white')) +
  annotate(x=c(HEATT)-5,y=0.35,label=c(letters[length(HEATT):1]),
           vjust=1,geom="label", colour = "black", fontface = 2, size = 3.5,
           label.padding = unit(0.18, "lines"),
           label.r = unit(0.11, "lines"),
           label.size = NA)

# generate plots 

grid.arrange(a,b, ncol = 2)


# alternative used for plot 2 - normalised to all pyrite-bearing 'rocks',  - R indicates stacked

framboidsS <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1)/(framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1),nodules_bins$V2))
nodulesS <- as.data.frame(cbind((nodules_bins$V1)/(framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1),nodules_bins$V2))
veinsS <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1+veins_bins$V1)/(framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1),nodules_bins$V2))
undifS <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1)/(framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1),nodules_bins$V2))

#### output - Figure 1B ###

c <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 542), expand = c(0,0)) +
  geom_vline(xintercept = c(Extras))  +
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.45)) +
  geom_stepribbon(data = framboidsS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "red", colour = "black") +
  geom_stepribbon(data = nodulesS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "blue", colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.11,0,0.11,0), "cm"),
        axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        plot.tag.position = c(0.01,0.99),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size = 12),
        panel.background = element_rect(fill = 'white')) +
  ylab("Proportion of\nframboid and\nnodule-bearing\npyritic \n(meta)sedimentary\nrocks") +
  labs(tag = "B")

d <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(542, Top), expand = c(0,0))  +
  geom_vline(xintercept = c(Extras)) +
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.45)) +
  geom_stepribbon(data = framboidsS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "red", colour = "black") +
  geom_stepribbon(data = nodulesS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "blue", colour = "black") +
  geom_step(data = framboidsS, aes(V2-0.5, V1)) + 
  geom_step(data = nodulesS, aes(V2-0.5, V1)) + 
  #geom_step(data = undifS, aes(V2-0.5, V1)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.11,0.15,0.11,-0.16), "cm"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = 'white'))

#Extras, ME, OAEs, 

# generate plots 

grid.arrange(c,d, ncol = 2)

#### output - Figure 1C ###

# wilkin-framboids

wilk.framboidsS.all <- as.data.frame(cbind((wilkin_comp_bins$V1)/sediments_bins$V1,nodules_bins$V2))
wilk.framboidsS.one <- as.data.frame(cbind((wilkin_comp.one_bins$V1)/sediments_bins$V1,nodules_bins$V2))
wilk.framboidsS.three <- as.data.frame(cbind((wilkin_comp.three_bins$V1)/sediments_bins$V1,nodules_bins$V2))
wilk.framboidsS.five <- as.data.frame(cbind((wilkin_comp.five_bins$V1)/sediments_bins$V1,nodules_bins$V2))

e <-  ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 542), expand = c(0,0.11)) +
  geom_vline(xintercept = c(Extras)) + 
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.25)) +
  #geom_stepribbon(data = wilk.framboidsS.all, aes(V2-0.5, ymin = 0, ymax = V1, fill = "Max"), colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.five, aes(V2-0.5, ymin = 0, ymax = V1, fill = "5"), colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.three, aes(V2-0.5, ymin = 0, ymax = V1, fill = "3"), colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.one, aes(V2-0.5, ymin = 0, ymax = V1, fill = "1"), colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0,0,0), "cm"),
        axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        plot.tag.position = c(0.01,0.99),
        panel.background = element_rect(fill = 'white'),
        axis.text.y=element_text(size = 12),
        axis.text.x=element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = c(0.48, 0.9),
        legend.background = element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.key.size = unit(0.7, 'cm')) +
  ylab("Proportion of\nWilkin-pyrite\n(meta)sedimentary\nrocks") +
  coord_geo(dat = list("periods", "eras"), xlim = c(Bottom, 542), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Cm", "Msr")) +
  labs(tag = "C") +
  scale_fill_manual(values = c("Max" = "lightsalmon",
                               "5" = "tomato1",
                               "3" = "orangered2",
                               "1" = "orangered4"),
                    name = 'Mentions of (meta)sedimentary rocks \nin "Wilkin" pyrite publications', guide = 'legend', 
                    labels = c('Max','5', '3', '1'))


f <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(542, Top), expand = c(0,0.11)) +
  geom_vline(xintercept = c(Extras))  +
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  #geom_stepribbon(data = wilk.framboidsS.all, aes(V2-0.5, ymin = 0, ymax = V1), fill = "lightsalmon", colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.five, aes(V2-0.5, ymin = 0, ymax = V1), fill = "tomato1", colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.three, aes(V2-0.5, ymin = 0, ymax = V1), fill = "orangered2", colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.one, aes(V2-0.5, ymin = 0, ymax = V1), fill = "orangered4", colour = "black") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0.15,0,-0.16), "cm"),
        panel.background = element_rect(fill = 'white'),
        axis.text.x=element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0,0.25), expand = c(0,0)) +
  coord_geo(dat = list("periods", "eras"), xlim = c(542, Top),
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE), size = 3)

grid.arrange(a,b,c,d,e,f,ncol = 2, bottom = textGrob("Age (Ma)"))



# optional data export

out <- cbind(framboids_bins, nodules_bins, veins_bins, pyrite_undif_bins, sediments_bins, wilkin_comp.one_bins,wilkin_comp.three_bins,wilkin_comp.five_bins,wilkin_comp_bins)
out <- out[,c(1,3,5,7,9,11,13,15,17:18)]
names(out)[1] <- "framboids"
names(out)[2] <- "nodules"
names(out)[3] <- "min"
names(out)[4] <- "undif"
names(out)[5] <- "seds"
names(out)[6] <- "Wilkin.one"
names(out)[7] <- "Wilkin.three"
names(out)[8] <- "Wilkin.five"
names(out)[9] <- "Wilkin.all"
names(out)[10] <- "Age"

write.csv(out, "xdd_binned_results.csv") 

# block end

###### PART 2 - this section is to generate an output without propagation of units ##### 
# Retain sed.list generated in part 1
# sed.list is derived the assumption that a strat package can be considered
# 'sedimentary' if the majority of units within that package are sedimentary, if
# no other information is available (i.e., from concepts)

# block start

# import file

data_p1 <-  read.csv(file = "data_part1_comp.csv", row.names = 1)

# optional compilation of target words and summary
phrases <- as.data.frame(as.array(data_p1$target_word))
list.phrases <- as.data.frame(summary(data_p1$target_word))

# compile framboid mentions, including in 'phrase' for pyrite undif hits

framboids <- grepl("framboid", data_p1$target_word, ignore.case=TRUE) |
  grepl("framboid", data_p1$phrase, ignore.case=TRUE)

framboids <- subset(data_p1, framboids)

# check mentions of pyrite framboid in concept description

framboids1 <- grepl("framboid", data_p1$other, ignore.case=TRUE)
framboids1 <- subset(data_p1, framboids1)
framboids <- rbind(framboids, framboids1)

# some duplicate results_id should not be removed - these are present for strat packages spanning the Precambrian-Phanerozoic, e.g., Tindir Gp
# therefore need to subdivide before removing duplicates associated with multiple hits

cutoff <- 541 # Ma 
phanerozoic_increment <- 1 # in Ma
precambrian_increment <- 10 # in Ma

framboids1 <- subset(framboids, t_age2 > cutoff & b_age2 > cutoff)
framboids1 <- framboids1[!duplicated(framboids1$strat_name_id),]

framboids2 <- subset(framboids, t_age2 < cutoff & b_age2 < cutoff)
framboids2 <- framboids2[!duplicated(framboids2$strat_name_id),]

framboids3 <- subset(framboids, t_age2 < cutoff & b_age2 == cutoff)
framboids3 <- framboids3[!duplicated(framboids3$strat_name_id),]

framboids4 <- subset(framboids, t_age2 == cutoff & b_age2 > cutoff)
framboids4 <- framboids4[!duplicated(framboids4$strat_name_id),]

framboids <- rbind(framboids1, framboids2, framboids3, framboids4)

# subset so only sediments extracted - use sed.list from part 1 - based on unit lith descriptions
# sed.list is defined using 'rocks'

framboids <- safe_right_join(sed.list, framboids, by = "strat_name_id", conflict = coalesce)

sediments <-  grepl("TRUE", framboids$lith2) |
  (grepl(paste(rocks, collapse = "|"), framboids$other, ignore.case=TRUE) |
  grepl(paste(rocks, collapse = "|"), framboids$strat_phrase_root, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), framboids$phrase, ignore.case=TRUE) | 
  grepl(paste(rocks, collapse = "|"), framboids$strat_flag, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), framboids$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), framboids$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), framboids$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), framboids$phrase, ignore.case=TRUE))

framboids <- subset(framboids, sediments) # subset to 'rocks' of interest

framboids <- select(framboids, -c(lith2))

# block end

# repeat for nodules

# block start

nodules <- grepl("nodul", data_p1$target_word, ignore.case=TRUE) | 
  grepl("concretion", data_p1$target_word, ignore.case=TRUE) 

nodules <- subset(data_p1, nodules)

# include nodule mentions pulled from lexicons (concepts)

nodules1 <- grepl("pyrite nodul", data_p1$other, ignore.case=TRUE) | grepl("pyrite concretion", data_p1$other, ignore.case=TRUE) |
  grepl("pyritic nodul", data_p1$other, ignore.case=TRUE) | grepl("pyritic concretion", data_p1$other, ignore.case=TRUE) |
  grepl("nodules of pyrite", data_p1$other, ignore.case=TRUE) | grepl("concretions of pyrite", data_p1$other, ignore.case=TRUE) | 
  grepl("nodular pyrite", data_p1$other, ignore.case=TRUE) | grepl("concretionary pyrite", data_p1$other, ignore.case=TRUE)

nodules1 <- subset(data_p1, nodules1)
nodules <- rbind(nodules, nodules1)

nodules1 <- subset(nodules, t_age2 > cutoff & b_age2 > cutoff)
nodules1 <- nodules1[!duplicated(nodules1$strat_name_id),]

nodules2 <- subset(nodules, t_age2 < cutoff & b_age2 < cutoff)
nodules2 <- nodules2[!duplicated(nodules2$strat_name_id),]

nodules3 <- subset(nodules, t_age2 < cutoff & b_age2 == cutoff)
nodules3 <- nodules3[!duplicated(nodules3$strat_name_id),]

nodules4 <- subset(nodules, t_age2 == cutoff & b_age2 > cutoff)
nodules4 <- nodules4[!duplicated(nodules4$strat_name_id),]

nodules <- rbind(nodules1, nodules2, nodules3, nodules4)

# subset so only sediments extracted - use sed.list from part 1 - based on unit lith descriptions
# sed.list is defined using 'rocks'

nodules <- safe_right_join(sed.list, nodules, by = "strat_name_id", conflict = coalesce)

sediments <- grepl("TRUE", nodules$lith2) |
  (grepl(paste(rocks, collapse = "|"), nodules$other, ignore.case=TRUE) |
     grepl(paste(rocks, collapse = "|"), nodules$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks, collapse = "|"), nodules$phrase, ignore.case=TRUE) | 
     grepl(paste(rocks, collapse = "|"), nodules$strat_flag, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), nodules$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), nodules$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), nodules$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), nodules$phrase, ignore.case=TRUE))

nodules <- subset(nodules, sediments) # subset to 'rocks' of interest

nodules <- select(nodules, -c(lith2))

# block end

# now to process undifferentiated pyrite mentions

# block start

# pyritic strat concepts pulled from lexicons (undifferentiated)
pyritic_strat <- grepl("pyrit", data_p1$other, ignore.case=TRUE)
pyritic_strat <- subset(data_p1, pyritic_strat)
false <- !grepl("non-pyrit", pyritic_strat$other, ignore.case=TRUE)
pyritic_strat <- subset(pyritic_strat, false)

all0 <- rbind(framboids, nodules)

all0 <- all0[!duplicated(all0$strat_name_id),]

pyritic_strat <- safe_right_join(all0, pyritic_strat, by = "strat_name_id", conflict = coalesce)
pyritic_strat <- pyritic_strat[is.na(pyritic_strat$target_word),]

pyrite <- grepl("\\<pyrite\\>", data_p1$target_word, ignore.case=TRUE) | grepl("\\<pyritic\\>", data_p1$target_word, ignore.case=TRUE)

pyrite_undif <- subset(data_p1, pyrite)

pyrite_undif <- rbind(pyritic_strat, pyrite_undif)

# subset so only sediments extracted - use sed.list from part 1 - based on unit lith descriptions
# sed.list is defined using 'rocks'

pyrite_undif <- safe_right_join(sed.list, pyrite_undif, by = "strat_name_id", conflict = coalesce)

# in theory lith2 is the only requirement here (if defined by ALL sediments), but other terms added as a precaution

sediments <- grepl("TRUE", pyrite_undif$lith2) |
  (grepl(paste(rocks, collapse = "|"), pyrite_undif$other, ignore.case=TRUE) |
     grepl(paste(rocks, collapse = "|"), pyrite_undif$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks, collapse = "|"), pyrite_undif$phrase, ignore.case=TRUE) | 
     grepl(paste(rocks, collapse = "|"), pyrite_undif$strat_flag, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), pyrite_undif$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), pyrite_undif$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), pyrite_undif$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), pyrite_undif$phrase, ignore.case=TRUE))

pyrite_undif <- subset(pyrite_undif, sediments)

##### carry forward only strat packages not already containing framboid or nodule phrases

pyrite_undif <- anti_join(pyrite_undif, all0, by = "strat_name_id") 

# drop lith2 as not needed from this point onwards

pyrite_undif <- select(pyrite_undif, -c(lith2))

pyrite_undif1 <- subset(pyrite_undif, t_age2 > cutoff & b_age2 > cutoff)
pyrite_undif1 <- pyrite_undif1[!duplicated(pyrite_undif1$strat_name_id),]

pyrite_undif2 <- subset(pyrite_undif, t_age2 < cutoff & b_age2 < cutoff)
pyrite_undif2 <- pyrite_undif2[!duplicated(pyrite_undif2$strat_name_id),]

pyrite_undif3 <- subset(pyrite_undif, t_age2 < cutoff & b_age2 == cutoff)
pyrite_undif3 <- pyrite_undif3[!duplicated(pyrite_undif3$strat_name_id),]

pyrite_undif4 <- subset(pyrite_undif, t_age2 == cutoff & b_age2 > cutoff)
pyrite_undif4 <- pyrite_undif4[!duplicated(pyrite_undif4$strat_name_id),]

pyrite_undif <- rbind(pyrite_undif1, pyrite_undif2, pyrite_undif3, pyrite_undif4)

# alternative for pyrite undif taking only pyrite extracted via xdd

pyrite_xdd <- data_p1

pyrite_xdd <- pyrite_xdd[!is.na(pyrite_xdd$result_id),]

pyrite_xdd <- safe_right_join(sed.list, pyrite_xdd, by = "strat_name_id", conflict = coalesce)

sediments <- grepl("TRUE", pyrite_xdd$lith2) |
  (grepl(paste(rocks, collapse = "|"), pyrite_xdd$other, ignore.case=TRUE) |
     grepl(paste(rocks, collapse = "|"), pyrite_xdd$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks, collapse = "|"), pyrite_xdd$phrase, ignore.case=TRUE) | 
     grepl(paste(rocks, collapse = "|"), pyrite_xdd$strat_flag, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), pyrite_xdd$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), pyrite_xdd$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), pyrite_xdd$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), pyrite_xdd$phrase, ignore.case=TRUE))

pyrite_xdd <- subset(pyrite_xdd, sediments)

##### carry forward only strat packages not already containing framboid or nodule phrases

pyrite_xdd <- anti_join(pyrite_xdd, all0, by = "strat_name_id") 

# drop lith2 as not needed from this point onwards

pyrite_xdd <- select(pyrite_xdd, -c(lith2))

pyrite_xdd1 <- subset(pyrite_xdd, t_age2 > cutoff & b_age2 > cutoff)
pyrite_xdd1 <- pyrite_xdd1[!duplicated(pyrite_xdd1$strat_name_id),]

pyrite_xdd2 <- subset(pyrite_xdd, t_age2 < cutoff & b_age2 < cutoff)
pyrite_xdd2 <- pyrite_xdd2[!duplicated(pyrite_xdd2$strat_name_id),]

pyrite_xdd3 <- subset(pyrite_xdd, t_age2 < cutoff & b_age2 == cutoff)
pyrite_xdd3 <- pyrite_xdd3[!duplicated(pyrite_xdd3$strat_name_id),]

pyrite_xdd4 <- subset(pyrite_xdd, t_age2 == cutoff & b_age2 > cutoff)
pyrite_xdd4 <- pyrite_xdd4[!duplicated(pyrite_xdd4$strat_name_id),]

pyrite_xdd <- rbind(pyrite_xdd1, pyrite_xdd2, pyrite_xdd3, pyrite_xdd4)

# veins/mineralisation extraction as in part 1

veins <- grepl("vein", pyrite_undif$phrase, ignore.case=TRUE)  |  
  grepl("mineralisation", pyrite_undif$phrase, ignore.case=TRUE)

veins <- subset(pyrite_undif, veins)

# block end

# block start

# melt datasets - the same workflow as part 1, but with different column ranges

framboids <- melt(framboids, id.vars = c(1:44), na.rm=TRUE)
framboids$value <- as.numeric(framboids$value)
framboids <- subset(framboids, variable != "b_age2")

nodules <- melt(nodules, id.vars = c(1:44), na.rm=TRUE)
nodules$value <- as.numeric(nodules$value)
nodules <- subset(nodules, variable != "b_age2")

pyrite_undif <- melt(pyrite_undif, id.vars = c(1:44), na.rm=TRUE)
pyrite_undif$value <- as.numeric(pyrite_undif$value)
pyrite_undif <- subset(pyrite_undif, variable != "b_age2")

pyrite_xdd <- melt(pyrite_xdd, id.vars = c(1:44), na.rm=TRUE)
pyrite_xdd$value <- as.numeric(pyrite_xdd$value)
pyrite_xdd <- subset(pyrite_xdd, variable != "b_age2")

veins <- melt(veins, id.vars = c(1:44), na.rm=TRUE)
veins$value <- as.numeric(veins$value)
veins <- subset(veins, variable != "b_age2")

# block end

# bin the datasets

# block start

phanerozoic_increment <- 1 # in Ma
precambrian_increment <- 10 # in MA

# Phanerozoic bins

framboids_bins1 <- hist(framboids$value[framboids$value >= 0 & framboids$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
framboids_bins1 <- as.data.frame(cbind(framboids_bins1$counts, framboids_bins1$breaks))
framboids_bins1 <- framboids_bins1[1:542,]

# Precambrian bins
framboids_bins2 <- hist(framboids$value[framboids$value >= cutoff & framboids$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
framboids_bins2 <- as.data.frame(cbind(framboids_bins2$counts, framboids_bins2$breaks))

# Bind
framboids_bins <- rbind(framboids_bins1,framboids_bins2)

framboids_bins <- framboids_bins[c(1:541,543:889),] # remove duplicate 541

# Phanerozoic bins
nodules_bins1 <- hist(nodules$value[nodules$value >= 0 & nodules$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
nodules_bins1 <- as.data.frame(cbind(nodules_bins1$counts, nodules_bins1$breaks))
nodules_bins1 <- nodules_bins1[1:542,]

# Precambrian bins
nodules_bins2 <- hist(nodules$value[nodules$value >= cutoff & nodules$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
nodules_bins2 <- as.data.frame(cbind(nodules_bins2$counts, nodules_bins2$breaks))

# Bind
nodules_bins <- rbind(nodules_bins1,nodules_bins2)

nodules_bins <- nodules_bins[c(1:541,543:889),] # remove duplicate 541

# undifferentiated pyrite mentions

# Phanerozoic bins
pyrite_undif_bins1 <- hist(pyrite_undif$value[pyrite_undif$value >= 0 & pyrite_undif$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
pyrite_undif_bins1 <- as.data.frame(cbind(pyrite_undif_bins1$counts, pyrite_undif_bins1$breaks))
pyrite_undif_bins1 <- pyrite_undif_bins1[1:542,]

# Precambrian bins
pyrite_undif_bins2 <- hist(pyrite_undif$value[pyrite_undif$value >= cutoff & pyrite_undif$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
pyrite_undif_bins2 <- as.data.frame(cbind(pyrite_undif_bins2$counts, pyrite_undif_bins2$breaks))

# Bind
pyrite_undif_bins <- rbind(pyrite_undif_bins1,pyrite_undif_bins2)

pyrite_undif_bins <- pyrite_undif_bins[c(1:541,543:889),] # remove duplicate 541

# undifferentiated pyrite mentions - xdd only

# Phanerozoic bins
pyrite_xdd_bins1 <- hist(pyrite_xdd$value[pyrite_xdd$value >= 0 & pyrite_xdd$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
pyrite_xdd_bins1 <- as.data.frame(cbind(pyrite_xdd_bins1$counts, pyrite_xdd_bins1$breaks))
pyrite_xdd_bins1 <- pyrite_xdd_bins1[1:542,]

# Precambrian bins
pyrite_xdd_bins2 <- hist(pyrite_xdd$value[pyrite_xdd$value >= cutoff & pyrite_xdd$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
pyrite_xdd_bins2 <- as.data.frame(cbind(pyrite_xdd_bins2$counts, pyrite_xdd_bins2$breaks))

# Bind
pyrite_xdd_bins <- rbind(pyrite_xdd_bins1,pyrite_xdd_bins2)

pyrite_xdd_bins <- pyrite_xdd_bins[c(1:541,543:889),] # remove duplicate 541

# undifferentiated pyrite mentions - associated with veins and mineralisation mentions

# Phanerozoic bins
veins_bins1 <- hist(veins$value[veins$value >= 0 & veins$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
veins_bins1 <- as.data.frame(cbind(veins_bins1$counts, veins_bins1$breaks))
veins_bins1 <- veins_bins1[1:542,]

# Precambrian bins
veins_bins2 <- hist(veins$value[veins$value >= cutoff & veins$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
veins_bins2 <- as.data.frame(cbind(veins_bins2$counts, veins_bins2$breaks))

# Bind
veins_bins <- rbind(veins_bins1,veins_bins2)

veins_bins <- veins_bins[c(1:541,543:889),] # remove duplicate 541

# block end

# now to normalise to sedimentary packages - including nodules and framboids if not identified as 'sedimentary'

# block start

# use sed info from lith for N. American units, upgraded to strat package level, using Mode as defined in part 1

data.sed <- safe_right_join(sed.list, data_p1, by = "strat_name_id", conflict = coalesce)

sediments <- grepl("TRUE", data.sed$lith2) |
  (grepl(paste(rocks, collapse = "|"), data.sed$other, ignore.case=TRUE) |
     grepl(paste(rocks, collapse = "|"), data.sed$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks, collapse = "|"), data.sed$phrase, ignore.case=TRUE) | 
     grepl(paste(rocks, collapse = "|"), data.sed$strat_flag, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), data.sed$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), data.sed$strat_phrase_root, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), data.sed$strat_flag, ignore.case=TRUE) |
     grepl(paste(rocks2, collapse = "|"), data.sed$phrase, ignore.case=TRUE))

sediments <- subset(data.sed, sediments)

sediments <- select(sediments, -c(lith2))

sediments1 <- subset(sediments, t_age2 > cutoff & b_age2 > cutoff)
sediments1 <- sediments1[!duplicated(sediments1$strat_name_id),]

sediments2 <- subset(sediments, t_age2 < cutoff & b_age2 < cutoff)
sediments2 <- sediments2[!duplicated(sediments2$strat_name_id),]

sediments3 <- subset(sediments, t_age2 < cutoff & b_age2 == cutoff)
sediments3 <- sediments3[!duplicated(sediments3$strat_name_id),]

sediments4 <- subset(sediments, t_age2 == cutoff & b_age2 > cutoff)
sediments4 <- sediments4[!duplicated(sediments4$strat_name_id),]

sediments <- rbind(sediments1, sediments2, sediments3, sediments4)

sediments <- melt(sediments, id.vars = c(1:44), na.rm=TRUE)
sediments$value <- as.numeric(sediments$value)

sediments <- subset(sediments, variable != "b_age2")

sediments_bins1 <- hist(sediments$value[sediments$value >= 0 & sediments$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
sediments_bins1 <- as.data.frame(cbind(sediments_bins1$counts, sediments_bins1$breaks))
sediments_bins1 <- sediments_bins1[1:542,]

# Precambrian bins
sediments_bins2 <- hist(sediments$value[sediments$value >= cutoff & sediments$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
sediments_bins2 <- as.data.frame(cbind(sediments_bins2$counts, sediments_bins2$breaks))

# Bind
sediments_bins <- rbind(sediments_bins1,sediments_bins2)

sediments_bins <- sediments_bins[c(1:541,543:889),] # remove duplicate 541

# block end

### import Wilkin-framboids ##

wilkin_comp <-  read.csv(file = "wilkin_framboids_strat.csv", row.names = 1)

sediments <- (grepl(paste(rocks, collapse = "|"), wilkin_comp$other, ignore.case=TRUE) | 
                grepl(paste(rocks, collapse = "|"), wilkin_comp$strat_name, ignore.case=TRUE) | 
                grepl(paste(rocks, collapse = "|"), wilkin_comp$highlight, ignore.case=TRUE))  &
  (grepl(paste(rocks2, collapse = "|"), wilkin_comp$other, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), wilkin_comp$strat_name, ignore.case=TRUE) | 
     grepl(paste(rocks2, collapse = "|"), wilkin_comp$highlight, ignore.case=TRUE))

wilkin_comp <- subset(wilkin_comp, sediments) # subset to 'rocks' of interest

wilkin_comp <- melt(wilkin_comp, id.vars = c(1:48), na.rm=TRUE)
wilkin_comp$value <- as.numeric(wilkin_comp$value)
wilkin_comp <- subset(wilkin_comp, variable != "b_age2")

# all strat names mentioned in Wilkin-pyrite papers

# Phanerozoic bins

wilkin_comp_bins1 <- hist(wilkin_comp$value[wilkin_comp$value >= 0 & wilkin_comp$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
wilkin_comp_bins1 <- as.data.frame(cbind(wilkin_comp_bins1$counts, wilkin_comp_bins1$breaks))
wilkin_comp_bins1 <- wilkin_comp_bins1[1:542,]

# Precambrian bins
wilkin_comp_bins2 <- hist(wilkin_comp$value[wilkin_comp$value >= cutoff & wilkin_comp$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
wilkin_comp_bins2 <- as.data.frame(cbind(wilkin_comp_bins2$counts, wilkin_comp_bins2$breaks))

# Bind
wilkin_comp_bins <- rbind(wilkin_comp_bins1,wilkin_comp_bins2)

wilkin_comp_bins <- wilkin_comp_bins[c(1:541,543:889),] # remove duplicate 541

# cull to 1 strat names mentioned in Wilkin-pyrite papers

wilkin_comp.one <- subset(wilkin_comp, count <= 1)

# Phanerozoic bins

wilkin_comp.one_bins1 <- hist(wilkin_comp.one$value[wilkin_comp.one$value >= 0 & wilkin_comp.one$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
wilkin_comp.one_bins1 <- as.data.frame(cbind(wilkin_comp.one_bins1$counts, wilkin_comp.one_bins1$breaks))
wilkin_comp.one_bins1 <- wilkin_comp.one_bins1[1:542,]

# Precambrian bins
wilkin_comp.one_bins2 <- hist(wilkin_comp.one$value[wilkin_comp.one$value >= cutoff & wilkin_comp.one$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
wilkin_comp.one_bins2 <- as.data.frame(cbind(wilkin_comp.one_bins2$counts, wilkin_comp.one_bins2$breaks))

# Bind
wilkin_comp.one_bins <- rbind(wilkin_comp.one_bins1,wilkin_comp.one_bins2)

wilkin_comp.one_bins <- wilkin_comp.one_bins[c(1:541,543:889),] # remove duplicate 541

# cull to 3 strat names mentioned in Wilkin-pyrite papers

wilkin_comp.three <- subset(wilkin_comp, count <= 3)

# Phanerozoic bins

wilkin_comp.three_bins1 <- hist(wilkin_comp.three$value[wilkin_comp.three$value >= 0 & wilkin_comp.three$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
wilkin_comp.three_bins1 <- as.data.frame(cbind(wilkin_comp.three_bins1$counts, wilkin_comp.three_bins1$breaks))
wilkin_comp.three_bins1 <- wilkin_comp.three_bins1[1:542,]

# Precambrian bins
wilkin_comp.three_bins2 <- hist(wilkin_comp.three$value[wilkin_comp.three$value >= cutoff & wilkin_comp.three$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
wilkin_comp.three_bins2 <- as.data.frame(cbind(wilkin_comp.three_bins2$counts, wilkin_comp.three_bins2$breaks))

# Bind
wilkin_comp.three_bins <- rbind(wilkin_comp.three_bins1,wilkin_comp.three_bins2)

wilkin_comp.three_bins <- wilkin_comp.three_bins[c(1:541,543:889),] # remove duplicate 541

# cull to 5 strat names mentioned in Wilkin-pyrite papers

wilkin_comp.five <- subset(wilkin_comp, count <= 5)

# Phanerozoic bins

wilkin_comp.five_bins1 <- hist(wilkin_comp.five$value[wilkin_comp.five$value >= 0 & wilkin_comp.five$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
wilkin_comp.five_bins1 <- as.data.frame(cbind(wilkin_comp.five_bins1$counts, wilkin_comp.five_bins1$breaks))
wilkin_comp.five_bins1 <- wilkin_comp.five_bins1[1:542,]

# Precambrian bins
wilkin_comp.five_bins2 <- hist(wilkin_comp.five$value[wilkin_comp.five$value >= cutoff & wilkin_comp.five$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
wilkin_comp.five_bins2 <- as.data.frame(cbind(wilkin_comp.five_bins2$counts, wilkin_comp.five_bins2$breaks))

# Bind
wilkin_comp.five_bins <- rbind(wilkin_comp.five_bins1,wilkin_comp.five_bins2)

wilkin_comp.five_bins <- wilkin_comp.five_bins[c(1:541,543:889),] # remove duplicate 541

# now for normalisation to sedimentary bins, as stacked output

# normalisation as in part 1

# block start


# optional additional time points for plotting

HEATT <- c(56, 66, 93, 116, 183, 200, 251, 359, 372, 383, 444, 514, 541) 
ME <- c(445, 372, 252, 201, 65)
OAEs <- c(183, 120, 111, 93) # from Jenkyns 2010
PETM <- 55.8 # from Jenkyns 2010
Sturt <- c(716, 663) # Sturtian glaciation from ???
Es <- c(517, 502, 405, 393, 388, 382, 359, 330, 249, 240, 230, 220, 188, 145)
Extras <- c(2500, 1600, 1000, 850, 635, 541, 485.4, 443.8, 419.2, 358.9, 298.9, 
            251.9, 201.3, 145, 66, 23.03, 2.58) # Chronostrat divisions
supercontinents <- c(320, 170, 900, 700, 1600, 1400) # from Li et al. 2019 Precambrian Research
carb.intervals <- c(541, 465, 372, 323, 265, 227, 164, 133)
lows <- c(323, 299, 201, 170) # carb.intervals and lows from Riding et al. 2019
glacials <- cbind(c(2360, 2400,2450,717,650,467,330,315,307,305,300,34),
                  c(2380,2440,2460,660,635,444,317,311,308,306,288,0))
glacials <- as.data.frame(glacials)
BIF <- cbind(c(3000,1000, 541), c(1800, 542, 527))
BIF <- as.data.frame(BIF)
HEBS <- cbind(c(130,310, 335,330,362,388,490,505,2050), c(89, 305, 336, 323,359,383,491,506,2060))
HEBS <- as.data.frame(HEBS)

SEDEX <- read.csv("SEDEX.csv")

MVT <- subset(SEDEX, prevtype == "MVT")
SEDEX <- subset(SEDEX, prevtype == "SEDEX")

# plots

Top <- -0 # top age for plot (in Ma), set at -0.5 in order to centre bins
Bottom <- 3000 # bottom age for plot (in Ma)
phanerozoic_increment <- 1
precambrian_increment <- 10

# use below for plot 1 - S indicates stacked

framboidsS <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1)/sediments_bins$V1,nodules_bins$V2))
nodulesS <- as.data.frame(cbind((nodules_bins$V1)/sediments_bins$V1,nodules_bins$V2))
veinsS <- as.data.frame(cbind((veins_bins$V1+framboids_bins$V1+nodules_bins$V1)/sediments_bins$V1,veins_bins$V2))
undifS <- as.data.frame(cbind((pyrite_undif_bins$V1+framboids_bins$V1+nodules_bins$V1)/sediments_bins$V1,nodules_bins$V2))

#### output - Figure 1A ###

a <- ggplot() + theme_bw() +
  geom_linerange(data = glacials, aes(xmin= V1, xmax = V2, y = 0.38), colour = "deepskyblue", size = 2) +
  annotate(x=c(3000,2850,2350,2100, 1900),y=0.4,label= c("BIF", "Major glacials", "SEDEX", "MVT", "HEBS"),
           vjust=0,hjust = 0, geom="text", colour = c("orange", "deepskyblue", "black", "grey50", "red"), fontface = 2, size = 3.5) +
  geom_linerange(data = BIF, aes(xmin= V1, xmax = V2, y = 0.39), colour = "orange", size = 2) +
  geom_point(data = SEDEX, aes(agemy, y = 0.36)) +
  geom_point(data = MVT, aes(agemy, y = 0.37), colour = "grey50") +
  geom_point(data = HEBS, aes(V1, y = 0.37), colour = "red") +
  geom_point(data = HEBS, aes(V2, y = 0.37), colour = "red") +
  geom_linerange(data = HEBS, aes(xmax= V1, xmin = V2, y = 0.37), colour = "red", size = 2) +
  scale_x_reverse(limits = c(Bottom, 542), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0, 0.35), clip="off") +
  geom_vline(xintercept = c(Extras)) +
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  geom_stepribbon(data = undifS, aes(V2-0.5, ymin = 0, ymax = V1, fill = "Undifferentiated pyrite"), colour = "black") +
  geom_stepribbon(data = framboidsS, aes(V2-0.5, ymin = 0, ymax = V1, fill = "Framboids"), colour = "black") +
  geom_stepribbon(data = nodulesS, aes(V2-0.5, ymin = 0, ymax = V1, fill = "Nodules/concretions"), colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.95,0,0.11,0), "cm"),
        axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        plot.tag.position = c(0.01,0.99),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = 'white'),
        axis.text.y=element_text(size = 12),
        legend.position = c(0.44, 0.9),
        legend.background = element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.key.size = unit(0.7, 'cm')) +
  ylab("Proportion of\npyrite-bearing\n(meta)sedimentary\nrocks") +
  labs(tag = "A")  +
  annotate(x=c(HEATT)-5,y=0.35,label=c(letters[length(HEATT):1]),
           vjust=1,geom="label", colour = "black", fontface = 2, size = 3.5,
           label.padding = unit(0.18, "lines"),
           label.r = unit(0.11, "lines"),
           label.size = NA) +
  scale_fill_manual(values = c("Undifferentiated pyrite" = "grey90",
                               "Framboids" = "red",
                               "Nodules/concretions" = "blue"),
                    name = NULL, guide = 'legend', 
                    labels = c('Undifferentiated pyrite','Framboids', 'Nodules/concretions'))

b <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(542, Top), expand = c(0,0)) +
  geom_linerange(data = BIF, aes(xmin= V1, xmax = V2, y = 0.39), colour = "orange", size = 2) +
  geom_linerange(data = glacials, aes(xmin= V1, xmax = V2, y = 0.38), colour = "deepskyblue", size = 2) +
  geom_point(data = SEDEX, aes(agemy, y = 0.36)) +
  geom_point(data = MVT, aes(agemy, y = 0.37), colour = "grey50") +
  geom_point(data = HEBS, aes(V1, y = 0.37), colour = "red") +
  geom_point(data = HEBS, aes(V2, y = 0.37), colour = "red") +
  geom_linerange(data = HEBS, aes(xmax= V1, xmin = V2, y = 0.37), colour = "red", size = 2) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0, 0.35), clip="off") +
  geom_vline(xintercept = c(Extras)) +
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  geom_stepribbon(data = undifS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "grey90", colour = "black") +
  geom_stepribbon(data = framboidsS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "red", colour = "black") +
  geom_stepribbon(data = nodulesS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "blue", colour = "black") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.95,0.15,0.11,-0.16), "cm"),
        panel.background = element_rect(fill = 'white')) +
  annotate(x=c(HEATT)-5,y=0.35,label=c(letters[length(HEATT):1]),
           vjust=1,geom="label", colour = "black", fontface = 2, size = 3.5,
           label.padding = unit(0.18, "lines"),
           label.r = unit(0.11, "lines"),
           label.size = NA)

# generate plots 

grid.arrange(a,b, ncol = 2)


# alternative used for plot 2 - normalised to all pyrite-bearing 'rocks',  - R indicates stacked

framboidsS <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1)/(framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1),nodules_bins$V2))
nodulesS <- as.data.frame(cbind((nodules_bins$V1)/(framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1),nodules_bins$V2))
veinsS <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1+veins_bins$V1)/(framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1),nodules_bins$V2))
undifS <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1)/(framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1),nodules_bins$V2))

#### output - Figure 1B ###

c <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 542), expand = c(0,0)) +
  geom_vline(xintercept = c(Extras))  +
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.45)) +
  geom_stepribbon(data = framboidsS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "red", colour = "black") +
  geom_stepribbon(data = nodulesS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "blue", colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.11,0,0.11,0), "cm"),
        axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        plot.tag.position = c(0.01,0.99),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size = 12),
        panel.background = element_rect(fill = 'white')) +
  ylab("Proportion of\nframboid and\nnodule-bearing\npyritic \n(meta)sedimentary\nrocks") +
  labs(tag = "B")

d <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(542, Top), expand = c(0,0))  +
  geom_vline(xintercept = c(Extras)) +
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.45)) +
  geom_stepribbon(data = framboidsS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "red", colour = "black") +
  geom_stepribbon(data = nodulesS, aes(V2-0.5, ymin = 0, ymax = V1), fill = "blue", colour = "black") +
  geom_step(data = framboidsS, aes(V2-0.5, V1)) + 
  geom_step(data = nodulesS, aes(V2-0.5, V1)) + 
  #geom_step(data = undifS, aes(V2-0.5, V1)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.11,0.15,0.11,-0.16), "cm"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = 'white'))

#Extras, ME, OAEs, 

# generate plots 

grid.arrange(c,d, ncol = 2)

#### output - Figure 1C ###

# wilkin-framboids

wilk.framboidsS.all <- as.data.frame(cbind((wilkin_comp_bins$V1)/sediments_bins$V1,nodules_bins$V2))
wilk.framboidsS.one <- as.data.frame(cbind((wilkin_comp.one_bins$V1)/sediments_bins$V1,nodules_bins$V2))
wilk.framboidsS.three <- as.data.frame(cbind((wilkin_comp.three_bins$V1)/sediments_bins$V1,nodules_bins$V2))
wilk.framboidsS.five <- as.data.frame(cbind((wilkin_comp.five_bins$V1)/sediments_bins$V1,nodules_bins$V2))

e <-  ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 542), expand = c(0,0.11)) +
  geom_vline(xintercept = c(Extras)) + 
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.1)) +
  #geom_stepribbon(data = wilk.framboidsS.all, aes(V2-0.5, ymin = 0, ymax = V1, fill = "Max"), colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.five, aes(V2-0.5, ymin = 0, ymax = V1, fill = "5"), colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.three, aes(V2-0.5, ymin = 0, ymax = V1, fill = "3"), colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.one, aes(V2-0.5, ymin = 0, ymax = V1, fill = "1"), colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0,0,0), "cm"),
        axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        plot.tag.position = c(0.01,0.99),
        panel.background = element_rect(fill = 'white'),
        axis.text.y=element_text(size = 12),
        axis.text.x=element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = c(0.48, 0.9),
        legend.background = element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.key.size = unit(0.7, 'cm')) +
  ylab("Proportion of\nWilkin-pyrite\n(meta)sedimentary\nrocks") +
  coord_geo(dat = list("periods", "eras"), xlim = c(Bottom, 542), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Cm", "Msr")) +
  labs(tag = "C") +
  scale_fill_manual(values = c("Max" = "lightsalmon",
                               "5" = "tomato1",
                               "3" = "orangered2",
                               "1" = "orangered4"),
                    name = 'Mentions of (meta)sedimentary rocks \nin "Wilkin" pyrite publications', guide = 'legend', 
                    labels = c('Max','5', '3', '1'))


f <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(542, Top), expand = c(0,0.11)) +
  geom_vline(xintercept = c(Extras))  +
  geom_vline(xintercept = c(HEATT), colour = "white", linetype = "dashed", size = 0.5) +
  #geom_stepribbon(data = wilk.framboidsS.all, aes(V2-0.5, ymin = 0, ymax = V1), fill = "lightsalmon", colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.five, aes(V2-0.5, ymin = 0, ymax = V1), fill = "tomato1", colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.three, aes(V2-0.5, ymin = 0, ymax = V1), fill = "orangered2", colour = "black") +
  geom_stepribbon(data = wilk.framboidsS.one, aes(V2-0.5, ymin = 0, ymax = V1), fill = "orangered4", colour = "black") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0.15,0,-0.16), "cm"),
        panel.background = element_rect(fill = 'white'),
        axis.text.x=element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0,0.1), expand = c(0,0)) +
  coord_geo(dat = list("periods", "eras"), xlim = c(542, Top),
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE), size = 3)

grid.arrange(a,b,c,d,e,f,ncol = 2, bottom = textGrob("Age (Ma)"))



# optional data export

out <- cbind(framboids_bins, nodules_bins, veins_bins, pyrite_undif_bins, sediments_bins, wilkin_comp.one_bins,wilkin_comp.three_bins,wilkin_comp.five_bins,wilkin_comp_bins)
out <- out[,c(1,3,5,7,9,11,13,15,17:18)]
names(out)[1] <- "framboids"
names(out)[2] <- "nodules"
names(out)[3] <- "min"
names(out)[4] <- "undif"
names(out)[5] <- "seds"
names(out)[6] <- "Wilkin.one"
names(out)[7] <- "Wilkin.three"
names(out)[8] <- "Wilkin.five"
names(out)[9] <- "Wilkin.all"
names(out)[10] <- "Age"

write.csv(out, "xdd_binned_results_strat_only.csv") 

## END ##

#### Supplementary Materials derived from PART 1 ####

## correlation chart - Fig. S8-S10 ##

# re-run PART 1 up to melt of the framboid and nodule datasets

framboids$type <- "Framboids"
nodules$type <- "Nodules"

length(unique(framboids$result_id))+length(unique(nodules$result_id))

framboids.2 <- framboids[,c("result_id", "docid", 
                            "strat_name_id", "strat_name_long",
                            "unit_id", "t_age", "b_age", 
                            "target_word", "phrase", "type", 
                            "other", "lith", "environ","clat","clng")] 

framboids.2 <- framboids.2[!duplicated(framboids.2[,c("strat_name_id", "unit_id")]),]

nodules.2 <- nodules[,c("result_id", "docid", 
                        "strat_name_id", "strat_name_long",
                        "unit_id", "t_age", "b_age", 
                        "target_word", "phrase", "type", 
                        "other", "lith", "environ","clat","clng")] 

nodules.2 <- nodules.2[!duplicated(nodules.2[,c("strat_name_id", "unit_id")]),]

# units join

combined.units <- full_join(framboids.2, nodules.2, by = "unit_id")
combined.units <- combined.units[!is.na(combined.units$unit_id),]

combined.units1 <- subset(combined.units, type.x == "Framboids" & type.y == "Nodules")
combined.units1$type <- "Both"

combined.units2 <- subset(combined.units, type.x == "Framboids" & is.na(type.y) == TRUE)
combined.units2$type <- "Framboids"

combined.units3 <- subset(combined.units, type.y == "Nodules" & is.na(type.x) == TRUE)
combined.units3$type <- "Nodules"

combined.units <- rbind(combined.units1, combined.units2, combined.units3)

# strat names join

combined.strat <- full_join(framboids.2, nodules.2, by = "strat_name_id")
combined.strat <- combined.strat[is.na(combined.strat$unit_id.x),]
combined.strat <- combined.strat[is.na(combined.strat$unit_id.y),]

combined.strat1 <- subset(combined.strat, type.x == "Framboids" & type.y == "Nodules")
combined.strat1$type <- "Both"

combined.strat2 <- subset(combined.strat, type.x == "Framboids" & is.na(type.y) == TRUE)
combined.strat2$type <- "Framboids"

combined.strat3 <- subset(combined.strat, type.y == "Nodules" & is.na(type.x) == TRUE)
combined.strat3$type <- "Nodules"

combined.strat <- rbind(combined.strat1, combined.strat2, combined.strat3)

names(combined.units)[3] <- "strat_name_id"
names(combined.strat)[5] <- "unit_id"

combined.units <- combined.units[,-18]
combined.strat <- combined.strat[,-19]

combined <- rbind(combined.units, combined.strat)

combined.x <- combined[!is.na(combined$t_age.x),]

combined.x$t_age <- combined.x$t_age.x
combined.x$b_age <- combined.x$b_age.x
combined.x$strat_name_long <- combined.x$strat_name_long.x
combined.x$lith <- combined.x$lith.x
combined.x$environ <- combined.x$environ.x
combined.x$other <- combined.x$other.x
combined.x$clat <- combined.x$clat.x
combined.x$clng <- combined.x$clng.x

combined.y <- combined[is.na(combined$t_age.x),]

combined.y$t_age <- combined.y$t_age.y
combined.y$b_age <- combined.y$b_age.y
combined.y$strat_name_long <- combined.y$strat_name_long.y
combined.y$lith <- combined.y$lith.y
combined.y$environ <- combined.y$environ.y
combined.y$other <- combined.y$other.y
combined.y$clat <- combined.y$clat.y
combined.y$clng <- combined.y$clng.y

combined <- rbind(combined.x, combined.y)

#combined <- combined[,c(-4,-6,-7,-10,-11,-12,-13,-17,-18,-22,-23,-24)]

ggplot(combined, aes(t_age, reorder(strat_name_long, t_age), group = strat_name_long)) +
  geom_linerange(aes(xmin = b_age, xmax = t_age, colour = type)) + theme_bw() +
  scale_x_reverse(limits = c(3000, 0), breaks = seq(0,3000, by = 100))

# export combined file for gplates figure

write.csv(combined, "pyrite_combined.csv")

## lithology metrics - estimates - Fig. S6A ##

# % pure mudstone

combined.new.units <- rbind(framboids, nodules)
combined.new.units <- combined.new.units[!duplicated(combined.new.units[,c("strat_name_id", "unit_id")]),]
combined.new.units <- combined.new.units[,1:86]

muds <- grepl(paste(mud_rocks, collapse = "|"), combined.new.units$lith, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new.units$environ, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new.units$other, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new.units$strat_name_long, ignore.case=TRUE)

muds.interbeds <- subset(combined.new.units, muds)

muds.interbeds.total <- (100/nrow(combined.new.units))*nrow(muds.interbeds) # % mudstones - approximation includes subordinate beds

manual <- muds.interbeds[sample(nrow(muds.interbeds),53), ] # manual assessment - check for approx. proportion of mudstone-dominated units/packages

others <- grepl(paste(other_rocks, collapse = "|"), muds.interbeds$lith, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$environ, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$other, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$strat_name_long, ignore.case=TRUE)

muds.pure <- subset(muds.interbeds, !others)

muds.total <- (100/nrow(combined.new.units))*nrow(muds.pure) # % mudstones - approximation for pure mudstones

others <- grepl(paste(other_rocks, collapse = "|"), muds.interbeds$lith, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$environ, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$other, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$strat_name_long, ignore.case=TRUE)

muds.sands <- subset(muds.interbeds, others)

muds.sands.total <- (100/nrow(combined.new.units))*nrow(muds.sands) # % mudstones interbedded with sands


# % pure sand

sand <- grepl("sand", combined.new.units$lith, ignore.case=TRUE) |
  grepl("sand", combined.new.units$environ, ignore.case=TRUE) |
  grepl("sand", combined.new.units$other, ignore.case=TRUE) |
  grepl("sand", combined.new.units$strat_name_long, ignore.case=TRUE) |
  grepl("conglomerate", combined.new.units$lith, ignore.case=TRUE) |
  grepl("conglomerate", combined.new.units$environ, ignore.case=TRUE) |
  grepl("conglomerate", combined.new.units$other, ignore.case=TRUE) |
  grepl("conglomerate", combined.new.units$strat_name_long, ignore.case=TRUE) |
  grepl("breccia", combined.new.units$lith, ignore.case=TRUE) |
  grepl("breccia", combined.new.units$environ, ignore.case=TRUE) |
  grepl("breccia", combined.new.units$other, ignore.case=TRUE) |
  grepl("breccia", combined.new.units$strat_name_long, ignore.case=TRUE) |
  grepl("quartzite", combined.new.units$lith, ignore.case=TRUE) |
  grepl("quartzite", combined.new.units$environ, ignore.case=TRUE) |
  grepl("quartzite", combined.new.units$other, ignore.case=TRUE) |
  grepl("quartzite", combined.new.units$strat_name_long, ignore.case=TRUE)

sand.interbeds <- subset(combined.new.units, sand)

sand.interbeds.total <- (100/nrow(combined.new.units))*nrow(sand.interbeds) # % sandstones - approximation includes subordinate beds

muds <- grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$lith, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$environ, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$other, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$strat_name_long, ignore.case=TRUE)

sand.interbeds.not.muds <- subset(sand.interbeds, !muds)

sand.interbeds.not.muds.total <- (100/nrow(combined.new.units))*nrow(sand.interbeds.not.muds) # % sandstones - approximation includes subordinate beds, excluding mudstones


others <- grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$lith, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$environ, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$other, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl("lime", sand.interbeds$lith, ignore.case=TRUE) |
  grepl("lime", sand.interbeds$environ, ignore.case=TRUE) |
  grepl("lime", sand.interbeds$other, ignore.case=TRUE) |
  grepl("lime", sand.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl("basalt", sand.interbeds$lith, ignore.case=TRUE) |
  grepl("basalt", sand.interbeds$environ, ignore.case=TRUE) |
  grepl("basalt", sand.interbeds$other, ignore.case=TRUE) |
  grepl("basalt", sand.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl("dolomite", sand.interbeds$lith, ignore.case=TRUE) |
  grepl("dolomite", sand.interbeds$environ, ignore.case=TRUE) |
  grepl("dolomite", sand.interbeds$other, ignore.case=TRUE) |
  grepl("dolomite", sand.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl("chert", sand.interbeds$lith, ignore.case=TRUE) |
  grepl("chert", sand.interbeds$environ, ignore.case=TRUE) |
  grepl("chert", sand.interbeds$other, ignore.case=TRUE) |
  grepl("chert", sand.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl("volcaniclastic", sand.interbeds$lith, ignore.case=TRUE) |
  grepl("volcaniclastic", sand.interbeds$environ, ignore.case=TRUE) |
  grepl("volcaniclastic", sand.interbeds$other, ignore.case=TRUE) |
  grepl("volcaniclastic", sand.interbeds$strat_name_long, ignore.case=TRUE)



sand.pure <- subset(sand.interbeds, !others)

sand.total <- (100/nrow(combined.new.units))*nrow(sand.pure) # % sandstones - approximation for pure sandstones & conglomerates

# % pure carbonate and dolomite

carbs <- grepl("lime", combined.new.units$lith, ignore.case=TRUE) |
  grepl("lime", combined.new.units$environ, ignore.case=TRUE) |
  grepl("lime", combined.new.units$other, ignore.case=TRUE) |
  grepl("lime", combined.new.units$strat_name_long, ignore.case=TRUE) |
  grepl("dolomite", combined.new.units$lith, ignore.case=TRUE) |
  grepl("dolomite", combined.new.units$environ, ignore.case=TRUE) |
  grepl("dolomite", combined.new.units$other, ignore.case=TRUE) |
  grepl("dolomite", combined.new.units$strat_name_long, ignore.case=TRUE) |
  grepl("chert", combined.new.units$lith, ignore.case=TRUE) |
  grepl("chert", combined.new.units$environ, ignore.case=TRUE) |
  grepl("chert", combined.new.units$other, ignore.case=TRUE) |
  grepl("chert", combined.new.units$strat_name_long, ignore.case=TRUE)

carbs.interbeds <- subset(combined.new.units, carbs)

carbs.interbeds.total <- (100/nrow(combined.new.units))*nrow(carbs.interbeds) # % carbonates - approximation includes subordinate beds

others <- grepl(paste(mud_rocks, collapse = "|"), carbs.interbeds$lith, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), carbs.interbeds$environ, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), carbs.interbeds$other, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), carbs.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl("basalt", carbs.interbeds$lith, ignore.case=TRUE) |
  grepl("basalt", carbs.interbeds$environ, ignore.case=TRUE) |
  grepl("basalt", carbs.interbeds$other, ignore.case=TRUE) |
  grepl("basalt", carbs.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl("volcaniclastic", carbs.interbeds$lith, ignore.case=TRUE) |
  grepl("volcaniclastic", carbs.interbeds$environ, ignore.case=TRUE) |
  grepl("volcaniclastic", carbs.interbeds$other, ignore.case=TRUE) |
  grepl("volcaniclastic", carbs.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl("sand", carbs.interbeds$lith, ignore.case=TRUE) |
  grepl("sand", carbs.interbeds$environ, ignore.case=TRUE) |
  grepl("sand", carbs.interbeds$other, ignore.case=TRUE) |
  grepl("sand", carbs.interbeds$strat_name_long, ignore.case=TRUE)  |
  grepl("conglomerate", carbs.interbeds$lith, ignore.case=TRUE) |
  grepl("conglomerate", carbs.interbeds$environ, ignore.case=TRUE) |
  grepl("conglomerate", carbs.interbeds$other, ignore.case=TRUE) |
  grepl("conglomerate", carbs.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl("breccia", carbs.interbeds$lith, ignore.case=TRUE) |
  grepl("breccia", carbs.interbeds$environ, ignore.case=TRUE) |
  grepl("breccia", carbs.interbeds$other, ignore.case=TRUE) |
  grepl("breccia", carbs.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl("quartzite", carbs.interbeds$lith, ignore.case=TRUE) |
  grepl("quartzite", carbs.interbeds$environ, ignore.case=TRUE) |
  grepl("quartzite", carbs.interbeds$other, ignore.case=TRUE) |
  grepl("quartzite", carbs.interbeds$strat_name_long, ignore.case=TRUE)


carbs.pure <- subset(carbs.interbeds, !others)

carb.total <- (100/nrow(combined.new.units))*nrow(carbs.pure) # % carbonates, dolomites, cherts - approximation for pure

100 - carb.total - sand.total - muds.total

# summary liths

mud.dominant.interbeds <- (muds.interbeds.total-muds.total)*0.72 # 72% based on manual assessment of interbedded lithologies
mud.subordinate.interbeds <- (muds.interbeds.total-muds.total)*0.28 # based on manual assessment of interbedded lithologies

remaining <- 100-sum(muds.total,mud.dominant.interbeds,mud.subordinate.interbeds,sand.total, carb.total)

summary.liths <- as.data.frame(rbind(muds.total,mud.dominant.interbeds,mud.subordinate.interbeds,sand.total, carb.total, remaining))
summary.liths$liths <- rownames(summary.liths)
summary.liths$group <- "A"

ggplot(summary.liths, aes(group, V1, fill = reorder(liths, V1))) + geom_col() # Fig. S5A

# source summary - Fig. S6B

unique(combined.new.units$author)

sources <- table(combined.new.units$author, useNA = "ifany")
sources <- as.data.frame(sources)
sources$group <- "A"

ggplot(sources, aes(group, Freq, fill = reorder(Var1, Freq))) + geom_col(position = "fill")

100/sum(sources$Freq)*127

## counts - Fig. S4 ##

a <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = out, aes(Age-0.5, framboids), colour = "red")

b <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = out, aes(Age-0.5, framboids), colour = "red") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

c <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = out, aes(Age-0.5, nodules), colour = "blue")

d <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = out, aes(Age-0.5, nodules), colour = "blue") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

e <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = out, aes(Age-0.5, undif), colour = "grey70")

f <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = out, aes(Age-0.5, undif), colour = "grey70") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

g <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = out, aes(Age-0.5, seds), colour = "black")

h <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = out, aes(Age-0.5, seds), colour = "black") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


grid.arrange(a,b,c,d,e,f,g,h, ncol = 2)

# 'mineralisation' and evaporite plots - Fig. S7D-E

# run this block for sediment normalisation

evap <- (grepl("evap", data_p2$lith, ignore.case=TRUE) | 
           grepl("evap", data_p2$other, ignore.case=TRUE) | 
           grepl("evap", data_p2$environ, ignore.case=TRUE) |
           grepl("evap", data_p2$strat_phrase_root, ignore.case=TRUE) | 
           grepl("evap", data_p2$phrase, ignore.case=TRUE) | 
           grepl("evap", data_p2$strat_flag, ignore.case=TRUE) |
           grepl("evap", data_p2$environ, ignore.case=TRUE))  &
  (grepl("evap", data_p2$lith, ignore.case=TRUE) | 
     grepl("evap", data_p2$other, ignore.case=TRUE) | 
     grepl("evap", data_p2$strat_phrase_root, ignore.case=TRUE) | 
     grepl("evap", data_p2$strat_flag, ignore.case=TRUE) |
     grepl("evap", data_p2$phrase, ignore.case=TRUE) |
     grepl("evap", data_p2$environ, ignore.case=TRUE))

evap <- subset(data_p2, evap)

evap1 <- subset(evap, t_age2 > cutoff & b_age2 > cutoff)
evap2 <- evap1[!duplicated(evap1$unit_id), ]
evap3 <- evap1[is.na(evap1$unit_id),]
evap2 <- evap2[!is.na(evap2$unit_id),]
evap3 <- evap3[!duplicated(evap3$strat_name_id),]
evap1 <- rbind(evap2, evap3)

evap2 <- subset(evap, t_age2 < cutoff & b_age2 < cutoff)
evap3 <- evap2[!duplicated(evap2$unit_id), ]
evap4 <- evap2[is.na(evap2$unit_id),]
evap3 <- evap3[!is.na(evap3$unit_id),]
evap4 <- evap4[!duplicated(evap4$strat_name_id),]
evap2 <- rbind(evap3, evap4)

evap3 <- subset(evap, t_age2 < cutoff & b_age2 == cutoff)
evap4 <- evap3[!duplicated(evap3$unit_id), ]
evap5 <- evap3[is.na(evap3$unit_id),]
evap4 <- evap4[!is.na(evap4$unit_id),]
evap5 <- evap5[!duplicated(evap5$strat_name_id),]
evap3 <- rbind(evap4, evap5)

evap4 <- subset(evap, t_age2 == cutoff & b_age2 > cutoff)
evap5 <- evap4[!duplicated(evap4$unit_id), ]
evap6 <- evap4[is.na(evap4$unit_id),]
evap5 <- evap5[!is.na(evap5$unit_id),]
evap6 <- evap6[!duplicated(evap6$strat_name_id),]
evap4 <- rbind(evap5, evap6)

evap <- rbind(evap1, evap2, evap3, evap4)

evap <- melt(evap, id.vars = c(1:83,626), na.rm=TRUE)
evap$value <- as.numeric(evap$value)
evap <- subset(evap, variable != "b_age2")

# Phanerozoic bins

evap_bins1 <- hist(evap$value[evap$value >= 0 & evap$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
evap_bins1 <- as.data.frame(cbind(evap_bins1$counts, evap_bins1$breaks))
evap_bins1 <- evap_bins1[1:542,]

# Precambrian bins
evap_bins2 <- hist(evap$value[evap$value >= cutoff & evap$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
evap_bins2 <- as.data.frame(cbind(evap_bins2$counts, evap_bins2$breaks))

# Bind
evap_bins <- rbind(evap_bins1,evap_bins2)
evap_bins <- evap_bins[c(1:541,543:889),] # remove duplicate 541

# block end

veinsS <- as.data.frame(cbind((veins_bins$V1)/(framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1),nodules_bins$V2))
veins2S <- as.data.frame(cbind((veins2_bins$V1)/(framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1),nodules_bins$V2))
evapS <- as.data.frame(cbind((evap_bins$V1)/(sediments_bins$V1),nodules_bins$V2))

# note veins2S (including framboid & nodule-bearing rocks) and veinsS (excluding framboid and nodule-bearing pyritic rocks) are
# produce nearly identitical results. VeinsS is shown on Fig. S7D.

a <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = veins2S, aes(V2-0.5, V1), colour = "red") +
  geom_step(data = veinsS, aes(V2-0.5, V1), colour = "black")

b <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = veins2S, aes(V2-0.5, V1), colour = "red") +
  geom_step(data = veinsS, aes(V2-0.5, V1), colour = "black") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

c <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = evapS, aes(V2-0.5, V1), colour = "blue")

d <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = evapS, aes(V2-0.5, V1), colour = "blue") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

grid.arrange(a,b,c,d, ncol = 2)

# block end

##### Supplementary Materials derived from PART 2 ####

## lithology metrics - estimates - Fig. S6C ##

## utilise units lith descriptions where possible by using dcast

# run PART 2 up to melt of framboids and nodule datasets

unit.liths <- combined.new.units[,c(7,40)]

test <- grepl("Strom", framboids$target_word)
test <- subset(framboids, test)

unit.liths <- cbind(unit.liths, ave(unit.liths$strat_name_id, unit.liths$strat_name_id, FUN=seq_along))
names(unit.liths)[3] <- "ID"
unit.liths$liths <- "lith"
unit.liths$liths <- paste(unit.liths$liths, unit.liths$ID)

unit.liths <- dcast(unit.liths, strat_name_id ~ liths, value.var = "lith")

# % pure mudstone

combined.new <- rbind(framboids, nodules)
combined.new <- combined.new[!duplicated(combined.new[,c("strat_name_id")]),]
combined.new <- combined.new[,1:46]

combined.new <- left_join(combined.new, unit.liths, by = "strat_name_id")

muds <- grepl(paste(mud_rocks, collapse = "|"), combined.new$other, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$strat_name_long, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 1`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 2`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 3`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 4`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 5`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 6`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 7`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 8`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 9`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 10`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 11`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 12`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 13`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 14`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 15`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 16`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 17`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 18`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 19`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 20`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 21`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 22`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), combined.new$`lith 23`, ignore.case=TRUE)


muds.interbeds <- subset(combined.new, muds)

muds.interbeds.total <- (100/nrow(combined.new))*nrow(muds.interbeds) # % mudstones - approximation includes subordinate beds


others <- grepl(paste(other_rocks, collapse = "|"), muds.interbeds$other, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 1`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 2`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 3`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 4`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 5`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 6`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 7`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 8`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 9`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 10`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 11`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 12`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 13`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 14`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 15`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 16`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 17`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 18`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 19`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 20`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 21`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 22`, ignore.case=TRUE) |
  grepl(paste(other_rocks, collapse = "|"), muds.interbeds$`lith 23`, ignore.case=TRUE)

muds.pure <- subset(muds.interbeds, !others)

muds.total <- (100/nrow(combined.new))*nrow(muds.pure) # % mudstones - approximation for pure mudstones

# % pure sand

targets <- c("sand", "conglomerate", "breccia", "quarzite")

sand <- grepl(paste(targets, collapse = "|"), combined.new$other, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$strat_name_long, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 1`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 2`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 3`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 4`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 5`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 6`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 7`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 8`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 9`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 10`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 11`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 12`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 13`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 14`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 15`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 16`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 17`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 18`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 19`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 20`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 21`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 22`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), combined.new$`lith 23`, ignore.case=TRUE)

sand.interbeds <- subset(combined.new, sand)

sand.interbeds.total <- (100/nrow(combined.new))*nrow(sand.interbeds) # % sandstones - approximation includes subordinate beds

muds <- grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$other, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 1`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 2`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 3`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 4`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 5`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 6`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 7`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 8`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 9`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 10`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 11`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 12`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 13`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 14`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 15`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 16`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 17`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 18`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 19`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 20`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 21`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 22`, ignore.case=TRUE) |
  grepl(paste(mud_rocks, collapse = "|"), sand.interbeds$`lith 23`, ignore.case=TRUE)

sand.interbeds.not.muds <- subset(sand.interbeds, !muds)

sand.interbeds.not.muds.total <- (100/nrow(combined.new))*nrow(sand.interbeds.not.muds) # % sandstones - approximation includes subordinate beds, excluding mudstones

other_targets <- c("lime", "volcaniclastic", "chert", "dolomite")

others <- grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$other, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$strat_name_long, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 1`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 2`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 3`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 4`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 5`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 6`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 7`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 8`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 9`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 10`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 11`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 12`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 13`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 14`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 15`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 16`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 17`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 18`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 19`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 20`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 21`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 22`, ignore.case=TRUE) |
  grepl(paste(other_targets, collapse = "|"), sand.interbeds.not.muds$`lith 23`, ignore.case=TRUE)


sand.pure <- subset(sand.interbeds.not.muds, !others)

sand.total <- (100/nrow(combined.new))*nrow(sand.pure) # % sandstones - approximation for pure sandstones & conglomerates

# % pure carbonate and dolomite

bios <- c("lime", "chert", "dolomite")

bios <- grepl(paste(bios, collapse = "|"), combined.new$other, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$strat_name_long, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 1`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 2`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 3`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 4`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 5`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 6`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 7`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 8`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 9`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 10`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 11`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 12`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 13`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 14`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 15`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 16`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 17`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 18`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 19`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 20`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 21`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 22`, ignore.case=TRUE) |
  grepl(paste(bios, collapse = "|"), combined.new$`lith 23`, ignore.case=TRUE)

carbs.interbeds <- subset(combined.new, bios)

carbs.interbeds.total <- (100/nrow(combined.new))*nrow(carbs.interbeds) # % carbonates - approximation includes subordinate beds

targets <- c(mud_rocks, "volcaniclastic", "sand", "conglomerate", "breccia", "quarzite")

others <-  grepl(paste(targets, collapse = "|"), carbs.interbeds$other, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$strat_name_long, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 1`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 2`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 3`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 4`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 5`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 6`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 7`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 8`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 9`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 10`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 11`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 12`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 13`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 14`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 15`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 16`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 17`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 18`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 19`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 20`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 21`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 22`, ignore.case=TRUE) |
  grepl(paste(targets, collapse = "|"), carbs.interbeds$`lith 23`, ignore.case=TRUE)

carbs.pure <- subset(carbs.interbeds, !others)

carb.total <- (100/nrow(combined.new))*nrow(carbs.pure) # % carbonates, dolomites, cherts - approximation for pure

100 - carb.total - sand.total - muds.total

# summary liths

mud.dominant.interbeds <- (muds.interbeds.total-muds.total)*0.72 # based on manual assessment of interbedded mudstone packages
mud.subordinate.interbeds <- (muds.interbeds.total-muds.total)*0.28 # based on manual assessment of interbedded mudstone packages

remaining <- 100-sum(muds.total,mud.dominant.interbeds,mud.subordinate.interbeds,sand.total, carb.total)

summary.liths <- as.data.frame(rbind(muds.total,mud.dominant.interbeds,mud.subordinate.interbeds,sand.total, carb.total, remaining))
summary.liths$liths <- rownames(summary.liths)
summary.liths$group <- "A"

ggplot(summary.liths, aes(group, V1, fill = reorder(liths, V1))) + geom_col()

# source summary - Fig. S6D

unique(combined.new$author)

sources <- table(combined.new$author, useNA = "ifany")
sources <- as.data.frame(sources)
sources$group <- "A"

ggplot(sources, aes(group, Freq, fill = reorder(Var1, Freq))) + geom_col(position = "fill")

(100/sum(sources$Freq))*3

## counts - Fig. S5 ##

out$framboids
out$Age

n <- 2

a <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = out, aes(Age-0.5, framboids), colour = "red")

b <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = out, aes(Age-0.5, framboids), colour = "red") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

c <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = out, aes(Age-0.5, nodules), colour = "blue")

d <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = out, aes(Age-0.5, nodules), colour = "blue") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

e <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = out, aes(Age-0.5, undif), colour = "grey70")

f <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = out, aes(Age-0.5, undif), colour = "grey70") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

g <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = out, aes(Age-0.5, Wilkin.one), colour = "grey70") +
  geom_step(data = out, aes(Age-0.5, Wilkin.three), colour = "grey50") +
  geom_step(data = out, aes(Age-0.5, Wilkin.five), colour = "grey30") +
  geom_step(data = out, aes(Age-0.5, Wilkin.all), colour = "grey10")


h <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = out, aes(Age-0.5, Wilkin.one), colour = "grey70") +
  geom_step(data = out, aes(Age-0.5, Wilkin.three), colour = "grey50") +
  geom_step(data = out, aes(Age-0.5, Wilkin.five), colour = "grey30") +
  geom_step(data = out, aes(Age-0.5, Wilkin.all), colour = "grey10") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

i <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  geom_step(data = out, aes(Age-0.5, seds), colour = "black")

j <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  geom_step(data = out, aes(Age-0.5, seds), colour = "black") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


grid.arrange(a,b,c,d,e,f,g,h,i,j, ncol = 2)

# END #


