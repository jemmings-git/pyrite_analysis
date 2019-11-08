##### This R script was written by Joe Emmings and Jo Walsh (British Geological Survey) ####
##### This script is designed to plot GeoDeepDive extractions
##### interfaced with the Macrostrat database

##### for for further development - 08/11/19 - explore subdividing the analysis into 'shallow water' and 'deep water' components
##### perhaps using environment definitions from
##### envs <- macrostrat_data("envs.json", "https://macrostrat.org/api/defs/environments?all")
##### as the basis?


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


rm(list=ls())

source('macrostrat_data.R')

project_home <- 'N:/Data/xGDD/analysis'
tryCatch({
  setwd(project_home)
}, error = function(err) { 
  if (dir.exists('./data')) {
    setwd('./data') }
}
)

###### PART 1 - using 'output 2' ('p2') - composite analysis of units and strat packages #####


data_p2 <-  read.csv(file = "data_part2_comp.csv", row.names = 1)

# generate complete list of sedimentary units for part 2 - 
# this is done by propagating unit lith info to strat IDs
# so, for any given strat package, if the majority of units within that package
# are 'sedimentary', then the overall package is considered sedimentary
# this info is then joined to the other sediment definitions derived from concepts descriptions

# run this block below - this is required for normalisation to all sediments
# Further subdivision, such as 'marine', etc. is not presently implemented


data_p2$lith2 <- grepl(paste(sedimentary_rocks, collapse = "|"), data_p2$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), data_p2$other, ignore.case=TRUE)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# keep this file for part 2

sed.list <- aggregate(lith2 ~ strat_name_id, data_p2, Mode)

# block end

# optional compilation of target words and summary
phrases <- as.data.frame(as.array(data_p2$target_word))
list <- as.data.frame(summary(data_p2$target_word))

# compile framboid mentions

# block start

framboids <- grepl("framboid", data_p2$target_word, ignore.case=TRUE)

framboids <- subset(data_p2, framboids)

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

framboids1 <- subset(framboids, t_age > cutoff & b_age > cutoff)
framboids2 <- framboids1[!duplicated(framboids1$unit_id), ]
framboids3 <- framboids1[is.na(framboids1$unit_id),]
framboids2 <- framboids2[!is.na(framboids2$unit_id),] 

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), framboids2$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), framboids2$other, ignore.case=TRUE)

framboids2 <- subset(framboids2, sediments) #  subset only sed units
framboids3 <- framboids3[!duplicated(framboids3$strat_name_id),]
framboids1 <- rbind(framboids2, framboids3)

framboids2 <- subset(framboids, t_age < cutoff & b_age < cutoff)
framboids3 <- framboids2[!duplicated(framboids2$unit_id), ]
framboids4 <- framboids2[is.na(framboids2$unit_id),]
framboids3 <- framboids3[!is.na(framboids3$unit_id),]

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), framboids3$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), framboids3$other, ignore.case=TRUE)

framboids3 <- subset(framboids3, sediments) #  subset only sed units
framboids4 <- framboids4[!duplicated(framboids4$strat_name_id),]
framboids2 <- rbind(framboids3, framboids4)

framboids3 <- subset(framboids, t_age < cutoff & b_age == cutoff)
framboids4 <- framboids3[!duplicated(framboids3$unit_id), ]
framboids5 <- framboids3[is.na(framboids3$unit_id),]
framboids4 <- framboids4[!is.na(framboids4$unit_id),]

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), framboids4$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), framboids4$other, ignore.case=TRUE)

framboids4 <- subset(framboids4, sediments) #  subset only sed units
framboids5 <- framboids5[!duplicated(framboids5$strat_name_id),]
framboids3 <- rbind(framboids4, framboids5)

framboids4 <- subset(framboids, t_age == cutoff & b_age > cutoff)
framboids5 <- framboids4[!duplicated(framboids4$unit_id), ]
framboids6 <- framboids4[is.na(framboids4$unit_id),]
framboids5 <- framboids5[!is.na(framboids5$unit_id),]

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), framboids5$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), framboids5$other, ignore.case=TRUE)

framboids5 <- subset(framboids5, sediments) #  subset only sed units
framboids6 <- framboids6[!duplicated(framboids6$strat_name_id),]
framboids4 <- rbind(framboids5, framboids6)

framboids <- rbind(framboids1, framboids2, framboids3, framboids4)

# block end

# repeat the same process for any other 'similar' phrases - in this case 'pyrite nodules'

# block start

nodules <- grepl("nodul", data_p2$target_word, ignore.case=TRUE) | grepl("concretion", data_p2$target_word, ignore.case=TRUE) 

nodules <- subset(data_p2, nodules)

# include nodule mentions pulled from lexicons (concepts)

nodules1 <- grepl("pyrite nodul", data_p2$other, ignore.case=TRUE) | grepl("pyrite concretion", data_p2$other, ignore.case=TRUE) |
  grepl("pyritic nodul", data_p2$other, ignore.case=TRUE) | grepl("pyritic concretion", data_p2$other, ignore.case=TRUE) |
  grepl("nodules of pyrite", data_p2$other, ignore.case=TRUE) | grepl("concretions of pyrite", data_p2$other, ignore.case=TRUE) | 
  grepl("nodular pyrite", data_p2$other, ignore.case=TRUE) | grepl("concretionary pyrite", data_p2$other, ignore.case=TRUE)

nodules1 <- subset(data_p2, nodules1)
nodules <- rbind(nodules, nodules1)

nodules1 <- subset(nodules, t_age > cutoff & b_age > cutoff)
nodules2 <- nodules1[!duplicated(nodules1$unit_id), ]
nodules3 <- nodules1[is.na(nodules1$unit_id),]
nodules2 <- nodules2[!is.na(nodules2$unit_id),]

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), nodules2$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), nodules2$other, ignore.case=TRUE)

nodules2 <- subset(nodules2, sediments) #  subset only sed units
nodules3 <- nodules3[!duplicated(nodules3$strat_name_id),]
nodules1 <- rbind(nodules2, nodules3)

nodules2 <- subset(nodules, t_age < cutoff & b_age < cutoff)
nodules3 <- nodules2[!duplicated(nodules2$unit_id), ]
nodules4 <- nodules2[is.na(nodules2$unit_id),]
nodules3 <- nodules3[!is.na(nodules3$unit_id),]

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), nodules3$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), nodules3$other, ignore.case=TRUE)

nodules3 <- subset(nodules3, sediments) #  subset only sed units
nodules4 <- nodules4[!duplicated(nodules4$strat_name_id),]
nodules2 <- rbind(nodules3, nodules4)

nodules3 <- subset(nodules, t_age < cutoff & b_age == cutoff)
nodules4 <- nodules3[!duplicated(nodules3$unit_id), ]
nodules5 <- nodules3[is.na(nodules3$unit_id),]
nodules4 <- nodules4[!is.na(nodules4$unit_id),]

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), nodules4$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), nodules4$other, ignore.case=TRUE)

nodules4 <- subset(nodules4, sediments) #  subset only sed units
nodules5 <- nodules5[!duplicated(nodules5$strat_name_id),]
nodules3 <- rbind(nodules4, nodules5)

nodules4 <- subset(nodules, t_age == cutoff & b_age > cutoff)
nodules5 <- nodules4[!duplicated(nodules4$unit_id), ]
nodules6 <- nodules4[is.na(nodules4$unit_id),]
nodules5 <- nodules5[!is.na(nodules5$unit_id),]

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), nodules5$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), nodules5$other, ignore.case=TRUE)

nodules5 <- subset(nodules5, sediments) #  subset only sed units
nodules6 <- nodules6[!duplicated(nodules6$strat_name_id),]
nodules4 <- rbind(nodules5, nodules6)

nodules <- rbind(nodules1, nodules2, nodules3, nodules4)

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

all <- rbind(framboids, nodules)

all1 <- all[!duplicated(all$unit_id), ]
all2 <- all[is.na(all$unit_id),]
all1 <- all1[!is.na(all1$unit_id),]
all2 <- all2[!duplicated(all2$strat_name_id),]
all <- rbind(all1, all2)

# right join framboids+nodules to the pyritic strat record

pyritic_strat <- safe_right_join(all, pyritic_strat, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)
pyritic_strat <- pyritic_strat[is.na(pyritic_strat$target_word),]

pyrite_undif <- grepl("\\<pyrite\\>", data_p2$target_word, ignore.case=TRUE) | grepl("\\<pyritic\\>", data_p2$target_word, ignore.case=TRUE)

pyrite_undif <- subset(data_p2, pyrite_undif)

pyrite_undif <- rbind(pyritic_strat, pyrite_undif)

# use anti_join to again carry forward only packages not already identified as containing framboids or nodules

pyrite_undif <- anti_join(pyrite_undif, all, by = c("strat_name_id", "strat_name_id")) 

# subset so only sediments extracted

# run this line if ALL sediments - note slightly different ordering of subset compared to framboids and nodules
# this is because ALL framboids and nodules are assumed to be sedimentary, UNLESS the strat ID is matched to
# a non-sedimentary unit. Whereas undifferentiated pyrite is included ONLY if explicitly
# linked to a sedimentary unit - this is because framboids + nodules form primarily in 
# a 'sedimentary' environment whereas 'pyrite' in general is present in many settings

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), pyrite_undif$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), pyrite_undif$other, ignore.case=TRUE)

pyrite_undif <- subset(pyrite_undif, sediments)

pyrite_undif1 <- subset(pyrite_undif, t_age > cutoff & b_age > cutoff)
pyrite_undif2 <- pyrite_undif1[!duplicated(pyrite_undif1$unit_id), ]
pyrite_undif3 <- pyrite_undif1[is.na(pyrite_undif1$unit_id),]
pyrite_undif2 <- pyrite_undif2[!is.na(pyrite_undif2$unit_id),]
pyrite_undif3 <- pyrite_undif3[!duplicated(pyrite_undif3$strat_name_id),]
pyrite_undif1 <- rbind(pyrite_undif2, pyrite_undif3)

pyrite_undif2 <- subset(pyrite_undif, t_age < cutoff & b_age < cutoff)
pyrite_undif3 <- pyrite_undif2[!duplicated(pyrite_undif2$unit_id), ]
pyrite_undif4 <- pyrite_undif2[is.na(pyrite_undif2$unit_id),]
pyrite_undif3 <- pyrite_undif3[!is.na(pyrite_undif3$unit_id),]
pyrite_undif4 <- pyrite_undif4[!duplicated(pyrite_undif4$strat_name_id),]
pyrite_undif2 <- rbind(pyrite_undif3, pyrite_undif4)

pyrite_undif3 <- subset(pyrite_undif, t_age < cutoff & b_age == cutoff)
pyrite_undif4 <- pyrite_undif3[!duplicated(pyrite_undif3$unit_id), ]
pyrite_undif5 <- pyrite_undif3[is.na(pyrite_undif3$unit_id),]
pyrite_undif4 <- pyrite_undif4[!is.na(pyrite_undif4$unit_id),]
pyrite_undif5 <- pyrite_undif5[!duplicated(pyrite_undif5$strat_name_id),]
pyrite_undif3 <- rbind(pyrite_undif4, pyrite_undif5)

pyrite_undif4 <- subset(pyrite_undif, t_age == cutoff & b_age > cutoff)
pyrite_undif5 <- pyrite_undif4[!duplicated(pyrite_undif4$unit_id), ]
pyrite_undif6 <- pyrite_undif4[is.na(pyrite_undif4$unit_id),]
pyrite_undif5 <- pyrite_undif5[!is.na(pyrite_undif5$unit_id),]
pyrite_undif6 <- pyrite_undif6[!duplicated(pyrite_undif6$strat_name_id),]
pyrite_undif4 <- rbind(pyrite_undif5, pyrite_undif6)

pyrite_undif <- rbind(pyrite_undif1, pyrite_undif2, pyrite_undif3, pyrite_undif4)

# extract nearby mentions of veins or mineralisation

veins <- grepl("vein", pyrite_undif$phrase, ignore.case=TRUE)  |  grepl("mineralisation", pyrite_undif$phrase, ignore.case=TRUE)

veins <- subset(pyrite_undif, veins)

# block end

# the following section melts all datasets (i.e., convert wide to long format)

# block start

# for example framboids

framboids <- melt(framboids, id.vars = c(1:20, 23:83), na.rm=TRUE)
framboids$value <- as.numeric(framboids$value)

# remove b_age - otherwise end up with 2x in the final bin

framboids <- subset(framboids, variable != "b_age")

# repeat for the others

nodules <- melt(nodules, id.vars = c(1:20, 23:83), na.rm=TRUE)
nodules$value <- as.numeric(nodules$value)
nodules <- subset(nodules, variable != "b_age")

pyrite_undif <- melt(pyrite_undif, id.vars = c(1:20, 23:83), na.rm=TRUE)
pyrite_undif$value <- as.numeric(pyrite_undif$value)
pyrite_undif <- subset(pyrite_undif, variable != "b_age")

veins <- melt(veins, id.vars = c(1:20, 23:83), na.rm=TRUE)
veins$value <- as.numeric(veins$value)
veins <- subset(veins, variable != "b_age")

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

framboids_bins <- framboids_bins[!duplicated(framboids_bins$V2),] # remove second 541

# repeat for the others

# nodules

nodules_bins1 <- hist(nodules$value[nodules$value >= 0 & nodules$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
nodules_bins1 <- as.data.frame(cbind(nodules_bins1$counts, nodules_bins1$breaks))
nodules_bins1 <- nodules_bins1[1:542,]

nodules_bins2 <- hist(nodules$value[nodules$value >= cutoff & nodules$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
nodules_bins2 <- as.data.frame(cbind(nodules_bins2$counts, nodules_bins2$breaks))

nodules_bins <- rbind(nodules_bins1,nodules_bins2)
nodules_bins <- nodules_bins[!duplicated(nodules_bins$V2),] # remove second 541

# undifferentiated pyrite mentions

pyrite_undif_bins1 <- hist(pyrite_undif$value[pyrite_undif$value >= 0 & pyrite_undif$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
pyrite_undif_bins1 <- as.data.frame(cbind(pyrite_undif_bins1$counts, pyrite_undif_bins1$breaks))
pyrite_undif_bins1 <- pyrite_undif_bins1[1:542,]

pyrite_undif_bins2 <- hist(pyrite_undif$value[pyrite_undif$value >= cutoff & pyrite_undif$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
pyrite_undif_bins2 <- as.data.frame(cbind(pyrite_undif_bins2$counts, pyrite_undif_bins2$breaks))

pyrite_undif_bins <- rbind(pyrite_undif_bins1,pyrite_undif_bins2)
pyrite_undif_bins <- pyrite_undif_bins[!duplicated(pyrite_undif_bins$V2),] # remove second 541

# undifferentiated pyrite mentions - associated with veins and mineralisation mentions

veins_bins1 <- hist(veins$value[veins$value >= 0 & veins$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
veins_bins1 <- as.data.frame(cbind(veins_bins1$counts, veins_bins1$breaks))
veins_bins1 <- veins_bins1[1:542,]

veins_bins2 <- hist(veins$value[veins$value >= cutoff & veins$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
veins_bins2 <- as.data.frame(cbind(veins_bins2$counts, veins_bins2$breaks))

veins_bins <- rbind(veins_bins1,veins_bins2)
veins_bins <- veins_bins[!duplicated(veins_bins$V2),] # remove second 541

# block end

# now to normalise to sedimentary packages - including nodules and framboids if not identified as 'sedimentary'

# block start

# run this block if ALL sediments - includes framboid + nodule mentions (in order to scale normalisation between 0 and 1)

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), data_p2$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), data_p2$other, ignore.case=TRUE) |
  grepl("framboid", data_p2$target_word, ignore.case=TRUE) | 
  grepl("nodul", data_p2$target_word, ignore.case=TRUE) | 
  grepl("concretion", data_p2$target_word, ignore.case=TRUE) |
  grepl("framboid", data_p2$other, ignore.case=TRUE) | 
  grepl("pyrite nodul", data_p2$other, ignore.case=TRUE) | 
  grepl("pyrite concretion", data_p2$other, ignore.case=TRUE) |
  grepl("pyritic nodul", data_p2$other, ignore.case=TRUE) | 
  grepl("pyritic concretion", data_p2$other, ignore.case=TRUE) |
  grepl("nodules of pyrite", data_p2$other, ignore.case=TRUE) | 
  grepl("concretions of pyrite", data_p2$other, ignore.case=TRUE) | 
  grepl("nodular pyrite", data_p2$other, ignore.case=TRUE) | 
  grepl("concretionary pyrite", data_p2$other, ignore.case=TRUE)


sediments <- subset(data_p2, sediments)

### updated 04/11 to include units

sediments1 <- subset(sediments, t_age > cutoff & b_age > cutoff)
sediments2 <- sediments1[!duplicated(sediments1$unit_id), ]
sediments3 <- sediments1[is.na(sediments1$unit_id),]
sediments2 <- sediments2[!is.na(sediments2$unit_id),]
sediments3 <- sediments3[!duplicated(sediments3$strat_name_id),]
sediments1 <- rbind(sediments2, sediments3)

sediments2 <- subset(sediments, t_age < cutoff & b_age < cutoff)
sediments3 <- sediments2[!duplicated(sediments2$unit_id), ]
sediments4 <- sediments2[is.na(sediments2$unit_id),]
sediments3 <- sediments3[!is.na(sediments3$unit_id),]
sediments4 <- sediments4[!duplicated(sediments4$strat_name_id),]
sediments2 <- rbind(sediments3, sediments4)

sediments3 <- subset(sediments, t_age < cutoff & b_age == cutoff)
sediments4 <- sediments3[!duplicated(sediments3$unit_id), ]
sediments5 <- sediments3[is.na(sediments3$unit_id),]
sediments4 <- sediments4[!is.na(sediments4$unit_id),]
sediments5 <- sediments5[!duplicated(sediments5$strat_name_id),]
sediments3 <- rbind(sediments4, sediments5)

sediments4 <- subset(sediments, t_age == cutoff & b_age > cutoff)
sediments5 <- sediments4[!duplicated(sediments4$unit_id), ]
sediments6 <- sediments4[is.na(sediments4$unit_id),]
sediments5 <- sediments5[!is.na(sediments5$unit_id),]
sediments6 <- sediments6[!duplicated(sediments6$strat_name_id),]
sediments4 <- rbind(sediments5, sediments6)

sediments <- rbind(sediments1, sediments2, sediments3, sediments4)

sediments <- melt(sediments, id.vars = c(1:20, 23:83), na.rm=TRUE)
sediments$value <- as.numeric(sediments$value)

sediments <- subset(sediments, variable != "b_age")

# Phanerozoic bins

sediments_bins1 <- hist(sediments$value[sediments$value >= 0 & sediments$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
sediments_bins1 <- as.data.frame(cbind(sediments_bins1$counts, sediments_bins1$breaks))
sediments_bins1 <- sediments_bins1[1:542,]

# Precambrian bins
sediments_bins2 <- hist(sediments$value[sediments$value >= cutoff & sediments$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
sediments_bins2 <- as.data.frame(cbind(sediments_bins2$counts, sediments_bins2$breaks))

# Bind
sediments_bins <- rbind(sediments_bins1,sediments_bins2)
sediments_bins <- sediments_bins[!duplicated(sediments_bins$V2),] # remove second 541

# block end

# now for normalisation to sedimentary bins, as stacked output

# block start

# remove vein counts from pyrite_undif

pyrite_undif_bins$V1 <- pyrite_undif_bins$V1-veins_bins$V1

# ratios

framboidsR <- as.data.frame(cbind((framboids_bins$V1+veins_bins$V1)/sediments_bins$V1,nodules_bins$V2))
nodulesR <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1+veins_bins$V1)/sediments_bins$V1,nodules_bins$V2))
veinsR <- as.data.frame(cbind((veins_bins$V1)/sediments_bins$V1,veins_bins$V2))
all <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1)/sediments_bins$V1,nodules_bins$V2))
ratio <- as.data.frame(cbind((framboidsR$V1)/(framboidsR$V1+nodulesR$V1),nodules_bins$V2))
ratio$V1[is.nan(ratio$V1)] <- NA
frams_nods <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1)/sediments_bins$V1,nodules_bins$V2))


# optional generation of additional age points for plotting

ME <- c(445, 372, 252, 201, 65)
OAEs <- c(183, 120, 111, 93) # from Jenkyns 2010
PETM <- 55.8 # from Jenkyns 2010
Sturt <- c(716, 663) # Sturtian glaciation from ???
Es <- c(517, 502, 405, 393, 388, 382, 359, 330, 249, 240, 230, 220, 188, 145)
Extras <- c(2500, 1600, 1000, 720, 635, 541, 485.4, 443.8, 419.2, 358.9, 298.9, 
            251.9, 201.3, 145, 66, 23.03, 2.58) # Chronostrat divisions
supercontinents <- c(320, 170, 900, 700, 1600, 1400) # from Li et al. 2019 Precambrian Research
carb.intervals <- c(541, 465, 372, 323, 265, 227, 164, 133)
lows <- c(323, 299, 201, 170) # carb.intervals and lows from Riding et al. 2019

# block end

# plot the results

# block start

Top <- -0 # top age for plot (in Ma), set at -0.5 in order to centre bins
Bottom <- 3500 # bottom age for plot (in Ma)
phanerozoic_increment <- 1
precambrian_increment <- 10


a <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  scale_y_continuous(limits = c(0,0.4)) +
  geom_stepribbon(data = all, aes(V2-0.5, ymin = 0, ymax = V1), fill = "grey") +
  geom_stepribbon(data = nodulesR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "blue") +
  geom_stepribbon(data = framboidsR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "red") +
  geom_stepribbon(data = veinsR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "green") +
  geom_step(data = all, aes(V2-0.5, V1)) +
  geom_step(data = nodulesR, aes(V2-0.5, V1)) + 
  geom_step(data = veinsR, aes(V2-0.5, V1)) + 
  geom_step(data = framboidsR, aes(V2-0.5, V1)) #+ 
#geom_vline(xintercept = c(Extras))

b <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  scale_y_continuous(limits = c(0,0.4)) +
  geom_stepribbon(data = all, aes(V2-0.5, ymin = 0, ymax = V1), fill = "grey") +
  geom_stepribbon(data = nodulesR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "blue") +
  geom_stepribbon(data = framboidsR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "red") +
  geom_stepribbon(data = veinsR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "green") +
  geom_step(data = framboidsR, aes(V2-0.5, V1)) + 
  geom_step(data = nodulesR, aes(V2-0.5, V1)) + 
  geom_step(data = veinsR, aes(V2-0.5, V1)) + 
  geom_step(data = all, aes(V2-0.5, V1)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) #+ 
#geom_vline(xintercept = ME) + 
#geom_vline(xintercept = c(Es, OAEs))

#Extras, ME, OAEs, 

# generate plots

grid.arrange(a,b, ncol = 2)

# block end

##### Part 1 follow-on - fragmentation index and convex hull facets for each period ####
## in development 05/11/19 ##

# import geological period info

periods <- read.delim("periods.txt")

periods<- periods[seq(dim(periods)[1],1),]

breaks <- periods$Breaks

periods$Period <- factor(periods$Period, levels=unique(periods$Period))
periods <- periods$Period

# import zaffos file - https://github.com/UW-Macrostrat/PNAS_201702297/blob/master/FinalData/ContinuousTimeSeries.csv

zaffos <- read.delim("Zaffos_et_al.txt", sep=",")

zaffos <- zaffos[,1:2]

names(zaffos)<- c("V2", "fragmentation")

names(veinsR)<- c("veinsR", "V2")
names(nodulesR)<- c("nodulesR", "V2")
names(framboidsR)<- c("framboidsR", "V2")
names(all)<- c("all", "V2")
names(frams_nods)<- c("frams_nods", "V2")

zaffos <- left_join(veinsR, zaffos, by = c("V2", "V2"))
zaffos <- left_join(nodulesR, zaffos, by = c("V2", "V2"))
zaffos <- left_join(framboidsR, zaffos, by = c("V2", "V2"))
zaffos <- left_join(all, zaffos, by = c("V2", "V2"))
zaffos <- left_join(frams_nods, zaffos, by = c("V2", "V2"))

test <- zaffos$fragmentation

# plot - not sure if it is worth further investigating? doesn't appear to be any key relationship

ggplot(zaffos) + geom_line(aes(V2, fragmentation)) +
  stat_smooth(aes(V2, fragmentation), colour = "black", method=lm, formula = y ~ poly(x,6)) +
  theme_bw() + scale_x_continuous(limits = c(0,541))

# one option might be to explore D() i.e., differentiating the fragmentation curve

# superimposed onto strat_only plot 

p <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  scale_y_continuous(limits = c(0,0.3)) +
  geom_stepribbon(data = all, aes(V2-0.5, ymin = 0, ymax = all), fill = "grey") +
  geom_stepribbon(data = nodulesR, aes(V2-0.5, ymin = 0, ymax = nodulesR), fill = "blue") +
  geom_stepribbon(data = framboidsR, aes(V2-0.5, ymin = 0, ymax = framboidsR), fill = "red") +
  geom_stepribbon(data = veinsR, aes(V2-0.5, ymin = 0, ymax = veinsR), fill = "green") +
  geom_step(data = framboidsR, aes(V2-0.5, framboidsR)) + 
  geom_step(data = nodulesR, aes(V2-0.5, nodulesR)) + 
  geom_step(data = veinsR, aes(V2-0.5, veinsR)) + 
  geom_step(data = all, aes(V2-0.5, all)) +
  geom_line(data = zaffos, aes(V2, fragmentation/3)) +
  scale_y_continuous(sec.axis = sec_axis(~.), limits = c(0, 0.35))

print(p)

# facets by period - in this example plotting the ratio of framboids/nodules versus pyritic vs. non-pyritic packages and units

pyrites <- rbind(framboids, nodules, pyrite_undif)

length(unique(pyrites$docid)) # metrics - no. of documents
length(unique(pyrites$strat_name_id)) # metrics - no. of strat packages
length(unique(pyrites$unit_id)) # metrics - no. of units
length(unique(pyrites$result_id)) # metrics - no. of target phrases

pyrites1 <- subset(pyrites, value > cutoff)
pyrites2 <- pyrites1[!duplicated(pyrites1$unit_id), ]
pyrites3 <- pyrites1[is.na(pyrites1$unit_id),]
pyrites2 <- pyrites2[!is.na(pyrites2$unit_id),]
pyrites3 <- pyrites3[!duplicated(pyrites3$strat_name_id),]
pyrites1 <- rbind(pyrites2, pyrites3)

pyrites2 <- subset(pyrites, value < cutoff)
pyrites3 <- pyrites2[!duplicated(pyrites2$unit_id), ]
pyrites4 <- pyrites2[is.na(pyrites2$unit_id),]
pyrites3 <- pyrites3[!is.na(pyrites3$unit_id),]
pyrites4 <- pyrites4[!duplicated(pyrites4$strat_name_id),]
pyrites2 <- rbind(pyrites3, pyrites4)

pyrites <- rbind(pyrites1, pyrites2)

# fix for Phanerozoic vs. Precambrian below

pyrites1 <- subset(pyrites, value <= 541)
pyrites1$V2 <- round(pyrites1$value, digits = 0)

pyrites2 <- subset(pyrites, value > 541)
pyrites2$value <- pyrites2$value-1
pyrites2$V2 <- round(pyrites2$value, digits = -1)
pyrites2$value <- pyrites2$value+1
pyrites2$V2 <- pyrites2$V2+1

pyrites <- rbind(pyrites1, pyrites2)

# use frams_nods from this point onwards (i.e, framboids+nodules/ sed packages)
# or 'all' is (framboids+nodules+pyrite undif/sed packages)

test <- left_join(pyrites, ratio, by = c("V2", "V2"))
test <- left_join(test, all, by = c("V2", "V2"))
test <- left_join(test, frams_nods, by = c("V2", "V2"))

test$period <- cut(test$V2, 
                   breaks=c(-Inf, 2500, 541, 485, 444, 419, 359, 299, 252, 201, 145, 66, 23, 3, Inf),
                   labels=periods)

#unique(test$econ)

test <- test[!is.infinite(test$V1),]
test <- test[!is.na(test$V1),]
test <- test[!is.infinite(test$frams_nods),]
test <- test[!is.na(test$frams_nods),]

hull_test <- test %>%
  group_by(strat_name_id) %>%
  slice(chull(V1, frams_nods))

# optional removal of quaternary & archaean

#hull_test <- hull_test[which(hull_test$period != "Quaternary"),]
#hull_test <- hull_test[which(hull_test$period != "Archaean"),]

# x axis - ratio of framboids:nodules, y-axis proportion of framboids+nodules normalised to sed. packages

ggplot(hull_test, aes(V1, frams_nods)) + geom_polygon(aes(group = strat_name_id, fill = period), alpha = 0.5) +
  theme_bw() + 
  #scale_y_log10(limits = c(0.003, 0.3)) +
  #scale_x_continuous(limits = c(0.1, 0.5)) +
  geom_polygon(aes(group = strat_name_id), colour = "black", alpha = 0) +
  facet_wrap(~period) + geom_hline(yintercept = c(0.02, 0.06)) +
  geom_vline(xintercept = 0.35) +
  geom_point(aes(fill = period), colour = "black", pch = 21)

# or for y axis as all pyrite mentions (including undifferentiated)

test <- test[!is.infinite(test$V1),]
test <- test[!is.na(test$V1),]
test <- test[!is.infinite(test$all),]
test <- test[!is.na(test$all),]

hull_test <- test %>%
  group_by(strat_name_id) %>%
  slice(chull(V1, all))

# optional removal of quaternary

#hull_test <- hull_test[which(hull_test$period != "Quaternary"),]
#hull_test <- hull_test[which(hull_test$period != "Archaean"),]

# x axis - ratio of framboids:nodules, y-axis proportion of framboids+nodules normalised to sed. packages

ggplot(hull_test, aes(V1, all)) + geom_polygon(aes(group = strat_name_id, fill = period), alpha = 0.5) +
  theme_bw() + 
  #scale_y_log10(limits = c(0.003, 0.3)) +
  #scale_x_continuous(limits = c(0.1, 0.5)) +
  geom_polygon(aes(group = strat_name_id), colour = "black", alpha = 0) +
  facet_wrap(~period) + geom_hline(yintercept = c(0.15, 0.25)) +
  geom_vline(xintercept = 0.35) +
  geom_point(aes(fill = period), colour = "black", pch = 21)


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
list <- as.data.frame(summary(data_p1$target_word))

# compile framboid mentions

framboids <- grepl("framboid", data_p1$target_word, ignore.case=TRUE)

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

framboids1 <- subset(framboids, t_age > cutoff & b_age > cutoff)
framboids1 <- framboids1[!duplicated(framboids1$strat_name_id),]

framboids2 <- subset(framboids, t_age < cutoff & b_age < cutoff)
framboids2 <- framboids2[!duplicated(framboids2$strat_name_id),]

framboids3 <- subset(framboids, t_age < cutoff & b_age == cutoff)
framboids3 <- framboids3[!duplicated(framboids3$strat_name_id),]

framboids4 <- subset(framboids, t_age == cutoff & b_age > cutoff)
framboids4 <- framboids4[!duplicated(framboids4$strat_name_id),]

framboids <- rbind(framboids1, framboids2, framboids3, framboids4)

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

nodules1 <- subset(nodules, t_age > cutoff & b_age > cutoff)
nodules1 <- nodules1[!duplicated(nodules1$strat_name_id),]

nodules2 <- subset(nodules, t_age < cutoff & b_age < cutoff)
nodules2 <- nodules2[!duplicated(nodules2$strat_name_id),]

nodules3 <- subset(nodules, t_age < cutoff & b_age == cutoff)
nodules3 <- nodules3[!duplicated(nodules3$strat_name_id),]

nodules4 <- subset(nodules, t_age == cutoff & b_age > cutoff)
nodules4 <- nodules4[!duplicated(nodules4$strat_name_id),]

nodules <- rbind(nodules1, nodules2, nodules3, nodules4)


# block end

# now to process undifferentiated pyrite mentions

# block start

# pyritic strat concepts pulled from lexicons (undifferentiated)
pyritic_strat <- grepl("pyrit", data_p1$other, ignore.case=TRUE)
pyritic_strat <- subset(data_p1, pyritic_strat)
false <- !grepl("non-pyrit", pyritic_strat$other, ignore.case=TRUE)
pyritic_strat <- subset(pyritic_strat, false)

all <- rbind(framboids, nodules)

all <- all[!duplicated(all$strat_name_id),]

# 01/11/19 note -
# something which could be done, would be to upgrade pyritic units to strat packages using Mode
# in the same way as 'sediments' via the use of sed.list..

pyritic_strat <- safe_right_join(all, pyritic_strat, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)
pyritic_strat <- pyritic_strat[is.na(pyritic_strat$target_word),]

pyrite <- grepl("\\<pyrite\\>", data_p1$target_word, ignore.case=TRUE) | grepl("\\<pyritic\\>", data_p1$target_word, ignore.case=TRUE)

pyrite_undif <- subset(data_p1, pyrite)

pyrite_undif <- rbind(pyritic_strat, pyrite_undif)

# subset so only sediments extracted - use sed.list from part 3 - based on unit lith descriptions
# sed.list is presently defined as all sediments

pyrite_undif <- safe_right_join(sed.list, pyrite_undif, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)

# in theory lith2 is the only requirement here (if defined by ALL sediments), but other terms added as a precaution

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), pyrite_undif$other, ignore.case=TRUE) | 
  grepl("TRUE", pyrite_undif$lith2) 

pyrite_undif <- subset(pyrite_undif, sediments)

##### carry forward only strat packages not already containing framboid or nodule phrases

pyrite_undif <- anti_join(pyrite_undif, all, by = c("strat_name_id", "strat_name_id")) 

# drop lith2 as not needed from this point onwards

pyrite_undif <- select(pyrite_undif, -c(lith2))

pyrite_undif1 <- subset(pyrite_undif, t_age > cutoff & b_age > cutoff)
pyrite_undif1 <- pyrite_undif1[!duplicated(pyrite_undif1$strat_name_id),]

pyrite_undif2 <- subset(pyrite_undif, t_age < cutoff & b_age < cutoff)
pyrite_undif2 <- pyrite_undif2[!duplicated(pyrite_undif2$strat_name_id),]

pyrite_undif3 <- subset(pyrite_undif, t_age < cutoff & b_age == cutoff)
pyrite_undif3 <- pyrite_undif3[!duplicated(pyrite_undif3$strat_name_id),]

pyrite_undif4 <- subset(pyrite_undif, t_age == cutoff & b_age > cutoff)
pyrite_undif4 <- pyrite_undif4[!duplicated(pyrite_undif4$strat_name_id),]

pyrite_undif <- rbind(pyrite_undif1, pyrite_undif2, pyrite_undif3, pyrite_undif4)

# veins/mineralisation extraction as in part 1

veins <- grepl("vein", pyrite_undif$phrase, ignore.case=TRUE)  |  
  grepl("mineralisation", pyrite_undif$phrase, ignore.case=TRUE)

veins <- subset(pyrite_undif, veins)

# block end

# block start

# melt datasets - the same workflow as part 1, but with different column ranges

framboids <- melt(framboids, id.vars = c(1:26, 29:44), na.rm=TRUE)
framboids$value <- as.numeric(framboids$value)
framboids <- subset(framboids, variable != "b_age")

nodules <- melt(nodules, id.vars = c(1:26, 29:44), na.rm=TRUE)
nodules$value <- as.numeric(nodules$value)
nodules <- subset(nodules, variable != "b_age")

pyrite_undif <- melt(pyrite_undif, id.vars = c(1:26, 29:44), na.rm=TRUE)
pyrite_undif$value <- as.numeric(pyrite_undif$value)
pyrite_undif <- subset(pyrite_undif, variable != "b_age")

veins <- melt(veins, id.vars = c(1:26, 29:44), na.rm=TRUE)
veins$value <- as.numeric(veins$value)
veins <- subset(veins, variable != "b_age")

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

framboids_bins <- framboids_bins[!duplicated(framboids_bins$V2),] # remove second 541

# Phanerozoic bins
nodules_bins1 <- hist(nodules$value[nodules$value >= 0 & nodules$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
nodules_bins1 <- as.data.frame(cbind(nodules_bins1$counts, nodules_bins1$breaks))
nodules_bins1 <- nodules_bins1[1:542,]

# Precambrian bins
nodules_bins2 <- hist(nodules$value[nodules$value >= cutoff & nodules$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
nodules_bins2 <- as.data.frame(cbind(nodules_bins2$counts, nodules_bins2$breaks))

# Bind
nodules_bins <- rbind(nodules_bins1,nodules_bins2)

nodules_bins <- nodules_bins[!duplicated(nodules_bins$V2),] # remove second 541

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

pyrite_undif_bins <- pyrite_undif_bins[!duplicated(pyrite_undif_bins$V2),] # remove second 541

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

veins_bins <- veins_bins[!duplicated(veins_bins$V2),] # remove second 541

# block end

# now to normalise to sedimentary packages - including nodules and framboids if not identified as 'sedimentary'

# block start

# use sed info from lith for N. American units, upgraded to strat package level, using Mode as defined in part 1

data.sed <- safe_right_join(sed.list, data_p1, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)

sediments <- grepl(paste(sedimentary_rocks, collapse = "|"), data.sed$lith, ignore.case=TRUE) | 
  grepl(paste(sedimentary_rocks, collapse = "|"), data.sed$other, ignore.case=TRUE) |
  grepl("framboid", data.sed$target_word, ignore.case=TRUE) | 
  grepl("nodul", data.sed$target_word, ignore.case=TRUE) | 
  grepl("concretion", data.sed$target_word, ignore.case=TRUE) |
  grepl("framboid", data.sed$other, ignore.case=TRUE) | 
  grepl("pyrite nodul", data.sed$other, ignore.case=TRUE) | 
  grepl("pyritic nodul", data.sed$other, ignore.case=TRUE) | 
  grepl("nodules of pyrite", data.sed$other, ignore.case=TRUE) | 
  grepl("pyrite concretion", data.sed$other, ignore.case=TRUE) |
  grepl("concretionary pyrite", data.sed$other, ignore.case=TRUE) |
  grepl("concretions of pyrite", data.sed$other, ignore.case=TRUE) |
  grepl("TRUE", data.sed$lith2) 

sediments <- subset(data.sed, sediments)

sediments <- select(sediments, -c(lith2))

sediments1 <- subset(sediments, t_age > cutoff & b_age > cutoff)
sediments1 <- sediments1[!duplicated(sediments1$strat_name_id),]

sediments2 <- subset(sediments, t_age < cutoff & b_age < cutoff)
sediments2 <- sediments2[!duplicated(sediments2$strat_name_id),]

sediments3 <- subset(sediments, t_age < cutoff & b_age == cutoff)
sediments3 <- sediments3[!duplicated(sediments3$strat_name_id),]

sediments4 <- subset(sediments, t_age == cutoff & b_age > cutoff)
sediments4 <- sediments4[!duplicated(sediments4$strat_name_id),]

sediments <- rbind(sediments1, sediments2, sediments3, sediments4)

sediments <- melt(sediments, id.vars = c(1:26, 29:44), na.rm=TRUE)
sediments$value <- as.numeric(sediments$value)

sediments <- subset(sediments, variable != "b_age")

sediments_bins1 <- hist(sediments$value[sediments$value >= 0 & sediments$value < cutoff+1], breaks = seq(0, cutoff+1, by = phanerozoic_increment))
sediments_bins1 <- as.data.frame(cbind(sediments_bins1$counts, sediments_bins1$breaks))
sediments_bins1 <- sediments_bins1[1:542,]

# Precambrian bins
sediments_bins2 <- hist(sediments$value[sediments$value >= cutoff & sediments$value <= 4001], breaks = seq(cutoff, 4001, by = precambrian_increment))
sediments_bins2 <- as.data.frame(cbind(sediments_bins2$counts, sediments_bins2$breaks))

# Bind
sediments_bins <- rbind(sediments_bins1,sediments_bins2)

sediments_bins <- sediments_bins[!duplicated(sediments_bins$V2),] # remove second 541

# block end

# normalisation as in part 1

# block start

# remove vein counts from pyrite_undif

pyrite_undif_bins$V1 <- pyrite_undif_bins$V1-veins_bins$V1

# ratios

framboidsR <- as.data.frame(cbind((framboids_bins$V1+veins_bins$V1)/sediments_bins$V1,nodules_bins$V2))
nodulesR <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1+veins_bins$V1)/sediments_bins$V1,nodules_bins$V2))
veinsR <- as.data.frame(cbind((veins_bins$V1)/sediments_bins$V1,veins_bins$V2))
all <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1+pyrite_undif_bins$V1)/sediments_bins$V1,nodules_bins$V2))
ratio <- as.data.frame(cbind((framboidsR$V1)/(framboidsR$V1+nodulesR$V1),nodules_bins$V2))
ratio$V1[is.nan(ratio$V1)] <- NA
frams_nods <- as.data.frame(cbind((framboids_bins$V1+nodules_bins$V1)/sediments_bins$V1,nodules_bins$V2))


# optional additional time points for plotting

ME <- c(445, 372, 252, 201, 65)
OAEs <- c(183, 120, 111, 93) # from Jenkyns 2010
PETM <- 55.8 # from Jenkyns 2010
Sturt <- c(716, 663) # Sturtian glaciation from ???
Es <- c(517, 502, 405, 393, 388, 382, 359, 330, 249, 240, 230, 220, 188, 145)
Extras <- c(2500, 1600, 1000, 720, 635, 541, 485.4, 443.8, 419.2, 358.9, 298.9, 
            251.9, 201.3, 145, 66, 23.03, 2.58) # Chronostrat divisions
supercontinents <- c(320, 170, 900, 700, 1600, 1400) # from Li et al. 2019 Precambrian Research
carb.intervals <- c(541, 465, 372, 323, 265, 227, 164, 133)
lows <- c(323, 299, 201, 170) # carb.intervals and lows from Riding et al. 2019

# plots

Top <- -0 # top age for plot (in Ma), set at -0.5 in order to centre bins
Bottom <- 3500 # bottom age for plot (in Ma)
phanerozoic_increment <- 1
precambrian_increment <- 10


a <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(Bottom, 541)) +
  scale_y_continuous(limits = c(0,0.2)) +
  geom_stepribbon(data = all, aes(V2-0.5, ymin = 0, ymax = V1), fill = "grey") +
  geom_stepribbon(data = nodulesR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "blue") +
  geom_stepribbon(data = framboidsR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "red") +
  geom_stepribbon(data = veinsR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "green") +
  geom_step(data = all, aes(V2-0.5, V1)) +
  geom_step(data = nodulesR, aes(V2-0.5, V1)) + 
  geom_step(data = veinsR, aes(V2-0.5, V1)) + 
  geom_step(data = framboidsR, aes(V2-0.5, V1)) #+ 
#geom_vline(xintercept = c(Extras))

b <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  scale_y_continuous(limits = c(0,0.2)) +
  geom_stepribbon(data = all, aes(V2-0.5, ymin = 0, ymax = V1), fill = "grey") +
  geom_stepribbon(data = nodulesR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "blue") +
  geom_stepribbon(data = framboidsR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "red") +
  geom_stepribbon(data = veinsR, aes(V2-0.5, ymin = 0, ymax = V1), fill = "green") +
  geom_step(data = framboidsR, aes(V2-0.5, V1)) + 
  geom_step(data = nodulesR, aes(V2-0.5, V1)) + 
  geom_step(data = veinsR, aes(V2-0.5, V1)) + 
  geom_step(data = all, aes(V2-0.5, V1)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) #+ 
#geom_vline(xintercept = ME) #+ 
#geom_vline(xintercept = c(Es, OAEs))

#Extras, ME, OAEs, 

# generate plots

b <- b + geom_vline(xintercept = c(Extras))

grid.arrange(a,b, ncol = 2)

# block end

##### Part 2 follow-on - fragmentation index and convex hull facets for each period ####
## in development 05/11/19 ##

# import geological period info

periods <- read.delim("periods.txt")

periods<- periods[seq(dim(periods)[1],1),]

breaks <- periods$Breaks

periods$Period <- factor(periods$Period, levels=unique(periods$Period))
periods <- periods$Period

# import zaffos file - https://github.com/UW-Macrostrat/PNAS_201702297/blob/master/FinalData/ContinuousTimeSeries.csv

zaffos <- read.delim("Zaffos_et_al.txt", sep=",")

zaffos <- zaffos[,1:2]

names(zaffos)<- c("V2", "fragmentation")

names(veinsR)<- c("veinsR", "V2")
names(nodulesR)<- c("nodulesR", "V2")
names(framboidsR)<- c("framboidsR", "V2")
names(all)<- c("all", "V2")
names(frams_nods)<- c("frams_nods", "V2")


zaffos <- left_join(veinsR, zaffos, by = c("V2", "V2"))
zaffos <- left_join(nodulesR, zaffos, by = c("V2", "V2"))
zaffos <- left_join(framboidsR, zaffos, by = c("V2", "V2"))
zaffos <- left_join(all, zaffos, by = c("V2", "V2"))
zaffos <- left_join(frams_nods, zaffos, by = c("V2", "V2"))

test <- zaffos$fragmentation

# plot - not sure if it is worth further investigating? doesn't appear to be any key relationship

ggplot(zaffos) + geom_line(aes(V2, fragmentation)) +
  stat_smooth(aes(V2, fragmentation), colour = "black", method=lm, formula = y ~ poly(x,6)) +
  theme_bw() + scale_x_continuous(limits = c(0,541))

# one option might be to explore D() i.e., differentiating the fragmentation curve

# superimposed onto strat_only plot 

p <- ggplot() + theme_bw() +
  scale_x_reverse(limits = c(541, Top)) +
  scale_y_continuous(limits = c(0,0.3)) +
  geom_stepribbon(data = all, aes(V2-0.5, ymin = 0, ymax = all), fill = "grey") +
  geom_stepribbon(data = nodulesR, aes(V2-0.5, ymin = 0, ymax = nodulesR), fill = "blue") +
  geom_stepribbon(data = framboidsR, aes(V2-0.5, ymin = 0, ymax = framboidsR), fill = "red") +
  geom_stepribbon(data = veinsR, aes(V2-0.5, ymin = 0, ymax = veinsR), fill = "green") +
  geom_step(data = framboidsR, aes(V2-0.5, framboidsR)) + 
  geom_step(data = nodulesR, aes(V2-0.5, nodulesR)) + 
  geom_step(data = veinsR, aes(V2-0.5, veinsR)) + 
  geom_step(data = all, aes(V2-0.5, all)) +
  geom_line(data = zaffos, aes(V2, fragmentation/3)) +
  scale_y_continuous(sec.axis = sec_axis(~.), limits = c(0, 0.35))

print(p)

# facets by period - in this example plotting the ratio of framboids/nodules versus pyritic vs. non-pyritic packages and units

pyrites <- rbind(framboids, nodules, pyrite_undif)

pyrites1 <- subset(pyrites, value > cutoff)
pyrites1 <- pyrites1[!duplicated(pyrites1$strat_name_id),]

pyrites2 <- pyrites2[!duplicated(pyrites2$strat_name_id),]

pyrites <- rbind(pyrites1, pyrites2)

# fix for Phanerozoic vs. Precambrian below

pyrites1 <- subset(pyrites, value <= 541)
pyrites1$V2 <- round(pyrites1$value, digits = 0)

pyrites2 <- subset(pyrites, value > 541)
pyrites2$value <- pyrites2$value-1
pyrites2$V2 <- round(pyrites2$value, digits = -1)
pyrites2$value <- pyrites2$value+1
pyrites2$V2 <- pyrites2$V2+1

pyrites <- rbind(pyrites1, pyrites2)

test <- left_join(pyrites, ratio, by = c("V2", "V2"))
test <- left_join(test, all, by = c("V2", "V2"))
test <- left_join(test, frams_nods, by = c("V2", "V2"))

test$period <- cut(test$V2, 
                   breaks=c(-Inf, 2500, 541, 485, 444, 419, 359, 299, 252, 201, 145, 66, 23, 3, Inf),
                   labels=periods)

#unique(test$econ)

# use frams_nods from this point onwards (i.e, framboids+nodules/ sed packages)
# or 'all' is (framboids+nodules+pyrite undif/sed packages)

library(tidyverse)

test <- test[!is.infinite(test$V1),]
test <- test[!is.na(test$V1),]
test <- test[!is.infinite(test$frams_nods),]
test <- test[!is.na(test$frams_nods),]

hull_test <- test %>%
  group_by(strat_name_id) %>%
  slice(chull(V1, frams_nods))

# optional removal of quaternary & archaean

#hull_test <- hull_test[which(hull_test$period != "Quaternary"),]
#hull_test <- hull_test[which(hull_test$period != "Archaean"),]

ggplot(hull_test, aes(V1, frams_nods)) + geom_polygon(aes(group = strat_name_id, fill = period), alpha = 0.5) +
  theme_bw() + 
  #scale_y_log10(limits = c(0.003, 0.3)) +
  #scale_x_continuous(limits = c(0.1, 0.5)) +
  geom_polygon(aes(group = strat_name_id), colour = "black", alpha = 0) +
  facet_wrap(~period) + geom_hline(yintercept = c(0.02, 0.06)) +
  geom_vline(xintercept = 0.35) +
  geom_point(aes(fill = period), colour = "black", pch = 21)

# or for all pyrite

test <- test[!is.infinite(test$V1),]
test <- test[!is.na(test$V1),]
test <- test[!is.infinite(test$all),]
test <- test[!is.na(test$all),]

hull_test <- test %>%
  group_by(strat_name_id) %>%
  slice(chull(V1, all))

# optional removal of quaternary & archaean

#hull_test <- hull_test[which(hull_test$period != "Quaternary"),]
#hull_test <- hull_test[which(hull_test$period != "Archaean"),]

ggplot(hull_test, aes(V1, all)) + geom_polygon(aes(group = strat_name_id, fill = period), alpha = 0.5) +
  theme_bw() + 
  #scale_y_log10(limits = c(0.003, 0.3)) +
  #scale_x_continuous(limits = c(0.1, 0.5)) +
  geom_polygon(aes(group = strat_name_id), colour = "black", alpha = 0) +
  facet_wrap(~period) + geom_hline(yintercept = c(0.1, 0.15)) +
  geom_vline(xintercept = 0.35) +
  geom_point(aes(fill = period), colour = "black", pch = 21)

##### END #####
