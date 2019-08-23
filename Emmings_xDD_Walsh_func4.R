##### This R script was written by Joe Emmings (British Geological Survey) ####
##### This script is designed to manipulate and process GeoDeepDive extractions
##### interfaced with the Macrostrat database

##### PART 1 - Load and merge the GDD extractions and Macrostrat database ####

source("data_preparation.R")

# install the following libraries (and dependencies)

library(dplyr) # plyr is not needed (note loading plyr after dplyr will prevent execution of PART 2)
library(tidyr)
library(readr)
library(jsonlite)
library(reshape2)
library(ggplot2)
library(devtools) # run once
#devtools::install_github("adibender/pammtools") # provides geom_stepribbon for ggplot2
library(pammtools)
library(gridExtra)

# clear working environment if re-running this script from Part 1 in the same R session
rm(list=ls())

data <- import_and_strat_align()
data <- subset_and_bin(data)

# Optionally save and/or load datafile
# export_data(data)

#input datafile
#data <- read.delim("data_comp2.txt")

data <- simplify_lith_env(data)

# optional export of the merged data file
# write.table(data, "data5.txt", sep="\t")

##### PART 4 - Analysis ####

# summarise target words (in a meaningful way)
supersaturation_keyphrases <- c("framboidal pyrite",
                                "pyrite framboid",
                                "pyrite framboids")
supersaturation <- subset(data, tolower(target_word) %in% supersaturation_keyphrases)

# We are picking up more data by checking for all cases (7653 vs 7335 this run)
advection_keyphrases <- c("pyrite nodule",
                          "pyrite nodules",
                          "pyritic nodule",
                          "pyritic nodules",
                          "nodular pyrite",
                          "nodules of pyrite",
                          "pyrite concretion",
                          "pyrite concretions",
                          "pyritic concretions",
                          "concretionary pyrite",
                          "concretions of pyrite")
advection <- subset(data, tolower(target_word) %in% advection_keyphrases)

# example binned data for 'supersaturation' and 'advection'

phanerozoic_increment <- 1 # in Ma
precambrian_increment <- 10 # in MA

# Phanerozoic bins
supersaturation_bins1 <- hist(supersaturation$value[supersaturation$value >= -0.5 & supersaturation$value <= 540.5], breaks = seq(-0.5, 540.5, by = phanerozoic_increment))
supersaturation_bins1 <- as.data.frame(cbind(supersaturation_bins1$counts, supersaturation_bins1$breaks))

# Precambrian bins
supersaturation_bins2 <- hist(supersaturation$value[supersaturation$value > 540.5 & supersaturation$value <= 2500.5], breaks = seq(540.5, 2500.5, by = precambrian_increment))
supersaturation_bins2 <- as.data.frame(cbind(supersaturation_bins2$counts, supersaturation_bins2$breaks))

# Bind
supersaturation_bins <- rbind(supersaturation_bins1,supersaturation_bins2)

# Phanerozoic bins
advection_bins1 <- hist(advection$value[advection$value >= -0.5 & advection$value <= 540.5], breaks = seq(-0.5, 540.5, by = phanerozoic_increment))
advection_bins1 <- as.data.frame(cbind(advection_bins1$counts, advection_bins1$breaks))

# Precambrian bins
advection_bins2 <- hist(advection$value[advection$value >= 540.5 & advection$value <= 2500.5], breaks = seq(540.5, 2500.5, by = precambrian_increment))
advection_bins2 <- as.data.frame(cbind(advection_bins2$counts, advection_bins2$breaks))

# Bind
advection_bins <- rbind(advection_bins1,advection_bins2)

##### PART 4b - normalisation now corrects for additional strat_name_IDs in 'strat' with missing environment of deposition info ####

# this introduces assumptions that the features of interest are primarily sedimentary and 'marine'
# it is possible to subset 'data' on env and lith but this will reduce n

# this also means the sedimentary and marine binned data is now dependent on target word(s)
# therefore this section must be run and adapted for each analysis

marine <- subset(data, tolower(target_word) %in% c(supersaturation_keyphrases, advection_keyphrases)
              | envCODE == "Marine")

sedimentary <- subset(data, tolower(target_word) %in% c(supersaturation_keyphrases, advection_keyphrases)
                      | lithCODE == "Sedimentary")

# binned data for 'marine' and 'sedimentary' units - this is a requirement for normalisation

# Phanerozoic bins
marine_bins1 <- hist(marine$value[marine$value >= -0.5 & marine$value <= 540.5], breaks = seq(-0.5, 540.5, by = phanerozoic_increment))
marine_bins1 <- as.data.frame(cbind(marine_bins1$counts, marine_bins1$breaks))

# Precambrian bins
marine_bins2 <- hist(marine$value[marine$value >= 540.5 & marine$value <= 2500.5], breaks = seq(540.5, 2500.5, by = precambrian_increment))
marine_bins2 <- as.data.frame(cbind(marine_bins2$counts, marine_bins2$breaks))

# Bind
marine_bins <- rbind(marine_bins1,marine_bins2)

# Phanerozoic bins
sedimentary_bins1 <- hist(sedimentary$value[sedimentary$value >= -0.5 & sedimentary$value <= 540.5], breaks = seq(-0.5, 540.5, by = phanerozoic_increment))
sedimentary_bins1 <- as.data.frame(cbind(sedimentary_bins1$counts, sedimentary_bins1$breaks))

# Precambrian bins
sedimentary_bins2 <- hist(sedimentary$value[sedimentary$value >= 540.5 & sedimentary$value <= 2500.5], breaks = seq(540.5, 2500.5, by = precambrian_increment))
sedimentary_bins2 <- as.data.frame(cbind(sedimentary_bins2$counts, sedimentary_bins2$breaks))

# Bind
sedimentary_bins <- rbind(sedimentary_bins1,sedimentary_bins2)

# Ratio of 'supersaturation' to 'advection'
pyrite_ratio <- as.data.frame(cbind(supersaturation_bins$V1/advection_bins$V1,advection_bins$V2))

# Proportion of 'supersaturation' (framboids) to 'advection' (nodules/concretions)
pyrite_prop <- as.data.frame(cbind(supersaturation_bins$V1/(supersaturation_bins$V1+advection_bins$V1),advection_bins$V2))

# normalised to sedimentary units
pyrite_sed_prop <- as.data.frame(cbind((supersaturation_bins$V1+advection_bins$V1)/sedimentary_bins$V1,advection_bins$V2))

# normalised to marine units
pyrite_marine_prop <- as.data.frame(cbind((supersaturation_bins$V1+advection_bins$V1)/marine_bins$V1,advection_bins$V2))

# fraction of supersaturation vs. advection plotted on top of the total normalised data
framboid_marine_prop <- as.data.frame(cbind((pyrite_prop$V1*pyrite_marine_prop$V1),supersaturation_bins$V2))

# replace NaN with 0 (required for correct plotting on ggplot)
framboid_marine_prop$V1[is.nan(framboid_marine_prop$V1)] <- 0

# or can do normalised to sedimentary units
framboid_sed_prop <- as.data.frame(cbind((pyrite_prop$V1*pyrite_sed_prop$V1),advection_bins$V2))
framboid_sed_prop$V1[is.nan(framboid_sed_prop$V1)] <- 0

# mass extinctions

ME <- c(445, 372, 252, 201, 65)
OAEs <- c(183, 120, 111, 93) # from Jenkyns 2010
PETM <- 55.8 # from Jenkyns 2010

# possible anoxic events

# step plot - marine
a <- ggplot(pyrite_marine_prop, aes(V2, V1)) + theme_bw() +
  geom_stepribbon(aes(ymin = 0, ymax = V1), fill = "gray85") +
  scale_x_reverse(limits = c(2500, 540)) +
  scale_y_continuous(limits = c(0,0.75)) + geom_step() +
  geom_stepribbon(data = framboid_marine_prop, aes(ymin = 0, ymax = V1), fill = "red") +
  geom_step(data = framboid_marine_prop)

b <- ggplot(pyrite_marine_prop, aes(V2, V1)) + theme_bw() +
  geom_stepribbon(aes(ymin = 0, ymax = V1), fill = "gray85") +
  scale_x_reverse(limits = c(541, -0.5)) +
  scale_y_continuous(limits = c(0,0.75)) + geom_step() +
  geom_stepribbon(data = framboid_marine_prop, aes(ymin = 0, ymax = V1), fill = "red") +
  geom_step(data = framboid_marine_prop) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_vline(xintercept = ME, colour = 'blue') +
  geom_vline(xintercept = OAEs, colour = 'red') +
  geom_vline(xintercept = PETM, colour = 'black')

grid.arrange(a,b, ncol = 2)

# step plot - sediments
a <- ggplot(pyrite_sed_prop, aes(V2, V1)) + theme_bw() +
  geom_stepribbon(aes(ymin = 0, ymax = V1), fill = "gray85") +
  scale_x_reverse(limits = c(2500, 540)) +
  scale_y_continuous(limits = c(0,0.75)) + geom_step() +
  geom_stepribbon(data = framboid_sed_prop, aes(ymin = 0, ymax = V1), fill = "red") +
  geom_step(data = framboid_sed_prop)

b <- ggplot(pyrite_sed_prop, aes(V2, V1)) + theme_bw() +
  geom_stepribbon(aes(ymin = 0, ymax = V1), fill = "gray85") +
  scale_x_reverse(limits = c(541, -0.5)) +
  scale_y_continuous(limits = c(0,0.75)) + geom_step() +
  geom_stepribbon(data = framboid_sed_prop, aes(ymin = 0, ymax = V1), fill = "red") +
  geom_step(data = framboid_sed_prop) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_vline(xintercept = ME, colour = 'blue') +
  geom_vline(xintercept = OAEs, colour = 'red') +
  geom_vline(xintercept = PETM, colour = 'black')

grid.arrange(a,b, ncol = 2)

