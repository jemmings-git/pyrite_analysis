# normalising to marine units or sedimentary units yields very similar results
# this is not surprising because pyrite precipitation requires a source of S
# this is dominantly SO4 in seawater, which is reduced (can be biotic or abiotic) in the
# water column or during diagenesis
# therefore non-marine sediments contain trace pyrite (or entirely lack pyrite)
# this builds confidence that the designation of 'marine' in the Macrostrat database is robust
source('data_preparation.R')

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

###### PART 5 generic target script for running in one go (normalised to 'marine' units) ####
# this example attempts to replicate the Peters et al. stromatolites workflow

# select phrases of interest in grepl function below, this can be done quickly by appraising the 'phrases' dataframe

# unique(phrases)

# quick replication of stromatolites demo below - not perfect
# includes subset for marine settings (Phanerozoic) - otherwise the output is skewed by
# freshwater data such as the Palaeocene-Eocene Green River Fm.
# this method is consistent with Peters et al. 2017

target <- grepl("stromat", data$target_word) | grepl("Stromat", data$target_word) |
  grepl("thromb", data$target_word) | grepl("Thromb", data$target_word)

target <- subset(data, target)

target1 <- subset(target, value < 541)
target1 <- subset(target1, envCODE == "Marine")
target2 <- subset(target, value > 541)

target <- rbind(target1, target2)

# remove mentions of 'non' i.e., non-stromatolitic

false <- !grepl("non", target$target_word) & !grepl("Non", target$target_word)

target <- subset(target, false)

marine <- grepl("stromat", data$target_word) | grepl("Stromat", data$target_word) |
  grepl("thromb", data$target_word) | grepl("Thromb", data$target_word) |
  data$envCODE == "Marine"

marine <- subset(data, marine)

marine <- subset(marine, false)

sedimentary <- grepl("stromat", data$target_word) | grepl("Stromat", data$target_word) |
  grepl("thromb", data$target_word) | grepl("Thromb", data$target_word) |
  data$envCODE == "Sedimentary"

sedimentary <- subset(data, sedimentary)

sedimentary <- subset(sedimentary, false)

Top <- -0 # top age for plot (in Ma), set at -0.5 in order to centre bins
Bottom <- 2500 # bottom age for plot (in Ma)
phanerozoic_increment <- 1
precambrian_increment <- 10

# bin calculation

# Phanerozoic bins
marine_bins1 <- hist(marine$value[marine$value >= -0.5 & marine$value <= 540.5], breaks = seq(-0.5, 540.5, by = phanerozoic_increment))
marine_bins1 <- as.data.frame(cbind(marine_bins1$counts, marine_bins1$breaks))

# Precambrian bins
marine_bins2 <- hist(marine$value[marine$value >= 540.5 & marine$value <= 2500.5], breaks = seq(540.5, 2500.5, by = precambrian_increment))
marine_bins2 <- as.data.frame(cbind(marine_bins2$counts, marine_bins2$breaks))

# Bind
marine_bins <- rbind(marine_bins1,marine_bins2)

# Ignore V2 = 0 (anomaly)
marine_bins_se <- marine_bins
marine_bins_se[1,1] <- NA

# Standard error
marine_bins$se <- (median((marine_bins_se$V1-median(marine_bins_se$V1, na.rm = TRUE))^2, na.rm = TRUE))^0.5 # combined se calc

# Phanerozoic bins
sedimentary_bins1 <- hist(sedimentary$value[sedimentary$value >= -0.5 & sedimentary$value <= 540.5], breaks = seq(-0.5, 540.5, by = phanerozoic_increment))
sedimentary_bins1 <- as.data.frame(cbind(sedimentary_bins1$counts, sedimentary_bins1$breaks))

# Precambrian bins
sedimentary_bins2 <- hist(sedimentary$value[sedimentary$value >= 540.5 & sedimentary$value <= 2500.5], breaks = seq(540.5, 2500.5, by = precambrian_increment))
sedimentary_bins2 <- as.data.frame(cbind(sedimentary_bins2$counts, sedimentary_bins2$breaks))

# Bind
sedimentary_bins <- rbind(sedimentary_bins1,sedimentary_bins2)

# Ignore V2 = 0 (anomaly)
sedimentary_bins_se <- sedimentary_bins
sedimentary_bins_se[1,1] <- NA

# Standard error - provisional - this needs some more thought
sedimentary_bins$se <- (median((sedimentary_bins_se$V1-median(sedimentary_bins_se$V1, na.rm = TRUE))^2, na.rm = TRUE))^0.5 # combined se calc

target <- subset(data, target)

# may need to consider seperate se for precambrian & phanerozoic?

# Phanerozoic bins
target_bins1 <- hist(target$value[target$value >= -0.5 & target$value <= 540.5], breaks = seq(-0.5, 540.5, by = phanerozoic_increment))
target_bins1 <- as.data.frame(cbind(target_bins1$counts, target_bins1$breaks))

# Precambrian bins
target_bins2 <- hist(target$value[target$value >= 540.5 & target$value <= 2500.5], breaks = seq(540.5, 2500.5, by = precambrian_increment))
target_bins2 <- as.data.frame(cbind(target_bins2$counts, target_bins2$breaks))

# Bind
target_bins <- rbind(target_bins1,target_bins2)

# Ignore V2 = 0 (anomaly)
target_bins_se <- target_bins
target_bins_se[1,1] <- NA

# Standard error
target_bins$se <- (median((target_bins_se$V1-median(target_bins_se$V1, na.rm = TRUE))^2, na.rm = TRUE))^0.5 # combined se calc

# propagate standard errors in fractional quadrature

standard_error_quad <- ((((100/target_bins$V1)*target_bins$se)^2)+(((100/marine_bins$V1)*marine_bins$se)^2))^0.5


target_marine_prop <- as.data.frame(cbind((target_bins$V1)/marine_bins$V1,
                                          target_bins$V2,
                                          ((target_bins$V1)/marine_bins$V1)*(standard_error_quad/100)))

target_marine_prop$topse <- target_marine_prop$V1+target_marine_prop$V3
target_marine_prop$basese <- target_marine_prop$V1-target_marine_prop$V3

target_marine_prop$topse[target_marine_prop$topse > 1] <- 1
target_marine_prop$basese[target_marine_prop$basese < 0] <- 0

# Output (Peters et al. style step)

a <- ggplot(target_marine_prop, aes(V2, V1)) + theme_bw() +
  geom_stepribbon(aes(ymin = basese, ymax = topse), fill = "gray85") +
  scale_x_reverse(limits = c(Bottom, 540)) +
  scale_y_continuous(limits = c(0,1)) +
  geom_step(aes(V2, topse), colour = 'gray') +
  geom_step(aes(V2, basese), colour = 'gray') + geom_step()

b <- ggplot(target_marine_prop, aes(V2, V1)) + theme_bw() +
  geom_stepribbon(aes(ymin = basese, ymax = topse), fill = "gray85") +
  scale_x_reverse(limits = c(541, Top)) +
  scale_y_continuous(limits = c(0,1)) +
  geom_step(aes(V2, topse), colour = 'gray') +
  geom_step(aes(V2, basese), colour = 'gray') + geom_step() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

grid.arrange(a,b, ncol = 2)


# note it might be possible to map the international stratigraphic chart onto ggplot using
# https://macrostrat.org/api/defs/intervals?all&format=json

# optional - calculate number of stromatolite entries matched to Macrostrat (as a comparison to Peters et al. 2017)

#overall number of documents containing "stromatolitic"-terms with strat_name_id
extracts <- read.csv("data/results.csv")
extracts <- subset(extracts, grepl("stromat", extracts$target_word) | grepl("Stromat", extracts$target_word) | grepl("thromb", extracts$target_word) | grepl("Thromb", extracts$target_word))
length(unique(extracts$docid))

# number of documents containing "stromatolitic"-terms matched to Macrostrat database:
stromatolites <- subset(data, grepl("stromat", data$target_word) | grepl("Stromat", data$target_word) | grepl("thromb", data$target_word) | grepl("Thromb", data$target_word))
stromatolites <- stromatolites[!is.na(stromatolites$strat_name_id),]
length(unique(as.factor(stromatolites$docid)))


###### END #####
