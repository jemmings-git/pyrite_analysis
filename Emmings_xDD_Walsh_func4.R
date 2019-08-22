##### This R script was written by Joe Emmings (British Geological Survey) ####
##### This script is designed to manipulate and process GeoDeepDive extractions
##### interfaced with the Macrostrat database

##### PART 1 - Load and merge the GDD extractions and Macrostrat database ####

# install the following libraries (and dependencies)

library(dplyr) # plyr is not needed (note loading plyr after dplyr will prevent execution of PART 2)
library(tidyr)
library(readr)
library(jsonlite)
library(reshape2)
library(ggplot2)
library(devtools) # run once
devtools::install_github("moodymudskipper/safejoin") # needed for coalesce join # run once
library(safejoin)
devtools::install_github("adibender/pammtools") # provides geom_stepribbon for ggplot2
library(pammtools)
library(gridExtra)

# clear working environment if re-running this script from Part 1 in the same R session
rm(list=ls())

# set working directory
project_home <- paste(getwd(), "data", sep="/")
setwd(project_home)

#import GDD extractions

extracts <- read_csv("results.csv")

# remove unresolved strat_name_id hits

extracts <- extracts[!grepl(pattern = "\\~", extracts$strat_name_id),]

# import TWO Macrostrat database files - download and keep local copies
macrostrat_data <- function(filename, data_url){
  if (!file.exists(filename)) {
    download.file(data_url, filename, method="auto")
  }
  data <- fromJSON(filename)
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

# OR OPTION 2 - merged datasets without propagation of 'units'
#strat_units <- safe_right_join(units, strat, by = c("strat_name_id", "strat_name_id"), conflict = ~.y)
#strat_units <- strat_units[!duplicated(strat_units[c("strat_name_id", "t_age", "b_age")]),] # remove duplicate IDs but only for those with matching age ranges
# produces a more 'blocky' and apparently lower resolution output

extracts$strat_name_id  = as.character(extracts$strat_name_id)
data <- safe_right_join(extracts, strat_units, by = c("strat_name_id", "strat_name_id"), conflict = coalesce)

# optional export of the merged data file
# write.table(data, "data5.txt", sep="\t")

# optional remove extracts and strat from working environment
# rm(extracts, strat)

##### PART 2 - Prepare dataset for normalisation to Macrostrat units #####
##### - Addition of Phanerozoic time bins between top and base ages ####
##### subset dataframe into phanerozoic & precambrian (ultimately bins = 1 and 10, respectively)

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

data$lith  = as.character(data$lith)
data$environ  = as.character(data$environ)
data$econ  = as.character(data$econ)
data$measure = as.character(data$measure)
data$units_above = as.character(data$units_above)
data$units_below = as.character(data$units_below)
data$refs = as.character(data$refs)

# export the data
write.table(data, "data_comp2.txt", sep="\t")

# clear working environment
rm(list=ls())


##### PART 3 Simplify datafile, simplify lith and environ variables ####

#input datafile
data <- read.delim("data_comp2.txt")

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

# normalising to marine units or sedimentary units yields very similar results
# this is not surprising because pyrite precipitation requires a source of S
# this is dominantly SO4 in seawater, which is reduced (can be biotic or abiotic) in the
# water column or during diagenesis
# therefore non-marine sediments contain trace pyrite (or entirely lack pyrite)
# this builds confidence that the designation of 'marine' in the Macrostrat database is robust

###### PART 5 generic target script for running in one go (normalised to 'marine' units) ####
# this example attempts to replicate the Peters et al. stromatolites workflow

# select phrases of interest in grepl function below, this can be done quickly by appraising the 'phrases' dataframe

unique(phrases)

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
extracts <- read.delim("results.txt")
extracts <- subset(extracts, grepl("stromat", extracts$target_word) | grepl("Stromat", extracts$target_word) | grepl("thromb", extracts$target_word) | grepl("Thromb", extracts$target_word))
length(unique(extracts$docid))

# number of documents containing "stromatolitic"-terms matched to Macrostrat database:
stromatolites <- subset(data, grepl("stromat", data$target_word) | grepl("Stromat", data$target_word) | grepl("thromb", data$target_word) | grepl("Thromb", data$target_word))
stromatolites <- stromatolites[!is.na(stromatolites$strat_name_id),]
length(unique(as.factor(stromatolites$docid)))


###### END #####
