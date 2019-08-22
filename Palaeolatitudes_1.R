# This is the palaelatitude binning code preserved for future reference, won't run as-is

# subset palaeolatitude data (G-plates model) - see PART 6

plates <- subset(data, select = c(strat_name_id, target_word, unit_name, t_age, b_age, strat_name, t_plat, t_plng, b_plat, b_plng, lith, environ))
plates <- plates[!is.na(plates$t_plat),]
plates$lith  = as.character(plates$lith)
plates$environ  = as.character(plates$environ)
write.table(plates, "plates.txt", sep="\t")

#### PART 6 - G-plates palaeolatitude - method 1 - plotted as polygons ####

# re-import selected file

plates <- read.delim("plates.txt")

plates$envCODE <- gsub("TRUE", "Non-Marine", grepl("non-marine", plates$environ))
plates$envCODE <- gsub("FALSE", "Marine", plates$envCODE)
plates$envCODE[is.na(plates$unit_name)] <- NA
plates$lithCODE <- gsub("TRUE", "Sedimentary", grepl("sedimentary", plates$lith))
plates$lithCODE <- gsub("FALSE", "Non-Sed", plates$lithCODE)
plates$lithCODE[is.na(plates$unit_name)] <- NA

# select targets - here all target words are selected in one go, but ultimately it is probably
# appopriate to plot polygons for concretions, framboids, and marine units, etc. seperately

plates_target <- subset(plates, target_word == "framboidal pyrite"
                        | target_word == "Framboidal pyrite"
                        | target_word == "Framboidal Pyrite"
                        | target_word == "pyrite framboids"
                        | target_word == "pyrite framboid"
                        | target_word == "Pyrite framboids"
                        | target_word == "Pyrite framboid"
                        | target_word == "Pyrite Framboids"
                        | target_word == "Pyrite Framboid"
                        | target_word == "PYRITE FRAMBOID"
                        | target_word == "PYRITE FRAMBOIDS"
                        | target_word == "pyrite nodules"
                        | target_word == "pyrite concretions"
                        | target_word == "Pyrite nodules"
                        | target_word == "Pyrite nodule"
                        | target_word == "nodular pyrite"
                        | target_word == "concretionary pyrite"
                        | target_word == "pyritic nodules"
                        | target_word == "pyritic concretions"
                        | target_word == "concretions of pyrite"
                        | target_word == "pyritic nodule"
                        | target_word == "nodules of pyrite"
                        | target_word == "Pyrite nodule"
                        | target_word == "Pyrite concretions"
                        | target_word == "Nodular pyrite"
                        | target_word == "Pyritic nodules"
                        | target_word == "pyrite concretion"
                        | target_word == "Pyrite Nodules"
                        | target_word == "Pyrite concretion"
                        | target_word == "PYRITE NODULES"
                        | target_word == "PYRITE NODULE"
                        | envCODE == "Marine")

plates_target <- plates_target[unique(plates_target$strat_name_id),]

colnames(plates_target)

# melt age
plates_target2 <- melt(plates_target, id.vars = c("strat_name_id", "target_word", 
                                                  "unit_name", "strat_name", "t_plat", "t_plng",
                                                  "b_plat", "b_plng", "lith", "environ",
                                                  "envCODE","lithCODE"),
                       variable.name = "age_var", value.name = "age")


plates_target3 <- subset(plates_target2, age_var != "b_age")

# melt plat for t_age

plates_target3 <- melt(plates_target3, id.vars = c("strat_name_id", "target_word", 
                                                   "unit_name", "strat_name", "t_plng",
                                                   "b_plng", "lith", "environ",
                                                   "envCODE","lithCODE", "age_var", "age"),
                       variable.name = "plat_var", value.name = "plat")

plates_target3 <- plates_target3[!is.na(plates_target3$plat),]

plates_target4 <- subset(plates_target2, age_var == "b_age")

# melt plat for b_age

plates_target4 <- melt(plates_target4, id.vars = c("strat_name_id", "target_word", 
                                                   "unit_name", "strat_name", "t_plng",
                                                   "b_plng", "lith", "environ",
                                                   "envCODE","lithCODE", "age_var", "age"),
                       variable.name = "plat_var", value.name = "plat")

# combine 

plates_target4 <- plates_target4[!is.na(plates_target4$plat),]
plates_target4 <- plates_target4[order(as.character(plates_target4$plat_var)),]
plates_target5 <- rbind(plates_target3, plates_target4)
plates_target6 <- plates_target5[order(plates_target5$unit_name),]
plates_target6 <- unique(plates_target6[c("strat_name_id", "age_var", "plat_var")])
plates_target5$Row.names <- as.factor(row.names(plates_target5))
plates_target6$Row.names <- as.factor(row.names(plates_target6))
plates_target <- left_join(plates_target6, plates_target5, by = "Row.names")

plates_target <- plates_target[order(plates_target$unit_name),]
plates_target <- plates_target[,5:18]

# output

ggplot(plates_target) + geom_polygon(aes(age, plat, group = strat_name_id.y), alpha =0.5) +
  scale_x_reverse(limits = c(2500,0)) + theme_bw()

# plotted with framboids & concretions seperately

plates_framboids <- subset(plates_target, target_word == "framboidal pyrite"
                           | target_word == "Framboidal pyrite"
                           | target_word == "Framboidal Pyrite"
                           | target_word == "pyrite framboids"
                           | target_word == "pyrite framboid"
                           | target_word == "Pyrite framboids"
                           | target_word == "Pyrite framboid"
                           | target_word == "Pyrite Framboids"
                           | target_word == "Pyrite Framboid"
                           | target_word == "PYRITE FRAMBOID"
                           | target_word == "PYRITE FRAMBOIDS")

plates_nods <- subset(plates_target, target_word == "pyrite nodules"
                      | target_word == "pyrite concretions"
                      | target_word == "Pyrite nodules"
                      | target_word == "Pyrite nodule"
                      | target_word == "nodular pyrite"
                      | target_word == "concretionary pyrite"
                      | target_word == "pyritic nodules"
                      | target_word == "pyritic concretions"
                      | target_word == "concretions of pyrite"
                      | target_word == "pyritic nodule"
                      | target_word == "nodules of pyrite"
                      | target_word == "Pyrite nodule"
                      | target_word == "Pyrite concretions"
                      | target_word == "Nodular pyrite"
                      | target_word == "Pyritic nodules"
                      | target_word == "pyrite concretion"
                      | target_word == "Pyrite Nodules"
                      | target_word == "Pyrite concretion"
                      | target_word == "PYRITE NODULES")

# combined plot - but it doesn't work very well, because few strat units with target mentions contain plat info

ggplot(plates_target) + geom_polygon(aes(age, plat, group = strat_name_id.y), alpha =0.5) +
  scale_x_reverse(limits = c(2500,0)) + theme_bw() +
  geom_polygon(data = plates_framboids, aes(age, plat, group = strat_name_id.y), alpha =0.5, fill = 'red') +
  geom_polygon(data = plates_nods, aes(age, plat, group = strat_name_id.y), alpha =0.5, fill = 'blue')



##### PART 6 - G-plates palaeolatitude - method 2 - plotted as rasters ####

# age split

plates1 <- data.frame(plates) # duplicate for testing

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

plates1_precambrian <- subset(plates1, b_age > cutoff)
plates1_precambrian <- add_all_bins(plates1_precambrian, precambrian_increment)

plates1_phanerozoic <- subset(plates1, b_age < cutoff)
plates1_phanerozoic <- add_all_bins(plates1_phanerozoic, phanerozoic_increment)

# Bind the Phanerozoic and Precambrian datasets

plates1 <- bind_rows(plates1_phanerozoic, plates1_precambrian)

plates1 <- melt(plates1, id.vars = c("strat_name_id", "target_word","unit_name", "strat_name", 
                                     "t_plat", "t_plng", "b_plat", "b_plng", 
                                     "lith", "environ", "envCODE", "lithCODE"), na.rm=TRUE,
                variable.name = "age_var", value.name = "age")

plates1$age <- as.numeric(plates1$age)

# plat split

cutoff <- 541 # Ma 
phanerozoic_increment <- 1 # in degrees
precambrian_increment <- 1 # in degrees

plates1$b_plat <- plates1$b_plat + 100 # transform to avoid negative values
plates1$t_plat <- plates1$t_plat + 100 # transform to avoid negative values


add_all_bins <- function(df, increment) {
  ncols <- ceiling((max(df$b_plat - df$t_plat))/increment)
  
  df <- df %>% mutate(bin1 = ifelse(b_plat-increment > t_plat+increment, t_plat+increment, NA))
  for(n in 2:(ncols-1)){
    next_bin <- paste('bin', n, sep='')
    bin <- paste0('bin', (n-1), sep='')
    # See https://stackoverflow.com/a/49311813/4319767 for why this syntax 
    df <- df %>% mutate(!!next_bin := ifelse(UQ(rlang::sym(bin)) < b_plat-increment, UQ(rlang::sym(bin))+increment, NA))
  }  
  return(df)
}

plates1_precambrian <- subset(plates1, age > cutoff)
plates1_precambrian <- add_all_bins(plates1_precambrian, precambrian_increment)

plates1_phanerozoic <- subset(plates1, age < cutoff)
plates1_phanerozoic <- add_all_bins(plates1_phanerozoic, phanerozoic_increment)

# Bind the Phanerozoic and Precambrian datasets

plates1 <- bind_rows(plates1_phanerozoic, plates1_precambrian)

plates1 <- melt(plates1, id.vars = c("strat_name_id", "target_word","unit_name", "strat_name", 
                                     "t_plng", "b_plng", "age_var", "age",
                                     "lith", "environ", "envCODE", "lithCODE"), na.rm=TRUE,
                variable.name = "plat_var", value.name = "plat")

plates1$plat <- as.numeric(plates1$plat)
plates1$plat <- plates1$plat - 100 # reverse transform

Top <- 0 # top age for plot (in Ma)
Bottom <- 2500 # bottom age for plot (in Ma)

# density plots

a <- ggplot(plates1, aes(age, plat)) + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  theme_bw() + scale_fill_viridis_c(option="inferno") +
  scale_x_reverse(limits = c(Bottom, 540)) +
  scale_y_continuous(limits = c(-90,90)) + theme(legend.position = "none")

b <- ggplot(plates1, aes(age, plat)) + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  theme_bw() + scale_fill_viridis_c(option="inferno") +
  scale_x_reverse(limits = c(540, Top)) +
  scale_y_continuous(limits = c(-90,90)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

grid.arrange(a,b, ncol = 2)

# alternative - 2d binned plots

a <- ggplot(plates1, aes(age, plat)) + geom_bin2d(binwidth = 10) +
  theme_bw() + scale_fill_viridis_c(option="inferno", trans = "log10", limits = c(1, 1000)) +
  scale_x_reverse(limits = c(Bottom, 540)) +
  scale_y_continuous(limits = c(-90,90)) + theme(legend.position = "none")

b <- ggplot(plates1, aes(age, plat)) + geom_bin2d(binwidth = 1) +
  theme_bw() + scale_fill_viridis_c(option="inferno", trans = "log10", limits = c(1, 1000)) +
  scale_x_reverse(limits = c(540, Top)) +
  scale_y_continuous(limits = c(-90,90)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

grid.arrange(a,b, ncol = 2)

# can now subset as per. method 1 if desired