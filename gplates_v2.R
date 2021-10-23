
## GPlates analysis ##

library(chronosphere)
library(rgdal)
library(gridExtra)
library(plyr)
library(dplyr)
library(ggplot2)
library(ncdf4)
library(raster)
library(reshape2)
library(scales)
library(ggnewscale)
library(scatterpie)
library(geoR)
library(ggrepel)
library(grid)
library(tidyverse)
library(deeptime)
library(ggfittext)
library(gstat)
library(geosphere)
library(RColorBrewer)

#out <- datasets(dat = NULL)
#oneDat <- datasets("paleomap")

###### Palaeo-bathymetry (part 1) #####

dem <- fetch(dat="paleomap", var="dem", res=1)

dem.ages <- row.names(as.data.frame(dem@index)) # list available dem ages

#out <- dem[c(1,8,seq(13, length(dem), by = 2))]

dem1 <- chronosphere::as.list(dem)
#dem <- dem1

# optional downsample

factor <- 1

dem1 <- lapply(dem1, function(d) raster::aggregate(d, fact=factor)) #  downsample

# an interesting idea - could plot standard deviation of bathymetry
# sensu Scotese restriction metric - see PDF download
# try moving standard deviation across matrix, select a window

dem1 <- lapply(dem1, function(d) raster::as.matrix(d))
dem1 <- lapply(dem1, function(d) as.data.frame(d))
names(dem1) <- dem.ages

dem2 <- ldply(dem1, data.frame)
names(dem2)[1] <- "Interval"

colnames(dem2) <- c("Interval", seq(-180,180,by=1*factor))
dem2$lat <- (seq(-90,90,by =1*factor))*(-1)

dem3 <- reshape2::melt(dem2, id.vars = c("Interval","lat"))
names(dem3)[3] <- "long"

dem3$Interval <- as.character(dem3$Interval) 
dem3$long <- as.character(dem3$long) 
dem3$long <- as.numeric(dem3$long) 

dem3$zone <- as.numeric(as.character(dem3$Interval))
dem3$zone2 <- cut(dem3$zone,c(540,510,480,455,420,385,380,365,330,300,280,250,100,80))
dem3 <- dem3[!is.na(dem3$zone2),]

dem4 <- dem3 %>% 
  group_by(zone2,lat,long) %>% 
  summarise(value = min(value)) # flatten using mean or min topo

dem4 <- subset(dem4, zone2 != "(100,250]")

targets <- c(90,265,290,315,350,375,385,400,430,470,495,525) # select key time slices

#land2 <- subset(dem3, value > 0) # option 2
#land2 <- land2[land2$Interval %in% targets, ]

# projected coastlines (optional)

coasts <- reconstruct("coastlines", age= c(0,targets)) #c(seq(0,750, by = 50),

coasts1 <- lapply(coasts, function(d) fortify(d))

coasts1 <- ldply(coasts1, data.frame)
names(coasts1)[1] <- "Interval"
coasts1$Interval <- as.numeric(coasts1$Interval)

coasts1$zone2 <- cut(coasts1$Interval,c(540,510,480,455,420,385,380,365,330,300,280,250,100,80))
#coasts1 <- coasts1[!is.na(coasts1$zone2),]

# view bathymetry min or mean across time intervals of interest

ggplot(dem4) +
  geom_tile(aes(x = long, y = lat, fill = value)) +
  facet_wrap(~zone2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  #geom_polygon(data = coasts1, aes(x = long, y = lat, group = group), colour = "black", fill = NA) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value)))
  #scale_fill_distiller(limits = c(-2000,0), palette = "RdBu", oob = squish)

## import Mukherjee and Large pyrites with clusters 1-5

setwd("D:/Pyrite_backup/pyrite_analysis-master/pyrite_analysis-master")

pyrites <- read.csv("pyrite_stats_with_long_lats.csv")

# subset < 750 Ma

pyrites <- subset(pyrites, Age.Ma <= 750)
pyrites <- pyrites[!is.na(pyrites$Lat),]

pyrites <- pyrites[,-1]

modern.coast <- subset(coasts1, Interval == 0)

# plot pyrites dataset on present day map

ggplot(pyrites) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_polygon(data = modern.coast, aes(x = long, y = lat, group = group), colour = "black", fill = "white") +
  geom_point(aes(Long, Lat), colour = "red")

# test scatterpie

pyrites2 <- pyrites[,-4]

pyrites2 <- reshape2::dcast(pyrites2, Location+Age.Ma+Lat+Long+Phrase+Source+Notes~clusters, value.var = "perc")

pyrites2[c("Type.1","Type.2","Type.3","Type.4","Type.5")][is.na(pyrites2[c("Type.1","Type.2","Type.3","Type.4","Type.5")])] <- 0

ggplot(pyrites2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_polygon(data = modern.coast, aes(x = long, y = lat, group = group), colour = "black", fill = "white") +
  #geom_point(aes(Long, Lat), colour = "red") +
  geom_scatterpie(aes(Long, Lat, r= 3), data = pyrites2,
                  cols=c("Type.1","Type.2","Type.3","Type.4","Type.5"), color= NA)

pyrites2$Long <- round(pyrites2$Long, digits = 1)
pyrites2$Lat <- round(pyrites2$Lat, digits = 1)

pyrites3 <- pyrites2[1:100,]
pyrites4 <- pyrites2[101:150,]
pyrites5 <- pyrites2[151:nrow(pyrites2),]

coords1 <- reconstruct(pyrites3[, c("Long", "Lat"), drop=FALSE], age= round(pyrites3$Age.Ma/5)*5, enumerate=FALSE, verbose=FALSE)
coords2 <- reconstruct(pyrites4[, c("Long", "Lat"), drop=FALSE], age= round(pyrites4$Age.Ma/5)*5, enumerate=FALSE, verbose=FALSE)
coords3 <- reconstruct(pyrites5[, c("Long", "Lat"), drop=FALSE], age= round(pyrites5$Age.Ma/5)*5, enumerate=FALSE, verbose=FALSE)

pyrites3 <- cbind(pyrites3, coords1)
pyrites4 <- cbind(pyrites4, coords2)
pyrites5 <- cbind(pyrites5, coords3)

pyrites2 <- rbind(pyrites3,pyrites4,pyrites5)

names(pyrites2)[13:14] <- c("Pal.Long","Pal.Lat")

pyrites2$zone2 <- cut(pyrites2$Age.Ma,c(540,510,480,455,420,385,380,365,330,300,280,250,100,80))
pyrites2 <- pyrites2[!is.na(pyrites2$zone2),]
pyrites2 <- pyrites2[!is.na(pyrites2$Pal.Long),]

pyrites2 <- subset(pyrites2, zone2 != "(100,250]")

pyrites2$group <- 1:nrow(pyrites2)
pyrites2$group <- as.character(pyrites2$group)

# get repel coordinates

unique(pyrites2$zone2)

# optional group pyrites by location - as presented in manuscript

# to plot by location AND age bin use 'pyrites.all' and Pal.Long.jit, Pal.Lat.jit

pyrites2 <- pyrites2 %>%
  group_by(Location, zone2) %>%
  summarise(Type.1 = mean(Type.1),Type.2 = mean(Type.2),
            Type.3 = mean(Type.3),Type.4 = mean(Type.4),
            Type.5 = mean(Type.5), Pal.Long = mean(Pal.Long),
            Pal.Lat = mean(Pal.Lat))

pyrites2$group <- 1:nrow(pyrites2)
pyrites2$group <- as.character(pyrites2$group)

# use names(pyrites2.sub1)[17:18] (and so on) for pyrites.all
# use pyrites2.sub1 <- pyrites2.sub1[-19] (and so on) for pyrites.all


##### (80,100] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(80,100]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(0.5, "lines")) #+
  #geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub1 <- cbind(pyrites2.sub,out)
names(pyrites2.sub1)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub1$Pal.Long.jit <- (pyrites2.sub1$Pal.Long.jit*360)-180
pyrites2.sub1$Pal.Lat.jit <- (pyrites2.sub1$Pal.Lat.jit*180)-90
pyrites2.sub1 <- pyrites2.sub1[-13]

##### (100,250] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(100,250]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub2 <- cbind(pyrites2.sub,out)
names(pyrites2.sub2)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub2$Pal.Long.jit <- (pyrites2.sub2$Pal.Long.jit*360)-180
pyrites2.sub2$Pal.Lat.jit <- (pyrites2.sub2$Pal.Lat.jit*180)-90
pyrites2.sub2 <- pyrites2.sub2[-13]

##### (250,280] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(250,280]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub3 <- cbind(pyrites2.sub,out)
names(pyrites2.sub3)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub3$Pal.Long.jit <- (pyrites2.sub3$Pal.Long.jit*360)-180
pyrites2.sub3$Pal.Lat.jit <- (pyrites2.sub3$Pal.Lat.jit*180)-90
pyrites2.sub3 <- pyrites2.sub3[-13]

##### (280,300] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(280,300]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub4 <- cbind(pyrites2.sub,out)
names(pyrites2.sub4)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub4$Pal.Long.jit <- (pyrites2.sub4$Pal.Long.jit*360)-180
pyrites2.sub4$Pal.Lat.jit <- (pyrites2.sub4$Pal.Lat.jit*180)-90
pyrites2.sub4 <- pyrites2.sub4[-13]

##### (300,330] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(300,330]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub5 <- cbind(pyrites2.sub,out)
names(pyrites2.sub5)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub5$Pal.Long.jit <- (pyrites2.sub5$Pal.Long.jit*360)-180
pyrites2.sub5$Pal.Lat.jit <- (pyrites2.sub5$Pal.Lat.jit*180)-90
pyrites2.sub5 <- pyrites2.sub5[-13]

##### (330,365] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(330,365]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)

out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub6 <- cbind(pyrites2.sub,out)
names(pyrites2.sub6)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub6$Pal.Long.jit <- (pyrites2.sub6$Pal.Long.jit*360)-180
pyrites2.sub6$Pal.Lat.jit <- (pyrites2.sub6$Pal.Lat.jit*180)-90
pyrites2.sub6 <- pyrites2.sub6[-13]

##### (365,380] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(365,380]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub7 <- cbind(pyrites2.sub,out)
names(pyrites2.sub7)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub7$Pal.Long.jit <- (pyrites2.sub7$Pal.Long.jit*360)-180
pyrites2.sub7$Pal.Lat.jit <- (pyrites2.sub7$Pal.Lat.jit*180)-90
pyrites2.sub7 <- pyrites2.sub7[-13]

##### (380,385] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(380,385]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub8 <- cbind(pyrites2.sub,out)
names(pyrites2.sub8)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub8$Pal.Long.jit <- (pyrites2.sub8$Pal.Long.jit*360)-180
pyrites2.sub8$Pal.Lat.jit <- (pyrites2.sub8$Pal.Lat.jit*180)-90
pyrites2.sub8 <- pyrites2.sub8[-13]

##### (385,420] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(385,420]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub9 <- cbind(pyrites2.sub,out)
names(pyrites2.sub9)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub9$Pal.Long.jit <- (pyrites2.sub9$Pal.Long.jit*360)-180
pyrites2.sub9$Pal.Lat.jit <- (pyrites2.sub9$Pal.Lat.jit*180)-90
pyrites2.sub9 <- pyrites2.sub9[-13]

##### (420,455] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(420,455]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub10 <- cbind(pyrites2.sub,out)
names(pyrites2.sub10)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub10$Pal.Long.jit <- (pyrites2.sub10$Pal.Long.jit*360)-180
pyrites2.sub10$Pal.Lat.jit <- (pyrites2.sub10$Pal.Lat.jit*180)-90
pyrites2.sub10 <- pyrites2.sub10[-13]

##### (455,480] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(455,480]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub11 <- cbind(pyrites2.sub,out)
names(pyrites2.sub11)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub11$Pal.Long.jit <- (pyrites2.sub11$Pal.Long.jit*360)-180
pyrites2.sub11$Pal.Lat.jit <- (pyrites2.sub11$Pal.Lat.jit*180)-90
pyrites2.sub11 <- pyrites2.sub11[-13]

##### (480,510] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(480,510]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(1.5, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub12 <- cbind(pyrites2.sub,out)
names(pyrites2.sub12)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub12$Pal.Long.jit <- (pyrites2.sub12$Pal.Long.jit*360)-180
pyrites2.sub12$Pal.Lat.jit <- (pyrites2.sub12$Pal.Lat.jit*180)-90
pyrites2.sub12 <- pyrites2.sub12[-13]

##### (510,540] #####

pyrites2.sub <- subset(pyrites2, zone2 == "(510,540]")

p1 <- ggplot(dem4) +
  #geom_tile(aes(x = long, y = lat, fill = value)) +
  #facet_wrap(~zone2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "colorbar", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites2, aes(Pal.Long, Pal.Lat)) +
  #geom_point(data = out, aes((V1*360)-180, (V2*180)-90), colour = "red") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01)) +
  geom_text_repel(data = pyrites2.sub, aes(x = Pal.Long, y = Pal.Lat, label = group), max.overlaps = 1000, force_pull = 1, force = 5, size = 10,
                  min.segment.length = 0.01, box.padding = unit(0.8, "lines")) #+
#geom_text_repel(data = test, aes((V1*360)-180, (V2*180)-90, label = row.names(test)), colour = "red")

p1

grid.force()

kids <- grid.get("textrepeltree", grep = TRUE)
kids <- kids$children

for (i in 1:142) {
  x2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["x1"]]
  y2 <-kids[[paste("segmentrepelgrob",i, sep = "")]][[".ORIGINAL"]][["y1"]]
  assign(paste("segment_x2",i, sep = ""), x2)
  assign(paste("segment_y2",i, sep = ""), y2)
  rm(surface)
  output.x2 <- do.call(list, mget(ls(pattern="segment_x2[0-9]")))
  output.y2 <- do.call(list, mget(ls(pattern="segment_y2[0-9]")))
}

df_x <- as.data.frame(do.call("rbind", output.x2))
df_y <- as.data.frame(do.call("rbind", output.y2))

out <- cbind(df_x,df_y)
names(out)[2] <- "V2"

out$group <- row.names(out)
out$group <- gsub("segment_x2","",out$group)
out$group <- as.numeric(out$group)
out <- out[order(out$group),]

pyrites2.sub13 <- cbind(pyrites2.sub,out)
names(pyrites2.sub13)[11:12] <- c("Pal.Long.jit","Pal.Lat.jit")
pyrites2.sub13$Pal.Long.jit <- (pyrites2.sub13$Pal.Long.jit*360)-180
pyrites2.sub13$Pal.Lat.jit <- (pyrites2.sub13$Pal.Lat.jit*180)-90
pyrites2.sub13 <- pyrites2.sub13[-13]

pyrites.all <- rbind(pyrites2.sub1, pyrites2.sub2, pyrites2.sub3, pyrites2.sub4, pyrites2.sub5,
                     pyrites2.sub6, pyrites2.sub7, pyrites2.sub8, pyrites2.sub9, pyrites2.sub10,
                     pyrites2.sub11, pyrites2.sub12, pyrites2.sub13)

##### plot together ####

# import xdd combined outputs #

pyrites.xdd <- read.csv("pyrite_combined_manual.csv")

pyrites.xdd <- pyrites.xdd[!is.na(pyrites.xdd$clat),]

pyrites.xdd <- subset(pyrites.xdd, b_age < 600)

pyrites.xdd$av_age <- round((pyrites.xdd$t_age+(pyrites.xdd$b_age-pyrites.xdd$t_age)/2)/5)*5

pyrites.xdd0 <- pyrites.xdd[1:50,]
pyrites.xdd1 <- pyrites.xdd[51:100,]
pyrites.xdd2 <- pyrites.xdd[101:150,]
pyrites.xdd3 <- pyrites.xdd[151:200,]
pyrites.xdd4 <- pyrites.xdd[201:250,]
pyrites.xdd5 <- pyrites.xdd[251:300,]
pyrites.xdd6 <- pyrites.xdd[301:350,]
pyrites.xdd7 <- pyrites.xdd[351:400,]
pyrites.xdd8 <- pyrites.xdd[401:450,]
pyrites.xdd9 <- pyrites.xdd[451:500,]
pyrites.xdd10 <- pyrites.xdd[501:550,]
pyrites.xdd11 <- pyrites.xdd[551:nrow(pyrites.xdd),]

coords0 <- reconstruct(pyrites.xdd0[, c("clng", "clat"), drop=FALSE], age= c(pyrites.xdd0$av_age), enumerate=FALSE, verbose=FALSE)
pyrites.xdd0 <- cbind(pyrites.xdd0,coords0)
names(pyrites.xdd0)[40:41] <- c("clng.pal","clat.pal")

coords1 <- reconstruct(pyrites.xdd1[, c("clng", "clat"), drop=FALSE], age= c(pyrites.xdd1$av_age), enumerate=FALSE, verbose=FALSE)
pyrites.xdd1 <- cbind(pyrites.xdd1,coords1)
names(pyrites.xdd1)[40:41] <- c("clng.pal","clat.pal")

coords2 <- reconstruct(pyrites.xdd2[, c("clng", "clat"), drop=FALSE], age= pyrites.xdd2$av_age, enumerate=FALSE, verbose=FALSE)
pyrites.xdd2 <- cbind(pyrites.xdd2,coords2)
names(pyrites.xdd2)[40:41] <- c("clng.pal","clat.pal")

coords3 <- reconstruct(pyrites.xdd3[, c("clng", "clat"), drop=FALSE], age= pyrites.xdd3$av_age, enumerate=FALSE, verbose=FALSE)
pyrites.xdd3 <- cbind(pyrites.xdd3,coords3)
names(pyrites.xdd3)[40:41] <- c("clng.pal","clat.pal")

coords4 <- reconstruct(pyrites.xdd4[, c("clng", "clat"), drop=FALSE], age= pyrites.xdd4$av_age, enumerate=FALSE, verbose=FALSE)
pyrites.xdd4 <- cbind(pyrites.xdd4,coords4)
names(pyrites.xdd4)[40:41] <- c("clng.pal","clat.pal")

coords5 <- reconstruct(pyrites.xdd5[, c("clng", "clat"), drop=FALSE], age= pyrites.xdd5$av_age, enumerate=FALSE, verbose=FALSE)
pyrites.xdd5 <- cbind(pyrites.xdd5,coords5)
names(pyrites.xdd5)[40:41] <- c("clng.pal","clat.pal")

coords6 <- reconstruct(pyrites.xdd6[, c("clng", "clat"), drop=FALSE], age= pyrites.xdd6$av_age, enumerate=FALSE, verbose=FALSE)
pyrites.xdd6 <- cbind(pyrites.xdd6,coords6)
names(pyrites.xdd6)[40:41] <- c("clng.pal","clat.pal")

coords7 <- reconstruct(pyrites.xdd7[, c("clng", "clat"), drop=FALSE], age= pyrites.xdd7$av_age, enumerate=FALSE, verbose=FALSE)
pyrites.xdd7 <- cbind(pyrites.xdd7,coords7)
names(pyrites.xdd7)[40:41] <- c("clng.pal","clat.pal")

coords8 <- reconstruct(pyrites.xdd8[, c("clng", "clat"), drop=FALSE], age= pyrites.xdd8$av_age, enumerate=FALSE, verbose=FALSE)
pyrites.xdd8 <- cbind(pyrites.xdd8,coords8)
names(pyrites.xdd8)[40:41] <- c("clng.pal","clat.pal")

coords9 <- reconstruct(pyrites.xdd9[, c("clng", "clat"), drop=FALSE], age= pyrites.xdd9$av_age, enumerate=FALSE, verbose=FALSE)
pyrites.xdd9 <- cbind(pyrites.xdd9,coords9)
names(pyrites.xdd9)[40:41] <- c("clng.pal","clat.pal")

coords10 <- reconstruct(pyrites.xdd10[, c("clng", "clat"), drop=FALSE], age= pyrites.xdd10$av_age, enumerate=FALSE, verbose=FALSE)
pyrites.xdd10 <- cbind(pyrites.xdd10,coords10)
names(pyrites.xdd10)[40:41] <- c("clng.pal","clat.pal")

coords11 <- reconstruct(pyrites.xdd11[, c("clng", "clat"), drop=FALSE], age= pyrites.xdd11$av_age, enumerate=FALSE, verbose=FALSE)
pyrites.xdd11 <- cbind(pyrites.xdd11,coords11)
names(pyrites.xdd11)[40:41] <- c("clng.pal","clat.pal")

pyrites.xdd.av <- rbind(pyrites.xdd0, pyrites.xdd1, pyrites.xdd2, pyrites.xdd3, pyrites.xdd4, pyrites.xdd5,
                     pyrites.xdd6, pyrites.xdd7, pyrites.xdd8, pyrites.xdd9, pyrites.xdd10,
                     pyrites.xdd11)

pyrites.xdd.av$zone2 <- cut(pyrites.xdd.av$av_age,c(540,510,480,455,420,385,380,365,330,300,280,250,100,80))
pyrites.xdd.av <- pyrites.xdd.av[!is.na(pyrites.xdd.av$zone2),]

pyrites.xdd.av <- subset(pyrites.xdd.av, zone2 != "(100,250]")

## range charts ##

## base ages ##

coords0 <- reconstruct(pyrites.xdd0[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd0$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd0 <- cbind(pyrites.xdd0,coords0)
names(pyrites.xdd0)[42:43] <- c("clng.pal_base","clat.pal_base")

coords1 <- reconstruct(pyrites.xdd1[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd1$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd1 <- cbind(pyrites.xdd1,coords1)
names(pyrites.xdd1)[42:43] <- c("clng.pal_base","clat.pal_base")

coords2 <- reconstruct(pyrites.xdd2[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd2$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd2 <- cbind(pyrites.xdd2,coords2)
names(pyrites.xdd2)[42:43] <- c("clng.pal_base","clat.pal_base")

coords3 <- reconstruct(pyrites.xdd3[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd3$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd3 <- cbind(pyrites.xdd3,coords3)
names(pyrites.xdd3)[42:43] <- c("clng.pal_base","clat.pal_base")

coords4 <- reconstruct(pyrites.xdd4[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd4$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd4 <- cbind(pyrites.xdd4,coords4)
names(pyrites.xdd4)[42:43] <- c("clng.pal_base","clat.pal_base")

coords5 <- reconstruct(pyrites.xdd5[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd5$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd5 <- cbind(pyrites.xdd5,coords5)
names(pyrites.xdd5)[42:43] <- c("clng.pal_base","clat.pal_base")

coords6 <- reconstruct(pyrites.xdd6[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd6$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd6 <- cbind(pyrites.xdd6,coords6)
names(pyrites.xdd6)[42:43] <- c("clng.pal_base","clat.pal_base")

coords7 <- reconstruct(pyrites.xdd7[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd7$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd7 <- cbind(pyrites.xdd7,coords7)
names(pyrites.xdd7)[42:43] <- c("clng.pal_base","clat.pal_base")

coords8 <- reconstruct(pyrites.xdd8[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd8$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd8 <- cbind(pyrites.xdd8,coords8)
names(pyrites.xdd8)[42:43] <- c("clng.pal_base","clat.pal_base")

coords9 <- reconstruct(pyrites.xdd9[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd9$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd9 <- cbind(pyrites.xdd9,coords9)
names(pyrites.xdd9)[42:43] <- c("clng.pal_base","clat.pal_base")

coords10 <- reconstruct(pyrites.xdd10[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd10$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd10 <- cbind(pyrites.xdd10,coords10)
names(pyrites.xdd10)[42:43] <- c("clng.pal_base","clat.pal_base")

coords11 <- reconstruct(pyrites.xdd11[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd11$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd11 <- cbind(pyrites.xdd11,coords11)
names(pyrites.xdd11)[42:43] <- c("clng.pal_base","clat.pal_base")

coords12 <- reconstruct(pyrites.xdd12[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd12$b_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd12 <- cbind(pyrites.xdd12,coords12)
names(pyrites.xdd12)[42:43] <- c("clng.pal_base","clat.pal_base")

## top ages ##

coords0 <- reconstruct(pyrites.xdd0[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd0$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd0 <- cbind(pyrites.xdd0,coords0)
names(pyrites.xdd0)[44:45] <- c("clng.pal_top","clat.pal_top")

coords1 <- reconstruct(pyrites.xdd1[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd1$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd1 <- cbind(pyrites.xdd1,coords1)
names(pyrites.xdd1)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords2 <- reconstruct(pyrites.xdd2[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd2$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd2 <- cbind(pyrites.xdd2,coords2)
names(pyrites.xdd2)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords3 <- reconstruct(pyrites.xdd3[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd3$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd3 <- cbind(pyrites.xdd3,coords3)
names(pyrites.xdd3)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords4 <- reconstruct(pyrites.xdd4[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd4$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd4 <- cbind(pyrites.xdd4,coords4)
names(pyrites.xdd4)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords5 <- reconstruct(pyrites.xdd5[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd5$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd5 <- cbind(pyrites.xdd5,coords5)
names(pyrites.xdd5)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords6 <- reconstruct(pyrites.xdd6[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd6$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd6 <- cbind(pyrites.xdd6,coords6)
names(pyrites.xdd6)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords7 <- reconstruct(pyrites.xdd7[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd7$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd7 <- cbind(pyrites.xdd7,coords7)
names(pyrites.xdd7)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords8 <- reconstruct(pyrites.xdd8[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd8$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd8 <- cbind(pyrites.xdd8,coords8)
names(pyrites.xdd8)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords9 <- reconstruct(pyrites.xdd9[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd9$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd9 <- cbind(pyrites.xdd9,coords9)
names(pyrites.xdd9)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords10 <- reconstruct(pyrites.xdd10[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd10$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd10 <- cbind(pyrites.xdd10,coords10)
names(pyrites.xdd10)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords11 <- reconstruct(pyrites.xdd11[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd11$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd11 <- cbind(pyrites.xdd11,coords11)
names(pyrites.xdd11)[44:45] <-  c("clng.pal_top","clat.pal_top")

coords12 <- reconstruct(pyrites.xdd12[, c("clng", "clat"), drop=FALSE], age= c(round(pyrites.xdd12$t_age/5)*5), enumerate=FALSE, verbose=FALSE)
pyrites.xdd12 <- cbind(pyrites.xdd12,coords12)
names(pyrites.xdd12)[44:45] <-  c("clng.pal_top","clat.pal_top")

pyrites.xdd.all <- rbind(pyrites.xdd0, pyrites.xdd1, pyrites.xdd2, pyrites.xdd3, pyrites.xdd4, pyrites.xdd5,
                        pyrites.xdd6, pyrites.xdd7, pyrites.xdd8, pyrites.xdd9, pyrites.xdd10,
                        pyrites.xdd11)

names(pyrites.xdd.all)[30] <- "type"

# plot pyrite samples onto palaeolat panel

pyrites.samples <- pyrites[,-4]

pyrites.samples <- reshape2::dcast(pyrites.samples, Location+Age.Ma+Lat+Long+Phrase+Source+Notes~clusters, value.var = "perc")

pyrites.samples[c("Type.1","Type.2","Type.3","Type.4","Type.5")][is.na(pyrites.samples[c("Type.1","Type.2","Type.3","Type.4","Type.5")])] <- 0

pyrites.samples$Long <- round(pyrites.samples$Long, digits = 1)
pyrites.samples$Lat <- round(pyrites.samples$Lat, digits = 1)

pyrites3 <- pyrites.samples[1:100,]
pyrites4 <- pyrites.samples[101:150,]
pyrites5 <- pyrites.samples[151:nrow(pyrites.samples),]

coords1 <- reconstruct(pyrites3[, c("Long", "Lat"), drop=FALSE], age= round(pyrites3$Age.Ma), enumerate=FALSE, verbose=FALSE)
coords2 <- reconstruct(pyrites4[, c("Long", "Lat"), drop=FALSE], age= round(pyrites4$Age.Ma), enumerate=FALSE, verbose=FALSE)
coords3 <- reconstruct(pyrites5[, c("Long", "Lat"), drop=FALSE], age= round(pyrites5$Age.Ma), enumerate=FALSE, verbose=FALSE)

pyrites3 <- cbind(pyrites3, coords1)
pyrites4 <- cbind(pyrites4, coords2)
pyrites5 <- cbind(pyrites5, coords3)

pyrites.samples <- rbind(pyrites3,pyrites4,pyrites5)

names(pyrites.samples)[13:14] <- c("Pal.Long","Pal.Lat")


a <- ggplot(pyrites.xdd.all) + 
  geom_segment(aes(x = t_age, xend = b_age, y=clat.pal_top, yend = clat.pal_base, colour = type), size = 1) +
  scale_x_reverse(limits = c(543,-1), expand = c(0,0)) + theme_bw() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(breaks = c(-90,-45,0,45,90), limits = c(-90,90), expand = c(0, 0)) +
  xlab("Age (Ma)") + ylab("Palaeolatitude (°)") + 
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_scatterpie(aes(Age.Ma, Pal.Lat, r = 3), colour = "black", data = pyrites.samples, 
                  cols=c("Type.5","Type.2","Type.1","Type.3","Type.4")) +
  geom_vline(xintercept = c(540,510,480,455,420,385,380,365,330,300,280,250,100,80)) +
  geom_rect(aes(xmin = 300, xmax = 420, ymax= 45, ymin = -45), fill = NA, colour = "black", size = 1)

b <- ggplot(pyrites.xdd.all) + 
  geom_segment(aes(x = t_age, xend = b_age, y=clat.pal_top, yend = clat.pal_base, colour = type), size = 1) +
  scale_x_reverse(limits = c(422,300), expand = c(0,0)) + theme_bw() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(breaks = c(-45,-22.5,0,22.5,45), limits = c(-45,45), expand = c(0, 0)) +
  xlab("Age (Ma)") + ylab("Palaeolatitude (°)") + 
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_scatterpie(aes(Age.Ma, Pal.Lat, r = 1.5), colour = "black", data = pyrites.samples, 
                  cols=c("Type.5","Type.2","Type.1","Type.3","Type.4")) +
  geom_vline(xintercept = c(540,510,480,455,420,385,380,365,330,300,280,250,100,80))


dem4$zone2 <- as.character(dem4$zone2)

dem4$zone2[dem4$zone2 == "(80,100]"] <- "80-100 Ma"
dem4$zone2[dem4$zone2 == "(250,280]"] <- "250-280 Ma"
dem4$zone2[dem4$zone2 == "(280,300]"] <- "280-300 Ma"
dem4$zone2[dem4$zone2 == "(300,330]"] <- "300-330 Ma"
dem4$zone2[dem4$zone2 == "(330,365]"] <- "330-365 Ma"
dem4$zone2[dem4$zone2 == "(365,380]"] <- "365-380 Ma"
dem4$zone2[dem4$zone2 == "(380,385]"] <- "380-385 Ma"
dem4$zone2[dem4$zone2 == "(385,420]"] <- "385-420 Ma"
dem4$zone2[dem4$zone2 == "(420,455]"] <- "420-455 Ma"
dem4$zone2[dem4$zone2 == "(455,480]"] <- "455-480 Ma"
dem4$zone2[dem4$zone2 == "(480,510]"] <- "480-510 Ma"
dem4$zone2[dem4$zone2 == "(510,540]"] <- "510-540 Ma"

dem4$zone2 <- factor(dem4$zone2,
                               levels=c("80-100 Ma","250-280 Ma","280-300 Ma", "300-330 Ma", "330-365 Ma", "365-380 Ma", 
                                        "380-385 Ma", "385-420 Ma", "420-455 Ma", "455-480 Ma", "480-510 Ma", "510-540 Ma"))


pyrites.xdd.av$zone2 <- as.character(pyrites.xdd.av$zone2)

pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(80,100]"] <- "80-100 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(250,280]"] <- "250-280 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(280,300]"] <- "280-300 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(300,330]"] <- "300-330 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(330,365]"] <- "330-365 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(365,380]"] <- "365-380 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(380,385]"] <- "380-385 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(385,420]"] <- "385-420 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(420,455]"] <- "420-455 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(455,480]"] <- "455-480 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(480,510]"] <- "480-510 Ma"
pyrites.xdd.av$zone2[pyrites.xdd.av$zone2 == "(510,540]"] <- "510-540 Ma"

pyrites.xdd.av$zone2 <- factor(pyrites.xdd.av$zone2,
                            levels=c("80-100 Ma","250-280 Ma","280-300 Ma", "300-330 Ma", "330-365 Ma", "365-380 Ma", 
                                     "380-385 Ma", "385-420 Ma", "420-455 Ma", "455-480 Ma", "480-510 Ma", "510-540 Ma"))



pyrites.all$zone2 <- as.character(pyrites.all$zone2)

pyrites.all$zone2[pyrites.all$zone2 == "(80,100]"] <- "80-100 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(250,280]"] <- "250-280 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(280,300]"] <- "280-300 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(300,330]"] <- "300-330 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(330,365]"] <- "330-365 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(365,380]"] <- "365-380 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(380,385]"] <- "380-385 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(385,420]"] <- "385-420 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(420,455]"] <- "420-455 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(455,480]"] <- "455-480 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(480,510]"] <- "480-510 Ma"
pyrites.all$zone2[pyrites.all$zone2 == "(510,540]"] <- "510-540 Ma"

pyrites.all$zone2 <- factor(pyrites.all$zone2,
                            levels=c("80-100 Ma","250-280 Ma","280-300 Ma", "300-330 Ma", "330-365 Ma", "365-380 Ma", 
                                     "380-385 Ma", "385-420 Ma", "420-455 Ma", "455-480 Ma", "480-510 Ma", "510-540 Ma"))
names(pyrites.all)[3:7] <- c(1:5)

c <- ggplot(dem4) +
  geom_tile(aes(x = long, y = lat, fill = value)) +
  scale_fill_gradientn(colours = c("dodgerblue4","dodgerblue3","dodgerblue2","deepskyblue1","deepskyblue1","turquoise1","burlywood2","grey70","grey50","grey50"), 
                       values = rescale(c(min(dem3$value),-2000,-1000,-500,-200,0,200,1000,max(dem3$value))),
                       guide = "none", limits=c(min(dem3$value),max(dem3$value))) +
  new_scale("fill") +
  facet_wrap(~fct_rev(zone2)) +
  theme_bw() +
  geom_segment(data = pyrites.all, aes(x = Pal.Long, y = Pal.Lat, xend = Pal.Long.jit, yend = Pal.Lat.jit, group = group)) +
  geom_scatterpie(aes(Pal.Long.jit, Pal.Lat.jit, group = group, r = 10), colour = "black", data = pyrites.all, 
                  cols=c("5","2","1","3","4")) +
  scale_fill_discrete(name = "Pyrite types") +
  new_scale("fill") +
  coord_fixed(ratio = 1) +
  geom_point(data = pyrites.xdd.av, aes(clng.pal, clat.pal, shape = type, fill = type),size = 2, colour = "black") +
  scale_fill_manual(values = c("#F8766D","#00BA38","#619CFF"), name = "Pyrite-bearing rocks") +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01), breaks = c(-180,0,180)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01), breaks = c(-90,0,90)) +
  scale_shape_manual(values = c(21:23), name = "Pyrite-bearing rocks") + 
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Palaeolongitude (°)") + ylab("Palaeolatitude (°)")

lay <- rbind(c(1),
             c(1),
             c(2),
             c(2),
             c(3),
             c(3),
             c(3),
             c(3))

grid.arrange(a,b,c, layout_matrix = lay) # alternative fig

##### main figure (Fig 6) ####

phan.micro.amb <- read.csv("pyrite_redox_ratio_downsampled_phanerozoic.csv")
phan_results.target1.melt.bins <- read.csv("SGP_Fe_redox_fractions_downsampled_phanerozoic.csv")
phan_results.target1.melt.bins <- reshape2::melt(phan_results.target1.melt.bins, id.vars = c("X","edge","span","mid"))
phan_results.target2.melt.bins <- read.csv("SGP_TOCP_redox_fractions_downsampled_phanerozoic.csv")
phan_results.target2.melt.bins <- reshape2::melt(phan_results.target2.melt.bins, id.vars = c("X","edge","span","mid"))

SGP_fractions <- read.csv("SGP_Fe_ferruginous_fractions.csv")
SGP_fractions$age <- round(SGP_fractions$age)
SGP_fractions$long <- round(SGP_fractions$long, digits = 1)
SGP_fractions$lat <- round(SGP_fractions$lat, digits = 1)
SGP_fractions <- subset(SGP_fractions, age < 540)

row.names(SGP_fractions) <- NULL

SGP_fractions0 <- SGP_fractions[1:50,]
SGP_fractions1 <- SGP_fractions[51:100,]
SGP_fractions2 <- SGP_fractions[101:150,]
SGP_fractions3 <- SGP_fractions[151:200,]
SGP_fractions4 <- SGP_fractions[201:250,]
SGP_fractions5 <- SGP_fractions[251:300,]
SGP_fractions6 <- SGP_fractions[301:350,]
SGP_fractions7 <- SGP_fractions[351:400,]
SGP_fractions8 <- SGP_fractions[401:nrow(SGP_fractions),]

coords0 <- reconstruct(SGP_fractions0[, c("long","lat"), drop=FALSE], age= c(round(SGP_fractions0$age/5)*5), enumerate=FALSE, verbose=FALSE)
SGP_fractions0 <- cbind(SGP_fractions0,coords0)
names(SGP_fractions0)[8:9] <- c("long.pal","lat.pal")

coords1 <- reconstruct(SGP_fractions1[, c("long","lat"), drop=FALSE], age= round(SGP_fractions1$age/5)*5, enumerate=FALSE, verbose=FALSE)
SGP_fractions1 <- cbind(SGP_fractions1,coords1)
names(SGP_fractions1)[8:9] <- c("long.pal","lat.pal")

coords2 <- reconstruct(SGP_fractions2[, c("long","lat"), drop=FALSE], age= round(SGP_fractions2$age/5)*5, enumerate=FALSE, verbose=FALSE)
SGP_fractions2 <- cbind(SGP_fractions2,coords2)
names(SGP_fractions2)[8:9] <- c("long.pal","lat.pal")

coords3 <- reconstruct(SGP_fractions3[, c("long","lat"), drop=FALSE], age= round(SGP_fractions3$age/5)*5, enumerate=FALSE, verbose=FALSE)
SGP_fractions3 <- cbind(SGP_fractions3,coords3)
names(SGP_fractions3)[8:9] <- c("long.pal","lat.pal")

coords4 <- reconstruct(SGP_fractions4[, c("long","lat"), drop=FALSE], age= round(SGP_fractions4$age/5)*5, enumerate=FALSE, verbose=FALSE)
SGP_fractions4 <- cbind(SGP_fractions4,coords4)
names(SGP_fractions4)[8:9] <- c("long.pal","lat.pal")

coords5 <- reconstruct(SGP_fractions5[, c("long","lat"), drop=FALSE], age= round(SGP_fractions5$age/5)*5, enumerate=FALSE, verbose=FALSE)
SGP_fractions5 <- cbind(SGP_fractions5,coords5)
names(SGP_fractions5)[8:9] <- c("long.pal","lat.pal")

coords6 <- reconstruct(SGP_fractions6[, c("long","lat"), drop=FALSE], age= round(SGP_fractions6$age/5)*5, enumerate=FALSE, verbose=FALSE)
SGP_fractions6 <- cbind(SGP_fractions6,coords6)
names(SGP_fractions6)[8:9] <- c("long.pal","lat.pal")

coords7 <- reconstruct(SGP_fractions7[, c("long","lat"), drop=FALSE], age= round(SGP_fractions7$age/5)*5, enumerate=FALSE, verbose=FALSE)
SGP_fractions7 <- cbind(SGP_fractions7,coords7)
names(SGP_fractions7)[8:9] <- c("long.pal","lat.pal")

coords8 <- reconstruct(SGP_fractions8[, c("long","lat"), drop=FALSE], age= round(SGP_fractions8$age/5)*5, enumerate=FALSE, verbose=FALSE)
SGP_fractions8 <- cbind(SGP_fractions8,coords8)
names(SGP_fractions8)[8:9] <- c("long.pal","lat.pal")

SGP_fractions.all <- rbind(SGP_fractions0, SGP_fractions1, SGP_fractions2, SGP_fractions3, SGP_fractions4, SGP_fractions5,
                         SGP_fractions6, SGP_fractions7, SGP_fractions8)

names(SGP_fractions.all)[c(4,8:9)] <- c("Age.Ma","Pal.Long","Pal.Lat")

a <- ggplot() +
  #geom_stepribbon(data = phan.micro.amb,
  #aes(x = edge, ymin = lower, ymax = upper), fill = "grey50") +
  geom_step(data = phan.micro.amb, aes(y =ratio, x = edge),
            direction = "hv", colour = "black", size = 0.75) +
  theme_bw() + scale_x_reverse(limits = c(550,0), expand = c(0,0)) +
  geom_step(data = subset(phan_results.target2.melt.bins, variable == "ferruginous"), 
            aes(y = value, x = edge, group = variable), 
            position = "stack", direction = "hv", size = 0.75, colour = "red") +
  geom_step(data = subset(phan_results.target1.melt.bins, variable == "ferruginous"), 
            aes(y = value, x = edge, group = variable), 
            position = "stack", direction = "hv", size = 0.75, colour = "blue") +
  geom_vline(xintercept = c(540,510,480,455,420,385,380,365,330,300,280,250,100,80)) +
  ylab("Fraction ferruginous") +
  coord_cartesian(ylim = c(0,1), clip = "off") +
  labs(tag = "A") +
  theme(plot.tag.position = c(0.01,0.99))

combined <- rbind(pyrites.samples[,c("Age.Ma","Pal.Long","Pal.Lat")],SGP_fractions.all[,c("Age.Ma","Pal.Long","Pal.Lat")])
combined <- combined[complete.cases(combined),]

# range plot

pyrites.samples$fraction <- (pyrites.samples$Type.5+pyrites.samples$Type.2+pyrites.samples$Type.1)/(pyrites.samples$Type.5+pyrites.samples$Type.2+pyrites.samples$Type.1+pyrites.samples$Type.3+pyrites.samples$Type.4)
names(SGP_fractions.all)[7] <- "fraction"

b <- ggplot(pyrites.xdd.all) + 
  #stat_density_2d(data = combined, aes(Age.Ma, Pal.Lat, fill = ..density.., alpha = ..density..),
                  #h = NULL, adjust = 2, n = 54, geom = "raster", contour = FALSE) +
  #geom_bin2d(data = combined, aes(Age.Ma, Pal.Lat, fill = ..density..), bins = 54) +
  #scale_fill_distiller(palette = "Spectral",direction = -1, limits = c(0.25,0,75), oob = squish) +
  geom_segment(aes(x = t_age, xend = b_age, y=clat.pal_top, yend = clat.pal_base, colour = type), size = 1) +
  scale_colour_manual(values = c("green","black","#3288BD")) +
  scale_x_reverse(limits = c(540,1), expand = c(0,0)) + theme_bw() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(breaks = c(-90,-45,0,45,90), limits = c(-90,90), expand = c(0, 0)) +
  xlab("Age (Ma)") + ylab("Palaeolatitude (°)") + 
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.tag.position = c(0.01,0.99)) +
  geom_point(aes(Age.Ma, Pal.Lat, fill = fraction), colour = "black", data = pyrites.samples, pch = 21, size = 2) +
  geom_point(aes(Age.Ma, Pal.Lat, fill = fraction), colour = "black", data = SGP_fractions.all, pch = 22, size = 2) +
  geom_vline(xintercept = c(540,510,480,455,420,385,380,365,330,300,280,250,100,80)) +
  scale_fill_distiller(palette = "Spectral", limits = c(0.25,0.75), oob = squish) +
  labs(tag = "B") 

grid.arrange(a,b,ncol = 1)

# pyrite TE spatial proximities and weights (sensu Mehra et al. 2021) - for insertion into Fig 4

XY <- as.matrix(pyrites.samples[,c(4,3)])

degree_bin <- 0.5
t_t <- as.data.frame(apply(XY, 1, FUN=function(X) distHaversine(X, XY)))/degree_bin
test <- 1/((t_t^2)+1)
test <- as.data.frame(test)
prox_t <- rowSums(test)/2
pyrites.samples$prox2 <- prox_t
pyrites.samples$weights2  <- 1/((prox_t*median(2/prox_t))+1)

pyrites.samples.labs <- pyrites.samples[!duplicated(pyrites.samples$Location),]
pyrites.samples.labs$index <- seq(1,nrow(pyrites.samples.labs), by = 1)

# Fig. 4F overlay

ggplot(pyrites.samples, aes(Age.Ma, prox2)) + geom_point() +
  scale_x_reverse(limits = c(540,0)) + theme_bw() #+  geom_text_repel(data = pyrites.samples.labs, aes(Age.Ma, prox2, label = index), size = 2, max.overlaps = 2000)

# Fig. S23

ggplot(pyrites.samples, aes(Age.Ma, prox2)) + geom_point() +
  scale_x_reverse(limits = c(540,0)) + theme_bw() + 
  scale_y_continuous(limits = c(0,30)) +
  geom_text_repel(data = pyrites.samples.labs, aes(Age.Ma, prox2, label = index), size = 5, max.overlaps = 2000) +
  coord_geo(dat = list("periods", "eras"), xlim = c(540, 0), size = 5, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr")) +
  ylab("Spatial proximity") +
  xlab("Age (Ma)")

#write.csv(pyrites.samples.labs, "pyrites.samples.labs.csv")


# idw - fraction ferruginous

pyrites.xdd.all2 <- subset(pyrites.xdd.all, type == "Nodules" | type == "Both")
pyrites.xdd.all2 <- pyrites.xdd.all2[,c(31,32,37:45)]
pyrites.xdd.all2 <- reshape2::melt(pyrites.xdd.all2, id.vars = c(1:5,7,9,11))
pyrites.xdd.all2 <- pyrites.xdd.all2[,-9]
names(pyrites.xdd.all2)[9] <- "Pal.Long"

pyrites.xdd.all2 <- reshape2::melt(pyrites.xdd.all2, id.vars = c(1:5,9))
pyrites.xdd.all2 <- pyrites.xdd.all2[,-7]
names(pyrites.xdd.all2)[7] <- "Pal.Lat"

pyrites.xdd.all2$fraction <- 0.25
pyrites.xdd.all2 <- pyrites.xdd.all2[,-c(3:4)]
pyrites.xdd.all2 <- reshape2::melt(pyrites.xdd.all2, id.vars = c(4:6))
pyrites.xdd.all2 <- pyrites.xdd.all2[,-4]
names(pyrites.xdd.all2)[4] <- "Age.Ma"

combined2 <- rbind(pyrites.samples[,c("Age.Ma","Pal.Long","Pal.Lat","fraction")],
                  SGP_fractions.all[,c("Age.Ma","Pal.Long","Pal.Lat","fraction")],
                  pyrites.xdd.all2[,c("Age.Ma","Pal.Long","Pal.Lat","fraction")])

combined2 <- combined2[complete.cases(combined2),]

combined2 <- subset(combined2, Age.Ma <= 540)

# age weights for composite dataset (sensu Mehra et al. 2021) - use to help balance idw

age_bin <- 10 # span_age parameter
t_t <- as.data.frame(outer(combined2$Age.Ma, combined2$Age.Ma, `-`))/age_bin
t_t <- ((t_t^2)^0.5) # positive distance matrix
test <- 1/((t_t^2)+1)
test <- as.data.frame(test)
prox_t <- rowSums(test)
combined2$prox <- prox_t
combined2$weights  <- 1/((prox_t*median(2/prox_t))+1)

# spatial weights (sensu Mehra et al. 2021) - use to help balance idw

XY <- as.matrix(combined2[,2:3])

degree_bin <- 0.5
t_t <- as.data.frame(apply(XY, 1, FUN=function(X) distHaversine(X, XY)))/degree_bin
test <- 1/((t_t^2)+1)
test <- as.data.frame(test)
prox_t <- rowSums(test)/2
combined2$prox2 <- prox_t
combined2$weights2  <- 1/((prox_t*median(2/prox_t))+1)

combined2$prox_total <- combined2$prox+combined2$prox2
combined2$weight_total <- 1/((combined2$prox_total*median(3/combined2$prox_total))+1)

# visualise weights

xxx <- ggplot(combined2, aes(weight_total)) + geom_histogram() +
  geom_vline(xintercept = 0.2, colour = "red") +
  theme_bw()


combined2$zone2 <- cut(combined2$Age.Ma,c(540,510,480,455,420,385,380,365,330,300,280,250,100,80))
combined2$zone2 <- as.character(combined2$zone2)

combined2$zone2[combined2$zone2 == "(80,100]"] <- "80-100 Ma"
combined2$zone2[combined2$zone2 == "(250,280]"] <- "250-280 Ma"
combined2$zone2[combined2$zone2 == "(280,300]"] <- "280-300 Ma"
combined2$zone2[combined2$zone2 == "(300,330]"] <- "300-330 Ma"
combined2$zone2[combined2$zone2 == "(330,365]"] <- "330-365 Ma"
combined2$zone2[combined2$zone2 == "(365,380]"] <- "365-380 Ma"
combined2$zone2[combined2$zone2 == "(380,385]"] <- "380-385 Ma"
combined2$zone2[combined2$zone2 == "(385,420]"] <- "385-420 Ma"
combined2$zone2[combined2$zone2 == "(420,455]"] <- "420-455 Ma"
combined2$zone2[combined2$zone2 == "(455,480]"] <- "455-480 Ma"
combined2$zone2[combined2$zone2 == "(480,510]"] <- "480-510 Ma"
combined2$zone2[combined2$zone2 == "(510,540]"] <- "510-540 Ma"

combined2$zone2 <- factor(combined2$zone2,
                            levels=c("80-100 Ma","250-280 Ma","280-300 Ma", "300-330 Ma", "330-365 Ma", "365-380 Ma", 
                                     "380-385 Ma", "385-420 Ma", "420-455 Ma", "455-480 Ma", "480-510 Ma", "510-540 Ma"))

# idw 

x.range <- as.integer(c(0,540))
y.range <- as.integer(c(-90,90))

bin <- 10

grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=bin), y=seq(from=y.range[1], to=y.range[2], by=bin))
coordinates(grd) <- ~ x+y
gridded(grd) <- TRUE

slices <- slice_sample(combined2, n = 1000000, replace = TRUE, weight_by = combined2$weight_total)

slices[,2:3] <- jitter2d(slices[,2:3], max = 0.01)

coords <- slices

coordinates(coords) = ~Age.Ma+Pal.Lat

fraction.idw <- idw(formula=fraction ~ 1, locations=coords, newdata=grd, idp = 2)

fraction.idw.output = as.data.frame(fraction.idw)
names(fraction.idw.output)[1:3]<-c("Age.Ma","Pal.Lat","fraction")

yyy <- ggplot(slices, aes(weight_total)) + geom_histogram() +
  geom_vline(xintercept = 0.2, colour = "red") +
  theme_bw()

#grid.arrange(xxx,yyy, ncol = 2) # Fig. S24

c <- ggplot(pyrites.xdd.all) + 
  geom_tile(data=fraction.idw.output,aes(x=Age.Ma, y = Pal.Lat, fill=fraction)) +
  #geom_contour_filled(data = fraction.idw.output,  
                    #aes(x=Age.Ma, y = Pal.Lat, z=fraction)) +
  #geom_segment(aes(x = t_age, xend = b_age, y=clat.pal_top, yend = clat.pal_base, colour = type), size = 1) +
  scale_x_reverse(limits = c(540,1), expand = c(0,0)) + theme_bw() +
  scale_fill_distiller(palette = "Spectral", limits = c(0.25,0.75), oob = squish) +
  #scale_fill_gradientn(colours = c("blue","blue", "white","white","red", "red"),
                       #values = scales::rescale(c(0,0.1, 0.4, 0.6,0.9,1))) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(breaks = c(-90,-45,0,45,90), limits = c(-90,90), expand = c(0, 0)) +
  xlab("Age (Ma)") + ylab("Palaeolatitude (°)") + 
  #theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #geom_point(aes(Age.Ma, Pal.Lat), colour = "black", data = pyrites.samples) +
  #geom_point(aes(Age.Ma, Pal.Lat), colour = "black", data = SGP_fractions.all, pch = 21) +
  geom_vline(xintercept = c(540,510,480,455,420,385,380,365,330,300,280,250,100,80)) +
  theme(legend.position = c(0.95, 0.22), plot.tag.position = c(0.01,0.99)) +
  coord_geo(dat = list("periods", "eras"), xlim = c(540, 0), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr")) +
  labs(tag = "C") 

grid.arrange(a,b,c,ncol=1)

# idw by group

intervals <- as.data.frame(unique(combined2$zone2))
intervals <- as.data.frame(intervals[!is.na(intervals)])
names(intervals)[1] <- "interval"

fraction.idw.output2 <- data.frame()

for (i in 1:nrow(intervals)) {
  x.range <- as.integer(c(-180,180))
  y.range <- as.integer(c(-90,90))
  bin <- 10
  grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=bin), y=seq(from=y.range[1], to=y.range[2], by=bin))
  coordinates(grd) <- ~ x+y
  gridded(grd) <- TRUE
  slices.sub <- subset(slices, zone2 == paste(intervals[i,]))
  coords <- slices.sub
  coordinates(coords) = ~Pal.Long+Pal.Lat
  fraction.idw <- idw(formula=fraction ~ 1, locations=coords, newdata=grd, idp = 2)
  fraction.idw.output = as.data.frame(fraction.idw)
  names(fraction.idw.output)[1:3]<-c("Pal.Long","Pal.Lat","fraction")
  fraction.idw.output <- fraction.idw.output[,-4]
  test <- reshape2::dcast(fraction.idw.output, Pal.Lat~Pal.Long)
  test.lats <- test[,1]
  test.longs <- colnames(test)
  test <- test[,-1]
  test <- as.matrix(test)
  test <- raster(test, xmn = -180, xmx = 180, ymn = -90, ymx = 90)
  combined2.sub <- subset(combined2, zone2 == paste(intervals[i,]))
  out <- raster::extract(test, combined2.sub[,2:3], buffer=2500000, cellnumbers = T)
  xy <- xyFromCell(test, cell = do.call(rbind, out)[,1])
  xy[,2] <- xy[,2]*(-1)
  xy <- SpatialPoints(xy)
  test <- raster::mask(test, xy, invers  = FALSE)
  test <- as.matrix(test)
  test <- cbind(test.lats, test)
  colnames(test) <- test.longs
  test <- as.data.frame(test)
  fraction.idw.output <- reshape2::melt(test, id.vars = "Pal.Lat")
  names(fraction.idw.output)[1:3]<-c("Pal.Lat","Pal.Long","fraction")
  fraction.idw.output$Pal.Long <- as.numeric(as.character(fraction.idw.output$Pal.Long))
  fraction.idw.output$zone2 <- paste(intervals[i,])
  fraction.idw.output2 <- rbind(fraction.idw.output2, fraction.idw.output)
}

# cut each dataframe by 'zone' (i.e., selected flattened interval)

# SGP (Fe speciation)

combined3 <- subset(combined2, zone2 != "(100,250]")

SGP_fractions.all2 <- SGP_fractions.all

SGP_fractions.all2$zone2 <- cut(SGP_fractions.all2$Age.Ma,c(540,510,480,455,420,385,380,365,330,300,280,250,100,80))
SGP_fractions.all2$zone2 <- as.character(SGP_fractions.all2$zone2)

SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(80,100]"] <- "80-100 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(250,280]"] <- "250-280 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(280,300]"] <- "280-300 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(300,330]"] <- "300-330 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(330,365]"] <- "330-365 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(365,380]"] <- "365-380 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(380,385]"] <- "380-385 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(385,420]"] <- "385-420 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(420,455]"] <- "420-455 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(455,480]"] <- "455-480 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(480,510]"] <- "480-510 Ma"
SGP_fractions.all2$zone2[SGP_fractions.all2$zone2 == "(510,540]"] <- "510-540 Ma"

SGP_fractions.all2$zone2 <- factor(SGP_fractions.all2$zone2,
                          levels=c("80-100 Ma","250-280 Ma","280-300 Ma", "300-330 Ma", "330-365 Ma", "365-380 Ma", 
                                   "380-385 Ma", "385-420 Ma", "420-455 Ma", "455-480 Ma", "480-510 Ma", "510-540 Ma"))

SGP_fractions.all2 <- subset(SGP_fractions.all2, zone2 != "(100,250]")

# Pyrite TEs

pyrites.samples2 <- pyrites.samples

pyrites.samples2$zone2 <- cut(pyrites.samples2$Age.Ma,c(540,510,480,455,420,385,380,365,330,300,280,250,100,80))
pyrites.samples2$zone2 <- as.character(pyrites.samples2$zone2)

pyrites.samples2$zone2[pyrites.samples2$zone2 == "(80,100]"] <- "80-100 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(250,280]"] <- "250-280 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(280,300]"] <- "280-300 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(300,330]"] <- "300-330 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(330,365]"] <- "330-365 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(365,380]"] <- "365-380 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(380,385]"] <- "380-385 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(385,420]"] <- "385-420 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(420,455]"] <- "420-455 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(455,480]"] <- "455-480 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(480,510]"] <- "480-510 Ma"
pyrites.samples2$zone2[pyrites.samples2$zone2 == "(510,540]"] <- "510-540 Ma"

pyrites.samples2$zone2 <- factor(pyrites.samples2$zone2,
                                   levels=c("80-100 Ma","250-280 Ma","280-300 Ma", "300-330 Ma", "330-365 Ma", "365-380 Ma", 
                                            "380-385 Ma", "385-420 Ma", "420-455 Ma", "455-480 Ma", "480-510 Ma", "510-540 Ma"))

pyrites.samples2 <- subset(pyrites.samples2, zone2 != "(100,250]")

# xDD

pyrites.xdd.all3 <- pyrites.xdd.all2

pyrites.xdd.all3$zone2 <- cut(pyrites.xdd.all3$Age.Ma,c(540,510,480,455,420,385,380,365,330,300,280,250,100,80))
pyrites.xdd.all3$zone2 <- as.character(pyrites.xdd.all3$zone2)

pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(80,100]"] <- "80-100 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(250,280]"] <- "250-280 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(280,300]"] <- "280-300 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(300,330]"] <- "300-330 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(330,365]"] <- "330-365 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(365,380]"] <- "365-380 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(380,385]"] <- "380-385 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(385,420]"] <- "385-420 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(420,455]"] <- "420-455 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(455,480]"] <- "455-480 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(480,510]"] <- "480-510 Ma"
pyrites.xdd.all3$zone2[pyrites.xdd.all3$zone2 == "(510,540]"] <- "510-540 Ma"

pyrites.xdd.all3$zone2 <- factor(pyrites.xdd.all3$zone2,
                                 levels=c("80-100 Ma","250-280 Ma","280-300 Ma", "300-330 Ma", "330-365 Ma", "365-380 Ma", 
                                          "380-385 Ma", "385-420 Ma", "420-455 Ma", "455-480 Ma", "480-510 Ma", "510-540 Ma"))

pyrites.xdd.all3 <- subset(pyrites.xdd.all3, zone2 != "(100,250]")

land1 <- subset(dem4, value > 100) 

d <- ggplot(land1) +
  geom_tile(data = fraction.idw.output2, aes(x = Pal.Long, y = Pal.Lat, fill = fraction)) +
  scale_fill_distiller(palette = "Spectral", limits = c(0.25,0.75), oob = squish) +
  new_scale_fill() +
  #stat_density_2d(data = combined3, aes(Pal.Long, Pal.Lat, fill = (..density../..n..)*100),
                  #h = NULL, adjust = 3, geom = "raster", contour = FALSE, n = 80) +
  geom_tile(aes(x = long, y = lat), fill = "white") +
  facet_wrap(~fct_rev(zone2)) +
  theme_bw() +
  #scale_fill_gradientn(colours = c("grey80","grey80",NA,NA),
                       #values = scales::rescale(c(0,0.0000001,0.000001, 0.01))) +
  geom_point(data = SGP_fractions.all2, aes(Pal.Long, Pal.Lat), pch = 16, size = 0.2) +
  geom_point(data = pyrites.samples2, aes(Pal.Long, Pal.Lat), pch = 16, size = 0.2) +
  geom_point(data = pyrites.xdd.all3, aes(Pal.Long, Pal.Lat), pch = 16, size = 0.2) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits = c(-180,180), expand = c(0.01, 0.01), breaks = c(-180,0,180)) + 
  scale_y_continuous(limits = c(-90,90), expand = c(0.01, 0.01), breaks = c(-90,0,90)) +
  xlab("Palaeolongitude (°)") + ylab("Palaeolatitude (°)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.tag.position = c(0.01,0.99)) +
  labs(tag = "D") 

lay <- rbind(c(1),
             c(1),
             c(2),
             c(2),
             c(3),
             c(3),
             c(4),
             c(4),
             c(4),
             c(4))

grid.arrange(a,b,c,d, layout_matrix = lay) # main fig (Fig. 6)

###### end #####
