##### This R script was written by Joe Emmings (British Geological Survey) ####
##### This R script implements data manipulation and statistical analysis #####
##### as in Figures 2-6 and Figs. S11-S19 #####

##### KEY REFERENCES ####

# Figure 2 utilises data from the Sedimentary Geochemistry and Palaeoenvironments Project (SGP) http://sgp-search.io/
# See Mehra et al. 2021, GSA Today, 31, https://doi.org/10.1130/GSATG484A.1. and Farrell et al., Geobiology (in review) (2021).

# Figures 3-6 utilise pyrite trace element data accessible here: https://doi.org/10.1130/GEOL.S.12456332
# The data derive from:
# Large et al., Earth and Planetary Science Letters 389, 209-220 (2014).
#	Large et al., Gondwana Research 28, 1282-1293 (2015).
# Large et al., Earth and Planetary Science Letters 428, 139-150 (2015).
# Large et al., Mineralium Deposita 54, 485-506 (2019).
# Mukherjee and Large, Geology 48, 1018-1022 (2020).

# install packages

library(ggplot2)
library(reshape2)
library(ggbiplot)
library(BBmisc)
library(gridExtra)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(useful)
library(compositions)
library(data.table)
library(boot)
library(fANCOVA)
library(colorspace)

# set working directory

setwd("N:/Papers/Pyrite_geochemistry_revisited")

##### Figure 1 - text mining outputs ####

# install.R, preparation.R, analysis.R

##### Figure 2 - SGP violin plots #####

# import redox stage definitions derived from Figure 1

zones <- read.delim("redox_zones.txt")
zone.names <- zones$Zone
zones <- zones$Break

##### Figure 2 ####

SGP <- read.csv("SGP2.csv") # import SGP siliciclastic sediments in focal area

# the data are accessible via the following API call
#{"type":"samples","filters":{"country":["North America","United Kingdom","United States","Canada","Australia","New Zealand"],"lithology_type":["siliciclastic","carbonate","organic"],"lithology_class":["sedimentary"]},"show":["fe","fe_carb","fe_ox","fe_mag","fe_py","tic","toc","tot_c","del_13c_carb","del_13c_org","tmax","s2","s1","s3","su","s_py","s_org","del_34s_py","del_34s_cas","del_34s_gyp","del_34s_obs","del_34s_bulk","n","del_15n","alu","ars","cu","mo","mn","ni","p","u","v","zr","interpreted_age","fe_hr"]}zones <- read.delim("redox_zones2.txt") # use v2 for more conservative zone A def

# add redox stage definitions to SGP data

SGP$zone <- cut(SGP$interpreted.age,c(zones))
levels(SGP$zone) <- c(paste(rev(zone.names)))
SGP[,ncol(SGP)][is.na(SGP$zone)] <- "Meso"

# subset for anoxic facies
SGP.sub <- subset(SGP, FeHR/Fe..wt.. >= 0.38 & Fe..wt.. >= 0.5)

a <- ggplot(SGP.sub, aes(x=zone, y=Fe.py..wt../FeHR, fill=zone)) + 
  geom_violin(scale = "width", draw_quantiles = c(0.25, 0.75)) +
  theme_bw() + 
  scale_x_discrete(limits = rev(levels(SGP$zone))) +
  stat_summary(fun="median", geom="point", aes(group = zone)) +
  stat_summary(fun="median", geom="line", group = 1) +
  scale_y_continuous(limits = c(0,1)) + geom_hline(yintercept = c(0.7,0.8))

b <- ggplot(SGP.sub, aes(x=zone, y=TOC..wt../P..ppm., fill=zone)) + 
  geom_violin(scale = "width", draw_quantiles = c(0.25, 0.75)) +
  theme_bw() + 
  scale_x_discrete(limits = rev(levels(SGP$zone))) +
  stat_summary(fun="median", geom="point", aes(group = zone)) +
  stat_summary(fun="median", geom="line", group = 1) +
  scale_y_log10()

grid.arrange(a,b, ncol= 1)

# tidy using graphics package

##### Figure 3A - Hierarchical cluster analysis (HCA) ####

# download the pyrite trace element dataset of Mukherjee and Large (2020) (CC-BY-NC)
# https://doi.org/10.1130/GEOL.S.12456332 

pyrites <- read.delim("Mukherjee_and_Large2.txt") # manually culled for shales in the focal area

# add redox stages derived from Fig. 1

pyrites$zone <- cut(pyrites$Age.Ma,c(zones))
levels(pyrites$zone) <- c(paste(rev(zone.names)))
pyrites[,21][is.na(pyrites$zone)] <- "Meso"

HCA <- pyrites[4:20]

# clr trans
HCA2 <- clr(HCA)  # replace clr with log transform for Fig. S11
HCA2 <- as.data.frame(HCA2) 
HCA3 <- as.data.frame(t(HCA2))

HCA2 <- BBmisc::normalize(HCA2, method = "standardize") # scale and centre for HCA matrix

# cluster
hr <- hclust(dist(HCA2, method = "euclidean"), method = "ward.D2")
hc <- hclust(dist(HCA3, method = "euclidean"), method = "ward.D2")

#plot(hr) # display dendrogram
#rect.hclust(hr, k=5, border="red") 
clusters <- cutree(hr, k=6)

#plot(hc) # display dendrogram
#rect.hclust(hc, k=4, border="red") 

pyrites <- cbind(pyrites, clusters)

col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))

# output
labels <- pyrites[4:22]
my_group <- as.numeric(as.factor(substr(labels$zone, 1 , 2))) # select zone or cluster
colSide <- brewer.pal(6, "Set3")[my_group] # select 10 or 6

# Fig. 3A

heatmap(as.matrix(HCA2), Rowv = as.dendrogram(hr), Colv = as.dendrogram(hc), 
        col = col, RowSideColors=colSide, labRow = labels$zone)

# for colour scale
pheatmap(as.matrix(HCA2), cluster_rows = hr, 
         cluster_cols = hc)

# tidy using graphics package

##### Figure 3B - Principal Component Analysis (PCA) ####

# select variables
myvars <- c("Se", "Zn", "Cd", "Mo", "Tl", "Bi", "Te", "Mn", "Cu", "Au", "Ag", "Sb") # selected elements - include Co, Ni, As, Pb for Fig. S12
pca.dat <- pyrites[myvars]
pca.dat <- cbind(pyrites[,1:3],pca.dat, pyrites[,21:22])
pca.dat <- pca.dat[complete.cases(pca.dat[,]),]

# clr transform
pca.dat[,4:(ncol(pca.dat)-2)] <- clr(pca.dat[,4:(ncol(pca.dat)-2)]) # replace clr with log transform for Fig. S11

labels <- pca.dat[,2]
age <- pca.dat[,3]
groups <- pca.dat[,ncol(pca.dat)-1]
HCA.clusters <- pca.dat[,ncol(pca.dat)]
HCA.clusters2 <- as.character(HCA.clusters)

pca.dat <- pca.dat[,4:(ncol(pca.dat)-2)]

# compute principal components
PCAdata <- prcomp(pca.dat,
                  center = TRUE,
                  scale. = FALSE) # covariance matrix

summary(PCAdata)
plot(PCAdata, type = "l") # scree plot

# print PCA biplot (Fig 3B)

g <- ggbiplot(PCAdata, obs.scale = 1, var.scale = 1, 
              ellipse = TRUE, ellipse.prob = 0.68,
              circle = FALSE, choices = c(1,2), groups = groups) # select clusters or groups
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top') + theme_bw() + geom_point(aes(shape = groups, fill = groups)) +
  scale_fill_manual(values = c("black", "red", "orange", "blue", "green", "grey")) +
  scale_shape_manual(values = c(3,21:25,4))

print(g) 

# tidy using graphics package

##### Figure 3C - total concentration boxplots #####

pyrites$totals <- rowSums(pyrites[,c("Mn","Co","Ni","Cu","Zn","As","Se","Mo","Ag","Cd","Sb","Te","Pt","Au","Tl","Pb","Bi" )])

# boxplots

pyrites3 <- subset(pyrites, clusters != "Cluster High-Au") # remove High-Au group

# print - Fig. 3C

ggplot(pyrites3, aes(clusters, totals)) + 
  geom_boxplot(aes(group = clusters), notch = TRUE, outlier.shape = NA, coef = 0) + theme_bw() +
  stat_summary(fun="median", geom="point", aes(group = clusters)) +
  coord_cartesian(ylim=c(0,12000))

##### Figure 4 - PCA polar coordinate histograms #####

# extract PCA biplot cartersian coordinates
coords <- as.data.frame(PCAdata[["x"]])
coords <- coords[,1:2]
coords2 <- cbind(coords, groups)

# convert cartesian to polar coords

poles <- cart2pol(coords2$PC1, coords2$PC2, degrees = TRUE)
poles <- cbind(poles, groups, labels, HCA.clusters, age)

# meaningful cluster names

poles$HCA.clusters[poles$HCA.clusters == 1] <- "Cluster 1,6" 
poles$HCA.clusters[poles$HCA.clusters == 2] <- "Cluster 5"
poles$HCA.clusters[poles$HCA.clusters == 3] <- "Cluster 2" 
poles$HCA.clusters[poles$HCA.clusters == 4] <- "Cluster 4" 
poles$HCA.clusters[poles$HCA.clusters == 5] <- "Cluster High-Au" 
poles$HCA.clusters[poles$HCA.clusters == 6] <- "Cluster 3" 

poles$theta <- poles$theta+4.78 # rotate (cluster 1 centred on Cd)

poles1 <- subset(poles, theta <= 360)
poles2 <- subset(poles, theta > 360)

poles2$theta <- poles2$theta-360

poles <- rbind(poles1, poles2)

# weighted histograms

# by HCA cluster (Fig. 4, left hand side)

ggplot(subset(poles, HCA.clusters != "Cluster High-Au"), aes(theta)) + 
  geom_histogram(binwidth = 15, aes(weight = r, fill = HCA.clusters), center = 7.5) + 
  theme_bw() + coord_flip() +
  facet_wrap(~HCA.clusters, ncol = 6, scale="free_x") + scale_x_reverse() +
  geom_vline(xintercept = c(0,45,135,225,270,315)) +
  theme(legend.position = "none") #+ geom_density(aes(y=..density..*20000)) # optional density line

# by redox stage (Fig. 4, right hand side)

ggplot(poles, aes(theta)) + 
  geom_histogram(binwidth = 15, aes(weight = r, fill = groups), center = 7.5) + 
  theme_bw() + coord_flip() +
  facet_wrap(~groups, ncol = 6, scale="free_x") + scale_x_reverse() +
  geom_vline(xintercept = c(0,45,135,225,270,315))

# coloured by location (Figs. S12-S16)
ggplot(poles, aes(theta)) + 
  geom_histogram(binwidth = 15, aes(weight = r, fill = groups), center = 7.5) + 
  theme_bw() + coord_flip() +
  facet_wrap(~labels, ncol = 6, scale="free_x") + scale_x_reverse() +
  geom_vline(xintercept = c(0,45,135,225,270,315))

# tidy all plots using graphics package

##### Figure 5 ####

# as factor - i.e., free and collapsed scale
poles$age2 <- as.factor(poles$age)

Periods <- as.factor(c(2500, 1600, 1000, 720, 635, 541, 485.4, 443.8, 419.2, 358.9, 298.9, 
                       251.9, 201.3, 145, 66, 23.03, 2.58))

quartile1 <- function(x){
  out <- quantile(x, probs = c(0.25))
  names(out) <- c("ymin")
  return(out) 
}

quartile2 <- function(x){
  out <- quantile(x, probs = c(0.75))
  names(out) <- c("ymin")
  return(out) 
}

# age difference bar (top part of Fig. 5A)

diffs <-  age
diffs <- as.data.frame(diffs[order(diffs)])

names(diffs)[1] <- "Age"

diffs <- as.data.frame(diffs[!duplicated(diffs$Age),])

diffs_df <- diffs[-1,] - diffs[-nrow(diffs),]
diffs_df <- as.data.frame(diffs_df)

diffs <- cbind(diffs[1:201,],diffs_df)
names(diffs)[1] <- "Age"
names(diffs)[2] <- "Diff"

diffs$Age <- as.factor(diffs$Age)

ggplot(diffs, aes(Age, Diff)) + geom_line(group = 1) + theme_bw() + scale_y_log10()

# for main panel on Fig. 5A

xdd <- read.csv("xdd_binned_results.csv") # import text mining outputs - use write.csv following Figure 1B output in analysis.R

xdd1 <- data.table(xdd)
poles1 <- data.table(poles)

setkey(xdd1, Age)
setkey(poles1, age)

out <- xdd1[poles1, roll = "nearest"] # join to nearest datum

ggplot(poles, aes(age2, theta)) + 
  scale_y_reverse() + 
  scale_x_discrete(limits = rev(levels(poles$age2)), guide = guide_axis(n.dodge=2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  stat_summary(fun=quartile1, geom="line", size = 0.75, group = 1, colour = "red") +
  stat_summary(fun=quartile2, geom="line", size = 0.75, group = 1, colour = "red") +
  stat_summary(fun=median, geom="line", size = 0.75, group = 1) +
  geom_hline(yintercept = c(0,45,135,225,270,315)) +
  geom_step(data = out, aes(x = age2, y = (nodules+framboids)/(framboids+nodules+undif)*800), 
            group = 1, colour = "red")

# Boxplots by unit - Fig. S17

Periods <- c(2500, 1600, 1000, 720, 635, 541, 485.4, 443.8, 419.2, 358.9, 298.9, 
                      251.9, 201.3, 145, 66, 23.03, 2.58) # add Period boundaries (ICS)

HEATT <- c(56, 66, 93, 116, 183, 200, 251, 359, 372, 383, 444, 514, 542) 

poles.sub1 <- subset(poles, age > 681)
poles.sub1$bin <- cut(poles.sub1$age, seq(from = 681, to = 3520, by = 20)) # 20 Ma bins

a <- ggplot(poles.sub1, aes(age, theta)) + 
  geom_boxplot(aes(group = bin, fill = groups, colour = groups), width = 5, size = 1) + 
  scale_y_reverse(sec.axis = sec_axis(~./800)) + scale_x_reverse(limits = c(3520, 676)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  geom_hline(yintercept = c(0,45,135,225,270,315)) +
  geom_step(data = xdd, aes(x = Age, y = (nodules)/(framboids+nodules+undif)*800), 
            group = 1, colour = "blue") +
  geom_step(data = xdd, aes(x = Age, y = (nodules+framboids)/(framboids+nodules+undif)*800), 
             group = 1, colour = "red") +
  geom_vline(xintercept = c(Periods)) +  geom_vline(xintercept = c(HEATT), colour = "green") +
  geom_vline(xintercept = c(zones), colour = "orange") +
  scale_fill_manual(values=c("#f564e3")) + scale_colour_manual(values=c("#f564e3"))

poles.sub2 <- subset(poles, age <= 681)
poles.sub2$bin <- cut(poles.sub2$age, seq(from = 0, to = 681, by = 5)) # 5 Ma bins

b <- ggplot(poles.sub2, aes(age, theta)) + 
  geom_boxplot(aes(group = bin, fill = groups, colour = groups), width = 1, size = 1) + 
  scale_y_reverse(sec.axis = sec_axis(~./800)) + scale_x_reverse(limits = c(686, 305)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  geom_hline(yintercept = c(0,45,135,225,270,315)) +
  geom_step(data = xdd, aes(x = Age, y = (nodules)/(framboids+nodules+undif)*800), 
            group = 1, colour = "blue") +
  geom_step(data = xdd, aes(x = Age, y = (nodules+framboids)/(framboids+nodules+undif)*800), 
            group = 1, colour = "red") +
  geom_vline(xintercept = c(Periods)) +  geom_vline(xintercept = c(HEATT), colour = "green") +
  geom_vline(xintercept = c(zones), colour = "orange") +
  scale_fill_manual(values=c("#00ba38","#00bfc4", "#619cff")) + 
  scale_colour_manual(values=c("#00ba38","#00bfc4", "#619cff"))

c <- ggplot(poles.sub2, aes(age, theta)) + 
  geom_boxplot(aes(group = bin, fill = groups, colour = groups), width = 1, size = 1) + 
  scale_y_reverse(sec.axis = sec_axis(~./800)) + scale_x_reverse(limits = c(315, -1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  geom_hline(yintercept = c(0,45,135,225,270,315)) +
  geom_step(data = xdd, aes(x = Age, y = (nodules)/(framboids+nodules+undif)*800), 
            group = 1, colour = "blue") +
  geom_step(data = xdd, aes(x = Age, y = (nodules+framboids)/(framboids+nodules+undif)*800), 
            group = 1, colour = "red") +
  geom_vline(xintercept = c(Periods)) +  geom_vline(xintercept = c(HEATT), colour = "green") +
  geom_vline(xintercept = c(zones), colour = "orange") +
  scale_fill_manual(values=c("#f8766d", "#b79f00","#00ba38")) + 
  scale_colour_manual(values=c("#f8766d", "#b79f00","#00ba38"))

grid.arrange(a,b,c,ncol=1)

##### Figure 6A - loess predictions #####

# age weighting - see Mehra et al. 2021. Curation and Analysis of Global Sedimentary Geochemical Data to Inform Earth History. GSA Today, 31. https://doi.org/10.1130/GSATG484A.1

age_bin <- 100 # span_age parameter

t_t <- as.data.frame(outer(poles$age, poles$age, `-`))/age_bin

t_t <- ((t_t^2)^0.5) # positive distance matrix

test <- 1/((t_t^2)+1)

test <- as.data.frame(test)

prox_t <- rowSums(test)

# generate two smooths, Precambrian and Phanerozoic

boot.R <- 10 # select number to bootstrap - final version boot n = 10000

# include all data (including outliers)
poles2 <- poles
prox_t2 <- prox_t

# Cross Validated (CV) loess span

# predict ideal spans (robust)

poles2.sub1 <- subset(poles2, age <= 540)
poles2.sub2 <- subset(poles2, age > 540)
prox_t22 <- cbind(poles2, prox_t2)
prox_t2.sub1 <- subset(prox_t22, age <= 540)
prox_t2.sub1 <- prox_t2.sub1[,10]
prox_t2.sub2 <- subset(prox_t22, age > 540)
prox_t2.sub2 <- prox_t2.sub2[,10]

loess.predict <- loess.as(poles2.sub1$age, poles2.sub1$theta, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, weights = 1/(prox_t2.sub1)) # , weights = 1/(prox_t2.sub1)
loess.predict[["pars"]][["span"]] # optimum span = 0.05 for Phanerozoic

loess.predict <- loess.as(poles2.sub2$age, poles2.sub2$theta, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, weights = 1/(prox_t2.sub2)) # , weights = 1/(prox_t2.sub2)
loess.predict[["pars"]][["span"]] # optimum span = 0.07 for Precambrian

# Phanerozoic loess boot

boot_fn <- function(poles2, indices) {
  d <- poles2[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(theta ~ age, d, span = 0.06, degree = 1, weights = 1/(prox_t2)) # remove weights here if replicating Fig. S19B 
  predict(loess_fit, data.frame(age = seq(0,539, 1)), se = T)$fit
}

loess_boot <- boot(poles2, R = boot.R, statistic = boot_fn) # for Fig. S19B, add the following: , weights = 1/((prox_t2*median(2/prox_t2))+1) # ca. 1 in 5

# 95% confidence intervals and median smooth
conf_97.5_A <- apply(loess_boot$t, 2, function(x) quantile(x, .975, na.rm = TRUE))
conf_50_A <- apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE))
conf_2.5_A <- apply(loess_boot$t, 2, function(x) quantile(x, .025, na.rm = TRUE))

conf_A <- as.data.frame(cbind(conf_2.5_A, conf_50_A, conf_97.5_A))
conf_A$age <- seq(0,539, 1)

outA <- as.data.frame(loess_boot$t)
names(outA)[1:ncol(outA)] <- seq(0,539, 1)
outA$group <- seq(1,(nrow(outA)), 1)

outA <- reshape2::melt(outA, id.vars = "group")
names(outA)[2:3] <- c("age", "theta")

outA$age2 <- rep(seq(0,539, 1),each = boot.R)

# Precambrian loess boot

boot_fn <- function(poles2, indices) {
  d <- poles2[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(theta ~ age, d, span = 0.07, degree = 1, weights = 1/(prox_t2)) # remove weights here if replicating Fig. S19B
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

loess_boot <- boot(poles2, R = boot.R, statistic = boot_fn) # for Fig. S19B, add the following: , weights = 1/((prox_t2*median(2/prox_t2))+1) # ca. 1 in 5

# 95% confidence intervals and median smooth
conf_97.5_B <- apply(loess_boot$t, 2, function(x) quantile(x, .975, na.rm = TRUE))
conf_50_B <- apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE))
conf_2.5_B <- apply(loess_boot$t, 2, function(x) quantile(x, .025, na.rm = TRUE))

conf_B <- as.data.frame(cbind(conf_2.5_B, conf_50_B, conf_97.5_B))
conf_B$age <- seq(540,4000, 10)

outB <- as.data.frame(loess_boot$t)
names(outB)[1:ncol(outB)] <- seq(540,4000, 10)
outB$group <- seq(1,(nrow(outB)), 1)

outB <- reshape2::melt(outB, id.vars = "group")
names(outB)[2:3] <- c("age", "theta")

outB$age2 <- rep(seq(540,4000, 10),each = boot.R)

# print - Fig. 5B

Periods <- c(2500, 1600, 1000, 720, 635, 541, 485.4, 443.8, 419.2, 358.9, 298.9, 
             251.9, 201.3, 145, 66, 23.03, 2.58) # add Period boundaries (ICS)

a <- ggplot() + 
  scale_x_reverse(limits = c(3500, 540)) +
  scale_y_reverse(limits = c(360,-30)) +
  stat_density_2d(data = poles2, aes(age, theta, fill = (..density..)^0.5), 
                  geom = "raster", contour = FALSE, adjust = 1/2) + 
  theme_bw() +
  geom_ribbon(data = conf_A, aes(x = age, ymin = conf_2.5_A, ymax = conf_97.5_A), colour = "red", fill = NA) +
  geom_ribbon(data = conf_B, aes(x = age, ymin = conf_2.5_B, ymax = conf_97.5_B), colour = "red", fill = NA) +
  geom_line(data = conf_A, aes(age, conf_50_A), colour = "red", size = 1) +
  geom_line(data = conf_B, aes(age, conf_50_B), colour = "red", size = 1) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  scale_fill_gradientn(colours = colorspace::sequential_hcl(6, palette = "Oslo", rev = FALSE)) +
  geom_hline(yintercept = c(0,45,135,225,270,315)) +
  geom_vline(xintercept = c(Periods))

b <- ggplot() + 
  scale_x_reverse(limits = c(541, -1)) +
  scale_y_reverse(limits = c(360,-30)) +
  stat_density_2d(data = poles2, aes(age, theta, fill = (..density..)^0.5), 
                  geom = "raster", contour = FALSE, adjust = 1/2) + 
  theme_bw() +
  geom_ribbon(data = conf_A, aes(x = age, ymin = conf_2.5_A, ymax = conf_97.5_A), colour = "red", fill = NA) +
  geom_ribbon(data = conf_B, aes(x = age, ymin = conf_2.5_B, ymax = conf_97.5_B), colour = "red", fill = NA) +
  geom_line(data = conf_A, aes(age, conf_50_A), colour = "red", size = 1) +
  geom_line(data = conf_B, aes(age, conf_50_B), colour = "red", size = 1) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  scale_fill_gradientn(colours = colorspace::sequential_hcl(6, palette = "Oslo", rev = FALSE)) +
  geom_hline(yintercept = c(0,45,135,225,270,315)) +
  geom_vline(xintercept = c(Periods)) + 
  geom_vline(xintercept = c(zones), colour = "red")

grid.arrange(a,b, ncol = 2)

##### Figure 6B - Interquartile range (IQR) loess #####

# compute IQR by age
theta.IQR <- poles %>%
  group_by(age, groups) %>%
  summarise(IQR = IQR(theta, na.rm = TRUE), age = age) # collapse age or not?

age_bin <- 100 # span_age parameter 

t_t <- as.data.frame(outer(theta.IQR$age, theta.IQR$age, `-`))/age_bin

t_t <- ((t_t^2)^0.5) # positive distance matrix

test <- 1/((t_t^2)+1)

test <- as.data.frame(test)

prox_t <- rowSums(test)

poles2 <- theta.IQR
prox_t2 <- prox_t

# Cross Validated (CV) loess span selection
poles2.sub1 <- subset(poles2, age <= 540)
poles2.sub2 <- subset(poles2, age > 540)

prox_t22 <- cbind(poles2, prox_t2)
prox_t2.sub1 <- subset(prox_t22, age <= 540)
prox_t2.sub1 <- prox_t2.sub1$...4
prox_t2.sub2 <- subset(prox_t22, age > 540)
prox_t2.sub2 <- prox_t2.sub2$...4

loess.predict <- loess.as(poles2.sub1$age, poles2.sub1$IQR, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, weights = 1/(prox_t2.sub1)) # 
loess.predict[["pars"]][["span"]] # 0.05 is Phanerozoic optimum

loess.predict <- loess.as(poles2.sub2$age, poles2.sub2$IQR, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, weights = 1/(prox_t2.sub2)) # 
loess.predict[["pars"]][["span"]] # 0.11 is Precambrian optimum 

boot.R <- 10 # select number to bootstrap n = 10000 for finalised

# Phanerozoic boot

boot_fn <- function(poles2, indices) {
  d <- poles2[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(IQR ~ age, d, span = 0.05, degree = 1, weights = 1/(prox_t2)) 
  predict(loess_fit, data.frame(age = seq(0,539, 1)), se = T)$fit
}

loess_boot <- boot(poles2, R = boot.R, statistic = boot_fn) 

# confidence intervals and median smooth
conf_97.5_A <- apply(loess_boot$t, 2, function(x) quantile(x, .975, na.rm = TRUE))
conf_50_A <- apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE))
conf_2.5_A <- apply(loess_boot$t, 2, function(x) quantile(x, .025, na.rm = TRUE))

conf_A <- as.data.frame(cbind(conf_2.5_A, conf_50_A, conf_97.5_A))
conf_A$age <- seq(0,539, 1)

outA <- as.data.frame(loess_boot$t)
names(outA)[1:ncol(outA)] <- seq(0,539, 1)
outA$group <- seq(1,(nrow(outA)), 1)

outA <- reshape2::melt(outA, id.vars = "group")
names(outA)[2:3] <- c("age", "IQR")

outA$age2 <- rep(seq(0,539, 1),each = boot.R)

# Precambrian boot

boot_fn <- function(poles2, indices) {
  d <- poles2[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(IQR ~ age, d, span = 0.11, degree = 1, weights = 1/(prox_t2)) 
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

loess_boot <- boot(poles2, R = boot.R, statistic = boot_fn) 

# confidence intervals and median smooth
conf_97.5_B <- apply(loess_boot$t, 2, function(x) quantile(x, .975, na.rm = TRUE))
conf_50_B <- apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE))
conf_2.5_B <- apply(loess_boot$t, 2, function(x) quantile(x, .025, na.rm = TRUE))

conf_B <- as.data.frame(cbind(conf_2.5_B, conf_50_B, conf_97.5_B))
conf_B$age <- seq(540,4000, 10)

outB <- as.data.frame(loess_boot$t)
names(outB)[1:ncol(outB)] <- seq(540,4000, 10)
outB$group <- seq(1,(nrow(outB)), 1)

outB <- reshape2::melt(outB, id.vars = "group")
names(outB)[2:3] <- c("age", "IQR")

outB$age2 <- rep(seq(540,4000, 10),each = boot.R)

# print - Fig. 5C

a <- ggplot() + 
  scale_x_reverse(limits = c(3500, 540)) +
  scale_y_continuous(limits = c(0,360)) +
  theme_bw() +
  geom_ribbon(data = conf_A, aes(x = age, ymin = conf_2.5_A, ymax = conf_97.5_A), colour = "red", fill = NA) +
  geom_ribbon(data = conf_B, aes(x = age, ymin = conf_2.5_B, ymax = conf_97.5_B), colour = "red", fill = NA) +
  geom_line(data = conf_A, aes(age, conf_50_A), colour = "red", size = 1) +
  geom_line(data = conf_B, aes(age, conf_50_B), colour = "red", size = 1) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  scale_fill_gradientn(colours = colorspace::sequential_hcl(6, palette = "Oslo", rev = FALSE)) +
  geom_hline(yintercept = median(poles2$IQR)) +
  geom_vline(xintercept = c(Periods))

b <- ggplot() + 
  scale_x_reverse(limits = c(541, -1)) +
  scale_y_continuous(limits = c(0,360)) +
  theme_bw() +
  geom_ribbon(data = conf_A, aes(x = age, ymin = conf_2.5_A, ymax = conf_97.5_A), colour = "red", fill = NA) +
  geom_ribbon(data = conf_B, aes(x = age, ymin = conf_2.5_B, ymax = conf_97.5_B), colour = "red", fill = NA) +
  geom_line(data = conf_A, aes(age, conf_50_A), colour = "red", size = 1) +
  geom_line(data = conf_B, aes(age, conf_50_B), colour = "red", size = 1) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  scale_fill_gradientn(colours = colorspace::sequential_hcl(6, palette = "Oslo", rev = FALSE)) +
  geom_hline(yintercept = median(poles2$IQR)) +
  geom_vline(xintercept = c(Periods)) + geom_vline(xintercept = c(zones), colour = "red")

grid.arrange(a,b, ncol = 2)

##### Fig. S18 - approximation of PCA biplot using element ratios #####

pyrites$ratio1 <- (pyrites$Mo+(pyrites$Cd*100)+(pyrites$Se*20))/(pyrites$Bi+pyrites$Te+pyrites$Au)
pyrites$ratio2 <- (pyrites$Sb+(pyrites$Ag*100))/(pyrites$Mn/100)

pyrites2 <- pyrites %>%
  group_by(clusters) %>%
  summarise(median1 = median(ratio1, na.rm = TRUE), IQR10.1 = quantile(ratio1, probs = 0.10, na.rm = TRUE), IQR90.1 = quantile(ratio1, probs = 0.90, na.rm = TRUE),
            median2 = median(ratio2, na.rm = TRUE), IQR10.2 = quantile(ratio2, probs = 0.10, na.rm = TRUE), IQR90.2 = quantile(ratio2, probs = 0.90, na.rm = TRUE))

pyrites2 <- subset(pyrites2, clusters != "Cluster High-Au") 

ggplot(pyrites2, aes(median1, median2)) + 
  stat_ellipse(data = subset(pyrites3, Location == "Cariaco basin"), 
               aes((Mo+(Cd*100)+(Se*20))/(Bi+Te+Au), (Sb+(Ag*100))/(Mn/100)),type = "t", level = 0.8, colour = "black") +
  geom_point(aes(fill = clusters), pch = 21) + theme_bw() +
  stat_ellipse(data = pyrites3, aes((Mo+(Cd*100)+(Se*20))/(Bi+Te+Au), (Sb+(Ag*100))/(Mn/100), 
                                      colour = clusters), type = "t", level = 0.8) +
  geom_errorbar(aes(ymin = IQR10.2, ymax = IQR90.2, colour = clusters), width = 0) +
  geom_errorbarh(aes(xmin = IQR10.1, xmax = IQR90.1, colour = clusters), height = 0) +
  scale_y_log10() + scale_x_log10()
