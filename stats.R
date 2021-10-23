##### This R script was written by Joe Emmings (British Geological Survey) ####
##### This R script implements data manipulation and statistical analysis #####
##### as in Figures 2-6A #####
##### please run 'gplates.R' to reproduce Figures 6B-D #####

# Each code block (for each figure) should be run consecutively #

##### KEY REFERENCES ####

# Figures 2-6 utilise pyrite trace element data accessible here: https://doi.org/10.1130/GEOL.S.12456332
# The data derive from:
# Large et al., Earth and Planetary Science Letters 389, 209-220 (2014).
# Large et al., Gondwana Research 28, 1282-1293 (2015).
# Large et al., Earth and Planetary Science Letters 428, 139-150 (2015).
# Large et al., Mineralium Deposita 54, 485-506 (2019).
# Mukherjee and Large, Geology 48, 1018-1022 (2020).

# Figures 5-6 utilise Phase 1 data from the Sedimentary Geochemistry and Palaeoenvironments Project (SGP) http://sgp-search.io/
# The citable reference for the Phase 1 data product is:
# Farrell et al., The Sedimentary Geochemistry and Palaeoenvironments Project. 788 Geobiology 00, 1-12 (2021).

# This code also sources the TOC/P data reported in:
# Reinhard et al., Evolution of the global phosphorus cycle. Nature 541, 386-389 (2017).

###### Packages #####

library(ggplot2)
library(reshape2)
library(ggbiplot)
library(BBmisc)
library(gridExtra)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(compositions)
library(data.table)
library(boot)
library(fANCOVA)
library(colorspace)
library(jsonlite)
library(httr)
library(rio)
library(msir)
library(caret)
library(zCompositions)
library(pammtools)
library(randomForest)
library(deeptime)
library(ggfittext)
library(stringr)
library(kernlab)

# set working directory

setwd(...)

##### Figure 1 - text mining outputs ####

# install.R, preparation.R, analysis.R

##### Figure 2A - Hierarchical cluster analysis (HCA) ####

# download the pyrite trace element dataset of Mukherjee and Large (2020) (CC-BY-NC)
# https://doi.org/10.1130/GEOL.S.12456332 

# import

source_pyrites <- 'https://gsapubs.figshare.com/ndownloader/files/23030852'
pyrites <- rio::import(source_pyrites)

# clean 
colnames(pyrites)[1:37] <- pyrites[1,1:37]
colnames(pyrites)[c(1,3)] <- c("Sample", "Age.Ma")
pyrites <- pyrites[-c(1,ncol(pyrites),ncol(pyrites)-1),-c(12, 24:37)]
names(pyrites) <- gsub(x = names(pyrites), pattern = "_Py", replacement = "")
pyrites <- pyrites[-c(nrow(pyrites),nrow(pyrites)-1),]
pyrites <- data.frame(lapply(pyrites, function(x) {
  gsub("<", "-", x)
  }))
pyrites[,3:22] <- sapply(pyrites[,3:22], as.numeric) # coerce to numeric
pyrites <- cbind(rownames(pyrites),pyrites)
names(pyrites)[1] <- "row"
pyrites <- reshape2::melt(pyrites, id.vars = c("row","Sample", "Location", "Age.Ma"))
pyrites1 <- subset(pyrites, value < 0)
pyrites1$value <- abs(pyrites1$value*(2/3)) # 2/3 value (sensu Boogaart and Tolosana-Delgado)
pyrites2 <- subset(pyrites, value >= 0)
pyrites <- rbind(pyrites1,pyrites2)
pyrites <- reshape2::dcast(pyrites, row+Sample+Location+Age.Ma~variable, value.var = "value")
pyrites$Pb <- pyrites$X206Pb+pyrites$X207Pb+pyrites$X208Pb
pyrites.all <- pyrites[,-c(1,20:22)]

# cull to focal area

sub <- grepl("Russia", pyrites.all$Location, ignore.case=TRUE) | 
  grepl("China", pyrites.all$Location, ignore.case=TRUE) |
  grepl("Ukraine", pyrites.all$Location, ignore.case=TRUE) |
  grepl("Germany", pyrites.all$Location, ignore.case=TRUE) |
  grepl("Turkey", pyrites.all$Location, ignore.case=TRUE) |
  grepl("Malaysia", pyrites.all$Location, ignore.case=TRUE) |
  grepl("Botswana", pyrites.all$Location, ignore.case=TRUE) |
  grepl("Africa", pyrites.all$Location, ignore.case=TRUE) |
  grepl("India", pyrites.all$Location, ignore.case=TRUE) |
  grepl("Sweden", pyrites.all$Location, ignore.case=TRUE) |
  grepl("South Afr", pyrites.all$Location, ignore.case=TRUE) |
  grepl("Zambia", pyrites.all$Location, ignore.case=TRUE) |
  grepl("Finland", pyrites.all$Location, ignore.case=TRUE)

pyrites.cull <- subset(pyrites.all, !sub)

# options
#pyrites <- pyrites.cull # culled dataset (not deployed)
pyrites <- pyrites.all # all samples

# explore missing (NA) values
sapply(pyrites, function(y) sum(length(which(is.na(y)))))

# assume any remaining 0 values (limited to Au) are not true zero but missing data (NA values)

pyrites[pyrites == 0] <- NA 

# impute missing values
imp <- as.matrix(pyrites[,c(4:20)])
imp <- multRepl(imp, imp.missing = TRUE, label = NA, closure=10^6)

pyrites <- cbind(pyrites[,1:3], imp)

# replace any NA values in Sample and Location fields
pyrites$Sample[is.na(pyrites$Sample)] <- "Unknown"
pyrites$Location[is.na(pyrites$Location)] <- "Unknown"

#  remove any remaining rows with NA (i.e., where age = NA)
pyrites <- pyrites[complete.cases(pyrites[,]),]

HCA <- pyrites[4:(ncol(pyrites))]

# clr trans
HCA2 <- clr(HCA)  # Note this step plus embedded dist in hclust below is equivalent to hclust of distance matrix of acomp class
HCA2 <- as.data.frame(HCA2) 

# experimental

HCA3 <- variation(acomp(HCA)) # sensu. Boogaart and Tolosana-Delgado

# cluster
hr <- hclust(dist(HCA2, method = "euclidean"), method = "ward.D2")
hc <- hclust(dist(HCA3, method = "euclidean"), method = "ward.D2") # ward.D2

plot(hr) # display dendrogram
rect.hclust(hr, k=5, border="red") 
clusters <- cutree(hr, k=5)

# meaningful cluster names

clusters[clusters == 1] <- "Type.4" 
clusters[clusters == 2] <- "Type.1"
clusters[clusters == 3] <- "Type.5" 
clusters[clusters == 4] <- "Type.2" 
clusters[clusters == 5] <- "Type.3" 

plot(hc) # display dendrogram
#rect.hclust(hc, k=4, border="red") 

pyrites$clusters <- clusters

col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))

# output
labels <- pyrites[4:21]
my_group <- as.numeric(as.factor(substr(labels$cluster, 1 , 2))) # select zone or cluster
colSide <- brewer.pal(6, "Set3")[my_group] # select 10 or 6

# Fig. 2A

HCA4 <- BBmisc::normalize(HCA2, method = "standardize") # scale and centre for HCA matrix visualization

heatmap(as.matrix(HCA4), Rowv = as.dendrogram(hr), Colv = as.dendrogram(hc), 
        col = col, labRow = labels$cluster)

# for colour scale
pheatmap(as.matrix(HCA4), cluster_rows = hr, 
         cluster_cols = hc)

# tidy using graphics package

##### Figure 2B - Principal Component Analysis (PCA) ####

pyrites$totals <- rowSums(pyrites[,c("Mn","Co","Ni","Cu","Zn","As","Se","Mo","Ag","Cd","Sb","Te","Au","Tl","Pb","Bi","Pt" )], na.rm = TRUE)
pyrites$clusters <- factor(pyrites$clusters, levels = c("Type.5", "Type.2", "Type.1", "Type.3", "Type.4"))

# option to select variables
#myvars <- c("Co", "Se", "Zn", "Cd", "Mo", "Tl", "Bi", "Mn") # selected elements
#pca.dat <- pyrites[myvars]

pca.dat <- pyrites[c(4:20)] # select all

pca.dat <- cbind(pyrites[,1:3],pca.dat, pyrites[,21])

# clr transform
pca.dat[,4:(ncol(pca.dat)-1)] <- clr(pca.dat[,4:(ncol(pca.dat)-1)]) 

labels <- pca.dat[,2]
age <- pca.dat[,3]
HCA.clusters <- pca.dat[,ncol(pca.dat)]
HCA.clusters2 <- as.character(HCA.clusters)

pca.dat <- pca.dat[,4:(ncol(pca.dat)-1)]

# compute principal components
PCAdata <- prcomp(pca.dat,
                  center = TRUE,
                  scale. = FALSE) # covariance matrix

summary(PCAdata)
plot(PCAdata, type = "l") # scree plot shows elbow after PC3

# print PCA biplot (Fig 2B)

g <- ggbiplot(PCAdata, obs.scale = 1, var.scale = 1, 
              ellipse = TRUE, ellipse.prob = 0.68,
              circle = FALSE, choices = c(1,2), groups = HCA.clusters2) # select clusters
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top') + theme_bw()
print(g) 

# tidy using graphics package

##### Figure 2C - total concentration boxplots #####

ggplot(pyrites, aes(clusters, totals)) + 
  geom_boxplot(aes(group = clusters), notch = TRUE, outlier.shape = NA, coef = 0, fill= 'grey') + theme_bw() +
  stat_summary(fun="median", geom="point", aes(group = clusters), size = 3) +
  coord_cartesian(ylim=c(0,12000))

# alternative violin plots

ggplot(pyrites, aes(clusters, totals)) + 
  geom_violin(aes(group = clusters), draw_quantiles = c(0.25,0.75), trim = TRUE, fill = 'grey') + theme_bw() +
  stat_summary(fun="median", geom="point", aes(group = clusters), size = 3) +
  stat_summary(fun="median", geom="line", group = 1) +
  scale_y_log10()

##### Figure 4 #####
##### extract PC1-7 scores ####

# extract cartesian coords

coords <- as.data.frame(PCAdata[["x"]])
coords.all <- coords[,1:7] # includes PCs 1-7

coords.all <- cbind(coords.all, labels, HCA.clusters, age)

###### prediction interval bootstraps up to 2sd #####

## in order to adopt a fixed span (Fig. S21-S22)
# the user should replace 'span_upper <- loess.predict[["pars"]][["span"]]' and 'span_lower <- loess.predict[["pars"]][["span"]]'
# with 'span_upper <- loess.predict[["pars"]][["span"]]' and 'span_lower <- loess.predict[["pars"]][["span"]]'

coords.all2 <- coords.all

# may want to set.seed for reproducibility

speed <- 10 # Use speed = 10 for full version. Use speed = 100 for quicker preliminary results

# construct normal distribution for sampling

n <- 10000
norm <- rnorm(n)
out <- hist(norm, breaks = c(seq(-4.3,4.3, by = 0.2)))
bins <- as.data.frame(out[[1]])
bins <- round(bins+0.1, digits = 2)
freqs <- as.data.frame(out[[2]])
out <- as.data.frame(cbind(bins[1:nrow(bins)-1,],freqs))
names(out)[1:2] <- c("sd","count")

##### Phanerozoic PC1 #####

coords.all_Phan <- subset(coords.all2, age < 541)

# 2 sigma

nsig <- 2

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_2sd_phan_upper <- as.data.frame(cbind(PC1_2sd_phan$x, PC1_2sd_phan$upper))
names(PC1_2sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_2sd_phan_upper$age, PC1_2sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_2sd_phan_upper, indices) {
  d <- PC1_2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_2sd_phan_upper_50)[1] <- "PC1"

PC1_2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_2sd_phan_upper_all)[1:ncol(PC1_2sd_phan_upper_all)] <- seq(0,540, 1)
PC1_2sd_phan_upper_all$group <- seq(1,(nrow(PC1_2sd_phan_upper_all)), 1)

PC1_2sd_phan_upper_all <- reshape2::melt(PC1_2sd_phan_upper_all, id.vars = "group")
names(PC1_2sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_2sd_phan_lower <- as.data.frame(cbind(PC1_2sd_phan$x, PC1_2sd_phan$lower))
names(PC1_2sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_2sd_phan_lower$age, PC1_2sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_2sd_phan_lower, indices) {
  d <- PC1_2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_2sd_phan_lower_50)[1] <- "PC1"

PC1_2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_2sd_phan_lower_all)[1:ncol(PC1_2sd_phan_lower_all)] <- seq(0,540, 1)
PC1_2sd_phan_lower_all$group <- seq(1,(nrow(PC1_2sd_phan_lower_all)), 1)

PC1_2sd_phan_lower_all <- reshape2::melt(PC1_2sd_phan_lower_all, id.vars = "group")
names(PC1_2sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.8sd_phan_upper <- as.data.frame(cbind(PC1_1.8sd_phan$x, PC1_1.8sd_phan$upper))
names(PC1_1.8sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.8sd_phan_upper$age, PC1_1.8sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_1.8sd_phan_upper, indices) {
  d <- PC1_1.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_1.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_1.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_1.8sd_phan_upper_50)[1] <- "PC1"

PC1_1.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_1.8sd_phan_upper_all)[1:ncol(PC1_1.8sd_phan_upper_all)] <- seq(0,540, 1)
PC1_1.8sd_phan_upper_all$group <- seq(1,(nrow(PC1_1.8sd_phan_upper_all)), 1)

PC1_1.8sd_phan_upper_all <- reshape2::melt(PC1_1.8sd_phan_upper_all, id.vars = "group")
names(PC1_1.8sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_1.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.8sd_phan_lower <- as.data.frame(cbind(PC1_1.8sd_phan$x, PC1_1.8sd_phan$lower))
names(PC1_1.8sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.8sd_phan_lower$age, PC1_1.8sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_1.8sd_phan_lower, indices) {
  d <- PC1_1.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_1.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_1.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_1.8sd_phan_lower_50)[1] <- "PC1"

PC1_1.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_1.8sd_phan_lower_all)[1:ncol(PC1_1.8sd_phan_lower_all)] <- seq(0,540, 1)
PC1_1.8sd_phan_lower_all$group <- seq(1,(nrow(PC1_1.8sd_phan_lower_all)), 1)

PC1_1.8sd_phan_lower_all <- reshape2::melt(PC1_1.8sd_phan_lower_all, id.vars = "group")
names(PC1_1.8sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_1.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.6sd_phan_upper <- as.data.frame(cbind(PC1_1.6sd_phan$x, PC1_1.6sd_phan$upper))
names(PC1_1.6sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.6sd_phan_upper$age, PC1_1.6sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_1.6sd_phan_upper, indices) {
  d <- PC1_1.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_1.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_1.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_1.6sd_phan_upper_50)[1] <- "PC1"

PC1_1.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_1.6sd_phan_upper_all)[1:ncol(PC1_1.6sd_phan_upper_all)] <- seq(0,540, 1)
PC1_1.6sd_phan_upper_all$group <- seq(1,(nrow(PC1_1.6sd_phan_upper_all)), 1)

PC1_1.6sd_phan_upper_all <- reshape2::melt(PC1_1.6sd_phan_upper_all, id.vars = "group")
names(PC1_1.6sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_1.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.6sd_phan_lower <- as.data.frame(cbind(PC1_1.6sd_phan$x, PC1_1.6sd_phan$lower))
names(PC1_1.6sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.6sd_phan_lower$age, PC1_1.6sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_1.6sd_phan_lower, indices) {
  d <- PC1_1.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_1.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_1.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_1.6sd_phan_lower_50)[1] <- "PC1"

PC1_1.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_1.6sd_phan_lower_all)[1:ncol(PC1_1.6sd_phan_lower_all)] <- seq(0,540, 1)
PC1_1.6sd_phan_lower_all$group <- seq(1,(nrow(PC1_1.6sd_phan_lower_all)), 1)

PC1_1.6sd_phan_lower_all <- reshape2::melt(PC1_1.6sd_phan_lower_all, id.vars = "group")
names(PC1_1.6sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_1.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.4sd_phan_upper <- as.data.frame(cbind(PC1_1.4sd_phan$x, PC1_1.4sd_phan$upper))
names(PC1_1.4sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.4sd_phan_upper$age, PC1_1.4sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_1.4sd_phan_upper, indices) {
  d <- PC1_1.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_1.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_1.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_1.4sd_phan_upper_50)[1] <- "PC1"

PC1_1.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_1.4sd_phan_upper_all)[1:ncol(PC1_1.4sd_phan_upper_all)] <- seq(0,540, 1)
PC1_1.4sd_phan_upper_all$group <- seq(1,(nrow(PC1_1.4sd_phan_upper_all)), 1)

PC1_1.4sd_phan_upper_all <- reshape2::melt(PC1_1.4sd_phan_upper_all, id.vars = "group")
names(PC1_1.4sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_1.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.4sd_phan_lower <- as.data.frame(cbind(PC1_1.4sd_phan$x, PC1_1.4sd_phan$lower))
names(PC1_1.4sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.4sd_phan_lower$age, PC1_1.4sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_1.4sd_phan_lower, indices) {
  d <- PC1_1.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_1.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_1.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_1.4sd_phan_lower_50)[1] <- "PC1"

PC1_1.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_1.4sd_phan_lower_all)[1:ncol(PC1_1.4sd_phan_lower_all)] <- seq(0,540, 1)
PC1_1.4sd_phan_lower_all$group <- seq(1,(nrow(PC1_1.4sd_phan_lower_all)), 1)

PC1_1.4sd_phan_lower_all <- reshape2::melt(PC1_1.4sd_phan_lower_all, id.vars = "group")
names(PC1_1.4sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_1.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.2sd_phan_upper <- as.data.frame(cbind(PC1_1.2sd_phan$x, PC1_1.2sd_phan$upper))
names(PC1_1.2sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.2sd_phan_upper$age, PC1_1.2sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_1.2sd_phan_upper, indices) {
  d <- PC1_1.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_1.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_1.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_1.2sd_phan_upper_50)[1] <- "PC1"

PC1_1.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_1.2sd_phan_upper_all)[1:ncol(PC1_1.2sd_phan_upper_all)] <- seq(0,540, 1)
PC1_1.2sd_phan_upper_all$group <- seq(1,(nrow(PC1_1.2sd_phan_upper_all)), 1)

PC1_1.2sd_phan_upper_all <- reshape2::melt(PC1_1.2sd_phan_upper_all, id.vars = "group")
names(PC1_1.2sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_1.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.2sd_phan_lower <- as.data.frame(cbind(PC1_1.2sd_phan$x, PC1_1.2sd_phan$lower))
names(PC1_1.2sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.2sd_phan_lower$age, PC1_1.2sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_1.2sd_phan_lower, indices) {
  d <- PC1_1.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_1.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_1.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_1.2sd_phan_lower_50)[1] <- "PC1"

PC1_1.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_1.2sd_phan_lower_all)[1:ncol(PC1_1.2sd_phan_lower_all)] <- seq(0,540, 1)
PC1_1.2sd_phan_lower_all$group <- seq(1,(nrow(PC1_1.2sd_phan_lower_all)), 1)

PC1_1.2sd_phan_lower_all <- reshape2::melt(PC1_1.2sd_phan_lower_all, id.vars = "group")
names(PC1_1.2sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_1.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.0sd_phan_upper <- as.data.frame(cbind(PC1_1.0sd_phan$x, PC1_1.0sd_phan$upper))
names(PC1_1.0sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.0sd_phan_upper$age, PC1_1.0sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_1.0sd_phan_upper, indices) {
  d <- PC1_1.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_1.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_1.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_1.0sd_phan_upper_50)[1] <- "PC1"

PC1_1.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_1.0sd_phan_upper_all)[1:ncol(PC1_1.0sd_phan_upper_all)] <- seq(0,540, 1)
PC1_1.0sd_phan_upper_all$group <- seq(1,(nrow(PC1_1.0sd_phan_upper_all)), 1)

PC1_1.0sd_phan_upper_all <- reshape2::melt(PC1_1.0sd_phan_upper_all, id.vars = "group")
names(PC1_1.0sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_1.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.0sd_phan_lower <- as.data.frame(cbind(PC1_1.0sd_phan$x, PC1_1.0sd_phan$lower))
names(PC1_1.0sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.0sd_phan_lower$age, PC1_1.0sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_1.0sd_phan_lower, indices) {
  d <- PC1_1.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_1.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_1.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_1.0sd_phan_lower_50)[1] <- "PC1"

PC1_1.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_1.0sd_phan_lower_all)[1:ncol(PC1_1.0sd_phan_lower_all)] <- seq(0,540, 1)
PC1_1.0sd_phan_lower_all$group <- seq(1,(nrow(PC1_1.0sd_phan_lower_all)), 1)

PC1_1.0sd_phan_lower_all <- reshape2::melt(PC1_1.0sd_phan_lower_all, id.vars = "group")
names(PC1_1.0sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_1.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.8sd_phan_upper <- as.data.frame(cbind(PC1_0.8sd_phan$x, PC1_0.8sd_phan$upper))
names(PC1_0.8sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.8sd_phan_upper$age, PC1_0.8sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_0.8sd_phan_upper, indices) {
  d <- PC1_0.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_0.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_0.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_0.8sd_phan_upper_50)[1] <- "PC1"

PC1_0.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_0.8sd_phan_upper_all)[1:ncol(PC1_0.8sd_phan_upper_all)] <- seq(0,540, 1)
PC1_0.8sd_phan_upper_all$group <- seq(1,(nrow(PC1_0.8sd_phan_upper_all)), 1)

PC1_0.8sd_phan_upper_all <- reshape2::melt(PC1_0.8sd_phan_upper_all, id.vars = "group")
names(PC1_0.8sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_0.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.8sd_phan_lower <- as.data.frame(cbind(PC1_0.8sd_phan$x, PC1_0.8sd_phan$lower))
names(PC1_0.8sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.8sd_phan_lower$age, PC1_0.8sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_0.8sd_phan_lower, indices) {
  d <- PC1_0.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_0.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_0.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_0.8sd_phan_lower_50)[1] <- "PC1"

PC1_0.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_0.8sd_phan_lower_all)[1:ncol(PC1_0.8sd_phan_lower_all)] <- seq(0,540, 1)
PC1_0.8sd_phan_lower_all$group <- seq(1,(nrow(PC1_0.8sd_phan_lower_all)), 1)

PC1_0.8sd_phan_lower_all <- reshape2::melt(PC1_0.8sd_phan_lower_all, id.vars = "group")
names(PC1_0.8sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_0.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.6sd_phan_upper <- as.data.frame(cbind(PC1_0.6sd_phan$x, PC1_0.6sd_phan$upper))
names(PC1_0.6sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.6sd_phan_upper$age, PC1_0.6sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_0.6sd_phan_upper, indices) {
  d <- PC1_0.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_0.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_0.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_0.6sd_phan_upper_50)[1] <- "PC1"

PC1_0.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_0.6sd_phan_upper_all)[1:ncol(PC1_0.6sd_phan_upper_all)] <- seq(0,540, 1)
PC1_0.6sd_phan_upper_all$group <- seq(1,(nrow(PC1_0.6sd_phan_upper_all)), 1)

PC1_0.6sd_phan_upper_all <- reshape2::melt(PC1_0.6sd_phan_upper_all, id.vars = "group")
names(PC1_0.6sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_0.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.6sd_phan_lower <- as.data.frame(cbind(PC1_0.6sd_phan$x, PC1_0.6sd_phan$lower))
names(PC1_0.6sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.6sd_phan_lower$age, PC1_0.6sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_0.6sd_phan_lower, indices) {
  d <- PC1_0.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_0.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_0.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_0.6sd_phan_lower_50)[1] <- "PC1"

PC1_0.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_0.6sd_phan_lower_all)[1:ncol(PC1_0.6sd_phan_lower_all)] <- seq(0,540, 1)
PC1_0.6sd_phan_lower_all$group <- seq(1,(nrow(PC1_0.6sd_phan_lower_all)), 1)

PC1_0.6sd_phan_lower_all <- reshape2::melt(PC1_0.6sd_phan_lower_all, id.vars = "group")
names(PC1_0.6sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_0.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.4sd_phan_upper <- as.data.frame(cbind(PC1_0.4sd_phan$x, PC1_0.4sd_phan$upper))
names(PC1_0.4sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.4sd_phan_upper$age, PC1_0.4sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_0.4sd_phan_upper, indices) {
  d <- PC1_0.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_0.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_0.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_0.4sd_phan_upper_50)[1] <- "PC1"

PC1_0.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_0.4sd_phan_upper_all)[1:ncol(PC1_0.4sd_phan_upper_all)] <- seq(0,540, 1)
PC1_0.4sd_phan_upper_all$group <- seq(1,(nrow(PC1_0.4sd_phan_upper_all)), 1)

PC1_0.4sd_phan_upper_all <- reshape2::melt(PC1_0.4sd_phan_upper_all, id.vars = "group")
names(PC1_0.4sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_0.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.4sd_phan_lower <- as.data.frame(cbind(PC1_0.4sd_phan$x, PC1_0.4sd_phan$lower))
names(PC1_0.4sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.4sd_phan_lower$age, PC1_0.4sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_0.4sd_phan_lower, indices) {
  d <- PC1_0.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_0.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_0.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_0.4sd_phan_lower_50)[1] <- "PC1"

PC1_0.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_0.4sd_phan_lower_all)[1:ncol(PC1_0.4sd_phan_lower_all)] <- seq(0,540, 1)
PC1_0.4sd_phan_lower_all$group <- seq(1,(nrow(PC1_0.4sd_phan_lower_all)), 1)

PC1_0.4sd_phan_lower_all <- reshape2::melt(PC1_0.4sd_phan_lower_all, id.vars = "group")
names(PC1_0.4sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_0.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.2sd_phan_upper <- as.data.frame(cbind(PC1_0.2sd_phan$x, PC1_0.2sd_phan$upper))
names(PC1_0.2sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.2sd_phan_upper$age, PC1_0.2sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_0.2sd_phan_upper, indices) {
  d <- PC1_0.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_0.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_0.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_0.2sd_phan_upper_50)[1] <- "PC1"

PC1_0.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_0.2sd_phan_upper_all)[1:ncol(PC1_0.2sd_phan_upper_all)] <- seq(0,540, 1)
PC1_0.2sd_phan_upper_all$group <- seq(1,(nrow(PC1_0.2sd_phan_upper_all)), 1)

PC1_0.2sd_phan_upper_all <- reshape2::melt(PC1_0.2sd_phan_upper_all, id.vars = "group")
names(PC1_0.2sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_0.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.2sd_phan_lower <- as.data.frame(cbind(PC1_0.2sd_phan$x, PC1_0.2sd_phan$lower))
names(PC1_0.2sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.2sd_phan_lower$age, PC1_0.2sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_0.2sd_phan_lower, indices) {
  d <- PC1_0.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_0.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_0.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_0.2sd_phan_lower_50)[1] <- "PC1"

PC1_0.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_0.2sd_phan_lower_all)[1:ncol(PC1_0.2sd_phan_lower_all)] <- seq(0,540, 1)
PC1_0.2sd_phan_lower_all$group <- seq(1,(nrow(PC1_0.2sd_phan_lower_all)), 1)

PC1_0.2sd_phan_lower_all <- reshape2::melt(PC1_0.2sd_phan_lower_all, id.vars = "group")
names(PC1_0.2sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_0.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.0sd_phan_upper <- as.data.frame(cbind(PC1_0.0sd_phan$x, PC1_0.0sd_phan$upper))
names(PC1_0.0sd_phan_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.0sd_phan_upper$age, PC1_0.0sd_phan_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_0.0sd_phan_upper, indices) {
  d <- PC1_0.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_0.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_0.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC1_0.0sd_phan_upper_50)[1] <- "PC1"

PC1_0.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC1_0.0sd_phan_upper_all)[1:ncol(PC1_0.0sd_phan_upper_all)] <- seq(0,540, 1)
PC1_0.0sd_phan_upper_all$group <- seq(1,(nrow(PC1_0.0sd_phan_upper_all)), 1)

PC1_0.0sd_phan_upper_all <- reshape2::melt(PC1_0.0sd_phan_upper_all, id.vars = "group")
names(PC1_0.0sd_phan_upper_all)[2:3] <- c("age", "PC1")

PC1_0.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.0sd_phan_lower <- as.data.frame(cbind(PC1_0.0sd_phan$x, PC1_0.0sd_phan$lower))
names(PC1_0.0sd_phan_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.0sd_phan_lower$age, PC1_0.0sd_phan_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_0.0sd_phan_lower, indices) {
  d <- PC1_0.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_0.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_0.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC1_0.0sd_phan_lower_50)[1] <- "PC1"

PC1_0.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC1_0.0sd_phan_lower_all)[1:ncol(PC1_0.0sd_phan_lower_all)] <- seq(0,540, 1)
PC1_0.0sd_phan_lower_all$group <- seq(1,(nrow(PC1_0.0sd_phan_lower_all)), 1)

PC1_0.0sd_phan_lower_all <- reshape2::melt(PC1_0.0sd_phan_lower_all, id.vars = "group")
names(PC1_0.0sd_phan_lower_all)[2:3] <- c("age", "PC1")

PC1_0.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

a <- ggplot() + 
  #geom_line(data = PC1_2sd_phan_upper_all, aes(age2,PC1, group = group)) +
  #geom_line(data = PC1_2sd_phan_lower_all, aes(age2,PC1, group = group)) +
  geom_line(data = PC1_2sd_phan_upper_50, aes(age,PC1), colour = 'blue') +
  geom_line(data = PC1_2sd_phan_lower_50, aes(age,PC1), colour = 'blue') +
  geom_line(data = PC1_1.6sd_phan_upper_50, aes(age,PC1), colour = 'green') +
  geom_line(data = PC1_1.6sd_phan_lower_50, aes(age,PC1), colour = 'green') +
  geom_line(data = PC1_1.4sd_phan_upper_50, aes(age,PC1), colour = 'orange') +
  geom_line(data = PC1_1.4sd_phan_lower_50, aes(age,PC1), colour = 'orange') +
  geom_line(data = PC1_1.2sd_phan_upper_50, aes(age,PC1), colour = 'pink') +
  geom_line(data = PC1_1.2sd_phan_lower_50, aes(age,PC1), colour = 'pink') +
  geom_line(data = PC1_1.0sd_phan_upper_50, aes(age,PC1), colour = 'brown') +
  geom_line(data = PC1_1.0sd_phan_lower_50, aes(age,PC1), colour = 'brown') +
  geom_line(data = PC1_0.8sd_phan_upper_50, aes(age,PC1), colour = 'yellow') +
  geom_line(data = PC1_0.8sd_phan_lower_50, aes(age,PC1), colour = 'yellow') +
  geom_line(data = PC1_0.6sd_phan_upper_50, aes(age,PC1), colour = 'red') +
  geom_line(data = PC1_0.6sd_phan_lower_50, aes(age,PC1), colour = 'red') +
  geom_line(data = PC1_0.4sd_phan_upper_50, aes(age,PC1), colour = 'blue') +
  geom_line(data = PC1_0.4sd_phan_lower_50, aes(age,PC1), colour = 'blue') +
  geom_line(data = PC1_0.2sd_phan_upper_50, aes(age,PC1), colour = 'green') +
  geom_line(data = PC1_0.2sd_phan_lower_50, aes(age,PC1), colour = 'green') +
  geom_line(data = PC1_0.0sd_phan_upper_50, aes(age,PC1), colour = 'orange') +
  geom_line(data = PC1_0.0sd_phan_lower_50, aes(age,PC1), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC1)) + 
  #geom_point(data = PC1_2sd_phan_upper, aes(age, PC1), colour = 'red') +
  scale_x_reverse(limits = c(540,0)) + theme_bw()

##### Precambrian PC1 #####

coords.all_pre <- subset(coords.all2, age > 541)

# 2 sigma

nsig <- 2

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_2sd_pre_upper <- as.data.frame(cbind(PC1_2sd_pre$x, PC1_2sd_pre$upper))
names(PC1_2sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_2sd_pre_upper$age, PC1_2sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[1], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_2sd_pre_upper, indices) {
  d <- PC1_2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_2sd_pre_upper_50)[1] <- "PC1"

PC1_2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_2sd_pre_upper_all)[1:ncol(PC1_2sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_2sd_pre_upper_all$group <- seq(1,(nrow(PC1_2sd_pre_upper_all)), 1)

PC1_2sd_pre_upper_all <- reshape2::melt(PC1_2sd_pre_upper_all, id.vars = "group")
names(PC1_2sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_2sd_pre_lower <- as.data.frame(cbind(PC1_2sd_pre$x, PC1_2sd_pre$lower))
names(PC1_2sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_2sd_pre_lower$age, PC1_2sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_2sd_pre_lower, indices) {
  d <- PC1_2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_2sd_pre_lower_50)[1] <- "PC1"

PC1_2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_2sd_pre_lower_all)[1:ncol(PC1_2sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_2sd_pre_lower_all$group <- seq(1,(nrow(PC1_2sd_pre_lower_all)), 1)

PC1_2sd_pre_lower_all <- reshape2::melt(PC1_2sd_pre_lower_all, id.vars = "group")
names(PC1_2sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.8sd_pre_upper <- as.data.frame(cbind(PC1_1.8sd_pre$x, PC1_1.8sd_pre$upper))
names(PC1_1.8sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.8sd_pre_upper$age, PC1_1.8sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_1.8sd_pre_upper, indices) {
  d <- PC1_1.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_1.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_1.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_1.8sd_pre_upper_50)[1] <- "PC1"

PC1_1.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_1.8sd_pre_upper_all)[1:ncol(PC1_1.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_1.8sd_pre_upper_all$group <- seq(1,(nrow(PC1_1.8sd_pre_upper_all)), 1)

PC1_1.8sd_pre_upper_all <- reshape2::melt(PC1_1.8sd_pre_upper_all, id.vars = "group")
names(PC1_1.8sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_1.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.8sd_pre_lower <- as.data.frame(cbind(PC1_1.8sd_pre$x, PC1_1.8sd_pre$lower))
names(PC1_1.8sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.8sd_pre_lower$age, PC1_1.8sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_1.8sd_pre_lower, indices) {
  d <- PC1_1.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_1.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_1.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_1.8sd_pre_lower_50)[1] <- "PC1"

PC1_1.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_1.8sd_pre_lower_all)[1:ncol(PC1_1.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_1.8sd_pre_lower_all$group <- seq(1,(nrow(PC1_1.8sd_pre_lower_all)), 1)

PC1_1.8sd_pre_lower_all <- reshape2::melt(PC1_1.8sd_pre_lower_all, id.vars = "group")
names(PC1_1.8sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_1.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.6sd_pre_upper <- as.data.frame(cbind(PC1_1.6sd_pre$x, PC1_1.6sd_pre$upper))
names(PC1_1.6sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.6sd_pre_upper$age, PC1_1.6sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_1.6sd_pre_upper, indices) {
  d <- PC1_1.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_1.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_1.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_1.6sd_pre_upper_50)[1] <- "PC1"

PC1_1.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_1.6sd_pre_upper_all)[1:ncol(PC1_1.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_1.6sd_pre_upper_all$group <- seq(1,(nrow(PC1_1.6sd_pre_upper_all)), 1)

PC1_1.6sd_pre_upper_all <- reshape2::melt(PC1_1.6sd_pre_upper_all, id.vars = "group")
names(PC1_1.6sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_1.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.6sd_pre_lower <- as.data.frame(cbind(PC1_1.6sd_pre$x, PC1_1.6sd_pre$lower))
names(PC1_1.6sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.6sd_pre_lower$age, PC1_1.6sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_1.6sd_pre_lower, indices) {
  d <- PC1_1.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_1.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_1.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_1.6sd_pre_lower_50)[1] <- "PC1"

PC1_1.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_1.6sd_pre_lower_all)[1:ncol(PC1_1.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_1.6sd_pre_lower_all$group <- seq(1,(nrow(PC1_1.6sd_pre_lower_all)), 1)

PC1_1.6sd_pre_lower_all <- reshape2::melt(PC1_1.6sd_pre_lower_all, id.vars = "group")
names(PC1_1.6sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_1.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.4sd_pre_upper <- as.data.frame(cbind(PC1_1.4sd_pre$x, PC1_1.4sd_pre$upper))
names(PC1_1.4sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.4sd_pre_upper$age, PC1_1.4sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_1.4sd_pre_upper, indices) {
  d <- PC1_1.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_1.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_1.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_1.4sd_pre_upper_50)[1] <- "PC1"

PC1_1.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_1.4sd_pre_upper_all)[1:ncol(PC1_1.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_1.4sd_pre_upper_all$group <- seq(1,(nrow(PC1_1.4sd_pre_upper_all)), 1)

PC1_1.4sd_pre_upper_all <- reshape2::melt(PC1_1.4sd_pre_upper_all, id.vars = "group")
names(PC1_1.4sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_1.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.4sd_pre_lower <- as.data.frame(cbind(PC1_1.4sd_pre$x, PC1_1.4sd_pre$lower))
names(PC1_1.4sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.4sd_pre_lower$age, PC1_1.4sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_1.4sd_pre_lower, indices) {
  d <- PC1_1.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_1.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_1.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_1.4sd_pre_lower_50)[1] <- "PC1"

PC1_1.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_1.4sd_pre_lower_all)[1:ncol(PC1_1.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_1.4sd_pre_lower_all$group <- seq(1,(nrow(PC1_1.4sd_pre_lower_all)), 1)

PC1_1.4sd_pre_lower_all <- reshape2::melt(PC1_1.4sd_pre_lower_all, id.vars = "group")
names(PC1_1.4sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_1.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.2sd_pre_upper <- as.data.frame(cbind(PC1_1.2sd_pre$x, PC1_1.2sd_pre$upper))
names(PC1_1.2sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.2sd_pre_upper$age, PC1_1.2sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_1.2sd_pre_upper, indices) {
  d <- PC1_1.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_1.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_1.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_1.2sd_pre_upper_50)[1] <- "PC1"

PC1_1.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_1.2sd_pre_upper_all)[1:ncol(PC1_1.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_1.2sd_pre_upper_all$group <- seq(1,(nrow(PC1_1.2sd_pre_upper_all)), 1)

PC1_1.2sd_pre_upper_all <- reshape2::melt(PC1_1.2sd_pre_upper_all, id.vars = "group")
names(PC1_1.2sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_1.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.2sd_pre_lower <- as.data.frame(cbind(PC1_1.2sd_pre$x, PC1_1.2sd_pre$lower))
names(PC1_1.2sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.2sd_pre_lower$age, PC1_1.2sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_1.2sd_pre_lower, indices) {
  d <- PC1_1.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_1.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_1.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_1.2sd_pre_lower_50)[1] <- "PC1"

PC1_1.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_1.2sd_pre_lower_all)[1:ncol(PC1_1.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_1.2sd_pre_lower_all$group <- seq(1,(nrow(PC1_1.2sd_pre_lower_all)), 1)

PC1_1.2sd_pre_lower_all <- reshape2::melt(PC1_1.2sd_pre_lower_all, id.vars = "group")
names(PC1_1.2sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_1.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.0sd_pre_upper <- as.data.frame(cbind(PC1_1.0sd_pre$x, PC1_1.0sd_pre$upper))
names(PC1_1.0sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.0sd_pre_upper$age, PC1_1.0sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_1.0sd_pre_upper, indices) {
  d <- PC1_1.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_1.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_1.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_1.0sd_pre_upper_50)[1] <- "PC1"

PC1_1.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_1.0sd_pre_upper_all)[1:ncol(PC1_1.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_1.0sd_pre_upper_all$group <- seq(1,(nrow(PC1_1.0sd_pre_upper_all)), 1)

PC1_1.0sd_pre_upper_all <- reshape2::melt(PC1_1.0sd_pre_upper_all, id.vars = "group")
names(PC1_1.0sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_1.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_1.0sd_pre_lower <- as.data.frame(cbind(PC1_1.0sd_pre$x, PC1_1.0sd_pre$lower))
names(PC1_1.0sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_1.0sd_pre_lower$age, PC1_1.0sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_1.0sd_pre_lower, indices) {
  d <- PC1_1.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_1.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_1.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_1.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_1.0sd_pre_lower_50)[1] <- "PC1"

PC1_1.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_1.0sd_pre_lower_all)[1:ncol(PC1_1.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_1.0sd_pre_lower_all$group <- seq(1,(nrow(PC1_1.0sd_pre_lower_all)), 1)

PC1_1.0sd_pre_lower_all <- reshape2::melt(PC1_1.0sd_pre_lower_all, id.vars = "group")
names(PC1_1.0sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_1.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.8sd_pre_upper <- as.data.frame(cbind(PC1_0.8sd_pre$x, PC1_0.8sd_pre$upper))
names(PC1_0.8sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.8sd_pre_upper$age, PC1_0.8sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_0.8sd_pre_upper, indices) {
  d <- PC1_0.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_0.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_0.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_0.8sd_pre_upper_50)[1] <- "PC1"

PC1_0.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_0.8sd_pre_upper_all)[1:ncol(PC1_0.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_0.8sd_pre_upper_all$group <- seq(1,(nrow(PC1_0.8sd_pre_upper_all)), 1)

PC1_0.8sd_pre_upper_all <- reshape2::melt(PC1_0.8sd_pre_upper_all, id.vars = "group")
names(PC1_0.8sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_0.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.8sd_pre_lower <- as.data.frame(cbind(PC1_0.8sd_pre$x, PC1_0.8sd_pre$lower))
names(PC1_0.8sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.8sd_pre_lower$age, PC1_0.8sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_0.8sd_pre_lower, indices) {
  d <- PC1_0.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_0.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_0.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_0.8sd_pre_lower_50)[1] <- "PC1"

PC1_0.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_0.8sd_pre_lower_all)[1:ncol(PC1_0.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_0.8sd_pre_lower_all$group <- seq(1,(nrow(PC1_0.8sd_pre_lower_all)), 1)

PC1_0.8sd_pre_lower_all <- reshape2::melt(PC1_0.8sd_pre_lower_all, id.vars = "group")
names(PC1_0.8sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_0.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.6sd_pre_upper <- as.data.frame(cbind(PC1_0.6sd_pre$x, PC1_0.6sd_pre$upper))
names(PC1_0.6sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.6sd_pre_upper$age, PC1_0.6sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_0.6sd_pre_upper, indices) {
  d <- PC1_0.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_0.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_0.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_0.6sd_pre_upper_50)[1] <- "PC1"

PC1_0.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_0.6sd_pre_upper_all)[1:ncol(PC1_0.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_0.6sd_pre_upper_all$group <- seq(1,(nrow(PC1_0.6sd_pre_upper_all)), 1)

PC1_0.6sd_pre_upper_all <- reshape2::melt(PC1_0.6sd_pre_upper_all, id.vars = "group")
names(PC1_0.6sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_0.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.6sd_pre_lower <- as.data.frame(cbind(PC1_0.6sd_pre$x, PC1_0.6sd_pre$lower))
names(PC1_0.6sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.6sd_pre_lower$age, PC1_0.6sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_0.6sd_pre_lower, indices) {
  d <- PC1_0.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_0.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_0.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_0.6sd_pre_lower_50)[1] <- "PC1"

PC1_0.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_0.6sd_pre_lower_all)[1:ncol(PC1_0.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_0.6sd_pre_lower_all$group <- seq(1,(nrow(PC1_0.6sd_pre_lower_all)), 1)

PC1_0.6sd_pre_lower_all <- reshape2::melt(PC1_0.6sd_pre_lower_all, id.vars = "group")
names(PC1_0.6sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_0.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.4sd_pre_upper <- as.data.frame(cbind(PC1_0.4sd_pre$x, PC1_0.4sd_pre$upper))
names(PC1_0.4sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.4sd_pre_upper$age, PC1_0.4sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_0.4sd_pre_upper, indices) {
  d <- PC1_0.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_0.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_0.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_0.4sd_pre_upper_50)[1] <- "PC1"

PC1_0.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_0.4sd_pre_upper_all)[1:ncol(PC1_0.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_0.4sd_pre_upper_all$group <- seq(1,(nrow(PC1_0.4sd_pre_upper_all)), 1)

PC1_0.4sd_pre_upper_all <- reshape2::melt(PC1_0.4sd_pre_upper_all, id.vars = "group")
names(PC1_0.4sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_0.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.4sd_pre_lower <- as.data.frame(cbind(PC1_0.4sd_pre$x, PC1_0.4sd_pre$lower))
names(PC1_0.4sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.4sd_pre_lower$age, PC1_0.4sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_0.4sd_pre_lower, indices) {
  d <- PC1_0.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_0.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_0.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_0.4sd_pre_lower_50)[1] <- "PC1"

PC1_0.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_0.4sd_pre_lower_all)[1:ncol(PC1_0.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_0.4sd_pre_lower_all$group <- seq(1,(nrow(PC1_0.4sd_pre_lower_all)), 1)

PC1_0.4sd_pre_lower_all <- reshape2::melt(PC1_0.4sd_pre_lower_all, id.vars = "group")
names(PC1_0.4sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_0.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.2sd_pre_upper <- as.data.frame(cbind(PC1_0.2sd_pre$x, PC1_0.2sd_pre$upper))
names(PC1_0.2sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.2sd_pre_upper$age, PC1_0.2sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_0.2sd_pre_upper, indices) {
  d <- PC1_0.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_0.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_0.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_0.2sd_pre_upper_50)[1] <- "PC1"

PC1_0.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_0.2sd_pre_upper_all)[1:ncol(PC1_0.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_0.2sd_pre_upper_all$group <- seq(1,(nrow(PC1_0.2sd_pre_upper_all)), 1)

PC1_0.2sd_pre_upper_all <- reshape2::melt(PC1_0.2sd_pre_upper_all, id.vars = "group")
names(PC1_0.2sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_0.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.2sd_pre_lower <- as.data.frame(cbind(PC1_0.2sd_pre$x, PC1_0.2sd_pre$lower))
names(PC1_0.2sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.2sd_pre_lower$age, PC1_0.2sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_0.2sd_pre_lower, indices) {
  d <- PC1_0.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_0.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_0.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_0.2sd_pre_lower_50)[1] <- "PC1"

PC1_0.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_0.2sd_pre_lower_all)[1:ncol(PC1_0.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_0.2sd_pre_lower_all$group <- seq(1,(nrow(PC1_0.2sd_pre_lower_all)), 1)

PC1_0.2sd_pre_lower_all <- reshape2::melt(PC1_0.2sd_pre_lower_all, id.vars = "group")
names(PC1_0.2sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_0.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC1 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.0sd_pre_upper <- as.data.frame(cbind(PC1_0.0sd_pre$x, PC1_0.0sd_pre$upper))
names(PC1_0.0sd_pre_upper)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.0sd_pre_upper$age, PC1_0.0sd_pre_upper$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC1_0.0sd_pre_upper, indices) {
  d <- PC1_0.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC1_0.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC1_0.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC1_0.0sd_pre_upper_50)[1] <- "PC1"

PC1_0.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC1_0.0sd_pre_upper_all)[1:ncol(PC1_0.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC1_0.0sd_pre_upper_all$group <- seq(1,(nrow(PC1_0.0sd_pre_upper_all)), 1)

PC1_0.0sd_pre_upper_all <- reshape2::melt(PC1_0.0sd_pre_upper_all, id.vars = "group")
names(PC1_0.0sd_pre_upper_all)[2:3] <- c("age", "PC1")

PC1_0.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC1 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC1_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC1, degree=1, nsigma = nsig, span = span)
PC1_0.0sd_pre_lower <- as.data.frame(cbind(PC1_0.0sd_pre$x, PC1_0.0sd_pre$lower))
names(PC1_0.0sd_pre_lower)[1:2] <- c("age","PC1")

# CV spans
loess.predict <- loess.as(PC1_0.0sd_pre_lower$age, PC1_0.0sd_pre_lower$PC1, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC1_0.0sd_pre_lower, indices) {
  d <- PC1_0.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC1 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC1_0.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC1_0.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC1_0.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC1_0.0sd_pre_lower_50)[1] <- "PC1"

PC1_0.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC1_0.0sd_pre_lower_all)[1:ncol(PC1_0.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC1_0.0sd_pre_lower_all$group <- seq(1,(nrow(PC1_0.0sd_pre_lower_all)), 1)

PC1_0.0sd_pre_lower_all <- reshape2::melt(PC1_0.0sd_pre_lower_all, id.vars = "group")
names(PC1_0.0sd_pre_lower_all)[2:3] <- c("age", "PC1")

PC1_0.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

b <- ggplot() + 
  #geom_line(data = PC1_2sd_pre_upper_all, aes(age2,PC1, group = group)) +
  #geom_line(data = PC1_2sd_pre_lower_all, aes(age2,PC1, group = group)) +
  geom_line(data = PC1_2sd_pre_upper_50, aes(age,PC1), colour = 'blue') +
  geom_line(data = PC1_2sd_pre_lower_50, aes(age,PC1), colour = 'blue') +
  geom_line(data = PC1_1.6sd_pre_upper_50, aes(age,PC1), colour = 'green') +
  geom_line(data = PC1_1.6sd_pre_lower_50, aes(age,PC1), colour = 'green') +
  geom_line(data = PC1_1.4sd_pre_upper_50, aes(age,PC1), colour = 'orange') +
  geom_line(data = PC1_1.4sd_pre_lower_50, aes(age,PC1), colour = 'orange') +
  geom_line(data = PC1_1.2sd_pre_upper_50, aes(age,PC1), colour = 'pink') +
  geom_line(data = PC1_1.2sd_pre_lower_50, aes(age,PC1), colour = 'pink') +
  geom_line(data = PC1_1.0sd_pre_upper_50, aes(age,PC1), colour = 'brown') +
  geom_line(data = PC1_1.0sd_pre_lower_50, aes(age,PC1), colour = 'brown') +
  geom_line(data = PC1_0.8sd_pre_upper_50, aes(age,PC1), colour = 'yellow') +
  geom_line(data = PC1_0.8sd_pre_lower_50, aes(age,PC1), colour = 'yellow') +
  geom_line(data = PC1_0.6sd_pre_upper_50, aes(age,PC1), colour = 'red') +
  geom_line(data = PC1_0.6sd_pre_lower_50, aes(age,PC1), colour = 'red') +
  geom_line(data = PC1_0.4sd_pre_upper_50, aes(age,PC1), colour = 'blue') +
  geom_line(data = PC1_0.4sd_pre_lower_50, aes(age,PC1), colour = 'blue') +
  geom_line(data = PC1_0.2sd_pre_upper_50, aes(age,PC1), colour = 'green') +
  geom_line(data = PC1_0.2sd_pre_lower_50, aes(age,PC1), colour = 'green') +
  geom_line(data = PC1_0.0sd_pre_upper_50, aes(age,PC1), colour = 'orange') +
  geom_line(data = PC1_0.0sd_pre_lower_50, aes(age,PC1), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC1)) + 
  #geom_point(data = PC1_2sd_pre_upper, aes(age, PC1), colour = 'red') +
  scale_x_reverse(limits = c(4000,540)) + theme_bw()

grid.arrange(b,a,ncol=2)


##### Phanerozoic PC2 #####

coords.all_Phan <- subset(coords.all2, age < 541)

# 2 sigma

nsig <- 2

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_2sd_phan_upper <- as.data.frame(cbind(PC2_2sd_phan$x, PC2_2sd_phan$upper))
names(PC2_2sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_2sd_phan_upper$age, PC2_2sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_2sd_phan_upper, indices) {
  d <- PC2_2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_2sd_phan_upper_50)[1] <- "PC2"

PC2_2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_2sd_phan_upper_all)[1:ncol(PC2_2sd_phan_upper_all)] <- seq(0,540, 1)
PC2_2sd_phan_upper_all$group <- seq(1,(nrow(PC2_2sd_phan_upper_all)), 1)

PC2_2sd_phan_upper_all <- reshape2::melt(PC2_2sd_phan_upper_all, id.vars = "group")
names(PC2_2sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_2sd_phan_lower <- as.data.frame(cbind(PC2_2sd_phan$x, PC2_2sd_phan$lower))
names(PC2_2sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_2sd_phan_lower$age, PC2_2sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_2sd_phan_lower, indices) {
  d <- PC2_2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_2sd_phan_lower_50)[1] <- "PC2"

PC2_2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_2sd_phan_lower_all)[1:ncol(PC2_2sd_phan_lower_all)] <- seq(0,540, 1)
PC2_2sd_phan_lower_all$group <- seq(1,(nrow(PC2_2sd_phan_lower_all)), 1)

PC2_2sd_phan_lower_all <- reshape2::melt(PC2_2sd_phan_lower_all, id.vars = "group")
names(PC2_2sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.8sd_phan_upper <- as.data.frame(cbind(PC2_1.8sd_phan$x, PC2_1.8sd_phan$upper))
names(PC2_1.8sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.8sd_phan_upper$age, PC2_1.8sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_1.8sd_phan_upper, indices) {
  d <- PC2_1.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_1.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_1.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_1.8sd_phan_upper_50)[1] <- "PC2"

PC2_1.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_1.8sd_phan_upper_all)[1:ncol(PC2_1.8sd_phan_upper_all)] <- seq(0,540, 1)
PC2_1.8sd_phan_upper_all$group <- seq(1,(nrow(PC2_1.8sd_phan_upper_all)), 1)

PC2_1.8sd_phan_upper_all <- reshape2::melt(PC2_1.8sd_phan_upper_all, id.vars = "group")
names(PC2_1.8sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_1.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.8sd_phan_lower <- as.data.frame(cbind(PC2_1.8sd_phan$x, PC2_1.8sd_phan$lower))
names(PC2_1.8sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.8sd_phan_lower$age, PC2_1.8sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_1.8sd_phan_lower, indices) {
  d <- PC2_1.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_1.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_1.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_1.8sd_phan_lower_50)[1] <- "PC2"

PC2_1.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_1.8sd_phan_lower_all)[1:ncol(PC2_1.8sd_phan_lower_all)] <- seq(0,540, 1)
PC2_1.8sd_phan_lower_all$group <- seq(1,(nrow(PC2_1.8sd_phan_lower_all)), 1)

PC2_1.8sd_phan_lower_all <- reshape2::melt(PC2_1.8sd_phan_lower_all, id.vars = "group")
names(PC2_1.8sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_1.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.6sd_phan_upper <- as.data.frame(cbind(PC2_1.6sd_phan$x, PC2_1.6sd_phan$upper))
names(PC2_1.6sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.6sd_phan_upper$age, PC2_1.6sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_1.6sd_phan_upper, indices) {
  d <- PC2_1.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_1.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_1.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_1.6sd_phan_upper_50)[1] <- "PC2"

PC2_1.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_1.6sd_phan_upper_all)[1:ncol(PC2_1.6sd_phan_upper_all)] <- seq(0,540, 1)
PC2_1.6sd_phan_upper_all$group <- seq(1,(nrow(PC2_1.6sd_phan_upper_all)), 1)

PC2_1.6sd_phan_upper_all <- reshape2::melt(PC2_1.6sd_phan_upper_all, id.vars = "group")
names(PC2_1.6sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_1.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.6sd_phan_lower <- as.data.frame(cbind(PC2_1.6sd_phan$x, PC2_1.6sd_phan$lower))
names(PC2_1.6sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.6sd_phan_lower$age, PC2_1.6sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_1.6sd_phan_lower, indices) {
  d <- PC2_1.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_1.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_1.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_1.6sd_phan_lower_50)[1] <- "PC2"

PC2_1.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_1.6sd_phan_lower_all)[1:ncol(PC2_1.6sd_phan_lower_all)] <- seq(0,540, 1)
PC2_1.6sd_phan_lower_all$group <- seq(1,(nrow(PC2_1.6sd_phan_lower_all)), 1)

PC2_1.6sd_phan_lower_all <- reshape2::melt(PC2_1.6sd_phan_lower_all, id.vars = "group")
names(PC2_1.6sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_1.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.4sd_phan_upper <- as.data.frame(cbind(PC2_1.4sd_phan$x, PC2_1.4sd_phan$upper))
names(PC2_1.4sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.4sd_phan_upper$age, PC2_1.4sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_1.4sd_phan_upper, indices) {
  d <- PC2_1.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_1.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_1.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_1.4sd_phan_upper_50)[1] <- "PC2"

PC2_1.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_1.4sd_phan_upper_all)[1:ncol(PC2_1.4sd_phan_upper_all)] <- seq(0,540, 1)
PC2_1.4sd_phan_upper_all$group <- seq(1,(nrow(PC2_1.4sd_phan_upper_all)), 1)

PC2_1.4sd_phan_upper_all <- reshape2::melt(PC2_1.4sd_phan_upper_all, id.vars = "group")
names(PC2_1.4sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_1.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.4sd_phan_lower <- as.data.frame(cbind(PC2_1.4sd_phan$x, PC2_1.4sd_phan$lower))
names(PC2_1.4sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.4sd_phan_lower$age, PC2_1.4sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_1.4sd_phan_lower, indices) {
  d <- PC2_1.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_1.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_1.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_1.4sd_phan_lower_50)[1] <- "PC2"

PC2_1.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_1.4sd_phan_lower_all)[1:ncol(PC2_1.4sd_phan_lower_all)] <- seq(0,540, 1)
PC2_1.4sd_phan_lower_all$group <- seq(1,(nrow(PC2_1.4sd_phan_lower_all)), 1)

PC2_1.4sd_phan_lower_all <- reshape2::melt(PC2_1.4sd_phan_lower_all, id.vars = "group")
names(PC2_1.4sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_1.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.2sd_phan_upper <- as.data.frame(cbind(PC2_1.2sd_phan$x, PC2_1.2sd_phan$upper))
names(PC2_1.2sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.2sd_phan_upper$age, PC2_1.2sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_1.2sd_phan_upper, indices) {
  d <- PC2_1.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_1.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_1.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_1.2sd_phan_upper_50)[1] <- "PC2"

PC2_1.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_1.2sd_phan_upper_all)[1:ncol(PC2_1.2sd_phan_upper_all)] <- seq(0,540, 1)
PC2_1.2sd_phan_upper_all$group <- seq(1,(nrow(PC2_1.2sd_phan_upper_all)), 1)

PC2_1.2sd_phan_upper_all <- reshape2::melt(PC2_1.2sd_phan_upper_all, id.vars = "group")
names(PC2_1.2sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_1.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.2sd_phan_lower <- as.data.frame(cbind(PC2_1.2sd_phan$x, PC2_1.2sd_phan$lower))
names(PC2_1.2sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.2sd_phan_lower$age, PC2_1.2sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_1.2sd_phan_lower, indices) {
  d <- PC2_1.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_1.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_1.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_1.2sd_phan_lower_50)[1] <- "PC2"

PC2_1.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_1.2sd_phan_lower_all)[1:ncol(PC2_1.2sd_phan_lower_all)] <- seq(0,540, 1)
PC2_1.2sd_phan_lower_all$group <- seq(1,(nrow(PC2_1.2sd_phan_lower_all)), 1)

PC2_1.2sd_phan_lower_all <- reshape2::melt(PC2_1.2sd_phan_lower_all, id.vars = "group")
names(PC2_1.2sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_1.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.0sd_phan_upper <- as.data.frame(cbind(PC2_1.0sd_phan$x, PC2_1.0sd_phan$upper))
names(PC2_1.0sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.0sd_phan_upper$age, PC2_1.0sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_1.0sd_phan_upper, indices) {
  d <- PC2_1.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_1.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_1.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_1.0sd_phan_upper_50)[1] <- "PC2"

PC2_1.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_1.0sd_phan_upper_all)[1:ncol(PC2_1.0sd_phan_upper_all)] <- seq(0,540, 1)
PC2_1.0sd_phan_upper_all$group <- seq(1,(nrow(PC2_1.0sd_phan_upper_all)), 1)

PC2_1.0sd_phan_upper_all <- reshape2::melt(PC2_1.0sd_phan_upper_all, id.vars = "group")
names(PC2_1.0sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_1.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.0sd_phan_lower <- as.data.frame(cbind(PC2_1.0sd_phan$x, PC2_1.0sd_phan$lower))
names(PC2_1.0sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.0sd_phan_lower$age, PC2_1.0sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_1.0sd_phan_lower, indices) {
  d <- PC2_1.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_1.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_1.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_1.0sd_phan_lower_50)[1] <- "PC2"

PC2_1.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_1.0sd_phan_lower_all)[1:ncol(PC2_1.0sd_phan_lower_all)] <- seq(0,540, 1)
PC2_1.0sd_phan_lower_all$group <- seq(1,(nrow(PC2_1.0sd_phan_lower_all)), 1)

PC2_1.0sd_phan_lower_all <- reshape2::melt(PC2_1.0sd_phan_lower_all, id.vars = "group")
names(PC2_1.0sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_1.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.8sd_phan_upper <- as.data.frame(cbind(PC2_0.8sd_phan$x, PC2_0.8sd_phan$upper))
names(PC2_0.8sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.8sd_phan_upper$age, PC2_0.8sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_0.8sd_phan_upper, indices) {
  d <- PC2_0.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_0.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_0.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_0.8sd_phan_upper_50)[1] <- "PC2"

PC2_0.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_0.8sd_phan_upper_all)[1:ncol(PC2_0.8sd_phan_upper_all)] <- seq(0,540, 1)
PC2_0.8sd_phan_upper_all$group <- seq(1,(nrow(PC2_0.8sd_phan_upper_all)), 1)

PC2_0.8sd_phan_upper_all <- reshape2::melt(PC2_0.8sd_phan_upper_all, id.vars = "group")
names(PC2_0.8sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_0.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.8sd_phan_lower <- as.data.frame(cbind(PC2_0.8sd_phan$x, PC2_0.8sd_phan$lower))
names(PC2_0.8sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.8sd_phan_lower$age, PC2_0.8sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_0.8sd_phan_lower, indices) {
  d <- PC2_0.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_0.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_0.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_0.8sd_phan_lower_50)[1] <- "PC2"

PC2_0.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_0.8sd_phan_lower_all)[1:ncol(PC2_0.8sd_phan_lower_all)] <- seq(0,540, 1)
PC2_0.8sd_phan_lower_all$group <- seq(1,(nrow(PC2_0.8sd_phan_lower_all)), 1)

PC2_0.8sd_phan_lower_all <- reshape2::melt(PC2_0.8sd_phan_lower_all, id.vars = "group")
names(PC2_0.8sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_0.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.6sd_phan_upper <- as.data.frame(cbind(PC2_0.6sd_phan$x, PC2_0.6sd_phan$upper))
names(PC2_0.6sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.6sd_phan_upper$age, PC2_0.6sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_0.6sd_phan_upper, indices) {
  d <- PC2_0.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_0.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_0.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_0.6sd_phan_upper_50)[1] <- "PC2"

PC2_0.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_0.6sd_phan_upper_all)[1:ncol(PC2_0.6sd_phan_upper_all)] <- seq(0,540, 1)
PC2_0.6sd_phan_upper_all$group <- seq(1,(nrow(PC2_0.6sd_phan_upper_all)), 1)

PC2_0.6sd_phan_upper_all <- reshape2::melt(PC2_0.6sd_phan_upper_all, id.vars = "group")
names(PC2_0.6sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_0.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.6sd_phan_lower <- as.data.frame(cbind(PC2_0.6sd_phan$x, PC2_0.6sd_phan$lower))
names(PC2_0.6sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.6sd_phan_lower$age, PC2_0.6sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_0.6sd_phan_lower, indices) {
  d <- PC2_0.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_0.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_0.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_0.6sd_phan_lower_50)[1] <- "PC2"

PC2_0.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_0.6sd_phan_lower_all)[1:ncol(PC2_0.6sd_phan_lower_all)] <- seq(0,540, 1)
PC2_0.6sd_phan_lower_all$group <- seq(1,(nrow(PC2_0.6sd_phan_lower_all)), 1)

PC2_0.6sd_phan_lower_all <- reshape2::melt(PC2_0.6sd_phan_lower_all, id.vars = "group")
names(PC2_0.6sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_0.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.4sd_phan_upper <- as.data.frame(cbind(PC2_0.4sd_phan$x, PC2_0.4sd_phan$upper))
names(PC2_0.4sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.4sd_phan_upper$age, PC2_0.4sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_0.4sd_phan_upper, indices) {
  d <- PC2_0.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_0.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_0.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_0.4sd_phan_upper_50)[1] <- "PC2"

PC2_0.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_0.4sd_phan_upper_all)[1:ncol(PC2_0.4sd_phan_upper_all)] <- seq(0,540, 1)
PC2_0.4sd_phan_upper_all$group <- seq(1,(nrow(PC2_0.4sd_phan_upper_all)), 1)

PC2_0.4sd_phan_upper_all <- reshape2::melt(PC2_0.4sd_phan_upper_all, id.vars = "group")
names(PC2_0.4sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_0.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.4sd_phan_lower <- as.data.frame(cbind(PC2_0.4sd_phan$x, PC2_0.4sd_phan$lower))
names(PC2_0.4sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.4sd_phan_lower$age, PC2_0.4sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_0.4sd_phan_lower, indices) {
  d <- PC2_0.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_0.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_0.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_0.4sd_phan_lower_50)[1] <- "PC2"

PC2_0.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_0.4sd_phan_lower_all)[1:ncol(PC2_0.4sd_phan_lower_all)] <- seq(0,540, 1)
PC2_0.4sd_phan_lower_all$group <- seq(1,(nrow(PC2_0.4sd_phan_lower_all)), 1)

PC2_0.4sd_phan_lower_all <- reshape2::melt(PC2_0.4sd_phan_lower_all, id.vars = "group")
names(PC2_0.4sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_0.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.2sd_phan_upper <- as.data.frame(cbind(PC2_0.2sd_phan$x, PC2_0.2sd_phan$upper))
names(PC2_0.2sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.2sd_phan_upper$age, PC2_0.2sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_0.2sd_phan_upper, indices) {
  d <- PC2_0.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_0.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_0.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_0.2sd_phan_upper_50)[1] <- "PC2"

PC2_0.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_0.2sd_phan_upper_all)[1:ncol(PC2_0.2sd_phan_upper_all)] <- seq(0,540, 1)
PC2_0.2sd_phan_upper_all$group <- seq(1,(nrow(PC2_0.2sd_phan_upper_all)), 1)

PC2_0.2sd_phan_upper_all <- reshape2::melt(PC2_0.2sd_phan_upper_all, id.vars = "group")
names(PC2_0.2sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_0.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.2sd_phan_lower <- as.data.frame(cbind(PC2_0.2sd_phan$x, PC2_0.2sd_phan$lower))
names(PC2_0.2sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.2sd_phan_lower$age, PC2_0.2sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_0.2sd_phan_lower, indices) {
  d <- PC2_0.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_0.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_0.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_0.2sd_phan_lower_50)[1] <- "PC2"

PC2_0.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_0.2sd_phan_lower_all)[1:ncol(PC2_0.2sd_phan_lower_all)] <- seq(0,540, 1)
PC2_0.2sd_phan_lower_all$group <- seq(1,(nrow(PC2_0.2sd_phan_lower_all)), 1)

PC2_0.2sd_phan_lower_all <- reshape2::melt(PC2_0.2sd_phan_lower_all, id.vars = "group")
names(PC2_0.2sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_0.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.0sd_phan_upper <- as.data.frame(cbind(PC2_0.0sd_phan$x, PC2_0.0sd_phan$upper))
names(PC2_0.0sd_phan_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.0sd_phan_upper$age, PC2_0.0sd_phan_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_0.0sd_phan_upper, indices) {
  d <- PC2_0.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_0.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_0.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC2_0.0sd_phan_upper_50)[1] <- "PC2"

PC2_0.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC2_0.0sd_phan_upper_all)[1:ncol(PC2_0.0sd_phan_upper_all)] <- seq(0,540, 1)
PC2_0.0sd_phan_upper_all$group <- seq(1,(nrow(PC2_0.0sd_phan_upper_all)), 1)

PC2_0.0sd_phan_upper_all <- reshape2::melt(PC2_0.0sd_phan_upper_all, id.vars = "group")
names(PC2_0.0sd_phan_upper_all)[2:3] <- c("age", "PC2")

PC2_0.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.0sd_phan_lower <- as.data.frame(cbind(PC2_0.0sd_phan$x, PC2_0.0sd_phan$lower))
names(PC2_0.0sd_phan_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.0sd_phan_lower$age, PC2_0.0sd_phan_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_0.0sd_phan_lower, indices) {
  d <- PC2_0.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_0.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_0.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC2_0.0sd_phan_lower_50)[1] <- "PC2"

PC2_0.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC2_0.0sd_phan_lower_all)[1:ncol(PC2_0.0sd_phan_lower_all)] <- seq(0,540, 1)
PC2_0.0sd_phan_lower_all$group <- seq(1,(nrow(PC2_0.0sd_phan_lower_all)), 1)

PC2_0.0sd_phan_lower_all <- reshape2::melt(PC2_0.0sd_phan_lower_all, id.vars = "group")
names(PC2_0.0sd_phan_lower_all)[2:3] <- c("age", "PC2")

PC2_0.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

a <- ggplot() + 
  #geom_line(data = PC2_2sd_phan_upper_all, aes(age2,PC2, group = group)) +
  #geom_line(data = PC2_2sd_phan_lower_all, aes(age2,PC2, group = group)) +
  geom_line(data = PC2_2sd_phan_upper_50, aes(age,PC2), colour = 'blue') +
  geom_line(data = PC2_2sd_phan_lower_50, aes(age,PC2), colour = 'blue') +
  geom_line(data = PC2_1.6sd_phan_upper_50, aes(age,PC2), colour = 'green') +
  geom_line(data = PC2_1.6sd_phan_lower_50, aes(age,PC2), colour = 'green') +
  geom_line(data = PC2_1.4sd_phan_upper_50, aes(age,PC2), colour = 'orange') +
  geom_line(data = PC2_1.4sd_phan_lower_50, aes(age,PC2), colour = 'orange') +
  geom_line(data = PC2_1.2sd_phan_upper_50, aes(age,PC2), colour = 'pink') +
  geom_line(data = PC2_1.2sd_phan_lower_50, aes(age,PC2), colour = 'pink') +
  geom_line(data = PC2_1.0sd_phan_upper_50, aes(age,PC2), colour = 'brown') +
  geom_line(data = PC2_1.0sd_phan_lower_50, aes(age,PC2), colour = 'brown') +
  geom_line(data = PC2_0.8sd_phan_upper_50, aes(age,PC2), colour = 'yellow') +
  geom_line(data = PC2_0.8sd_phan_lower_50, aes(age,PC2), colour = 'yellow') +
  geom_line(data = PC2_0.6sd_phan_upper_50, aes(age,PC2), colour = 'red') +
  geom_line(data = PC2_0.6sd_phan_lower_50, aes(age,PC2), colour = 'red') +
  geom_line(data = PC2_0.4sd_phan_upper_50, aes(age,PC2), colour = 'blue') +
  geom_line(data = PC2_0.4sd_phan_lower_50, aes(age,PC2), colour = 'blue') +
  geom_line(data = PC2_0.2sd_phan_upper_50, aes(age,PC2), colour = 'green') +
  geom_line(data = PC2_0.2sd_phan_lower_50, aes(age,PC2), colour = 'green') +
  geom_line(data = PC2_0.0sd_phan_upper_50, aes(age,PC2), colour = 'orange') +
  geom_line(data = PC2_0.0sd_phan_lower_50, aes(age,PC2), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC2)) + 
  #geom_point(data = PC2_2sd_phan_upper, aes(age, PC2), colour = 'red') +
  scale_x_reverse(limits = c(540,0)) + theme_bw()

##### Precambrian PC2 #####

coords.all_pre <- subset(coords.all2, age > 541)

# 2 sigma

nsig <- 2

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_2sd_pre_upper <- as.data.frame(cbind(PC2_2sd_pre$x, PC2_2sd_pre$upper))
names(PC2_2sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_2sd_pre_upper$age, PC2_2sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[1], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_2sd_pre_upper, indices) {
  d <- PC2_2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_2sd_pre_upper_50)[1] <- "PC2"

PC2_2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_2sd_pre_upper_all)[1:ncol(PC2_2sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_2sd_pre_upper_all$group <- seq(1,(nrow(PC2_2sd_pre_upper_all)), 1)

PC2_2sd_pre_upper_all <- reshape2::melt(PC2_2sd_pre_upper_all, id.vars = "group")
names(PC2_2sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_2sd_pre_lower <- as.data.frame(cbind(PC2_2sd_pre$x, PC2_2sd_pre$lower))
names(PC2_2sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_2sd_pre_lower$age, PC2_2sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_2sd_pre_lower, indices) {
  d <- PC2_2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_2sd_pre_lower_50)[1] <- "PC2"

PC2_2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_2sd_pre_lower_all)[1:ncol(PC2_2sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_2sd_pre_lower_all$group <- seq(1,(nrow(PC2_2sd_pre_lower_all)), 1)

PC2_2sd_pre_lower_all <- reshape2::melt(PC2_2sd_pre_lower_all, id.vars = "group")
names(PC2_2sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.8sd_pre_upper <- as.data.frame(cbind(PC2_1.8sd_pre$x, PC2_1.8sd_pre$upper))
names(PC2_1.8sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.8sd_pre_upper$age, PC2_1.8sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_1.8sd_pre_upper, indices) {
  d <- PC2_1.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_1.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_1.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_1.8sd_pre_upper_50)[1] <- "PC2"

PC2_1.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_1.8sd_pre_upper_all)[1:ncol(PC2_1.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_1.8sd_pre_upper_all$group <- seq(1,(nrow(PC2_1.8sd_pre_upper_all)), 1)

PC2_1.8sd_pre_upper_all <- reshape2::melt(PC2_1.8sd_pre_upper_all, id.vars = "group")
names(PC2_1.8sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_1.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.8sd_pre_lower <- as.data.frame(cbind(PC2_1.8sd_pre$x, PC2_1.8sd_pre$lower))
names(PC2_1.8sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.8sd_pre_lower$age, PC2_1.8sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_1.8sd_pre_lower, indices) {
  d <- PC2_1.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_1.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_1.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_1.8sd_pre_lower_50)[1] <- "PC2"

PC2_1.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_1.8sd_pre_lower_all)[1:ncol(PC2_1.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_1.8sd_pre_lower_all$group <- seq(1,(nrow(PC2_1.8sd_pre_lower_all)), 1)

PC2_1.8sd_pre_lower_all <- reshape2::melt(PC2_1.8sd_pre_lower_all, id.vars = "group")
names(PC2_1.8sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_1.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.6sd_pre_upper <- as.data.frame(cbind(PC2_1.6sd_pre$x, PC2_1.6sd_pre$upper))
names(PC2_1.6sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.6sd_pre_upper$age, PC2_1.6sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_1.6sd_pre_upper, indices) {
  d <- PC2_1.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_1.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_1.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_1.6sd_pre_upper_50)[1] <- "PC2"

PC2_1.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_1.6sd_pre_upper_all)[1:ncol(PC2_1.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_1.6sd_pre_upper_all$group <- seq(1,(nrow(PC2_1.6sd_pre_upper_all)), 1)

PC2_1.6sd_pre_upper_all <- reshape2::melt(PC2_1.6sd_pre_upper_all, id.vars = "group")
names(PC2_1.6sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_1.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.6sd_pre_lower <- as.data.frame(cbind(PC2_1.6sd_pre$x, PC2_1.6sd_pre$lower))
names(PC2_1.6sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.6sd_pre_lower$age, PC2_1.6sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_1.6sd_pre_lower, indices) {
  d <- PC2_1.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_1.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_1.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_1.6sd_pre_lower_50)[1] <- "PC2"

PC2_1.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_1.6sd_pre_lower_all)[1:ncol(PC2_1.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_1.6sd_pre_lower_all$group <- seq(1,(nrow(PC2_1.6sd_pre_lower_all)), 1)

PC2_1.6sd_pre_lower_all <- reshape2::melt(PC2_1.6sd_pre_lower_all, id.vars = "group")
names(PC2_1.6sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_1.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.4sd_pre_upper <- as.data.frame(cbind(PC2_1.4sd_pre$x, PC2_1.4sd_pre$upper))
names(PC2_1.4sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.4sd_pre_upper$age, PC2_1.4sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_1.4sd_pre_upper, indices) {
  d <- PC2_1.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_1.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_1.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_1.4sd_pre_upper_50)[1] <- "PC2"

PC2_1.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_1.4sd_pre_upper_all)[1:ncol(PC2_1.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_1.4sd_pre_upper_all$group <- seq(1,(nrow(PC2_1.4sd_pre_upper_all)), 1)

PC2_1.4sd_pre_upper_all <- reshape2::melt(PC2_1.4sd_pre_upper_all, id.vars = "group")
names(PC2_1.4sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_1.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.4sd_pre_lower <- as.data.frame(cbind(PC2_1.4sd_pre$x, PC2_1.4sd_pre$lower))
names(PC2_1.4sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.4sd_pre_lower$age, PC2_1.4sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_1.4sd_pre_lower, indices) {
  d <- PC2_1.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_1.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_1.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_1.4sd_pre_lower_50)[1] <- "PC2"

PC2_1.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_1.4sd_pre_lower_all)[1:ncol(PC2_1.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_1.4sd_pre_lower_all$group <- seq(1,(nrow(PC2_1.4sd_pre_lower_all)), 1)

PC2_1.4sd_pre_lower_all <- reshape2::melt(PC2_1.4sd_pre_lower_all, id.vars = "group")
names(PC2_1.4sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_1.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.2sd_pre_upper <- as.data.frame(cbind(PC2_1.2sd_pre$x, PC2_1.2sd_pre$upper))
names(PC2_1.2sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.2sd_pre_upper$age, PC2_1.2sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_1.2sd_pre_upper, indices) {
  d <- PC2_1.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_1.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_1.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_1.2sd_pre_upper_50)[1] <- "PC2"

PC2_1.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_1.2sd_pre_upper_all)[1:ncol(PC2_1.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_1.2sd_pre_upper_all$group <- seq(1,(nrow(PC2_1.2sd_pre_upper_all)), 1)

PC2_1.2sd_pre_upper_all <- reshape2::melt(PC2_1.2sd_pre_upper_all, id.vars = "group")
names(PC2_1.2sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_1.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.2sd_pre_lower <- as.data.frame(cbind(PC2_1.2sd_pre$x, PC2_1.2sd_pre$lower))
names(PC2_1.2sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.2sd_pre_lower$age, PC2_1.2sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_1.2sd_pre_lower, indices) {
  d <- PC2_1.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_1.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_1.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_1.2sd_pre_lower_50)[1] <- "PC2"

PC2_1.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_1.2sd_pre_lower_all)[1:ncol(PC2_1.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_1.2sd_pre_lower_all$group <- seq(1,(nrow(PC2_1.2sd_pre_lower_all)), 1)

PC2_1.2sd_pre_lower_all <- reshape2::melt(PC2_1.2sd_pre_lower_all, id.vars = "group")
names(PC2_1.2sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_1.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.0sd_pre_upper <- as.data.frame(cbind(PC2_1.0sd_pre$x, PC2_1.0sd_pre$upper))
names(PC2_1.0sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.0sd_pre_upper$age, PC2_1.0sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_1.0sd_pre_upper, indices) {
  d <- PC2_1.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_1.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_1.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_1.0sd_pre_upper_50)[1] <- "PC2"

PC2_1.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_1.0sd_pre_upper_all)[1:ncol(PC2_1.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_1.0sd_pre_upper_all$group <- seq(1,(nrow(PC2_1.0sd_pre_upper_all)), 1)

PC2_1.0sd_pre_upper_all <- reshape2::melt(PC2_1.0sd_pre_upper_all, id.vars = "group")
names(PC2_1.0sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_1.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_1.0sd_pre_lower <- as.data.frame(cbind(PC2_1.0sd_pre$x, PC2_1.0sd_pre$lower))
names(PC2_1.0sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_1.0sd_pre_lower$age, PC2_1.0sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_1.0sd_pre_lower, indices) {
  d <- PC2_1.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_1.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_1.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_1.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_1.0sd_pre_lower_50)[1] <- "PC2"

PC2_1.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_1.0sd_pre_lower_all)[1:ncol(PC2_1.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_1.0sd_pre_lower_all$group <- seq(1,(nrow(PC2_1.0sd_pre_lower_all)), 1)

PC2_1.0sd_pre_lower_all <- reshape2::melt(PC2_1.0sd_pre_lower_all, id.vars = "group")
names(PC2_1.0sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_1.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.8sd_pre_upper <- as.data.frame(cbind(PC2_0.8sd_pre$x, PC2_0.8sd_pre$upper))
names(PC2_0.8sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.8sd_pre_upper$age, PC2_0.8sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_0.8sd_pre_upper, indices) {
  d <- PC2_0.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_0.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_0.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_0.8sd_pre_upper_50)[1] <- "PC2"

PC2_0.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_0.8sd_pre_upper_all)[1:ncol(PC2_0.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_0.8sd_pre_upper_all$group <- seq(1,(nrow(PC2_0.8sd_pre_upper_all)), 1)

PC2_0.8sd_pre_upper_all <- reshape2::melt(PC2_0.8sd_pre_upper_all, id.vars = "group")
names(PC2_0.8sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_0.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.8sd_pre_lower <- as.data.frame(cbind(PC2_0.8sd_pre$x, PC2_0.8sd_pre$lower))
names(PC2_0.8sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.8sd_pre_lower$age, PC2_0.8sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_0.8sd_pre_lower, indices) {
  d <- PC2_0.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_0.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_0.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_0.8sd_pre_lower_50)[1] <- "PC2"

PC2_0.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_0.8sd_pre_lower_all)[1:ncol(PC2_0.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_0.8sd_pre_lower_all$group <- seq(1,(nrow(PC2_0.8sd_pre_lower_all)), 1)

PC2_0.8sd_pre_lower_all <- reshape2::melt(PC2_0.8sd_pre_lower_all, id.vars = "group")
names(PC2_0.8sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_0.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.6sd_pre_upper <- as.data.frame(cbind(PC2_0.6sd_pre$x, PC2_0.6sd_pre$upper))
names(PC2_0.6sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.6sd_pre_upper$age, PC2_0.6sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_0.6sd_pre_upper, indices) {
  d <- PC2_0.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_0.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_0.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_0.6sd_pre_upper_50)[1] <- "PC2"

PC2_0.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_0.6sd_pre_upper_all)[1:ncol(PC2_0.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_0.6sd_pre_upper_all$group <- seq(1,(nrow(PC2_0.6sd_pre_upper_all)), 1)

PC2_0.6sd_pre_upper_all <- reshape2::melt(PC2_0.6sd_pre_upper_all, id.vars = "group")
names(PC2_0.6sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_0.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.6sd_pre_lower <- as.data.frame(cbind(PC2_0.6sd_pre$x, PC2_0.6sd_pre$lower))
names(PC2_0.6sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.6sd_pre_lower$age, PC2_0.6sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_0.6sd_pre_lower, indices) {
  d <- PC2_0.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_0.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_0.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_0.6sd_pre_lower_50)[1] <- "PC2"

PC2_0.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_0.6sd_pre_lower_all)[1:ncol(PC2_0.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_0.6sd_pre_lower_all$group <- seq(1,(nrow(PC2_0.6sd_pre_lower_all)), 1)

PC2_0.6sd_pre_lower_all <- reshape2::melt(PC2_0.6sd_pre_lower_all, id.vars = "group")
names(PC2_0.6sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_0.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.4sd_pre_upper <- as.data.frame(cbind(PC2_0.4sd_pre$x, PC2_0.4sd_pre$upper))
names(PC2_0.4sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.4sd_pre_upper$age, PC2_0.4sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_0.4sd_pre_upper, indices) {
  d <- PC2_0.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_0.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_0.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_0.4sd_pre_upper_50)[1] <- "PC2"

PC2_0.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_0.4sd_pre_upper_all)[1:ncol(PC2_0.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_0.4sd_pre_upper_all$group <- seq(1,(nrow(PC2_0.4sd_pre_upper_all)), 1)

PC2_0.4sd_pre_upper_all <- reshape2::melt(PC2_0.4sd_pre_upper_all, id.vars = "group")
names(PC2_0.4sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_0.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.4sd_pre_lower <- as.data.frame(cbind(PC2_0.4sd_pre$x, PC2_0.4sd_pre$lower))
names(PC2_0.4sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.4sd_pre_lower$age, PC2_0.4sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_0.4sd_pre_lower, indices) {
  d <- PC2_0.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_0.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_0.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_0.4sd_pre_lower_50)[1] <- "PC2"

PC2_0.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_0.4sd_pre_lower_all)[1:ncol(PC2_0.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_0.4sd_pre_lower_all$group <- seq(1,(nrow(PC2_0.4sd_pre_lower_all)), 1)

PC2_0.4sd_pre_lower_all <- reshape2::melt(PC2_0.4sd_pre_lower_all, id.vars = "group")
names(PC2_0.4sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_0.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.2sd_pre_upper <- as.data.frame(cbind(PC2_0.2sd_pre$x, PC2_0.2sd_pre$upper))
names(PC2_0.2sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.2sd_pre_upper$age, PC2_0.2sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_0.2sd_pre_upper, indices) {
  d <- PC2_0.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_0.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_0.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_0.2sd_pre_upper_50)[1] <- "PC2"

PC2_0.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_0.2sd_pre_upper_all)[1:ncol(PC2_0.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_0.2sd_pre_upper_all$group <- seq(1,(nrow(PC2_0.2sd_pre_upper_all)), 1)

PC2_0.2sd_pre_upper_all <- reshape2::melt(PC2_0.2sd_pre_upper_all, id.vars = "group")
names(PC2_0.2sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_0.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.2sd_pre_lower <- as.data.frame(cbind(PC2_0.2sd_pre$x, PC2_0.2sd_pre$lower))
names(PC2_0.2sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.2sd_pre_lower$age, PC2_0.2sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_0.2sd_pre_lower, indices) {
  d <- PC2_0.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_0.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_0.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_0.2sd_pre_lower_50)[1] <- "PC2"

PC2_0.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_0.2sd_pre_lower_all)[1:ncol(PC2_0.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_0.2sd_pre_lower_all$group <- seq(1,(nrow(PC2_0.2sd_pre_lower_all)), 1)

PC2_0.2sd_pre_lower_all <- reshape2::melt(PC2_0.2sd_pre_lower_all, id.vars = "group")
names(PC2_0.2sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_0.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC2 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.0sd_pre_upper <- as.data.frame(cbind(PC2_0.0sd_pre$x, PC2_0.0sd_pre$upper))
names(PC2_0.0sd_pre_upper)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.0sd_pre_upper$age, PC2_0.0sd_pre_upper$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC2_0.0sd_pre_upper, indices) {
  d <- PC2_0.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC2_0.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC2_0.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC2_0.0sd_pre_upper_50)[1] <- "PC2"

PC2_0.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC2_0.0sd_pre_upper_all)[1:ncol(PC2_0.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC2_0.0sd_pre_upper_all$group <- seq(1,(nrow(PC2_0.0sd_pre_upper_all)), 1)

PC2_0.0sd_pre_upper_all <- reshape2::melt(PC2_0.0sd_pre_upper_all, id.vars = "group")
names(PC2_0.0sd_pre_upper_all)[2:3] <- c("age", "PC2")

PC2_0.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC2 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC2_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC2, degree=1, nsigma = nsig, span = span)
PC2_0.0sd_pre_lower <- as.data.frame(cbind(PC2_0.0sd_pre$x, PC2_0.0sd_pre$lower))
names(PC2_0.0sd_pre_lower)[1:2] <- c("age","PC2")

# CV spans
loess.predict <- loess.as(PC2_0.0sd_pre_lower$age, PC2_0.0sd_pre_lower$PC2, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC2_0.0sd_pre_lower, indices) {
  d <- PC2_0.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC2 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC2_0.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC2_0.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC2_0.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC2_0.0sd_pre_lower_50)[1] <- "PC2"

PC2_0.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC2_0.0sd_pre_lower_all)[1:ncol(PC2_0.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC2_0.0sd_pre_lower_all$group <- seq(1,(nrow(PC2_0.0sd_pre_lower_all)), 1)

PC2_0.0sd_pre_lower_all <- reshape2::melt(PC2_0.0sd_pre_lower_all, id.vars = "group")
names(PC2_0.0sd_pre_lower_all)[2:3] <- c("age", "PC2")

PC2_0.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

b <- ggplot() + 
  #geom_line(data = PC2_2sd_pre_upper_all, aes(age2,PC2, group = group)) +
  #geom_line(data = PC2_2sd_pre_lower_all, aes(age2,PC2, group = group)) +
  geom_line(data = PC2_2sd_pre_upper_50, aes(age,PC2), colour = 'blue') +
  geom_line(data = PC2_2sd_pre_lower_50, aes(age,PC2), colour = 'blue') +
  geom_line(data = PC2_1.6sd_pre_upper_50, aes(age,PC2), colour = 'green') +
  geom_line(data = PC2_1.6sd_pre_lower_50, aes(age,PC2), colour = 'green') +
  geom_line(data = PC2_1.4sd_pre_upper_50, aes(age,PC2), colour = 'orange') +
  geom_line(data = PC2_1.4sd_pre_lower_50, aes(age,PC2), colour = 'orange') +
  geom_line(data = PC2_1.2sd_pre_upper_50, aes(age,PC2), colour = 'pink') +
  geom_line(data = PC2_1.2sd_pre_lower_50, aes(age,PC2), colour = 'pink') +
  geom_line(data = PC2_1.0sd_pre_upper_50, aes(age,PC2), colour = 'brown') +
  geom_line(data = PC2_1.0sd_pre_lower_50, aes(age,PC2), colour = 'brown') +
  geom_line(data = PC2_0.8sd_pre_upper_50, aes(age,PC2), colour = 'yellow') +
  geom_line(data = PC2_0.8sd_pre_lower_50, aes(age,PC2), colour = 'yellow') +
  geom_line(data = PC2_0.6sd_pre_upper_50, aes(age,PC2), colour = 'red') +
  geom_line(data = PC2_0.6sd_pre_lower_50, aes(age,PC2), colour = 'red') +
  geom_line(data = PC2_0.4sd_pre_upper_50, aes(age,PC2), colour = 'blue') +
  geom_line(data = PC2_0.4sd_pre_lower_50, aes(age,PC2), colour = 'blue') +
  geom_line(data = PC2_0.2sd_pre_upper_50, aes(age,PC2), colour = 'green') +
  geom_line(data = PC2_0.2sd_pre_lower_50, aes(age,PC2), colour = 'green') +
  geom_line(data = PC2_0.0sd_pre_upper_50, aes(age,PC2), colour = 'orange') +
  geom_line(data = PC2_0.0sd_pre_lower_50, aes(age,PC2), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC2)) + 
  #geom_point(data = PC2_2sd_pre_upper, aes(age, PC2), colour = 'red') +
  scale_x_reverse(limits = c(4000,540)) + theme_bw()

grid.arrange(b,a,ncol=2)




##### Phanerozoic PC3 #####

coords.all_Phan <- subset(coords.all2, age < 541)

# 2 sigma

nsig <- 2

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_2sd_phan_upper <- as.data.frame(cbind(PC3_2sd_phan$x, PC3_2sd_phan$upper))
names(PC3_2sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_2sd_phan_upper$age, PC3_2sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_2sd_phan_upper, indices) {
  d <- PC3_2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_2sd_phan_upper_50)[1] <- "PC3"

PC3_2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_2sd_phan_upper_all)[1:ncol(PC3_2sd_phan_upper_all)] <- seq(0,540, 1)
PC3_2sd_phan_upper_all$group <- seq(1,(nrow(PC3_2sd_phan_upper_all)), 1)

PC3_2sd_phan_upper_all <- reshape2::melt(PC3_2sd_phan_upper_all, id.vars = "group")
names(PC3_2sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_2sd_phan_lower <- as.data.frame(cbind(PC3_2sd_phan$x, PC3_2sd_phan$lower))
names(PC3_2sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_2sd_phan_lower$age, PC3_2sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_2sd_phan_lower, indices) {
  d <- PC3_2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_2sd_phan_lower_50)[1] <- "PC3"

PC3_2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_2sd_phan_lower_all)[1:ncol(PC3_2sd_phan_lower_all)] <- seq(0,540, 1)
PC3_2sd_phan_lower_all$group <- seq(1,(nrow(PC3_2sd_phan_lower_all)), 1)

PC3_2sd_phan_lower_all <- reshape2::melt(PC3_2sd_phan_lower_all, id.vars = "group")
names(PC3_2sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.8sd_phan_upper <- as.data.frame(cbind(PC3_1.8sd_phan$x, PC3_1.8sd_phan$upper))
names(PC3_1.8sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.8sd_phan_upper$age, PC3_1.8sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_1.8sd_phan_upper, indices) {
  d <- PC3_1.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_1.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_1.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_1.8sd_phan_upper_50)[1] <- "PC3"

PC3_1.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_1.8sd_phan_upper_all)[1:ncol(PC3_1.8sd_phan_upper_all)] <- seq(0,540, 1)
PC3_1.8sd_phan_upper_all$group <- seq(1,(nrow(PC3_1.8sd_phan_upper_all)), 1)

PC3_1.8sd_phan_upper_all <- reshape2::melt(PC3_1.8sd_phan_upper_all, id.vars = "group")
names(PC3_1.8sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_1.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.8sd_phan_lower <- as.data.frame(cbind(PC3_1.8sd_phan$x, PC3_1.8sd_phan$lower))
names(PC3_1.8sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.8sd_phan_lower$age, PC3_1.8sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_1.8sd_phan_lower, indices) {
  d <- PC3_1.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_1.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_1.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_1.8sd_phan_lower_50)[1] <- "PC3"

PC3_1.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_1.8sd_phan_lower_all)[1:ncol(PC3_1.8sd_phan_lower_all)] <- seq(0,540, 1)
PC3_1.8sd_phan_lower_all$group <- seq(1,(nrow(PC3_1.8sd_phan_lower_all)), 1)

PC3_1.8sd_phan_lower_all <- reshape2::melt(PC3_1.8sd_phan_lower_all, id.vars = "group")
names(PC3_1.8sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_1.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.6sd_phan_upper <- as.data.frame(cbind(PC3_1.6sd_phan$x, PC3_1.6sd_phan$upper))
names(PC3_1.6sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.6sd_phan_upper$age, PC3_1.6sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_1.6sd_phan_upper, indices) {
  d <- PC3_1.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_1.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_1.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_1.6sd_phan_upper_50)[1] <- "PC3"

PC3_1.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_1.6sd_phan_upper_all)[1:ncol(PC3_1.6sd_phan_upper_all)] <- seq(0,540, 1)
PC3_1.6sd_phan_upper_all$group <- seq(1,(nrow(PC3_1.6sd_phan_upper_all)), 1)

PC3_1.6sd_phan_upper_all <- reshape2::melt(PC3_1.6sd_phan_upper_all, id.vars = "group")
names(PC3_1.6sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_1.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.6sd_phan_lower <- as.data.frame(cbind(PC3_1.6sd_phan$x, PC3_1.6sd_phan$lower))
names(PC3_1.6sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.6sd_phan_lower$age, PC3_1.6sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_1.6sd_phan_lower, indices) {
  d <- PC3_1.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_1.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_1.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_1.6sd_phan_lower_50)[1] <- "PC3"

PC3_1.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_1.6sd_phan_lower_all)[1:ncol(PC3_1.6sd_phan_lower_all)] <- seq(0,540, 1)
PC3_1.6sd_phan_lower_all$group <- seq(1,(nrow(PC3_1.6sd_phan_lower_all)), 1)

PC3_1.6sd_phan_lower_all <- reshape2::melt(PC3_1.6sd_phan_lower_all, id.vars = "group")
names(PC3_1.6sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_1.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.4sd_phan_upper <- as.data.frame(cbind(PC3_1.4sd_phan$x, PC3_1.4sd_phan$upper))
names(PC3_1.4sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.4sd_phan_upper$age, PC3_1.4sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_1.4sd_phan_upper, indices) {
  d <- PC3_1.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_1.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_1.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_1.4sd_phan_upper_50)[1] <- "PC3"

PC3_1.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_1.4sd_phan_upper_all)[1:ncol(PC3_1.4sd_phan_upper_all)] <- seq(0,540, 1)
PC3_1.4sd_phan_upper_all$group <- seq(1,(nrow(PC3_1.4sd_phan_upper_all)), 1)

PC3_1.4sd_phan_upper_all <- reshape2::melt(PC3_1.4sd_phan_upper_all, id.vars = "group")
names(PC3_1.4sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_1.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.4sd_phan_lower <- as.data.frame(cbind(PC3_1.4sd_phan$x, PC3_1.4sd_phan$lower))
names(PC3_1.4sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.4sd_phan_lower$age, PC3_1.4sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_1.4sd_phan_lower, indices) {
  d <- PC3_1.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_1.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_1.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_1.4sd_phan_lower_50)[1] <- "PC3"

PC3_1.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_1.4sd_phan_lower_all)[1:ncol(PC3_1.4sd_phan_lower_all)] <- seq(0,540, 1)
PC3_1.4sd_phan_lower_all$group <- seq(1,(nrow(PC3_1.4sd_phan_lower_all)), 1)

PC3_1.4sd_phan_lower_all <- reshape2::melt(PC3_1.4sd_phan_lower_all, id.vars = "group")
names(PC3_1.4sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_1.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.2sd_phan_upper <- as.data.frame(cbind(PC3_1.2sd_phan$x, PC3_1.2sd_phan$upper))
names(PC3_1.2sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.2sd_phan_upper$age, PC3_1.2sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_1.2sd_phan_upper, indices) {
  d <- PC3_1.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_1.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_1.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_1.2sd_phan_upper_50)[1] <- "PC3"

PC3_1.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_1.2sd_phan_upper_all)[1:ncol(PC3_1.2sd_phan_upper_all)] <- seq(0,540, 1)
PC3_1.2sd_phan_upper_all$group <- seq(1,(nrow(PC3_1.2sd_phan_upper_all)), 1)

PC3_1.2sd_phan_upper_all <- reshape2::melt(PC3_1.2sd_phan_upper_all, id.vars = "group")
names(PC3_1.2sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_1.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.2sd_phan_lower <- as.data.frame(cbind(PC3_1.2sd_phan$x, PC3_1.2sd_phan$lower))
names(PC3_1.2sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.2sd_phan_lower$age, PC3_1.2sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_1.2sd_phan_lower, indices) {
  d <- PC3_1.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_1.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_1.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_1.2sd_phan_lower_50)[1] <- "PC3"

PC3_1.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_1.2sd_phan_lower_all)[1:ncol(PC3_1.2sd_phan_lower_all)] <- seq(0,540, 1)
PC3_1.2sd_phan_lower_all$group <- seq(1,(nrow(PC3_1.2sd_phan_lower_all)), 1)

PC3_1.2sd_phan_lower_all <- reshape2::melt(PC3_1.2sd_phan_lower_all, id.vars = "group")
names(PC3_1.2sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_1.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.0sd_phan_upper <- as.data.frame(cbind(PC3_1.0sd_phan$x, PC3_1.0sd_phan$upper))
names(PC3_1.0sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.0sd_phan_upper$age, PC3_1.0sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_1.0sd_phan_upper, indices) {
  d <- PC3_1.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_1.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_1.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_1.0sd_phan_upper_50)[1] <- "PC3"

PC3_1.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_1.0sd_phan_upper_all)[1:ncol(PC3_1.0sd_phan_upper_all)] <- seq(0,540, 1)
PC3_1.0sd_phan_upper_all$group <- seq(1,(nrow(PC3_1.0sd_phan_upper_all)), 1)

PC3_1.0sd_phan_upper_all <- reshape2::melt(PC3_1.0sd_phan_upper_all, id.vars = "group")
names(PC3_1.0sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_1.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.0sd_phan_lower <- as.data.frame(cbind(PC3_1.0sd_phan$x, PC3_1.0sd_phan$lower))
names(PC3_1.0sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.0sd_phan_lower$age, PC3_1.0sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_1.0sd_phan_lower, indices) {
  d <- PC3_1.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_1.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_1.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_1.0sd_phan_lower_50)[1] <- "PC3"

PC3_1.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_1.0sd_phan_lower_all)[1:ncol(PC3_1.0sd_phan_lower_all)] <- seq(0,540, 1)
PC3_1.0sd_phan_lower_all$group <- seq(1,(nrow(PC3_1.0sd_phan_lower_all)), 1)

PC3_1.0sd_phan_lower_all <- reshape2::melt(PC3_1.0sd_phan_lower_all, id.vars = "group")
names(PC3_1.0sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_1.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.8sd_phan_upper <- as.data.frame(cbind(PC3_0.8sd_phan$x, PC3_0.8sd_phan$upper))
names(PC3_0.8sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.8sd_phan_upper$age, PC3_0.8sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_0.8sd_phan_upper, indices) {
  d <- PC3_0.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_0.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_0.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_0.8sd_phan_upper_50)[1] <- "PC3"

PC3_0.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_0.8sd_phan_upper_all)[1:ncol(PC3_0.8sd_phan_upper_all)] <- seq(0,540, 1)
PC3_0.8sd_phan_upper_all$group <- seq(1,(nrow(PC3_0.8sd_phan_upper_all)), 1)

PC3_0.8sd_phan_upper_all <- reshape2::melt(PC3_0.8sd_phan_upper_all, id.vars = "group")
names(PC3_0.8sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_0.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.8sd_phan_lower <- as.data.frame(cbind(PC3_0.8sd_phan$x, PC3_0.8sd_phan$lower))
names(PC3_0.8sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.8sd_phan_lower$age, PC3_0.8sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_0.8sd_phan_lower, indices) {
  d <- PC3_0.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_0.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_0.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_0.8sd_phan_lower_50)[1] <- "PC3"

PC3_0.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_0.8sd_phan_lower_all)[1:ncol(PC3_0.8sd_phan_lower_all)] <- seq(0,540, 1)
PC3_0.8sd_phan_lower_all$group <- seq(1,(nrow(PC3_0.8sd_phan_lower_all)), 1)

PC3_0.8sd_phan_lower_all <- reshape2::melt(PC3_0.8sd_phan_lower_all, id.vars = "group")
names(PC3_0.8sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_0.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.6sd_phan_upper <- as.data.frame(cbind(PC3_0.6sd_phan$x, PC3_0.6sd_phan$upper))
names(PC3_0.6sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.6sd_phan_upper$age, PC3_0.6sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_0.6sd_phan_upper, indices) {
  d <- PC3_0.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_0.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_0.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_0.6sd_phan_upper_50)[1] <- "PC3"

PC3_0.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_0.6sd_phan_upper_all)[1:ncol(PC3_0.6sd_phan_upper_all)] <- seq(0,540, 1)
PC3_0.6sd_phan_upper_all$group <- seq(1,(nrow(PC3_0.6sd_phan_upper_all)), 1)

PC3_0.6sd_phan_upper_all <- reshape2::melt(PC3_0.6sd_phan_upper_all, id.vars = "group")
names(PC3_0.6sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_0.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.6sd_phan_lower <- as.data.frame(cbind(PC3_0.6sd_phan$x, PC3_0.6sd_phan$lower))
names(PC3_0.6sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.6sd_phan_lower$age, PC3_0.6sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_0.6sd_phan_lower, indices) {
  d <- PC3_0.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_0.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_0.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_0.6sd_phan_lower_50)[1] <- "PC3"

PC3_0.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_0.6sd_phan_lower_all)[1:ncol(PC3_0.6sd_phan_lower_all)] <- seq(0,540, 1)
PC3_0.6sd_phan_lower_all$group <- seq(1,(nrow(PC3_0.6sd_phan_lower_all)), 1)

PC3_0.6sd_phan_lower_all <- reshape2::melt(PC3_0.6sd_phan_lower_all, id.vars = "group")
names(PC3_0.6sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_0.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.4sd_phan_upper <- as.data.frame(cbind(PC3_0.4sd_phan$x, PC3_0.4sd_phan$upper))
names(PC3_0.4sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.4sd_phan_upper$age, PC3_0.4sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_0.4sd_phan_upper, indices) {
  d <- PC3_0.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_0.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_0.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_0.4sd_phan_upper_50)[1] <- "PC3"

PC3_0.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_0.4sd_phan_upper_all)[1:ncol(PC3_0.4sd_phan_upper_all)] <- seq(0,540, 1)
PC3_0.4sd_phan_upper_all$group <- seq(1,(nrow(PC3_0.4sd_phan_upper_all)), 1)

PC3_0.4sd_phan_upper_all <- reshape2::melt(PC3_0.4sd_phan_upper_all, id.vars = "group")
names(PC3_0.4sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_0.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.4sd_phan_lower <- as.data.frame(cbind(PC3_0.4sd_phan$x, PC3_0.4sd_phan$lower))
names(PC3_0.4sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.4sd_phan_lower$age, PC3_0.4sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_0.4sd_phan_lower, indices) {
  d <- PC3_0.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_0.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_0.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_0.4sd_phan_lower_50)[1] <- "PC3"

PC3_0.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_0.4sd_phan_lower_all)[1:ncol(PC3_0.4sd_phan_lower_all)] <- seq(0,540, 1)
PC3_0.4sd_phan_lower_all$group <- seq(1,(nrow(PC3_0.4sd_phan_lower_all)), 1)

PC3_0.4sd_phan_lower_all <- reshape2::melt(PC3_0.4sd_phan_lower_all, id.vars = "group")
names(PC3_0.4sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_0.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.2sd_phan_upper <- as.data.frame(cbind(PC3_0.2sd_phan$x, PC3_0.2sd_phan$upper))
names(PC3_0.2sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.2sd_phan_upper$age, PC3_0.2sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_0.2sd_phan_upper, indices) {
  d <- PC3_0.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_0.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_0.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_0.2sd_phan_upper_50)[1] <- "PC3"

PC3_0.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_0.2sd_phan_upper_all)[1:ncol(PC3_0.2sd_phan_upper_all)] <- seq(0,540, 1)
PC3_0.2sd_phan_upper_all$group <- seq(1,(nrow(PC3_0.2sd_phan_upper_all)), 1)

PC3_0.2sd_phan_upper_all <- reshape2::melt(PC3_0.2sd_phan_upper_all, id.vars = "group")
names(PC3_0.2sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_0.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.2sd_phan_lower <- as.data.frame(cbind(PC3_0.2sd_phan$x, PC3_0.2sd_phan$lower))
names(PC3_0.2sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.2sd_phan_lower$age, PC3_0.2sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_0.2sd_phan_lower, indices) {
  d <- PC3_0.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_0.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_0.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_0.2sd_phan_lower_50)[1] <- "PC3"

PC3_0.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_0.2sd_phan_lower_all)[1:ncol(PC3_0.2sd_phan_lower_all)] <- seq(0,540, 1)
PC3_0.2sd_phan_lower_all$group <- seq(1,(nrow(PC3_0.2sd_phan_lower_all)), 1)

PC3_0.2sd_phan_lower_all <- reshape2::melt(PC3_0.2sd_phan_lower_all, id.vars = "group")
names(PC3_0.2sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_0.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.0sd_phan_upper <- as.data.frame(cbind(PC3_0.0sd_phan$x, PC3_0.0sd_phan$upper))
names(PC3_0.0sd_phan_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.0sd_phan_upper$age, PC3_0.0sd_phan_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_0.0sd_phan_upper, indices) {
  d <- PC3_0.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_0.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_0.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC3_0.0sd_phan_upper_50)[1] <- "PC3"

PC3_0.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC3_0.0sd_phan_upper_all)[1:ncol(PC3_0.0sd_phan_upper_all)] <- seq(0,540, 1)
PC3_0.0sd_phan_upper_all$group <- seq(1,(nrow(PC3_0.0sd_phan_upper_all)), 1)

PC3_0.0sd_phan_upper_all <- reshape2::melt(PC3_0.0sd_phan_upper_all, id.vars = "group")
names(PC3_0.0sd_phan_upper_all)[2:3] <- c("age", "PC3")

PC3_0.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.0sd_phan_lower <- as.data.frame(cbind(PC3_0.0sd_phan$x, PC3_0.0sd_phan$lower))
names(PC3_0.0sd_phan_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.0sd_phan_lower$age, PC3_0.0sd_phan_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_0.0sd_phan_lower, indices) {
  d <- PC3_0.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_0.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_0.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC3_0.0sd_phan_lower_50)[1] <- "PC3"

PC3_0.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC3_0.0sd_phan_lower_all)[1:ncol(PC3_0.0sd_phan_lower_all)] <- seq(0,540, 1)
PC3_0.0sd_phan_lower_all$group <- seq(1,(nrow(PC3_0.0sd_phan_lower_all)), 1)

PC3_0.0sd_phan_lower_all <- reshape2::melt(PC3_0.0sd_phan_lower_all, id.vars = "group")
names(PC3_0.0sd_phan_lower_all)[2:3] <- c("age", "PC3")

PC3_0.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

a <- ggplot() + 
  #geom_line(data = PC3_2sd_phan_upper_all, aes(age2,PC3, group = group)) +
  #geom_line(data = PC3_2sd_phan_lower_all, aes(age2,PC3, group = group)) +
  geom_line(data = PC3_2sd_phan_upper_50, aes(age,PC3), colour = 'blue') +
  geom_line(data = PC3_2sd_phan_lower_50, aes(age,PC3), colour = 'blue') +
  geom_line(data = PC3_1.6sd_phan_upper_50, aes(age,PC3), colour = 'green') +
  geom_line(data = PC3_1.6sd_phan_lower_50, aes(age,PC3), colour = 'green') +
  geom_line(data = PC3_1.4sd_phan_upper_50, aes(age,PC3), colour = 'orange') +
  geom_line(data = PC3_1.4sd_phan_lower_50, aes(age,PC3), colour = 'orange') +
  geom_line(data = PC3_1.2sd_phan_upper_50, aes(age,PC3), colour = 'pink') +
  geom_line(data = PC3_1.2sd_phan_lower_50, aes(age,PC3), colour = 'pink') +
  geom_line(data = PC3_1.0sd_phan_upper_50, aes(age,PC3), colour = 'brown') +
  geom_line(data = PC3_1.0sd_phan_lower_50, aes(age,PC3), colour = 'brown') +
  geom_line(data = PC3_0.8sd_phan_upper_50, aes(age,PC3), colour = 'yellow') +
  geom_line(data = PC3_0.8sd_phan_lower_50, aes(age,PC3), colour = 'yellow') +
  geom_line(data = PC3_0.6sd_phan_upper_50, aes(age,PC3), colour = 'red') +
  geom_line(data = PC3_0.6sd_phan_lower_50, aes(age,PC3), colour = 'red') +
  geom_line(data = PC3_0.4sd_phan_upper_50, aes(age,PC3), colour = 'blue') +
  geom_line(data = PC3_0.4sd_phan_lower_50, aes(age,PC3), colour = 'blue') +
  geom_line(data = PC3_0.2sd_phan_upper_50, aes(age,PC3), colour = 'green') +
  geom_line(data = PC3_0.2sd_phan_lower_50, aes(age,PC3), colour = 'green') +
  geom_line(data = PC3_0.0sd_phan_upper_50, aes(age,PC3), colour = 'orange') +
  geom_line(data = PC3_0.0sd_phan_lower_50, aes(age,PC3), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC3)) + 
  #geom_point(data = PC3_2sd_phan_upper, aes(age, PC3), colour = 'red') +
  scale_x_reverse(limits = c(540,0)) + theme_bw()

##### Precambrian PC3 #####

coords.all_pre <- subset(coords.all2, age > 541)

# 2 sigma

nsig <- 2

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_2sd_pre_upper <- as.data.frame(cbind(PC3_2sd_pre$x, PC3_2sd_pre$upper))
names(PC3_2sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_2sd_pre_upper$age, PC3_2sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[1], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_2sd_pre_upper, indices) {
  d <- PC3_2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_2sd_pre_upper_50)[1] <- "PC3"

PC3_2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_2sd_pre_upper_all)[1:ncol(PC3_2sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_2sd_pre_upper_all$group <- seq(1,(nrow(PC3_2sd_pre_upper_all)), 1)

PC3_2sd_pre_upper_all <- reshape2::melt(PC3_2sd_pre_upper_all, id.vars = "group")
names(PC3_2sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_2sd_pre_lower <- as.data.frame(cbind(PC3_2sd_pre$x, PC3_2sd_pre$lower))
names(PC3_2sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_2sd_pre_lower$age, PC3_2sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_2sd_pre_lower, indices) {
  d <- PC3_2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_2sd_pre_lower_50)[1] <- "PC3"

PC3_2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_2sd_pre_lower_all)[1:ncol(PC3_2sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_2sd_pre_lower_all$group <- seq(1,(nrow(PC3_2sd_pre_lower_all)), 1)

PC3_2sd_pre_lower_all <- reshape2::melt(PC3_2sd_pre_lower_all, id.vars = "group")
names(PC3_2sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.8sd_pre_upper <- as.data.frame(cbind(PC3_1.8sd_pre$x, PC3_1.8sd_pre$upper))
names(PC3_1.8sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.8sd_pre_upper$age, PC3_1.8sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_1.8sd_pre_upper, indices) {
  d <- PC3_1.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_1.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_1.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_1.8sd_pre_upper_50)[1] <- "PC3"

PC3_1.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_1.8sd_pre_upper_all)[1:ncol(PC3_1.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_1.8sd_pre_upper_all$group <- seq(1,(nrow(PC3_1.8sd_pre_upper_all)), 1)

PC3_1.8sd_pre_upper_all <- reshape2::melt(PC3_1.8sd_pre_upper_all, id.vars = "group")
names(PC3_1.8sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_1.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.8sd_pre_lower <- as.data.frame(cbind(PC3_1.8sd_pre$x, PC3_1.8sd_pre$lower))
names(PC3_1.8sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.8sd_pre_lower$age, PC3_1.8sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_1.8sd_pre_lower, indices) {
  d <- PC3_1.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_1.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_1.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_1.8sd_pre_lower_50)[1] <- "PC3"

PC3_1.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_1.8sd_pre_lower_all)[1:ncol(PC3_1.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_1.8sd_pre_lower_all$group <- seq(1,(nrow(PC3_1.8sd_pre_lower_all)), 1)

PC3_1.8sd_pre_lower_all <- reshape2::melt(PC3_1.8sd_pre_lower_all, id.vars = "group")
names(PC3_1.8sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_1.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.6sd_pre_upper <- as.data.frame(cbind(PC3_1.6sd_pre$x, PC3_1.6sd_pre$upper))
names(PC3_1.6sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.6sd_pre_upper$age, PC3_1.6sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_1.6sd_pre_upper, indices) {
  d <- PC3_1.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_1.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_1.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_1.6sd_pre_upper_50)[1] <- "PC3"

PC3_1.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_1.6sd_pre_upper_all)[1:ncol(PC3_1.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_1.6sd_pre_upper_all$group <- seq(1,(nrow(PC3_1.6sd_pre_upper_all)), 1)

PC3_1.6sd_pre_upper_all <- reshape2::melt(PC3_1.6sd_pre_upper_all, id.vars = "group")
names(PC3_1.6sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_1.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.6sd_pre_lower <- as.data.frame(cbind(PC3_1.6sd_pre$x, PC3_1.6sd_pre$lower))
names(PC3_1.6sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.6sd_pre_lower$age, PC3_1.6sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_1.6sd_pre_lower, indices) {
  d <- PC3_1.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_1.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_1.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_1.6sd_pre_lower_50)[1] <- "PC3"

PC3_1.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_1.6sd_pre_lower_all)[1:ncol(PC3_1.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_1.6sd_pre_lower_all$group <- seq(1,(nrow(PC3_1.6sd_pre_lower_all)), 1)

PC3_1.6sd_pre_lower_all <- reshape2::melt(PC3_1.6sd_pre_lower_all, id.vars = "group")
names(PC3_1.6sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_1.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.4sd_pre_upper <- as.data.frame(cbind(PC3_1.4sd_pre$x, PC3_1.4sd_pre$upper))
names(PC3_1.4sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.4sd_pre_upper$age, PC3_1.4sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_1.4sd_pre_upper, indices) {
  d <- PC3_1.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_1.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_1.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_1.4sd_pre_upper_50)[1] <- "PC3"

PC3_1.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_1.4sd_pre_upper_all)[1:ncol(PC3_1.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_1.4sd_pre_upper_all$group <- seq(1,(nrow(PC3_1.4sd_pre_upper_all)), 1)

PC3_1.4sd_pre_upper_all <- reshape2::melt(PC3_1.4sd_pre_upper_all, id.vars = "group")
names(PC3_1.4sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_1.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.4sd_pre_lower <- as.data.frame(cbind(PC3_1.4sd_pre$x, PC3_1.4sd_pre$lower))
names(PC3_1.4sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.4sd_pre_lower$age, PC3_1.4sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_1.4sd_pre_lower, indices) {
  d <- PC3_1.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_1.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_1.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_1.4sd_pre_lower_50)[1] <- "PC3"

PC3_1.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_1.4sd_pre_lower_all)[1:ncol(PC3_1.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_1.4sd_pre_lower_all$group <- seq(1,(nrow(PC3_1.4sd_pre_lower_all)), 1)

PC3_1.4sd_pre_lower_all <- reshape2::melt(PC3_1.4sd_pre_lower_all, id.vars = "group")
names(PC3_1.4sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_1.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.2sd_pre_upper <- as.data.frame(cbind(PC3_1.2sd_pre$x, PC3_1.2sd_pre$upper))
names(PC3_1.2sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.2sd_pre_upper$age, PC3_1.2sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_1.2sd_pre_upper, indices) {
  d <- PC3_1.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_1.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_1.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_1.2sd_pre_upper_50)[1] <- "PC3"

PC3_1.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_1.2sd_pre_upper_all)[1:ncol(PC3_1.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_1.2sd_pre_upper_all$group <- seq(1,(nrow(PC3_1.2sd_pre_upper_all)), 1)

PC3_1.2sd_pre_upper_all <- reshape2::melt(PC3_1.2sd_pre_upper_all, id.vars = "group")
names(PC3_1.2sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_1.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.2sd_pre_lower <- as.data.frame(cbind(PC3_1.2sd_pre$x, PC3_1.2sd_pre$lower))
names(PC3_1.2sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.2sd_pre_lower$age, PC3_1.2sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_1.2sd_pre_lower, indices) {
  d <- PC3_1.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_1.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_1.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_1.2sd_pre_lower_50)[1] <- "PC3"

PC3_1.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_1.2sd_pre_lower_all)[1:ncol(PC3_1.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_1.2sd_pre_lower_all$group <- seq(1,(nrow(PC3_1.2sd_pre_lower_all)), 1)

PC3_1.2sd_pre_lower_all <- reshape2::melt(PC3_1.2sd_pre_lower_all, id.vars = "group")
names(PC3_1.2sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_1.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.0sd_pre_upper <- as.data.frame(cbind(PC3_1.0sd_pre$x, PC3_1.0sd_pre$upper))
names(PC3_1.0sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.0sd_pre_upper$age, PC3_1.0sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_1.0sd_pre_upper, indices) {
  d <- PC3_1.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_1.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_1.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_1.0sd_pre_upper_50)[1] <- "PC3"

PC3_1.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_1.0sd_pre_upper_all)[1:ncol(PC3_1.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_1.0sd_pre_upper_all$group <- seq(1,(nrow(PC3_1.0sd_pre_upper_all)), 1)

PC3_1.0sd_pre_upper_all <- reshape2::melt(PC3_1.0sd_pre_upper_all, id.vars = "group")
names(PC3_1.0sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_1.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_1.0sd_pre_lower <- as.data.frame(cbind(PC3_1.0sd_pre$x, PC3_1.0sd_pre$lower))
names(PC3_1.0sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_1.0sd_pre_lower$age, PC3_1.0sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_1.0sd_pre_lower, indices) {
  d <- PC3_1.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_1.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_1.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_1.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_1.0sd_pre_lower_50)[1] <- "PC3"

PC3_1.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_1.0sd_pre_lower_all)[1:ncol(PC3_1.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_1.0sd_pre_lower_all$group <- seq(1,(nrow(PC3_1.0sd_pre_lower_all)), 1)

PC3_1.0sd_pre_lower_all <- reshape2::melt(PC3_1.0sd_pre_lower_all, id.vars = "group")
names(PC3_1.0sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_1.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.8sd_pre_upper <- as.data.frame(cbind(PC3_0.8sd_pre$x, PC3_0.8sd_pre$upper))
names(PC3_0.8sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.8sd_pre_upper$age, PC3_0.8sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_0.8sd_pre_upper, indices) {
  d <- PC3_0.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_0.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_0.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_0.8sd_pre_upper_50)[1] <- "PC3"

PC3_0.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_0.8sd_pre_upper_all)[1:ncol(PC3_0.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_0.8sd_pre_upper_all$group <- seq(1,(nrow(PC3_0.8sd_pre_upper_all)), 1)

PC3_0.8sd_pre_upper_all <- reshape2::melt(PC3_0.8sd_pre_upper_all, id.vars = "group")
names(PC3_0.8sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_0.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.8sd_pre_lower <- as.data.frame(cbind(PC3_0.8sd_pre$x, PC3_0.8sd_pre$lower))
names(PC3_0.8sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.8sd_pre_lower$age, PC3_0.8sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_0.8sd_pre_lower, indices) {
  d <- PC3_0.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_0.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_0.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_0.8sd_pre_lower_50)[1] <- "PC3"

PC3_0.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_0.8sd_pre_lower_all)[1:ncol(PC3_0.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_0.8sd_pre_lower_all$group <- seq(1,(nrow(PC3_0.8sd_pre_lower_all)), 1)

PC3_0.8sd_pre_lower_all <- reshape2::melt(PC3_0.8sd_pre_lower_all, id.vars = "group")
names(PC3_0.8sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_0.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.6sd_pre_upper <- as.data.frame(cbind(PC3_0.6sd_pre$x, PC3_0.6sd_pre$upper))
names(PC3_0.6sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.6sd_pre_upper$age, PC3_0.6sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_0.6sd_pre_upper, indices) {
  d <- PC3_0.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_0.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_0.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_0.6sd_pre_upper_50)[1] <- "PC3"

PC3_0.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_0.6sd_pre_upper_all)[1:ncol(PC3_0.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_0.6sd_pre_upper_all$group <- seq(1,(nrow(PC3_0.6sd_pre_upper_all)), 1)

PC3_0.6sd_pre_upper_all <- reshape2::melt(PC3_0.6sd_pre_upper_all, id.vars = "group")
names(PC3_0.6sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_0.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.6sd_pre_lower <- as.data.frame(cbind(PC3_0.6sd_pre$x, PC3_0.6sd_pre$lower))
names(PC3_0.6sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.6sd_pre_lower$age, PC3_0.6sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_0.6sd_pre_lower, indices) {
  d <- PC3_0.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_0.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_0.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_0.6sd_pre_lower_50)[1] <- "PC3"

PC3_0.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_0.6sd_pre_lower_all)[1:ncol(PC3_0.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_0.6sd_pre_lower_all$group <- seq(1,(nrow(PC3_0.6sd_pre_lower_all)), 1)

PC3_0.6sd_pre_lower_all <- reshape2::melt(PC3_0.6sd_pre_lower_all, id.vars = "group")
names(PC3_0.6sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_0.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.4sd_pre_upper <- as.data.frame(cbind(PC3_0.4sd_pre$x, PC3_0.4sd_pre$upper))
names(PC3_0.4sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.4sd_pre_upper$age, PC3_0.4sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_0.4sd_pre_upper, indices) {
  d <- PC3_0.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_0.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_0.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_0.4sd_pre_upper_50)[1] <- "PC3"

PC3_0.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_0.4sd_pre_upper_all)[1:ncol(PC3_0.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_0.4sd_pre_upper_all$group <- seq(1,(nrow(PC3_0.4sd_pre_upper_all)), 1)

PC3_0.4sd_pre_upper_all <- reshape2::melt(PC3_0.4sd_pre_upper_all, id.vars = "group")
names(PC3_0.4sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_0.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.4sd_pre_lower <- as.data.frame(cbind(PC3_0.4sd_pre$x, PC3_0.4sd_pre$lower))
names(PC3_0.4sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.4sd_pre_lower$age, PC3_0.4sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_0.4sd_pre_lower, indices) {
  d <- PC3_0.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_0.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_0.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_0.4sd_pre_lower_50)[1] <- "PC3"

PC3_0.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_0.4sd_pre_lower_all)[1:ncol(PC3_0.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_0.4sd_pre_lower_all$group <- seq(1,(nrow(PC3_0.4sd_pre_lower_all)), 1)

PC3_0.4sd_pre_lower_all <- reshape2::melt(PC3_0.4sd_pre_lower_all, id.vars = "group")
names(PC3_0.4sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_0.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.2sd_pre_upper <- as.data.frame(cbind(PC3_0.2sd_pre$x, PC3_0.2sd_pre$upper))
names(PC3_0.2sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.2sd_pre_upper$age, PC3_0.2sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_0.2sd_pre_upper, indices) {
  d <- PC3_0.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_0.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_0.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_0.2sd_pre_upper_50)[1] <- "PC3"

PC3_0.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_0.2sd_pre_upper_all)[1:ncol(PC3_0.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_0.2sd_pre_upper_all$group <- seq(1,(nrow(PC3_0.2sd_pre_upper_all)), 1)

PC3_0.2sd_pre_upper_all <- reshape2::melt(PC3_0.2sd_pre_upper_all, id.vars = "group")
names(PC3_0.2sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_0.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.2sd_pre_lower <- as.data.frame(cbind(PC3_0.2sd_pre$x, PC3_0.2sd_pre$lower))
names(PC3_0.2sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.2sd_pre_lower$age, PC3_0.2sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_0.2sd_pre_lower, indices) {
  d <- PC3_0.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_0.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_0.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_0.2sd_pre_lower_50)[1] <- "PC3"

PC3_0.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_0.2sd_pre_lower_all)[1:ncol(PC3_0.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_0.2sd_pre_lower_all$group <- seq(1,(nrow(PC3_0.2sd_pre_lower_all)), 1)

PC3_0.2sd_pre_lower_all <- reshape2::melt(PC3_0.2sd_pre_lower_all, id.vars = "group")
names(PC3_0.2sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_0.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC3 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.0sd_pre_upper <- as.data.frame(cbind(PC3_0.0sd_pre$x, PC3_0.0sd_pre$upper))
names(PC3_0.0sd_pre_upper)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.0sd_pre_upper$age, PC3_0.0sd_pre_upper$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC3_0.0sd_pre_upper, indices) {
  d <- PC3_0.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC3_0.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC3_0.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC3_0.0sd_pre_upper_50)[1] <- "PC3"

PC3_0.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC3_0.0sd_pre_upper_all)[1:ncol(PC3_0.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC3_0.0sd_pre_upper_all$group <- seq(1,(nrow(PC3_0.0sd_pre_upper_all)), 1)

PC3_0.0sd_pre_upper_all <- reshape2::melt(PC3_0.0sd_pre_upper_all, id.vars = "group")
names(PC3_0.0sd_pre_upper_all)[2:3] <- c("age", "PC3")

PC3_0.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC3 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC3_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC3, degree=1, nsigma = nsig, span = span)
PC3_0.0sd_pre_lower <- as.data.frame(cbind(PC3_0.0sd_pre$x, PC3_0.0sd_pre$lower))
names(PC3_0.0sd_pre_lower)[1:2] <- c("age","PC3")

# CV spans
loess.predict <- loess.as(PC3_0.0sd_pre_lower$age, PC3_0.0sd_pre_lower$PC3, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC3_0.0sd_pre_lower, indices) {
  d <- PC3_0.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC3 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC3_0.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC3_0.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC3_0.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC3_0.0sd_pre_lower_50)[1] <- "PC3"

PC3_0.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC3_0.0sd_pre_lower_all)[1:ncol(PC3_0.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC3_0.0sd_pre_lower_all$group <- seq(1,(nrow(PC3_0.0sd_pre_lower_all)), 1)

PC3_0.0sd_pre_lower_all <- reshape2::melt(PC3_0.0sd_pre_lower_all, id.vars = "group")
names(PC3_0.0sd_pre_lower_all)[2:3] <- c("age", "PC3")

PC3_0.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

b <- ggplot() + 
  #geom_line(data = PC3_2sd_pre_upper_all, aes(age2,PC3, group = group)) +
  #geom_line(data = PC3_2sd_pre_lower_all, aes(age2,PC3, group = group)) +
  geom_line(data = PC3_2sd_pre_upper_50, aes(age,PC3), colour = 'blue') +
  geom_line(data = PC3_2sd_pre_lower_50, aes(age,PC3), colour = 'blue') +
  geom_line(data = PC3_1.6sd_pre_upper_50, aes(age,PC3), colour = 'green') +
  geom_line(data = PC3_1.6sd_pre_lower_50, aes(age,PC3), colour = 'green') +
  geom_line(data = PC3_1.4sd_pre_upper_50, aes(age,PC3), colour = 'orange') +
  geom_line(data = PC3_1.4sd_pre_lower_50, aes(age,PC3), colour = 'orange') +
  geom_line(data = PC3_1.2sd_pre_upper_50, aes(age,PC3), colour = 'pink') +
  geom_line(data = PC3_1.2sd_pre_lower_50, aes(age,PC3), colour = 'pink') +
  geom_line(data = PC3_1.0sd_pre_upper_50, aes(age,PC3), colour = 'brown') +
  geom_line(data = PC3_1.0sd_pre_lower_50, aes(age,PC3), colour = 'brown') +
  geom_line(data = PC3_0.8sd_pre_upper_50, aes(age,PC3), colour = 'yellow') +
  geom_line(data = PC3_0.8sd_pre_lower_50, aes(age,PC3), colour = 'yellow') +
  geom_line(data = PC3_0.6sd_pre_upper_50, aes(age,PC3), colour = 'red') +
  geom_line(data = PC3_0.6sd_pre_lower_50, aes(age,PC3), colour = 'red') +
  geom_line(data = PC3_0.4sd_pre_upper_50, aes(age,PC3), colour = 'blue') +
  geom_line(data = PC3_0.4sd_pre_lower_50, aes(age,PC3), colour = 'blue') +
  geom_line(data = PC3_0.2sd_pre_upper_50, aes(age,PC3), colour = 'green') +
  geom_line(data = PC3_0.2sd_pre_lower_50, aes(age,PC3), colour = 'green') +
  geom_line(data = PC3_0.0sd_pre_upper_50, aes(age,PC3), colour = 'orange') +
  geom_line(data = PC3_0.0sd_pre_lower_50, aes(age,PC3), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC3)) + 
  #geom_point(data = PC3_2sd_pre_upper, aes(age, PC3), colour = 'red') +
  scale_x_reverse(limits = c(4000,540)) + theme_bw()

grid.arrange(b,a,ncol=2)





##### Phanerozoic PC4 #####

coords.all_Phan <- subset(coords.all2, age < 541)

# 2 sigma

nsig <- 2

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_2sd_phan_upper <- as.data.frame(cbind(PC4_2sd_phan$x, PC4_2sd_phan$upper))
names(PC4_2sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_2sd_phan_upper$age, PC4_2sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_2sd_phan_upper, indices) {
  d <- PC4_2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_2sd_phan_upper_50)[1] <- "PC4"

PC4_2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_2sd_phan_upper_all)[1:ncol(PC4_2sd_phan_upper_all)] <- seq(0,540, 1)
PC4_2sd_phan_upper_all$group <- seq(1,(nrow(PC4_2sd_phan_upper_all)), 1)

PC4_2sd_phan_upper_all <- reshape2::melt(PC4_2sd_phan_upper_all, id.vars = "group")
names(PC4_2sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_2sd_phan_lower <- as.data.frame(cbind(PC4_2sd_phan$x, PC4_2sd_phan$lower))
names(PC4_2sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_2sd_phan_lower$age, PC4_2sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_2sd_phan_lower, indices) {
  d <- PC4_2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_2sd_phan_lower_50)[1] <- "PC4"

PC4_2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_2sd_phan_lower_all)[1:ncol(PC4_2sd_phan_lower_all)] <- seq(0,540, 1)
PC4_2sd_phan_lower_all$group <- seq(1,(nrow(PC4_2sd_phan_lower_all)), 1)

PC4_2sd_phan_lower_all <- reshape2::melt(PC4_2sd_phan_lower_all, id.vars = "group")
names(PC4_2sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.8sd_phan_upper <- as.data.frame(cbind(PC4_1.8sd_phan$x, PC4_1.8sd_phan$upper))
names(PC4_1.8sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.8sd_phan_upper$age, PC4_1.8sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_1.8sd_phan_upper, indices) {
  d <- PC4_1.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_1.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_1.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_1.8sd_phan_upper_50)[1] <- "PC4"

PC4_1.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_1.8sd_phan_upper_all)[1:ncol(PC4_1.8sd_phan_upper_all)] <- seq(0,540, 1)
PC4_1.8sd_phan_upper_all$group <- seq(1,(nrow(PC4_1.8sd_phan_upper_all)), 1)

PC4_1.8sd_phan_upper_all <- reshape2::melt(PC4_1.8sd_phan_upper_all, id.vars = "group")
names(PC4_1.8sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_1.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.8sd_phan_lower <- as.data.frame(cbind(PC4_1.8sd_phan$x, PC4_1.8sd_phan$lower))
names(PC4_1.8sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.8sd_phan_lower$age, PC4_1.8sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_1.8sd_phan_lower, indices) {
  d <- PC4_1.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_1.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_1.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_1.8sd_phan_lower_50)[1] <- "PC4"

PC4_1.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_1.8sd_phan_lower_all)[1:ncol(PC4_1.8sd_phan_lower_all)] <- seq(0,540, 1)
PC4_1.8sd_phan_lower_all$group <- seq(1,(nrow(PC4_1.8sd_phan_lower_all)), 1)

PC4_1.8sd_phan_lower_all <- reshape2::melt(PC4_1.8sd_phan_lower_all, id.vars = "group")
names(PC4_1.8sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_1.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.6sd_phan_upper <- as.data.frame(cbind(PC4_1.6sd_phan$x, PC4_1.6sd_phan$upper))
names(PC4_1.6sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.6sd_phan_upper$age, PC4_1.6sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_1.6sd_phan_upper, indices) {
  d <- PC4_1.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_1.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_1.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_1.6sd_phan_upper_50)[1] <- "PC4"

PC4_1.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_1.6sd_phan_upper_all)[1:ncol(PC4_1.6sd_phan_upper_all)] <- seq(0,540, 1)
PC4_1.6sd_phan_upper_all$group <- seq(1,(nrow(PC4_1.6sd_phan_upper_all)), 1)

PC4_1.6sd_phan_upper_all <- reshape2::melt(PC4_1.6sd_phan_upper_all, id.vars = "group")
names(PC4_1.6sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_1.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.6sd_phan_lower <- as.data.frame(cbind(PC4_1.6sd_phan$x, PC4_1.6sd_phan$lower))
names(PC4_1.6sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.6sd_phan_lower$age, PC4_1.6sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_1.6sd_phan_lower, indices) {
  d <- PC4_1.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_1.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_1.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_1.6sd_phan_lower_50)[1] <- "PC4"

PC4_1.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_1.6sd_phan_lower_all)[1:ncol(PC4_1.6sd_phan_lower_all)] <- seq(0,540, 1)
PC4_1.6sd_phan_lower_all$group <- seq(1,(nrow(PC4_1.6sd_phan_lower_all)), 1)

PC4_1.6sd_phan_lower_all <- reshape2::melt(PC4_1.6sd_phan_lower_all, id.vars = "group")
names(PC4_1.6sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_1.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.4sd_phan_upper <- as.data.frame(cbind(PC4_1.4sd_phan$x, PC4_1.4sd_phan$upper))
names(PC4_1.4sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.4sd_phan_upper$age, PC4_1.4sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_1.4sd_phan_upper, indices) {
  d <- PC4_1.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_1.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_1.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_1.4sd_phan_upper_50)[1] <- "PC4"

PC4_1.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_1.4sd_phan_upper_all)[1:ncol(PC4_1.4sd_phan_upper_all)] <- seq(0,540, 1)
PC4_1.4sd_phan_upper_all$group <- seq(1,(nrow(PC4_1.4sd_phan_upper_all)), 1)

PC4_1.4sd_phan_upper_all <- reshape2::melt(PC4_1.4sd_phan_upper_all, id.vars = "group")
names(PC4_1.4sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_1.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.4sd_phan_lower <- as.data.frame(cbind(PC4_1.4sd_phan$x, PC4_1.4sd_phan$lower))
names(PC4_1.4sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.4sd_phan_lower$age, PC4_1.4sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_1.4sd_phan_lower, indices) {
  d <- PC4_1.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_1.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_1.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_1.4sd_phan_lower_50)[1] <- "PC4"

PC4_1.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_1.4sd_phan_lower_all)[1:ncol(PC4_1.4sd_phan_lower_all)] <- seq(0,540, 1)
PC4_1.4sd_phan_lower_all$group <- seq(1,(nrow(PC4_1.4sd_phan_lower_all)), 1)

PC4_1.4sd_phan_lower_all <- reshape2::melt(PC4_1.4sd_phan_lower_all, id.vars = "group")
names(PC4_1.4sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_1.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.2sd_phan_upper <- as.data.frame(cbind(PC4_1.2sd_phan$x, PC4_1.2sd_phan$upper))
names(PC4_1.2sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.2sd_phan_upper$age, PC4_1.2sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_1.2sd_phan_upper, indices) {
  d <- PC4_1.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_1.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_1.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_1.2sd_phan_upper_50)[1] <- "PC4"

PC4_1.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_1.2sd_phan_upper_all)[1:ncol(PC4_1.2sd_phan_upper_all)] <- seq(0,540, 1)
PC4_1.2sd_phan_upper_all$group <- seq(1,(nrow(PC4_1.2sd_phan_upper_all)), 1)

PC4_1.2sd_phan_upper_all <- reshape2::melt(PC4_1.2sd_phan_upper_all, id.vars = "group")
names(PC4_1.2sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_1.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.2sd_phan_lower <- as.data.frame(cbind(PC4_1.2sd_phan$x, PC4_1.2sd_phan$lower))
names(PC4_1.2sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.2sd_phan_lower$age, PC4_1.2sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_1.2sd_phan_lower, indices) {
  d <- PC4_1.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_1.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_1.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_1.2sd_phan_lower_50)[1] <- "PC4"

PC4_1.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_1.2sd_phan_lower_all)[1:ncol(PC4_1.2sd_phan_lower_all)] <- seq(0,540, 1)
PC4_1.2sd_phan_lower_all$group <- seq(1,(nrow(PC4_1.2sd_phan_lower_all)), 1)

PC4_1.2sd_phan_lower_all <- reshape2::melt(PC4_1.2sd_phan_lower_all, id.vars = "group")
names(PC4_1.2sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_1.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.0sd_phan_upper <- as.data.frame(cbind(PC4_1.0sd_phan$x, PC4_1.0sd_phan$upper))
names(PC4_1.0sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.0sd_phan_upper$age, PC4_1.0sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_1.0sd_phan_upper, indices) {
  d <- PC4_1.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_1.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_1.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_1.0sd_phan_upper_50)[1] <- "PC4"

PC4_1.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_1.0sd_phan_upper_all)[1:ncol(PC4_1.0sd_phan_upper_all)] <- seq(0,540, 1)
PC4_1.0sd_phan_upper_all$group <- seq(1,(nrow(PC4_1.0sd_phan_upper_all)), 1)

PC4_1.0sd_phan_upper_all <- reshape2::melt(PC4_1.0sd_phan_upper_all, id.vars = "group")
names(PC4_1.0sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_1.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.0sd_phan_lower <- as.data.frame(cbind(PC4_1.0sd_phan$x, PC4_1.0sd_phan$lower))
names(PC4_1.0sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.0sd_phan_lower$age, PC4_1.0sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_1.0sd_phan_lower, indices) {
  d <- PC4_1.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_1.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_1.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_1.0sd_phan_lower_50)[1] <- "PC4"

PC4_1.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_1.0sd_phan_lower_all)[1:ncol(PC4_1.0sd_phan_lower_all)] <- seq(0,540, 1)
PC4_1.0sd_phan_lower_all$group <- seq(1,(nrow(PC4_1.0sd_phan_lower_all)), 1)

PC4_1.0sd_phan_lower_all <- reshape2::melt(PC4_1.0sd_phan_lower_all, id.vars = "group")
names(PC4_1.0sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_1.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.8sd_phan_upper <- as.data.frame(cbind(PC4_0.8sd_phan$x, PC4_0.8sd_phan$upper))
names(PC4_0.8sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.8sd_phan_upper$age, PC4_0.8sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_0.8sd_phan_upper, indices) {
  d <- PC4_0.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_0.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_0.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_0.8sd_phan_upper_50)[1] <- "PC4"

PC4_0.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_0.8sd_phan_upper_all)[1:ncol(PC4_0.8sd_phan_upper_all)] <- seq(0,540, 1)
PC4_0.8sd_phan_upper_all$group <- seq(1,(nrow(PC4_0.8sd_phan_upper_all)), 1)

PC4_0.8sd_phan_upper_all <- reshape2::melt(PC4_0.8sd_phan_upper_all, id.vars = "group")
names(PC4_0.8sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_0.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.8sd_phan_lower <- as.data.frame(cbind(PC4_0.8sd_phan$x, PC4_0.8sd_phan$lower))
names(PC4_0.8sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.8sd_phan_lower$age, PC4_0.8sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_0.8sd_phan_lower, indices) {
  d <- PC4_0.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_0.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_0.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_0.8sd_phan_lower_50)[1] <- "PC4"

PC4_0.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_0.8sd_phan_lower_all)[1:ncol(PC4_0.8sd_phan_lower_all)] <- seq(0,540, 1)
PC4_0.8sd_phan_lower_all$group <- seq(1,(nrow(PC4_0.8sd_phan_lower_all)), 1)

PC4_0.8sd_phan_lower_all <- reshape2::melt(PC4_0.8sd_phan_lower_all, id.vars = "group")
names(PC4_0.8sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_0.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.6sd_phan_upper <- as.data.frame(cbind(PC4_0.6sd_phan$x, PC4_0.6sd_phan$upper))
names(PC4_0.6sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.6sd_phan_upper$age, PC4_0.6sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_0.6sd_phan_upper, indices) {
  d <- PC4_0.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_0.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_0.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_0.6sd_phan_upper_50)[1] <- "PC4"

PC4_0.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_0.6sd_phan_upper_all)[1:ncol(PC4_0.6sd_phan_upper_all)] <- seq(0,540, 1)
PC4_0.6sd_phan_upper_all$group <- seq(1,(nrow(PC4_0.6sd_phan_upper_all)), 1)

PC4_0.6sd_phan_upper_all <- reshape2::melt(PC4_0.6sd_phan_upper_all, id.vars = "group")
names(PC4_0.6sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_0.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.6sd_phan_lower <- as.data.frame(cbind(PC4_0.6sd_phan$x, PC4_0.6sd_phan$lower))
names(PC4_0.6sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.6sd_phan_lower$age, PC4_0.6sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_0.6sd_phan_lower, indices) {
  d <- PC4_0.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_0.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_0.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_0.6sd_phan_lower_50)[1] <- "PC4"

PC4_0.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_0.6sd_phan_lower_all)[1:ncol(PC4_0.6sd_phan_lower_all)] <- seq(0,540, 1)
PC4_0.6sd_phan_lower_all$group <- seq(1,(nrow(PC4_0.6sd_phan_lower_all)), 1)

PC4_0.6sd_phan_lower_all <- reshape2::melt(PC4_0.6sd_phan_lower_all, id.vars = "group")
names(PC4_0.6sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_0.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.4sd_phan_upper <- as.data.frame(cbind(PC4_0.4sd_phan$x, PC4_0.4sd_phan$upper))
names(PC4_0.4sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.4sd_phan_upper$age, PC4_0.4sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_0.4sd_phan_upper, indices) {
  d <- PC4_0.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_0.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_0.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_0.4sd_phan_upper_50)[1] <- "PC4"

PC4_0.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_0.4sd_phan_upper_all)[1:ncol(PC4_0.4sd_phan_upper_all)] <- seq(0,540, 1)
PC4_0.4sd_phan_upper_all$group <- seq(1,(nrow(PC4_0.4sd_phan_upper_all)), 1)

PC4_0.4sd_phan_upper_all <- reshape2::melt(PC4_0.4sd_phan_upper_all, id.vars = "group")
names(PC4_0.4sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_0.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.4sd_phan_lower <- as.data.frame(cbind(PC4_0.4sd_phan$x, PC4_0.4sd_phan$lower))
names(PC4_0.4sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.4sd_phan_lower$age, PC4_0.4sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_0.4sd_phan_lower, indices) {
  d <- PC4_0.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_0.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_0.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_0.4sd_phan_lower_50)[1] <- "PC4"

PC4_0.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_0.4sd_phan_lower_all)[1:ncol(PC4_0.4sd_phan_lower_all)] <- seq(0,540, 1)
PC4_0.4sd_phan_lower_all$group <- seq(1,(nrow(PC4_0.4sd_phan_lower_all)), 1)

PC4_0.4sd_phan_lower_all <- reshape2::melt(PC4_0.4sd_phan_lower_all, id.vars = "group")
names(PC4_0.4sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_0.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.2sd_phan_upper <- as.data.frame(cbind(PC4_0.2sd_phan$x, PC4_0.2sd_phan$upper))
names(PC4_0.2sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.2sd_phan_upper$age, PC4_0.2sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_0.2sd_phan_upper, indices) {
  d <- PC4_0.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_0.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_0.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_0.2sd_phan_upper_50)[1] <- "PC4"

PC4_0.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_0.2sd_phan_upper_all)[1:ncol(PC4_0.2sd_phan_upper_all)] <- seq(0,540, 1)
PC4_0.2sd_phan_upper_all$group <- seq(1,(nrow(PC4_0.2sd_phan_upper_all)), 1)

PC4_0.2sd_phan_upper_all <- reshape2::melt(PC4_0.2sd_phan_upper_all, id.vars = "group")
names(PC4_0.2sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_0.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.2sd_phan_lower <- as.data.frame(cbind(PC4_0.2sd_phan$x, PC4_0.2sd_phan$lower))
names(PC4_0.2sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.2sd_phan_lower$age, PC4_0.2sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_0.2sd_phan_lower, indices) {
  d <- PC4_0.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_0.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_0.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_0.2sd_phan_lower_50)[1] <- "PC4"

PC4_0.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_0.2sd_phan_lower_all)[1:ncol(PC4_0.2sd_phan_lower_all)] <- seq(0,540, 1)
PC4_0.2sd_phan_lower_all$group <- seq(1,(nrow(PC4_0.2sd_phan_lower_all)), 1)

PC4_0.2sd_phan_lower_all <- reshape2::melt(PC4_0.2sd_phan_lower_all, id.vars = "group")
names(PC4_0.2sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_0.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.0sd_phan_upper <- as.data.frame(cbind(PC4_0.0sd_phan$x, PC4_0.0sd_phan$upper))
names(PC4_0.0sd_phan_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.0sd_phan_upper$age, PC4_0.0sd_phan_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_0.0sd_phan_upper, indices) {
  d <- PC4_0.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_0.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_0.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC4_0.0sd_phan_upper_50)[1] <- "PC4"

PC4_0.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC4_0.0sd_phan_upper_all)[1:ncol(PC4_0.0sd_phan_upper_all)] <- seq(0,540, 1)
PC4_0.0sd_phan_upper_all$group <- seq(1,(nrow(PC4_0.0sd_phan_upper_all)), 1)

PC4_0.0sd_phan_upper_all <- reshape2::melt(PC4_0.0sd_phan_upper_all, id.vars = "group")
names(PC4_0.0sd_phan_upper_all)[2:3] <- c("age", "PC4")

PC4_0.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.0sd_phan_lower <- as.data.frame(cbind(PC4_0.0sd_phan$x, PC4_0.0sd_phan$lower))
names(PC4_0.0sd_phan_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.0sd_phan_lower$age, PC4_0.0sd_phan_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_0.0sd_phan_lower, indices) {
  d <- PC4_0.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_0.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_0.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC4_0.0sd_phan_lower_50)[1] <- "PC4"

PC4_0.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC4_0.0sd_phan_lower_all)[1:ncol(PC4_0.0sd_phan_lower_all)] <- seq(0,540, 1)
PC4_0.0sd_phan_lower_all$group <- seq(1,(nrow(PC4_0.0sd_phan_lower_all)), 1)

PC4_0.0sd_phan_lower_all <- reshape2::melt(PC4_0.0sd_phan_lower_all, id.vars = "group")
names(PC4_0.0sd_phan_lower_all)[2:3] <- c("age", "PC4")

PC4_0.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

a <- ggplot() + 
  #geom_line(data = PC4_2sd_phan_upper_all, aes(age2,PC4, group = group)) +
  #geom_line(data = PC4_2sd_phan_lower_all, aes(age2,PC4, group = group)) +
  geom_line(data = PC4_2sd_phan_upper_50, aes(age,PC4), colour = 'blue') +
  geom_line(data = PC4_2sd_phan_lower_50, aes(age,PC4), colour = 'blue') +
  geom_line(data = PC4_1.6sd_phan_upper_50, aes(age,PC4), colour = 'green') +
  geom_line(data = PC4_1.6sd_phan_lower_50, aes(age,PC4), colour = 'green') +
  geom_line(data = PC4_1.4sd_phan_upper_50, aes(age,PC4), colour = 'orange') +
  geom_line(data = PC4_1.4sd_phan_lower_50, aes(age,PC4), colour = 'orange') +
  geom_line(data = PC4_1.2sd_phan_upper_50, aes(age,PC4), colour = 'pink') +
  geom_line(data = PC4_1.2sd_phan_lower_50, aes(age,PC4), colour = 'pink') +
  geom_line(data = PC4_1.0sd_phan_upper_50, aes(age,PC4), colour = 'brown') +
  geom_line(data = PC4_1.0sd_phan_lower_50, aes(age,PC4), colour = 'brown') +
  geom_line(data = PC4_0.8sd_phan_upper_50, aes(age,PC4), colour = 'yellow') +
  geom_line(data = PC4_0.8sd_phan_lower_50, aes(age,PC4), colour = 'yellow') +
  geom_line(data = PC4_0.6sd_phan_upper_50, aes(age,PC4), colour = 'red') +
  geom_line(data = PC4_0.6sd_phan_lower_50, aes(age,PC4), colour = 'red') +
  geom_line(data = PC4_0.4sd_phan_upper_50, aes(age,PC4), colour = 'blue') +
  geom_line(data = PC4_0.4sd_phan_lower_50, aes(age,PC4), colour = 'blue') +
  geom_line(data = PC4_0.2sd_phan_upper_50, aes(age,PC4), colour = 'green') +
  geom_line(data = PC4_0.2sd_phan_lower_50, aes(age,PC4), colour = 'green') +
  geom_line(data = PC4_0.0sd_phan_upper_50, aes(age,PC4), colour = 'orange') +
  geom_line(data = PC4_0.0sd_phan_lower_50, aes(age,PC4), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC4)) + 
  #geom_point(data = PC4_2sd_phan_upper, aes(age, PC4), colour = 'red') +
  scale_x_reverse(limits = c(540,0)) + theme_bw()

##### Precambrian PC4 #####

coords.all_pre <- subset(coords.all2, age > 541)

# 2 sigma

nsig <- 2

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_2sd_pre_upper <- as.data.frame(cbind(PC4_2sd_pre$x, PC4_2sd_pre$upper))
names(PC4_2sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_2sd_pre_upper$age, PC4_2sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[1], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_2sd_pre_upper, indices) {
  d <- PC4_2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_2sd_pre_upper_50)[1] <- "PC4"

PC4_2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_2sd_pre_upper_all)[1:ncol(PC4_2sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_2sd_pre_upper_all$group <- seq(1,(nrow(PC4_2sd_pre_upper_all)), 1)

PC4_2sd_pre_upper_all <- reshape2::melt(PC4_2sd_pre_upper_all, id.vars = "group")
names(PC4_2sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_2sd_pre_lower <- as.data.frame(cbind(PC4_2sd_pre$x, PC4_2sd_pre$lower))
names(PC4_2sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_2sd_pre_lower$age, PC4_2sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_2sd_pre_lower, indices) {
  d <- PC4_2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_2sd_pre_lower_50)[1] <- "PC4"

PC4_2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_2sd_pre_lower_all)[1:ncol(PC4_2sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_2sd_pre_lower_all$group <- seq(1,(nrow(PC4_2sd_pre_lower_all)), 1)

PC4_2sd_pre_lower_all <- reshape2::melt(PC4_2sd_pre_lower_all, id.vars = "group")
names(PC4_2sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.8sd_pre_upper <- as.data.frame(cbind(PC4_1.8sd_pre$x, PC4_1.8sd_pre$upper))
names(PC4_1.8sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.8sd_pre_upper$age, PC4_1.8sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_1.8sd_pre_upper, indices) {
  d <- PC4_1.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_1.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_1.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_1.8sd_pre_upper_50)[1] <- "PC4"

PC4_1.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_1.8sd_pre_upper_all)[1:ncol(PC4_1.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_1.8sd_pre_upper_all$group <- seq(1,(nrow(PC4_1.8sd_pre_upper_all)), 1)

PC4_1.8sd_pre_upper_all <- reshape2::melt(PC4_1.8sd_pre_upper_all, id.vars = "group")
names(PC4_1.8sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_1.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.8sd_pre_lower <- as.data.frame(cbind(PC4_1.8sd_pre$x, PC4_1.8sd_pre$lower))
names(PC4_1.8sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.8sd_pre_lower$age, PC4_1.8sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_1.8sd_pre_lower, indices) {
  d <- PC4_1.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_1.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_1.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_1.8sd_pre_lower_50)[1] <- "PC4"

PC4_1.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_1.8sd_pre_lower_all)[1:ncol(PC4_1.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_1.8sd_pre_lower_all$group <- seq(1,(nrow(PC4_1.8sd_pre_lower_all)), 1)

PC4_1.8sd_pre_lower_all <- reshape2::melt(PC4_1.8sd_pre_lower_all, id.vars = "group")
names(PC4_1.8sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_1.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.6sd_pre_upper <- as.data.frame(cbind(PC4_1.6sd_pre$x, PC4_1.6sd_pre$upper))
names(PC4_1.6sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.6sd_pre_upper$age, PC4_1.6sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_1.6sd_pre_upper, indices) {
  d <- PC4_1.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_1.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_1.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_1.6sd_pre_upper_50)[1] <- "PC4"

PC4_1.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_1.6sd_pre_upper_all)[1:ncol(PC4_1.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_1.6sd_pre_upper_all$group <- seq(1,(nrow(PC4_1.6sd_pre_upper_all)), 1)

PC4_1.6sd_pre_upper_all <- reshape2::melt(PC4_1.6sd_pre_upper_all, id.vars = "group")
names(PC4_1.6sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_1.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.6sd_pre_lower <- as.data.frame(cbind(PC4_1.6sd_pre$x, PC4_1.6sd_pre$lower))
names(PC4_1.6sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.6sd_pre_lower$age, PC4_1.6sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_1.6sd_pre_lower, indices) {
  d <- PC4_1.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_1.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_1.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_1.6sd_pre_lower_50)[1] <- "PC4"

PC4_1.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_1.6sd_pre_lower_all)[1:ncol(PC4_1.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_1.6sd_pre_lower_all$group <- seq(1,(nrow(PC4_1.6sd_pre_lower_all)), 1)

PC4_1.6sd_pre_lower_all <- reshape2::melt(PC4_1.6sd_pre_lower_all, id.vars = "group")
names(PC4_1.6sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_1.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.4sd_pre_upper <- as.data.frame(cbind(PC4_1.4sd_pre$x, PC4_1.4sd_pre$upper))
names(PC4_1.4sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.4sd_pre_upper$age, PC4_1.4sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_1.4sd_pre_upper, indices) {
  d <- PC4_1.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_1.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_1.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_1.4sd_pre_upper_50)[1] <- "PC4"

PC4_1.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_1.4sd_pre_upper_all)[1:ncol(PC4_1.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_1.4sd_pre_upper_all$group <- seq(1,(nrow(PC4_1.4sd_pre_upper_all)), 1)

PC4_1.4sd_pre_upper_all <- reshape2::melt(PC4_1.4sd_pre_upper_all, id.vars = "group")
names(PC4_1.4sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_1.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.4sd_pre_lower <- as.data.frame(cbind(PC4_1.4sd_pre$x, PC4_1.4sd_pre$lower))
names(PC4_1.4sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.4sd_pre_lower$age, PC4_1.4sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_1.4sd_pre_lower, indices) {
  d <- PC4_1.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_1.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_1.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_1.4sd_pre_lower_50)[1] <- "PC4"

PC4_1.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_1.4sd_pre_lower_all)[1:ncol(PC4_1.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_1.4sd_pre_lower_all$group <- seq(1,(nrow(PC4_1.4sd_pre_lower_all)), 1)

PC4_1.4sd_pre_lower_all <- reshape2::melt(PC4_1.4sd_pre_lower_all, id.vars = "group")
names(PC4_1.4sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_1.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.2sd_pre_upper <- as.data.frame(cbind(PC4_1.2sd_pre$x, PC4_1.2sd_pre$upper))
names(PC4_1.2sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.2sd_pre_upper$age, PC4_1.2sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_1.2sd_pre_upper, indices) {
  d <- PC4_1.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_1.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_1.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_1.2sd_pre_upper_50)[1] <- "PC4"

PC4_1.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_1.2sd_pre_upper_all)[1:ncol(PC4_1.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_1.2sd_pre_upper_all$group <- seq(1,(nrow(PC4_1.2sd_pre_upper_all)), 1)

PC4_1.2sd_pre_upper_all <- reshape2::melt(PC4_1.2sd_pre_upper_all, id.vars = "group")
names(PC4_1.2sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_1.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.2sd_pre_lower <- as.data.frame(cbind(PC4_1.2sd_pre$x, PC4_1.2sd_pre$lower))
names(PC4_1.2sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.2sd_pre_lower$age, PC4_1.2sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_1.2sd_pre_lower, indices) {
  d <- PC4_1.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_1.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_1.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_1.2sd_pre_lower_50)[1] <- "PC4"

PC4_1.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_1.2sd_pre_lower_all)[1:ncol(PC4_1.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_1.2sd_pre_lower_all$group <- seq(1,(nrow(PC4_1.2sd_pre_lower_all)), 1)

PC4_1.2sd_pre_lower_all <- reshape2::melt(PC4_1.2sd_pre_lower_all, id.vars = "group")
names(PC4_1.2sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_1.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.0sd_pre_upper <- as.data.frame(cbind(PC4_1.0sd_pre$x, PC4_1.0sd_pre$upper))
names(PC4_1.0sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.0sd_pre_upper$age, PC4_1.0sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_1.0sd_pre_upper, indices) {
  d <- PC4_1.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_1.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_1.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_1.0sd_pre_upper_50)[1] <- "PC4"

PC4_1.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_1.0sd_pre_upper_all)[1:ncol(PC4_1.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_1.0sd_pre_upper_all$group <- seq(1,(nrow(PC4_1.0sd_pre_upper_all)), 1)

PC4_1.0sd_pre_upper_all <- reshape2::melt(PC4_1.0sd_pre_upper_all, id.vars = "group")
names(PC4_1.0sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_1.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_1.0sd_pre_lower <- as.data.frame(cbind(PC4_1.0sd_pre$x, PC4_1.0sd_pre$lower))
names(PC4_1.0sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_1.0sd_pre_lower$age, PC4_1.0sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_1.0sd_pre_lower, indices) {
  d <- PC4_1.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_1.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_1.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_1.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_1.0sd_pre_lower_50)[1] <- "PC4"

PC4_1.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_1.0sd_pre_lower_all)[1:ncol(PC4_1.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_1.0sd_pre_lower_all$group <- seq(1,(nrow(PC4_1.0sd_pre_lower_all)), 1)

PC4_1.0sd_pre_lower_all <- reshape2::melt(PC4_1.0sd_pre_lower_all, id.vars = "group")
names(PC4_1.0sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_1.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.8sd_pre_upper <- as.data.frame(cbind(PC4_0.8sd_pre$x, PC4_0.8sd_pre$upper))
names(PC4_0.8sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.8sd_pre_upper$age, PC4_0.8sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_0.8sd_pre_upper, indices) {
  d <- PC4_0.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_0.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_0.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_0.8sd_pre_upper_50)[1] <- "PC4"

PC4_0.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_0.8sd_pre_upper_all)[1:ncol(PC4_0.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_0.8sd_pre_upper_all$group <- seq(1,(nrow(PC4_0.8sd_pre_upper_all)), 1)

PC4_0.8sd_pre_upper_all <- reshape2::melt(PC4_0.8sd_pre_upper_all, id.vars = "group")
names(PC4_0.8sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_0.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.8sd_pre_lower <- as.data.frame(cbind(PC4_0.8sd_pre$x, PC4_0.8sd_pre$lower))
names(PC4_0.8sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.8sd_pre_lower$age, PC4_0.8sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_0.8sd_pre_lower, indices) {
  d <- PC4_0.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_0.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_0.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_0.8sd_pre_lower_50)[1] <- "PC4"

PC4_0.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_0.8sd_pre_lower_all)[1:ncol(PC4_0.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_0.8sd_pre_lower_all$group <- seq(1,(nrow(PC4_0.8sd_pre_lower_all)), 1)

PC4_0.8sd_pre_lower_all <- reshape2::melt(PC4_0.8sd_pre_lower_all, id.vars = "group")
names(PC4_0.8sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_0.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.6sd_pre_upper <- as.data.frame(cbind(PC4_0.6sd_pre$x, PC4_0.6sd_pre$upper))
names(PC4_0.6sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.6sd_pre_upper$age, PC4_0.6sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_0.6sd_pre_upper, indices) {
  d <- PC4_0.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_0.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_0.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_0.6sd_pre_upper_50)[1] <- "PC4"

PC4_0.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_0.6sd_pre_upper_all)[1:ncol(PC4_0.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_0.6sd_pre_upper_all$group <- seq(1,(nrow(PC4_0.6sd_pre_upper_all)), 1)

PC4_0.6sd_pre_upper_all <- reshape2::melt(PC4_0.6sd_pre_upper_all, id.vars = "group")
names(PC4_0.6sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_0.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.6sd_pre_lower <- as.data.frame(cbind(PC4_0.6sd_pre$x, PC4_0.6sd_pre$lower))
names(PC4_0.6sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.6sd_pre_lower$age, PC4_0.6sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_0.6sd_pre_lower, indices) {
  d <- PC4_0.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_0.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_0.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_0.6sd_pre_lower_50)[1] <- "PC4"

PC4_0.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_0.6sd_pre_lower_all)[1:ncol(PC4_0.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_0.6sd_pre_lower_all$group <- seq(1,(nrow(PC4_0.6sd_pre_lower_all)), 1)

PC4_0.6sd_pre_lower_all <- reshape2::melt(PC4_0.6sd_pre_lower_all, id.vars = "group")
names(PC4_0.6sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_0.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.4sd_pre_upper <- as.data.frame(cbind(PC4_0.4sd_pre$x, PC4_0.4sd_pre$upper))
names(PC4_0.4sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.4sd_pre_upper$age, PC4_0.4sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_0.4sd_pre_upper, indices) {
  d <- PC4_0.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_0.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_0.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_0.4sd_pre_upper_50)[1] <- "PC4"

PC4_0.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_0.4sd_pre_upper_all)[1:ncol(PC4_0.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_0.4sd_pre_upper_all$group <- seq(1,(nrow(PC4_0.4sd_pre_upper_all)), 1)

PC4_0.4sd_pre_upper_all <- reshape2::melt(PC4_0.4sd_pre_upper_all, id.vars = "group")
names(PC4_0.4sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_0.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.4sd_pre_lower <- as.data.frame(cbind(PC4_0.4sd_pre$x, PC4_0.4sd_pre$lower))
names(PC4_0.4sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.4sd_pre_lower$age, PC4_0.4sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_0.4sd_pre_lower, indices) {
  d <- PC4_0.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_0.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_0.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_0.4sd_pre_lower_50)[1] <- "PC4"

PC4_0.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_0.4sd_pre_lower_all)[1:ncol(PC4_0.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_0.4sd_pre_lower_all$group <- seq(1,(nrow(PC4_0.4sd_pre_lower_all)), 1)

PC4_0.4sd_pre_lower_all <- reshape2::melt(PC4_0.4sd_pre_lower_all, id.vars = "group")
names(PC4_0.4sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_0.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.2sd_pre_upper <- as.data.frame(cbind(PC4_0.2sd_pre$x, PC4_0.2sd_pre$upper))
names(PC4_0.2sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.2sd_pre_upper$age, PC4_0.2sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_0.2sd_pre_upper, indices) {
  d <- PC4_0.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_0.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_0.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_0.2sd_pre_upper_50)[1] <- "PC4"

PC4_0.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_0.2sd_pre_upper_all)[1:ncol(PC4_0.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_0.2sd_pre_upper_all$group <- seq(1,(nrow(PC4_0.2sd_pre_upper_all)), 1)

PC4_0.2sd_pre_upper_all <- reshape2::melt(PC4_0.2sd_pre_upper_all, id.vars = "group")
names(PC4_0.2sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_0.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.2sd_pre_lower <- as.data.frame(cbind(PC4_0.2sd_pre$x, PC4_0.2sd_pre$lower))
names(PC4_0.2sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.2sd_pre_lower$age, PC4_0.2sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_0.2sd_pre_lower, indices) {
  d <- PC4_0.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_0.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_0.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_0.2sd_pre_lower_50)[1] <- "PC4"

PC4_0.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_0.2sd_pre_lower_all)[1:ncol(PC4_0.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_0.2sd_pre_lower_all$group <- seq(1,(nrow(PC4_0.2sd_pre_lower_all)), 1)

PC4_0.2sd_pre_lower_all <- reshape2::melt(PC4_0.2sd_pre_lower_all, id.vars = "group")
names(PC4_0.2sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_0.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC4 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.0sd_pre_upper <- as.data.frame(cbind(PC4_0.0sd_pre$x, PC4_0.0sd_pre$upper))
names(PC4_0.0sd_pre_upper)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.0sd_pre_upper$age, PC4_0.0sd_pre_upper$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC4_0.0sd_pre_upper, indices) {
  d <- PC4_0.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC4_0.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC4_0.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC4_0.0sd_pre_upper_50)[1] <- "PC4"

PC4_0.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC4_0.0sd_pre_upper_all)[1:ncol(PC4_0.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC4_0.0sd_pre_upper_all$group <- seq(1,(nrow(PC4_0.0sd_pre_upper_all)), 1)

PC4_0.0sd_pre_upper_all <- reshape2::melt(PC4_0.0sd_pre_upper_all, id.vars = "group")
names(PC4_0.0sd_pre_upper_all)[2:3] <- c("age", "PC4")

PC4_0.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC4 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC4_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC4, degree=1, nsigma = nsig, span = span)
PC4_0.0sd_pre_lower <- as.data.frame(cbind(PC4_0.0sd_pre$x, PC4_0.0sd_pre$lower))
names(PC4_0.0sd_pre_lower)[1:2] <- c("age","PC4")

# CV spans
loess.predict <- loess.as(PC4_0.0sd_pre_lower$age, PC4_0.0sd_pre_lower$PC4, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC4_0.0sd_pre_lower, indices) {
  d <- PC4_0.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC4 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC4_0.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC4_0.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC4_0.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC4_0.0sd_pre_lower_50)[1] <- "PC4"

PC4_0.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC4_0.0sd_pre_lower_all)[1:ncol(PC4_0.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC4_0.0sd_pre_lower_all$group <- seq(1,(nrow(PC4_0.0sd_pre_lower_all)), 1)

PC4_0.0sd_pre_lower_all <- reshape2::melt(PC4_0.0sd_pre_lower_all, id.vars = "group")
names(PC4_0.0sd_pre_lower_all)[2:3] <- c("age", "PC4")

PC4_0.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

b <- ggplot() + 
  #geom_line(data = PC4_2sd_pre_upper_all, aes(age2,PC4, group = group)) +
  #geom_line(data = PC4_2sd_pre_lower_all, aes(age2,PC4, group = group)) +
  geom_line(data = PC4_2sd_pre_upper_50, aes(age,PC4), colour = 'blue') +
  geom_line(data = PC4_2sd_pre_lower_50, aes(age,PC4), colour = 'blue') +
  geom_line(data = PC4_1.6sd_pre_upper_50, aes(age,PC4), colour = 'green') +
  geom_line(data = PC4_1.6sd_pre_lower_50, aes(age,PC4), colour = 'green') +
  geom_line(data = PC4_1.4sd_pre_upper_50, aes(age,PC4), colour = 'orange') +
  geom_line(data = PC4_1.4sd_pre_lower_50, aes(age,PC4), colour = 'orange') +
  geom_line(data = PC4_1.2sd_pre_upper_50, aes(age,PC4), colour = 'pink') +
  geom_line(data = PC4_1.2sd_pre_lower_50, aes(age,PC4), colour = 'pink') +
  geom_line(data = PC4_1.0sd_pre_upper_50, aes(age,PC4), colour = 'brown') +
  geom_line(data = PC4_1.0sd_pre_lower_50, aes(age,PC4), colour = 'brown') +
  geom_line(data = PC4_0.8sd_pre_upper_50, aes(age,PC4), colour = 'yellow') +
  geom_line(data = PC4_0.8sd_pre_lower_50, aes(age,PC4), colour = 'yellow') +
  geom_line(data = PC4_0.6sd_pre_upper_50, aes(age,PC4), colour = 'red') +
  geom_line(data = PC4_0.6sd_pre_lower_50, aes(age,PC4), colour = 'red') +
  geom_line(data = PC4_0.4sd_pre_upper_50, aes(age,PC4), colour = 'blue') +
  geom_line(data = PC4_0.4sd_pre_lower_50, aes(age,PC4), colour = 'blue') +
  geom_line(data = PC4_0.2sd_pre_upper_50, aes(age,PC4), colour = 'green') +
  geom_line(data = PC4_0.2sd_pre_lower_50, aes(age,PC4), colour = 'green') +
  geom_line(data = PC4_0.0sd_pre_upper_50, aes(age,PC4), colour = 'orange') +
  geom_line(data = PC4_0.0sd_pre_lower_50, aes(age,PC4), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC4)) + 
  #geom_point(data = PC4_2sd_pre_upper, aes(age, PC4), colour = 'red') +
  scale_x_reverse(limits = c(4000,540)) + theme_bw()

grid.arrange(b,a,ncol=2)








##### Phanerozoic PC5 #####

coords.all_Phan <- subset(coords.all2, age < 541)

# 2 sigma

nsig <- 2

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_2sd_phan_upper <- as.data.frame(cbind(PC5_2sd_phan$x, PC5_2sd_phan$upper))
names(PC5_2sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_2sd_phan_upper$age, PC5_2sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_2sd_phan_upper, indices) {
  d <- PC5_2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_2sd_phan_upper_50)[1] <- "PC5"

PC5_2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_2sd_phan_upper_all)[1:ncol(PC5_2sd_phan_upper_all)] <- seq(0,540, 1)
PC5_2sd_phan_upper_all$group <- seq(1,(nrow(PC5_2sd_phan_upper_all)), 1)

PC5_2sd_phan_upper_all <- reshape2::melt(PC5_2sd_phan_upper_all, id.vars = "group")
names(PC5_2sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_2sd_phan_lower <- as.data.frame(cbind(PC5_2sd_phan$x, PC5_2sd_phan$lower))
names(PC5_2sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_2sd_phan_lower$age, PC5_2sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_2sd_phan_lower, indices) {
  d <- PC5_2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_2sd_phan_lower_50)[1] <- "PC5"

PC5_2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_2sd_phan_lower_all)[1:ncol(PC5_2sd_phan_lower_all)] <- seq(0,540, 1)
PC5_2sd_phan_lower_all$group <- seq(1,(nrow(PC5_2sd_phan_lower_all)), 1)

PC5_2sd_phan_lower_all <- reshape2::melt(PC5_2sd_phan_lower_all, id.vars = "group")
names(PC5_2sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.8sd_phan_upper <- as.data.frame(cbind(PC5_1.8sd_phan$x, PC5_1.8sd_phan$upper))
names(PC5_1.8sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.8sd_phan_upper$age, PC5_1.8sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_1.8sd_phan_upper, indices) {
  d <- PC5_1.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_1.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_1.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_1.8sd_phan_upper_50)[1] <- "PC5"

PC5_1.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_1.8sd_phan_upper_all)[1:ncol(PC5_1.8sd_phan_upper_all)] <- seq(0,540, 1)
PC5_1.8sd_phan_upper_all$group <- seq(1,(nrow(PC5_1.8sd_phan_upper_all)), 1)

PC5_1.8sd_phan_upper_all <- reshape2::melt(PC5_1.8sd_phan_upper_all, id.vars = "group")
names(PC5_1.8sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_1.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.8sd_phan_lower <- as.data.frame(cbind(PC5_1.8sd_phan$x, PC5_1.8sd_phan$lower))
names(PC5_1.8sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.8sd_phan_lower$age, PC5_1.8sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_1.8sd_phan_lower, indices) {
  d <- PC5_1.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_1.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_1.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_1.8sd_phan_lower_50)[1] <- "PC5"

PC5_1.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_1.8sd_phan_lower_all)[1:ncol(PC5_1.8sd_phan_lower_all)] <- seq(0,540, 1)
PC5_1.8sd_phan_lower_all$group <- seq(1,(nrow(PC5_1.8sd_phan_lower_all)), 1)

PC5_1.8sd_phan_lower_all <- reshape2::melt(PC5_1.8sd_phan_lower_all, id.vars = "group")
names(PC5_1.8sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_1.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.6sd_phan_upper <- as.data.frame(cbind(PC5_1.6sd_phan$x, PC5_1.6sd_phan$upper))
names(PC5_1.6sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.6sd_phan_upper$age, PC5_1.6sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_1.6sd_phan_upper, indices) {
  d <- PC5_1.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_1.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_1.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_1.6sd_phan_upper_50)[1] <- "PC5"

PC5_1.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_1.6sd_phan_upper_all)[1:ncol(PC5_1.6sd_phan_upper_all)] <- seq(0,540, 1)
PC5_1.6sd_phan_upper_all$group <- seq(1,(nrow(PC5_1.6sd_phan_upper_all)), 1)

PC5_1.6sd_phan_upper_all <- reshape2::melt(PC5_1.6sd_phan_upper_all, id.vars = "group")
names(PC5_1.6sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_1.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.6sd_phan_lower <- as.data.frame(cbind(PC5_1.6sd_phan$x, PC5_1.6sd_phan$lower))
names(PC5_1.6sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.6sd_phan_lower$age, PC5_1.6sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_1.6sd_phan_lower, indices) {
  d <- PC5_1.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_1.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_1.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_1.6sd_phan_lower_50)[1] <- "PC5"

PC5_1.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_1.6sd_phan_lower_all)[1:ncol(PC5_1.6sd_phan_lower_all)] <- seq(0,540, 1)
PC5_1.6sd_phan_lower_all$group <- seq(1,(nrow(PC5_1.6sd_phan_lower_all)), 1)

PC5_1.6sd_phan_lower_all <- reshape2::melt(PC5_1.6sd_phan_lower_all, id.vars = "group")
names(PC5_1.6sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_1.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.4sd_phan_upper <- as.data.frame(cbind(PC5_1.4sd_phan$x, PC5_1.4sd_phan$upper))
names(PC5_1.4sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.4sd_phan_upper$age, PC5_1.4sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_1.4sd_phan_upper, indices) {
  d <- PC5_1.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_1.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_1.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_1.4sd_phan_upper_50)[1] <- "PC5"

PC5_1.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_1.4sd_phan_upper_all)[1:ncol(PC5_1.4sd_phan_upper_all)] <- seq(0,540, 1)
PC5_1.4sd_phan_upper_all$group <- seq(1,(nrow(PC5_1.4sd_phan_upper_all)), 1)

PC5_1.4sd_phan_upper_all <- reshape2::melt(PC5_1.4sd_phan_upper_all, id.vars = "group")
names(PC5_1.4sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_1.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.4sd_phan_lower <- as.data.frame(cbind(PC5_1.4sd_phan$x, PC5_1.4sd_phan$lower))
names(PC5_1.4sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.4sd_phan_lower$age, PC5_1.4sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_1.4sd_phan_lower, indices) {
  d <- PC5_1.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_1.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_1.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_1.4sd_phan_lower_50)[1] <- "PC5"

PC5_1.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_1.4sd_phan_lower_all)[1:ncol(PC5_1.4sd_phan_lower_all)] <- seq(0,540, 1)
PC5_1.4sd_phan_lower_all$group <- seq(1,(nrow(PC5_1.4sd_phan_lower_all)), 1)

PC5_1.4sd_phan_lower_all <- reshape2::melt(PC5_1.4sd_phan_lower_all, id.vars = "group")
names(PC5_1.4sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_1.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.2sd_phan_upper <- as.data.frame(cbind(PC5_1.2sd_phan$x, PC5_1.2sd_phan$upper))
names(PC5_1.2sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.2sd_phan_upper$age, PC5_1.2sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_1.2sd_phan_upper, indices) {
  d <- PC5_1.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_1.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_1.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_1.2sd_phan_upper_50)[1] <- "PC5"

PC5_1.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_1.2sd_phan_upper_all)[1:ncol(PC5_1.2sd_phan_upper_all)] <- seq(0,540, 1)
PC5_1.2sd_phan_upper_all$group <- seq(1,(nrow(PC5_1.2sd_phan_upper_all)), 1)

PC5_1.2sd_phan_upper_all <- reshape2::melt(PC5_1.2sd_phan_upper_all, id.vars = "group")
names(PC5_1.2sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_1.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.2sd_phan_lower <- as.data.frame(cbind(PC5_1.2sd_phan$x, PC5_1.2sd_phan$lower))
names(PC5_1.2sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.2sd_phan_lower$age, PC5_1.2sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_1.2sd_phan_lower, indices) {
  d <- PC5_1.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_1.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_1.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_1.2sd_phan_lower_50)[1] <- "PC5"

PC5_1.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_1.2sd_phan_lower_all)[1:ncol(PC5_1.2sd_phan_lower_all)] <- seq(0,540, 1)
PC5_1.2sd_phan_lower_all$group <- seq(1,(nrow(PC5_1.2sd_phan_lower_all)), 1)

PC5_1.2sd_phan_lower_all <- reshape2::melt(PC5_1.2sd_phan_lower_all, id.vars = "group")
names(PC5_1.2sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_1.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.0sd_phan_upper <- as.data.frame(cbind(PC5_1.0sd_phan$x, PC5_1.0sd_phan$upper))
names(PC5_1.0sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.0sd_phan_upper$age, PC5_1.0sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_1.0sd_phan_upper, indices) {
  d <- PC5_1.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_1.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_1.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_1.0sd_phan_upper_50)[1] <- "PC5"

PC5_1.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_1.0sd_phan_upper_all)[1:ncol(PC5_1.0sd_phan_upper_all)] <- seq(0,540, 1)
PC5_1.0sd_phan_upper_all$group <- seq(1,(nrow(PC5_1.0sd_phan_upper_all)), 1)

PC5_1.0sd_phan_upper_all <- reshape2::melt(PC5_1.0sd_phan_upper_all, id.vars = "group")
names(PC5_1.0sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_1.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.0sd_phan_lower <- as.data.frame(cbind(PC5_1.0sd_phan$x, PC5_1.0sd_phan$lower))
names(PC5_1.0sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.0sd_phan_lower$age, PC5_1.0sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_1.0sd_phan_lower, indices) {
  d <- PC5_1.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_1.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_1.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_1.0sd_phan_lower_50)[1] <- "PC5"

PC5_1.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_1.0sd_phan_lower_all)[1:ncol(PC5_1.0sd_phan_lower_all)] <- seq(0,540, 1)
PC5_1.0sd_phan_lower_all$group <- seq(1,(nrow(PC5_1.0sd_phan_lower_all)), 1)

PC5_1.0sd_phan_lower_all <- reshape2::melt(PC5_1.0sd_phan_lower_all, id.vars = "group")
names(PC5_1.0sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_1.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.8sd_phan_upper <- as.data.frame(cbind(PC5_0.8sd_phan$x, PC5_0.8sd_phan$upper))
names(PC5_0.8sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.8sd_phan_upper$age, PC5_0.8sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_0.8sd_phan_upper, indices) {
  d <- PC5_0.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_0.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_0.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_0.8sd_phan_upper_50)[1] <- "PC5"

PC5_0.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_0.8sd_phan_upper_all)[1:ncol(PC5_0.8sd_phan_upper_all)] <- seq(0,540, 1)
PC5_0.8sd_phan_upper_all$group <- seq(1,(nrow(PC5_0.8sd_phan_upper_all)), 1)

PC5_0.8sd_phan_upper_all <- reshape2::melt(PC5_0.8sd_phan_upper_all, id.vars = "group")
names(PC5_0.8sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_0.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.8sd_phan_lower <- as.data.frame(cbind(PC5_0.8sd_phan$x, PC5_0.8sd_phan$lower))
names(PC5_0.8sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.8sd_phan_lower$age, PC5_0.8sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_0.8sd_phan_lower, indices) {
  d <- PC5_0.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_0.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_0.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_0.8sd_phan_lower_50)[1] <- "PC5"

PC5_0.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_0.8sd_phan_lower_all)[1:ncol(PC5_0.8sd_phan_lower_all)] <- seq(0,540, 1)
PC5_0.8sd_phan_lower_all$group <- seq(1,(nrow(PC5_0.8sd_phan_lower_all)), 1)

PC5_0.8sd_phan_lower_all <- reshape2::melt(PC5_0.8sd_phan_lower_all, id.vars = "group")
names(PC5_0.8sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_0.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.6sd_phan_upper <- as.data.frame(cbind(PC5_0.6sd_phan$x, PC5_0.6sd_phan$upper))
names(PC5_0.6sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.6sd_phan_upper$age, PC5_0.6sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_0.6sd_phan_upper, indices) {
  d <- PC5_0.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_0.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_0.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_0.6sd_phan_upper_50)[1] <- "PC5"

PC5_0.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_0.6sd_phan_upper_all)[1:ncol(PC5_0.6sd_phan_upper_all)] <- seq(0,540, 1)
PC5_0.6sd_phan_upper_all$group <- seq(1,(nrow(PC5_0.6sd_phan_upper_all)), 1)

PC5_0.6sd_phan_upper_all <- reshape2::melt(PC5_0.6sd_phan_upper_all, id.vars = "group")
names(PC5_0.6sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_0.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.6sd_phan_lower <- as.data.frame(cbind(PC5_0.6sd_phan$x, PC5_0.6sd_phan$lower))
names(PC5_0.6sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.6sd_phan_lower$age, PC5_0.6sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_0.6sd_phan_lower, indices) {
  d <- PC5_0.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_0.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_0.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_0.6sd_phan_lower_50)[1] <- "PC5"

PC5_0.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_0.6sd_phan_lower_all)[1:ncol(PC5_0.6sd_phan_lower_all)] <- seq(0,540, 1)
PC5_0.6sd_phan_lower_all$group <- seq(1,(nrow(PC5_0.6sd_phan_lower_all)), 1)

PC5_0.6sd_phan_lower_all <- reshape2::melt(PC5_0.6sd_phan_lower_all, id.vars = "group")
names(PC5_0.6sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_0.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.4sd_phan_upper <- as.data.frame(cbind(PC5_0.4sd_phan$x, PC5_0.4sd_phan$upper))
names(PC5_0.4sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.4sd_phan_upper$age, PC5_0.4sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_0.4sd_phan_upper, indices) {
  d <- PC5_0.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_0.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_0.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_0.4sd_phan_upper_50)[1] <- "PC5"

PC5_0.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_0.4sd_phan_upper_all)[1:ncol(PC5_0.4sd_phan_upper_all)] <- seq(0,540, 1)
PC5_0.4sd_phan_upper_all$group <- seq(1,(nrow(PC5_0.4sd_phan_upper_all)), 1)

PC5_0.4sd_phan_upper_all <- reshape2::melt(PC5_0.4sd_phan_upper_all, id.vars = "group")
names(PC5_0.4sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_0.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.4sd_phan_lower <- as.data.frame(cbind(PC5_0.4sd_phan$x, PC5_0.4sd_phan$lower))
names(PC5_0.4sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.4sd_phan_lower$age, PC5_0.4sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_0.4sd_phan_lower, indices) {
  d <- PC5_0.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_0.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_0.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_0.4sd_phan_lower_50)[1] <- "PC5"

PC5_0.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_0.4sd_phan_lower_all)[1:ncol(PC5_0.4sd_phan_lower_all)] <- seq(0,540, 1)
PC5_0.4sd_phan_lower_all$group <- seq(1,(nrow(PC5_0.4sd_phan_lower_all)), 1)

PC5_0.4sd_phan_lower_all <- reshape2::melt(PC5_0.4sd_phan_lower_all, id.vars = "group")
names(PC5_0.4sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_0.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.2sd_phan_upper <- as.data.frame(cbind(PC5_0.2sd_phan$x, PC5_0.2sd_phan$upper))
names(PC5_0.2sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.2sd_phan_upper$age, PC5_0.2sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_0.2sd_phan_upper, indices) {
  d <- PC5_0.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_0.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_0.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_0.2sd_phan_upper_50)[1] <- "PC5"

PC5_0.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_0.2sd_phan_upper_all)[1:ncol(PC5_0.2sd_phan_upper_all)] <- seq(0,540, 1)
PC5_0.2sd_phan_upper_all$group <- seq(1,(nrow(PC5_0.2sd_phan_upper_all)), 1)

PC5_0.2sd_phan_upper_all <- reshape2::melt(PC5_0.2sd_phan_upper_all, id.vars = "group")
names(PC5_0.2sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_0.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.2sd_phan_lower <- as.data.frame(cbind(PC5_0.2sd_phan$x, PC5_0.2sd_phan$lower))
names(PC5_0.2sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.2sd_phan_lower$age, PC5_0.2sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_0.2sd_phan_lower, indices) {
  d <- PC5_0.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_0.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_0.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_0.2sd_phan_lower_50)[1] <- "PC5"

PC5_0.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_0.2sd_phan_lower_all)[1:ncol(PC5_0.2sd_phan_lower_all)] <- seq(0,540, 1)
PC5_0.2sd_phan_lower_all$group <- seq(1,(nrow(PC5_0.2sd_phan_lower_all)), 1)

PC5_0.2sd_phan_lower_all <- reshape2::melt(PC5_0.2sd_phan_lower_all, id.vars = "group")
names(PC5_0.2sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_0.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.0sd_phan_upper <- as.data.frame(cbind(PC5_0.0sd_phan$x, PC5_0.0sd_phan$upper))
names(PC5_0.0sd_phan_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.0sd_phan_upper$age, PC5_0.0sd_phan_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_0.0sd_phan_upper, indices) {
  d <- PC5_0.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_0.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_0.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC5_0.0sd_phan_upper_50)[1] <- "PC5"

PC5_0.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC5_0.0sd_phan_upper_all)[1:ncol(PC5_0.0sd_phan_upper_all)] <- seq(0,540, 1)
PC5_0.0sd_phan_upper_all$group <- seq(1,(nrow(PC5_0.0sd_phan_upper_all)), 1)

PC5_0.0sd_phan_upper_all <- reshape2::melt(PC5_0.0sd_phan_upper_all, id.vars = "group")
names(PC5_0.0sd_phan_upper_all)[2:3] <- c("age", "PC5")

PC5_0.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.0sd_phan_lower <- as.data.frame(cbind(PC5_0.0sd_phan$x, PC5_0.0sd_phan$lower))
names(PC5_0.0sd_phan_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.0sd_phan_lower$age, PC5_0.0sd_phan_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_0.0sd_phan_lower, indices) {
  d <- PC5_0.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_0.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_0.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC5_0.0sd_phan_lower_50)[1] <- "PC5"

PC5_0.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC5_0.0sd_phan_lower_all)[1:ncol(PC5_0.0sd_phan_lower_all)] <- seq(0,540, 1)
PC5_0.0sd_phan_lower_all$group <- seq(1,(nrow(PC5_0.0sd_phan_lower_all)), 1)

PC5_0.0sd_phan_lower_all <- reshape2::melt(PC5_0.0sd_phan_lower_all, id.vars = "group")
names(PC5_0.0sd_phan_lower_all)[2:3] <- c("age", "PC5")

PC5_0.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

a <- ggplot() + 
  #geom_line(data = PC5_2sd_phan_upper_all, aes(age2,PC5, group = group)) +
  #geom_line(data = PC5_2sd_phan_lower_all, aes(age2,PC5, group = group)) +
  geom_line(data = PC5_2sd_phan_upper_50, aes(age,PC5), colour = 'blue') +
  geom_line(data = PC5_2sd_phan_lower_50, aes(age,PC5), colour = 'blue') +
  geom_line(data = PC5_1.6sd_phan_upper_50, aes(age,PC5), colour = 'green') +
  geom_line(data = PC5_1.6sd_phan_lower_50, aes(age,PC5), colour = 'green') +
  geom_line(data = PC5_1.4sd_phan_upper_50, aes(age,PC5), colour = 'orange') +
  geom_line(data = PC5_1.4sd_phan_lower_50, aes(age,PC5), colour = 'orange') +
  geom_line(data = PC5_1.2sd_phan_upper_50, aes(age,PC5), colour = 'pink') +
  geom_line(data = PC5_1.2sd_phan_lower_50, aes(age,PC5), colour = 'pink') +
  geom_line(data = PC5_1.0sd_phan_upper_50, aes(age,PC5), colour = 'brown') +
  geom_line(data = PC5_1.0sd_phan_lower_50, aes(age,PC5), colour = 'brown') +
  geom_line(data = PC5_0.8sd_phan_upper_50, aes(age,PC5), colour = 'yellow') +
  geom_line(data = PC5_0.8sd_phan_lower_50, aes(age,PC5), colour = 'yellow') +
  geom_line(data = PC5_0.6sd_phan_upper_50, aes(age,PC5), colour = 'red') +
  geom_line(data = PC5_0.6sd_phan_lower_50, aes(age,PC5), colour = 'red') +
  geom_line(data = PC5_0.4sd_phan_upper_50, aes(age,PC5), colour = 'blue') +
  geom_line(data = PC5_0.4sd_phan_lower_50, aes(age,PC5), colour = 'blue') +
  geom_line(data = PC5_0.2sd_phan_upper_50, aes(age,PC5), colour = 'green') +
  geom_line(data = PC5_0.2sd_phan_lower_50, aes(age,PC5), colour = 'green') +
  geom_line(data = PC5_0.0sd_phan_upper_50, aes(age,PC5), colour = 'orange') +
  geom_line(data = PC5_0.0sd_phan_lower_50, aes(age,PC5), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC5)) + 
  #geom_point(data = PC5_2sd_phan_upper, aes(age, PC5), colour = 'red') +
  scale_x_reverse(limits = c(540,0)) + theme_bw()

##### Precambrian PC5 #####

coords.all_pre <- subset(coords.all2, age > 541)

# 2 sigma

nsig <- 2

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_2sd_pre_upper <- as.data.frame(cbind(PC5_2sd_pre$x, PC5_2sd_pre$upper))
names(PC5_2sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_2sd_pre_upper$age, PC5_2sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[1], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_2sd_pre_upper, indices) {
  d <- PC5_2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_2sd_pre_upper_50)[1] <- "PC5"

PC5_2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_2sd_pre_upper_all)[1:ncol(PC5_2sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_2sd_pre_upper_all$group <- seq(1,(nrow(PC5_2sd_pre_upper_all)), 1)

PC5_2sd_pre_upper_all <- reshape2::melt(PC5_2sd_pre_upper_all, id.vars = "group")
names(PC5_2sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_2sd_pre_lower <- as.data.frame(cbind(PC5_2sd_pre$x, PC5_2sd_pre$lower))
names(PC5_2sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_2sd_pre_lower$age, PC5_2sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_2sd_pre_lower, indices) {
  d <- PC5_2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_2sd_pre_lower_50)[1] <- "PC5"

PC5_2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_2sd_pre_lower_all)[1:ncol(PC5_2sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_2sd_pre_lower_all$group <- seq(1,(nrow(PC5_2sd_pre_lower_all)), 1)

PC5_2sd_pre_lower_all <- reshape2::melt(PC5_2sd_pre_lower_all, id.vars = "group")
names(PC5_2sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.8sd_pre_upper <- as.data.frame(cbind(PC5_1.8sd_pre$x, PC5_1.8sd_pre$upper))
names(PC5_1.8sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.8sd_pre_upper$age, PC5_1.8sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_1.8sd_pre_upper, indices) {
  d <- PC5_1.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_1.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_1.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_1.8sd_pre_upper_50)[1] <- "PC5"

PC5_1.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_1.8sd_pre_upper_all)[1:ncol(PC5_1.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_1.8sd_pre_upper_all$group <- seq(1,(nrow(PC5_1.8sd_pre_upper_all)), 1)

PC5_1.8sd_pre_upper_all <- reshape2::melt(PC5_1.8sd_pre_upper_all, id.vars = "group")
names(PC5_1.8sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_1.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.8sd_pre_lower <- as.data.frame(cbind(PC5_1.8sd_pre$x, PC5_1.8sd_pre$lower))
names(PC5_1.8sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.8sd_pre_lower$age, PC5_1.8sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_1.8sd_pre_lower, indices) {
  d <- PC5_1.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_1.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_1.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_1.8sd_pre_lower_50)[1] <- "PC5"

PC5_1.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_1.8sd_pre_lower_all)[1:ncol(PC5_1.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_1.8sd_pre_lower_all$group <- seq(1,(nrow(PC5_1.8sd_pre_lower_all)), 1)

PC5_1.8sd_pre_lower_all <- reshape2::melt(PC5_1.8sd_pre_lower_all, id.vars = "group")
names(PC5_1.8sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_1.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.6sd_pre_upper <- as.data.frame(cbind(PC5_1.6sd_pre$x, PC5_1.6sd_pre$upper))
names(PC5_1.6sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.6sd_pre_upper$age, PC5_1.6sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_1.6sd_pre_upper, indices) {
  d <- PC5_1.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_1.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_1.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_1.6sd_pre_upper_50)[1] <- "PC5"

PC5_1.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_1.6sd_pre_upper_all)[1:ncol(PC5_1.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_1.6sd_pre_upper_all$group <- seq(1,(nrow(PC5_1.6sd_pre_upper_all)), 1)

PC5_1.6sd_pre_upper_all <- reshape2::melt(PC5_1.6sd_pre_upper_all, id.vars = "group")
names(PC5_1.6sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_1.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.6sd_pre_lower <- as.data.frame(cbind(PC5_1.6sd_pre$x, PC5_1.6sd_pre$lower))
names(PC5_1.6sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.6sd_pre_lower$age, PC5_1.6sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_1.6sd_pre_lower, indices) {
  d <- PC5_1.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_1.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_1.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_1.6sd_pre_lower_50)[1] <- "PC5"

PC5_1.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_1.6sd_pre_lower_all)[1:ncol(PC5_1.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_1.6sd_pre_lower_all$group <- seq(1,(nrow(PC5_1.6sd_pre_lower_all)), 1)

PC5_1.6sd_pre_lower_all <- reshape2::melt(PC5_1.6sd_pre_lower_all, id.vars = "group")
names(PC5_1.6sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_1.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.4sd_pre_upper <- as.data.frame(cbind(PC5_1.4sd_pre$x, PC5_1.4sd_pre$upper))
names(PC5_1.4sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.4sd_pre_upper$age, PC5_1.4sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_1.4sd_pre_upper, indices) {
  d <- PC5_1.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_1.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_1.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_1.4sd_pre_upper_50)[1] <- "PC5"

PC5_1.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_1.4sd_pre_upper_all)[1:ncol(PC5_1.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_1.4sd_pre_upper_all$group <- seq(1,(nrow(PC5_1.4sd_pre_upper_all)), 1)

PC5_1.4sd_pre_upper_all <- reshape2::melt(PC5_1.4sd_pre_upper_all, id.vars = "group")
names(PC5_1.4sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_1.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.4sd_pre_lower <- as.data.frame(cbind(PC5_1.4sd_pre$x, PC5_1.4sd_pre$lower))
names(PC5_1.4sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.4sd_pre_lower$age, PC5_1.4sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_1.4sd_pre_lower, indices) {
  d <- PC5_1.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_1.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_1.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_1.4sd_pre_lower_50)[1] <- "PC5"

PC5_1.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_1.4sd_pre_lower_all)[1:ncol(PC5_1.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_1.4sd_pre_lower_all$group <- seq(1,(nrow(PC5_1.4sd_pre_lower_all)), 1)

PC5_1.4sd_pre_lower_all <- reshape2::melt(PC5_1.4sd_pre_lower_all, id.vars = "group")
names(PC5_1.4sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_1.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.2sd_pre_upper <- as.data.frame(cbind(PC5_1.2sd_pre$x, PC5_1.2sd_pre$upper))
names(PC5_1.2sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.2sd_pre_upper$age, PC5_1.2sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_1.2sd_pre_upper, indices) {
  d <- PC5_1.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_1.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_1.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_1.2sd_pre_upper_50)[1] <- "PC5"

PC5_1.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_1.2sd_pre_upper_all)[1:ncol(PC5_1.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_1.2sd_pre_upper_all$group <- seq(1,(nrow(PC5_1.2sd_pre_upper_all)), 1)

PC5_1.2sd_pre_upper_all <- reshape2::melt(PC5_1.2sd_pre_upper_all, id.vars = "group")
names(PC5_1.2sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_1.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.2sd_pre_lower <- as.data.frame(cbind(PC5_1.2sd_pre$x, PC5_1.2sd_pre$lower))
names(PC5_1.2sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.2sd_pre_lower$age, PC5_1.2sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_1.2sd_pre_lower, indices) {
  d <- PC5_1.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_1.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_1.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_1.2sd_pre_lower_50)[1] <- "PC5"

PC5_1.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_1.2sd_pre_lower_all)[1:ncol(PC5_1.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_1.2sd_pre_lower_all$group <- seq(1,(nrow(PC5_1.2sd_pre_lower_all)), 1)

PC5_1.2sd_pre_lower_all <- reshape2::melt(PC5_1.2sd_pre_lower_all, id.vars = "group")
names(PC5_1.2sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_1.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.0sd_pre_upper <- as.data.frame(cbind(PC5_1.0sd_pre$x, PC5_1.0sd_pre$upper))
names(PC5_1.0sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.0sd_pre_upper$age, PC5_1.0sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_1.0sd_pre_upper, indices) {
  d <- PC5_1.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_1.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_1.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_1.0sd_pre_upper_50)[1] <- "PC5"

PC5_1.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_1.0sd_pre_upper_all)[1:ncol(PC5_1.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_1.0sd_pre_upper_all$group <- seq(1,(nrow(PC5_1.0sd_pre_upper_all)), 1)

PC5_1.0sd_pre_upper_all <- reshape2::melt(PC5_1.0sd_pre_upper_all, id.vars = "group")
names(PC5_1.0sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_1.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_1.0sd_pre_lower <- as.data.frame(cbind(PC5_1.0sd_pre$x, PC5_1.0sd_pre$lower))
names(PC5_1.0sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_1.0sd_pre_lower$age, PC5_1.0sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_1.0sd_pre_lower, indices) {
  d <- PC5_1.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_1.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_1.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_1.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_1.0sd_pre_lower_50)[1] <- "PC5"

PC5_1.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_1.0sd_pre_lower_all)[1:ncol(PC5_1.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_1.0sd_pre_lower_all$group <- seq(1,(nrow(PC5_1.0sd_pre_lower_all)), 1)

PC5_1.0sd_pre_lower_all <- reshape2::melt(PC5_1.0sd_pre_lower_all, id.vars = "group")
names(PC5_1.0sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_1.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.8sd_pre_upper <- as.data.frame(cbind(PC5_0.8sd_pre$x, PC5_0.8sd_pre$upper))
names(PC5_0.8sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.8sd_pre_upper$age, PC5_0.8sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_0.8sd_pre_upper, indices) {
  d <- PC5_0.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_0.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_0.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_0.8sd_pre_upper_50)[1] <- "PC5"

PC5_0.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_0.8sd_pre_upper_all)[1:ncol(PC5_0.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_0.8sd_pre_upper_all$group <- seq(1,(nrow(PC5_0.8sd_pre_upper_all)), 1)

PC5_0.8sd_pre_upper_all <- reshape2::melt(PC5_0.8sd_pre_upper_all, id.vars = "group")
names(PC5_0.8sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_0.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.8sd_pre_lower <- as.data.frame(cbind(PC5_0.8sd_pre$x, PC5_0.8sd_pre$lower))
names(PC5_0.8sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.8sd_pre_lower$age, PC5_0.8sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_0.8sd_pre_lower, indices) {
  d <- PC5_0.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_0.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_0.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_0.8sd_pre_lower_50)[1] <- "PC5"

PC5_0.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_0.8sd_pre_lower_all)[1:ncol(PC5_0.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_0.8sd_pre_lower_all$group <- seq(1,(nrow(PC5_0.8sd_pre_lower_all)), 1)

PC5_0.8sd_pre_lower_all <- reshape2::melt(PC5_0.8sd_pre_lower_all, id.vars = "group")
names(PC5_0.8sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_0.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.6sd_pre_upper <- as.data.frame(cbind(PC5_0.6sd_pre$x, PC5_0.6sd_pre$upper))
names(PC5_0.6sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.6sd_pre_upper$age, PC5_0.6sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_0.6sd_pre_upper, indices) {
  d <- PC5_0.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_0.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_0.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_0.6sd_pre_upper_50)[1] <- "PC5"

PC5_0.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_0.6sd_pre_upper_all)[1:ncol(PC5_0.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_0.6sd_pre_upper_all$group <- seq(1,(nrow(PC5_0.6sd_pre_upper_all)), 1)

PC5_0.6sd_pre_upper_all <- reshape2::melt(PC5_0.6sd_pre_upper_all, id.vars = "group")
names(PC5_0.6sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_0.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.6sd_pre_lower <- as.data.frame(cbind(PC5_0.6sd_pre$x, PC5_0.6sd_pre$lower))
names(PC5_0.6sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.6sd_pre_lower$age, PC5_0.6sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_0.6sd_pre_lower, indices) {
  d <- PC5_0.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_0.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_0.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_0.6sd_pre_lower_50)[1] <- "PC5"

PC5_0.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_0.6sd_pre_lower_all)[1:ncol(PC5_0.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_0.6sd_pre_lower_all$group <- seq(1,(nrow(PC5_0.6sd_pre_lower_all)), 1)

PC5_0.6sd_pre_lower_all <- reshape2::melt(PC5_0.6sd_pre_lower_all, id.vars = "group")
names(PC5_0.6sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_0.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.4sd_pre_upper <- as.data.frame(cbind(PC5_0.4sd_pre$x, PC5_0.4sd_pre$upper))
names(PC5_0.4sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.4sd_pre_upper$age, PC5_0.4sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_0.4sd_pre_upper, indices) {
  d <- PC5_0.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_0.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_0.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_0.4sd_pre_upper_50)[1] <- "PC5"

PC5_0.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_0.4sd_pre_upper_all)[1:ncol(PC5_0.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_0.4sd_pre_upper_all$group <- seq(1,(nrow(PC5_0.4sd_pre_upper_all)), 1)

PC5_0.4sd_pre_upper_all <- reshape2::melt(PC5_0.4sd_pre_upper_all, id.vars = "group")
names(PC5_0.4sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_0.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.4sd_pre_lower <- as.data.frame(cbind(PC5_0.4sd_pre$x, PC5_0.4sd_pre$lower))
names(PC5_0.4sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.4sd_pre_lower$age, PC5_0.4sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_0.4sd_pre_lower, indices) {
  d <- PC5_0.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_0.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_0.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_0.4sd_pre_lower_50)[1] <- "PC5"

PC5_0.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_0.4sd_pre_lower_all)[1:ncol(PC5_0.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_0.4sd_pre_lower_all$group <- seq(1,(nrow(PC5_0.4sd_pre_lower_all)), 1)

PC5_0.4sd_pre_lower_all <- reshape2::melt(PC5_0.4sd_pre_lower_all, id.vars = "group")
names(PC5_0.4sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_0.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.2sd_pre_upper <- as.data.frame(cbind(PC5_0.2sd_pre$x, PC5_0.2sd_pre$upper))
names(PC5_0.2sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.2sd_pre_upper$age, PC5_0.2sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_0.2sd_pre_upper, indices) {
  d <- PC5_0.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_0.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_0.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_0.2sd_pre_upper_50)[1] <- "PC5"

PC5_0.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_0.2sd_pre_upper_all)[1:ncol(PC5_0.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_0.2sd_pre_upper_all$group <- seq(1,(nrow(PC5_0.2sd_pre_upper_all)), 1)

PC5_0.2sd_pre_upper_all <- reshape2::melt(PC5_0.2sd_pre_upper_all, id.vars = "group")
names(PC5_0.2sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_0.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.2sd_pre_lower <- as.data.frame(cbind(PC5_0.2sd_pre$x, PC5_0.2sd_pre$lower))
names(PC5_0.2sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.2sd_pre_lower$age, PC5_0.2sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_0.2sd_pre_lower, indices) {
  d <- PC5_0.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_0.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_0.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_0.2sd_pre_lower_50)[1] <- "PC5"

PC5_0.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_0.2sd_pre_lower_all)[1:ncol(PC5_0.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_0.2sd_pre_lower_all$group <- seq(1,(nrow(PC5_0.2sd_pre_lower_all)), 1)

PC5_0.2sd_pre_lower_all <- reshape2::melt(PC5_0.2sd_pre_lower_all, id.vars = "group")
names(PC5_0.2sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_0.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC5 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.0sd_pre_upper <- as.data.frame(cbind(PC5_0.0sd_pre$x, PC5_0.0sd_pre$upper))
names(PC5_0.0sd_pre_upper)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.0sd_pre_upper$age, PC5_0.0sd_pre_upper$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC5_0.0sd_pre_upper, indices) {
  d <- PC5_0.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC5_0.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC5_0.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC5_0.0sd_pre_upper_50)[1] <- "PC5"

PC5_0.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC5_0.0sd_pre_upper_all)[1:ncol(PC5_0.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC5_0.0sd_pre_upper_all$group <- seq(1,(nrow(PC5_0.0sd_pre_upper_all)), 1)

PC5_0.0sd_pre_upper_all <- reshape2::melt(PC5_0.0sd_pre_upper_all, id.vars = "group")
names(PC5_0.0sd_pre_upper_all)[2:3] <- c("age", "PC5")

PC5_0.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC5 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC5_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC5, degree=1, nsigma = nsig, span = span)
PC5_0.0sd_pre_lower <- as.data.frame(cbind(PC5_0.0sd_pre$x, PC5_0.0sd_pre$lower))
names(PC5_0.0sd_pre_lower)[1:2] <- c("age","PC5")

# CV spans
loess.predict <- loess.as(PC5_0.0sd_pre_lower$age, PC5_0.0sd_pre_lower$PC5, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC5_0.0sd_pre_lower, indices) {
  d <- PC5_0.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC5 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC5_0.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC5_0.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC5_0.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC5_0.0sd_pre_lower_50)[1] <- "PC5"

PC5_0.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC5_0.0sd_pre_lower_all)[1:ncol(PC5_0.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC5_0.0sd_pre_lower_all$group <- seq(1,(nrow(PC5_0.0sd_pre_lower_all)), 1)

PC5_0.0sd_pre_lower_all <- reshape2::melt(PC5_0.0sd_pre_lower_all, id.vars = "group")
names(PC5_0.0sd_pre_lower_all)[2:3] <- c("age", "PC5")

PC5_0.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

b <- ggplot() + 
  #geom_line(data = PC5_2sd_pre_upper_all, aes(age2,PC5, group = group)) +
  #geom_line(data = PC5_2sd_pre_lower_all, aes(age2,PC5, group = group)) +
  geom_line(data = PC5_2sd_pre_upper_50, aes(age,PC5), colour = 'blue') +
  geom_line(data = PC5_2sd_pre_lower_50, aes(age,PC5), colour = 'blue') +
  geom_line(data = PC5_1.6sd_pre_upper_50, aes(age,PC5), colour = 'green') +
  geom_line(data = PC5_1.6sd_pre_lower_50, aes(age,PC5), colour = 'green') +
  geom_line(data = PC5_1.4sd_pre_upper_50, aes(age,PC5), colour = 'orange') +
  geom_line(data = PC5_1.4sd_pre_lower_50, aes(age,PC5), colour = 'orange') +
  geom_line(data = PC5_1.2sd_pre_upper_50, aes(age,PC5), colour = 'pink') +
  geom_line(data = PC5_1.2sd_pre_lower_50, aes(age,PC5), colour = 'pink') +
  geom_line(data = PC5_1.0sd_pre_upper_50, aes(age,PC5), colour = 'brown') +
  geom_line(data = PC5_1.0sd_pre_lower_50, aes(age,PC5), colour = 'brown') +
  geom_line(data = PC5_0.8sd_pre_upper_50, aes(age,PC5), colour = 'yellow') +
  geom_line(data = PC5_0.8sd_pre_lower_50, aes(age,PC5), colour = 'yellow') +
  geom_line(data = PC5_0.6sd_pre_upper_50, aes(age,PC5), colour = 'red') +
  geom_line(data = PC5_0.6sd_pre_lower_50, aes(age,PC5), colour = 'red') +
  geom_line(data = PC5_0.4sd_pre_upper_50, aes(age,PC5), colour = 'blue') +
  geom_line(data = PC5_0.4sd_pre_lower_50, aes(age,PC5), colour = 'blue') +
  geom_line(data = PC5_0.2sd_pre_upper_50, aes(age,PC5), colour = 'green') +
  geom_line(data = PC5_0.2sd_pre_lower_50, aes(age,PC5), colour = 'green') +
  geom_line(data = PC5_0.0sd_pre_upper_50, aes(age,PC5), colour = 'orange') +
  geom_line(data = PC5_0.0sd_pre_lower_50, aes(age,PC5), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC5)) + 
  #geom_point(data = PC5_2sd_pre_upper, aes(age, PC5), colour = 'red') +
  scale_x_reverse(limits = c(4000,540)) + theme_bw()

grid.arrange(b,a,ncol=2)










##### Phanerozoic PC6 #####

coords.all_Phan <- subset(coords.all2, age < 541)

# 2 sigma

nsig <- 2

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_2sd_phan_upper <- as.data.frame(cbind(PC6_2sd_phan$x, PC6_2sd_phan$upper))
names(PC6_2sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_2sd_phan_upper$age, PC6_2sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_2sd_phan_upper, indices) {
  d <- PC6_2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_2sd_phan_upper_50)[1] <- "PC6"

PC6_2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_2sd_phan_upper_all)[1:ncol(PC6_2sd_phan_upper_all)] <- seq(0,540, 1)
PC6_2sd_phan_upper_all$group <- seq(1,(nrow(PC6_2sd_phan_upper_all)), 1)

PC6_2sd_phan_upper_all <- reshape2::melt(PC6_2sd_phan_upper_all, id.vars = "group")
names(PC6_2sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_2sd_phan_lower <- as.data.frame(cbind(PC6_2sd_phan$x, PC6_2sd_phan$lower))
names(PC6_2sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_2sd_phan_lower$age, PC6_2sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_2sd_phan_lower, indices) {
  d <- PC6_2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_2sd_phan_lower_50)[1] <- "PC6"

PC6_2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_2sd_phan_lower_all)[1:ncol(PC6_2sd_phan_lower_all)] <- seq(0,540, 1)
PC6_2sd_phan_lower_all$group <- seq(1,(nrow(PC6_2sd_phan_lower_all)), 1)

PC6_2sd_phan_lower_all <- reshape2::melt(PC6_2sd_phan_lower_all, id.vars = "group")
names(PC6_2sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.8sd_phan_upper <- as.data.frame(cbind(PC6_1.8sd_phan$x, PC6_1.8sd_phan$upper))
names(PC6_1.8sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.8sd_phan_upper$age, PC6_1.8sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_1.8sd_phan_upper, indices) {
  d <- PC6_1.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_1.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_1.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_1.8sd_phan_upper_50)[1] <- "PC6"

PC6_1.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_1.8sd_phan_upper_all)[1:ncol(PC6_1.8sd_phan_upper_all)] <- seq(0,540, 1)
PC6_1.8sd_phan_upper_all$group <- seq(1,(nrow(PC6_1.8sd_phan_upper_all)), 1)

PC6_1.8sd_phan_upper_all <- reshape2::melt(PC6_1.8sd_phan_upper_all, id.vars = "group")
names(PC6_1.8sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_1.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.8sd_phan_lower <- as.data.frame(cbind(PC6_1.8sd_phan$x, PC6_1.8sd_phan$lower))
names(PC6_1.8sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.8sd_phan_lower$age, PC6_1.8sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_1.8sd_phan_lower, indices) {
  d <- PC6_1.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_1.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_1.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_1.8sd_phan_lower_50)[1] <- "PC6"

PC6_1.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_1.8sd_phan_lower_all)[1:ncol(PC6_1.8sd_phan_lower_all)] <- seq(0,540, 1)
PC6_1.8sd_phan_lower_all$group <- seq(1,(nrow(PC6_1.8sd_phan_lower_all)), 1)

PC6_1.8sd_phan_lower_all <- reshape2::melt(PC6_1.8sd_phan_lower_all, id.vars = "group")
names(PC6_1.8sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_1.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.6sd_phan_upper <- as.data.frame(cbind(PC6_1.6sd_phan$x, PC6_1.6sd_phan$upper))
names(PC6_1.6sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.6sd_phan_upper$age, PC6_1.6sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_1.6sd_phan_upper, indices) {
  d <- PC6_1.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_1.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_1.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_1.6sd_phan_upper_50)[1] <- "PC6"

PC6_1.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_1.6sd_phan_upper_all)[1:ncol(PC6_1.6sd_phan_upper_all)] <- seq(0,540, 1)
PC6_1.6sd_phan_upper_all$group <- seq(1,(nrow(PC6_1.6sd_phan_upper_all)), 1)

PC6_1.6sd_phan_upper_all <- reshape2::melt(PC6_1.6sd_phan_upper_all, id.vars = "group")
names(PC6_1.6sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_1.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.6sd_phan_lower <- as.data.frame(cbind(PC6_1.6sd_phan$x, PC6_1.6sd_phan$lower))
names(PC6_1.6sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.6sd_phan_lower$age, PC6_1.6sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_1.6sd_phan_lower, indices) {
  d <- PC6_1.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_1.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_1.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_1.6sd_phan_lower_50)[1] <- "PC6"

PC6_1.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_1.6sd_phan_lower_all)[1:ncol(PC6_1.6sd_phan_lower_all)] <- seq(0,540, 1)
PC6_1.6sd_phan_lower_all$group <- seq(1,(nrow(PC6_1.6sd_phan_lower_all)), 1)

PC6_1.6sd_phan_lower_all <- reshape2::melt(PC6_1.6sd_phan_lower_all, id.vars = "group")
names(PC6_1.6sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_1.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.4sd_phan_upper <- as.data.frame(cbind(PC6_1.4sd_phan$x, PC6_1.4sd_phan$upper))
names(PC6_1.4sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.4sd_phan_upper$age, PC6_1.4sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_1.4sd_phan_upper, indices) {
  d <- PC6_1.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_1.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_1.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_1.4sd_phan_upper_50)[1] <- "PC6"

PC6_1.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_1.4sd_phan_upper_all)[1:ncol(PC6_1.4sd_phan_upper_all)] <- seq(0,540, 1)
PC6_1.4sd_phan_upper_all$group <- seq(1,(nrow(PC6_1.4sd_phan_upper_all)), 1)

PC6_1.4sd_phan_upper_all <- reshape2::melt(PC6_1.4sd_phan_upper_all, id.vars = "group")
names(PC6_1.4sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_1.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.4sd_phan_lower <- as.data.frame(cbind(PC6_1.4sd_phan$x, PC6_1.4sd_phan$lower))
names(PC6_1.4sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.4sd_phan_lower$age, PC6_1.4sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_1.4sd_phan_lower, indices) {
  d <- PC6_1.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_1.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_1.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_1.4sd_phan_lower_50)[1] <- "PC6"

PC6_1.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_1.4sd_phan_lower_all)[1:ncol(PC6_1.4sd_phan_lower_all)] <- seq(0,540, 1)
PC6_1.4sd_phan_lower_all$group <- seq(1,(nrow(PC6_1.4sd_phan_lower_all)), 1)

PC6_1.4sd_phan_lower_all <- reshape2::melt(PC6_1.4sd_phan_lower_all, id.vars = "group")
names(PC6_1.4sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_1.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.2sd_phan_upper <- as.data.frame(cbind(PC6_1.2sd_phan$x, PC6_1.2sd_phan$upper))
names(PC6_1.2sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.2sd_phan_upper$age, PC6_1.2sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_1.2sd_phan_upper, indices) {
  d <- PC6_1.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_1.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_1.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_1.2sd_phan_upper_50)[1] <- "PC6"

PC6_1.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_1.2sd_phan_upper_all)[1:ncol(PC6_1.2sd_phan_upper_all)] <- seq(0,540, 1)
PC6_1.2sd_phan_upper_all$group <- seq(1,(nrow(PC6_1.2sd_phan_upper_all)), 1)

PC6_1.2sd_phan_upper_all <- reshape2::melt(PC6_1.2sd_phan_upper_all, id.vars = "group")
names(PC6_1.2sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_1.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.2sd_phan_lower <- as.data.frame(cbind(PC6_1.2sd_phan$x, PC6_1.2sd_phan$lower))
names(PC6_1.2sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.2sd_phan_lower$age, PC6_1.2sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_1.2sd_phan_lower, indices) {
  d <- PC6_1.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_1.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_1.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_1.2sd_phan_lower_50)[1] <- "PC6"

PC6_1.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_1.2sd_phan_lower_all)[1:ncol(PC6_1.2sd_phan_lower_all)] <- seq(0,540, 1)
PC6_1.2sd_phan_lower_all$group <- seq(1,(nrow(PC6_1.2sd_phan_lower_all)), 1)

PC6_1.2sd_phan_lower_all <- reshape2::melt(PC6_1.2sd_phan_lower_all, id.vars = "group")
names(PC6_1.2sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_1.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.0sd_phan_upper <- as.data.frame(cbind(PC6_1.0sd_phan$x, PC6_1.0sd_phan$upper))
names(PC6_1.0sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.0sd_phan_upper$age, PC6_1.0sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_1.0sd_phan_upper, indices) {
  d <- PC6_1.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_1.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_1.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_1.0sd_phan_upper_50)[1] <- "PC6"

PC6_1.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_1.0sd_phan_upper_all)[1:ncol(PC6_1.0sd_phan_upper_all)] <- seq(0,540, 1)
PC6_1.0sd_phan_upper_all$group <- seq(1,(nrow(PC6_1.0sd_phan_upper_all)), 1)

PC6_1.0sd_phan_upper_all <- reshape2::melt(PC6_1.0sd_phan_upper_all, id.vars = "group")
names(PC6_1.0sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_1.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.0sd_phan_lower <- as.data.frame(cbind(PC6_1.0sd_phan$x, PC6_1.0sd_phan$lower))
names(PC6_1.0sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.0sd_phan_lower$age, PC6_1.0sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_1.0sd_phan_lower, indices) {
  d <- PC6_1.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_1.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_1.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_1.0sd_phan_lower_50)[1] <- "PC6"

PC6_1.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_1.0sd_phan_lower_all)[1:ncol(PC6_1.0sd_phan_lower_all)] <- seq(0,540, 1)
PC6_1.0sd_phan_lower_all$group <- seq(1,(nrow(PC6_1.0sd_phan_lower_all)), 1)

PC6_1.0sd_phan_lower_all <- reshape2::melt(PC6_1.0sd_phan_lower_all, id.vars = "group")
names(PC6_1.0sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_1.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.8sd_phan_upper <- as.data.frame(cbind(PC6_0.8sd_phan$x, PC6_0.8sd_phan$upper))
names(PC6_0.8sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.8sd_phan_upper$age, PC6_0.8sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_0.8sd_phan_upper, indices) {
  d <- PC6_0.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_0.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_0.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_0.8sd_phan_upper_50)[1] <- "PC6"

PC6_0.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_0.8sd_phan_upper_all)[1:ncol(PC6_0.8sd_phan_upper_all)] <- seq(0,540, 1)
PC6_0.8sd_phan_upper_all$group <- seq(1,(nrow(PC6_0.8sd_phan_upper_all)), 1)

PC6_0.8sd_phan_upper_all <- reshape2::melt(PC6_0.8sd_phan_upper_all, id.vars = "group")
names(PC6_0.8sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_0.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.8sd_phan_lower <- as.data.frame(cbind(PC6_0.8sd_phan$x, PC6_0.8sd_phan$lower))
names(PC6_0.8sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.8sd_phan_lower$age, PC6_0.8sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_0.8sd_phan_lower, indices) {
  d <- PC6_0.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_0.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_0.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_0.8sd_phan_lower_50)[1] <- "PC6"

PC6_0.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_0.8sd_phan_lower_all)[1:ncol(PC6_0.8sd_phan_lower_all)] <- seq(0,540, 1)
PC6_0.8sd_phan_lower_all$group <- seq(1,(nrow(PC6_0.8sd_phan_lower_all)), 1)

PC6_0.8sd_phan_lower_all <- reshape2::melt(PC6_0.8sd_phan_lower_all, id.vars = "group")
names(PC6_0.8sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_0.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.6sd_phan_upper <- as.data.frame(cbind(PC6_0.6sd_phan$x, PC6_0.6sd_phan$upper))
names(PC6_0.6sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.6sd_phan_upper$age, PC6_0.6sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_0.6sd_phan_upper, indices) {
  d <- PC6_0.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_0.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_0.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_0.6sd_phan_upper_50)[1] <- "PC6"

PC6_0.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_0.6sd_phan_upper_all)[1:ncol(PC6_0.6sd_phan_upper_all)] <- seq(0,540, 1)
PC6_0.6sd_phan_upper_all$group <- seq(1,(nrow(PC6_0.6sd_phan_upper_all)), 1)

PC6_0.6sd_phan_upper_all <- reshape2::melt(PC6_0.6sd_phan_upper_all, id.vars = "group")
names(PC6_0.6sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_0.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.6sd_phan_lower <- as.data.frame(cbind(PC6_0.6sd_phan$x, PC6_0.6sd_phan$lower))
names(PC6_0.6sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.6sd_phan_lower$age, PC6_0.6sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_0.6sd_phan_lower, indices) {
  d <- PC6_0.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_0.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_0.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_0.6sd_phan_lower_50)[1] <- "PC6"

PC6_0.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_0.6sd_phan_lower_all)[1:ncol(PC6_0.6sd_phan_lower_all)] <- seq(0,540, 1)
PC6_0.6sd_phan_lower_all$group <- seq(1,(nrow(PC6_0.6sd_phan_lower_all)), 1)

PC6_0.6sd_phan_lower_all <- reshape2::melt(PC6_0.6sd_phan_lower_all, id.vars = "group")
names(PC6_0.6sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_0.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.4sd_phan_upper <- as.data.frame(cbind(PC6_0.4sd_phan$x, PC6_0.4sd_phan$upper))
names(PC6_0.4sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.4sd_phan_upper$age, PC6_0.4sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_0.4sd_phan_upper, indices) {
  d <- PC6_0.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_0.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_0.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_0.4sd_phan_upper_50)[1] <- "PC6"

PC6_0.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_0.4sd_phan_upper_all)[1:ncol(PC6_0.4sd_phan_upper_all)] <- seq(0,540, 1)
PC6_0.4sd_phan_upper_all$group <- seq(1,(nrow(PC6_0.4sd_phan_upper_all)), 1)

PC6_0.4sd_phan_upper_all <- reshape2::melt(PC6_0.4sd_phan_upper_all, id.vars = "group")
names(PC6_0.4sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_0.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.4sd_phan_lower <- as.data.frame(cbind(PC6_0.4sd_phan$x, PC6_0.4sd_phan$lower))
names(PC6_0.4sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.4sd_phan_lower$age, PC6_0.4sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_0.4sd_phan_lower, indices) {
  d <- PC6_0.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_0.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_0.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_0.4sd_phan_lower_50)[1] <- "PC6"

PC6_0.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_0.4sd_phan_lower_all)[1:ncol(PC6_0.4sd_phan_lower_all)] <- seq(0,540, 1)
PC6_0.4sd_phan_lower_all$group <- seq(1,(nrow(PC6_0.4sd_phan_lower_all)), 1)

PC6_0.4sd_phan_lower_all <- reshape2::melt(PC6_0.4sd_phan_lower_all, id.vars = "group")
names(PC6_0.4sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_0.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.2sd_phan_upper <- as.data.frame(cbind(PC6_0.2sd_phan$x, PC6_0.2sd_phan$upper))
names(PC6_0.2sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.2sd_phan_upper$age, PC6_0.2sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_0.2sd_phan_upper, indices) {
  d <- PC6_0.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_0.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_0.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_0.2sd_phan_upper_50)[1] <- "PC6"

PC6_0.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_0.2sd_phan_upper_all)[1:ncol(PC6_0.2sd_phan_upper_all)] <- seq(0,540, 1)
PC6_0.2sd_phan_upper_all$group <- seq(1,(nrow(PC6_0.2sd_phan_upper_all)), 1)

PC6_0.2sd_phan_upper_all <- reshape2::melt(PC6_0.2sd_phan_upper_all, id.vars = "group")
names(PC6_0.2sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_0.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.2sd_phan_lower <- as.data.frame(cbind(PC6_0.2sd_phan$x, PC6_0.2sd_phan$lower))
names(PC6_0.2sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.2sd_phan_lower$age, PC6_0.2sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_0.2sd_phan_lower, indices) {
  d <- PC6_0.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_0.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_0.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_0.2sd_phan_lower_50)[1] <- "PC6"

PC6_0.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_0.2sd_phan_lower_all)[1:ncol(PC6_0.2sd_phan_lower_all)] <- seq(0,540, 1)
PC6_0.2sd_phan_lower_all$group <- seq(1,(nrow(PC6_0.2sd_phan_lower_all)), 1)

PC6_0.2sd_phan_lower_all <- reshape2::melt(PC6_0.2sd_phan_lower_all, id.vars = "group")
names(PC6_0.2sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_0.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.0sd_phan_upper <- as.data.frame(cbind(PC6_0.0sd_phan$x, PC6_0.0sd_phan$upper))
names(PC6_0.0sd_phan_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.0sd_phan_upper$age, PC6_0.0sd_phan_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_0.0sd_phan_upper, indices) {
  d <- PC6_0.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_0.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_0.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC6_0.0sd_phan_upper_50)[1] <- "PC6"

PC6_0.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC6_0.0sd_phan_upper_all)[1:ncol(PC6_0.0sd_phan_upper_all)] <- seq(0,540, 1)
PC6_0.0sd_phan_upper_all$group <- seq(1,(nrow(PC6_0.0sd_phan_upper_all)), 1)

PC6_0.0sd_phan_upper_all <- reshape2::melt(PC6_0.0sd_phan_upper_all, id.vars = "group")
names(PC6_0.0sd_phan_upper_all)[2:3] <- c("age", "PC6")

PC6_0.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.0sd_phan_lower <- as.data.frame(cbind(PC6_0.0sd_phan$x, PC6_0.0sd_phan$lower))
names(PC6_0.0sd_phan_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.0sd_phan_lower$age, PC6_0.0sd_phan_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_0.0sd_phan_lower, indices) {
  d <- PC6_0.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_0.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_0.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC6_0.0sd_phan_lower_50)[1] <- "PC6"

PC6_0.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC6_0.0sd_phan_lower_all)[1:ncol(PC6_0.0sd_phan_lower_all)] <- seq(0,540, 1)
PC6_0.0sd_phan_lower_all$group <- seq(1,(nrow(PC6_0.0sd_phan_lower_all)), 1)

PC6_0.0sd_phan_lower_all <- reshape2::melt(PC6_0.0sd_phan_lower_all, id.vars = "group")
names(PC6_0.0sd_phan_lower_all)[2:3] <- c("age", "PC6")

PC6_0.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

a <- ggplot() + 
  #geom_line(data = PC6_2sd_phan_upper_all, aes(age2,PC6, group = group)) +
  #geom_line(data = PC6_2sd_phan_lower_all, aes(age2,PC6, group = group)) +
  geom_line(data = PC6_2sd_phan_upper_50, aes(age,PC6), colour = 'blue') +
  geom_line(data = PC6_2sd_phan_lower_50, aes(age,PC6), colour = 'blue') +
  geom_line(data = PC6_1.6sd_phan_upper_50, aes(age,PC6), colour = 'green') +
  geom_line(data = PC6_1.6sd_phan_lower_50, aes(age,PC6), colour = 'green') +
  geom_line(data = PC6_1.4sd_phan_upper_50, aes(age,PC6), colour = 'orange') +
  geom_line(data = PC6_1.4sd_phan_lower_50, aes(age,PC6), colour = 'orange') +
  geom_line(data = PC6_1.2sd_phan_upper_50, aes(age,PC6), colour = 'pink') +
  geom_line(data = PC6_1.2sd_phan_lower_50, aes(age,PC6), colour = 'pink') +
  geom_line(data = PC6_1.0sd_phan_upper_50, aes(age,PC6), colour = 'brown') +
  geom_line(data = PC6_1.0sd_phan_lower_50, aes(age,PC6), colour = 'brown') +
  geom_line(data = PC6_0.8sd_phan_upper_50, aes(age,PC6), colour = 'yellow') +
  geom_line(data = PC6_0.8sd_phan_lower_50, aes(age,PC6), colour = 'yellow') +
  geom_line(data = PC6_0.6sd_phan_upper_50, aes(age,PC6), colour = 'red') +
  geom_line(data = PC6_0.6sd_phan_lower_50, aes(age,PC6), colour = 'red') +
  geom_line(data = PC6_0.4sd_phan_upper_50, aes(age,PC6), colour = 'blue') +
  geom_line(data = PC6_0.4sd_phan_lower_50, aes(age,PC6), colour = 'blue') +
  geom_line(data = PC6_0.2sd_phan_upper_50, aes(age,PC6), colour = 'green') +
  geom_line(data = PC6_0.2sd_phan_lower_50, aes(age,PC6), colour = 'green') +
  geom_line(data = PC6_0.0sd_phan_upper_50, aes(age,PC6), colour = 'orange') +
  geom_line(data = PC6_0.0sd_phan_lower_50, aes(age,PC6), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC6)) + 
  #geom_point(data = PC6_2sd_phan_upper, aes(age, PC6), colour = 'red') +
  scale_x_reverse(limits = c(540,0)) + theme_bw()

##### Precambrian PC6 #####

coords.all_pre <- subset(coords.all2, age > 541)

# 2 sigma

nsig <- 2

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_2sd_pre_upper <- as.data.frame(cbind(PC6_2sd_pre$x, PC6_2sd_pre$upper))
names(PC6_2sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_2sd_pre_upper$age, PC6_2sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[1], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_2sd_pre_upper, indices) {
  d <- PC6_2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_2sd_pre_upper_50)[1] <- "PC6"

PC6_2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_2sd_pre_upper_all)[1:ncol(PC6_2sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_2sd_pre_upper_all$group <- seq(1,(nrow(PC6_2sd_pre_upper_all)), 1)

PC6_2sd_pre_upper_all <- reshape2::melt(PC6_2sd_pre_upper_all, id.vars = "group")
names(PC6_2sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_2sd_pre_lower <- as.data.frame(cbind(PC6_2sd_pre$x, PC6_2sd_pre$lower))
names(PC6_2sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_2sd_pre_lower$age, PC6_2sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_2sd_pre_lower, indices) {
  d <- PC6_2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_2sd_pre_lower_50)[1] <- "PC6"

PC6_2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_2sd_pre_lower_all)[1:ncol(PC6_2sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_2sd_pre_lower_all$group <- seq(1,(nrow(PC6_2sd_pre_lower_all)), 1)

PC6_2sd_pre_lower_all <- reshape2::melt(PC6_2sd_pre_lower_all, id.vars = "group")
names(PC6_2sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.8sd_pre_upper <- as.data.frame(cbind(PC6_1.8sd_pre$x, PC6_1.8sd_pre$upper))
names(PC6_1.8sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.8sd_pre_upper$age, PC6_1.8sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_1.8sd_pre_upper, indices) {
  d <- PC6_1.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_1.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_1.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_1.8sd_pre_upper_50)[1] <- "PC6"

PC6_1.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_1.8sd_pre_upper_all)[1:ncol(PC6_1.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_1.8sd_pre_upper_all$group <- seq(1,(nrow(PC6_1.8sd_pre_upper_all)), 1)

PC6_1.8sd_pre_upper_all <- reshape2::melt(PC6_1.8sd_pre_upper_all, id.vars = "group")
names(PC6_1.8sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_1.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.8sd_pre_lower <- as.data.frame(cbind(PC6_1.8sd_pre$x, PC6_1.8sd_pre$lower))
names(PC6_1.8sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.8sd_pre_lower$age, PC6_1.8sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_1.8sd_pre_lower, indices) {
  d <- PC6_1.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_1.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_1.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_1.8sd_pre_lower_50)[1] <- "PC6"

PC6_1.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_1.8sd_pre_lower_all)[1:ncol(PC6_1.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_1.8sd_pre_lower_all$group <- seq(1,(nrow(PC6_1.8sd_pre_lower_all)), 1)

PC6_1.8sd_pre_lower_all <- reshape2::melt(PC6_1.8sd_pre_lower_all, id.vars = "group")
names(PC6_1.8sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_1.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.6sd_pre_upper <- as.data.frame(cbind(PC6_1.6sd_pre$x, PC6_1.6sd_pre$upper))
names(PC6_1.6sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.6sd_pre_upper$age, PC6_1.6sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_1.6sd_pre_upper, indices) {
  d <- PC6_1.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_1.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_1.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_1.6sd_pre_upper_50)[1] <- "PC6"

PC6_1.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_1.6sd_pre_upper_all)[1:ncol(PC6_1.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_1.6sd_pre_upper_all$group <- seq(1,(nrow(PC6_1.6sd_pre_upper_all)), 1)

PC6_1.6sd_pre_upper_all <- reshape2::melt(PC6_1.6sd_pre_upper_all, id.vars = "group")
names(PC6_1.6sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_1.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.6sd_pre_lower <- as.data.frame(cbind(PC6_1.6sd_pre$x, PC6_1.6sd_pre$lower))
names(PC6_1.6sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.6sd_pre_lower$age, PC6_1.6sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_1.6sd_pre_lower, indices) {
  d <- PC6_1.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_1.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_1.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_1.6sd_pre_lower_50)[1] <- "PC6"

PC6_1.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_1.6sd_pre_lower_all)[1:ncol(PC6_1.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_1.6sd_pre_lower_all$group <- seq(1,(nrow(PC6_1.6sd_pre_lower_all)), 1)

PC6_1.6sd_pre_lower_all <- reshape2::melt(PC6_1.6sd_pre_lower_all, id.vars = "group")
names(PC6_1.6sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_1.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.4sd_pre_upper <- as.data.frame(cbind(PC6_1.4sd_pre$x, PC6_1.4sd_pre$upper))
names(PC6_1.4sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.4sd_pre_upper$age, PC6_1.4sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_1.4sd_pre_upper, indices) {
  d <- PC6_1.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_1.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_1.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_1.4sd_pre_upper_50)[1] <- "PC6"

PC6_1.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_1.4sd_pre_upper_all)[1:ncol(PC6_1.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_1.4sd_pre_upper_all$group <- seq(1,(nrow(PC6_1.4sd_pre_upper_all)), 1)

PC6_1.4sd_pre_upper_all <- reshape2::melt(PC6_1.4sd_pre_upper_all, id.vars = "group")
names(PC6_1.4sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_1.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.4sd_pre_lower <- as.data.frame(cbind(PC6_1.4sd_pre$x, PC6_1.4sd_pre$lower))
names(PC6_1.4sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.4sd_pre_lower$age, PC6_1.4sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_1.4sd_pre_lower, indices) {
  d <- PC6_1.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_1.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_1.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_1.4sd_pre_lower_50)[1] <- "PC6"

PC6_1.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_1.4sd_pre_lower_all)[1:ncol(PC6_1.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_1.4sd_pre_lower_all$group <- seq(1,(nrow(PC6_1.4sd_pre_lower_all)), 1)

PC6_1.4sd_pre_lower_all <- reshape2::melt(PC6_1.4sd_pre_lower_all, id.vars = "group")
names(PC6_1.4sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_1.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.2sd_pre_upper <- as.data.frame(cbind(PC6_1.2sd_pre$x, PC6_1.2sd_pre$upper))
names(PC6_1.2sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.2sd_pre_upper$age, PC6_1.2sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_1.2sd_pre_upper, indices) {
  d <- PC6_1.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_1.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_1.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_1.2sd_pre_upper_50)[1] <- "PC6"

PC6_1.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_1.2sd_pre_upper_all)[1:ncol(PC6_1.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_1.2sd_pre_upper_all$group <- seq(1,(nrow(PC6_1.2sd_pre_upper_all)), 1)

PC6_1.2sd_pre_upper_all <- reshape2::melt(PC6_1.2sd_pre_upper_all, id.vars = "group")
names(PC6_1.2sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_1.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.2sd_pre_lower <- as.data.frame(cbind(PC6_1.2sd_pre$x, PC6_1.2sd_pre$lower))
names(PC6_1.2sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.2sd_pre_lower$age, PC6_1.2sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_1.2sd_pre_lower, indices) {
  d <- PC6_1.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_1.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_1.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_1.2sd_pre_lower_50)[1] <- "PC6"

PC6_1.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_1.2sd_pre_lower_all)[1:ncol(PC6_1.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_1.2sd_pre_lower_all$group <- seq(1,(nrow(PC6_1.2sd_pre_lower_all)), 1)

PC6_1.2sd_pre_lower_all <- reshape2::melt(PC6_1.2sd_pre_lower_all, id.vars = "group")
names(PC6_1.2sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_1.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.0sd_pre_upper <- as.data.frame(cbind(PC6_1.0sd_pre$x, PC6_1.0sd_pre$upper))
names(PC6_1.0sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.0sd_pre_upper$age, PC6_1.0sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_1.0sd_pre_upper, indices) {
  d <- PC6_1.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_1.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_1.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_1.0sd_pre_upper_50)[1] <- "PC6"

PC6_1.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_1.0sd_pre_upper_all)[1:ncol(PC6_1.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_1.0sd_pre_upper_all$group <- seq(1,(nrow(PC6_1.0sd_pre_upper_all)), 1)

PC6_1.0sd_pre_upper_all <- reshape2::melt(PC6_1.0sd_pre_upper_all, id.vars = "group")
names(PC6_1.0sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_1.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_1.0sd_pre_lower <- as.data.frame(cbind(PC6_1.0sd_pre$x, PC6_1.0sd_pre$lower))
names(PC6_1.0sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_1.0sd_pre_lower$age, PC6_1.0sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_1.0sd_pre_lower, indices) {
  d <- PC6_1.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_1.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_1.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_1.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_1.0sd_pre_lower_50)[1] <- "PC6"

PC6_1.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_1.0sd_pre_lower_all)[1:ncol(PC6_1.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_1.0sd_pre_lower_all$group <- seq(1,(nrow(PC6_1.0sd_pre_lower_all)), 1)

PC6_1.0sd_pre_lower_all <- reshape2::melt(PC6_1.0sd_pre_lower_all, id.vars = "group")
names(PC6_1.0sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_1.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.8sd_pre_upper <- as.data.frame(cbind(PC6_0.8sd_pre$x, PC6_0.8sd_pre$upper))
names(PC6_0.8sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.8sd_pre_upper$age, PC6_0.8sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_0.8sd_pre_upper, indices) {
  d <- PC6_0.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_0.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_0.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_0.8sd_pre_upper_50)[1] <- "PC6"

PC6_0.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_0.8sd_pre_upper_all)[1:ncol(PC6_0.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_0.8sd_pre_upper_all$group <- seq(1,(nrow(PC6_0.8sd_pre_upper_all)), 1)

PC6_0.8sd_pre_upper_all <- reshape2::melt(PC6_0.8sd_pre_upper_all, id.vars = "group")
names(PC6_0.8sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_0.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.8sd_pre_lower <- as.data.frame(cbind(PC6_0.8sd_pre$x, PC6_0.8sd_pre$lower))
names(PC6_0.8sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.8sd_pre_lower$age, PC6_0.8sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_0.8sd_pre_lower, indices) {
  d <- PC6_0.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_0.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_0.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_0.8sd_pre_lower_50)[1] <- "PC6"

PC6_0.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_0.8sd_pre_lower_all)[1:ncol(PC6_0.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_0.8sd_pre_lower_all$group <- seq(1,(nrow(PC6_0.8sd_pre_lower_all)), 1)

PC6_0.8sd_pre_lower_all <- reshape2::melt(PC6_0.8sd_pre_lower_all, id.vars = "group")
names(PC6_0.8sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_0.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.6sd_pre_upper <- as.data.frame(cbind(PC6_0.6sd_pre$x, PC6_0.6sd_pre$upper))
names(PC6_0.6sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.6sd_pre_upper$age, PC6_0.6sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_0.6sd_pre_upper, indices) {
  d <- PC6_0.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_0.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_0.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_0.6sd_pre_upper_50)[1] <- "PC6"

PC6_0.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_0.6sd_pre_upper_all)[1:ncol(PC6_0.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_0.6sd_pre_upper_all$group <- seq(1,(nrow(PC6_0.6sd_pre_upper_all)), 1)

PC6_0.6sd_pre_upper_all <- reshape2::melt(PC6_0.6sd_pre_upper_all, id.vars = "group")
names(PC6_0.6sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_0.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.6sd_pre_lower <- as.data.frame(cbind(PC6_0.6sd_pre$x, PC6_0.6sd_pre$lower))
names(PC6_0.6sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.6sd_pre_lower$age, PC6_0.6sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_0.6sd_pre_lower, indices) {
  d <- PC6_0.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_0.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_0.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_0.6sd_pre_lower_50)[1] <- "PC6"

PC6_0.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_0.6sd_pre_lower_all)[1:ncol(PC6_0.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_0.6sd_pre_lower_all$group <- seq(1,(nrow(PC6_0.6sd_pre_lower_all)), 1)

PC6_0.6sd_pre_lower_all <- reshape2::melt(PC6_0.6sd_pre_lower_all, id.vars = "group")
names(PC6_0.6sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_0.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.4sd_pre_upper <- as.data.frame(cbind(PC6_0.4sd_pre$x, PC6_0.4sd_pre$upper))
names(PC6_0.4sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.4sd_pre_upper$age, PC6_0.4sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_0.4sd_pre_upper, indices) {
  d <- PC6_0.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_0.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_0.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_0.4sd_pre_upper_50)[1] <- "PC6"

PC6_0.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_0.4sd_pre_upper_all)[1:ncol(PC6_0.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_0.4sd_pre_upper_all$group <- seq(1,(nrow(PC6_0.4sd_pre_upper_all)), 1)

PC6_0.4sd_pre_upper_all <- reshape2::melt(PC6_0.4sd_pre_upper_all, id.vars = "group")
names(PC6_0.4sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_0.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.4sd_pre_lower <- as.data.frame(cbind(PC6_0.4sd_pre$x, PC6_0.4sd_pre$lower))
names(PC6_0.4sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.4sd_pre_lower$age, PC6_0.4sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_0.4sd_pre_lower, indices) {
  d <- PC6_0.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_0.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_0.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_0.4sd_pre_lower_50)[1] <- "PC6"

PC6_0.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_0.4sd_pre_lower_all)[1:ncol(PC6_0.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_0.4sd_pre_lower_all$group <- seq(1,(nrow(PC6_0.4sd_pre_lower_all)), 1)

PC6_0.4sd_pre_lower_all <- reshape2::melt(PC6_0.4sd_pre_lower_all, id.vars = "group")
names(PC6_0.4sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_0.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.2sd_pre_upper <- as.data.frame(cbind(PC6_0.2sd_pre$x, PC6_0.2sd_pre$upper))
names(PC6_0.2sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.2sd_pre_upper$age, PC6_0.2sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_0.2sd_pre_upper, indices) {
  d <- PC6_0.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_0.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_0.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_0.2sd_pre_upper_50)[1] <- "PC6"

PC6_0.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_0.2sd_pre_upper_all)[1:ncol(PC6_0.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_0.2sd_pre_upper_all$group <- seq(1,(nrow(PC6_0.2sd_pre_upper_all)), 1)

PC6_0.2sd_pre_upper_all <- reshape2::melt(PC6_0.2sd_pre_upper_all, id.vars = "group")
names(PC6_0.2sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_0.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.2sd_pre_lower <- as.data.frame(cbind(PC6_0.2sd_pre$x, PC6_0.2sd_pre$lower))
names(PC6_0.2sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.2sd_pre_lower$age, PC6_0.2sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_0.2sd_pre_lower, indices) {
  d <- PC6_0.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_0.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_0.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_0.2sd_pre_lower_50)[1] <- "PC6"

PC6_0.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_0.2sd_pre_lower_all)[1:ncol(PC6_0.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_0.2sd_pre_lower_all$group <- seq(1,(nrow(PC6_0.2sd_pre_lower_all)), 1)

PC6_0.2sd_pre_lower_all <- reshape2::melt(PC6_0.2sd_pre_lower_all, id.vars = "group")
names(PC6_0.2sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_0.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC6 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.0sd_pre_upper <- as.data.frame(cbind(PC6_0.0sd_pre$x, PC6_0.0sd_pre$upper))
names(PC6_0.0sd_pre_upper)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.0sd_pre_upper$age, PC6_0.0sd_pre_upper$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC6_0.0sd_pre_upper, indices) {
  d <- PC6_0.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC6_0.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC6_0.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC6_0.0sd_pre_upper_50)[1] <- "PC6"

PC6_0.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC6_0.0sd_pre_upper_all)[1:ncol(PC6_0.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC6_0.0sd_pre_upper_all$group <- seq(1,(nrow(PC6_0.0sd_pre_upper_all)), 1)

PC6_0.0sd_pre_upper_all <- reshape2::melt(PC6_0.0sd_pre_upper_all, id.vars = "group")
names(PC6_0.0sd_pre_upper_all)[2:3] <- c("age", "PC6")

PC6_0.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC6 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC6_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC6, degree=1, nsigma = nsig, span = span)
PC6_0.0sd_pre_lower <- as.data.frame(cbind(PC6_0.0sd_pre$x, PC6_0.0sd_pre$lower))
names(PC6_0.0sd_pre_lower)[1:2] <- c("age","PC6")

# CV spans
loess.predict <- loess.as(PC6_0.0sd_pre_lower$age, PC6_0.0sd_pre_lower$PC6, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC6_0.0sd_pre_lower, indices) {
  d <- PC6_0.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC6 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC6_0.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC6_0.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC6_0.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC6_0.0sd_pre_lower_50)[1] <- "PC6"

PC6_0.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC6_0.0sd_pre_lower_all)[1:ncol(PC6_0.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC6_0.0sd_pre_lower_all$group <- seq(1,(nrow(PC6_0.0sd_pre_lower_all)), 1)

PC6_0.0sd_pre_lower_all <- reshape2::melt(PC6_0.0sd_pre_lower_all, id.vars = "group")
names(PC6_0.0sd_pre_lower_all)[2:3] <- c("age", "PC6")

PC6_0.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

b <- ggplot() + 
  #geom_line(data = PC6_2sd_pre_upper_all, aes(age2,PC6, group = group)) +
  #geom_line(data = PC6_2sd_pre_lower_all, aes(age2,PC6, group = group)) +
  geom_line(data = PC6_2sd_pre_upper_50, aes(age,PC6), colour = 'blue') +
  geom_line(data = PC6_2sd_pre_lower_50, aes(age,PC6), colour = 'blue') +
  geom_line(data = PC6_1.6sd_pre_upper_50, aes(age,PC6), colour = 'green') +
  geom_line(data = PC6_1.6sd_pre_lower_50, aes(age,PC6), colour = 'green') +
  geom_line(data = PC6_1.4sd_pre_upper_50, aes(age,PC6), colour = 'orange') +
  geom_line(data = PC6_1.4sd_pre_lower_50, aes(age,PC6), colour = 'orange') +
  geom_line(data = PC6_1.2sd_pre_upper_50, aes(age,PC6), colour = 'pink') +
  geom_line(data = PC6_1.2sd_pre_lower_50, aes(age,PC6), colour = 'pink') +
  geom_line(data = PC6_1.0sd_pre_upper_50, aes(age,PC6), colour = 'brown') +
  geom_line(data = PC6_1.0sd_pre_lower_50, aes(age,PC6), colour = 'brown') +
  geom_line(data = PC6_0.8sd_pre_upper_50, aes(age,PC6), colour = 'yellow') +
  geom_line(data = PC6_0.8sd_pre_lower_50, aes(age,PC6), colour = 'yellow') +
  geom_line(data = PC6_0.6sd_pre_upper_50, aes(age,PC6), colour = 'red') +
  geom_line(data = PC6_0.6sd_pre_lower_50, aes(age,PC6), colour = 'red') +
  geom_line(data = PC6_0.4sd_pre_upper_50, aes(age,PC6), colour = 'blue') +
  geom_line(data = PC6_0.4sd_pre_lower_50, aes(age,PC6), colour = 'blue') +
  geom_line(data = PC6_0.2sd_pre_upper_50, aes(age,PC6), colour = 'green') +
  geom_line(data = PC6_0.2sd_pre_lower_50, aes(age,PC6), colour = 'green') +
  geom_line(data = PC6_0.0sd_pre_upper_50, aes(age,PC6), colour = 'orange') +
  geom_line(data = PC6_0.0sd_pre_lower_50, aes(age,PC6), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC6)) + 
  #geom_point(data = PC6_2sd_pre_upper, aes(age, PC6), colour = 'red') +
  scale_x_reverse(limits = c(4000,540)) + theme_bw()

grid.arrange(b,a,ncol=2)


##### Phanerozoic PC7 #####

coords.all_Phan <- subset(coords.all2, age < 541)

# 2 sigma

nsig <- 2

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_2sd_phan_upper <- as.data.frame(cbind(PC7_2sd_phan$x, PC7_2sd_phan$upper))
names(PC7_2sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_2sd_phan_upper$age, PC7_2sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_2sd_phan_upper, indices) {
  d <- PC7_2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_2sd_phan_upper_50)[1] <- "PC7"

PC7_2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_2sd_phan_upper_all)[1:ncol(PC7_2sd_phan_upper_all)] <- seq(0,540, 1)
PC7_2sd_phan_upper_all$group <- seq(1,(nrow(PC7_2sd_phan_upper_all)), 1)

PC7_2sd_phan_upper_all <- reshape2::melt(PC7_2sd_phan_upper_all, id.vars = "group")
names(PC7_2sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_2sd_phan_lower <- as.data.frame(cbind(PC7_2sd_phan$x, PC7_2sd_phan$lower))
names(PC7_2sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_2sd_phan_lower$age, PC7_2sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_2sd_phan_lower, indices) {
  d <- PC7_2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_2sd_phan_lower_50)[1] <- "PC7"

PC7_2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_2sd_phan_lower_all)[1:ncol(PC7_2sd_phan_lower_all)] <- seq(0,540, 1)
PC7_2sd_phan_lower_all$group <- seq(1,(nrow(PC7_2sd_phan_lower_all)), 1)

PC7_2sd_phan_lower_all <- reshape2::melt(PC7_2sd_phan_lower_all, id.vars = "group")
names(PC7_2sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.8sd_phan_upper <- as.data.frame(cbind(PC7_1.8sd_phan$x, PC7_1.8sd_phan$upper))
names(PC7_1.8sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.8sd_phan_upper$age, PC7_1.8sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_1.8sd_phan_upper, indices) {
  d <- PC7_1.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_1.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_1.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_1.8sd_phan_upper_50)[1] <- "PC7"

PC7_1.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_1.8sd_phan_upper_all)[1:ncol(PC7_1.8sd_phan_upper_all)] <- seq(0,540, 1)
PC7_1.8sd_phan_upper_all$group <- seq(1,(nrow(PC7_1.8sd_phan_upper_all)), 1)

PC7_1.8sd_phan_upper_all <- reshape2::melt(PC7_1.8sd_phan_upper_all, id.vars = "group")
names(PC7_1.8sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_1.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.8sd_phan_lower <- as.data.frame(cbind(PC7_1.8sd_phan$x, PC7_1.8sd_phan$lower))
names(PC7_1.8sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.8sd_phan_lower$age, PC7_1.8sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_1.8sd_phan_lower, indices) {
  d <- PC7_1.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_1.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_1.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_1.8sd_phan_lower_50)[1] <- "PC7"

PC7_1.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_1.8sd_phan_lower_all)[1:ncol(PC7_1.8sd_phan_lower_all)] <- seq(0,540, 1)
PC7_1.8sd_phan_lower_all$group <- seq(1,(nrow(PC7_1.8sd_phan_lower_all)), 1)

PC7_1.8sd_phan_lower_all <- reshape2::melt(PC7_1.8sd_phan_lower_all, id.vars = "group")
names(PC7_1.8sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_1.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.6sd_phan_upper <- as.data.frame(cbind(PC7_1.6sd_phan$x, PC7_1.6sd_phan$upper))
names(PC7_1.6sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.6sd_phan_upper$age, PC7_1.6sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_1.6sd_phan_upper, indices) {
  d <- PC7_1.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_1.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_1.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_1.6sd_phan_upper_50)[1] <- "PC7"

PC7_1.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_1.6sd_phan_upper_all)[1:ncol(PC7_1.6sd_phan_upper_all)] <- seq(0,540, 1)
PC7_1.6sd_phan_upper_all$group <- seq(1,(nrow(PC7_1.6sd_phan_upper_all)), 1)

PC7_1.6sd_phan_upper_all <- reshape2::melt(PC7_1.6sd_phan_upper_all, id.vars = "group")
names(PC7_1.6sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_1.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.6sd_phan_lower <- as.data.frame(cbind(PC7_1.6sd_phan$x, PC7_1.6sd_phan$lower))
names(PC7_1.6sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.6sd_phan_lower$age, PC7_1.6sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_1.6sd_phan_lower, indices) {
  d <- PC7_1.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_1.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_1.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_1.6sd_phan_lower_50)[1] <- "PC7"

PC7_1.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_1.6sd_phan_lower_all)[1:ncol(PC7_1.6sd_phan_lower_all)] <- seq(0,540, 1)
PC7_1.6sd_phan_lower_all$group <- seq(1,(nrow(PC7_1.6sd_phan_lower_all)), 1)

PC7_1.6sd_phan_lower_all <- reshape2::melt(PC7_1.6sd_phan_lower_all, id.vars = "group")
names(PC7_1.6sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_1.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.4sd_phan_upper <- as.data.frame(cbind(PC7_1.4sd_phan$x, PC7_1.4sd_phan$upper))
names(PC7_1.4sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.4sd_phan_upper$age, PC7_1.4sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_1.4sd_phan_upper, indices) {
  d <- PC7_1.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_1.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_1.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_1.4sd_phan_upper_50)[1] <- "PC7"

PC7_1.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_1.4sd_phan_upper_all)[1:ncol(PC7_1.4sd_phan_upper_all)] <- seq(0,540, 1)
PC7_1.4sd_phan_upper_all$group <- seq(1,(nrow(PC7_1.4sd_phan_upper_all)), 1)

PC7_1.4sd_phan_upper_all <- reshape2::melt(PC7_1.4sd_phan_upper_all, id.vars = "group")
names(PC7_1.4sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_1.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.4sd_phan_lower <- as.data.frame(cbind(PC7_1.4sd_phan$x, PC7_1.4sd_phan$lower))
names(PC7_1.4sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.4sd_phan_lower$age, PC7_1.4sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_1.4sd_phan_lower, indices) {
  d <- PC7_1.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_1.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_1.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_1.4sd_phan_lower_50)[1] <- "PC7"

PC7_1.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_1.4sd_phan_lower_all)[1:ncol(PC7_1.4sd_phan_lower_all)] <- seq(0,540, 1)
PC7_1.4sd_phan_lower_all$group <- seq(1,(nrow(PC7_1.4sd_phan_lower_all)), 1)

PC7_1.4sd_phan_lower_all <- reshape2::melt(PC7_1.4sd_phan_lower_all, id.vars = "group")
names(PC7_1.4sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_1.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.2sd_phan_upper <- as.data.frame(cbind(PC7_1.2sd_phan$x, PC7_1.2sd_phan$upper))
names(PC7_1.2sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.2sd_phan_upper$age, PC7_1.2sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_1.2sd_phan_upper, indices) {
  d <- PC7_1.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_1.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_1.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_1.2sd_phan_upper_50)[1] <- "PC7"

PC7_1.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_1.2sd_phan_upper_all)[1:ncol(PC7_1.2sd_phan_upper_all)] <- seq(0,540, 1)
PC7_1.2sd_phan_upper_all$group <- seq(1,(nrow(PC7_1.2sd_phan_upper_all)), 1)

PC7_1.2sd_phan_upper_all <- reshape2::melt(PC7_1.2sd_phan_upper_all, id.vars = "group")
names(PC7_1.2sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_1.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.2sd_phan_lower <- as.data.frame(cbind(PC7_1.2sd_phan$x, PC7_1.2sd_phan$lower))
names(PC7_1.2sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.2sd_phan_lower$age, PC7_1.2sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_1.2sd_phan_lower, indices) {
  d <- PC7_1.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_1.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_1.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_1.2sd_phan_lower_50)[1] <- "PC7"

PC7_1.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_1.2sd_phan_lower_all)[1:ncol(PC7_1.2sd_phan_lower_all)] <- seq(0,540, 1)
PC7_1.2sd_phan_lower_all$group <- seq(1,(nrow(PC7_1.2sd_phan_lower_all)), 1)

PC7_1.2sd_phan_lower_all <- reshape2::melt(PC7_1.2sd_phan_lower_all, id.vars = "group")
names(PC7_1.2sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_1.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.0sd_phan_upper <- as.data.frame(cbind(PC7_1.0sd_phan$x, PC7_1.0sd_phan$upper))
names(PC7_1.0sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.0sd_phan_upper$age, PC7_1.0sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_1.0sd_phan_upper, indices) {
  d <- PC7_1.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_1.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_1.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_1.0sd_phan_upper_50)[1] <- "PC7"

PC7_1.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_1.0sd_phan_upper_all)[1:ncol(PC7_1.0sd_phan_upper_all)] <- seq(0,540, 1)
PC7_1.0sd_phan_upper_all$group <- seq(1,(nrow(PC7_1.0sd_phan_upper_all)), 1)

PC7_1.0sd_phan_upper_all <- reshape2::melt(PC7_1.0sd_phan_upper_all, id.vars = "group")
names(PC7_1.0sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_1.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.0sd_phan_lower <- as.data.frame(cbind(PC7_1.0sd_phan$x, PC7_1.0sd_phan$lower))
names(PC7_1.0sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.0sd_phan_lower$age, PC7_1.0sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_1.0sd_phan_lower, indices) {
  d <- PC7_1.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_1.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_1.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_1.0sd_phan_lower_50)[1] <- "PC7"

PC7_1.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_1.0sd_phan_lower_all)[1:ncol(PC7_1.0sd_phan_lower_all)] <- seq(0,540, 1)
PC7_1.0sd_phan_lower_all$group <- seq(1,(nrow(PC7_1.0sd_phan_lower_all)), 1)

PC7_1.0sd_phan_lower_all <- reshape2::melt(PC7_1.0sd_phan_lower_all, id.vars = "group")
names(PC7_1.0sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_1.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.8sd_phan_upper <- as.data.frame(cbind(PC7_0.8sd_phan$x, PC7_0.8sd_phan$upper))
names(PC7_0.8sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.8sd_phan_upper$age, PC7_0.8sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_0.8sd_phan_upper, indices) {
  d <- PC7_0.8sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_0.8sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_0.8sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.8sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_0.8sd_phan_upper_50)[1] <- "PC7"

PC7_0.8sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_0.8sd_phan_upper_all)[1:ncol(PC7_0.8sd_phan_upper_all)] <- seq(0,540, 1)
PC7_0.8sd_phan_upper_all$group <- seq(1,(nrow(PC7_0.8sd_phan_upper_all)), 1)

PC7_0.8sd_phan_upper_all <- reshape2::melt(PC7_0.8sd_phan_upper_all, id.vars = "group")
names(PC7_0.8sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_0.8sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.8sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.8sd_phan_lower <- as.data.frame(cbind(PC7_0.8sd_phan$x, PC7_0.8sd_phan$lower))
names(PC7_0.8sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.8sd_phan_lower$age, PC7_0.8sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_0.8sd_phan_lower, indices) {
  d <- PC7_0.8sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_0.8sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_0.8sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.8sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_0.8sd_phan_lower_50)[1] <- "PC7"

PC7_0.8sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_0.8sd_phan_lower_all)[1:ncol(PC7_0.8sd_phan_lower_all)] <- seq(0,540, 1)
PC7_0.8sd_phan_lower_all$group <- seq(1,(nrow(PC7_0.8sd_phan_lower_all)), 1)

PC7_0.8sd_phan_lower_all <- reshape2::melt(PC7_0.8sd_phan_lower_all, id.vars = "group")
names(PC7_0.8sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_0.8sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.6sd_phan_upper <- as.data.frame(cbind(PC7_0.6sd_phan$x, PC7_0.6sd_phan$upper))
names(PC7_0.6sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.6sd_phan_upper$age, PC7_0.6sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_0.6sd_phan_upper, indices) {
  d <- PC7_0.6sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_0.6sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_0.6sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.6sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_0.6sd_phan_upper_50)[1] <- "PC7"

PC7_0.6sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_0.6sd_phan_upper_all)[1:ncol(PC7_0.6sd_phan_upper_all)] <- seq(0,540, 1)
PC7_0.6sd_phan_upper_all$group <- seq(1,(nrow(PC7_0.6sd_phan_upper_all)), 1)

PC7_0.6sd_phan_upper_all <- reshape2::melt(PC7_0.6sd_phan_upper_all, id.vars = "group")
names(PC7_0.6sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_0.6sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.6sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.6sd_phan_lower <- as.data.frame(cbind(PC7_0.6sd_phan$x, PC7_0.6sd_phan$lower))
names(PC7_0.6sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.6sd_phan_lower$age, PC7_0.6sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_0.6sd_phan_lower, indices) {
  d <- PC7_0.6sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_0.6sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_0.6sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.6sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_0.6sd_phan_lower_50)[1] <- "PC7"

PC7_0.6sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_0.6sd_phan_lower_all)[1:ncol(PC7_0.6sd_phan_lower_all)] <- seq(0,540, 1)
PC7_0.6sd_phan_lower_all$group <- seq(1,(nrow(PC7_0.6sd_phan_lower_all)), 1)

PC7_0.6sd_phan_lower_all <- reshape2::melt(PC7_0.6sd_phan_lower_all, id.vars = "group")
names(PC7_0.6sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_0.6sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.4sd_phan_upper <- as.data.frame(cbind(PC7_0.4sd_phan$x, PC7_0.4sd_phan$upper))
names(PC7_0.4sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.4sd_phan_upper$age, PC7_0.4sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_0.4sd_phan_upper, indices) {
  d <- PC7_0.4sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_0.4sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_0.4sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.4sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_0.4sd_phan_upper_50)[1] <- "PC7"

PC7_0.4sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_0.4sd_phan_upper_all)[1:ncol(PC7_0.4sd_phan_upper_all)] <- seq(0,540, 1)
PC7_0.4sd_phan_upper_all$group <- seq(1,(nrow(PC7_0.4sd_phan_upper_all)), 1)

PC7_0.4sd_phan_upper_all <- reshape2::melt(PC7_0.4sd_phan_upper_all, id.vars = "group")
names(PC7_0.4sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_0.4sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.4sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.4sd_phan_lower <- as.data.frame(cbind(PC7_0.4sd_phan$x, PC7_0.4sd_phan$lower))
names(PC7_0.4sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.4sd_phan_lower$age, PC7_0.4sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_0.4sd_phan_lower, indices) {
  d <- PC7_0.4sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_0.4sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_0.4sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.4sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_0.4sd_phan_lower_50)[1] <- "PC7"

PC7_0.4sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_0.4sd_phan_lower_all)[1:ncol(PC7_0.4sd_phan_lower_all)] <- seq(0,540, 1)
PC7_0.4sd_phan_lower_all$group <- seq(1,(nrow(PC7_0.4sd_phan_lower_all)), 1)

PC7_0.4sd_phan_lower_all <- reshape2::melt(PC7_0.4sd_phan_lower_all, id.vars = "group")
names(PC7_0.4sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_0.4sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.2sd_phan_upper <- as.data.frame(cbind(PC7_0.2sd_phan$x, PC7_0.2sd_phan$upper))
names(PC7_0.2sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.2sd_phan_upper$age, PC7_0.2sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_0.2sd_phan_upper, indices) {
  d <- PC7_0.2sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_0.2sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_0.2sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.2sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_0.2sd_phan_upper_50)[1] <- "PC7"

PC7_0.2sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_0.2sd_phan_upper_all)[1:ncol(PC7_0.2sd_phan_upper_all)] <- seq(0,540, 1)
PC7_0.2sd_phan_upper_all$group <- seq(1,(nrow(PC7_0.2sd_phan_upper_all)), 1)

PC7_0.2sd_phan_upper_all <- reshape2::melt(PC7_0.2sd_phan_upper_all, id.vars = "group")
names(PC7_0.2sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_0.2sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.2sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.2sd_phan_lower <- as.data.frame(cbind(PC7_0.2sd_phan$x, PC7_0.2sd_phan$lower))
names(PC7_0.2sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.2sd_phan_lower$age, PC7_0.2sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_0.2sd_phan_lower, indices) {
  d <- PC7_0.2sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_0.2sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_0.2sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.2sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_0.2sd_phan_lower_50)[1] <- "PC7"

PC7_0.2sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_0.2sd_phan_lower_all)[1:ncol(PC7_0.2sd_phan_lower_all)] <- seq(0,540, 1)
PC7_0.2sd_phan_lower_all$group <- seq(1,(nrow(PC7_0.2sd_phan_lower_all)), 1)

PC7_0.2sd_phan_lower_all <- reshape2::melt(PC7_0.2sd_phan_lower_all, id.vars = "group")
names(PC7_0.2sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_0.2sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.0sd_phan_upper <- as.data.frame(cbind(PC7_0.0sd_phan$x, PC7_0.0sd_phan$upper))
names(PC7_0.0sd_phan_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.0sd_phan_upper$age, PC7_0.0sd_phan_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_0.0sd_phan_upper, indices) {
  d <- PC7_0.0sd_phan_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_0.0sd_phan_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_0.0sd_phan_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.0sd_phan_upper_50$age <- seq(0,540, 1)
names(PC7_0.0sd_phan_upper_50)[1] <- "PC7"

PC7_0.0sd_phan_upper_all <- as.data.frame(loess_boot$t)
names(PC7_0.0sd_phan_upper_all)[1:ncol(PC7_0.0sd_phan_upper_all)] <- seq(0,540, 1)
PC7_0.0sd_phan_upper_all$group <- seq(1,(nrow(PC7_0.0sd_phan_upper_all)), 1)

PC7_0.0sd_phan_upper_all <- reshape2::melt(PC7_0.0sd_phan_upper_all, id.vars = "group")
names(PC7_0.0sd_phan_upper_all)[2:3] <- c("age", "PC7")

PC7_0.0sd_phan_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_Phan$age, coords.all_Phan$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.0sd_phan <- msir::loess.sd(coords.all_Phan$age, coords.all_Phan$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.0sd_phan_lower <- as.data.frame(cbind(PC7_0.0sd_phan$x, PC7_0.0sd_phan$lower))
names(PC7_0.0sd_phan_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.0sd_phan_lower$age, PC7_0.0sd_phan_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_0.0sd_phan_lower, indices) {
  d <- PC7_0.0sd_phan_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_0.0sd_phan_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_0.0sd_phan_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.0sd_phan_lower_50$age <- seq(0,540, 1)
names(PC7_0.0sd_phan_lower_50)[1] <- "PC7"

PC7_0.0sd_phan_lower_all <- as.data.frame(loess_boot$t)
names(PC7_0.0sd_phan_lower_all)[1:ncol(PC7_0.0sd_phan_lower_all)] <- seq(0,540, 1)
PC7_0.0sd_phan_lower_all$group <- seq(1,(nrow(PC7_0.0sd_phan_lower_all)), 1)

PC7_0.0sd_phan_lower_all <- reshape2::melt(PC7_0.0sd_phan_lower_all, id.vars = "group")
names(PC7_0.0sd_phan_lower_all)[2:3] <- c("age", "PC7")

PC7_0.0sd_phan_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

a <- ggplot() + 
  #geom_line(data = PC7_2sd_phan_upper_all, aes(age2,PC7, group = group)) +
  #geom_line(data = PC7_2sd_phan_lower_all, aes(age2,PC7, group = group)) +
  geom_line(data = PC7_2sd_phan_upper_50, aes(age,PC7), colour = 'blue') +
  geom_line(data = PC7_2sd_phan_lower_50, aes(age,PC7), colour = 'blue') +
  geom_line(data = PC7_1.6sd_phan_upper_50, aes(age,PC7), colour = 'green') +
  geom_line(data = PC7_1.6sd_phan_lower_50, aes(age,PC7), colour = 'green') +
  geom_line(data = PC7_1.4sd_phan_upper_50, aes(age,PC7), colour = 'orange') +
  geom_line(data = PC7_1.4sd_phan_lower_50, aes(age,PC7), colour = 'orange') +
  geom_line(data = PC7_1.2sd_phan_upper_50, aes(age,PC7), colour = 'pink') +
  geom_line(data = PC7_1.2sd_phan_lower_50, aes(age,PC7), colour = 'pink') +
  geom_line(data = PC7_1.0sd_phan_upper_50, aes(age,PC7), colour = 'brown') +
  geom_line(data = PC7_1.0sd_phan_lower_50, aes(age,PC7), colour = 'brown') +
  geom_line(data = PC7_0.8sd_phan_upper_50, aes(age,PC7), colour = 'yellow') +
  geom_line(data = PC7_0.8sd_phan_lower_50, aes(age,PC7), colour = 'yellow') +
  geom_line(data = PC7_0.6sd_phan_upper_50, aes(age,PC7), colour = 'red') +
  geom_line(data = PC7_0.6sd_phan_lower_50, aes(age,PC7), colour = 'red') +
  geom_line(data = PC7_0.4sd_phan_upper_50, aes(age,PC7), colour = 'blue') +
  geom_line(data = PC7_0.4sd_phan_lower_50, aes(age,PC7), colour = 'blue') +
  geom_line(data = PC7_0.2sd_phan_upper_50, aes(age,PC7), colour = 'green') +
  geom_line(data = PC7_0.2sd_phan_lower_50, aes(age,PC7), colour = 'green') +
  geom_line(data = PC7_0.0sd_phan_upper_50, aes(age,PC7), colour = 'orange') +
  geom_line(data = PC7_0.0sd_phan_lower_50, aes(age,PC7), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC7)) + 
  #geom_point(data = PC7_2sd_phan_upper, aes(age, PC7), colour = 'red') +
  scale_x_reverse(limits = c(540,0)) + theme_bw()

##### Precambrian PC7 #####

coords.all_pre <- subset(coords.all2, age > 541)

# 2 sigma

nsig <- 2

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_2sd_pre_upper <- as.data.frame(cbind(PC7_2sd_pre$x, PC7_2sd_pre$upper))
names(PC7_2sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_2sd_pre_upper$age, PC7_2sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[1], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_2sd_pre_upper, indices) {
  d <- PC7_2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_2sd_pre_upper_50)[1] <- "PC7"

PC7_2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_2sd_pre_upper_all)[1:ncol(PC7_2sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_2sd_pre_upper_all$group <- seq(1,(nrow(PC7_2sd_pre_upper_all)), 1)

PC7_2sd_pre_upper_all <- reshape2::melt(PC7_2sd_pre_upper_all, id.vars = "group")
names(PC7_2sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_2sd_pre_lower <- as.data.frame(cbind(PC7_2sd_pre$x, PC7_2sd_pre$lower))
names(PC7_2sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_2sd_pre_lower$age, PC7_2sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_2sd_pre_lower, indices) {
  d <- PC7_2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_2sd_pre_lower_50)[1] <- "PC7"

PC7_2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_2sd_pre_lower_all)[1:ncol(PC7_2sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_2sd_pre_lower_all$group <- seq(1,(nrow(PC7_2sd_pre_lower_all)), 1)

PC7_2sd_pre_lower_all <- reshape2::melt(PC7_2sd_pre_lower_all, id.vars = "group")
names(PC7_2sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.8 sigma

nsig <- 1.8

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.8sd_pre_upper <- as.data.frame(cbind(PC7_1.8sd_pre$x, PC7_1.8sd_pre$upper))
names(PC7_1.8sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.8sd_pre_upper$age, PC7_1.8sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_1.8sd_pre_upper, indices) {
  d <- PC7_1.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_1.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_1.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_1.8sd_pre_upper_50)[1] <- "PC7"

PC7_1.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_1.8sd_pre_upper_all)[1:ncol(PC7_1.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_1.8sd_pre_upper_all$group <- seq(1,(nrow(PC7_1.8sd_pre_upper_all)), 1)

PC7_1.8sd_pre_upper_all <- reshape2::melt(PC7_1.8sd_pre_upper_all, id.vars = "group")
names(PC7_1.8sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_1.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.8sd_pre_lower <- as.data.frame(cbind(PC7_1.8sd_pre$x, PC7_1.8sd_pre$lower))
names(PC7_1.8sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.8sd_pre_lower$age, PC7_1.8sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_1.8sd_pre_lower, indices) {
  d <- PC7_1.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_1.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_1.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_1.8sd_pre_lower_50)[1] <- "PC7"

PC7_1.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_1.8sd_pre_lower_all)[1:ncol(PC7_1.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_1.8sd_pre_lower_all$group <- seq(1,(nrow(PC7_1.8sd_pre_lower_all)), 1)

PC7_1.8sd_pre_lower_all <- reshape2::melt(PC7_1.8sd_pre_lower_all, id.vars = "group")
names(PC7_1.8sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_1.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.6 sigma

nsig <- 1.6

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.6sd_pre_upper <- as.data.frame(cbind(PC7_1.6sd_pre$x, PC7_1.6sd_pre$upper))
names(PC7_1.6sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.6sd_pre_upper$age, PC7_1.6sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_1.6sd_pre_upper, indices) {
  d <- PC7_1.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_1.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_1.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_1.6sd_pre_upper_50)[1] <- "PC7"

PC7_1.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_1.6sd_pre_upper_all)[1:ncol(PC7_1.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_1.6sd_pre_upper_all$group <- seq(1,(nrow(PC7_1.6sd_pre_upper_all)), 1)

PC7_1.6sd_pre_upper_all <- reshape2::melt(PC7_1.6sd_pre_upper_all, id.vars = "group")
names(PC7_1.6sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_1.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.6sd_pre_lower <- as.data.frame(cbind(PC7_1.6sd_pre$x, PC7_1.6sd_pre$lower))
names(PC7_1.6sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.6sd_pre_lower$age, PC7_1.6sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_1.6sd_pre_lower, indices) {
  d <- PC7_1.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_1.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_1.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_1.6sd_pre_lower_50)[1] <- "PC7"

PC7_1.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_1.6sd_pre_lower_all)[1:ncol(PC7_1.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_1.6sd_pre_lower_all$group <- seq(1,(nrow(PC7_1.6sd_pre_lower_all)), 1)

PC7_1.6sd_pre_lower_all <- reshape2::melt(PC7_1.6sd_pre_lower_all, id.vars = "group")
names(PC7_1.6sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_1.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.4 sigma

nsig <- 1.4

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.4sd_pre_upper <- as.data.frame(cbind(PC7_1.4sd_pre$x, PC7_1.4sd_pre$upper))
names(PC7_1.4sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.4sd_pre_upper$age, PC7_1.4sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_1.4sd_pre_upper, indices) {
  d <- PC7_1.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_1.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_1.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_1.4sd_pre_upper_50)[1] <- "PC7"

PC7_1.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_1.4sd_pre_upper_all)[1:ncol(PC7_1.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_1.4sd_pre_upper_all$group <- seq(1,(nrow(PC7_1.4sd_pre_upper_all)), 1)

PC7_1.4sd_pre_upper_all <- reshape2::melt(PC7_1.4sd_pre_upper_all, id.vars = "group")
names(PC7_1.4sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_1.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.4sd_pre_lower <- as.data.frame(cbind(PC7_1.4sd_pre$x, PC7_1.4sd_pre$lower))
names(PC7_1.4sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.4sd_pre_lower$age, PC7_1.4sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_1.4sd_pre_lower, indices) {
  d <- PC7_1.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_1.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_1.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_1.4sd_pre_lower_50)[1] <- "PC7"

PC7_1.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_1.4sd_pre_lower_all)[1:ncol(PC7_1.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_1.4sd_pre_lower_all$group <- seq(1,(nrow(PC7_1.4sd_pre_lower_all)), 1)

PC7_1.4sd_pre_lower_all <- reshape2::melt(PC7_1.4sd_pre_lower_all, id.vars = "group")
names(PC7_1.4sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_1.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.2 sigma

nsig <- 1.2

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.2sd_pre_upper <- as.data.frame(cbind(PC7_1.2sd_pre$x, PC7_1.2sd_pre$upper))
names(PC7_1.2sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.2sd_pre_upper$age, PC7_1.2sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_1.2sd_pre_upper, indices) {
  d <- PC7_1.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_1.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_1.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_1.2sd_pre_upper_50)[1] <- "PC7"

PC7_1.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_1.2sd_pre_upper_all)[1:ncol(PC7_1.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_1.2sd_pre_upper_all$group <- seq(1,(nrow(PC7_1.2sd_pre_upper_all)), 1)

PC7_1.2sd_pre_upper_all <- reshape2::melt(PC7_1.2sd_pre_upper_all, id.vars = "group")
names(PC7_1.2sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_1.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.2sd_pre_lower <- as.data.frame(cbind(PC7_1.2sd_pre$x, PC7_1.2sd_pre$lower))
names(PC7_1.2sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.2sd_pre_lower$age, PC7_1.2sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_1.2sd_pre_lower, indices) {
  d <- PC7_1.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_1.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_1.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_1.2sd_pre_lower_50)[1] <- "PC7"

PC7_1.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_1.2sd_pre_lower_all)[1:ncol(PC7_1.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_1.2sd_pre_lower_all$group <- seq(1,(nrow(PC7_1.2sd_pre_lower_all)), 1)

PC7_1.2sd_pre_lower_all <- reshape2::melt(PC7_1.2sd_pre_lower_all, id.vars = "group")
names(PC7_1.2sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_1.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.0 sigma

nsig <- 1.0

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.0sd_pre_upper <- as.data.frame(cbind(PC7_1.0sd_pre$x, PC7_1.0sd_pre$upper))
names(PC7_1.0sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.0sd_pre_upper$age, PC7_1.0sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_1.0sd_pre_upper, indices) {
  d <- PC7_1.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_1.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_1.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_1.0sd_pre_upper_50)[1] <- "PC7"

PC7_1.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_1.0sd_pre_upper_all)[1:ncol(PC7_1.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_1.0sd_pre_upper_all$group <- seq(1,(nrow(PC7_1.0sd_pre_upper_all)), 1)

PC7_1.0sd_pre_upper_all <- reshape2::melt(PC7_1.0sd_pre_upper_all, id.vars = "group")
names(PC7_1.0sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_1.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_1.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_1.0sd_pre_lower <- as.data.frame(cbind(PC7_1.0sd_pre$x, PC7_1.0sd_pre$lower))
names(PC7_1.0sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_1.0sd_pre_lower$age, PC7_1.0sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_1.0sd_pre_lower, indices) {
  d <- PC7_1.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_1.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_1.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_1.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_1.0sd_pre_lower_50)[1] <- "PC7"

PC7_1.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_1.0sd_pre_lower_all)[1:ncol(PC7_1.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_1.0sd_pre_lower_all$group <- seq(1,(nrow(PC7_1.0sd_pre_lower_all)), 1)

PC7_1.0sd_pre_lower_all <- reshape2::melt(PC7_1.0sd_pre_lower_all, id.vars = "group")
names(PC7_1.0sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_1.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.8 sigma

nsig <- 0.8

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.8sd_pre_upper <- as.data.frame(cbind(PC7_0.8sd_pre$x, PC7_0.8sd_pre$upper))
names(PC7_0.8sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.8sd_pre_upper$age, PC7_0.8sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_0.8sd_pre_upper, indices) {
  d <- PC7_0.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_0.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_0.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_0.8sd_pre_upper_50)[1] <- "PC7"

PC7_0.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_0.8sd_pre_upper_all)[1:ncol(PC7_0.8sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_0.8sd_pre_upper_all$group <- seq(1,(nrow(PC7_0.8sd_pre_upper_all)), 1)

PC7_0.8sd_pre_upper_all <- reshape2::melt(PC7_0.8sd_pre_upper_all, id.vars = "group")
names(PC7_0.8sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_0.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.8sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.8sd_pre_lower <- as.data.frame(cbind(PC7_0.8sd_pre$x, PC7_0.8sd_pre$lower))
names(PC7_0.8sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.8sd_pre_lower$age, PC7_0.8sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_0.8sd_pre_lower, indices) {
  d <- PC7_0.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_0.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_0.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_0.8sd_pre_lower_50)[1] <- "PC7"

PC7_0.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_0.8sd_pre_lower_all)[1:ncol(PC7_0.8sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_0.8sd_pre_lower_all$group <- seq(1,(nrow(PC7_0.8sd_pre_lower_all)), 1)

PC7_0.8sd_pre_lower_all <- reshape2::melt(PC7_0.8sd_pre_lower_all, id.vars = "group")
names(PC7_0.8sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_0.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.6 sigma

nsig <- 0.6

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.6sd_pre_upper <- as.data.frame(cbind(PC7_0.6sd_pre$x, PC7_0.6sd_pre$upper))
names(PC7_0.6sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.6sd_pre_upper$age, PC7_0.6sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_0.6sd_pre_upper, indices) {
  d <- PC7_0.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_0.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_0.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_0.6sd_pre_upper_50)[1] <- "PC7"

PC7_0.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_0.6sd_pre_upper_all)[1:ncol(PC7_0.6sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_0.6sd_pre_upper_all$group <- seq(1,(nrow(PC7_0.6sd_pre_upper_all)), 1)

PC7_0.6sd_pre_upper_all <- reshape2::melt(PC7_0.6sd_pre_upper_all, id.vars = "group")
names(PC7_0.6sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_0.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.6sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.6sd_pre_lower <- as.data.frame(cbind(PC7_0.6sd_pre$x, PC7_0.6sd_pre$lower))
names(PC7_0.6sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.6sd_pre_lower$age, PC7_0.6sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_0.6sd_pre_lower, indices) {
  d <- PC7_0.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_0.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_0.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_0.6sd_pre_lower_50)[1] <- "PC7"

PC7_0.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_0.6sd_pre_lower_all)[1:ncol(PC7_0.6sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_0.6sd_pre_lower_all$group <- seq(1,(nrow(PC7_0.6sd_pre_lower_all)), 1)

PC7_0.6sd_pre_lower_all <- reshape2::melt(PC7_0.6sd_pre_lower_all, id.vars = "group")
names(PC7_0.6sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_0.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.4 sigma

nsig <- 0.4

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.4sd_pre_upper <- as.data.frame(cbind(PC7_0.4sd_pre$x, PC7_0.4sd_pre$upper))
names(PC7_0.4sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.4sd_pre_upper$age, PC7_0.4sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_0.4sd_pre_upper, indices) {
  d <- PC7_0.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_0.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_0.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_0.4sd_pre_upper_50)[1] <- "PC7"

PC7_0.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_0.4sd_pre_upper_all)[1:ncol(PC7_0.4sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_0.4sd_pre_upper_all$group <- seq(1,(nrow(PC7_0.4sd_pre_upper_all)), 1)

PC7_0.4sd_pre_upper_all <- reshape2::melt(PC7_0.4sd_pre_upper_all, id.vars = "group")
names(PC7_0.4sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_0.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.4sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.4sd_pre_lower <- as.data.frame(cbind(PC7_0.4sd_pre$x, PC7_0.4sd_pre$lower))
names(PC7_0.4sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.4sd_pre_lower$age, PC7_0.4sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_0.4sd_pre_lower, indices) {
  d <- PC7_0.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_0.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_0.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_0.4sd_pre_lower_50)[1] <- "PC7"

PC7_0.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_0.4sd_pre_lower_all)[1:ncol(PC7_0.4sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_0.4sd_pre_lower_all$group <- seq(1,(nrow(PC7_0.4sd_pre_lower_all)), 1)

PC7_0.4sd_pre_lower_all <- reshape2::melt(PC7_0.4sd_pre_lower_all, id.vars = "group")
names(PC7_0.4sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_0.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.2 sigma

nsig <- 0.2

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.2sd_pre_upper <- as.data.frame(cbind(PC7_0.2sd_pre$x, PC7_0.2sd_pre$upper))
names(PC7_0.2sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.2sd_pre_upper$age, PC7_0.2sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_0.2sd_pre_upper, indices) {
  d <- PC7_0.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_0.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_0.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_0.2sd_pre_upper_50)[1] <- "PC7"

PC7_0.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_0.2sd_pre_upper_all)[1:ncol(PC7_0.2sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_0.2sd_pre_upper_all$group <- seq(1,(nrow(PC7_0.2sd_pre_upper_all)), 1)

PC7_0.2sd_pre_upper_all <- reshape2::melt(PC7_0.2sd_pre_upper_all, id.vars = "group")
names(PC7_0.2sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_0.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.2sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.2sd_pre_lower <- as.data.frame(cbind(PC7_0.2sd_pre$x, PC7_0.2sd_pre$lower))
names(PC7_0.2sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.2sd_pre_lower$age, PC7_0.2sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_0.2sd_pre_lower, indices) {
  d <- PC7_0.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_0.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_0.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_0.2sd_pre_lower_50)[1] <- "PC7"

PC7_0.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_0.2sd_pre_lower_all)[1:ncol(PC7_0.2sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_0.2sd_pre_lower_all$group <- seq(1,(nrow(PC7_0.2sd_pre_lower_all)), 1)

PC7_0.2sd_pre_lower_all <- reshape2::melt(PC7_0.2sd_pre_lower_all, id.vars = "group")
names(PC7_0.2sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_0.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.0 sigma

nsig <- 0.0

# PC7 upper

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.0sd_pre_upper <- as.data.frame(cbind(PC7_0.0sd_pre$x, PC7_0.0sd_pre$upper))
names(PC7_0.0sd_pre_upper)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.0sd_pre_upper$age, PC7_0.0sd_pre_upper$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- loess.predict[["pars"]][["span"]]

# boot upper
boot_fn <- function(PC7_0.0sd_pre_upper, indices) {
  d <- PC7_0.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_upper, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(PC7_0.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
PC7_0.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(PC7_0.0sd_pre_upper_50)[1] <- "PC7"

PC7_0.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(PC7_0.0sd_pre_upper_all)[1:ncol(PC7_0.0sd_pre_upper_all)] <- seq(540,4000, 10)
PC7_0.0sd_pre_upper_all$group <- seq(1,(nrow(PC7_0.0sd_pre_upper_all)), 1)

PC7_0.0sd_pre_upper_all <- reshape2::melt(PC7_0.0sd_pre_upper_all, id.vars = "group")
names(PC7_0.0sd_pre_upper_all)[2:3] <- c("age", "PC7")

PC7_0.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# PC7 lower

# CV span
loess.predict <- loess.as(coords.all_pre$age, coords.all_pre$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

PC7_0.0sd_pre <- msir::loess.sd(coords.all_pre$age, coords.all_pre$PC7, degree=1, nsigma = nsig, span = span)
PC7_0.0sd_pre_lower <- as.data.frame(cbind(PC7_0.0sd_pre$x, PC7_0.0sd_pre$lower))
names(PC7_0.0sd_pre_lower)[1:2] <- c("age","PC7")

# CV spans
loess.predict <- loess.as(PC7_0.0sd_pre_lower$age, PC7_0.0sd_pre_lower$PC7, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- loess.predict[["pars"]][["span"]]

# boot lower
boot_fn <- function(PC7_0.0sd_pre_lower, indices) {
  d <- PC7_0.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(PC7 ~ age, d, span = span_lower, degree = 1)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == (nsig*-1)])/speed

loess_boot <- boot(PC7_0.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# median smooth
PC7_0.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
PC7_0.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(PC7_0.0sd_pre_lower_50)[1] <- "PC7"

PC7_0.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(PC7_0.0sd_pre_lower_all)[1:ncol(PC7_0.0sd_pre_lower_all)] <- seq(540,4000, 10)
PC7_0.0sd_pre_lower_all$group <- seq(1,(nrow(PC7_0.0sd_pre_lower_all)), 1)

PC7_0.0sd_pre_lower_all <- reshape2::melt(PC7_0.0sd_pre_lower_all, id.vars = "group")
names(PC7_0.0sd_pre_lower_all)[2:3] <- c("age", "PC7")

PC7_0.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

b <- ggplot() + 
  #geom_line(data = PC7_2sd_pre_upper_all, aes(age2,PC7, group = group)) +
  #geom_line(data = PC7_2sd_pre_lower_all, aes(age2,PC7, group = group)) +
  geom_line(data = PC7_2sd_pre_upper_50, aes(age,PC7), colour = 'blue') +
  geom_line(data = PC7_2sd_pre_lower_50, aes(age,PC7), colour = 'blue') +
  geom_line(data = PC7_1.6sd_pre_upper_50, aes(age,PC7), colour = 'green') +
  geom_line(data = PC7_1.6sd_pre_lower_50, aes(age,PC7), colour = 'green') +
  geom_line(data = PC7_1.4sd_pre_upper_50, aes(age,PC7), colour = 'orange') +
  geom_line(data = PC7_1.4sd_pre_lower_50, aes(age,PC7), colour = 'orange') +
  geom_line(data = PC7_1.2sd_pre_upper_50, aes(age,PC7), colour = 'pink') +
  geom_line(data = PC7_1.2sd_pre_lower_50, aes(age,PC7), colour = 'pink') +
  geom_line(data = PC7_1.0sd_pre_upper_50, aes(age,PC7), colour = 'brown') +
  geom_line(data = PC7_1.0sd_pre_lower_50, aes(age,PC7), colour = 'brown') +
  geom_line(data = PC7_0.8sd_pre_upper_50, aes(age,PC7), colour = 'yellow') +
  geom_line(data = PC7_0.8sd_pre_lower_50, aes(age,PC7), colour = 'yellow') +
  geom_line(data = PC7_0.6sd_pre_upper_50, aes(age,PC7), colour = 'red') +
  geom_line(data = PC7_0.6sd_pre_lower_50, aes(age,PC7), colour = 'red') +
  geom_line(data = PC7_0.4sd_pre_upper_50, aes(age,PC7), colour = 'blue') +
  geom_line(data = PC7_0.4sd_pre_lower_50, aes(age,PC7), colour = 'blue') +
  geom_line(data = PC7_0.2sd_pre_upper_50, aes(age,PC7), colour = 'green') +
  geom_line(data = PC7_0.2sd_pre_lower_50, aes(age,PC7), colour = 'green') +
  geom_line(data = PC7_0.0sd_pre_upper_50, aes(age,PC7), colour = 'orange') +
  geom_line(data = PC7_0.0sd_pre_lower_50, aes(age,PC7), colour = 'orange') +
  geom_point(data = coords.all, aes(age, PC7)) + 
  #geom_point(data = PC7_2sd_pre_upper, aes(age, PC7), colour = 'red') +
  scale_x_reverse(limits = c(4000,540)) + theme_bw()

grid.arrange(b,a,ncol=2)

##### Median ribbons and raw plots ####

# median ribbons (not presented in manuscript)

# replace PC as necessary

phan_2sd <- cbind(PC3_2sd_phan_upper_50,PC3_2sd_phan_lower_50$PC3)
names(phan_2sd)[1:3] <- c("upper","age","lower")
phan_1.8sd <- cbind(PC3_1.8sd_phan_upper_50,PC3_1.8sd_phan_lower_50$PC3)
names(phan_1.8sd)[1:3] <- c("upper","age","lower")
phan_1.6sd <- cbind(PC3_1.6sd_phan_upper_50,PC3_1.6sd_phan_lower_50$PC3)
names(phan_1.6sd)[1:3] <- c("upper","age","lower")
phan_1.4sd <- cbind(PC3_1.4sd_phan_upper_50,PC3_1.4sd_phan_lower_50$PC3)
names(phan_1.4sd)[1:3] <- c("upper","age","lower")
phan_1.2sd <- cbind(PC3_1.2sd_phan_upper_50,PC3_1.2sd_phan_lower_50$PC3)
names(phan_1.2sd)[1:3] <- c("upper","age","lower")
phan_1.0sd <- cbind(PC3_1.0sd_phan_upper_50,PC3_1.0sd_phan_lower_50$PC3)
names(phan_1.0sd)[1:3] <- c("upper","age","lower")
phan_0.8sd <- cbind(PC3_0.8sd_phan_upper_50,PC3_0.8sd_phan_lower_50$PC3)
names(phan_0.8sd)[1:3] <- c("upper","age","lower")
phan_0.6sd <- cbind(PC3_0.6sd_phan_upper_50,PC3_0.6sd_phan_lower_50$PC3)
names(phan_0.6sd)[1:3] <- c("upper","age","lower")
phan_0.4sd <- cbind(PC3_0.4sd_phan_upper_50,PC3_0.4sd_phan_lower_50$PC3)
names(phan_0.4sd)[1:3] <- c("upper","age","lower")
phan_0.2sd <- cbind(PC3_0.2sd_phan_upper_50,PC3_0.2sd_phan_lower_50$PC3)
names(phan_0.2sd)[1:3] <- c("upper","age","lower")

pre_2sd <- cbind(PC3_2sd_pre_upper_50,PC3_2sd_pre_lower_50$PC3)
names(pre_2sd)[1:3] <- c("upper","age","lower")
pre_1.8sd <- cbind(PC3_1.8sd_pre_upper_50,PC3_1.8sd_pre_lower_50$PC3)
names(pre_1.8sd)[1:3] <- c("upper","age","lower")
pre_1.6sd <- cbind(PC3_1.6sd_pre_upper_50,PC3_1.6sd_pre_lower_50$PC3)
names(pre_1.6sd)[1:3] <- c("upper","age","lower")
pre_1.4sd <- cbind(PC3_1.4sd_pre_upper_50,PC3_1.4sd_pre_lower_50$PC3)
names(pre_1.4sd)[1:3] <- c("upper","age","lower")
pre_1.2sd <- cbind(PC3_1.2sd_pre_upper_50,PC3_1.2sd_pre_lower_50$PC3)
names(pre_1.2sd)[1:3] <- c("upper","age","lower")
pre_1.0sd <- cbind(PC3_1.0sd_pre_upper_50,PC3_1.0sd_pre_lower_50$PC3)
names(pre_1.0sd)[1:3] <- c("upper","age","lower")
pre_0.8sd <- cbind(PC3_0.8sd_pre_upper_50,PC3_0.8sd_pre_lower_50$PC3)
names(pre_0.8sd)[1:3] <- c("upper","age","lower")
pre_0.6sd <- cbind(PC3_0.6sd_pre_upper_50,PC3_0.6sd_pre_lower_50$PC3)
names(pre_0.6sd)[1:3] <- c("upper","age","lower")
pre_0.4sd <- cbind(PC3_0.4sd_pre_upper_50,PC3_0.4sd_pre_lower_50$PC3)
names(pre_0.4sd)[1:3] <- c("upper","age","lower")
pre_0.2sd <- cbind(PC3_0.2sd_pre_upper_50,PC3_0.2sd_pre_lower_50$PC3)
names(pre_0.2sd)[1:3] <- c("upper","age","lower")

a <- ggplot() +
  geom_ribbon(data = phan_2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_1.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_1.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_1.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_1.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_1.0sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_0.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_0.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_0.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_0.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_line(data = PC3_0.0sd_phan_lower_50, aes(age,PC3), colour = 'blue', size = 1) +
  geom_point(data = coords.all, aes(age, PC3), size = 0.5) + 
  scale_x_reverse(limits = c(540,0)) + theme_bw()

b <- ggplot() +
  geom_ribbon(data = pre_2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_1.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_1.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_1.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_1.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_1.0sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_0.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_0.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_0.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_0.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_line(data = PC3_0.0sd_pre_lower_50, aes(age,PC3), colour = 'blue', size = 1) +
  geom_point(data = coords.all, aes(age, PC3), size = 0.5) + 
  scale_x_reverse(limits = c(4000,540)) + theme_bw() + 
  scale_y_continuous(limits = c(-9.25,5.25))

grid.arrange(b,a,ncol=2)

# raw plots (select PC as necessary)

a <- ggplot() + 
  geom_line(data = PC2_2sd_phan_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_2sd_phan_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.6sd_phan_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.6sd_phan_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.4sd_phan_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.4sd_phan_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.2sd_phan_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.2sd_phan_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.0sd_phan_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.0sd_phan_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.8sd_phan_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.8sd_phan_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.6sd_phan_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.6sd_phan_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.4sd_phan_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.4sd_phan_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.2sd_phan_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.2sd_phan_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.0sd_phan_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.0sd_phan_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_point(data = coords.all, aes(age, PC2)) + 
  scale_x_reverse(limits = c(540,0)) + theme_bw()

b <- ggplot() + 
  geom_line(data = PC2_2sd_pre_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_2sd_pre_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.6sd_pre_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.6sd_pre_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.4sd_pre_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.4sd_pre_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.2sd_pre_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.2sd_pre_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.0sd_pre_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_1.0sd_pre_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.8sd_pre_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.8sd_pre_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.6sd_pre_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.6sd_pre_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.4sd_pre_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.4sd_pre_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.2sd_pre_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.2sd_pre_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.0sd_pre_upper_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_line(data = PC2_0.0sd_pre_lower_all, aes(age2,PC2, group = group), colour = 'blue') +
  geom_point(data = coords.all, aes(age, PC2)) + 
  scale_x_reverse(limits = c(3600,540)) + theme_bw()

grid.arrange(b,a,ncol = 2)

##### Bind outputs ####
# note 0.0sd lower and upper outputs are duplicates of the median line
# therefore 0.0sd outputs only bound once

# Phanerozoic
# PC1 Phanerozoic
PC1_phan_comp <- rbind(PC1_2sd_phan_upper_all,PC1_1.8sd_phan_upper_all,
                       PC1_1.6sd_phan_upper_all,PC1_1.4sd_phan_upper_all,
                       PC1_1.2sd_phan_upper_all,PC1_1.0sd_phan_upper_all,
                       PC1_0.8sd_phan_upper_all,PC1_0.6sd_phan_upper_all,
                       PC1_0.4sd_phan_upper_all,PC1_0.2sd_phan_upper_all,
                       PC1_2sd_phan_lower_all,PC1_1.8sd_phan_lower_all,
                       PC1_1.6sd_phan_lower_all,PC1_1.4sd_phan_lower_all,
                       PC1_1.2sd_phan_lower_all,PC1_1.0sd_phan_lower_all,
                       PC1_0.8sd_phan_lower_all,PC1_0.6sd_phan_lower_all,
                       PC1_0.4sd_phan_lower_all,PC1_0.2sd_phan_lower_all,
                       PC1_0.0sd_phan_lower_all)

# PC2 Phanerozoic
PC2_phan_comp <- rbind(PC2_2sd_phan_upper_all,PC2_1.8sd_phan_upper_all,
                       PC2_1.6sd_phan_upper_all,PC2_1.4sd_phan_upper_all,
                       PC2_1.2sd_phan_upper_all,PC2_1.0sd_phan_upper_all,
                       PC2_0.8sd_phan_upper_all,PC2_0.6sd_phan_upper_all,
                       PC2_0.4sd_phan_upper_all,PC2_0.2sd_phan_upper_all,
                       PC2_2sd_phan_lower_all,PC2_1.8sd_phan_lower_all,
                       PC2_1.6sd_phan_lower_all,PC2_1.4sd_phan_lower_all,
                       PC2_1.2sd_phan_lower_all,PC2_1.0sd_phan_lower_all,
                       PC2_0.8sd_phan_lower_all,PC2_0.6sd_phan_lower_all,
                       PC2_0.4sd_phan_lower_all,PC2_0.2sd_phan_lower_all,
                       PC2_0.0sd_phan_lower_all)

# PC3 Phanerozoic
PC3_phan_comp <- rbind(PC3_2sd_phan_upper_all,PC3_1.8sd_phan_upper_all,
                       PC3_1.6sd_phan_upper_all,PC3_1.4sd_phan_upper_all,
                       PC3_1.2sd_phan_upper_all,PC3_1.0sd_phan_upper_all,
                       PC3_0.8sd_phan_upper_all,PC3_0.6sd_phan_upper_all,
                       PC3_0.4sd_phan_upper_all,PC3_0.2sd_phan_upper_all,
                       PC3_2sd_phan_lower_all,PC3_1.8sd_phan_lower_all,
                       PC3_1.6sd_phan_lower_all,PC3_1.4sd_phan_lower_all,
                       PC3_1.2sd_phan_lower_all,PC3_1.0sd_phan_lower_all,
                       PC3_0.8sd_phan_lower_all,PC3_0.6sd_phan_lower_all,
                       PC3_0.4sd_phan_lower_all,PC3_0.2sd_phan_lower_all,
                       PC3_0.0sd_phan_lower_all)

# PC4 Phanerozoic
PC4_phan_comp <- rbind(PC4_2sd_phan_upper_all,PC4_1.8sd_phan_upper_all,
                       PC4_1.6sd_phan_upper_all,PC4_1.4sd_phan_upper_all,
                       PC4_1.2sd_phan_upper_all,PC4_1.0sd_phan_upper_all,
                       PC4_0.8sd_phan_upper_all,PC4_0.6sd_phan_upper_all,
                       PC4_0.4sd_phan_upper_all,PC4_0.2sd_phan_upper_all,
                       PC4_2sd_phan_lower_all,PC4_1.8sd_phan_lower_all,
                       PC4_1.6sd_phan_lower_all,PC4_1.4sd_phan_lower_all,
                       PC4_1.2sd_phan_lower_all,PC4_1.0sd_phan_lower_all,
                       PC4_0.8sd_phan_lower_all,PC4_0.6sd_phan_lower_all,
                       PC4_0.4sd_phan_lower_all,PC4_0.2sd_phan_lower_all,
                       PC4_0.0sd_phan_lower_all)

# PC5 Phanerozoic
PC5_phan_comp <- rbind(PC5_2sd_phan_upper_all,PC5_1.8sd_phan_upper_all,
                       PC5_1.6sd_phan_upper_all,PC5_1.4sd_phan_upper_all,
                       PC5_1.2sd_phan_upper_all,PC5_1.0sd_phan_upper_all,
                       PC5_0.8sd_phan_upper_all,PC5_0.6sd_phan_upper_all,
                       PC5_0.4sd_phan_upper_all,PC5_0.2sd_phan_upper_all,
                       PC5_2sd_phan_lower_all,PC5_1.8sd_phan_lower_all,
                       PC5_1.6sd_phan_lower_all,PC5_1.4sd_phan_lower_all,
                       PC5_1.2sd_phan_lower_all,PC5_1.0sd_phan_lower_all,
                       PC5_0.8sd_phan_lower_all,PC5_0.6sd_phan_lower_all,
                       PC5_0.4sd_phan_lower_all,PC5_0.2sd_phan_lower_all,
                       PC5_0.0sd_phan_lower_all)

# PC6 Phanerozoic
PC6_phan_comp <- rbind(PC6_2sd_phan_upper_all,PC6_1.8sd_phan_upper_all,
                       PC6_1.6sd_phan_upper_all,PC6_1.4sd_phan_upper_all,
                       PC6_1.2sd_phan_upper_all,PC6_1.0sd_phan_upper_all,
                       PC6_0.8sd_phan_upper_all,PC6_0.6sd_phan_upper_all,
                       PC6_0.4sd_phan_upper_all,PC6_0.2sd_phan_upper_all,
                       PC6_2sd_phan_lower_all,PC6_1.8sd_phan_lower_all,
                       PC6_1.6sd_phan_lower_all,PC6_1.4sd_phan_lower_all,
                       PC6_1.2sd_phan_lower_all,PC6_1.0sd_phan_lower_all,
                       PC6_0.8sd_phan_lower_all,PC6_0.6sd_phan_lower_all,
                       PC6_0.4sd_phan_lower_all,PC6_0.2sd_phan_lower_all,
                       PC6_0.0sd_phan_lower_all)

# PC7 Phanerozoic
PC7_phan_comp <- rbind(PC7_2sd_phan_upper_all,PC7_1.8sd_phan_upper_all,
                       PC7_1.6sd_phan_upper_all,PC7_1.4sd_phan_upper_all,
                       PC7_1.2sd_phan_upper_all,PC7_1.0sd_phan_upper_all,
                       PC7_0.8sd_phan_upper_all,PC7_0.6sd_phan_upper_all,
                       PC7_0.4sd_phan_upper_all,PC7_0.2sd_phan_upper_all,
                       PC7_2sd_phan_lower_all,PC7_1.8sd_phan_lower_all,
                       PC7_1.6sd_phan_lower_all,PC7_1.4sd_phan_lower_all,
                       PC7_1.2sd_phan_lower_all,PC7_1.0sd_phan_lower_all,
                       PC7_0.8sd_phan_lower_all,PC7_0.6sd_phan_lower_all,
                       PC7_0.4sd_phan_lower_all,PC7_0.2sd_phan_lower_all,
                       PC7_0.0sd_phan_lower_all)

phan_comp <- cbind(PC1_phan_comp,PC2_phan_comp,PC3_phan_comp,PC4_phan_comp,PC5_phan_comp,PC6_phan_comp,PC7_phan_comp)
phan_comp <- phan_comp[,c(1:3,7,11,15,19,23,27)]
names(phan_comp)[1:2] <- c("group","age")
phan_comp <- phan_comp[complete.cases(phan_comp[,]),] 

# Precambrian
# PC1 Phanerozoic
PC1_pre_comp <- rbind(PC1_2sd_pre_upper_all,PC1_1.8sd_pre_upper_all,
                      PC1_1.6sd_pre_upper_all,PC1_1.4sd_pre_upper_all,
                      PC1_1.2sd_pre_upper_all,PC1_1.0sd_pre_upper_all,
                      PC1_0.8sd_pre_upper_all,PC1_0.6sd_pre_upper_all,
                      PC1_0.4sd_pre_upper_all,PC1_0.2sd_pre_upper_all,
                      PC1_2sd_pre_lower_all,PC1_1.8sd_pre_lower_all,
                      PC1_1.6sd_pre_lower_all,PC1_1.4sd_pre_lower_all,
                      PC1_1.2sd_pre_lower_all,PC1_1.0sd_pre_lower_all,
                      PC1_0.8sd_pre_lower_all,PC1_0.6sd_pre_lower_all,
                      PC1_0.4sd_pre_lower_all,PC1_0.2sd_pre_lower_all,
                      PC1_0.0sd_pre_lower_all)

# PC2 Precambrian
PC2_pre_comp <- rbind(PC2_2sd_pre_upper_all,PC2_1.8sd_pre_upper_all,
                      PC2_1.6sd_pre_upper_all,PC2_1.4sd_pre_upper_all,
                      PC2_1.2sd_pre_upper_all,PC2_1.0sd_pre_upper_all,
                      PC2_0.8sd_pre_upper_all,PC2_0.6sd_pre_upper_all,
                      PC2_0.4sd_pre_upper_all,PC2_0.2sd_pre_upper_all,
                      PC2_2sd_pre_lower_all,PC2_1.8sd_pre_lower_all,
                      PC2_1.6sd_pre_lower_all,PC2_1.4sd_pre_lower_all,
                      PC2_1.2sd_pre_lower_all,PC2_1.0sd_pre_lower_all,
                      PC2_0.8sd_pre_lower_all,PC2_0.6sd_pre_lower_all,
                      PC2_0.4sd_pre_lower_all,PC2_0.2sd_pre_lower_all,
                      PC2_0.0sd_pre_lower_all)

# PC3 Precambrian
PC3_pre_comp <- rbind(PC3_2sd_pre_upper_all,PC3_1.8sd_pre_upper_all,
                      PC3_1.6sd_pre_upper_all,PC3_1.4sd_pre_upper_all,
                      PC3_1.2sd_pre_upper_all,PC3_1.0sd_pre_upper_all,
                      PC3_0.8sd_pre_upper_all,PC3_0.6sd_pre_upper_all,
                      PC3_0.4sd_pre_upper_all,PC3_0.2sd_pre_upper_all,
                      PC3_2sd_pre_lower_all,PC3_1.8sd_pre_lower_all,
                      PC3_1.6sd_pre_lower_all,PC3_1.4sd_pre_lower_all,
                      PC3_1.2sd_pre_lower_all,PC3_1.0sd_pre_lower_all,
                      PC3_0.8sd_pre_lower_all,PC3_0.6sd_pre_lower_all,
                      PC3_0.4sd_pre_lower_all,PC3_0.2sd_pre_lower_all,
                      PC3_0.0sd_pre_lower_all)

# PC4 Precambrian
PC4_pre_comp <- rbind(PC4_2sd_pre_upper_all,PC4_1.8sd_pre_upper_all,
                      PC4_1.6sd_pre_upper_all,PC4_1.4sd_pre_upper_all,
                      PC4_1.2sd_pre_upper_all,PC4_1.0sd_pre_upper_all,
                      PC4_0.8sd_pre_upper_all,PC4_0.6sd_pre_upper_all,
                      PC4_0.4sd_pre_upper_all,PC4_0.2sd_pre_upper_all,
                      PC4_2sd_pre_lower_all,PC4_1.8sd_pre_lower_all,
                      PC4_1.6sd_pre_lower_all,PC4_1.4sd_pre_lower_all,
                      PC4_1.2sd_pre_lower_all,PC4_1.0sd_pre_lower_all,
                      PC4_0.8sd_pre_lower_all,PC4_0.6sd_pre_lower_all,
                      PC4_0.4sd_pre_lower_all,PC4_0.2sd_pre_lower_all,
                      PC4_0.0sd_pre_lower_all)

# PC5 Precambrian
PC5_pre_comp <- rbind(PC5_2sd_pre_upper_all,PC5_1.8sd_pre_upper_all,
                      PC5_1.6sd_pre_upper_all,PC5_1.4sd_pre_upper_all,
                      PC5_1.2sd_pre_upper_all,PC5_1.0sd_pre_upper_all,
                      PC5_0.8sd_pre_upper_all,PC5_0.6sd_pre_upper_all,
                      PC5_0.4sd_pre_upper_all,PC5_0.2sd_pre_upper_all,
                      PC5_2sd_pre_lower_all,PC5_1.8sd_pre_lower_all,
                      PC5_1.6sd_pre_lower_all,PC5_1.4sd_pre_lower_all,
                      PC5_1.2sd_pre_lower_all,PC5_1.0sd_pre_lower_all,
                      PC5_0.8sd_pre_lower_all,PC5_0.6sd_pre_lower_all,
                      PC5_0.4sd_pre_lower_all,PC5_0.2sd_pre_lower_all,
                      PC5_0.0sd_pre_lower_all)

# PC6 Precambrian
PC6_pre_comp <- rbind(PC6_2sd_pre_upper_all,PC6_1.8sd_pre_upper_all,
                      PC6_1.6sd_pre_upper_all,PC6_1.4sd_pre_upper_all,
                      PC6_1.2sd_pre_upper_all,PC6_1.0sd_pre_upper_all,
                      PC6_0.8sd_pre_upper_all,PC6_0.6sd_pre_upper_all,
                      PC6_0.4sd_pre_upper_all,PC6_0.2sd_pre_upper_all,
                      PC6_2sd_pre_lower_all,PC6_1.8sd_pre_lower_all,
                      PC6_1.6sd_pre_lower_all,PC6_1.4sd_pre_lower_all,
                      PC6_1.2sd_pre_lower_all,PC6_1.0sd_pre_lower_all,
                      PC6_0.8sd_pre_lower_all,PC6_0.6sd_pre_lower_all,
                      PC6_0.4sd_pre_lower_all,PC6_0.2sd_pre_lower_all,
                      PC6_0.0sd_pre_lower_all)

# PC7 Precambrian
PC7_pre_comp <- rbind(PC7_2sd_pre_upper_all,PC7_1.8sd_pre_upper_all,
                      PC7_1.6sd_pre_upper_all,PC7_1.4sd_pre_upper_all,
                      PC7_1.2sd_pre_upper_all,PC7_1.0sd_pre_upper_all,
                      PC7_0.8sd_pre_upper_all,PC7_0.6sd_pre_upper_all,
                      PC7_0.4sd_pre_upper_all,PC7_0.2sd_pre_upper_all,
                      PC7_2sd_pre_lower_all,PC7_1.8sd_pre_lower_all,
                      PC7_1.6sd_pre_lower_all,PC7_1.4sd_pre_lower_all,
                      PC7_1.2sd_pre_lower_all,PC7_1.0sd_pre_lower_all,
                      PC7_0.8sd_pre_lower_all,PC7_0.6sd_pre_lower_all,
                      PC7_0.4sd_pre_lower_all,PC7_0.2sd_pre_lower_all,
                      PC7_0.0sd_pre_lower_all)


pre_comp <- cbind(PC1_pre_comp,PC2_pre_comp,PC3_pre_comp,PC4_pre_comp,PC5_pre_comp,PC6_pre_comp,PC7_pre_comp)
pre_comp <- pre_comp[,c(1:3,7,11,15,19,23,27)]
names(pre_comp)[1:2] <- c("group","age")
pre_comp <- pre_comp[complete.cases(pre_comp[,]),] 

# clean environment

rm(list= ls()[!(ls() %in% c('pre_comp', 'phan_comp', 'zones','zone.names', 'coords.all', 'coords.all2', 'n','speed'))])

###### Train model ######

coords.train <- coords.all2[,c(1:7,ncol(coords.all2)-1)]

validation_index <- createDataPartition(coords.train$HCA.clusters, p=0.80, list=FALSE)
validation <- coords.train[-validation_index,]
MLdataset <- coords.train[validation_index,]

control <- trainControl(method="cv", number=10)
metric <- "Accuracy"

# a) linear algorithms
set.seed(7)
fit.lda <- train(HCA.clusters~., data=MLdataset, method="lda", metric=metric, trControl=control)
# b) nonlinear algorithms
# CART
set.seed(7)
fit.cart <- train(HCA.clusters~., data=MLdataset, method="rpart", metric=metric, trControl=control)
# kNN
set.seed(7)
fit.knn <- train(HCA.clusters~., data=MLdataset, method="knn", metric=metric, trControl=control)
# c) advanced algorithms
# SVM
set.seed(7)
fit.svm <- train(HCA.clusters~., data=MLdataset, method="svmRadial", metric=metric, trControl=control)
# Random Forest
set.seed(7)
fit.rf <- train(HCA.clusters~., data=MLdataset, method="rf", metric=metric, trControl=control)

results <- resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf))
summary(results)

# compare models
dotplot(results)

fit <- fit.knn # select model

predictions <- predict(fit, validation) 
confusionMatrix(predictions, as.factor(validation$HCA.clusters))

##### Shuffle, deploy model and aggregate ####

s <- 100 # shuffle interations 

datalist = list()

for (i in 1:s) {
  dat <- phan_comp %>%
    group_by(age) %>%
    mutate(PC1 = sample(PC1, replace = FALSE),
           PC2 = sample(PC2, replace = FALSE),
           PC3 = sample(PC3, replace = FALSE),
           PC4 = sample(PC4, replace = FALSE),
           PC5 = sample(PC5, replace = FALSE),
           PC6 = sample(PC6, replace = FALSE),
           PC7 = sample(PC7, replace = FALSE))
  datML <- dat[,3:(ncol(dat))]
  dat_predictions <- predict(fit, datML)
  dat_predictions <- as.data.frame(dat_predictions)
  dat <- cbind(dat, dat_predictions)
  dat <- dat %>% 
    group_by(age) %>%
    mutate(Type.5 = mean(dat_predictions == "Type.5"),
           Type.1 = mean(dat_predictions == "Type.1"),
           Type.2 = mean(dat_predictions == "Type.2"),
           Type.3 = mean(dat_predictions == "Type.3"),
           Type.4 = mean(dat_predictions == "Type.4"))
  dat <- dat[!duplicated(dat[,2]),]
  dat <- dat[,c(2,11:ncol(dat))] 
  dat$age <- as.numeric(as.character(dat$age))
  dat$i <- i
  datalist[[i]] <- dat
}

phan_results = do.call(rbind, datalist)
#write.csv(phan_results, "phanerozoic_pyrite_types.csv") # optional export

# melt
phan_results.melt <- reshape2::melt(phan_results, id.vars = c("age","i"))

datalist = list()

for (i in 1:s) {
  dat <- pre_comp %>%
    group_by(age) %>%
    mutate(PC1 = sample(PC1, replace = FALSE),
           PC2 = sample(PC2, replace = FALSE),
           PC3 = sample(PC3, replace = FALSE),
           PC4 = sample(PC4, replace = FALSE),
           PC5 = sample(PC5, replace = FALSE),
           PC6 = sample(PC6, replace = FALSE),
           PC7 = sample(PC7, replace = FALSE))
  datML <- dat[,3:(ncol(dat))]
  dat_predictions <- predict(fit, datML)
  dat_predictions <- as.data.frame(dat_predictions)
  dat <- cbind(dat, dat_predictions)
  dat <- dat %>% 
    group_by(age) %>%
    mutate(Type.5 = mean(dat_predictions == "Type.5"),
           Type.1 = mean(dat_predictions == "Type.1"),
           Type.2 = mean(dat_predictions == "Type.2"),
           Type.3 = mean(dat_predictions == "Type.3"),
           Type.4 = mean(dat_predictions == "Type.4"))
  dat <- dat[!duplicated(dat[,2]),]
  dat <- dat[,c(2,11:ncol(dat))] 
  dat$age <- as.numeric(as.character(dat$age))
  dat <- as.data.frame(dat)
  dat$i <- i
  datalist[[i]] <- dat
}

pre_results = do.call(rbind, datalist)
#write.csv(pre_results, "precambrian_pyrite_types.csv") # optional export

# melt
pre_results.melt <- reshape2::melt(pre_results, id.vars = c("age","i"))

##### Downsample, export and plot results ####

# Phanerozoic age differences
diffs <-  coords.all$age
diffs <- as.data.frame(diffs[order(diffs)])
names(diffs)[1] <- "Age"
diffs <- as.data.frame(diffs[!duplicated(diffs$Age),])
diffs_df <- round(diffs[-1,]) - round(diffs[-nrow(diffs),])
diffs_df <- as.data.frame(diffs_df)
diffs <- cbind(diffs[1:nrow(diffs)-1,],diffs_df)
names(diffs)[1] <- "Age"
names(diffs)[2] <- "Diff"
diffs <- subset(diffs, Diff > 0)
diffs$Diff2 <- diffs$Diff/2
diffs$edge <- round(diffs$Diff2+diffs$Age)
diffs$span <- c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4])
diffs$mid <- (c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4]))/2
diffs$mid <- diffs$edge - diffs$mid
diffs <- subset(diffs, span > 0)
diffs$bins <- cut(diffs$mid, c(0, unique(diffs$edge)), include.lowest = TRUE)
diffs.phan <- subset(diffs, Age < 541)

# Precambrian
diffs <-  c(coords.all$age,3640) # dummy base age added
diffs <- as.data.frame(diffs[order(diffs)])
names(diffs)[1] <- "Age"
diffs <- as.data.frame(diffs[!duplicated(diffs$Age),])
diffs <- round(diffs, digits = -1)
diffs_df <- round(diffs[-1,], digits = -1) - round(diffs[-nrow(diffs),], digits = -1)
diffs_df <- as.data.frame(diffs_df)
diffs <- cbind(diffs[1:nrow(diffs)-1,],diffs_df)
names(diffs)[1] <- "Age"
names(diffs)[2] <- "Diff"
diffs <- subset(diffs, Diff > 0)
diffs$Diff2 <- diffs$Diff/2
diffs$edge <- round(diffs$Diff2+diffs$Age)
diffs$span <- c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4])
diffs$mid <- (c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4]))/2
diffs$mid <- diffs$edge - diffs$mid
diffs <- subset(diffs, span > 0)
diffs$bins <- cut(diffs$mid, c(540, unique(diffs$edge)), include.lowest = TRUE)
diffs.pre <- subset(diffs, Age >= 540)

# Phanerozoic
phan_results.melt.bins <- phan_results.melt
phan_results.melt.bins$bins <- cut(phan_results.melt.bins$age, c(0, unique(diffs.phan$edge)), include.lowest = TRUE)

phan_results.melt.bins <- phan_results.melt.bins %>% 
  group_by(variable,bins) %>% 
  mutate(mean = mean(value),
         min = mean(value)-(sd(value)*2),
         max = mean(value)+(sd(value)*2))

phan_results.melt.bins[,6:8][phan_results.melt.bins[,6:8] < 0] <- 0
phan_results.melt.bins[,6:8][phan_results.melt.bins[,6:8] > 1] <- 1
phan_results.melt.bins <- phan_results.melt.bins[!duplicated(phan_results.melt.bins[c(3,5)]),]
phan_results.melt.bins <- phan_results.melt.bins[,-c(2,4)]
phan_results.melt.bins <- reshape2::melt(phan_results.melt.bins, id.vars = c("age","variable","bins"))
names(phan_results.melt.bins)[2] <- "type"
phan_results.bins <- reshape2::dcast(phan_results.melt.bins, age + bins + variable ~ type)
phan_results.bins <- join(diffs.phan, phan_results.bins, by = "bins")
phan_results.bins <- phan_results.bins[,c(4:6,9:ncol(phan_results.bins))]
phan_results.bins <- rbind(phan_results.bins[1:3,],phan_results.bins)
phan_results.bins[1:3,3] <- NA
phan_results.bins[1:3,1] <- 0
phan_results.melt.bins <- reshape2::melt(phan_results.bins, id.vars = c("edge","span", "mid","variable"))
names(phan_results.melt.bins)[5] <- "type"
write.csv(phan_results.bins, "phanerozoic_pyrite_types_downsampled.csv") # optional export

# Precambrian
pre_results.melt.bins <- pre_results.melt
pre_results.melt.bins$bins <- cut(pre_results.melt.bins$age, c(540, unique(diffs.pre$edge)), include.lowest = TRUE)

pre_results.melt.bins <- pre_results.melt.bins %>% 
  group_by(variable,bins) %>% 
  mutate(mean = mean(value),
         min = mean(value)-(sd(value)*2),
         max = mean(value)+(sd(value)*2))

pre_results.melt.bins[,6:8][pre_results.melt.bins[,6:8] < 0] <- 0
pre_results.melt.bins[,6:8][pre_results.melt.bins[,6:8] > 1] <- 1
pre_results.melt.bins <- pre_results.melt.bins[!duplicated(pre_results.melt.bins[c(3,5)]),]
pre_results.melt.bins <- pre_results.melt.bins[,-c(2,4)]
pre_results.melt.bins <- reshape2::melt(pre_results.melt.bins, id.vars = c("age","variable","bins"))
names(pre_results.melt.bins)[2] <- "type"
pre_results.bins <- reshape2::dcast(pre_results.melt.bins, age + bins + variable ~ type)
pre_results.bins <- join(diffs.pre, pre_results.bins, by = "bins")
pre_results.bins <- pre_results.bins[,c(4:6,9:ncol(pre_results.bins))]
pre_results.bins <- rbind(pre_results.bins[1,],pre_results.bins[1,],pre_results.bins)
pre_results.bins[1:3,4:9] <- pre_results.bins[4:6,4:9]
pre_results.melt.bins <- reshape2::melt(pre_results.bins, id.vars = c("edge","span", "mid", "variable"))
names(pre_results.melt.bins)[5] <- "type"

write.csv(pre_results.bins, "precambrian_pyrite_types_downsampled.csv") # optional export

# optional import results (to avoid re-running the above code every time)

#phan_results.bins <- read.csv("phanerozoic_pyrite_types_downsampled.csv")
#phan_results.bins <- phan_results.bins[,-1]
#names(phan_results.bins)[5:9] <- c("Type.5", "Type.1","Type.2", "Type.3", "Type.4")
#phan_results.melt.bins <- reshape2::melt(phan_results.bins, id.vars = c("edge","span", "mid", "variable"))
#names(phan_results.melt.bins)[5] <- "type"

#pre_results.bins <- read.csv("precambrian_pyrite_types_downsampled.csv")
#pre_results.bins <- pre_results.bins[,-1]
#names(pre_results.bins)[5:9] <- c("Type.5", "Type.1","Type.2", "Type.3", "Type.4")
#pre_results.melt.bins <- reshape2::melt(pre_results.bins, id.vars = c("edge","span", "mid", "variable"))
#names(pre_results.melt.bins)[5] <- "type"

HEATT <- c(56, 66, 93, 116, 183, 200, 251, 359, 372, 383, 444, 514, 542) 

xdd <- read.csv("xdd_binned_results.csv") # import text mining outputs

# print

phan_results.melt.bins$type <- factor(phan_results.melt.bins$type, levels = c("Type.5", "Type.2", "Type.1", "Type.3", "Type.4"))
pre_results.melt.bins$type <- factor(pre_results.melt.bins$type, levels = c("Type.5", "Type.2", "Type.1", "Type.3", "Type.4"))

phan_results.melt.bins.mean <- subset(phan_results.melt.bins, variable == "mean")
pre_results.melt.bins.mean <- subset(pre_results.melt.bins, variable == "mean")

C <- ggplot() +
  geom_bar(data = phan_results.melt.bins.mean, aes(y = value, x = mid, fill = type, width = span),
           stat="identity") +
  theme_bw() + scale_x_reverse(limits = c(550,0)) +
  geom_point(aes(x = unique(coords.all$age), y = 1.02)) +
  geom_step(data = phan_results.melt.bins.mean, aes(y = value, x = edge, group = type), 
            position = "stack", direction = "hv") +
  geom_vline(xintercept = HEATT, colour = "red") + 
  theme(legend.position = "none",
        axis.title.x = element_blank()) + 
  geom_vline(xintercept = 540) +
  coord_geo(dat = list("periods", "eras"), xlim = c(550, 0), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c(""))

D <- ggplot() +
  geom_bar(data = pre_results.melt.bins.mean, aes(y = value, x = mid, fill = type, width = span),
           stat="identity") +
  theme_bw() + scale_x_reverse(limits = c(3600,540)) +
  geom_point(aes(x = unique(coords.all$age), y = 1.02)) +
  geom_step(data = pre_results.melt.bins.mean, aes(y = value, x = edge, group = type), 
            position = "stack", direction = "hv") +
  geom_vline(xintercept = HEATT, colour = "red") + 
  theme(legend.position = "none",
        axis.title.x = element_blank()) + 
  coord_geo(dat = list("periods", "eras"), xlim = c(3600, 540), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c(""))

grid.arrange(D,C,ncol=2)

### scaled to xdd output ##

phan_results.bins.mean <- subset(phan_results.bins, variable == "mean")
pre_results.bins.mean <- subset(pre_results.bins, variable == "mean")

xdd_phan <- subset(xdd, Age < 540)
xdd_phan.dt <- data.table(xdd_phan , key = "Age" )
phan_bins.dt <- data.table(phan_results.bins.mean, key = "mid")
phan_out <- as.data.frame(phan_bins.dt[xdd_phan.dt, list(Age, framboids, nodules, undif,Type.5,Type.2,Type.1,Type.3,Type.4) , roll = "nearest" ]) # DType.3,PType.3,DECOUPLED-EUX
phan_out$ratio <- (phan_out$framboids+phan_out$nodules)/(phan_out$framboids+phan_out$undif+phan_out$nodules)
phan_out <- phan_out[,-c(2:4)]
phan_out <- reshape2::melt(phan_out, id.vars = c("Age","ratio"))
phan_out$value <- phan_out$ratio*phan_out$value # scaled output (to Fig. 1B)

xdd_pre <- subset(xdd, Age >= 540 & Age < 3600)
xdd_pre.dt <- data.table(xdd_pre , key = "Age" )
pre_bins.dt <- data.table(pre_results.bins.mean, key = "mid")
pre_out <- as.data.frame(pre_bins.dt[xdd_pre.dt, list(Age, framboids, nodules, undif,Type.5,Type.2,Type.1,Type.3,Type.4) , roll = "nearest" ]) #DType.3,PType.3,DECOUPLED-EUX
pre_out$ratio <- (pre_out$framboids+pre_out$nodules)/(pre_out$framboids+pre_out$undif+pre_out$nodules)
pre_out <- pre_out[,-c(2:4)]
pre_out <- reshape2::melt(pre_out, id.vars = c("Age","ratio"))
pre_out$value <- pre_out$ratio*pre_out$value # scaled output (to Fig. 1B)

scaleFUN <- function(x) sprintf("%.2f", x)

A <- ggplot() +
  geom_bar(data = phan_out, aes(y = value, x = Age, fill = variable),
           stat="identity") +
  theme_bw() + scale_x_reverse(limits = c(550,0)) +
  geom_point(aes(x = unique(coords.all$age), y = -0.01)) +
  geom_step(data = phan_out, aes(y = value, x = Age, group = variable), 
            position = "stack", direction = "hv") +
  geom_vline(xintercept = HEATT, colour = "red") +
  theme(legend.position = "none",
        axis.title.x = element_blank()) + 
  scale_y_continuous(limits = c(0,0.4), labels = scaleFUN) +
  coord_geo(dat = list("periods", "eras"), xlim = c(550, 0), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("t","t"), abbrv = list(TRUE, TRUE),
            skip = c(""))

B <- ggplot() +
  geom_bar(data = pre_out, aes(y = value, x = Age, fill = variable),
           stat="identity", width = 10) +
  theme_bw() + scale_x_reverse(limits = c(3600,540)) +
  geom_point(aes(x = unique(coords.all$age), y = -0.01)) +
  geom_step(data = pre_out, aes(y = value, x = Age, group = variable), 
            position = "stack", direction = "mid") +
  geom_vline(xintercept = HEATT, colour = "red") +
  theme(legend.position = "none",
        axis.title.x = element_blank()) + 
  scale_y_continuous(limits = c(0,0.4), labels = scaleFUN) +
  coord_geo(dat = list("periods", "eras"), xlim = c(3600, 540), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("t","t"), abbrv = list(TRUE, TRUE),
            skip = c(""))

# Fig 4A-B

grid.arrange(B,A,D,C,ncol=2)

# uncertainty estimates per pyrite type

phan_results.melt.bins.error <- subset(phan_results.melt.bins, variable == "min" | variable == "max")
phan_results.bins.error <- reshape2::dcast(phan_results.melt.bins.error, edge+span+mid+type~variable)

pre_results.melt.bins.error <- subset(pre_results.melt.bins, variable == "min" | variable == "max")
pre_results.bins.error <- reshape2::dcast(pre_results.melt.bins.error, edge+span+mid+type~variable)

# Type.5

a <- ggplot() +
  geom_stepribbon(data = subset(phan_results.bins.error, type == "Type.5"),
                  aes(x = edge, ymin = min, ymax = max), fill = "grey50") +
  geom_step(data = phan_results.bins.mean, aes(y = Type.5, x = edge),
            direction = "hv") +
  theme_bw() + scale_x_reverse(limits = c(550,0)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Age Ma") +
  coord_geo(dat = list("periods", "eras"), xlim = c(540, 0), size = 3, 
                            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
                            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
                            skip = c("Msr"))

b <- ggplot() +
  geom_stepribbon(data = subset(pre_results.bins.error, type == "Type.5"),
                  aes(x = edge, ymin = min, ymax = max), fill = "grey50") +
  geom_step(data = pre_results.bins.mean, aes(y = Type.5, x = edge),
            direction = "hv") +
  theme_bw() + scale_x_reverse(limits = c(4000, 540)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Age Ma") +
  coord_geo(dat = list("periods", "eras"), xlim = c(4000, 540), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr"))

grid.arrange(b,a,ncol=2)

# Type.2

a <- ggplot() +
  geom_stepribbon(data = subset(phan_results.bins.error, type == "Type.2"),
                  aes(x = edge, ymin = min, ymax = max), fill = "grey50") +
  geom_step(data = phan_results.bins.mean, aes(y = Type.2, x = edge),
            direction = "hv") +
  theme_bw() + scale_x_reverse(limits = c(550,0)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Age Ma") +
  coord_geo(dat = list("periods", "eras"), xlim = c(550, 0), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr"))

b <- ggplot() +
  geom_stepribbon(data = subset(pre_results.bins.error, type == "Type.2"),
                  aes(x = edge, ymin = min, ymax = max), fill = "grey50") +
  geom_step(data = pre_results.bins.mean, aes(y = Type.2, x = edge),
            direction = "hv") +
  theme_bw() + scale_x_reverse(limits = c(4000, 540)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Age Ma") +
  coord_geo(dat = list("periods", "eras"), xlim = c(4000, 540), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr"))

grid.arrange(b,a,ncol=2)

# Type.1

a <- ggplot() +
  geom_stepribbon(data = subset(phan_results.bins.error, type == "Type.1"),
                  aes(x = edge, ymin = min, ymax = max), fill = "grey50") +
  geom_step(data = phan_results.bins.mean, aes(y = Type.1, x = edge),
            direction = "hv") +
  theme_bw() + scale_x_reverse(limits = c(550,0)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Age Ma") +
  coord_geo(dat = list("periods", "eras"), xlim = c(540, 0), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr"))


b <- ggplot() +
  geom_stepribbon(data = subset(pre_results.bins.error, type == "Type.1"),
                  aes(x = edge, ymin = min, ymax = max), fill = "grey50") +
  geom_step(data = pre_results.bins.mean, aes(y = Type.1, x = edge),
            direction = "hv") +
  theme_bw() + scale_x_reverse(limits = c(4000, 540)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Age Ma") +
  coord_geo(dat = list("periods", "eras"), xlim = c(4000, 540), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr"))


grid.arrange(b,a,ncol=2)

# Type.3

a <- ggplot() +
  geom_stepribbon(data = subset(phan_results.bins.error, type == "Type.3"),
                  aes(x = edge, ymin = min, ymax = max), fill = "grey50") +
  geom_step(data = phan_results.bins.mean, aes(y = Type.3, x = edge),
            direction = "hv") +
  theme_bw() + scale_x_reverse(limits = c(550,0)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Age Ma") +
  coord_geo(dat = list("periods", "eras"), xlim = c(550, 0), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr"))


b <- ggplot() +
  geom_stepribbon(data = subset(pre_results.bins.error, type == "Type.3"),
                  aes(x = edge, ymin = min, ymax = max), fill = "grey50") +
  geom_step(data = pre_results.bins.mean, aes(y = Type.3, x = edge),
            direction = "hv") +
  theme_bw() + scale_x_reverse(limits = c(4000, 540)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Age Ma") +
  coord_geo(dat = list("periods", "eras"), xlim = c(4000, 540), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr"))


grid.arrange(b,a,ncol=2)

# Type.4

a <- ggplot() +
  geom_stepribbon(data = subset(phan_results.bins.error, type == "Type.4"),
                  aes(x = edge, ymin = min, ymax = max), fill = "grey50") +
  geom_step(data = phan_results.bins.mean, aes(y = Type.4, x = edge),
            direction = "hv") +
  theme_bw() + scale_x_reverse(limits = c(550,0)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Age Ma") +
  coord_geo(dat = list("periods", "eras"), xlim = c(550, 0), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr"))

b <- ggplot() +
  geom_stepribbon(data = subset(pre_results.bins.error, type == "Type.4"),
                  aes(x = edge, ymin = min, ymax = max), fill = "grey50") +
  geom_step(data = pre_results.bins.mean, aes(y = Type.4, x = edge),
            direction = "hv") +
  theme_bw() + scale_x_reverse(limits = c(4000, 540)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Age Ma") +
  coord_geo(dat = list("periods", "eras"), xlim = c(4000, 540), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr"))

grid.arrange(b,a,ncol=2)

# binned ratios with uncertainty estimates

# Type.5 error
phan_results.melt.bins.error.POM <- subset(phan_results.melt.bins, type == "Type.5")
phan_results.melt.bins.error.POM <- reshape2::dcast(phan_results.melt.bins.error.POM,
                                                    edge+span+mid+type~variable)
phan_results.melt.bins.error.POM$error <- (phan_results.melt.bins.error.POM$max-phan_results.melt.bins.error.POM$mean)
pre_results.melt.bins.error.POM <- subset(pre_results.melt.bins, type == "Type.5")
pre_results.melt.bins.error.POM <- reshape2::dcast(pre_results.melt.bins.error.POM,
                                                   edge+span+mid+type~variable)
pre_results.melt.bins.error.POM$error <- (pre_results.melt.bins.error.POM$max-pre_results.melt.bins.error.POM$mean)

# Type.2 error
phan_results.melt.bins.error.DOM.OX <- subset(phan_results.melt.bins, type == "Type.2")
phan_results.melt.bins.error.DOM.OX <- reshape2::dcast(phan_results.melt.bins.error.DOM.OX,
                                                       edge+span+mid+type~variable)
phan_results.melt.bins.error.DOM.OX$error <- (phan_results.melt.bins.error.DOM.OX$max-phan_results.melt.bins.error.DOM.OX$mean)
pre_results.melt.bins.error.DOM.OX <- subset(pre_results.melt.bins, type == "Type.2")
pre_results.melt.bins.error.DOM.OX <- reshape2::dcast(pre_results.melt.bins.error.DOM.OX,
                                                      edge+span+mid+type~variable)
pre_results.melt.bins.error.DOM.OX$error <- (pre_results.melt.bins.error.DOM.OX$max-pre_results.melt.bins.error.DOM.OX$mean)

# Type.1 error
phan_results.melt.bins.error.DOM.DFe <- subset(phan_results.melt.bins, type == "Type.1")
phan_results.melt.bins.error.DOM.DFe <- reshape2::dcast(phan_results.melt.bins.error.DOM.DFe,
                                                        edge+span+mid+type~variable)
phan_results.melt.bins.error.DOM.DFe$error <- (phan_results.melt.bins.error.DOM.DFe$max-phan_results.melt.bins.error.DOM.DFe$mean)
pre_results.melt.bins.error.DOM.DFe <- subset(pre_results.melt.bins, type == "Type.1")
pre_results.melt.bins.error.DOM.DFe <- reshape2::dcast(pre_results.melt.bins.error.DOM.DFe,
                                                       edge+span+mid+type~variable)
pre_results.melt.bins.error.DOM.DFe$error <- (pre_results.melt.bins.error.DOM.DFe$max-pre_results.melt.bins.error.DOM.DFe$mean)

# Type.3 error
phan_results.melt.bins.error.OM.EUX <- subset(phan_results.melt.bins, type == "Type.3")
phan_results.melt.bins.error.OM.EUX <- reshape2::dcast(phan_results.melt.bins.error.OM.EUX,
                                                       edge+span+mid+type~variable)
phan_results.melt.bins.error.OM.EUX$error <- (phan_results.melt.bins.error.OM.EUX$max-phan_results.melt.bins.error.OM.EUX$mean)
pre_results.melt.bins.error.OM.EUX <- subset(pre_results.melt.bins, type == "Type.3")
pre_results.melt.bins.error.OM.EUX <- reshape2::dcast(pre_results.melt.bins.error.OM.EUX,
                                                      edge+span+mid+type~variable)
pre_results.melt.bins.error.OM.EUX$error <- (pre_results.melt.bins.error.OM.EUX$max-pre_results.melt.bins.error.OM.EUX$mean)

# Type.4 error
phan_results.melt.bins.error.OX.EUX <- subset(phan_results.melt.bins, type == "Type.4")
phan_results.melt.bins.error.OX.EUX <- reshape2::dcast(phan_results.melt.bins.error.OX.EUX,
                                                       edge+span+mid+type~variable)
phan_results.melt.bins.error.OX.EUX$error <- (phan_results.melt.bins.error.OX.EUX$max-phan_results.melt.bins.error.OX.EUX$mean)
pre_results.melt.bins.error.OX.EUX <- subset(pre_results.melt.bins, type == "Type.4")
pre_results.melt.bins.error.OX.EUX <- reshape2::dcast(pre_results.melt.bins.error.OX.EUX,
                                                      edge+span+mid+type~variable)
pre_results.melt.bins.error.OX.EUX$error <- (pre_results.melt.bins.error.OX.EUX$max-pre_results.melt.bins.error.OX.EUX$mean)

# POM vs DOM - Fig. 4C

# error in quadrature
phan.POM.DOM <- phan_results.melt.bins.error.POM[,1:3]
phan.POM.DOM$POM <- phan_results.melt.bins.error.POM$mean
phan.POM.DOM$DOM <- phan_results.melt.bins.error.DOM.OX$mean+phan_results.melt.bins.error.DOM.DFe$mean
phan.POM.DOM$ratio <- phan.POM.DOM$POM/(phan.POM.DOM$POM+phan.POM.DOM$DOM)
phan.POM.DOM$den <- ((phan_results.melt.bins.error.POM$error)^2+(phan_results.melt.bins.error.DOM.OX$error)^2+(phan_results.melt.bins.error.DOM.DFe$error)^2)^0.5
phan.POM.DOM$num <- (phan_results.melt.bins.error.POM$error)
phan.POM.DOM$quad <- ((phan.POM.DOM$num/phan.POM.DOM$POM)^2+(phan.POM.DOM$den/(phan.POM.DOM$DOM+phan.POM.DOM$POM))^2)^0.5
phan.POM.DOM$quad <- (phan.POM.DOM$ratio)*(phan.POM.DOM$quad)
phan.POM.DOM$upper <- phan.POM.DOM$ratio+phan.POM.DOM$quad
phan.POM.DOM$lower <- phan.POM.DOM$ratio-phan.POM.DOM$quad
phan.POM.DOM[,10:11][phan.POM.DOM[,10:11] < 0] <- 0
phan.POM.DOM[,10:11][phan.POM.DOM[,10:11] > 1] <- 1

#error in quad
pre.POM.DOM <- pre_results.melt.bins.error.POM[,1:3]
pre.POM.DOM$POM <- pre_results.melt.bins.error.POM$mean
pre.POM.DOM$DOM <- pre_results.melt.bins.error.DOM.OX$mean+pre_results.melt.bins.error.DOM.DFe$mean
pre.POM.DOM$ratio <- pre.POM.DOM$POM/(pre.POM.DOM$POM+pre.POM.DOM$DOM)
pre.POM.DOM$den <- ((pre_results.melt.bins.error.POM$error)^2+(pre_results.melt.bins.error.DOM.OX$error)^2+(pre_results.melt.bins.error.DOM.DFe$error)^2)^0.5
pre.POM.DOM$num <- (pre_results.melt.bins.error.POM$error)
pre.POM.DOM$quad <- ((pre.POM.DOM$num/pre.POM.DOM$POM)^2+(pre.POM.DOM$den/(pre.POM.DOM$DOM+pre.POM.DOM$POM))^2)^0.5
pre.POM.DOM$quad <- (pre.POM.DOM$ratio)*(pre.POM.DOM$quad)
pre.POM.DOM$upper <- pre.POM.DOM$ratio+pre.POM.DOM$quad
pre.POM.DOM$lower <- pre.POM.DOM$ratio-pre.POM.DOM$quad
pre.POM.DOM[,10:11][pre.POM.DOM[,10:11] < 0] <- 0
pre.POM.DOM[,10:11][pre.POM.DOM[,10:11] > 1] <- 1

# POM vs. DOM

a1 <- ggplot() +
  geom_stepribbon(data = phan.POM.DOM,
                  aes(x = edge, ymin = lower, ymax = upper), fill = "grey50") +
  geom_step(data = phan.POM.DOM, aes(y =ratio, x = edge),
            direction = "hv", colour = "blue") +
  theme_bw() + 
  scale_x_reverse(limits = c(550,0), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1)) 

b1 <- ggplot() +
  geom_stepribbon(data = pre.POM.DOM,
                  aes(x = edge, ymin = lower, ymax = upper), fill = "grey50") +
  geom_step(data = pre.POM.DOM, aes(y =ratio, x = edge),
            direction = "hv", colour = "blue") +
  theme_bw() + 
  scale_x_reverse(limits = c(4000, 540), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1))

grid.arrange(b1,a1,ncol=2)

# Microenvironment vs. ambient - Fig. 4D

# error in quadrature
phan.micro.amb <- phan_results.melt.bins.error.POM[,1:3]
phan.micro.amb$micro <- phan_results.melt.bins.error.POM$mean+phan_results.melt.bins.error.DOM.OX$mean+phan_results.melt.bins.error.DOM.DFe$mean
phan.micro.amb$amb <- phan_results.melt.bins.error.OX.EUX$mean+phan_results.melt.bins.error.OM.EUX$mean
phan.micro.amb$ratio <- phan.micro.amb$micro/(phan.micro.amb$micro+phan.micro.amb$amb)
phan.micro.amb$den <- ((phan_results.melt.bins.error.OX.EUX$error)^2+(phan_results.melt.bins.error.OM.EUX$error)^2+(phan_results.melt.bins.error.POM$error)^2+(phan_results.melt.bins.error.DOM.OX$error)^2+(phan_results.melt.bins.error.DOM.DFe$error)^2)^0.5
phan.micro.amb$num <- ((phan_results.melt.bins.error.POM$error)^2+(phan_results.melt.bins.error.DOM.OX$error)^2+(phan_results.melt.bins.error.DOM.DFe$error)^2)^0.5
phan.micro.amb$quad <- ((phan.micro.amb$num/phan.micro.amb$micro)^2+(phan.micro.amb$den/(phan.micro.amb$amb+phan.micro.amb$micro))^2)^0.5
phan.micro.amb$quad <- (phan.micro.amb$ratio)*(phan.micro.amb$quad)
phan.micro.amb$upper <- phan.micro.amb$ratio+phan.micro.amb$quad
phan.micro.amb$lower <- phan.micro.amb$ratio-phan.micro.amb$quad
phan.micro.amb[,10:11][phan.micro.amb[,10:11] < 0] <- 0
phan.micro.amb[,10:11][phan.micro.amb[,10:11] > 1] <- 1

#error in quad
pre.micro.amb <- pre_results.melt.bins.error.POM[,1:3]
pre.micro.amb$micro <- pre_results.melt.bins.error.POM$mean+pre_results.melt.bins.error.DOM.OX$mean+pre_results.melt.bins.error.DOM.DFe$mean
pre.micro.amb$amb <- pre_results.melt.bins.error.OX.EUX$mean+pre_results.melt.bins.error.OM.EUX$mean
pre.micro.amb$ratio <- pre.micro.amb$micro/(pre.micro.amb$micro+pre.micro.amb$amb)
pre.micro.amb$den <- ((pre_results.melt.bins.error.OX.EUX$error)^2+(pre_results.melt.bins.error.OM.EUX$error)^2+(pre_results.melt.bins.error.POM$error)^2+(pre_results.melt.bins.error.DOM.OX$error)^2+(pre_results.melt.bins.error.DOM.DFe$error)^2)^0.5
pre.micro.amb$num <- ((pre_results.melt.bins.error.POM$error)^2+(pre_results.melt.bins.error.DOM.OX$error)^2+(pre_results.melt.bins.error.DOM.DFe$error)^2)^0.5
pre.micro.amb$quad <- ((pre.micro.amb$num/pre.micro.amb$micro)^2+(pre.micro.amb$den/(pre.micro.amb$amb+pre.micro.amb$micro))^2)^0.5
pre.micro.amb$quad <- (pre.micro.amb$ratio)*(pre.micro.amb$quad)
pre.micro.amb$upper <- pre.micro.amb$ratio+pre.micro.amb$quad
pre.micro.amb$lower <- pre.micro.amb$ratio-pre.micro.amb$quad
pre.micro.amb[,10:11][pre.micro.amb[,10:11] < 0] <- 0
pre.micro.amb[,10:11][pre.micro.amb[,10:11] > 1] <- 1

a2 <- ggplot() +
  geom_stepribbon(data = phan.micro.amb,
                  aes(x = edge, ymin = lower, ymax = upper), fill = "grey50") +
  geom_step(data = phan.micro.amb, aes(y =ratio, x = edge),
            direction = "hv", colour = "blue") +
  theme_bw() +
  scale_x_reverse(limits = c(550,0), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1))

b2 <- ggplot() +
  geom_stepribbon(data = pre.micro.amb,
                  aes(x = edge, ymin = lower, ymax = upper), fill = "grey50") +
  geom_step(data = pre.micro.amb, aes(y =ratio, x = edge),
            direction = "hv", colour = "blue") +
  theme_bw() +
  scale_x_reverse(limits = c(4000,540), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1))

grid.arrange(b2,a2,ncol=2)

# Fig. 4E

# Fe-OX (Type 5 + Type 2) vs. DFe (Type 1) microenvironments

# error in quadrature
phan.OX.DFe <- phan_results.melt.bins.error.POM[,1:3]
phan.OX.DFe$OX <- phan_results.melt.bins.error.DOM.OX$mean+phan_results.melt.bins.error.POM$mean
phan.OX.DFe$DFe <- phan_results.melt.bins.error.DOM.DFe$mean
phan.OX.DFe$ratio <- phan.OX.DFe$OX/(phan.OX.DFe$DFe+phan.OX.DFe$OX)
phan.OX.DFe$den <- ((phan_results.melt.bins.error.DOM.OX$error)^2+(phan_results.melt.bins.error.DOM.DFe$error)^2+(phan_results.melt.bins.error.POM$error)^2)^0.5
phan.OX.DFe$num <- ((phan_results.melt.bins.error.DOM.OX$error)^2+(phan_results.melt.bins.error.POM$error)^2)^0.5
phan.OX.DFe$quad <- ((phan.OX.DFe$num/phan.OX.DFe$OX)^2+(phan.OX.DFe$den/(phan.OX.DFe$DFe+phan.OX.DFe$OX))^2)^0.5
phan.OX.DFe$quad <- (phan.OX.DFe$ratio)*(phan.OX.DFe$quad)
phan.OX.DFe$upper <- phan.OX.DFe$ratio+phan.OX.DFe$quad
phan.OX.DFe$lower <- phan.OX.DFe$ratio-phan.OX.DFe$quad
phan.OX.DFe[,10:11][phan.OX.DFe[,10:11] < 0] <- 0
phan.OX.DFe[,10:11][phan.OX.DFe[,10:11] > 1] <- 1

pre.OX.DFe <- pre_results.melt.bins.error.POM[,1:3]
pre.OX.DFe$OX <- pre_results.melt.bins.error.DOM.OX$mean+pre_results.melt.bins.error.POM$mean
pre.OX.DFe$DFe <- pre_results.melt.bins.error.DOM.DFe$mean
pre.OX.DFe$ratio <- pre.OX.DFe$OX/(pre.OX.DFe$DFe+pre.OX.DFe$OX)
pre.OX.DFe$den <- ((pre_results.melt.bins.error.DOM.OX$error)^2+(pre_results.melt.bins.error.DOM.DFe$error)^2+(pre_results.melt.bins.error.POM$error)^2)^0.5
pre.OX.DFe$num <- ((pre_results.melt.bins.error.DOM.OX$error)^2+(pre_results.melt.bins.error.POM$error)^2)^0.5
pre.OX.DFe$quad <- ((pre.OX.DFe$num/pre.OX.DFe$OX)^2+(pre.OX.DFe$den/(pre.OX.DFe$DFe+pre.OX.DFe$OX))^2)^0.5
pre.OX.DFe$quad <- (pre.OX.DFe$ratio)*(pre.OX.DFe$quad)
pre.OX.DFe$upper <- pre.OX.DFe$ratio+pre.OX.DFe$quad
pre.OX.DFe$lower <- pre.OX.DFe$ratio-pre.OX.DFe$quad
pre.OX.DFe[,10:11][pre.OX.DFe[,10:11] < 0] <- 0
pre.OX.DFe[,10:11][pre.OX.DFe[,10:11] > 1] <- 1

# Fe-OX (Type 5 + Type 2) vs. DFe (Type 1)

a3 <- ggplot() +
  geom_stepribbon(data = phan.OX.DFe,
                  aes(x = edge, ymin = lower, ymax = upper), fill = "grey50") +
  geom_step(data = phan.OX.DFe, aes(y =ratio, x = edge),
            direction = "hv", colour = "blue") +
  theme_bw() + 
  scale_x_reverse(limits = c(550,0), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1)) 

b3 <- ggplot() +
  geom_stepribbon(data = pre.OX.DFe,
                  aes(x = edge, ymin = lower, ymax = upper), fill = "grey50") +
  geom_step(data = pre.OX.DFe, aes(y =ratio, x = edge),
            direction = "hv", colour = "blue") +
  theme_bw() + 
  scale_x_reverse(limits = c(4000,540), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1)) 


grid.arrange(b3,a3,ncol=2)

# Figure 4 - tidy in professional graphics package (e.g., Adobe Illustrator)

grid.arrange(B,A,D,C,b1,a1,b2,a2,b3,a3,ncol=2)

###### SGP Phase 1 analysis - used in Figs. 5-6 ####

##### import SGP data 

# call SGP API - not culled to focal area # 

body = '{"type":"samples","filters":{"lithology_class":["sedimentary","metamorphic"]},"show":["fe","mn","fe_carb","fe_ox","fe_mag","fe_py","fe_py_fe_hr","fe_hr_fe_t","toc","mo","u","p","ni","v","interpreted_age","fe_hr","analysis_ref_long", "coord_lat","coord_long"]}'

r <- POST("https://sgp-search.io/api/v1/post", body=body, 
          httr::add_headers(`accept` = 'application/json'), 
          httr::content_type('application/json'))

out <- content(r, "text")

SGP <- jsonlite::fromJSON(out)

# write SGP refs
SGP_refs <- subset(SGP, select = `analysis reference (long)`)
SGP_refs <- as.data.frame(SGP_refs[!duplicated(SGP_refs),])
SGP_refs <- as.data.frame(SGP_refs[complete.cases(SGP_refs),])
names(SGP_refs)[1] <- "SGP reference"
SGP_refs <- as.data.frame(SGP_refs[order(SGP_refs$`SGP reference`),])
#write.csv(SGP_refs, "SGP_refs.csv") # export SGP reference list

SGP[,3:20] <- sapply(SGP[,3:20], as.numeric)

ppm <- dplyr::select(SGP,contains(" (ppm)")) 
ppb <- dplyr::select(SGP,contains(" (ppb)"))
ppb <- ppb/1000
wt <- dplyr::select(SGP,contains(" (wt%)"))
wt <- wt*10000

SGP2 <- cbind(SGP[,c(1:5,21)],ppm,ppb,wt)

names(SGP2) <- gsub(x = names(SGP2), pattern = " (ppm)", replacement = "", fixed = TRUE)
names(SGP2) <- gsub(x = names(SGP2), pattern = " (ppb)", replacement = "", fixed = TRUE)
names(SGP2) <- gsub(x = names(SGP2), pattern = " (wt%)", replacement = "", fixed = TRUE)
names(SGP2) <- gsub(x = names(SGP2), pattern = "-", replacement = "", fixed = TRUE)
names(SGP2)[5] <- "age"

# define intervals based on pyrite type distributions #

zones <- c(3400,3100,3000,2800,2100,1700,1500,1400,900,800,740,690,655,610,532,514,480,456,419,383,372,367,359,330,299,280,250,201,54,0)

SGP2$zones <- cut(SGP2$age, zones)

get_midpoint <- function(cut_label) {
  mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ","))))
}

SGP2$zones2 <- sapply(SGP2$zones, get_midpoint)

###### SGP - Fepy/FeHR ####

# optional anoxic subset

SGP2 <- subset(SGP2, (Fecarb+Feox+Femag+Fepy)/Fe > 0.38 & Fe > (0.5*10000))

SGP2$target1 <- SGP2$Fepy/(SGP2$Fecarb+SGP2$Feox+SGP2$Femag+SGP2$Fepy)

SGP2_phan <- subset(SGP2, age <= 540)
SGP2_phan1 <- subset(SGP2_phan, age > 150)
SGP2_phan2 <- subset(SGP2_phan, age < 200)
SGP2_pre <- subset(SGP2, age > 540)

# build normal distribution

speed <- 10 # speed = 10 for full version, but can increase as above.

# construct normal distribution for sampling

n <- 10000
norm <- rnorm(n)
out <- hist(norm, breaks = c(seq(-4.3,4.3, by = 0.2)))
bins <- as.data.frame(out[[1]])
bins <- round(bins+0.1, digits = 2)
freqs <- as.data.frame(out[[2]])
out <- as.data.frame(cbind(bins[1:nrow(bins)-1,],freqs))
names(out)[1:2] <- c("sd","count")

# phanerozoic

# span_upper <- loess.predict[["pars"]][["span"]]
# span_upper <- loess.predict[["pars"]][["span"]]

# 2sd

nsig <- 2

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_2sd_phan1_upper <- as.data.frame(cbind(target1_2sd_phan1$x, target1_2sd_phan1$upper))
names(target1_2sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_2sd_phan1_upper$age, target1_2sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_2sd_phan1_upper, indices) {
  d <- target1_2sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_2sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_2sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_2sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_2sd_phan1_upper_50)[1] <- "target1"

target1_2sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_2sd_phan1_upper_all)[1:ncol(target1_2sd_phan1_upper_all)] <- seq(0,540, 1)
target1_2sd_phan1_upper_all$group <- seq(1,(nrow(target1_2sd_phan1_upper_all)), 1)

target1_2sd_phan1_upper_all <- reshape2::melt(target1_2sd_phan1_upper_all, id.vars = "group")
names(target1_2sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_2sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_2sd_phan1_lower <- as.data.frame(cbind(target1_2sd_phan1$x, target1_2sd_phan1$lower))
names(target1_2sd_phan1_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_2sd_phan1_lower$age, target1_2sd_phan1_lower$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_2sd_phan1_lower, indices) {
  d <- target1_2sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_2sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_2sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_2sd_phan1_lower_50$age <- seq(0,540, 1)
names(target1_2sd_phan1_lower_50)[1] <- "target1"

target1_2sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target1_2sd_phan1_lower_all)[1:ncol(target1_2sd_phan1_lower_all)] <- seq(0,540, 1)
target1_2sd_phan1_lower_all$group <- seq(1,(nrow(target1_2sd_phan1_lower_all)), 1)

target1_2sd_phan1_lower_all <- reshape2::melt(target1_2sd_phan1_lower_all, id.vars = "group")
names(target1_2sd_phan1_lower_all)[2:3] <- c("age", "target1")

target1_2sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8sd

nsig <- 1.8

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.8sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_1.8sd_phan1_upper <- as.data.frame(cbind(target1_1.8sd_phan1$x, target1_1.8sd_phan1$upper))
names(target1_1.8sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.8sd_phan1_upper$age, target1_1.8sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_1.8sd_phan1_upper, indices) {
  d <- target1_1.8sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.8sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.8sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.8sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_1.8sd_phan1_upper_50)[1] <- "target1"

target1_1.8sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.8sd_phan1_upper_all)[1:ncol(target1_1.8sd_phan1_upper_all)] <- seq(0,540, 1)
target1_1.8sd_phan1_upper_all$group <- seq(1,(nrow(target1_1.8sd_phan1_upper_all)), 1)

target1_1.8sd_phan1_upper_all <- reshape2::melt(target1_1.8sd_phan1_upper_all, id.vars = "group")
names(target1_1.8sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_1.8sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.8sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_1.8sd_phan1_lower <- as.data.frame(cbind(target1_1.8sd_phan1$x, target1_1.8sd_phan1$lower))
names(target1_1.8sd_phan1_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.8sd_phan1_lower$age, target1_1.8sd_phan1_lower$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_1.8sd_phan1_lower, indices) {
  d <- target1_1.8sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.8sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.8sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.8sd_phan1_lower_50$age <- seq(0,540, 1)
names(target1_1.8sd_phan1_lower_50)[1] <- "target1"

target1_1.8sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.8sd_phan1_lower_all)[1:ncol(target1_1.8sd_phan1_lower_all)] <- seq(0,540, 1)
target1_1.8sd_phan1_lower_all$group <- seq(1,(nrow(target1_1.8sd_phan1_lower_all)), 1)

target1_1.8sd_phan1_lower_all <- reshape2::melt(target1_1.8sd_phan1_lower_all, id.vars = "group")
names(target1_1.8sd_phan1_lower_all)[2:3] <- c("age", "target1")

target1_1.8sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6sd

nsig <- 1.6

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.6sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_1.6sd_phan1_upper <- as.data.frame(cbind(target1_1.6sd_phan1$x, target1_1.6sd_phan1$upper))
names(target1_1.6sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.6sd_phan1_upper$age, target1_1.6sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_1.6sd_phan1_upper, indices) {
  d <- target1_1.6sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.6sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.6sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.6sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_1.6sd_phan1_upper_50)[1] <- "target1"

target1_1.6sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.6sd_phan1_upper_all)[1:ncol(target1_1.6sd_phan1_upper_all)] <- seq(0,540, 1)
target1_1.6sd_phan1_upper_all$group <- seq(1,(nrow(target1_1.6sd_phan1_upper_all)), 1)

target1_1.6sd_phan1_upper_all <- reshape2::melt(target1_1.6sd_phan1_upper_all, id.vars = "group")
names(target1_1.6sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_1.6sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.6sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_1.6sd_phan1_lower <- as.data.frame(cbind(target1_1.6sd_phan1$x, target1_1.6sd_phan1$lower))
names(target1_1.6sd_phan1_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.6sd_phan1_lower$age, target1_1.6sd_phan1_lower$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_1.6sd_phan1_lower, indices) {
  d <- target1_1.6sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.6sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.6sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.6sd_phan1_lower_50$age <- seq(0,540, 1)
names(target1_1.6sd_phan1_lower_50)[1] <- "target1"

target1_1.6sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.6sd_phan1_lower_all)[1:ncol(target1_1.6sd_phan1_lower_all)] <- seq(0,540, 1)
target1_1.6sd_phan1_lower_all$group <- seq(1,(nrow(target1_1.6sd_phan1_lower_all)), 1)

target1_1.6sd_phan1_lower_all <- reshape2::melt(target1_1.6sd_phan1_lower_all, id.vars = "group")
names(target1_1.6sd_phan1_lower_all)[2:3] <- c("age", "target1")

target1_1.6sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4sd

nsig <- 1.4

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.4sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_1.4sd_phan1_upper <- as.data.frame(cbind(target1_1.4sd_phan1$x, target1_1.4sd_phan1$upper))
names(target1_1.4sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.4sd_phan1_upper$age, target1_1.4sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_1.4sd_phan1_upper, indices) {
  d <- target1_1.4sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.4sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.4sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.4sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_1.4sd_phan1_upper_50)[1] <- "target1"

target1_1.4sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.4sd_phan1_upper_all)[1:ncol(target1_1.4sd_phan1_upper_all)] <- seq(0,540, 1)
target1_1.4sd_phan1_upper_all$group <- seq(1,(nrow(target1_1.4sd_phan1_upper_all)), 1)

target1_1.4sd_phan1_upper_all <- reshape2::melt(target1_1.4sd_phan1_upper_all, id.vars = "group")
names(target1_1.4sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_1.4sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.4sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_1.4sd_phan1_lower <- as.data.frame(cbind(target1_1.4sd_phan1$x, target1_1.4sd_phan1$lower))
names(target1_1.4sd_phan1_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.4sd_phan1_lower$age, target1_1.4sd_phan1_lower$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_1.4sd_phan1_lower, indices) {
  d <- target1_1.4sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.4sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.4sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.4sd_phan1_lower_50$age <- seq(0,540, 1)
names(target1_1.4sd_phan1_lower_50)[1] <- "target1"

target1_1.4sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.4sd_phan1_lower_all)[1:ncol(target1_1.4sd_phan1_lower_all)] <- seq(0,540, 1)
target1_1.4sd_phan1_lower_all$group <- seq(1,(nrow(target1_1.4sd_phan1_lower_all)), 1)

target1_1.4sd_phan1_lower_all <- reshape2::melt(target1_1.4sd_phan1_lower_all, id.vars = "group")
names(target1_1.4sd_phan1_lower_all)[2:3] <- c("age", "target1")

target1_1.4sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2sd

nsig <- 1.2

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_1.2sd_phan1_upper <- as.data.frame(cbind(target1_1.2sd_phan1$x, target1_1.2sd_phan1$upper))
names(target1_1.2sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.2sd_phan1_upper$age, target1_1.2sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_1.2sd_phan1_upper, indices) {
  d <- target1_1.2sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.2sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.2sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.2sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_1.2sd_phan1_upper_50)[1] <- "target1"

target1_1.2sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.2sd_phan1_upper_all)[1:ncol(target1_1.2sd_phan1_upper_all)] <- seq(0,540, 1)
target1_1.2sd_phan1_upper_all$group <- seq(1,(nrow(target1_1.2sd_phan1_upper_all)), 1)

target1_1.2sd_phan1_upper_all <- reshape2::melt(target1_1.2sd_phan1_upper_all, id.vars = "group")
names(target1_1.2sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_1.2sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_1.2sd_phan1_lower <- as.data.frame(cbind(target1_1.2sd_phan1$x, target1_1.2sd_phan1$lower))
names(target1_1.2sd_phan1_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.2sd_phan1_lower$age, target1_1.2sd_phan1_lower$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_1.2sd_phan1_lower, indices) {
  d <- target1_1.2sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.2sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.2sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.2sd_phan1_lower_50$age <- seq(0,540, 1)
names(target1_1.2sd_phan1_lower_50)[1] <- "target1"

target1_1.2sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.2sd_phan1_lower_all)[1:ncol(target1_1.2sd_phan1_lower_all)] <- seq(0,540, 1)
target1_1.2sd_phan1_lower_all$group <- seq(1,(nrow(target1_1.2sd_phan1_lower_all)), 1)

target1_1.2sd_phan1_lower_all <- reshape2::melt(target1_1.2sd_phan1_lower_all, id.vars = "group")
names(target1_1.2sd_phan1_lower_all)[2:3] <- c("age", "target1")

target1_1.2sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0sd

nsig <- 1.0

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.0sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_1.0sd_phan1_upper <- as.data.frame(cbind(target1_1.0sd_phan1$x, target1_1.0sd_phan1$upper))
names(target1_1.0sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.0sd_phan1_upper$age, target1_1.0sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_1.0sd_phan1_upper, indices) {
  d <- target1_1.0sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.0sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.0sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.0sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_1.0sd_phan1_upper_50)[1] <- "target1"

target1_1.0sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.0sd_phan1_upper_all)[1:ncol(target1_1.0sd_phan1_upper_all)] <- seq(0,540, 1)
target1_1.0sd_phan1_upper_all$group <- seq(1,(nrow(target1_1.0sd_phan1_upper_all)), 1)

target1_1.0sd_phan1_upper_all <- reshape2::melt(target1_1.0sd_phan1_upper_all, id.vars = "group")
names(target1_1.0sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_1.0sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.0sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_1.0sd_phan1_lower <- as.data.frame(cbind(target1_1.0sd_phan1$x, target1_1.0sd_phan1$lower))
names(target1_1.0sd_phan1_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.0sd_phan1_lower$age, target1_1.0sd_phan1_lower$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_1.0sd_phan1_lower, indices) {
  d <- target1_1.0sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.0sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.0sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.0sd_phan1_lower_50$age <- seq(0,540, 1)
names(target1_1.0sd_phan1_lower_50)[1] <- "target1"

target1_1.0sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.0sd_phan1_lower_all)[1:ncol(target1_1.0sd_phan1_lower_all)] <- seq(0,540, 1)
target1_1.0sd_phan1_lower_all$group <- seq(1,(nrow(target1_1.0sd_phan1_lower_all)), 1)

target1_1.0sd_phan1_lower_all <- reshape2::melt(target1_1.0sd_phan1_lower_all, id.vars = "group")
names(target1_1.0sd_phan1_lower_all)[2:3] <- c("age", "target1")

target1_1.0sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8sd

nsig <- 0.8

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.8sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_0.8sd_phan1_upper <- as.data.frame(cbind(target1_0.8sd_phan1$x, target1_0.8sd_phan1$upper))
names(target1_0.8sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.8sd_phan1_upper$age, target1_0.8sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_0.8sd_phan1_upper, indices) {
  d <- target1_0.8sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.8sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.8sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.8sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_0.8sd_phan1_upper_50)[1] <- "target1"

target1_0.8sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.8sd_phan1_upper_all)[1:ncol(target1_0.8sd_phan1_upper_all)] <- seq(0,540, 1)
target1_0.8sd_phan1_upper_all$group <- seq(1,(nrow(target1_0.8sd_phan1_upper_all)), 1)

target1_0.8sd_phan1_upper_all <- reshape2::melt(target1_0.8sd_phan1_upper_all, id.vars = "group")
names(target1_0.8sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_0.8sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.8sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_0.8sd_phan1_lower <- as.data.frame(cbind(target1_0.8sd_phan1$x, target1_0.8sd_phan1$lower))
names(target1_0.8sd_phan1_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.8sd_phan1_lower$age, target1_0.8sd_phan1_lower$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_0.8sd_phan1_lower, indices) {
  d <- target1_0.8sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.8sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.8sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.8sd_phan1_lower_50$age <- seq(0,540, 1)
names(target1_0.8sd_phan1_lower_50)[1] <- "target1"

target1_0.8sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.8sd_phan1_lower_all)[1:ncol(target1_0.8sd_phan1_lower_all)] <- seq(0,540, 1)
target1_0.8sd_phan1_lower_all$group <- seq(1,(nrow(target1_0.8sd_phan1_lower_all)), 1)

target1_0.8sd_phan1_lower_all <- reshape2::melt(target1_0.8sd_phan1_lower_all, id.vars = "group")
names(target1_0.8sd_phan1_lower_all)[2:3] <- c("age", "target1")

target1_0.8sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6sd

nsig <- 0.6

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.6sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_0.6sd_phan1_upper <- as.data.frame(cbind(target1_0.6sd_phan1$x, target1_0.6sd_phan1$upper))
names(target1_0.6sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.6sd_phan1_upper$age, target1_0.6sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_0.6sd_phan1_upper, indices) {
  d <- target1_0.6sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.6sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.6sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.6sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_0.6sd_phan1_upper_50)[1] <- "target1"

target1_0.6sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.6sd_phan1_upper_all)[1:ncol(target1_0.6sd_phan1_upper_all)] <- seq(0,540, 1)
target1_0.6sd_phan1_upper_all$group <- seq(1,(nrow(target1_0.6sd_phan1_upper_all)), 1)

target1_0.6sd_phan1_upper_all <- reshape2::melt(target1_0.6sd_phan1_upper_all, id.vars = "group")
names(target1_0.6sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_0.6sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.6sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_0.6sd_phan1_lower <- as.data.frame(cbind(target1_0.6sd_phan1$x, target1_0.6sd_phan1$lower))
names(target1_0.6sd_phan1_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.6sd_phan1_lower$age, target1_0.6sd_phan1_lower$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_0.6sd_phan1_lower, indices) {
  d <- target1_0.6sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.6sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.6sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.6sd_phan1_lower_50$age <- seq(0,540, 1)
names(target1_0.6sd_phan1_lower_50)[1] <- "target1"

target1_0.6sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.6sd_phan1_lower_all)[1:ncol(target1_0.6sd_phan1_lower_all)] <- seq(0,540, 1)
target1_0.6sd_phan1_lower_all$group <- seq(1,(nrow(target1_0.6sd_phan1_lower_all)), 1)

target1_0.6sd_phan1_lower_all <- reshape2::melt(target1_0.6sd_phan1_lower_all, id.vars = "group")
names(target1_0.6sd_phan1_lower_all)[2:3] <- c("age", "target1")

target1_0.6sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4sd

nsig <- 0.4

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.4sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_0.4sd_phan1_upper <- as.data.frame(cbind(target1_0.4sd_phan1$x, target1_0.4sd_phan1$upper))
names(target1_0.4sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.4sd_phan1_upper$age, target1_0.4sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_0.4sd_phan1_upper, indices) {
  d <- target1_0.4sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.4sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.4sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.4sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_0.4sd_phan1_upper_50)[1] <- "target1"

target1_0.4sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.4sd_phan1_upper_all)[1:ncol(target1_0.4sd_phan1_upper_all)] <- seq(0,540, 1)
target1_0.4sd_phan1_upper_all$group <- seq(1,(nrow(target1_0.4sd_phan1_upper_all)), 1)

target1_0.4sd_phan1_upper_all <- reshape2::melt(target1_0.4sd_phan1_upper_all, id.vars = "group")
names(target1_0.4sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_0.4sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.4sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_0.4sd_phan1_lower <- as.data.frame(cbind(target1_0.4sd_phan1$x, target1_0.4sd_phan1$lower))
names(target1_0.4sd_phan1_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.4sd_phan1_lower$age, target1_0.4sd_phan1_lower$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_0.4sd_phan1_lower, indices) {
  d <- target1_0.4sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.4sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.4sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.4sd_phan1_lower_50$age <- seq(0,540, 1)
names(target1_0.4sd_phan1_lower_50)[1] <- "target1"

target1_0.4sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.4sd_phan1_lower_all)[1:ncol(target1_0.4sd_phan1_lower_all)] <- seq(0,540, 1)
target1_0.4sd_phan1_lower_all$group <- seq(1,(nrow(target1_0.4sd_phan1_lower_all)), 1)

target1_0.4sd_phan1_lower_all <- reshape2::melt(target1_0.4sd_phan1_lower_all, id.vars = "group")
names(target1_0.4sd_phan1_lower_all)[2:3] <- c("age", "target1")

target1_0.4sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2sd

nsig <- 0.2

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_0.2sd_phan1_upper <- as.data.frame(cbind(target1_0.2sd_phan1$x, target1_0.2sd_phan1$upper))
names(target1_0.2sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.2sd_phan1_upper$age, target1_0.2sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_0.2sd_phan1_upper, indices) {
  d <- target1_0.2sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.2sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.2sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.2sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_0.2sd_phan1_upper_50)[1] <- "target1"

target1_0.2sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.2sd_phan1_upper_all)[1:ncol(target1_0.2sd_phan1_upper_all)] <- seq(0,540, 1)
target1_0.2sd_phan1_upper_all$group <- seq(1,(nrow(target1_0.2sd_phan1_upper_all)), 1)

target1_0.2sd_phan1_upper_all <- reshape2::melt(target1_0.2sd_phan1_upper_all, id.vars = "group")
names(target1_0.2sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_0.2sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_0.2sd_phan1_lower <- as.data.frame(cbind(target1_0.2sd_phan1$x, target1_0.2sd_phan1$lower))
names(target1_0.2sd_phan1_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.2sd_phan1_lower$age, target1_0.2sd_phan1_lower$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_0.2sd_phan1_lower, indices) {
  d <- target1_0.2sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.2sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.2sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.2sd_phan1_lower_50$age <- seq(0,540, 1)
names(target1_0.2sd_phan1_lower_50)[1] <- "target1"

target1_0.2sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.2sd_phan1_lower_all)[1:ncol(target1_0.2sd_phan1_lower_all)] <- seq(0,540, 1)
target1_0.2sd_phan1_lower_all$group <- seq(1,(nrow(target1_0.2sd_phan1_lower_all)), 1)

target1_0.2sd_phan1_lower_all <- reshape2::melt(target1_0.2sd_phan1_lower_all, id.vars = "group")
names(target1_0.2sd_phan1_lower_all)[2:3] <- c("age", "target1")

target1_0.2sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0sd (once)

nsig <- 0

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target1, degree=1, nsigma = nsig, span = span)
target1_0sd_phan1_upper <- as.data.frame(cbind(target1_0sd_phan1$x, target1_0sd_phan1$upper))
names(target1_0sd_phan1_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0sd_phan1_upper$age, target1_0sd_phan1_upper$target1, degree=1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_0sd_phan1_upper, indices) {
  d <- target1_0sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=1)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0sd_phan1_upper_50$age <- seq(0,540, 1)
names(target1_0sd_phan1_upper_50)[1] <- "target1"

target1_0sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target1_0sd_phan1_upper_all)[1:ncol(target1_0sd_phan1_upper_all)] <- seq(0,540, 1)
target1_0sd_phan1_upper_all$group <- seq(1,(nrow(target1_0sd_phan1_upper_all)), 1)

target1_0sd_phan1_upper_all <- reshape2::melt(target1_0sd_phan1_upper_all, id.vars = "group")
names(target1_0sd_phan1_upper_all)[2:3] <- c("age", "target1")

target1_0sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# Mesozoic-Cenozoic

# phanerozoic

# span_upper <- loess.predict[["pars"]][["span"]]
# span_upper <- loess.predict[["pars"]][["span"]]

# 2sd

nsig <- 2

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_2sd_phan2_upper <- as.data.frame(cbind(target1_2sd_phan2$x, target1_2sd_phan2$upper))
names(target1_2sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_2sd_phan2_upper$age, target1_2sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_2sd_phan2_upper, indices) {
  d <- target1_2sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_2sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_2sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_2sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_2sd_phan2_upper_50)[1] <- "target1"

target1_2sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_2sd_phan2_upper_all)[1:ncol(target1_2sd_phan2_upper_all)] <- seq(0,540, 1)
target1_2sd_phan2_upper_all$group <- seq(1,(nrow(target1_2sd_phan2_upper_all)), 1)

target1_2sd_phan2_upper_all <- reshape2::melt(target1_2sd_phan2_upper_all, id.vars = "group")
names(target1_2sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_2sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_2sd_phan2_lower <- as.data.frame(cbind(target1_2sd_phan2$x, target1_2sd_phan2$lower))
names(target1_2sd_phan2_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_2sd_phan2_lower$age, target1_2sd_phan2_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- 1

# boot lower
boot_fn <- function(target1_2sd_phan2_lower, indices) {
  d <- target1_2sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_2sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_2sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_2sd_phan2_lower_50$age <- seq(0,540, 1)
names(target1_2sd_phan2_lower_50)[1] <- "target1"

target1_2sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target1_2sd_phan2_lower_all)[1:ncol(target1_2sd_phan2_lower_all)] <- seq(0,540, 1)
target1_2sd_phan2_lower_all$group <- seq(1,(nrow(target1_2sd_phan2_lower_all)), 1)

target1_2sd_phan2_lower_all <- reshape2::melt(target1_2sd_phan2_lower_all, id.vars = "group")
names(target1_2sd_phan2_lower_all)[2:3] <- c("age", "target1")

target1_2sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8sd

nsig <- 1.8

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.8sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_1.8sd_phan2_upper <- as.data.frame(cbind(target1_1.8sd_phan2$x, target1_1.8sd_phan2$upper))
names(target1_1.8sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.8sd_phan2_upper$age, target1_1.8sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_1.8sd_phan2_upper, indices) {
  d <- target1_1.8sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.8sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.8sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.8sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_1.8sd_phan2_upper_50)[1] <- "target1"

target1_1.8sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.8sd_phan2_upper_all)[1:ncol(target1_1.8sd_phan2_upper_all)] <- seq(0,540, 1)
target1_1.8sd_phan2_upper_all$group <- seq(1,(nrow(target1_1.8sd_phan2_upper_all)), 1)

target1_1.8sd_phan2_upper_all <- reshape2::melt(target1_1.8sd_phan2_upper_all, id.vars = "group")
names(target1_1.8sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_1.8sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.8sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_1.8sd_phan2_lower <- as.data.frame(cbind(target1_1.8sd_phan2$x, target1_1.8sd_phan2$lower))
names(target1_1.8sd_phan2_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.8sd_phan2_lower$age, target1_1.8sd_phan2_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- 1

# boot lower
boot_fn <- function(target1_1.8sd_phan2_lower, indices) {
  d <- target1_1.8sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.8sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.8sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.8sd_phan2_lower_50$age <- seq(0,540, 1)
names(target1_1.8sd_phan2_lower_50)[1] <- "target1"

target1_1.8sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.8sd_phan2_lower_all)[1:ncol(target1_1.8sd_phan2_lower_all)] <- seq(0,540, 1)
target1_1.8sd_phan2_lower_all$group <- seq(1,(nrow(target1_1.8sd_phan2_lower_all)), 1)

target1_1.8sd_phan2_lower_all <- reshape2::melt(target1_1.8sd_phan2_lower_all, id.vars = "group")
names(target1_1.8sd_phan2_lower_all)[2:3] <- c("age", "target1")

target1_1.8sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6sd

nsig <- 1.6

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.6sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_1.6sd_phan2_upper <- as.data.frame(cbind(target1_1.6sd_phan2$x, target1_1.6sd_phan2$upper))
names(target1_1.6sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.6sd_phan2_upper$age, target1_1.6sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_1.6sd_phan2_upper, indices) {
  d <- target1_1.6sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.6sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.6sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.6sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_1.6sd_phan2_upper_50)[1] <- "target1"

target1_1.6sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.6sd_phan2_upper_all)[1:ncol(target1_1.6sd_phan2_upper_all)] <- seq(0,540, 1)
target1_1.6sd_phan2_upper_all$group <- seq(1,(nrow(target1_1.6sd_phan2_upper_all)), 1)

target1_1.6sd_phan2_upper_all <- reshape2::melt(target1_1.6sd_phan2_upper_all, id.vars = "group")
names(target1_1.6sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_1.6sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.6sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_1.6sd_phan2_lower <- as.data.frame(cbind(target1_1.6sd_phan2$x, target1_1.6sd_phan2$lower))
names(target1_1.6sd_phan2_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.6sd_phan2_lower$age, target1_1.6sd_phan2_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- 1

# boot lower
boot_fn <- function(target1_1.6sd_phan2_lower, indices) {
  d <- target1_1.6sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.6sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.6sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.6sd_phan2_lower_50$age <- seq(0,540, 1)
names(target1_1.6sd_phan2_lower_50)[1] <- "target1"

target1_1.6sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.6sd_phan2_lower_all)[1:ncol(target1_1.6sd_phan2_lower_all)] <- seq(0,540, 1)
target1_1.6sd_phan2_lower_all$group <- seq(1,(nrow(target1_1.6sd_phan2_lower_all)), 1)

target1_1.6sd_phan2_lower_all <- reshape2::melt(target1_1.6sd_phan2_lower_all, id.vars = "group")
names(target1_1.6sd_phan2_lower_all)[2:3] <- c("age", "target1")

target1_1.6sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4sd

nsig <- 1.4

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.4sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_1.4sd_phan2_upper <- as.data.frame(cbind(target1_1.4sd_phan2$x, target1_1.4sd_phan2$upper))
names(target1_1.4sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.4sd_phan2_upper$age, target1_1.4sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_1.4sd_phan2_upper, indices) {
  d <- target1_1.4sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.4sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.4sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.4sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_1.4sd_phan2_upper_50)[1] <- "target1"

target1_1.4sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.4sd_phan2_upper_all)[1:ncol(target1_1.4sd_phan2_upper_all)] <- seq(0,540, 1)
target1_1.4sd_phan2_upper_all$group <- seq(1,(nrow(target1_1.4sd_phan2_upper_all)), 1)

target1_1.4sd_phan2_upper_all <- reshape2::melt(target1_1.4sd_phan2_upper_all, id.vars = "group")
names(target1_1.4sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_1.4sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.4sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_1.4sd_phan2_lower <- as.data.frame(cbind(target1_1.4sd_phan2$x, target1_1.4sd_phan2$lower))
names(target1_1.4sd_phan2_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.4sd_phan2_lower$age, target1_1.4sd_phan2_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- 1

# boot lower
boot_fn <- function(target1_1.4sd_phan2_lower, indices) {
  d <- target1_1.4sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.4sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.4sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.4sd_phan2_lower_50$age <- seq(0,540, 1)
names(target1_1.4sd_phan2_lower_50)[1] <- "target1"

target1_1.4sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.4sd_phan2_lower_all)[1:ncol(target1_1.4sd_phan2_lower_all)] <- seq(0,540, 1)
target1_1.4sd_phan2_lower_all$group <- seq(1,(nrow(target1_1.4sd_phan2_lower_all)), 1)

target1_1.4sd_phan2_lower_all <- reshape2::melt(target1_1.4sd_phan2_lower_all, id.vars = "group")
names(target1_1.4sd_phan2_lower_all)[2:3] <- c("age", "target1")

target1_1.4sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2sd

nsig <- 1.2

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_1.2sd_phan2_upper <- as.data.frame(cbind(target1_1.2sd_phan2$x, target1_1.2sd_phan2$upper))
names(target1_1.2sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.2sd_phan2_upper$age, target1_1.2sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_1.2sd_phan2_upper, indices) {
  d <- target1_1.2sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.2sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.2sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.2sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_1.2sd_phan2_upper_50)[1] <- "target1"

target1_1.2sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.2sd_phan2_upper_all)[1:ncol(target1_1.2sd_phan2_upper_all)] <- seq(0,540, 1)
target1_1.2sd_phan2_upper_all$group <- seq(1,(nrow(target1_1.2sd_phan2_upper_all)), 1)

target1_1.2sd_phan2_upper_all <- reshape2::melt(target1_1.2sd_phan2_upper_all, id.vars = "group")
names(target1_1.2sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_1.2sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_1.2sd_phan2_lower <- as.data.frame(cbind(target1_1.2sd_phan2$x, target1_1.2sd_phan2$lower))
names(target1_1.2sd_phan2_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.2sd_phan2_lower$age, target1_1.2sd_phan2_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- 1

# boot lower
boot_fn <- function(target1_1.2sd_phan2_lower, indices) {
  d <- target1_1.2sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.2sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.2sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.2sd_phan2_lower_50$age <- seq(0,540, 1)
names(target1_1.2sd_phan2_lower_50)[1] <- "target1"

target1_1.2sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.2sd_phan2_lower_all)[1:ncol(target1_1.2sd_phan2_lower_all)] <- seq(0,540, 1)
target1_1.2sd_phan2_lower_all$group <- seq(1,(nrow(target1_1.2sd_phan2_lower_all)), 1)

target1_1.2sd_phan2_lower_all <- reshape2::melt(target1_1.2sd_phan2_lower_all, id.vars = "group")
names(target1_1.2sd_phan2_lower_all)[2:3] <- c("age", "target1")

target1_1.2sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0sd

nsig <- 1.0

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.0sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_1.0sd_phan2_upper <- as.data.frame(cbind(target1_1.0sd_phan2$x, target1_1.0sd_phan2$upper))
names(target1_1.0sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.0sd_phan2_upper$age, target1_1.0sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_1.0sd_phan2_upper, indices) {
  d <- target1_1.0sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.0sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.0sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.0sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_1.0sd_phan2_upper_50)[1] <- "target1"

target1_1.0sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.0sd_phan2_upper_all)[1:ncol(target1_1.0sd_phan2_upper_all)] <- seq(0,540, 1)
target1_1.0sd_phan2_upper_all$group <- seq(1,(nrow(target1_1.0sd_phan2_upper_all)), 1)

target1_1.0sd_phan2_upper_all <- reshape2::melt(target1_1.0sd_phan2_upper_all, id.vars = "group")
names(target1_1.0sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_1.0sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.0sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_1.0sd_phan2_lower <- as.data.frame(cbind(target1_1.0sd_phan2$x, target1_1.0sd_phan2$lower))
names(target1_1.0sd_phan2_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.0sd_phan2_lower$age, target1_1.0sd_phan2_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- 1

# boot lower
boot_fn <- function(target1_1.0sd_phan2_lower, indices) {
  d <- target1_1.0sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.0sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.0sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.0sd_phan2_lower_50$age <- seq(0,540, 1)
names(target1_1.0sd_phan2_lower_50)[1] <- "target1"

target1_1.0sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.0sd_phan2_lower_all)[1:ncol(target1_1.0sd_phan2_lower_all)] <- seq(0,540, 1)
target1_1.0sd_phan2_lower_all$group <- seq(1,(nrow(target1_1.0sd_phan2_lower_all)), 1)

target1_1.0sd_phan2_lower_all <- reshape2::melt(target1_1.0sd_phan2_lower_all, id.vars = "group")
names(target1_1.0sd_phan2_lower_all)[2:3] <- c("age", "target1")

target1_1.0sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8sd

nsig <- 0.8

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.8sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_0.8sd_phan2_upper <- as.data.frame(cbind(target1_0.8sd_phan2$x, target1_0.8sd_phan2$upper))
names(target1_0.8sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.8sd_phan2_upper$age, target1_0.8sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_0.8sd_phan2_upper, indices) {
  d <- target1_0.8sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.8sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.8sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.8sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_0.8sd_phan2_upper_50)[1] <- "target1"

target1_0.8sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.8sd_phan2_upper_all)[1:ncol(target1_0.8sd_phan2_upper_all)] <- seq(0,540, 1)
target1_0.8sd_phan2_upper_all$group <- seq(1,(nrow(target1_0.8sd_phan2_upper_all)), 1)

target1_0.8sd_phan2_upper_all <- reshape2::melt(target1_0.8sd_phan2_upper_all, id.vars = "group")
names(target1_0.8sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_0.8sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.8sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_0.8sd_phan2_lower <- as.data.frame(cbind(target1_0.8sd_phan2$x, target1_0.8sd_phan2$lower))
names(target1_0.8sd_phan2_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.8sd_phan2_lower$age, target1_0.8sd_phan2_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- 1

# boot lower
boot_fn <- function(target1_0.8sd_phan2_lower, indices) {
  d <- target1_0.8sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.8sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.8sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.8sd_phan2_lower_50$age <- seq(0,540, 1)
names(target1_0.8sd_phan2_lower_50)[1] <- "target1"

target1_0.8sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.8sd_phan2_lower_all)[1:ncol(target1_0.8sd_phan2_lower_all)] <- seq(0,540, 1)
target1_0.8sd_phan2_lower_all$group <- seq(1,(nrow(target1_0.8sd_phan2_lower_all)), 1)

target1_0.8sd_phan2_lower_all <- reshape2::melt(target1_0.8sd_phan2_lower_all, id.vars = "group")
names(target1_0.8sd_phan2_lower_all)[2:3] <- c("age", "target1")

target1_0.8sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6sd

nsig <- 0.6

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.6sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_0.6sd_phan2_upper <- as.data.frame(cbind(target1_0.6sd_phan2$x, target1_0.6sd_phan2$upper))
names(target1_0.6sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.6sd_phan2_upper$age, target1_0.6sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_0.6sd_phan2_upper, indices) {
  d <- target1_0.6sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.6sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.6sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.6sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_0.6sd_phan2_upper_50)[1] <- "target1"

target1_0.6sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.6sd_phan2_upper_all)[1:ncol(target1_0.6sd_phan2_upper_all)] <- seq(0,540, 1)
target1_0.6sd_phan2_upper_all$group <- seq(1,(nrow(target1_0.6sd_phan2_upper_all)), 1)

target1_0.6sd_phan2_upper_all <- reshape2::melt(target1_0.6sd_phan2_upper_all, id.vars = "group")
names(target1_0.6sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_0.6sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.6sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_0.6sd_phan2_lower <- as.data.frame(cbind(target1_0.6sd_phan2$x, target1_0.6sd_phan2$lower))
names(target1_0.6sd_phan2_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.6sd_phan2_lower$age, target1_0.6sd_phan2_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- 1

# boot lower
boot_fn <- function(target1_0.6sd_phan2_lower, indices) {
  d <- target1_0.6sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.6sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.6sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.6sd_phan2_lower_50$age <- seq(0,540, 1)
names(target1_0.6sd_phan2_lower_50)[1] <- "target1"

target1_0.6sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.6sd_phan2_lower_all)[1:ncol(target1_0.6sd_phan2_lower_all)] <- seq(0,540, 1)
target1_0.6sd_phan2_lower_all$group <- seq(1,(nrow(target1_0.6sd_phan2_lower_all)), 1)

target1_0.6sd_phan2_lower_all <- reshape2::melt(target1_0.6sd_phan2_lower_all, id.vars = "group")
names(target1_0.6sd_phan2_lower_all)[2:3] <- c("age", "target1")

target1_0.6sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4sd

nsig <- 0.4

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.4sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_0.4sd_phan2_upper <- as.data.frame(cbind(target1_0.4sd_phan2$x, target1_0.4sd_phan2$upper))
names(target1_0.4sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.4sd_phan2_upper$age, target1_0.4sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_0.4sd_phan2_upper, indices) {
  d <- target1_0.4sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.4sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.4sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.4sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_0.4sd_phan2_upper_50)[1] <- "target1"

target1_0.4sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.4sd_phan2_upper_all)[1:ncol(target1_0.4sd_phan2_upper_all)] <- seq(0,540, 1)
target1_0.4sd_phan2_upper_all$group <- seq(1,(nrow(target1_0.4sd_phan2_upper_all)), 1)

target1_0.4sd_phan2_upper_all <- reshape2::melt(target1_0.4sd_phan2_upper_all, id.vars = "group")
names(target1_0.4sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_0.4sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.4sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_0.4sd_phan2_lower <- as.data.frame(cbind(target1_0.4sd_phan2$x, target1_0.4sd_phan2$lower))
names(target1_0.4sd_phan2_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.4sd_phan2_lower$age, target1_0.4sd_phan2_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- 1

# boot lower
boot_fn <- function(target1_0.4sd_phan2_lower, indices) {
  d <- target1_0.4sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.4sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.4sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.4sd_phan2_lower_50$age <- seq(0,540, 1)
names(target1_0.4sd_phan2_lower_50)[1] <- "target1"

target1_0.4sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.4sd_phan2_lower_all)[1:ncol(target1_0.4sd_phan2_lower_all)] <- seq(0,540, 1)
target1_0.4sd_phan2_lower_all$group <- seq(1,(nrow(target1_0.4sd_phan2_lower_all)), 1)

target1_0.4sd_phan2_lower_all <- reshape2::melt(target1_0.4sd_phan2_lower_all, id.vars = "group")
names(target1_0.4sd_phan2_lower_all)[2:3] <- c("age", "target1")

target1_0.4sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2sd

nsig <- 0.2

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_0.2sd_phan2_upper <- as.data.frame(cbind(target1_0.2sd_phan2$x, target1_0.2sd_phan2$upper))
names(target1_0.2sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.2sd_phan2_upper$age, target1_0.2sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_0.2sd_phan2_upper, indices) {
  d <- target1_0.2sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.2sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.2sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.2sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_0.2sd_phan2_upper_50)[1] <- "target1"

target1_0.2sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.2sd_phan2_upper_all)[1:ncol(target1_0.2sd_phan2_upper_all)] <- seq(0,540, 1)
target1_0.2sd_phan2_upper_all$group <- seq(1,(nrow(target1_0.2sd_phan2_upper_all)), 1)

target1_0.2sd_phan2_upper_all <- reshape2::melt(target1_0.2sd_phan2_upper_all, id.vars = "group")
names(target1_0.2sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_0.2sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_0.2sd_phan2_lower <- as.data.frame(cbind(target1_0.2sd_phan2$x, target1_0.2sd_phan2$lower))
names(target1_0.2sd_phan2_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.2sd_phan2_lower$age, target1_0.2sd_phan2_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- 1

# boot lower
boot_fn <- function(target1_0.2sd_phan2_lower, indices) {
  d <- target1_0.2sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.2sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.2sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.2sd_phan2_lower_50$age <- seq(0,540, 1)
names(target1_0.2sd_phan2_lower_50)[1] <- "target1"

target1_0.2sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.2sd_phan2_lower_all)[1:ncol(target1_0.2sd_phan2_lower_all)] <- seq(0,540, 1)
target1_0.2sd_phan2_lower_all$group <- seq(1,(nrow(target1_0.2sd_phan2_lower_all)), 1)

target1_0.2sd_phan2_lower_all <- reshape2::melt(target1_0.2sd_phan2_lower_all, id.vars = "group")
names(target1_0.2sd_phan2_lower_all)[2:3] <- c("age", "target1")

target1_0.2sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0sd (once)

nsig <- 0

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target1, degree=0, nsigma = nsig, span = span)
target1_0sd_phan2_upper <- as.data.frame(cbind(target1_0sd_phan2$x, target1_0sd_phan2$upper))
names(target1_0sd_phan2_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0sd_phan2_upper$age, target1_0sd_phan2_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- 1

# boot upper
boot_fn <- function(target1_0sd_phan2_upper, indices) {
  d <- target1_0sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0sd_phan2_upper_50$age <- seq(0,540, 1)
names(target1_0sd_phan2_upper_50)[1] <- "target1"

target1_0sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target1_0sd_phan2_upper_all)[1:ncol(target1_0sd_phan2_upper_all)] <- seq(0,540, 1)
target1_0sd_phan2_upper_all$group <- seq(1,(nrow(target1_0sd_phan2_upper_all)), 1)

target1_0sd_phan2_upper_all <- reshape2::melt(target1_0sd_phan2_upper_all, id.vars = "group")
names(target1_0sd_phan2_upper_all)[2:3] <- c("age", "target1")

target1_0sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# precambrian


# 2sd

nsig <- 2

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_2sd_pre_upper <- as.data.frame(cbind(target1_2sd_pre$x, target1_2sd_pre$upper))
names(target1_2sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_2sd_pre_upper$age, target1_2sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_2sd_pre_upper, indices) {
  d <- target1_2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_2sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_2sd_pre_upper_50)[1] <- "target1"

target1_2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_2sd_pre_upper_all)[1:ncol(target1_2sd_pre_upper_all)] <- seq(540,4000, 10)
target1_2sd_pre_upper_all$group <- seq(1,(nrow(target1_2sd_pre_upper_all)), 1)

target1_2sd_pre_upper_all <- reshape2::melt(target1_2sd_pre_upper_all, id.vars = "group")
names(target1_2sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_2sd_pre_lower <- as.data.frame(cbind(target1_2sd_pre$x, target1_2sd_pre$lower))
names(target1_2sd_pre_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_2sd_pre_lower$age, target1_2sd_pre_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_2sd_pre_lower, indices) {
  d <- target1_2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_2sd_pre_lower_50$age <- seq(540,4000, 10)
names(target1_2sd_pre_lower_50)[1] <- "target1"

target1_2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target1_2sd_pre_lower_all)[1:ncol(target1_2sd_pre_lower_all)] <- seq(540,4000, 10)
target1_2sd_pre_lower_all$group <- seq(1,(nrow(target1_2sd_pre_lower_all)), 1)

target1_2sd_pre_lower_all <- reshape2::melt(target1_2sd_pre_lower_all, id.vars = "group")
names(target1_2sd_pre_lower_all)[2:3] <- c("age", "target1")

target1_2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.8sd

nsig <- 1.8

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.8sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_1.8sd_pre_upper <- as.data.frame(cbind(target1_1.8sd_pre$x, target1_1.8sd_pre$upper))
names(target1_1.8sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.8sd_pre_upper$age, target1_1.8sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_1.8sd_pre_upper, indices) {
  d <- target1_1.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_1.8sd_pre_upper_50)[1] <- "target1"

target1_1.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.8sd_pre_upper_all)[1:ncol(target1_1.8sd_pre_upper_all)] <- seq(540,4000, 10)
target1_1.8sd_pre_upper_all$group <- seq(1,(nrow(target1_1.8sd_pre_upper_all)), 1)

target1_1.8sd_pre_upper_all <- reshape2::melt(target1_1.8sd_pre_upper_all, id.vars = "group")
names(target1_1.8sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_1.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.8sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_1.8sd_pre_lower <- as.data.frame(cbind(target1_1.8sd_pre$x, target1_1.8sd_pre$lower))
names(target1_1.8sd_pre_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.8sd_pre_lower$age, target1_1.8sd_pre_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_1.8sd_pre_lower, indices) {
  d <- target1_1.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(target1_1.8sd_pre_lower_50)[1] <- "target1"

target1_1.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.8sd_pre_lower_all)[1:ncol(target1_1.8sd_pre_lower_all)] <- seq(540,4000, 10)
target1_1.8sd_pre_lower_all$group <- seq(1,(nrow(target1_1.8sd_pre_lower_all)), 1)

target1_1.8sd_pre_lower_all <- reshape2::melt(target1_1.8sd_pre_lower_all, id.vars = "group")
names(target1_1.8sd_pre_lower_all)[2:3] <- c("age", "target1")

target1_1.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.6sd

nsig <- 1.6

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.6sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_1.6sd_pre_upper <- as.data.frame(cbind(target1_1.6sd_pre$x, target1_1.6sd_pre$upper))
names(target1_1.6sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.6sd_pre_upper$age, target1_1.6sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_1.6sd_pre_upper, indices) {
  d <- target1_1.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_1.6sd_pre_upper_50)[1] <- "target1"

target1_1.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.6sd_pre_upper_all)[1:ncol(target1_1.6sd_pre_upper_all)] <- seq(540,4000, 10)
target1_1.6sd_pre_upper_all$group <- seq(1,(nrow(target1_1.6sd_pre_upper_all)), 1)

target1_1.6sd_pre_upper_all <- reshape2::melt(target1_1.6sd_pre_upper_all, id.vars = "group")
names(target1_1.6sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_1.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.6sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_1.6sd_pre_lower <- as.data.frame(cbind(target1_1.6sd_pre$x, target1_1.6sd_pre$lower))
names(target1_1.6sd_pre_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.6sd_pre_lower$age, target1_1.6sd_pre_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_1.6sd_pre_lower, indices) {
  d <- target1_1.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(target1_1.6sd_pre_lower_50)[1] <- "target1"

target1_1.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.6sd_pre_lower_all)[1:ncol(target1_1.6sd_pre_lower_all)] <- seq(540,4000, 10)
target1_1.6sd_pre_lower_all$group <- seq(1,(nrow(target1_1.6sd_pre_lower_all)), 1)

target1_1.6sd_pre_lower_all <- reshape2::melt(target1_1.6sd_pre_lower_all, id.vars = "group")
names(target1_1.6sd_pre_lower_all)[2:3] <- c("age", "target1")

target1_1.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.4sd

nsig <- 1.4

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.4sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_1.4sd_pre_upper <- as.data.frame(cbind(target1_1.4sd_pre$x, target1_1.4sd_pre$upper))
names(target1_1.4sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.4sd_pre_upper$age, target1_1.4sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_1.4sd_pre_upper, indices) {
  d <- target1_1.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_1.4sd_pre_upper_50)[1] <- "target1"

target1_1.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.4sd_pre_upper_all)[1:ncol(target1_1.4sd_pre_upper_all)] <- seq(540,4000, 10)
target1_1.4sd_pre_upper_all$group <- seq(1,(nrow(target1_1.4sd_pre_upper_all)), 1)

target1_1.4sd_pre_upper_all <- reshape2::melt(target1_1.4sd_pre_upper_all, id.vars = "group")
names(target1_1.4sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_1.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.4sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_1.4sd_pre_lower <- as.data.frame(cbind(target1_1.4sd_pre$x, target1_1.4sd_pre$lower))
names(target1_1.4sd_pre_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.4sd_pre_lower$age, target1_1.4sd_pre_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_1.4sd_pre_lower, indices) {
  d <- target1_1.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(target1_1.4sd_pre_lower_50)[1] <- "target1"

target1_1.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.4sd_pre_lower_all)[1:ncol(target1_1.4sd_pre_lower_all)] <- seq(540,4000, 10)
target1_1.4sd_pre_lower_all$group <- seq(1,(nrow(target1_1.4sd_pre_lower_all)), 1)

target1_1.4sd_pre_lower_all <- reshape2::melt(target1_1.4sd_pre_lower_all, id.vars = "group")
names(target1_1.4sd_pre_lower_all)[2:3] <- c("age", "target1")

target1_1.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.2sd

nsig <- 1.2

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_1.2sd_pre_upper <- as.data.frame(cbind(target1_1.2sd_pre$x, target1_1.2sd_pre$upper))
names(target1_1.2sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.2sd_pre_upper$age, target1_1.2sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_1.2sd_pre_upper, indices) {
  d <- target1_1.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_1.2sd_pre_upper_50)[1] <- "target1"

target1_1.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.2sd_pre_upper_all)[1:ncol(target1_1.2sd_pre_upper_all)] <- seq(540,4000, 10)
target1_1.2sd_pre_upper_all$group <- seq(1,(nrow(target1_1.2sd_pre_upper_all)), 1)

target1_1.2sd_pre_upper_all <- reshape2::melt(target1_1.2sd_pre_upper_all, id.vars = "group")
names(target1_1.2sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_1.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_1.2sd_pre_lower <- as.data.frame(cbind(target1_1.2sd_pre$x, target1_1.2sd_pre$lower))
names(target1_1.2sd_pre_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.2sd_pre_lower$age, target1_1.2sd_pre_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_1.2sd_pre_lower, indices) {
  d <- target1_1.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(target1_1.2sd_pre_lower_50)[1] <- "target1"

target1_1.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.2sd_pre_lower_all)[1:ncol(target1_1.2sd_pre_lower_all)] <- seq(540,4000, 10)
target1_1.2sd_pre_lower_all$group <- seq(1,(nrow(target1_1.2sd_pre_lower_all)), 1)

target1_1.2sd_pre_lower_all <- reshape2::melt(target1_1.2sd_pre_lower_all, id.vars = "group")
names(target1_1.2sd_pre_lower_all)[2:3] <- c("age", "target1")

target1_1.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.0sd

nsig <- 1.0

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.0sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_1.0sd_pre_upper <- as.data.frame(cbind(target1_1.0sd_pre$x, target1_1.0sd_pre$upper))
names(target1_1.0sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.0sd_pre_upper$age, target1_1.0sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_1.0sd_pre_upper, indices) {
  d <- target1_1.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_1.0sd_pre_upper_50)[1] <- "target1"

target1_1.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_1.0sd_pre_upper_all)[1:ncol(target1_1.0sd_pre_upper_all)] <- seq(540,4000, 10)
target1_1.0sd_pre_upper_all$group <- seq(1,(nrow(target1_1.0sd_pre_upper_all)), 1)

target1_1.0sd_pre_upper_all <- reshape2::melt(target1_1.0sd_pre_upper_all, id.vars = "group")
names(target1_1.0sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_1.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_1.0sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_1.0sd_pre_lower <- as.data.frame(cbind(target1_1.0sd_pre$x, target1_1.0sd_pre$lower))
names(target1_1.0sd_pre_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_1.0sd_pre_lower$age, target1_1.0sd_pre_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_1.0sd_pre_lower, indices) {
  d <- target1_1.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_1.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_1.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_1.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(target1_1.0sd_pre_lower_50)[1] <- "target1"

target1_1.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target1_1.0sd_pre_lower_all)[1:ncol(target1_1.0sd_pre_lower_all)] <- seq(540,4000, 10)
target1_1.0sd_pre_lower_all$group <- seq(1,(nrow(target1_1.0sd_pre_lower_all)), 1)

target1_1.0sd_pre_lower_all <- reshape2::melt(target1_1.0sd_pre_lower_all, id.vars = "group")
names(target1_1.0sd_pre_lower_all)[2:3] <- c("age", "target1")

target1_1.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.8sd

nsig <- 0.8

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.8sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_0.8sd_pre_upper <- as.data.frame(cbind(target1_0.8sd_pre$x, target1_0.8sd_pre$upper))
names(target1_0.8sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.8sd_pre_upper$age, target1_0.8sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_0.8sd_pre_upper, indices) {
  d <- target1_0.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_0.8sd_pre_upper_50)[1] <- "target1"

target1_0.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.8sd_pre_upper_all)[1:ncol(target1_0.8sd_pre_upper_all)] <- seq(540,4000, 10)
target1_0.8sd_pre_upper_all$group <- seq(1,(nrow(target1_0.8sd_pre_upper_all)), 1)

target1_0.8sd_pre_upper_all <- reshape2::melt(target1_0.8sd_pre_upper_all, id.vars = "group")
names(target1_0.8sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_0.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.8sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_0.8sd_pre_lower <- as.data.frame(cbind(target1_0.8sd_pre$x, target1_0.8sd_pre$lower))
names(target1_0.8sd_pre_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.8sd_pre_lower$age, target1_0.8sd_pre_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_0.8sd_pre_lower, indices) {
  d <- target1_0.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(target1_0.8sd_pre_lower_50)[1] <- "target1"

target1_0.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.8sd_pre_lower_all)[1:ncol(target1_0.8sd_pre_lower_all)] <- seq(540,4000, 10)
target1_0.8sd_pre_lower_all$group <- seq(1,(nrow(target1_0.8sd_pre_lower_all)), 1)

target1_0.8sd_pre_lower_all <- reshape2::melt(target1_0.8sd_pre_lower_all, id.vars = "group")
names(target1_0.8sd_pre_lower_all)[2:3] <- c("age", "target1")

target1_0.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.6sd

nsig <- 0.6

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.6sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_0.6sd_pre_upper <- as.data.frame(cbind(target1_0.6sd_pre$x, target1_0.6sd_pre$upper))
names(target1_0.6sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.6sd_pre_upper$age, target1_0.6sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_0.6sd_pre_upper, indices) {
  d <- target1_0.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_0.6sd_pre_upper_50)[1] <- "target1"

target1_0.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.6sd_pre_upper_all)[1:ncol(target1_0.6sd_pre_upper_all)] <- seq(540,4000, 10)
target1_0.6sd_pre_upper_all$group <- seq(1,(nrow(target1_0.6sd_pre_upper_all)), 1)

target1_0.6sd_pre_upper_all <- reshape2::melt(target1_0.6sd_pre_upper_all, id.vars = "group")
names(target1_0.6sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_0.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.6sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_0.6sd_pre_lower <- as.data.frame(cbind(target1_0.6sd_pre$x, target1_0.6sd_pre$lower))
names(target1_0.6sd_pre_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.6sd_pre_lower$age, target1_0.6sd_pre_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_0.6sd_pre_lower, indices) {
  d <- target1_0.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(target1_0.6sd_pre_lower_50)[1] <- "target1"

target1_0.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.6sd_pre_lower_all)[1:ncol(target1_0.6sd_pre_lower_all)] <- seq(540,4000, 10)
target1_0.6sd_pre_lower_all$group <- seq(1,(nrow(target1_0.6sd_pre_lower_all)), 1)

target1_0.6sd_pre_lower_all <- reshape2::melt(target1_0.6sd_pre_lower_all, id.vars = "group")
names(target1_0.6sd_pre_lower_all)[2:3] <- c("age", "target1")

target1_0.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.4sd

nsig <- 0.4

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.4sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_0.4sd_pre_upper <- as.data.frame(cbind(target1_0.4sd_pre$x, target1_0.4sd_pre$upper))
names(target1_0.4sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.4sd_pre_upper$age, target1_0.4sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_0.4sd_pre_upper, indices) {
  d <- target1_0.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_0.4sd_pre_upper_50)[1] <- "target1"

target1_0.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.4sd_pre_upper_all)[1:ncol(target1_0.4sd_pre_upper_all)] <- seq(540,4000, 10)
target1_0.4sd_pre_upper_all$group <- seq(1,(nrow(target1_0.4sd_pre_upper_all)), 1)

target1_0.4sd_pre_upper_all <- reshape2::melt(target1_0.4sd_pre_upper_all, id.vars = "group")
names(target1_0.4sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_0.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.4sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_0.4sd_pre_lower <- as.data.frame(cbind(target1_0.4sd_pre$x, target1_0.4sd_pre$lower))
names(target1_0.4sd_pre_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.4sd_pre_lower$age, target1_0.4sd_pre_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_0.4sd_pre_lower, indices) {
  d <- target1_0.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(target1_0.4sd_pre_lower_50)[1] <- "target1"

target1_0.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.4sd_pre_lower_all)[1:ncol(target1_0.4sd_pre_lower_all)] <- seq(540,4000, 10)
target1_0.4sd_pre_lower_all$group <- seq(1,(nrow(target1_0.4sd_pre_lower_all)), 1)

target1_0.4sd_pre_lower_all <- reshape2::melt(target1_0.4sd_pre_lower_all, id.vars = "group")
names(target1_0.4sd_pre_lower_all)[2:3] <- c("age", "target1")

target1_0.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.2sd

nsig <- 0.2

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_0.2sd_pre_upper <- as.data.frame(cbind(target1_0.2sd_pre$x, target1_0.2sd_pre$upper))
names(target1_0.2sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.2sd_pre_upper$age, target1_0.2sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_0.2sd_pre_upper, indices) {
  d <- target1_0.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_0.2sd_pre_upper_50)[1] <- "target1"

target1_0.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_0.2sd_pre_upper_all)[1:ncol(target1_0.2sd_pre_upper_all)] <- seq(540,4000, 10)
target1_0.2sd_pre_upper_all$group <- seq(1,(nrow(target1_0.2sd_pre_upper_all)), 1)

target1_0.2sd_pre_upper_all <- reshape2::melt(target1_0.2sd_pre_upper_all, id.vars = "group")
names(target1_0.2sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_0.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target1 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0.2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_0.2sd_pre_lower <- as.data.frame(cbind(target1_0.2sd_pre$x, target1_0.2sd_pre$lower))
names(target1_0.2sd_pre_lower)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0.2sd_pre_lower$age, target1_0.2sd_pre_lower$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target1_0.2sd_pre_lower, indices) {
  d <- target1_0.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(target1_0.2sd_pre_lower_50)[1] <- "target1"

target1_0.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target1_0.2sd_pre_lower_all)[1:ncol(target1_0.2sd_pre_lower_all)] <- seq(540,4000, 10)
target1_0.2sd_pre_lower_all$group <- seq(1,(nrow(target1_0.2sd_pre_lower_all)), 1)

target1_0.2sd_pre_lower_all <- reshape2::melt(target1_0.2sd_pre_lower_all, id.vars = "group")
names(target1_0.2sd_pre_lower_all)[2:3] <- c("age", "target1")

target1_0.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0sd (once)

nsig <- 0

# target1 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target1_0sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target1, degree=0, nsigma = nsig, span = span)
target1_0sd_pre_upper <- as.data.frame(cbind(target1_0sd_pre$x, target1_0sd_pre$upper))
names(target1_0sd_pre_upper)[1:2] <- c("age","target1")

# CV spans
loess.predict <- loess.as(target1_0sd_pre_upper$age, target1_0sd_pre_upper$target1, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target1_0sd_pre_upper, indices) {
  d <- target1_0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target1 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target1_0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target1_0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target1_0sd_pre_upper_50$age <- seq(540,4000, 10)
names(target1_0sd_pre_upper_50)[1] <- "target1"

target1_0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target1_0sd_pre_upper_all)[1:ncol(target1_0sd_pre_upper_all)] <- seq(540,4000, 10)
target1_0sd_pre_upper_all$group <- seq(1,(nrow(target1_0sd_pre_upper_all)), 1)

target1_0sd_pre_upper_all <- reshape2::melt(target1_0sd_pre_upper_all, id.vars = "group")
names(target1_0sd_pre_upper_all)[2:3] <- c("age", "target1")

target1_0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))


x <- ggplot() + 
  geom_line(data = target1_2sd_pre_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_2sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.8sd_pre_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.8sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.6sd_pre_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.6sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.4sd_pre_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.4sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.2sd_pre_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.2sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.0sd_pre_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.0sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.8sd_pre_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.8sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.6sd_pre_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.6sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.4sd_pre_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.4sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.2sd_pre_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.2sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0sd_pre_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_point(data = SGP2, aes(age, target1)) + 
  scale_x_reverse(limits = c(4000,540)) + theme_bw() +
  coord_cartesian(ylim = c(0, 1))


y <- ggplot() + 
  geom_line(data = target1_2sd_phan1_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_2sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.8sd_phan1_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.8sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.6sd_phan1_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.6sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.4sd_phan1_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.4sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.2sd_phan1_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.2sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.0sd_phan1_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.0sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.8sd_phan1_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.8sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.6sd_phan1_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.6sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.4sd_phan1_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.4sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.2sd_phan1_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.2sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0sd_phan1_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_2sd_phan2_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_2sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.8sd_phan2_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.8sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.6sd_phan2_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.6sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.4sd_phan2_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.4sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.2sd_phan2_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.2sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.0sd_phan2_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_1.0sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.8sd_phan2_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.8sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.6sd_phan2_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.6sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.4sd_phan2_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.4sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.2sd_phan2_lower_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0.2sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_line(data = target1_0sd_phan2_upper_all, aes(age2,target1, group = group), colour = 'blue') +
  geom_point(data = SGP2, aes(age, target1)) + 
  scale_x_reverse(limits = c(540,0)) + theme_bw() +
  coord_cartesian(ylim = c(0, 1))

grid.arrange(x,y,ncol = 2)

##### ribbons ####

##### Median ribbons and raw plots ####

phan1_2sd <- cbind(target1_2sd_phan1_upper_50,target1_2sd_phan1_lower_50$target1)
names(phan1_2sd)[1:3] <- c("upper","age","lower")
phan2_2sd <- cbind(target1_2sd_phan2_upper_50,target1_2sd_phan2_lower_50$target1)
names(phan2_2sd)[1:3] <- c("upper","age","lower")
phan_2sd <- rbind(phan1_2sd,phan2_2sd)
phan_2sd <- aggregate(cbind(lower,upper) ~ age, phan_2sd, mean )

phan1_1.8sd <- cbind(target1_1.8sd_phan1_upper_50,target1_1.8sd_phan1_lower_50$target1)
names(phan1_1.8sd)[1:3] <- c("upper","age","lower")
phan2_1.8sd <- cbind(target1_1.8sd_phan2_upper_50,target1_1.8sd_phan2_lower_50$target1)
names(phan2_1.8sd)[1:3] <- c("upper","age","lower")
phan_1.8sd <- rbind(phan1_1.8sd,phan2_1.8sd)
phan_1.8sd <- aggregate(cbind(lower,upper) ~ age, phan_1.8sd, mean )

phan1_1.6sd <- cbind(target1_1.6sd_phan1_upper_50,target1_1.6sd_phan1_lower_50$target1)
names(phan1_1.6sd)[1:3] <- c("upper","age","lower")
phan2_1.6sd <- cbind(target1_1.6sd_phan2_upper_50,target1_1.6sd_phan2_lower_50$target1)
names(phan2_1.6sd)[1:3] <- c("upper","age","lower")
phan_1.6sd <- rbind(phan1_1.6sd,phan2_1.6sd)
phan_1.6sd <- aggregate(cbind(lower,upper) ~ age, phan_1.6sd, mean )

phan1_1.4sd <- cbind(target1_1.4sd_phan1_upper_50,target1_1.4sd_phan1_lower_50$target1)
names(phan1_1.4sd)[1:3] <- c("upper","age","lower")
phan2_1.4sd <- cbind(target1_1.4sd_phan2_upper_50,target1_1.4sd_phan2_lower_50$target1)
names(phan2_1.4sd)[1:3] <- c("upper","age","lower")
phan_1.4sd <- rbind(phan1_1.4sd,phan2_1.4sd)
phan_1.4sd <- aggregate(cbind(lower,upper) ~ age, phan_1.4sd, mean )

phan1_1.2sd <- cbind(target1_1.2sd_phan1_upper_50,target1_1.2sd_phan1_lower_50$target1)
names(phan1_1.2sd)[1:3] <- c("upper","age","lower")
phan2_1.2sd <- cbind(target1_1.2sd_phan2_upper_50,target1_1.2sd_phan2_lower_50$target1)
names(phan2_1.2sd)[1:3] <- c("upper","age","lower")
phan_1.2sd <- rbind(phan1_1.2sd,phan2_1.2sd)
phan_1.2sd <- aggregate(cbind(lower,upper) ~ age, phan_1.2sd, mean )

phan1_1.0sd <- cbind(target1_1.0sd_phan1_upper_50,target1_1.0sd_phan1_lower_50$target1)
names(phan1_1.0sd)[1:3] <- c("upper","age","lower")
phan2_1.0sd <- cbind(target1_1.0sd_phan2_upper_50,target1_1.0sd_phan2_lower_50$target1)
names(phan2_1.0sd)[1:3] <- c("upper","age","lower")
phan_1.0sd <- rbind(phan1_1.0sd,phan2_1.0sd)
phan_1.0sd <- aggregate(cbind(lower,upper) ~ age, phan_1.0sd, mean )

phan1_0.8sd <- cbind(target1_0.8sd_phan1_upper_50,target1_0.8sd_phan1_lower_50$target1)
names(phan1_0.8sd)[1:3] <- c("upper","age","lower")
phan2_0.8sd <- cbind(target1_0.8sd_phan2_upper_50,target1_0.8sd_phan2_lower_50$target1)
names(phan2_0.8sd)[1:3] <- c("upper","age","lower")
phan_0.8sd <- rbind(phan1_0.8sd,phan2_0.8sd)
phan_0.8sd <- aggregate(cbind(lower,upper) ~ age, phan_0.8sd, mean )

phan1_0.6sd <- cbind(target1_0.6sd_phan1_upper_50,target1_0.6sd_phan1_lower_50$target1)
names(phan1_0.6sd)[1:3] <- c("upper","age","lower")
phan2_0.6sd <- cbind(target1_0.6sd_phan2_upper_50,target1_0.6sd_phan2_lower_50$target1)
names(phan2_0.6sd)[1:3] <- c("upper","age","lower")
phan_0.6sd <- rbind(phan1_0.6sd,phan2_0.6sd)
phan_0.6sd <- aggregate(cbind(lower,upper) ~ age, phan_0.6sd, mean )

phan1_0.4sd <- cbind(target1_0.4sd_phan1_upper_50,target1_0.4sd_phan1_lower_50$target1)
names(phan1_0.4sd)[1:3] <- c("upper","age","lower")
phan2_0.4sd <- cbind(target1_0.4sd_phan2_upper_50,target1_0.4sd_phan2_lower_50$target1)
names(phan2_0.4sd)[1:3] <- c("upper","age","lower")
phan_0.4sd <- rbind(phan1_0.4sd,phan2_0.4sd)
phan_0.4sd <- aggregate(cbind(lower,upper) ~ age, phan_0.4sd, mean )

phan1_0.2sd <- cbind(target1_0.2sd_phan1_upper_50,target1_0.2sd_phan1_lower_50$target1)
names(phan1_0.2sd)[1:3] <- c("upper","age","lower")
phan2_0.2sd <- cbind(target1_0.2sd_phan2_upper_50,target1_0.2sd_phan2_lower_50$target1)
names(phan2_0.2sd)[1:3] <- c("upper","age","lower")
phan_0.2sd <- rbind(phan1_0.2sd,phan2_0.2sd)
phan_0.2sd <- aggregate(cbind(lower,upper) ~ age, phan_0.2sd, mean )

phan1_0sd <- cbind(target1_0sd_phan1_upper_50)
names(phan1_0sd)[1:2] <- c("upper","age")
phan2_0sd <- cbind(target1_0sd_phan2_upper_50)
names(phan2_0sd)[1:2] <- c("upper","age")
phan_0sd <- rbind(phan1_0sd,phan2_0sd)
phan_0sd <- aggregate(cbind(upper) ~ age, phan_0sd, mean )

pre_2sd <- cbind(target1_2sd_pre_upper_50,target1_2sd_pre_lower_50$target1)
names(pre_2sd)[1:3] <- c("upper","age","lower")
pre_1.8sd <- cbind(target1_1.8sd_pre_upper_50,target1_1.8sd_pre_lower_50$target1)
names(pre_1.8sd)[1:3] <- c("upper","age","lower")
pre_1.6sd <- cbind(target1_1.6sd_pre_upper_50,target1_1.6sd_pre_lower_50$target1)
names(pre_1.6sd)[1:3] <- c("upper","age","lower")
pre_1.4sd <- cbind(target1_1.4sd_pre_upper_50,target1_1.4sd_pre_lower_50$target1)
names(pre_1.4sd)[1:3] <- c("upper","age","lower")
pre_1.2sd <- cbind(target1_1.2sd_pre_upper_50,target1_1.2sd_pre_lower_50$target1)
names(pre_1.2sd)[1:3] <- c("upper","age","lower")
pre_1.0sd <- cbind(target1_1.0sd_pre_upper_50,target1_1.0sd_pre_lower_50$target1)
names(pre_1.0sd)[1:3] <- c("upper","age","lower")
pre_0.8sd <- cbind(target1_0.8sd_pre_upper_50,target1_0.8sd_pre_lower_50$target1)
names(pre_0.8sd)[1:3] <- c("upper","age","lower")
pre_0.6sd <- cbind(target1_0.6sd_pre_upper_50,target1_0.6sd_pre_lower_50$target1)
names(pre_0.6sd)[1:3] <- c("upper","age","lower")
pre_0.4sd <- cbind(target1_0.4sd_pre_upper_50,target1_0.4sd_pre_lower_50$target1)
names(pre_0.4sd)[1:3] <- c("upper","age","lower")
pre_0.2sd <- cbind(target1_0.2sd_pre_upper_50,target1_0.2sd_pre_lower_50$target1)
names(pre_0.2sd)[1:3] <- c("upper","age","lower")

a <- ggplot() +
  geom_ribbon(data = phan_2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "blue") +
  #geom_ribbon(data = phan_1.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_1.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_1.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_1.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_1.0sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "blue") +
  #geom_ribbon(data = phan_0.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_0.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_0.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_0.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_line(data = phan_0sd, aes(age,upper), colour = 'blue', size = 1) +
  geom_point(data = SGP2, aes(age, target1), size = 0.25) + 
  scale_x_reverse(limits = c(540,0), expand = c(0,0)) + theme_bw() +
  coord_cartesian(ylim = c(0, 1), expand = c(0,0)) +
  geom_hline(yintercept = c(0.7,0.8), colour = "black") +
  ylab("Fepy/FeHR")

b <- ggplot() +
  geom_ribbon(data = pre_2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "blue") +
  #geom_ribbon(data = pre_1.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_1.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_1.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_1.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_1.0sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "blue") +
  #geom_ribbon(data = pre_0.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_0.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_0.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_0.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_line(data = target1_0sd_pre_upper_50, aes(age,target1), colour = 'blue', size = 1) +
  geom_point(data = SGP2, aes(age, target1), size = 0.25) + 
  scale_x_reverse(limits = c(1800,540), expand = c(0,0)) + theme_bw() +
  coord_cartesian(ylim = c(0, 1), expand = c(0,0)) +
  geom_hline(yintercept = c(0.7,0.8), colour = "black") +
  ylab("Fepy/FeHR") +
  labs(tag = "A") +
  theme(plot.tag.position = c(0.01,0.99))

grid.arrange(b,a,ncol=2)

##### Downsample, export and plot results ####

# Phanerozoic age differences
diffs <-  SGP2$age
diffs <- as.data.frame(diffs[order(diffs)])
names(diffs)[1] <- "Age"
diffs <- as.data.frame(diffs[!duplicated(diffs$Age),])
diffs_df <- round(diffs[-1,]) - round(diffs[-nrow(diffs),])
diffs_df <- as.data.frame(diffs_df)
diffs <- cbind(diffs[1:nrow(diffs)-1,],diffs_df)
names(diffs)[1] <- "Age"
names(diffs)[2] <- "Diff"
diffs <- subset(diffs, Diff > 0)
diffs$Diff2 <- diffs$Diff/2
diffs$edge <- round(diffs$Diff2+diffs$Age)
diffs$span <- c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4])
diffs$mid <- (c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4]))/2
diffs$mid <- diffs$edge - diffs$mid
diffs <- subset(diffs, span > 0)
diffs$bins <- cut(diffs$mid, c(0, unique(diffs$edge)), include.lowest = TRUE)
diffs.phan <- subset(diffs, Age < 541)

# Precambrian
diffs <-  c(SGP2$age,2000) # dummy base age added
diffs <- as.data.frame(diffs[order(diffs)])
names(diffs)[1] <- "Age"
diffs <- as.data.frame(diffs[!duplicated(diffs$Age),])
diffs <- round(diffs, digits = -1)
diffs_df <- round(diffs[-1,], digits = -1) - round(diffs[-nrow(diffs),], digits = -1)
diffs_df <- as.data.frame(diffs_df)
diffs <- cbind(diffs[1:nrow(diffs)-1,],diffs_df)
names(diffs)[1] <- "Age"
names(diffs)[2] <- "Diff"
diffs <- subset(diffs, Diff > 0)
diffs$Diff2 <- diffs$Diff/2
diffs$edge <- round(diffs$Diff2+diffs$Age)
diffs$span <- c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4])
diffs$mid <- (c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4]))/2
diffs$mid <- diffs$edge - diffs$mid
diffs <- subset(diffs, span > 0)
diffs$bins <- cut(diffs$mid, c(540, unique(diffs$edge)), include.lowest = TRUE)
diffs.pre <- subset(diffs, Age >= 540)

# Phanerozoic
target1_phan_comp <- rbind(target1_2sd_phan1_upper_all,target1_1.8sd_phan1_upper_all,
                           target1_1.6sd_phan1_upper_all,target1_1.4sd_phan1_upper_all,
                           target1_1.2sd_phan1_upper_all,target1_1.0sd_phan1_upper_all,
                           target1_0.8sd_phan1_upper_all,target1_0.6sd_phan1_upper_all,
                           target1_0.4sd_phan1_upper_all,target1_0.2sd_phan1_upper_all,
                           target1_2sd_phan1_lower_all,target1_1.8sd_phan1_lower_all,
                           target1_1.6sd_phan1_lower_all,target1_1.4sd_phan1_lower_all,
                           target1_1.2sd_phan1_lower_all,target1_1.0sd_phan1_lower_all,
                           target1_0.8sd_phan1_lower_all,target1_0.6sd_phan1_lower_all,
                           target1_0.4sd_phan1_lower_all,target1_0.2sd_phan1_lower_all,
                           target1_0sd_phan1_upper_all,
                           target1_2sd_phan2_upper_all,target1_1.8sd_phan2_upper_all,
                           target1_1.6sd_phan2_upper_all,target1_1.4sd_phan2_upper_all,
                           target1_1.2sd_phan2_upper_all,target1_1.0sd_phan2_upper_all,
                           target1_0.8sd_phan2_upper_all,target1_0.6sd_phan2_upper_all,
                           target1_0.4sd_phan2_upper_all,target1_0.2sd_phan2_upper_all,
                           target1_2sd_phan2_lower_all,target1_1.8sd_phan2_lower_all,
                           target1_1.6sd_phan2_lower_all,target1_1.4sd_phan2_lower_all,
                           target1_1.2sd_phan2_lower_all,target1_1.0sd_phan2_lower_all,
                           target1_0.8sd_phan2_lower_all,target1_0.6sd_phan2_lower_all,
                           target1_0.4sd_phan2_lower_all,target1_0.2sd_phan2_lower_all,
                           target1_0sd_phan2_upper_all)

target1_phan_comp <- target1_phan_comp[complete.cases(target1_phan_comp[,]),] 
target1_phan_comp$target1[target1_phan_comp$target1 > 1] <- 1
target1_phan_comp$target1[target1_phan_comp$target1 < 0] <- 0

target1_phan_comp$redox <- cut(target1_phan_comp$target1, c(0,0.7,0.8,1), include.lowest = TRUE,
                              labels = c("ferruginous", "possible euxinic","euxinic"))

target1_phan_comp.bins <- target1_phan_comp
target1_phan_comp.bins$bins <- cut(target1_phan_comp.bins$age2, c(0, unique(diffs.phan$edge)), include.lowest = TRUE)

target1_phan_comp.bins <- target1_phan_comp.bins %>% 
  group_by(bins, redox) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

target1_phan_comp.bins <- target1_phan_comp.bins[,-3]
target1_phan_comp.bins <- reshape2::dcast(target1_phan_comp.bins, bins~redox)
target1_phan_comp.bins[is.na(target1_phan_comp.bins)] <- 0

phan_results.target1.bins <- join(diffs.phan, target1_phan_comp.bins, by = "bins")
phan_results.target1.bins <- phan_results.target1.bins[,c(4:6,8:ncol(phan_results.target1.bins))]
phan_results.target1.bins <- rbind(phan_results.target1.bins[1,],phan_results.target1.bins)
phan_results.target1.bins[1,3] <- NA
phan_results.target1.bins[1,1] <- 0
phan_results.target1.melt.bins <- reshape2::melt(phan_results.target1.bins, id.vars = c("edge","span", "mid"))

#write.csv(phan_results.target1.bins, "SGP_Fe_redox_fractions_downsampled_phanerozoic.csv") # optional export

# Precambrian
target1_pre_comp <- rbind(target1_2sd_pre_upper_all,target1_1.8sd_pre_upper_all,
                           target1_1.6sd_pre_upper_all,target1_1.4sd_pre_upper_all,
                           target1_1.2sd_pre_upper_all,target1_1.0sd_pre_upper_all,
                           target1_0.8sd_pre_upper_all,target1_0.6sd_pre_upper_all,
                           target1_0.4sd_pre_upper_all,target1_0.2sd_pre_upper_all,
                           target1_2sd_pre_lower_all,target1_1.8sd_pre_lower_all,
                           target1_1.6sd_pre_lower_all,target1_1.4sd_pre_lower_all,
                           target1_1.2sd_pre_lower_all,target1_1.0sd_pre_lower_all,
                           target1_0.8sd_pre_lower_all,target1_0.6sd_pre_lower_all,
                           target1_0.4sd_pre_lower_all,target1_0.2sd_pre_lower_all,
                           target1_0sd_pre_upper_all)

target1_pre_comp <- target1_pre_comp[complete.cases(target1_pre_comp[,]),] 
target1_pre_comp$target1[target1_pre_comp$target1 > 1] <- 1
target1_pre_comp$target1[target1_pre_comp$target1 < 0] <- 0

target1_pre_comp$redox <- cut(target1_pre_comp$target1, c(0,0.7,0.8,1), include.lowest = TRUE,
                               labels = c("ferruginous", "possible euxinic","euxinic"))

target1_pre_comp.bins <- target1_pre_comp
target1_pre_comp.bins$bins <- cut(target1_pre_comp.bins$age2, c(0, unique(diffs.pre$edge)), include.lowest = TRUE)

target1_pre_comp.bins <- target1_pre_comp.bins %>% 
  group_by(bins, redox) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

target1_pre_comp.bins <- target1_pre_comp.bins[,-3]
target1_pre_comp.bins <- reshape2::dcast(target1_pre_comp.bins, bins~redox)
target1_pre_comp.bins[is.na(target1_pre_comp.bins)] <- 0

pre_results.target1.bins <- join(diffs.pre, target1_pre_comp.bins, by = "bins")
pre_results.target1.bins <- pre_results.target1.bins[,c(4:6,8:ncol(pre_results.target1.bins))]
pre_results.target1.bins <- rbind(pre_results.target1.bins[1,],pre_results.target1.bins)
pre_results.target1.bins[1,3] <- NA
pre_results.target1.bins[2,4:6] <- pre_results.target1.bins[3,4:6]
pre_results.target1.bins <- pre_results.target1.bins[-1,]
pre_results.target1.melt.bins <- reshape2::melt(pre_results.target1.bins, id.vars = c("edge","span", "mid"))

#write.csv(pre_results.target1.bins, "SGP_Fe_redox_fractions_downsampled_precambrian.csv") # optional export

my_colours <- RColorBrewer::brewer.pal(11, "RdGy")[-c(1,3,4,5,6,7,10,11)]

c <- ggplot() +
  geom_bar(data = phan_results.target1.melt.bins, aes(y = value, x = mid, fill = variable, width = span),
           stat="identity") +
  theme_bw() + scale_x_reverse(limits = c(540,0)) +
  geom_step(data = phan_results.target1.melt.bins, aes(y = value, x = edge, group = variable), 
            position = "stack", direction = "hv") +
  geom_step(data = subset(phan_results.target1.melt.bins, variable == "euxinic" | variable == "possible euxinic" ), aes(y = value, x = edge, group = variable), 
            position = "stack", direction = "hv") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1), expand = c(0,0)) +
  ylab("Fraction euxinic") +
  scale_fill_manual(values = my_colours)

d <- ggplot() +
  geom_bar(data = pre_results.target1.melt.bins, aes(y = value, x = mid, fill = variable, width = span),
           stat="identity") +
  theme_bw() + scale_x_reverse(limits = c(1800,540)) +
  geom_step(data = pre_results.target1.melt.bins, aes(y = value, x = edge, group = variable), 
            position = "stack", direction = "hv") +
  geom_step(data = subset(pre_results.target1.melt.bins, variable == "euxinic" | variable == "possible euxinic" ), aes(y = value, x = edge, group = variable), 
            position = "stack", direction = "hv") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1), expand = c(0,0)) +
  ylab("Fraction euxinic") +
  scale_fill_manual(values = my_colours)

grid.arrange(b,a,d,c,ncol=2)

###### SGP - TOC/P ####

source_Reinhard <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fnature20772/MediaObjects/41586_2017_BFnature20772_MOESM300_ESM.xlsx'
Reinhard1 <- rio::import(source_Reinhard, which  = 8) # V > 300

Reinhard1$`TOC (wt%)` <- as.numeric(Reinhard1$`TOC (wt%)`)

Reinhard1 <- Reinhard1[!is.na(Reinhard1$`TOC (wt%)`),]
Reinhard1 <- Reinhard1[!is.na(Reinhard1$`P (wt%)`),]

Reinhard1$target2 <- log10(Reinhard1$`TOC (wt%)`/Reinhard1$`P (wt%)`)
Reinhard1 <- Reinhard1[!is.na(Reinhard1$target2),]
Reinhard1 <- Reinhard1[!is.infinite(Reinhard1$target2),]
names(Reinhard1)[1] <- "age"
Reinhard1 <- subset(Reinhard1, select = c(age, target2, Source))
Reinhard1$source <- "Reinhard"

source_Reinhard <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fnature20772/MediaObjects/41586_2017_BFnature20772_MOESM300_ESM.xlsx'
Reinhard2 <- rio::import(source_Reinhard, which  = 1) # all Precambrian

names(Reinhard2)[2] <- "age"
Reinhard2 <- subset(Reinhard2, age > 540)

Reinhard2$`TOC (wt%)` <- as.numeric(Reinhard2$`TOC (wt%)`)

Reinhard2 <- Reinhard2[!is.na(Reinhard2$`TOC (wt%)`),]
Reinhard2 <- Reinhard2[!is.na(Reinhard2$`P (wt%)`),]

Reinhard2$target2 <- log10(Reinhard2$`TOC (wt%)`/Reinhard2$`P (wt%)`)
Reinhard2 <- Reinhard2[!is.na(Reinhard2$target2),]
Reinhard2 <- Reinhard2[!is.infinite(Reinhard2$target2),]
Reinhard2 <- subset(Reinhard2, select = c(age, target2, Source))
Reinhard2$source <- "Reinhard"

Reinhard <- rbind(Reinhard1, Reinhard2)

# SGP 

SGP2 <- cbind(SGP[,c(1:5,21)],ppm,ppb,wt)

names(SGP2) <- gsub(x = names(SGP2), pattern = " (ppm)", replacement = "", fixed = TRUE)
names(SGP2) <- gsub(x = names(SGP2), pattern = " (ppb)", replacement = "", fixed = TRUE)
names(SGP2) <- gsub(x = names(SGP2), pattern = " (wt%)", replacement = "", fixed = TRUE)
names(SGP2) <- gsub(x = names(SGP2), pattern = "-", replacement = "", fixed = TRUE)
names(SGP2)[5] <- "age"

# define intervals based on pyrite type distributions #

zones <- c(3400,3100,3000,2800,2100,1700,1500,1400,900,800,740,690,655,610,532,514,480,456,419,383,372,367,359,330,299,280,250,201,54,0)

SGP2$zones <- cut(SGP2$age, zones)

get_midpoint <- function(cut_label) {
  mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ","))))
}

SGP2$zones2 <- sapply(SGP2$zones, get_midpoint)

# optional anoxic subset

SGP2 <- subset(SGP2, (Fecarb+Feox+Femag+Fepy)/Fe > 0.38 & Fe > (0.5*10000))

SGP2$target2 <- log10(SGP2$TOC/(SGP2$P))

SGP2 <- SGP2[!is.na(SGP2$target2),]
SGP2 <- SGP2[!is.infinite(SGP2$target2),]

SGP2 <- subset(SGP2, select = c(age, target2, `analysis reference (long)`))
SGP2$source <- "SGP"

# check for duplicates

Reinhard_names <- Reinhard
Reinhard_names$first <- word(Reinhard_names$Source, 1)

SGP_names <- SGP2
SGP_names$first <- word(SGP_names$`analysis reference (long)`, 1)
SGP_names$first <- gsub(",","",SGP_names$first)

check <- inner_join(SGP_names, Reinhard_names, by = "first")

# remove duplicates

Reinhard <- subset(Reinhard, Source != "Brumsack 1991")

# combine Reinhard and SGP TOC/P

all <- rbind(SGP2[,-3], Reinhard[,-3])
all <- all[complete.cases(all), ]

# SGP focus for regression - see main text for explanation

SGP2_phan <- subset(SGP2, age <= 540)
SGP2_phan1 <- subset(SGP2_phan, age > 150)
SGP2_phan2 <- subset(SGP2_phan, age < 200)
SGP2_pre <- subset(SGP2, age > 540)


# build normal distribution

speed <- 10 # speed = 10 for full version

# construct normal distribution for sampling

n <- 10000
norm <- rnorm(n)
out <- hist(norm, breaks = c(seq(-4.3,4.3, by = 0.2)))
bins <- as.data.frame(out[[1]])
bins <- round(bins+0.1, digits = 2)
freqs <- as.data.frame(out[[2]])
out <- as.data.frame(cbind(bins[1:nrow(bins)-1,],freqs))
names(out)[1:2] <- c("sd","count")

# phanerozoic

# span_upper <- loess.predict[["pars"]][["span"]]
# span_upper <- loess.predict[["pars"]][["span"]]

# 2sd

nsig <- 2

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_2sd_phan1_upper <- as.data.frame(cbind(target2_2sd_phan1$x, target2_2sd_phan1$upper))
names(target2_2sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_2sd_phan1_upper$age, target2_2sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_2sd_phan1_upper, indices) {
  d <- target2_2sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_2sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_2sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_2sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_2sd_phan1_upper_50)[1] <- "target2"

target2_2sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_2sd_phan1_upper_all)[1:ncol(target2_2sd_phan1_upper_all)] <- seq(0,540, 1)
target2_2sd_phan1_upper_all$group <- seq(1,(nrow(target2_2sd_phan1_upper_all)), 1)

target2_2sd_phan1_upper_all <- reshape2::melt(target2_2sd_phan1_upper_all, id.vars = "group")
names(target2_2sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_2sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_2sd_phan1_lower <- as.data.frame(cbind(target2_2sd_phan1$x, target2_2sd_phan1$lower))
names(target2_2sd_phan1_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_2sd_phan1_lower$age, target2_2sd_phan1_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_2sd_phan1_lower, indices) {
  d <- target2_2sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_2sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_2sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_2sd_phan1_lower_50$age <- seq(0,540, 1)
names(target2_2sd_phan1_lower_50)[1] <- "target2"

target2_2sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target2_2sd_phan1_lower_all)[1:ncol(target2_2sd_phan1_lower_all)] <- seq(0,540, 1)
target2_2sd_phan1_lower_all$group <- seq(1,(nrow(target2_2sd_phan1_lower_all)), 1)

target2_2sd_phan1_lower_all <- reshape2::melt(target2_2sd_phan1_lower_all, id.vars = "group")
names(target2_2sd_phan1_lower_all)[2:3] <- c("age", "target2")

target2_2sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8sd

nsig <- 1.8

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.8sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_1.8sd_phan1_upper <- as.data.frame(cbind(target2_1.8sd_phan1$x, target2_1.8sd_phan1$upper))
names(target2_1.8sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.8sd_phan1_upper$age, target2_1.8sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.8sd_phan1_upper, indices) {
  d <- target2_1.8sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.8sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.8sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.8sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_1.8sd_phan1_upper_50)[1] <- "target2"

target2_1.8sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.8sd_phan1_upper_all)[1:ncol(target2_1.8sd_phan1_upper_all)] <- seq(0,540, 1)
target2_1.8sd_phan1_upper_all$group <- seq(1,(nrow(target2_1.8sd_phan1_upper_all)), 1)

target2_1.8sd_phan1_upper_all <- reshape2::melt(target2_1.8sd_phan1_upper_all, id.vars = "group")
names(target2_1.8sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_1.8sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.8sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_1.8sd_phan1_lower <- as.data.frame(cbind(target2_1.8sd_phan1$x, target2_1.8sd_phan1$lower))
names(target2_1.8sd_phan1_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.8sd_phan1_lower$age, target2_1.8sd_phan1_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.8sd_phan1_lower, indices) {
  d <- target2_1.8sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.8sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.8sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.8sd_phan1_lower_50$age <- seq(0,540, 1)
names(target2_1.8sd_phan1_lower_50)[1] <- "target2"

target2_1.8sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.8sd_phan1_lower_all)[1:ncol(target2_1.8sd_phan1_lower_all)] <- seq(0,540, 1)
target2_1.8sd_phan1_lower_all$group <- seq(1,(nrow(target2_1.8sd_phan1_lower_all)), 1)

target2_1.8sd_phan1_lower_all <- reshape2::melt(target2_1.8sd_phan1_lower_all, id.vars = "group")
names(target2_1.8sd_phan1_lower_all)[2:3] <- c("age", "target2")

target2_1.8sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6sd

nsig <- 1.6

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.6sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_1.6sd_phan1_upper <- as.data.frame(cbind(target2_1.6sd_phan1$x, target2_1.6sd_phan1$upper))
names(target2_1.6sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.6sd_phan1_upper$age, target2_1.6sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.6sd_phan1_upper, indices) {
  d <- target2_1.6sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.6sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.6sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.6sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_1.6sd_phan1_upper_50)[1] <- "target2"

target2_1.6sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.6sd_phan1_upper_all)[1:ncol(target2_1.6sd_phan1_upper_all)] <- seq(0,540, 1)
target2_1.6sd_phan1_upper_all$group <- seq(1,(nrow(target2_1.6sd_phan1_upper_all)), 1)

target2_1.6sd_phan1_upper_all <- reshape2::melt(target2_1.6sd_phan1_upper_all, id.vars = "group")
names(target2_1.6sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_1.6sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.6sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_1.6sd_phan1_lower <- as.data.frame(cbind(target2_1.6sd_phan1$x, target2_1.6sd_phan1$lower))
names(target2_1.6sd_phan1_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.6sd_phan1_lower$age, target2_1.6sd_phan1_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.6sd_phan1_lower, indices) {
  d <- target2_1.6sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.6sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.6sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.6sd_phan1_lower_50$age <- seq(0,540, 1)
names(target2_1.6sd_phan1_lower_50)[1] <- "target2"

target2_1.6sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.6sd_phan1_lower_all)[1:ncol(target2_1.6sd_phan1_lower_all)] <- seq(0,540, 1)
target2_1.6sd_phan1_lower_all$group <- seq(1,(nrow(target2_1.6sd_phan1_lower_all)), 1)

target2_1.6sd_phan1_lower_all <- reshape2::melt(target2_1.6sd_phan1_lower_all, id.vars = "group")
names(target2_1.6sd_phan1_lower_all)[2:3] <- c("age", "target2")

target2_1.6sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4sd

nsig <- 1.4

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.4sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_1.4sd_phan1_upper <- as.data.frame(cbind(target2_1.4sd_phan1$x, target2_1.4sd_phan1$upper))
names(target2_1.4sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.4sd_phan1_upper$age, target2_1.4sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.4sd_phan1_upper, indices) {
  d <- target2_1.4sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.4sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.4sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.4sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_1.4sd_phan1_upper_50)[1] <- "target2"

target2_1.4sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.4sd_phan1_upper_all)[1:ncol(target2_1.4sd_phan1_upper_all)] <- seq(0,540, 1)
target2_1.4sd_phan1_upper_all$group <- seq(1,(nrow(target2_1.4sd_phan1_upper_all)), 1)

target2_1.4sd_phan1_upper_all <- reshape2::melt(target2_1.4sd_phan1_upper_all, id.vars = "group")
names(target2_1.4sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_1.4sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.4sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_1.4sd_phan1_lower <- as.data.frame(cbind(target2_1.4sd_phan1$x, target2_1.4sd_phan1$lower))
names(target2_1.4sd_phan1_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.4sd_phan1_lower$age, target2_1.4sd_phan1_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.4sd_phan1_lower, indices) {
  d <- target2_1.4sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.4sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.4sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.4sd_phan1_lower_50$age <- seq(0,540, 1)
names(target2_1.4sd_phan1_lower_50)[1] <- "target2"

target2_1.4sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.4sd_phan1_lower_all)[1:ncol(target2_1.4sd_phan1_lower_all)] <- seq(0,540, 1)
target2_1.4sd_phan1_lower_all$group <- seq(1,(nrow(target2_1.4sd_phan1_lower_all)), 1)

target2_1.4sd_phan1_lower_all <- reshape2::melt(target2_1.4sd_phan1_lower_all, id.vars = "group")
names(target2_1.4sd_phan1_lower_all)[2:3] <- c("age", "target2")

target2_1.4sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2sd

nsig <- 1.2

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_1.2sd_phan1_upper <- as.data.frame(cbind(target2_1.2sd_phan1$x, target2_1.2sd_phan1$upper))
names(target2_1.2sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.2sd_phan1_upper$age, target2_1.2sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.2sd_phan1_upper, indices) {
  d <- target2_1.2sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.2sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.2sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.2sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_1.2sd_phan1_upper_50)[1] <- "target2"

target2_1.2sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.2sd_phan1_upper_all)[1:ncol(target2_1.2sd_phan1_upper_all)] <- seq(0,540, 1)
target2_1.2sd_phan1_upper_all$group <- seq(1,(nrow(target2_1.2sd_phan1_upper_all)), 1)

target2_1.2sd_phan1_upper_all <- reshape2::melt(target2_1.2sd_phan1_upper_all, id.vars = "group")
names(target2_1.2sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_1.2sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_1.2sd_phan1_lower <- as.data.frame(cbind(target2_1.2sd_phan1$x, target2_1.2sd_phan1$lower))
names(target2_1.2sd_phan1_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.2sd_phan1_lower$age, target2_1.2sd_phan1_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.2sd_phan1_lower, indices) {
  d <- target2_1.2sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.2sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.2sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.2sd_phan1_lower_50$age <- seq(0,540, 1)
names(target2_1.2sd_phan1_lower_50)[1] <- "target2"

target2_1.2sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.2sd_phan1_lower_all)[1:ncol(target2_1.2sd_phan1_lower_all)] <- seq(0,540, 1)
target2_1.2sd_phan1_lower_all$group <- seq(1,(nrow(target2_1.2sd_phan1_lower_all)), 1)

target2_1.2sd_phan1_lower_all <- reshape2::melt(target2_1.2sd_phan1_lower_all, id.vars = "group")
names(target2_1.2sd_phan1_lower_all)[2:3] <- c("age", "target2")

target2_1.2sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0sd

nsig <- 1.0

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.0sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_1.0sd_phan1_upper <- as.data.frame(cbind(target2_1.0sd_phan1$x, target2_1.0sd_phan1$upper))
names(target2_1.0sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.0sd_phan1_upper$age, target2_1.0sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.0sd_phan1_upper, indices) {
  d <- target2_1.0sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.0sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.0sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.0sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_1.0sd_phan1_upper_50)[1] <- "target2"

target2_1.0sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.0sd_phan1_upper_all)[1:ncol(target2_1.0sd_phan1_upper_all)] <- seq(0,540, 1)
target2_1.0sd_phan1_upper_all$group <- seq(1,(nrow(target2_1.0sd_phan1_upper_all)), 1)

target2_1.0sd_phan1_upper_all <- reshape2::melt(target2_1.0sd_phan1_upper_all, id.vars = "group")
names(target2_1.0sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_1.0sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.0sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_1.0sd_phan1_lower <- as.data.frame(cbind(target2_1.0sd_phan1$x, target2_1.0sd_phan1$lower))
names(target2_1.0sd_phan1_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.0sd_phan1_lower$age, target2_1.0sd_phan1_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.0sd_phan1_lower, indices) {
  d <- target2_1.0sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.0sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.0sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.0sd_phan1_lower_50$age <- seq(0,540, 1)
names(target2_1.0sd_phan1_lower_50)[1] <- "target2"

target2_1.0sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.0sd_phan1_lower_all)[1:ncol(target2_1.0sd_phan1_lower_all)] <- seq(0,540, 1)
target2_1.0sd_phan1_lower_all$group <- seq(1,(nrow(target2_1.0sd_phan1_lower_all)), 1)

target2_1.0sd_phan1_lower_all <- reshape2::melt(target2_1.0sd_phan1_lower_all, id.vars = "group")
names(target2_1.0sd_phan1_lower_all)[2:3] <- c("age", "target2")

target2_1.0sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8sd

nsig <- 0.8

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.8sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_0.8sd_phan1_upper <- as.data.frame(cbind(target2_0.8sd_phan1$x, target2_0.8sd_phan1$upper))
names(target2_0.8sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.8sd_phan1_upper$age, target2_0.8sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.8sd_phan1_upper, indices) {
  d <- target2_0.8sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.8sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.8sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.8sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_0.8sd_phan1_upper_50)[1] <- "target2"

target2_0.8sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.8sd_phan1_upper_all)[1:ncol(target2_0.8sd_phan1_upper_all)] <- seq(0,540, 1)
target2_0.8sd_phan1_upper_all$group <- seq(1,(nrow(target2_0.8sd_phan1_upper_all)), 1)

target2_0.8sd_phan1_upper_all <- reshape2::melt(target2_0.8sd_phan1_upper_all, id.vars = "group")
names(target2_0.8sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_0.8sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.8sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_0.8sd_phan1_lower <- as.data.frame(cbind(target2_0.8sd_phan1$x, target2_0.8sd_phan1$lower))
names(target2_0.8sd_phan1_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.8sd_phan1_lower$age, target2_0.8sd_phan1_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.8sd_phan1_lower, indices) {
  d <- target2_0.8sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.8sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.8sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.8sd_phan1_lower_50$age <- seq(0,540, 1)
names(target2_0.8sd_phan1_lower_50)[1] <- "target2"

target2_0.8sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.8sd_phan1_lower_all)[1:ncol(target2_0.8sd_phan1_lower_all)] <- seq(0,540, 1)
target2_0.8sd_phan1_lower_all$group <- seq(1,(nrow(target2_0.8sd_phan1_lower_all)), 1)

target2_0.8sd_phan1_lower_all <- reshape2::melt(target2_0.8sd_phan1_lower_all, id.vars = "group")
names(target2_0.8sd_phan1_lower_all)[2:3] <- c("age", "target2")

target2_0.8sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6sd

nsig <- 0.6

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.6sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_0.6sd_phan1_upper <- as.data.frame(cbind(target2_0.6sd_phan1$x, target2_0.6sd_phan1$upper))
names(target2_0.6sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.6sd_phan1_upper$age, target2_0.6sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.6sd_phan1_upper, indices) {
  d <- target2_0.6sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.6sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.6sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.6sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_0.6sd_phan1_upper_50)[1] <- "target2"

target2_0.6sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.6sd_phan1_upper_all)[1:ncol(target2_0.6sd_phan1_upper_all)] <- seq(0,540, 1)
target2_0.6sd_phan1_upper_all$group <- seq(1,(nrow(target2_0.6sd_phan1_upper_all)), 1)

target2_0.6sd_phan1_upper_all <- reshape2::melt(target2_0.6sd_phan1_upper_all, id.vars = "group")
names(target2_0.6sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_0.6sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.6sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_0.6sd_phan1_lower <- as.data.frame(cbind(target2_0.6sd_phan1$x, target2_0.6sd_phan1$lower))
names(target2_0.6sd_phan1_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.6sd_phan1_lower$age, target2_0.6sd_phan1_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.6sd_phan1_lower, indices) {
  d <- target2_0.6sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.6sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.6sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.6sd_phan1_lower_50$age <- seq(0,540, 1)
names(target2_0.6sd_phan1_lower_50)[1] <- "target2"

target2_0.6sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.6sd_phan1_lower_all)[1:ncol(target2_0.6sd_phan1_lower_all)] <- seq(0,540, 1)
target2_0.6sd_phan1_lower_all$group <- seq(1,(nrow(target2_0.6sd_phan1_lower_all)), 1)

target2_0.6sd_phan1_lower_all <- reshape2::melt(target2_0.6sd_phan1_lower_all, id.vars = "group")
names(target2_0.6sd_phan1_lower_all)[2:3] <- c("age", "target2")

target2_0.6sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4sd

nsig <- 0.4

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.4sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_0.4sd_phan1_upper <- as.data.frame(cbind(target2_0.4sd_phan1$x, target2_0.4sd_phan1$upper))
names(target2_0.4sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.4sd_phan1_upper$age, target2_0.4sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.4sd_phan1_upper, indices) {
  d <- target2_0.4sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.4sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.4sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.4sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_0.4sd_phan1_upper_50)[1] <- "target2"

target2_0.4sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.4sd_phan1_upper_all)[1:ncol(target2_0.4sd_phan1_upper_all)] <- seq(0,540, 1)
target2_0.4sd_phan1_upper_all$group <- seq(1,(nrow(target2_0.4sd_phan1_upper_all)), 1)

target2_0.4sd_phan1_upper_all <- reshape2::melt(target2_0.4sd_phan1_upper_all, id.vars = "group")
names(target2_0.4sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_0.4sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.4sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_0.4sd_phan1_lower <- as.data.frame(cbind(target2_0.4sd_phan1$x, target2_0.4sd_phan1$lower))
names(target2_0.4sd_phan1_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.4sd_phan1_lower$age, target2_0.4sd_phan1_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.4sd_phan1_lower, indices) {
  d <- target2_0.4sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.4sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.4sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.4sd_phan1_lower_50$age <- seq(0,540, 1)
names(target2_0.4sd_phan1_lower_50)[1] <- "target2"

target2_0.4sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.4sd_phan1_lower_all)[1:ncol(target2_0.4sd_phan1_lower_all)] <- seq(0,540, 1)
target2_0.4sd_phan1_lower_all$group <- seq(1,(nrow(target2_0.4sd_phan1_lower_all)), 1)

target2_0.4sd_phan1_lower_all <- reshape2::melt(target2_0.4sd_phan1_lower_all, id.vars = "group")
names(target2_0.4sd_phan1_lower_all)[2:3] <- c("age", "target2")

target2_0.4sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2sd

nsig <- 0.2

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_0.2sd_phan1_upper <- as.data.frame(cbind(target2_0.2sd_phan1$x, target2_0.2sd_phan1$upper))
names(target2_0.2sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.2sd_phan1_upper$age, target2_0.2sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.2sd_phan1_upper, indices) {
  d <- target2_0.2sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.2sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.2sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.2sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_0.2sd_phan1_upper_50)[1] <- "target2"

target2_0.2sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.2sd_phan1_upper_all)[1:ncol(target2_0.2sd_phan1_upper_all)] <- seq(0,540, 1)
target2_0.2sd_phan1_upper_all$group <- seq(1,(nrow(target2_0.2sd_phan1_upper_all)), 1)

target2_0.2sd_phan1_upper_all <- reshape2::melt(target2_0.2sd_phan1_upper_all, id.vars = "group")
names(target2_0.2sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_0.2sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.2sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_0.2sd_phan1_lower <- as.data.frame(cbind(target2_0.2sd_phan1$x, target2_0.2sd_phan1$lower))
names(target2_0.2sd_phan1_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.2sd_phan1_lower$age, target2_0.2sd_phan1_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.2sd_phan1_lower, indices) {
  d <- target2_0.2sd_phan1_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.2sd_phan1_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.2sd_phan1_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.2sd_phan1_lower_50$age <- seq(0,540, 1)
names(target2_0.2sd_phan1_lower_50)[1] <- "target2"

target2_0.2sd_phan1_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.2sd_phan1_lower_all)[1:ncol(target2_0.2sd_phan1_lower_all)] <- seq(0,540, 1)
target2_0.2sd_phan1_lower_all$group <- seq(1,(nrow(target2_0.2sd_phan1_lower_all)), 1)

target2_0.2sd_phan1_lower_all <- reshape2::melt(target2_0.2sd_phan1_lower_all, id.vars = "group")
names(target2_0.2sd_phan1_lower_all)[2:3] <- c("age", "target2")

target2_0.2sd_phan1_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0sd (once)

nsig <- 0

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan1$age, SGP2_phan1$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0sd_phan1 <- msir::loess.sd(SGP2_phan1$age, SGP2_phan1$target2, degree=0, nsigma = nsig, span = span)
target2_0sd_phan1_upper <- as.data.frame(cbind(target2_0sd_phan1$x, target2_0sd_phan1$upper))
names(target2_0sd_phan1_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0sd_phan1_upper$age, target2_0sd_phan1_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0sd_phan1_upper, indices) {
  d <- target2_0sd_phan1_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0sd_phan1_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0sd_phan1_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0sd_phan1_upper_50$age <- seq(0,540, 1)
names(target2_0sd_phan1_upper_50)[1] <- "target2"

target2_0sd_phan1_upper_all <- as.data.frame(loess_boot$t)
names(target2_0sd_phan1_upper_all)[1:ncol(target2_0sd_phan1_upper_all)] <- seq(0,540, 1)
target2_0sd_phan1_upper_all$group <- seq(1,(nrow(target2_0sd_phan1_upper_all)), 1)

target2_0sd_phan1_upper_all <- reshape2::melt(target2_0sd_phan1_upper_all, id.vars = "group")
names(target2_0sd_phan1_upper_all)[2:3] <- c("age", "target2")

target2_0sd_phan1_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# Mesozoic-Cenozoic

# phanerozoic

# span_upper <- loess.predict[["pars"]][["span"]]
# span_upper <- loess.predict[["pars"]][["span"]]

# 2sd

nsig <- 2

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_2sd_phan2_upper <- as.data.frame(cbind(target2_2sd_phan2$x, target2_2sd_phan2$upper))
names(target2_2sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_2sd_phan2_upper$age, target2_2sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_2sd_phan2_upper, indices) {
  d <- target2_2sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_2sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_2sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_2sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_2sd_phan2_upper_50)[1] <- "target2"

target2_2sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_2sd_phan2_upper_all)[1:ncol(target2_2sd_phan2_upper_all)] <- seq(0,540, 1)
target2_2sd_phan2_upper_all$group <- seq(1,(nrow(target2_2sd_phan2_upper_all)), 1)

target2_2sd_phan2_upper_all <- reshape2::melt(target2_2sd_phan2_upper_all, id.vars = "group")
names(target2_2sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_2sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_2sd_phan2_lower <- as.data.frame(cbind(target2_2sd_phan2$x, target2_2sd_phan2$lower))
names(target2_2sd_phan2_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_2sd_phan2_lower$age, target2_2sd_phan2_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_2sd_phan2_lower, indices) {
  d <- target2_2sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_2sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_2sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_2sd_phan2_lower_50$age <- seq(0,540, 1)
names(target2_2sd_phan2_lower_50)[1] <- "target2"

target2_2sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target2_2sd_phan2_lower_all)[1:ncol(target2_2sd_phan2_lower_all)] <- seq(0,540, 1)
target2_2sd_phan2_lower_all$group <- seq(1,(nrow(target2_2sd_phan2_lower_all)), 1)

target2_2sd_phan2_lower_all <- reshape2::melt(target2_2sd_phan2_lower_all, id.vars = "group")
names(target2_2sd_phan2_lower_all)[2:3] <- c("age", "target2")

target2_2sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.8sd

nsig <- 1.8

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.8sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_1.8sd_phan2_upper <- as.data.frame(cbind(target2_1.8sd_phan2$x, target2_1.8sd_phan2$upper))
names(target2_1.8sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.8sd_phan2_upper$age, target2_1.8sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.8sd_phan2_upper, indices) {
  d <- target2_1.8sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.8sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.8sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.8sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_1.8sd_phan2_upper_50)[1] <- "target2"

target2_1.8sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.8sd_phan2_upper_all)[1:ncol(target2_1.8sd_phan2_upper_all)] <- seq(0,540, 1)
target2_1.8sd_phan2_upper_all$group <- seq(1,(nrow(target2_1.8sd_phan2_upper_all)), 1)

target2_1.8sd_phan2_upper_all <- reshape2::melt(target2_1.8sd_phan2_upper_all, id.vars = "group")
names(target2_1.8sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_1.8sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.8sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_1.8sd_phan2_lower <- as.data.frame(cbind(target2_1.8sd_phan2$x, target2_1.8sd_phan2$lower))
names(target2_1.8sd_phan2_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.8sd_phan2_lower$age, target2_1.8sd_phan2_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.8sd_phan2_lower, indices) {
  d <- target2_1.8sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.8sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.8sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.8sd_phan2_lower_50$age <- seq(0,540, 1)
names(target2_1.8sd_phan2_lower_50)[1] <- "target2"

target2_1.8sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.8sd_phan2_lower_all)[1:ncol(target2_1.8sd_phan2_lower_all)] <- seq(0,540, 1)
target2_1.8sd_phan2_lower_all$group <- seq(1,(nrow(target2_1.8sd_phan2_lower_all)), 1)

target2_1.8sd_phan2_lower_all <- reshape2::melt(target2_1.8sd_phan2_lower_all, id.vars = "group")
names(target2_1.8sd_phan2_lower_all)[2:3] <- c("age", "target2")

target2_1.8sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.6sd

nsig <- 1.6

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.6sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_1.6sd_phan2_upper <- as.data.frame(cbind(target2_1.6sd_phan2$x, target2_1.6sd_phan2$upper))
names(target2_1.6sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.6sd_phan2_upper$age, target2_1.6sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.6sd_phan2_upper, indices) {
  d <- target2_1.6sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.6sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.6sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.6sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_1.6sd_phan2_upper_50)[1] <- "target2"

target2_1.6sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.6sd_phan2_upper_all)[1:ncol(target2_1.6sd_phan2_upper_all)] <- seq(0,540, 1)
target2_1.6sd_phan2_upper_all$group <- seq(1,(nrow(target2_1.6sd_phan2_upper_all)), 1)

target2_1.6sd_phan2_upper_all <- reshape2::melt(target2_1.6sd_phan2_upper_all, id.vars = "group")
names(target2_1.6sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_1.6sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.6sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_1.6sd_phan2_lower <- as.data.frame(cbind(target2_1.6sd_phan2$x, target2_1.6sd_phan2$lower))
names(target2_1.6sd_phan2_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.6sd_phan2_lower$age, target2_1.6sd_phan2_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.6sd_phan2_lower, indices) {
  d <- target2_1.6sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.6sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.6sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.6sd_phan2_lower_50$age <- seq(0,540, 1)
names(target2_1.6sd_phan2_lower_50)[1] <- "target2"

target2_1.6sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.6sd_phan2_lower_all)[1:ncol(target2_1.6sd_phan2_lower_all)] <- seq(0,540, 1)
target2_1.6sd_phan2_lower_all$group <- seq(1,(nrow(target2_1.6sd_phan2_lower_all)), 1)

target2_1.6sd_phan2_lower_all <- reshape2::melt(target2_1.6sd_phan2_lower_all, id.vars = "group")
names(target2_1.6sd_phan2_lower_all)[2:3] <- c("age", "target2")

target2_1.6sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.4sd

nsig <- 1.4

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.4sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_1.4sd_phan2_upper <- as.data.frame(cbind(target2_1.4sd_phan2$x, target2_1.4sd_phan2$upper))
names(target2_1.4sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.4sd_phan2_upper$age, target2_1.4sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.4sd_phan2_upper, indices) {
  d <- target2_1.4sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.4sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.4sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.4sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_1.4sd_phan2_upper_50)[1] <- "target2"

target2_1.4sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.4sd_phan2_upper_all)[1:ncol(target2_1.4sd_phan2_upper_all)] <- seq(0,540, 1)
target2_1.4sd_phan2_upper_all$group <- seq(1,(nrow(target2_1.4sd_phan2_upper_all)), 1)

target2_1.4sd_phan2_upper_all <- reshape2::melt(target2_1.4sd_phan2_upper_all, id.vars = "group")
names(target2_1.4sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_1.4sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.4sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_1.4sd_phan2_lower <- as.data.frame(cbind(target2_1.4sd_phan2$x, target2_1.4sd_phan2$lower))
names(target2_1.4sd_phan2_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.4sd_phan2_lower$age, target2_1.4sd_phan2_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.4sd_phan2_lower, indices) {
  d <- target2_1.4sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.4sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.4sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.4sd_phan2_lower_50$age <- seq(0,540, 1)
names(target2_1.4sd_phan2_lower_50)[1] <- "target2"

target2_1.4sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.4sd_phan2_lower_all)[1:ncol(target2_1.4sd_phan2_lower_all)] <- seq(0,540, 1)
target2_1.4sd_phan2_lower_all$group <- seq(1,(nrow(target2_1.4sd_phan2_lower_all)), 1)

target2_1.4sd_phan2_lower_all <- reshape2::melt(target2_1.4sd_phan2_lower_all, id.vars = "group")
names(target2_1.4sd_phan2_lower_all)[2:3] <- c("age", "target2")

target2_1.4sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.2sd

nsig <- 1.2

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_1.2sd_phan2_upper <- as.data.frame(cbind(target2_1.2sd_phan2$x, target2_1.2sd_phan2$upper))
names(target2_1.2sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.2sd_phan2_upper$age, target2_1.2sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.2sd_phan2_upper, indices) {
  d <- target2_1.2sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.2sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.2sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.2sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_1.2sd_phan2_upper_50)[1] <- "target2"

target2_1.2sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.2sd_phan2_upper_all)[1:ncol(target2_1.2sd_phan2_upper_all)] <- seq(0,540, 1)
target2_1.2sd_phan2_upper_all$group <- seq(1,(nrow(target2_1.2sd_phan2_upper_all)), 1)

target2_1.2sd_phan2_upper_all <- reshape2::melt(target2_1.2sd_phan2_upper_all, id.vars = "group")
names(target2_1.2sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_1.2sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_1.2sd_phan2_lower <- as.data.frame(cbind(target2_1.2sd_phan2$x, target2_1.2sd_phan2$lower))
names(target2_1.2sd_phan2_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.2sd_phan2_lower$age, target2_1.2sd_phan2_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.2sd_phan2_lower, indices) {
  d <- target2_1.2sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.2sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.2sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.2sd_phan2_lower_50$age <- seq(0,540, 1)
names(target2_1.2sd_phan2_lower_50)[1] <- "target2"

target2_1.2sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.2sd_phan2_lower_all)[1:ncol(target2_1.2sd_phan2_lower_all)] <- seq(0,540, 1)
target2_1.2sd_phan2_lower_all$group <- seq(1,(nrow(target2_1.2sd_phan2_lower_all)), 1)

target2_1.2sd_phan2_lower_all <- reshape2::melt(target2_1.2sd_phan2_lower_all, id.vars = "group")
names(target2_1.2sd_phan2_lower_all)[2:3] <- c("age", "target2")

target2_1.2sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 1.0sd

nsig <- 1.0

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.0sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_1.0sd_phan2_upper <- as.data.frame(cbind(target2_1.0sd_phan2$x, target2_1.0sd_phan2$upper))
names(target2_1.0sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.0sd_phan2_upper$age, target2_1.0sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.0sd_phan2_upper, indices) {
  d <- target2_1.0sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.0sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.0sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.0sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_1.0sd_phan2_upper_50)[1] <- "target2"

target2_1.0sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.0sd_phan2_upper_all)[1:ncol(target2_1.0sd_phan2_upper_all)] <- seq(0,540, 1)
target2_1.0sd_phan2_upper_all$group <- seq(1,(nrow(target2_1.0sd_phan2_upper_all)), 1)

target2_1.0sd_phan2_upper_all <- reshape2::melt(target2_1.0sd_phan2_upper_all, id.vars = "group")
names(target2_1.0sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_1.0sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.0sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_1.0sd_phan2_lower <- as.data.frame(cbind(target2_1.0sd_phan2$x, target2_1.0sd_phan2$lower))
names(target2_1.0sd_phan2_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.0sd_phan2_lower$age, target2_1.0sd_phan2_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.0sd_phan2_lower, indices) {
  d <- target2_1.0sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.0sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.0sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.0sd_phan2_lower_50$age <- seq(0,540, 1)
names(target2_1.0sd_phan2_lower_50)[1] <- "target2"

target2_1.0sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.0sd_phan2_lower_all)[1:ncol(target2_1.0sd_phan2_lower_all)] <- seq(0,540, 1)
target2_1.0sd_phan2_lower_all$group <- seq(1,(nrow(target2_1.0sd_phan2_lower_all)), 1)

target2_1.0sd_phan2_lower_all <- reshape2::melt(target2_1.0sd_phan2_lower_all, id.vars = "group")
names(target2_1.0sd_phan2_lower_all)[2:3] <- c("age", "target2")

target2_1.0sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.8sd

nsig <- 0.8

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.8sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_0.8sd_phan2_upper <- as.data.frame(cbind(target2_0.8sd_phan2$x, target2_0.8sd_phan2$upper))
names(target2_0.8sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.8sd_phan2_upper$age, target2_0.8sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.8sd_phan2_upper, indices) {
  d <- target2_0.8sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.8sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.8sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.8sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_0.8sd_phan2_upper_50)[1] <- "target2"

target2_0.8sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.8sd_phan2_upper_all)[1:ncol(target2_0.8sd_phan2_upper_all)] <- seq(0,540, 1)
target2_0.8sd_phan2_upper_all$group <- seq(1,(nrow(target2_0.8sd_phan2_upper_all)), 1)

target2_0.8sd_phan2_upper_all <- reshape2::melt(target2_0.8sd_phan2_upper_all, id.vars = "group")
names(target2_0.8sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_0.8sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.8sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_0.8sd_phan2_lower <- as.data.frame(cbind(target2_0.8sd_phan2$x, target2_0.8sd_phan2$lower))
names(target2_0.8sd_phan2_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.8sd_phan2_lower$age, target2_0.8sd_phan2_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.8sd_phan2_lower, indices) {
  d <- target2_0.8sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.8sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.8sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.8sd_phan2_lower_50$age <- seq(0,540, 1)
names(target2_0.8sd_phan2_lower_50)[1] <- "target2"

target2_0.8sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.8sd_phan2_lower_all)[1:ncol(target2_0.8sd_phan2_lower_all)] <- seq(0,540, 1)
target2_0.8sd_phan2_lower_all$group <- seq(1,(nrow(target2_0.8sd_phan2_lower_all)), 1)

target2_0.8sd_phan2_lower_all <- reshape2::melt(target2_0.8sd_phan2_lower_all, id.vars = "group")
names(target2_0.8sd_phan2_lower_all)[2:3] <- c("age", "target2")

target2_0.8sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.6sd

nsig <- 0.6

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.6sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_0.6sd_phan2_upper <- as.data.frame(cbind(target2_0.6sd_phan2$x, target2_0.6sd_phan2$upper))
names(target2_0.6sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.6sd_phan2_upper$age, target2_0.6sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.6sd_phan2_upper, indices) {
  d <- target2_0.6sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.6sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.6sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.6sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_0.6sd_phan2_upper_50)[1] <- "target2"

target2_0.6sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.6sd_phan2_upper_all)[1:ncol(target2_0.6sd_phan2_upper_all)] <- seq(0,540, 1)
target2_0.6sd_phan2_upper_all$group <- seq(1,(nrow(target2_0.6sd_phan2_upper_all)), 1)

target2_0.6sd_phan2_upper_all <- reshape2::melt(target2_0.6sd_phan2_upper_all, id.vars = "group")
names(target2_0.6sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_0.6sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.6sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_0.6sd_phan2_lower <- as.data.frame(cbind(target2_0.6sd_phan2$x, target2_0.6sd_phan2$lower))
names(target2_0.6sd_phan2_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.6sd_phan2_lower$age, target2_0.6sd_phan2_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.6sd_phan2_lower, indices) {
  d <- target2_0.6sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.6sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.6sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.6sd_phan2_lower_50$age <- seq(0,540, 1)
names(target2_0.6sd_phan2_lower_50)[1] <- "target2"

target2_0.6sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.6sd_phan2_lower_all)[1:ncol(target2_0.6sd_phan2_lower_all)] <- seq(0,540, 1)
target2_0.6sd_phan2_lower_all$group <- seq(1,(nrow(target2_0.6sd_phan2_lower_all)), 1)

target2_0.6sd_phan2_lower_all <- reshape2::melt(target2_0.6sd_phan2_lower_all, id.vars = "group")
names(target2_0.6sd_phan2_lower_all)[2:3] <- c("age", "target2")

target2_0.6sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.4sd

nsig <- 0.4

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.4sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_0.4sd_phan2_upper <- as.data.frame(cbind(target2_0.4sd_phan2$x, target2_0.4sd_phan2$upper))
names(target2_0.4sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.4sd_phan2_upper$age, target2_0.4sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.4sd_phan2_upper, indices) {
  d <- target2_0.4sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.4sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.4sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.4sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_0.4sd_phan2_upper_50)[1] <- "target2"

target2_0.4sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.4sd_phan2_upper_all)[1:ncol(target2_0.4sd_phan2_upper_all)] <- seq(0,540, 1)
target2_0.4sd_phan2_upper_all$group <- seq(1,(nrow(target2_0.4sd_phan2_upper_all)), 1)

target2_0.4sd_phan2_upper_all <- reshape2::melt(target2_0.4sd_phan2_upper_all, id.vars = "group")
names(target2_0.4sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_0.4sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.4sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_0.4sd_phan2_lower <- as.data.frame(cbind(target2_0.4sd_phan2$x, target2_0.4sd_phan2$lower))
names(target2_0.4sd_phan2_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.4sd_phan2_lower$age, target2_0.4sd_phan2_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.4sd_phan2_lower, indices) {
  d <- target2_0.4sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.4sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.4sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.4sd_phan2_lower_50$age <- seq(0,540, 1)
names(target2_0.4sd_phan2_lower_50)[1] <- "target2"

target2_0.4sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.4sd_phan2_lower_all)[1:ncol(target2_0.4sd_phan2_lower_all)] <- seq(0,540, 1)
target2_0.4sd_phan2_lower_all$group <- seq(1,(nrow(target2_0.4sd_phan2_lower_all)), 1)

target2_0.4sd_phan2_lower_all <- reshape2::melt(target2_0.4sd_phan2_lower_all, id.vars = "group")
names(target2_0.4sd_phan2_lower_all)[2:3] <- c("age", "target2")

target2_0.4sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0.2sd

nsig <- 0.2

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_0.2sd_phan2_upper <- as.data.frame(cbind(target2_0.2sd_phan2$x, target2_0.2sd_phan2$upper))
names(target2_0.2sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.2sd_phan2_upper$age, target2_0.2sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.2sd_phan2_upper, indices) {
  d <- target2_0.2sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.2sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.2sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.2sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_0.2sd_phan2_upper_50)[1] <- "target2"

target2_0.2sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.2sd_phan2_upper_all)[1:ncol(target2_0.2sd_phan2_upper_all)] <- seq(0,540, 1)
target2_0.2sd_phan2_upper_all$group <- seq(1,(nrow(target2_0.2sd_phan2_upper_all)), 1)

target2_0.2sd_phan2_upper_all <- reshape2::melt(target2_0.2sd_phan2_upper_all, id.vars = "group")
names(target2_0.2sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_0.2sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.2sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_0.2sd_phan2_lower <- as.data.frame(cbind(target2_0.2sd_phan2$x, target2_0.2sd_phan2$lower))
names(target2_0.2sd_phan2_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.2sd_phan2_lower$age, target2_0.2sd_phan2_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.2sd_phan2_lower, indices) {
  d <- target2_0.2sd_phan2_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.2sd_phan2_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.2sd_phan2_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.2sd_phan2_lower_50$age <- seq(0,540, 1)
names(target2_0.2sd_phan2_lower_50)[1] <- "target2"

target2_0.2sd_phan2_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.2sd_phan2_lower_all)[1:ncol(target2_0.2sd_phan2_lower_all)] <- seq(0,540, 1)
target2_0.2sd_phan2_lower_all$group <- seq(1,(nrow(target2_0.2sd_phan2_lower_all)), 1)

target2_0.2sd_phan2_lower_all <- reshape2::melt(target2_0.2sd_phan2_lower_all, id.vars = "group")
names(target2_0.2sd_phan2_lower_all)[2:3] <- c("age", "target2")

target2_0.2sd_phan2_lower_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))

# 0sd (once)

nsig <- 0

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_phan2$age, SGP2_phan2$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0sd_phan2 <- msir::loess.sd(SGP2_phan2$age, SGP2_phan2$target2, degree=0, nsigma = nsig, span = span)
target2_0sd_phan2_upper <- as.data.frame(cbind(target2_0sd_phan2$x, target2_0sd_phan2$upper))
names(target2_0sd_phan2_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0sd_phan2_upper$age, target2_0sd_phan2_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0sd_phan2_upper, indices) {
  d <- target2_0sd_phan2_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(0,540, 1)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0sd_phan2_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0sd_phan2_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0sd_phan2_upper_50$age <- seq(0,540, 1)
names(target2_0sd_phan2_upper_50)[1] <- "target2"

target2_0sd_phan2_upper_all <- as.data.frame(loess_boot$t)
names(target2_0sd_phan2_upper_all)[1:ncol(target2_0sd_phan2_upper_all)] <- seq(0,540, 1)
target2_0sd_phan2_upper_all$group <- seq(1,(nrow(target2_0sd_phan2_upper_all)), 1)

target2_0sd_phan2_upper_all <- reshape2::melt(target2_0sd_phan2_upper_all, id.vars = "group")
names(target2_0sd_phan2_upper_all)[2:3] <- c("age", "target2")

target2_0sd_phan2_upper_all$age2 <- rep(seq(0,540, 1),each = round(boot.R))


# precambrian

# 2sd

nsig <- 2

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_2sd_pre_upper <- as.data.frame(cbind(target2_2sd_pre$x, target2_2sd_pre$upper))
names(target2_2sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_2sd_pre_upper$age, target2_2sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_2sd_pre_upper, indices) {
  d <- target2_2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_2sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_2sd_pre_upper_50)[1] <- "target2"

target2_2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_2sd_pre_upper_all)[1:ncol(target2_2sd_pre_upper_all)] <- seq(540,4000, 10)
target2_2sd_pre_upper_all$group <- seq(1,(nrow(target2_2sd_pre_upper_all)), 1)

target2_2sd_pre_upper_all <- reshape2::melt(target2_2sd_pre_upper_all, id.vars = "group")
names(target2_2sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_2sd_pre_lower <- as.data.frame(cbind(target2_2sd_pre$x, target2_2sd_pre$lower))
names(target2_2sd_pre_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_2sd_pre_lower$age, target2_2sd_pre_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_2sd_pre_lower, indices) {
  d <- target2_2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_2sd_pre_lower_50$age <- seq(540,4000, 10)
names(target2_2sd_pre_lower_50)[1] <- "target2"

target2_2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target2_2sd_pre_lower_all)[1:ncol(target2_2sd_pre_lower_all)] <- seq(540,4000, 10)
target2_2sd_pre_lower_all$group <- seq(1,(nrow(target2_2sd_pre_lower_all)), 1)

target2_2sd_pre_lower_all <- reshape2::melt(target2_2sd_pre_lower_all, id.vars = "group")
names(target2_2sd_pre_lower_all)[2:3] <- c("age", "target2")

target2_2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.8sd

nsig <- 1.8

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.8sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_1.8sd_pre_upper <- as.data.frame(cbind(target2_1.8sd_pre$x, target2_1.8sd_pre$upper))
names(target2_1.8sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.8sd_pre_upper$age, target2_1.8sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.8sd_pre_upper, indices) {
  d <- target2_1.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_1.8sd_pre_upper_50)[1] <- "target2"

target2_1.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.8sd_pre_upper_all)[1:ncol(target2_1.8sd_pre_upper_all)] <- seq(540,4000, 10)
target2_1.8sd_pre_upper_all$group <- seq(1,(nrow(target2_1.8sd_pre_upper_all)), 1)

target2_1.8sd_pre_upper_all <- reshape2::melt(target2_1.8sd_pre_upper_all, id.vars = "group")
names(target2_1.8sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_1.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.8sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_1.8sd_pre_lower <- as.data.frame(cbind(target2_1.8sd_pre$x, target2_1.8sd_pre$lower))
names(target2_1.8sd_pre_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.8sd_pre_lower$age, target2_1.8sd_pre_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.8sd_pre_lower, indices) {
  d <- target2_1.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(target2_1.8sd_pre_lower_50)[1] <- "target2"

target2_1.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.8sd_pre_lower_all)[1:ncol(target2_1.8sd_pre_lower_all)] <- seq(540,4000, 10)
target2_1.8sd_pre_lower_all$group <- seq(1,(nrow(target2_1.8sd_pre_lower_all)), 1)

target2_1.8sd_pre_lower_all <- reshape2::melt(target2_1.8sd_pre_lower_all, id.vars = "group")
names(target2_1.8sd_pre_lower_all)[2:3] <- c("age", "target2")

target2_1.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.6sd

nsig <- 1.6

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.6sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_1.6sd_pre_upper <- as.data.frame(cbind(target2_1.6sd_pre$x, target2_1.6sd_pre$upper))
names(target2_1.6sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.6sd_pre_upper$age, target2_1.6sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.6sd_pre_upper, indices) {
  d <- target2_1.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_1.6sd_pre_upper_50)[1] <- "target2"

target2_1.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.6sd_pre_upper_all)[1:ncol(target2_1.6sd_pre_upper_all)] <- seq(540,4000, 10)
target2_1.6sd_pre_upper_all$group <- seq(1,(nrow(target2_1.6sd_pre_upper_all)), 1)

target2_1.6sd_pre_upper_all <- reshape2::melt(target2_1.6sd_pre_upper_all, id.vars = "group")
names(target2_1.6sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_1.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.6sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_1.6sd_pre_lower <- as.data.frame(cbind(target2_1.6sd_pre$x, target2_1.6sd_pre$lower))
names(target2_1.6sd_pre_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.6sd_pre_lower$age, target2_1.6sd_pre_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.6sd_pre_lower, indices) {
  d <- target2_1.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(target2_1.6sd_pre_lower_50)[1] <- "target2"

target2_1.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.6sd_pre_lower_all)[1:ncol(target2_1.6sd_pre_lower_all)] <- seq(540,4000, 10)
target2_1.6sd_pre_lower_all$group <- seq(1,(nrow(target2_1.6sd_pre_lower_all)), 1)

target2_1.6sd_pre_lower_all <- reshape2::melt(target2_1.6sd_pre_lower_all, id.vars = "group")
names(target2_1.6sd_pre_lower_all)[2:3] <- c("age", "target2")

target2_1.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.4sd

nsig <- 1.4

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.4sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_1.4sd_pre_upper <- as.data.frame(cbind(target2_1.4sd_pre$x, target2_1.4sd_pre$upper))
names(target2_1.4sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.4sd_pre_upper$age, target2_1.4sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.4sd_pre_upper, indices) {
  d <- target2_1.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_1.4sd_pre_upper_50)[1] <- "target2"

target2_1.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.4sd_pre_upper_all)[1:ncol(target2_1.4sd_pre_upper_all)] <- seq(540,4000, 10)
target2_1.4sd_pre_upper_all$group <- seq(1,(nrow(target2_1.4sd_pre_upper_all)), 1)

target2_1.4sd_pre_upper_all <- reshape2::melt(target2_1.4sd_pre_upper_all, id.vars = "group")
names(target2_1.4sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_1.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.4sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_1.4sd_pre_lower <- as.data.frame(cbind(target2_1.4sd_pre$x, target2_1.4sd_pre$lower))
names(target2_1.4sd_pre_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.4sd_pre_lower$age, target2_1.4sd_pre_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.4sd_pre_lower, indices) {
  d <- target2_1.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(target2_1.4sd_pre_lower_50)[1] <- "target2"

target2_1.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.4sd_pre_lower_all)[1:ncol(target2_1.4sd_pre_lower_all)] <- seq(540,4000, 10)
target2_1.4sd_pre_lower_all$group <- seq(1,(nrow(target2_1.4sd_pre_lower_all)), 1)

target2_1.4sd_pre_lower_all <- reshape2::melt(target2_1.4sd_pre_lower_all, id.vars = "group")
names(target2_1.4sd_pre_lower_all)[2:3] <- c("age", "target2")

target2_1.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.2sd

nsig <- 1.2

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_1.2sd_pre_upper <- as.data.frame(cbind(target2_1.2sd_pre$x, target2_1.2sd_pre$upper))
names(target2_1.2sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.2sd_pre_upper$age, target2_1.2sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.2sd_pre_upper, indices) {
  d <- target2_1.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_1.2sd_pre_upper_50)[1] <- "target2"

target2_1.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.2sd_pre_upper_all)[1:ncol(target2_1.2sd_pre_upper_all)] <- seq(540,4000, 10)
target2_1.2sd_pre_upper_all$group <- seq(1,(nrow(target2_1.2sd_pre_upper_all)), 1)

target2_1.2sd_pre_upper_all <- reshape2::melt(target2_1.2sd_pre_upper_all, id.vars = "group")
names(target2_1.2sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_1.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_1.2sd_pre_lower <- as.data.frame(cbind(target2_1.2sd_pre$x, target2_1.2sd_pre$lower))
names(target2_1.2sd_pre_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.2sd_pre_lower$age, target2_1.2sd_pre_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.2sd_pre_lower, indices) {
  d <- target2_1.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(target2_1.2sd_pre_lower_50)[1] <- "target2"

target2_1.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.2sd_pre_lower_all)[1:ncol(target2_1.2sd_pre_lower_all)] <- seq(540,4000, 10)
target2_1.2sd_pre_lower_all$group <- seq(1,(nrow(target2_1.2sd_pre_lower_all)), 1)

target2_1.2sd_pre_lower_all <- reshape2::melt(target2_1.2sd_pre_lower_all, id.vars = "group")
names(target2_1.2sd_pre_lower_all)[2:3] <- c("age", "target2")

target2_1.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 1.0sd

nsig <- 1.0

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.0sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_1.0sd_pre_upper <- as.data.frame(cbind(target2_1.0sd_pre$x, target2_1.0sd_pre$upper))
names(target2_1.0sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.0sd_pre_upper$age, target2_1.0sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_1.0sd_pre_upper, indices) {
  d <- target2_1.0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.0sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_1.0sd_pre_upper_50)[1] <- "target2"

target2_1.0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_1.0sd_pre_upper_all)[1:ncol(target2_1.0sd_pre_upper_all)] <- seq(540,4000, 10)
target2_1.0sd_pre_upper_all$group <- seq(1,(nrow(target2_1.0sd_pre_upper_all)), 1)

target2_1.0sd_pre_upper_all <- reshape2::melt(target2_1.0sd_pre_upper_all, id.vars = "group")
names(target2_1.0sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_1.0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_1.0sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_1.0sd_pre_lower <- as.data.frame(cbind(target2_1.0sd_pre$x, target2_1.0sd_pre$lower))
names(target2_1.0sd_pre_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_1.0sd_pre_lower$age, target2_1.0sd_pre_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_1.0sd_pre_lower, indices) {
  d <- target2_1.0sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_1.0sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_1.0sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_1.0sd_pre_lower_50$age <- seq(540,4000, 10)
names(target2_1.0sd_pre_lower_50)[1] <- "target2"

target2_1.0sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target2_1.0sd_pre_lower_all)[1:ncol(target2_1.0sd_pre_lower_all)] <- seq(540,4000, 10)
target2_1.0sd_pre_lower_all$group <- seq(1,(nrow(target2_1.0sd_pre_lower_all)), 1)

target2_1.0sd_pre_lower_all <- reshape2::melt(target2_1.0sd_pre_lower_all, id.vars = "group")
names(target2_1.0sd_pre_lower_all)[2:3] <- c("age", "target2")

target2_1.0sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.8sd

nsig <- 0.8

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.8sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_0.8sd_pre_upper <- as.data.frame(cbind(target2_0.8sd_pre$x, target2_0.8sd_pre$upper))
names(target2_0.8sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.8sd_pre_upper$age, target2_0.8sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.8sd_pre_upper, indices) {
  d <- target2_0.8sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.8sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.8sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.8sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_0.8sd_pre_upper_50)[1] <- "target2"

target2_0.8sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.8sd_pre_upper_all)[1:ncol(target2_0.8sd_pre_upper_all)] <- seq(540,4000, 10)
target2_0.8sd_pre_upper_all$group <- seq(1,(nrow(target2_0.8sd_pre_upper_all)), 1)

target2_0.8sd_pre_upper_all <- reshape2::melt(target2_0.8sd_pre_upper_all, id.vars = "group")
names(target2_0.8sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_0.8sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.8sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_0.8sd_pre_lower <- as.data.frame(cbind(target2_0.8sd_pre$x, target2_0.8sd_pre$lower))
names(target2_0.8sd_pre_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.8sd_pre_lower$age, target2_0.8sd_pre_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.8sd_pre_lower, indices) {
  d <- target2_0.8sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.8sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.8sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.8sd_pre_lower_50$age <- seq(540,4000, 10)
names(target2_0.8sd_pre_lower_50)[1] <- "target2"

target2_0.8sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.8sd_pre_lower_all)[1:ncol(target2_0.8sd_pre_lower_all)] <- seq(540,4000, 10)
target2_0.8sd_pre_lower_all$group <- seq(1,(nrow(target2_0.8sd_pre_lower_all)), 1)

target2_0.8sd_pre_lower_all <- reshape2::melt(target2_0.8sd_pre_lower_all, id.vars = "group")
names(target2_0.8sd_pre_lower_all)[2:3] <- c("age", "target2")

target2_0.8sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.6sd

nsig <- 0.6

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.6sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_0.6sd_pre_upper <- as.data.frame(cbind(target2_0.6sd_pre$x, target2_0.6sd_pre$upper))
names(target2_0.6sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.6sd_pre_upper$age, target2_0.6sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.6sd_pre_upper, indices) {
  d <- target2_0.6sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.6sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.6sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.6sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_0.6sd_pre_upper_50)[1] <- "target2"

target2_0.6sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.6sd_pre_upper_all)[1:ncol(target2_0.6sd_pre_upper_all)] <- seq(540,4000, 10)
target2_0.6sd_pre_upper_all$group <- seq(1,(nrow(target2_0.6sd_pre_upper_all)), 1)

target2_0.6sd_pre_upper_all <- reshape2::melt(target2_0.6sd_pre_upper_all, id.vars = "group")
names(target2_0.6sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_0.6sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.6sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_0.6sd_pre_lower <- as.data.frame(cbind(target2_0.6sd_pre$x, target2_0.6sd_pre$lower))
names(target2_0.6sd_pre_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.6sd_pre_lower$age, target2_0.6sd_pre_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.6sd_pre_lower, indices) {
  d <- target2_0.6sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.6sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.6sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.6sd_pre_lower_50$age <- seq(540,4000, 10)
names(target2_0.6sd_pre_lower_50)[1] <- "target2"

target2_0.6sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.6sd_pre_lower_all)[1:ncol(target2_0.6sd_pre_lower_all)] <- seq(540,4000, 10)
target2_0.6sd_pre_lower_all$group <- seq(1,(nrow(target2_0.6sd_pre_lower_all)), 1)

target2_0.6sd_pre_lower_all <- reshape2::melt(target2_0.6sd_pre_lower_all, id.vars = "group")
names(target2_0.6sd_pre_lower_all)[2:3] <- c("age", "target2")

target2_0.6sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.4sd

nsig <- 0.4

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.4sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_0.4sd_pre_upper <- as.data.frame(cbind(target2_0.4sd_pre$x, target2_0.4sd_pre$upper))
names(target2_0.4sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.4sd_pre_upper$age, target2_0.4sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.4sd_pre_upper, indices) {
  d <- target2_0.4sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.4sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.4sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.4sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_0.4sd_pre_upper_50)[1] <- "target2"

target2_0.4sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.4sd_pre_upper_all)[1:ncol(target2_0.4sd_pre_upper_all)] <- seq(540,4000, 10)
target2_0.4sd_pre_upper_all$group <- seq(1,(nrow(target2_0.4sd_pre_upper_all)), 1)

target2_0.4sd_pre_upper_all <- reshape2::melt(target2_0.4sd_pre_upper_all, id.vars = "group")
names(target2_0.4sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_0.4sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.4sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_0.4sd_pre_lower <- as.data.frame(cbind(target2_0.4sd_pre$x, target2_0.4sd_pre$lower))
names(target2_0.4sd_pre_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.4sd_pre_lower$age, target2_0.4sd_pre_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.4sd_pre_lower, indices) {
  d <- target2_0.4sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.4sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.4sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.4sd_pre_lower_50$age <- seq(540,4000, 10)
names(target2_0.4sd_pre_lower_50)[1] <- "target2"

target2_0.4sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.4sd_pre_lower_all)[1:ncol(target2_0.4sd_pre_lower_all)] <- seq(540,4000, 10)
target2_0.4sd_pre_lower_all$group <- seq(1,(nrow(target2_0.4sd_pre_lower_all)), 1)

target2_0.4sd_pre_lower_all <- reshape2::melt(target2_0.4sd_pre_lower_all, id.vars = "group")
names(target2_0.4sd_pre_lower_all)[2:3] <- c("age", "target2")

target2_0.4sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0.2sd

nsig <- 0.2

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_0.2sd_pre_upper <- as.data.frame(cbind(target2_0.2sd_pre$x, target2_0.2sd_pre$upper))
names(target2_0.2sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.2sd_pre_upper$age, target2_0.2sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0.2sd_pre_upper, indices) {
  d <- target2_0.2sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.2sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.2sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.2sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_0.2sd_pre_upper_50)[1] <- "target2"

target2_0.2sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_0.2sd_pre_upper_all)[1:ncol(target2_0.2sd_pre_upper_all)] <- seq(540,4000, 10)
target2_0.2sd_pre_upper_all$group <- seq(1,(nrow(target2_0.2sd_pre_upper_all)), 1)

target2_0.2sd_pre_upper_all <- reshape2::melt(target2_0.2sd_pre_upper_all, id.vars = "group")
names(target2_0.2sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_0.2sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# target2 lower

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0.2sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_0.2sd_pre_lower <- as.data.frame(cbind(target2_0.2sd_pre$x, target2_0.2sd_pre$lower))
names(target2_0.2sd_pre_lower)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0.2sd_pre_lower$age, target2_0.2sd_pre_lower$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_lower <- span

# boot lower
boot_fn <- function(target2_0.2sd_pre_lower, indices) {
  d <- target2_0.2sd_pre_lower[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_lower, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0.2sd_pre_lower, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0.2sd_pre_lower_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0.2sd_pre_lower_50$age <- seq(540,4000, 10)
names(target2_0.2sd_pre_lower_50)[1] <- "target2"

target2_0.2sd_pre_lower_all <- as.data.frame(loess_boot$t)
names(target2_0.2sd_pre_lower_all)[1:ncol(target2_0.2sd_pre_lower_all)] <- seq(540,4000, 10)
target2_0.2sd_pre_lower_all$group <- seq(1,(nrow(target2_0.2sd_pre_lower_all)), 1)

target2_0.2sd_pre_lower_all <- reshape2::melt(target2_0.2sd_pre_lower_all, id.vars = "group")
names(target2_0.2sd_pre_lower_all)[2:3] <- c("age", "target2")

target2_0.2sd_pre_lower_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

# 0sd (once)

nsig <- 0

# target2 upper

# CV span
loess.predict <- loess.as(SGP2_pre$age, SGP2_pre$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span <- loess.predict[["pars"]][["span"]]

target2_0sd_pre <- msir::loess.sd(SGP2_pre$age, SGP2_pre$target2, degree=0, nsigma = nsig, span = span)
target2_0sd_pre_upper <- as.data.frame(cbind(target2_0sd_pre$x, target2_0sd_pre$upper))
names(target2_0sd_pre_upper)[1:2] <- c("age","target2")

# CV spans
loess.predict <- loess.as(target2_0sd_pre_upper$age, target2_0sd_pre_upper$target2, degree=0, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, family = "symmetric") 
span_upper <- span

# boot upper
boot_fn <- function(target2_0sd_pre_upper, indices) {
  d <- target2_0sd_pre_upper[indices, ]
  d <- d[order(d$age),]
  loess_fit <- loess(target2 ~ age, d, span = span_upper, degree=0)
  predict(loess_fit, data.frame(age = seq(540,4000, 10)), se = T)$fit
}

boot.R <- (out$count[out$sd == nsig])/speed

loess_boot <- boot(target2_0sd_pre_upper, R = round(boot.R), statistic = boot_fn)  

# 95% confidence intervals and median smooth
target2_0sd_pre_upper_50 <- as.data.frame(apply(loess_boot$t, 2, function(x) quantile(x, .50, na.rm = TRUE)))
target2_0sd_pre_upper_50$age <- seq(540,4000, 10)
names(target2_0sd_pre_upper_50)[1] <- "target2"

target2_0sd_pre_upper_all <- as.data.frame(loess_boot$t)
names(target2_0sd_pre_upper_all)[1:ncol(target2_0sd_pre_upper_all)] <- seq(540,4000, 10)
target2_0sd_pre_upper_all$group <- seq(1,(nrow(target2_0sd_pre_upper_all)), 1)

target2_0sd_pre_upper_all <- reshape2::melt(target2_0sd_pre_upper_all, id.vars = "group")
names(target2_0sd_pre_upper_all)[2:3] <- c("age", "target2")

target2_0sd_pre_upper_all$age2 <- rep(seq(540,4000, 10),each = round(boot.R))

x <- ggplot() + 
  geom_line(data = target2_2sd_pre_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_2sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.8sd_pre_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.8sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.6sd_pre_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.6sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.4sd_pre_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.4sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.2sd_pre_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.2sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.0sd_pre_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.0sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.8sd_pre_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.8sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.6sd_pre_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.6sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.4sd_pre_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.4sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.2sd_pre_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.2sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0sd_pre_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_point(data = all, aes(age, target2, shape = source), fill = "white") + 
  scale_shape_manual(values = c(22,16)) +
  scale_x_reverse(limits = c(1800,540)) + theme_bw()


y <- ggplot() + 
  geom_line(data = target2_2sd_phan1_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_2sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.8sd_phan1_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.8sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.6sd_phan1_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.6sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.4sd_phan1_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.4sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.2sd_phan1_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.2sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.0sd_phan1_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.0sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.8sd_phan1_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.8sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.6sd_phan1_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.6sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.4sd_phan1_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.4sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.2sd_phan1_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.2sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0sd_phan1_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_2sd_phan2_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_2sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.8sd_phan2_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.8sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.6sd_phan2_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.6sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.4sd_phan2_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.4sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.2sd_phan2_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.2sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.0sd_phan2_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_1.0sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.8sd_phan2_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.8sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.6sd_phan2_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.6sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.4sd_phan2_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.4sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.2sd_phan2_lower_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0.2sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_line(data = target2_0sd_phan2_upper_all, aes(age2,target2, group = group), colour = 'blue') +
  geom_point(data = all, aes(age, target2, shape = source), fill = "white") + 
  scale_shape_manual(values = c(22,16)) +
  scale_x_reverse(limits = c(540,0)) + theme_bw()

grid.arrange(x,y,ncol = 2)

##### ribbons ####

##### Median ribbons and raw plots ####

phan1_2sd <- cbind(target2_2sd_phan1_upper_50,target2_2sd_phan1_lower_50$target2)
names(phan1_2sd)[1:3] <- c("upper","age","lower")
phan2_2sd <- cbind(target2_2sd_phan2_upper_50,target2_2sd_phan2_lower_50$target2)
names(phan2_2sd)[1:3] <- c("upper","age","lower")
phan_2sd <- rbind(phan1_2sd,phan2_2sd)
phan_2sd <- aggregate(cbind(lower,upper) ~ age, phan_2sd, mean )

phan1_1.8sd <- cbind(target2_1.8sd_phan1_upper_50,target2_1.8sd_phan1_lower_50$target2)
names(phan1_1.8sd)[1:3] <- c("upper","age","lower")
phan2_1.8sd <- cbind(target2_1.8sd_phan2_upper_50,target2_1.8sd_phan2_lower_50$target2)
names(phan2_1.8sd)[1:3] <- c("upper","age","lower")
phan_1.8sd <- rbind(phan1_1.8sd,phan2_1.8sd)
phan_1.8sd <- aggregate(cbind(lower,upper) ~ age, phan_1.8sd, mean )

phan1_1.6sd <- cbind(target2_1.6sd_phan1_upper_50,target2_1.6sd_phan1_lower_50$target2)
names(phan1_1.6sd)[1:3] <- c("upper","age","lower")
phan2_1.6sd <- cbind(target2_1.6sd_phan2_upper_50,target2_1.6sd_phan2_lower_50$target2)
names(phan2_1.6sd)[1:3] <- c("upper","age","lower")
phan_1.6sd <- rbind(phan1_1.6sd,phan2_1.6sd)
phan_1.6sd <- aggregate(cbind(lower,upper) ~ age, phan_1.6sd, mean )

phan1_1.4sd <- cbind(target2_1.4sd_phan1_upper_50,target2_1.4sd_phan1_lower_50$target2)
names(phan1_1.4sd)[1:3] <- c("upper","age","lower")
phan2_1.4sd <- cbind(target2_1.4sd_phan2_upper_50,target2_1.4sd_phan2_lower_50$target2)
names(phan2_1.4sd)[1:3] <- c("upper","age","lower")
phan_1.4sd <- rbind(phan1_1.4sd,phan2_1.4sd)
phan_1.4sd <- aggregate(cbind(lower,upper) ~ age, phan_1.4sd, mean )

phan1_1.2sd <- cbind(target2_1.2sd_phan1_upper_50,target2_1.2sd_phan1_lower_50$target2)
names(phan1_1.2sd)[1:3] <- c("upper","age","lower")
phan2_1.2sd <- cbind(target2_1.2sd_phan2_upper_50,target2_1.2sd_phan2_lower_50$target2)
names(phan2_1.2sd)[1:3] <- c("upper","age","lower")
phan_1.2sd <- rbind(phan1_1.2sd,phan2_1.2sd)
phan_1.2sd <- aggregate(cbind(lower,upper) ~ age, phan_1.2sd, mean )

phan1_1.0sd <- cbind(target2_1.0sd_phan1_upper_50,target2_1.0sd_phan1_lower_50$target2)
names(phan1_1.0sd)[1:3] <- c("upper","age","lower")
phan2_1.0sd <- cbind(target2_1.0sd_phan2_upper_50,target2_1.0sd_phan2_lower_50$target2)
names(phan2_1.0sd)[1:3] <- c("upper","age","lower")
phan_1.0sd <- rbind(phan1_1.0sd,phan2_1.0sd)
phan_1.0sd <- aggregate(cbind(lower,upper) ~ age, phan_1.0sd, mean )

phan1_0.8sd <- cbind(target2_0.8sd_phan1_upper_50,target2_0.8sd_phan1_lower_50$target2)
names(phan1_0.8sd)[1:3] <- c("upper","age","lower")
phan2_0.8sd <- cbind(target2_0.8sd_phan2_upper_50,target2_0.8sd_phan2_lower_50$target2)
names(phan2_0.8sd)[1:3] <- c("upper","age","lower")
phan_0.8sd <- rbind(phan1_0.8sd,phan2_0.8sd)
phan_0.8sd <- aggregate(cbind(lower,upper) ~ age, phan_0.8sd, mean )

phan1_0.6sd <- cbind(target2_0.6sd_phan1_upper_50,target2_0.6sd_phan1_lower_50$target2)
names(phan1_0.6sd)[1:3] <- c("upper","age","lower")
phan2_0.6sd <- cbind(target2_0.6sd_phan2_upper_50,target2_0.6sd_phan2_lower_50$target2)
names(phan2_0.6sd)[1:3] <- c("upper","age","lower")
phan_0.6sd <- rbind(phan1_0.6sd,phan2_0.6sd)
phan_0.6sd <- aggregate(cbind(lower,upper) ~ age, phan_0.6sd, mean )

phan1_0.4sd <- cbind(target2_0.4sd_phan1_upper_50,target2_0.4sd_phan1_lower_50$target2)
names(phan1_0.4sd)[1:3] <- c("upper","age","lower")
phan2_0.4sd <- cbind(target2_0.4sd_phan2_upper_50,target2_0.4sd_phan2_lower_50$target2)
names(phan2_0.4sd)[1:3] <- c("upper","age","lower")
phan_0.4sd <- rbind(phan1_0.4sd,phan2_0.4sd)
phan_0.4sd <- aggregate(cbind(lower,upper) ~ age, phan_0.4sd, mean )

phan1_0.2sd <- cbind(target2_0.2sd_phan1_upper_50,target2_0.2sd_phan1_lower_50$target2)
names(phan1_0.2sd)[1:3] <- c("upper","age","lower")
phan2_0.2sd <- cbind(target2_0.2sd_phan2_upper_50,target2_0.2sd_phan2_lower_50$target2)
names(phan2_0.2sd)[1:3] <- c("upper","age","lower")
phan_0.2sd <- rbind(phan1_0.2sd,phan2_0.2sd)
phan_0.2sd <- aggregate(cbind(lower,upper) ~ age, phan_0.2sd, mean )

phan1_0sd <- cbind(target2_0sd_phan1_upper_50)
names(phan1_0sd)[1:2] <- c("upper","age")
phan2_0sd <- cbind(target2_0sd_phan2_upper_50)
names(phan2_0sd)[1:2] <- c("upper","age")
phan_0sd <- rbind(phan1_0sd,phan2_0sd)
phan_0sd <- aggregate(cbind(upper) ~ age, phan_0sd, mean )

pre_2sd <- cbind(target2_2sd_pre_upper_50,target2_2sd_pre_lower_50$target2)
names(pre_2sd)[1:3] <- c("upper","age","lower")
pre_1.8sd <- cbind(target2_1.8sd_pre_upper_50,target2_1.8sd_pre_lower_50$target2)
names(pre_1.8sd)[1:3] <- c("upper","age","lower")
pre_1.6sd <- cbind(target2_1.6sd_pre_upper_50,target2_1.6sd_pre_lower_50$target2)
names(pre_1.6sd)[1:3] <- c("upper","age","lower")
pre_1.4sd <- cbind(target2_1.4sd_pre_upper_50,target2_1.4sd_pre_lower_50$target2)
names(pre_1.4sd)[1:3] <- c("upper","age","lower")
pre_1.2sd <- cbind(target2_1.2sd_pre_upper_50,target2_1.2sd_pre_lower_50$target2)
names(pre_1.2sd)[1:3] <- c("upper","age","lower")
pre_1.0sd <- cbind(target2_1.0sd_pre_upper_50,target2_1.0sd_pre_lower_50$target2)
names(pre_1.0sd)[1:3] <- c("upper","age","lower")
pre_0.8sd <- cbind(target2_0.8sd_pre_upper_50,target2_0.8sd_pre_lower_50$target2)
names(pre_0.8sd)[1:3] <- c("upper","age","lower")
pre_0.6sd <- cbind(target2_0.6sd_pre_upper_50,target2_0.6sd_pre_lower_50$target2)
names(pre_0.6sd)[1:3] <- c("upper","age","lower")
pre_0.4sd <- cbind(target2_0.4sd_pre_upper_50,target2_0.4sd_pre_lower_50$target2)
names(pre_0.4sd)[1:3] <- c("upper","age","lower")
pre_0.2sd <- cbind(target2_0.2sd_pre_upper_50,target2_0.2sd_pre_lower_50$target2)
names(pre_0.2sd)[1:3] <- c("upper","age","lower")

scaleFUN <- function(x) sprintf("%.1f", x)

e <- ggplot() +
  geom_ribbon(data = phan_2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_1.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_1.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_1.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_1.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = phan_1.0sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_0.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_0.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_0.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = phan_0.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_line(data = phan_0sd, aes(age,upper), colour = 'red', size = 1) +
  geom_point(data = Reinhard, aes(age, target2), fill = "white", size = 0.75, pch = 22) +
  geom_point(data = SGP2, aes(age, target2), size = 0.25) +
  scale_x_reverse(limits = c(540,0)) + theme_bw() +
  coord_cartesian(ylim = c(-2,4),expand = c(0,0)) +
  geom_hline(yintercept = c(1.8), colour = "black") +
  ylab("log10(TOC/P)") + scale_y_continuous(labels=scaleFUN) +
  theme(legend.position = "none")

f <- ggplot() +
  geom_ribbon(data = pre_2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_1.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_1.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_1.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_1.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_ribbon(data = pre_1.0sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_0.8sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_0.6sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_0.4sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  #geom_ribbon(data = pre_0.2sd, aes(age, ymin = lower, ymax = upper), alpha = 0.2, colour = "red") +
  geom_line(data = target2_0sd_pre_upper_50, aes(age,target2), colour = 'red', size = 1) +
  geom_point(data = Reinhard, aes(age, target2), fill = "white", size = 0.75, pch = 22) +
  geom_point(data = SGP2, aes(age, target2), size = 0.25) +
  scale_x_reverse(limits = c(1800,540)) + theme_bw() +
  coord_cartesian(ylim = c(-2,4),expand = c(0,0)) +
  geom_hline(yintercept = c(1.8), colour = "black") +
  ylab("log10(TOC/P)") + scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  theme(plot.tag.position = c(0.01,0.99), legend.position = "none")

grid.arrange(f,e,ncol=2)

##### Downsample, export and plot results ####

# Phanerozoic age differences
diffs <-  SGP2$age
diffs <- as.data.frame(diffs[order(diffs)])
names(diffs)[1] <- "Age"
diffs <- as.data.frame(diffs[!duplicated(diffs$Age),])
diffs_df <- round(diffs[-1,]) - round(diffs[-nrow(diffs),])
diffs_df <- as.data.frame(diffs_df)
diffs <- cbind(diffs[1:nrow(diffs)-1,],diffs_df)
names(diffs)[1] <- "Age"
names(diffs)[2] <- "Diff"
diffs <- subset(diffs, Diff > 0)
diffs$Diff2 <- diffs$Diff/2
diffs$edge <- round(diffs$Diff2+diffs$Age)
diffs$span <- c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4])
diffs$mid <- (c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4]))/2
diffs$mid <- diffs$edge - diffs$mid
diffs <- subset(diffs, span > 0)
diffs$bins <- cut(diffs$mid, c(0, unique(diffs$edge)), include.lowest = TRUE)
diffs.phan <- subset(diffs, Age < 541)

# Precambrian
diffs <-  c(SGP2$age,2000) # dummy base age added
diffs <- as.data.frame(diffs[order(diffs)])
names(diffs)[1] <- "Age"
diffs <- as.data.frame(diffs[!duplicated(diffs$Age),])
diffs <- round(diffs, digits = -1)
diffs_df <- round(diffs[-1,], digits = -1) - round(diffs[-nrow(diffs),], digits = -1)
diffs_df <- as.data.frame(diffs_df)
diffs <- cbind(diffs[1:nrow(diffs)-1,],diffs_df)
names(diffs)[1] <- "Age"
names(diffs)[2] <- "Diff"
diffs <- subset(diffs, Diff > 0)
diffs$Diff2 <- diffs$Diff/2
diffs$edge <- round(diffs$Diff2+diffs$Age)
diffs$span <- c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4])
diffs$mid <- (c(diffs[1,4], diffs[-1,4] - diffs[-nrow(diffs),4]))/2
diffs$mid <- diffs$edge - diffs$mid
diffs <- subset(diffs, span > 0)
diffs$bins <- cut(diffs$mid, c(540, unique(diffs$edge)), include.lowest = TRUE)
diffs.pre <- subset(diffs, Age >= 540)

# Phanerozoic
target2_phan_comp <- rbind(target2_2sd_phan1_upper_all,target2_1.8sd_phan1_upper_all,
                           target2_1.6sd_phan1_upper_all,target2_1.4sd_phan1_upper_all,
                           target2_1.2sd_phan1_upper_all,target2_1.0sd_phan1_upper_all,
                           target2_0.8sd_phan1_upper_all,target2_0.6sd_phan1_upper_all,
                           target2_0.4sd_phan1_upper_all,target2_0.2sd_phan1_upper_all,
                           target2_2sd_phan1_lower_all,target2_1.8sd_phan1_lower_all,
                           target2_1.6sd_phan1_lower_all,target2_1.4sd_phan1_lower_all,
                           target2_1.2sd_phan1_lower_all,target2_1.0sd_phan1_lower_all,
                           target2_0.8sd_phan1_lower_all,target2_0.6sd_phan1_lower_all,
                           target2_0.4sd_phan1_lower_all,target2_0.2sd_phan1_lower_all,
                           target2_0sd_phan1_upper_all,
                           target2_2sd_phan2_upper_all,target2_1.8sd_phan2_upper_all,
                           target2_1.6sd_phan2_upper_all,target2_1.4sd_phan2_upper_all,
                           target2_1.2sd_phan2_upper_all,target2_1.0sd_phan2_upper_all,
                           target2_0.8sd_phan2_upper_all,target2_0.6sd_phan2_upper_all,
                           target2_0.4sd_phan2_upper_all,target2_0.2sd_phan2_upper_all,
                           target2_2sd_phan2_lower_all,target2_1.8sd_phan2_lower_all,
                           target2_1.6sd_phan2_lower_all,target2_1.4sd_phan2_lower_all,
                           target2_1.2sd_phan2_lower_all,target2_1.0sd_phan2_lower_all,
                           target2_0.8sd_phan2_lower_all,target2_0.6sd_phan2_lower_all,
                           target2_0.4sd_phan2_lower_all,target2_0.2sd_phan2_lower_all,
                           target2_0sd_phan2_upper_all)

target2_phan_comp <- target2_phan_comp[complete.cases(target2_phan_comp[,]),] 

target2_phan_comp$redox <- cut(target2_phan_comp$target2, c(-4,1.8,8), include.lowest = TRUE,
                               labels = c("ferruginous","euxinic"))

target2_phan_comp.bins <- target2_phan_comp
target2_phan_comp.bins$bins <- cut(target2_phan_comp.bins$age2, c(0, unique(diffs.phan$edge)), include.lowest = TRUE)

target2_phan_comp.bins <- target2_phan_comp.bins %>% 
  group_by(bins, redox) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

target2_phan_comp.bins <- target2_phan_comp.bins[,-3]
target2_phan_comp.bins <- reshape2::dcast(target2_phan_comp.bins, bins~redox)
target2_phan_comp.bins[is.na(target2_phan_comp.bins)] <- 0

phan_results.target2.bins <- join(diffs.phan, target2_phan_comp.bins, by = "bins")
phan_results.target2.bins <- phan_results.target2.bins[,c(4:6,8:ncol(phan_results.target2.bins))]
phan_results.target2.bins <- rbind(phan_results.target2.bins[1,],phan_results.target2.bins)
phan_results.target2.bins[1,3] <- NA
phan_results.target2.bins[1,1] <- 0
phan_results.target2.melt.bins <- reshape2::melt(phan_results.target2.bins, id.vars = c("edge","span", "mid"))

#write.csv(phan_results.target2.bins, "SGP_TOCP_redox_fractions_downsampled_phanerozoic.csv") # optional export

# Precambrian
target2_pre_comp <- rbind(target2_2sd_pre_upper_all,target2_1.8sd_pre_upper_all,
                          target2_1.6sd_pre_upper_all,target2_1.4sd_pre_upper_all,
                          target2_1.2sd_pre_upper_all,target2_1.0sd_pre_upper_all,
                          target2_0.8sd_pre_upper_all,target2_0.6sd_pre_upper_all,
                          target2_0.4sd_pre_upper_all,target2_0.2sd_pre_upper_all,
                          target2_2sd_pre_lower_all,target2_1.8sd_pre_lower_all,
                          target2_1.6sd_pre_lower_all,target2_1.4sd_pre_lower_all,
                          target2_1.2sd_pre_lower_all,target2_1.0sd_pre_lower_all,
                          target2_0.8sd_pre_lower_all,target2_0.6sd_pre_lower_all,
                          target2_0.4sd_pre_lower_all,target2_0.2sd_pre_lower_all,
                          target2_0sd_pre_upper_all)

target2_pre_comp <- target2_pre_comp[complete.cases(target2_pre_comp[,]),] 

target2_pre_comp$redox <- cut(target2_pre_comp$target2, c(-2,1.8,8), include.lowest = TRUE,
                              labels = c("ferruginous","euxinic"))

target2_pre_comp.bins <- target2_pre_comp
target2_pre_comp.bins$bins <- cut(target2_pre_comp.bins$age2, c(0, unique(diffs.pre$edge)), include.lowest = TRUE)

target2_pre_comp.bins <- target2_pre_comp.bins %>% 
  group_by(bins, redox) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

target2_pre_comp.bins <- target2_pre_comp.bins[,-3]
target2_pre_comp.bins <- reshape2::dcast(target2_pre_comp.bins, bins~redox)
target2_pre_comp.bins[is.na(target2_pre_comp.bins)] <- 0

pre_results.target2.bins <- join(diffs.pre, target2_pre_comp.bins, by = "bins")
pre_results.target2.bins <- pre_results.target2.bins[,c(4:6,8:ncol(pre_results.target2.bins))]
pre_results.target2.bins <- rbind(pre_results.target2.bins[1,],pre_results.target2.bins)
pre_results.target2.bins[1,3] <- NA
pre_results.target2.bins[2,4:5] <- pre_results.target2.bins[3,4:5]
pre_results.target2.bins <- pre_results.target2.bins[-1,]
pre_results.target2.melt.bins <- reshape2::melt(pre_results.target2.bins, id.vars = c("edge","span", "mid"))

#write.csv(pre_results.target2.bins, "SGP_TOCP_redox_fractions_downsampled_precambrian.csv") # optional export

my_colours <- RColorBrewer::brewer.pal(11, "RdGy")[-c(1,3,4,5,6,7,8,10,11)]

g <- ggplot() +
  geom_bar(data = phan_results.target2.melt.bins, aes(y = value, x = mid, fill = variable, width = span),
           stat="identity") +
  theme_bw() + scale_x_reverse(limits = c(540,0), expand = c(0,0)) +
  geom_step(data = phan_results.target2.melt.bins, aes(y = value, x = edge, group = variable), 
            position = "stack", direction = "hv") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1), expand = c(0,0)) +
  ylab("Fraction euxinic") +
  scale_fill_manual(values = my_colours)

h <- ggplot() +
  geom_bar(data = pre_results.target2.melt.bins, aes(y = value, x = mid, fill = variable, width = span),
           stat="identity") +
  theme_bw() + scale_x_reverse(limits = c(1800,540), expand = c(0,0)) +
  geom_step(data = pre_results.target2.melt.bins, aes(y = value, x = edge, group = variable), 
            position = "stack", direction = "hv") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1), expand = c(0,0)) +
  ylab("Fraction euxinic")+
  scale_fill_manual(values = my_colours)

grid.arrange(b,a,d,c,f,e,h,g,ncol=2) 

# import pyrites (assuming prior export from above)

phan_results.pyrite.bins <- read.csv("phanerozoic_pyrite_types_downsampled.csv")
phan_results.pyrite.bins <- phan_results.pyrite.bins[,-1]
names(phan_results.pyrite.bins)[5:9] <- c("Type.5", "Type.1","Type.2", "Type.3", "Type.4")
phan_results.pyrite.melt.bins <- reshape2::melt(phan_results.pyrite.bins, id.vars = c("edge","span", "mid", "variable"))
names(phan_results.pyrite.melt.bins)[5] <- "type"

pre_results.pyrite.bins <- read.csv("precambrian_pyrite_types_downsampled.csv")
pre_results.pyrite.bins <- pre_results.pyrite.bins[,-1]
names(pre_results.pyrite.bins)[5:9] <- c("Type.5", "Type.1","Type.2", "Type.3", "Type.4")
pre_results.pyrite.melt.bins <- reshape2::melt(pre_results.pyrite.bins, id.vars = c("edge","span", "mid", "variable"))
names(pre_results.pyrite.melt.bins)[5] <- "type"

HEATT <- c(56, 66, 93, 116, 183, 200, 251, 359, 372, 383, 444, 514, 542) 

xdd <- read.csv("xdd_binned_results.csv") # import text mining outputs

# print

phan_results.pyrite.melt.bins$type <- factor(phan_results.pyrite.melt.bins$type, levels = c("Type.5", "Type.2", "Type.1", "Type.3", "Type.4"))
pre_results.pyrite.melt.bins$type <- factor(pre_results.pyrite.melt.bins$type, levels = c("Type.5", "Type.2", "Type.1", "Type.3", "Type.4"))

phan_results.melt.bins.mean <- subset(phan_results.pyrite.melt.bins, variable == "mean")
pre_results.melt.bins.mean <- subset(pre_results.pyrite.melt.bins, variable == "mean")

# binned ratios with uncertainty estimates

# Type.5 error
phan_results.pyrite.melt.bins.error.POM <- subset(phan_results.pyrite.melt.bins, type == "Type.5")
phan_results.pyrite.melt.bins.error.POM <- reshape2::dcast(phan_results.pyrite.melt.bins.error.POM,
                                                    edge+span+mid+type~variable)
phan_results.pyrite.melt.bins.error.POM$error <- (phan_results.pyrite.melt.bins.error.POM$max-phan_results.pyrite.melt.bins.error.POM$mean)
pre_results.pyrite.melt.bins.error.POM <- subset(pre_results.pyrite.melt.bins, type == "Type.5")
pre_results.pyrite.melt.bins.error.POM <- reshape2::dcast(pre_results.pyrite.melt.bins.error.POM,
                                                   edge+span+mid+type~variable)
pre_results.pyrite.melt.bins.error.POM$error <- (pre_results.pyrite.melt.bins.error.POM$max-pre_results.pyrite.melt.bins.error.POM$mean)

# Type.2 error
phan_results.pyrite.melt.bins.error.DOM.OX <- subset(phan_results.pyrite.melt.bins, type == "Type.2")
phan_results.pyrite.melt.bins.error.DOM.OX <- reshape2::dcast(phan_results.pyrite.melt.bins.error.DOM.OX,
                                                       edge+span+mid+type~variable)
phan_results.pyrite.melt.bins.error.DOM.OX$error <- (phan_results.pyrite.melt.bins.error.DOM.OX$max-phan_results.pyrite.melt.bins.error.DOM.OX$mean)
pre_results.pyrite.melt.bins.error.DOM.OX <- subset(pre_results.pyrite.melt.bins, type == "Type.2")
pre_results.pyrite.melt.bins.error.DOM.OX <- reshape2::dcast(pre_results.pyrite.melt.bins.error.DOM.OX,
                                                      edge+span+mid+type~variable)
pre_results.pyrite.melt.bins.error.DOM.OX$error <- (pre_results.pyrite.melt.bins.error.DOM.OX$max-pre_results.pyrite.melt.bins.error.DOM.OX$mean)

# Type.1 error
phan_results.pyrite.melt.bins.error.DOM.DFe <- subset(phan_results.pyrite.melt.bins, type == "Type.1")
phan_results.pyrite.melt.bins.error.DOM.DFe <- reshape2::dcast(phan_results.pyrite.melt.bins.error.DOM.DFe,
                                                        edge+span+mid+type~variable)
phan_results.pyrite.melt.bins.error.DOM.DFe$error <- (phan_results.pyrite.melt.bins.error.DOM.DFe$max-phan_results.pyrite.melt.bins.error.DOM.DFe$mean)
pre_results.pyrite.melt.bins.error.DOM.DFe <- subset(pre_results.pyrite.melt.bins, type == "Type.1")
pre_results.pyrite.melt.bins.error.DOM.DFe <- reshape2::dcast(pre_results.pyrite.melt.bins.error.DOM.DFe,
                                                       edge+span+mid+type~variable)
pre_results.pyrite.melt.bins.error.DOM.DFe$error <- (pre_results.pyrite.melt.bins.error.DOM.DFe$max-pre_results.pyrite.melt.bins.error.DOM.DFe$mean)

# Type.3 error
phan_results.pyrite.melt.bins.error.OM.EUX <- subset(phan_results.pyrite.melt.bins, type == "Type.3")
phan_results.pyrite.melt.bins.error.OM.EUX <- reshape2::dcast(phan_results.pyrite.melt.bins.error.OM.EUX,
                                                       edge+span+mid+type~variable)
phan_results.pyrite.melt.bins.error.OM.EUX$error <- (phan_results.pyrite.melt.bins.error.OM.EUX$max-phan_results.pyrite.melt.bins.error.OM.EUX$mean)
pre_results.pyrite.melt.bins.error.OM.EUX <- subset(pre_results.pyrite.melt.bins, type == "Type.3")
pre_results.pyrite.melt.bins.error.OM.EUX <- reshape2::dcast(pre_results.pyrite.melt.bins.error.OM.EUX,
                                                      edge+span+mid+type~variable)
pre_results.pyrite.melt.bins.error.OM.EUX$error <- (pre_results.pyrite.melt.bins.error.OM.EUX$max-pre_results.pyrite.melt.bins.error.OM.EUX$mean)

# Type.4 error
phan_results.pyrite.melt.bins.error.OX.EUX <- subset(phan_results.pyrite.melt.bins, type == "Type.4")
phan_results.pyrite.melt.bins.error.OX.EUX <- reshape2::dcast(phan_results.pyrite.melt.bins.error.OX.EUX,
                                                       edge+span+mid+type~variable)
phan_results.pyrite.melt.bins.error.OX.EUX$error <- (phan_results.pyrite.melt.bins.error.OX.EUX$max-phan_results.pyrite.melt.bins.error.OX.EUX$mean)
pre_results.pyrite.melt.bins.error.OX.EUX <- subset(pre_results.pyrite.melt.bins, type == "Type.4")
pre_results.pyrite.melt.bins.error.OX.EUX <- reshape2::dcast(pre_results.pyrite.melt.bins.error.OX.EUX,
                                                      edge+span+mid+type~variable)
pre_results.pyrite.melt.bins.error.OX.EUX$error <- (pre_results.pyrite.melt.bins.error.OX.EUX$max-pre_results.pyrite.melt.bins.error.OX.EUX$mean)

# Microenvironment vs. ambient - Fig. 4D

# error in quadrature
phan.micro.amb <- phan_results.pyrite.melt.bins.error.POM[,1:3]
phan.micro.amb$micro <- phan_results.pyrite.melt.bins.error.POM$mean+phan_results.pyrite.melt.bins.error.DOM.OX$mean+phan_results.pyrite.melt.bins.error.DOM.DFe$mean
phan.micro.amb$amb <- phan_results.pyrite.melt.bins.error.OX.EUX$mean+phan_results.pyrite.melt.bins.error.OM.EUX$mean
phan.micro.amb$ratio <- phan.micro.amb$micro/(phan.micro.amb$micro+phan.micro.amb$amb)
phan.micro.amb$den <- ((phan_results.pyrite.melt.bins.error.OX.EUX$error)^2+(phan_results.pyrite.melt.bins.error.OM.EUX$error)^2+(phan_results.pyrite.melt.bins.error.POM$error)^2+(phan_results.pyrite.melt.bins.error.DOM.OX$error)^2+(phan_results.pyrite.melt.bins.error.DOM.DFe$error)^2)^0.5
phan.micro.amb$num <- ((phan_results.pyrite.melt.bins.error.POM$error)^2+(phan_results.pyrite.melt.bins.error.DOM.OX$error)^2+(phan_results.pyrite.melt.bins.error.DOM.DFe$error)^2)^0.5
phan.micro.amb$quad <- ((phan.micro.amb$num/phan.micro.amb$micro)^2+(phan.micro.amb$den/(phan.micro.amb$amb+phan.micro.amb$micro))^2)^0.5
phan.micro.amb$quad <- (phan.micro.amb$ratio)*(phan.micro.amb$quad)
phan.micro.amb$upper <- phan.micro.amb$ratio+phan.micro.amb$quad
phan.micro.amb$lower <- phan.micro.amb$ratio-phan.micro.amb$quad
phan.micro.amb[,10:11][phan.micro.amb[,10:11] < 0] <- 0
phan.micro.amb[,10:11][phan.micro.amb[,10:11] > 1] <- 1

#error in quad
pre.micro.amb <- pre_results.pyrite.melt.bins.error.POM[,1:3]
pre.micro.amb$micro <- pre_results.pyrite.melt.bins.error.POM$mean+pre_results.pyrite.melt.bins.error.DOM.OX$mean+pre_results.pyrite.melt.bins.error.DOM.DFe$mean
pre.micro.amb$amb <- pre_results.pyrite.melt.bins.error.OX.EUX$mean+pre_results.pyrite.melt.bins.error.OM.EUX$mean
pre.micro.amb$ratio <- pre.micro.amb$micro/(pre.micro.amb$micro+pre.micro.amb$amb)
pre.micro.amb$den <- ((pre_results.pyrite.melt.bins.error.OX.EUX$error)^2+(pre_results.pyrite.melt.bins.error.OM.EUX$error)^2+(pre_results.pyrite.melt.bins.error.POM$error)^2+(pre_results.pyrite.melt.bins.error.DOM.OX$error)^2+(pre_results.pyrite.melt.bins.error.DOM.DFe$error)^2)^0.5
pre.micro.amb$num <- ((pre_results.pyrite.melt.bins.error.POM$error)^2+(pre_results.pyrite.melt.bins.error.DOM.OX$error)^2+(pre_results.pyrite.melt.bins.error.DOM.DFe$error)^2)^0.5
pre.micro.amb$quad <- ((pre.micro.amb$num/pre.micro.amb$micro)^2+(pre.micro.amb$den/(pre.micro.amb$amb+pre.micro.amb$micro))^2)^0.5
pre.micro.amb$quad <- (pre.micro.amb$ratio)*(pre.micro.amb$quad)
pre.micro.amb$upper <- pre.micro.amb$ratio+pre.micro.amb$quad
pre.micro.amb$lower <- pre.micro.amb$ratio-pre.micro.amb$quad
pre.micro.amb[,10:11][pre.micro.amb[,10:11] < 0] <- 0
pre.micro.amb[,10:11][pre.micro.amb[,10:11] > 1] <- 1

# optional export

#write.csv(phan.micro.amb, "pyrite_redox_ratio_downsampled_phanerozoic.csv")
#write.csv(pre.micro.amb, "pyrite_redox_ratio_downsampled_precambrian.csv")

A <- ggplot() +
  geom_bar(data = phan_results.melt.bins.mean, aes(y = value, x = mid, fill = type, width = span),
           stat="identity") +
  theme_bw() + scale_x_reverse(limits = c(550,0)) +
  geom_step(data = phan_results.melt.bins.mean, aes(y = value, x = edge, group = type), 
            position = "stack", direction = "hv") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0,1), expand = c(0,0))


B <- ggplot() +
  geom_bar(data = pre_results.melt.bins.mean, aes(y = value, x = mid, fill = type, width = span),
           stat="identity") +
  theme_bw() + scale_x_reverse(limits = c(1800,540)) +
  geom_step(data = pre_results.melt.bins.mean, aes(y = value, x = edge, group = type), 
            position = "stack", direction = "hv") +
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0,1), expand = c(0,0)) +
  labs(tag = "X") +
  theme(plot.tag.position = c(0.01,0.99))


C <- ggplot() +
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
  ylab("Fraction ferruginous") +
  coord_cartesian(ylim = c(0,1), clip = "off") +
  coord_geo(dat = list("periods", "eras"), xlim = c(540, 0), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr"))
  
D <- ggplot() +
  #geom_stepribbon(data = pre.micro.amb,
                  #aes(x = edge, ymin = lower, ymax = upper), fill = "grey50") +
  geom_step(data = pre.micro.amb, aes(y =ratio, x = edge),
            direction = "hv", colour = "black", size = 0.75) +
  theme_bw() + scale_x_reverse(limits = c(1800,540), expand = c(0,0)) +
  geom_step(data = subset(pre_results.target2.melt.bins, variable == "ferruginous"), 
            aes(y = value, x = edge, group = variable), 
            position = "stack", direction = "hv", size = 0.75, colour = "red") +
  geom_step(data = subset(pre_results.target1.melt.bins, variable == "ferruginous"), 
            aes(y = value, x = edge, group = variable), 
            position = "stack", direction = "hv", size = 0.75, colour = "blue") +
  ylab("Fraction ferruginous") +
  coord_cartesian(ylim = c(0,1), clip = "off") +
  coord_geo(dat = list("periods", "eras"), xlim = c(1800, 540), size = 3, 
            height = list(unit(0.75, "lines"), unit(0.75, "lines")),
            pos = list("b", "b"), abbrv = list(TRUE, TRUE),
            skip = c("Msr")) +
  labs(tag = "C") +
  theme(plot.tag.position = c(0.01,0.99))


grid.arrange(b,a,f,e,D,C,ncol=2)

# export SGP ferruginous % for palaeo-mapping #

body = '{"type":"samples","filters":{"lithology_class":["sedimentary","metamorphic"]},"show":["fe","mn","fe_carb","fe_ox","fe_mag","fe_py","fe_py_fe_hr","fe_hr_fe_t","toc","mo","u","p","ni","v","interpreted_age","fe_hr","analysis_ref_long", "coord_lat","coord_long", "strat_name"]}'

r <- POST("https://sgp-search.io/api/v1/post", body=body, 
          httr::add_headers(`accept` = 'application/json'), 
          httr::content_type('application/json'))

out <- content(r, "text")

SGP <- jsonlite::fromJSON(out)

# write SGP refs
SGP_refs <- subset(SGP, select = `analysis reference (long)`)
SGP_refs <- as.data.frame(SGP_refs[!duplicated(SGP_refs),])
SGP_refs <- as.data.frame(SGP_refs[complete.cases(SGP_refs),])
names(SGP_refs)[1] <- "SGP reference"
SGP_refs <- as.data.frame(SGP_refs[order(SGP_refs$`SGP reference`),])
#write.csv(SGP_refs, "SGP_refs.csv")

SGP[,c(3:5,7:21)] <- sapply(SGP[,c(3:5,7:21)], as.numeric)

ppm <- dplyr::select(SGP,contains(" (ppm)")) 
ppb <- dplyr::select(SGP,contains(" (ppb)"))
ppb <- ppb/1000
wt <- dplyr::select(SGP,contains(" (wt%)"))
wt <- wt*10000

SGP2 <- cbind(SGP[,c(1:6,22)],ppm,ppb,wt)

names(SGP2) <- gsub(x = names(SGP2), pattern = " (ppm)", replacement = "", fixed = TRUE)
names(SGP2) <- gsub(x = names(SGP2), pattern = " (ppb)", replacement = "", fixed = TRUE)
names(SGP2) <- gsub(x = names(SGP2), pattern = " (wt%)", replacement = "", fixed = TRUE)
names(SGP2) <- gsub(x = names(SGP2), pattern = "-", replacement = "", fixed = TRUE)
names(SGP2)[3:6] <- c("lat","long","age","strat")

SGP2$lat <- round(SGP2$lat, digits = 1)
SGP2$long <- round(SGP2$long, digits = 1)
SGP2$age <- round(SGP2$age, digits = 0)

# define intervals based on pyrite type distributions #

zones <- c(3400,3100,3000,2800,2100,1700,1500,1400,900,800,740,690,655,610,532,514,480,456,419,383,372,367,359,330,299,280,250,201,54,0)

SGP2$zones <- cut(SGP2$age, zones)

get_midpoint <- function(cut_label) {
  mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ","))))
}

SGP2$zones2 <- sapply(SGP2$zones, get_midpoint)

###### SGP - Fepy/FeHR ####

# anoxic subset

SGP2 <- subset(SGP2, (Fecarb+Feox+Femag+Fepy)/Fe > 0.38 & Fe > (0.5*10000))

SGP3 <- subset(SGP2, (Fepy)/(Fecarb+Feox+Femag+Fepy) >= 0.7)
SGP3$Fe_redox <- "Euxinic"

SGP4 <- subset(SGP2, (Fepy)/(Fecarb+Feox+Femag+Fepy) < 0.7)
SGP4$Fe_redox <- "Ferruginous"

SGP_Fe <- rbind(SGP3, SGP4)

SGP_Fe <- SGP_Fe %>% 
  group_by(lat, long, age, Fe_redox) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

SGP_Fe1 <- subset(SGP_Fe, Fe_redox != "Euxinic")
SGP_Fe2 <- subset(SGP_Fe, Fe_redox == "Euxinic" & perc == 1)
SGP_Fe2$Fe_redox <- "Ferruginous"
SGP_Fe2$perc <- 0

SGP_Fe <- rbind(SGP_Fe1, SGP_Fe2)

SGP_Fe <- SGP_Fe[complete.cases(SGP_Fe),]

write.csv(SGP_Fe, "SGP_Fe_ferruginous_fractions.csv") # export - excluding strat names

# optional with strat names

SGP_Fe <- rbind(SGP3, SGP4)

SGP_Fe <- SGP_Fe %>% 
  group_by(lat, long, strat, age, Fe_redox) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

SGP_Fe1 <- subset(SGP_Fe, Fe_redox != "Euxinic")
SGP_Fe2 <- subset(SGP_Fe, Fe_redox == "Euxinic" & perc == 1)
SGP_Fe2$Fe_redox <- "Ferruginous"
SGP_Fe2$perc <- 0

SGP_Fe <- rbind(SGP_Fe1, SGP_Fe2)

SGP_Fe <- SGP_Fe[complete.cases(SGP_Fe),]

#write.csv(SGP_Fe, "SGP_Fe_ferruginous_fractions_strat.csv") # export

######### END ########



