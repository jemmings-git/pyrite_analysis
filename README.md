# Experiment to analyse the GeoDeepDive library for pyrite types through geological time, and statistical analysis of pyrite trace element geochemistry

## Steps

At present we recommend running this in RStudio but aim to supply a docker image with which to reproduce the analysis

 * `install.R` - install the project's dependencies
 * `preparation.R` - collect and link related data from Macrostrat, do binning on the dataframe, save the results
 * `analysis.R` - produce a series of analytical plots illustrating the observable ocean redox events as presented on Figure 1
 * `stats.R`- conduct statistics and data manipulation as presented on Figures 2-5 and supplementary materials
