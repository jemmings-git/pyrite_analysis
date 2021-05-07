# Analyse the xDD library for pyrite types through geological time (Figure 1), and conduct statistical analysis of pyrite trace element geochemistry (Figures 2-5)

## Steps

At present we recommend running this in RStudio

 * `install.R` - install the project's dependencies
 * `preparation.R` - collect and link related data from Macrostrat, do binning on the dataframe, save the results
 * `analysis.R` - produce a series of analytical plots illustrating the observable pyrite distributions as presented on Figure 1. This code is split into two parts, dealing with Figs 1A-B and 1C-D.
 * `stats.R`- conduct statistics and data manipulation as presented on Figures 2-6 and supplementary materials. Each code chunk (for each figure) should be run consecutively.

`xdd_binned_results.csv` - the results of text mining as presented on Figures 1A-1B

`xdd_binned_results_strat_only.csv` - the results of text mining as presented on Figures 1C-1D

`redox_zones.txt` - anoxic temporal stages A-E derived from interpretation of the text mining outputs

`pyrite_refs.txt` - list of refences in the xDD library containing at least one mention of framboidal, nodular, concretionary or undifferentiated pyrite matched to a (meta)sedimentary rock in the Macrostrat database. This is an expanded version of the reference list contained in the Supplementary Materials.

`SGP_refs.txt` - list of references used to construct Figure 2
