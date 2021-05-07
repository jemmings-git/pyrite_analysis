# Analyse the xDD library for pyrite types through geological time (Figure 1), and conduct statistical analysis of pyrite trace element geochemistry (Figures 2-5)

We recommend implementation in RStudio.

## Steps following deployment of the pyrite machine reading application (https://github.com/jemmings-git/pyrite_app):

 * `install.R` - install the project's dependencies
 * `preparation.R` - collect and link related data from Macrostrat, do binning on the dataframe, save the results
 * `analysis.R` - produce a series of analytical plots illustrating the observable pyrite distributions as presented on Figure 1. This code is split into two parts, dealing with Figs 1A-B and 1C-D.
 * `stats.R`- conduct statistics and data manipulation as presented on Figures 2-5 and supplementary materials. Each code block (for each figure) should be run consecutively. SGP and pyrite trace element data are downloaded, culled and processed within R. This script reproduces Figures 2-5.

`xdd_binned_results.csv` - the results of text mining as presented on Figures 1A-1B

`xdd_binned_results_strat_only.csv` - the results of text mining as presented on Figures 1C-1D

`redox_zones.txt` - anoxic temporal stages A-E derived from interpretation of the text mining outputs

`pyrite_refs.txt` - list of refences in the xDD library containing at least one mention of framboidal, nodular, concretionary or undifferentiated pyrite matched to a (meta)sedimentary rock in the Macrostrat database. This is an expanded version of the reference list contained in the Supplementary Materials.

`SGP_refs.txt` - list of references used to construct Figure 2 (Sedimentary Geochemistry and Palaeoenvironments Project, SGP: http://sgp-search.io/)

`precambrian_pyrite_types_downsampled.csv` - Results of statistical analysis, machine learning, downsampling and aggregation of Precambrian pyrite trace element data (as presented on Figure 5)

`phanerozoic_pyrite_types_downsampled.csv` - Results of statistical analysis, machine learning, downsampling and aggregation of Phanerozoic pyrite trace element data (as presented on Figure 5)

# Acknowledgements

## Figure 2 utilises Phase 1 data from the Sedimentary Geochemistry and Palaeoenvironments Project (SGP) http://sgp-search.io/
The citable references for the Phase 1 data product are:
Mehra et al. 2021, GSA Today, 31, https://doi.org/10.1130/GSATG484A.1.
Farrell et al., Geobiology (in review) (2021).

## Figures 3-5 utilise pyrite trace element data accessible here: https://doi.org/10.1130/GEOL.S.12456332
The data derive from:
Large et al., Earth and Planetary Science Letters 389, 209-220 (2014).
Large et al., Gondwana Research 28, 1282-1293 (2015).
Large et al., Earth and Planetary Science Letters 428, 139-150 (2015).
Large et al., Mineralium Deposita 54, 485-506 (2019).
Mukherjee and Large, Geology 48, 1018-1022 (2020).

