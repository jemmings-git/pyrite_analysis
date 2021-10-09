# Sedimentary pyrite mega-analysis
### We analysed the xDD library for mentions of pyrite types through geological time (Figure 1). Then we conducted statistical analysis, machine learning and interpolation of pyrite and bulk sediment geochemistry in order to map modes of anoxia through geological time and space (Figures 2-6).

### This code accompanies a manuscript submitted to *Science Advances* in May 2021
#### Title: Pyrite mega-analysis reveals modes of anoxia through geological time  
#### Authors: Joseph F. Emmings<sup>1,2</sup>*, Simon W. Poulton<sup>3</sup>, Joanna Walsh<sup>4,5</sup>, Kathryn A. Leeming<sup>1</sup>, Ian Ross<sup>6</sup>, Shanan Peters<sup>7</sup>
#### Affiliations 
<sup>1</sup>British Geological Survey, Keyworth, Nottingham, NG12 5GG, UK.  
<sup>2</sup>School of Geography, Geology and the Environment, University of Leicester, Leicester, LE1 7RH, UK.  
<sup>3</sup>School of Earth and Environment, University of Leeds, Leeds, LS2 9JT, UK.  
<sup>4</sup>Lyell Centre, British Geological Survey, Riccarton, Edinburgh EH14 4AS, UK.  
<sup>5</sup>Ordnance Survey, Explorer House, Adanac Drive, Southampton, SO16 0AS, UK.  
<sup>6</sup>Department of Computer Sciences, University of Wisconsin–Madison, Madison, Wisconsin 53706, USA.   
<sup>7</sup>Department of Geoscience, University of Wisconsin–Madison, Madison, Wisconsin 53706, USA.   
\* Present address: CGG, Tyn Y Coed, LLanrhos, Llandudno, LL30 1SA.

#### We recommend implementation in RStudio.

## Steps following deployment of the [pyrite machine reading application](https://github.com/jemmings-git/pyrite_app):

#### R scripts:
 * `install.R` - install the project's dependencies
 * `macrostrat_data.R` - download Macrostrat definitions
 * `preparation.R` - [import xDD pyrite results](https://geodeepdive.org/app_output/jemmings_with_pyrite_24Oct2019.zip), collect and link related data from Macrostrat, do binning on the dataframe, save the results
 * `wilkin_snippets.R` - process Wilkin-pyrite mentions (Fig 1C) via the [xDD snippets API](https://xdd.wisc.edu/api/snippets?term=Wilkin,framboid&full_results=true&inclusive=true&clean&known_terms=stratigraphic_names)
 * `analysis.R` - produce a series of analytical plots illustrating the observable xDD pyrite distributions as presented on Figure 1. This code is split into two parts, dealing with Figs 1A-C and Supplementary figures 
 * `stats.R`- conduct data manipulation, statistical analysis and machine learning as presented on Figures 2-5 and supplementary materials. Each code block (for each figure) should be run consecutively. Pyrite trace element and SGP data are downloaded and analysed. This script reproduces Figures 2-5.
 * `gplates.R` - reconstruct palaeogeographic coordinates and interpolate results (in time and space), via the [GPlates Web Service](https://gws.gplates.org/) interfaced with the [*chronosphere* R package](https://github.com/chronosphere-portal/r_package/). This script reproduces Figures 6-7.

#### Datasets & results:

 * `xdd_binned_results.csv` - the results of text mining as presented on Figures 1A-1B  
 * `xdd_binned_results_strat_only.csv` - the results of text mining as presented on supplementary xDD figures  
 * `pyrite_refs.txt` - list of refences in the xDD library containing at least one mention of framboidal, nodular, concretionary or undifferentiated pyrite matched to a (meta)sedimentary rock in the Macrostrat database. This is an expanded version of the reference list contained in the Supplementary Materials.
 * `SGP_refs.csv` - list of references underlying the Fe-speciation and TOC/P data stored in the [Sedimentary Geochemistry and Palaeoenvironments Project database, SGP](http://sgp-search.io/)  
 * `precambrian_pyrite_types_downsampled.csv` - Results of statistical analysis, machine learning, downsampling and aggregation of [Precambrian pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332) (as presented on Figure 4)
 * `phanerozoic_pyrite_types_downsampled.csv` - Results of statistical analysis, machine learning, downsampling and aggregation of [Phanerozoic pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332) (as presented on Figure 4)
 * `precambrian_pyrite_types_downsampled_singleCV.csv` - Results of statistical analysis, machine learning, downsampling and aggregation of [Precambrian pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332), using a single coherent cross-validated loess span (as presented in supplementary materials)  
 * `phanerozoic_pyrite_types_downsampled_singleCV.csv` - Results of statistical analysis, machine learning, downsampling and aggregation of [Phanerozoic pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332), using a single coherent cross-validated loess span (as presented in supplementary materials)  
 * `SEDEX.csv` - an open access USGS dataset that is not stored in this repository. [Singer *et al*., 2009](https://pubs.usgs.gov/of/2009/1252/). [Download](https://pubs.usgs.gov/of/2009/1252/SedZn-PbEX2009.xls). BIFs and glaciations are hard-coded.
 * 
 * `pyrite_stats_with_long_lats.csv` - aggregated pyrite type fractions derived from [pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332) by age and present day geographic coordinates (added manually using [Google Earth](https://earth.google.com/web/)). Used in figures 5-6.

# Acknowledgements

### Figure 1 is derived from the [GeoDeepDive](https://geodeepdive.org) (now known as xDD) library and machine reading system
As the basis of our pyrite application we [forked](https://github.com/jemmings-git/pyrite_app) the [original text mining application](https://github.com/UW-Macrostrat/stromatolites_demo) used in [Peters, Husson and Wilcots 2017, *Geology*](http://doi.org/10.1130/G38931.1). 

### Figure 2 utilises Phase 1 data from the [Sedimentary Geochemistry and Palaeoenvironments Project (SGP)](http://sgp-search.io/)
The citable references for the Phase 1 data product are:

Mehra et al. 2021, GSA Today, 31, https://doi.org/10.1130/GSATG484A.1.

Farrell et al., Geobiology (in review) (2021).

### Figures 3-5 utilise [previously published pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332) : https://doi.org/10.1130/GEOL.S.12456332
The data derive from:

Large et al., Earth and Planetary Science Letters 389, 209-220 (2014).

Large et al., Gondwana Research 28, 1282-1293 (2015).

Large et al., Earth and Planetary Science Letters 428, 139-150 (2015).

Large et al., Mineralium Deposita 54, 485-506 (2019).

Mukherjee and Large, Geology 48, 1018-1022 (2020).

