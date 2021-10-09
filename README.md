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
 * `wilkin_snippets.json` - the aggregated Wilkin-pyrite xDD snippets results (reproduced in this repository for ease). Used in **Figure 1C & supplementary xDD figures**
 *  `SEDEX.csv` - an open access USGS dataset that is not stored in this repository. [Singer *et al*., 2009](https://pubs.usgs.gov/of/2009/1252/). [Download](https://pubs.usgs.gov/of/2009/1252/SedZn-PbEX2009.xls). BIFs and glaciations are hard-coded. These data are plotted on **Figures 1 & 4**
 * `xdd_binned_results.csv` - the results of text mining as presented on **Figures 1A-1B**  
 * `xdd_binned_results_strat_only.csv` - the results of text mining as presented on **supplementary xDD figures**  
 * `precambrian_pyrite_types_downsampled.csv` - Results of statistical analysis, machine learning, downsampling and aggregation of [Precambrian pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332) (as presented on **Figure 4**)
 * `phanerozoic_pyrite_types_downsampled.csv` - Results of statistical analysis, machine learning, downsampling and aggregation of [Phanerozoic pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332) (as presented on **Figure 4**)
 * `precambrian_pyrite_types_downsampled_singleCV.csv` - Results of statistical analysis, machine learning, downsampling and aggregation of [Precambrian pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332), using a single coherent cross-validated loess span (as presented in the **supplementary materials**)  
 * `phanerozoic_pyrite_types_downsampled_singleCV.csv` - Results of statistical analysis, machine learning, downsampling and aggregation of [Phanerozoic pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332), using a single coherent cross-validated loess span (as presented in the **supplementary materials**)  
 * `pyrite_redox_ratio_downsampled_phanerozoic.csv` - Phanerozoic age-interpolated ratio of microenvironment vs ambient euxinic conditions derived from pyrite trace element geochemistry. As presented on **Figure 4D, 5C, 6A**
 * `pyrite_redox_ratio_downsampled_precambrian.csv`- Precambrian age-interpolated ratio of microenvironment vs ambient euxinic conditions derived from pyrite trace element geochemistry. As presented on **Figure 4D, 5C, 6A**
 * `SGP_Fe_redox_fractions_downsampled_phanerozoic.csv` - Phanerozoic age-interpolated fraction of Phanerozoic ferruginous conditions based on the key Fe speciation thresholds (see main text). A derivative of the SGP dataset shown in **Figure 5A**. As presented on **Figure 5C, 6A**
 * `SGP_Fe_redox_fractions_downsampled_precambrian.csv` - Precambrian age-interpolated fraction of Precambrian ferruginous conditions based on the key Fe speciation thresholds (see main text). A derivative of the SGP dataset shown in **Figure 5A**. As presented on **Figure 5C, 6A**
 * `SGP_TOCP_redox_fractions_downsampled_phanerozoic.csv` - Phanerozoic age-interpolated fraction of Phanerozoic ferruginous conditions based on a key TOC/P threshold (see main text). A derivative of the SGP dataset shown in **Figure 5B**. As presented on **Figure 5C, 6A**
 * `SGP_TOCP_redox_fractions_downsampled_precambrian.csv` - Precambrian age-interpolated fraction of Precambrian ferruginous conditions based on a key TOC/P threshold (see main text). A derivative of the SGP dataset shown in **Figure 5B**. As presented on **Figure 5C, 6A**
 * `pyrite_combined_manual.csv` - present day geographic coordinates (clat, clng) for the xDD pyrite framboid and/or nodule-bearing rock record. Regarding *Phanerozoic* xDD pyrite-bearing rocks outside the Macrostrat focal area, present day geographic coordinates are estimates based on the type sections reported in the [British Geological Survey](https://www.bgs.ac.uk/technologies/the-bgs-lexicon-of-named-rock-units/) and [Geoscience Australia](https://asud.ga.gov.au/) stratigraphic lexicons). These xDD data were converted to a quasi-ferruginous fraction, and subsequently interpolated using reconstructed palaeogeographic coordinates. Presented on **Figure 6**
 * `pyrite_stats_with_long_lats.csv` - the fraction of ferruginous vs. euxinic conditions (microenvironment vs. ambient) on an age (1 Ma bin) and unit-by-unit basis for [pyrite trace element analyses](https://doi.org/10.1130/GEOL.S.12456332) reported by [Mukherjee and Large, 2021](https://doi.org/10.1130/GEOL.S.12456332). Present day geographic coordinates were added manually by inspection of the original dataset sample descriptions, cross-referenced to coordinates extracted primarily from [Google Earth](https://earth.google.com/web/), and rarely SGP from entries. These data were interpolated using reconstructed palaeogeographic coordinates. Presented on **Figure 6**
 * `SGP_Fe_ferruginous_fractions.csv` - the fraction of ferruginous vs. euxinic conditions on an age (1 Ma bin) and long-lat basis for [SGP Fe-speciation data](http://sgp-search.io/). These data were interpolated using reconstructed palaeogeographic coordinates. Presented on **Figure 6**

#### Reference lists:
 * `pyrite_refs.txt` - list of refences in the xDD library containing at least one mention of framboidal, nodular, concretionary or undifferentiated pyrite matched to a (meta)sedimentary rock in the Macrostrat database. This is an expanded version of the reference list contained in the Supplementary Materials.
 * `SGP_refs.csv` - list of references underlying the Fe-speciation and TOC/P data stored in the [Sedimentary Geochemistry and Palaeoenvironments Project database, SGP](http://sgp-search.io/)  

# Acknowledgements

### Figure 1 is derived from the [GeoDeepDive](https://geodeepdive.org) (now known as xDD) library and machine reading system
As the basis of our pyrite application we [forked](https://github.com/jemmings-git/pyrite_app) the [original text mining application](https://github.com/UW-Macrostrat/stromatolites_demo) used in [Peters, Husson and Wilcots 2017, *Geology*](http://doi.org/10.1130/G38931.1). 

### Figures 2-7 utilise [previously published pyrite trace element data](https://doi.org/10.1130/GEOL.S.12456332):  

[Large et al., 2014. *Earth and Planetary Science Letters*, 389, 209-220.](http://dx.doi.org/10.1016/j.epsl.2013.12.020)  
[Large et al., 2015. *Gondwana Research*, 28, 1282-1293.](http://dx.doi.org/10.1016/j.gr.2015.06.004)  
[Large et al., 2015. *Earth and Planetary Science Letters*, 428, 139-150.](http://dx.doi.org/10.1016/j.epsl.2015.07.026)  
[Large et al., 2019. *Mineralium Deposita*, 54, 485-506.](https://doi.org/10.1007/s00126-019-00873-9)  
[Mukherjee and Large, 2020. *Geology*, 48, 1018-1022.](https://doi.org/10.1130/G47890.1)  

### Figures 5-7 utilise Fe-speciation and TOC/P data from the [Sedimentary Geochemistry and Palaeoenvironments Project (SGP)](http://sgp-search.io/)
The citable references for the Phase 1 data product are:

[Mehra et al., 2021. *GSA Today*, 31](https://doi.org/10.1130/GSATG484A.1.)  
[Farrell et al., 2021. *Geobiology*, 00, 1– 12.](https://doi.org/10.1111/gbi.12462)  

### Reconstructions ###

Palaeogeographic reconstructions derive from the [GPlates Web Service](https://gws.gplates.org/) interfaced with the [*chronosphere* R package](https://github.com/chronosphere-portal/r_package/). The reconstructions derive from the model of [Scotese, C. R. 2016. Tutorial: PALEOMAP PaleoAtlas for GPlates and the PaleoData Plotter Program](https://www.earthbyte.org/paleomap-paleoatlas-for-gplates/). Palaeo-digital elevation models derive from [Scotese, C. R, and Wright, N. 2018. PALEOMAP Paleodigital Elevation Models (PaleoDEMS) for the Phanerozoic](https://www.earthbyte.org/paleodem-resource-scotese-and-wright-2018/).

Chronostratigraphic timescales were plotted in R using the [*deeptime* R package](https://github.com/willgearty/deeptime).

Please refer to the **Supplementary materials** for the complete R package reference list.

