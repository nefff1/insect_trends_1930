# Codes and data to analyse trends of saproxylic beetle and butterfly species richness between 1930 and 2021 in Switzerland

[![](https://zenodo.org/badge/DOI/10.5281/zenodo.17256301.svg)](https://doi.org/10.5281/zenodo.17256301)

This repository contains codes and additional data which were used in the analyses for the following manuscript:

Neff F, Bollmann K, Chittaro Y, Gossner MM, Herzog F, Korner-Nievergelt F, Litsios G, Martínez-Núñez C, Moretti M, Rey E, Sanchez A, Knop E. **Ninety-year trends reveal sharpest insect declines mid-20th century.**

In this study, based on species records data, mean annual occupancy was estimated for 595 saproxylic beetle species and 216 buttefly species (Papilionoidea, incl. Zygaenidae moths) in six biogeographic zones in Switzerland. The data cover the years 1930–2021 (occupancy estimates were determined per two-year interval; 1930/31, 1932/33, etc.). Trends in saproxylic beetle and butterfly species richness were reconstructed for the biogeographic zones and the whole of Switzerland. Richness trends were related to a set of environmental variables (land use, climate) and to species traits (size, specialisation, temperature niche).

The following R code files are included in the folder *R_Code*:

-   **R_Main_analyses.R**: Contains all the main analyses necessary to reproduce figures and numbers in the manuscript. Based on processed data (occupancy estimates, drivers, traits).

-   **R_Drivers.R**: Contains codes used to calculate the driver variables from the raw data

-   **R_Traits.R**: Contains codes used to calculate the trait variables from the raw data

-   **R_Occupancy_detection_models_sapro.R**: Contains the codes used to fit the occupancy-detection models for the saproxylic beetle species. The models were fitted on an HPC cluster.

-   **R_Occupancy_detection_models_butter.R**: Contains the codes used to fit the occupancy-detection models for the butterfly species. The models were fitted on an HPC cluster.

-   **f_occ_det.R**: Helper function used to fit occupancy-detection models in Stan (called by the previous files).

The following Stan code file is included in the folder *Stan_Code*:

-   **Stan_occ_det_cmdstan.stan**: Stan model code for the occupancy-detection models (called by the respective R code files)

The sources of external data are indicated in the code files. Occupancy estimates are available from Zenodo (<https://doi.org/10.5281/zenodo.17255265>). The following additional data files necessary to reproduce the analyses are available from the folder *Data*:

-   **splist_sapro.txt**: Text file containing the names of the 595 saproxylic beetle species analysed

-   **splist_butter.txt**: Text file containing the names of the 216 butterfly species analysed

-   **Traits/Traits_saproxylic_beetles.csv**: Categorical trait data for all saproxylic beetle species. Each row is a species. Contains the following columns:

    -   Name_std: species name

    -   size_cat: body size category (categorical: small, medium, large)

    -   stenotopy: habitat specialisation (categorical: eurytopic, oligotopic, stenotopic)

    -   food_spec: food specialisation (categorical: polyphagous, oligophagous, monophagous)

    -   Tniche: temperature niche (categorical: cold, intermediate, warm)

-   **Traits/Traits_butterflies.csv**: Categorical trait data for all butterfly species. Each row is a species. Contains the following columns:

    -   Name_std: species name

    -   size_cat: body size category (categorical: small, medium, large)

    -   stenotopy: habitat specialisation (categorical: eurytopic, oligotopic, stenotopic)

    -   food_spec: food specialisation (categorical: polyphagous, oligophagous, monophagous)

    -   Tniche: temperature niche (categorical: cold, intermediate, warm)

    -   voltinism: number of generations per year (categorical: 0.5, 1, 2, 3)

    -   hibernation: overwintering stage (categorical: egg, egg/larva, larva, larva/pupa, pupa, adult)

-   **Drivers/Drivers.csv**: Processed data on environmental variables. Each row contains data for one biogeographic zone and interval. Contains the following columns:

    -   zone: biogeographic zone

    -   year_mean: Mean year of the interval

    -   length_interval: Length of the intervals, either 8 (eight years) or 12 (twelve years)

    -   forest_area_change: estimated change in forest area (real number)

    -   forest_area_abs: estimated absolute forest area (proportion) (real number)

    -   wood_harvest_int_change: estimated change in wood harvest intensity (real number)

    -   after_storm: factor denoting whether the interval follows a storm event (boolean)

    -   human_pop_change: estimated change in human population density (real number)

    -   mechanisation_abs: estimated absolute mechanisation degree (real number)

    -   mechanisation_period: factor denoting whether the intervals falls within the major mechanisation period (boolean)

    -   grassland_area_change: estimated change in grassland area (real number)

    -   grassland_area_abs: estimated absolute grassland area (proportion) (real number)

    -   temperature_change: estimated mean annual temperature change (real number)

    -   temperature_abs: estimate absolute mean annual temperature anomaly (real number)

-   **Drivers/Census_data_raw.csv**: Raw data from different censuses, which were used to calculate some of the environmental variables. At the level of Swiss cantons. Contains the following columns:

    -   year: year of the data

    -   canton: name of the canton

    -   value: value of the censused variable (real number)

    -   parameter: name of the censused variable. One of the following:

        -   "forest: area" (forest area in ha)

        -   "forest: wood harvest" (wood harvest in m3)

        -   "grassland area" (grassland area in ha)

        -   "human population" (human population)

        -   "tractors" (number of tractors)

    -   forest_owner: for forest data, denotes whether data is on publicly or privately owned forests

    -   source: source of the data

-   **Other/GBIF/d_taxonomy_GBIF_sapro.csv**: File used to match GBIF taxon keys to taxonomic names used in the present analyses. Used during temperature niche calculation of saproxylic beetles.
