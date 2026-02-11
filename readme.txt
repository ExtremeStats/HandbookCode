NAME OF DATASET 

waveheight.RData

(The data are also available from the R package DATAstudio)

DATA DESCRIPTION 

The data considered are significant wave heights $H_s$ extracted from the opensource RESOURCECODE wave database developed at Ifremer. It consist in high-resolution regional hindcasts, including  the evolution in time (hourly) and space of the directional spectrum and several integrated wave parameters. The modelled area extends from 12°W to 13.5°E longitude, and from 36°N to 63°N latitude and covers the period from 1994 to 2020. The data have been downloaded using the R package resourcecode. The following pre-processing steps have been performed:
1. Select a spatial aera extending from 5.25°W to 4.25°E longitude, and from 46.6°N to 47.5°N latitude
2. Sample 100 locations within this area (using the seed 111 in R)
3. At this 100 locations, extract only daily (i.e. observation at noon) $H_s$ from January to March (included) between 1995 to 2015. This is achieved using the function get_parameters from the package resourcecode.

The final .RData file contains the following objects:
- hs_dat: a dataframe containing  observations at each location (n=1895 observations, d=100 locations)
- Locations: a dataframe containing the coordinates of each location (Longitude, Latitude and Node of the RESOURCECODE grid)


CONTACT INFORMATION

Data and toolbox maintainer: Nicolas raillard,  nicolas.raillard@ifremer.fr
Website: https://resourcecode.ifremer.fr/

ACKNOWLEDGMENTS

The ResourceCODE project has received support under the framework of the OCEANERA-NET COFUND project, with funding provided by national/regional sources and co-funding by the European Union's Horizon 2020 research and innovation programme. 

REFERENCES

Accensi Mickael, Alday Gonzalez Matias Felipe, Maisondieu Christophe (2022). RESOURCECODE - Resource Characterization to Reduce the Cost of Energy through Coordinated Data Enterprise. Database user manual. WP3| Database development. https://doi.org/10.13155/86306

Raillard Nicolas (2024). resourcecode: Access to the RESOURCECODE hindcast database. R package version 0.0.1, https://nraillard.github.io/resourcecode/
