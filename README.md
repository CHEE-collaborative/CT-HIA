## About
This repository holds the data processing, simulation, and analysis code for the research project, "Simulating desegregation through affordable housing development: an environmental health impact assessment of Connecticut zoning law" by Saira Prasanth, Nire Oloyede, Xuezhixing Zhang, Kai Chen, and Daniel Carrión, 2024. A preprint is available through medRxiv here: https://doi.org/10.1101/2024.02.13.24302645

## Installation and Data Downloads
First, install R version 4.3.2. R is available for free download from the R Project for Statistical Computing here: https://www.r-project.org/

Then, download a copy of the repository to your device and open the project, "CT HIA.Rproj", in R or R Studio. The `renv` package will be installed, if it is not currently installed on your device. To load package versions identical to those used in this analysis (as recorded in the lockfile, "renv.lock"), run `renv::restore()`. This will only change versions of packages loaded within this project, and package versions for any other projects will remain unchanged.

Below are instructions to reproduce the datasets used in this analysis.

To retrieve population data from the U.S. Census Bureau using the `tidycensus` package, you will need to create a Census API key. You can request a key by following the instructions at this link: https://api.census.gov/data/key_signup.html

To retrieve income limits from the U.S. Department of Housing and Urban Development (HUD) using the `hudr` package, you will need to create a HUD account and access token. You can create an access token by following the instructions at this link: https://www.huduser.gov/portal/dataset/fmr-api.html

The following datasets must be downloaded manually from their respective host sites. This will require you to have a free NASA Earthdata account. If you do not have an account, you can create one by following the directions outlined here: https://urs.earthdata.nasa.gov/documentation/for_users/how_to_register 

The first dataset to download manually is annual 0.0083-degree resolution nitrogen dioxide (NO2) levels, from NASA's Health and Air Quality Applied Sciences Team (HAQAST) Nitrogen Dioxide Surface-Level Annual Average Concentrations V1 product.

+ Visit https://disc.gsfc.nasa.gov/datasets/SFC_NITROGEN_DIOXIDE_CONC_1/summary
+ Click "Get Data" and use the following parameters to download one netCDF file for the main analysis:
    - Download Method: Get Original Files
    - Refine Date Range: 2019-01-01 to 2019-12-31
    - File Format: netCDF
+ Click "Get Data". 
+ Click the link with the netCDF file name. You may be prompted to log in with your NASA Earthdata account. The downloaded file should be named: "SurfaceNO2_0.0083deg_2019.nc"; place this in the "data" folder.
+ Download a second netCDF file using the previous three steps but instead setting the "Refine Date Range" parameter to the following date range: 2018-01-01 to 2018-12-31. This is for the sensitivity analysis. The downloaded file should be named: "SurfaceNO2_0.0083deg_2018.nc"; place this in the "data" folder.

The second dataset to download manually is the Gridded Population of the World (GPW) Version 4, Population Count, 2020, with 30-second resolution from the Columbia University Center for International Earth Science Information Network (CIESIN) and NASA Socioeconomic Data and Applications Center (SEDAC).

+ Visit https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-rev11/data-download
+ You may be prompted to log in with your NASA Earthdata account.
+ Set the following parameters using the drop-down menus:
    - Temporal: Single Year
    - FileFormat: GeoTiff
    - Resolution: 30 Second (approx. 1km)
+ Select the checkbox for Year 2020.
+ Click "Create Download". The downloaded .zip folder should be named: "gpw-v4-population-count-rev11_2020_30_sec_tif"; place this in the "data" folder.

The remaining datasets used in this analysis are stored in this repository in the "data" folder. These datasets are:
+ Monthly 1-km x 1-km normalized difference vegetation index (NDVI) in 2019 and 2018, retrieved from NASA's Terra Moderate Resolution Imaging Spectroradiometer (MODIS) Vegetation Indices Monthly (MOD13A3) Version 6.1, available at: https://lpdaac.usgs.gov/products/mod13a3v061/
    - There are 24 total GeoTIFF files containing NDVI data in this folder (1 for each month in 2019 and 2018).
+ Statewide annual and monthly (May-September) crude mortality rates by race and ethnicity in 2019 and 2018 in Connecticut, retrieved from the Centers for Disease Control and Prevention Wide-ranging ONline Data for Epidemiologic Research (CDC WONDER) Underlying Cause of Death by Single-Race Categories database, available at: https://wonder.cdc.gov/ucd-icd10-expanded.html
    - There are 12 total text files containing mortality data in this folder (annual and monthly rates, stratified by race, by ethnicity, and by both race and ethnicity, in 2019 and 2018)
    
In aggregate, after all datasets have been downloaded to the "data" folder, the size of this folder will be approximately 7.67 GB.

## Running the Analysis
Open the `_targets.R` and `Functions.R` scripts. Then, in the R console, run the following lines (also located at the beginning of `_targets.R`) to load the packages required to define the targets pipeline and to execute the pipeline:
```
library(targets)
library(tarchetypes)
library(future)
library(crew)
tar_make()
```
You can cancel the execution at any point by pressing the escape key.

The results shared in the corresponding publication are based on 1,000 repeats of the simulation. The number of repeats is determined by the value of `n` in line 18 of `_targets.R` and is currently set to 10 to allow users to run a smaller, quicker version of the simulation by default. You can assign any value to `n` that is greater than 1 to test these, and you can set `n <- 1000` to replicate the full simulation.

## Viewing Results

Each table and figure has a corresponding "target," which can be read using `tar_read()` and loaded into your R environment using `tar_load()`. For example, the main tables and figures can be viewed in R by running the following lines:
```
tar_read(exhibit1)
tar_read(exhibit2)
tar_read(exhibit3)
tar_read(exhibit4)
```
Supplementary tables and figures can also be read this way, by supplying the corresponding target names (e.g. "exhibitA13").

Additionally, after running the pipeline, all main or supplementary figures and tables wil be saved in the "output" folder created in your local repository.
