
# download data from URL
download_data <- function() {
  
  # set longer timeout
  options(timeout = 7200)

  # PM2.5
  if(!file.exists(here("data", "2019_pm25_daily_average.txt.gz"))){download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2019_pm25_daily_average.txt.gz", destfile = here("data", "2019_pm25_daily_average.txt.gz"))}
  if(!file.exists(here("data", "2018_pm25_daily_average.txt.gz"))){download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2018_pm25_daily_average.txt.gz", destfile = here("data", "2018_pm25_daily_average.txt.gz"))}

  # O3
  if(!file.exists(here("data", "2019_ozone_daily_8hour_maximum.txt.gz"))){download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2019_ozone_daily_8hour_maximum.txt.gz", destfile = here("data", "2019_ozone_daily_8hour_maximum.txt.gz"))}
  if(!file.exists(here("data", "2018_ozone_daily_8hour_maximum.txt.gz"))){download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2018_ozone_daily_8hour_maximum.txt.gz", destfile = here("data", "2018_ozone_daily_8hour_maximum.txt.gz"))}

  # noise
  if(!file.exists(here("data", "CONUS_road_noise_2020.zip"))){download.file("https://www.bts.gov/bts-net-storage/CONUS_road_noise_2020.zip", destfile = here("data", "CONUS_road_noise_2020.zip"))}
  if(!file.exists(here("data", "CONUS_road_noise_2018.zip"))){download.file("https://www.bts.gov/bts-net-storage/CONUS_road_noise_2018.zip", destfile = here("data", "CONUS_road_noise_2018.zip"))}

  # housing
  if(!file.exists(here("data", "Affordable_Housing_by_Town_2011-2022.csv"))){download.file("https://data.ct.gov/resource/3udy-56vi.csv", destfile = here("data", "Affordable_Housing_by_Town_2011-2022.csv"))}

  downloads <- 1
  return(downloads)
}

# load township polygons
Load_CT_townships <- function() {
  
  options(tigris_use_cache = TRUE)
  
  # get township polygons
  CT_townships <- county_subdivisions(state = "CT", year = 2019) %>% 
    filter(NAME != "County subdivisions not defined") %>% 
    dplyr::select(GEOID, NAME, geometry)
  
  return(CT_townships)
}

get_exposures_town <- function(year = "2019", CT_townships, downloads) {
  
  # unzip GPW data
  unzip(here("data", "gpw-v4-population-count-rev11_2020_30_sec_tif.zip"), files = "gpw_v4_population_count_rev11_2020_30_sec.tif", exdir = here("data"))
  
  # load GPW for population weighting
  gpw <- rast(here("data", "gpw_v4_population_count_rev11_2020_30_sec.tif")) %>% 
    # change CRS
    project(crs(CT_townships)) %>%
    # match extent
    crop(ext(CT_townships))
  
  # read in PM2.5 daily average by census tract from EPA FAQSD
  pm_tract <- read.delim(gzfile(here("data", paste0(year, "_pm25_daily_average.txt.gz"))), sep = ",") %>%
    ## include only Connecticut tracts - FIPS starting with 9
    filter(str_detect(`FIPS`, "^9")) %>% 
    ## annual average by tract
    group_by(FIPS) %>% 
    summarize(pm_mean = mean(pm25_daily_average.ug.m3.)) %>% 
    ## format FIPS to match GEOID format used by tidycensus get_acs
    mutate(GEOID = str_pad(FIPS, 11, side = "left", pad = "0"))
  
  # summarize mean tract-level ozone estimates by town
  pm_town <- tracts(state = "CT", year = year) %>%
    # join tract-level exposures (excluding area of Long Island Sound)
    right_join(pm_tract) %>% 
    # calculate areal average by town
    st_interpolate_aw(x = .["pm_mean"], to = CT_townships, extensive = FALSE, keep_NA = TRUE) %>% 
    # add town names based on row order
    mutate(NAME = CT_townships$NAME) %>% 
    st_drop_geometry()
  
  # read in ozone daily 8-hour maximum by census tract from EPA FAQSD (external file)
  o3_tract <- read.delim(gzfile(here("data", paste0(year, "_ozone_daily_8hour_maximum.txt.gz"))), sep = ",") %>%
    ## include only Connecticut tracts - FIPS starting with 9
    filter(str_detect(`FIPS`, "^9")) %>% 
    ## annual average by tract
    group_by(FIPS) %>% 
    summarize(o3_mean = mean(ozone_daily_8hour_maximum.ppb.)) %>% 
    ## format FIPS to match GEOID format used by tidycensus get_acs
    mutate(GEOID = str_pad(FIPS, 11, side = "left", pad = "0"))
  
  # summarize mean tract-level ozone estimates by town
  o3_town <- tracts(state = "CT", year = year) %>%
    # join tract-level exposures (excluding area of Long Island Sound)
    right_join(o3_tract) %>% 
    # calculate areal average by town
    st_interpolate_aw(x = .["o3_mean"], to = CT_townships, extensive = FALSE) %>% 
    # add town names
    mutate(NAME = CT_townships$NAME) %>% 
    st_drop_geometry()
  
  # read in NO2 data
  no2_r <- rast(here("data", paste0("SurfaceNO2_0.0083deg_", year, ".nc"))) %>%
    project(crs(CT_townships)) %>%
    # match extent
    crop(ext(CT_townships))
  
  # resample GPW to match grid cells
  gpw_no2 <- resample(gpw, no2_r)
  
  # extract raster values and join to town polygon
  no2_vals <- exact_extract(no2_r, CT_townships, weights = gpw_no2, fun = "weighted_mean")
  
  # write function to find saturation VP: t>0 formula
  calc_saturationVP <- function(t) exp(34.494 - (4924.99/(t + 237.1)))/(t + 105)^1.57  
  
  # assign heat season dates (5/1 - 9/30)
  heat_season <- seq.Date(as.Date(paste0(year, "/5/1")), as.Date(paste0(year, "/9/30")), "days") %>% 
    format("%m/%d")
  
  # pull temp data from Daymet
  download_daymet_ncss(location = c(42.05, -73.73, 40.98, -71.79), start = year, end = year, param = "tmax", path = here("data"))
  
  #pull vapor pressure data from Daymet
  download_daymet_ncss(location = c(42.05, -73.73, 40.98, -71.79), start = year, end = year, param = "vp", path = here("data"))
  
  # use spatraster to adjust temp data crs
  temp <- rast(here("data", paste0("tmax_daily_", year, "_ncss.nc"))) %>% 
    project(crs(CT_townships)) %>%
    subset(121:273) %>% 
    crop(ext(CT_townships))
  
  # use spatraster to adjust vp data crs (ADJUST)
  vp <- rast(here("data", paste0("vp_daily_", year, "_ncss.nc"))) %>% 
    project(crs(CT_townships)) %>% 
    subset(121:273) %>% 
    crop(ext(CT_townships))
  
  # initialize raster stacks for saturation VP and relative humidity
  satvp <- temp
  rh <- temp
  heatindex <- temp
  
  for(i in 1:length(temp[1])) {
    satvp[[i]] <- temp[[i]] %>% calc_saturationVP()
    rh[[i]] <- vp[[i]]/satvp[[i]]
    
    for(j in length(temp[[i]])) {
      heatindex[[i]][j] <- heat.index(t = temp[[i]][j], rh = rh[[i]][j], temperature.metric = "celsius", output.metric = "celsius", round = 2)
    }
  }
  
  gpw_heat <- terra::resample(gpw, heatindex)
  
  # replace NA with 0
  gpw_heat[] <- ifelse(is.na(gpw_heat[]), 0, gpw_heat[])
  
  # use exactextractr to average heat index over townships
  heatindex_town <- purrr::map(1:153, ~ exact_extract(heatindex[[.x]], CT_townships, weights = gpw_heat, "weighted_mean"))
  
  # initialize data frame for daily heat index by town
  heat_town <- CT_townships %>% st_drop_geometry()
  
  # append to CT town shapefile
  for(i in 1:length(heatindex_town)) {
    heat_town[, 2+i] <- heatindex_town[[i]]
  }
  
  heat_town <- heat_town %>% 
    rename_with(~ heat_season[1:153], .cols = 3:155)
  
  heat_town_avg <- heat_town %>% 
    # pivot to long format
    pivot_longer(cols = 3:155, names_to = "date", values_to = "heatindex") %>% 
    # average daily max heat index across 5/1-9/30 (for plotting/summarizing annual exposures)
    group_by(GEOID, NAME) %>% 
    summarize(heat_annual_mean = mean(heatindex)) %>% 
    ungroup()
  
  # list NDVI file names
  ndvi_files <- list.files(path = here("data"), pattern = paste0("MOD13A3.061__1_km_monthly_NDVI_doy", year)) %>%
    as.list()
  
  # create raster for each NDVI file
  ndvi <- purrr::map(ndvi_files, ~ crop(project(rast(here("data", .x)), crs(CT_townships)), ext(CT_townships)))
  
  # resample GPW to match NDVI grid cells
  gpw_ndvi <- resample(gpw, ndvi[[1]])
  
  # replace NA with 0
  gpw_ndvi[] <- ifelse(is.na(gpw_ndvi[]), 0, gpw_ndvi[])
  
  # extract average NDVI over townships
  ndvi_town_vals <- map(ndvi, ~ exact_extract(.x, CT_townships, weights = gpw_ndvi, fun = "weighted_mean"))
  
  # initialize data frame for average NDVI by town
  ndvi_town <- CT_townships
  
  # append to CT town shapefile
  for(i in 1:length(ndvi_town_vals)) {
    ndvi_town[, 3+i] <- ndvi_town_vals[[i]]
  }
  
  # pivot to long format
  ndvi_town <- ndvi_town %>%
    pivot_longer(cols = 4:15, names_to = "month", values_to = "ndvi") %>% 
    group_by(NAME, geometry) %>% 
    summarize(ndvi_mean = mean(ndvi)) %>% 
    ungroup() %>% 
    st_drop_geometry()
  
  # load noise data in 2018
  noise_file_names <- paste0("CT_road_noise_", 2018, ".tif")
  
  # unzip data
  unzip(here("data", "CONUS_road_noise_2018.zip"), exdir = here("data"))
  unzip(here("data", "CONUS_road_noise_2020.zip"), exdir = here("data"))
  
  noise_data <- rast(here("data", "CONUS_road_noise_2018", "State_rasters", noise_file_names)) %>% 
    project(crs(CT_townships)) %>% 
    crop(ext(CT_townships))
  
  # resample gpw to match noise grid cells
  gpw_noise <- resample(gpw, noise_data)
  # replace NA with 0
  gpw_noise[] <- ifelse(is.na(gpw_noise[]), 0, gpw_noise[])
  
  # extract noise values over townships
  noise_2018 <- exact_extract(noise_data, CT_townships, weights = gpw_noise, fun = "weighted_mean")
  
  # load noise data in 2018 and 2020
  noise_file_names <- paste0("CT_road_noise_", 2020, ".tif")
  noise_data <- rast(here("data", "CONUS_road_noise_2020", "State_rasters", noise_file_names)) %>% 
    project(crs(CT_townships)) %>% 
    crop(ext(CT_townships))
  
  # resample gpw to match noise grid cells
  gpw_noise <- resample(gpw, noise_data)
  # replace NA with 0
  gpw_noise[] <- ifelse(is.na(gpw_noise[]), 0, gpw_noise[])
  
  # extract noise values over townships
  noise_2020 <- exact_extract(noise_data, CT_townships, weights = gpw_noise, fun = "weighted_mean")
  
  # use linear imputed values for 2019
  noise_town <- data.frame(NAME = CT_townships$NAME, year = year) %>% 
    mutate(noise_mean = if_else(year == "2018", noise_2018, (noise_2018 + noise_2020)/2))
  
  # merge all exposures
  exposures_town <- CT_townships %>% 
    left_join(pm_town, by = "NAME") %>% 
    left_join(o3_town, by = "NAME") %>% 
    mutate(no2_mean = no2_vals) %>% 
    left_join(heat_town_avg, by = c("NAME", "GEOID")) %>% 
    left_join(ndvi_town, by = "NAME") %>% 
    left_join(noise_town, by = "NAME") %>% 
    left_join(heat_town, by = c("NAME", "GEOID"))
  
  return(exposures_town)
}

allocate_housing <- function(year_h = "2019", downloads) {
  
  # set 10% benchmark
  pct = 10
  
  # read housing data from CT Open Data
  housing_town <- read.csv(here("data", "Affordable_Housing_by_Town_2011-2022.csv")) %>% 
    filter(year == as.numeric(year_h)) %>% 
    mutate(at_least_10 = cut(percent_assisted, breaks = c(0,9.99999,100),
                             labels = c("No", "Yes"))) %>% 
    rename("NAME" = "town") %>% 
    # calculate new units needed to meet % affordable benchmark (0 if already met)
    mutate(new_units = case_when(percent_assisted >= pct ~ 0, percent_assisted < pct ~  (0.01*pct*X_2010_census-total_assisted)/(1-0.01*pct))) %>% 
    mutate(new_units = ceiling(new_units))
  
  return(housing_town)
  
}

get_mortality <- function(year = "2019") {
  
  # read statewide all-cause mortality rates from CDC WONDER by race/ethnicity
  mortality_black_white <- read.delim(here("data", paste0("allcause_", year, "_race_ethnicity.txt")), sep = "\t") %>% 
    filter(Single.Race.6 %in% c("White", "Black or African American"), Hispanic.Origin == "Not Hispanic or Latino")
  
  mortality_asian <- read.delim(here("data", paste0("allcause_", year, "_race.txt")), sep = "\t") %>% 
    filter(Single.Race.6 == "Asian") %>% 
    mutate(Population = as.character(Population))
  
  mortality_hispanic <- read.delim(here("data", paste0("allcause_", year, "_ethnicity.txt")), sep = "\t") %>% 
    filter(Hispanic.Origin == "Hispanic or Latino")
  
  mortality <- bind_rows(mortality_black_white, mortality_asian, mortality_hispanic) %>%
    mutate(race_ethnicity = case_match(Single.Race.6, "Black or African American" ~ "black_nh", "White" ~ "white_nh", "Asian" ~ "asian", NA ~ "hispanic")) %>% 
    dplyr::select(race_ethnicity, mr_pre = Crude.Rate, Population) %>% 
    mutate(mr_pre = as.numeric(mr_pre))
  
  return(mortality)
  
}

get_mortality_daily <- function(year = "2019", mortality) {
  
  # load monthly mortality counts to estimate daily deaths
  mortality_monthly_black_white <- read.delim(here("data", paste0("allcause_", year, "_race_ethnicity_month.txt")), sep = "\t") %>% 
    filter(Single.Race.6 %in% c("White", "Black or African American"), Hispanic.Origin == "Not Hispanic or Latino")
  
  mortality_monthly_asian <- read.delim(here("data", paste0("allcause_", year, "_race_month.txt")), sep = "\t")%>% 
    filter(Single.Race.6 == "Asian")
  
  mortality_monthly_hispanic <- read.delim(here("data", paste0("allcause_", year, "_ethnicity_month.txt")), sep = "\t") %>% 
    filter(Hispanic.Origin == "Hispanic or Latino") %>% 
    mutate(Deaths = as.numeric(Deaths))
  
  mortality_daily <- bind_rows(mortality_monthly_black_white, mortality_monthly_asian, mortality_monthly_hispanic) %>%
    mutate(month = str_sub(Month, start = 1, end = 3)) %>% 
    mutate(days = case_when(month %in% c("May", "Jul", "Aug") ~ 31, month %in% c("Jun", "Sep") ~ 30)) %>% 
    mutate(race_ethnicity = case_match(Single.Race.6, "Black or African American" ~ "black_nh", "White" ~ "white_nh", "Asian" ~ "asian", NA ~ "hispanic")) %>% 
    dplyr::select(race_ethnicity, month, Month.Code, Deaths, days) %>%
    left_join(mortality, by = "race_ethnicity") %>% 
    mutate(daily_rate_per100000 = Deaths/as.numeric(Population)*100000/days)
  
  return(mortality_daily)
}

# pull 3-person household low-income limits by town
Load_AMI_limits <- function(year, hud_key = hud_key) {

  # assign HUD access token
  hud_key <- hud_key
  
  # get FIPS code by town
  town_codes <- get_hud_fmr_listcounties(hud_key = hud_key, stateid = "CT")
  
  # retrieve low-income limits by town in blocks (cannot retrieve more than 100 rows at once)
  Sys.sleep(60)
  AMI_limits_a <- purrr::map(town_codes$fips_code[1:100], ~ get_hud_il_data(entityid = .x, yr = year, hud_key = hud_key) %>% .$low %>% .$il80_p3) %>% unlist()
  
  Sys.sleep(60)
  AMI_limits_b <- purrr::map(town_codes$fips_code[101:169], ~ get_hud_il_data(entityid = .x, yr = year, hud_key = hud_key) %>% .$low %>% .$il80_p3) %>% unlist()
  
  # combine limits
  AMI_limits <- data.frame(GEOID = town_codes$fips_code, NAME = sub(" town", "", x = town_codes$town_name), limit = c(AMI_limits_a, AMI_limits_b))
  
  return(AMI_limits)
}

# load population data and clean population variables
get_low_income_pop <- function(year, AMI_limits){
  
  # list ACS variables with income brackets by race/ethnicity
  incomes_vars <- c("B19001B", "B19001D", "B19001H", "B19001I")
  
  # list ethnoracial group labels
  race_eth <- c("black_nh", "asian", "white_nh", "hispanic")
  
  # list variable suffixes
  variable <- c("002", "003", "004", "005", "006", "007", "008", "009", "010", "011", "012", "013", "014", "015", "016", "017")
  
  # list maximum income for each bracket
  max_income <- c(9999, 14999, 19999, 24999, 29999, 34999, 39999, 44999, 49999, 59999, 74999, 99999, 124999, 149999, 199999, 200000)
  
  # join maximum incomes for each bracket with variable suffixes
  income_brackets <- data.frame(variable, max_income)
  
  # get population by race and race by ethnicity in each town
  eth_by_race <- get_decennial(geography = "county subdivision", state = "CT", variables = c("P1_004N",  "P2_006N"), year = 2020, output = "wide") %>% 
    mutate(black_nh_pct = P2_006N/P1_004N)
  
  # get average household size by town
  household_size <- get_acs(geography = "county subdivision", state = "CT", variables = "B25010_001", year = as.numeric(year), output = "tidy", survey = "acs5") %>% 
    filter(!str_detect(NAME, "County subdivisions not defined"))
  
  # get households in each income bracket by town for each race/ethnicity group (wide table format)
  incomes_race <- map(incomes_vars, ~ get_acs(geography = "county subdivision", state = "CT", table = .x, year = as.numeric(year), output = "tidy", survey = "acs5")) %>% 
    data.frame() %>%
    dplyr::select(GEOID, variable, black = estimate, asian = estimate.1, white_nh = estimate.2, hispanic = estimate.3) %>%
    # join ethnicity by race estimates for each town
    left_join(eth_by_race, by = "GEOID") %>%
    # estimate non-Hispanic Black population using race by ethnicity per town
    mutate(black_nh = round(black*black_nh_pct, 0)) %>% 
    # filter out "County subdivisions not defined"
    filter(!str_detect(NAME, "County subdivisions not defined")) %>% 
    mutate(variable = sub(".*_", "", variable)) %>% 
    dplyr::select(GEOID, variable, black_nh, white_nh, asian, hispanic)
  
  # merge township data
  low_income <- incomes_race %>% 
    left_join(AMI_limits, by = "GEOID") %>% 
    # inner join to filter out "001" variable (total)
    inner_join(income_brackets, by = "variable") %>%
    # join average household size by town
    left_join(household_size, by = "GEOID") %>% 
    # make table long
    pivot_longer(cols = c(black_nh, white_nh, asian, hispanic), names_to = "race_ethnicity") %>% 
    # clean town name
    mutate(NAME = sub(" town.*", "", NAME.x))
  
  low_income_pop <- low_income %>%
    # replace estimates above income limit with 0
    mutate(value_li = case_when(max_income <= limit ~ value, max_income > limit ~ 0)) %>% 
    group_by(NAME, race_ethnicity, max_income) %>% 
    # sum total low-income households in each race/ethnicity group and low-income bracket by town
    summarize(pop_pre = sum(value_li), avg_size = first(estimate))
  
  return(low_income_pop)
}

get_town_polygon_sf <- function(CT_townships, housing_town, low_income_pop) {
  
  town_polygon <- CT_townships %>% 
    # join housing data by town - using inner join to exclude area of Long Island Sound
    inner_join(housing_town, by = "NAME")
  
  # list ethnoracial group labels
  race_eth <- c("black_nh", "asian", "white_nh", "hispanic")
  
  # join population with low income for each race/ethnicity by town
  town_polygon_list <- map(race_eth, ~ left_join(town_polygon, filter(low_income_pop, race_ethnicity == .x), by = "NAME"))
  
  # bind shapefiles by race/ethnicity
  town_polygon_sf <- bind_rows(town_polygon_list) %>% 
    filter(is.na(race_ethnicity) == FALSE)
  
  return(town_polygon_sf)
}

get_weights <- function(CT_townships) {
  
  # unzip GPW data
  unzip(here("data", "gpw-v4-population-count-rev11_2020_30_sec_tif.zip"), files = "gpw_v4_population_count_rev11_2020_30_sec.tif", exdir = here("data"))
  
  gpw <- read_stars(here("data", "gpw_v4_population_count_rev11_2020_30_sec.tif"))
  
  CT_towns <- CT_townships %>% 
    st_transform(st_crs(gpw))
  
  gpw_CT <- st_crop(gpw, CT_towns)
  
  centroids <- st_as_sf(gpw_CT) %>% 
    st_centroid() %>%
    rename(pop = "gpw_v4_population_count_rev11_2020_30_sec.tif") %>%
    st_join(., CT_towns, join = st_intersects) %>%
    st_transform(., crs = 2234) %>%
    bind_cols(., st_coordinates(.)) %>%
    st_drop_geometry() %>%
    group_by(NAME) %>%
    summarise(pop_centroid_x = weighted.mean(x = X, w = pop),
              pop_centroid_y = weighted.mean(x = Y, w = pop)) %>%
    st_as_sf(coords = c("pop_centroid_x", "pop_centroid_y"), crs = 2234) %>% 
    st_transform(crs = crs(CT_townships))
  
  # assign weights to each pair of towns equal to inverse distance squared
  ## list every town pair
  weights <- data.frame(NAME = rep(centroids$NAME, times = nrow(centroids))) %>% 
    mutate(moved_to = rep(centroids$NAME, times = rep(nrow(centroids), nrow(centroids)))) %>%
    ## create columns for pre- and post-move town centroids
    left_join(centroids, by = "NAME") %>% 
    rename(moved_from = NAME, NAME = moved_to, centr_from = geometry) %>% 
    left_join(centroids, by = "NAME") %>% 
    rename(moved_to = NAME, centr_to = geometry) %>% 
    ## calculate distances between centroids
    mutate(distance = st_distance(centr_from, centr_to, by_element = TRUE)) %>% 
    mutate(distance = as.numeric(distance))
  
  return(weights)
  
}

get_weights_post <- function(weights, CT_townships, housing_town, low_income_pop, AMI_limits) {
  
  # initialize data frame for same-town pairs
  dist_0 <- weights %>% 
    # filter for same-town pairs (staying in town)
    filter(moved_from == moved_to) %>% 
    dplyr::select(moved_from, moved_to)
  
  # replace 0 distance with estimated town radius
  for(i in 1:nrow(dist_0)) {
    dist_0$replace[i] <- CT_townships %>%
      filter(NAME == dist_0$moved_from[i]) %>% 
      # retrieve convex hull
      st_convex_hull() %>% 
      # cast to multipoint
      st_cast("MULTIPOINT") %>%
      # cast to point
      st_cast("POINT") %>%
      # calculate distances between points
      st_distance() %>%
      #find max distance between any two points on town convex hull and divide by two
      max(na.rm = TRUE)/2
  }
  
  # select columns from low_income_pop to join households by pre-move town with weights
  families_pre_town <- low_income_pop %>% 
    dplyr::select(moved_from = NAME, race_ethnicity, max_income, pop_pre, avg_size) 
  
  # select columns from housing_town to join new units with weights by post-move town
  dev_town <- housing_town %>% 
    left_join(AMI_limits, by = "NAME") %>% 
    dplyr::select(moved_to = NAME, new_units, limit) 
  
  weights_pre <- weights %>%
    # join pre-move town households
    left_join(families_pre_town, by = "moved_from") %>% 
    # join post-move new units
    left_join(dev_town, by = "moved_to") %>% 
    # join new distances to replace distances of 0 for same-town pairs
    left_join(dist_0, by = c("moved_from", "moved_to")) %>% 
    # replace zero-length distances with replacements (half of max distance w/in town)
    mutate(dist_new = if_else(distance == 0, replace, distance)) %>% 
    # filter for only eligible moves
    filter(max_income <= limit) %>% 
    # collapse income brackets
    group_by(moved_from, moved_to, race_ethnicity) %>% 
    summarize(pop_pre = sum(pop_pre), avg_size = first(avg_size), new_units = first(new_units), dist_new = first(dist_new)) %>% 
    ungroup() %>% 
    # group by town pre-move and group
    group_by(moved_from, race_ethnicity) %>% 
    # round up units to integer
    mutate(prop_units = new_units/sum(new_units)) %>%
    # take inverse of distance squared
    mutate(inv_sq = 1/dist_new^2) %>% 
    # calculate crude penalty using families (pre-move town), proportion of new units (post-move town), and inverse distance squared (town pair)
    mutate(penalty = prop_units*inv_sq) %>% 
    # turn penalty into weighted probability that sums to 1 for any given starting town and group
    mutate(lambda = penalty/sum(penalty)) %>% 
    ungroup() 
  
  weights_post <- weights_pre %>%  
    # arrange by starting town
    arrange(moved_from) %>% 
    # create column for population moved in each town pair
    mutate(pop_moved = 0)
  
  return(weights_post)
}

Merge_data <- function(weights_post, mortality, exposures_town) {
  
  # annual average PM2.5, O3, NO2, NDVI, noise
  exposures_town_annual <- exposures_town[, -10:-162] %>% 
    st_drop_geometry()
  
  post_sim_1 <- weights_post %>% 
    # collapse income brackets
    group_by(NAME = moved_from, moved_to, race_ethnicity) %>% 
    summarize(pop_moved = sum(pop_moved), avg_size = first(avg_size)) %>% 
    ungroup() %>% 
    # join baseline all-cause mortality rates (per 100,000, currently using statewide)
    left_join(mortality, by = "race_ethnicity") %>% 
    # join annual average exposures pre-move
    left_join(exposures_town_annual, by = "NAME") %>% 
    rename(pm_pre = pm_mean, o3_pre = o3_mean, no2_pre = no2_mean, ndvi_pre = ndvi_mean, noise_pre = noise_mean, moved_from = NAME, NAME = moved_to) %>% 
    # join annual average exposures post-move
    left_join(exposures_town_annual, by = "NAME") %>% 
    rename(pm_post = pm_mean, o3_post = o3_mean, no2_post = no2_mean, ndvi_post = ndvi_mean, noise_post = noise_mean, moved_to = NAME) %>% 
    ## calculate difference (post-pre) in each annual exposure
    mutate(pm_diff = pm_post - pm_pre, o3_diff = o3_post - o3_pre, no2_diff = no2_post - no2_pre, ndvi_diff = ndvi_post - ndvi_pre, noise_diff = noise_post - noise_pre) %>% 
    ## calculate term to multiply by attributable fraction during simulation
    mutate(factor_term = avg_size/100000*mr_pre)
  
  return(post_sim_1)
}

# process heat data
get_post_heat_1 <- function(year = "2019", weights_post, exposures_town, mortality_daily) {
  
  # daily maximum heat index
  exposures_town_daily <- exposures_town[, -3:-9] %>% 
    # pivot to long format
    pivot_longer(cols = 3:155, names_to = "date", values_to = "heatindex")
  
  post_heat_1 <- weights_post %>% 
    # collapse income brackets
    rename(NAME = moved_from) %>% 
    group_by(NAME, moved_to, race_ethnicity) %>% 
    summarize(pop_moved = sum(pop_moved), avg_size = first(avg_size)) %>% 
    ungroup() %>% 
    # join pre-move heat index by town for each day
    full_join(exposures_town_daily, by = "NAME") %>% 
    rename(moved_from = NAME, NAME = moved_to, heatindex_pre = heatindex) %>% 
    # join post-move heat index by town and day
    left_join(exposures_town_daily, by = c("NAME", "date")) %>% 
    rename(moved_to = NAME) %>% 
    # convert pre/post heat index to Fahrenheit
    mutate(heatindex_pre = heatindex_pre*9/5+32, heatindex = heatindex*9/5+32) %>% 
    mutate(heatindex = if_else(heatindex_pre >= 75 & heatindex < 75, 75, heatindex)) %>%
    mutate(heatindex_pre = if_else(heatindex_pre < 75 & heatindex >= 75, 75, heatindex_pre)) %>% 
    mutate(fahr_diff = if_else(heatindex_pre < 75 & heatindex < 75, 0, (heatindex - heatindex_pre))) %>% 
    # pull and reformat month to join baseline all-cause daily mortality rates
    mutate(Month.Code = str_c(year, "/", substr(date, 1, 2), sep = "")) %>% 
    # join daily mortality rates (per 100,000, currently using statewide)
    left_join(mortality_daily, by = c("race_ethnicity", "Month.Code")) %>% 
    group_by(race_ethnicity, moved_from, moved_to) %>% 
    summarize(fahr_diff_times_rate = sum(fahr_diff*daily_rate_per100000), pop_moved = first(pop_moved), avg_size = first(avg_size)) %>% 
    mutate(factor_heat = avg_size/100000*(fahr_diff_times_rate/5))
  
  return(post_heat_1)
}

# simulation
simulate_cf <- function(weights_post, post_sim_1, post_heat_1, n = 2) {
  # assign dose-response relationships / health impact functions for all-cause mortality
  beta_pm_mean <- log(1.12) # beta per 10 μg/m3 increase, Pope et al., 2019
  SE_pm <- (beta_pm_mean - log(1.08))/1.96 # standard error
  
  beta_o3_mean <- log(1.02) # beta per 10 ppb inc., Turner et al., 2016, adjusted for PM2.5, NO2
  SE_o3 <- (beta_o3_mean - log(1.01))/1.96
  
  beta_no2_mean <- log(1.06) # beta per 10 ppb increase, Huang et al., 2021
  SE_no2 <- (beta_no2_mean - log(1.04))/1.96
  
  rr_heat_mean <- 1.014 # rate ratio per 5 F increase (lag 0), Wellenius et al., 2017
  SE_heat <- (rr_heat_mean - 1.004)/1.96
  
  beta_ndvi_mean <- log(0.96) # hazard ratio per 0.1 NDVI increase, Rojas-Rueda et al., 2019
  SE_ndvi <- (beta_ndvi_mean - log(0.94))/1.96
  
  beta_noise_mean <- log(1.05) # hazard ratio per 10 dB increase, Hao et al., 2022
  SE_noise <- (beta_noise_mean - log(1.02))/1.96
  
  simulate_once <- function() {
    
    # create vector of race/ethnic groups present in summary data
    race_eth_2 <- weights_post %>% pull(race_ethnicity) %>% unique()
    # initialize list for weights and town pairs for each racial/ethnic grouping
    weights_by_race <- list()
    
    # simulate movement of households into new housing
    for(i in 1:length(race_eth_2)) {
      ## filter town pairs and weights for each racial/ethnic grouping
      weights_by_race[[i]] <- weights_post %>% filter(race_ethnicity == race_eth_2[i])
      ## draw iteratively using binomial distribution from pool of households by starting town
      for(j in 0:168) {
        weights_by_race[[i]]$pop_moved[169*j+1] <- rbinom(1, size = weights_by_race[[i]]$pop_pre[169*j+1], prob = weights_by_race[[i]]$lambda[169*j+1])
        for(k in 2:169) {
          weights_by_race[[i]]$pop_moved[169*j+k] <- rbinom(1, size = max(0, weights_by_race[[i]]$pop_pre[169*j+k] - sum(weights_by_race[[i]]$pop_moved[(169*j+1):(169*j+k)])), prob = weights_by_race[[i]]$lambda[169*j+k])
        }
      }
    }
    
    # unnest list to data frame
    weights_post_1 <- bind_rows(weights_by_race) %>% 
      dplyr::select(moved_from, moved_to, race_ethnicity, eligible = pop_pre, pop_moved)
    
    # summarize families moved out by town
    moved_out <- weights_post_1 %>% 
      group_by(moved_from, race_ethnicity) %>%
      summarize(pop_pre = max(eligible), pop_out = sum(pop_moved)) %>% 
      rename(NAME = moved_from)
    
    # summarize families moved in by town
    moved_in <- weights_post_1 %>% 
      group_by(moved_to, race_ethnicity) %>% 
      summarize(pop_in = sum(pop_moved)) %>% 
      rename(NAME = moved_to)
    
    # join moved_out and moved_in by town
    pop_town_moved <- moved_out %>% 
      left_join(moved_in, by = c("NAME", "race_ethnicity")) %>% 
      #subtract families moved out and add families moved in by town
      mutate(pop_post = pop_pre + pop_in - pop_out) %>% 
      arrange(race_ethnicity, NAME) %>% 
      pull(pop_post)
    
    # draw random coefficient or RR from Normal distribution
    beta_pm <- rnorm(1, mean = beta_pm_mean, sd = SE_pm)
    beta_o3 <- rnorm(1, mean = beta_o3_mean, sd = SE_o3)
    beta_no2 <- rnorm(1, mean = beta_no2_mean, sd = SE_no2)
    rr_heat <- rnorm(1, mean = rr_heat_mean, sd = SE_heat)
    beta_ndvi <- rnorm(1, mean = beta_ndvi_mean, sd = SE_ndvi)
    beta_noise <- rnorm(1, mean = beta_noise_mean, sd = SE_noise)
    
    # estimate change in exposure and incidence
    post_sim <- post_sim_1 %>%
      # join population moved by town pair
      left_join(weights_post_1, by = c("moved_from", "moved_to", "race_ethnicity")) %>% 
      # Calculate deaths averted by multiplying by attributable fraction
      mutate(averted_pm = pop_moved.y*factor_term*(1-exp(pm_diff/10*beta_pm))) %>% 
      mutate(averted_o3 = pop_moved.y*factor_term*(1-exp(o3_diff/10*beta_o3))) %>% 
      mutate(averted_no2 = pop_moved.y*factor_term*(1-exp(no2_diff/10*beta_no2))) %>% 
      mutate(averted_ndvi = pop_moved.y*factor_term*(1-exp(ndvi_diff/0.1*beta_ndvi))) %>% 
      mutate(averted_noise = pop_moved.y*factor_term*(1-exp(noise_diff/10*beta_noise)))
    
    post_heat <- post_heat_1 %>% 
      left_join(weights_post_1, by = c("moved_from", "moved_to", "race_ethnicity")) %>% 
      mutate(averted_heat = pop_moved.y*factor_heat*(1-rr_heat)) %>% 
      group_by(moved_from, moved_to, race_ethnicity) %>% 
      summarize(averted_heat = sum(averted_heat), pop_moved = first(pop_moved.y))
    
    # summarize deaths averted by ethnoracial group
    averted_by_race <- post_sim %>%
      # join deaths averted from heat
      left_join(post_heat, by = c("moved_to", "moved_from", "race_ethnicity")) %>% 
      group_by(race_ethnicity) %>% 
      # sum deaths averted within each ethnoracial group
      summarize(across(c(averted_pm, averted_o3, averted_no2, averted_heat, averted_ndvi, averted_noise), ~ sum(.x))) %>% 
      pivot_longer(values_to = "averted", names_to = "exposure", cols = starts_with("averted")) %>%
      # arrange alphabetically
      arrange(race_ethnicity, exposure) %>% 
      # store as vector
      pull(averted)
    
    total_number_moved <- weights_post_1 %>% 
      filter(moved_from != moved_to) %>% 
      pull(pop_moved) %>% 
      sum()
    
    averted_by_race_pop_town_moved <- c(averted_by_race, pop_town_moved, total_number_moved)
    
    return(averted_by_race_pop_town_moved)
  }
  
  # set seed for reproducible samples  
  set.seed(99814)
  
  # repeat simulation x times
  all_sims <- future_replicate(n, simulate_once(), simplify = "array")
  
  return(all_sims)
  
}

exhibit3 <- function(all_sims) {
  
  # assign labels
  group_labels <- c("Asian", "Non-Hispanic Black", "Hispanic/Latino", "Non-Hispanic White")
  exposure_labels <- c("Heat", "NDVI", "NO2", "Noise", "O3", "PM2.5")
  
  # average estimates across all simulations, by pollutant and ethnoracial group
  all_sims_avg <- all_sims %>%
    unlist() %>%
    .[1:24,] %>% 
    data.frame(sample = ., race_ethnicity = rep(group_labels, each = 6), exposure = rep(exposure_labels, times = 4)) %>% 
    # convert to long format
    pivot_longer(cols = starts_with("sample"), names_to = "sample", names_prefix = "sample.", values_to = "averted") %>% 
    group_by(race_ethnicity, exposure) %>%
    # calculate mean and ME for confidence interval bounds
    summarize(averted_mean = mean(averted), averted_ME = sd(averted)*1.96) %>% 
    mutate(cl = averted_mean - averted_ME, cu = averted_mean + averted_ME)
  
  # average estimates across all simulations, by pollutant (collapsing ethnoracial groups)
  all_sims_totals <- all_sims %>% 
    unlist() %>%
    .[1:24,] %>% 
    data.frame(sample = ., race_ethnicity = rep(group_labels, each = 6), exposure = rep(exposure_labels, times = 4)) %>% 
    # convert to long format
    pivot_longer(cols = starts_with("sample"), names_to = "sample", names_prefix = "sample.", values_to = "averted") %>%
    group_by(exposure, sample) %>% 
    # sum deaths averted across ethnoracial groups from each pollutant for each sample
    summarize(averted = sum(averted)) %>% 
    ungroup() %>% 
    group_by(exposure) %>%
    # calculate mean and ME for confidence interval bounds
    summarize(averted_mean = mean(averted), averted_ME = sd(averted)*1.96) %>% 
    mutate(cl = averted_mean - averted_ME, cu = averted_mean + averted_ME)
  
  # total deaths averted across groups
  total_tbl <- all_sims_totals %>% 
    mutate(race_ethnicity = "Total")
  
  # deaths averted by ethnoracial group
  stratified_tbl <- all_sims_avg %>%
    bind_rows(total_tbl) %>% 
    mutate(across(averted_mean:cu, ~ round(.x, 0))) 
  
  # format final table of deaths averted by group and totals
  exhibit3 <- stratified_tbl %>% 
    mutate(averted = paste0(averted_mean, " (", cl, ", ", cu, ")")) %>%
    dplyr::select(`Ethnoracial Group` = race_ethnicity, exposure, averted) %>% 
    pivot_wider(names_from = exposure, values_from = averted)

  return(exhibit3)
}

get_avg_moved <- function(all_sims, town_polygon_sf) {
  
  # average number of households moved to new town
  avg_number_moved <- all_sims %>% 
    unlist() %>% 
    .[701,] %>% 
    mean()
  
  # average proportion of baseline households that moved to new town
  avg_prop_moved <- avg_number_moved/sum(town_polygon_sf$pop_pre)
  
  avg_moved <- c(avg_number_moved, avg_prop_moved)
  
  return(avg_moved)
}

exhibit4 <- function(town_polygon_sf, all_sims) {
  
  # assign labels
  group_labels <- c("Asian", "Non-Hispanic Black", "Hispanic/Latino", "Non-Hispanic White")
  exposure_labels <- c("Heat", "NDVI", "NO2", "Noise", "O3", "PM2.5")
  
  # average estimates across all simulations, by pollutant and ethnoracial group
  all_sims_avg <- all_sims %>%
    unlist() %>%
    .[1:24,] %>% 
    data.frame(sample = ., race_ethnicity = rep(group_labels, each = 6), exposure = rep(exposure_labels, times = 4)) %>% 
    # convert to long format
    pivot_longer(cols = starts_with("sample"), names_to = "sample", names_prefix = "sample.", values_to = "averted") %>% 
    group_by(race_ethnicity, exposure) %>%
    # calculate mean and ME for confidence interval bounds
    summarize(averted_mean = mean(averted), averted_ME = sd(averted)*1.96) %>% 
    mutate(cl = averted_mean - averted_ME, cu = averted_mean + averted_ME)
  
  # average estimates across all simulations, by pollutant (collapsing ethnoracial groups)
  all_sims_totals <- all_sims %>% 
    unlist() %>%
    .[1:24,] %>% 
    data.frame(sample = ., race_ethnicity = rep(group_labels, each = 6), exposure = rep(exposure_labels, times = 4)) %>% 
    # convert to long format
    pivot_longer(cols = starts_with("sample"), names_to = "sample", names_prefix = "sample.", values_to = "averted") %>%
    group_by(exposure, sample) %>% 
    # sum deaths averted across ethnoracial groups from each pollutant for each sample
    summarize(averted = sum(averted)) %>% 
    ungroup() %>% 
    group_by(exposure) %>%
    # calculate mean and ME for confidence interval bounds
    summarize(averted_mean = mean(averted), averted_ME = sd(averted)*1.96) %>% 
    mutate(cl = averted_mean - averted_ME, cu = averted_mean + averted_ME)
  
  # total deaths averted across groups
  total_tbl <- all_sims_totals %>% 
    mutate(race_ethnicity = "Total")  
  
  # collapse income brackets - sf with low-income population by town and ethnoracial group
  town_polygon_sf_li <- town_polygon_sf %>%
    st_drop_geometry() %>% 
    group_by(GEOID, NAME, race_ethnicity) %>% 
    summarize(pop_pre = sum(pop_pre), avg_size = first(avg_size))
  
  # population by race/ethnicity
  pop_by_ethnoracial <- town_polygon_sf_li %>%
    mutate(pop_li = pop_pre*avg_size) %>% 
    group_by(race_ethnicity) %>% 
    summarize(population_li = sum(pop_li)) %>% 
    data.frame() %>% 
    adorn_totals()
  
  # rates of deaths averted
  exhibit4 <- all_sims_avg %>%
    bind_rows(total_tbl) %>% 
    mutate(race_ethnicity = case_match(race_ethnicity, "Asian" ~ "asian", "Hispanic/Latino" ~ "hispanic", "Non-Hispanic Black" ~ "black_nh", "Non-Hispanic White" ~ "white_nh", "Total" ~ "Total")) %>% 
    left_join(pop_by_ethnoracial, by = "race_ethnicity") %>% 
    mutate(across(.cols = c(averted_mean, cl, cu), ~ round(.x/population_li*100000, 2))) %>% 
    mutate(averted_rate = paste0(averted_mean, " (", cl, ", ", cu, ")")) %>%
    dplyr::select(`Ethnoracial Group` = race_ethnicity, exposure, averted_rate) %>% 
    pivot_wider(names_from = exposure, values_from = averted_rate) %>% 
    mutate(`Ethnoracial Group` = case_match(`Ethnoracial Group`, "asian" ~ "Asian", "hispanic" ~ "Hispanic/Latino", "black_nh" ~ "Non-Hispanic Black", "white_nh" ~ "Non-Hispanic White", "Total" ~ "Total"))

  return(exhibit4)
  
}

tableS2 <- function(town_polygon_sf, all_sims, exposures_town) {
  
  # collapse max_income
  town_polygon_sf_2 <- town_polygon_sf %>%
    group_by(GEOID, NAME, race_ethnicity, percent_assisted, at_least_10, new_units) %>% 
    summarize(pop_pre = sum(pop_pre)) %>% 
    ungroup()
  
  # average population moved over all simulation runs
  sim_moved_avg <- all_sims %>%
    unlist() %>%
    .[25:700,] %>% 
    rowMeans()
  
  # join average population moved to town polygon
  sim_moved_totals <- town_polygon_sf_2 %>% 
    arrange(race_ethnicity, NAME) %>% 
    mutate(pop_post = round(sim_moved_avg, 0))
  
  # collapse ethnoracial groupings
  sim_moved_totals_2 <- sim_moved_totals %>% 
    group_by(GEOID, NAME, percent_assisted, at_least_10, new_units) %>% 
    summarize(pop_pre = sum(pop_pre), pop_post = sum(pop_post))
  
  # join pollutants by town
  pollutants_prepost <- sim_moved_totals_2 %>%
    left_join(st_drop_geometry(exposures_town), by = "NAME") 
  
  # make vectors to calculate mean and percentiles - pre-move
  exposures_pre <- data.frame(pm = rep(x = pollutants_prepost$pm_mean, times = pollutants_prepost$pop_pre)) %>% 
    mutate(o3 = rep(x = pollutants_prepost$o3_mean, times = pollutants_prepost$pop_pre)) %>% 
    mutate(no2 = rep(x = pollutants_prepost$no2_mean, times = pollutants_prepost$pop_pre)) %>% 
    mutate(heat = rep(x = pollutants_prepost$heat_annual_mean, times = pollutants_prepost$pop_pre)) %>% 
    mutate(ndvi = rep(x = pollutants_prepost$ndvi_mean, times = pollutants_prepost$pop_pre)) %>% 
    mutate(noise = rep(x = pollutants_prepost$noise_mean, times = pollutants_prepost$pop_pre)) 
  
  exposure_pre_summary <- exposures_pre %>% 
    summarize(across(.cols = c(pm, o3, no2, ndvi, heat, noise), list(mean = mean, q25 = function(x) stats::quantile(x, probs = 0.25), q75 = function(x) stats::quantile(x, probs = 0.75)))) %>% 
    mutate(when = "pre")
  
  # make vectors to calculate mean and percentiles - post-move
  exposures_post <- data.frame(pm = rep(x = pollutants_prepost$pm_mean, times = pollutants_prepost$pop_post)) %>% 
    mutate(o3 = rep(x = pollutants_prepost$o3_mean, times = pollutants_prepost$pop_post)) %>% 
    mutate(no2 = rep(x = pollutants_prepost$no2_mean, times = pollutants_prepost$pop_post)) %>% 
    mutate(heat = rep(x = pollutants_prepost$heat_annual_mean, times = pollutants_prepost$pop_post)) %>% 
    mutate(ndvi = rep(x = pollutants_prepost$ndvi_mean, times = pollutants_prepost$pop_post)) %>%
    mutate(noise = rep(x = pollutants_prepost$noise_mean, times = pollutants_prepost$pop_post)) 
  
  exposure_post_summary <- exposures_post %>%
    summarize(across(.cols = c(pm, o3, no2, ndvi, heat, noise), list(mean = mean, q25 = function(x) stats::quantile(x, probs = 0.25), q75 = function(x) stats::quantile(x, probs = 0.75)))) %>% 
    mutate(when = "post")
  
  # bind rows
  pollutants_prepost_tbl <- bind_rows(exposure_pre_summary, exposure_post_summary) %>% 
    mutate(across(pm_mean:noise_q75, ~ round(.x, 2))) %>% 
    pivot_longer(cols = pm_mean:noise_q75) %>% 
    pivot_wider(names_from = when) %>% 
    mutate(exposure = sub(x = name, pattern = regex("_.+"), replacement = "")) %>% 
    mutate(quartile = sub(x = name, pattern = regex(".+_"), replacement = "")) %>%
    dplyr::select(-name) %>% 
    pivot_wider(names_from = quartile, values_from = c("pre", "post")) %>% 
    mutate(pre = paste0(pre_mean, " (", pre_q25, ", ", pre_q75, ")"), post = paste0(post_mean, " (", post_q25, ", ", post_q75, ")")) %>% 
    dplyr::select(exposure, pre, post)
    
  return(pollutants_prepost_tbl)
}

tableS1 <- function(all_sims, town_polygon_sf) {
  
  # average population moved over all simulation runs
  sim_moved_avg <- all_sims %>%
    unlist() %>%
    .[25:700,] %>% 
    rowMeans()
  
  race_town <- get_decennial(geography = "county subdivision", state = "CT", variable = c("P2_006N", "P2_005N", "P2_002N", "P1_006N"), year = 2020, output = "wide") %>% 
    rename(black_nh = P2_006N, white_nh = P2_005N, hispanic = P2_002N, asian = P1_006N) %>%
    filter(!str_detect(NAME, "County subdivisions not defined")) %>% 
    mutate(county = sub(pattern = ".*town, ", replacement = "", x = NAME), .after = NAME) %>% 
    mutate(county = sub(pattern = " County, Connecticut", replacement = "", x = county))
  
  D_multi_pre <- map(unique(race_town$county), ~ filter(race_town, county == .x) %>% 
                       dplyr::select(black_nh, white_nh, hispanic, asian) %>% 
                       DMulti()) %>% 
    unlist() %>% 
    data.frame(D_multi_pre = ., NAME = unique(race_town$county))
  
  pop_change_town <- town_polygon_sf %>%
    group_by(GEOID, NAME, race_ethnicity, percent_assisted, at_least_10, new_units) %>% 
    summarize(pop_pre = sum(pop_pre)) %>% 
    ungroup() %>% 
    arrange(race_ethnicity, NAME) %>% 
    mutate(pop_post = round(sim_moved_avg, 0)) %>% 
    st_drop_geometry() %>% 
    mutate(pop_change = pop_post - pop_pre) %>% 
    dplyr::select(GEOID, race_ethnicity, pop_change) %>% 
    pivot_wider(names_from = "race_ethnicity", values_from = "pop_change") %>% 
    left_join(race_town, by = "GEOID") %>% 
    mutate(asian_post = asian.x + asian.y, black_nh_post = black_nh.x + black_nh.y, hispanic_post = hispanic.x + hispanic.y, white_nh_post = white_nh.x + white_nh.y) %>% 
    dplyr::select(GEOID, NAME, asian_post, black_nh_post, hispanic_post, white_nh_post) %>% 
    mutate(county = sub(pattern = ".*town, ", replacement = "", x = NAME), .after = NAME) %>% 
    mutate(county = sub(pattern = " County, Connecticut", replacement = "", x = county))
  
  D_multi_post <- map(unique(pop_change_town$county), ~ filter(pop_change_town, county == .x) %>%
                        dplyr::select(asian_post, black_nh_post, hispanic_post, white_nh_post) %>%
                        DMulti()) %>% 
    unlist() %>% 
    data.frame(D_multi_post = ., NAME = unique(pop_change_town$county))
  
  D_multi_sf <- get_acs(state = "CT", geography = "county", variables = "B01001_001", geometry = TRUE, year = 2019) %>%
    mutate(NAME = sub(pattern = " County, Connecticut", replacement = "", x = NAME)) %>% 
    left_join(D_multi_pre, by = "NAME") %>% 
    left_join(D_multi_post, by = "NAME")
    
  # for Table S1
  D_pre_post_tbl <- D_multi_sf %>% 
    st_drop_geometry() %>% 
    dplyr::select(NAME, D_multi_pre, D_multi_post) %>% 
    mutate(across(c(D_multi_pre, D_multi_post), ~ round(.x, 2)))
  
  return(D_pre_post_tbl)
  
}

figureS7 <- function(tableS1) {
  
  D_multi_sf <- get_acs(state = "CT", geography = "county", variables = "B01001_001", geometry = TRUE, year = 2019) %>%
    mutate(NAME = sub(pattern = " County, Connecticut", replacement = "", x = NAME)) %>% 
    left_join(tableS1, by = "NAME")
  
  D_pre <- ggplot() + geom_sf(data = D_multi_sf, aes(fill = D_multi_pre)) + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 0.6)) + labs(fill = "Multigroup D", title = "A. Index of dissimilarity in 2019") + theme_minimal() + theme(plot.title = element_text(size = 10))
  
  legend_D_multi <- get_legend(D_pre) %>% as_ggplot()
  
  D_pre <- D_pre + theme(legend.position = "none")
  
  D_post <- ggplot() + geom_sf(data = D_multi_sf, aes(fill = D_multi_post)) + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 0.6)) + labs(fill = "Multigroup D", title = "B. Index of dissimilarity post-simulations") + theme_minimal() + theme(legend.position = "none") + theme(legend.position = "none") + theme(plot.title = element_text(size = 10))
  
  D_pre_post_map <- grid.arrange(D_pre, D_post, legend_D_multi, nrow = 1, widths = c(3,3,1))
  
  return(D_pre_post_map)
}  

figureS1 <- function(exposures_town) {
  
  # map 2019 annual mean PM2.5 by town
  pm_2019_map <- ggplot() + geom_sf(data = exposures_town, aes(fill = pm_mean)) + scale_fill_viridis_c(option = "magma", direction = -1) + labs(fill = "PM2.5 (μg/m3)") + theme_minimal()
  
  # map 2019 annual mean ozone by town
  o3_2019_map <- ggplot() + geom_sf(data = exposures_town, aes(fill = o3_mean)) + scale_fill_viridis_c(option = "magma", direction = -1) + labs(fill = "O3 (ppb)") + theme_minimal()
  
  # map 2019 annual mean NO2 by town
  no2_2019_map <- ggplot() + geom_sf(data = exposures_town, aes(fill = no2_mean)) + scale_fill_viridis_c(option = "magma", direction = -1) + labs(fill = "NO2 (ppb)") + theme_minimal()
  
  # map 2019 mean heat index during 5/1-9/30 by town
  heat_2019_map <- ggplot() + geom_sf(data = exposures_town, aes(fill = heat_annual_mean)) + scale_fill_viridis_c(option = "magma", direction = -1) + labs(fill = "Heat Index (°C)") + theme_minimal()
  
  # map 2019 annual mean NDVI by town
  ndvi_2019_map <- ggplot() + geom_sf(data = exposures_town, aes(fill = ndvi_mean)) + scale_fill_viridis_c(option = "magma", direction = -1) + labs(fill = "NDVI") + theme_minimal()
  
  # map 2019 annual mean noise by town
  noise_2019_map <- ggplot() + geom_sf(data = exposures_town, aes(fill = noise_mean)) + scale_fill_viridis_c(option = "magma", direction = -1) + labs(fill = "Noise (dB)") + theme_minimal()
  
  figureS1 <- list(pm_2019_map, o3_2019_map, no2_2019_map, heat_2019_map, ndvi_2019_map, noise_2019_map)
  
  return(figureS1)
  
}

figureS3_S6 <- function(all_sims, town_polygon_sf) {
  
  # collapse max_income
  town_polygon_sf_2 <- town_polygon_sf %>%
    group_by(GEOID, NAME, race_ethnicity, percent_assisted, at_least_10, new_units) %>% 
    summarize(pop_pre = sum(pop_pre)) %>% 
    ungroup()
  
  # average population moved over all simulation runs
  sim_moved_avg <- all_sims %>%
    unlist() %>%
    .[25:700,] %>% 
    rowMeans()
  
  # join average population moved to town polygon
  sim_moved_totals <- town_polygon_sf_2 %>% 
    arrange(race_ethnicity, NAME) %>% 
    mutate(pop_post = round(sim_moved_avg, 0))
  
  # round up maximum population in any one town for each group to set upper limit for pre/post map legends
  find_map_limit <- function(group) {
    pre_max <- sim_moved_totals %>% filter(race_ethnicity == group) %>% pull(pop_pre) %>% max() 
    post_max <- sim_moved_totals %>% filter(race_ethnicity == group) %>% pull(pop_post) %>% max() 
    
    max_limit <- (max(pre_max, post_max)/100) %>% ceiling()*100
    return(max_limit)
  }
  
  asian_max <- find_map_limit(group = "asian")
  hispanic_max <- find_map_limit(group = "hispanic")
  black_nh_max <- find_map_limit(group = "black_nh")
  white_nh_max <- find_map_limit(group = "white_nh")
  
  # Asian households
  asian_pre_map <- ggplot() + geom_sf(data = town_polygon_sf_2[town_polygon_sf_2$race_ethnicity == "asian",], aes(fill = pop_pre)) + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, asian_max)) + labs(fill = "Households") + theme_minimal() # + theme(legend.position = "None",)
  
  asian_post_map <- ggplot() + geom_sf(data = sim_moved_totals[sim_moved_totals$race_ethnicity == "asian",], aes(fill = pop_post)) + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, asian_max)) + labs(fill = "Households") + theme_minimal() # + theme(legend.direction = "horizontal", legend.position = "bottom")
  
  # asian_pre_post_map <- grid.arrange(asian_pre_map, asian_post_map, nrow = 2)
  
  # Non-Hispanic Black households
  black_nh_pre_map <- ggplot() + geom_sf(data = town_polygon_sf_2[town_polygon_sf_2$race_ethnicity == "black_nh",], aes(fill = pop_pre)) + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, black_nh_max)) + labs(fill = "Households") + theme_minimal()
  
  black_nh_post_map <- ggplot() + geom_sf(data = sim_moved_totals[sim_moved_totals$race_ethnicity == "black_nh",], aes(fill = pop_post)) + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, black_nh_max)) + labs(fill = "Households") + theme_minimal()
  
  # black_nh_pre_post_map <- grid.arrange(black_nh_pre_map, black_nh_post_map, nrow = 1)
  
  # Hispanic/Latino households
  hispanic_pre_map <- ggplot() + geom_sf(data = town_polygon_sf_2[town_polygon_sf_2$race_ethnicity == "hispanic",], aes(fill = pop_pre)) + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, hispanic_max)) + labs(fill = "Households") + theme_minimal()
  
  hispanic_post_map <- ggplot() + geom_sf(data = sim_moved_totals[sim_moved_totals$race_ethnicity == "hispanic",], aes(fill = pop_post)) + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, hispanic_max)) + labs(fill = "Households") + theme_minimal()
  
  # hispanic_pre_post_map <- grid.arrange(hispanic_pre_map, hispanic_post_map, nrow = 1)
  
  # Non-Hispanic White households
  white_nh_pre_map <- ggplot() + geom_sf(data = town_polygon_sf_2[town_polygon_sf_2$race_ethnicity == "white_nh",], aes(fill = pop_pre)) + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, white_nh_max)) + labs(fill = "Households") + theme_minimal()
  
  white_nh_post_map <- ggplot() + geom_sf(data = sim_moved_totals[sim_moved_totals$race_ethnicity == "white_nh",], aes(fill = pop_post)) + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, white_nh_max)) + labs(fill = "Households") + theme_minimal()
  
  # white_nh_pre_post_map <- grid.arrange(white_nh_pre_map, white_nh_post_map, nrow = 1)
  
  figureS3_S6 <- list(asian_pre_map, asian_post_map, black_nh_pre_map, black_nh_post_map, hispanic_pre_map, hispanic_post_map, white_nh_pre_map, white_nh_post_map)
  
  return(figureS3_S6)
  
}

sim_moved_totals_2 <- function(all_sims, town_polygon_sf) {
  
  # collapse max_income
  town_polygon_sf_2 <- town_polygon_sf %>%
    group_by(GEOID, NAME, race_ethnicity, percent_assisted, at_least_10, new_units) %>% 
    summarize(pop_pre = sum(pop_pre)) %>% 
    ungroup()
  
  # average population moved over all simulation runs
  sim_moved_avg <- all_sims %>%
    unlist() %>%
    .[25:700,] %>% 
    rowMeans()
  
  # join average population moved to town polygon
  sim_moved_totals <- town_polygon_sf_2 %>% 
    arrange(race_ethnicity, NAME) %>% 
    mutate(pop_post = round(sim_moved_avg, 0))
  
  # collapse ethnoracial groupings
  sim_moved_totals_2 <- sim_moved_totals %>% 
    group_by(GEOID, NAME, percent_assisted, at_least_10, new_units) %>% 
    summarize(pop_pre = sum(pop_pre), pop_post = sum(pop_post)) %>% 
    mutate(cat_pre = case_when(pop_pre %in% 0:5499 ~ "<5,500", 
                               pop_pre %in% 5500:10999 ~ "5,500-10,999",
                               pop_pre %in% 11000:16499 ~ "11,000-16,499", 
                               pop_pre %in% 16500:21999 ~ "16,500-21,999", 
                               pop_pre %in% 22000:27499 ~ "22,000-27,499",
                               pop_pre %in% 27500:33000 ~ "27,500-32,999")) %>% 
    mutate(cat_post = case_when(pop_post %in% 0:5499 ~ "<5,500", 
                               pop_post %in% 5500:10999 ~ "5,500-10,999",
                               pop_post %in% 11000:16499 ~ "11,000-16,499", 
                               pop_post %in% 16500:21999 ~ "16,500-21,999", 
                               pop_post %in% 22000:27499 ~ "22,000-27,499",
                               pop_post %in% 27500:33000 ~ "27,500-32,999"))
  
  return(sim_moved_totals_2)
  
}

exhibit1 <- function(sim_moved_totals_2) {
  
  exhibit1 <- ggplot() + geom_sf(data = sim_moved_totals_2, aes(fill = cat_pre)) + scale_colour_brewer(limits = c("<5,500", "5,500-10,999", "11,000-16,499", "16,500-21,999", "22,000-27,499", "27,500-32,999"), 
    type = "seq", palette = "YlOrBr", aesthetics = "fill") + labs(fill = "Households", tag = "A.") + theme_minimal() + theme(plot.subtitle = element_text(size = 8))
  
  return(exhibit1)
  
}

exhibit1_data <- function(sim_moved_totals_2) {
  
  # spreadsheet for plotting data
  exhibit1_data <- sim_moved_totals_2 %>% 
    dplyr::select(GEOID, NAME, `Households`= pop_pre, `Range` = cat_pre) %>% 
    st_drop_geometry()
  
  return(exhibit1_data)
  
}

exhibit2 <- function(sim_moved_totals_2) {
  
  exhibit2 <- ggplot() + geom_sf(data = sim_moved_totals_2, aes(fill = cat_post)) + scale_colour_brewer(limits = c("<5,500", "5,500-10,999", "11,000-16,499", "16,500-21,999", "22,000-27,499", "27,500-32,999"), 
    type = "seq", palette = "YlOrBr", aesthetics = "fill") + labs(fill = "Households", tag = "B.") + theme_minimal() + theme(plot.subtitle = element_text(size = 8))
  
  return(exhibit2)
  
}

exhibit2_data <- function(sim_moved_totals_2) {
  
  # spreadsheet for plotting data
  exhibit2_data <- sim_moved_totals_2 %>% 
    dplyr::select(GEOID, NAME, `Households` = pop_post, `Range` = cat_post) %>% 
    st_drop_geometry()
  
  return(exhibit2_data)
  
}

figureS2 <- function(town_polygon_sf) {
  
  # collapse max_income
  town_polygon_sf_2 <- town_polygon_sf %>%
    group_by(GEOID, NAME, race_ethnicity, percent_assisted, at_least_10, new_units) %>% 
    summarize(pop_pre = sum(pop_pre)) %>% 
    ungroup()
  
  # map percent of affordable housing units by town
  percent_affordable <- ggplot() + geom_sf(data = town_polygon_sf_2, aes(fill = percent_assisted)) + scale_fill_viridis_c(option = "magma", direction = -1) + labs(fill = "Percent Affordable") + theme_minimal()
  
  # map 8-30g exemptions by town
  exemptions_town <- ggplot() + geom_sf(data = town_polygon_sf_2, aes(fill = at_least_10)) + labs(fill = "At Least 10%?") + scale_colour_brewer(type = "div", aesthetics = "fill", palette = "Set2", direction = -1) + theme_minimal()
  
  figureS2 <- list(percent_affordable, exemptions_town)
  
  return(figureS2)
  
}

simulate_cf_s1 <- function(weights, housing_town, low_income_pop, AMI_limits, post_sim_1, post_heat_1, CT_townships, n = 2) {
  
    # initialize data frame for same-town pairs
    dist_0 <- weights %>% 
      # filter for same-town pairs (staying in town)
      filter(moved_from == moved_to) %>% 
      dplyr::select(moved_from, moved_to)
    
    # replace 0 distance with estimated town radius
    for(i in 1:nrow(dist_0)) {
      dist_0$replace[i] <- CT_townships %>%
        filter(NAME == dist_0$moved_from[i]) %>% 
        # retrieve convex hull
        st_convex_hull() %>% 
        # cast to multipoint
        st_cast("MULTIPOINT") %>%
        # cast to point
        st_cast("POINT") %>%
        # calculate distances between points
        st_distance() %>%
        #find max distance between any two points on town convex hull and divide by two
        max(na.rm = TRUE)/2
    }
    
    # select columns from low_income_pop to join households by pre-move town with weights
    families_pre_town <- low_income_pop %>% 
      dplyr::select(moved_from = NAME, race_ethnicity, max_income, pop_pre, avg_size) 
    
    # select columns from housing_town to join new units with weights by post-move town
    dev_town <- housing_town %>% 
      left_join(AMI_limits, by = "NAME") %>% 
      dplyr::select(moved_to = NAME, new_units, limit) 
    
    weights_post_s1 <- weights %>%
      # join pre-move town households
      left_join(families_pre_town, by = "moved_from") %>% 
      # join post-move new units
      left_join(dev_town, by = "moved_to") %>% 
      # join new distances to replace distances of 0 for same-town pairs
      left_join(dist_0, by = c("moved_from", "moved_to")) %>% 
      # replace zero-length distances with replacements (half of max distance w/in town)
      mutate(dist_new = if_else(distance == 0, replace, distance)) %>% 
      # filter for only eligible moves
      filter(max_income <= limit) %>% 
      # collapse income brackets
      group_by(moved_from, moved_to, race_ethnicity) %>% 
      summarize(pop_pre = sum(pop_pre), avg_size = first(avg_size), new_units = first(new_units), dist_new = first(dist_new)) %>% 
      ungroup() %>% 
      # group by town pre-move and group
      group_by(moved_from, race_ethnicity) %>% 
      # round up units to integer
      mutate(prop_units = new_units/sum(new_units)) %>%
      # take inverse of distance squared - not using in this sensitivity analysis
      mutate(inv_sq = 1/dist_new^2) %>% 
      # turn penalty (no distance weighting) into weighted probability that sums to 1 for any starting town and group
      mutate(lambda = prop_units/sum(prop_units)) %>% 
      ungroup() %>%  
      # arrange by starting town
      arrange(moved_from) %>% 
      # create column for population moved in each town pair
      mutate(pop_moved = 0)
  
  # set seed for reproducible samples  
  set.seed(99814)
  
  # repeat simulation x times
  all_sims_s1 <- simulate_cf(weights_post = weights_post_s1, post_sim_1, post_heat_1, n)
  
  return(all_sims_s1)
}

simulate_cf_18 <- function(exposures_town_18, CT_townships, downloads, n) {
  
  # replace housing data
  housing_town_18 <- allocate_housing(year_h = "2018", downloads = downloads)
  
  # replace mortality data
  mortality_18 <- get_mortality(year = "2018")
  mortality_daily_18 <- get_mortality_daily(year = "2018", mortality = mortality_18)
  
  # replace AMI limits
  AMI_limits_18 <- Load_AMI_limits(year = "2018", hud_key = hud_key)
  
  # replace population data and average household sizes
  low_income_pop_18 <- get_low_income_pop(year = "2018", AMI_limits = AMI_limits_18)
  
  # re-run simulation
  town_polygon_sf_18 <- get_town_polygon_sf(CT_townships = CT_townships, housing_town = housing_town_18, low_income_pop = low_income_pop_18)
  
  weights_18 <- get_weights(CT_townships = CT_townships)
  
  weights_post_18 <- get_weights_post(weights = weights_18, CT_townships = CT_townships, housing_town = housing_town_18, low_income_pop = low_income_pop_18, AMI_limits = AMI_limits_18)
  
  post_sim_1_18 <- Merge_data(weights_post = weights_post_18, mortality = mortality_18, exposures_town = exposures_town_18)
  
  post_heat_1_18 <- get_post_heat_1(year = "2018", weights_post = weights_post_18, exposures_town = exposures_town_18, mortality_daily = mortality_daily_18)
  
  # set seed for reproducible samples
  set.seed(99814)
  
  # repeat simulation x times
  all_sims_18 <- simulate_cf(weights_post = weights_post_18, post_sim_1 = post_sim_1_18, post_heat_1 = post_heat_1_18, n = n)

  return(all_sims_18)

}

save_exhibits <- function(exhibit1, exhibit2, exhibit1_data, exhibit2_data, exhibit3, exhibit4, figureS1, figureS2, figureS3_S6, tableS1, figureS7, tableS2, exhibitA12, exhibitA13, avg_moved, avg_moved_s1, avg_moved_18) {
  
  if(!dir.exists(here("output"))){dir.create(here("output"))}
  
  # exhibit 1 and 2 - map of population, pre/post-simulation and underlying data
  ## exhibit 1
  ggsave(plot = exhibit1, here("output", "exhibit1.pdf"), width = 4.5, height = 3, units = "in")
  write.csv(exhibit1_data, here("output", "exhibit1_data.csv"))
  
  ## exhibit 2
  ggsave(plot = exhibit2, here("output", "exhibit2.pdf"), width = 4.5, height = 3, units = "in")
  write.csv(exhibit2_data, here("output", "exhibit2_data.csv"))
  
  ## exhibit 1 and 2, combined
  exhibit1_nolegend <- exhibit1 + theme(legend.position = "none")
  exhibit2_nolegend <- exhibit2 + theme(legend.position = "none")
  legend <- as_ggplot(get_legend(exhibit1))
  
  exhibit1_2 <- grid.arrange(exhibit1_nolegend, exhibit2_nolegend, legend, nrow = 1, widths = c(3,3,1))
  ggsave(plot = exhibit1_2, here("output", "exhibit1_2.pdf"), width = 10, height = 3, units = "in")
  
  # exhibit 3 - deaths averted
  write.csv(exhibit3, file = here("output", "exhibit3.csv"))
  
  # exhibit 4 - rates of deaths averted
  write.csv(exhibit4, here("output", "exhibit4.csv"))
  
  # exhibits A1a-f - map of baseline exposures by town
  ggsave(figureS1[[1]], filename = here("output", "exhibitA1a.png"), width = 4.5, height = 3, units = "in")
  ggsave(figureS1[[2]], filename = here("output", "exhibitA1b.png"), width = 4.5, height = 3, units = "in")
  ggsave(figureS1[[3]], filename = here("output", "exhibitA1c.png"), width = 4.5, height = 3, units = "in")
  ggsave(figureS1[[4]], filename = here("output", "exhibitA1d.png"), width = 4.5, height = 3, units = "in")
  ggsave(figureS1[[5]], filename = here("output", "exhibitA1e.png"), width = 4.5, height = 3, units = "in")
  ggsave(figureS1[[6]], filename = here("output", "exhibitA1f.png"), width = 4.5, height = 3, units = "in")
  
  # exhibit A2a-b - baseline percent affordable and >= 10% by town
  ggsave(figureS2[[1]], filename = here("output", "exhibitA2a.png"), width = 4.5, height = 3, units = "in")
  ggsave(figureS2[[2]], filename = here("output", "exhibitA2b.png"), width = 4.5, height = 3, units = "in")
  
  # exhibit A3 and exhibit A4 are text
  
  # exhibit A5a-b - map of Asian low-income households, pre/post-simulation
  ggsave(figureS3_S6[[1]], filename = here("output", "exhibitA5a.png"), width = 4.5, height = 3, units = "in")
  ggsave(figureS3_S6[[2]], filename = here("output", "exhibitA5b.png"), width = 4.5, height = 3, units = "in")
  
  # exhibit A6 - map of non-Hispanic Black low-income households, pre/post-simulation
  ggsave(figureS3_S6[[3]], filename = here("output", "exhibitA6a.png"), width = 4.5, height = 3, units = "in")
  ggsave(figureS3_S6[[4]], filename = here("output", "exhibitA6b.png"), width = 4.5, height = 3, units = "in")
  
  # exhibit A7 - map of Latino low-income households, pre/post-simulation
  ggsave(figureS3_S6[[5]], filename = here("output", "exhibitA7a.png"), width = 4.5, height = 3, units = "in")
  ggsave(figureS3_S6[[6]], filename = here("output", "exhibitA7b.png"), width = 4.5, height = 3, units = "in")
  
  # exhibit A8 - map of non-Hispanic White low-income households, pre/post-simulation
  ggsave(figureS3_S6[[7]], filename = here("output", "exhibitA8a.png"), width = 4.5, height = 3, units = "in")
  ggsave(figureS3_S6[[8]], filename = here("output", "exhibitA8b.png"), width = 4.5, height = 3, units = "in")
  
  # exhibit A9 - dissimilarity index by county, pre/post table
  write.csv(tableS1, here("output", "exhibitA9.csv"))
  
  # exhibit A10 - dissimilarity index by county, pre/post map
  ggsave(figureS7, filename = here("output", "exhibitA10.png"), width = 10, height = 3, units = "in")
  
  # exhibit A11 - exposures pre/post table
  write.csv(tableS2, here("output", "exhibitA11.csv"))
  
  # exhibit A12 - sensitivity analysis, no inverse distance weighting penalty
  write.csv(exhibitA12, here("output", "exhibitA12.csv"))
  
  # exhibit A13 - sensitivity analysis, using 2018 data
  write.csv(exhibitA13, here("output", "exhibitA13.csv"))
  
  # average moved - main analysis
  write.csv(avg_moved, here("output", "avg_moved.csv"))
  
  # average moved - sensitivity analysis 1
  write.csv(avg_moved_s1, here("output", "avg_moved_s1.csv"))
  
  # average moved - sensitivity analysis 2
  write.csv(avg_moved_18, here("output", "avg_moved_s2.csv"))
  
}
