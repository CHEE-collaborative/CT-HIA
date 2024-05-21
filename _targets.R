
library(targets)
library(tarchetypes)
library(future)

source("Functions.R")

plan(multisession, workers = availableCores() - 1)

options(tigris_use_cache = TRUE)

################ CUSTOMIZE THE FOLLOWING BEFORE RUNNING tar_make() #############
# required - define your Census API key here
census_key <- ""

# required - define your HUD API key here
hud_key <- ""

# optional - replace 10 with desired number of simulation repeats (1,000 for full simulation)
n <- 10

# optional - set to TRUE if you would like to skip sensitivity analyses
options(SKIP_SUPPL_ANALYSIS = FALSE)
################################################################################

tar_option_set(
  packages = c(
    "tidyverse",
    "tidycensus",
    "here",
    "sp",
    "sf",
    "readxl",
    "ncdf4",
    "daymetr",
    "tidyr",
    "tigris",
    "weathermetrics",
    "exactextractr",
    "terra",
    "stars",
    "humidity",
    "hudr",
    "janitor",
    "gridExtra",
    "future",
    "future.apply",
    "OasisR",
    "ggpubr",
    "osrm"
  ),
  format = 'rds', 
  workspace_on_error = TRUE)

list(
  tar_target(install_census_key, {
    census_api_key(key = census_key)
  }),
  tar_target(downloads, {
    download_data()
  }),
  tar_target(CT_townships, {
    Load_CT_townships()
  }),
  tar_target(exposures_town, {
    get_exposures_town(year = "2019", CT_townships = CT_townships, downloads = downloads)
  }),
  tar_target(housing_town, {
    allocate_housing(year_h = "2019", downloads = downloads)
  }),
  tar_target(mortality, {
    get_mortality(year = "2019")
  }),
  tar_target(mortality_daily, {
    get_mortality_daily(year = "2019", mortality = mortality)
  }),
  tar_target(AMI_limits, {
    Load_AMI_limits(year = "2019", hud_key = hud_key)
  }),
  tar_target(low_income_pop, {
    get_low_income_pop(year = "2019", AMI_limits = AMI_limits)
  }),
  tar_target(town_polygon_sf, {
    get_town_polygon_sf(CT_townships = CT_townships,
                        housing_town = housing_town,
                        low_income_pop = low_income_pop)
  }),
  tar_target(weights, {
    get_weights(CT_townships = CT_townships)
  }),
  tar_target(weights_post, {
    get_weights_post(weights = weights,
                     CT_townships = CT_townships,
                     housing_town = housing_town,
                     low_income_pop = low_income_pop,
                     AMI_limits = AMI_limits)
  }),
  tar_target(post_sim_1, {
    Merge_data(weights_post = weights_post, 
               mortality = mortality, 
               exposures_town = exposures_town)
  }),
  tar_target(post_heat_1, {
    get_post_heat_1(year = "2019",
                    weights_post = weights_post, 
                    exposures_town = exposures_town, 
                    mortality_daily = mortality_daily)
  }),
  # Simulation
  tar_target(all_sims, {
    simulate_cf(weights_post = weights_post, 
                post_sim_1 = post_sim_1, 
                post_heat_1 = post_heat_1,
                n = n)
  }),
  # Tables and Figures
  tar_target(sim_moved_totals_2, {
    get_sim_moved_totals_2(all_sims = all_sims, town_polygon_sf = town_polygon_sf)
  }),
  tar_target(figure2A, {
    get_figure2A(sim_moved_totals_2 = sim_moved_totals_2)
  }),
  tar_target(figure2A_data, {
    get_figure2A_data(sim_moved_totals_2 = sim_moved_totals_2)
  }),
  tar_target(figure2B, {
    get_figure2B(sim_moved_totals_2 = sim_moved_totals_2)
  }),
  tar_target(figure2B_data, {
    get_figure2B_data(sim_moved_totals_2 = sim_moved_totals_2)
  }),
  tar_target(table1, {
    get_table1(all_sims = all_sims)
  }),
  tar_target(avg_moved, {
    get_avg_moved(all_sims = all_sims, town_polygon_sf = town_polygon_sf)
  }),
  tar_target(table2, {
    get_table2(town_polygon_sf = town_polygon_sf, all_sims = all_sims)
  }),
  tar_target(tableA2, {
    get_tableA2(all_sims = all_sims, town_polygon_sf = town_polygon_sf)
  }),
  tar_target(figureA8, {
    get_figureA8(town_polygon_sf = town_polygon_sf, all_sims = all_sims, exposures_town = exposures_town)
  }),
  tar_target(tableA3, {
    get_tableA3(town_polygon_sf = town_polygon_sf, all_sims = all_sims, exposures_town = exposures_town)
  }),
  tar_target(figureA1, {
    get_figureA1(exposures_town = exposures_town)
  }),
  tar_target(figureA2, {
    get_figureA2(town_polygon_sf = town_polygon_sf)
  }),
  tar_target(figureA3_A6, {
    get_figureA3_A6(all_sims = all_sims, town_polygon_sf = town_polygon_sf)
  }),
  tar_target(figureA7, {
    get_figureA7(tableA2 = tableA2)
  }),
  tar_target(figureA9, {
    get_figureA9(CT_townships = CT_townships, AMI_limits = AMI_limits)
  }),
  tar_target(number_AMI_restricted, {
    get_number_AMI_restricted(low_income_pop = low_income_pop, AMI_limits = AMI_limits)
  }),
  tar_target(car_usage_pre_post, {
    get_car_usage_pre_post(sim_moved_totals_2 = sim_moved_totals_2)
  }),
  tar_target(tableA1, {
    get_tableA1(all_sims = all_sims, town_polygon_sf = town_polygon_sf)
  }),
  # Sensitivity Analysis 1
  tar_target(all_sims_s1, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    simulate_cf_s1(weights = weights, 
                   housing_town = housing_town, 
                   low_income_pop = low_income_pop, 
                   AMI_limits = AMI_limits, 
                   post_sim_1 = post_sim_1, 
                   post_heat_1 = post_heat_1, 
                   CT_townships = CT_townships, 
                   n = n)
  }),
  tar_target(tableA4, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_table1(all_sims = all_sims_s1)
  }),
  tar_target(avg_moved_s1, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_avg_moved(all_sims = all_sims_s1, town_polygon_sf = town_polygon_sf)
  }),
  # Sensitivity Analysis 2
  tar_target(weights_drive, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_weights_drive(CT_townships = CT_townships)
  }),
  tar_target(all_sims_drive, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    simulate_cf_drive(low_income_pop = low_income_pop,
                      housing_town = housing_town,
                      weights = weights_drive, 
                      AMI_limits = AMI_limits,
                      post_sim_1 = post_sim_1,
                      post_heat_1 = post_heat_1,
                      n = n)
  }),
  tar_target(tableA5, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_table1(all_sims = all_sims_drive)
  }),
  tar_target(avg_moved_drive, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_avg_moved(all_sims = all_sims_drive, town_polygon_sf = town_polygon_sf)
  }),
  # Sensitivity Analysis 3
  tar_target(exposures_town_18, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_exposures_town(year = "2018", CT_townships = CT_townships, downloads = downloads)
  }),
  tar_target(all_sims_18, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    simulate_cf_18(exposures_town_18 = exposures_town_18, CT_townships = CT_townships, downloads = downloads, n = n)
  }),
  tar_target(tableA6, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_table1(all_sims = all_sims_18)
  }),
  tar_target(avg_moved_18, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_avg_moved(all_sims = all_sims_18, town_polygon_sf = town_polygon_sf)
  }),
  # Save results
  tar_target(save_exhibits_all, {
    save_exhibits(figure2A = figure2A,
                  figure2A_data = figure2A_data,
                  figure2B = figure2B,
                  figure2B_data = figure2B_data,
                  table1 = table1,
                  table2 = table2,
                  figureA1 = figureA1,
                  figureA2 = figureA2,
                  figureA3_A6 = figureA3_A6,
                  tableA1 = tableA1,
                  tableA2 = tableA2,
                  figureA7 = figureA7,
                  figureA8 = figureA8,
                  tableA3 = tableA3,
                  tableA4 = tableA4,
                  tableA5 = tableA5,
                  tableA6 = tableA6,
                  figureA9 = figureA9,
                  avg_moved = avg_moved,
                  avg_moved_s1 = avg_moved_s1,
                  avg_moved_18 = avg_moved_18,
                  avg_moved_drive = avg_moved_drive,
                  car_usage_pre_post = car_usage_pre_post,
                  number_AMI_restricted = number_AMI_restricted)
  })
)

