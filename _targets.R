
library(targets)
library(tarchetypes)
library(future)
library(crew)

source("Functions.R")
options(tigris_use_cache = TRUE)

################ CUSTOMIZE THE FOLLOWING BEFORE RUNNING tar_make() #############
# required - define your Census API key here
census_key <- ""

# required - define your HUD API key here
hud_key <- ""

# optional - replace 10 with desired number of simulation repeats (1,000 for full simulation)
n <- 10

# optional - replace 1 with desired number of parallel workers (if greater than availableCores()-1, will default to availableCores()-1)
max_workers <- 1

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
    "ggpubr"
  ),
  format = 'rds', 
  workspace_on_error = TRUE,
  controller = crew_controller_local(workers = min(availableCores(constraints = "connections", omit = 1), max_workers))
)

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
    sim_moved_totals_2(all_sims = all_sims, town_polygon_sf = town_polygon_sf)
  }),
  tar_target(exhibit1, {
    exhibit1(sim_moved_totals_2 = sim_moved_totals_2)
  }),
  tar_target(exhibit1_data, {
    exhibit1_data(sim_moved_totals_2 = sim_moved_totals_2)
  }),
  tar_target(exhibit2, {
    exhibit2(sim_moved_totals_2 = sim_moved_totals_2)
  }),
  tar_target(exhibit2_data, {
    exhibit2_data(sim_moved_totals_2 = sim_moved_totals_2)
  }),
  tar_target(exhibit3, {
    exhibit3(all_sims = all_sims)
  }),
  tar_target(avg_moved, {
    get_avg_moved(all_sims = all_sims, town_polygon_sf = town_polygon_sf)
  }),
  tar_target(exhibit4, {
    exhibit4(town_polygon_sf = town_polygon_sf, all_sims = all_sims)
  }),
  tar_target(tableS1, {
    tableS1(all_sims = all_sims, town_polygon_sf = town_polygon_sf)
  }),
  tar_target(tableS2, {
    tableS2(town_polygon_sf = town_polygon_sf, all_sims = all_sims, exposures_town = exposures_town)
  }),
  tar_target(figureS1, {
    figureS1(exposures_town = exposures_town)
  }),
  tar_target(figureS2, {
    figureS2(town_polygon_sf = town_polygon_sf)
  }),
  tar_target(figureS3_S6, {
    figureS3_S6(all_sims = all_sims, town_polygon_sf = town_polygon_sf)
  }),
  tar_target(figureS7, {
    figureS7(tableS1 = tableS1)
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
  tar_target(exhibitA12, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    exhibit3(all_sims = all_sims_s1)
    }),
  tar_target(avg_moved_s1, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_avg_moved(all_sims = all_sims_s1, town_polygon_sf = town_polygon_sf)
    }),
  # Sensitivity Analysis 2
  tar_target(exposures_town_18, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_exposures_town(year = "2018", CT_townships = CT_townships, downloads = downloads)
    }),
  tar_target(all_sims_18, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    simulate_cf_18(exposures_town_18 = exposures_town_18, CT_townships = CT_townships, downloads = downloads, n = n)
    }),
  tar_target(exhibitA13, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    exhibit3(all_sims = all_sims_18)
    }),
  tar_target(avg_moved_18, {
    if (isTRUE(getOption("SKIP_SUPPL_ANALYSIS"))) return(NULL)
    get_avg_moved(all_sims = all_sims_18, town_polygon_sf = town_polygon_sf)
    }),
  # Save results
  tar_target(save_exhibits, {
    save_exhibits(exhibit1 = exhibit1,
                  exhibit1_data = exhibit1_data,
                  exhibit2 = exhibit2,
                  exhibit2_data = exhibit2_data,
                  exhibit3 = exhibit3,
                  exhibit4 = exhibit4,
                  figureS1 = figureS1,
                  figureS2 = figureS2,
                  figureS3_S6 = figureS3_S6,
                  tableS1 = tableS1,
                  figureS7 = figureS7,
                  tableS2 = tableS2,
                  exhibitA12 = exhibitA12,
                  exhibitA13 = exhibitA13,
                  avg_moved = avg_moved,
                  avg_moved_s1 = avg_moved_s1,
                  avg_moved_18 = avg_moved_18)
  })
)

