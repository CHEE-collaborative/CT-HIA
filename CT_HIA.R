
# CT HIA
# SP 5.31.23

# CT Townships & Census Tracts - EPA FAQSD pollution data

library(tidyverse)
library(sf)
library(here)
library(tidycensus)

# PM2.5
# read FAQSD and filter for July 4, 2019

pm <- read.delim(gzfile("2019_pm25_daily_average.txt.gz"), sep = ",") %>%
  filter(str_detect(`FIPS`, "^9")) %>% 
  filter(`Date` == "2019/07/04") %>% 
  mutate(GEOID = str_pad(FIPS, 11, side = "left", pad = "0"))

# get CT census tracts and join FAQSD to 2019 shapefile
# census_api_key()

options(tigris_use_cache = TRUE)

ct_tracts <- get_acs(
  state = "CT",
  geography = "tract",
  variables = "B02001_001",
  geometry = TRUE,
  year = 2019) %>% 
  inner_join(., pm, by = "GEOID")

# get CT town lines and polygon

town_lines <- st_read("CT Town Lines/Town_Lines.shp")

# plot air pollution by census tract

ggplot() + geom_sf(data = ct_tracts, aes(fill = pm25_daily_average.ug.m3.)) + 
  scale_fill_viridis_c(option = "magma", direction = -1) + geom_sf(data = town_lines)

# Ozone
# Temperature



# Affordable housing
# read in affordable housing by town

housing <- read.csv("Affordable_Housing_by_Town_2011-2022.csv") %>% 
  filter(`Year` == 2019) %>% 
  mutate("At least 10%?" = cut(`Percent.Affordable`, breaks = c(0,9.99999,100),
  labels = c("No", "Yes"))) %>% 
  rename("TOWN_NAME" = "Town")

# get CT town polygon

ct_units <- st_read("CT Vicinity Town Polygon/CT_Vicinity_Town_Polygon.shp") %>%
  filter(STATE_COD == "CT") %>% 
  inner_join(., housing, by = "TOWN_NAME")

# plot affordable housing units by town and 8-30g exemptions

ggplot() + geom_sf(data = ct_units, aes(fill = `Percent.Affordable`)) + 
  scale_fill_viridis_c(option = "magma", direction = -1)

ggplot() + geom_sf(data = ct_units, aes(fill = `At least 10%?`))

# Income
# plot income/SES

ggplot() + geom_sf(data = ct_tracts, aes(fill = estimate)) + 
  scale_fill_viridis_c(option = "magma", direction = -1) + geom_sf(data = town_lines)

# get low-income population by tract

# plot low-income population by tract

# Baseline incidence




