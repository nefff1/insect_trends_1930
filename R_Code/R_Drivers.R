################################################################################.
# SETUP ------------------------------------ ###################################
################################################################################.

# R version 4.3.3

# load packages ################################################################
################################################################################.

library(tidyverse); theme_set(theme_classic()) #2.0.0
library(sf) #1.0-15
library(data.table) #1.15.4
library(zoo) #1.8-12
library(ncdf4) #1.23
library(raster) #3.6-26
library(parallel)
library(lubridate) #1.9.3
library(sfheaders) #0.4.4

select <- dplyr::select

# set global parameters ########################################################
################################################################################.

v_zones <- c(Jura = "Jura",
             Plateau = "Plateau", 
             NorthernAlps = "Northern Alps",
             CentralAlps = "Central Alps",
             SouthernAlps = "Southern Alps",
             HighAlps = "High Alps")


v_textsize <- c(plotlabel = 9,
                axis.title = 8, axis.text = 7, 
                legend.title = 8, legend.text = 7,
                additional.text = 6) 

# read data ####################################################################
################################################################################.

# Census data ------------------------------------------------------------------.
# (see table for sources)
d_censuses <- read.csv2("Data/Drivers/Census_data_raw.csv")

# Temperature ------------------------------------------------------------------.
# reconstructed gridded temperature data from MeteoSwiss (www.meteoswiss.admin.ch)
nc_temperature <- nc_open("Data/Drivers/TrecabsM1901_ch02.lonlat_190101010000_202112010000.nc")

# biogeographic zones ----------------------------------------------------------.
# derived from elevation data and biogeographic regions of Switzerland
# both are available form https://data.geo.admin.ch (ch.bafu.biogeographische_regionen, ch.swisstopo.swissalti3d)
sf_zones <- st_read("Data/Other/biogeographic_zones.gpkg")

# size of the different zones
d_prop_zones <- sf_zones |> 
  as.data.frame() |> 
  filter(in_cscf) |> 
  group_by(zone) |> 
  summarise(n = n())

# cantons (spatial) ------------------------------------------------------------.
# Swiss boundaries layer 
# available from https://data.geo.admin.ch/ (ch.swisstopo.swissboundaries3d-kanton-flaeche.fill)
# version of 2022-05
sf_cantons <- st_read("Data/Other/swissBOUNDARIES3D_1_5_LV95_LN02.gpkg",
                      layer = "TLM_KANTONSGEBIET") |> 
  rename(NAME = "name")

# Swiss land-use statistics (Arealstatistik) -----------------------------------.
# available from https://data.geo.admin.ch/ (ch.bfs.arealstatistik)
d_area <- fread("Data/Other/ag-b-00.03-37-area-all-csv.csv")
d_area <- d_area %>% 
  rename(X = E_COORD,
         Y = N_COORD)

sf_area <- d_area |> #interpret as points
  mutate(X = X + 50, # set point in the middle of the quadrat (this is imporant later when joining the zone squares!!)
         Y = Y + 50) |> 
  st_as_sf(coords = c("X", "Y"), crs = 2056)

################################################################################.
# CALCULATE DRIVERS  -------------- ############################################
################################################################################.

length_interval <- 8 # years
length_interval2 <- 12 # years

# define consecutive time intervals of certain length (for insect data). 
# Double-years, always last observation overlapping (to include the gap as well)
year_insect <- seq(1930, 2020, 2)

year_insect_start <- seq(min(year_insect), max(year_insect),
                         length_interval - 2)
year_insect_start <-
  year_insect_start[year_insect_start <= max(year_insect) + 1 - (length_interval - 1)]

year_insect_start2 <- seq(min(year_insect), max(year_insect),
                          length_interval2 - 2)
year_insect_start2 <-
  year_insect_start2[year_insect_start2 <= max(year_insect) + 1 - (length_interval2 - 1)]


# Human population change ######################################################
################################################################################.

d_pop <- d_censuses |> 
  filter(parameter == "human population")

# deal with Jura ---------------------------------------------------------------.

d_prop_Jura <-
  d_pop |> 
  group_by(year) |> 
  reframe(prop_Jura = value[canton == "Jura"]/
            (value[canton == "Bern"] + value[canton == "Jura"])) |> 
  filter(year == 1971)


d_pop_Jura <- d_pop |> 
  filter(canton == "Bern", year < 1971) |> 
  mutate(value = value * d_prop_Jura$prop_Jura,
         canton = "Jura")

d_pop <- d_pop |> 
  bind_rows(d_pop_Jura) |> 
  group_by(year) |> 
  mutate(value = ifelse(canton == "Bern" & year < 1971,
                        value - value[canton == "Jura"],
                        value)) |> 
  ungroup() 

# distribute to biogeographic zones --------------------------------------------.

sf_comb <-
  sf_area |> 
  filter(LU85_10 == 100) |> # where buildings are 
  st_join(sf_cantons |> 
            select(canton = NAME)) |> 
  st_join(sf_zones) |> 
  filter(!is.na(canton)) # some squares at edge of Switzerland

d_comb_smry <- sf_comb |> 
  as.data.frame() |> 
  group_by(canton, zone, in_cscf) |> 
  summarise(n_points = n(),
            .groups = "drop") |> 
  group_by(canton) |> 
  mutate(perc = n_points / sum(n_points)) |> 
  ungroup()

d_pop_zone <-
  d_pop |> 
  left_join(d_comb_smry, by = "canton", relationship = "many-to-many") |> 
  filter(in_cscf) |>
  group_by(zone, year) |> 
  summarise(value = sum(value * perc),
            .groups = "drop")

# d_pop_zone |> 
#   ggplot(aes(x = year, y = value)) +
#   geom_point() +
#   stat_smooth() +
#   facet_grid(~ zone)

# zonal summary ----------------------------------------------------------------.

d_pop_zone_inter <- expand_grid(year = seq(min(d_pop_zone$year), max(d_pop_zone$year), 1),
                                zone = unique(d_pop_zone$zone)) |>   
  full_join(d_pop_zone, by = c("year", "zone")) %>%
  group_by(zone) |> 
  mutate(value = na.approx(value)) |> 
  ungroup()


d_LM <-
  d_pop_zone_inter |> 
  select(-c(value)) |> 
  left_join(d_pop_zone_inter |> 
              left_join(d_prop_zones, by = "zone") |> 
              mutate(value_rel = value / (n * 25)) |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  filter(year %in% year_insect_start, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  summarise(lm_coef = lm(value_rel ~ year_end)$coefficients[2],
            lm_intercept = lm(value_rel ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval - 1)/2,
         pred_start = lm_intercept + year_start * lm_coef,
         pred_end = lm_intercept + (year_start + length_interval - 1) * lm_coef) 

d_LM2 <-
  d_pop_zone_inter |> 
  select(-c(value)) |> 
  left_join(d_pop_zone_inter |> 
              left_join(d_prop_zones, by = "zone") |> 
              mutate(value_rel = value / (n * 25)) |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  filter(year %in% year_insect_start2, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval2 & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  summarise(lm_coef = lm(value_rel ~ year_end)$coefficients[2],
            lm_intercept = lm(value_rel ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval2 - 1)/2,
         pred_start = lm_intercept + year_start * lm_coef,
         pred_end = lm_intercept + (year_start + length_interval2 - 1) * lm_coef) 


d_drivers <- d_LM |> 
  select(zone, year_mean, human_pop_change = lm_coef) |> 
  mutate(length_interval = length_interval) |> 
  bind_rows(d_LM2 |> 
              select(zone, year_mean, human_pop_change = lm_coef) |> 
              mutate(length_interval = length_interval2))

# plotting ---------------------------------------------------------------------.

d_pop_zone_inter |>
  left_join(d_pop_zone |> 
              select(-value) |> 
              mutate(Interpolated = F),
            by = c("zone", "year")) |> 
  left_join(d_prop_zones, by = "zone") |> 
  mutate(value_rel = value / (n * 25),
         Interpolated = ifelse(is.na(Interpolated), T, Interpolated),
         zone = factor(v_zones[zone], level = v_zones)) |> 
  filter(year >= 1930) |> 
  ggplot(aes(x = year, y = value_rel)) +
  geom_tile(data = d_drivers |> 
              filter(length_interval == 8) |> 
              mutate(human_pop_change = human_pop_change / 
                       (2 * sd(human_pop_change)),
                     zone = factor(v_zones[zone], level = v_zones)), 
            aes(x = year_mean, y = 1, fill = human_pop_change),
            height = Inf, width = length_interval - 2,
            colour = "grey20", alpha= .8) +
  geom_point(aes(size = Interpolated)) +
  geom_segment(data = d_LM |> 
                 mutate(zone = factor(v_zones[zone], level = v_zones),
                        pred_start = pred_start,
                        pred_end = pred_end), 
               aes(x = year_start, xend = year_start + length_interval - 1,
                   y = pred_start, yend = pred_end), 
               linewidth = .75, alpha = .7, colour = "sienna1") +
  facet_wrap(~ zone) +
  scale_fill_gradient2(low = "#003c30", high = "#543005",
                       name = "Std.\nhuman population\ndensity change") +
  labs(y = "Human population density\n(individuals / km2)", x = "Year") +
  scale_x_continuous(breaks = c(1950, 1975, 2000)) +
  scale_size_manual(values = c(.75, .25)) +
  theme(axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        strip.text = element_text(size = v_textsize["axis.title"]),
        legend.title = element_text(size = v_textsize["legend.title"]),
        legend.text = element_text(size = v_textsize["legend.text"]))

ggsave("Output/Figures/Drivers_human_pop.jpeg", width = 10, height = 5.5)

# Forest: area #################################################################
################################################################################.

d_forestarea <- d_censuses |> 
  filter(parameter == "forest: area")

# deal with Jura ---------------------------------------------------------------.

d_prop_Jura <- d_forestarea |> 
  group_by(forest_owner, year) |> 
  reframe(prop_Jura = value[canton == "Jura"]/
            (value[canton == "Bern"] + value[canton == "Jura"])) |> 
  filter(year == 1979)

d_forestarea_Jura <- d_forestarea |> 
  filter(canton == "Bern", year < 1979) |> 
  left_join(d_prop_Jura |> select(-year),
            by = "forest_owner") |> 
  mutate(value = value * prop_Jura,
         canton = "Jura") |> 
  select(-prop_Jura)


d_forestarea <- d_forestarea |> 
  bind_rows(d_forestarea_Jura) |> 
  group_by(year, forest_owner) |> 
  mutate(value = ifelse(canton == "Bern" & year < 1979,
                        value - value[canton == "Jura"],
                        value)) |> 
  ungroup() 

# distribute to biogeographic zones --------------------------------------------.

sf_comb <-
  sf_area |> 
  filter(LU85_10 == 300) |> # where forest is
  st_join(sf_cantons |> 
            select(canton = NAME)) |> 
  st_join(sf_zones) |> 
  filter(!is.na(canton)) # some squares at edge of Switzerland

d_comb_smry <- sf_comb |> 
  as.data.frame() |> 
  group_by(canton, zone, in_cscf) |> 
  summarise(n_points = n(),
            .groups = "drop") |> 
  group_by(canton) |> 
  mutate(perc = n_points / sum(n_points)) |> 
  ungroup()



d_forestarea_zone <- d_forestarea |> 
  left_join(d_comb_smry, by = "canton", relationship = "many-to-many") |> 
  filter(!is.na(perc)) |> # SBB, etc
  filter(in_cscf) |>
  group_by(zone, year, forest_owner) |> 
  summarise(value = sum(value * perc),
            .groups = "drop")

# deal with year without private data 
d_odds_priv <-
  d_forestarea_zone |> 
  group_by(zone, year) |> 
  reframe(odds_private = value[forest_owner == "private"] / 
            value[forest_owner == "public"]) |> 
  filter(year %in% c(1950, 1955)) |> #first years with all data
  group_by(zone) |> 
  summarise(odds_private = mean(odds_private),
            .groups = "drop")

d_addon_private <- d_forestarea_zone |> 
  filter(year < 1950) |> 
  left_join(d_odds_priv, by = "zone") |> 
  mutate(value = value * odds_private,
         forest_owner = "private") |> 
  select(-odds_private)

d_forestarea_zone <- d_forestarea_zone |> 
  bind_rows(d_addon_private)


# summarise everything up
d_forestarea_zone <- d_forestarea_zone |> 
  group_by(zone, year) |> 
  summarise(value = sum(value),
            .groups = "drop")


# reconstruct for missing years 1951-1954, 1956, 1957, 1959

d_addon1 <- d_forestarea_zone |> 
  filter(year %in% c(1949, 1955)) |> # 1950 is somewhat of an outlier, thus take 1949 as starting year
  group_by(zone) |> 
  group_modify(~data.frame(year = 1951:1954,
                           value = predict(lm(value ~ year, data = .),
                                           data.frame(year = 1951:1954))))

d_addon2 <- d_forestarea_zone |> 
  filter(year %in% c(1955, 1958)) |> # 1950 is somewhat of an outlier, thus take 1949 as starting year
  group_by(zone) |> 
  group_modify(~data.frame(year = 1956:1957,
                           value = predict(lm(value ~ year, data = .),
                                           data.frame(year = 1956:1957))))

d_addon3 <- d_forestarea_zone |> 
  filter(year %in% c(1958, 1960)) |> # 1950 is somewhat of an outlier, thus take 1949 as starting year
  group_by(zone) |> 
  group_modify(~data.frame(year = 1959,
                           value = predict(lm(value ~ year, data = .),
                                           data.frame(year = 1959))))


d_forestarea_zone <- d_forestarea_zone |> 
  bind_rows(d_addon1) |>
  bind_rows(d_addon2) |>
  bind_rows(d_addon3)

# zonal summary ----------------------------------------------------------------.

d_LM <-
  d_forestarea_zone |> 
  select(-value) |> 
  left_join(d_forestarea_zone |> 
              left_join(d_prop_zones, by = "zone") |> 
              mutate(value_rel = value / (n * 25)) |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  filter(year %in% year_insect_start, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  mutate(year_end = year_end - year_start - (length_interval - 1) / 2) |> # center year such that intercept refers to estimated mean
  summarise(lm_coef = lm(value_rel ~ year_end)$coefficients[2],
            lm_intercept = lm(value_rel ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval - 1)/2,
         pred_start = lm_intercept - (length_interval - 1)/2 * lm_coef,
         pred_end = lm_intercept + (length_interval - 1)/2 * lm_coef) 

d_LM2 <-
  d_forestarea_zone |> 
  select(-value) |> 
  left_join(d_forestarea_zone |> 
              left_join(d_prop_zones, by = "zone") |> 
              mutate(value_rel = value / (n * 25)) |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  filter(year %in% year_insect_start2, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval2 & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  mutate(year_end = year_end - year_start - (length_interval2 - 1) / 2) |> # center year such that intercept refers to estimated mean
  summarise(lm_coef = lm(value_rel ~ year_end)$coefficients[2],
            lm_intercept = lm(value_rel ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval2 - 1)/2,
         pred_start = lm_intercept - (length_interval2 - 1)/2 * lm_coef,
         pred_end = lm_intercept + (length_interval2 - 1)/2 * lm_coef) 

d_drivers <-
  d_drivers |> 
  full_join(d_LM |> 
              mutate(length_interval = length_interval) |> 
              select(zone, year_mean, length_interval, 
                     forest_area_change = lm_coef) |> 
              bind_rows(d_LM2 |> 
                          mutate(length_interval = length_interval2) |> 
                          select(zone, year_mean, length_interval, 
                                 forest_area_change = lm_coef)) |> 
              mutate(forest_area_change = sign(forest_area_change) * 
                       (abs(forest_area_change) ^ (1/3))),
            by = c("zone", "year_mean", "length_interval"), relationship = "one-to-one") |>
  full_join(d_LM |> 
              select(zone, year_mean, forest_area_abs = lm_intercept) |> 
              mutate(length_interval = length_interval) |> 
              bind_rows(d_LM2 |> 
                          select(zone, year_mean, forest_area_abs = lm_intercept) |> 
                          mutate(length_interval = length_interval2)),
            by = c("zone", "year_mean", "length_interval"), relationship = "one-to-one")

# plotting ---------------------------------------------------------------------.

d_forestarea_zone |>
  left_join(d_prop_zones, by = "zone") |> 
  mutate(value_rel = value / (25 * n),
         zone = factor(v_zones[zone], level = v_zones),
         Interpolated = year %in% c(1951:1954, 1956:1957, 1959)) |> 
  ggplot(aes(x = year, y = value_rel)) +
  geom_tile(data = d_drivers |> 
              filter(length_interval == 8) |> 
              mutate(forest_area_change = forest_area_change / 
                       (2 * sd(forest_area_change)),
                     zone = factor(v_zones[zone], level = v_zones)), 
            aes(x = year_mean, y = 1, fill = forest_area_change),
            height = Inf, width = length_interval - 2,
            colour = "grey20", alpha = .8) +
  geom_point(aes(size = Interpolated)) +
  geom_segment(data = d_LM |> 
                 mutate(zone = factor(v_zones[zone], level = v_zones),
                        pred_start = pred_start,
                        pred_end = pred_end), 
               aes(x = year_start, xend = year_start + length_interval - 1,
                   y = pred_start, yend = pred_end), 
               linewidth = .75, alpha = .7, colour = "sienna1") +
  facet_wrap(~ zone) +
  scale_fill_gradient2(low = "#003c30", high = "#543005",
                       name = "Std.\nforest area\nchange") +
  scale_size_manual(values = c(.5, 1.5)) +
  labs(y = "Forest area (%)", x = "Year") +
  scale_x_continuous(breaks = c(1950, 1975, 2000)) +
  scale_size_manual(values = c(.75, .25)) +
  theme(axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        strip.text = element_text(size = v_textsize["axis.title"]),
        legend.title = element_text(size = v_textsize["legend.title"]),
        legend.text = element_text(size = v_textsize["legend.text"]))


ggsave("Output/Figures/Drivers_forest_area.jpeg", width = 10, height = 5.5)

# Forest: wood harvest #########################################################
################################################################################.


d_woodharvest <- d_censuses |> 
  filter(parameter == "forest: wood harvest")

# deal with Jura ---------------------------------------------------------------.

d_prop_Jura <-
  d_woodharvest |> 
  group_by(forest_owner, year) |> 
  reframe(prop_Jura = value[canton == "Jura"]/
            (value[canton == "Bern"] + value[canton == "Jura"])) |> 
  filter(year %in% 1979:1984) |> # more volatile than area, so take mean of several years
  group_by(forest_owner) |> 
  summarise(prop_Jura = mean(prop_Jura),
            .groups = "drop")

d_woodharvest_Jura <- d_woodharvest |> 
  filter(canton == "Bern", year < 1979) |> 
  left_join(d_prop_Jura,
            by = "forest_owner") |> 
  mutate(value = value * prop_Jura,
         canton = "Jura") |> 
  select(-prop_Jura)

d_woodharvest <- d_woodharvest |> 
  bind_rows(d_woodharvest_Jura) |> 
  group_by(year, forest_owner) |> 
  mutate(value = ifelse(canton == "Bern" & year < 1979,
                        value - value[canton == "Jura"],
                        value)) |> 
  ungroup() 

# distribute to biogeographic zones --------------------------------------------.

sf_comb <-
  sf_area |> 
  filter(LU85_10 == 300) |>  # where forest is
  st_join(sf_cantons |> 
            select(canton = NAME)) |> 
  st_join(sf_zones) |> 
  filter(!is.na(canton)) # some squares at edge of Switzerland

d_comb_smry <- sf_comb |> 
  as.data.frame() |> 
  group_by(canton, zone, in_cscf) |> 
  summarise(n_points = n(),
            .groups = "drop") |> 
  group_by(canton) |> 
  mutate(perc = n_points / sum(n_points)) |> 
  ungroup()

# now distribute (holz) --------------------------------------------------------.

d_woodharvest_zone <-
  d_woodharvest |> 
  left_join(d_comb_smry, by = "canton", relationship = "many-to-many") |> 
  filter(!is.na(perc)) |> # SBB, etc
  filter(in_cscf) |>
  group_by(zone, year, forest_owner) |> 
  summarise(value = sum(value * perc))


# deal with year without private data 

d_woodharvest_odds_priv <-
  d_woodharvest_zone |> 
  group_by(zone, year) |> 
  reframe(odds_private = value[forest_owner == "private"] / 
            value[forest_owner == "public"]) |>  
  filter(year %in% 1950:1955) |> #first years with all data (more volatile than area)
  group_by(zone) |> 
  summarise(odds_private = mean(odds_private),
            .groups = "drop")

d_addon_private <- d_woodharvest_zone |> 
  filter(year < 1950) |> 
  left_join(d_woodharvest_odds_priv, by = "zone") |> 
  mutate(value = value * odds_private,
         forest_owner = "private") |> 
  select(-odds_private)

d_woodharvest_zone <- d_woodharvest_zone |> 
  bind_rows(d_addon_private)

# summarise everything up
d_woodharvest_zone <- d_woodharvest_zone |> 
  group_by(zone, year) |> 
  summarise(value = sum(value),
            .groups = "drop")

# zonal summary ----------------------------------------------------------------.

d_intensity_zone <- d_forestarea_zone |> 
  inner_join(d_woodharvest_zone, by = c("zone", "year"),
             suffix = c(".forestarea", ".woodharvest")) |> 
  mutate(intensity = value.woodharvest / value.forestarea)

d_intensity_zone |> 
  ggplot(aes(x = year, y = intensity)) +
  geom_point() +
  facet_wrap(~ zone)


# major storm events as evident from the plots (sometimes two outlier years apparent)
# outliers related to WW II not included
d_storm_outliers <-  
  data.frame(zone = "NorthernAlps", year = 1990) |> # Viviane
  add_row(zone = "NorthernAlps", year = 2000) |> # Lothar
  add_row(zone = "NorthernAlps", year = 2001) |> # Lothar
  add_row(zone = "HighAlps", year = 1990) |> # Viviane
  add_row(zone = "HighAlps", year = 1991) |> # Viviane
  add_row(zone = "HighAlps", year = 2000) |> # Lothar
  add_row(zone = "Jura", year = 2000) |> # Lothar
  add_row(zone = "Plateau", year = 1962) |> # unknown
  add_row(zone = "Plateau", year = 1967) |> # unknown
  add_row(zone = "Plateau", year = 1990) |> # Viviane
  add_row(zone = "Plateau", year = 2000) |> # Lothar
  add_row(zone = "Plateau", year = 2001) |> # Lothar
  add_row(zone = "CentralAlps", year = 1990) |>  # Viviane
  add_row(zone = "CentralAlps", year = 1991) # Viviane

d_intensity_zone_robust <- d_intensity_zone |> 
  left_join(d_storm_outliers |> 
              mutate(storm = T),
            by = c("zone", "year")) |>
  filter(is.na(storm)) |> 
  select(-storm)

d_LM <-
  d_intensity_zone_robust |> 
  select(-c(value.forestarea, value.woodharvest, intensity)) |> 
  left_join(d_intensity_zone_robust |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  # trick to include missing data years
  mutate(year = ifelse(year == 1989 & zone %in% c("NorthernAlps", "HighAlps", "Plateau", "CentralAlps"), 
                       1990, year)) |> 
  filter(year %in% year_insect_start, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  summarise(lm_coef = lm(intensity ~ year_end)$coefficients[2],
            lm_intercept = lm(intensity ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval - 1)/2,
         pred_start = lm_intercept + year_start * lm_coef,
         pred_end = lm_intercept + (year_start + length_interval - 1) * lm_coef) 

d_LM2 <-
  d_intensity_zone_robust |> 
  select(-c(value.forestarea, value.woodharvest, intensity)) |> 
  left_join(d_intensity_zone_robust |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  # trick to include missing data years
  mutate(year = ifelse(year == 1989 & zone %in% c("NorthernAlps", "HighAlps", "Plateau", "CentralAlps"), 
                       1990, year),
         year = ifelse(year == 1999 & zone %in% c("NorthernAlps", "HighAlps", "Jura", "Plateau"), 
                       2000, year)) |> 
  filter(year %in% year_insect_start2, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval2 & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  summarise(lm_coef = lm(intensity ~ year_end)$coefficients[2],
            lm_intercept = lm(intensity ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval2 - 1)/2,
         pred_start = lm_intercept + year_start * lm_coef,
         pred_end = lm_intercept + (year_start + length_interval2 - 1) * lm_coef) 


d_stormperiods <- d_LM |> 
  left_join(d_storm_outliers, by = "zone",
            relationship = "many-to-many") |> 
  mutate(year_diff = year_start - year,
         stormgroups = case_when(year %in% c(1990, 1991) ~ 1990.5,
                                 year %in% c(2000, 2001) ~ 2000.5,
                                 .default = year)) |> 
  group_by(zone, stormgroups) |> 
  filter(year_diff >= 0) |> 
  filter(year_diff == min(year_diff)) |> 
  ungroup() |> 
  select(zone, year_mean) |> 
  mutate(stormperiod = T) |> 
  distinct()

d_stormperiods2 <- d_LM2 |> 
  left_join(d_storm_outliers, by = "zone",
            relationship = "many-to-many") |> 
  mutate(year_diff = year_start - year,
         stormgroups = case_when(year %in% c(1990, 1991) ~ 1990.5,
                                 year %in% c(2000, 2001) ~ 2000.5,
                                 .default = year)) |> 
  group_by(zone, stormgroups) |> 
  filter(year_diff >= 0) |> 
  filter(year_diff == min(year_diff)) |> 
  ungroup() |> 
  select(zone, year_mean) |> 
  mutate(stormperiod = T) |> 
  distinct()


d_drivers <- d_drivers |> 
  full_join(d_LM |> 
              select(zone, year_mean, wood_harvest_int_change = lm_coef) |> 
              mutate(length_interval = length_interval) |> 
              bind_rows(d_LM2 |> 
                          select(zone, year_mean, wood_harvest_int_change = lm_coef) |> 
                          mutate(length_interval = length_interval2)),
            by = c("zone", "year_mean", "length_interval"), relationship = "one-to-one") |> 
  full_join(d_LM |> 
              left_join(d_stormperiods, by = c("zone", "year_mean")) |> 
              mutate(after_storm = !is.na(stormperiod)) |> 
              select(zone, year_mean, after_storm) |> 
              mutate(length_interval = length_interval) |> 
              bind_rows(d_LM2 |> 
                          left_join(d_stormperiods2, by = c("zone", "year_mean")) |> 
                          mutate(after_storm = !is.na(stormperiod)) |> 
                          select(zone, year_mean, after_storm) |> 
                          mutate(length_interval = length_interval2)),
            by = c("zone", "year_mean", "length_interval"), relationship = "one-to-one")

# plotting ---------------------------------------------------------------------.

d_intensity_zone |>
  left_join(d_storm_outliers |> mutate(Storm = T),
            by = join_by(zone, year)) |> 
  mutate(Storm = ifelse(is.na(Storm), F, Storm),
         zone = factor(v_zones[zone], level = v_zones)) |> 
  ggplot(aes(x = year, y = intensity)) +
  geom_tile(data = d_drivers |> 
              filter(length_interval == 8) |> 
              mutate(wood_harvest_int_change = (wood_harvest_int_change - mean(wood_harvest_int_change)) / 
                       (2 * sd(wood_harvest_int_change)),
                     zone = factor(v_zones[zone], level = v_zones)), 
            aes(x = year_mean, y = 1, fill = wood_harvest_int_change),
            height = Inf, width = length_interval - 2,
            colour = "grey20", alpha = .8) +
  geom_point(aes(size = Storm)) +
  geom_segment(data = d_LM |> 
                 mutate(zone = factor(v_zones[zone], level = v_zones)), 
               aes(x = year_start, xend = year_start + length_interval - 1,
                   y = pred_start, yend = pred_end), 
               linewidth = .75, alpha = .7, colour = "sienna1") +
  geom_tile(data = d_stormperiods |> 
              mutate(zone = factor(v_zones[zone], level = v_zones)),
            aes(x = year_mean, y = min(d_intensity_zone$intensity) + 1 / 2),
            fill = "red2", 
            height = 1,
            width = length_interval - 2,
            colour = NA) +
  facet_wrap(~ zone) +
  scale_fill_gradient2(low = "#003c30", high = "#543005",
                       name = "Std.\nwood harvest\nintensity change") +
  scale_size_manual(values = c(.5, 1.5)) +
  labs(y = "Wood harvest intensity\n(m3 / hectare)", x = "Year") +
  scale_x_continuous(breaks = c(1950, 1975, 2000)) +
  theme(axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        strip.text = element_text(size = v_textsize["axis.title"]),
        legend.title = element_text(size = v_textsize["legend.title"]),
        legend.text = element_text(size = v_textsize["legend.text"]))

ggsave("Output/Figures/Drivers_woodharvest.jpeg", width = 10, height = 5.5)

# Tractors (mechanisation) #####################################################
################################################################################.

d_tractors <- d_censuses |> 
  filter(parameter == "tractors")

# deal with Jura ---------------------------------------------------------------.

d_prop_Jura <-
  d_tractors |> 
  group_by(year) |> 
  reframe(prop_Jura = value[canton == "Jura"]/
            (value[canton == "Bern"] + value[canton == "Jura"])) |> 
  filter(year == min(year))


d_tractors_Jura <- d_tractors |> 
  filter(canton == "Bern", year < 1980) |> 
  mutate(value = value * d_prop_Jura$prop_Jura, # assuming constant proportion since the very beginning...
         canton = "Jura")


d_tractors <- d_tractors |> 
  bind_rows(d_tractors_Jura) |> 
  group_by(year) |> 
  mutate(value = ifelse(canton == "Bern" & year < 1980,
                        value - value[canton == "Jura"],
                        value)) |> 
  ungroup() 

# distribute to biogeographic zones --------------------------------------------.

sf_comb <-
  sf_area |> 
  filter(LU85_10 == 220) |> # agricultural areas in which tractors are used
  st_join(sf_cantons |> 
            select(canton = NAME)) |> 
  st_join(sf_zones) |> 
  filter(!is.na(canton)) # some squares at edge of Switzerland

d_comb_smry <- sf_comb |> 
  as.data.frame() |> 
  group_by(canton, zone, in_cscf) |> 
  summarise(n_points = n(),
            .groups = "drop") |> 
  group_by(canton) |> 
  mutate(perc = n_points / sum(n_points)) |> 
  ungroup()

d_tractors_zone <-
  d_tractors |> 
  left_join(d_comb_smry, by = "canton", relationship = "many-to-many") |> 
  filter(in_cscf) |>
  group_by(zone, year) |> 
  summarise(value = sum(value * perc),
            .groups = "drop")

# zonal summary ----------------------------------------------------------------.

d_tractors_zone_inter <- expand_grid(year = seq(min(d_tractors_zone$year), max(d_tractors_zone$year), 1),
                                     zone = unique(d_tractors_zone$zone)) |>   
  full_join(d_tractors_zone, by = c("year", "zone")) %>%
  group_by(zone) |> 
  mutate(value = na.approx(value)) |> 
  ungroup()

d_LM <-
  d_tractors_zone_inter |> 
  select(-c(value)) |> 
  left_join(d_tractors_zone_inter |> 
              left_join(d_prop_zones, by = "zone") |> 
              mutate(value_rel = value / (n * 25)) |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  filter(year %in% year_insect_start, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  mutate(year_end = year_end - year_start - (length_interval - 1) / 2) |> # center year such that intercept refers to estimated mean
  summarise(lm_coef = lm(value_rel ~ year_end)$coefficients[2],
            lm_intercept = lm(value_rel ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval - 1)/2,
         pred_start = lm_intercept - (length_interval - 1) / 2 * lm_coef,
         pred_end = lm_intercept + (length_interval - 1) / 2 * lm_coef) 

d_LM2 <-
  d_tractors_zone_inter |> 
  select(-c(value)) |> 
  left_join(d_tractors_zone_inter |> 
              left_join(d_prop_zones, by = "zone") |> 
              mutate(value_rel = value / (n * 25)) |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  filter(year %in% year_insect_start2, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval2 & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  mutate(year_end = year_end - year_start - (length_interval2 - 1) / 2) |> # center year such that intercept refers to estimated mean
  summarise(lm_coef = lm(value_rel ~ year_end)$coefficients[2],
            lm_intercept = lm(value_rel ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval2 - 1)/2,
         pred_start = lm_intercept - (length_interval - 1) / 2 * lm_coef,
         pred_end = lm_intercept + (length_interval - 1) / 2 * lm_coef) 

d_drivers <-
  d_drivers |> 
  full_join(d_LM |> 
              select(zone, year_mean, mechanisation_abs = lm_intercept) |> 
              mutate(length_interval = length_interval) |> 
              bind_rows(d_LM2 |> 
                          select(zone, year_mean, mechanisation_abs = lm_intercept) |> 
                          mutate(length_interval = length_interval2)),
            by = c("zone", "year_mean", "length_interval"), relationship = "one-to-one") |> 
  full_join(d_LM |> 
              select(zone, year_mean) |> 
              mutate(mechanisation_period = year_mean >= 1950 & year_mean <= 1990 &
                       zone %in% c("Plateau", "NorthernAlps", "Jura"),
                     length_interval = length_interval) |> 
              bind_rows(d_LM2 |> 
                          select(zone, year_mean) |> 
                          mutate(mechanisation_period = year_mean >= 1950 & year_mean <= 1990 &
                                   zone %in% c("Plateau", "NorthernAlps", "Jura"),
                                 length_interval = length_interval2)),
            by = c("zone", "year_mean", "length_interval"), relationship = "one-to-one")


# plotting ---------------------------------------------------------------------.

d_tractors_zone_inter |>
  left_join(d_tractors_zone |> 
              select(-value) |> 
              mutate(Interpolated = F),
            by = c("zone", "year")) |> 
  left_join(d_prop_zones, by = "zone") |> 
  mutate(value_rel = value / (n * 25),
         Interpolated = ifelse(is.na(Interpolated), T, Interpolated),
         zone = factor(v_zones[zone], level = v_zones)) |> 
  filter(year >= 1930) |> 
  ggplot(aes(x = year, y = value_rel)) +
  geom_tile(data = d_drivers |>
              filter(length_interval == 8) |> 
              mutate(zone = factor(v_zones[zone], level = v_zones)),
            aes(x = year_mean, y = 1, fill = mechanisation_period),
            height = Inf, width = length_interval - 2,
            colour = "grey20", alpha= .8) +
  geom_point(aes(size = Interpolated)) +
  geom_segment(data = d_LM |> 
                 mutate(zone = factor(v_zones[zone], level = v_zones),
                        pred_start = pred_start,
                        pred_end = pred_end), 
               aes(x = year_start, xend = year_start + length_interval - 1,
                   y = pred_start, yend = pred_end), 
               linewidth = .75, alpha = .7, colour = "sienna1") +
  facet_wrap(~ zone) +
  scale_fill_manual(values = c("white", "turquoise1"),
                    name = "Mechanisation period") +
  labs(y = "Mechanisation\n(tractors / km2)", x = "Year") +
  scale_x_continuous(breaks = c(1950, 1975, 2000)) +
  scale_size_manual(values = c(.75, .25)) +
  theme(axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        strip.text = element_text(size = v_textsize["axis.title"]),
        legend.title = element_text(size = v_textsize["legend.title"]),
        legend.text = element_text(size = v_textsize["legend.text"]))

ggsave("Output/Figures/Drivers_mechanisation.jpeg", width = 10, height = 5.5)

# Grassland area ###############################################################
################################################################################.

d_GL <- d_censuses |> 
  filter(parameter == "grassland area")

# deal with Jura ---------------------------------------------------------------.
d_prop_Jura <- d_GL |> 
  filter(year == 1975,
         canton %in% c("Bern", "Jura")) |> 
  reframe(prop_Jura = value[canton == "Jura"]/
            (value[canton == "Bern"] + value[canton == "Jura"]))


d_GL_Jura <- d_GL |> 
  filter(canton == "Bern", year < 1975) |> 
  mutate(value = value * d_prop_Jura$prop_Jura, # assuming constant proportion since the very beginning...
         canton = "Jura")


d_GL <- d_GL |> 
  bind_rows(d_GL_Jura) |> 
  group_by(year) |> 
  mutate(value = ifelse(canton == "Bern" & year < 1975,
                        value - value[canton == "Jura"],
                        value)) |> 
  ungroup() 

# distribute to biogeographic zones --------------------------------------------.

sf_comb <-
  sf_area |> 
  filter(LU85_10 == 220) |> # categories containing grasslands (temporary and permanent)
  st_join(sf_cantons |> 
            select(canton = NAME)) |> 
  st_join(sf_zones) |> 
  filter(!is.na(canton)) # some squares at edge of Switzerland


d_comb_smry <- sf_comb |> 
  as.data.frame() |> 
  group_by(canton, zone, in_cscf) |> 
  summarise(n_points = n(),
            .groups = "drop") |> 
  group_by(canton) |> 
  mutate(perc = n_points / sum(n_points)) |> 
  ungroup()

d_GL_zone <-
  d_GL |> 
  left_join(d_comb_smry, by = "canton", relationship = "many-to-many") |> 
  filter(in_cscf) |>
  group_by(zone, year) |> 
  summarise(value = sum(value * perc),
            .groups = "drop")

# interpolate missing values
d_GL_zone_inter <- expand_grid(year = seq(min(d_GL_zone$year), max(d_GL_zone$year), 1),
                               zone = unique(d_GL_zone$zone)) |>   
  full_join(d_GL_zone, by = c("year", "zone")) %>%
  group_by(zone) |> 
  mutate(value = na.approx(value)) |> 
  ungroup()


# add summering pastures -------------------------------------------------------.

d_summering_lt <-  sf_area |>
  filter(LU85_10 %in% c(240) |
           LU97_10 %in% c(240) |
           LU09_10 %in% c(240) |
           LU18_10 %in% c(240)) |> 
  st_join(sf_cantons |> 
            select(canton = NAME)) |> 
  st_join(sf_zones) |> 
  filter(!is.na(canton)) |> # some squares at edge of Switzerland 
  as.data.frame() |> 
  select(zone, in_cscf, RELI, 
         FJ_85 = FJ85, FJ_97 = FJ97, FJ_09 = FJ09, FJ_18 = FJ18, 
         LU_85 = LU85_10, LU_97 = LU97_10, LU_09 = LU09_10, LU_18 = LU18_10) |> 
  pivot_longer(-c(zone, in_cscf, RELI),
               names_pattern = "(.*)_(.*)",
               names_to = c(".value", "repetition"))

d_summering_smry <- d_summering_lt |> 
  filter(in_cscf, LU == 240) |> 
  group_by(zone, repetition) |> 
  summarise(year = mean(FJ),
            hectares = n(),
            .groups = "drop")

# interpolate
d_sommerung_zone_inter <- d_GL_zone_inter |> 
  select(zone, year) |> 
  full_join(d_summering_smry,
            by = c("zone", "year")) |> 
  arrange(zone, year) |> 
  group_by(zone) |> 
  mutate(hectares = approx(year, hectares, year)$y,
         # first values for all before AS 1985
         hectares = ifelse(year < min(year[!is.na(repetition)]),
                           hectares[year == min(year[!is.na(repetition)])],
                           hectares),
         # last value for all after AS 2018
         hectares = ifelse(year > max(year[!is.na(repetition)]),
                           hectares[year == max(year[!is.na(repetition)])],
                           hectares)) |> 
  ungroup() |> 
  filter(is.na(repetition)) |> 
  select(-repetition)

d_GL_zone_inter_comb <- d_GL_zone_inter |> 
  mutate(type = "BFS") |> 
  bind_rows(d_sommerung_zone_inter |> 
              rename(value = hectares) |> 
              mutate(type = "AS"))

# zonal summary ----------------------------------------------------------------.

d_LM <-
  d_GL_zone_inter_comb |> 
  select(zone, year) |> 
  distinct() |> 
  left_join(d_GL_zone_inter_comb |> 
              group_by(zone, year) |> 
              summarise(value = sum(value),
                        .groups = "drop") |> 
              left_join(d_prop_zones, by = "zone") |> 
              mutate(value_rel = value / (n  * 25)) |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  filter(year %in% year_insect_start, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  mutate(year_end = year_end - year_start - (length_interval - 1) / 2) |> # center year such that intercept refers to estimated mean
  summarise(lm_coef = lm(value_rel ~ year_end)$coefficients[2],
            lm_intercept = lm(value_rel ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval - 1)/2,
         pred_start = lm_intercept - (length_interval - 1) / 2 * lm_coef,
         pred_end = lm_intercept + (length_interval - 1) / 2 * lm_coef) 

d_LM2 <-
  d_GL_zone_inter_comb |> 
  select(zone, year) |> 
  distinct() |> 
  left_join(d_GL_zone_inter_comb |> 
              group_by(zone, year) |> 
              summarise(value = sum(value),
                        .groups = "drop") |> 
              left_join(d_prop_zones, by = "zone") |> 
              mutate(value_rel = value / (n  * 25)) |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  filter(year %in% year_insect_start2, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval2 & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  mutate(year_end = year_end - year_start - (length_interval2 - 1) / 2) |> # center year such that intercept refers to estimated mean
  summarise(lm_coef = lm(value_rel ~ year_end)$coefficients[2],
            lm_intercept = lm(value_rel ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval2 - 1)/2,
         pred_start = lm_intercept - (length_interval2 - 1) / 2 * lm_coef,
         pred_end = lm_intercept + (length_interval2 - 1) / 2 * lm_coef) 

d_drivers <- d_drivers |> 
  full_join(d_LM |> 
              mutate(length_interval = length_interval) |> 
              select(zone, year_mean, length_interval, 
                     grassland_area_change = lm_coef) |> 
              bind_rows(d_LM2 |> 
                          mutate(length_interval = length_interval2) |> 
                          select(zone, year_mean, length_interval, 
                                 grassland_area_change = lm_coef)),
            by = c("zone", "year_mean", "length_interval"), relationship = "one-to-one") |> 
  full_join(d_LM |> 
              select(zone, year_mean, grassland_area_abs = lm_intercept) |> 
              mutate(length_interval = length_interval) |> 
              bind_rows(d_LM2 |> 
                          select(zone, year_mean, grassland_area_abs = lm_intercept) |> 
                          mutate(length_interval = length_interval2)),
            by = c("zone", "year_mean", "length_interval"), relationship = "one-to-one")

# plotting ---------------------------------------------------------------------.

d_GL_zone_inter |> 
  left_join(d_GL_zone |> 
              select(-value) |> 
              mutate(Interpolated = F),
            by = c("zone", "year")) |> 
  mutate(Interpolated = ifelse(is.na(Interpolated), T, Interpolated)) |> 
  left_join(d_prop_zones, by = "zone") |> 
  mutate(value_rel = value / (25 * n),
         zone = factor(v_zones[zone], level = v_zones),
         source = "Managed grasslands\n(agricultural statistics)") |>
  filter(year >= 1930) |> 
  ggplot(aes(x = year, y = value_rel)) +
  geom_tile(data = d_drivers |> 
              filter(length_interval == 8) |> 
              mutate(grassland_area_change = grassland_area_change / 
                       (2 * sd(grassland_area_change)),
                     zone = factor(v_zones[zone], level = v_zones)), 
            aes(x = year_mean, y = 1, fill = grassland_area_change),
            height = Inf, width = length_interval - 2,
            colour = "grey20", alpha = .8) +
  geom_point(aes(size = Interpolated, colour = source)) +
  geom_point(data = d_sommerung_zone_inter |> 
               left_join(d_prop_zones, by = "zone") |> 
               mutate(value_rel = hectares / (25 * n),
                      zone = factor(v_zones[zone], level = v_zones),
                      Interpolated = T,
                      source = "Summering pastures\n(land-use statistics)"),
             aes(colour = source, size = Interpolated)) +
  geom_point(data = d_summering_smry |> 
               left_join(d_prop_zones, by = "zone") |> 
               mutate(value_rel = hectares / (25 * n),
                      zone = factor(v_zones[zone], level = v_zones),
                      Interpolated = F,
                      source = "Summering pastures\n(land-use statistics)"),
             aes(colour = source, size = Interpolated)) +
  geom_segment(data = d_LM |> 
                 mutate(zone = factor(v_zones[zone], level = v_zones),
                        pred_start = pred_start,
                        pred_end = pred_end), 
               aes(x = year_start, xend = year_start + length_interval - 1,
                   y = pred_start, yend = pred_end), 
               linewidth = .75, alpha = .7, colour = "sienna1") +
  facet_wrap(~ zone) +
  scale_fill_gradient2(low = "#003c30", high = "#543005",
                       name = "Std.\ngrassland value\nchange") +
  labs(y = "Grassland value (%)", x = "Year") +
  scale_x_continuous(breaks = c(1950, 1975, 2000)) +
  scale_size_manual(values = c(.75, .15), name = "Inter-/Extrapolated") +
  scale_colour_manual(values = c("mediumpurple", "mediumpurple4"),
                      name = "Source") +
  theme(axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        strip.text = element_text(size = v_textsize["axis.title"]),
        legend.title = element_text(size = v_textsize["legend.title"]),
        legend.text = element_text(size = v_textsize["legend.text"]))

ggsave("Output/Figures/Drivers_grassland_area.jpeg", width = 10, height = 5.5)

# Average yearly temperature ###################################################
################################################################################.

# read data --------------------------------------------------------------------.

# extract longitude & latitude
lon <- ncvar_get(nc_temperature, "lon")
lat <- ncvar_get(nc_temperature, "lat", verbose = F)

# extract time
t <- ncvar_get(nc_temperature, "time")
t <- as.Date("1800-01-01") %m+% months(t)

# extract data (T, etc)
ndvi_array <- ncvar_get(nc_temperature, "TrecabsM1901") # store the data in a 3-dimensional array

# loop through months
f_append <- \(i){
  
  ndvi_slice <- ndvi_array[, , i] 
  
  data.frame(lon = rep(lon, ncol(ndvi_slice)),
             lat =  rep(lat, each = nrow(ndvi_slice)),
             Trec = as.vector(ndvi_slice),
             date = t[i]) |> 
    filter(!is.na(Trec))
  
}

cl <- makeCluster(6)
clusterExport(cl, c("ndvi_array", "lon", "lat", "t"))
clusterEvalQ(cl, {library(tidyverse)})
d_trec <- parLapplyLB(cl, 1:length(t), f_append) |> 
  bind_rows()
stopCluster(cl)


d_trec <- d_trec %>% 
  mutate(month = month(date),
         year = year(date))

# add LV95n coordinates
d_trec[, c("CX", "CY")] <- sf_project(from = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                                      to = "+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=2600000 +y_0=1200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs",
                                      cbind(d_trec$lon, d_trec$lat))

# Summarise data per year and zone ---------------------------------------------.

d_trec_agg <- d_trec %>% 
  mutate(n_days = days_in_month(date)) %>% 
  group_by(lon, lat, CX, CY, year) %>% 
  summarise(Trec_m = sum(Trec * n_days) / sum(n_days),
            .groups = "drop")


# create sf polygons
radius_x <- 1600 / 2
radius_y <- 2320 / 2
radius_x <- mean(unique(diff(lon))) / 2
radius_y <- mean(unique(diff(lat))) / 2


sf_trec_agg_sub <-
  d_trec_agg %>% 
  filter(year == 1901) %>% # dummy year
  mutate(ID = paste(lon, lat),
         xlow = lon - radius_x,
         ylow = lat - radius_y,
         xhigh = lon + radius_x,
         yhigh = lat + radius_y,
         X1 = xlow, Y1 = ylow,
         X2 = xlow, Y2 = yhigh,
         X3 = xhigh, Y3 = yhigh,
         X4 = xhigh, Y4 = ylow,
         X5 = xlow, Y5 = ylow) %>% 
  pivot_longer(c(X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5)) %>% 
  mutate(dim = substr(name, 1, 1),
         nr = substr(name, 2, 2)) %>% 
  select(-name) %>% 
  spread(dim, value) %>% 
  sf_polygon(polygon_id = "ID", x = "X", y = "Y") 

st_crs(sf_trec_agg_sub) <- 4326 # define coordinate reference system
sf_trec_agg_sub$Trec_m <- d_trec_agg$Trec_m[d_trec_agg$year == 1901]

sf_trec_agg_sub <- st_transform(sf_trec_agg_sub, crs = 2056)

sf_trec_agg_sub <- sf_trec_agg_sub %>%
  st_join(sf_zones %>%
            group_by(zone, in_cscf) %>%
            summarise(.groups = "drop"), 
          join = st_intersects, largest = T)

# add to overall dataset
d_trec_agg <- d_trec_agg %>% 
  mutate(ID = paste(lon, lat)) %>% 
  left_join(sf_trec_agg_sub %>% 
              as.data.frame() %>% 
              select(ID, zone, in_cscf),
            by = "ID") 

d_trec_zone_yr <- d_trec_agg %>% 
  filter(!is.na(zone),
         in_cscf) %>% 
  group_by(zone, year) %>% 
  summarise(Trec_m = mean(Trec_m),
            .groups = "drop") 

# zonal summary ----------------------------------------------------------------.

d_LM <-
  d_trec_zone_yr |>
  select(-c(Trec_m)) |> 
  left_join(d_trec_zone_yr |> 
              group_by(zone) |> 
              mutate(Trec_m_anom = Trec_m - mean(Trec_m[year %in% 1901:2000])) |> 
              ungroup() |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  filter(year %in% year_insect_start, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  mutate(year_end = year_end - year_start - (length_interval - 1) / 2) |> # center year such that intercept refers to estimated mean
  summarise(lm_coef = lm(Trec_m_anom ~ year_end)$coefficients[2],
            lm_intercept = lm(Trec_m_anom ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval - 1)/2,
         pred_start = lm_intercept -+ (length_interval - 1) / 2 * lm_coef,
         pred_end = lm_intercept + (length_interval - 1) / 2 * lm_coef) 

d_LM2 <-
  d_trec_zone_yr |> 
  select(-c(Trec_m)) |> 
  left_join(d_trec_zone_yr |> 
              group_by(zone) |> 
              mutate(Trec_m_anom = Trec_m - mean(Trec_m[year %in% 1901:2000])) |> 
              ungroup() |> 
              rename(year_end = year),
            by = "zone", relationship = "many-to-many") |> 
  filter(year %in% year_insect_start2, # reduce to consecutive intervals (sharing always two observations because of double years for insects)
         year_end < year + length_interval2 & 
           year_end >= year) |>
  group_by(zone, year_start = year) %>% 
  mutate(year_end = year_end - year_start - (length_interval2 - 1) / 2) |> # center year such that intercept refers to estimated mean
  summarise(lm_coef = lm(Trec_m_anom ~ year_end)$coefficients[2],
            lm_intercept = lm(Trec_m_anom ~ year_end)$coefficients[1],
            .groups = "drop") |> 
  mutate(year_mean = year_start + (length_interval2 - 1)/2,
         pred_start = lm_intercept - (length_interval2 - 1) / 2 * lm_coef,
         pred_end = lm_intercept + (length_interval2 - 1) / 2 * lm_coef) 

d_drivers <- d_drivers |> 
  full_join(d_LM |> 
              mutate(length_interval = length_interval) |> 
              select(zone, year_mean, length_interval, 
                     temperature_change = lm_coef) |> 
              bind_rows(d_LM2 |> 
                          mutate(length_interval = length_interval2) |> 
                          select(zone, year_mean, length_interval, 
                                 temperature_change = lm_coef)),
            by = c("zone", "year_mean", "length_interval"), relationship = "one-to-one") |> 
  full_join(d_LM |> 
              select(zone, year_mean, temperature_abs = lm_intercept) |> 
              mutate(length_interval = length_interval) |> 
              bind_rows(d_LM2 |> 
                          select(zone, year_mean, temperature_abs = lm_intercept) |> 
                          mutate(length_interval = length_interval2)),
            by = c("zone", "year_mean", "length_interval"), relationship = "one-to-one")

# plotting ---------------------------------------------------------------------.

d_trec_zone_yr |>
  group_by(zone) |> 
  mutate(Trec_m_anom = Trec_m - mean(Trec_m[year %in% 1901:2000])) |> 
  ungroup() |> 
  filter(year >= 1930) |> 
  mutate(zone = factor(v_zones[zone], level = v_zones)) |> 
  ggplot(aes(x = year, y = Trec_m_anom)) +
  geom_tile(data = d_drivers |> 
              filter(length_interval == 8) |> 
              mutate(temperature_change = temperature_change / 
                       (2 * sd(temperature_change)),
                     zone = factor(v_zones[zone], level = v_zones)), 
            aes(x = year_mean, y = 1, fill = temperature_change),
            height = Inf, width = length_interval - 2,
            colour = "grey20", alpha= .8) +
  geom_point(size = .75) +
  geom_segment(data = d_LM |> 
                 mutate(zone = factor(v_zones[zone], level = v_zones),
                        pred_start = pred_start,
                        pred_end = pred_end), 
               aes(x = year_start, xend = year_start + length_interval - 1,
                   y = pred_start, yend = pred_end), 
               linewidth = .75, alpha = .7, colour = "sienna1") +
  facet_wrap(~ zone) +
  scale_fill_gradient2(low = "#003c30", high = "#543005",
                       name = "Std.\nmean temperature\nchange") +
  labs(y = "Mean temperature anomaly (1901–2000)\n(°C)", x = "Year") +
  scale_x_continuous(breaks = c(1950, 1975, 2000)) +
  theme(axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        strip.text = element_text(size = v_textsize["axis.title"]),
        legend.title = element_text(size = v_textsize["legend.title"]),
        legend.text = element_text(size = v_textsize["legend.text"]))

ggsave("Output/Figures/Drivers_temperature.jpeg", width = 10, height = 5.5)

################################################################################.
# EXPORT DRIVERS DATASET  ------- ##############################################
################################################################################.

d_drivers <- d_drivers |>
  select(zone, year_mean, length_interval, everything()) |> 
  arrange(zone, year_mean, length_interval)
  
d_drivers |> 
  fwrite("Data/Drivers/Drivers.csv")

