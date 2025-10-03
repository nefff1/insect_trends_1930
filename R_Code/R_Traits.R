# setup ########################################################################
################################################################################.

# R version 4.3.3

library(tidyverse) #2.0.0
library(raster) #3.6-26
library(sf) #1.0-15
library(stars) #0.6-6
library(giscoR) #0.6.0
library(data.table) #1.15.4

select <- dplyr::select

sel_species_sapro <- readLines("Data/splist_sapro.txt")
sel_species_butter <- readLines("Data/splist_butter.txt")

# temperature niches saproxylic beetles ########################################
################################################################################.

# data on Europe ---------------------------------------------------------------.

# European raster CGRS:

# CGRS grid available from https://www.eea.europa.eu/data-and-maps/data/common-european-chorological-grid-reference-system-cgrs
sf_cgrs <- st_read(dsn = "Data/Other/cgrs_grid/cgrs_grid.shp")
sf_cgrs <- sf_cgrs[-11267, ] 

# European countries (NUTS):

sf_europe <- gisco_get_countries(region = "Europe")

# to cut off overseas areas
cut_pol <-  st_sfc(st_polygon(list(cbind(c(-5, 40, 45, -40, -5),
                                         c(25, 30, 75, 65, 25)))),
                   crs = 4326) %>%
  st_transform(crs = st_crs(sf_europe)) # harmonize coordinate system


# cut off overseas areas and other countries
sf_europe <- st_intersection(sf_europe, cut_pol)
sf_europe <- sf_europe %>%
  filter(!NAME_ENGL %in% c("Russian Federation", "Belarus", "Ukraine", "Iceland",
                           "Svalbard and Jan Mayen", "Faroes")) %>% 
  st_transform(crs = st_crs(sf_cgrs))

# Europen subset of CGRS grid
sf_cgrs_europe <- sf_cgrs %>%
  st_filter(sf_europe)

# read and process Wordclim data -----------------------------------------------.

# data available from https://www.worldclim.org/data/worldclim21.html

# mean number of days in the included months
d_days <- data.frame()
for (mn_i in 1:12){
  out <- c()
  for (yr_i in 1970:2020){
    out <- c(out, lubridate::days_in_month(as.Date(paste(yr_i, mn_i, 1, sep = "-"))))
  }
  d_days <- data.frame(month = mn_i, days = mean(out)) %>%
    bind_rows(d_days, .)
}


# read through the twelve monthly files
dir <- "Data/Other/Worldclim/wc2.1_2.5m_tavg"
cut <- st_buffer(st_union(sf_cgrs_europe), 100000)

r_wordclim <- raster(list.files(dir, full.names = T)[1]) # first month
r_wordclim <- crop(r_wordclim, extent(st_bbox(cut)))
r_wordclim <- r_wordclim * d_days$days[1]

# add other months
for (i in 2:12){
  file <- list.files(dir, full.names = T)[i]
  
  r_target <- raster(file)
  r_target <- crop(r_target, extent(st_bbox(cut)))
  
  r_wordclim <- r_wordclim + r_target *  d_days$days[i]
}

r_wordclim <- r_wordclim / sum(d_days$days) # mean temperature


# convert to sf
sf_wordclim <- st_as_sf(st_as_stars(r_wordclim))
sf_wordclim <- st_transform(sf_wordclim, crs = st_crs(sf_cgrs_europe))

names(sf_wordclim)[1] <- "Tniche"

# Average T over CGRS grid:

sf_wordclim_cgrs <- st_join(sf_wordclim, sf_cgrs_europe,
                            join = st_within) # choose only those completely within CGRS grid. Pragmatic and time-saving choice

sf_wordclim_cgrs <- sf_wordclim_cgrs %>%
  filter(!is.na(CGRSNAME)) # points outside grid / between two grid cells

# mean across grid
sf_wordclim_cgrs <- sf_wordclim_cgrs %>%
  group_by(CGRSNAME) %>%
  summarise(Tniche = mean(Tniche),
            n = n())

# add to CGRS data
d_cgrs_europe <- sf_cgrs_europe %>%
  as.data.frame() %>%
  left_join(sf_wordclim_cgrs, by = "CGRSNAME")

sf_cgrs_europe$Tniche <- d_cgrs_europe$Tniche
sf_cgrs_europe <- sf_cgrs_europe %>%
  filter(!is.na(Tniche))


# read GBIF data ---------------------------------------------------------------.
# available from https://doi.org/10.15468/dl.mak332

d_gbif_full <- fread("Data/Other/GBIF/0041029-240906103802322/occurrence.csv")
d_std <- fread("Data/Other/GBIF/d_taxonomy_GBIF_sapro.csv")
# use taxonomy of present study
d_gbif_full <- d_gbif_full |> 
  left_join(d_std,
            join_by(taxonKey))

# list of species to be analysed
specieslist <- sel_species_sapro[sel_species_sapro %in% d_gbif_full$species_std] # no data for two species


sf_europe_buff <- st_buffer(sf_europe, 20000)

d_climatic_niche_sapro <- data.frame()
for (sp_i in specieslist){
  
  d_gbif <- d_gbif_full %>%
    filter(species_std == sp_i)
  
  # transform to sf object
  sf_gbif <- d_gbif %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
             crs = "epsg:4121") %>%
    st_transform(crs = st_crs(sf_cgrs_europe))
  
  # exclude non-european observations
  sel <- st_intersects(sf_gbif, sf_europe_buff, sparse = T)
  sel <- which(unlist(lapply(sel, length)) > 0)
  sf_gbif <- sf_gbif[sel, ]
  
  sf_gbif <- st_join(sf_gbif, sf_cgrs_europe)
  
  out <- sf_gbif %>%
    as.data.frame() %>%
    select(CGRSNAME, Tniche) %>%
    distinct() %>%
    summarise(Tniche_mean = mean(Tniche, na.rm = T),
              Tniche_sd = sd(Tniche, na.rm = T),
              Tniche_q05 = quantile(Tniche, .05, na.rm = T),
              Tniche_q95 = quantile(Tniche, .95, na.rm = T))
  
  out$species_std <- sp_i
  
  if (class(out) == "data.frame"){
    d_climatic_niche_sapro <- d_climatic_niche_sapro %>%
      bind_rows(out)
  }
  
  
  
}

d_climatic_niche_sapro <- d_climatic_niche_sapro %>%
  arrange(species_std)

# temperature niches butterflies ###############################################
################################################################################.

# use data generated in Neff et al. 2022 Nature Communications

d_climatic_niche_butter <- readRDS("Data/Traits/RAW/Neff_et_al_2022_NatComms/d_climatic_niche.rds")
d_climatic_niche_butter <- d_climatic_niche_butter |> 
  mutate(species = case_when(species == "Melitaea athalia aggr"~ "Melitaea athalia aggr.", 
                             species %in% c("Colias alfacariensis", "Colias hyale") ~ "Colias hyale aggr.",
                             species %in% c("Leptidea juvernica", "Leptidea sinapis") ~ "Leptidea sinapis aggr.",
                             species == "Colias crocea"~ "Colias croceus",
                             species == "Erebia albergana"~ "Erebia alberganus",
                             species == "Erebia montana"~ "Erebia montanus",
                             species == "Eumedonia eumedon"~ "Aricia eumedon",
                             .default = species)) |> 
  group_by(species) |> # aggregates (2)
  summarise(Tniche_mean = mean(Tniche_mean),
            .groups = "drop")


# other traits saproxylic beetles ##############################################
################################################################################.

# trait data collection from several sources (see manuscript)
d_traits_sapro_raw <- fread('Data/Traits/RAW/beetle_traits_RAW.csv')
d_std_sapro <- fread('Data/Traits/RAW/beetle_taxonomy.csv')

d_traits_sapro <-
  d_traits_sapro_raw |> 
  mutate(SpeciesID = paste(trimws(Genre), trimws(Espèce)),
         SpeciesID = ifelse(SpeciesID == "Abraeus globosus (=perpusillus)",
                            "Abraeus perpusillus", SpeciesID)) |> 
  left_join(d_std_sapro, by = c(SpeciesID = "Name_original")) |> 
  select(Name_std, Family, Genus, Species,
         Saproxylique, Saproxylique.facultatif, `Taille.(mm)`,
         Guildes.trophiques, Microhabitat, Taille.préférentielle.du.bois.mort,
         Saproxylation.préférentielle.du.bois.mort, Héliophilie, Hygrophilie,
         Essences.hôtes.préférentielles, `Groupe.d'essences.hôtes`) |> 
  # for aggregates with unclear value, take the one from the more abundant species (based on info fauna data)
  mutate(Hygrophilie = ifelse(Name_std == "Ampedus nigrinus aggr.",
                              "Meso [Mésophile/indifférent]", Hygrophilie),
         Hygrophilie = ifelse(Name_std == "Ampedus pomorum aggr.",
                              "Meso [Mésophile/indifférent]", Hygrophilie)) |> 
  # manual combination
  mutate(`Groupe.d'essences.hôtes` = ifelse(Name_std == "Ampedus nigrinus aggr.",
                                            "FeuRes [essences hôtes feuillues et résineuses]", `Groupe.d'essences.hôtes`),
         `Groupe.d'essences.hôtes` = ifelse(Name_std == "Ampedus pomorum aggr.",
                                            "Feu(Res) [essences hôtes surtout feuillues, parfois résineuses]", `Groupe.d'essences.hôtes`)) |> 
  group_by(Name_std, Family, Genus, Species) |>
  summarise(across(`Taille.(mm)`, ~ mean(.)),
            across(c(Saproxylique, Saproxylique.facultatif,
                     Guildes.trophiques, Hygrophilie,
                     Saproxylation.préférentielle.du.bois.mort,
                     `Groupe.d'essences.hôtes`), ~ unique(.)),
            across(c(Taille.préférentielle.du.bois.mort, Héliophilie),
                   ~ paste(sort(unique(trimws(strsplit(., "/")[[1]]))),
                           collapse = "/")),
            across(c(Microhabitat),
                   ~ paste(sort(unique(trimws(strsplit(., ",")[[1]]))),
                           collapse = ",")),
            across(c(Essences.hôtes.préférentielles),
                   ~ ifelse(any(is.na(.)), .[!is.na(.)],
                            paste(sort(unique(gsub("*", "", trimws(strsplit(., ",")[[1]])))),
                                  collapse = ","))),
            .groups = "drop") 


# standardise and translate naming
d_traits_sapro <-
  d_traits_sapro |> 
  select(-c(Saproxylique, Saproxylique.facultatif)) |>
  rename(size = `Taille.(mm)`,
         feeding_guild = Guildes.trophiques,
         hygrophilia = Hygrophilie,
         heliophilia = Héliophilie,
         tree_group = `Groupe.d'essences.hôtes`,
         diameter = Taille.préférentielle.du.bois.mort,
         microhabitat = Microhabitat,
         hosts = Essences.hôtes.préférentielles) |>
  rowwise() |>
  mutate(feeding_guild = paste(sort(c(ifelse(grepl("Zoo ", feeding_guild), "zoophagous", NA),
                                      ifelse(grepl("Xyl ", feeding_guild), "xylophagous", NA),
                                      ifelse(grepl("Sxy ", feeding_guild), "saproxylophagous", NA),
                                      ifelse(grepl("Sap ", feeding_guild), "saporphagous", NA),
                                      ifelse(grepl("Myc ", feeding_guild), "mycetophagous", NA),
                                      ifelse(grepl("Rhizo ", feeding_guild), "rhiziphagous", NA))),
                               collapse = "/")) |>
  ungroup() |>
  mutate(hygrophilia = recode(hygrophilia,
                              "Meso [Mésophile/indifférent]" = "mesophilous",
                              "Xer [Xérophile]" = "xerophilous",
                              "Hyg [hygrophile]" = "hygrophilous"),
         heliophilia = recode(heliophilia,
                              "Hel [Héliophile]" = "heliophilous",
                              "Hel [Héliophile]/Sci [Sciaphile]" = "unspecific",
                              "Sci [Sciaphile]" = "sciaphilous"),
         tree_group = recode(tree_group,
                             "Feu (essences hôtes feuillues)" = "deciduous",
                             "Feu(Res) [essences hôtes surtout feuillues, parfois résineuses]" = "mostly_deciduous",
                             "FeuRes [essences hôtes feuillues et résineuses]" = "both",
                             "Res(Feu) [essences hôtes surtout résineuses, parfois feuillues]" = "mostly_coniferous",
                             "Res [Essences hôtes résineuses]" = "coniferous"),
         diameter = recode(diameter,
                           "PBM [bois mort de petit diamètre <10 cm]" = "thin",
                           "GBM [bois mort de gros diamètre >40 cm]" = "thick",
                           "GBM [bois mort de gros diamètre >40 cm]/PBM [bois mort de petit diamètre <10 cm]" = "unspecific")) |>
  # create a stenotopy measure based on habitat requirements
  # (the more specialised requirements, the higher the stenotopy)
  rowwise() |> 
  mutate(n_specific = sum(hygrophilia %in% c("hygrophilous", "xerophilous"),
                          heliophilia %in% c("sciaphilous", "heliophilous"),
                          tree_group %in% c("coniferous", "deciduous"),
                          diameter %in% c("thin", "thick"))) |> 
  ungroup() |> 
  mutate(stenotopy = case_when(n_specific == 4 ~ "stenotopic",
                               n_specific == 3 ~ "oligotopic",
                               n_specific < 3 ~ "eurytopic")) |> 
  select(-n_specific) |> 
  full_join(d_climatic_niche_sapro, by = c("Name_std" = "species_std"), relationship = "one-to-one") |> 
  filter(Name_std %in% sel_species_sapro) |>  # Cis submicans was not included in the analyses (too few years)
  # categorize continuous traits
  mutate(size_cat = cut(log(size), # take log to get more equal group size
                        breaks = 3, labels = c("small", "medium", "large")),
         Tniche = cut(Tniche_mean, breaks = 3, labels = c("cold", "intermediate", "warm")))


# food specialisation ----------------------------------------------------------.

d_german_sl <- fread("Data/Traits/RAW/GermanSL_edit.csv") # plant name standardisation
v_monotypic <- c(Abies = "Abies alba",
                 Aesculus = "Aesculus hippocastanum",
                 Ailanthus = "Ailanthus altissima",
                 Carpinus = "Carpinus betulus",
                 Castanea = "Castanea sativa",
                 Corylus = "Corylus avellana",
                 Fagus = "Fagus sylvatica",
                 Ficus = "Ficus carica",
                 Frangula = "Frangula alnus",
                 Hedera = "Hedera helix",
                 Ilex = "Ilex aquifolium",
                 Juglans = "Juglans regia",
                 Larix = "Larix decidua",
                 Ligustrum = "Ligustrum vulgare",
                 Mespilus = "Mespilus germanica",
                 Myricaria = "Myricaria germanica",
                 Ostrya = "Ostrya carpinifolia",
                 Picea = "Picea abies",
                 Pseudotsuga = "Pseudotsuga menziesii",
                 Robinia = "Robinia pseudoacacia",
                 Spartium = "Spartium junceum",
                 Syringa = "Syringa vulgaris",
                 Taxus = "Taxus baccata",
                 Ulex = "Ulex europaeus",
                 Viscum = "Viscum album")

d_foodspec_sapro <- d_traits_sapro_raw |> 
  mutate(SpeciesID = paste(trimws(Genre), trimws(Espèce)),
         SpeciesID = ifelse(SpeciesID == "Abraeus globosus (=perpusillus)",
                            "Abraeus perpusillus", SpeciesID)) |> 
  left_join(d_std_sapro, by = c(SpeciesID = "Name_original")) |> 
  select(Name_std, Hosts = Essences.hôtes.préférentielles) |>
  mutate(Hosts = gsub(", ", ",", Hosts),
         Hosts = gsub("\\(aussi ", ",", Hosts),
         Hosts = gsub("\\*", "", Hosts),
         Hosts = gsub("\\(que Pinus selon Bouget\\)", "", Hosts),
         Hosts = gsub("Salix Populus", "Salix,Populus", Hosts),
         Hosts = gsub("\\)", "", Hosts),
         Hosts = trimws(Hosts)) |> 
  separate_rows(Hosts, sep = ",") |> 
  mutate(Hosts = trimws(Hosts),
         Hosts = gsub(" spp.", "", Hosts),
         Hosts = recode(Hosts,
                        Caprinus = "Carpinus",
                        "Castanea sativa." = "Castanea sativa",
                        "Malus sylvetris" = "Malus sylvestris",
                        "Ulmus campestris" = "Ulmus glabra",
                        "Myrica" = "Myricaria",
                        "Pinus sylvetris" = "Pinus sylvestris",
                        "Populus sp." = "Populus",
                        "Quercus spp" = "Quercus",
                        "Sarothamnus" = "Cytisus",
                        "Taxinus" = "Taxus",
                        "Tillus" = "Tilia",
                        "pyrus" = "Pyrus",
                        "robinia" = "Robinia",
                        "Labrunum" = "Laburnum")) |> 
  # exclude non relevant and faulty cases:
  filter(!Hosts %in% c("Baccharis", "Eucalyptus", "Magnolia", "Nerium", 
                       "Tamarix", "Tsuga", "Zelkova", 
                       "dendron", "spp.", "", 
                       "Prunus malus"),
         !is.na(Hosts)) |> 
  rowwise() |> 
  mutate(Hosts = ifelse(Hosts %in% names(v_monotypic),
                        v_monotypic[Hosts], Hosts)) |> 
  ungroup() |> 
  left_join(d_german_sl |> select(Name_original, Name_std_host = Name_std, Rank, Family, Genus, Species), 
            by = c(Hosts = "Name_original")) |> 
  group_by(Name_std) |> 
  summarise(n_family = n_distinct(Family),
            n_genus = n_distinct(Genus),
            n_species = n_distinct(Species),
            Hosts = paste(sort(unique(Hosts)), collapse = "/"),
            .groups = "drop") |> 
  mutate(food_spec = case_when(n_family > 1 ~ "polyphagous",
                               n_genus > 1 ~ "oligophagous",
                               n_species > 1 ~ "monophagous",
                               n_species == 1 ~ "monophagous")) |> 
  filter(Name_std %in% sel_species_sapro)  

d_traits_sapro <- d_traits_sapro |> 
  full_join(d_foodspec_sapro |> select(Name_std, food_spec), 
            by = "Name_std", relationship = "one-to-one")

d_traits_sapro |> 
  select(Name_std, size_cat, stenotopy, food_spec, Tniche) |> 
  fwrite("Data/Traits/Traits_saproxylic_beetles.csv")

# other traits butterflies #####################################################
################################################################################.

d_std_butter <- fread('Data/Traits/RAW/butterfly_taxonomy.csv')
  
# collect data from different sources (see manuscript)
d_traits_middleton <- fread('Data/Traits/RAW/European_&_Maghreb_Butterfly_Trait_data_v1.2.csv')
d_traits_middleton <- d_traits_middleton |> 
  mutate(Taxa.name = gsub("_", " ", Taxa.name)) |> 
  left_join(d_std_butter, by = c("Taxa.name" = "Name_original"))
d_traits_cook <- fread('Data/Traits/RAW/Cook_et_al_ecological_traits.csv')
d_traits_cook <- d_traits_cook |> 
  left_join(d_std_butter, by = c("scientific_name" = "Name_original"))

# Klaiber et al (2017) data:
d_traits_butter <- fread('Data/Traits/RAW/Fauna_Indicativa_Tagfalter_et_Zygaenidae2022.csv')

d_traits_butter <- d_traits_butter |> 
  mutate(Taxa.Unterart = ifelse(Taxa.Unterart == "", NA, Taxa.Unterart),
         SpeciesID = paste(Taxa.Gattung, Taxa.Art),
         SpeciesID = ifelse(!is.na(Taxa.Unterart),
                            paste(SpeciesID, Taxa.Unterart),
                            SpeciesID)) |> 
  left_join(d_std_butter, 
            by = c(SpeciesID = "Name_original"))


d_traits_butter <- d_traits_butter |> 
  filter(Name_std %in% sel_species_butter) |>
  rename(dry_meadow = `TWW Kennartstatus.Kennartstatus`,
         stenotopy = Biotopbindung.Kürzel,
         food_spec = `Phagie der Raupe.Kürzel`,
         hibernation = `Überwinterung.Stadium, welches überwintert`,
         voltinism = `Voltinismus.Generationen pro Jahr in der Schweiz`,
         wing_length = `Morphologie.Flügel.Länge (mm)`) |> 
  mutate(stenotopy = recode(stenotopy,
                            E = "eurytopic",
                            O = "oligotopic",
                            S = "stenotopic"),
         food_spec = recode(food_spec,
                            M = "monophagous",
                            NO = "narrowly oligophagous",
                            O = "oligophagous",
                            P = "polyphagous"),
          hibernation = recode(hibernation,
                              Imago = "adult",
                              Puppe = "pupa",
                              Raupe = "larva",
                              Ei = "egg",
                              "Raupe in Eihülle" = "egg",
                              "(Raupe)" = "larva",
                              "Raupe, Raupe in Eihülle" = "egg/larva",
                              "Raupe, Puppe" = "larva/pupa",
                              "Ei, Raupe" = "egg/larva"),
         hibernation = ifelse(hibernation %in% c("NA", "keine"), NA, hibernation)) |> 
  select(Name_std, SpeciesID, dry_meadow, stenotopy, food_spec,
         hibernation, voltinism, wing_length) |> 
  mutate(across(everything(), ~ ifelse(. == "NA", NA, .)),
         wing_length = as.numeric(wing_length)) |> 
  mutate(wing_length = ifelse(Name_std %in% c("Parnassius phoebus", "Apatura iris"), NA, wing_length)) |>  # apparent mistakes in data
  left_join(d_traits_middleton |> 
              group_by(Name_std) |> 
              summarise(wing_length_midd = mean((FoL_var_male_average + FoL_var_female_average) / 2,
                                                na.rm = T),
                        .groups = "drop"),
            by = "Name_std", relationship = "many-to-one") |> 
  mutate(wing_length = ifelse(is.na(wing_length), wing_length_midd, wing_length)) |> 
  select(-wing_length_midd) |> 
  left_join(d_traits_cook |> 
              mutate(wing_length_cook = ifelse(Family == "Zygaenidae",
                                               (forewing_minimum + forewing_maximum) / 2,
                                               (forewing_minimum + forewing_maximum) / 4)) |> # weirdly doubled values...
              group_by(Name_std) |> 
              summarise(wing_length_cook = mean(wing_length_cook, na.rm = T),
                        .groups = "drop"),
            by = "Name_std", relationship = "many-to-one") |> 
  mutate(wing_length = ifelse(is.na(wing_length), wing_length_cook, wing_length)) |> 
  select(-wing_length_cook) |> 
  left_join(d_climatic_niche_butter, by = c("Name_std" = "species"),
            relationship = "many-to-one") |> 
  group_by(Name_std) |> 
  summarise(across(-SpeciesID,
                   ~ifelse(is.numeric(.), mean(.), paste(unique(.), collapse = "__")))) |> 
  mutate(stenotopy = recode(stenotopy, 
                            oligotopic__stenotopic = "stenotopic"), # Phengaris alcon rebeli more common
         food_spec = recode(food_spec,
                            "narrowly oligophagous__monophagous" = "monophagous")) |>  # Phengaris alcon rebeli more common
  mutate(food_spec = ifelse(food_spec == "narrowly oligophagous", "monophagous", food_spec)) |> 
  mutate((across(everything(), ~ ifelse(. == "NA", NA, .))))

# categorise categorical traits:
d_traits_butter <- d_traits_butter |> 
  mutate(size_cat = cut(log(wing_length), # take log to get more equal group size
                              breaks = 3, labels = c("small", "medium", "large")),
         Tniche = cut(Tniche_mean, breaks = 3, labels = c("cold", "intermediate", "warm")))

d_traits_butter |>
  select(Name_std, size_cat, stenotopy, food_spec, Tniche, voltinism, hibernation) |> 
  fwrite("Data/Traits/Traits_butterflies.csv")
