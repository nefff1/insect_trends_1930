# Initialise system ############################################################
################################################################################.

# R version 4.0.2

rm(list = ls()); graphics.off()
Sys.setlocale('LC_ALL','C.UTF-8')
options(max.print = 500)

library(tidyverse) #1.3.1
library(cmdstanr); set_cmdstan_path(path = '/home/neff/cmdstan-2.27.0') #0.4.0
library(parallel)
library(data.table) #1.14.0
library(bayestestR) #0.10.0
library(posterior) #0.1.5

# Rhat function (edited from the rstan package function 'Rhat')
f_Rhat <- function(sims) {
  
  f_tomatrix <- function(obj_draws){
    matrix(as.numeric(obj_draws), ncol = 8, byrow = F) # 4 chains split in two
  }
  
  f_rhat_rfun <- function(sims) {
    chains <- ncol(sims)
    n_samples <- nrow(sims)
    chain_mean <- numeric(chains)
    chain_var <- numeric(chains)
    for (i in seq_len(chains)) {
      chain_mean[i] <- mean(sims[, i])
      chain_var[i] <- var(sims[, i])
    }
    var_between <- n_samples * var(chain_mean)
    var_within <- mean(chain_var)
    sqrt((var_between/var_within + n_samples - 1)/n_samples)
  }
  
  f_z_scale <- function(x){
    S <- length(x)
    r <- rank(x, ties.method = 'average')
    z <- qnorm((r - 1/2)/S)
    z[is.na(x)] <- NA
    if (!is.null(dim(x))) {
      z <- array(z, dim = dim(x), dimnames = dimnames(x))
    }
    z
  }
  
  bulk_rhat <- f_rhat_rfun(f_z_scale(f_tomatrix(sims)))
  sims_folded <- abs(sims - median(sims))
  tail_rhat <- f_rhat_rfun(f_z_scale(f_tomatrix(sims_folded)))
  max(bulk_rhat, tail_rhat)
}


source('R_Code/f_occ_det.R')
stan_mod <- cmdstan_model('Stan_Code/Stan_occ_det_cmdstan.stan')

options(mc.cores = 4)

# information on the model run (from cluster)
sp_index_i <- as.numeric(commandArgs(trailingOnly = TRUE)[1]) # species index
iter <- as.numeric(commandArgs(trailingOnly = TRUE)[2]) # number of iterations to run


d_records <- fread('Data/RAW_occdet/d_records_butter.csv')

sel_species <- read_lines("Data/splist_butter.txt")

sp_i <- sel_species[sp_index_i]

# projects with restricted species focus
d_rest_projects <- fread('Data/RAW_occdet/Projects_restricted_focus_butter.csv')
d_other_projects <- fread("Data/RAW_occdet/Projects_other_focus_butter.csv")

# expert persons / projects
d_experts_butter <- fread('Data/RAW_occdet/d_experts_butter.csv')
d_expertprojs_butter <- fread('Data/RAW_occdet/d_expertprojs_butter.csv')

# biogeographic regions data
d_biogeo <- fread('Data/RAW_occdet/Biogeographische_regionen_five_km2.csv')
d_biogeo <- d_biogeo %>% 
  mutate(five_km2_ID = paste(CX_five_km2, CY_five_km2, sep = ' / ')) |> 
  mutate(biogeo6 = recode(biogeo6,
                          Jura = "Jura",
                          Mittelland = "Plateau",
                          Alpennordflanke = "NorthernAlps",
                          "Westliche Zentralalpen" = "WesternCentralAlps",
                          "Östliche Zentralalpen" = "EasternCentralAlps",
                          Alpensüdflanke = "SouthernAlps"))

# height above sealevel data
d_height <- fread('Data/RAW_occdet/mean_height_km2.csv')
d_height <- d_height %>% 
  mutate(CX_five_km2 = floor(CX_km2 / 5000) * 5000,
         CY_five_km2 = floor(CY_km2 / 5000) * 5000) %>% 
  group_by(CX_five_km2, CY_five_km2) %>% 
  summarise(height = sum(Z * n) / sum(n),
            .groups = "drop") %>% 
  mutate(five_km2_ID = paste(CX_five_km2, CY_five_km2, sep = ' / '))

d_zones <- d_biogeo %>% 
  left_join(d_height, by = c("CX_five_km2", "CY_five_km2", "five_km2_ID")) %>% 
  mutate(zone = ifelse(height > 1400 | biogeo12 == "Engadin", "HighAlps",
                       ifelse(grepl("CentralAlps", biogeo6),
                              "CentralAlps", biogeo6)))

d_sites <- expand.grid(five_km2_ID = unique(d_records$five_km2_ID),
                       two_A = seq(min(d_records$two_A), max(d_records$two_A), 2)) %>% 
  arrange(five_km2_ID, two_A)

d_visits <-
  d_cscf %>% 
  mutate(observer = ifelse(PROJET == "", LEG, PROJET)) %>% 
  group_by(observer) %>% 
  mutate(n_visit = n_distinct(visit_ID)) %>% 
  ungroup() %>% 
  left_join(d_rest_projects, by = 'PROJET') %>%
  rowwise() %>%
  filter((!Name_std %in% strsplit(focal_species, ' \\| ')[[1]]) |
           any(strsplit(focal_species, ' \\| ')[[1]] == sp_i)) %>%
  ungroup() %>%
  group_by(five_km2_ID, two_A, visit_ID) %>% 
  summarise(pres = ifelse(any(Name_std == sp_i), 1, 0),
            list_length = length(unique(Name_std)),
            list_length_cat = cut(list_length, breaks = c(0, 1, 3, Inf),
                                  labels = c('single', 'short', 'long')),
            # determine source category
            source = ifelse(sp_i %in% unlist(strsplit(focal_species, ' \\| ')),
                            'Project_targeted',
                            ifelse(any(PROJET %in% c('LRPAP')), 
                                   'Project_RL',
                                   ifelse(any(PROJET != '') & !any(PROJET %in% d_other_projects$PROJET),
                                          ifelse(PROJET %in% d_expertprojs_butter$PROJET, 'Project expert',
                                                 'Project'),
                                          ifelse(any(LEG %in% d_experts_butter$LEG),
                                                 'CitSc_expert', 'CitSc')))),
            museum_record = as.numeric(all(MUS != "")),
            observer = ifelse(unique(n_visit) >= 15, unique(observer), "OTHER"), 
            .groups = 'drop') %>% 
  mutate(five_km2_ID = as.factor(five_km2_ID),
         fact_two_A = as.factor(two_A),
         source = as.factor(source),
         museum_record = as.factor(museum_record)) %>% 
  arrange(five_km2_ID, two_A, visit_ID)


# add environmental data #######################################################
################################################################################.

d_sites <- d_sites %>% 
  left_join(d_biogeo %>% 
              select(five_km2_ID, biogeo6, biogeo12), by = 'five_km2_ID') %>% 
  left_join(d_zones %>% 
              select(five_km2_ID, zone), by = c("five_km2_ID")) %>% 
  mutate_at(vars(biogeo6, biogeo12, zone), as.factor) %>% 
  left_join(d_height %>% 
              filter(five_km2_ID %in% d_sites$five_km2_ID) %>% 
              select(five_km2_ID, height) %>% 
              mutate(height_z = as.numeric(scale(height))), by = 'five_km2_ID')

# subset based on zones ########################################################
################################################################################.

# only consider zones, where species was reported at least once

sel_zone <- d_visits %>%
  filter(pres == 1) %>%
  left_join(d_sites, by = 'five_km2_ID') %>%
  select(zone) %>%
  distinct()

d_sites <- d_sites %>%
  filter(zone %in% sel_zone$zone) %>%
  droplevels()

d_visits <- d_visits %>%
  filter(five_km2_ID %in% d_sites$five_km2_ID)

d_sites <- d_sites %>% 
  arrange(five_km2_ID, two_A)

d_visits <- d_visits %>% 
  arrange(five_km2_ID, two_A)

print(paste('N sites =', nrow(d_sites), '| N visits =', nrow(d_visits)))

# run MCMC chains ##############################################################
################################################################################.

f_occ_det(d_sites = d_sites, d_visits = d_visits, 
          formula_occ = ~ height_z + I(height_z^2) + (1 | biogeo12) + (1 | five_km2_ID),
          formula_det = ~ list_length_cat + source + museum_record + (1 | fact_two_A) + (1 | observer),
          var_site = 'five_km2_ID', var_year = 'two_A', 
          var_presabs = 'pres',
          var_area = 'zone',
          stan_mod = stan_mod,
          iter_warmup = floor(iter / 2), iter_sampling = iter - floor(iter / 2),
          chains = 4,
          output_dir = 'Output',
          output_basename = gsub(' ', '_', sp_i))

# post-sampling processing #####################################################
################################################################################.

# move to enclosing folder

dir.create(paste0('Output/Done/', gsub(' ', '_', sp_i)))

system(paste0('mv Output/', gsub(' ', '_', sp_i), '-1.csv',
              ' Output/Done/', gsub(' ', '_', sp_i)))
system(paste0('mv Output/', gsub(' ', '_', sp_i), '-2.csv',
              ' Output/Done/', gsub(' ', '_', sp_i)))
system(paste0('mv Output/', gsub(' ', '_', sp_i), '-3.csv',
              ' Output/Done/', gsub(' ', '_', sp_i)))
system(paste0('mv Output/', gsub(' ', '_', sp_i), '-4.csv',
              ' Output/Done/', gsub(' ', '_', sp_i)))

# create folder for final data -------------------------------------------------.

dir <- paste0('Output/occdet/Butter_', iter, '/')
dir.create(paste0(dir, gsub(' ', '_', sp_i)))

# export site IDs for matching -------------------------------------------------.
d_sites %>% 
  select(five_km2_ID, two_A) %>% 
  fwrite(paste0(dir, gsub(' ', '_', sp_i), '/SITEIDS_', gsub(' ', '_', sp_i), '.csv'))


# extract parameter estimates and write to file --------------------------------

files <- list.files(paste0('/home/neff/INSECT/Output/Done/', gsub(' ', '_', sp_i)),
                    full.names = T)
files <- files[!grepl("SITEIDS_", files)]

fit <- read_cmdstan_csv(files,
                        variables = c('mu_o', 'alpha_year', 'sigma_a_year',
                                      'sigma_a_ro',
                                      'beta_fo',
                                      'mu_d', 'alpha_rd', 'sigma_a_rd',
                                      'beta_fd', 'z'), 
                        format = 'draws_list')


#  estimates of a selection of parameters --------------------------------------.

pars <- c('z', 'alpha_year', 'beta_fo', 'mu_o')

out1 <- lapply(fit$post_warmup_draws, function(x) x[, grepl(paste(paste0('^', pars), 
                                                                  collapse = '|'), 
                                                            names(x)), drop = F]) %>% 
  bind_rows() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("var")

names(out1)[-1] <- paste0('I', 1:(ncol(out1) - 1))
fwrite(out1, paste0(dir, '/', gsub(' ', '_', sp_i), '/RAW_', gsub(' ', '_', sp_i), '.csv'))

# aggregate data per zone and two-year interval --------------------------------.


d_mcmc_z <- out1 %>% 
  filter(substr(var, 1, 2) == "z[") %>% 
  column_to_rownames("var")

d_sites <- d_sites %>% 
  select(five_km2_ID, two_A) %>% 
  left_join(d_zones %>% 
              select(five_km2_ID, biogeo6, biogeo12, zone), by = "five_km2_ID")

d_z_intmean_reg <- data.frame()
for (i in unique(d_sites$two_A)){
  for (j in unique(d_sites$zone)) {
    sel <- which(d_sites$two_A == i & d_sites$zone == j)
    
    d_z_intmean_reg <- data.frame(two_A = i, zone = j, iter = 1:ncol(d_mcmc_z),
                                  z_mean = apply(d_mcmc_z[sel, ], 2, mean)) %>%
      bind_rows(d_z_intmean_reg, .)
  }
}


out_agg <- d_z_intmean_reg |> 
  mutate(iter = paste0("occ_m_i", formatC(iter, width = 4, flag = 0)),
         Name_std = sp_i) |> 
  pivot_wider(values_from = z_mean, names_from = iter) |> 
  left_join(d_sites |> 
              group_by(zone, two_A) |> 
              summarise(n_squares = n_distinct(five_km2_ID),
                        .groups = "drop"),
            join_by(zone, two_A),
            relationship = "one-to-one") |> 
  select(Name_std, two_A, zone, n_squares, everything())

fwrite(out_agg,
       paste0(dir_data, "Data/Occupancy_estimates/occ_means_1930_butter.csv"),
       append = T)


# Rhat values for yearly mean occupancy values ---------------------------------.

if (iter < 500) { # doesn't run for high numbers of iterations
  
  f_diag <- function(dat, n_chains = 4){
    mat <- dat %>% 
      as.numeric() %>% 
      matrix(ncol = n_chains, byrow = F)
    
    data.frame(Rhat = f_Rhat(mat),
               ess_bulk = ess_bulk(mat),
               ess_tail = ess_tail(mat))
  }
  
  
  d_mcmc_z <- out1 %>% 
    filter(substr(var, 1, 2) == "z[") %>% 
    column_to_rownames("var")
  
  d_diag <- data.frame()
  for (i in unique(d_sites$two_A)){
    
    sel <- which(d_sites$two_A == i)
    
    diag <- apply(d_mcmc_z[sel, ], 2, mean) %>% 
      f_diag
    
    d_diag <- data.frame(two_A = i,
                         diag) %>% 
      bind_rows(d_diag, .)
  }
  
  
  f_diag_zone <- function(j){
    out <- data.frame()
    for (i in unique(d_sites$two_A)){
      sel <- which(d_sites$two_A == i & d_sites$zone == j)
      
      diag <- apply(d_mcmc_z[sel, ], 2, mean) %>% 
        f_diag()
      
      out <- data.frame(two_A = i,
                        zone = j,
                        diag) %>% 
        bind_rows(out, .)
    }
    
    out
  }
  
  d_diag_zone <- mclapply(unique(d_sites$zone), f_diag_zone, mc.cores = 4) %>% 
    bind_rows()
  
  d_diag <- d_diag %>% 
    bind_rows(d_diag_zone) %>% 
    mutate(Name_std = sp_i)
  
  fwrite(d_diag, paste0(dir, '/', gsub(' ', '_', sp_i), '/DIAGZ_', gsub(' ', '_', sp_i), '.csv'))
}

