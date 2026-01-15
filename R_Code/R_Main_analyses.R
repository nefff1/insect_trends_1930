################################################################################.
# SETUP -------------------- ###################################################
################################################################################.

# R version 4.3.3

# load packages ################################################################
################################################################################.

library(data.table) #1.15.4
library(tidyverse); theme_set(theme_classic()) #2.0.0
library(parallel)
library(bayestestR) #0.16.0
library(DHARMa) #0.4.6
library(ggh4x) #0.2.8
library(brms) #2.21.0
library(ggpubr) #0.6.0
library(cowplot) #1.1.3
library(grid)

# set global parameters ########################################################
################################################################################.

# insect group colours and labels (for plotting)
v_groupcols <- c(Sapro = "#043565", Butter = "#EB4B98")
v_grouplabels <- c(Sapro = "Saproxylic beetles", Butter = "Butterflies")

# names and order of biogeographic zones
v_zones <- c(Jura = "Jura",
             Plateau = "Plateau", 
             NorthernAlps = "Northern Alps",
             CentralAlps = "Central Alps",
             SouthernAlps = "Southern Alps",
             HighAlps = "High Alps")

# same including all of Switzerland
v_zones_v2 <- c(Switzerland = "Overall",
                Jura = "Jura",
                Plateau = "Plateau", 
                NorthernAlps = "Northern Alps",
                CentralAlps = "Central Alps",
                SouthernAlps = "Southern Alps",
                HighAlps = "High Alps")


# main traits saproxylic beetles
sel_traits_sapro <- c(Size = "size_cat", 
                      "Habitat specialisation" = "stenotopy",
                      "Food specialisation" = "food_spec",
                      "Temperature niche" = "Tniche")

# main traits butterflies
sel_traits_butter <- c(Size = "size_cat", 
                      "Habitat specialisation" = "stenotopy",
                      "Food specialisation" = "food_spec",
                      "Temperature niche" = "Tniche",
                      "Overw. stage" = "hibernation",
                      Voltinism = "voltinism")

# traitvalues for all main traits (both groups)
d_traitsel_comb <- data.frame(trait = "Size",
                              traitvalue = c("small", "medium", "large"),
                              traitvalue_short = c("small", "medium", "large")) |> 
  bind_rows(data.frame(trait = "Habitat specialisation",
                       traitvalue = c("stenotopic", "oligotopic", "eurytopic"),
                       traitvalue_short = c("stenotopic", "oligotopic", "eurytopic"))) |> 
  bind_rows(data.frame(trait = "Food specialisation",
                       traitvalue = c("monophagous", "oligophagous", "polyphagous"),
                       traitvalue_short = c("monoph.", "oligoph.", "polyph."))) |> 
  bind_rows(data.frame(trait = "Temperature niche",
                       traitvalue = c("cold", "intermediate", "warm"),
                       traitvalue_short = c("cold", "interm.", "warm")))

# traitvalues for all traits (butterflies)
d_traitsel_butter <-  data.frame(trait = "Overw. stage",
                                 traitvalue = c("egg", "egg/larva", "larva", "pupa", "adult"),
                                 traitvalue_short = c("egg", "egg/larva", "larva", "pupa", "adult")) |> 
  bind_rows(data.frame(trait = "Voltinism",
                       traitvalue = c("0.5", "1", "2", "3"),
                       traitvalue_short = c("0.5", "1", "2", "3"))) |> 
  bind_rows(data.frame(trait = "Size",
                       traitvalue = c("small", "medium", "large"),
                       traitvalue_short = c("small", "medium", "large"))) |> 
  bind_rows(data.frame(trait = "Habitat specialisation",
                       traitvalue = c("stenotopic", "oligotopic", "eurytopic"),
                       traitvalue_short = c("stenotopic", "oligotopic", "eurytopic"))) |> 
  bind_rows(data.frame(trait = "Food specialisation",
                       traitvalue = c("monophagous", "oligophagous", "polyphagous"),
                       traitvalue_short = c("monoph.", "oligoph.", "polyph."))) |> 
  bind_rows(data.frame(trait = "Temperature niche",
                       traitvalue = c("cold", "intermediate", "warm"),
                       traitvalue_short = c("cold", "interm.", "warm")))

# order and (long) naming of environmental variables
v_vars_spec <- c("Temperature abs." = "temperature_abs",
                 "Temperature change" = "temperature_change",
                 "Human population density change" = "human_pop_change",
                 "Mechanisation abs." = "mechanisation_abs",
                 "Mechanisation period" = "mechanisation_period",
                 "Grassland area abs." = "grassland_area_abs",
                 "Grassland area change" = "grassland_area_change",
                 "Wood harvest intensity change" = "wood_harvest_int_change",
                 "Storm aftermath" = "after_storm",
                 "Forest area abs." = "forest_area_abs",
                 "Forest area change" = "forest_area_change")

# order and (short) naming of environmental variables
v_vars_spec_short <- c("Temp." = "temperature_abs",
                       "Δ\u00A0Temp." = "temperature_change",
                       "Δ\u00A0Human pop.\u00A0dens." = "human_pop_change",
                       "Mechanis." = "mechanisation_abs",
                       "Mechanis. period" = "mechanisation_period",
                       "Grassl. area" = "grassland_area_abs",
                       "Δ\u00A0Grassl. area" = "grassland_area_change",
                       "Δ\u00A0W.\u00A0harvest intensity" = "wood_harvest_int_change",
                       "Storm aftermath" = "after_storm",
                       "Forest area" = "forest_area_abs",
                       "Δ\u00A0Forest area" = "forest_area_change")

n_sim <- 5000 # number of Monte Carlo simulations

n_iter <- 20000 # number of iterations in regression models

# interval length at which trends and drivers are related
length_interval <- 8 # years (standard length)
length_interval2 <- 12 # years (for sensitivity analyses)

# text size (for plotting)
v_textsize <- c(plotlabel = 9,
                axis.title = 8, axis.text = 7, 
                legend.title = 8, legend.text = 7,
                additional.text = 6) 


# define functions #############################################################
################################################################################.

# runs (and saves) Monte Carlo richness simulations for a single species
f_ric_sim <- function(sp_i){
  
  # filter mean occupancy estimates for the target species
  x <- d_occ_means |> filter(Name_std == sp_i) |> select(-n_squares)
  
  # draw iteration number (n_sim times)
  d_sim <- data.frame(draw = seq_len(n_sim),
                      iter = sample(seq_len(sum(grepl("occ_m_i", names(x)))),
                                    n_sim,
                                    replace = T))
  # add data for the simulations
  d_sim <- d_sim |> 
    left_join(x |> 
                pivot_longer(contains("occ_m_i"),
                             names_to = "iter",
                             values_to = "occ_m") |> 
                mutate(iter = as.integer(gsub("occ_m_i", "", iter))),
              join_by(iter),
              relationship = "many-to-many") |> 
    select(-iter)
  
  # output the data
  fwrite(d_sim,
         paste0("Data/tmp/sric_sim/sric_sim_", gsub(" |\\.", "_", sp_i, ".csv")))
}

# functionruns a linear model and outputs the coefficients and pvalues
f_lmcoefs<- function(formula){
  mod <- lm(formula)
  lm_coef <- mod$coefficients[2]
  pval <- summary(mod)$coefficients[2, 4]
  data.frame(lm_coef, pval)
}

# function to predict the change in trends when a variable (var_i) is increased 
# by a certain amount (default to 0.5 SD) (while keeping all other variables constant)
# based on regression models relating trends and drivers
f_pred_overall <- function(group_i, var_i, diff_i = NA, ci = .95){
  
  if (is.numeric(d_drivers_raw[, var_i])){ # for numeric environmental variables
    sd_i <- d_drivers_raw |> 
      filter(length_interval == sel_interval) |> 
      select(sym(var_i)) |> 
      deframe() |> 
      sd()
    
    if (is.na(diff_i)){ # 0.5 sd by default
      level1_i <- 0
      level2_i <- 0.5
      
      diff_i <- .5 * sd_i
    } else {
      
      diff_z <- diff_i / sd_i
      
      level1_i <- 0
      level2_i <- diff_z
    }
  } else { # for non-numeric environmental variables (which are all T/F here)
    level1_i <- FALSE
    level2_i <- TRUE
  } 
  
  if (group_i == "butter"){ # model data frame and model for butterflies
    d_mod_i <- d_mod_butter
    mod_i <- mod_main_butter
  } else if (group_i == "sapro"){ # model data frame and model for saproxylic beetles
    d_mod_i <- d_mod_sapro
    mod_i <- mod_main_sapro
  }
  
  
  d_new <-
    d_mod_i |> 
    # use mean values for all numeric, non-target variables (to get marginal effects)
    select(!!! syms(unname(v_vars_spec))) |> 
    select(- sym(var_i)) |> 
    summarise_all(~ ifelse(is.numeric(.), mean(.), .[1])) |> 
    # add target variable at two levels
    expand_grid(!! sym(var_i) := c(level1_i, level2_i)) |>
    # add missing model variables which are not of interest
    mutate(year_mean = NA, zone = NA)
  
  # posterior predictions for the generated dataframe
  m_pred <- posterior_epred(mod_i,
                            newdata = d_new) 
  
  # calculate means and 95% credible intervals
  c_out <- c(mean(apply(m_pred, 1, diff)),
             ci(apply(m_pred, 1, diff), ci = ci)[-1] |> as.numeric())
  
  names(c_out) <- c("mean", paste0("lower", ci), paste0("upper", ci))
  
  # add difference at original variable scale for numeric variables (i.e., what does 0.5SD relate to)
  if (is.numeric(d_drivers_raw[, var_i])){
    c_out <- c(c_out, diff_i)
    names(c_out)[4] <- "diff_original"
  }
  
  c_out
}

# function that run a regression model relating richness trends to environmental variables and their interactions with traitvalues
# the trait (e.g. body size) and the group (e.g. saproxylic beetles) can be specified
f_traitmod <- function(trait_i, group_i){
  # prepare list for outputs
  out <- list()
  out$plotsdata <- list()
  
  # load model dataset and trait specifications
  if (group_i == "sapro"){
    d_mod_trait <- d_mod_sapro_trait
    d_traitsel <- d_traitsel_comb
  } else if (group_i == "butter"){
    d_mod_trait <- d_mod_butter_trait
    d_traitsel <- d_traitsel_butter
  }
  
  # model dataframe
  d_target <- d_mod_trait |> 
    filter(trait == trait_i,
           # remove cases with no traitvalues or very rare traitvalues (i.e. larva/pupa)
           !is.na(traitvalue),
           traitvalue != "larva/pupa",
           traitvalue != "") |> 
    rowwise() |> 
    # convert to short traitvalue format
    mutate(traitvalue = d_traitsel$traitvalue_short[d_traitsel$traitvalue == traitvalue]) |> 
    ungroup() |> 
    # convert to factor
    mutate(traitvalue = factor(traitvalue, 
                               levels = d_traitsel$traitvalue_short[d_traitsel$trait == names(sel_traits_butter[sel_traits_butter == trait_i])]))
  
  
  
  set.seed(92) # set seed for MCMC sampling
  mod_target <- d_target |>
    # model formula:
    (\(x) brm(mean ~ (
      temperature_abs + temperature_change  +
        human_pop_change +
        mechanisation_abs + mechanisation_period +
        grassland_area_abs + grassland_area_change +
        forest_area_abs + forest_area_change  + 
        wood_harvest_int_change + after_storm
    ) * traitvalue  +
      (1 | year_mean) + (1 | zone) + (1 | year_mean:zone),
    # priors:
    prior = prior(normal(0,5), class = b) +
      prior(normal(0,5), class = Intercept) +
      prior(cauchy(0,1), class = sd) +
      prior(cauchy(0,25), class = sigma), 
    data = x,
    iter = n_iter, # number of iterations
    family = student(),
    file = paste0("Models/mod_",
                  group_i, "_trait_", trait_i, "_",
                  sel_interval, "yrs_", 
                  n_iter, "iters.rds")))()
  
  # create dataframe for posterior model prediction of overall effect of traitvalues
  d_new_overall <- d_target |> 
    ungroup() |> 
    # take mean values for all env. variables
    summarise(across(c(temperature_abs, temperature_change,
                       human_pop_change,
                       mechanisation_abs,
                       grassland_area_abs, grassland_area_change,
                       wood_harvest_int_change, 
                       forest_area_abs, forest_area_change),
                     ~ mean(.))) |> 
    # for factors, take the (two) different valeus
    expand_grid(mechanisation_period = sort(unique(d_target$mechanisation_period))) |> 
    expand_grid(after_storm = sort(unique(d_target$after_storm)))
  
  # add all traitvalues
  d_new_overall <- expand_grid(d_new_overall, traitvalue = levels(d_target$traitvalue))
  
  # posterior prediction
  d_pred_overall_raw <- posterior_epred(mod_target, d_new_overall,
                                        re_formula = NA) |> 
    as.data.frame()
  
  names(d_pred_overall_raw) <- d_new_overall$traitvalue
  
  # take mean of factors / logicals / etc.
  d_pred_overall <- lapply(unique(names(d_pred_overall_raw)),
                           function(x) tibble(!! sym(x) := apply(d_pred_overall_raw[, names(d_pred_overall_raw) == x], 
                                                                 1, mean))) |> 
    bind_cols()
  
  # summaries of posterior predictions for different traitvalues
  d_pred_overall <- d_pred_overall |> 
    pivot_longer(everything(), names_to = "traitvalue", values_to = "pred") |> 
    group_by(traitvalue) |> 
    summarise(ci(pred, .95) |> 
                rename(lower95 = CI_low,
                       upper95 = CI_high),
              ci(pred, .9) |> 
                rename(lower90 = CI_low,
                       upper90 = CI_high),
              ci(pred, .8) |> 
                rename(lower80 = CI_low,
                       upper80 = CI_high),
              mean = mean(pred),
              .groups = "drop") |> 
    select(-CI) |> 
    mutate(traitvalue = factor(traitvalue, levels = levels(d_target$traitvalue)))
  
  out$plotsdata[["overall"]] <- d_pred_overall
  
  # now iterate through all env. variables and change one after the other to get posterior predictions for marginal effects
  for (var_i in v_vars_spec){
    
    if (d_target[, var_i, drop = T] |> is.numeric()){ # for numeric variables
      # prepare dataframe for predictions
      d_new <-
        tibble(!! sym(var_i) := seq(min(d_target[, var_i]), max(d_target[, var_i]), length.out = 50)) |> # 50 values from min to max
        expand_grid(d_target |> 
                      ungroup() |> 
                      summarise(across(c(temperature_abs, temperature_change,
                                         human_pop_change,
                                         mechanisation_abs,
                                         grassland_area_abs, grassland_area_change,
                                         wood_harvest_int_change, 
                                         forest_area_abs, forest_area_change),
                                       ~ mean(.))) |> 
                      expand_grid(mechanisation_period = sort(unique(d_target$mechanisation_period))) |> 
                      expand_grid(after_storm = sort(unique(d_target$after_storm))) |> 
                      select(- !! sym(var_i)))
      
      d_new <- expand_grid(d_new, traitvalue = levels(d_target$traitvalue)) |> 
        mutate(traitvalue = factor(traitvalue, levels = levels(d_target$traitvalue)))
    } else { # for factor variables
      # prepare dataframe for predictions
      d_new <-
        tibble(!! sym(var_i) := unique(sort(d_target[, var_i, drop = T]))) |> 
        expand_grid(d_target |> 
                      ungroup() |> 
                      summarise(across(c(temperature_abs, temperature_change,
                                         human_pop_change,
                                         mechanisation_abs,
                                         grassland_area_abs, grassland_area_change,
                                         wood_harvest_int_change, 
                                         forest_area_abs, forest_area_change),
                                       ~ mean(.))) |> 
                      expand_grid(mechanisation_period = sort(unique(d_target$mechanisation_period))) |> 
                      expand_grid(after_storm = sort(unique(d_target$after_storm))) |> 
                      select(- !! sym(var_i)) |> 
                      distinct())
      
      d_new <- expand_grid(d_new, traitvalue = levels(d_target$traitvalue)) |> 
        mutate(traitvalue = factor(traitvalue, levels = levels(d_target$traitvalue)))
    }
    
    d_new <- d_new |> 
      mutate(ID = paste(!! sym(var_i), traitvalue)) |> 
      arrange(ID)
    
    # calculate posterior prediction
    m_pred_raw <- posterior_epred(mod_target, d_new,
                                  re_formula = NA) 
    
    # take mean of factors / logicals / etc.
    m_pred <- sapply(unique(d_new$ID),
                     function(x) matrix(apply(m_pred_raw[, d_new$ID == x, drop = F], 
                                              1, mean), ncol = 1)) |> 
      unname()
    
    d_new_cond <- d_new |> 
      select(!!sym(var_i), traitvalue) |> 
      distinct()
    
    # caluclate means and credible intervals
    d_new_cond$pred_mean <- apply(m_pred, 2, mean)
    d_new_cond$pred_upper <- apply(m_pred, 2, function(x) ci(x, .95)$CI_high)
    d_new_cond$pred_lower <- apply(m_pred, 2, function(x) ci(x, .95)$CI_low)
    
    # save output
    out$plotsdata[[var_i]]$d_new_cond <- d_new_cond
    
    if (d_target[, var_i, drop = T] |> is.numeric()){ # for numeric variables
      
      # prediction at minimum value
      d_pred_min <-
        d_new |> 
        filter(!! sym(var_i) == min(!! sym(var_i))) |> 
        group_by(ID) |> 
        group_map(~ data.frame(pred = apply(posterior_epred(mod_target, .,
                                                            re_formula = NA), 1, mean),
                               traitvalue = unique(.$traitvalue))) |> 
        bind_rows()
      
      d_smry_min <- d_pred_min |> 
        group_by(traitvalue) |> 
        summarise(ci(pred, .95) |> 
                    rename(lower95 = CI_low,
                           upper95 = CI_high),
                  ci(pred, .9) |> 
                    rename(lower90 = CI_low,
                           upper90 = CI_high),
                  ci(pred, .8) |> 
                    rename(lower80 = CI_low,
                           upper80 = CI_high),
                  mean = mean(pred),
                  .groups = "drop") |> 
        select(-CI) 
      
      # prediction at maximum value
      d_pred_max <-
        d_new |> 
        filter(!! sym(var_i) == max(!! sym(var_i))) |> 
        group_by(ID) |> 
        group_map(~ data.frame(pred = apply(posterior_epred(mod_target, .,
                                                            re_formula = NA), 1, mean),
                               traitvalue = unique(.$traitvalue))) |> 
        bind_rows()
      
      d_smry_max <- d_pred_max |> 
        group_by(traitvalue) |> 
        summarise(ci(pred, .95) |> 
                    rename(lower95 = CI_low,
                           upper95 = CI_high),
                  ci(pred, .9) |> 
                    rename(lower90 = CI_low,
                           upper90 = CI_high),
                  ci(pred, .8) |> 
                    rename(lower80 = CI_low,
                           upper80 = CI_high),
                  mean = mean(pred),
                  .groups = "drop") |> 
        select(-CI) 
      
      # collate minimum and maximum value predictions
      d_smry <- d_smry_min |> 
        mutate(what = "min") |> 
        bind_rows(d_smry_max |> 
                    mutate(what = "max"))
      
      # save to output list
      out$plotsdata[[var_i]]$d_smry <- d_smry
      
    } else { # for factor variables
      
      
      d_smry <- data.frame()
      for (value_i in unique(d_new[, var_i, drop = T])){ # prediction for each factor variable
        d_pred_i <-
          d_new |> 
          filter(!! sym(var_i) == value_i) |> 
          group_by(ID) |> 
          group_map(~ data.frame(pred = apply(posterior_epred(mod_target, .,
                                                              re_formula = NA), 1, mean),
                                 traitvalue = unique(.$traitvalue))) |> 
          bind_rows()
        
        d_smry <- d_pred_i |> 
          group_by(traitvalue) |> 
          summarise(ci(pred, .95) |> 
                      rename(lower95 = CI_low,
                             upper95 = CI_high),
                    ci(pred, .9) |> 
                      rename(lower90 = CI_low,
                             upper90 = CI_high),
                    ci(pred, .8) |> 
                      rename(lower80 = CI_low,
                             upper80 = CI_high),
                    mean = mean(pred),
                    .groups = "drop") |> 
          select(-CI) |> 
          mutate(what = value_i) |> 
          (\(x) bind_rows(d_smry, x))()
        
      }
      
      # save to output list
      out$plotsdata[[var_i]]$d_smry <- d_smry
      
    }
    
  }
  
  # output
  out
}

# function rearranges output from f_traitmod for plotting of trait-variable-regression models (Extended data figures)
# trait and output of f_traitmod need to be specified
f_comb_plotsdata <- function(trait_i, l_traitmod_i){
  
  # extract data necessary to plot all the different panels (i.e. overall effects, different environmental variables)
  plotsdata_i <- l_traitmod_i[[trait_i]]$plotsdata
  
  
  # function to get data for lines connecting minimum- and maximum-level predictions
  f_linedata <- function(x){
    x$var <- names(x)[1]
    
    x$type = ifelse(n_distinct(x[, 1]) <= 2, "factor", "cont.")
    
    if (x$type[1] == "cont."){
      mean_raw <- mean(d_drivers_raw[d_drivers_raw$length_interval == sel_interval, names(x)[1]]) # mean of raw values of the variable (not scaled)
      sd_raw <- sd(d_drivers_raw[d_drivers_raw$length_interval == sel_interval, names(x)[1]]) # sd of raw variables of the variable (not scaled)
      
      x[, 1] <- x[, 1] * 2 * sd_raw + mean_raw # transformation of the scaled data (backtransform to original scale)
    }
    
    x <- x[x[, 1, drop = T] %in% range(x[, 1]), ] # restrict to min/max
    names(x)[1] <- "x"
    
    
    
    x <- x |> 
      mutate(traitvalue_num = as.numeric(traitvalue),
             traitvalue_chr = as.character(traitvalue),
             diff_x = diff(range(x)))
    
    x
    
  }
  
  out <- list()
  
  # etract the data for the overall plots of the different traits (i.e. without env. variables)
  # rangedata contains the data for the vertical lineranges (i.e. credible intervals, means)
  out$rangedata <- plotsdata_i$overall |> 
    mutate(traitvalue_num = as.numeric(traitvalue), # traitvalue as numeric value
           traitvalue_chr = as.character(traitvalue), # traitvalue as character value
           x = traitvalue_num, # x axis value
           diff_x = n_distinct(x) - 1, # number of gaps between unique x values
           var = "overall") # variable (for facetting)
  
  # start to compile the data for the plots for each env. variable
  # linedata contains data for the lines connecting minimum- and maximum-level predictions
  out$linedata <- data.frame()
  for (var_i in names(plotsdata_i)[-1]){ # loop through env. variables
    
    if (plotsdata_i[[var_i]]$d_new_cond[, 1, drop = T] |> is.numeric()){ # d_new_cond contains the predictions along the range of numeric values for var_i
      mean_raw <- mean(d_drivers_raw[d_drivers_raw$length_interval == sel_interval, var_i]) # mean of raw environmental variable
      sd_raw <- sd(d_drivers_raw[d_drivers_raw$length_interval == sel_interval, var_i]) # sf of raw environmental variable
      
      out$rangedata <- plotsdata_i[[var_i]]$d_smry |> 
        mutate(traitvalue_num = as.numeric(traitvalue), # traitvalue as numeric value
               traitvalue_chr = as.character(traitvalue), # traitvalue as character value
               x = ifelse(what == "min", # x axis value
                          min(plotsdata_i[[var_i]]$d_new_cond[, 1]),
                          max(plotsdata_i[[var_i]]$d_new_cond[, 1])),
               x = x * 2 * sd_raw + mean_raw, # backtransform to original scale of env. variable
               diff_x = diff(range(x)),  # number of gaps between unique x values
               var = var_i) |>  # variable (for facetting)
        select(-what) |> 
        (\(x) bind_rows(out$rangedata, x))()
      
    } else { # for factors
      out$rangedata <- plotsdata_i[[var_i]]$d_smry |> 
        mutate(traitvalue_num = as.numeric(traitvalue), # traitvalue as numeric value
               traitvalue_chr = as.character(traitvalue), # traitvalue as character value
               x = as.numeric(what), # x axis value
               diff_x = n_distinct(x) - 1, # number of gaps between unique x values
               var = var_i) |>  # variable (for facetting)
        select(-what) |> 
        (\(x) bind_rows(out$rangedata, x))()
    }
    
    
    
    out$linedata <- f_linedata(plotsdata_i[[var_i]]$d_new_cond) |> 
      (\(x) bind_rows(out$linedata, x))()
    
  }
  
  # put together output
  out$rangedata$trait <- trait_i
  out$linedata$trait <- trait_i
  
  out
}

# function to calculate the percentage point and proportional change in richness per specified time interval
# output already summarised to means and credible intervals
f_sric_change <- function(data, length_interval_i = 10){
  
  data |> 
    group_by(draw) |> 
    mutate(sric = sric / sric[two_A == 1930]) |> # value relative to the 1930 richness value
    ungroup() |> 
    # self-join for being able to get differences for different time intervals
    inner_join(data |> 
                group_by(draw) |> 
                mutate(sric = sric / sric[two_A == 1930],
                       two_A = two_A - length_interval_i) |> # shift year values back certain years such that observations that are apart by the defined amount of years are joined
                ungroup(),
              join_by(two_A, draw),
              suffix = c(".start", ".end")) |> 
    mutate(diff_abs = (sric.end - sric.start) * 100, # absolute difference in percentage points
           diff_rel = diff_abs / sric.start, # difference relative to starting condition
           year_end = two_A + length_interval_i) |> 
    select(year_start = two_A, year_end, draw, diff_abs, diff_rel) |> 
    # now summarise across the MCMC iterations
    group_by(year_start, year_end) |> 
    summarise(diff_abs_mean = mean(diff_abs),
              ci(diff_abs, .95) |> 
                rename(diff_abs_lower = CI_low,
                       diff_abs_upper = CI_high),
              diff_rel_mean = mean(diff_rel),
              ci(diff_rel, .95) |> 
                rename(diff_rel_lower = CI_low,
                       diff_rel_upper = CI_high),
              .groups = "drop") |> 
    mutate(across(contains("diff"), ~round(., 1))) |> # round values
    select(year_start, year_end, CI, everything())

}

# load data ####################################################################
################################################################################.

# environmental (change) variables ---------------------------------------------.

d_drivers_raw <- fread("Data/Drivers/Drivers.csv") |> 
  as.data.frame()

# species lists ----------------------------------------------------------------.

splist_sapro <- read_lines("Data/splist_sapro.txt")
splist_butter <- read_lines("Data/splist_butter.txt")

# species traits ---------------------------------------------------------------.

d_traits_sapro <- fread("Data/Traits/Traits_saproxylic_beetles.csv")
d_traits_butter <- fread("Data/Traits/Traits_butterflies.csv") |> 
  mutate(voltinism = as.character(voltinism))

# mean occupancy values --------------------------------------------------------.

# available from Zenodo: https://doi.org/10.5281/zenodo.17255265
d_occ_means <- fread("Data/Occupancy_estimates/occ_means_1930_sapro.csv") |> 
  bind_rows(fread("Data/Occupancy_estimates/occ_means_1930_butter.csv"))|> 
  as.data.frame()

# number of 5x5 km squares included in the saproxylic-beetle occupancy detection models (per zone)
d_nsquares_sapro <- d_occ_means |> 
  filter(Name_std %in% splist_sapro) |> 
  select(zone, n_squares) |> 
  distinct()

# number of 5x5 km squares included in the butterfly occupancy detection models (per zone)
d_nsquares_butter <- d_occ_means |> 
  filter(Name_std %in% splist_butter) |> 
  select(zone, n_squares) |> 
  distinct()

################################################################################.
# PROCESS DATA ------------- ###################################################
################################################################################.

# run richness simulations #####################################################
################################################################################.

# Monte Carlo simulations of richness
# Basically draws occupancy estimates for each species, which can then be aggregated into richness estimates

# run parallelised to speed up

# run for saproxylic beetles ---------------------------------------------------.
set.seed(121)
cl <- makeCluster(8)
clusterExport(cl, c("d_occ_means", "n_sim"))
clusterEvalQ(cl, {library(tidyverse); library(data.table)})
parLapplyLB(cl, splist_sapro, f_ric_sim)
stopCluster(cl)

# run for butterflies ----------------------------------------------------------.
set.seed(175)
cl <- makeCluster(8)
clusterExport(cl, c("d_occ_means", "n_sim"))
clusterEvalQ(cl, {library(tidyverse); library(data.table)})
parLapplyLB(cl, splist_butter, f_ric_sim)
stopCluster(cl)

# calculate richness ###########################################################
################################################################################.

# based on the previous richness simulations, aggregate the richness estimates per Monte Carlo simulation

# saproxylic beetles -----------------------------------------------------------.

# prepare dataframes
d_ric_sapro <- data.frame()
d_ric_reg_sapro <- data.frame()
d_ric_trait_sapro <- data.frame()
d_ric_trait_reg_sapro <- data.frame()

for (sp_i in splist_sapro){ # run through all species
  d_sim <- fread(paste0("Data/tmp/sric_sim/sric_sim_", gsub(" |\\.", "_", sp_i, ".csv"))) # read simulation data
  
  # Switzerland-wide mean occupancy of this species (zones weighted by number of squares analysed)
  d_sim_CH <-
    d_nsquares_sapro |> 
    expand_grid(two_A = seq(1930, 2020, 2),
                Name_std = sp_i,
                draw = seq_len(n_sim)) |> 
    left_join(d_sim,
              join_by(Name_std, zone, two_A, draw),
              relationship = "one-to-many") |> 
    mutate(occ_m = ifelse(is.na(occ_m), 0, occ_m)) |> 
    group_by(Name_std, two_A, draw) |>  
    summarise(sric = sum(occ_m * n_squares) / sum(n_squares),
              .groups = "drop") 
  
  # add to output and summarise (iteratively)
  d_ric_sapro <- d_ric_sapro |> 
    bind_rows(d_sim_CH |> 
                select(-Name_std) |> 
                mutate(count = 1)) |> 
    group_by(two_A, draw) |> 
    summarise(sric = sum(sric),
              count = sum(count),
              .groups = "drop")
  
  # add to zonal output and summarise (iteratively)
  d_ric_reg_sapro <- d_ric_reg_sapro |> 
    bind_rows(d_sim |> 
                mutate(count = 1) |> 
                rename(sric = occ_m) |> 
                select(-Name_std)) |> 
    group_by(zone, two_A, draw) |> 
    summarise(sric = sum(sric),
              count = sum(count),
              .groups = "drop")
  
  # add to trait-grouped output and summarise (iteratively)
  d_ric_trait_sapro <- d_ric_trait_sapro |> 
      bind_rows(d_sim_CH |> 
                  left_join(d_traits_sapro |> 
                              pivot_longer(-Name_std, names_to = "trait", values_to = "traitvalue"),
                            join_by(Name_std),
                            relationship = "many-to-many") |> 
                  mutate(count = 1) |> 
                  select(-Name_std)) |> 
      group_by(two_A, trait, traitvalue, draw) |> 
      summarise(sric = sum(sric),
                count = sum(count),
                .groups = "drop")
  
  # add to zonal and trait-grouped output and summarise (iteratively)
  d_ric_trait_reg_sapro <- d_ric_trait_reg_sapro |> 
    bind_rows(d_sim |> 
                left_join(d_traits_sapro |> 
                            pivot_longer(-Name_std, names_to = "trait", values_to = "traitvalue"),
                          join_by(Name_std),
                          relationship = "many-to-many") |> 
                mutate(count = 1) |> 
                rename(sric = occ_m) |> 
                select(-Name_std)) |> 
    group_by(two_A, zone, trait, traitvalue, draw) |> 
    summarise(sric = sum(sric),
              count = sum(count),
              .groups = "drop")
}

# save for later processing
saveRDS(d_ric_sapro, "Data/tmp/sric_agg/d_ric_sapro.rds")
saveRDS(d_ric_reg_sapro, "Data/tmp/sric_agg/d_ric_reg_sapro.rds")
saveRDS(d_ric_trait_sapro, "Data/tmp/sric_agg/d_ric_trait_sapro.rds")
saveRDS(d_ric_trait_reg_sapro, "Data/tmp/sric_agg/d_ric_trait_reg_sapro.rds")

# butterflies ------------------------------------------------------------------.

# prepare dataframes
d_ric_butter <- data.frame()
d_ric_reg_butter <- data.frame()
d_ric_trait_butter <- data.frame()
d_ric_trait_reg_butter <- data.frame()

for (sp_i in splist_butter){
  d_sim <- fread(paste0("Data/tmp/sric_sim/sric_sim_", gsub(" |\\.", "_", sp_i, ".csv"))) # read simulation data
  
  # Switzerland-wide mean occupancy of this species (zones weighted by number of squares analysed)
  d_sim_CH <-
    d_nsquares_butter |> 
    expand_grid(two_A = seq(1930, 2020, 2),
                Name_std = sp_i,
                draw = seq_len(n_sim)) |> 
    left_join(d_sim,
              join_by(Name_std, zone, two_A, draw),
              relationship = "one-to-many") |> 
    mutate(occ_m = ifelse(is.na(occ_m), 0, occ_m)) |> 
    group_by(Name_std, two_A, draw) |>  
    summarise(sric = sum(occ_m * n_squares) / sum(n_squares),
              .groups = "drop") 
  
  # add to output and summarise (iteratively)
  d_ric_butter <-
    d_ric_butter |> 
    bind_rows(d_sim_CH |> 
                select(-Name_std) |> 
                mutate(count = 1)) |> 
    group_by(two_A, draw) |> 
    summarise(sric = sum(sric),
              count = sum(count),
              .groups = "drop")
  
  # add to zonal output and summarise (iteratively)
  d_ric_reg_butter <- d_ric_reg_butter |> 
    bind_rows(d_sim |> 
                mutate(count = 1) |> 
                rename(sric = occ_m) |> 
                select(-Name_std)) |> 
    group_by(zone, two_A, draw) |> 
    summarise(sric = sum(sric),
              count = sum(count),
              .groups = "drop")
  
  # add to trait-grouped output and summarise (iteratively)
  d_ric_trait_butter <- d_ric_trait_butter |>
    bind_rows(d_sim_CH |>
                left_join(d_traits_butter |>
                            pivot_longer(-Name_std, names_to = "trait", values_to = "traitvalue"),
                          join_by(Name_std),
                          relationship = "many-to-many") |>
                mutate(count = 1) |>
                select(-Name_std)) |>
    group_by(two_A, trait, traitvalue, draw) |>
    summarise(sric = sum(sric),
              count = sum(count),
              .groups = "drop")

  # add to zonal and trait-grouped output and summarise (iteratively)
  d_ric_trait_reg_butter <- d_ric_trait_reg_butter |>
    bind_rows(d_sim |>
                left_join(d_traits_butter |>
                            pivot_longer(-Name_std, names_to = "trait", values_to = "traitvalue"),
                          join_by(Name_std),
                          relationship = "many-to-many") |>
                mutate(count = 1) |>
                rename(sric = occ_m) |>
                select(-Name_std)) |>
    group_by(two_A, zone, trait, traitvalue, draw) |>
    summarise(sric = sum(sric),
              count = sum(count),
              .groups = "drop")
}

# save for later processing
saveRDS(d_ric_butter, "Data/tmp/sric_agg/d_ric_butter.rds")
saveRDS(d_ric_reg_butter, "Data/tmp/sric_agg/d_ric_reg_butter.rds")
saveRDS(d_ric_trait_butter, "Data/tmp/sric_agg/d_ric_trait_butter.rds")
saveRDS(d_ric_trait_reg_butter, "Data/tmp/sric_agg/d_ric_trait_reg_butter.rds")

# aggregate richness ###########################################################
################################################################################.

# so far, richness estimates are at the level of Monte Carlo simulations
# now, aggregate them to get means and credible intervals

# saproxylic beetles -----------------------------------------------------------.

d_ric_sapro <- readRDS("Data/tmp/sric_agg/d_ric_sapro.rds")
d_ric_reg_sapro <- readRDS("Data/tmp/sric_agg/d_ric_reg_sapro.rds")
d_ric_trait_sapro <- readRDS("Data/tmp/sric_agg/d_ric_trait_sapro.rds")
d_ric_trait_reg_sapro <- readRDS("Data/tmp/sric_agg/d_ric_trait_reg_sapro.rds")

d_ric_mean_sapro <- d_ric_sapro |> 
  group_by(two_A) |> 
  summarise(mean = mean(sric),
            lower = ci(sric)$CI_low,
            upper = ci(sric)$CI_high,
            .groups = "drop")

d_ric_reg_mean_sapro <- d_ric_reg_sapro |> 
  group_by(two_A, zone) |> 
  summarise(mean = mean(sric),
            lower = ci(sric)$CI_low,
            upper = ci(sric)$CI_high,
            .groups = "drop")

d_ric_trait_mean_sapro <- d_ric_trait_sapro |> 
  group_by(two_A, trait, traitvalue) |> 
  summarise(mean = mean(sric),
            lower = ci(sric)$CI_low,
            upper = ci(sric)$CI_high,
            .groups = "drop")

d_ric_trait_reg_mean_sapro <- d_ric_trait_reg_sapro |> 
  group_by(two_A, zone, trait, traitvalue) |> 
  summarise(mean = mean(sric),
            lower = ci(sric)$CI_low,
            upper = ci(sric)$CI_high,
            .groups = "drop")

# butterflies ------------------------------------------------------------------.

d_ric_butter <- readRDS("Data/tmp/sric_agg/d_ric_butter.rds")
d_ric_reg_butter <- readRDS("Data/tmp/sric_agg/d_ric_reg_butter.rds")
d_ric_trait_butter <- readRDS("Data/tmp/sric_agg/d_ric_trait_butter.rds")
d_ric_trait_reg_butter <- readRDS("Data/tmp/sric_agg/d_ric_trait_reg_butter.rds")

d_ric_mean_butter <- d_ric_butter |> 
  group_by(two_A) |> 
  summarise(mean = mean(sric),
            lower = ci(sric)$CI_low,
            upper = ci(sric)$CI_high,
            .groups = "drop")

d_ric_reg_mean_butter <- d_ric_reg_butter |> 
  group_by(two_A, zone) |> 
  summarise(mean = mean(sric),
            lower = ci(sric)$CI_low,
            upper = ci(sric)$CI_high,
            .groups = "drop")

d_ric_trait_mean_butter <- d_ric_trait_butter |> 
  group_by(two_A, trait, traitvalue) |> 
  summarise(mean = mean(sric),
            lower = ci(sric)$CI_low,
            upper = ci(sric)$CI_high,
            .groups = "drop")

d_ric_trait_reg_mean_butter <- d_ric_trait_reg_butter |> 
  group_by(two_A, zone, trait, traitvalue) |> 
  summarise(mean = mean(sric),
            lower = ci(sric)$CI_low,
            upper = ci(sric)$CI_high,
            .groups = "drop")



# calculate richness trends ####################################################
################################################################################.

# for each interval (either 8 years or 12 years), calculate the richness trends with linear models
# calculations at the level of single Monte Carlo simulations
# in the end, summarise to means and credible intervals

# saproxylic beetles -----------------------------------------------------------.

d_ric_sapro <- d_ric_sapro |> 
  # richness relative to the 1930 value
  mutate(sric_rel = sric / d_ric_mean_sapro$mean[d_ric_mean_sapro$two_A == 1930])

d_ric_reg_sapro <- d_ric_reg_sapro |>
  # richness relative to the 1930 value
  mutate(sric_rel = sric / d_ric_mean_sapro$mean[d_ric_mean_sapro == 1930])

# Switzerland wide, 8 years
d_lm_ric_interval_sapro <- 
  d_ric_sapro |> 
  select(draw, year_start = two_A) |> 
  left_join(d_ric_sapro,
            join_by(draw), 
            relationship = "many-to-many") |> 
  filter(two_A < year_start + length_interval &  # means length_interval/2 points per regression
           two_A >= year_start) |> 
  mutate(year_mean = year_start + (length_interval - 1)/2) |> 
  group_by(draw, year_mean) %>% 
  filter(n() == floor(length_interval / 2)) |> # exclude cases with not enough observations (at the end)
  summarise(f_lmcoefs(sric_rel ~ two_A),
            .groups = "drop") |> 
  group_by(year_mean) |> 
  summarise(mean  = mean(lm_coef),
            lower = ci(lm_coef)$CI_low,
            upper = ci(lm_coef)$CI_high,
            .groups = "drop")

# Switzerland wide, 12 years
d_lm_ric_interval2_sapro <- 
  d_ric_sapro |> 
  select(draw, year_start = two_A) |> 
  left_join(d_ric_sapro,
            join_by(draw), 
            relationship = "many-to-many") |> 
  filter(two_A < year_start + length_interval2 &  # means length_interval/2 points per regression
           two_A >= year_start) |> 
  mutate(year_mean = year_start + (length_interval2 - 1)/2) |> 
  group_by(draw, year_mean) %>% 
  filter(n() == floor(length_interval2 / 2)) |> # exclude cases with not enough observations (at the end)
  summarise(f_lmcoefs(sric_rel ~ two_A),
            .groups = "drop") |> 
  group_by(year_mean) |> 
  summarise(mean  = mean(lm_coef),
            lower = ci(lm_coef)$CI_low,
            upper = ci(lm_coef)$CI_high,
            .groups = "drop")

# Zonal, 8 years
d_lm_ric_interval_reg_sapro <-
  d_ric_reg_sapro |> 
  select(draw, zone, year_start = two_A) |> 
  left_join(d_ric_reg_sapro,
            join_by(draw, zone), 
            relationship = "many-to-many") |> 
  filter(two_A < year_start + length_interval &  # means length_interval/2 points per regression
           two_A >= year_start) |> 
  mutate(year_mean = year_start + (length_interval - 1)/2) |> 
  group_by(draw, zone, year_mean) %>% 
  filter(n() == floor(length_interval / 2)) |> # exclude cases with not enough observations (at the end)
  summarise(f_lmcoefs(sric_rel ~ two_A),
            .groups = "drop") |>  
  group_by(zone, year_mean) |> 
  summarise(mean  = mean(lm_coef),
            lower = ci(lm_coef)$CI_low,
            upper = ci(lm_coef)$CI_high,
            .groups = "drop")

# Zonal, 12 years
d_lm_ric_interval2_reg_sapro <-
  d_ric_reg_sapro |> 
  select(draw, zone, year_start = two_A) |> 
  left_join(d_ric_reg_sapro,
            join_by(draw, zone), 
            relationship = "many-to-many") |> 
  filter(two_A < year_start + length_interval2 &  # means length_interval/2 points per regression
           two_A >= year_start) |> 
  mutate(year_mean = year_start + (length_interval2 - 1)/2) |> 
  group_by(draw, zone, year_mean) %>% 
  filter(n() == floor(length_interval2 / 2)) |> # exclude cases with not enough observations (at the end)
  summarise(f_lmcoefs(sric_rel ~ two_A),
            .groups = "drop") |>  
  group_by(zone, year_mean) |> 
  summarise(mean  = mean(lm_coef),
            lower = ci(lm_coef)$CI_low,
            upper = ci(lm_coef)$CI_high,
            .groups = "drop")

# for all the different traitvalues
d_lm_ric_interval_trait_reg_sapro <- data.frame()
d_lm_ric_interval2_trait_reg_sapro <- data.frame()
for (trait_i in sel_traits_sapro){
  d_ric_i <- d_ric_trait_sapro |> 
    filter(trait == trait_i)
  
  d_ric_reg_i <- d_ric_trait_reg_sapro |> 
    filter(trait == trait_i)
  
  d_ric_reg_i <- d_ric_reg_i |> 
    filter(traitvalue != "" & !is.na(traitvalue)) |> 
    left_join(d_ric_i |> 
                filter(two_A == 1930) |> 
                group_by(traitvalue) |> 
                summarise(sric_baseline = mean(sric)),
              join_by(traitvalue)) |> 
    mutate(sric_rel = sric / sric_baseline)
  
  d_lm_ric_interval_trait_reg_sapro <- d_ric_reg_i |> 
    select(traitvalue, draw, zone, year_start = two_A) |> 
    left_join(d_ric_reg_i,
              join_by(draw, zone, traitvalue), 
              relationship = "many-to-many") |> 
    filter(two_A < year_start + length_interval &  # means length_interval/2 points per regression
             two_A >= year_start) |> 
    mutate(year_mean = year_start + (length_interval - 1)/2) |> 
    group_by(trait, traitvalue, draw, zone, year_mean) %>% 
    filter(n() == floor(length_interval / 2)) |> # exclude cases with not enough observations (at the end)
    summarise(f_lmcoefs(sric_rel ~ two_A),
              .groups = "drop") |>  
    group_by(trait, traitvalue, zone, year_mean) |> 
    summarise(mean  = mean(lm_coef),
              lower = ci(lm_coef)$CI_low,
              upper = ci(lm_coef)$CI_high,
              .groups = "drop") |> 
    (\(x) bind_rows(d_lm_ric_interval_trait_reg_sapro, x))()
  
  d_lm_ric_interval2_trait_reg_sapro <- d_ric_reg_i |> 
    select(traitvalue, draw, zone, year_start = two_A) |> 
    left_join(d_ric_reg_i,
              join_by(draw, zone, traitvalue), 
              relationship = "many-to-many") |> 
    filter(two_A < year_start + length_interval2 &  # means length_interval/2 points per regression
             two_A >= year_start) |> 
    mutate(year_mean = year_start + (length_interval2 - 1)/2) |> 
    group_by(trait, traitvalue, draw, zone, year_mean) %>% 
    filter(n() == floor(length_interval2 / 2)) |> # exclude cases with not enough observations (at the end)
    summarise(f_lmcoefs(sric_rel ~ two_A),
              .groups = "drop") |>  
    group_by(trait, traitvalue, zone, year_mean) |> 
    summarise(mean  = mean(lm_coef),
              lower = ci(lm_coef)$CI_low,
              upper = ci(lm_coef)$CI_high,
              .groups = "drop") |> 
    (\(x) bind_rows(d_lm_ric_interval2_trait_reg_sapro, x))()
  
}


saveRDS(d_lm_ric_interval_sapro, "Data/tmp/sric_agg/d_lm_ric_interval_sapro.rds")
saveRDS(d_lm_ric_interval2_sapro, "Data/tmp/sric_agg/d_lm_ric_interval2_sapro.rds")
saveRDS(d_lm_ric_interval_reg_sapro, "Data/tmp/sric_agg/d_lm_ric_interval_reg_sapro.rds")
saveRDS(d_lm_ric_interval2_reg_sapro, "Data/tmp/sric_agg/d_lm_ric_interval2_reg_sapro.rds")
saveRDS(d_lm_ric_interval_trait_reg_sapro, "Data/tmp/sric_agg/d_lm_ric_interval_trait_reg_sapro.rds")
saveRDS(d_lm_ric_interval2_trait_reg_sapro, "Data/tmp/sric_agg/d_lm_ric_interval2_trait_reg_sapro.rds")


# butterflies ------------------------------------------------------------------.

d_ric_butter <- d_ric_butter |> 
  # richness relative to the 1930 value
  mutate(sric_rel = sric / d_ric_mean_butter$mean[d_ric_mean_butter$two_A == 1930])

d_ric_reg_butter <- d_ric_reg_butter |>
  # richness relative to the 1930 value
  mutate(sric_rel = sric / d_ric_mean_butter$mean[d_ric_mean_butter == 1930])

# Switzerland wide, 8 years
d_lm_ric_interval_butter <- 
  d_ric_butter |> 
  select(draw, year_start = two_A) |> 
  left_join(d_ric_butter,
            join_by(draw), 
            relationship = "many-to-many") |> 
  filter(two_A < year_start + length_interval &  # means length_interval/2 points per regression
           two_A >= year_start) |> 
  mutate(year_mean = year_start + (length_interval - 1)/2) |> 
  group_by(draw, year_mean) %>% 
  filter(n() == floor(length_interval / 2)) |> # exclude cases with not enough observations (at the end)
  summarise(f_lmcoefs(sric_rel ~ two_A),
            .groups = "drop") |> 
  group_by(year_mean) |> 
  summarise(mean  = mean(lm_coef),
            lower = ci(lm_coef)$CI_low,
            upper = ci(lm_coef)$CI_high,
            .groups = "drop")

# Switzerland wide, 12 years
d_lm_ric_interval2_butter <- 
  d_ric_butter |> 
  select(draw, year_start = two_A) |> 
  left_join(d_ric_butter,
            join_by(draw), 
            relationship = "many-to-many") |> 
  filter(two_A < year_start + length_interval2 &  # means length_interval/2 points per regression
           two_A >= year_start) |> 
  mutate(year_mean = year_start + (length_interval2 - 1)/2) |> 
  group_by(draw, year_mean) %>% 
  filter(n() == floor(length_interval2 / 2)) |> # exclude cases with not enough observations (at the end)
  summarise(f_lmcoefs(sric_rel ~ two_A),
            .groups = "drop") |> 
  group_by(year_mean) |> 
  summarise(mean  = mean(lm_coef),
            lower = ci(lm_coef)$CI_low,
            upper = ci(lm_coef)$CI_high,
            .groups = "drop")

# Zonal, 8 years
d_lm_ric_interval_reg_butter <-
  d_ric_reg_butter |> 
  select(draw, zone, year_start = two_A) |> 
  left_join(d_ric_reg_butter,
            join_by(draw, zone), 
            relationship = "many-to-many") |> 
  filter(two_A < year_start + length_interval &  # means length_interval/2 points per regression
           two_A >= year_start) |> 
  mutate(year_mean = year_start + (length_interval - 1)/2) |> 
  group_by(draw, zone, year_mean) %>% 
  filter(n() == floor(length_interval / 2)) |> # exclude cases with not enough observations (at the end)
  summarise(f_lmcoefs(sric_rel ~ two_A),
            .groups = "drop") |>  
  group_by(zone, year_mean) |> 
  summarise(mean  = mean(lm_coef),
            lower = ci(lm_coef)$CI_low,
            upper = ci(lm_coef)$CI_high,
            .groups = "drop")

# Zonal, 12 years
d_lm_ric_interval2_reg_butter <-
  d_ric_reg_butter |> 
  select(draw, zone, year_start = two_A) |> 
  left_join(d_ric_reg_butter,
            join_by(draw, zone), 
            relationship = "many-to-many") |> 
  filter(two_A < year_start + length_interval2 &  # means length_interval/2 points per regression
           two_A >= year_start) |> 
  mutate(year_mean = year_start + (length_interval2 - 1)/2) |> 
  group_by(draw, zone, year_mean) %>% 
  filter(n() == floor(length_interval2 / 2)) |> # exclude cases with not enough observations (at the end)
  summarise(f_lmcoefs(sric_rel ~ two_A),
            .groups = "drop") |>  
  group_by(zone, year_mean) |> 
  summarise(mean  = mean(lm_coef),
            lower = ci(lm_coef)$CI_low,
            upper = ci(lm_coef)$CI_high,
            .groups = "drop")

# for all the different traitvalues
d_lm_ric_interval_trait_reg_butter <- data.frame()
d_lm_ric_interval2_trait_reg_butter <- data.frame()
for (trait_i in sel_traits_butter){
  d_ric_i <- d_ric_trait_butter |> 
    filter(trait == trait_i)
  
  d_ric_reg_i <- d_ric_trait_reg_butter |> 
    filter(trait == trait_i)
  
  d_ric_reg_i <- d_ric_reg_i |> 
    filter(traitvalue != "" & !is.na(traitvalue)) |> 
    left_join(d_ric_i |> 
                filter(two_A == 1930) |> 
                group_by(traitvalue) |> 
                summarise(sric_baseline = mean(sric)),
              join_by(traitvalue)) |> 
    mutate(sric_rel = sric / sric_baseline)
  
  d_lm_ric_interval_trait_reg_butter <- d_ric_reg_i |> 
    select(traitvalue, draw, zone, year_start = two_A) |> 
    left_join(d_ric_reg_i,
              join_by(draw, zone, traitvalue), 
              relationship = "many-to-many") |> 
    filter(two_A < year_start + length_interval &  # means length_interval/2 points per regression
             two_A >= year_start) |> 
    mutate(year_mean = year_start + (length_interval - 1)/2) |> 
    group_by(trait, traitvalue, draw, zone, year_mean) %>% 
    filter(n() == floor(length_interval / 2)) |> # exclude cases with not enough observations (at the end)
    summarise(f_lmcoefs(sric_rel ~ two_A),
              .groups = "drop") |>  
    group_by(trait, traitvalue, zone, year_mean) |> 
    summarise(mean  = mean(lm_coef),
              lower = ci(lm_coef)$CI_low,
              upper = ci(lm_coef)$CI_high,
              .groups = "drop") |> 
    (\(x) bind_rows(d_lm_ric_interval_trait_reg_butter, x))()
  
  d_lm_ric_interval2_trait_reg_butter <- d_ric_reg_i |> 
    select(traitvalue, draw, zone, year_start = two_A) |> 
    left_join(d_ric_reg_i,
              join_by(draw, zone, traitvalue), 
              relationship = "many-to-many") |> 
    filter(two_A < year_start + length_interval2 &  # means length_interval/2 points per regression
             two_A >= year_start) |> 
    mutate(year_mean = year_start + (length_interval2 - 1)/2) |> 
    group_by(trait, traitvalue, draw, zone, year_mean) %>% 
    filter(n() == floor(length_interval2 / 2)) |> # exclude cases with not enough observations (at the end)
    summarise(f_lmcoefs(sric_rel ~ two_A),
              .groups = "drop") |>  
    group_by(trait, traitvalue, zone, year_mean) |> 
    summarise(mean  = mean(lm_coef),
              lower = ci(lm_coef)$CI_low,
              upper = ci(lm_coef)$CI_high,
              .groups = "drop") |> 
    (\(x) bind_rows(d_lm_ric_interval2_trait_reg_butter, x))()
  
}

# save outputs for later processing
saveRDS(d_lm_ric_interval_butter, "Data/tmp/sric_agg/d_lm_ric_interval_butter.rds")
saveRDS(d_lm_ric_interval2_butter, "Data/tmp/sric_agg/d_lm_ric_interval2_butter.rds")
saveRDS(d_lm_ric_interval_reg_butter, "Data/tmp/sric_agg/d_lm_ric_interval_reg_butter.rds")
saveRDS(d_lm_ric_interval2_reg_butter, "Data/tmp/sric_agg/d_lm_ric_interval2_reg_butter.rds")
saveRDS(d_lm_ric_interval_trait_reg_butter, "Data/tmp/sric_agg/d_lm_ric_interval_trait_reg_butter.rds")
saveRDS(d_lm_ric_interval2_trait_reg_butter, "Data/tmp/sric_agg/d_lm_ric_interval2_trait_reg_butter.rds")

# load processed data ##########################################################
################################################################################.

# richness trends saproxylic beetles -------------------------------------------.

d_lm_ric_interval_sapro <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval_sapro.rds")
d_lm_ric_interval2_sapro <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval2_sapro.rds")
d_lm_ric_interval_reg_sapro <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval_reg_sapro.rds")
d_lm_ric_interval2_reg_sapro <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval2_reg_sapro.rds")
d_lm_ric_interval_trait_reg_sapro <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval_trait_reg_sapro.rds")
d_lm_ric_interval2_trait_reg_sapro <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval2_trait_reg_sapro.rds")


# richness trends butterflies --------------------------------------------------.

d_lm_ric_interval_butter <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval_butter.rds")
d_lm_ric_interval2_butter <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval2_butter.rds")
d_lm_ric_interval_reg_butter <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval_reg_butter.rds")
d_lm_ric_interval2_reg_butter <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval2_reg_butter.rds")
d_lm_ric_interval_trait_reg_butter <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval_trait_reg_butter.rds")
d_lm_ric_interval2_trait_reg_butter <- readRDS("Data/tmp/sric_agg/d_lm_ric_interval2_trait_reg_butter.rds")


################################################################################.
# DRIVER MODELS ---------- #####################################################
################################################################################.

# in this section, the models relating richness trends and environmental variables are fitted

# select one of the two interval lengths (length of interval in analyses)

sel_interval <- length_interval # main analyses
# sel_interval <- length_interval2 # sensitivity analyses

# prepare data for modelling ###################################################
################################################################################.

# scale all numeric variables
d_drivers <- d_drivers_raw |> 
  group_by(length_interval) |> 
  mutate(across(where(is.numeric) & !c(year_mean, length_interval) | 
                  mechanisation_abs | human_pop_change, 
                ~ (. - mean(.)) / (2 * sd(.)))) |> 
  ungroup()

# saproxylic beetles -----------------------------------------------------------.

# define model data frames
# join trend data and environmental variables
# one dataset not including traits, another including traits
if (sel_interval == 8){
  d_mod_sapro <- d_drivers |> 
    filter(length_interval == sel_interval) |> 
    select(-length_interval) |> 
    left_join(d_lm_ric_interval_reg_sapro, by = c("zone", "year_mean"),
              relationship = "one-to-one")
  
  
  d_mod_sapro_trait <- d_drivers |> 
    filter(length_interval == sel_interval) |> 
    select(-length_interval) |> 
    left_join(d_lm_ric_interval_trait_reg_sapro, by = c("zone", "year_mean"),
              relationship = "one-to-many")
} else if (sel_interval == 12){
  d_mod_sapro <- d_drivers |> 
    filter(length_interval == sel_interval) |> 
    select(-length_interval) |> 
    left_join(d_lm_ric_interval2_reg_sapro, by = c("zone", "year_mean"),
              relationship = "one-to-one")
  
  d_mod_sapro_trait <- d_drivers |> 
    filter(length_interval == sel_interval) |> 
    select(-length_interval) |> 
    left_join(d_lm_ric_interval2_trait_reg_sapro, by = c("zone", "year_mean"),
              relationship = "one-to-many")
}

# butterflies ------------------------------------------------------------------.

# define model data frames
# join trend data and environmental variables
# one dataset not including traits, another including traits
if (sel_interval == 8){
  d_mod_butter <- d_drivers |> 
    filter(length_interval == sel_interval) |> 
    select(-length_interval)  |> 
    left_join(d_lm_ric_interval_reg_butter, by = c("zone", "year_mean"),
              relationship = "one-to-one")
  
  d_mod_butter_trait <- d_drivers |> 
    filter(length_interval == sel_interval) |> 
    select(-length_interval) |> 
    left_join(d_lm_ric_interval_trait_reg_butter, 
              by = c("zone", "year_mean"),
              relationship = "one-to-many") |> 
    mutate(traitvalue = ifelse(traitvalue == "", NA, traitvalue))
} else if (sel_interval == 12){
  d_mod_butter <- d_drivers |> 
    filter(length_interval == sel_interval) |> 
    select(-length_interval)  |> 
    left_join(d_lm_ric_interval2_reg_butter, by = c("zone", "year_mean"),
              relationship = "one-to-one")
  
  d_mod_butter_trait <- d_drivers |> 
    filter(length_interval == sel_interval) |> 
    select(-length_interval) |> 
    left_join(d_lm_ric_interval_trait_reg_butter, 
              by = c("zone", "year_mean"),
              relationship = "one-to-many") |> 
    mutate(traitvalue = ifelse(traitvalue == "", NA, traitvalue))
}

# check correlations ###########################################################
################################################################################.

# check correlations between environmental variables

d_cor_8yrs <- d_drivers |> 
  filter(length_interval == 8) |>
  select(forest_area_change:temperature_abs) |> 
  select(-mechanisation_period, -after_storm) |> 
  cor() |> 
  as.data.frame() %>% 
  rownames_to_column("var1") %>% 
  pivot_longer(-var1, names_to = "var2", values_to = "cor") %>% 
  mutate(cor = ifelse(var1 != var2, cor, NA),
         var1 = factor(var1, levels = unique(var1)),
         var2 = factor(var2, levels = rev(levels(var1)))) %>% 
  mutate(cor_str = formatC(cor, digits = 2, format = "f"))


d_cor_8yrs %>% 
  ggplot(aes(x = var1, y = var2, fill = cor)) +
  geom_tile() +
  geom_text(data = \(x) filter(x, !is.na(cor)), aes(label = cor_str)) +
  scale_fill_gradient2(limits = c(-1, 1), name = "Correlation",
                       na.value = "grey70") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  coord_fixed() 

d_cor_12yrs <- d_drivers |> 
  filter(length_interval == 12) |>
  select(forest_area_change:temperature_abs) |> 
  select(-mechanisation_period, -after_storm) |> 
  cor() |> 
  as.data.frame() %>% 
  rownames_to_column("var1") %>% 
  pivot_longer(-var1, names_to = "var2", values_to = "cor") %>% 
  mutate(cor = ifelse(var1 != var2, cor, NA),
         var1 = factor(var1, levels = unique(var1)),
         var2 = factor(var2, levels = rev(levels(var1)))) %>% 
  mutate(cor_str = formatC(cor, digits = 2, format = "f"))


d_cor_12yrs %>% 
  ggplot(aes(x = var1, y = var2, fill = cor)) +
  geom_tile() +
  geom_text(data = \(x) filter(x, !is.na(cor)), aes(label = cor_str)) +
  scale_fill_gradient2(limits = c(-1, 1), name = "Correlation",
                       na.value = "grey70") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  coord_fixed() 

# overall models ###############################################################
################################################################################.

# fit the overarching models (not including trait data) per insect group

# saproxylic beetles -----------------------------------------------------------.

set.seed(143)
mod_main_sapro <- d_mod_sapro |>
  (\(x) brm(mean ~
              temperature_abs + temperature_change  +
              human_pop_change +
              mechanisation_abs + mechanisation_period +
              grassland_area_abs + grassland_area_change +
              forest_area_abs + forest_area_change  +
              wood_harvest_int_change + after_storm +
              (1 | year_mean) + (1 | zone),
            prior = prior(normal(0,5), class = b) +
              prior(normal(0,5), class = Intercept) +
              prior(cauchy(0,1), class = sd) +
              prior(cauchy(0,25), class = sigma),
            data = x,
            iter = n_iter,
            family = student(),
            file = paste0("Models/mod_sapro_",
                          sel_interval, "yrs_",
                          n_iter, "iters.rds")))()


# Model checks
# pp_check(mod_main_sapro)
# model.check <- createDHARMa(
#   simulatedResponse = t(posterior_predict(mod_main_sapro)),
#   observedResponse = d_mod_sapro$mean,
#   fittedPredictedResponse = apply(t(posterior_epred(mod_main_sapro)), 1, mean),
#   integerResponse = FALSE)
# 
# plot(model.check)
# plot(model.check, form = d_mod_sapro$temperature_abs)
# plot(model.check, form = d_mod_sapro$temperature_change)
# plot(model.check, form = d_mod_sapro$human_pop_change)
# plot(model.check, form = d_mod_sapro$forest_area_abs)
# plot(model.check, form = d_mod_sapro$forest_area_change)
# plot(model.check, form = d_mod_sapro$wood_harvest_int_change)
# plot(model.check, form = d_mod_sapro$after_storm)

# summarise model outputs (for later plotting)
d_smry_sapro <- as.data.frame(mod_main_sapro) |> 
  select(contains("b_"),
         -b_Intercept) |> 
  rename_all(~ gsub("b_|TRUE", "", .)) |> 
  pivot_longer(everything(), names_to = "var", values_to = "value") |> 
  group_by(var) |> 
  summarise(ci(value, .95) |> 
              rename(lower95 = CI_low,
                     upper95 = CI_high),
            ci(value, .9) |> 
              rename(lower90 = CI_low,
                     upper90 = CI_high),
            ci(value, .8) |> 
              rename(lower80 = CI_low,
                     upper80 = CI_high),
            mean = mean(value),
            .groups = "drop") |> 
  select(-CI)


# butterflies ------------------------------------------------------------------.

set.seed(1211)
mod_main_butter <- d_mod_butter |>
  (\(x) brm(mean ~ 
              temperature_abs + temperature_change  +
              human_pop_change +
              mechanisation_abs + mechanisation_period +
              grassland_area_abs + grassland_area_change +
              forest_area_abs + forest_area_change  + 
              wood_harvest_int_change + after_storm +
              (1 | year_mean) + (1 | zone),
            prior = prior(normal(0,5), class = b) +
              prior(normal(0,5), class = Intercept) +
              prior(cauchy(0,1), class = sd) +
              prior(cauchy(0,25), class = sigma),
            data = x,
            iter = n_iter,
            family = student(),
            file = paste0("Models/mod_butter_", 
                          sel_interval, "yrs_", 
                          n_iter, "iters.rds")))()

# Model checks
# pp_check(mod_main_butter)
# model.check <- createDHARMa(
#   simulatedResponse = t(posterior_predict(mod_main_butter)),
#   observedResponse = d_mod_butter$mean,
#   fittedPredictedResponse = apply(t(posterior_epred(mod_main_butter)), 1, mean),
#   integerResponse = FALSE)
# 
# plot(model.check)
# plot(model.check, form = d_mod_butter$temperature_abs)
# plot(model.check, form = d_mod_butter$temperature_change)
# plot(model.check, form = d_mod_butter$human_pop_change)
# plot(model.check, form = d_mod_butter$mechanisation_abs)
# plot(model.check, form = d_mod_butter$mechanisation_period)
# plot(model.check, form = d_mod_butter$grassland_area_abs)
# plot(model.check, form = d_mod_butter$grassland_area_change)

# summarise model outputs (for later plotting)
d_smry_butter <- as.data.frame(mod_main_butter) |> 
  select(contains("b_"),
         -b_Intercept) |> 
  rename_all(~ gsub("b_|TRUE", "", .)) |> 
  pivot_longer(everything(), names_to = "var", values_to = "value") |> 
  group_by(var) |> 
  summarise(ci(value, .95) |> 
              rename(lower95 = CI_low,
                     upper95 = CI_high),
            ci(value, .9) |> 
              rename(lower90 = CI_low,
                     upper90 = CI_high),
            ci(value, .8) |> 
              rename(lower80 = CI_low,
                     upper80 = CI_high),
            mean = mean(value),
            .groups = "drop") |> 
  select(-CI)

# trait models (incl. plotting) ################################################
################################################################################.

# fit the models per insect group and trait, including the interaction of environmental variables and traitvalues
# then, plot the outputs (Extended data figures)

# saproxylic beetles -----------------------------------------------------------.

# run models and predictions
l_traitmod_sapro <- lapply(sel_traits_sapro,
                           f_traitmod, group_i = "sapro")
names(l_traitmod_sapro) <- sel_traits_sapro

# rearrange data for plotting
l_comb_plotsdata <- lapply(sel_traits_sapro, f_comb_plotsdata,
                           l_traitmod_sapro)

l_plotsdata_sapro <- list()
# join all the range data (i.e. data for vertical CI ranges and means)
l_plotsdata_sapro$rangedata <- lapply(l_comb_plotsdata, function(x) x$rangedata) |> 
  bind_rows() |> 
  mutate(across(lower95:mean, ~ . * sel_interval * 100)) # change to percentage point change per interval
# join all the line data (i.e. data for lines connecting minimum- and maximum-level predictions)
l_plotsdata_sapro$linedata <- lapply(l_comb_plotsdata, function(x) x$linedata) |> 
  bind_rows() |> 
  mutate(across(pred_mean:pred_lower, ~ . * sel_interval * 100)) # change to percentage point change per interval

# rearrange data along x axis and other smaller modifications
l_plotsdata_sapro <- lapply(l_plotsdata_sapro,
                            function(x) x |> 
                              group_by(trait, var) |> 
                              mutate(x = ifelse(var != "overall",
                                                x + (traitvalue_num - median(traitvalue_num)) * diff_x * 0.15, # position dodge
                                                x)) |> 
                              ungroup() |> 
                              rowwise() |> 
                              mutate(trait = names(sel_traits_sapro)[sel_traits_sapro == trait],
                                     var = ifelse(var != "overall", 
                                                  names(v_vars_spec_short)[v_vars_spec_short == var],
                                                  var)) |> 
                              ungroup() |> 
                              mutate(trait = factor(trait,
                                                    levels = unique(names(sel_traits_sapro))),
                                     var = factor(var, levels = c("overall",
                                                                  names(v_vars_spec_short)))))

# get x axis ranges
l_xscales_raw <- l_plotsdata_sapro$rangedata |> 
  group_by(var) |> 
  summarise(xmin = min(x) - 0.15 * unique(diff_x),
            xmax = max(x) + 0.15 * unique(diff_x),
            .groups = "drop") |> 
  (\(x) split(x, x$var))()

# convert to x axis argument
l_xscales <- l_xscales_raw |> 
  map(function(x) scale_x_continuous(limits = c(x$xmin, x$xmax)))

# different x scales for factors
l_xscales[c("Mechanis. period", "Storm aftermath")] <- 
  l_xscales_raw[c("Mechanis. period", "Storm aftermath")] |> 
  map(function(x) scale_x_continuous(limits = c(x$xmin, x$xmax), breaks = c(0, 1),
                                     labels = c("no", "yes")))

# different x axis scales for overall plot
l_xscales[c("overall")] <- 
  l_xscales_raw[c("overall")] |> 
  map(function(x) scale_x_continuous(limits = c(x$xmin, x$xmax), breaks = 999))


# strip texts for environmental variables
v_labeller <- c("", gsub(" ", "\n", names(v_vars_spec_short)), unique(names(sel_traits_sapro)))
names(v_labeller) <- c("overall", names(v_vars_spec_short), unique(names(sel_traits_sapro)))


pgrid <-
  l_plotsdata_sapro$rangedata |> 
  ggplot(aes(x = x, xend = x,
             colour = as.factor(traitvalue_num))) +
  # add different grey backgrounds for different trais
  geom_rect(data = function(x) x |> 
              select(trait, var) |> 
              distinct() |> 
              filter(trait %in% unique(trait)[seq(1, 9, 2)]),
            inherit.aes = F,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_rect(data = function(x) x |> 
              select(trait, var) |> 
              distinct() |> 
              filter(trait %in% unique(trait)[seq(2, 10, 2)]),
            inherit.aes = F,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey75") +
  # add border for overall plot
  geom_rect(data = function(x) x |> 
              select(trait, var) |> 
              distinct() |> 
              filter(var == "overall"),
            inherit.aes = F,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, colour = 1, fill = NA) +
  # add zero line
  geom_hline(yintercept = 0, colour = "grey20") +
  # add lines connecting minimum- and maximum-level predictions
  geom_line(data = l_plotsdata_sapro$linedata,
            aes(y = pred_mean, linetype = type), size = .5) +
  # 95% CI (vertical line)
  geom_segment(aes(y = lower95, yend = upper95),
               linewidth = .25) +
  # 90% CI (vertical line)
  geom_segment(aes(y = lower90, yend = upper90),
               linewidth = .75) +
  # 80% CI (vertical line)
  geom_segment(aes(y = lower80, yend = upper80),
               linewidth = 1.5) +
  # add mean prediction
  geom_point(aes(y = mean, fill = as.factor(traitvalue_num)),
             size = 2, pch = 21, colour = "grey35") +
  # plot adjustments
  scale_linetype_manual(values = c(cont. = 1, factor = 11)) +
  scale_y_continuous(guide = guide_axis(position = "right")) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  # facetting
  facet_grid(trait ~ var, scales = "free", switch = "both",
             labeller = as_labeller(v_labeller)) +
  # add the different x axes to the different facets
  facetted_pos_scales(x = l_xscales)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text.y = element_text(size = v_textsize["axis.title"]),
        strip.text.x = element_text(size = v_textsize["axis.title"],
                                    angle = 90, vjust = .5, hjust = 1),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        strip.placement = "outside",
        strip.background.x = element_blank()) +
  labs(y = paste0("Species richness trend\n(percentage point change per ", sel_interval, " years)"))


# adjust width of first column
gt <- ggplot_gtable(ggplot_build(pgrid))
gt$widths[8] <- .7 * gt$widths[8]

# export plot
cairo_pdf(paste0("Output/Figures/lm_ric_mod_trait_sapro_", sel_interval, "yrs.pdf"),
          width = 180 / 25.4, height = 170 / 25.4, fallback_resolution = 400) #170mm is max height
grid.draw(gt)

# add legend insets
for (i in seq_len(n_distinct(l_plotsdata_sapro$rangedata$trait))){
  trait_i <- unique((l_plotsdata_sapro$rangedata$trait))[i]
  
  fill_i <- ifelse(i %% 2 == 0, "grey75", "grey90")
  
  p_inset <-
    d_traitsel_comb |> 
    filter(trait == trait_i) |> 
    mutate(traitvalue_short = factor(traitvalue_short, levels = unique(traitvalue_short))) |> 
    ggplot(aes(colour = traitvalue_short, x = 1, y = 1)) +
    geom_line(size = 1, group = 1) +
    theme_void() +
    scale_colour_manual(values = RColorBrewer::brewer.pal(3, "Dark2"), name = "") +
    theme(legend.position = "inside",
          legend.position.inside = c(0, 1),
          legend.justification.inside = c(0, .67),
          legend.key.height = unit(.2, "cm"),
          legend.key.width = unit(.4, "cm"),
          legend.key.spacing = unit(.05, "cm"),
          plot.background = element_rect(fill = fill_i, colour = 1),
          legend.text = element_text(size = v_textsize["legend.text"]))
  
  
  inset_vp <- viewport(x = 0.085, y = .972 - (i - 1)/4.72, width = 0.092, height = 0.05)
  pushViewport(inset_vp)
  inset_grob <- ggplotGrob(p_inset)
  grid.draw(inset_grob)
  popViewport()
}

dev.off()

# butterflies ------------------------------------------------------------------.

# run models and predictions
l_traitmod_butter <- lapply(sel_traits_butter,
                            f_traitmod, group_i = "butter")
names(l_traitmod_butter) <- sel_traits_butter

# rearrange data for plotting
l_comb_plotsdata <- lapply(sel_traits_butter, f_comb_plotsdata,
                           l_traitmod_butter)

l_plotsdata_butter <- list()
# join all the range data (i.e. data for vertical CI ranges and means)
l_plotsdata_butter$rangedata <- lapply(l_comb_plotsdata, function(x) x$rangedata) |> 
  bind_rows() |> 
  mutate(across(lower95:mean, ~ . * sel_interval * 100)) # change to percentage point change per interval
# join all the line data (i.e. data for lines connecting minimum- and maximum-level predictions)
l_plotsdata_butter$linedata <- lapply(l_comb_plotsdata, function(x) x$linedata) |> 
  bind_rows() |> 
  mutate(across(pred_mean:pred_lower, ~ . * sel_interval * 100)) # change to percentage point change per interval

# rearrange data along x axis and other smaller modifications
l_plotsdata_butter <- lapply(l_plotsdata_butter,
                             function(x) x |> 
                               filter(trait %in% sel_traits_butter[1:4]) |> 
                               group_by(trait, var) |> 
                               mutate(x = ifelse(var != "overall",
                                                 x + (traitvalue_num - median(traitvalue_num)) * diff_x * 0.15, # position dodge
                                                 x)) |> 
                               ungroup() |> 
                               rowwise() |> 
                               mutate(trait = names(sel_traits_butter)[sel_traits_butter == trait],
                                      var = ifelse(var != "overall", 
                                                   names(v_vars_spec_short)[v_vars_spec_short == var],
                                                   var)) |> 
                               ungroup() |> 
                               mutate(trait = factor(trait,
                                                     levels = unique(names(sel_traits_butter[1:4]))),
                                      var = factor(var, levels = c("overall",
                                                                   names(v_vars_spec_short)))))

# get x axis ranges
l_xscales_raw <- l_plotsdata_butter$rangedata |> 
  group_by(var) |> 
  summarise(xmin = min(x) - 0.15 * unique(diff_x),
            xmax = max(x) + 0.15 * unique(diff_x),
            .groups = "drop") |> 
  (\(x) split(x, x$var))()

# convert to x axis argument
l_xscales <- l_xscales_raw |> 
  map(function(x) scale_x_continuous(limits = c(x$xmin, x$xmax)))

# different x scales for factors
l_xscales[c("Mechanis. period", "Storm aftermath")] <- 
  l_xscales_raw[c("Mechanis. period", "Storm aftermath")] |> 
  map(function(x) scale_x_continuous(limits = c(x$xmin, x$xmax), breaks = c(0, 1),
                                     labels = c("no", "yes")))

# different x axis scales for overall plot
l_xscales[c("overall")] <- 
  l_xscales_raw[c("overall")] |> 
  map(function(x) scale_x_continuous(limits = c(x$xmin, x$xmax), breaks = 999))


# strip texts for environmental variables
v_labeller <- c("", gsub(" ", "\n", names(v_vars_spec_short)), unique(names(sel_traits_butter[1:4])))
names(v_labeller) <- c("overall", names(v_vars_spec_short), unique(names(sel_traits_butter[1:4])))


pgrid <-
  l_plotsdata_butter$rangedata |> 
  ggplot(aes(x = x, xend = x,
             colour = as.factor(traitvalue_num))) +
  # add different grey backgrounds for different trais
  geom_rect(data = function(x) x |> 
              select(trait, var) |> 
              distinct() |> 
              filter(trait %in% unique(trait)[seq(1, 9, 2)]),
            inherit.aes = F,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_rect(data = function(x) x |> 
              select(trait, var) |> 
              distinct() |> 
              filter(trait %in% unique(trait)[seq(2, 10, 2)]),
            inherit.aes = F,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey75") +
  # add border for overall plot
  geom_rect(data = function(x) x |> 
              select(trait, var) |> 
              distinct() |> 
              filter(var == "overall"),
            inherit.aes = F,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, colour = 1, fill = NA) +
  # add zero line
  geom_hline(yintercept = 0, colour = "grey20") +
  # add lines connecting minimum- and maximum-level predictions
  geom_line(data = l_plotsdata_butter$linedata,
            aes(y = pred_mean, linetype = type), size = .5) +
  # 95% CI (vertical line)
  geom_segment(aes(y = lower95, yend = upper95),
               linewidth = .25) +
  # 90% CI (vertical line)
  geom_segment(aes(y = lower90, yend = upper90),
               linewidth = .75) +
  # 80% CI (vertical line)
  geom_segment(aes(y = lower80, yend = upper80),
               linewidth = 1.5) +
  # add mean prediction
  geom_point(aes(y = mean, fill = as.factor(traitvalue_num)),
             size = 2, pch = 21, colour = "grey35") +
  # plot adjustments
  scale_linetype_manual(values = c(cont. = 1, factor = 11)) +
  scale_y_continuous(guide = guide_axis(position = "right")) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  # facetting
  facet_grid(trait ~ var, scales = "free", switch = "both",
             labeller = as_labeller(v_labeller)) +
  # add the different x axes to the different facets
  facetted_pos_scales(x = l_xscales)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text.y = element_text(size = v_textsize["axis.title"]),
        strip.text.x = element_text(size = v_textsize["axis.title"],
                                    angle = 90, vjust = .5, hjust = 1),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        strip.placement = "outside",
        strip.background.x = element_blank()) +
  labs(y = paste0("Species richness trend\n(percentage point change per ", sel_interval, " years)"))


# adjust width of first column
gt <- ggplot_gtable(ggplot_build(pgrid))
gt$widths[8] <- .7 * gt$widths[8]

# export plot
cairo_pdf(paste0("Output/Figures/lm_ric_mod_trait_butter_", sel_interval, "yrs.pdf"),
          width = 180 / 25.4, height = 170 / 25.4, fallback_resolution = 400) #170mm is max height
grid.draw(gt)

# add legend insets
for (i in seq_len(n_distinct(l_plotsdata_butter$rangedata$trait))){
  trait_i <- unique((l_plotsdata_butter$rangedata$trait))[i]
  
  fill_i <- ifelse(i %% 2 == 0, "grey75", "grey90")
  
  p_inset <-
    d_traitsel_comb |> 
    filter(trait == trait_i) |> 
    mutate(traitvalue_short = factor(traitvalue_short, levels = unique(traitvalue_short))) |> 
    ggplot(aes(colour = traitvalue_short, x = 1, y = 1)) +
    geom_line(size = 1, group = 1) +
    theme_void() +
    scale_colour_manual(values = RColorBrewer::brewer.pal(3, "Dark2"), name = "") +
    theme(legend.position = "inside",
          legend.position.inside = c(0, 1),
          legend.justification.inside = c(0, .67),
          legend.key.height = unit(.2, "cm"),
          legend.key.width = unit(.4, "cm"),
          legend.key.spacing = unit(.05, "cm"),
          plot.background = element_rect(fill = fill_i, colour = 1),
          legend.text = element_text(size = v_textsize["legend.text"]))
  
  
  inset_vp <- viewport(x = 0.087, y = .972 - (i - 1)/4.72, width = 0.095, height = 0.05)
  pushViewport(inset_vp)
  inset_grob <- ggplotGrob(p_inset)
  grid.draw(inset_grob)
  popViewport()
}

dev.off()

################################################################################.
# TEXT OUTPUT -------------- ###################################################
################################################################################.

# in this section, calculate numbers which are included in the main text of the manuscript

# richness change ##############################################################
################################################################################.

# calculate the change in richness (percentage points, proportional) for different insect groups and time periods

# saproxylic beetles, per decade, whole of Switzerland
f_sric_change(d_ric_sapro,
              length_interval_i = 8)|> 
  filter(year_start %in% seq(1930, 2010, 10))

# butterflies, per decade, whole of Switzerland
f_sric_change(d_ric_butter,
              length_interval_i = 8)|> 
  filter(year_start %in% seq(1930, 2010, 10))

# saproxylic beetles, whole study period, whole of Switzerland
f_sric_change(d_ric_sapro,
              length_interval_i = 90)

# butterflies, whole study period, whole of Switzerland
f_sric_change(d_ric_butter,
              length_interval_i = 90)

# butterflies, whole study period, per zone
d_ric_reg_butter |> 
  group_by(zone) |> 
  group_map(~ tibble(zone = unique(.$zone),
                         f_sric_change(., length_interval_i = 90)),
            .keep = T) |> 
  bind_rows()


# saproxylic beetles, last 4 decades only
f_sric_change(d_ric_sapro,
              length_interval_i = 40) |> 
  filter(year_end == 2020)

# butterflies, last 4 decades only
f_sric_change(d_ric_butter,
              length_interval_i = 40) |> 
  filter(year_end == 2020)

# butterflies, last 2 decades only
f_sric_change(d_ric_butter,
              length_interval_i = 20) |> 
  filter(year_end == 2020)

# Traits -----------------------------------------------------------------------.

# richness changes for different traits / traitvalues

# body size
d_ric_trait_sapro |> 
  filter(trait == "size_cat",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  group_map(~ tibble(traitvalue = unique(.$traitvalue),
                     f_sric_change(., length_interval_i = 90)),
            .keep = T) |> 
  bind_rows()

d_ric_trait_butter |> 
  filter(trait == "size_cat",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  group_map(~ tibble(traitvalue = unique(.$traitvalue),
                     f_sric_change(., length_interval_i = 90)),
            .keep = T) |> 
  bind_rows()

# habitat specialisation
d_ric_trait_sapro |> 
  filter(trait == "stenotopy",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  group_map(~ tibble(traitvalue = unique(.$traitvalue),
                     f_sric_change(., length_interval_i = 90)),
            .keep = T) |> 
  bind_rows()

d_ric_trait_butter |> 
  filter(trait == "stenotopy",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  group_map(~ tibble(traitvalue = unique(.$traitvalue),
                     f_sric_change(., length_interval_i = 90)),
            .keep = T) |> 
  bind_rows()

# food specialisation
d_ric_trait_sapro |> 
  filter(trait == "food_spec",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  group_map(~ tibble(traitvalue = unique(.$traitvalue),
                     f_sric_change(., length_interval_i = 90)),
            .keep = T) |> 
  bind_rows()

d_ric_trait_butter |> 
  filter(trait == "food_spec",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  group_map(~ tibble(traitvalue = unique(.$traitvalue),
                     f_sric_change(., length_interval_i = 90)),
            .keep = T) |> 
  bind_rows()

# temperature niche
d_ric_trait_sapro |> 
  filter(trait == "Tniche",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  group_map(~ tibble(traitvalue = unique(.$traitvalue),
                     f_sric_change(., length_interval_i = 50)),
            .keep = T) |> 
  bind_rows() |> 
  filter(year_start == 1930)

d_ric_trait_butter |> 
  filter(trait == "Tniche",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  group_map(~ tibble(traitvalue = unique(.$traitvalue),
                     f_sric_change(., length_interval_i = 50)),
            .keep = T) |> 
  bind_rows() |> 
  filter(year_start == 1930)

d_ric_trait_sapro |> 
  filter(trait == "Tniche",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  group_map(~ tibble(traitvalue = unique(.$traitvalue),
                     f_sric_change(., length_interval_i = 90)),
            .keep = T) |> 
  bind_rows()

d_ric_trait_butter |> 
  filter(trait == "Tniche",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  group_map(~ tibble(traitvalue = unique(.$traitvalue),
                     f_sric_change(., length_interval_i = 90)),
            .keep = T) |> 
  bind_rows()


# change since 1980 for different temperature niche groups
d_ric_trait_mean_sapro |> 
  filter(trait == "Tniche",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  mutate(across(c(mean, lower, upper), ~ (1 - . / mean[two_A == 1930])) * 100) |> 
  filter(two_A %in% c(1980, max(two_A)))

d_ric_trait_mean_butter |> 
  filter(trait == "Tniche",
         traitvalue != "") |> 
  group_by(traitvalue) |> 
  mutate(across(c(mean, lower, upper), ~ (1 - . / mean[two_A == 1930])) * 100) |> 
  filter(two_A %in% c(1980, max(two_A)))

# driver models ################################################################
################################################################################.

# quantify the effect sizes from the regression models relating richness trends and environmental variabels

f_pred_overall("sapro", "mechanisation_abs")
(f_pred_overall("sapro", "mechanisation_abs") * 8 * 100) |> signif(2) # percentage point change per 8-year interval

f_pred_overall("sapro", "temperature_abs")
(f_pred_overall("sapro", "temperature_abs") * 8 * 100) |> signif(2) # percentage point change per 8-year interval

f_pred_overall("sapro", "temperature_abs", 2)
(f_pred_overall("sapro", "temperature_abs", 2) * 8 * 100) |> signif(2) # percentage point change per 8-year interval

f_pred_overall("sapro", "after_storm")
(f_pred_overall("sapro", "after_storm") * 8 * 100) |> signif(2) # percentage point change per 8-year interval


f_pred_overall("butter", "temperature_abs")
(f_pred_overall("butter", "temperature_abs") * 8 * 100) |> signif(2) # percentage point change per 8-year interval

f_pred_overall("butter", "temperature_abs", 2)
(f_pred_overall("butter", "temperature_abs", 2) * 8 * 100) |> signif(2) # percentage point change per 8-year interval

f_pred_overall("butter", "mechanisation_abs")
(f_pred_overall("butter", "mechanisation_abs") * 8 * 100) |> signif(2) # percentage point change per 8-year interval

f_pred_overall("butter", "mechanisation_period")
(f_pred_overall("butter", "mechanisation_period") * 8 * 100) |> signif(2) # percentage point change per 8-year interval

f_pred_overall("butter", "after_storm")
(f_pred_overall("butter", "after_storm") * 8 * 100) |> signif(2) # percentage point change per 8-year interval


################################################################################.
# PLOTTING -------------------- ################################################
################################################################################.

# remaining plots are created in this section

# richness trajectories ########################################################
################################################################################.

# richness trajectories for both insect groups and for whole of Switzerland as well as the different biogeographic zones

d_ric_mean_sapro %>% 
  mutate(group = "Sapro") %>% 
  bind_rows(d_ric_mean_butter %>%
              mutate(group = "Butter")) %>%
  mutate(two_A = as.integer(two_A)) %>% 
  group_by(group) %>% 
  # convert to richness relative to the 1930 value
  mutate(upper = upper / mean[two_A == 1930],
         lower = lower / mean[two_A == 1930],
         mean = mean / mean[two_A == 1930]) %>%
  ungroup() %>% 
  mutate(zone = "Switzerland") |> 
  bind_rows(d_ric_reg_mean_sapro %>% 
              mutate(group = "Sapro") %>% 
              bind_rows(d_ric_reg_mean_butter %>%
                          mutate(group = "Butter")) %>%
              # convert to richness relative to the 1930 value (Switzerland-wide)
              left_join(d_ric_mean_sapro %>% 
                          filter(two_A == 1930) %>% 
                          mutate(group = "Sapro") %>% 
                          select(group, sric_ref = mean) %>% 
                          bind_rows(d_ric_mean_butter %>% 
                                      filter(two_A == 1930) %>%
                                      mutate(group = "Butter") %>%
                                      select(group, sric_ref = mean)),
                        by = "group") %>% 
              mutate(mean = mean / sric_ref,
                     upper = upper / sric_ref,
                     lower = lower / sric_ref)) |> 
  mutate(zone = factor(v_zones_v2[zone], levels = v_zones_v2),
         group = factor(group, levels = names(v_grouplabels))) |> 
  ggplot(aes(x = two_A, y = mean, col = group)) +
  # highlight the Switzerland-wide panel
  geom_rect(data = function(x) x |> filter(zone == "Overall"),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "navajowhite", colour = NA) +
  # decade-wise vertical lines
  geom_vline(xintercept = seq(1930, 2020, 10), colour = "grey80") +
  # 100% horizontal line
  geom_hline(yintercept = 1) +
  # credible intervals
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group),
              alpha = .5, col = NA) +
  # mean line
  geom_line(size = 1) +
  scale_colour_manual(values = v_groupcols, name = "Group", labels = v_grouplabels) +
  scale_fill_manual(values = v_groupcols, name = "Group", labels = v_grouplabels)  +
  xlab("Year") +
  ylab("Relative richness (%)") +
  scale_y_continuous(labels = \(x) x * 100) + # convert to percentage
  facet_grid(. ~ zone) +
  coord_cartesian(xlim = c(1930, 2020)) +
  theme(legend.position = "top",
        legend.background = element_rect(fill = "grey80"),
        legend.box.margin = margin(-8, 0, -7, 0),
        legend.key.size = unit(3, "mm"),
        legend.key.spacing.x = unit(13, "mm"),
        axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        legend.title = element_text(size = v_textsize["legend.title"], 
                                    margin = margin(r = 17, 
                                                    t = 1.1, b = 1,
                                                    unit = "mm")),
        legend.text = element_text(size = v_textsize["legend.text"]),
        strip.text = element_text(size = v_textsize["axis.title"]),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Output/Figures/sric_trend.pdf", 
       width = 180, height = 55, unit = "mm", dpi = 400)

# richness trajectories per trait ##############################################
################################################################################.

# main (shared) traits ---------------------------------------------------------.

# main figure including the four traits shared between saproxylic beetles and butterflies

l_plots_sric_trait <- list()
for (trait_i in c("Size", "Habitat specialisation", "Food specialisation", # loop through traits
                  "Temperature niche")){
  
  d_sric_target <- 
    # bind and standardise beetle and butterfly data for specific trait
    d_ric_trait_mean_sapro %>% 
    filter(trait == sel_traits_sapro[trait_i]) |> 
    mutate(group = v_grouplabels["Sapro"]) |> 
    bind_rows(d_ric_trait_mean_butter %>% 
                filter(trait == sel_traits_butter[trait_i]) |> 
                mutate(group = v_grouplabels["Butter"])) |> 
    # exclude very rare or missing traitvalues
    filter(!is.na(traitvalue),
           traitvalue != "larva/pupa",
           traitvalue != "") |> 
    rowwise() |> 
    mutate(traitvalue = d_traitsel_comb$traitvalue_short[d_traitsel_comb$traitvalue == traitvalue]) |> 
    ungroup() |> 
    # preparations for plotting
    mutate(two_A = as.integer(two_A),
           traitvalue = factor(traitvalue, 
                               levels = d_traitsel_comb$traitvalue_short[d_traitsel_comb$trait == trait_i]),
           trait = trait_i,
           group = factor(group, levels = v_grouplabels))
  
  
  d_plot <- d_sric_target |> 
    # convert richness to relative richness (relative to 1930 value for that trait and traitvalue)
    group_by(group, trait, traitvalue) %>% 
    mutate(upper = upper / mean[two_A == 1930],
           lower = lower / mean[two_A == 1930],
           mean = mean / mean[two_A == 1930]) %>%
    ungroup()
  
  # range of y values (used for plotting)
  yrange <- c(min(d_plot$lower),
              max(d_plot$upper))

  # get nice breaks for the y axis range
  break_names <- get_breaks(n = 5)(x = yrange)
  break_names <- break_names[break_names > 0]
  # convert underlying breaks to values between 0 and 1 (for later plotting because of two different plots in one)
  breaks <- (break_names - yrange[1]) / diff(yrange)
  
  # prepate facetting
  strips <- strip_nested(
    text_x = list(element_text(), element_blank()),
    background_x = list(element_rect(), element_blank()),
    by_layer_x = T
  )
  
  # parameters for bar plot inlet
  width_bar <- 15 # in years (trend is shown for 90 years)
  spacing_bar <- 15 # in years (trend is shown for 90 years)
  
  l_plots_sric_trait[[trait_i]] <-
    d_plot |> 
    # convert y values to the range 0 to 1
    mutate(across(c(upper, lower, mean), ~ (. - yrange[1]) / diff(yrange))) |>
    ggplot(aes(x = two_A, y = mean, col = traitvalue)) +
    # vertical lines for decades
    geom_vline(xintercept = seq(1930, 2020, 10), colour = "grey80") +
    # 100% relative richness horizontal line
    annotate(geom = "segment", x = -Inf, xend = 2021, 
             y = (1 - yrange[1]) / diff(yrange), 
             yend = (1 - yrange[1]) / diff(yrange)) +
    # 95% CIn ribbon
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = traitvalue),
                alpha = .5, col = NA) +
    # mean line
    geom_line(size = 1) +
    # bar plot inlet at the right side of the plot
    geom_area(data = d_sric_target |>
                group_by(group, two_A) |>
                # proportion of total richness in a certain two-year interval
                mutate(mean = mean / sum(mean)) |>
                ungroup() |>
                # convert year to x axis values at the right side of the main plot
                mutate(two_A = (two_A - 1930) * width_bar/90 + 2020 + spacing_bar),
              aes(fill = traitvalue), colour = NA) +
    # x axis line which is interrupted between the two subplots
    annotate(geom = "segment", x = -Inf, xend = 2021, 
             y = -.05,  yend = -.05) +
    annotate(geom = "segment", x = 2020 + spacing_bar - 2.5, xend = Inf, 
             y = -.05,  yend = -.05) +
    scale_colour_brewer(palette = "Dark2", name = NULL)+
    scale_fill_brewer(palette = "Dark2", name = NULL)  +
    xlab("Year") +
    ylab("Rel. richness (%)") +
    # breaks for y axis which were determined above
    # second y axis for bar plot inlet
    scale_y_continuous(breaks = breaks,
                       labels = break_names * 100,
                       sec.axis = sec_axis( ~ ., breaks = c(0, 0.5, 1),
                                            name = "Ov. proportion")) +
    # x axis breaks and labels for interrupted axis
    scale_x_continuous(breaks = c(1950, 1980, 2010, 
                                  2020 + spacing_bar + width_bar /2),
                       labels = c(1950, 1980, 2010, "1930-\n2021")) +
    # to not clip off the manually added x axis line
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    # facetting
    facet_nested_wrap(~ trait + group,
                      strip = strips) +
    theme(axis.line.x = element_blank(), # manually added (interrupted)
          legend.background = element_rect(fill = "grey80"),
          strip.background = element_blank(),
          panel.spacing = unit(3, "mm"),
          legend.key.size = unit(3, "mm"),
          legend.box.margin = margin(-10, -8, -10, -7),
          axis.title = element_text(size = v_textsize["axis.title"]),
          axis.text = element_text(size = v_textsize["axis.text"]),
          legend.text = element_text(size = v_textsize["legend.text"]),
          strip.text = element_text(size = v_textsize["axis.title"]),
          axis.text.x = element_text(angle = 45, hjust = 1))
}
# only include x axis title and text in the lowest plot
l_plots_sric_trait[c(1:3)] <- lapply(l_plots_sric_trait[c(1:3)], function(x) x +
                                       theme(axis.title.x = element_blank(),
                                             axis.text.x = element_blank()))

# put together plots for different plots
plot_grid(plotlist = l_plots_sric_trait,
          ncol = 1, align = "v",
          rel_heights = c(1, 1, 1, 
                          ((170-13.4)/4+13.4) / ((170-13.4)/4)), # 1.34 cm is the height of title and axis text (x)
          labels = letters[1:4],
          label_size = v_textsize["plotlabel"]) + 
  # draw vertical separation lines between panels
  draw_line(x = c(0.4, 0.4), y = c(0.085, 0.085 + 0.17), size = .75) +
  draw_line(x = c(0.4, 0.4), y = c(0.085 + 0.232, 0.085 + 0.232 + 0.17), size = .75)  +
  draw_line(x = c(0.4, 0.4), y = c(0.085 + 0.232 * 2, 0.085 + 0.232 * 2 + 0.17), size = .75) +
  draw_line(x = c(0.4, 0.4), y = c(0.085 + 0.232 * 3, 0.085 + 0.232 * 3 + 0.17), size = .75) 

# save figure
ggsave("Output/Figures/sric_trend_trait_both.pdf", 
       width = 90, height = 170, unit = "mm", dpi = 400)

# additional traits (supplement) -----------------------------------------------.

# for the two additional butterfly traits (hibernation stage, voltinism)

l_plots_sric_trait_butter_suppl_sel <- list()
for (trait_i in c("hibernation", "voltinism")){ # loop through the two traits
  d_plot <-  d_ric_trait_mean_butter %>% 
    # restrict to looped trait, exclude rare and empty traitvalues
    filter(trait == trait_i,
           !is.na(traitvalue),
           traitvalue != "larva/pupa",
           traitvalue != "") |> 
    mutate(two_A = as.integer(two_A)) %>% 
    # convert to relative richness (relative to 1930 values)
    group_by(traitvalue) %>% 
    mutate(upper = upper / mean[two_A == 1930],
           lower = lower / mean[two_A == 1930],
           mean = mean / mean[two_A == 1930]) %>%
    ungroup() |> 
    mutate(traitvalue = factor(traitvalue, levels = d_traitsel_butter$traitvalue[d_traitsel_butter$trait == names(sel_traits_butter[sel_traits_butter == trait_i])]),
           trait = names(sel_traits_butter[sel_traits_butter == trait_i]))
  
  # y axis range
  yrange <- c(min(d_plot$lower),  
              max(d_plot$upper))
  
  
  # get nice breaks for the y axis range
  break_names <- get_breaks(n = 5)(x = yrange)
  break_names <- break_names[break_names > 0]
  # convert underlying breaks to values between 0 and 1 (for later plotting because of two different plots in one)
  breaks <- (break_names - yrange[1]) / diff(yrange)
  
  p <-
    d_plot |> 
    # convert y values to the range 0 to 1
    mutate(across(c(upper, lower, mean), ~ (. - yrange[1]) / diff(yrange))) |>
    ggplot(aes(x = two_A, y = mean, col = traitvalue)) +
    # 100% relative richness horizontal line
    annotate(geom = "segment", x = -Inf, xend = 2021, 
             y = (1 - yrange[1]) / diff(yrange), 
             yend = (1 - yrange[1]) / diff(yrange)) +
    # vertical lines for decades
    geom_vline(xintercept = seq(1930, 2020, 10), colour = "grey80") +
    # 95% CI ribbon
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = traitvalue),
                alpha = .5, col = NA) +
    # mean line
    geom_line(size = 1) +
    # bar plot inlets
    geom_area(data = d_ric_trait_mean_butter %>% 
                filter(trait == trait_i,
                       !is.na(traitvalue),
                       traitvalue != "larva/pupa",
                       traitvalue != "") |> 
                # proportion of total richness in a certain two-year interval
                group_by(two_A) |> 
                mutate(mean = mean / sum(mean)) |> 
                ungroup() |> 
                # convert year to x axis values at the right side of the main plot
                mutate(two_A = (as.integer(two_A) - 1930) / 10 + 2029,
                       traitvalue = factor(traitvalue, levels = d_traitsel_butter$traitvalue[d_traitsel_butter$trait == names(sel_traits_butter[sel_traits_butter == trait_i])]),
                       trait = names(sel_traits_butter[sel_traits_butter == trait_i])),
              aes(fill = traitvalue),
              colour = NA) +
    # x axis line which is interrupted between the two subplots
    annotate(geom = "segment", x = -Inf, xend = 2021, 
             y = -.05,  yend = -.05) +
    annotate(geom = "segment", x = 2026.5, xend = Inf, 
             y = -.05,  yend = -.05) +
    scale_colour_brewer(palette = "Dark2", name = NULL)+
    scale_fill_brewer(palette = "Dark2", name = NULL)  +
    xlab("Year") +
    ylab("Rel. richness (%)") +
    # breaks for y axis which were determined above
    # second y axis for bar plot inlet
    scale_y_continuous(breaks = breaks,
                       labels = break_names * 100,
                       sec.axis = sec_axis( ~ ., breaks = c(0, 0.5, 1),
                                            name = "Ov. proportion")) +
    # x axis breaks and labels for interrupted axis
    scale_x_continuous(breaks = c(1950, 1980, 2010, 2033.5),
                       labels = c(1950, 1980, 2010, "1930-\n2021")) +
    # to not clip off the manually added x axis line
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    facet_grid(~ trait) +
    theme(axis.line.x = element_blank(), # manually added (interrupted)
          legend.position = "bottom",
          legend.background = element_rect(fill = "grey80"),
          panel.spacing = unit(3, "mm"),
          legend.key.size = unit(3, "mm"),
          legend.box.margin = margin(-10, -8, -10, -7),
          axis.title = element_text(size = v_textsize["axis.title"]),
          axis.text = element_text(size = v_textsize["axis.text"]),
          legend.title = element_text(size = v_textsize["legend.title"]),
          legend.text = element_text(size = v_textsize["legend.text"]),
          strip.text = element_text(size = v_textsize["axis.title"]),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # define number of rows in the fill legend
  if (nlevels(d_plot$traitvalue) > 4){
    p <-  p +
      guides(color = guide_legend(nrow = 3, byrow = T))
  } else if (nlevels(d_plot$traitvalue) > 2){
    p <-  p +
      guides(color = guide_legend(nrow = 2, byrow = T))
  } else {
    p <-  p +
      guides(color = guide_legend(nrow = 1, byrow = T))
  }
  
  
  l_plots_sric_trait_butter_suppl_sel[[trait_i]] <-
    p
  
}

# combine the plots for the two traits
plot_grid(plotlist = l_plots_sric_trait_butter_suppl_sel,
          ncol = 2, align = "hv", rel_heights = c(1, 1),
          labels = letters[1:2],
          label_size = v_textsize["plotlabel"])

# save figure
ggsave("Output/Figures/sric_trend_trait_butter_suppl_sel.jpeg", 
       width = 90, height = 70, unit = "mm", dpi = 400)

# overall driver models ########################################################
################################################################################.

# plot model estimates for the overarching regression model relating richness trends and environmental variables

p1 <- 
  # combine beetle and butterfly data
  d_smry_sapro |> 
  mutate(group = "Sapro") |> 
  bind_rows(d_smry_butter |> 
              mutate(group = "Butter")) |> 
  mutate(across(lower95:mean, ~ . * sel_interval * 100), # convert to percentage point change per (8-/12-year) interval
         var = factor(var, levels = rev(v_vars_spec_short), # add environmental factor ordering and labels for plotting
                      labels = rev(names(v_vars_spec_short))),
         group = factor(group, levels = names(v_grouplabels))) |>
  ggplot(aes(y = var, yend = var, colour = group)) +
  geom_vline(xintercept = 0, lty = 2) + # zeroline
  geom_segment(aes(x = lower95, xend = upper95), #95% CI
               linewidth = .25) +
  geom_segment(aes(x = lower90, xend = upper90), #90% CI
               linewidth = .75) +
  geom_segment(aes(x = lower80, xend = upper80), #80% CI
               linewidth = 1.5) +
  geom_point(aes(x = mean), # mean
             size = 2) +
  facet_grid(~ group, scales = "free_x", # facet by insect group
             labeller = as_labeller(sapply(v_grouplabels, # add decent group labels
                                           function(x) paste0("\n\n\n", x)))) +
  scale_colour_manual(values = v_groupcols) +
  # Nicer x axis breaks
  scale_x_continuous(breaks = case_when(sel_interval == 8 ~ c(-4, -2, 0, 2, 4), 
                                        sel_interval == 12 ~ c(-10, -5, 0, 5, 10))) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        strip.text = element_text(size = v_textsize["axis.title"])) +
  labs(x = paste0("Species richness trend change\n(percentage point change per ",
                  sel_interval, " years)"))

# save figure
cairo_pdf(paste0("Output/Figures/lm_ric_mod_both_",
                 sel_interval, "yrs.pdf"),
          width = 90 / 25.4, height = 90 / 25.4, fallback_resolution  = 400)
p1
dev.off()
