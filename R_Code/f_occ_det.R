# function to prepare data for occupancy-detection models and to fit the models

f_occ_det <- function(d_sites, # data frame containing all the site-year-combinations included in the analysis
                      d_visits, # data frame containing all the vitis included in the analysis
                      formula_occ, # formula for the occurrence probability model
                      formula_det, # formula for the detection probability model
                      var_site, # name of the variable denoting the site
                      var_year, # name of the variable denoting the year
                      var_presabs, # name of the variable denoting the presence-absence of the species
                      var_area = NULL, # name of an optional variable defining regions for which separate random walks are fitted
                      stan_mod, # Stan model to be fitted
                      ...) {
  
  # arrange site data frame along sites and years
  d_sites <- d_sites %>% 
    arrange(!! sym(var_site), !! sym(var_year))
  # arrange visits along sites and years (imported that both data frames have the same ordering!!)
  d_visits <- d_visits %>% 
    arrange(!! sym(var_site), !! sym(var_year))
  
  P <- nrow(d_sites) # number of site x year combinations
  
  # determine number of visits for each Q*Y combination:
  V <- d_visits %>% 
    mutate(visit = 1) %>% 
    full_join(d_sites %>% select(!! sym(var_site), !! sym(var_year)),
              by = c(var_site, var_year)) %>% 
    group_by(!! sym(var_site), !! sym(var_year)) %>% 
    summarise(V = sum(visit, na.rm = T),
              .groups = "drop") %>% 
    arrange(!! sym(var_site), !! sym(var_year)) %>% 
    select(V) %>% deframe()
  
  # get terms of occupancy formula
  terms_occ <- attr(terms(formula_occ), "term.labels")
  terms_fixed_occ <- terms_occ[!grepl("\\|", terms_occ)]
  terms_random_occ <- terms_occ[grepl("\\|", terms_occ)]
  terms_random_occ <- substr(terms_random_occ, 2, nchar(terms_random_occ))
  terms_random_occ <- gsub(" |\\|", "", terms_random_occ)
  
  # year effect occupancy (random walk)
  ff <- formula(paste("~ ", var_year))
  x_year <- model.matrix(ff , data = d_sites %>% 
                           mutate_at(vars(!! sym(var_year)), 
                                     ~as.numeric(as.factor(.))))
  x_year <- x_year[, -1, drop = F] # exclude intercept
  
  # if separate random walks per region are run
  if (!is.null(var_area)){
    ff <- formula(paste("~ ", var_area))
    x_area <- model.matrix(ff , data = d_sites %>% 
                             mutate_at(vars(!! sym(var_area)), 
                                       ~as.numeric(as.factor(.))))
    x_area <- x_area[, -1, drop = F] # exclude intercept
  } else { # if only one global random run is used
    x_area <- matrix(rep(1, nrow(d_sites)), ncol = 1)
  }
  
  
  # model matrix fixed effects occupancy
  if (length(terms_fixed_occ) > 0) {
    ff <- formula(paste("~ ", paste(terms_fixed_occ, collapse = " + ")))
  } else {
    ff <- formula(~ 1)
  }
  x_fo <-  model.matrix(ff, data = d_sites)
  x_fo <- x_fo[, -1, drop = F] # exclude intercept
  
  # model matrix random effects occupancy
  if (length(terms_random_occ) > 0){
    ff <- formula(paste("~ ", paste(terms_random_occ, collapse = " + ")))
  } else {
    ff <- formula(~ 1)
  }
  x_ro <- model.matrix(ff, data = d_sites %>% 
                         mutate_at(vars(!!! syms(terms_random_occ)), 
                                   ~as.numeric(as.factor(.))))
  x_ro <- x_ro[, -1, drop = F] # exclude intercept
  
  # get terms of detectability formula
  terms_det <- attr(terms(formula_det), "term.labels")
  terms_fixed_det <- terms_det[!grepl("\\|", terms_det)]
  terms_random_det <- terms_det[grepl("\\|", terms_det)]
  terms_random_det <- substr(terms_random_det, 2, nchar(terms_random_det))
  terms_random_det <- gsub(" |\\|", "", terms_random_det)
  
  # model matrix fixed effects detectability
  if (length(terms_fixed_det) > 0){
    ff <- formula(paste("~ ", paste(terms_fixed_det, collapse = " + ")))
  } else {
    ff <- formula(~ 1)
  }
  x_fd <-  model.matrix(ff, data = d_visits)
  x_fd <- x_fd[, -1, drop = F] # exclude intercept
  
  # model matrix random effects detupancy
  if (length(terms_random_det) > 0){
    ff <- formula(paste("~ ", paste(terms_random_det, collapse = " + ")))
  } else {
    ff <- formula(~ 1)
  }
  x_rd <- model.matrix(ff, data =   d_visits %>% 
                         mutate_at(vars(!!! syms(terms_random_det)), 
                                   ~as.numeric(as.factor(.))))
  x_rd <- x_rd[, -1, drop = F] # exclude intercept
  
  # determine for each site-year combination, if there was an occurrence registered for the target species
  m_occ <- d_visits %>% 
    full_join(d_sites %>% select(!! sym(var_site), !! sym(var_year)),
              by = c(var_site, var_year)) %>% 
    group_by(!! sym(var_site), !! sym(var_year)) %>% 
    summarise(occ = sum(!! sym(var_presabs), na.rm = T),
              not_occ = sum(!! sym(var_presabs) == 0, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(occ = ifelse(occ > 0, 1, occ),
           not_occ = ifelse(occ > 0, 0, not_occ), # only 1 if not occupied at all
           not_occ = ifelse(not_occ > 0, 1, not_occ)) %>% 
    arrange(!! sym(var_site), !! sym(var_year))  # make sure the order is correct
  
  # which site-year combination are occupied / not occupied
  occ <- which(m_occ$occ == 1)
  not_occ <- which(m_occ$not_occ == 1)
  
  # in a vector containing all the visits, where does a new site-year combination start? And where does it end?
  # create two vectors for this
  start_i <- c(0, cumsum(V)) + 1
  start_i <- start_i[-length(start_i)]
  end_i <- cumsum(V)

  # add occupied / not-occupied data to the visits dataframe
  m_occ_visit <- d_visits %>% 
    left_join(m_occ,
              by = c(var_site, var_year)) %>% 
    arrange(!! sym(var_site), !! sym(var_year)) # make sure the order is correct
  
  # which are visits of site-year combinations, in which the species was recorded at all?
  occ_visit <- which(m_occ_visit$occ == 1)

  # number of years (or whatever time unit used)
  L_year <- max(x_year)
  # number of regions with separate random walks
  L_area <- max(x_area)
  
  # number of levels of the different random terms in the occurrence probability model
  L_ro <- apply(x_ro, 2, max)
  dim(L_ro) <- ncol(x_ro)
  
  # Stan model fits for all entries in the L_ro matrix, even if they are not defined. This is only problematic, if they are later interpreted
  if (length(unique(L_ro)) > 1) warning("Not all random varibles in the occupancy model have the same number of levels. But all returned alpha_ro values will have values. Make sure to exclude the non-defined combinations!")
  
  # number of levels of the different random terms in the detection probability model
  L_rd <- apply(x_rd, 2, max)
  dim(L_rd) <- ncol(x_rd)
  
  # Stan model fits for all entries in the L_rd matrix, even if they are not defined. This is only problematic, if they are later interpreted
  if (length(unique(L_rd)) > 1) warning("Not all random varibles in the detectability model have the same number of levels. But all returned alpha_rd values will have values. Make sure to exclude the non-defined combinations!")
  
  # prepare data list
  l_data <- list(P = P, # number of unique site x year combinations
                 V = V, # number of visits per site x year combination (vector)
                 N = nrow(d_visits), # total number of visits
                 OCC = length(occ), # number of occupied sites
                 OCC_N = length(occ_visit), # number of occupied sites
                 NOCC = length(not_occ), # number of not occupied sites (but visited)
                 x_year = x_year, # year of observation
                 x_area = x_area, # region of observation (for random walks)
                 x_fo = x_fo, # fixed variables in the occurrence probability model
                 x_ro = x_ro, # random variables in the occurrence probability model (integer values!)
                 x_fd = x_fd, # fixed variables in the detection probability model
                 x_rd = x_rd, # random variables in the detection probability model (integer values!)
                 K_fo = ncol(x_fo), # number of fixed variables in the occurrence probability model
                 K_ro = ncol(x_ro), # number of random variables in the occurrence probability model
                 K_fd = ncol(x_fd), # number of fixed  variables in the detection probability model
                 K_rd = ncol(x_rd), # number of random  variables in the detection probability model
                 L_year = L_year, # number of years (or whatever time unit used)
                 L_area = L_area, # number of regions for which separate random walks are estimated
                 L_ro = L_ro, # number of levels of random variables in the occurrence probability model
                 L_rd = L_rd,# number of levels of random variabless in the detection probability model
                 y = unlist(d_visits[, var_presabs]), # observation (pres-abs)
                 occ = occ, # vector of occupied site-year combinations
                 not_occ = not_occ, # vector not occupied site-year combinations
                 occ_visit = occ_visit, # vector of occupied site-year combinations (visit level)
                 start_i = start_i, # start of vector entries per site-year combination
                 end_i = end_i) # end of vector entries per site-year combination
  
  # run model with Cmdstanr or rstan
  if (any(class(stan_mod) == "CmdStanModel")){
    stan_mod$sample(data = l_data, ...)
  } else if (any(class(stan_mod) == "stanmodel")){
    sampling(object = stan_mod, data = l_data,
             ...)
  } else {
    warning("stan_mod not of class CmdStanModel or stanmodel")
  }
}
