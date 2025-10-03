f_occ_det <- function(d_sites, d_visits, 
                      formula_occ, formula_det,
                      var_site, var_year, var_presabs, 
                      var_area = NULL,
                      stan_mod, ...) {
  
  d_sites <- d_sites %>% 
    arrange(!! sym(var_site), !! sym(var_year))
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
    arrange(!! sym(var_site), !! sym(var_year)) 
  
  
  occ <- which(m_occ$occ == 1)
  not_occ <- which(m_occ$not_occ == 1)
  
  start_i <- c(0, cumsum(V)) + 1
  start_i <- start_i[-length(start_i)]
  end_i <- cumsum(V)
  

  m_occ_visit <- d_visits %>% 
    left_join(m_occ,
              by = c(var_site, var_year)) %>% 
    arrange(!! sym(var_site), !! sym(var_year))
  
  occ_visit <- which(m_occ_visit$occ == 1)

  
  L_year <- max(x_year)
  L_area <- max(x_area)
  
  L_ro <- apply(x_ro, 2, max)
  dim(L_ro) <- ncol(x_ro)
  
  if (length(unique(L_ro)) > 1) warning("Not all random varibles in the occupancy model have the same number of levels. But all returned alpha_ro values will have values. Make sure to exclude the non-defined combinations!")
  
  L_rd <- apply(x_rd, 2, max)
  dim(L_rd) <- ncol(x_rd)
  
  if (length(unique(L_rd)) > 1) warning("Not all random varibles in the detectability model have the same number of levels. But all returned alpha_rd values will have values. Make sure to exclude the non-defined combinations!")
  
  l_data <- list(P=P,
                 V = V,
                 N = nrow(d_visits),
                 OCC = length(occ),
                 OCC_N = length(occ_visit),
                 NOCC = length(not_occ),
                 x_year = x_year,
                 x_area = x_area,
                 x_fo = x_fo,
                 x_ro = x_ro,
                 x_fd = x_fd,
                 x_rd = x_rd,
                 K_fo = ncol(x_fo),
                 K_ro = ncol(x_ro),
                 K_fd = ncol(x_fd),
                 K_rd = ncol(x_rd),
                 L_year = L_year,
                 L_area = L_area,
                 L_ro = L_ro,
                 L_rd = L_rd,
                 y = unlist(d_visits[, var_presabs]),
                 occ = occ, 
                 not_occ = not_occ,
                 occ_visit = occ_visit, 
                 start_i = start_i,
                 end_i = end_i)
  
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
