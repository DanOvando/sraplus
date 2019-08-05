#' Run sraplus
#'
#' @param driors a list of driors passed from sraplus::format_driors
#' @param include_fit logical indicating whether to return the fitted object
#' @param seed seed for model runs
#' @param plim cutoff (in units of B/K) for hockey stick PT function
#' @param model the name of the sraplus TMB version to be run
#' @param fit_catches logical indicating whether catches should be fit or passed
#' @param randos random effects when passing to TMB
#' @param draws the number of SIR samples to run
#' @param engine one of 'sir','stan', or 'tmb'
#' @param cores number of cores for stan fits
#' @param chains number of chains for stan fits
#' @param n_keep the number of SIR samples to keep
#'
#' @return a fitted sraplus object
#' @export
#'
#'
fit_sraplus <- function(driors,
                        include_fit = TRUE,
                        seed = 42,
                        plim = 0.2,
                        model = "sraplus_tmb",
                        fit_catches = TRUE,
                        randos = "uc_proc_errors",
                        draws = 1e6,
                        n_keep = 2000,
                        engine = "sir",
                        cores = 4,
                        chains = 1,
                        cleanup = FALSE,
                        max_treedepth = 10,
                        adapt_delta = 0.8,
                        estimate_m = FALSE,
                        estimate_qslope = FALSE,
                        estimate_proc_error = TRUE
) {
  knockout <-
    list() #parameters to knockout from TMB estimation using TMB::map
  
  index_years <- driors$index_years
  
  effort_years <- driors$effort_years
  
  time <- length(driors$catch)
  
  if (all(is.na(driors$effort))) {
    index_t = driors$index
    
  } else {
    index_t = rep(0, length(driors$effort))
  }
  
  sra_data <- list(
    catch_t = driors$catch,
    index_t = index_t,
    effort_t = driors$effort,
    index_years = which(driors$index_years %in% driors$years),
    log_r_prior = log(driors$growth_rate),
    log_r_cv = driors$growth_rate_cv,
    log_init_dep_prior = log(driors$initial_b),
    log_init_dep_cv = driors$initial_b_cv,
    log_final_dep_prior = log(driors$terminal_b),
    log_final_dep_cv = driors$terminal_b_cv,
    time = time,
    fit_index = as.numeric(!all(is.na(driors$index)) |
                             !all(is.na(driors$effort))),
    calc_cpue =  as.numeric(!all(is.na(driors$effort))),
    use_u_prior = as.numeric(!all(is.na(driors$u_v_umsy))),
    u_years = which(driors$u_years %in% driors$years),
    u_priors = driors$u_v_umsy,
    u_cv = driors$u_cv,
    plim = plim,
    sigma_proc_prior = driors$sigma_r,
    sigma_proc_prior_cv = driors$sigma_r_cv,
    ref_type = ifelse(driors$ref_type == "k", 0, 1),
    use_final = !is.na(driors$terminal_b),
    use_final_u = as.numeric(!all(is.na(
      driors$log_final_u
    ))),
    log_final_u = driors$log_final_u,
    log_final_u_cv = driors$log_final_u_cv,
    use_init =  !is.na(driors$initial_b),
    sigma_u = driors$u_cv,
    log_k_guess = log(10 * max(driors$catch)),
    # f_cv = driors$f_cv,
    # q_slope = driors$q_slope,
    eps = 1e-3,
    nat_m = driors$m
  )
  
  k_guess <- log(10 * max(driors$catch))
  
  
  if (sra_data$fit_index == 1 & sra_data$calc_cpue == 0) {
    q_guess = pmin(1e-2, median((sra_data$index_t / sra_data$catch_t[sra_data$index_years])))
    
    
  } else if (sra_data$calc_cpue == 1) {
    q_guess = pmin(1e-2, median(0.2 / sra_data$effort_t))
    
  } else {
    q_guess <- 1e-2
    
    
  }
  
  sra_data$log_q_guess = log(q_guess)
  
  inits <- list(
    log_k = log(10 * max(driors$catch)),
    log_r = log(driors$growth_rate),
    # q = q_guess,
    log_q = log(q_guess),
    log_sigma_obs = log(0.2),
    log_init_dep = log(1),
    log_sigma_proc = log(0.01),
    uc_proc_errors = rep(0, time - 1),
    log_m = log(2),
    q_slope = ifelse(estimate_qslope == TRUE && sra_data$calc_cpue == 1,0.025,0)
  )
  
  
  if (sra_data$fit_index == 0) {
    knockout$log_q <- NA
    # knockout$q <- NA
    
    knockout$log_sigma_obs <- NA
    
    
  }
  
  # knockout$log_init_dep = NA
  
  if (sra_data$fit_index == 0 & sra_data$use_u_prior == 0) {
    knockout$log_sigma_proc <- NA
    
    knockout$uc_proc_errors <- rep(NA, time - 1)
    
    inits$log_sigma_proc <- log(1e-6)
    
    inits$uc_proc_errors <- rep(0, time - 1)
    
    randos <- NULL
    # randos <- "inv_f_t"
    
    # randos <- "log_f_t"
  }
  
  if (estimate_m == FALSE) {
    knockout$log_m <- NA
  }
  
  if (!(estimate_qslope == TRUE && sra_data$calc_cpue == 1)){
    
    knockout$q_slope <- NA
    
  }
  
  if (estimate_proc_error == FALSE){
    
    knockout$uc_proc_errors <- rep(NA, time - 1)
    
    knockout$log_sigma_proc <- NA
    
    inits$log_sigma_proc <- log(1e-6)
    
    randos <- NULL
    
  }
  
  knockout <- purrr::map(knockout, as.factor)
  
  
  # fit models
  
  if ((sra_data$fit_index == 0 &
       sra_data$use_u_prior == 0) | engine == "sir") {
    sra_fit <- sraplus::sraplus(
      catches = sra_data$catch_t,
      r = pmax(
        0.01,
        rnorm(draws, driors$growth_rate, driors$growth_rate_cv)
      ),
      m = runif(draws, 0.2, 6),
      init_dep = exp(
        rnorm(
          draws,
          sra_data$log_init_dep_prior,
          sra_data$log_init_dep_cv
        )
      ),
      k = runif(
        draws,
        1.15 * max(sra_data$catch_t),
        50 * max(sra_data$catch_t)
      ),
      sigma_procs = runif(draws, 0, 0.15),
      draws = draws,
      log_final_ref = ifelse(
        is.na(sra_data$log_final_dep_prior),
        0.5,
        sra_data$log_final_dep_prior
      ),
      sigma_dep = ifelse(
        is.na(sra_data$log_final_dep_prior),
        1,
        sra_data$log_final_dep_cv
      ),
      u_prior = sra_data$use_u_prior,
      u_priors = sra_data$u_priors,
      u_years = sra_data$u_years,
      sigma_u = sra_data$u_cv,
      ref_type = sra_data$ref_type,
      n_keep = n_keep,
      drawdex = 0:(draws - 1),
      qs = runif(draws, 1e-9, 1e-1),
      index_t = sra_data$index_t,
      index_years = sra_data$index_years,
      fit_index = sra_data$fit_index,
      sigma_obs = exp(rnorm(draws, log(0.2), 0.1)),
      plim = plim,
      use_final_u = sra_data$use_final_u,
      log_final_u = sra_data$log_final_u,
      log_final_u_cv =  sra_data$log_final_u_cv
    )
    
    
    keepers <- sra_fit$keepers
    
    outs <- stringr::str_detect(names(sra_fit), "_t")
    
    sra_fit$b_t[, keepers] -> a
    
    tidy_fits <-
      purrr::map_df(
        purrr::keep(sra_fit, outs),
        ~ as.data.frame(.x[, keepers]) %>% dplyr::mutate(year = 1:nrow(.)) %>% tidyr::gather(draw, value, -year),
        keepers = keepers,
        .id = "variable"
      ) %>%
      dplyr::mutate(draw = stringr::str_replace_all(draw, "\\D", "") %>% as.numeric())
    
    out <- tidy_fits %>%
      dplyr::group_by(year, variable) %>%
      dplyr::summarise(
        mean = mean(value),
        sd = sd(value),
        lower = quantile(value, 0.1),
        upper = quantile(value, 0.9)
      ) %>%
      dplyr::ungroup()
    
    out$variable <-
      dplyr::case_when(
        out$variable == "b_bmsy_t" ~ "b_div_bmsy",
        out$variable == "b_t" ~ "b",
        out$variable == "c_msy_t" ~ "c_div_msy",
        out$variable == "dep_t" ~ "depletion",
        out$variable == "u_umsy_t" ~ "u_div_umsy",
        TRUE ~ out$variable
      )
    out <- list(results = out,
                fit = tidy_fits)
    
    out$results$year <- out$results$year - 1 + min(driors$years)
    
    out$fit$year <- out$fit$year - 1 + min(driors$years)
    
    # index_tests <- sra_fit$index_hat[, sra_fit$keepers]
    
    
  }
  else if (engine == "stan") {
    if (sra_data$use_u_prior == 0) {
      udist <- NULL
      
    } else {
      udist <- 2
      
    }
    
    lower_k <- log(1.25 * max(driors$catch))
    
    upper_k <- log(50 * max(driors$catch))
    
    if (fit_catches == TRUE) {
      # inits$inv_f_t = rep(-2, time - 1)
      
      # inits$log_f_t = rep(-2, time - 1)
      
      
    }
    
    
    sraplus::get_tmb_model(model_name = model)
    
    sra_model <-
      TMB::MakeADFun(
        data = sra_data,
        parameters = inits,
        DLL = model,
        random = randos,
        silent = TRUE,
        inner.control = list(maxit = 1e3),
        hessian = FALSE,
        map = knockout
      )
    
    lower = rep(-Inf, length(sra_model$par)) %>%
      purrr::set_names(names(sra_model$par))
    
    # lower['q'] <- 1e-10
    lower['log_k'] <- lower_k
    
    upper = rep(Inf, length(sra_model$par)) %>%
      purrr::set_names(names(sra_model$par))
    
    upper['log_k'] <- upper_k
    
    upper["log_init_dep"] <- log(1.5)
    
    upper['log_q'] <- log(1)
    set.seed(seed)
    fit <-
      tmbstan::tmbstan(
        sra_model,
        lower = lower,
        upper = upper,
        cores = cores,
        chains = chains,
        iter = n_keep,
        init = purrr::map(chains,  ~ purrr::map(inits, jitter), inits = inits),
        control = list(max_treedepth = max_treedepth,
                       adapt_delta = adapt_delta)
      )
    
    draws = tidybayes::tidy_draws(fit) %>%
      tidyr::nest(-.chain, .iteration,-.draw,-.iteration)
    
    draws <- draws %>%
      dplyr::mutate(
        pars = purrr::map(
          data,
          get_posterior,
          inits = inits,
          sra_data = sra_data,
          model = model,
          randos = randos,
          knockout = knockout
        )
      ) %>%
      dplyr::select(-data)
    
    draws <- draws %>%
      dplyr::mutate(stack = purrr::map(pars, stack_stan)) %>%
      dplyr::select(-pars) %>%
      tidyr::unnest()
    
    draws <- draws %>%
      dplyr::group_by(variable, .draw) %>%
      dplyr::mutate(year = seq_along(value))
    
    
    draws$variable <-
      dplyr::case_when(
        draws$variable == "log_b" ~ "log_b_div_bmsy",
        draws$variable == "log_bt" ~ "log_b",
        draws$variable == "log_dep" ~ "log_depletion",
        draws$variable == "log_u" ~ "log_u_div_umsy",
        TRUE ~ draws$variable
      )
    
    
    logs <- draws %>%
      dplyr::ungroup() %>%
      dplyr::filter(stringr::str_detect(variable, "log_")) %>%
      dplyr::mutate(value = exp(value)) %>%
      dplyr::mutate(variable = stringr::str_remove_all(variable, "log_"))
    
    draws <- draws %>%
      dplyr::filter(!stringr::str_detect(variable, "log_")) %>%
      dplyr::bind_rows(logs)
    
    out <- logs %>%
      dplyr::group_by(variable, year) %>%
      dplyr::summarise(
        mean = mean(value),
        sd = sd(value),
        lower = quantile(value, 0.1),
        upper = quantile(value, 0.9)
      ) %>%
      dplyr::ungroup()
    
    if (include_fit == FALSE) {
      fit = NA
    }
    
    out <- list(results = out,
                fit = fit)
    
    rm(sra_model)
    
    rm(fit)
    
    # dyn.unload(TMB::dynlib(file.path("tmb", model)))
    
    
  }
  else {
    # fit TMB model
    
    if (sra_data$use_u_prior == 0) {
      udist <- NULL
      
    } else {
      udist <- 2
      
    }
    
    lower_k <- log(1.25 * max(driors$catch))
    
    upper_k <- log(50 * max(driors$catch))
    
    if (fit_catches == TRUE) {
      # inits$inv_f_t = rep(-2, time - 1)
      # inits$log_f_t = rep(-2, time - 1)
      
    }
    
    sraplus::get_tmb_model(model_name = model)
    sra_model <-
      TMB::MakeADFun(
        data = sra_data,
        parameters = inits,
        DLL = model,
        random = randos,
        silent = TRUE,
        inner.control = list(maxit = 1e6),
        hessian = TRUE,
        map = knockout
      )
    
    lower = rep(-Inf, length(sra_model$par)) %>%
      purrr::set_names(names(sra_model$par))
    
    # lower['q'] <- 1e-8
    
    lower['log_k'] <- lower_k
    
    upper = rep(Inf, length(sra_model$par)) %>%
      purrr::set_names(names(sra_model$par))
    
    upper['log_k'] <- upper_k
    
    upper['log_q'] <- 0
    
    upper["log_init_dep"] <- log(1.5)
    
    set.seed(seed)
    fit <- TMBhelper::fit_tmb(
      sra_model,
      fn = sra_model$fn,
      gr = sra_model$gr,
      newtonsteps = 6,
      lower = lower,
      upper = upper,
      getsd = FALSE,
      control = list(
        eval.max = 1e3,
        iter.max = 1e3,
        rel.tol = 1e-10
      )
    )
    
    
    if (fit$max_gradient > 1e-3) {
      fit <- TMBhelper::fit_tmb(
        sra_model,
        fn = sra_model$fn,
        gr = sra_model$gr,
        newtonsteps = 10,
        lower = lower,
        upper = upper,
        getsd = FALSE,
        control = list(
          eval.max = 1e3,
          iter.max = 1e3,
          rel.tol = 1e-10
        )
      )
      
    }
    
    fit_save <- sra_model
    
    fit_report <- fit_save$report()
    
    fit <- TMB::sdreport(fit_save, bias.correct = ifelse(is.null(randos), FALSE,TRUE))
    
    out <-
      dplyr::tibble(
        variable = names(fit$value),
        mean = fit$value,
        sd = fit$sd
      )
    
    out$variable <-
      dplyr::case_when(
        out$variable == "log_b" ~ "log_b_div_bmsy",
        out$variable == "log_bt" ~ "log_b",
        out$variable == "log_dep" ~ "log_depletion",
        out$variable == "log_u" ~ "log_u_div_umsy",
        TRUE ~ out$variable
      )
    
    
    out <- out %>%
      dplyr::mutate(lower = mean - 1.96 * sd,
                    upper = mean + 1.96 * sd)
    
    logs <- out %>%
      dplyr::filter(stringr::str_detect(variable, "log_")) %>%
      dplyr::mutate(
        mean = exp(mean),
        lower = exp(lower),
        upper = exp(upper),
        variable = stringr::str_remove_all(variable, "log_")
      )
    
    out <- out %>%
      dplyr::bind_rows(logs)
    
    if (include_fit == FALSE) {
      fit = NA
    }
    
    out <- list(results = out,
                fit = fit)
    
    rm(sra_model)
    
    rm(fit)
    
    # dyn.unload(TMB::dynlib(file.path("tmb", model)))
    
    
  } # close else
  
  if (cleanup == TRUE) {
    unlink(file.path(getwd(), "tmb"), recursive = TRUE)
  }
  
  # ldll <- getLoadedDLLs()
  #
  # if (any(names(ldll) == model)){
  #
  # temp <- ldll[names(ldll) == model][[1]]
  #
  # dyn.unload(temp[['path']])
  #
  # }
  #
  # put years back
  
  
  return(out)
  
  
}