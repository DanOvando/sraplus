#' Run sraplus
#'
#' The main function for taking an object created by \code{format_driors} and fitting a
#' model using \code{sraplus}
#'
#' @param driors a list of driors passed from sraplus::format_driors
#' @param include_fit logical indicating whether to return the fitted object
#' @param seed seed for model runs
#' @param plim cutoff (in units of B/K) for hockey stick PT function
#' @param model_name the name of the sraplus TMB version to be run, defaults to "sraplus_tmb"
#' @param randos random effects when passing to TMB
#' @param draws the number of SIR samples to run
#' @param engine one of 'sir','stan', or 'tmb'
#' @param cores number of cores for stan fits
#' @param chains number of chains for stan fits
#' @param n_keep the number of SIR samples to keep
#' @param cleanup logical indicating whether to remove the compiled TMB model after running
#' @param max_treedepth max_treedepth for models fit using stan
#' @param adapt_delta adap_delta for models fit using stan
#' @param estimate_shape logical indicating whether to estimate the shape parameter of the Pella-Tomlinson model.
#' If FALSE shape parameter is held at the initial value, either default or supplied
#' @param estimate_qslope logical indicating whether to estimate a slope parameter for q in CPUE fitting.
#' If FALSE q_slope is held at the initial value, either default or supplied. if TRUE, \code{estimate_proc_error} should be set to FALSE
#' @param estimate_proc_error logical indicating whether to estimate process errors.
#' If FALSE process errors are not included in the model
#' @param ci confidence/credible interval range for summaries
#'
#' @return a fitted sraplus object
#' @export
#' @examples
#' \dontrun{
#'
#' catch_only_driors <- sraplus::format_driors(
#' taxa = example_taxa,
#' catch = cod$catch,
#' years = cod$year,
#' use_heuristics = TRUE
#' )
#'
#'  catch_only_fit <- sraplus::fit_sraplus(driors = catch_only_driors,
#'engine = "sir",
#'draws = 1e6,
#'n_keep = 2000)
#'
#' }
fit_sraplus <- function(driors,
                        include_fit = TRUE,
                        seed = 42,
                        plim = 0.2,
                        model_name = "sraplus_tmb",
                        randos = "log_proc_errors",
                        draws = 1e6,
                        n_keep = 2000,
                        engine = "sir",
                        cores = 4,
                        chains = 1,
                        cleanup = FALSE,
                        max_treedepth = 10,
                        adapt_delta = 0.8,
                        estimate_shape = FALSE,
                        estimate_qslope = FALSE,
                        estimate_proc_error = TRUE,
                        estimate_k = FALSE,
                        learn_rate = 1e-3,
                        marginalize_q = FALSE,
                        use_baranov = TRUE,
                        include_m = FALSE,
                        ci = 0.89,
                        try_again = FALSE,
                        eps = 1e-6,
                        max_time = Inf,
                        eval.max = 200,
                        iter.max = 150,
                        rel.tol = 1e-10,
                        loopnum = 1,
                        newtonsteps = 1) {
  
  if (max_time < Inf){
    
    setTimeLimit(elapsed = max_time)
    
  }
  
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
    index_years = which(driors$years %in% driors$index_years),
    log_r_prior = log(driors$growth_rate_prior),
    log_r_cv = driors$growth_rate_prior_cv,
    log_init_dep_prior = log(driors$initial_state),
    log_init_dep_cv = driors$initial_state_cv,
    log_final_dep_prior = log(driors$terminal_state),
    log_final_dep_cv = driors$terminal_state_cv,
    time = time,
    fit_index = as.numeric(!all(is.na(driors$index)) |
                             !all(is.na(driors$effort))),
    calc_cpue =  as.numeric(!all(is.na(driors$effort))),
    use_baranov = as.numeric(use_baranov),
    use_u_prior = as.numeric(!all(is.na(driors$u_v_umsy))),
    u_years = which(driors$years %in% driors$u_years),
    u_priors = driors$u_v_umsy,
    u_cv = driors$u_cv,
    plim = plim,
    sigma_proc_prior = driors$sigma_r_prior,
    sigma_proc_prior_cv = driors$sigma_r_prior_cv,
    b_ref_type = ifelse(driors$b_ref_type == "k", 0, 1),
    f_ref_type = ifelse(driors$f_ref_type == "f", 0, 1),
    use_final_state = !is.na(driors$terminal_state),
    use_final_u = as.numeric(!all(is.na(
      driors$log_final_u
    ))),
    log_final_u = driors$log_final_u,
    log_final_u_cv = driors$log_final_u_cv,
    use_init =  !is.na(driors$initial_state),
    sigma_u = driors$u_cv,
    k_prior = driors$k_prior,
    k_prior_cv = driors$k_prior_cv,
    q_slope_prior = driors$q_slope_prior,
    q_slope_cv = driors$q_slope_prior_cv,
    eps = eps,
    nat_m = ifelse(include_m, driors$m, 0),
    shape_prior = driors$shape_prior,
    shape_cv = driors$shape_prior_cv,
    sigma_obs_prior = driors$sigma_obs_prior,
    sigma_obs_prior_cv = driors$sigma_obs_prior_cv,
    marginalize_q = marginalize_q,
    est_k = estimate_k,
    learn_rate = learn_rate
  )
  
  
  if (sra_data$fit_index == 1 & sra_data$calc_cpue == 0) {
    q_prior = pmin(1e-2, median((sra_data$index_t / sra_data$catch_t[sra_data$index_years])))
    
    
  } else if (sra_data$calc_cpue == 1) {
    q_prior = pmin(1e-2, median(0.2 / sra_data$effort_t))
    
  } else {
    q_prior <- 1e-2
    
    
  }
  
  if (estimate_qslope == TRUE & estimate_proc_error == TRUE) {
    warning(
      "Trying to estiamte both qslope and process error, are you sure you want to do that? We recommend setting either estimate_qslope = FALSE or estimate_proc_error = FALSE"
    )
    
  }
  
  sra_data$q_prior = q_prior
  
  sra_data$q_prior_cv = driors$q_prior_cv
  
  inits <- list(
    log_anchor = ifelse(estimate_k, log(driors$k_prior),log(0.42)),
    log_r = log(driors$growth_rate_prior),
    # q = q_guess,
    log_q = log(q_prior),
    log_sigma_obs = log(0.2),
    log_init_dep = log(1),
    log_sigma_proc = log(0.01),
    log_proc_errors = rep(0, time - 1),
    log_shape = log(driors$shape_prior),
    log_q_slope = log(
      ifelse(
        estimate_qslope == TRUE &&
          sra_data$calc_cpue == 1,
        0.025,
        driors$q_slope_prior
      )
    )
  )
  
  
  if (sra_data$fit_index == 0) {
    knockout$log_q <- NA
    # knockout$q <- NA
    
    knockout$log_sigma_obs <- NA
    
    
  }
  
  # knockout$log_init_dep = NA
  
  if (sra_data$fit_index == 0 & sra_data$use_u_prior == 0) {
    knockout$log_sigma_proc <- NA
    
    knockout$log_proc_errors <- rep(NA, time - 1)
    
    inits$log_sigma_proc <- log(1e-6)
    
    inits$log_proc_errors <- rep(0, time - 1)
    
    randos <- NULL
    # randos <- "inv_f_t"
    
    # randos <- "log_f_t"
  }
  
  if (estimate_shape == FALSE) {
    knockout$log_shape <- NA
  }
  
  if (!(estimate_qslope == TRUE && sra_data$calc_cpue == 1)) {
    knockout$log_q_slope <- NA
    
  }
  
  if (sra_data$calc_cpue == TRUE & marginalize_q == 1) {
    sra_data$marginalize_q <-  0
    
    warning(
      "You can't calculate CPUE and marginalize q, defaulting to calculating CPUE. Either set marginalize_q = 0 or calc_cpue = 0"
    )
    
  }
  
  
  if (sra_data$calc_cpue == FALSE & marginalize_q == 1) {
    knockout$log_q <- NA
    
  }
  
  
  if (estimate_proc_error == FALSE) {
    knockout$log_proc_errors <- rep(NA, time - 1)
    
    knockout$log_sigma_proc <- NA
    
    inits$log_sigma_proc <- log(1e-6)
    
    randos <- NULL
    
  }
  
  knockout <- purrr::map(knockout, as.factor)
  
  # fit models
  if (sra_data$fit_index == 0 | engine == "sir") {
    if (engine != "sir") {
      warning("You tried to fit a model with nothing in the likelihood - using SIR instead")
      engine <- sir
    }
    
    if (sra_data$use_final_state == 0 &
        sra_data$use_final_u == 0 &
        sra_data$use_u_prior == 0 & is.na(driors$terminal_state)) {
      stop(
        "Trying to run SIR without priors on final status or fishing mortality! Specify one or both of these"
      )
    }
    sra_fit <- sraplus::sraplus(
      catches = sra_data$catch_t,
      r = pmax(
        0.01,
        rnorm(
          draws,
          driors$growth_rate_prior,
          driors$growth_rate_prior_cv
        )
      ),
      m = runif(
        draws,
        ifelse(estimate_shape, 0.2, sra_data$shape_prior),
        ifelse(estimate_shape, 6, sra_data$shape_prior)
      ),
      init_dep = exp(
        rnorm(
          draws,
          sra_data$log_init_dep_prior,
          sra_data$log_init_dep_cv
        )
      ),
      k = rlnorm(
        draws,
        log(10* max(sra_data$catch_t)),
        2
      ),
      sigma_procs = runif(draws, 0, ifelse(estimate_proc_error, 0.2, 0)),
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
      use_u_prior = sra_data$use_u_prior,
      u_priors = sra_data$u_priors,
      u_years = sra_data$u_years,
      sigma_u = sra_data$u_cv,
      b_ref_type = sra_data$b_ref_type,
      f_ref_type = sra_data$f_ref_type,
      n_keep = n_keep,
      drawdex = 0:(draws - 1),
      qs = runif(draws, 1e-9, 1e-1),
      index_t = sra_data$index_t,
      index_years = sra_data$index_years,
      fit_index = sra_data$fit_index,
      sigma_obs = exp(rnorm(draws, log(0.2), 0.1)),
      plim = plim,
      use_final_u = sra_data$use_final_u,
      use_final_state = sra_data$use_final_state,
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
    
    draw_names <- data.frame(draw = 1:length(keepers), draw_id = keepers)
    
    tidy_fits <- tidy_fits %>% 
      dplyr::left_join(draw_names, by = "draw")
    
    static_outs <- !stringr::str_detect(names(sra_fit), "_t") & names(sra_fit) != "keepers"
    
    
    tidy_static_fits <-   purrr::map_df(
      purrr::keep(sra_fit, static_outs),
      ~ data.frame(value = .x[keepers]) %>% dplyr::mutate(year = 1, draw = 1:nrow(.)),
      keepers = keepers,
      .id = "variable"
    )
    out <- tidy_fits %>%
      dplyr::bind_rows(tidy_static_fits) %>% 
      dplyr::group_by(year, variable) %>%
      dplyr::summarise(
        mean = mean(value),
        sd = sd(value),
        lower = quantile(value, (1 - ci) / 2),
        upper = quantile(value, 1 - (1 - ci) / 2)
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
  else {
    if (sra_data$use_u_prior == 0) {
      udist <- NULL
      
    } else {
      udist <- 2
      
    }
    
    if (estimate_k){
    lower_anchor <- log(1.25 * max(driors$catch))
    
    upper_anchor <- log(50 * max(driors$catch))
    } else {
      lower_anchor = log(1e-3)
      
      upper_anchor = log(0.9);
      
    }
    
    sraplus::get_tmb_model(model_name = model_name)
    sra_model <-
      TMB::MakeADFun(
        data = sra_data,
        parameters = inits,
        DLL = model_name,
        random = randos,
        silent = TRUE,
        inner.control = list(maxit = 1e6),
        hessian = FALSE,
        map = knockout
      )
    
    
    lower = rep(-Inf, length(sra_model$par)) %>%
      purrr::set_names(names(sra_model$par))
    
    # lower['q'] <- 1e-10
    lower['log_anchor'] <- lower_anchor
    
    upper = rep(Inf, length(sra_model$par)) %>%
      purrr::set_names(names(sra_model$par))
    
    upper['log_anchor'] <- upper_anchor
    
    upper["log_init_dep"] <- log(1.5)
    
    if (marginalize_q == 0 & sra_data$fit_index == 1) {
      upper['log_q'] <- log(1)
    }
    if (estimate_qslope == TRUE) {
      lower['log_q_slope'] <- -Inf
      
      upper['log_q_slope'] <- log(.1)
    }
    
    if (estimate_shape == TRUE) {
      upper['log_shape'] <- log(4)
    }
    
    if (engine == "stan") {
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
        group_by(.chain, .iteration, .draw) %>%
        tidyr::nest()
      
      draws <- draws %>%
        dplyr::mutate(
          pars = purrr::map(
            data,
            purrr::quietly(get_posterior),
            inits = inits,
            sra_data = sra_data,
            model_name = model_name,
            randos = randos,
            knockout = knockout
          )
        ) %>%
        dplyr::select(-data)
      
      draws$pars <- purrr::map(draws$pars, "result")
      
      draws <- draws %>%
        dplyr::mutate(stack = purrr::map(pars, stack_stan)) %>%
        dplyr::select(-pars) %>%
        tidyr::unnest(col = stack)
      
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
      
      out <- draws %>%
        dplyr::group_by(variable, year) %>%
        dplyr::summarise(
          mean = mean(value),
          sd = sd(value),
          lower = quantile(value, (1 - ci) / 2),
          upper = quantile(value, 1 - (1 - ci) / 2)
        ) %>%
        dplyr::ungroup()
      
      if (include_fit == FALSE) {
        fit = NA
      }
      
      out <- list(results = out,
                  fit = fit)
      
      out$results$year <- out$results$year - 1 + min(driors$years)
      
      rm(sra_model)
      
      rm(fit)
      
      # dyn.unload(TMB::dynlib(file.path("tmb", model)))
      
      
    } else if (engine == "tmb") {
      set.seed(seed)
      
      fit <- TMBhelper::fit_tmb(
        sra_model,
        fn = sra_model$fn,
        gr = sra_model$gr,
        loopnum = loopnum,
        newtonsteps = newtonsteps,
        lower = lower,
        upper = upper,
        getsd = FALSE,
        control = list(
          eval.max = eval.max,
          iter.max = iter.max,
          rel.tol = rel.tol
        )
      )
      
      if (fit$max_gradient > 1e-3 & try_again == TRUE) {
        fit <- TMBhelper::fit_tmb(
          sra_model,
          fn = sra_model$fn,
          gr = sra_model$gr,
          newtonsteps = newtonsteps * 2,
          lower = lower,
          upper = upper,
          getsd = FALSE,
          control = list(
            eval.max = eval.max * 2,
            iter.max = iter.max * 2,
            rel.tol = rel.tol
          )
        )
        
      }
      
      fit_save <- sra_model
      
      fit_report <- fit_save$report()
      
      fit_diagnostics <- fit$diagnostics
      
      fit_diagnostic_message <- fit$Convergence_check
      
      fit <-
        TMB::sdreport(fit_save, bias.correct = ifelse(is.null(randos), FALSE, TRUE))
      
      fit$diagnostics <- fit_diagnostics
      
      fit$fit_diagnostic_message <- fit_diagnostic_message
      
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
        dplyr::mutate(lower = mean - qnorm(1 - (1 - ci) / 2) * sd,
                      upper = mean +  qnorm(1 - (1 - ci) / 2) * sd) %>% 
        dplyr::group_by(variable) %>% 
        dplyr::mutate(year = seq_along(mean)) %>% 
        ungroup()
      
      
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
      
      out$results$year <- out$results$year - 1 + min(driors$years)
      
      
      rm(sra_model)
      
      rm(fit)
      
    }
    
    
  } # close tmb or stan else
  
  
  if (cleanup == TRUE) {
    model_path <-
      file.path(getwd(),
                paste(model_name, utils::packageVersion("sraplus"), sep = '_v'))
    
    # dyn.unload(TMB::dynlib(file.path(model_path, model_name)))
    
    unlink(model_path, recursive = TRUE)
  }
  
  out$engine <- engine
  
  return(out)
  
  
}