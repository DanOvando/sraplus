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
#' @param estimate_k
#' @param estimate_f
#' @param learn_rate
#' @param analytical_q
#' @param use_baranov
#' @param include_m
#' @param try_again
#' @param eps
#' @param max_time
#' @param eval.max
#' @param iter.max
#' @param rel.tol
#' @param loopnum
#' @param newtonsteps
#' @param tune_prior_predictive
#' @param refresh
#' @param ...
#'
#' @return a fitted sraplus object
#' @export

fit_sraplus <- function(driors,
                        include_fit = TRUE,
                        seed = 42,
                        plim = 0.05,
                        model_name = "sraplus_tmb",
                        randos = "log_proc_errors",
                        draws = 1e5,
                        n_keep = 2000,
                        engine = "sir",
                        cores = 4,
                        chains = 1,
                        cleanup = FALSE,
                        max_treedepth = 10,
                        adapt_delta = 0.8,
                        estimate_shape = FALSE,
                        estimate_qslope = FALSE,
                        estimate_q = TRUE,
                        estimate_proc_error = TRUE,
                        estimate_initial_state = TRUE,
                        estimate_k = TRUE,
                        estimate_f  = FALSE,
                        learn_rate = 1e-3,
                        analytical_q = FALSE,
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
                        newtonsteps = 1,
                        tune_prior_predictive = TRUE,
                        index_fit_tuner = "sir",
                        refresh = 250,
                        log_bias_correct = TRUE,
                        workers = workers,
                        thin_draws = FALSE,
                        thin_rate = 0.5,
                        ...) {
  opts <- list(...)
  
  if (max_time < Inf) {
    setTimeLimit(elapsed = max_time * 60, transient = TRUE)
    
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

  if (log_bias_correct){
    log_terminal_u <- driors$log_terminal_u - driors$log_terminal_u_cv^2/2
      
  } else {
    
    log_terminal_u <- driors$log_terminal_u
  }
  
  sra_data <- list(
    catch_t = driors$catch,
    index_t = index_t,
    effort_t = driors$effort,
    index_years = which(driors$years %in% driors$index_years),
    log_r_prior = log(driors$growth_rate_prior),
    log_r_cv = driors$growth_rate_prior_cv,
    log_init_dep_prior = ifelse(log_bias_correct,log(driors$initial_state) - driors$initial_state_cv^2/2 ,log(driors$initial_state)),
    log_init_dep_cv = driors$initial_state_cv,
    log_terminal_dep_prior = ifelse(log_bias_correct, log(driors$terminal_state) - driors$terminal_state_cv^2/2,log(driors$terminal_state)),
    log_terminal_dep_cv = driors$terminal_state_cv,
    time = time,
    fit_index = as.numeric(!all(is.na(driors$index)) |
                             !all(is.na(driors$effort))),
    calc_cpue =  as.numeric(!all(is.na(driors$effort))),
    use_baranov = as.numeric(use_baranov),
    use_u_prior = as.numeric(!all(is.na(driors$u))),
    u_years = which(driors$years %in% driors$u_years),
    u_priors = driors$u,
    plim = plim,
    sigma_ratio_prior = driors$sigma_ratio_prior,
    sigma_ratio_prior_cv = driors$sigma_ratio_prior_cv,
    b_ref_type = ifelse(driors$b_ref_type == "k", 0, 1),
    f_ref_type = ifelse(driors$f_ref_type == "f", 0, 1),
    use_terminal_state = !is.na(driors$terminal_state),
    use_terminal_u = as.numeric(!all(is.na(
      driors$log_terminal_u
    ))),
    log_terminal_u = log_terminal_u,
    log_terminal_u_cv = driors$log_terminal_u_cv,
    use_init =  !is.na(driors$initial_state),
    u_cv = driors$u_cv,
    sigma_u = driors$u_cv,
    log_k_prior = log(driors$k_prior),
    log_k_prior_cv = driors$k_prior_cv,
    q_slope_prior = driors$q_slope_prior,
    q_slope_cv = driors$q_slope_prior_cv,
    eps = eps,
    nat_m = ifelse(include_m, driors$m, 0),
    shape_prior = driors$shape_prior,
    shape_cv = driors$shape_prior_cv,
    sigma_obs_prior = driors$sigma_obs_prior,
    sigma_obs_prior_cv = driors$sigma_obs_prior_cv,
    analytical_q = analytical_q,
    est_k = estimate_k,
    estimate_initial_state = estimate_initial_state,
    estimate_proc_error = estimate_proc_error,
    estimate_shape = estimate_shape,
    estimate_qslope = estimate_qslope,
    estimate_f = estimate_f,
    f_prior_form = driors$f_prior_form,
    learn_rate = learn_rate
  )
  
  
  if (sra_data$fit_index == 1 & sra_data$calc_cpue == 0 & is.na(driors$q_prior)) {
    q_prior = pmin(1e-2, median((sra_data$index_t / sra_data$catch_t[sra_data$index_years])))
    
    
  } else if (sra_data$calc_cpue == 1) {
    q_prior = pmin(1e-2, median(0.2 / sra_data$effort_t))
    
  } else {
    q_prior <- driors$q_prior
  }
  
  if (estimate_qslope == TRUE & estimate_proc_error == TRUE) {
    warning(
      "Trying to estiamte both qslope and process error, are you sure you want to do that? We recommend setting either estimate_qslope = FALSE or estimate_proc_error = FALSE"
    )
    
  }
  
  sra_data$q_prior = q_prior
  
  sra_data$q_prior_cv = driors$q_prior_cv
  
  inits <- list(
    log_anchor = ifelse(estimate_k, log(driors$k_prior), log(0.42)),
    log_r = log(driors$growth_rate_prior),
    # q = q_guess,
    log_q = log(q_prior),
    # sigma_obs = (driors$sigma_obs_prior),
    log_sigma_obs = log(driors$sigma_obs_prior),
    log_init_dep = log(0.99),
    log_sigma_ratio = log(driors$sigma_ratio_prior + 1e-6),
    log_proc_errors = rep(0, time - 1),
    log_shape = log(driors$shape_prior),
    log_f_t = rep(log(.2), time),
    log_q_slope = log(
      ifelse(
        estimate_qslope == TRUE &&
          sra_data$calc_cpue == 1,
        0.025,
        driors$q_slope_prior + 1e-6
      )
    )
  )
  
  if (estimate_initial_state == 0) {
    knockout$log_init_dep <- NA
  }
  
  if (estimate_f == 0) {
    knockout$log_f_t <- rep(NA, time)
  }
  
  
  if (estimate_q == FALSE){
    
    knockout$log_q <- NA
    
  }
  
  if (sra_data$fit_index == 0) {
    knockout$log_q <- NA
    # knockout$q <- NA
    
    knockout$log_sigma_obs <- NA
    
    # knockout$sigma_obs <- NA
    
  }
  
  # knockout$log_init_dep = NA
  
  if (sra_data$fit_index == 0 & sra_data$use_u_prior == 0) {
    knockout$log_sigma_ratio <- NA
    
    knockout$log_proc_errors <- rep(NA, time - 1)
    
    inits$log_sigma_ratio <- log(1e-6)
    
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
  
  if (sra_data$calc_cpue == TRUE & analytical_q == 1) {
    sra_data$analytical_q <-  0
    
    warning(
      "You can't calculate CPUE and marginalize q, defaulting to calculating CPUE. Either set analytical_q = 0 or calc_cpue = 0"
    )
    
  }
  
  
  if (sra_data$calc_cpue == FALSE & analytical_q == 1) {
    knockout$log_q <- NA
    
  }
  
  
  if (estimate_proc_error == FALSE) {
    knockout$log_proc_errors <- rep(NA, time - 1)
    
    knockout$log_sigma_ratio <- NA
    
    inits$log_sigma_ratio <- log(1e-6)
    
    randos <- NULL
    
  }
  
  knockout <- purrr::map(knockout, as.factor)
  # fit models
  if (sra_data$fit_index == 0 |
      engine == "sir" | index_fit_tuner == "sir") {
    if (engine != "sir" & index_fit_tuner != "sir") {
      warning("You tried to fit a model with nothing in the likelihood - using SIR instead")
      engine <- "sir"
    }
    
    if (sra_data$use_terminal_state == 0 &
        sra_data$use_terminal_u == 0 &
        sra_data$use_u_prior == 0 & is.na(driors$terminal_state)) {
      # stop(
      #   "Trying to run SIR without priors on terminal status or fishing mortality! Specify one or both of these"
      # )
    }
    
    if (estimate_k == TRUE) {
      anchors <- rlnorm(draws,
                        log(10 * max(sra_data$catch_t)),
                        2)
      
      
      
    } else {
      anchors <- rlnorm(draws,
                        log(0.5),
                        0.5)
      
      if (estimate_proc_error == FALSE) {
        anchors = pmin(anchors, 1)
      }
      
    }
    
    init_dep <- exp(
      truncnorm::rtruncnorm(
        draws,
        b = ifelse(estimate_proc_error, log(1.1), log(1)),
        mean = sra_data$log_init_dep_prior,
        sd = sra_data$log_init_dep_cv
      )
    )
    
    if (estimate_proc_error == FALSE) {
      init_dep = pmin(init_dep, 1)
      
    }
    
    if (engine != "sir") {
      # if engine is not sir, then this is being used for prior-predictive tuning
      
      sir_index <- NA
      
      sir_effort <- NA
      
      sir_fit_index <- 0
      
    } else {
      sir_fit_index <- sra_data$fit_index
      
      sir_effort <- sra_data$effort_t
      
      sir_index <- sra_data$index_t
    }
    sra_fit <- sraplus::sraplus(
      catches = sra_data$catch_t,
      rs = pmax(0.005,
                rlnorm(
                  draws,
                  log(driors$growth_rate_prior),
                  driors$growth_rate_prior_cv
                )),
      ms = runif(
        draws,
        ifelse(estimate_shape, 0.2, sra_data$shape_prior),
        ifelse(estimate_shape, 6, sra_data$shape_prior)
      ),
      init_dep = init_dep,
      anchors = anchors,
      sigma_procs = runif(draws, 0, ifelse(estimate_proc_error, 0.2, 0)),
      draws = draws,
      log_terminal_ref = ifelse(
        is.na(sra_data$log_terminal_dep_prior),
        0.5,
        sra_data$log_terminal_dep_prior
      ),
      sigma_dep = ifelse(
        is.na(sra_data$log_terminal_dep_prior),
        1,
        sra_data$log_terminal_dep_cv
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
      index_t = sir_index,
      index_years = sra_data$index_years,
      fit_index = sir_fit_index,
      sigma_obs = exp(rnorm(draws, log(0.2), 0.1)),
      plim = plim,
      use_terminal_u = sra_data$use_terminal_u,
      use_terminal_state = sra_data$use_terminal_state,
      log_terminal_u = sra_data$log_terminal_u,
      log_terminal_u_cv =  sra_data$log_terminal_u_cv,
      estimate_k = estimate_k,
      learn_rate = learn_rate
    )
    
    
    if (tune_prior_predictive == TRUE) {
      # tune samples from SIR such that "posterior" (bernouli filtered) of say terminal depletion matches prior
      # if we state that catches are not data, out prior should not change through the SIR
      
      if (any(!is.na(driors$log_terminal_u))) {
        if (driors$f_ref_type == "fmsy") {
          state <- sra_fit$u_umsy_t[nrow(sra_fit$dep_t), ]
          
        } else {
          state <- sra_fit$u_t[nrow(sra_fit$dep_t), ]
          
          
        }
        
        
        
      } else if (driors$b_ref_type == "k") {
        state <- sra_fit$dep_t[nrow(sra_fit$dep_t), ]
        
      } else {
        state <- sra_fit$b_bmsy_t[nrow(sra_fit$dep_t), ]
        
      }
      
      
      state_breaks <- seq(0, 10, by = .01)
      
      # state_bins <-
      #   cut(state_breaks,
      #       state_breaks,
      #       include.lowest = FALSE,
      #       right = FALSE)
      
      # edge_p <-
      #   pnorm(log(state_breaks),
      #         log(driors$terminal_state),
      #         driors$terminal_state_cv)
      # 
      # p_bin <- dplyr::lead(edge_p) - (edge_p)
      # 
      # bin_frame <- data.frame(bin = state_bins, p_bin = p_bin) %>%
      #   dplyr::mutate(bin = as.character(bin))
      # 
      
        bin_frame <-
        data.frame(state = state, likelihood = sra_fit$likelihood) %>%
        dplyr::mutate(bin = as.character(
          cut(
            state,
            state_breaks,
            include.lowest = FALSE,
            right = FALSE
          )
        )) %>%
        dplyr::group_by(bin) %>%
        dplyr::summarise(p_bin = mean(likelihood, na.rm = TRUE))
      
      #   browser()
      #   
      # bin_frame %>%
      #   ggplot(aes(bin, p_bin)) +
      #   geom_col() +
      #   coord_flip()
      #
      draws <- data.frame(state = state) %>%
        dplyr::mutate(bin = as.character(
          cut(
            state,
            state_breaks,
            include.lowest = FALSE,
            right = FALSE
          )
        )) %>%
        dplyr::left_join(bin_frame, by = "bin") %>%
        dplyr::group_by(bin) %>%
        dplyr::mutate(weight = unique(p_bin) / length(p_bin)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(index = 1:nrow(.))
    
      sra_fit$keepers <-
        sample(draws$index,
               n_keep,
               replace = TRUE,
               prob = draws$weight)
      
    }
    
    keepers <- sra_fit$keepers
    
    outs <- stringr::str_detect(names(sra_fit), "_t")
    
    #     sra_fit$b_t[, keepers] -> a
    # browser()
    
    tidy_fits <-
      purrr::map_df(
        purrr::keep(sra_fit, outs),
        ~ as.data.frame(.x[, keepers]) %>% dplyr::mutate(year = 1:nrow(.)) %>% tidyr::gather(draw, value, -year),
        keepers = keepers,
        .id = "variable"
      ) %>%
      dplyr::mutate(draw = stringr::str_replace_all(draw, "\\D", "") %>% as.numeric())
    
    draw_names <-
      data.frame(draw = 1:length(keepers), draw_id = keepers)
    
    tidy_fits <- tidy_fits %>%
      dplyr::left_join(draw_names, by = "draw")
    
    static_outs <-
      !stringr::str_detect(names(sra_fit), "_t") &
      names(sra_fit) != "keepers"
    
    
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
  
  if (engine != "sir") {
    if (sra_data$use_u_prior == 0) {
      udist <- NULL
      
    } else {
      udist <- 2
      
    }
    
    sraplus::get_tmb_model(model_name = model_name)
    
    lower_init_dep <- log(1e-3)
    
    if (estimate_k) {
      if (index_fit_tuner == "sir") {
        implied_dep_prior <-
          sra_fit$dep_t[length(driors$years), sra_fit$keepers]
        
        implied_k_prior <- log(sra_fit$k[sra_fit$keepers])
        
        implied_r_prior <- log(sra_fit$r[sra_fit$keepers])
        
        implied_r_prior <- log(sra_fit$r[sra_fit$keepers])
        if (driors$b_ref_type == "k") {
          implied_init_dep_prior <- log(sra_fit$dep_t[1, sra_fit$keepers])
        } else {
          implied_init_dep_prior <- log(sra_fit$b_bmsy_t[1, sra_fit$keepers])
        }
        
        sra_data$log_init_dep_prior <- mean(implied_init_dep_prior)
        
        sra_data$log_init_dep_prior_cv <- sd(implied_init_dep_prior)
        
        sra_data$log_r_prior <- mean(implied_r_prior)
        
        sra_data$log_r_cv <- sd(implied_r_prior)
        
        sra_data$log_k_prior <- mean(implied_k_prior)
        
        sra_data$log_k_prior_cv <- sd(implied_k_prior)
        
        lower_anchor <-  0.9 * min(implied_k_prior)
        
        inits$log_anchor <- median(implied_k_prior)
        
        driors$log_k_prior <- mean(implied_k_prior)
        
        upper_anchor <-   1.2 * max(implied_k_prior)
        
        # lower_init_dep <-  0.8 * min(implied_init_dep_prior)
        
      } else {
        temp_inits <- inits
        
        temp_inits$log_r <- log(0.0275)
        
        lks <- seq(1, log(50 * sum(driors$catch)), length.out = 50)
        
        lks <- 200
        pens <- NA
        
        it <- 2000
        itframe <-
          data.frame(
            log_anchor = runif(it, min = log(1), max = log(50 * max(
              driors$catch
            ))),
            log_r =  runif(it, min = log(0.01), max = log(1)),
            pens = NA,
            terminal_dep = NA
          )
        
        
        for (i in 1:it) {
          temp_inits$log_anchor <- itframe$log_anchor[i]
          temp_inits$log_r <- itframe$log_r[i]
          sra_model <-
            TMB::MakeADFun(
              data = sra_data,
              parameters = temp_inits,
              DLL = model_name,
              random = randos,
              silent = TRUE,
              inner.control = list(maxit = 1e6),
              hessian = FALSE,
              map = knockout
            )
          
          sra_model$report() -> a
          
          itframe$pens[i] <- a$pen
          
          itframe$terminal_dep[i] <- a$dep_t[length(a$dep_t)]
          
        }
        # browser()
        #
        # itframe %>%
        #   ggplot(aes((log_r), (log_anchor), color = pens == 0)) +
        #   geom_point()
        
        #
        # lower_anchor <- log(1.25 * max(driors$catch))
        lower_anchor <-
          0.9 * min(itframe$log_anchor[(itframe$pens == 0)])
        
        inits$log_anchor <- log(10 * exp(lower_anchor))
        
        driors$log_k_prior <-
          median(itframe$log_anchor[itframe$pens == 0 &
                                      itframe$terminal_dep < 0.9])
        
        upper_anchor <-
          max(itframe$log_anchor[itframe$pens == 0 &
                                   itframe$terminal_dep < 0.95])
        
      }
      # lower_anchor <- -Inf
      #
      # upper_anchor <- Inf
    } else {
      lower_anchor = log(1e-3)
      
      upper_anchor = log(0.9)
      
      
    }
    
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
    
    # upper["log_init_dep"] <- log(1.5)
    
    # if (estimate_proc_error == FALSE & estimate_initial_state == TRUE){
    
    if (estimate_initial_state == TRUE) {
      upper["log_init_dep"] <- log(1.05)
      
      
      # lower["log_init_dep"] <- lower_init_dep
      
    }
    
    if (sra_data$fit_index == 1) {
      # lower["sigma_obs"] <- 0
    }
    
    if (analytical_q == 0 & sra_data$fit_index == 1 & estimate_q == TRUE) {
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
                         adapt_delta = adapt_delta),
          refresh = refresh
        )
      
      draws = tidybayes::tidy_draws(fit) %>%
        dplyr::group_by(.chain, .iteration, .draw) %>%
        tidyr::nest() %>% 
        dplyr::ungroup()
      
      if (thin_draws){
        
        draws <- draws %>% 
          dplyr::group_by(.chain) %>% 
          dplyr::slice_sample(prop = thin_rate)
        
        
      }

      # a <- Sys.time()

      if (.Platform$OS.type == "windows"){
      
      stacked_draws <- foreach::foreach(i = 1:nrow(draws)) %do% {
      
      qgp <-   purrr::quietly(sraplus::get_posterior)
        
      out <- qgp(draws = draws$data[[i]], inits = inits, sra_data = sra_data,     model_name = model_name,
                 randos = randos,
                 knockout = knockout)
      
      pars <- names(out$result)
      
      flatstack <- purrr::map2_df(pars,out$result, ~data.frame(variable = .x, value = .y))
      
      
      }
      
      } else {
        
        doParallel::registerDoParallel(cores = workers)
        
        stacked_draws <- foreach::foreach(i = 1:nrow(draws)) %dopar% {
          
          qgp <-   purrr::quietly(sraplus::get_posterior)
          
          out <- qgp(draws = draws$data[[i]], inits = inits, sra_data = sra_data,     model_name = model_name,
                     randos = randos,
                     knockout = knockout)
          
          pars <- names(out$result)
          
          flatstack <- purrr::map2_df(pars,out$result, ~data.frame(variable = .x, value = .y))
          
          
        }
        
        
      }
      draws$stack <- stacked_draws
      # d <- Sys.time() - a
      
      #       draws <- draws %>%
      #   dplyr::mutate(
      #     pars = purrr_map(
      #       data,
      #       purrr::quietly(get_posterior),
      #       inits = inits,
      #       sra_data = sra_data,
      #       model_name = model_name,
      #       randos = randos,
      #       knockout = knockout
      #     )
      #   ) %>%
      #   dplyr::select(-data)
      # 
      # draws$pars <- purrr::map(draws$pars, "result")
      # 
      
      draws <- draws %>%
        # dplyr::mutate(stack = purrr::map(pars, stack_stan)) %>%
        # dplyr::select(-pars) %>%
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
      
      
      # logs <- draws %>%
      #   dplyr::ungroup() %>%
      #   dplyr::filter(stringr::str_detect(variable, "log_")) %>%
      #   dplyr::mutate(value = exp(value)) %>%
      #   dplyr::mutate(variable = stringr::str_remove_all(variable, "log_"))
      #
      # draws <- draws %>%
      #   dplyr::filter(!stringr::str_detect(variable, "log_")) %>%
      #   dplyr::bind_rows(logs)
      
      out <- draws %>%
        dplyr::group_by(variable, year) %>%
        dplyr::summarise(
          mean = mean(value),
          sd = sd(value),
          lower = quantile(value, (1 - ci) / 2),
          upper = quantile(value, 1 - (1 - ci) / 2)
        ) %>%
        dplyr::ungroup()
      
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
        dplyr::ungroup()
      
      
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