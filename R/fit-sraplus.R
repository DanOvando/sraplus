#' Run sraplus
#'
#' @param driors a list of driors passed from sraplus::format_driors
#' @param include_fit logical indicating whether to return the fitted object
#' @param seed seed for model runs
#' @param plim minimum for hockey stick in pt model
#' @param model the name of the sraplus TMB version to be run
#' @param fit_catches logicatl indicating whether catches should be fit or passed
#' @param randos random effects when passing to TMB
#' @param draws the number of SIR samples to run
#' @param n_keep the number of SIR samples to keep
#' @param use_sir logical indicating whether to use SIR or TMB
#'
#' @return a fitted sraplus object
#' @export
#'
#' @examples
#' #' \dontrun{
#' sum("a")
#' }
fit_sraplus <- function(driors,
                        include_fit = TRUE,
                        seed = 42,
                        plim = 0.2,
                        model = "sraplus_v2",
                        fit_catches = TRUE,
                        randos = "uc_proc_errors",
                        draws = 1e6,
                        n_keep = 2000,
                        use_sir = FALSE) {
  
  
  knockout <- list()
  
  index_years <- driors$index_years
  
  effort_years <- driors$effort_years
  
  time <- length(driors$catch)
  
  sra_data <- list(
    catch_t = driors$catch,
    index_t = driors$index,
    index_years = driors$index_years,
    log_r_prior = log(driors$growth_rate),
    log_r_cv = driors$growth_rate_cv,
    log_init_dep_prior = log(driors$initial_b),
    log_init_dep_cv = driors$initial_b_cv,
    log_final_dep_prior = log(driors$terminal_b),
    log_final_dep_cv = driors$terminal_b_cv,
    time = time,
    fit_index = as.numeric(!all(is.na(driors$index))),
    use_u_prior = as.numeric(!all(is.na(driors$u_v_umsy))),
    u_years = driors$u_years,
    u_priors = driors$u_v_umsy,
    u_cv = driors$u_cv,
    plim = plim,
    sigma_proc_prior = driors$sigma_r/2,
    sigma_proc_prior_cv = driors$sigma_r_cv,
    ref_type = ifelse(driors$ref_type == "k", 0,1),
    use_final = !is.na(driors$terminal_b),
    use_final_u = as.numeric(!all(is.na(driors$log_final_u))),
    log_final_u = driors$log_final_u,
    log_final_u_cv = driors$log_final_u_cv,
    use_init =  !is.na(driors$initial_b),
    sigma_u = driors$u_cv,
    log_k_guess = log(10 * max(driors$catch)),
    f_cv = driors$f_cv
  )
  
  k_guess <- log(10* max(driors$catch))
  
  inits <- list(
    log_k = log(10* max(driors$catch)),
    log_r = log(driors$growth_rate),
    log_q = log(1e-3),
    log_sigma_obs = log(0.2),
    log_init_dep = log(1),
    log_sigma_proc = log(0.01),
    uc_proc_errors = rep(0, time - 1),
    log_m = log(2))
  
  if ((sra_data$fit_index == 0 & sra_data$use_u_prior == 0) | use_sir == TRUE) {
    
    sra_fit <- sraplus::sraplus(
      catches = sra_data$catch_t,
      r = pmax(0.01,rnorm(draws,driors$growth_rate,driors$growth_rate_cv)),
      m = runif(draws, 0.2, 6),
      init_dep = exp(rnorm(draws, sra_data$log_init_dep_prior, sra_data$log_init_dep_cv)),
      k = runif(draws, 1.15 * max(sra_data$catch_t), 50 * max(sra_data$catch_t)),
      sigma_procs = runif(draws, 0, 0.15),
      draws = draws,
      log_final_ref = ifelse(is.na(sra_data$log_final_dep_prior),0.5,sra_data$log_final_dep_prior),
      sigma_dep = ifelse(is.na(sra_data$log_final_dep_prior),1,sra_data$log_final_dep_cv),
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
    
    outs <- stringr::str_detect(names(sra_fit),"_t")
    
    tidy_fits <-
      purrr::map_df(
        purrr::keep(sra_fit, outs),
        ~ dplyr::tibble(.x[, keepers]) %>% dplyr::mutate(year = 1:nrow(.)) %>% tidyr::gather(draw, value, -year),
        keepers = keepers,
        .id = "variable"
      ) %>%
      dplyr::mutate(draw = stringr::str_replace_all(draw,"\\D","") %>% as.numeric())
    
    out <- tidy_fits %>%
      dplyr::group_by(year, variable) %>%
      dplyr::summarise(mean = mean(value),
                sd = sd(value),
                lower = quantile(value, 0.1),
                upper = quantile(value, 0.9)) %>%
      dplyr::ungroup()
    
    
    out <- list(results = out,
                fit = tidy_fits)
    
   
    # index_tests <- sra_fit$index_hat[, sra_fit$keepers]

  }
  
  else {
    
    if (sra_data$use_u_prior == 0){
      
      udist <- NULL
      
    } else {
      
      udist <- 2
      
    }
    
    # pen <-
    #   list(
    #     bstatus.mean =  sra_data$log_final_dep_prior * 2,
    #     bstatus.sd = sra_data$log_final_dep_cv * 2,
    #     bstatus.dist = 2,
    #     ustatus.mean = log(last(sra_data$u_priors)),
    #     ustatus.sd = sra_data$u_cv,
    #     ustatus.dist = udist ,
    #     initial.mean = sra_data$log_init_dep_prior * 2 * 2,
    #     initial.sd = sra_data$log_init_dep_cv,
    #     initial.dist = 2
    #   )
    
    # driors$scientific_name <- Taxon <- c(Class="Actinopterygii", Order="Perciformes",
    #            Family="Scombridae", Genus="Thunnus", Species="albacares")
    #
    # start <- Sys.time()
    # fit <- sraplus::run.SIR(nrep=20000, Catch=sra_data$catch_t, Taxon = driors$scientific_name, penalties=pen,
    #                         years=seq_along(sra_data$catch_t))
    #
    # old_priors <- pen
    #
    # old_priors$ustatus.mean <- old_priors$fstatus.mean
    #
    # old_priors$ustatus.sd <- old_priors$fstatus.sd
    #
    # old_priors$ustatus.dist <- old_priors$fstatus.dist
    #
    # old_priors$fstatus.dist <- NULL
    #
    # start <- Sys.time()
    #
    # set.seed(42)
    # oldfit <-
    #   colesraplus::run.SIR(
    #     nrep = 500000,
    #     Catch = sra_data$catch_t,
    #     Taxon = driors$scientific_name,
    #     penalties = old_priors,
    #     years = seq_along(sra_data$catch_t)
    #   )
    # Sys.time() - start
    #
    # sraplus::plot_fit( oldfit)
    
    
    
    
    
    lower_k <- log(1.25 * max(driors$catch))
    
    upper_k <- log(50 * max(driors$catch))
    
    if (fit_catches == TRUE){
      
      inits$inv_f_t = rep(-2, time - 1);
      
    }
    
    
    if (sra_data$fit_index == 0){
      
      knockout$log_q <- NA
      
      knockout$log_sigma_obs <- NA
      
    }
    
    # knockout$log_init_dep = NA
    
    if (sra_data$fit_index == 0 & sra_data$use_u_prior == 0){
      
      knockout$log_sigma_proc <- NA
      
      knockout$uc_proc_errors <- rep(NA, time - 1)
      
      inits$log_sigma_proc <- -1e6
      
      inits$uc_proc_errors <- rep(0, time - 1)
      
      randos <- NULL
      
      randos <- "inv_f_t"
      
    }
    
    knockout <- purrr::map(knockout, as.factor)
    
    sra_model <-
      TMB::MakeADFun(
        data = sra_data,
        parameters = inits,
        DLL = model,
        random = randos,
        silent = TRUE,
        inner.control=list(maxit=1e3),
        hessian=FALSE,
        map = knockout
      )
    
    lower = rep(-Inf,length(sra_model$par)) %>%
      set_names(names(sra_model$par))
    
    lower['log_k'] <- lower_k
    
    upper = rep(Inf,length(sra_model$par)) %>%
      set_names(names(sra_model$par))
    
    upper['log_k'] <- upper_k
    
    upper["log_init_dep"] <- log(1.5)
    
    set.seed(seed)
    
    fit <- Optimize(
      sra_model,
      fn = sra_model$fn,
      gr = sra_model$gr,
      newtonsteps = 3,
      lower = lower,
      upper = upper,
      getsd = FALSE,
      control = list(
        eval.max = 1e3,
        iter.max = 1e3,
        rel.tol = 1e-10
      )
    )
    
    
    if (fit$max_gradient > 1e-3){
      
      fit <- Optimize(
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
    
    
    # wtf <- fit_save$report()
    #
    # wtf$wtf
    #
    # wtf$dep_t
    #
    # fit_report$value["log_dep"]
    #
    # foo <- out %>% filter(variable == "log_dep")
    #
    # plot(wtf$dep_t, exp(foo$mean))
    
    #
    # plot(wtf$dep_t, wtf$growth_t)
    #
    # plot(wtf$index_hat_
    
    fit_report <- fit_save$report()
    
    fit_sd_report <- sdreport(fit_save, bias.correct = TRUE)
    
    out <-
      data_frame(
        variable = names(fit_sd_report$value),
        mean = fit_sd_report$value,
        sd = fit_sd_report$sd
      )
    
    # SERIOUS hack. sdreport ignores if statements inside population model, so need to put the report values back in
    # true_names <- names(fit_report)
    #
    # for (i in true_names){
    #
    #   if (i %in% out$variable){
    #
    #     out$mean[out$variable == i] <- fit_report[[i]]
    #
    #     # out$mean[out$variable == i] %>% View()
    #   }
    #
    # } # close hack
    
    out <- out %>%
      mutate(lower = mean - 1.96 * sd,
             upper = mean + 1.96 * sd)
    
    if (include_fit == FALSE){
      fit = NA
    }
    
    out <- list(results = out,
                fit = fit)
    
    rm(sra_model)
    
    rm(fit)
    
  } # close else
  
  
  return(out)
  
  
}