#' plot_prior_posterior
#'
#' @param fit 
#' @param driors 
#' @param draws 
#' @param prob 
#'
#' @return
#' @export
#'
get_prior_posterior <- function(fit, driors,
                                 draws = 1000,
                                 prob  = 0.9){
  
  
 
  timeseries <- dplyr::tibble(year = driors$years,
                              catch = driors$catch)
  
  fitseries <- as.data.frame(matrix(NA, nrow = 0, ncol = 5))
  
  colnames(fitseries) <- c("metric","year", "mean", 'lower','upper')
  any_fits <- TRUE
  if (any(!is.na(driors$index))) {
    index_frame <- dplyr::tibble(year = driors$index_years,
                                 index = driors$index)
    timeseries <- timeseries %>%
      dplyr::left_join(index_frame, by = "year")
    
    index_fit <- fit$results %>%
      dplyr::filter(variable == "index_hat_t") %>%
      dplyr::select(mean, lower, upper) %>%
      dplyr::mutate(metric = "index",
                    year = driors$years) %>%
      dplyr::filter(year %in% driors$index_years) %>%
      dplyr::mutate(meanmean = mean(mean),
             sdmean = sd(mean)) %>%
      dplyr::mutate(
        mean = (mean - meanmean) / sdmean,
        lower = (lower - meanmean) / sdmean,
        upper = (upper - meanmean) / sdmean
      ) %>%
      dplyr::select(metric, year, mean, lower, upper)
    
    fitseries <- rbind(fitseries, index_fit)
    
  }
  
  if (any(!is.na(driors$u))) {
    u_frame <- dplyr::tibble(year = driors$u_years,
                             u = driors$u)
    timeseries <- timeseries %>%
      dplyr::left_join(u_frame, by = "year")
    
    u_fit <- fit$results %>% 
      dplyr::filter(variable == "u_div_umsy") %>% 
      dplyr::select(mean,lower, upper) %>% 
      dplyr::mutate(metric = "u",
                    year = driors$years) %>% 
      # dplyr::filter(year %in% driors$u_years) %>% 
      dplyr::mutate(
        meanmean = mean(mean),
        sdmean = sd(mean)) %>% 
      dplyr::mutate(mean = (mean - meanmean) / sdmean,
                    lower = (lower - meanmean) / sdmean,
                    upper = (upper - meanmean) / sdmean) %>% 
      dplyr::select(metric, year, mean, lower, upper)
    
    fitseries <- rbind(fitseries, u_fit)
    
  }
  
  if (any(!is.na(driors$effort))) {
    e_frame <- dplyr::tibble(year = driors$effort_years,
                             effort = driors$effort)
    timeseries <- timeseries %>%
      dplyr::left_join(e_frame, by = "year") %>% 
      dplyr::mutate(cpue = catch / effort)
    
    cpue_fit <- fit$results %>% 
      dplyr::filter(variable == "index_hat_t") %>% 
      dplyr::select(mean,lower, upper) %>% 
      dplyr::mutate(metric = "cpue",
                    year = driors$years) %>% 
      dplyr::filter(year %in% driors$effort_years) %>% 
      dplyr::mutate(
        meanmean = mean(mean),
        sdmean = sd(mean)) %>% 
      dplyr::mutate(mean = (mean - meanmean) / sdmean,
                    lower = (lower - meanmean) / sdmean,
                    upper = (upper - meanmean) / sdmean) %>% 
      dplyr::select(metric, year, mean, lower, upper)
    
    fitseries <- rbind(fitseries, cpue_fit)
    
    
  }
  
  if (nrow(fitseries) == 0){
    
    any_fits <- FALSE
    
    fitseries <- as.data.frame(matrix(NA, nrow = nrow(timeseries), ncol = 3))
    
    colnames(fitseries) <- c("mean", 'lower','upper')
    
  }
  
  timeseries <- timeseries %>%
    tidyr::gather(metric, value,-year) %>% 
    dplyr::group_by(metric) %>% 
    dplyr::mutate(value = scale(value)) %>% 
    dplyr:: ungroup() %>% {
      if (any_fits){
        dplyr::left_join(.,fitseries, by = c("year","metric"))
      } else {
        dplyr::bind_cols(.,fitseries)
      }
    } 
  
 
  multiple_u <- FALSE
  
  if (length(driors$log_terminal_u) > 1) {
    multiple_u <- TRUE
    
    driors$log_terminal_u1 = driors$log_terminal_u[1]
    
    driors$log_terminal_u1_cv = driors$log_terminal_u_cv[1]
    
    driors$log_terminal_u2 = driors$log_terminal_u[2]
    
    driors$log_terminal_u2_cv = driors$log_terminal_u_cv[2]
    
    driors <-
      purrr::list_modify(driors,
                         "log_terminal_u" = NULL,
                         "log_terminal_u_cv" = NULL)
    
    driors$terminal_u1 = driors$terminal_u[1]
    
    driors$terminal_u1_cv = driors$terminal_u_cv[1]
    
    driors$terminal_u2 = driors$terminal_u[2]
    
    driors$terminal_u2_cv = driors$terminal_u_cv[2]
    
    driors <-
      purrr::list_modify(driors,
                         "terminal_u" = NULL,
                         "terminal_u_cv" = NULL)
    
  }
  
  
  vars <- names(driors)
  
  var_count <- table(stringr::str_replace_all(vars, "_cv", ""))
  plot_vars <- names(var_count)[var_count == 2]
  
  has_values <-
    purrr::map_lgl(plot_vars, ~ !all(is.na(driors[[.x]])))
  
  plot_vars <- plot_vars[has_values]
  
  
  
  foo <- function(var, driors, n = 1000) {
    if (stringr::str_detect(var, "log")) {
      sims <- rnorm(n, driors[[var]], driors[[paste0(var, "_cv")]])
      
      simquants <-
        quantile(sims, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2))
      
      sims <- sims[sims >= simquants[1] & sims <= simquants[2]]
      
    } else {
      sims <-
        exp(rnorm(n, log(pmax(driors[[var]], 1e-6)), driors[[paste0(var, "_cv")]]))
      
      simquants <-
        quantile(sims, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2))
      
      sims <- sims[sims >= simquants[1] & sims <= simquants[2]]
      
      
    }
  }
  
  priors <- dplyr::tibble(variable = plot_vars) %>%
    dplyr::mutate(draws = purrr::map(variable, foo, driors = driors)) %>%
    tidyr::unnest(cols = draws)
  
  # this is going to be very unelegant given differnces in naming
  
  fits <- fit$results %>% 
    dplyr::mutate(variable = dplyr::case_when(variable == "r" ~ "growth_rate",
                                              variable == "m" ~ "shape",
                                              TRUE ~ variable ))
  
  priors$variable <- stringr::str_replace_all(priors$variable,"_prior","") 
  
  priors <- priors %>% 
    dplyr::group_by(variable) %>% 
    dplyr::summarise(mean = mean(draws),
                     lower = quantile(draws,probs = c((1 - prob) / 2)),
                     upper = quantile(draws,probs = c(1 - (1 - prob) / 2))) %>% 
    dplyr::mutate(source = "Prior")

  if ("initial_state" %in% priors$variable) {
    
    if (driors$b_ref_type == "b") {
      
      initial_state <- fits[fits$variable == "b_div_bmsy",]
      
      initial_state <- as.data.frame(initial_state[1,])
      
      initial_state$variable <-  "initial_state"
      
      fits <- rbind(fits, (initial_state))
      
      
    } else if (driors$b_ref_type == "k"){
      
      initial_state <- fits[fits$variable == "depletion",]
      
      initial_state <- as.data.frame(initial_state[1,])
      
      initial_state$variable <-  "initial_state"
      
      fits <- rbind(fits, (initial_state))
      
      
    }
    
    
  }
  
  if ("terminal_state" %in% priors$variable){
    
    if (driors$b_ref_type == "b") {
      
      terminal_state <- fits[fits$variable == "b_div_bmsy",]
      
      terminal_state <- as.data.frame(terminal_state[nrow(terminal_state),])
      
      terminal_state$variable <-  "terminal_state"
      
      
      fits <- rbind(fits, (terminal_state))
      
      
    } else if (driors$b_ref_type == "k"){
      
      terminal_state <- fits[fits$variable == "depletion",]
      
      terminal_state <- as.data.frame(terminal_state[nrow(terminal_state),])
      
      terminal_state$variable <-  "terminal_state"
      
      fits <- rbind(fits, terminal_state)
      
      
    }
    
    
  } # close if terminal state
  
  if (any(stringr::str_detect(priors$variable, "terminal_u"))){
    
    if (driors$f_ref_type == "fmsy") {
      
      terminal_u <- fits[fits$variable == "u_div_umsy",]
      
      terminal_u <- as.data.frame(terminal_u[nrow(terminal_u),])
      
      terminal_u$variable <-  "terminal_u"
      
      if (multiple_u){
        
        terminal_u1 <- fits[fits$variable == "u_div_umsy",]
        
        terminal_u1 <- as.data.frame(terminal_u[nrow(terminal_u),])
        
        terminal_u1$variable <-  "terminal_u1"
        
        terminal_u2 <- fits[fits$variable == "u_div_umsy",]
        
        terminal_u2 <- as.data.frame(terminal_u[nrow(terminal_u),])
        
        terminal_u2$variable <-  "terminal_u2"
        
        terminal_u <- dplyr::bind_rows(terminal_u1 , terminal_u2)
        
      }
      
      fits <- rbind(fits, (terminal_u))
      
      
    } else if (driors$f_ref_type == "f"){
      
      terminal_u <- fits[fits$variable == "u_t",]
      
      terminal_u <- as.data.frame(terminal_u[nrow(terminal_u),])
      
      terminal_u$variable <-  "terminal_u"
      
      if (multiple_u){
        
        terminal_u1 <- fits[fits$variable == "u_t",]
        
        terminal_u1 <- as.data.frame(terminal_u[nrow(terminal_u),])
        
        terminal_u1$variable <-  "terminal_u1"
        
        terminal_u2 <- fits[fits$variable == "u_t",]
        
        terminal_u2 <- as.data.frame(terminal_u[nrow(terminal_u),])
        
        terminal_u2$variable <-  "terminal_u2"
        
        terminal_u <- dplyr::bind_rows(terminal_u1 , terminal_u2)
        
      }
      
      fits <- rbind(fits, terminal_u)
    }
    
    
  } # close if terminal state
  
  
posteriors <- fits %>% 
    dplyr::filter(variable %in% priors$variable) %>% 
    dplyr::mutate(source = "Posterior")
  
  prior_posterior <- posteriors[,colnames(priors)] %>% 
    rbind(priors %>% dplyr::filter(variable %in% posteriors$variable))
  
  
  return(list(prior_posterior = prior_posterior,
              fits = timeseries,
              any_fits = any_fits))  
  
  
}