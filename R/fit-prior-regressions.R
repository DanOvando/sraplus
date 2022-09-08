#' fir-sraplus-regressions
#'
#' @param split 
#' @param model_structure 
#' @param use_splits 
#' @param chains 
#' @param cores 
#' @param iter 
#' @param adapt_delta 
#' @param max_treedepth 
#' @param produce 
#' @param refresh 
#' @param post_draws 
#' @param ... 
#'
#' @return list with fitted object
#' @export
#'
fit_prior_regressions <-
  function(split,
           model_structure,
           use_splits = TRUE,
           chains = 4,
           cores = 4,
           iter = 2000,
           adapt_delta = adapt_delta,
           max_treedepth = max_treedepth,
           produce = "summary",
           refresh = 0,
           post_draws = 10,
           ...) {
    args <- list(...)
    
    
    # chains = 4
    #
    # cores = 4
    
    # ram_v_fmi <- ram_v_fmi %>%
    #   mutate(log_value = log(value))
    #
    #
    # training <- ram_v_fmi %>%
    #   filter(metric == "mean_b") %>%
    #   slice(1:100)
    #
    #
    # testing <- ram_v_fmi %>%
    #   filter(metric == "mean_b") %>%
    #   slice(101:173)
    #
    # model_structure <-
    #   "log_value ~ research + management + enforcement + socioeconomics +(1|country_rfmo)"
    
    if (use_splits == TRUE){
      
      training <- rsample::training(split)
      
      testing <- rsample::testing(split)
      
    } else {
      
      training <- split
      
      testing <- split
    }
    
    # has_nas <- map_lgl(training, ~any(is.na(.x))) | map_lgl(testing, ~any(is.na(.x)))
    
    
    model <- as.formula(model_structure)
    
    # if (str_detect(model_structure,"log")){
    #   
    #   training[,c("research","management","socioeconomics","enforcement")] <- training[,c("research","management","socioeconomics","enforcement")] + 1e-6
    #   
    #   testing[,c("research","management","socioeconomics","enforcement")] <- testing[,c("research","management","socioeconomics","enforcement")] + 1e-6
    #   
    #   
    # }
    
    if (str_detect(model_structure, "\\|")) {
      fit <-
        stan_glmer(
          model,
          data = training,
          chains = chains,
          cores = cores,
          refresh = refresh,
          iter = iter
        )
    } else if ((str_detect(model_structure, "s\\("))) {
      
      fit <-
        stan_gamm4(
          model,
          data = training,
          chains = chains,
          cores = cores,
          refresh = refresh,
          iter = iter
        )
      
    } else {
      fit <-
        stan_glm(
          model,
          data = training,
          chains = chains,
          cores = cores,
          refresh = refresh,
          iter = iter
        )
      
    }
    
    training_draws <-
      tidybayes::add_fitted_draws(model = fit,
                                  newdata = training,
                                  value = "fit_pred",
                                  n = post_draws) %>%
      mutate(type = "fitted") %>%
      ungroup()
    
    pp_training_draws <-
      tidybayes::add_predicted_draws(model = fit, newdata = training,n = post_draws) %>%
      ungroup()
    
    training_draws <- training_draws %>%
      mutate(pp_pred = pp_training_draws$.prediction) %>%
      gather(metric, prediction, fit_pred, pp_pred)
    
    pp_testing_draws <-
      tidybayes::add_predicted_draws(model = fit, newdata = testing, n = post_draws) %>%
      ungroup()
    # training_draws %>%
    #   ggplot(aes(log_value, prediction, color = metric)) +
    #   geom_point(alpha = 0.25)
    #
    #   ggplot(aes(factor(.row), prediction, fill = metric)) +
    #   geom_violin()
    
    if (produce == "summary") {
      training_summary <- training_draws %>%
        group_by(.row, metric) %>%
        summarise(mean_pred = mean(prediction),
                  observed = mean(log_value)) %>%
        spread(metric, mean_pred) %>%
        ungroup()
      
      testing_summary <- pp_testing_draws %>%
        group_by(.row) %>%
        summarise(pp_pred = mean(.prediction),
                  observed = mean(log_value)) %>%
        ungroup()
      
      out <- list(
        training_summary = training_summary,
        testing_summary = testing_summary,
        fit = fit
      )
      
    } else if (produce == "results") {
      out <- list(training_draws = training_draws,
                  pp_testing_draws = pp_testing_draws,
                  fit = fit)
      
    }
    
    return(out)
    
  }