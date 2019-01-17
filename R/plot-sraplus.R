#' Plot timeseries of estimates from multiple sraplus models
#'
#' @param ... any number of sraplus models
#'
#' @return a ggplot object
#' @export
#'
plot_sraplus <- function(...){
  
  fit_names<- names(list(...))
  
  fits <- list(...)  %>% 
    purrr::set_names(fit_names) 
  
  fitframe <- dplyr::tibble(fit = fit_names,
                            temp = fits) %>% 
    dplyr::mutate(results = purrr::map(temp,"results")) %>% 
    dplyr::select(-temp) %>% 
    tidyr::unnest()
  
  plotvars <- c("b_div_bmsy","c_div_msy","depletion","u_div_umsy")
  
 fitframe %>% 
   dplyr::filter(variable %in% plotvars) %>% 
   dplyr::group_by(variable,fit) %>% 
   dplyr::mutate(year = seq_along(mean)) %>% 
    dplyr::ungroup() %>% 
    ggplot() + 
    geom_ribbon(aes(year, ymin = lower, ymax = upper, fill = fit),
                size = 0.5, alpha = 0.5) +
    geom_line(aes(year, mean, color = fit),
              size = 1) +
    facet_wrap(~variable, scales = "free_y") + 
    sraplus::theme_sraplus() + 
    scale_y_continuous( name = "") +
    labs(x = "Time") +
    scale_fill_discrete(name = "Fit") + 
    scale_color_discrete(name = "Fit")
  
}