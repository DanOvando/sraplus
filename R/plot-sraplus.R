#' Plot timeseries of estimates from multiple sraplus models
#'
#' @param ... any number of sraplus models
#' @param fontsize font size for plots
#' @param years the years that the model fit applies to
#' @param plotvars variables to plot
#'
#' @return a ggplot object
#' @export
#'
plot_sraplus <- function(..., fontsize = 14,
                         years = NA,
                         plotvars = c("b_div_bmsy","c_div_msy","depletion","u_div_umsy")
){
  
  fit_names<- names(list(...))
  
  if (is.null(fit_names)){
    
    fit_names <- letters[seq_along(list(...))]
    
  }
  
  fits <- list(...)  %>% 
    purrr::set_names(fit_names) 
  
  fitframe <- dplyr::tibble(fit = fit_names,
                            temp = fits) %>% 
    dplyr::mutate(results = purrr::map(temp,"results")) %>% 
    dplyr::select(-temp) %>% 
    tidyr::unnest(cols = results)
  
  facet_labeller <- c(
    b_div_bmsy = "B/Bmsy",
    c_div_msy = "Catch/MSY",
    depletion = "Depletion",
    u_div_umsy = "U/Umsy"
  )
  
  fitframe %>% 
   dplyr::filter(variable %in% plotvars) %>% 
   dplyr::group_by(variable,fit) %>% {
   # dplyr::mutate(year = seq_along(mean)) %>% {
     if (!all(is.na(years))){
       dplyr::mutate(., year = years)
     } else {
       .
     }
   } %>% 
    dplyr::ungroup() %>% 
    ggplot2::ggplot() + 
   ggplot2::geom_ribbon(aes(year, ymin = lower, ymax = upper, fill = fit),
                size = 0.5, alpha = 0.5) +
   ggplot2::geom_line(aes(year, mean, color = fit),
              size = 1) +
   ggplot2::facet_wrap(~variable, scales = "free_y",
                       labeller = ggplot2::labeller(variable = facet_labeller)) + 
    sraplus::theme_sraplus(base_size = fontsize) + 
   ggplot2::scale_y_continuous( name = "", limits = c(0,NA)) +
    ggplot2::labs(x = "Year") +
   ggplot2::scale_fill_discrete(name = "Fit") + 
   ggplot2::scale_color_discrete(name = "Fit")
  
}