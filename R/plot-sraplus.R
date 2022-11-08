#' Plot timeseries of estimates from multiple sraplus models
#'
#' @param ... any number of sraplus models
#' @param fontsize font size for plots
#' @param years replacement vector of years for plotting if desired
#' @param plotvars variables to plot
#' @param kobe FALSE to return ribbons plot. TRUE to plot a Kobe plot
#' @param max_kobe_val the max value for the x and y axes of the Kobe plot
#'
#' @return a ggplot object
#' @export
#'
plot_sraplus <- function(..., fontsize = 14,
                         years = NA,
                         plotvars = c("b_div_bmsy","c_div_msy","depletion","u_div_umsy"),
                         kobe = FALSE,
                         max_kobe_val = 4
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
  
  if (kobe == FALSE){
  fitframe %>% 
   dplyr::filter(variable %in% plotvars) %>% 
   dplyr::group_by(variable,fit) %>% {
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
  } else {
    
    results <- fitframe %>% 
      dplyr::filter(variable %in% c("b_div_bmsy","u_div_umsy")) %>% 
      dplyr::group_by(variable,fit) %>% {
        if (!all(is.na(years))){
          dplyr::mutate(., year = years)
        } else {
          .
        }
      } %>% 
      dplyr::ungroup() 
    
    
    max_val <- max_kobe_val
    
    points <- results %>%
      dplyr::filter(variable %in% c("b_div_bmsy","u_div_umsy")) %>%
      dplyr::ungroup() %>%
      dplyr::select(year, fit,variable, mean) %>%
      dplyr::mutate(mean = pmin(mean, max_val)) %>%
      tidyr::pivot_wider(names_from = variable, values_from = mean)
    
    
    segments <- results %>%
      dplyr::filter(variable %in% c("b_div_bmsy","u_div_umsy")) %>%
      dplyr::ungroup() %>%
      dplyr::select(-sd) %>%
      tidyr::pivot_longer(c(lower, upper), names_to = "direction", values_to = "endpoint") %>%
      dplyr::mutate(
        radius = pmin(endpoint, max_val),
        mean = pmin(mean, max_val),
        angle = dplyr::case_when(
          variable == "b_div_bmsy" & direction == "upper" ~  0,
          variable == "b_div_bmsy" & direction == "lower" ~  pi,
          variable == "u_div_umsy" & direction == "upper" ~ 0.5 * pi,
          variable == "u_div_umsy" & direction == "lower" ~ -0.5 * pi
        ),
        radius = abs(endpoint - mean)) %>%
      dplyr::left_join(points, by = c("year", "fit"))
    
    
    segments %>%
      ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 1, linetype = 2) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2) +
      ggplot2::geom_spoke(aes(b_div_bmsy, u_div_umsy, angle = angle, radius = radius, color = year)) +
      ggplot2::geom_point(aes(b_div_bmsy, u_div_umsy, fill = year), size = 3, shape = 21) +
      ggplot2::geom_point(data = segments %>% dplyr::filter(year == max(year)),aes(b_div_bmsy, u_div_umsy, fill = year), size = 3, shape = 21, fill = "red") +
      ggplot2::facet_wrap(~fit) +
      ggplot2::scale_fill_viridis_c(name = "Year", option = "plasma", direction = -1) +
      ggplot2::scale_color_viridis_c(name = "Year",option = "plasma", direction = -1) +
      ggplot2::scale_x_continuous(name = "B/BMSY", limits = c(0, max_val), expand = ggplot2::expansion(mult = c(0, 0.05)),
                                  breaks = seq(0,max_val, by = 0.5)) +
      ggplot2::scale_y_continuous(name = "F/FMSY", limits = c(0, max_val), expand = ggplot2::expansion(mult = c(0, 0.05)),
                                  breaks = seq(0.5,max_val, by = 0.5)) +
      sraplus::theme_sraplus() 
  
    
    
    
    
  }
  
}