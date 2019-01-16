plot_driors <- function(driors) {
  timeseries <- dplyr::tibble(year = driors$years,
                              catch = driors$catch)
  
  if (any(!is.na(driors$index))) {
    index_frame <- dplyr::tibble(year = driors$index_years,
                                 index = driors$index)
    timeseries <- timeseries %>%
      dplyr::left_join(index_frame, by = "year")
    
  }
  
  if (any(!is.na(driors$u_v_umsy))) {
    u_frame <- dplyr::tibble(year = driors$u_v_umsy,
                             u_umsy = driors$u_years)
    timeseries <- timeseries %>%
      dplyr::left_join(u_frame, by = "year")
    
  }
  
  if (any(!is.na(driors$effort))) {
    e_frame <- dplyr::tibble(year = driors$effort,
                             effort = driors$effort_years)
    timeseries <- timeseries %>%
      dplyr::left_join(e_frame, by = "year")
    
  }
  
  
  
  
  
}