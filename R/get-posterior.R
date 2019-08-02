#' get posterior draws from tmbstan fit object
#'
#' @param draws parameter draws
#' @param inits initi values for default model
#' @param sra_data data
#' @param model model name
#' @param randos random effects
#' @param knockout parameter to knock out
#'
#' @return tibble of posterior draws
#' @export
#'
get_posterior <-  function(draws,
                           inits,
                           sra_data,
                           model,
                           randos,
                           knockout) {
  # draws <- a$data[[1]]
  
  draw_names <- names(draws)
  
  draw_names <-
    stringr::str_remove_all(draw_names, "\\d|__|\\[|\\]")
  
  draw_locs <- which(draw_names %in% names(inits))
  
  par_list <-
    dplyr::tibble(par = draw_names[draw_locs], vals = as.numeric(draws[draw_locs])) %>%
    tidyr::nest(-par)
  
  
  # SERIOUSLY annoying step to put back in things that you knocked out with map
  if (length(knockout) >= 1) {
    for (i in 1:length(knockout)) {
      temp <- names(knockout)[i]
      
      temp_entry <-
        dplyr::tibble(par = temp, data = list(dplyr::tibble(vals = inits[temp][[1]])))
      
      par_list <- par_list %>%
        dplyr::bind_rows(temp_entry)
      
    }
    
  }
  
  pars <- purrr::map(par_list$data, ~ dplyr::pull(.x[, 1])) %>%
    purrr::set_names(par_list$par)
  
  temp <-
    TMB::MakeADFun(
      data = sra_data,
      parameters = pars,
      DLL = model,
      random = randos,
      silent = TRUE,
      inner.control = list(maxit = 1e3),
      hessian = FALSE,
      map = knockout
    )
  
  full_pars <- temp$report()
  
  
}