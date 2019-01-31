#' stack reports from tmbstan fits into cohesive dataframe
#'
#' @param stack a stack of adreport outputs from one set of stan parameters
#'
#' @return a stacked set of parameters
#' @export
#'
stack_stan <- function(stack) {
  
  # stack <- draws$pars[[1]]
  
  pars <- names(stack)
  
  flatstack <- purrr::map2_df(pars,stack, ~dplyr::tibble(variable = .x, value = .y))
  
}