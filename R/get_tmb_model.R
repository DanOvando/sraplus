#' Compile and load TMB models for srapls
#'
#' @param model_name the name of the sraplus model to load
#'
#' @export
#'
get_tmb_model <- function(model_name = "sraplus_tmb" ){

invisible(TMB::compile(system.file("tmb", paste0(model_name,".cpp"), package = "sraplus")))

invisible(dyn.load(TMB::dynlib(
  gsub(".cpp","",system.file("tmb", paste0(model_name,".cpp"), package = "sraplus")))))

}
