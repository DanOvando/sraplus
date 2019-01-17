#' Compile and load TMB models for sraplus
#'
#' @param model_name the name of the sraplus model to load
#'
#' @export
#'
get_tmb_model <- function(model_name = "sraplus_tmb" ){

TMB::compile(system.file("tmb", paste0(model_name,".cpp"), package = "sraplus"))

dyn.load(TMB::dynlib(
  gsub(".cpp","",system.file("tmb", paste0(model_name,".cpp"), package = "sraplus"))))

}
