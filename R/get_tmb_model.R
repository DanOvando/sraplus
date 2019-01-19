#' Compile and load TMB models for sraplus
#'
#' @param model_name the name of the sraplus model to load
#'
#' @export
#'
get_tmb_model <- function(model_name = "sraplus_tmb") {
  if (!dir.exists(file.path(getwd(), "tmb"))) {
    dir.create("tmb")
    
  }
  
  file.copy(from = system.file("tmb", paste0(model_name, ".cpp"), package = "sraplus"),
            to =   file.path(getwd(), "tmb"),
            overwrite = FALSE)
  
  TMB::compile(file.path("tmb", paste0(model_name, ".cpp")))
  
  dyn.load(TMB::dynlib(file.path("tmb", model_name)))
  # gsub(".cpp","",system.file("tmb", paste0(model_name,".cpp"), package = "sraplus"))))
  
}
