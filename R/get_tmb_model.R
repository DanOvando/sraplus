#' Compile and load TMB models for sraplus
#'
#' @param model_name the name of the sraplus model to load
#'
#' @export
#'
get_tmb_model <- function(model_name = "sraplus_tmb") {
  model_path <-
    file.path(getwd(),
              paste(model_name, utils::packageVersion("sraplus"), sep = '_v'))
  
  # if (model_name %in% names(getLoadedDLLs())) {
  #   #if DLL is already loaded unload
  #   
  #   dyn.unload(TMB::dynlib(file.path(model_path, model_name)))
  #   
  # }
  # 
  if (!dir.exists(model_path)) {
    dir.create(model_path)
    
    file.copy(
      from = system.file("tmb", paste0(model_name, ".cpp"), package = "sraplus"),
      to =  model_path,
      overwrite = FALSE
    )
    
    TMB::compile(file.path(model_path, paste0(model_name, ".cpp")))
    
  }
  
  if (!model_name %in% names(getLoadedDLLs())) {
    # if DLL is already loaded unload
    # dyn.unload(TMB::dynlib(file.path(model_path, model_name)))
    dyn.load(TMB::dynlib(file.path(model_path, model_name)))
    
  }
  
}
