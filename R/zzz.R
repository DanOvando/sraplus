.onLoad <- function(libname, pkgname) {
  
  # if (dir.exists("tmb")){
  #   
  #   
  # }
  
  # dyn.load(TMB::dynlib(
  #   gsub(".cpp","",system.file("sraplus", "baranov.cpp", package = "sraplus"))))
  # 
}

# clean up DLLs on package unload
.onUnload <- function(libpath){
  
  library.dynam.unload("sraplus", libpath)
  
}