.onLoad <- function(libname, pkgname) {
  if (!XRJulia::findJulia(test = T)) {
    warning("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable.\n")
  }
  # XRJulia::RJulia()  # creates a Julia interface
  # print(XRJulia::RJulia())  # temporary
  # XRJulia::juliaUsing("ClobberingReload")
  # XRJulia::juliaCommand(paste0("sinclude(\"", system.file("julia", "sismonr.jl",
  #                                                         package = "sismonr"), "\")"))  # load the julia functions of the package in the new evaluator
}
