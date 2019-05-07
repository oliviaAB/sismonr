# .onLoad <- function(libname, pkgname) {
#   if (!XRJulia::findJulia(test = T)) {
#     warning("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable.\n")
#   }
#   testFile = system.file("julia", "testModules.jl", package = "sismonr")
#   julia_bin = XRJulia::findJulia()
#   if (.Platform$OS.type == "windows") cmd = paste0('"',julia_bin,'" ', testFile) else cmd = paste(julia_bin, "<", testFile)
#   hasModules = base::system(cmd, intern = FALSE)
#   # XRJulia::RJulia()  # creates a Julia interface
#   # print(XRJulia::RJulia())  # temporary
#   # XRJulia::juliaUsing("ClobberingReload")
#   # XRJulia::juliaCommand(paste0("sinclude(\"", system.file("julia", "sismonr.jl",
#   #                                                         package = "sismonr"), "\")"))  # load the julia functions of the package in the new evaluator
# }

.onAttach <- function(libname, pkgname) {
  if (!findJuliaNoError(test = T)) {
   warning("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable.\n")
  }else{
    ## Test that the correct version of Julia is installed
    julia_bin = XRJulia::findJulia()
    if (.Platform$OS.type == "windows") cmd = paste0('"', julia_bin, '" ', "-E \"VERSION < v\\\"1.0.0\\\"\"") else cmd = paste0(julia_bin, " ", "-E 'VERSION < v\"1.0.0\"'")
    if(base::system(cmd, intern = T, ignore.stdout = F)){
      warning("The current version of Julia is < to v1.0. Please install Julia v1.0 or later.")
    }else{
      testFile = system.file("julia", "testModules.jl", package = "sismonr")
      if (.Platform$OS.type == "windows") cmd = paste0('"',julia_bin,'" ', testFile) else cmd = paste(julia_bin, " ", testFile)
      packageStartupMessage("Checking if the required Julia modules are installed...")
      hasModules = base::system(cmd, intern = FALSE, ignore.stdout = T)
      packageStartupMessage("Done.")
    }
  }
}
