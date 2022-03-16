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

## Creating a new environment to store info about state of Julia modules
.pkgglobalenv <- new.env(parent=emptyenv())
assign("hasModules", NULL, envir = .pkgglobalenv) ## logical value, TRUE if all required Julia modules are available

.onLoad <- function(libname, pkgname) {
  ## If Julia is not installed
  if (!findJuliaNoError(test = T)) {
    packageStartupMessage("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable.\n")
  }else{

    ## Test that the correct version of Julia is installed
    min_version <- "1.6.5"
    check_version <- checkJuliaVersion(min_version)

    ## If we can't check the version
    if(check_version == -1){
      packageStartupMessage("Error when checking the existing Julia version. Please check that Julia is correctly installed.")
    }

    ## If the version is too old
    else if(check_version == 0){
      packageStartupMessage("The current version of Julia is < to v", min_version, ". Please install Julia v", min_version, " or later.")
    }

    ## If the version of Julia is ok, check whether the required Julia modules are installed
    else {
      checkJuliaModules(silence_warning = TRUE)
      hasModules <- get("hasModules", envir = .pkgglobalenv)
      if(!hasModules) packageStartupMessage("Required Julia modules not installed. Please run 'installJuliaModules()' to install them.")
    }
  }
}

# .onAttach <- function(libname, pkgname) {
#   if (!findJuliaNoError(test = T)) {
#    warning("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable.\n")
#   }else{
#     ## Test that the correct version of Julia is installed
#     julia_bin = XRJulia::findJulia()
#     if (.Platform$OS.type == "windows") cmd_args = "-E \"VERSION < v\\\"1.7.1\\\"\"" else cmd_args = "-E 'VERSION < v\"1.7.1\"'"
#     version_check = tryCatch(system2(julia_bin, args = cmd_args, stdout = TRUE, stderr = TRUE), warning = function(e) "", error = function(e) "")
#     if(version_check == ""){
#       warning("Error when checking the existing Julia version. Please check that Julia is correctly installed.")
#     }
#     else if(base::grepl("true", version_check, ignore.case = TRUE)){
#       warning("The current version of Julia is < to v1.7.1. Please install Julia v1.7.1 or later.")
#     }
#     else{
#       testFile = system.file("julia", "testModules.jl", package = "sismonr")
#       if (.Platform$OS.type == "windows") cmd = paste0('"',julia_bin,'" ', testFile) else cmd = paste(julia_bin, " ", testFile)
#       packageStartupMessage("Checking that the required Julia modules are installed...")
#       hasModules = base::system(cmd, intern = FALSE, ignore.stdout = T, ignore.stderr = T)
#       packageStartupMessage("Done.")
#     }
#   }
# }
