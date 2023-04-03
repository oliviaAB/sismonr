#' Find a Julia executable
#'
#' The function is almost identical to the \code{findJulia} function from the \code{XRJulia} package,
#' but prevent errors arising from looking for a Julia executable with the "which"/"where" function if Julia
#' is not present.
#'
#' @param test If TRUE, returns TRUE/FALSE depending on whether or not a Julia executable is found. If FALSE, returns the
#' path to the Julia executable if it exists.
#' @return TRUE/FALSE or the path to the Julia executable
#' @examples
#' \dontrun{
#' findJuliaNoError(test = T)
#' }
#' @export
findJuliaNoError = function(test = FALSE) {
  ## See if a location for the Julia executable has been specified
  ## environment variables JULIA_BIN or JULIA_SRC
  envvar <- Sys.getenv("JULIA_BIN")
  if(!nzchar(envvar)) {
    src <- Sys.getenv("JULIA_SRC")
    if(nzchar(src)) {
      ## try either the standard Julia source, or the directory above bin, etc.
      trybin <- paste(src, "usr","bin","julia", sep=.Platform$file.sep)
      if(file.exists(trybin))
        envvar <- trybin
      else {
        trybin <- paste(src, "julia", sep=.Platform$file.sep)
        if(file.exists(trybin))
          envvar <- trybin
      }
    }
  } # if none of these succeeds, `which julia` used to find an executable version
  if(!nzchar(envvar)) {
    command <-if (.Platform$OS.type == "windows") "where" else "which"
    envvar <- tryCatch(system2(command, "julia", stdout = TRUE, stderr = TRUE), warning = function(e) "")
    if(test)
      Sys.setenv(JULIA_BIN = envvar) # so next call finds this immediately
    else if(!nzchar(envvar))
      stop("No julia executable in search path and JULIA_BIN environment variable not set")
  }
  if(test)
    nzchar(envvar)
  else
    envvar
}

#' Checks whether the installed version of Julia is at least a specific version.
#'
#' Compares the installed Julia version to a specific version number.
#'
#' @param version Character, on the form "X.Y.Z" where X, Y and Z are integers,
#' representing the version to compare to.
#' @return
#' \itemize{
#' \item `-1`: if couldn't check the Julia version
#' \item `0`: if the current Julia version is lower than the specified `version`
#' \item `1`: if the current Julia version is at least the specified `version`
#' }
#' @export
checkJuliaVersion <- function(version){
  ## Find Julia executable
  julia_bin = XRJulia::findJulia()

  ## Test the version of Julia
  if (.Platform$OS.type == "windows") cmd_args = paste0("-E \"VERSION < v\\\"", version,"\\\"\"") else cmd_args = paste0("-E 'VERSION < v\"", version, "\"'")
  version_check = tryCatch(system2(julia_bin, args = cmd_args, stdout = TRUE, stderr = TRUE), warning = function(e) "", error = function(e) "")

  ## If we can't check the version
  if(version_check == ""){
    return(-1)
  }

  ## If the version is too old
  else if(base::grepl("true", version_check, ignore.case = TRUE)){
    return(0)
  }

  ## If the version is ok
  return(1)
}

#' Checking whether the required Julia modules are installed
#'
#' Checking whether the required Julia modules are installed.
#' @param silence_warning Logical, whether the warning issued if required
#' Julia modules are missing must be silenced. Only for testing. Default
#' value is FALSE (i.e. by default a warning will be issued if some of the
#' required modules are missing).
#' @return Nothing
#' @export
checkJuliaModules <- function(silence_warning = FALSE){

  ## Find Julia executable
  julia_bin = XRJulia::findJulia()

  ## Locate Julia script that tests the presence of modules
  testFile = system.file("julia", "testModules.jl", package = "sismonr")

  ## Run the script
  if (.Platform$OS.type == "windows") cmd = paste0('"',julia_bin,'" "', testFile, '"') else cmd = paste(julia_bin, " ", testFile)
  check_modules_res = base::system(cmd, intern = TRUE)

  ## The script returns an integer, giving the number of required Julia modules that are missing
  ## We save a logical indicating whether any of the required module is missing
  assign("hasModules", as.numeric(check_modules_res) == 0, envir = .pkgglobalenv)

  if((as.numeric(check_modules_res) > 0) & !silence_warning) warning("Required Julia modules not installed. Please run 'installJuliaModules()' to install them.")

  return(invisible(NULL))
}

#' Install required Julia modules
#'
#' Installed required Julia modules. Requires user permission before installing the modules.
#' @return Nothing.
#' @export
installJuliaModules <- function() {

  consent <- readline("This function is about to install required Julia modules. Do you want to proceed? [y/n] ")

  ## To make sure that the user is entering 'y' or 'n' (case insensitive)
  get_out <- !stringr::str_detect(consent, "(y|Y|n|N)")
  while(get_out){
    consent <- readline("ERROR: Expecting 'y' or 'n'. This function is about to install required Julia modules. Do you want to proceed? [y/n] ")
    get_out <- !stringr::str_detect(consent, "(y|Y|n|N)")
  }

  if(stringr::str_detect(consent, "(y|Y)")){
    julia_bin = XRJulia::findJulia()
    testFile = system.file("julia", "installModules.jl", package = "sismonr")

    if (.Platform$OS.type == "windows") cmd = paste0('"',julia_bin,'" "', testFile, '"') else cmd = paste(julia_bin, " ", testFile)
    check_modules_res = base::system(cmd, intern = FALSE)

    ## So that we know the required Julia modules are available
    assign("hasModules", TRUE, envir = .pkgglobalenv)

    message("Done.")
  }

  return(invisible(NULL))
}

#' Creates a new ready-to-use Julia evaluator.
#'
#' \code{newJuliaEvaluator} opens a new Julia evaluator and loads the required functions on it.
#'
#' If no port number is specified, the \code{RJulia} function chooses a port to open the evaluator.
#'
#' @param port An integer specifying the port to be used. Default \code{NULL}.
#' @return A Julia Evaluator from the XRJulia package.
#' @examples
#' \dontrun{
#' ev = newJuliaEvaluator()
#' }
#' @export
newJuliaEvaluator <- function(port = NULL) {
  if (is.null(port)) {
    excon = showConnections()[,"description"] ## existing open connections
    excon = excon[stringr::str_detect(excon, "->localhost:")]
    excon = as.integer(stringr::str_replace(excon, "->localhost:", ""))
    newport = ifelse(length(excon) == 0, 1118L, as.integer(max(excon)+1))
    ev <- RJulia(port = newport, .makeNew = T)  ## create a new julia evaluator
  } else {
    ev <- RJulia(port = as.integer(port), .makeNew = T)  ## create a new julia evaluator with a given port
  }

  ## Global variable (logical) that tells whether the required Julia modules are installed
  hasModules <- get("hasModules", envir = .pkgglobalenv)

  ## If no check has previously been done on whether the required Julia modules are available, do it now
  if(is.null(hasModules)){
    checkJuliaModules()
    hasModules <- get("hasModules", envir = .pkgglobalenv)
  }

  ## If some required Julia modules are not available, install them (will require the used to give their consent before installing anything)
  #if(!hasModules) installJuliaModules()

  ## Sourcing sismonr Julia functions
  juliaSource(system.file("julia", "network.jl", package = "sismonr"), evaluator = ev)  ## load the julia functions of the package in the new evaluator
  juliaSource(system.file("julia", "stochastic_model.jl", package = "sismonr"), evaluator = ev)  ## load the julia functions of the package in the new evaluator
  juliaSource(system.file("julia", "stochastic_simulation.jl", package = "sismonr"), evaluator = ev)  ## load the julia functions of the package in the new evaluator

  return(ev)
}

#' Returns the current Julia evaluator.
#'
#' Returns the current Julia evaluator; if none, starts a new one.
#'
#' \code{getJuliaEvaluator} is similar to the \code{XRJulia} function \code{RJulia},
#' but if no evaluator exists, creates a new one and loads \code{sismonr} Julia functions
#' on it.
#'
#' @return A Julia evaluator from \code{XRJulia} package.
#' @examples
#' \dontrun{
#' getJuliaEvaluator()
#' }
#' @export
getJuliaEvaluator <- function(){
  ev = RJulia() ## returns the current evaluator or creates one

  ## Global variable (logical) that tells whether the required Julia modules are installed
  hasModules <- get("hasModules", envir = .pkgglobalenv)

  ## If no check has previously been done on whether the required Julia modules are available, do it now
  if(is.null(hasModules)){
    checkJuliaModules()
    hasModules <- get("hasModules", envir = .pkgglobalenv)
  }
  ## If some required Julia modules are not available, install them (will require the used to give their consent before installing anything)
  if(!hasModules) installJuliaModules()

  ## Checking whether the sismonr Julia functions are already sourced on the evaluator (faster than sourcing them everytime)
  if(!juliaEval("@isdefined juliatest")){
    juliaSource(system.file("julia", "network.jl", package = "sismonr"), evaluator = ev)  ## load the julia functions of the package in the new evaluator
    juliaSource(system.file("julia", "stochastic_model.jl", package = "sismonr"), evaluator = ev)  ## load the julia functions of the package in the new evaluator
    juliaSource(system.file("julia", "stochastic_simulation.jl", package = "sismonr"), evaluator = ev)  ## load the julia functions of the package in the new evaluator
  }
  return(ev)
}

#' Closes a Julia evaluator.
#'
#' \code{removeJuliaEvaluator} closes a Julia evaluator from XRJulia package.
#'
#' @param ev A Julia evaluator from the XRJulia package.
#' @examples
#' \dontrun{
#' removeJuliaEvaluator(getJuliaEvaluator())
#' }
#' @export
removeJuliaEvaluator <- function(ev) {
  # ev$ServerQuit()
  # isdone = XR::rmInterface(ev)
  ev$finalize()
}

#' Check whether Julia is installed with correct version
#' and all required modules.
#'
#' Use for tests only.
check_julia <- function(){
  if (!findJuliaNoError(test = T)) {
    testthat::skip("Julia not installed.")
  }else{
    ## Test that the correct version of Julia is installed
    min_version <- "1.6.5"
    version_check <- checkJuliaVersion(min_version)
    if(version_check == -1){
      testthat::skip("Error when checking the existing Julia version. Please check that Julia is correctly installed.")
    }
    else if(version_check == 0){
      testthat::skip(paste0("Julia version < v", min_version, " installed, require v", min_version, " or later"))
    }

    ## Check whether required Julia modules are installed
    ## Global variable (logical) that tells whether the required Julia modules are installed
    hasModules <- get("hasModules", envir = .pkgglobalenv)

    ## If no check has previously been done on whether the required Julia modules are available, do it now
    if(is.null(hasModules)){
      checkJuliaModules(silence_warning = TRUE)
      hasModules <- get("hasModules", envir = .pkgglobalenv)
    }

    if(!hasModules) testthat::skip("Missing required Julia modules. Run `installJuliaModules()` to install them.")
  }
}
