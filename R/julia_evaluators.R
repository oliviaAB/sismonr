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
#' \donttest{
#' finJuliaNoError(test = T)
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

#' Creates a new ready-to-use Julia evaluator.
#'
#' \code{newJuliaEvaluator} opens a new Julia evaluator and loads the required functions on it.
#'
#' If no port number is specified, the \code{RJulia} function chooses a port to open the evaluator.
#'
#' @param port An integer specifying the port to be used. Default \code{NULL}.
#' @return A Julia Evaluator from the XRJulia package.
#' @examples
#' \donttest{
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
#  juliaUsing("ClobberingReload", evaluator = ev)
#  juliaCommand(paste0("sinclude(\"", system.file("julia", "sismonr.jl",
#                                                          package = "sismonr"), "\")"), evaluator = ev)  ## load the julia functions of the package in the new evaluator
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
#' \donttest{
#' getJuliaEvaluator()
#' }
#' @export
getJuliaEvaluator <- function(){
  ev = RJulia() ## returns the current evaluator or creates one
  if(!juliaEval("@isdefined juliatest")){ ## we want to know if the sismonr Julia functions are already sourced on the evaluator
                                                    ## (faster than sourcing them everytime)
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
#' \donttest{
#' removeJuliaEvaluator(getJuliaEvaluator())
#' }
#' @export
removeJuliaEvaluator <- function(ev) {
  # ev$ServerQuit()
  # isdone = XR::rmInterface(ev)
  ev$finalize()
}
