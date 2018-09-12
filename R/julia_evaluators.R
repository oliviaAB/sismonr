#' Creates a new ready-to-use Julia evaluator
#'
#' \code{newJuliaEvaluator} opens a new Julia evaluator and loads the required functions on it.
#'
#' If no port number is specified, the \code{RJulia} function chooses a port to open the evaluator.
#'
#' @param port An integer specifying the port to be used. Default \code{NULL}.
#' @return A Julia Evaluator from the XRJulia package.
#' @export
newJuliaEvaluator <- function(port = NULL) {
  if (is.null(port)) {
    ev <- RJulia(.makeNew = T)  ## create a new julia evaluator
  } else {
    ev <- RJulia(port = as.integer(port), .makeNew = T)  ## create a new julia evaluator with a given port
  }
  juliaUsing("ClobberingReload")
  juliaCommand(paste0("sinclude(\"", system.file("julia", "sismon.jl",
                                                          package = "sismon"), "\")"), evaluator = ev)  ## load the julia functions of the package in the new evaluator
  return(ev)
}

#' Returns the current Julia evaluator
#'
#' Returns the current Julia evaluator; if none, starts a new one
#'
#' \code{getJuliaEvaluator} is similar to the \code{XRJulia} function \code{RJulia},
#' but if no evaluator exists, creates a new one and load \code{sismon} Julia functions
#' on it.
#'
#'@return A Julia evalutator from \code{XRJulia} package.
#'@export
getJuliaEvaluator <- function(){
  ev = RJulia() ## returns the current evaluator or creates one
  if(!juliaEval("isdefined(:juliatest)")){ ## we want to know if the sismon Julia functions are already sourced on the evaluator
                                                    ## (faster than sourcing them everytime)
    juliaUsing("ClobberingReload", evaluator = ev)
    juliaCommand(paste0("sinclude(\"", system.file("julia", "sismon.jl",
                                                            package = "sismon"), "\")"), evaluator = ev)  # load the julia functions of the package in the new evaluator

  }
  return(ev)
}

#' Closes a Julia evaluator.
#'
#' \code{removeJuliaEvaluator} closes a Julia evaluator from XRJulia package.
#'
#' @param ev A Julia evaluator from the XRJulia package.
#' @export
removeJuliaEvaluator <- function(ev) {
  # ev$ServerQuit()
  # isdone = XR::rmInterface(ev)
  ev$finalize()
}
