#' Creates a new ready-to-use Julia evaluator
#'
#' \code{new_julia_evaluator} opens a new Julia evaluator and loads the required functions on it.
#'
#' If no port number is specified, the \code{RJulia} function chooses a port to open the evaluator.
#'
#' @param port An integer specifying the port to be used. Default \code{NULL}.
#' @return A Julia Evaluator from the XRJulia package.
#' @export
new_julia_evaluator <- function(port = NULL) {
  if (is.null(port)) {
    ev <- XRJulia::RJulia(.makeNew = T)  # create a new julia evaluator
  } else {
    ev <- XRJulia::RJulia(port = as.integer(port), .makeNew = T)  # create a new julia evaluator with a given port
  }
  XRJulia::juliaUsing("ClobberingReload")
  XRJulia::juliaCommand(paste0("sinclude(\"", system.file("julia", "sismon.jl",
                                                          package = "sismon"), "\")"), evaluator = ev)  # load the julia functions of the package in the new evaluator
  return(ev)
}


#' Closes a Julia evaluator.
#'
#' \code{remove_julia_evaluator} closes a Julia evaluator from XRJulia package.
#'
#' @param ev A Julia evaluator from the XRJulia package.
#' @export
remove_julia_evaluator <- function(ev) {
  # ev$ServerQuit()
  # isdone = XR::rmInterface(ev)
  ev$finalize()
}
