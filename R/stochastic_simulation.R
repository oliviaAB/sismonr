#' Tranforms a data-frame into a list.
#'
#' Transforms a data-frame into a list to be sent to Julia. The elements of the lists correspond to
#' the columns of the data-frame.
#'
#' @param mydf A data-frame.
#' @return A named list with elements corresponding to the columns of the input data-frame.
df2list = function(mydf){
  mylist = list()
  if(nrow(mydf)!=1){
    for(cols in colnames(mydf)){
      mylist[[cols]] = mydf[,cols]
    }
  }else{ ## If only 1 row in the data frame needs to convert the column element in a vector, otherwise Julia interprets each column as a single value
    for(cols in colnames(mydf)){
      mylist[[cols]] = list(mydf[,cols])
    }
  }
  return(mylist)
}

#' Creates a stochastic system from an in silico system.
#'
#' Creates a list of molecules, reactions and associated propensities to represent the in silico system.
#'
#' @param insilicosystem The in silico system (object of class \code{insilicosystem}, see \code{\link{createInSilicoSystem}}).
#' @param indargs An object of class \code{insilicoindividualargs} (i.e. a list with parameters for in silico individuals generation).
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' @export
createStochSystem = function(insilicosystem, indargs, ev = getJuliaEvaluator()){

  ## Creating the networks lists to be sent to Julia (converted to dictionaries in Julia)
  temp = names(insilicosystem$mosystem)
  for(t in temp){
    assign(t, df2list(insilicosystem[["mosystem"]][[t]]))
  }

  ## Creating the gene list to be sent to Julia (converted to dictionaries in Julia)
  assign("genes", df2list(insilicosystem[["genes"]]))

  ## If complexes and complexeskinetics are empty lists, we need to transform them into character(0)
  ## because Julia doesn't recognise empty lists
  complexes = switch((length(insilicosystem$mosystem$complexes) == 0) + 1, insilicosystem$mosystem$complexes, character(0))
  complexeskinetics = switch((length(insilicosystem$mosystem$complexeskinetics) == 0) + 1, insilicosystem$mosystem$complexeskinetics, character(0))

  message("Generating the stochastic system...")
  juliastochsystem = juliaCall("juliaCreateStochasticSystem",
                          genes, TCRN_edg, TLRN_edg, RDRN_edg, PDRN_edg, PTMRN_edg,
                          complexes, complexeskinetics,
                          as.integer(insilicosystem$sysargs$regcomplexes.size), indargs$gcnList, evaluator = ev)
  message("Done.")

  return(juliastochsystem)
}

