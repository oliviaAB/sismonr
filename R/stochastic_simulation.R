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
#' @param writefile Does the julia function write the species and reactions lists in a text file?
#' @param filepath If writefile = \code{TRUE}, path to the folder in the which the files will be created.
#' @param filename If writefile = \code{TRUE}, prefix of the files created to store the lists of species and reactions.
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' @export
createStochSystem = function(insilicosystem, indargs, writefile, filepath, filename, ev = getJuliaEvaluator()){

  ## Creating the networks lists to be sent to Julia (converted to dictionaries in Julia)
  temp = names(insilicosystem$mosystem)
  for(t in temp){
    assign(t, df2list(insilicosystem[["mosystem"]][[t]]))
  }

  ## Creating the gene list to be sent to Julia (converted to dictionaries in Julia)
  assign("genes", df2list(insilicosystem[["genes"]]))

  ## If complexes and complexeskinetics are empty lists, we need to transform them into character(0)
  ## because Julia doesn't recognise empty lists
  complexes = switch((length(insilicosystem$complexes) == 0) + 1, insilicosystem$complexes, character(0))
  complexeskinetics = switch((length(insilicosystem$complexeskinetics) == 0) + 1, insilicosystem$complexeskinetics, character(0))


  message("Generating the stochastic system...")
  juliastochsystem = juliaCall("juliaCreateStochasticSystem",
                          genes, TCRN_edg, TLRN_edg, RDRN_edg, PDRN_edg, PTMRN_edg,
                          complexes, complexeskinetics,
                          as.integer(insilicosystem$sysargs$regcomplexes.size), indargs$gcnList,
                          writefile, filepath, filename, evaluator = ev)
  message("Done.")

  return(juliastochsystem)
}

#' Calls the Julia simulation
#'
#' Calls the Julia function for simulating a stochastic system
#'
#' @param stochmodel A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' @param QTLeffects The list of QTL effects coefficients of the in silico individual to be simulated (see \code{\link{createIndividual}}).
#' @param InitVar The list of initial abundance variation coefficients of the in silico individual to be simulated (see \code{\link{createIndividual}}).
#' @param genes The data-frame of genes in the system.
#' @param simtime The amount of time to simulate the model (in seconds).
#' @param modelname The name of the model.
#' @param ntrialsPerInd The number of times the simulation must be replicated.
#' @param nepochs The number of times to record the state of the system during the simulation.
#' @param simalgorithm The name of the simulation algorithm to use in the Julia function \code{simulate} from the module \code{BioSimulator}.
#' Can be one of "Direct", "FirstReaction", "NextReaction", "OptimizedDirect", "TauLeaping", "StepAnticipation".
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return The result of the simulation (a data-frame).
#' @export
callJuliaStochasticSimulation = function(stochmodel, QTLeffects, InitVar, genes, simtime, modelname, ntrialsPerInd, nepochs, simalgorithm, ev = getJuliaEvaluator()){

  genesdf = df2list(genes)
  evXR = XR::getInterface(getClass("JuliaInterface"))
  expr = gettextf("%s(%s)","juliaStochasticSimulation", evXR$ServerArglist(stochmodel, QTLeffects, InitVar, genesdf,
                                                                           simtime, modelname = modelname, ntrials = ntrialsPerInd,
                                                                           nepochs = nepochs, simalgorithm = simalgorithm))
  key = evXR$ProxyName()
  cmd = jsonlite::toJSON(c("eval", expr, key, T))
  writeLines(cmd, evXR$connection)

  while(TRUE){
    value <- readLines(evXR$connection, 1)
    #print(value)
    if(length(value) == 0)
      Sys.sleep(1)
    else
      break
  }

  res = XR::valueFromServer(value, key, T, evXR)
  mycolnames = names(sort(unlist(res@fields$colindex@fields$lookup)))
  resdf = data.frame(matrix(unlist(res@fields$columns), ncol = length(res@fields$columns), dimnames = list(c(), mycolnames)))
  return(resdf)
}

#' Simulates a in-silico system
#'
#' Simulates (stochastically) the behavious of an in silico system over time
#'
#' @param insilicosystem The in silico system to be simulated (see \code{\link{createInSilicoSystem}}).
#' @param insilicopopulation The in silico insilicopopulation to be simulated (see \code{\link{createInSilicoPopulation}}).
#' @param simtime The amount of time to simulate the model (in seconds).
#' @param nepochs The number of times to record the state of the system during the simulation.
#' @param ntrialsPerInd The number of times the simulation must be replicated.
#' @param simalgorithm The name of the simulation algorithm to use in the Julia function \code{simulate} from the module \code{BioSimulator}.
#' Can be one of "Direct", "FirstReaction", "NextReaction", "OptimizedDirect", "TauLeaping", "StepAnticipation".
#' @param writefile Does the julia function write the species and reactions lists in a text file?
#' @param filepath If writefile = \code{TRUE}, path to the folder in the which the files will be created.
#' @param filename If writefile = \code{TRUE}, prefix of the files created to store the lists of species and reactions.
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return A list composed of:
#' \itemize{
#' \item \code{resTable}: A list where each element is the data-frame of simulated expression profiles for an individual in the in silico population.
#' \item \code{runningtime}: A vector of running time of all simulations
#' \item \code{stochmodel}: A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' }
#' @export
simulateInSilicoSystem = function(insilicosystem, insilicopopulation, simtime, nepochs = -1, ntrialsPerInd = 1, simalgorithm = "Direct", writefile = T, filepath = getwd(), filename = "simulation", ev = getJuliaEvaluator()){

  stochmodel = createStochSystem(insilicosystem, insilicopopulation$indargs, writefile, filepath, filename, ev = ev)
  cat("\n")

  ## Store the running time of each simulation
  runningtime = vector("numeric", length(insilicopopulation$individualsList)*ntrialsPerInd)
  ri = 1

  ## Set a progress bar
  cat("Starting simulations at", format(Sys.time(), usetz = T), "\n")
  progress = utils::txtProgressBar(min = 0, max = length(insilicopopulation$individualsList), style = 3)

  resTable = vector("list", length(insilicopopulation$individualsList))
  names(resTable) = names(insilicopopulation$individualsList)

  for(ind in names(insilicopopulation$individualsList)){
    tic()
    simJulia = callJuliaStochasticSimulation(stochmodel, insilicopopulation$individualsList[[ind]]$QTLeffects,
                                             insilicopopulation$individualsList[[ind]]$InitVar,
                                             insilicosystem$genes, simtime, modelname = ind, ntrials = ntrialsPerInd,
                                             nepochs = nepochs, simalgorithm = simalgorithm, ev)
    temp = toc(quiet = T)
    runningtime[ri]  = temp$toc - temp$tic
    utils::setTxtProgressBar(progress, ri)
    ri = ri + 1
    resTable[[ind]] = simJulia
  }

  cat("\nMean running time per simulation: ", mean(runningtime),"seconds. \n")
  return(list("resTable" = resTable, "runningtime" = runningtime, "stochmodel" = stochmodel))
}
