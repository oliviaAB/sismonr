#' Tranforms a data-frame into a list.
#'
#' Transforms a data-frame into a list. The elements of the list correspond to
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
#' @param writefile Does the julia function write the species and reactions lists in a text file?
#' @param filepath If writefile = \code{TRUE}, path to the folder in which the files will be created (default: current working directory).
#' @param filename If writefile = \code{TRUE}, prefix of the files created to store the lists of species and reactions (default: none).
#' @param verbose If TRUE (default), print messages to signal the start and finish of the function.
#' @param ev A Julia evaluator (for the XRJulia). If none provided select the current evaluator or create one if no evaluator exists.
#' @return A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5)
#' stochsys = createStochSystem(mysystem)
#' }
#' @export
createStochSystem = function(insilicosystem, writefile = F, filepath = NULL, filename = "simulation", verbose = T, ev = getJuliaEvaluator()){

  if(is.null(filepath)){
    writefile = F
    filepath = character(0)
  }

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


  if(verbose) message("Generating the stochastic system...")
  juliastochsystem = juliaCall("juliaCreateStochasticSystem",
                          genes, get("TCRN_edg"), get("TLRN_edg"), get("RDRN_edg"), get("PDRN_edg"), get("PTMRN_edg"),
                          complexes, complexeskinetics, insilicosystem$sysargs$gcnList,
                          writefile, filepath, filename, evaluator = ev)
  if(verbose) message("Done.")

  return(juliastochsystem)
}

#' Calls the Julia simulation function.
#'
#' Calls the Julia function for simulating a stochastic system. Should not be used by itself (this function is called by the wrapper functions \link{simulateInSilicoSystem} and \link{simulateParallelInSilicoSystem}).
#'
#' @param stochmodel A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' @param QTLeffects The list of QTL effects coefficients of the in silico individual to be simulated (see \code{\link{createIndividual}}).
#' @param InitAbundance The list of initial abundances of the molecules for the in silico individual to be simulated (see \code{\link{createIndividual}}).
#' @param genes The data-frame of genes in the system.
#' @param simtime Numeric. The amount of time to simulate the model (in seconds).
#' @param modelname String. The name of the model. Default value "MySimulation".
#' @param ntrials Integer. The number of times the simulation must be replicated.
#' @param nepochs Integer. The number of times to record the state of the system during the simulation.
#' @param simalgorithm String. The name of the simulation algorithm to use in the Julia function \code{simulate} from the module \code{BioSimulator}.
#' Possible values are: "Direct", "EnhancedDirect", "SortingDirect", "FirstReaction", "NextReaction", "TauLeapingDG2001", "TauLeapingDGLP2003", "StepAnticipation", "HybridSAL".
#' See \url{https://alanderos91.github.io/BioSimulator.jl/dev/man/algorithms/} for details about the algorithms.
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return The result of the simulation (a data-frame).
#' @export
callJuliaStochasticSimulation = function(stochmodel, QTLeffects, InitAbundance, genes, simtime, modelname = "MySimulation", ntrials, nepochs, simalgorithm, ev = getJuliaEvaluator()){

  genesdf = df2list(genes)
  evXR = XR::getInterface(getClass("JuliaInterface"))
  expr = gettextf("%s(%s)","juliaStochasticSimulation", evXR$ServerArglist(stochmodel, QTLeffects, InitAbundance, genesdf,
                                                                           simtime, modelname = modelname, ntrials = ntrials,
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

#' Simulates an in silico system.
#'
#' Simulates (stochastically) the behaviour of an in silico system over time,
#' i.e. the expression of the different genes.
#'
#' @param insilicosystem The in silico system to be simulated (see
#'   \code{\link{createInSilicoSystem}}).
#' @param insilicopopulation The in silico population to be simulated (see
#'   \code{\link{createInSilicoPopulation}}).
#' @param simtime The final time of the simulation (in seconds).
#' @param nepochs The number of times to record the state of the system during
#'   the simulation.
#' @param ntrials The number of times the simulation must be replicated (for
#'   each individual).
#' @param simalgorithm The name of the simulation algorithm to use in the Julia
#'   function \code{simulate} from the module \code{BioSimulator}. Possible
#'   values are: "Direct", "EnhancedDirect", "SortingDirect", "FirstReaction",
#'   "NextReaction", "TauLeapingDG2001", "TauLeapingDGLP2003",
#'   "StepAnticipation", "HybridSAL". See
#'   \url{https://alanderos91.github.io/BioSimulator.jl/dev/man/algorithms/} for
#'   details about the algorithms.
#' @param writefile Does the julia function write the species and reactions
#'   lists in a text file?
#' @param filepath If writefile = \code{TRUE}, path to the folder in which the
#'   files will be created (default: current working directory).
#' @param filename If writefile = \code{TRUE}, prefix of the files created to
#'   store the lists of species and reactions.
#' @param ev A Julia evaluator. If none provided select the current evaluator or
#'   create one if no evaluator exists.
#' @return A list composed of: \itemize{ \item \code{Simulation}: A data-frame
#'   with the simulated expression profiles of the genes for the different
#'   individuals in the in silico population. For gene i, "Ri" corresponds to
#'   the RNA form of the gene, "Pi" to the protein form of the gene. The suffix
#'   "GCNj" indicates that the molecule comes from the j-th allele of the gene.
#'   \item \code{runningtime}: A vector of running time of all runs of the
#'   simulation for each in silico individuals. \item \code{stochmodel}: A Julia
#'   proxy object to retrieve the stochastic system in the Julia evaluator. }
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none", ploidy = 2)
#' mypop = createInSilicoPopulation(1, mysystem)
#' sim = simulateInSilicoSystem(mysystem, mypop, simtime = 1000, ntrials = 10)
#' head(sim$Simulation)
#' ## Visualising the result
#' plotSimulation(sim$Simulation)
#' }
#' @export
simulateInSilicoSystem = function(insilicosystem, insilicopopulation, simtime, nepochs = -1, ntrials = 1, simalgorithm = "Direct", writefile = F, filepath = NULL, filename = "simulation", ev = getJuliaEvaluator()){

  if(is.null(filepath)) writefile = F
  stochmodel = createStochSystem(insilicosystem, writefile, filepath, filename, verbose = T, ev = ev)
  message("\n")

  ## Store the running time of each simulation
  runningtime = vector("numeric", length(insilicopopulation$individualsList))
  ri = 1

  ## Set a progress bar
  message("Starting simulations at ", format(Sys.time(), usetz = T), "\n")
  progress = utils::txtProgressBar(min = 0, max = length(insilicopopulation$individualsList), style = 3)

  resTable = vector("list", length(insilicopopulation$individualsList))
  names(resTable) = names(insilicopopulation$individualsList)

  for(ind in names(insilicopopulation$individualsList)){
    tic()
    simJulia = callJuliaStochasticSimulation(stochmodel, insilicopopulation$individualsList[[ind]]$QTLeffects,
                                             insilicopopulation$individualsList[[ind]]$InitAbundance,
                                             insilicosystem$genes, simtime, modelname = ind, ntrials = ntrials,
                                             nepochs = nepochs, simalgorithm = simalgorithm, ev)
    temp = toc(quiet = T)
    runningtime[ri]  = temp$toc - temp$tic
    utils::setTxtProgressBar(progress, ri)
    ri = ri + 1
    resTable[[ind]] = simJulia %>% mutate("Ind" = ind)
  }
  res = bind_rows(resTable)
  message("\nSimulations finished at ", format(Sys.time(), usetz = T), "\n")
  message("Mean running time per simulation: ", round(mean(runningtime), 3)," seconds. \n")
  return(list("Simulation" = res, "runningtime" = runningtime, "stochmodel" = stochmodel))
}


## function to start a Julia evaluator on a node of the cluster, given the port id (for parallel simulation)
startJuliaEvCluster = function(portid, insilicosystem){
  myev = newJuliaEvaluator(port = portid) ## start on the node a Julia evaluator with specified port number
  mystochmodel = createStochSystem(insilicosystem, writefile = F, verbose = F, ev = myev)
  return(list("myev" = myev, "mystochmodelvar" = mystochmodel@.Data)) ## return the Julia evaluator ID and the name on the Julia process of the stochmodel object
}

## function to run a simulation on a cluster (for parallel simulation)
simulateInCluster = function(i, indtosimulate, ntrialstosimulate, increment, individualsList, genes, simtime, nepochs, simalgorithm, infocores, no_cores, progress, maxprogress){
  myinfocore = infocores[[i - no_cores*(i-1)%/%no_cores]] ## get the infos of the corresponding cluster node
  myev = myinfocore$myev ## get the Julia evaluator corresponding to the current cluster node
  mystochmodel = juliaEval(paste0(myinfocore$mystochmodel), .get = F) ## get a proxy object corresponding to the stochmodel object on the Julia process using the variable name
  ind = indtosimulate[i]
  ntrialsclus = ntrialstosimulate[i]
  simJulia = callJuliaStochasticSimulation(mystochmodel, individualsList[[ind]]$QTLeffects,
                                           individualsList[[ind]]$InitAbundance,
                                           genes, simtime, modelname = ind, ntrials = ntrialsclus,
                                           nepochs = nepochs, simalgorithm = simalgorithm, myev)


  if(i %%no_cores == 1){
    utils::setTxtProgressBar(progress, min(i + (no_cores-1), maxprogress))
  }
  return(simJulia %>% mutate("trial" = !!sym("trial") + increment[i], "Ind" = ind))
}


#' Simulates an in silico system in parallel.
#'
#' Simulates (stochastically) the behaviour of an in silico system over time
#' using parallelisation, i.e. the expression of the different genes.
#'
#' @param insilicosystem The in silico system to be simulated (see
#'   \code{\link{createInSilicoSystem}}).
#' @param insilicopopulation The in silico population to be simulated (see
#'   \code{\link{createInSilicoPopulation}}).
#' @param simtime The final time of the simulation (in seconds).
#' @param nepochs The number of times to record the state of the system during
#'   the simulation.
#' @param ntrials The number of times the simulation must be replicated (for
#'   each individual).
#' @param simalgorithm The name of the simulation algorithm to use in the Julia
#'   function \code{simulate} from the module \code{BioSimulator}. Possible
#'   values are: "Direct", "EnhancedDirect", "SortingDirect", "FirstReaction",
#'   "NextReaction", "TauLeapingDG2001", "TauLeapingDGLP2003",
#'   "StepAnticipation", "HybridSAL". See
#'   \url{https://alanderos91.github.io/BioSimulator.jl/dev/man/algorithms/} for
#'   details about the algorithms.
#' @param writefile Does the julia function write the species and reactions
#'   lists in a text file?
#' @param filepath If writefile = \code{TRUE}, path to the folder in which the
#'   files will be created (default: current working directory).
#' @param filename If writefile = \code{TRUE}, prefix of the files created to
#'   store the lists of species and reactions.
#' @param no_cores The number of cores to use for the simulation. By default use
#'   the function \code{detectCores} from the \code{parallel} package to detect
#'   the number of available cores, and use this number - 1 for the simulation.
#' @param ev A Julia evaluator. If none provided select the current evaluator or
#'   create one if no evaluator exists.
#' @return A list composed of: \itemize{ \item \code{Simulation}: A data-frame
#'   with the simulated expression profiles of the genes for the different
#'   individuals in the in silico population. For gene i, "Ri" corresponds to
#'   the RNA form of the gene, "Pi" to the protein form of the gene. The suffix
#'   "GCNj" indicates that the molecule comes from the j-th allele of the gene.
#'   \item \code{runningtime}: The running time (elapsed seconds) of the
#'   parallel simulation (only 1 value). \item \code{stochmodel}: A Julia proxy
#'   object to retrieve the stochastic system in the Julia evaluator. }
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none", ploidy = 2)
#' mypop = createInSilicoPopulation(15, mysystem)
#' sim = simulateParallelInSilicoSystem(mysystem, mypop, 1000)
#' head(sim$Simulation)
#' ## Visualising the result
#' plotSimulation(sim$Simulation)
#' }
#' @export
simulateParallelInSilicoSystem = function(insilicosystem, insilicopopulation, simtime, nepochs = -1, ntrials = 1, simalgorithm = "Direct", writefile = F, filepath = NULL, filename = "simulation", no_cores = parallel::detectCores()-1, ev = getJuliaEvaluator()){

  if(is.null(filepath)) writefile = F

  stochmodel = createStochSystem(insilicosystem, writefile, filepath, filename, verbose = T, ev = ev)
  message("\n")


  mybaseport = ev$port ## Get the port of the current Julia evaluator
  portList = sapply(1:no_cores, sum, mybaseport) ## Assign to each core a port number, starting from 1+port number of the current evaluator

  mycluster = parallel::makeCluster(no_cores, outfile = "")

  parallel::clusterEvalQ(mycluster, library(XRJulia))
  parallel::clusterEvalQ(mycluster, library(utils))
  parallel::clusterExport(mycluster, "newJuliaEvaluator")
  parallel::clusterExport(mycluster, "createStochSystem")
  parallel::clusterExport(mycluster, "callJuliaStochasticSimulation")
  # parallel::clusterExport(mycluster, "stochmodel_string", envir = environment())
  # parallel::clusterExport(mycluster, "no_cores", envir = environment())

  ## Start a Julia evaluator on each cluster
  message("Starting Julia evaluators on each cluster node ... \n")
  infocores = parallel::clusterApply(mycluster, portList, startJuliaEvCluster, insilicosystem)
  message("Done.\n")
  # parallel::clusterExport(mycluster, "infocores", envir = environment())

  ## Split the workload for the different cluster nodes
  ## Case 1: there are more individuals than cluster nodes: in this case each node will be used to simulate a different individual
  ## Case 2: there are less individuals than cluster nodes: in this case each node will simulate approx ntrials/nnodes trials for each individual
  if(length(insilicopopulation$individualsList) >= no_cores){
    indtosimulate = names(insilicopopulation$individualsList)
    ntrialstosimulate = rep(ntrials, length(insilicopopulation$individualsList))
    increment = rep(0, length(insilicopopulation$individualsList)) ## no use here
  }else{
    split = rep(ntrials %/% no_cores, no_cores) ## vector with no_cores elements, giving the # of trials that each node has to run for an individual
    remain = ntrials %% no_cores ## we start by splitting equally the number of trials over the different cores
    split[min(1, remain):remain] = split[min(1, remain):remain] + 1 ## then if there is some remaining trial, we split them over the first cores
    split = split[split!=0] ## if ntrials<no_cores, we will have some cores that do nothing for an individual
    indtosimulate = rep(names(insilicopopulation$individualsList), each = length(split))
    ntrialstosimulate = rep(split, times = length(insilicopopulation$individualsList))
    increment = rep(c(0, cumsum(split)[-length(split)]), times = length(insilicopopulation$individualsList)) ## each node will increment the trials ID by the number provided in increment[i]
                                                                                                             ## so that if we have to simulate Ind1 10 times in total, and we split as following
                                                                                                             ## node 1 = 4 trials, node2 = 3 trials, node 3 = 3 trials
                                                                                                             ## node 1 will return trials 1 to 4, node 2: 5 to 7, node 3: 8 to 10
    }

  message("Starting simulations at ", format(Sys.time(), usetz = T), "\n")

  ## Create the progress bar
  maxprogress = length(indtosimulate)
  progress = utils::txtProgressBar(min = 0, max = maxprogress, style = 3)
  # parallel::clusterExport(mycluster, "progress", envir = environment())
  # parallel::clusterExport(mycluster, "maxprogress", envir = environment())

  ## Run the simulation on the cluster
  startsim = tic()
  resTable = parallel::clusterApply(mycluster, 1:length(indtosimulate), simulateInCluster, indtosimulate = indtosimulate, ntrialstosimulate = ntrialstosimulate, increment = increment, individualsList = insilicopopulation$individualsList, genes = insilicosystem$genes, simtime = simtime, nepochs = nepochs, simalgorithm = simalgorithm, infocores = infocores, no_cores = no_cores, progress = progress, maxprogress = maxprogress)
  stopsim = toc(quiet = T)$toc
  res = bind_rows(resTable)
  message("\nSimulations finished at ", format(Sys.time(), usetz = T), "\n")
  message("Running time of parallel simulations: ", round(stopsim - startsim, 3), " seconds\n")

  parallel::stopCluster(mycluster)

  return(list("Simulation" = res, "runningtime" = stopsim - startsim, "stochmodel" = stochmodel))
}


sumColAbundance = function(df, colsid){
  if(length(colsid) == 1){
    return(df[[colsid]])
  }else{
    return(rowSums(df[,colsid]))
  }
}

#' Merge the different allelic versions of the molecules.
#'
#' Merge (i.e. sum) the abundance of the different allelic versions of each molecule in the results of a simulation.
#'
#' @param df A dataframe with the abundance of the different molecules over time (from \code{\link{simulateInSilicoSystem}}
#' or \code{\link{simulateParallelInSilicoSystem}}).
#' @return A dataframe in which the abundance of the different allelic versions of the same molecule have been merged to give the abundance of the molecule (without distinction of the allele of origin).
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, empty = TRUE, ploidy = 2)
#' mypop = createInSilicoPopulation(1, mysystem)
#' sim = simulateInSilicoSystem(mysystem, mypop, 100)
#' head(sim$Simulation)
#' mergedAllelic = mergeAlleleAbundance(sim$Simulation)
#' head(mergedAllelic)
#' }
#' @export
mergeAlleleAbundance = function(df){
  mergeddf = df %>% select(!!sym("time"), !!sym("trial"), !!sym("Ind"))
  molsGCN = colnames(df)
  mols = stringr::str_replace_all(molsGCN, "GCN[[:digit:]]+", "")

  for(m in setdiff(unique(mols), c("time", "trial", "Ind"))){
    mergeddf[[m]] = sumColAbundance(df,which(mols == m))
  }

  return(mergeddf)
}

#' Merge the original and PTM versions of the proteins.
#'
#' Merge (i.e. sum) the abundance of the original and modified (PTM) versions of each protein in the results of a simulation.
#'
#' @param df A dataframe with the abundance of the different molecules over time (from \code{\link{simulateInSilicoSystem}}
#' or \code{\link{simulateParallelInSilicoSystem}}).
#' @return A dataframe in which the abundance of original and modified versions of a protein have been merged to give the abundance of the protein (without distinction of its post-translational modification state).
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, PC.p = 1, PC.PTM.p = 0.9, regcomplexes = "none", ploidy = 1)
#' mypop = createInSilicoPopulation(1, mysystem)
#' sim = simulateInSilicoSystem(mysystem, mypop, 100)
#' head(sim$Simulation)
#' mergedPTM = mergePTMAbundance(sim$Simulation)
#' head(mergedPTM)
#' }
#' @export
mergePTMAbundance = function(df){
  mergeddf = df %>% select(!!sym("time"), !!sym("trial"), !!sym("Ind"))
  molsPTM = colnames(df)
  mols = stringr::str_replace(molsPTM, "Pm", "P")

  for(m in setdiff(unique(mols), c("time", "trial", "Ind"))){
    mergeddf[[m]] = sumColAbundance(df,which(mols == m))
  }

  return(mergeddf)
}

#' Merge the free and in-complex versions of molecules.
#'
#' Merge (i.e. sum) the abundance of the free and in-complex versions of each molecule in the results of a simulation.
#'
#' @param df A dataframe with the abundance of the different molecules over time (from \code{\link{simulateInSilicoSystem}}
#' or \code{\link{simulateParallelInSilicoSystem}}).
#' @return A dataframe in which the abundance of free and in complex versions of a molecule have been merged to give the abundance of the molecule (without distinction of whether or not it is bound in a molecular complex).
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, PC.p = 1, PC.TC.p = 1, ploidy = 1)
#' mysystem = addComplex(mysystem, c(1, 2))
#' mypop = createInSilicoPopulation(1, mysystem)
#' sim = simulateInSilicoSystem(mysystem, mypop, 100)
#' head(sim$Simulation)
#' mergedComplex = mergeComplexesAbundance(sim$Simulation)
#' head(mergedComplex)
#' }
#' @export
mergeComplexesAbundance = function(df){
  molsComp = colnames(df)[stringr::str_detect(colnames(df), "^C")]
  mergeddf = df %>% select(-molsComp)

  for(i in molsComp){
    splt = stringr::str_split(i, "_")[[1]]
    splt = splt[!stringr::str_detect(splt, "^C")]
    for(j in unique(splt)){
      mergeddf[, j] = mergeddf[, j] + df[, i]*sum(splt == j)
    }
  }

  return(mergeddf)
}

chooseColourPalette = function(toplot, mergeAllele){

    if(!mergeAllele){
    ## we need to see how many different components exist in the system and how many allelic version of each exist
    compocols = data.frame(Components = unique(toplot$Components)) %>%
      mutate("isComplex" = stringr::str_detect(.data$Components, "^C")) %>%
      mutate("compoID" = stringr::str_replace_all(.data$Components, "GCN.+|_.+", ""),
             "Allele" = stringr::str_replace(.data$Components, "^[[:digit:]]+(?=GCN)|^[^_]+_{1}", "")) %>%
      dplyr::arrange(.data$compoID, .data$Allele)

    ncols = length(unique(compocols$compoID))
    #paletteCompo = randomcoloR::distinctColorPalette(ncols) ## create a colour for each component
    paletteCompo = grDevices::rainbow(ncols)
    names(paletteCompo) = sort(unique(compocols$compoID))
    palette = unlist(lapply(names(paletteCompo), function(x){ ## create a gradient (one colour for each allelic version of the component)
      nall = sum(compocols$compoID == x)
      return(rev(grDevices::colorRampPalette(c("#FFFFFF", paletteCompo[x]))(nall+1)[-1]))
    }))
    names(palette) = sort(compocols$Components)
  }else{
    #palette = randomcoloR::distinctColorPalette(length(unique(toplot$Components)))
    palette = grDevices::rainbow(length(unique(toplot$Components)))
    names(palette) = sort(unique(toplot$Components))
  }
  return(palette)
}

plotBase = function(toplot, palette, multitrials, yLogScale, ...){
  if(multitrials){
    simuplot = ggplot2::ggplot(toplot, aes_string(x = "time")) +
      ggplot2::geom_ribbon(aes_string(ymin = "LB", ymax = "UB", fill = "Components"), alpha = 0.5) +
      ggplot2::geom_line(aes_string(y = "mean", colour = "Components")) +
      ggplot2::scale_colour_manual(values = palette, breaks = names(palette), guide = guide_legend(ncol=10, byrow=T)) +
      ggplot2::scale_fill_manual(values = palette, breaks = names(palette), guide = guide_legend(ncol=10, byrow=T))
  }else{
    simuplot = ggplot2::ggplot(toplot, aes_string(x = "time", y = "Abundance", colour = "Components")) + ggplot2::geom_line() +
      ggplot2::scale_colour_manual(values = palette, breaks = names(palette), guide = guide_legend(ncol=10, byrow=T))
  }

  plotTitle = "Components absolute abundance"
  if(yLogScale){
    simuplot = simuplot + ggplot2::scale_y_log10()
    plotTitle = "log10(Components absolute abundance+0.5)"
  }

  simuplot = simuplot +
    ggplot2::facet_grid(Type~Ind, scales = "free_y") +
    ggplot2::xlab("Time (s)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.text = element_text(size = 7), strip.text = element_text(size = 10),
                   legend.spacing = grid::unit(5, "mm"), legend.key.size = grid::unit(10, "mm"),
                   legend.direction = "horizontal", ...) +
    ggplot2::ylab(plotTitle)
  return(simuplot)
}

#' Sort component names.
#'
#' Sorts the names of the components of an in silico system (for plotting or summary).
#'
#' Sort components into:
#' \itemize{
#' \item Genes (first) vs regulatory complexes
#' }
#' Then genes are sorted according to:
#' \itemize{
#' \item Gene ID (numeric)
#' \item Gene type (RNAs, then proteins, then modified proteins)
#' \item Allele
#' \item If the components don't have a type (e.g. see the legend resulting from \code{\link{plotSimulation}}),
#' non-modified vs modified
#' }
#' Then sort complexes according to:
#' \itemize{
#' \item Target reaction in the following order: regulators of transcription, translation, RNA decay,
#' protein decay, post-translational modification
#' \item Allele of their components
#' }
#'
#' @param componames A character vector, giving the names of the components.
#' @return A dataframe: first column: sorted names of the components, second column: is the component a
#' regulatory complex?, third column: is the component a modified protein?
#' @export
sortComponents = function(componames){

  componames = as.character(componames)
  ## To sort the genes
  genesDF = data.frame("Components" = componames, stringsAsFactors = F) %>%
    filter(!(stringr::str_detect(!!sym("Components"), "^C"))) %>%
    mutate("GeneID" = as.numeric(stringr::str_extract(!!sym("Components"),"(?<!GCN)\\d+")),
           "Allele" = as.numeric(stringr::str_extract(!!sym("Components"),"(?<=GCN)\\d+")),
           "Allele" = case_when(is.na(.data$Allele)~0,
                                !is.na(.data$Allele)~.data$Allele),
           "PTM" = case_when(stringr::str_detect(.data$Components, "^PTM")~1,
                             !stringr::str_detect(.data$Components, "^PTM")~0),
           "Type" = case_when(stringr::str_detect(.data$Components, "^\\d|^PTM")~0,
                              stringr::str_detect(.data$Components, "^R")~1,
                              stringr::str_detect(.data$Components, "^P\\d")~2,
                              stringr::str_detect(.data$Components, "^Pm")~3)) %>%
    dplyr::arrange(!!sym("GeneID"), !!sym("Type"), !!sym("Allele"), !!sym("PTM"))

  ## To sort the complexes when plotting the legend
  complexesDF = data.frame("Components" = componames, stringsAsFactors = F) %>%
    filter((stringr::str_detect(!!sym("Components"), "^C"))) %>%
    mutate("TargetReaction" = stringr::str_extract(!!sym("Components"),"(?<=^C)[[:alpha:]]{2,3}"),
           "TargetReaction" = case_when(.data$TargetReaction == "TC" ~ 1,
                                        .data$TargetReaction == "TL" ~ 2,
                                        .data$TargetReaction == "RD" ~ 3,
                                        .data$TargetReaction == "PD" ~ 4,
                                        .data$TargetReaction == "PTM" ~ 5),
           "ComplID" = as.numeric(stringr::str_extract(!!sym("Components"),"(?<=^C[[:alpha:]]{2,3})\\d+"))) %>%
    dplyr::arrange(!!sym("TargetReaction"), !!sym("ComplID"), !!sym("Components"))

  sorted = data.frame("Components" = c(genesDF$Components, complexesDF$Components),
                      "isComplex" = rep(c(F, T), c(nrow(genesDF), nrow(complexesDF))),
                      "isPTM" = c((genesDF$PTM == 1 | genesDF$Type == 3), rep(F, nrow(complexesDF))),stringsAsFactors = F)
  return(sorted)
}

plotLegendComponents = function(palette, nCompPerRow = 10, components){

  sortedComp = sortComponents(names(palette))
  sortedComp = cbind(sortedComp, data.frame("isPC" = sapply(sortedComp$Components, function(x){
    any(stringr::str_detect(components, paste0("^P", x, "(?=_|$)")))
  })))
  palette = palette[sortedComp$Components]

  nrows = length(palette) %/% nCompPerRow + (length(palette) %% nCompPerRow != 0)
  plots = vector("list", length = nrows)

  for(i in 1:nrows){
    rowpalette = palette[((i-1)*nCompPerRow+1):min(i*nCompPerRow, length(palette))]
    rowisComplex = sortedComp$isComplex[((i-1)*nCompPerRow+1):min(i*nCompPerRow, length(palette))]
    rowisPTM = sortedComp$isPTM[((i-1)*nCompPerRow+1):min(i*nCompPerRow, length(palette))]
    rowisPC = sortedComp$isPC[((i-1)*nCompPerRow+1):min(i*nCompPerRow, length(palette))]

    # namesComp = names(rowpalette)
    # namesComp[is.na(namesComp)] = ""
    # cols = unname(rowpalette)
    # cols[is.na(cols)] = "white"
    nbc = sum(rowisComplex)
    nbg = length(rowpalette) - nbc

    rowisPTM = rowisPTM[1:nbg]

    compoType = data.frame("isPC" = rowisPC[1:nbg], "isPTM" = rowisPTM)

    if(nbg > 0){
      legend_points = data.frame("Components" = rep(names(rowpalette)[1:nbg], sapply(1:nbg, function(x){
                                          if(compoType[x, "isPTM"] | !compoType[x, "isPC"]){
                                            1
                                          }else{2}
                                        })),
                                     "x" = rep(1:nbg, sapply(1:nbg, function(x){
                                       if(compoType[x, "isPTM"] | !compoType[x, "isPC"]){
                                         1
                                       }else{2}
                                     })),
                                     "y" = unlist(lapply(1:nbg, function(x){
                                       if(compoType[x, "isPTM"]){
                                         2
                                       }else if(compoType[x, "isPC"]){
                                           1:2
                                         }else{1}
                                       })))
    }else{
      legend_points = data.frame("Components" = character(),
                                 "x" = numeric(),
                                 "y" = numeric())
    }

    if(nbc > 0) legend_points = rbind(legend_points, data.frame("Components" = names(rowpalette)[(nbg+1):(nbg + nbc)],
                                     "x" = (1:nbc) + nbg,
                                     "y" = rep(3, nbc)))

    legend_text =  rbind(data.frame("Names" = names(rowpalette),
                                    "x" = 1:(nbg+nbc),
                                    "y" = 4.1,
                                    "angle" = 0,
                                    "size" = 3,
                                    "hjust" = 0.5, stringsAsFactors = F),
                         data.frame("Names" = c("RNAs", "Proteins", "Complexes"),
                                    "x" = 0.2,
                                    "y" = 1:3,
                                    "angle" = 0,
                                    "size" = 4,
                                    "hjust" = 1), stringsAsFactors = F)
    legend_lines = rbind(data.frame("x" = seq(0.5, length(rowpalette)+0.5, by = 1), "xend" = seq(0.5, length(rowpalette)+0.5, by = 1), "y" = 0.5, "yend" = 3.9),
                         data.frame("x" = 0.4, "xend" = length(rowpalette)+0.6, "y" = seq(0.5, 3.5, by = 1), "yend" = seq(0.5, 3.5, by = 1)))

    rowplot = ggplot2::ggplot() +
      ggplot2::geom_point(data = legend_points, size = 5, shape = 15, aes_string(x = "x", y = "y", group = "Components", colour = "Components"), show.legend = F) +
      ggplot2::geom_text(data = legend_text, aes_string(label = "Names", x = "x", y = "y", angle = "angle", size = "size", hjust = "hjust"), show.legend = F) +
      ggplot2::geom_segment(data = legend_lines, colour = "gray90", aes_string(x = "x", xend = "xend", y = "y", yend = "yend")) +
      ggplot2::scale_colour_manual(values = rowpalette, breaks = names(rowpalette)) +
      ggplot2::scale_size(range = c(2.5, 3.5), guide = F) +
      ggplot2::xlim(-1.6, nCompPerRow+1) + ggplot2::ylim(0.45, 4.9) +
      ggplot2::theme_void()

    plots[[i]] = rowplot

  }
  return(plots)
}

#' Plots the result of a simulation.
#'
#' Automatically plots the result of a simulation (i.e. the abundance of RNAs, proteins and complexes over time)
#' for the selected in silico individuals.
#'
#' If more than one trial is to be plotted, the mean abundance of each molecule over the different trials is plotted with a solid line,
#' and the min and max abundances represented as coloured areas around the mean.
#'
#' @param simdf The dataframe with the result of the simulation (see \code{\link{simulateInSilicoSystem}}).
#' @param molecules A vector of gene IDs (numeric or character) and/or complex IDs (e.g. CTC1) to be plotted.
#' @param inds A vector of in silico individual names for which to plot the expression profiles.
#' @param trials A vector of trials ID (= number) to use for the plot (see details).
#' @param timeMin Numeric. The minimum simulation time to plot. Default value set to the minimum time in the simulation.
#' @param timeMax Numeric. The maximum simulation time to plot. Default value set to the maximum time in the simulation.
#' @param mergeAllele Are the gene products originating from different alleles merged? Default TRUE. Also see \code{\link{mergeAlleleAbundance}}
#' @param mergePTM Are the modified and non-modified versions of the proteins merged? Default TRUE. Also see \code{\link{mergePTMAbundance}}
#' @param mergeComplexes Are the free and in complex gene products merged? Default FALSE. Also see \code{\link{mergeComplexesAbundance}}
#' @param yLogScale Plot the y-axis in log10-scale? If so, the abundance of each species at each time-point is increased by 1 to avoid zero values. Default TRUE.
#' @param nIndPerRow Positive integer, the number of individuals to plot per row. Default 3.
#' @param nCompPerRow Positive integer, the number of components to plot per row in the legend. Default 10.
#' @param ... Any additional parameter to be passed to \code{\link[ggplot2]{theme}} for the plot of each individual.
#' @return A plot from \code{\link[ggpubr]{ggarrange}}.
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none", ploidy = 2)
#' mypop = createInSilicoPopulation(15, mysystem)
#' sim = simulateInSilicoSystem(mysystem, mypop, 100, ntrials = 5)
#' plotSimulation(sim$Simulation,
#'  c(1, 2, 3),
#'  c("Ind1", "Ind2", "Ind3", "Ind4"),
#'  axis.title = element_text(color = "red"))
#' }
#' @export
plotSimulation = function(simdf, molecules = NULL, inds = unique(simdf$Ind), trials = unique(simdf$trial), timeMin = min(simdf$time), timeMax = max(simdf$time), mergeAllele = T, mergePTM = T, mergeComplexes = F, yLogScale  = T, nIndPerRow = 3, nCompPerRow = 10, ...){

  ## Select the requested inviduals (in principle we could only do this later but this avoids transforming the whole dataset
  ## if we want to plot only a couple of individuals)
  simind = simdf %>% dplyr::filter(!!sym("Ind") %in% inds) %>%
                     dplyr::filter(!!sym("trial") %in% trials) %>%
                     dplyr::filter(!!sym("time") >= timeMin & !!sym("time") <= timeMax)

  ## Are there mutliple trials?
  multitrials = length(unique(simind$trial)) > 1

  if(mergePTM)  simind = mergePTMAbundance(simind)
  if(mergeComplexes)  simind = mergeComplexesAbundance(simind)
  if(mergeAllele) simind = mergeAlleleAbundance(simind)

  ## First transformation into a long tibble
  toplot = simind %>%
    tidyr::gather(key = "ID", value = "Abundance", setdiff(names(simind), c("time", "trial", "Ind"))) %>%
    mutate("Type" = case_when(stringr::str_detect(.data$ID, "^R") ~ "RNAs",
                            stringr::str_detect(.data$ID, "^P") ~ "Proteins",
                            stringr::str_detect(.data$ID, "^C") ~ "Complexes"),
           "Components" = stringr::str_replace(.data$ID, "^R|^P", ""),
           "Components" = stringr::str_replace(.data$Components, "^m", "PTM"))

  ## Select the molecules to plot
  if(!is.null(molecules)){
    toplot = toplot %>%
      mutate("Molecule" = case_when(!!sym("Type") %in% c("RNAs", "Proteins") ~ stringr::str_extract(.data$Components, "(?<=^|(PTM))\\d+"),
                                    !!sym("Type") == "Complexes" ~ stringr::str_extract(.data$Components, "^C[[:alpha:]]{2,3}\\d+"))) %>%
      dplyr::filter(!!sym("Molecule") %in% molecules) %>%
      select(-!!sym("Molecule"))
  }

  ## If plot in log-scale, need to have non-zero abundances
  if(yLogScale) toplot = toplot %>% mutate("Abundance" = !!sym("Abundance") + 0.5)

  ## If multiple trials, summarise them with mean, min and max
  if(multitrials){
    toplot = toplot %>%
      dplyr::group_by(!!sym("Ind"), !!sym("time"), !!sym("Components"), !!sym("Type"), !!sym("ID")) %>%
      dplyr::summarise("mean" = mean(!!sym("Abundance")), "LB" = min(!!sym("Abundance")), "UB" = max(!!sym("Abundance")))
  }

  palette = chooseColourPalette(toplot, mergeAllele)

  # ## How many rows of plot given the nb of inds to plot and nb of inds per row
  # nbrows = c(rep(nIndsPerRow, length(inds) %/% nIndsPerRow), length(inds) %% nIndsPerRow)
  #
  # ## Create each row of plot
  # firstindex = 1
  # plots = vector("list", length = length(nbrows))
  # for(i in 1:length(nbrows)){
  #   indsrow = inds[firstindex:(firstindex + nbrows[i] - 1)]
  #   plots[[i]] = plotBase(filter(toplot, !!as.name("Ind") %in% indsrow), palette, multitrials, yLogScale)
  #   firstindex = firstindex + nbrows[i]
  # }

  ## Create a plot for each individual
  plots = vector("list", length = length(inds))
  for(i in 1:length(inds)){
    plots[[i]] = plotBase(filter(toplot, !!sym("Ind") == inds[i]), palette, multitrials, yLogScale, ...)
  }

  ## How many rows of plot given the nb of inds to plot and nb of inds per row
  nbrows = length(inds) %/% nIndPerRow + (length(inds) %% nIndPerRow != 0)

  simuplot = ggpubr::ggarrange(plotlist = plots, nrow = nbrows, ncol = min(nIndPerRow, length(inds)), common.legend = T, heights = rep(1, length(plots)), legend = "none")

  ## add the legend
  legendplots = plotLegendComponents(palette, nCompPerRow, unique(toplot$ID))
  finalplot = ggpubr::ggarrange(simuplot, plotlist = legendplots, nrow = length(legendplots) + 1, ncol = 1, heights = c(7, rep(1, length(legendplots))))


  ##simuplot
  return(finalplot)
}

plotBaseHM = function(toplot, yLogScale, VirPalOption, ...){

  simuplot = ggplot2::ggplot(toplot, aes_string(x = "time", y = "Components", fill = "mean")) + ggplot2::geom_tile()

  if(yLogScale){
    simuplot = simuplot + ggplot2::scale_fill_viridis_c(trans = "log10", option = VirPalOption, name = "log10(Components absolute abundance+0.5)")
  }else{
    simuplot = simuplot + ggplot2::scale_fill_viridis_c(option = VirPalOption, name = "Components absolute abundance")
  }

  simuplot = simuplot +
    ggplot2::facet_grid(Type~Ind, scales = "free_y") + #, space = "free_y"
    ggplot2::xlab("Time (s)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.text = element_text(size = 7), strip.text = element_text(size = 10),
                   legend.spacing = grid::unit(5, "mm"), legend.key.size = grid::unit(10, "mm"),
                   legend.direction = "horizontal", ...)

  return(simuplot)

}

#' Plots the result of a simulation as a heatmap.
#'
#' Automatically plots the result of a simulation for the selected in silico individuals as a heatmap.
#'
#' If more than one trial is to be plotted, the mean abundance of each molecule over the different trials is plotted with a solid line,
#' and the min and max abundances represented as coloured areas around the mean.
#'
#' @param simdf The dataframe with the result of the simulation (see \code{\link{simulateInSilicoSystem}}).
#' @param molecules A vector of gene IDs (numeric or character) and/or complex IDs (e.g. CTC1) to be plotted.
#' @param inds A vector of in silico individual names for which to plot the expression profiles.
#' @param trials A vector of trials ID (= number) to use for the plot (see details).
#' @param timeMin Numeric. The minimum simulation time to plot. Default value set to the minimum time in the simulation.
#' @param timeMax Numeric. The maximum simulation time to plot. Default value set to the maximum time in the simulation.
#' @param mergeAllele Are the gene products originating from different alleles merged? Default TRUE. Also see \code{\link{mergeAlleleAbundance}}
#' @param mergePTM Are the modified and non-modified versions of the proteins merged? Default TRUE. Also see \code{\link{mergePTMAbundance}}
#' @param mergeComplexes Are the free and in complex gene products merged? Default FALSE. Also see \code{\link{mergeComplexesAbundance}}
#' @param yLogScale Plot the y-axis in log10-scale? If so, the abundance of each species at each time-point is increased by 1 to avoid zero values. Default TRUE.
#' @param nIndPerRow Positive integer, the number of individuals to plot per row. Default 3.
#' @param VirPalOption String, palette name option to be passed to \code{\link[ggplot2]{scale_fill_viridis_c}}; can be one of "magma", "inferno", "plasma", "viridis" or "cividis". Default value is "plasma".
#' @param ... Any additional parameter to be passed to \code{\link[ggplot2]{theme}} for the plot of each individual.
#' @return A plot from \code{\link[ggpubr]{ggarrange}}.
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, ploidy = 2)
#' mypop = createInSilicoPopulation(10, mysystem)
#' sim = simulateInSilicoSystem(mysystem, mypop, 100, ntrials = 5)
#' plotHeatMap(sim$Simulation,
#'  c(1, 2, 3),
#'  c("Ind1", "Ind2", "Ind3", "Ind4"),
#'  axis.title = element_text(color = "red"))
#' }
#' @export
plotHeatMap = function(simdf, molecules = NULL, inds = unique(simdf$Ind), trials = unique(simdf$trial), timeMin = min(simdf$time), timeMax = max(simdf$time), mergeAllele = T, mergePTM = T, mergeComplexes = F, yLogScale  = T, nIndPerRow = 3, VirPalOption = "plasma", ...){
  simind = simdf %>% dplyr::filter(!!sym("Ind") %in% inds) %>%
                     dplyr::filter(!!sym("trial") %in% trials) %>%
                     dplyr::filter(!!sym("time") >= timeMin & !!sym("time") <= timeMax)

  if(mergePTM)  simind = mergePTMAbundance(simind)
  if(mergeComplexes)  simind = mergeComplexesAbundance(simind)
  if(mergeAllele) simind = mergeAlleleAbundance(simind)

  ## First transformation into a long tibble
  toplot = simind %>%
    tidyr::gather(key = "ID", value = "Abundance", setdiff(names(simind), c("time", "trial", "Ind"))) %>%
    mutate("Type" = case_when(stringr::str_detect(.data$ID, "^R") ~ "RNAs",
                              stringr::str_detect(.data$ID, "^P") ~ "Proteins",
                              stringr::str_detect(.data$ID, "^C") ~ "Complexes"),
           "Components" = stringr::str_replace(.data$ID, "^R|^P", ""),
           "Components" = stringr::str_replace(.data$Components, "^m", "PTM"))

  ## Select the molecules to plot
  if(!is.null(molecules)){
    toplot = toplot %>%
      mutate("Molecule" = case_when(!!sym("Type") %in% c("RNAs", "Proteins") ~ stringr::str_extract(.data$Components, "(?<=^|(PTM))\\d+"),
                                    !!sym("Type") == "Complexes" ~ stringr::str_extract(.data$Components, "^C[[:alpha:]]{2,3}\\d+"))) %>%
      dplyr::filter(!!sym("Molecule") %in% molecules) %>%
      select(-!!sym("Molecule"))
  }

  ## If plot in log-scale, need to have non-zero abundances
  if(yLogScale) toplot = toplot %>% mutate("Abundance" = !!sym("Abundance") + 0.5)

  ## Only plot the mean (even if only one trial to be plotted)
  toplot = toplot %>% dplyr::group_by(!!sym("Ind"), !!sym("time"), !!sym("Components"), !!sym("Type"), !!sym("ID")) %>%
           dplyr::summarise("mean" = mean(!!sym("Abundance")))


  ## Create a plot for each individual
  plots = vector("list", length = length(inds))
  for(i in 1:length(inds)){
    plots[[i]] = plotBaseHM(filter(toplot, !!sym("Ind") == inds[i]), yLogScale, VirPalOption, ...)
  }

  ## How many rows of plot given the nb of inds to plot and nb of inds per row
  nbrows = length(inds) %/% nIndPerRow + (length(inds) %% nIndPerRow != 0)

  simuplot = ggpubr::ggarrange(plotlist = plots, nrow = nbrows, ncol = min(nIndPerRow, length(inds)), common.legend = T, heights = rep(1, length(plots)), legend = "bottom")

  ##simuplot
  return(simuplot)
}

#' Returns a summary data-frame of a simulation.
#'
#' Returns a summary data-frame of a simulation giving the maximum average abundance of each component over the different trials and the average abundance of the components at the final time of the simulation, for the selected in silico individuals.
#'
#' @param simdf The data frame with the result of the simulation (see \code{\link{simulateInSilicoSystem}}).
#' @param inds A vector of in silico individual names for which to compute the summary values.
#' @param trials A vector of trials ID (= number) to use for the summary.
#' @param timeMin Numeric. The minimum simulation time to take into account. Default value set to the minimum time in the simulation.
#' @param timeMax Numeric. The maximum simulation time to take into account. Default value set to the maximum time in the simulation.
#' @param mergeAllele Are the gene products originating from different alleles merged? Default TRUE. Also see \code{\link{mergeAlleleAbundance}}
#' @param mergePTM Are the modified and non-modified versions of the proteins merged? Default TRUE. Also see \code{\link{mergePTMAbundance}}
#' @param mergeComplexes Are the free and in complex gene products merged? Default FALSE. Also see \code{\link{mergeComplexesAbundance}}
#' @param verbose If TRUE (default), print the individuals, trials, min and max time considered for the computation of the summary.
#' @return A data-frame giving for each component (rows) and each individual (columns) the max and final average abundance over the different trials.
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none", ploidy = 2)
#' mypop = createInSilicoPopulation(15, mysystem)
#' sim = simulateInSilicoSystem(mysystem, mypop, 100, ntrials = 5)
#' summariseSimulation(sim$Simulation, c("Ind1", "Ind2", "Ind3", "Ind4"))
#' }
#' @export
summariseSimulation = function(simdf, inds = unique(simdf$Ind), trials = unique(simdf$trial), timeMin = min(simdf$time), timeMax = max(simdf$time), mergeAllele = T, mergePTM = T, mergeComplexes = F, verbose = T){

  simind = simdf %>% dplyr::filter(!!sym("Ind") %in% inds) %>%
    dplyr::filter(!!sym("trial") %in% trials) %>%
    dplyr::filter(!!sym("time") >= timeMin & !!sym("time") <= timeMax)

  if(mergePTM)  simind = mergePTMAbundance(simind)
  if(mergeComplexes)  simind = mergeComplexesAbundance(simind)
  if(mergeAllele) simind = mergeAlleleAbundance(simind)

  sumtable = simind %>%
    tidyr::gather(key = "ID", value = "Abundance", setdiff(names(simind), c("time", "trial", "Ind"))) %>%
    dplyr::group_by(!!sym("Ind"), !!sym("ID"), !!sym("time")) %>%
    dplyr::summarise("mean" = mean(!!sym("Abundance")))

  sumtable_max = sumtable  %>%
    dplyr::summarise("Max" = max(!!sym("mean")))

  sumtable_final = sumtable %>%
    dplyr::filter(!!sym("time") == timeMax) %>%
    dplyr::rename("Final" = !!sym("mean")) %>%
    dplyr::select(-time)

  sumtable = dplyr::full_join(sumtable_max, sumtable_final, by = c("Ind", "ID")) %>%
    tidyr::gather(key = "Abundance", value = "val", "Max", "Final") %>%
    tidyr::spread(!!sym("Ind"), !!sym("val")) %>%
    dplyr::rename("Components" = !!sym("ID")) %>%
    dplyr::arrange(!!sym("Components"), dplyr::desc(!!sym("Abundance")))

  ## Sort values by component name
  sortedComp = sortComponents(unique(sumtable$Components))
  sumtable = dplyr::left_join(data.frame(Components =sortedComp$Components, stringsAsFactors = F),sumtable, by="Components")

  if(verbose){
    cat("--------------------------\n")
    cat("Summary of simulation for:\n")
    cat("--------------------------\n")
    cat("Individuals:", inds, "\n")
    cat("Trials:", trials, "\n")
    cat("Time:", timeMin, "s -", timeMax, "s\n")
    cat("--------------------------\n")
  }
  return(sumtable)
}


#' Samples the expected library size of individuals/samples
#'
#' Samples the expected library size of each individual/sample, accounting for potential lane size effects
#' (i.e. the impact of samples being processed on different lanes).
#'
#' The expected library size of each individual is sampled from a log-normal distribution. The mean of this distribution
#' depends on the lane on which the individual/sample is processed. By default, when \code{laneEffect = FALSE}, all samples
#' are assumed to be processed in a single batch. Thus their library size is sampled from a log-normal distribution with
#' identical mean (equal to \code{meanLogLibSize_lane}) and sd \code{sdLogLibSize_samples}. If \code{laneEffect = TRUE}, the
#' samples are assumed to be processed in \code{nLanes} batches, that each have a different mean log-library size. In this
#' case, the mean of the log-normal distribution for each lane is sampled from a normal distribution with mean
#' \code{meanLogLibSize_lane} and sd \code{sdLogLibSize_lane}. In turn, the expected library size of each individual/sample
#' is sampled form a log-normal distribution with the corresponding lane-dependent mean, and sd \code{sdLogLibSize_samples}.
#'
#' @param samples_list List of sample/individual names.
#' @param meanLogLibSize_lane Numeric. The mean of the log10 mean library size normal distribution (see Details). Default value of 7.
#' @param sdLogLibSize_lane Numeric. The sd of the log10 mean library size normal distribution (see Details). Default value of 0.5.
#' @param sdLogLibSize_samples Numeric. The sd of the log10 samples library size normal distribution (see Details). Default value of 0.2.
#' @param laneEffect Boolean. Are the samples processed on different lanes/batches? Default value is FALSE.
#' @param nLanes Numeric. How many lanes are there in the experiment? Automatically set to 1 if \code{laneEffect = F}. Default value is 2.
#' @return A list:
#' \itemize{
#' \item \code{lane}: the lane on which each sample is processed.
#' \item \code{expected_library_size}: the expected library size of each sample.
#' \item \code{lane_mean_library_size}: the mean library size of each lane.
#' }
#' @examples
#' samples_list = sapply(1:10, function(x){paste0("Ind", x)})
#' libsize = sampleLibrarySize(samples_list)
#' libsize = sampleLibrarySize(samples_list, laneEffect = TRUE, nLanes = 3)
#' @export
sampleLibrarySize = function(samples_list, meanLogLibSize_lane = 7, sdLogLibSize_lane = 0.5, sdLogLibSize_samples = 0.2, laneEffect = F, nLanes = 2){

  if(!laneEffect) nLanes = 1

  N = length(samples_list)

  ## Creating lane mean log-library size
  if(laneEffect){
    samplesLane = sample(1:nLanes, N, replace = T) ## to which lane belongs each sample
    lanesMean = rnorm(nLanes, meanLogLibSize_lane, sdLogLibSize_lane) ## what is the mean log-library size of each lane
    samplesMean = lanesMean[samplesLane] ## what is the mean log-library size for each sample
  }else{
    samplesLane = rep(1, N)
    lanesMean = c(meanLogLibSize_lane)
    samplesMean = rep(lanesMean, N)
  }

  ## Sampling the log10 of library size
  logLibSize = sapply(1:N, function(x){rnorm(1, mean = samplesMean[x], sd = sdLogLibSize_samples)})

  libSize = 10^logLibSize
  names(samplesLane) = samples_list
  names(libSize) = samples_list
  names(lanesMean) = 1:nLanes

  return(list("lane" = samplesLane, "expected_library_size" = libSize, "lane_mean_library_size" = 10^lanesMean))
}


#' Transforms a simulation time-point into RNA-seq-like data.
#'
#' Transforms a time-point of a simulation into RNA-seq-like data, i.e. simulates a read count for each RNA molecule per
#' individual (each individual is considered as a sample).
#'
#' The abundance of the RNA form of each gene at time \code{samplingTime} is extracted from the result of the simulation. If \code{mrnasOnly = TRUE}, non-coding RNAs are discarded. If the simulation contains
#' several trials per individual (see \code{\link{simulateInSilicoSystem}}), the abundance of each RNA is summed over the different trials for each individual. The abundance of
#' the RNAs in each individual are then transformed into "noisy proportions": for each individual, \code{tot_RNAs} times \code{propRnasSampled} RNAs molecules are sampled from the
#' total RNA molecules in the individual (where \code{tot_RNAs} is the total number of RNA molecules of a given individual). If \code{propRnasSampled} is 1, this step simply
#' corresponds to dividing the abundance of each RNA by the total number of RNAs for the individual. Otherwise, if \code{propRnasSampled} is less than 1, it introduces some
#' stochasticity in the "measurement" of RNAs as some low-expressed RNAs might not be represented. The noisy proportion of for each RNA is multiplied by the gene length (by default
#' set to 1 for all genes), to reproduce the bias of RNA-seq experiments in which longer genes get more reads. The proportions are re-scaled such that their sum over all RNAs for
#' each individual equals 1. The expected count of each RNA in each individual is then computed as the product of the corresponding noisy proportion and the library size of the individual.
#' The latter can be provided by the user (parameter \code{samplesLibSize}); otherwise it is sampled using the function \code{\link{sampleLibrarySize}}. Finally, the actual read count of each RNA for
#' each individual is sampled from a Poisson distribution with parameter lambda equal to the corresponding expected count.
#'
#' @param simdf The data-frame with the result of the simulation (see \code{\link{simulateInSilicoSystem}}).
#' @param insilicosystem The simulated in silico system (see \code{\link{createInSilicoSystem}}).
#' @param samplingTime Numeric. Time-point of the simulation to be transformed. By default, the maximum time of the simulation is used.
#' @param laneEffect Boolean. Are the samples processed on different lanes/batches (see \code{\link{sampleLibrarySize}})? Ignored if \code{samplesLibSize} is provided. Default value is FALSE.
#' @param nLanes Numeric. How many lanes are there in the experiment (see \code{\link{sampleLibrarySize}})? Automatically set to 1 if \code{laneEffect = F}. Ignored if \code{samplesLibSize} is provided. Default value is 2.
#' @param propRnasSampled Numeric. The proportion of molecules of RNAs that are sampled in each individual. Must be between 0 and 1. Default value is 0.9.
#' @param samplesLibSize Vector of expected library size for each individual/sample. If named, the names of the vector must correspond to the names of the individuals as specified
#' in the result of the simulation. If none provided, will be sampled from a log-normal distribution (see \code{\link{sampleLibrarySize}}). Default value is NULL.
#' @param genesLength Vector of gene length for each gene in the system. Its length must be equal to the number of protein-coding genes in the system (if \code{mrnasOnly} is TRUE)
#' or of genes in the system (if \code{mrnasOnly} is FALSE). If named, the names must correspond to the ID of the (protein-coding) genes in the system. If none provided, all genes
#' are assumed to be of length 1.
#' @param mrnasOnly Boolean. Are the noncoding RNAs to be discarded before the transformation? If TRUE, read counts will be returned only for protein-coding RNAs. Default value is TRUE.
#' @param mergeComplexes Boolean. Are the RNAs in complex accounted for in the read counts (i.e. are they detected by the RNA-seq experiment)? Default value is FALSE. See also \code{\link{mergeComplexesAbundance}}.
#' @param meanLogLibSize_lane Numeric. The mean of the log10 mean library size normal distribution (see \code{\link{sampleLibrarySize}}). Ignored if \code{samplesLibSize} is provided. Default value of 7.
#' @param sdLogLibSize_lane Numeric. The sd of the log10 mean library size normal distribution (see \code{\link{sampleLibrarySize}}). Ignored if \code{samplesLibSize} is provided. Default value of 0.5.
#' @param sdLogLibSize_samples Numeric. The sd of the log10 samples library size normal distribution (see \code{\link{sampleLibrarySize}}). Ignored if \code{samplesLibSize} is provided. Default value of 0.2.
#' @return A list:
#' \itemize{
#' \item \code{rnaSeqMatrix} A tibble giving for each RNA (column "Molecule") the observed read count in each individual (other columns, one per individual).
#' \item \code{samplesLibSize} The expected library size of each individual. May not be equal to the total read counts for the individual, as the actual counts are sampled from a Poisson distribution.
#' \item \code{genesLength} The length of each gene.
#' }
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none",
#'                                 ploidy = 2, PC.p = 1)
#' mypop = createInSilicoPopulation(10, mysystem)
#' sim = simulateInSilicoSystem(mysystem, mypop, simtime = 1000,
#'                              ntrials = 10, nepochs = 5)
#' rnaSeq = getRNAseqMatrix(sim$Simulation, mysystem, laneEffect = F)
#'
#' ## With a batch/lane effect on the library size of the samples
#' rnaSeq = getRNAseqMatrix(sim$Simulation, mysystem, laneEffect = T)
#'
#' ## Providing the library size of each sample/individual
#' libsize = rnorm(length(unique(sim$Simulation$Ind)), 1e7, 1e5)
#' names(libsize) = unique(sim$Simulation$Ind)
#' rnaSeq = getRNAseqMatrix(sim$Simulation, mysystem, samplesLibSize = libsize)
#'
#' ## Accounting for different gene lengths
#' genes_length = sample(1:200, nrow(mysystem$genes))
#' names(genes_length) = as.character(1:nrow(mysystem$genes))
#' rnaSeq = getRNAseqMatrix(sim$Simulation, mysystem,
#'                          samplesLibSize = libsize, genesLength = genes_length)
#' }
#' @export
getRNAseqMatrix = function(simdf, insilicosystem, samplingTime = max(simdf$time), laneEffect = F, nLanes = 2, propRnasSampled = 0.9, samplesLibSize = NULL, genesLength = NULL, mrnasOnly = T, mergeComplexes = F, meanLogLibSize_lane = 7, sdLogLibSize_lane = 0.5, sdLogLibSize_samples = 0.2){

  samples_list = unique(simdf$Ind)

  if(propRnasSampled <= 0 | propRnasSampled > 1) stop("propRnasSampled must be between 0 and 1.")

  ## Compute the library size of each sample/individual
  if(!is.null(samplesLibSize)){
    if(length(samplesLibSize) != length(samples_list)){
      stop("Vector samplesLibSize length must be equal to the number of simulated individuals.")
    }
    if(is.null(names(samplesLibSize))){
      names(samplesLibSize) = samples_list
    }else if(length(setdiff(samples_list, names(samplesLibSize))) > 0){
      stop("The names of the samplesLibSize vector must correspond to the names of the individuals in the simulation.")
    }
    samplesLibSize_list = list("expected_library_size" = samplesLibSize)
  }else{
    samplesLibSize_list = sampleLibrarySize(samples_list,
                                            meanLogLibSize_lane = meanLogLibSize_lane,
                                            sdLogLibSize_lane = sdLogLibSize_lane,
                                            sdLogLibSize_samples = sdLogLibSize_samples,
                                            laneEffect = laneEffect,
                                            nLanes = nLanes)
    samplesLibSize = samplesLibSize_list$expected_library_size
  }

  ## Get the molecules to retain (all RNAs if mrnasOnly = F, otherwise only keep protein-coding RNAs)
  molstokeep = sapply(insilicosystem$genes$id, function(x){paste0("R", x)})
  if(mrnasOnly) molstokeep = molstokeep[insilicosystem$genes$coding == "PC"]

  if(length(molstokeep) == 0) stop("No gene applicable for the transformation; try with mrnasOnly = F?")

  ## Get the length of each gene/RNA
  if(!is.null(genesLength)){
    if(length(genesLength) != length(molstokeep)){
      stop("Vector genesLength length must be equal to the number of protein-coding genes in the system (if mrnasOnly is TRUE) or of genes in the system (if mrnasOnly is FALSE).")
    }
    if(is.null(names(genesLength))){
      names(genesLength) = stringr::str_replace(molstokeep, "R", "")
    }else if(length(setdiff(stringr::str_replace(molstokeep, "R", ""), names(genesLength))) > 0){
      stop("The names of the genesLength vector must correspond to the IDs of the genes in the system.")
    }
  }else{
    genesLength = rep(1, length(molstokeep))
    names(genesLength) = stringr::str_replace(molstokeep, "R", "")
  }

  ## Merge RNAs abundance from RNAs in complexes if necessary
  if(mergeComplexes) simdf = mergeComplexesAbundance(simdf)

  ## Get the absolute abundance of the RNAs
  df = mergeAlleleAbundance(simdf) %>%
    dplyr::filter(time == samplingTime) %>% ## get the RNA abundance at the sampling time
    dplyr::select(-!!sym("time")) %>%
    tidyr::gather(key = !!sym("Molecule"), value = !!sym("Abundance"), -!!sym("trial"), -!!sym("Ind")) %>%
    dplyr::filter(!!sym("Molecule") %in% molstokeep) %>%
    dplyr::group_by(!!sym("Ind"), !!sym("Molecule")) %>%
    dplyr::summarise(abundance = sum(!!sym("Abundance")))

  ## Transform absolute abundance of RNAs into noisy proportions
  totRNAperInd = df %>%
    dplyr::summarise(tot_rna = sum(!!sym("abundance")))

  noisy_prop = Reduce(dplyr::bind_rows, lapply(unique(df$Ind), function(x){
    df_ind = dplyr::filter(df, !!sym("Ind") == x)
    sampledRNAs = sample(rep(df_ind$Molecule, times = df_ind$abundance),
                         size = dplyr::filter(totRNAperInd, !!sym("Ind") == x)$tot_rna * propRnasSampled,
                         replace = F)
    dplyr::tibble(Ind = x,
           Molecule = df_ind$Molecule,
           noisy_prop = sapply(!!sym("Molecule"), function(i){sum(sampledRNAs == i) / length(sampledRNAs)}))
  })) %>%
    dplyr::mutate(gene_length = unname(genesLength[stringr::str_replace(!!sym("Molecule"), "R", "")]),
                  noisy_prop = !!sym("noisy_prop") * !!sym("gene_length"))

  tot_prop = noisy_prop %>%
    dplyr::group_by(!!sym("Ind")) %>%
    dplyr::summarise(tot = sum(!!sym("noisy_prop")))

  temp = tot_prop$tot
  names(temp) = tot_prop$Ind

  noisy_prop$noisy_prop = noisy_prop$noisy_prop / temp[noisy_prop$Ind]

  ## Compute the expected count of each molecule for each individual/sample
  expected_counts = noisy_prop %>%
    dplyr::mutate(library_size = unname(samplesLibSize[!!sym("Ind")]),
                  expected_counts = !!sym("noisy_prop") * !!sym("library_size"))

  ## Sample the actual counts per molecule per individual/sample
  rnaseq_counts = sapply(expected_counts$expected_counts, function(x){stats::rpois(1, x)})

  res = dplyr::tibble(Ind = expected_counts$Ind,
                      Molecule = expected_counts$Molecule,
                      read_counts = rnaseq_counts) %>%
    tidyr::spread(key = !!sym("Ind"), value = !!sym("read_counts"))

  return(list("rnaSeqMatrix" = res,
              "samplesLibSize" = samplesLibSize_list,
              "genesLength" = genesLength))
}
