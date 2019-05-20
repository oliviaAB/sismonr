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
#' @param indargs An object of class \code{insilicoindividualargs} (i.e. a list with parameters for in silico individuals generation).
#' @param writefile Does the julia function write the species and reactions lists in a text file?
#' @param filepath If writefile = \code{TRUE}, path to the folder in which the files will be created (default: current working directory).
#' @param filename If writefile = \code{TRUE}, prefix of the files created to store the lists of species and reactions (default: none).
#' @param verbose If TRUE (default), print messages to signal the start and finish of the function.
#' @param ev A Julia evaluator (for the XRJulia). If none provided select the current evaluator or create one if no evaluator exists.
#' @return A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5)
#' indargs = insilicoindividualargs()
#' stochsys = createStochSystem(mysystem, indargs)
#' }
#' @export
createStochSystem = function(insilicosystem, indargs, writefile = F, filepath = NULL, filename = "simulation", verbose = T, ev = getJuliaEvaluator()){

  if(is.null(filepath)) writefile = F; filepath = character(0)

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
                          complexes, complexeskinetics, indargs$gcnList,
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
#' @param InitVar The list of initial abundance variation coefficients of the in silico individual to be simulated (see \code{\link{createIndividual}}).
#' @param genes The data-frame of genes in the system.
#' @param simtime Numeric. The amount of time to simulate the model (in seconds).
#' @param modelname String. The name of the model. Default value "MySimulation".
#' @param ntrials Integer. The number of times the simulation must be replicated.
#' @param nepochs Integer. The number of times to record the state of the system during the simulation.
#' @param simalgorithm String. The name of the simulation algorithm to use in the Julia function \code{simulate} from the module \code{BioSimulator}.
#' Can be one of "Direct", "FirstReaction", "NextReaction", "OptimizedDirect", "TauLeaping", "StepAnticipation".
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return The result of the simulation (a data-frame).
#' @export
callJuliaStochasticSimulation = function(stochmodel, QTLeffects, InitVar, genes, simtime, modelname = "MySimulation", ntrials, nepochs, simalgorithm, ev = getJuliaEvaluator()){

  genesdf = df2list(genes)
  evXR = XR::getInterface(getClass("JuliaInterface"))
  expr = gettextf("%s(%s)","juliaStochasticSimulation", evXR$ServerArglist(stochmodel, QTLeffects, InitVar, genesdf,
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
#' Simulates (stochastically) the behaviour of an in silico system over time, i.e. the expression of the different genes.
#'
#' @param insilicosystem The in silico system to be simulated (see \code{\link{createInSilicoSystem}}).
#' @param insilicopopulation The in silico population to be simulated (see \code{\link{createInSilicoPopulation}}).
#' @param simtime The final time of the simulation (in seconds).
#' @param nepochs The number of times to record the state of the system during the simulation.
#' @param ntrials The number of times the simulation must be replicated (for each individual).
#' @param simalgorithm The name of the simulation algorithm to use in the Julia function \code{simulate} from the module \code{BioSimulator}.
#' Can be one of "Direct", "FirstReaction", "NextReaction", "OptimizedDirect", "TauLeaping", "StepAnticipation".
#' @param writefile Does the julia function write the species and reactions lists in a text file?
#' @param filepath If writefile = \code{TRUE}, path to the folder in which the files will be created (default: current working directory).
#' @param filename If writefile = \code{TRUE}, prefix of the files created to store the lists of species and reactions.
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return A list composed of:
#' \itemize{
#' \item \code{Simulation}: A data-frame with the simulated expression profiles of the genes for the different individuals in the in silico population. For gene i, "Ri" corresponds to the
#' RNA form of the gene, "Pi" to the protein form of the gene. The suffix "GCNj" indicates that the molecule comes from the j-th allele of the gene.
#' \item \code{runningtime}: A vector of running time of all runs of the simulation for each in silico individuals.
#' \item \code{stochmodel}: A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' }
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none")
#' mypop = createInSilicoPopulation(1, mysystem, ploidy = 2)
#' sim = simulateInSilicoSystem(mysystem, mypop, simtime = 1000, ntrials = 10)
#' head(sim$Simulation)
#' ## Visualising the result
#' plotSimulation(sim$Simulation)
#' }
#' @export
simulateInSilicoSystem = function(insilicosystem, insilicopopulation, simtime, nepochs = -1, ntrials = 1, simalgorithm = "Direct", writefile = F, filepath = NULL, filename = "simulation", ev = getJuliaEvaluator()){

  if(is.null(filepath)) writefile = F
  stochmodel = createStochSystem(insilicosystem, insilicopopulation$indargs, writefile, filepath, filename, verbose = T, ev = ev)
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
                                             insilicopopulation$individualsList[[ind]]$InitVar,
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
startJuliaEvCluster = function(portid, insilicosystem, indargs, writefile, filepath, filename){
  myev = newJuliaEvaluator(port = portid) ## start on the node a Julia evaluator with specified port number
  mystochmodel = createStochSystem(insilicosystem, indargs, writefile, filepath, filename, verbose = F, ev = myev)
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
                                           individualsList[[ind]]$InitVar,
                                           genes, simtime, modelname = ind, ntrials = ntrialsclus,
                                           nepochs = nepochs, simalgorithm = simalgorithm, myev)


  if(i %%no_cores == 1){
    utils::setTxtProgressBar(progress, min(i + (no_cores-1), maxprogress))
  }
  return(simJulia %>% mutate("trial" = !!sym("trial") + increment[i], "Ind" = ind))
}


#' Simulates an in silico system in parallel.
#'
#' Simulates (stochastically) the behaviour of an in silico system over time using parallelisation, i.e. the expression of the different genes.
#'
#' @param insilicosystem The in silico system to be simulated (see \code{\link{createInSilicoSystem}}).
#' @param insilicopopulation The in silico population to be simulated (see \code{\link{createInSilicoPopulation}}).
#' @param simtime The final time of the simulation (in seconds).
#' @param nepochs The number of times to record the state of the system during the simulation.
#' @param ntrials The number of times the simulation must be replicated (for each individual).
#' @param simalgorithm The name of the simulation algorithm to use in the Julia function \code{simulate} from the module \code{BioSimulator}.
#' Can be one of "Direct", "FirstReaction", "NextReaction", "OptimizedDirect", "TauLeaping", "StepAnticipation".
#' @param writefile Does the julia function write the species and reactions lists in a text file?
#' @param filepath If writefile = \code{TRUE}, path to the folder in which the files will be created (default: current working directory).
#' @param filename If writefile = \code{TRUE}, prefix of the files created to store the lists of species and reactions.
#' @param no_cores The number of cores to use for the simulation. By default use the function \code{detectCores} from the \code{parallel}
#' package to detect the number of available cores, and use this number - 1 for the simulation.
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return A list composed of:
#' \itemize{
#' \item \code{Simulation}: A data-frame with the simulated expression profiles of the genes for the different individuals in the in silico population. For gene i, "Ri" corresponds to the
#' RNA form of the gene, "Pi" to the protein form of the gene. The suffix "GCNj" indicates that the molecule comes from the j-th allele of the gene.
#' \item \code{runningtime}: The running time (elapsed seconds) of the parallel simulation (only 1 value).
#' \item \code{stochmodel}: A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' }
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none")
#' mypop = createInSilicoPopulation(15, mysystem, ploidy = 2)
#' sim = simulateParallelInSilicoSystem(mysystem, mypop, 1000)
#' head(sim$Simulation)
#' ## Visualising the result
#' plotSimulation(sim$Simulation)
#' }
#' @export
simulateParallelInSilicoSystem= function(insilicosystem, insilicopopulation, simtime, nepochs = -1, ntrials = 1, simalgorithm = "Direct", writefile = F, filepath = NULL, filename = "simulation", no_cores = parallel::detectCores()-1, ev = getJuliaEvaluator()){

  if(is.null(filepath)) writefile = F

  stochmodel = createStochSystem(insilicosystem, insilicopopulation$indargs, writefile, filepath, filename, verbose = T, ev = ev)
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
  infocores = parallel::clusterApply(mycluster, portList, startJuliaEvCluster, insilicosystem, insilicopopulation$indargs, writefile, filepath, filename)
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
#' mysystem = createInSilicoSystem(G = 5, empty = TRUE)
#' mypop = createInSilicoPopulation(1, mysystem, ploidy = 2)
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
#' mysystem = createInSilicoSystem(G = 5, PC.p = 1, PC.PTM.p = 0.9, regcomplexes = "none")
#' mypop = createInSilicoPopulation(1, mysystem, ploidy = 1)
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
#' mysystem = createInSilicoSystem(G = 5, PC.p = 1, PC.TC.p = 1)
#' mysystem = addComplex(mysystem, c(1, 2))
#' mypop = createInSilicoPopulation(1, mysystem, ploidy = 1)
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
  if(yLogScale) simuplot = simuplot + ggplot2::scale_y_log10(); plotTitle = "log10(Components absolute abundance+1)"

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
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none")
#' mypop = createInSilicoPopulation(15, mysystem, ploidy = 2)
#' sim = simulateParallelInSilicoSystem(mysystem, mypop, 100, ntrials = 5)
#' plotSimulation(sim$Simulation, c("Ind1", "Ind2", "Ind3", "Ind4"),
#'  axis.title = element_text(color = "red"))
#' }
#' @export
plotSimulation = function(simdf, inds = unique(simdf$Ind), trials = unique(simdf$trial), timeMin = min(simdf$time), timeMax = max(simdf$time), mergeAllele = T, mergePTM = T, mergeComplexes = F, yLogScale  = T, nIndPerRow = 3, nCompPerRow = 10, ...){

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

  ## If plot in log-scale, need to have non-zero abundances
  if(yLogScale) toplot = toplot %>% mutate("Abundance" = !!sym("Abundance") + 1)

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
    simuplot = simuplot + ggplot2::scale_fill_viridis_c(trans = "log10", option = VirPalOption, name = "log10(Components absolute abundance+1)")
  }else{
    simuplot = simuplot + ggplot2::scale_fill_viridis_c(option = VirPalOption, name = "Components absolute abundance")
  }

  simuplot = simuplot +
    ggplot2::facet_grid(Type~Ind, scales = "free_y") +
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
#' mysystem = createInSilicoSystem(G = 5)
#' mypop = createInSilicoPopulation(10, mysystem, ploidy = 2)
#' sim = simulateParallelInSilicoSystem(mysystem, mypop, 100, ntrials = 5)
#' plotHeatMap(sim$Simulation, c("Ind1", "Ind2", "Ind3", "Ind4"),
#'  axis.title = element_text(color = "red"))
#' }
#' @export
plotHeatMap = function(simdf, inds = unique(simdf$Ind), trials = unique(simdf$trial), timeMin = min(simdf$time), timeMax = max(simdf$time), mergeAllele = T, mergePTM = T, mergeComplexes = F, yLogScale  = T, nIndPerRow = 3, VirPalOption = "plasma", ...){
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

  ## If plot in log-scale, need to have non-zero abundances
  if(yLogScale) toplot = toplot %>% mutate("Abundance" = !!sym("Abundance") + 1)

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
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none")
#' mypop = createInSilicoPopulation(15, mysystem, ploidy = 2)
#' sim = simulateParallelInSilicoSystem(mysystem, mypop, 100, ntrials = 5)
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
