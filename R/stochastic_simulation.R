#' Tranforms a data-frame into a list.
#'
#' Transforms a data-frame into a list. The elements of the lists correspond to
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
#' @param filepath If writefile = \code{TRUE}, path to the folder in the which the files will be created (default: current working directory).
#' @param filename If writefile = \code{TRUE}, prefix of the files created to store the lists of species and reactions (default: none).
#' @param ev A Julia evaluator (for the XRJulia). If none provided select the current evaluator or create one if no evaluator exists.
#' @return A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' @export
createStochSystem = function(insilicosystem, indargs, writefile, filepath = getwd(), filename = "simulation", ev = getJuliaEvaluator()){

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
                          complexes, complexeskinetics, indargs$gcnList,
                          writefile, filepath, filename, evaluator = ev)
  message("Done.")

  return(juliastochsystem)
}

#' Calls the Julia simulation
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

#' Simulates a in-silico system
#'
#' Simulates (stochastically) the behaviour of an in silico system over time.
#'
#' @param insilicosystem The in silico system to be simulated (see \code{\link{createInSilicoSystem}}).
#' @param insilicopopulation The in silico population to be simulated (see \code{\link{createInSilicoPopulation}}).
#' @param simtime The final time of the simulation (in seconds).
#' @param nepochs The number of times to record the state of the system during the simulation.
#' @param ntrials The number of times the simulation must be replicated.
#' @param simalgorithm The name of the simulation algorithm to use in the Julia function \code{simulate} from the module \code{BioSimulator}.
#' Can be one of "Direct", "FirstReaction", "NextReaction", "OptimizedDirect", "TauLeaping", "StepAnticipation".
#' @param writefile Does the julia function write the species and reactions lists in a text file?
#' @param filepath If writefile = \code{TRUE}, path to the folder in the which the files will be created (default: current working directory).
#' @param filename If writefile = \code{TRUE}, prefix of the files created to store the lists of species and reactions.
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return A list composed of:
#' \itemize{
#' \item \code{Simulation}: A data-frame with the simulated expression profiles of the genes for the different individuals in the in silico population. For gene i, "Ri" corresponds to the
#' RNA form of the gene, "Pi" to the protein form of the gene. "GCNj" corresponds to molecules (RNAs or proteins) originated from the j-th allele of the gene.
#' \item \code{runningtime}: A vector of running time of all simulations.
#' \item \code{stochmodel}: A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' }
#' @examples
#' \dontrun{
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none")
#' mypop = createInSilicoPopulation(1, mysystem, ploidy = 2)
#' sim = simulateInSilicoSystem(mysystem, mypop, 1000)
#' ## For visualisation
#' library(tidyr)
#' library(dplyr)
#' library(ggplot2)
#' library(RColorBrewer)
#' cols = brewer.pal(10, "Paired")
#' names(cols) = unlist(sapply(1:5, function(x){c(paste0(x, "GCN1"), paste0(x, "GCN2"))}))
#' toplot = sim$Simulation %>%
#' gather(key = "ID", value = "Abundance",
#'  setdiff(names(sim$Simulation), c("time", "trial", "Ind"))) %>%
#' mutate(Type = case_when(str_detect(ID, "^R") ~ "RNAs",
#'                          str_detect(ID, "^P") ~ "Proteins"),
#'         Components = stringr::str_replace(ID, "^R|^P", ""))
#' ggplot(toplot, aes(x = time, y = Abundance, colour = Components)) +
#' geom_line() + facet_grid(Type~., scales = "free_y") +
#' scale_colour_manual(values = cols, breaks = names(cols))
#' }
#' @export
simulateInSilicoSystem = function(insilicosystem, insilicopopulation, simtime, nepochs = -1, ntrials = 1, simalgorithm = "Direct", writefile = F, filepath = getwd(), filename = "simulation", ev = getJuliaEvaluator()){

  stochmodel = createStochSystem(insilicosystem, insilicopopulation$indargs, writefile, filepath, filename, ev = ev)
  message("\n")

  ## Store the running time of each simulation
  runningtime = vector("numeric", length(insilicopopulation$individualsList)*ntrials)
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

  message("\nMean running time per simulation: ", mean(runningtime),"seconds. \n")
  return(list("Simulation" = res, "runningtime" = runningtime, "stochmodel" = stochmodel))
}


## function to start a Julia evaluator on a node of the cluster, given the port id (for parallel simulation)
startJuliaEvCluster = function(portid){
  myev = newJuliaEvaluator(port = portid) ## start on the node a Julia evaluator with specified port number
  mystochmodel = juliaEval("eval(Meta.parse(%s))", stochmodel_string, .get = F, evaluator = myev) ## create in the Julia process the stochmodel object
  return(list("myev" = myev, "mystochmodelvar" = mystochmodel@.Data)) ## return the Julia evaluator ID and the name on the Julia process of the stochmodel object
}

## function to run a simulation on a cluster (for parallel simulation)
simulateInCluster = function(i, indtosimulate, ntrialstosimulate, increment, individualsList, genes, simtime, nepochs, simalgorithm){
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
  return(simJulia %>% mutate(trial = trial + increment[i], Ind = ind))
}


#' Simulates a in-silico system in parallel
#'
#' Simulates (stochastically) the behaviour of an in silico system over time using parallelisation.
#'
#' @param insilicosystem The in silico system to be simulated (see \code{\link{createInSilicoSystem}}).
#' @param insilicopopulation The in silico population to be simulated (see \code{\link{createInSilicoPopulation}}).
#' @param simtime The final time of the simulation (in seconds).
#' @param nepochs The number of times to record the state of the system during the simulation.
#' @param ntrials The number of times the simulation must be replicated.
#' @param simalgorithm The name of the simulation algorithm to use in the Julia function \code{simulate} from the module \code{BioSimulator}.
#' Can be one of "Direct", "FirstReaction", "NextReaction", "OptimizedDirect", "TauLeaping", "StepAnticipation".
#' @param writefile Does the julia function write the species and reactions lists in a text file?
#' @param filepath If writefile = \code{TRUE}, path to the folder in the which the files will be created (default: current working directory).
#' @param filename If writefile = \code{TRUE}, prefix of the files created to store the lists of species and reactions.
#' @param no_cores The number of cores to use for the simulation. By default use the function \code{detectCores} from the \code{parallel}
#' package to detect the number of available cores, and use this number-1 for the simulation.
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return A list composed of:
#' \itemize{
#' \item \code{Simulation}: A data-frame with the simulated expression profiles of the genes for the different individuals in the in silico population. For gene i, "Ri" corresponds to the
#' RNA form of the gene, "Pi" to the protein form of the gene. "GCNj" corresponds to molecules (RNAs or proteins) originated from the j-th allele of the gene.
#' \item \code{runningtime}: The running time (elapsed seconds) of the parallel simulation (only 1 value).
#' \item \code{stochmodel}: A Julia proxy object to retrieve the stochastic system in the Julia evaluator.
#' }
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5, regcomplexes = "none")
#' mypop = createInSilicoPopulation(15, mysystem, ploidy = 2)
#' sim = simulateParallelInSilicoSystem(mysystem, mypop, 1000)
#' }
#' @export
simulateParallelInSilicoSystem= function(insilicosystem, insilicopopulation, simtime, nepochs = -1, ntrials = 1, simalgorithm = "Direct", writefile = F, filepath = getwd(), filename = "simulation", no_cores = parallel::detectCores()-1, ev = getJuliaEvaluator()){

  stochmodel = createStochSystem(insilicosystem, insilicopopulation$indargs, writefile, filepath, filename, ev = ev)
  message("\n")

  stochmodel_string = juliaEval("string(%s)", stochmodel, evaluator = ev)

  mybaseport = ev$port ## Get the port of the current Julia evaluator
  portList = sapply(1:no_cores, sum, mybaseport) ## Assign to each core a port number, starting from 1+port number of the current evaluator

  mycluster = parallel::makeCluster(no_cores, outfile = "")

  parallel::clusterEvalQ(mycluster, library(XRJulia))
  parallel::clusterEvalQ(mycluster, library(utils))
  parallel::clusterExport(mycluster, "newJuliaEvaluator")
  parallel::clusterExport(mycluster, "callJuliaStochasticSimulation")
  parallel::clusterExport(mycluster, "stochmodel_string", envir = environment())
  parallel::clusterExport(mycluster, "no_cores", envir = environment())

  ## Start a Julia evaluator on each cluster
  message("Starting Julia evaluators on each cluster node ... \n")
  infocores = parallel::clusterApply(mycluster, portList, startJuliaEvCluster)
  message("Done.\n")
  parallel::clusterExport(mycluster, "infocores", envir = environment())

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
  parallel::clusterExport(mycluster, "progress", envir = environment())
  parallel::clusterExport(mycluster, "maxprogress", envir = environment())

  ## Run the simulation on the cluster
  startsim = tic()
  resTable = parallel::clusterApply(mycluster, 1:length(indtosimulate), simulateInCluster, indtosimulate = indtosimulate, ntrialstosimulate = ntrialstosimulate, increment = increment, individualsList = insilicopopulation$individualsList, genes = insilicosystem$genes, simtime = simtime, nepochs = nepochs, simalgorithm = simalgorithm)
  stopsim = toc(quiet = T)$toc
  res = bind_rows(resTable)
  message("\nRunning time of parallel simulations: ", stopsim - startsim, " seconds\n")

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
#' Merge (i.e. sum) the abundance of the different allelic versions of each molecule.
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
  mergeddf = df %>% select(time, trial, Ind)
  molsGCN = colnames(df)
  mols = stringr::str_replace_all(molsGCN, "GCN[[:digit:]]+", "")

  for(m in setdiff(unique(mols), c("time", "trial", "Ind"))){
    mergeddf[[m]] = sumColAbundance(df,which(mols == m))
  }

  return(mergeddf)
}

#' Merge the original and PTM versions of the proteins.
#'
#' Merge (i.e. sum) the abundance of the original and modified (PTM) versions of each protein.
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
  mergeddf = df %>% select(time, trial, Ind)
  molsPTM = colnames(df)
  mols = stringr::str_replace(molsPTM, "^Pm", "P")

  for(m in setdiff(unique(mols), c("time", "trial", "Ind"))){
    mergeddf[[m]] = sumColAbundance(df,which(mols == m))
  }

  return(mergeddf)
}

#' Merge the free and in-complex versions of molecules.
#'
#' Merge (i.e. sum) the abundance of the free and in-complex versions of each molecule.
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
    mergeddf[, splt[-1]] = mergeddf[, splt[-1]] + df[, i]
  }

  return(mergeddf)
}


# plotIndividual = function(simdf, ind, mergeAllele = F, mergePTM = F, mergeComplexes = F){
#
#   simind = simdf %>% filter(Ind == ind)
#
#   if(mergePTM)  simind = mergePTMAbundance(simind)
#   if(mergeComplexes)  simind = mergeComplexesAbundance(simind)
#   if(mergeAllele) simind = mergeAlleleAbundance(simind)
#
#   toplot = simind %>%
#     tidyr::gather(key = "ID", value = "Abundance", setdiff(names(simind), c("time", "trial", "Ind"))) %>%
#     mutate(Type = case_when(stringr::str_detect(ID, "^R") ~ "RNAs",
#                             stringr::str_detect(ID, "^P") ~ "Proteins",
#                             stringr::str_detect(ID, "^C") ~ "Complexes"),
#            Components = stringr::str_replace(ID, "^R|^P", ""))
#
#   ## GENERATE COLOUR PALETTE
#   if(!mergeAllele){
#     ## we need to see how many different components exist in the system and how many allelic version of each exist
#     compocols = data.frame(Components = unique(toplot$Components)) %>%
#       mutate(isComplex = stringr::str_detect(Components, "^C")) %>%
#       mutate(compoID = stringr::str_replace_all(Components, "GCN.+|_.+", ""),
#              Allele = stringr::str_replace(Components, "^[[:digit:]]+(?=GCN)|^[^_]+_{1}", "")) %>%
#       dplyr::arrange(compoID, Allele)
#
#     ncols = length(unique(compocols$compoID))
#     paletteCompo = randomcoloR::distinctColorPalette(ncols) ## create a colour for each component
#     names(paletteCompo) = unique(compocols$compoID)
#     palette = unlist(lapply(names(paletteCompo), function(x){ ## create a gradient (one colour for each allelic version of the component)
#       nall = sum(compocols$compoID == x)
#       return(grDevices::colorRampPalette(c("#FFFFFF", paletteCompo[x]))(nall+1)[-1])
#     }))
#     names(palette) = sort(compocols$Components)
#   }else{
#     palette = randomcoloR::distinctColorPalette(length(unique(toplot$Components)))
#     names(palette) = sort(unique(toplot$Components))
#   }
#
#
#   if(max(simind$trial) == 1){   ## If there is only one trial
#     simuplot = ggplot2::ggplot(toplot, aes(x = time, y = Abundance, colour = Components)) + ggplot2::geom_line() +
#       ggplot2::facet_grid(Type~., scales = "free_y") +
#       ggplot2::scale_colour_manual(values = palette, breaks = names(palette)) +
#       ggplot2::xlab("Time (s)") + ggplot2::ylab("Components absolute abundance") +
#       ggplot2::theme_minimal() +
#       ggplot2::theme(legend.text = element_text(size = 7), strip.text = element_text(size = 10))
#   }else{
#     toplotTrials = toplot %>%
#       dplyr::group_by(Ind, time, Components, Type, ID) %>%
#       dplyr::summarise("mean" = mean(Abundance), "LB" = min(Abundance), "UB" = max(Abundance))
#
#     simuplot = ggplot2::ggplot(toplotTrials, aes(x = time)) +
#       ggplot2::geom_ribbon(aes(ymin = LB, ymax = UB, fill = Components), alpha = 0.5) +
#       ggplot2::geom_line(aes(y = mean, colour = Components)) +
#       ggplot2::scale_colour_manual(values = palette, breaks = names(palette)) +
#       ggplot2::scale_fill_manual(values = palette, breaks = names(palette)) +
#       ggplot2::facet_grid(Type~., scales = "free_y") +
#       ggplot2::xlab("Time (s)") + ggplot2::ylab("Components absolute abundance") +
#       ggplot2::theme_minimal() +
#       ggplot2::theme(legend.text = element_text(size = 7), strip.text = element_text(size = 10))
#   }
#
#   print(simuplot)
#
# }
