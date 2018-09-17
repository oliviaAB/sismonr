## Generate the genes in the system and their attributes, according to the user parameters
## Inputs:
##  - nod: data frame created by the function createGenes
##  - sysargs: an object of class insiliosystemargs, i.e. list of all parameters for in silico network generation
## Outputs:
##    - nod: a data frame of genes (rows) and their attributes

#' Creates genes for the in silico system
#'
#' Generate the genes in the system and their attributes, according to the user parameters.
#'
#' @param sysargs An object of class \code{insilicosystemargs} (i.e. a list with parameters for in silico system generation).
#' @return A data frame of in silico genes. Attributes:
#' \itemize{
#' \item \code{id}: Integer, ID of the genes;
#' \item \code{coding}: coding status of the genes (either "PC" for protein-coding or "NC" for non-coding);
#' \item \code{TargetReaction}: the biological function of the genes ("TC": transcription regulator, "TL": translation regulator, "RD": RNA decay
#' regulator, "PD": protein decay regulator, "PTM": post-translational modification regulator, "MR": metabolic enzyme);
#' \item \code{PTMform}: Does the gene have a PTM form? "0" or "1" (here all "0", PTM form will be assigned later)
#' \item \code{Active form}: what is the active form of the gene? "R" for noncoding genes, "P" for protein-coding genes,
#' "Pm" for protein-coding genes with a PTM form;
#' \item \code{TCrate}: transcription rate of the genes;
#' \item \code{TLrate}: translation rate of the genes (0 for noncoding genes);
#' \item \code{RDrate}: RNA decay rate of the genes;
#' \item \code{PDrate}: Protein decay rate of the genes (0 for noncoding genes).
#' }
#' @export
createGenes = function(sysargs){

  genes = data.frame("id" = 1:sysargs[["G"]], "coding" = rep("", sysargs[["G"]]), "TargetReaction" = rep("", sysargs[["G"]]),  "PTMform" = rep("0", sysargs[["G"]]), "ActiveForm" = rep("", sysargs[["G"]]),
                   "TCrate" = rep(0,sysargs[["G"]]), "TLrate" = rep(0,sysargs[["G"]]), "RDrate" = rep(0,sysargs[["G"]]), "PDrate" = rep(0,sysargs[["G"]]), stringsAsFactors = F)
  rownames(genes) = genes$id

  ## Deciding gene status
  genes$coding = sample(c("PC", "NC"), sysargs[["G"]], prob = c(sysargs[["PC.p"]], sysargs[["NC.p"]]), replace = T)

  ## Deciding gene function (reaction to be regulated)
  genes$TargetReaction[genes$coding == "PC"] = sample(c("TC", "TL", "RD", "PD", "PTM", "MR"), sum(genes$coding == "PC"), prob = c(sysargs[["PC.TC.p"]], sysargs[["PC.TL.p"]], sysargs[["PC.RD.p"]], sysargs[["PC.PD.p"]], sysargs[["PC.PTM.p"]], sysargs[["PC.MR.p"]]), replace = T)
  genes$TargetReaction[genes$coding == "NC"] = sample(c("TC", "TL", "RD", "PD", "PTM"), sum(genes$coding == "NC"), prob = c(sysargs[["NC.TC.p"]], sysargs[["NC.TL.p"]], sysargs[["NC.RD.p"]], sysargs[["NC.PD.p"]], sysargs[["NC.PTM.p"]]), replace = T)

  ## In genes, state what is the active form of each gene, i.e. which form (i.e. RNA, protein, activated protein) is performing the regulation
  genes$ActiveForm[genes$coding == "NC"] = "R" ## noncoding genes act through their RNA
  genes$ActiveForm[genes$coding == "PC"] = "P" ## protein-coding genes act through their protein
#  genes$ActiveForm = sapply(1:nrow(genes), function(x){paste0(genes$ActiveForm[x], genes$id[x])})

  ## Sample the kinetic parameters of the genes

  ## Transcription rate: applicable to all genes
  genes$TCrate = sysargs[["basal_transcription_rate_samplingfct"]](sysargs[["G"]])

  ## RNA decay rate: applicable to all genes
  ## Sample the lifetime, the decay rate is defined as 1/lifetime
  genes$RDrate = 1/sysargs[["basal_RNAlifetime_samplingfct"]](sysargs[["G"]])

  ## Translation rate: applicable to protein-coding genes
  genes$TLrate[genes$coding == "PC"] = sysargs[["basal_translation_rate_samplingfct"]](sum(genes$coding == "PC"))

  ## Protein coding rate: applicable to protein-coding genes
  genes$PDrate[genes$coding == "PC"] = 1/sysargs[["basal_protlifetime_samplingfct"]](sum(genes$coding == "PC"))

  return(genes)
}

#' Creates an in silico regulatory network
#'
#' Creates an in silico regulatory network given a list of regulators and targets
#'
#' @param regsList A named list of length 2. Element "PC" (resp."NC") is a vector of gene ID of the protein-coding (resp. noncoding) regulators
#' for the network.
#' @param tarsList A named list of length 2. Element "PC" (resp."NC") is a vector of gene ID of the targets of protein-coding (resp. noncoding)
#' regulators.
#' @param reaction String. The ID of the reaction ("TC", "TL", "RD", "PD" or "PTM")
#' @param sysargs An object of class \code{insilicosystemargs} (i.e. a list with parameters for in silico system generation).
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return A list of two elements:
#' \itemize{
#' \item \code{edg} a data-frame of edges of the networks;
#' \item \code{complexes} a list of complexes composition (each element is named with the complex ID, the components are given as gene IDs).
#' }
#' The \code{edg} data-frame has the following variables:
#' \itemize{
#' \item \code{from} gene ID of the regulator, as a character;
#' \item \code{to} gene ID of the target, as an integer;
#' \item \code{TargetReaction} the ID of the reaction (as given by \code{reaction});
#' \item \code{RegSign} The sign of the reaction ("1" or "-1");
#' \item \code{RegBy} Is the regulator a protein-coding gene ("PC"), a noncoding gene ("NC") or a complex ("C")?
#' }
#' @export
createRegulatoryNetwork = function(regsList, tarsList, reaction, sysargs, ev = getJuliaEvaluator()){

  ## Call the julia function nwgeneration to generate the regulatory network
  juliaedg = juliaGet(juliaCall("juliaCreateNetwork", reaction, regsList[["PC"]], tarsList[["PC"]], sysargs[[paste(reaction, "PC", "indeg.distr", sep = ".")]],
                             sysargs[[paste(reaction, "PC", "outdeg.distr", sep = ".")]], sysargs[[paste(reaction, "PC", "outdeg.exp", sep = ".")]],
                             sysargs[[paste(reaction, "PC", "autoregproba", sep = ".")]], sysargs[[paste(reaction, "PC", "twonodesloop", sep = ".")]],
                             regsList[["NC"]], tarsList[["NC"]], sysargs[[paste(reaction, "NC", "indeg.distr", sep = ".")]],
                             sysargs[[paste(reaction, "NC", "outdeg.distr", sep = ".")]], sysargs[[paste(reaction, "NC", "outdeg.exp", sep = ".")]],
                             sysargs[[paste(reaction, "NC", "autoregproba", sep = ".")]], sysargs[[paste(reaction, "NC", "twonodesloop", sep = ".")]],
                             sysargs[["regcomplexes"]], sysargs[["regcomplexes.size"]], sysargs[["regcomplexes.p"]],
                             evaluator = ev), evaluator = ev)

  if(nrow(juliaedg$edg) == 0){
    edg = data.frame("from" = integer(), "to" = integer(), "TargetReaction" = character(), "RegSign" = character(), "RegBy" = character(), stringsAsFactors = F)
  } else{
    edg = data.frame("from" = unlist(juliaedg$edg[,1]), "to" = unlist(juliaedg$edg[,2]), "TargetReaction" = rep(reaction, nrow(juliaedg$edg)), "RegSign" = rep("", nrow(juliaedg$edg)), "RegBy" = unlist(juliaedg$edg[,3]), stringsAsFactors = F)

    ## Sample the sign of the edges
    if(reaction == "PTM"){ ## we want to make sure that for each gene targeted by PTM regulation at least one edge is positive (i.e. the target
                           ## gets transformed into its modified form)
      realtars = unique(edg$to) ## gives the id of target genes in the constructed network
      tarsedg = lapply(realtars, function(x){which(edg$to == x)}) ## for each target gene returns the rows id of juliaedg$edg corresponding to regulatory edges targeting the gene
      tarsedgpos = sapply(tarsedg, function(x){x[1]}) ## for each target gene returns the row id corresponding to its 1st regulatory edge
      tarsedgother = unlist(sapply(tarsedg, function(x){x[-1]})) ## for each target gene returns the row ids of the next (if existing) regulatory edges

      edg$RegSign[tarsedgpos] = "1" ## for each target gene, the first regulatory edge is positive, meaning that the regulator turns the protein into its modified form
      if(length(tarsedgother) != 0) edg$RegSign[tarsedgother] = sample(c("1","-1"), length(tarsedgother), prob = c(sysargs[["PTM.pos.p"]], 1 - sysargs[["PTM.pos.p"]]), replace = T)

    }else{
      edg$RegSign = sample(c("1","-1"), nrow(edg), prob = c(sysargs[[paste(reaction, "pos.p", sep = ".")]], 1 - sysargs[[paste(reaction, "pos.p", sep = ".")]]), replace = T)
    }
    edg = edg[order(edg$to),]
  }

  rownames(edg) = NULL
  complexes = lapply(juliaedg$complexes, unlist)
  return(list("edg" = edg, "complexes" = complexes))
}


#' Creates an in silico multi omic network
#'
#' Creates an in silico multi omic network from a data-frame of genes in the system.
#'
#' The \code{edg} data-frame represents the edges in a regulatory network (1 row = 1 edge). It contains the following parameters:
#' \itemize{
#' \item \code{from}: gene ID of the origin of the edge (character).
#' \item \code{to}: gene ID of the destination of the edge (integer).
#' \item \code{TargetReaction}: Type of regulation (ID of the controlled reaction: "TC", "TL", "RD", "PD" or "PTM").
#' \item \code{RegSign}: Sign of the reaction ("1" for positive regulation, "-1" for negative regulation).
#' \item \code{RegBy}: Type of the regulator: "PC" for protein-coding regulation, "NC" for noncoding regulator,
#' "C" for regulatory complex.
#' }
#'
#' @param genes A data-frame of genes in the system (see \code{\link{createGenes}}).
#' @param sysargs An object of class \code{insilicosystemargs} (i.e. a list with parameters for in silico system generation).
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return An in silico multi omic network, that is a list of:
#' \itemize{
#' \item \code{genes}: the modified data-frame of genes;
#' \item \code{edg}: A data-frame of edges in the multi-omic network (see details below).
#' \item \code{mosystem}: A list of the different regulation networks in the system and associated information.
#' \item \code{complexes}: a list of regulatory complexes composition. The names of the elements are the IDs of the complexes, and the
#' values are vectors of gene IDs constituting each regulatory complex.
#' \item \code{complexeskinetics}: a list of regulatory complexes kinetic parameters.
#' }
#' The \code{mosystem} list is composed of \code{TCRN_edg}, \code{TLRN_edg}, \code{RDRN_edg}, \code{PDRN_edg} and \code{PTMRN_edg}: data-frames of edges for the different
#' regulatory networks, with in addition to the usual fields in the edg data frame contains columns for kinetic parameters of the
#' regulation.
#' @export
createMultiOmicNetwork = function(genes, sysargs, ev = getJuliaEvaluator()){

  complexes = list()

  ## Define transcription regulatory network (TCRN) ----

  ## Identify transcription regulators in the system
  PCreg_id = genes$id[genes$coding == "PC" & genes$TargetReaction == "TC"] ## protein-coding regulators
  NCreg_id = genes$id[genes$coding == "NC" & genes$TargetReaction == "TC"] ## noncoding regulators

  ## Identify targets of transcription regulators in the system - here any gene
  PCtarget_id = genes$id
  NCtarget_id = genes$id
#  NCtarget_id = genes$id[genes$coding == "PC"]

  ## Construct the regulatory network
  TCRN = createRegulatoryNetwork(regsList = list("PC" = PCreg_id, "NC" = NCreg_id),
                                 tarsList = list("PC" = PCtarget_id, "NC" = NCtarget_id),
                                 reaction = "TC", sysargs = sysargs, ev = ev)
  TCRN_edg = TCRN[["edg"]]
  complexes = c(complexes, TCRN[["complexes"]])

  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameters for transcription regulation include the binding and unbinding rate of regulators to gene promoter,
  ##    and the fold change induced on transcription rate by a regulator bound to the promoter
  TCRN_edg = data.frame(TCRN_edg, "TCbindingrate" = sysargs[["TCbindingrate_samplingfct"]](nrow(TCRN_edg)),
                                  "TCunbindingrate" = sysargs[["TCunbindingrate_samplingfct"]](nrow(TCRN_edg)),
                                  "TCfoldchange" = rep(0, nrow(TCRN_edg)), stringsAsFactors = F)
  TCRN_edg$TCfoldchange[TCRN_edg$RegSign == "1"] = sysargs[["TCfoldchange_samplingfct"]](sum(TCRN_edg$RegSign == "1"))  ## Repressors induce a fold change of 0

  ## Define translation regulatory network (TLRN) ----

  ## Identify translation regulators in the system
  PCreg_id = genes$id[genes$coding == "PC" & genes$TargetReaction == "TL"] ## protein-coding regulators
  NCreg_id = genes$id[genes$coding == "NC" & genes$TargetReaction == "TL"] ## noncoding regulators

  ## Identify targets of translation regulators in the system - here any protein-coding gene
  PCtarget_id = genes$id[genes$coding == "PC"]
  NCtarget_id = genes$id[genes$coding == "PC"]

  ## Construct the regulatory network
  TLRN = createRegulatoryNetwork(regsList = list("PC" = PCreg_id, "NC" = NCreg_id),
                                 tarsList = list("PC" = PCtarget_id, "NC" = NCtarget_id),
                                 reaction = "TL", sysargs = sysargs, ev = ev)
  TLRN_edg = TLRN[["edg"]]
  complexes = c(complexes, TLRN[["complexes"]])

  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameters for translation regulation include the binding and unbinding rate of regulators to gene promoter,
  ##    and the fold change induced on translation rate by a regulator bound to the promoter
  TLRN_edg = data.frame(TLRN_edg, "TLbindingrate" = sysargs[["TLbindingrate_samplingfct"]](nrow(TLRN_edg)),
                        "TLunbindingrate" = sysargs[["TLunbindingrate_samplingfct"]](nrow(TLRN_edg)),
                        "TLfoldchange" = rep(0, nrow(TLRN_edg)), stringsAsFactors = F)
  TLRN_edg$TLfoldchange[TLRN_edg$RegSign == "1"] = sysargs[["TLfoldchange_samplingfct"]](sum(TLRN_edg$RegSign == "1"))  ## Repressors induce a fold change of 0


  ## Define RNA decay regulatory network (RDRN) ----

  ## Identify RNA decay regulators in the system
  PCreg_id = genes$id[genes$coding == "PC" & genes$TargetReaction == "RD"] ## protein-coding regulators
  NCreg_id = genes$id[genes$coding == "NC" & genes$TargetReaction == "RD"] ## noncoding regulators

  ## Identify targets of RNA decay regulators in the system - here any gene
  PCtarget_id = genes$id
  NCtarget_id = genes$id

  ## Construct the regulatory network
  RDRN = createRegulatoryNetwork(regsList = list("PC" = PCreg_id, "NC" = NCreg_id),
                                 tarsList = list("PC" = PCtarget_id, "NC" = NCtarget_id),
                                 reaction = "RD", sysargs = sysargs, ev = ev)
  RDRN_edg = RDRN[["edg"]]
  complexes = c(complexes, RDRN[["complexes"]])

  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameter for RNA decay regulation is the decay rate induced by the regulators
  RDRN_edg = data.frame(RDRN_edg, "RDregrate" = sysargs[["RDregrate_samplingfct"]](nrow(RDRN_edg)), stringsAsFactors = F)


  ## Define protein decay regulatory network (PDRN) ----

  ## Identify protein decay regulators in the system
  PCreg_id = genes$id[genes$coding == "PC" & genes$TargetReaction == "PD"] ## protein-coding regulators
  NCreg_id = genes$id[genes$coding == "NC" & genes$TargetReaction == "PD"] ## noncoding regulators

  ## Identify targets of protein decay regulators in the system - here any protein-coding gene
  PCtarget_id = genes$id[genes$coding == "PC"]
  NCtarget_id = genes$id[genes$coding == "PC"]

  ## Construct the regulatory network
  PDRN = createRegulatoryNetwork(regsList = list("PC" = PCreg_id, "NC" = NCreg_id),
                                 tarsList = list("PC" = PCtarget_id, "NC" = NCtarget_id),
                                 reaction = "PD", sysargs = sysargs, ev = ev)
  PDRN_edg = PDRN[["edg"]]
  complexes = c(complexes, PDRN[["complexes"]])

  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameter for protein decay regulation is the decay rate induced by the regulators
  PDRN_edg = data.frame(PDRN_edg, "PDregrate" = sysargs[["PDregrate_samplingfct"]](nrow(PDRN_edg)), stringsAsFactors = F)


  ## Define protein post-translational modification regulatory network (PTMRN) ----

  ## Identify protein post-translational modification regulators in the system
  PCreg_id = genes$id[genes$coding == "PC" & genes$TargetReaction == "PTM"] ## protein-coding regulators
  NCreg_id = genes$id[genes$coding == "NC" & genes$TargetReaction == "PTM"] ## noncoding regulators

  ## Identify targets of protein post-translational modification in the system - here any protein-coding gene
  PCtarget_id = genes$id[genes$coding == "PC"]
  NCtarget_id = genes$id[genes$coding == "PC"]

  ## Construct the regulatory network
  PTMRN = createRegulatoryNetwork(regsList = list("PC" = PCreg_id, "NC" = NCreg_id),
                                 tarsList = list("PC" = PCtarget_id, "NC" = NCtarget_id),
                                 reaction = "PTM", sysargs = sysargs, ev = ev)
  PTMRN_edg = PTMRN[["edg"]]
  complexes = c(complexes, PTMRN[["complexes"]])

  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameter for protein post-translational modification regulation is the decay rate induced by the regulators
  PTMRN_edg = data.frame(PTMRN_edg, "PTMregrate" = sysargs[["PTMregrate_samplingfct"]](nrow(PTMRN_edg)), stringsAsFactors = F)

  ## Uptade the PTM status of genes ----
  PTMtargets = unique(PTMRN_edg$to)
  genes[PTMtargets, "PTMform"] = "1"
  genes[PTMtargets, "ActiveForm"] = "Pm"
  genes$ActiveForm = sapply(1:nrow(genes), function(x){paste0(genes$ActiveForm[x], genes$id[x])})

  ## Define regulatory complexes kinetic parameters ----

  complexeskinetics = list()
  if(length(complexes)>0){
    formrates = sysargs[["complexesformationrate_samplingfct"]](length(complexes))
    dissrates = sysargs[["complexesdissociationrate_samplingfct"]](length(complexes))
    for(c in 1:length(complexes)){
      complexeskinetics[[names(complexes)[c]]] = list("formationrate" = formrates[c], "dissociationrate" = dissrates[c])
    }}


  ## Return the in silico system ----

  edg = do.call(rbind, list(TCRN_edg[, c("from", "to", "TargetReaction", "RegSign", "RegBy")],
                            TLRN_edg[, c("from", "to", "TargetReaction", "RegSign", "RegBy")],
                            RDRN_edg[, c("from", "to", "TargetReaction", "RegSign", "RegBy")],
                            PDRN_edg[, c("from", "to", "TargetReaction", "RegSign", "RegBy")],
                            PTMRN_edg[, c("from", "to", "TargetReaction", "RegSign", "RegBy")]))

  ## Return
  res = list("TCRN_edg" = TCRN_edg,
             "TLRN_edg" = TLRN_edg,
             "RDRN_edg" = RDRN_edg,
             "PDRN_edg" = PDRN_edg,
             "PTMRN_edg" = PTMRN_edg)

  return(list("mosystem" = res, "genes" = genes, "edg" = edg, "complexes" = complexes, "complexeskinetics" = complexeskinetics))

}

#' Creates an empty in silico multi omic network
#'
#' Creates an empty in silico multi omic network from a data-frame of genes in the system (no regulatory interactions) (cf
#' \code{\link{createMultiOmicNetwork}}).
#'
#' @param genes A data-frame of genes in the system (see \code{\link{createGenes}}).
#' @return An in silico multi omic network, that is a list of:
#' \itemize{
#' \item \code{genes}: the modified data-frame of genes;
#' \item \code{edg}: A data-frame of edges in the multi-omic network (see details below).
#' \item \code{mosystem}: A list of the different regulation networks in the system and associated information.
#' \item \code{complexes}: a list of regulatory complexes composition. The names of the elements are the IDs of the complexes, and the
#' values are vectors of gene IDs constituting each regulatory complex. Empty list.
#' \item \code{complexeskinetics}: a list of regulatory complexes kinetic parameters. Empty list.
#' }
#' The \code{mosystem} list is composed of \code{TCRN_edg}, \code{TLRN_edg}, \code{RDRN_edg}, \code{PDRN_edg} and \code{PTMRN_edg}: data-frames of edges for the different
#' regulatory networks, with in addition to the usual fields in the edg data frame contains columns for kinetic parameters of the
#' regulation. All empty.
#' @export
createEmptyMultiOmicNetwork = function(genes){

  edg = data.frame("from" = integer(), "to" = integer(), "TargetReaction" = character(), "RegSign" = character(), "RegBy" = character(), stringsAsFactors = F)
  res = list("TCRN_edg" = data.frame(edg,  "TCbindingrate" = numeric(), "TCunbindingrate" = numeric(), "TCfoldchange" = numeric(), stringsAsFactors = F),
             "TLRN_edg" = data.frame(edg,  "TLbindingrate" = numeric(), "TLunbindingrate" = numeric(), "TLfoldchange" = numeric(), stringsAsFactors = F),
             "RDRN_edg" = data.frame(edg, "RDregrate" = numeric(), stringsAsFactors = F),
             "PDRN_edg" = data.frame(edg, "PDregrate" = numeric(), stringsAsFactors = F),
             "PTMRN_edg" = data.frame(edg, "PTMregrate" = numeric(), stringsAsFactors = F))

  return(list("mosystem" = res, "genes" = genes, "edg" = edg, "complexes" = list(), "complexeskinetics" = list()))
}

#' Creates an in silico system.
#'
#' Creates an in silico system, i.e. the genes and the regulatory networks defining the system.
#'
#' @param empty Logical. Does the regulatory network is empty (= no regulation)? Default value is \code{FALSE}.
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return An object of class \code{insilicosystem}, that is a list composed of:
#' \itemize{
#' \item \code{genes}: a data-frame of genes (see \code{\link{createGenes}});
#' \item \code{edg}: a data-frame of edges in the regulatory network (see \code{\link{createMultiOmicNetwork}});
#' \item \code{mosystem}: a list defining the multi-omic regulatory network (see \code{\link{createMultiOmicNetwork}});
#' \item \code{sysargs}: An object of class \code{insilicosystemargs}; the parameters used to create the system.
#' }
#' @export
createInSilicoSystem = function(empty = F, ev = getJuliaEvaluator(), ...){

  sysargs = insilicosystemargs(...)
  genes = createGenes(sysargs)

  if(empty){
    monw = createEmptyMultiOmicNetwork(genes)
  }else{
    monw = createMultiOmicNetwork(genes, sysargs, ev)
  }

  value = c(monw, list("sysargs" = sysargs))
  attr(value, "class") = "insilicosystem"

  return(value)
}


#' Add an edge in the in silico system.
#'
#' Add an edge in the in silico system between specified genes.
#'
#' @param insilicosystem The in silico system (see \code{\link{createInSilicoSystem}}).
#' @param regID Integer. The ID of the regulator gene.
#' @param tarID Integer. The ID of the target gene.
#' @param targetreaction The type of regulation exerted ("TC", "TL", "RD", "PD" or "PTM").
#' @param regsign The sign of the regulation: either "1" (positive regulation) or "-1" (negative regulation).
#' @param kinetics Optional: list of kinetics parameters of the reaction. If not provided, they will be sampled
#' according to the distributions specified in \code{sysargs} (parameters of the simulation - an attribute of the in silico system).
#' @return The modified in silico system.
#' @export
addEdg = function(insilicosystem, regID, tarID, targetreaction, regsign, kinetics = NULL){

  ## Checking the input values ----
  if(class(insilicosystem) != "insilicosystem"){
    stop("Argument insilicosystem must be of class \"insilicosystem\".")
  }

  if(!(regID %in% insilicosystem$genes$id)){
    stop("Regulator id does not exist in the system.")
  }

  if(!(tarID %in% insilicosystem$genes$id)){
    stop("Target id does not exist in the system.")
  }

  if(!(targetreaction %in% c("TC", "TL", "RD", "PD", "PTM"))){
    stop("Target reaction unknown.")
  }

  if(!(regsign %in% c("1", "-1"))){
    stop("Interation sign must be either \"1\" or \"-1\".")
  }

  regby = dplyr::filter(insilicosystem$genes, id == as.character(regID))[1,"coding"]

  ## Checking if an edge already exists between the 2 genes ----
  if(nrow(dplyr::filter(insilicosystem$edg, from == regID & to == tarID)) != 0){
    stop(paste0("An edge already exists from gene ", regID, " to gene ", tarID, "."))
  }

  ## Adding the interaction in the edg data-frame ----
  insilicosystem$edg = dplyr::add_row(insilicosystem$edg, from = as.character(regID), to = as.integer(tarID),
                 TargetReaction = targetreaction, RegSign = regsign,
                 RegBy = regby)

  if(targetreaction == "TC"){
    if(is.null(kinetics)){ ## if no values are given for the kinetic parameters
      myTCbindingrate = insilicosystem$sysargs[["TCbindingrate_samplingfct"]](1)
      myTCunbindingrate = insilicosystem$sysargs[["TCunbindingrate_samplingfct"]](1)
      myTCfoldchange = insilicosystem$sysargs[["TCfoldchange_samplingfct"]](1) * (regsign == "1") ## if regsign = "-1" (repression) the fold change is 0
    }else if(length(setdiff(c("TCbindingrate", "TCunbindingrate", "TCfoldchange"), names(kinetics))) == 0){ ## if values are given for the appropriate kinetic parameters
      myTCbindingrate = kinetics$TCbindingrate
      myTCunbindingrate = kinetics$TCunbindingrate
      myTCfoldchange = kinetics$TCfoldchange
    }else{ ## if the kinetic vector does not provide values with appropriate kinetic parameter names
      stop("Vector kinetics does not provide the correct parameters: TCbindingrate, TCunbindingrate, TCfoldchange.")
    }
    insilicosystem$mosystem$TCRN_edg = dplyr::add_row(insilicosystem$mosystem$TCRN_edg, from = as.character(regID), to = as.integer(tarID),
                   TargetReaction = targetreaction, RegSign = regsign,
                   RegBy = regby,
                   TCbindingrate = myTCbindingrate, TCunbindingrate = myTCunbindingrate, TCfoldchange = myTCfoldchange)
  }

  if(targetreaction == "TL"){
    if(is.null(kinetics)){ ## if no values are given for the kinetic parameters
      myTLbindingrate = insilicosystem$sysargs[["TLbindingrate_samplingfct"]](1)
      myTLunbindingrate = insilicosystem$sysargs[["TLunbindingrate_samplingfct"]](1)
      myTLfoldchange = insilicosystem$sysargs[["TLfoldchange_samplingfct"]](1) * (regsign == "1") ## if regsign = "-1" (repression) the fold change is 0
    }else if(length(setdiff(c("TLbindingrate", "TLunbindingrate", "TLfoldchange"), names(kinetics))) == 0){ ## if values are given for the appropriate kinetic parameters
      myTLbindingrate = kinetics$TLbindingrate
      myTLunbindingrate = kinetics$TLunbindingrate
      myTLfoldchange = kinetics$TLfoldchange
    }else{ ## if the kinetic vector does not provide values with appropriate kinetic parameter names
      stop("Vector kinetics does not provide the correct parameters: TLbindingrate, TLunbindingrate, TLfoldchange.")
    }
    insilicosystem$mosystem$TLRN_edg = dplyr::add_row(insilicosystem$mosystem$TLRN_edg, from = as.character(regID), to = as.integer(tarID),
                   TargetReaction = targetreaction, RegSign = regsign,
                   RegBy = regby,
                   TLbindingrate = myTLbindingrate, TLunbindingrate = myTLunbindingrate, TLfoldchange = myTLfoldchange)
  }

  if(targetreaction == "RD"){
    if(is.null(kinetics)){ ## if no values are given for the kinetic parameters
      myRDregrate = insilicosystem$sysargs[["RDregrate_samplingfct"]](1)
    }else if(length(setdiff(c("RDbindingrate"), names(kinetics))) == 0){ ## if values are given for the appropriate kinetic parameters
      myRDregrate = kinetics$RDregrate
    }else{ ## if the kinetic vector does not provide values with appropriate kinetic parameter names
      stop("Vector kinetics does not provide the correct parameters: RDbindingrate.")
    }
    insilicosystem$mosystem$RDRN_edg = dplyr::add_row(insilicosystem$mosystem$RDRN_edg, from = as.character(regID), to = as.integer(tarID),
                   TargetReaction = targetreaction, RegSign = regsign,
                   RegBy = regby,
                   RDregrate = myRDregrate)
  }

  if(targetreaction == "PD"){
    if(is.null(kinetics)){ ## if no values are given for the kinetic parameters
      myPDregrate = insilicosystem$sysargs[["PDregrate_samplingfct"]](1)
    }else if(length(setdiff(c("PDbindingrate"), names(kinetics))) == 0){ ## if values are given for the appropriate kinetic parameters
      myPDregrate = kinetics$PDregrate
    }else{ ## if the kinetic vector does not provide values with appropriate kinetic parameter names
      stop("Vector kinetics does not provide the correct parameters: PDbindingrate.")
    }
    insilicosystem$mosystem$PDRN_edg = dplyr::add_row(insilicosystem$mosystem$PDRN_edg, from = as.character(regID), to = as.integer(tarID),
                   TargetReaction = targetreaction, RegSign = regsign,
                   RegBy = regby,
                   PDregrate = myPDregrate)
  }


  if(targetreaction == "PTM"){
    if(is.null(kinetics)){ ## if no values are given for the kinetic parameters
      myPTMregrate = insilicosystem$sysargs[["PTMregrate_samplingfct"]](1)
    }else if(length(setdiff(c("PTMbindingrate"), names(kinetics))) == 0){ ## if values are given for the appropriate kinetic parameters
      myPTMregrate = kinetics$PTMregrate
    }else{ ## if the kinetic vector does not provide values with appropriate kinetic parameter names
      stop("Vector kinetics does not provide the correct parameters: PTMbindingrate.")
    }
    insilicosystem$mosystem$PTMRN_edg = dplyr::add_row(insilicosystem$mosystem$PTMRN_edg, from = as.character(regID), to = as.integer(tarID),
                   TargetReaction = targetreaction, RegSign = regsign,
                   RegBy = regby,
                   PTMregrate = myPTMregrate)
  }

  return(insilicosystem)

}

#' Removes an edge from the in silico system.
#'
#' Add an edge in the in silico system between specified genes.
#'
#' @param insilicosystem The in silico system (see \code{\link{createInSilicoSystem}}).
#' @param regID Integer. The ID of the regulator gene.
#' @param tarID Integer. The ID of the target gene.
#' @return The modified in silico system.
#' @export
removeEdg = function(insilicosystem, regID, tarID){

  ## Checking the input values ----
  if(class(insilicosystem) != "insilicosystem"){
    stop("Argument insilicosystem must be of class \"insilicosystem\".")
  }

  if(!(regID %in% insilicosystem$genes$id)){
    stop("Regulator id does not exist in the system.")
  }

  if(!(tarID %in% insilicosystem$genes$id)){
    stop("Target id does not exist in the system.")
  }

  ## The row to remove ----
  theedg = dplyr::filter(insilicosystem$edg, from == regID & to == tarID)

  if(nrow(theedg) == 0){ ## if the edge doesn't exists, no need to remove it!
    message("No edge exists from gene ", regID, " to gene ", tarID,".", sep = "")
    return(insilicosystem)
  }else if(nrow(theedg) > 1){ ## if more than one edge exists between these genes, there is a problem
    stop("More than one edge in the system meets the criteria! There must be a mistake somewhere.")
  }else{ ## If no problem, remove the edge from the data frames edg and XXRN_edg
    targetreaction = theedg[1, "TargetReaction"]
    insilicosystem$edg = dplyr::filter(insilicosystem$edg, from != regID | to != tarID)
    insilicosystem$mosystem[[paste0(targetreaction, "RN_edg")]] = dplyr::filter(insilicosystem$mosystem[[paste0(targetreaction, "RN_edg")]], from != regID | to != tarID)
  }

  return(insilicosystem)
}
