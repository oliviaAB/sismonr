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
#' @param sysargs An object of class \code{\link{insilicosystemargs}} (i.e. a list with parameters for in silico system generation).
#' @return A data frame of in silico genes. Attributes:
#' \itemize{
#' \item \code{id}: Integer, ID of the genes;
#' \item \code{coding}: coding status of the genes (either "PC" for protein-coding or "NC" for non-coding). Sampled according to the parameter \code{PC.p} in \code{sysargs};
#' \item \code{TargetReaction}: the biological function of the genes ("TC": transcription regulator, "TL": translation regulator, "RD": RNA decay
#' regulator, "PD": protein decay regulator, "PTM": post-translational modification regulator, "MR": metabolic enzyme). Sampled according to the parameters \code{PC.TC.p}, etc for protein-coding genes or \code{NC.TC.p}, etc for noncoding genes, in \code{sysargs};
#' \item \code{PTMform}: Does the gene have a PTM form? "0" or "1" (here all "0", PTM form will be assigned later);
#' \item \code{Active form}: what is the active form of the gene? "R" for noncoding genes, "P" for protein-coding genes,
#' "Pm" for protein-coding genes with a PTM form;
#' \item \code{TCrate}: transcription rate of the genes. Sampled according to the parameter \code{basal_transcription_rate_samplingfct} in \code{sysargs};
#' \item \code{TLrate}: translation rate of the genes. Sampled according to the parameter \code{basal_translation_rate_samplingfct} in \code{sysargs} (0 for noncoding genes);
#' \item \code{RDrate}: RNA decay rate of the genes. Sampled according to the parameter \code{basal_RNAlifetime_rate_samplingfct} in \code{sysargs};
#' \item \code{PDrate}: Protein decay rate of the genes. Sampled according to the parameter \code{basal_protlifetime_rate_samplingfct} in \code{sysargs} (0 for noncoding genes).
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
#' @param regsList A named list of length 2. Element "PC" (resp."NC") is a vector of gene IDs of the protein-coding (resp. noncoding) regulators
#' for the network.
#' @param tarsList A named list of length 2. Element "PC" (resp."NC") is a vector of gene IDs of the potential targets of the protein-coding (resp. noncoding)
#' regulators.
#' @param reaction String. The ID of the reaction targeted by the interactions ("TC", "TL", "RD", "PD" or "PTM").
#' @param sysargs An object of class \code{\link{insilicosystemargs}} (i.e. a list with parameters for in silico system generation).
#' @param ev A Julia evaluator (for the XRJulia package). If none provided select the current evaluator or create one if no evaluator exists.
#' @return A list of two elements:
#' \itemize{
#' \item \code{edg} a data-frame of edges of the network with the following variables:
#' \itemize{
#' \item \code{from} gene ID of the regulator, as a character;
#' \item \code{to} gene ID of the target, as an integer;
#' \item \code{TargetReaction} the ID of the reaction (as given by \code{reaction});
#' \item \code{RegSign} The sign of the reaction ("1" or "-1");
#' \item \code{RegBy} Is the regulator a protein-coding gene ("PC"), a noncoding gene ("NC") or a complex ("C")?
#' };
#' \item \code{complexes} a list of complexes composition (each element is named with the complex ID, the components are given as gene IDs).
#' \item \code{complexesTargetReaction} a list defining which expression step the different regulatory complexes target (each element is named with the complex ID, the targeted reaction are given with a reaction ID, e.g. "TC" for transcription).
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
    edg = data.frame("from" = character(), "to" = character(), "TargetReaction" = character(), "RegSign" = character(), "RegBy" = character(), stringsAsFactors = F)
  } else{
    edg = data.frame("from" = sapply(unlist(juliaedg$edg[,1]), toString), "to" = sapply(unlist(juliaedg$edg[,2]), toString), "TargetReaction" = rep(reaction, nrow(juliaedg$edg)), "RegSign" = rep("", nrow(juliaedg$edg)), "RegBy" = unlist(juliaedg$edg[,3]), stringsAsFactors = F)

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
  complexes = lapply(complexes, as.character)
  complexesTargetReaction = as.list(rep(reaction, length(complexes)))
  names(complexesTargetReaction) = names(complexes)
  return(list("edg" = edg, "complexes" = complexes, "complexesTargetReaction" = complexesTargetReaction))
}


#' Creates an in silico multi omic network
#'
#' Creates an in silico multi omic network from a data-frame of genes in the system.
#'
#' @param genes A data-frame of genes in the system (see \code{\link{createGenes}}).
#' @param sysargs An object of class \code{\link{insilicosystemargs}} (i.e. a list with parameters for in silico system generation).
#' @param ev A Julia evaluator. If none provided select the current evaluator or create one if no evaluator exists.
#' @return An in silico multi omic network, that is a list of:
#' \itemize{
#' \item \code{genes}: the modified data-frame of genes;
#' \item \code{edg}: A data-frame of edges in the multi-omic network (1 row = 1 edge). It contains the following parameters:
#' \itemize{
#' \item \code{from}: gene ID of the origin of the edge (character).
#' \item \code{to}: gene ID of the destination of the edge (character).
#' \item \code{TargetReaction}: Type of regulation (ID of the controlled reaction: "TC", "TL", "RD", "PD" or "PTM").
#' \item \code{RegSign}: Sign of the reaction ("1" for positive regulation, "-1" for negative regulation).
#' \item \code{RegBy}: Type of the regulator: "PC" for protein-coding regulation, "NC" for noncoding regulator,
#' "C" for regulatory complex.
#' }
#' \item \code{mosystem}: A list of the different regulation networks in the system and associated information. Elements are \code{TCRN_edg}, \code{TLRN_edg}, \code{RDRN_edg}, \code{PDRN_edg} and \code{PTMRN_edg}: data-frames of edges for the different
#' regulatory networks, with in addition to the usual fields in the edg data frame contains columns for kinetic parameters of the
#' regulation.
#' \item \code{complexes}: a list of regulatory complexes composition. The names of the elements are the IDs of the complexes, and the
#' values are vectors of gene IDs constituting each regulatory complex.
#' \item \code{complexeskinetics}: a list of regulatory complexes kinetic parameters.
#' \item \code{complexesTargetReaction}: a list defining which expression step is targeted by each regulatory complex.
#' }
#' @export
createMultiOmicNetwork = function(genes, sysargs, ev = getJuliaEvaluator()){

  complexes = list()
  complexesTargetReaction = list()

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
  complexesTargetReaction = c(complexesTargetReaction, TCRN[["complexesTargetReaction"]])

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
  complexesTargetReaction = c(complexesTargetReaction, TLRN[["complexesTargetReaction"]])

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
  complexesTargetReaction = c(complexesTargetReaction, RDRN[["complexesTargetReaction"]])

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
  complexesTargetReaction = c(complexesTargetReaction, PDRN[["complexesTargetReaction"]])

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
  complexesTargetReaction = c(complexesTargetReaction, PTMRN[["complexesTargetReaction"]])

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

  return(list("mosystem" = res, "genes" = genes, "edg" = edg, "complexes" = complexes, "complexeskinetics" = complexeskinetics, "complexesTargetReaction" = complexesTargetReaction))

}

#' Creates an empty in silico multi-omic network
#'
#' Creates an empty in silico multi-omic network (i.e. no regulatory interactions) given a data-frame of genes (cf
#' \code{\link{createMultiOmicNetwork}}).
#'
#' @param genes A data-frame of genes in the system (see \code{\link{createGenes}}).
#' @return An in silico multi-omic network, that is a list of:
#' \itemize{
#' \item \code{genes}: the modified data-frame of genes;
#' \item \code{edg}: A data-frame of edges in the multi-omic network (1 row = 1 edge, here empty dataframe). It contains the following parameters:
#' \itemize{
#' \item \code{from}: gene ID of the origin of the edge (character).
#' \item \code{to}: gene ID of the destination of the edge (character).
#' \item \code{TargetReaction}: Type of regulation (ID of the controlled reaction: "TC", "TL", "RD", "PD" or "PTM").
#' \item \code{RegSign}: Sign of the reaction ("1" for positive regulation, "-1" for negative regulation).
#' \item \code{RegBy}: Type of the regulator: "PC" for protein-coding regulation, "NC" for noncoding regulator,
#' "C" for regulatory complex.
#' }
#' \item \code{mosystem}: A list of the different regulation networks in the system and associated information. Elements of the list are
#' \code{TCRN_edg}, \code{TLRN_edg}, \code{RDRN_edg}, \code{PDRN_edg} and \code{PTMRN_edg}: data-frames of edges for the different
#' regulatory networks, with in addition to the usual fields in the \code{edg} data frame contain columns for kinetic parameters of the
#' regulation. All empty.
#' \item \code{complexes}: a list of regulatory complexes composition. The names of the elements are the IDs of the complexes, and the
#' values are vectors of gene IDs constituting each regulatory complex. Empty list.
#' \item \code{complexeskinetics}: a list of regulatory complexes kinetic parameters. Empty list.
#' \item \code{complexesTargetReaction}: a list defining which expression step is targeted by each regulatory complex.
#' }
#' @export
createEmptyMultiOmicNetwork = function(genes){

  edg = data.frame("from" = character(), "to" = character(), "TargetReaction" = character(), "RegSign" = character(), "RegBy" = character(), stringsAsFactors = F)
  res = list("TCRN_edg" = data.frame(edg,  "TCbindingrate" = numeric(), "TCunbindingrate" = numeric(), "TCfoldchange" = numeric(), stringsAsFactors = F),
             "TLRN_edg" = data.frame(edg,  "TLbindingrate" = numeric(), "TLunbindingrate" = numeric(), "TLfoldchange" = numeric(), stringsAsFactors = F),
             "RDRN_edg" = data.frame(edg, "RDregrate" = numeric(), stringsAsFactors = F),
             "PDRN_edg" = data.frame(edg, "PDregrate" = numeric(), stringsAsFactors = F),
             "PTMRN_edg" = data.frame(edg, "PTMregrate" = numeric(), stringsAsFactors = F))

  genes$ActiveForm = sapply(1:nrow(genes), function(x){paste0(genes$ActiveForm[x], genes$id[x])})

  return(list("mosystem" = res, "genes" = genes, "edg" = edg, "complexes" = list(), "complexeskinetics" = list(), "complexesTargetReaction" = list()))
}

#' Creates an in silico system.
#'
#' Creates an in silico system, i.e. the genes and the regulatory networks defining the system.
#'
#' @param empty Logical. Does the regulatory network is empty (= no regulation)? Default value is \code{FALSE}.
#' @param ev A Julia evaluator (for the XRJulia package). If none provided select the current evaluator or create one if no evaluator exists.
#' @param ... Other arguments to be passed to the function \code{\link{insilicosystemargs}} (i.e. parameters for the generation of the in silico system).
#' @return An object of class \code{insilicosystem}, that is a list composed of:
#' \itemize{
#' \item \code{genes}: a data-frame of genes (see \code{\link{createGenes}});
#' \item \code{edg}: a data-frame of edges in the regulatory network (see \code{\link{createMultiOmicNetwork}});
#' \item \code{mosystem}: a list defining the multi-omic regulatory network (see \code{\link{createMultiOmicNetwork}});
#' \item \code{sysargs}: An object of class \code{insilicosystemargs}; the parameters used to create the system.
#' }
#' @examples
#' \donttest{
#' ## Creates an in silico system composed of 20 genes
#' mysystem1 = createInSilicoSystem(G = 20)
#' mysystem1$edg ## see all regulations in the system
#' mysystem1$mosystem$TCRN_edg ## see only regulations targeting transcription
#'
#' ## Creates an in silico systerm composed of 10 genes, all protein-coding
#' mysystem2 = createInSilicoSystem(G = 10, PC.p = 1)
#' mysystem2$genes
#'
#' ## Creates an in silico systerm composed of 5 genes,
#' ## all noncoding and all regulators of transcription
#' mysystem3 = createInSilicoSystem(G = 5, PC.p = 0, NC.TC.p = 1)
#' mysystem3$edg
#' mysystem3$TCRN_edg
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

#' Adds a gene in the in silico system.
#'
#' Adds a gene in the in silico system with specified parameters if provided, or with parameters sampled according to the system parameters.
#'
#' @param insilicosystem The in silico system (see \code{\link{createInSilicoSystem}}).
#' @param coding String. The coding status of the gene (either "PC" for protein-coding or "NC" for noncoding). If none provided, randomly chosen according to the
#' parameter \code{PC.p} provided in sysargs (see \code{\link{insilicosystemargs}}).
#' @param TargetReaction String. The biological function of the gene, i.e. the gene expression step targeted by the active product
#' of the gene. If none provided, randomly chosen according to the parameters \code{PC.TC.p}, etc or \code{NC.TC.p}, etc (depending on the coding status of the gene)
#' provided in sysargs (see \code{\link{insilicosystemargs}}).
#' @param TCrate Numeric. The transcription rate of the gene. If none provided, randomly chosen according to the
#' parameter \code{basal_transcription_rate_samplingfct} provided in sysargs (see \code{\link{insilicosystemargs}}).
#' @param TLrate Numeric. The translation rate of the gene. If none provided, randomly chosen according to the
#' parameter \code{basal_translation_rate_samplingfct} provided in sysargs (see \code{\link{insilicosystemargs}}).
#' @param RDrate Numeric. The RNA decay rate of the gene. If none provided, randomly chosen according to the
#' parameter \code{basal_RNAlifetime_samplingfct} provided in sysargs (see \code{\link{insilicosystemargs}}).
#' @param PDrate Numeric. The protein decay rate of the gene. If none provided, randomly chosen according to the
#' parameter \code{basal_protlifetime_samplingfct} provided in sysargs (see \code{\link{insilicosystemargs}}).
#' @return The modified in silico system.
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 5)
#' mysystem$genes
#' mysystem2 = addGene(mysystem, "PC", "TC", TCrate = 0.0001, TLrate = 0.001)
#' mysystem2$genes
#'
#' mysystem3 = addGene(mysystem2)
#' mysystem3$genes
#' }
#' @export
addGene = function(insilicosystem, coding = NULL, TargetReaction = NULL, TCrate = NULL, TLrate = NULL, RDrate = NULL, PDrate = NULL){

  ## Checking the input values ----
  if(class(insilicosystem) != "insilicosystem"){
    stop("Argument insilicosystem must be of class \"insilicosystem\".")
  }

  id = max(insilicosystem$genes$id) + 1


  if(!is.null(coding)){
    if(!(coding %in% c("PC", "NC"))){
      stop("\"coding\" argument must either be \"PC\" or \"NC\".")
    }
  }else{
    coding = sample(c("PC", "NC"), 1, prob = c(insilicosystem$sysargs[["PC.p"]], insilicosystem$sysargs[["NC.p"]]), replace = T)
  }

  ## If not provided, sampling kinetic parameters ----
  if(is.null(TCrate)){
    TCrate = insilicosystem$sysargs[["basal_transcription_rate_samplingfct"]](1)
  }

  if(is.null(RDrate)){
    RDrate = 1/insilicosystem$sysargs[["basal_RNAlifetime_samplingfct"]](1)
  }

  if(coding == "NC"){
    TLrate = 0.0
    PDrate = 0.0
    ActiveForm = paste0("R", id)

    ## Biological function of the gene
    if(is.null(TargetReaction)){
      TargetReaction = sample(c("TC", "TL", "RD", "PD", "PTM"), 1, prob = c(insilicosystem$sysargs[["NC.TC.p"]], insilicosystem$sysargs[["NC.TL.p"]],
                                                                            insilicosystem$sysargs[["NC.RD.p"]], insilicosystem$sysargs[["NC.PD.p"]],
                                                                            insilicosystem$sysargs[["NC.PTM.p"]]), replace = T)
    }else{
      if(!(TargetReaction %in% c("TC", "TL", "RD", "PD", "PTM"))){
        stop("Argument \"TargetReaction\" must be either \"TC\", \"TL\", \"RD\", \"PD\" or \"PTM\".")
      }
    }
  }else{ ## if coding == "PC"
    if(is.null(TLrate)){
      TLrate = insilicosystem$sysargs[["basal_translation_rate_samplingfct"]](1)
    }

    if(is.null(PDrate)){
      PDrate = 1/insilicosystem$sysargs[["basal_protlifetime_samplingfct"]](1)
    }
    ActiveForm = paste0("P", id)

    ## Biological function of the gene
    if(is.null(TargetReaction)){
      TargetReaction = sample(c("TC", "TL", "RD", "PD", "PTM", "MR"), 1, prob = c(insilicosystem$sysargs[["PC.TC.p"]], insilicosystem$sysargs[["PC.TL.p"]],
                                                                                  insilicosystem$sysargs[["PC.RD.p"]], insilicosystem$sysargs[["PC.PD.p"]],
                                                                                  insilicosystem$sysargs[["PC.PTM.p"]], insilicosystem$sysargs[["PC.MR.p"]]), replace = T)
    }else{
      if(!(TargetReaction %in% c("TC", "TL", "RD", "PD", "PTM", "MR"))){
        stop("Argument \"TargetReaction\" must be either \"TC\", \"TL\", \"RD\", \"PD\", \"PTM\" or \"MR\".")
      }
    }

  }

  PTMform = "0"

  insilicosystem$genes = dplyr::add_row(insilicosystem$genes, id = as.integer(id), coding = coding, TargetReaction = TargetReaction,
                                                              PTMform = PTMform, ActiveForm = ActiveForm, TCrate = TCrate,
                                                              TLrate = TLrate, RDrate = RDrate, PDrate = PDrate)

  insilicosystem$sysargs$G = insilicosystem$sysargs$G + 1

  return(insilicosystem)
}

# #' Removes a gene from the in silico system.
# #'
# #' Removes a gene from the in silico system. Any edge involving this gene is removed from the system,
# #' and the composition of the complexes comprising this gene are adjusted.
# #'
# #' @param insilicosystem The in silico system (see \code{\link{createInSilicoSystem}}).
# #' @param id Integer. The id of the gene to remove.
# #' @return The modified in silico system.
# #' @export
# removeGene = function(insilicosystem, id){
#   id2rem = as.integer(id) ## to avoid confusion when using dplyr::filter
#   id2rem.c = as.character(id)
#
#   ## Checking the input values ----
#   if(class(insilicosystem) != "insilicosystem"){
#     stop("Argument insilicosystem must be of class \"insilicosystem\".")
#   }
#
#   if(!(id2rem %in% insilicosystem$genes$id)){
#       stop("Gene ", id2rem, " is not in the system.")
#   }
#
#   ## Remove the gene from the gene list
#   targetreaction = paste0(insilicosystem$genes[insilicosystem$genes$id == id2rem, "TargetReaction"], "RN_edg")
#   insilicosystem$genes = dplyr::filter(insilicosystem$genes, id != id2rem)
#
#   ## Remove any edge involving the gene in the general edges list
#   insilicosystem$edg = dplyr::filter(insilicosystem$edg, from != id2rem.c & to != id2rem.c)
#
#   ## Remove any edge involving the gene in the specific multi-omic edges list
#   for(rn in names(insilicosystem$mosystem)){
#     insilicosystem$mosystem[[rn]] = dplyr::filter(insilicosystem$mosystem[[rn]], from != id2rem.c & to != id2rem.c)
#   }
#
#   ## If the gene is involved in a complex:
#   for(comp in names(insilicosystem$complexes)){
#     if(id2rem.c %in% insilicosystem$complexes[[comp]]){
#       newcompo = setdiff(insilicosystem$complexes[[comp]], id2rem.c )
#       if(length(newcompo) > 1){ ## if there are still at least 2 components in the complex, simply remove the gene from the complex component list
#         insilicosystem$complexes[[comp]] = newcompo
#       }else{
#         insilicosystem$complexes[[comp]] = NULL
#         insilicosystem$complexeskinetics[[comp]] = NULL
#
#         ## Replace all instances of the complex by the remaining component
#         newcoding = insilicosystem$genes[insilicosystem$genes$id == newcompo, "coding"]
#         rows2change = which(insilicosystem$edg$from == comp)
#         insilicosystem$edg[rows2change, "from"] = newcompo
#         insilicosystem$edg[rows2change, "RegBy"] = newcoding
#         rows2change = which(insilicosystem$mosystem[[targetreaction]]$from == comp)
#         insilicosystem$mosystem[[targetreaction]][rows2change, "from"] = newcompo
#         insilicosystem$mosystem[[targetreaction]][rows2change, "RegBy"] = newcoding
#       }
#     }
#   }
#   return(insilicosystem)
# }

#' Adds a regulatory complex in the in silico system.
#'
#' Adds a regulatory complex in the in silico system with specified parameters (if provided), or with parameters sampled according to the system parameters.
#'
#' @param insilicosystem The in silico system (see \code{\link{createInSilicoSystem}}).
#' @param compo An vector of integers, corresponding to the ID of the genes composing the complex. All genes composing the complex must
#' have the same biological function (i.e. same \code{TargetReaction} parameter).
#' @param formationrate The formation rate of the complex. If none provided, randomly chosen according to the
#' parameter \code{complexesformationrate_samplingfct} provided in \code{sysargs} (see \code{\link{insilicosystemargs}}).
#' @param dissociationrate The dissociation rate of the complex. If none provided, randomly chosen according to the
#' parameter \code{complexesdissociationrate_samplingfct} provided in \code{sysargs} (see \code{\link{insilicosystemargs}}).
#' @return Returns the modified in silico system.
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 10, PC.p = 1, PC.TC.p = 1)
#' mysystem$complexes ## no complex in the system
#' mysystem2 = addComplex(mysystem, c(1, 2, 3))
#' mysystem2$complexes
#' }
#' @export
addComplex = function(insilicosystem, compo, formationrate = NULL, dissociationrate = NULL){

  compo = as.character(compo) ## make sure the id of the components are strings

  compoG = compo[!stringr::str_detect(compo, "^C")] ## components of the complex that are gene products
  compoC = compo[stringr::str_detect(compo, "^C")] ## components of the complex that are regulatory complexes

  ## Checking the input values ----
  if(class(insilicosystem) != "insilicosystem"){
    stop("Argument insilicosystem must be of class \"insilicosystem\".")
  }

  if(!all(as.integer(compoG) %in% insilicosystem$genes$id)){
    stop("The components of the complex do not exist in the system.")
  }
  if(!all(compoC %in% names(insilicosystem$complexes))){
    stop("The components of the complex do not exist in the system.")
  }

  # if(length(compo) != insilicosystem$sysargs$regcomplexes.size){
  #  stop("Wrong number of components. The complex must be of size ", insilicosystem$sysargs$regcomplexes.size,".")
  # }

  targetreactions = c(insilicosystem$genes[which(insilicosystem$genes$id %in% as.integer(compoG)), "TargetReaction"],
                      unname(unlist(insilicosystem$complexesTargetReaction[compoC])))
  if(!all(targetreactions == targetreactions[1])){
    stop("The different components do not all have the same biological function.")
  }


  exnames = names(insilicosystem$complexes)[stringr::str_detect(names(insilicosystem$complexes), paste0("C", targetreactions[1]))] ## names of the complexes in the system targeting the same reaction/expression step
  if(length(exnames)>0){
    exnum = as.numeric(stringr::str_extract(exnames, "(\\d)+"))
  }else{
    exnum = 0
  }
  name = paste0("C", targetreactions[1], max(exnum)+1)


  if(is.null(formationrate)){
    formationrate = insilicosystem$sysargs[["complexesformationrate_samplingfct"]](1)
  }
  if(is.null(dissociationrate)){
    dissociationrate = insilicosystem$sysargs[["complexesdissociationrate_samplingfct"]](1)
  }

  insilicosystem$complexes[[name]] = compo

  insilicosystem$complexeskinetics[[name]] = list("formationrate" = formationrate, "dissociationrate" = dissociationrate)

  insilicosystem$complexesTargetReaction[[name]] = targetreactions[1]

  return(insilicosystem)
}

#' Removes a regulatory complex from the in silico system.
#'
#' Removes a regulatory complex from the in silico system. Any edge involving this complex is removed from the system.
#'
#' @param insilicosystem The in silico system (see \code{\link{createInSilicoSystem}}).
#' @param name String. The name of the regulatory complex to remove.
#' @return The modified in silico system.
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 10, PC.p = 1, PC.TC.p = 1, regcomplexes.p = 0.8)
#' mysystem$complexes
#' mysystem$edg
#' mysystem2 = removeComplex(mysystem, "CTC1")
#' mysystem2$complexes
#' mysystem2$edg
#' }
#' @export
removeComplex = function(insilicosystem, name){

  ## Checking the input values ----
  if(class(insilicosystem) != "insilicosystem"){
    stop("Argument insilicosystem must be of class \"insilicosystem\".")
  }

  if(!(name %in% names(insilicosystem$complexes))){
    stop("Complex ", name, " does not exist in the system.")
  }

  insilicosystem$complexes[[name]] = NULL
  insilicosystem$complexeskinetics[[name]] = NULL
  insilicosystem$complexesTargetReaction[[name]] = NULL

  ## Remove any edge involving the complex in the general edges list
  insilicosystem$edg = dplyr::filter_(insilicosystem$edg, "from" != name)

  ## Remove any edge involving the gene in the specific multi-omic edges list
  for(rn in names(insilicosystem$mosystem)){
    insilicosystem$mosystem[[rn]] = dplyr::filter_(insilicosystem$mosystem[[rn]], "from" != name)
  }

  return(insilicosystem)

}

#' Adds an edge in the in silico system's regulatory network.
#'
#' Adds an edge in the in silico system's regulatory network between specified genes.
#'
#' @param insilicosystem The in silico system (see \code{\link{createInSilicoSystem}}).
#' @param regID Integer. The ID of the regulator gene.
#' @param tarID Integer. The ID of the target gene.
#' @param regsign The sign of the regulation: either "1" (positive regulation) or "-1" (negative regulation). If none provided,
#' will be randomly chosen according to the parameter \code{TC.pos.p}, \code{TL.pos.p} or \code{PTM.pos.p} (depending on the type
#' of regulation - regulation of RNA or protein decay can only be negative) provided in sysargs (see \code{\link{insilicosystemargs}}).
#' @param kinetics Optional: named vector of kinetics parameters of the reaction. If none provided,
#' will be randomly chosen according to the parameters \code{[name of the param]_samplingfct} provided in sysargs (see \code{\link{insilicosystemargs}}). The parameters
#' to provide depend on the type of regulation (i.e. parameter \code{TargetReaction} of the regulator):
#' \itemize{
#' \item TargetReaction = "TC". The parameters to specify are "TCbindingrate", "TCunbindingrate" and "TCfoldchange";
#' \item TargetReaction = "TL". The parameters to specify are "TLbindingrate", "TLunbindingrate" and "TLfoldchange";
#' \item TargetReaction = "RD". The parameter to specify is "RDregrate";
#' \item TargetReaction = "PD". The parameter to specify is "PDregrate";
#' \item TargetReaction = "PTM". The parameter to specify is "PTMregrate".
#' }
#' @return The modified in silico system.
#' @examples
#' \donttest{
#' ## creates a system with no regulation
#' mysystem = createInSilicoSystem(G = 10, PC.p = 1, PC.TC.p = 1, empty = TRUE)
#' mysystem$edg
#' mysystem2 = addEdge(mysystem, 1, 2, regsign = "1",
#'  kinetics = c("TCbindingrate"= 0.01, "TCunbindingrate" = 0.1, "TCfoldchange" = 10))
#' ## check all existing interactions in the system (no kinetic parameters)
#' mysystem2$edg
#' ## check the interactions targeting transcription, with kinetic parameters
#' mysystem2$mosystem$TCRN_edg
#'
#' ## creates a system with no regulation
#' mysystem = createInSilicoSystem(G = 5, PC.p = 1, PC.PD.p = 1, empty = TRUE)
#' mysystem$edg
#' mysystem2 = addEdge(mysystem, 1, 2)
#' ## check all existing interactions in the system (no kinetic parameters)
#' mysystem2$edg
#' ## check the interactions targeting protein decay, with kinetic parameters
#' mysystem2$mosystem$PDRN_edg
#' }
#' @export
addEdge = function(insilicosystem, regID, tarID, regsign = NULL, kinetics = list()){
  regID = as.character(regID)
  tarID = as.character(tarID)
  ## Checking the input values ----
  if(class(insilicosystem) != "insilicosystem"){
    stop("Argument insilicosystem must be of class \"insilicosystem\".")
  }

  if(!(regID %in% as.character(insilicosystem$genes$id))){ ## if the regulator is not a gene product
    if(!(regID %in% names(insilicosystem$complexes))){
      stop("Regulator ", regID, "does not exist in the system.")
    }else{ ## if the regulator is a regulatory complex
      targetreaction = insilicosystem$complexesTargetReaction[[regID]]
      regby = "C"
    }
  }else{ ## if the regulator is a gene product
    targetreaction = insilicosystem$genes[insilicosystem$genes$id == as.integer(regID), "TargetReaction"]
    regby = dplyr::filter_(insilicosystem$genes, "id" == as.integer(regID))[1,"coding"]
  }

  if(grepl("^C", tarID)){ ## if the tarID given is a complex
    stop("Complexes cannot be regulated. Please provide a gene ID as target.")
  }
  if(!(as.integer(tarID) %in% insilicosystem$genes$id)){
    stop("Target ", tarID, " does not exist in the system.")
  }

  ## Checking if an edge already exists between the 2 genes ----
  if(nrow(dplyr::filter_(insilicosystem$edg, "from" == regID & "to" == tarID)) != 0){
    stop(paste0("An edge already exists from gene ", regID, " to gene ", tarID, "."))
  }

  abbr = c("TC" = "transcription (TC)", "TL" = "translation (TL)", "RD" = "RNA decay (RD)", "PD" = "protein decay (PD)", "PTM" = "post-translational modification (PTM)")

  codtar = insilicosystem$genes[insilicosystem$genes$id == as.integer(tarID), "coding"]
  if(codtar == "NC" & targetreaction %in% c("TL", "PD", "PTM")){
    stop("Target gene ", tarID, " is a noncoding gene. Cannot be regulated at the level of ", abbr[targetreaction], ".")
  }

  if(!is.null(regsign)){
    if(!(regsign %in% c("1", "-1"))){
      stop("Interation sign must be either \"1\" or \"-1\".")
    }else if(targetreaction %in% c("RD", "PD")){
      regsign = "-1"
    }
  }else{
    regsign = sample(c("1","-1"), 1, prob = c(insilicosystem$sysargs[[paste(targetreaction, "pos.p", sep = ".")]],
                                              1 - insilicosystem$sysargs[[paste(targetreaction, "pos.p", sep = ".")]]), replace = T)
  }

  ## which kinetic parameters must be created for the edge
  params = list("TC" = c("TCbindingrate" , "TCunbindingrate", "TCfoldchange"),
                "TL" = c("TLbindingrate" , "TLunbindingrate", "TLfoldchange"),
                "RD" = c("RDregrate"),
                "PD" = c("PDregrate"),
                "PTM" = c("PTMregrate"))

  ## Adding the interaction in the edg data-frame ----
  insilicosystem$edg = dplyr::add_row(insilicosystem$edg, from = regID, to = tarID,
                                      TargetReaction = targetreaction, RegSign = regsign, RegBy = regby)

  ## Retrieving/sampling the kinetic parameters
  kin2sample = setdiff(params[[targetreaction]], names(kinetics)) ## kinetic parameters not provided by the user
  kingiven = intersect(params[[targetreaction]], names(kinetics)) ## kinetic parameters provided by the user

  ## sampling the values for the parameters not provided
  kinsampled = sapply(kin2sample, function(i){
    insilicosystem$sysargs[[paste0(i,"_samplingfct")]](1)
  })

  kin = unlist(c(kinsampled, kinetics[kingiven]))
  names(kin) = c(kin2sample, kingiven)

  if(targetreaction %in% c("TC", "TL") & regsign == "-1"){
    kin[paste0(targetreaction, "foldchange")] = 0
  }

  ## add the edg and kinetic parameters to the adequate mutli-omic edge list
  edg2add = data.frame("from" = regID, "to" = tarID, "TargetReaction" = targetreaction,
              "RegSign" = regsign, "RegBy" = regby, t(kin), stringsAsFactors = F)

  insilicosystem$mosystem[[paste0(targetreaction, "RN_edg")]] = bind_rows(insilicosystem$mosystem[[paste0(targetreaction, "RN_edg")]],
                                                                          edg2add)

  return(insilicosystem)
}

#' Removes an edge from the in silico system.
#'
#' Removes an edge in the in silico system between specified genes.
#'
#' @param insilicosystem The in silico system (see \code{\link{createInSilicoSystem}}).
#' @param regID Integer or character. The ID of the regulator gene or the name of the regulatory complex.
#' @param tarID Integer or character. The ID of the target gene.
#' @return The modified in silico system.
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 10)
#' mysystem$edg
#' ## we'll remove the first edge
#' regToRemove = mysystem$edg$from[1]
#' tarToRemove = mysystem$edg$to[1]
#' mysystem2 = removeEdge(mysystem, regToRemove, tarToRemove)
#' mysystem2$edg
#' }
#' @export
removeEdge = function(insilicosystem, regID, tarID){

  regID = as.character(regID)
  tarID = as.character(tarID)

  ## Checking the input values ----
  if(class(insilicosystem) != "insilicosystem"){
    stop("Argument insilicosystem must be of class \"insilicosystem\".")
  }

  ## The row to remove ----
  theedg = dplyr::filter_(insilicosystem$edg, "from" == regID & "to" == tarID)

  if(nrow(theedg) == 0){ ## if the edge doesn't exists, no need to remove it!
    message("No edge exists from gene ", regID, " to gene ", tarID,".", sep = "")
    return(insilicosystem)
  }else if(nrow(theedg) > 1){ ## if more than one edge exists between these genes, there is a problem
    stop("More than one edge in the system meets the criterion! There must be a mistake somewhere.")
  }else{ ## If no problem, remove the edge from the data frames edg and XXRN_edg
    targetreaction = theedg[1, "TargetReaction"]
    insilicosystem$edg = dplyr::filter_(insilicosystem$edg, "from" != regID | "to" != tarID)
    insilicosystem$mosystem[[paste0(targetreaction, "RN_edg")]] = dplyr::filter_(insilicosystem$mosystem[[paste0(targetreaction, "RN_edg")]], "from" != regID | "to" != tarID)
  }

  return(insilicosystem)
}


# ## depends on package VisNetwork
# plotRegulatoryNetwork = function(insilicosystem){
#   mysystem = createInSilicoSystem(G = 15, PC.p = 0.5)
#   colsCS = c("PC" = "#4da6ff",  "NC" = "#ff3333", "Complexes" = "#808080")
#   nodes = mysystem$genes %>%
#     select(id, coding) %>%
#     dyplr::rename(group = coding) %>%
#     mutate(color = colsCS[group],
#            label = paste0("Gene ", id),
#            shadow = T)
#
#   for(i in names(mysystem$complexes)){
#     nodes = rbind(nodes, data.frame(id = i, group = "Complex", color = colsCS["Complexes"], label = paste0("Complex ", i), shadow = T, row.names = i))
#   }
#
#   colsGF = c("TC" = "#b30000", "TL" = "#0033cc", "RD" = "#ff6600", "PD" = "#33ccff", "PTM" = "#8000ff", "Complex formation" = "#808080")
#   edges = mysystem$edg %>%
#     select(-RegBy) %>%
#     mutate(dashes = RegSign == "-1",
#            color = colsGF[TargetReaction],
#            arrows = "to")
#
#   legendnodes = data.frame(label = c("Protein-coding genes", "Noncoding genes", "Complexes"),
#                            color = colsCS)
#
#   legendedges = data.frame(label = c("TC regulation", "TL regulation", "RD regulation", "PD regulation", "PTM regulation", "Complex formation", "Activation", "Repression"),
#                            color = c(colsGF, "black", "black"), dashes = rep(c(F, T), c(6, 2)))
#
#   for(i in names(mysystem$complexes)){
#     components = mysystem$complexes[[i]]
#     edges = rbind(edges, data.frame(from = components,
#                                     to = i,
#                                     TargetReaction = "Complex formation",
#                                     RegSign = "1",
#                                     dashes = F,
#                                     color = colsGF["Complex formation"],
#                                     arrows = "to"))
#   }
#
#   visNetwork::visNetwork(nodes, edges, main = paste(nrow(mysystem$genes), "genes,", nrow(mysystem$edg), "interactions", sep = " ")) %>%
#     visNetwork::visOptions(highlightNearest = T, , selectedBy = "group") %>%
#     visNetwork::visLegend(addEdges = legendedges, addNodes = legendnodes, useGroups = F)
# }
