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
#' \item \code{PDrate}: Protein decay rate of the genes (0 for noncoding genes);
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

#' Generates an in silico regulatory network
#'
#' Generates an in silico regulatory network given a list of regulators and targets
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
