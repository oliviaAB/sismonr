#' Create variants for genes in the system.
#'
#' Create variants that segregate in the in silico population for each gene in the system. Returns a list that gives for each gene in the system the existing variants segregating in the in silico population.
#'
#' @param genes A data frame of genes in the system (created by the function \code{\link{createGenes}}).
#' @param indargs An object of class \code{\link{insilicoindividualargs}} (i.e. a list with parameters for in silico individuals generation).
#' @return A list of size G (number of genes in the system) where each element is a matrix corresponding to the QTL effect
#' coefficients (rows) of the different variants (columns) segregating in the in silico population for the corresponding gene. A variant is defined by a set of QTL effect coefficients ("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDregrate",
#' "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregrate") that correspond to the impact of genetic
#' mutations carried by the variant on the different kinetic parameters of the gene.
#' @export
createVariants = function(genes, indargs){

  G = nrow(genes)

  variants = vector("list", G)
  names(variants) = genes$id

  if(indargs$ngenevariants > 1){  ## if the user requests more that one variant per gene
    for(i in genes$id){
      potentialqtls = 1:(5 + 5*(genes[i, "coding"] == "PC")) ## protein-coding genes can have mutations affecting more parameters compared to noncoding genes
      nchanges = sample(potentialqtls, indargs$ngenevariants - 1, replace = T) ## Select the number of "mutations" from the original allele for each variant (minimum 1 otherwise we would have several copies of the original allele)
      ## the 1st variant is the original allele - no mutation
      qtlchanges = unlist(sapply(1:(indargs$ngenevariants-1), function(x){ length(indargs$qtlnames)*x + sample(potentialqtls, nchanges[x], replace = F)})) ## sample which qtl are affected by mutations for each variant, and convert it into matrix coordinates
      ## For each gene, the variants are stocked in the form of a matrix, rows being the different qtls effect coefficients and columns being the different gene variants. Element i,j corresponds to the value of the QTL effect coefficient i for variant j
      ## The matrix is first filled with 1 (no effect, allele identical to the "original allele"), and 0 for QTL effect coefficients irrelevant for the given gene (e.g. coefficient affecting the translation rate for noncoding genes)
      temp = matrix(rep(c(1, 0), indargs$ngenevariants * c(length(potentialqtls), length(indargs$qtlnames) - length(potentialqtls))), byrow = T, nrow = length(indargs$qtlnames), ncol = indargs$ngenevariants, dimnames = list(indargs$qtlnames, 1:indargs$ngenevariants))
      temp[qtlchanges] = indargs$qtleffect_samplingfct(sum(nchanges)) ## for each QTL effect coefficient that has been selected to be different from 1, sample a value using the qtleffect_samplingfct function
      variants[[i]] = temp
    }
  }else{ ## if the user only requests 1 variant per gene, the variants have no mutations (QTL effect coefficients = 1)
    for(i in genes$id){
      variants[[i]] = matrix(rep(1, length(indargs$qtlnames)), byrow = T, nrow = length(indargs$qtlnames), ncol = 1, dimnames = list(indargs$qtlnames, 1))
    }
  }
  return(variants)
}

#' Create an In Silico individual
#'
#' Creates a in silico individual to be simulated (object of class \code{insilicoindividual}).
#'
#' @param variantsList A named list giving the variants segregating in the population for each gene (e.g. created by \code{\link{createVariants}}). Each element corresponds to one gene in the system (name of the element = gene ID).
#' Each element is a matrix, in which each column represents a variant of the gene segregating in the population. The rows represent the QTL effect coefficients of each variant
#' (i.e. the impact of each mutation the variant carries).
#' @param variantsFreq A named list giving for each gene the allelic frequency of each segregating variant. Each element corresponds to one gene in the system (name of the element = gene ID).
#' Each element is a vector, of length equal to the number of variants of the gene segregating in the population, giving the allele frequency of each of the variants.
#' @param indargs An object of class \code{\link{insilicoindividualargs}} (i.e. a list with parameters for in silico individuals generation).
#' @param sameInit Boolean. Does the individual have identical initial abundance to the rest of the population?
#' @return An object of class \code{insilicoindividual}, that is a list composed of:
#' \itemize{
#' \item \code{QTLeffects}: a list of the variants carried by the individual. 1st level of the list: the different "GCN" (Gene Copy Number),
#' that is the different alleles of the genes (as defined by the ploidy of the individual. A diploid will have GCN1 and GCN2); 2nd level: the different QTL effect coefficients. The elements
#' in this 2nd-level list are vectors of QTL effect coefficients for the different genes (coefficient for gene \code{i} at the \code{i}-th position
#' in the vector).
#' \item \code{haplotype}: data-frame (rows = genes, columns = Gene copy number) of the gene variants carried by the individual for each gene copy number (homologue).
#' \item \code{InitVar}: a list of the multiplicative coefficients to compute the initial abundance of the different molecules
#' (to be applied to the population mean initial abundance for the corresponding molecule). 1st level of the list: the different "GCN" (Gene Copy Number),
#' that is the different alleles of the genes (as defined by the ploidy of the individual. A diploid will have GCN1 and GCN2); 2nd level of the list:
#' vectors of the coefficients for the proteins ("P") and RNAs ("R") of the genes (coefficient for gene \code{i} at the \code{i}-th position in the vectors).
#' }
#' @export
createIndividual = function(variantsList, variantsFreq, indargs, sameInit = F){

  ## Create the QTL effect coefficients list ----

  G = length(variantsList)
  QTLeffects = vector("list", indargs$ploidy)
  names(QTLeffects) = indargs$gcnList

  ## individualvariants: data frame where rows are genes and columns are "copy number" ids i.e. each columns represent a homolog chromosom
  ## Element i, j in the data frame corresponds to the variant of gene i present in the homolog chromosom j
#  individualvariants = as.data.frame(matrix(sample(1:indargs$ngenevariants, G*indargs$ploidy, replace = T), nrow = G, ncol = indargs$ploidy))
  temp = vector("list", length = G)
  for(g in 1:G){
  temp[[g]] = sample(1:ncol(variantsList[[paste(g)]]), indargs$ploidy, replace = T, prob = variantsFreq[[paste(g)]])
  }
  individualvariants = as.data.frame(do.call(rbind, temp))
  names(individualvariants) = indargs$gcnList

  ## Work for each gene copy (here in the sense each homolog chromosome)
  for(gcn in indargs$gcnList){
    QTLeffects[[gcn]] = vector("list", length(indargs$qtlnames))
    names(QTLeffects[[gcn]]) = indargs$qtlnames
    for(q in indargs$qtlnames){
      for(g in 1:G){
        QTLeffects[[gcn]][[q]][g] = variantsList[[g]][q, individualvariants[g, gcn]]
      }
    }
  }


  ## Create the Initial abundance variation coefficients list ----

  InitVar = vector("list", indargs$ploidy)
  names(InitVar) = indargs$gcnList

  if(sameInit){ ## if sameInit = T, we want the initial abundance of each molecule to be equal to the default value (no variation between individuals)
    for(gcn in indargs$gcnList){
      InitVar[[gcn]] = list("R" = rep(1.0, G),
                            "P" = rep(1.0, G))
    }
  }else{
    for(gcn in indargs$gcnList){
      InitVar[[gcn]] = list("R" = indargs$initvar_samplingfct(G),
                            "P" = indargs$initvar_samplingfct(G))
    }
  }

  value = list("QTLeffects" = QTLeffects, "haplotype" = individualvariants, "InitVar" = InitVar)
  attr(value, "class") = "insilicoindividual"

  return(value)
}

#' Creates a population of In Silico individuals
#'
#' Creates a population of in silico individuals to be simulated.
#'
#' @param nInd Integer. The number of in silico individuals to create.
#' @param insilicosystem An \code{insilicosystem} object. The in silico system based on which which individuals are created. See \code{\link{createInSilicoSystem}}.
#' @param genvariants A named list giving the variants segregating in the population for each gene. Each element corresponds to one gene in the system (name of the element = gene ID).
#' Each element is a matrix, in which each column represents a variant of the gene segregating in the population. The rows represent the QTL effect coefficients of each variant
#' (i.e. the impact of each mutation the variant carries). If none provided, will be automatically generated by the function \code{\link{createVariants}}.
#' @param genvariants.freq A named list giving for each gene the allelic frequency of each segregating variant. Each element corresponds to one gene in the system (name of the element = gene ID).
#' Each element is a vector, of length equal to the number of variants of the gene segregating in the population, giving the allele frequency of each of the variants.
#' If none provided, it is assumed that all variants of a given gene have the same allelic frequency.
#' @param sameInit Logical. Is the initial abundance of the different molecules the same for all individuals in the population? Default value is \code{FALSE}.
#' @param ... Other arguments to be passed to the function \code{\link{insilicoindividualargs}} (i.e. parameters for the generation of the in silico individuals).
#' @return An object of class \code{insilicopopulation}, that is a list composed of:
#' \itemize{
#' \item \code{GenesVariants} A list of variants segregating in the population for each genes (see \code{\link{createVariants}}).
#' \item \code{individualsList} A list of In Silico individuals (i.e. objects of class \code{insilicoindividual}, see \code{\link{createIndividual}}).
#' \item \code{indargs} An object of class \code{\code{insilicoindividualargs}}; the parameters used to create the in silico individuals.
#' }
#' @examples
#' ## Creates the in silico system (with 6 genes)
#' mysystem = createInSilicoSystem(G = 6)
#' ## Creates a first population with 3 diploid individuals,
#' ## With 2 variants of each gene segregating in the population
#' mypop1 = createInSilicoPopulation(nInd = 3, mysystem, ploidy = 2, ngenevariants = 2)
#'
#' ## Creates an other population with 10 tetraploid individuals
#' mypop2 = createInSilicoPopulation(nInd = 10, mysystem, ploidy = 4)
#'
#' ## Creates a population with a given list of gene variants
#' mysystem = createInSilicoSystem(G = 3, PC.p = 1)
#'
#' ## We will create only 1 variant of gene 1, 3 variants of gene 2 and
#' ## 2 variants of gene 3
#' nbvariants = c(1, 3, 2)
#'
#' qtlnames = c("qtlTCrate", "qtlRDrate",
#'              "qtlTCregbind", "qtlRDregrate",
#'              "qtlactivity", "qtlTLrate",
#'              "qtlPDrate", "qtlTLregbind",
#'              "qtlPDregrate", "qtlPTMregrate")
#'
#' genvariants = lapply(nbvariants, function(x){
#'   matrix(1, nrow = length(qtlnames), ncol = x,
#'          dimnames = list(qtlnames, 1:x))
#' })
#' names(genvariants) = mysystem$genes$id
#'
#' ## the 2nd variant of gene 2 has a mutation reducing its transcription rate by 3
#' genvariants$`2`["qtlTCrate", 2] = 0.3
#' ## and the 3rd variant has an increased translation rate
#' genvariants$`2`["qtlTLrate", 2] = 1.5
#'
#' ## The 2nd variant of gene 3 has a mutation decreasing the activity of
#' ## its active product
#' genvariants$`3`["qtlactivity", 2] = 0.7
#'
#' ## Allelic frequency of each variant
#' genvariants.freq = list('1' = c(1),
#'                         '2' = c(0.6, 0.3, 0.1),
#'                         '3' = c(0.9, 0.1))
#'
#' mypop = createInSilicoPopulation(10, mysystem,
#'                                  genvariants = genvariants,
#'                                  genvariants.freq = genvariants.freq)
#' @export
createInSilicoPopulation = function(nInd, insilicosystem, genvariants = NULL, genvariants.freq = NULL, sameInit = F, ...){

  if(class(insilicosystem) != "insilicosystem") stop("Argument insilicosystem must be of class \"insilicosystem\".")

  indargs = insilicoindividualargs(...)

  if(is.null(genvariants)){
    genvariants = createVariants(insilicosystem$genes, indargs)
  }else{
    if(typeof(genvariants) != "list") stop("Argument genvariants must be a list of length G (number of genes in the in silico system)")
    if(length(genvariants) != nrow(insilicosystem$genes)) stop("Argument genvariants must be a list of length G (number of genes in the in silico system)")
    if(!all(sapply(genvariants, function(x){identical(rownames(x), indargs$qtlnames)}))) stop("The row names of the elements of genvariants must be c(\"", paste(indargs$qtlnames, collapse = "\", \""), "\").")
  }

  if(is.null(genvariants.freq)){
    genvariants.freq = lapply(genvariants, function(x){rep(1/ncol(x), ncol(x))})
    names(genvariants.freq) = names(genvariants)
  }else{
    if(typeof(genvariants.freq) != "list") stop("Argument genvariants.freq must be a list of length G (number of genes in the in silico system)")
    if(length(genvariants.freq) != nrow(insilicosystem$genes)) stop("Argument genvariants.freq must be a list of length G (number of genes in the in silico system)")
    if(!all(sapply(names(genvariants.freq), function(x){length(genvariants.freq[[x]]) == ncol(genvariants[[x]])}))) stop("Each element of genvariants.freq must have the same length as the number of columns of the corresponding element of genvariants.")
  }

  indnames = sapply(1:nInd, function(x){paste0("Ind", x)})
  individualsList = vector("list", nInd)
  names(individualsList) = indnames

  for(i in indnames){
    individualsList[[i]] = createIndividual(genvariants, genvariants.freq, indargs, sameInit = sameInit)
  }

  value = list("GenesVariants" = genvariants, "individualsList" = individualsList, "indargs" = indargs)
  attr(value, "class") = "insilicopopulation"

  return(value)
}

