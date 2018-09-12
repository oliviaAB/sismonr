#' Create variants for genes in the system.
#'
#' Create variants that segregate in the in silico population for each gene in the system.
#'
#' Returns a list that gives for each gene in the system the existing variants segregating in the in silico population.
#' A variant is defined by a set of QTL effect coefficients ("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDbindreg",
#' "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregbind") that correspond to the impact of genetic
#' mutations on the variant.
#'
#' @param genes A data frame of genes in the system (created by the function \code{\link{createGenes}}).
#' @param indargs An object of class \code{insilicoindividualargs} (i.e. a list with parameters for in silico individuals generation).
#' @return A list of size G (number of genes in the system) where each element is a matrix corresponding to the QTL effect
#' coefficients of the different variants for the corresponding gene.
#' @export
createVariants = function(genes, indargs){

  G = nrow(genes)

  variants = vector("list", G)
  names(variants) = genes$id

  if(indargs$ngenevariants > 1){  ## if the user requests more that one variant per gene
    for(i in genes$id){
      potentialqtls = 1:(5 + 4*(genes[i, "coding"] == "PC")) ## protein-coding genes can have mutations affecting more parameters compared to noncoding genes
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
#' @param variantsList List of segregating variants (created by \code{\link{createVariants}}). Each element is a matrix of the existing variants
#' for a given gene in the system.
#' @param indargs An object of class \code{insilicoindividualargs} (i.e. a list with parameters for in silico individuals generation).
#' @param sameInit Boolean. Does the individual have identical initial abundance to the rest of the population?
#' @return An object of class \code{insilicoindividual}.
#' @export
createIndividual = function(variantsList, indargs, sameInit = F){

  ## Create the QTL effect coefficients list ----

  G = length(variantsList)
  QTLeffects = vector("list", indargs$ploidy)
  names(QTLeffects) = indargs$gcnList

  ## individualvariants: data frame where rows are genes and columns are "copy number" ids i.e. each columns represent a homolog chromosom
  ## Element i, j in the data frame corresponds to the variant of gene i present in the homolog chromosom j
  individualvariants = as.data.frame(matrix(sample(1:indargs$ngenevariants, G*indargs$ploidy, replace = T), nrow = G, ncol = indargs$ploidy))
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
#' @param insilicosystem An \code{insilicosystem} object. The in silico system from which individuals are created.
#' @param indargs An object of class \code{insilicoindividualargs} (i.e. a list with parameters for in silico individuals generation).
#' @param sameInit Logical. Do the individuals in the population have the same initial abundance for the different molecules? Default value is \code{FALSE}.
#' @export
createPopulation = function(nInd, insilicosystem, indargs, sameInit = F){

  genvariants = createVariants(insilicosystem$genes, indargs)
  indnames = sapply(1:nInd, function(x){paste0("Ind", x)})
  individualsList = vector("list", nInd)
  names(individualsList) = indnames

  for(i in indnames){
    individualsList[[i]] = createIndividual(genvariants, indargs, sameInit = sameInit)
  }

  value = list("GenesVariants" = genvariants, "individualsList" = individualsList, "indargs" = indargs)
  attr(value, "class") = "insilicopopulation"

  return(value)
}

