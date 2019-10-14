#' Create variants for genes in the system.
#'
#' Create variants that segregate in the in silico population for each gene in the system.
#'
#' @param genes A data frame of genes in the system (created by the function \code{\link{createGenes}}).
#' @param indargs An object of class \code{\link{insilicoindividualargs}} (i.e. a list with parameters for in silico individuals generation).
#' @return A list of size G (number of genes in the system) where each element is a matrix corresponding to the QTL effect
#' coefficients (rows) of the different variants (columns) segregating in the in silico population for the corresponding gene. A variant is defined by a set of QTL effect coefficients ("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDregrate",
#' "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregrate") that correspond to the impact of genetic
#' mutations carried by the variant on the different kinetic parameters of the gene, as follow:
#' \itemize{
#' \item \code{qtlTCrate}: affects the basal transcription rate of the gene;
#' \item \code{qtlRDrate}: Affects the basal RNA decay rate of the gene;
#' \item \code{qtlTCregbind}: Affects the binding rate of the regulators of transcription on the gene's promoter (affects all transcription regulators targeting this gene);
#' \item \code{qtlRDregrate}: Affects the rate at which regulators of RNA decay encountering the RNAs of the gene trigger their degradation (affects all RNA decay regulators targeting this gene);
#' \item \code{qtlactivity}: Affects the activity of the active product of the gene. If the gene is encoding for a regulator of transcription or translation, this affects the binding rate of its active products (i.e. RNAs or proteins) to their
#' binding sites on their targets (affects the binding to all targets of the gene). If the gene encodes a regulator of RNA or protein decay or of protein post-translational modification, this affects the rate at which its active products
#' (i.e. RNAs or proteins) trigger the degradation/transformation of their targets (effect for all targets of the gene);
#' \item \code{qtlTLrate}: Affects the basal translation rate of the gene;
#' \item \code{qtlPDrate}: Affects the basal protein decay rate of the gene;
#' \item \code{qtlTLregbind}: Affects the binding rate of the regulators of translation on the gene's RNA binding sites (affects all translation regulators targeting this gene);
#' \item \code{qtlPDregrate}: Affects the rate at which regulators of protein decay encountering the proteins of the gene trigger their degradation (affects all protein decay regulators targeting this gene);
#' \item \code{qtlPTMregrate}: Affects the rate at which regulators of protein post-translational modification encountering the proteins of the gene trigger their modification (affects all protein post-translational modification regulators targeting this gene).
#'}
#' @examples
#' indargs = insilicoindividualargs()
#' genes = createGenes(insilicosystemargs(G = 5))
#' variants = createVariants(genes, indargs)
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

#' Creates an in silico individual.
#'
#' Creates a in silico individual to be simulated (object of class \code{insilicoindividual}).
#'
#' initialNoise: by default, the initial abundance of a molecule is equal to its steady state abundance in  the absence of any regulation
#' (e.g. for the RNA abundance of a gene, it is transcription rate / decay rate). If \code{initialNoise = TRUE}, instead the initial abundance of the
#' molecule will be sampled from a truncated Normal distribution of mean \code{SSabund} and SD \code{sqrt(SSabund)}, where \code{SSabund} is its
#' steady state abundance in the absence of any regulation, as specified above. The Normal distribution is truncated to only return positive values.
#'
#' @param insilicosystem An \code{insilicosystem} object. The in silico system based on which individuals are created. See \code{\link{createInSilicoSystem}}.
#' @param variantsList A named list giving the variants segregating in the population for each gene (e.g. created by \code{\link{createVariants}}). Each element corresponds to one gene in the system (name of the element = gene ID).
#' Each element is a matrix, in which each column represents a variant of the gene segregating in the population. The rows represent the QTL effect coefficients of each variant
#' (i.e. the impact of each mutation the variant carries).
#' @param variantsFreq A named list giving for each gene the allelic frequency of each segregating variant. Each element corresponds to one gene in the system (name of the element = gene ID).
#' Each element is a vector, of length equal to the number of variants of the gene segregating in the population, giving the allele frequency of each of the variants.
#' @param indargs An object of class \code{\link{insilicoindividualargs}} (i.e. a list with parameters for in silico individuals generation).
#' @param InitVar A list of the multiplicative coefficients to be applied to the initial abundance of the different molecules: elements "R" and "P" of the list giving the coefficients for the RNA
#' and protein form of the genes, respectively (coefficient for gene \code{i} at the \code{i}-th position in the vectors). If NULL, all coefficients set to 1.
#' @param initialNoise Logical. Is stochastic noise applied to the initial abundance of the different molecules? Default value is \code{TRUE} (see Details).
#'
#' @return An object of class \code{insilicoindividual}, that is a list composed of:
#' \itemize{
#' \item \code{QTLeffects}: a list of the variants carried by the individual. 1st level of the list: the different "GCN" (Gene Copy Number),
#' that is the different alleles of the genes (as defined by the ploidy of the individual: a diploid will have GCN1 and GCN2); 2nd level: the different QTL effect coefficients. The elements
#' in this 2nd-level list are vectors of QTL effect coefficients for the different genes (coefficient for gene \code{i} at the \code{i}-th position
#' in the vector).
#' \item \code{haplotype}: data-frame (rows = genes, columns = Gene copy number) giving the ID of the gene variant carried by the individual for each gene copy number (allele).
#' \item \code{InitAbundance}: A list of the initial abundance of the different molecules. 1st level of the list: the different "GCN" (Gene Copy Number),
#' that is the different alleles of the genes (as defined by the ploidy of the individual: a diploid will have GCN1 and GCN2); 2nd level of the list:
#' initial abundance of the protein ("P") and RNA ("R") form of the genes (coefficient for gene \code{i} at the \code{i}-th position in the vectors).
#' }
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 3, ploidy = 4)
#' indargs = insilicoindividualargs()
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
#' names(genvariants) = 1:length(nbvariants)
#'
#' ## the 2nd variant of gene 2 has a mutation reducing its transcription rate by 3
#' genvariants$`2`["qtlTCrate", 2] = 0.33
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
#' ## The third gene is not expressed at the beginning of the simulation
#' ## (its initial abundance is 0)
#' InitVar = list("R" = c(1, 1, 0), "P" = c(1, 1, 0))
#'
#' myind = createIndividual(mysystem, genvariants, genvariants.freq, indargs, InitVar = InitVar)
#' }
#' @export
createIndividual = function(insilicosystem, variantsList, variantsFreq, indargs, InitVar = NULL, initialNoise = TRUE){

  ## Create the QTL effect coefficients list ----

  G = length(variantsList)
  QTLeffects = vector("list", insilicosystem$sysargs$ploidy)
  names(QTLeffects) = insilicosystem$sysargs$gcnList

  ## individualvariants: data frame where rows are genes and columns are "copy number" ids i.e. each columns represent a homolog chromosom
  ## Element i, j in the data frame corresponds to the variant of gene i present in the homolog chromosom j
#  individualvariants = as.data.frame(matrix(sample(1:indargs$ngenevariants, G*insilicosystem$sysargs$ploidy, replace = T), nrow = G, ncol = insilicosystem$sysargs$ploidy))
  temp = vector("list", length = G)
  for(g in 1:G){
  temp[[g]] = sample(1:ncol(variantsList[[paste(g)]]), insilicosystem$sysargs$ploidy, replace = T, prob = variantsFreq[[paste(g)]])
  }
  individualvariants = as.data.frame(do.call(rbind, temp))
  names(individualvariants) = insilicosystem$sysargs$gcnList

  ## Work for each gene copy (here in the sense each homolog chromosome)
  for(gcn in insilicosystem$sysargs$gcnList){
    QTLeffects[[gcn]] = vector("list", length(indargs$qtlnames))
    names(QTLeffects[[gcn]]) = indargs$qtlnames
    for(q in indargs$qtlnames){
      for(g in 1:G){
        QTLeffects[[gcn]][[q]][g] = variantsList[[g]][q, individualvariants[g, gcn]]
      }
    }
  }


  ## Create the initial abundance list ----

  if(is.null(InitVar)){
    InitVar = list("R" = rep(1, G),
                   "P" = rep(1, G))
  }else{
      if(!all(names(InitVar) == c("R", "P"))){
        stop("Variable InitVar must be a list with elements \"R\" and \"P\".")
      }
      if(!all(sapply(InitVar, function(x){length(x) == G}))){
        stop("Each vector in InitVar must have a length equal to the number of genes in the system.")
      }
    }

  InitAbundance = vector("list", insilicosystem$sysargs$ploidy)
  names(InitAbundance) = insilicosystem$sysargs$gcnList

  for(gcn in insilicosystem$sysargs$gcnList){
    InitAbundance[[gcn]] = list("R" = rep(0, G),
                                "P" = rep(0, G))

    InitAbundance[[gcn]]$R = (insilicosystem$genes$TCrate * QTLeffects[[gcn]]$qtlTCrate) / (insilicosystem$genes$RDrate  * QTLeffects[[gcn]]$qtlRDrate)
    InitAbundance[[gcn]]$P = InitAbundance[[gcn]]$R * (insilicosystem$genes$TLrate * QTLeffects[[gcn]]$qtlTLrate) / (insilicosystem$genes$PDrate  * QTLeffects[[gcn]]$qtlPDrate)
    InitAbundance[[gcn]]$P[is.na(InitAbundance[[gcn]]$P)] = 0 ## true for non-coding genes, they don't have a protein form

    ## Apply the InitVar coefficients
    InitAbundance[[gcn]]$R = InitAbundance[[gcn]]$R * InitVar$R
    InitAbundance[[gcn]]$P = InitAbundance[[gcn]]$P * InitVar$P

    if(initialNoise){
      ## Applying initial noise: sample the actual initial abundance from a normal distribution
      ## with SD = sqrt(MEAN) so that there is not a lot of variation for molecules with small
      ## expected initial abundance - distribution truncated to 0 to only values (that can later be
      ## rounded to 0)
      InitAbundance[[gcn]]$R = sapply(InitAbundance[[gcn]]$R, function(x){
        if(x == 0) return(0)
        truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = x, sd = sqrt(x))
      })
      InitAbundance[[gcn]]$P = sapply(InitAbundance[[gcn]]$P, function(x){
        if(x == 0) return(0)
        truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = x, sd = sqrt(x))
      })
    }

    ## Round values because the initial abundances are discrete values
    InitAbundance[[gcn]]$R = round(InitAbundance[[gcn]]$R)
    InitAbundance[[gcn]]$P = round(InitAbundance[[gcn]]$P)
  }

  value = list("QTLeffects" = QTLeffects, "haplotype" = individualvariants, "InitAbundance" = InitAbundance)
  attr(value, "class") = "insilicoindividual"

  return(value)
}

#' Creates a population of in silico individuals.
#'
#' Creates a population of in silico individuals to be simulated.
#'
#' initialNoise: by default, the initial abundance of a molecule is equal to its steady state abundance in  the absence of any regulation
#' (e.g. for the RNA abundance of a gene, it is transcription rate / decay rate). If \code{initialNoise = TRUE}, instead the initial abundance of the
#' molecule will be sampled from a truncated Normal distribution of mean \code{SSabund} and SD \code{sqrt(SSabund)}, where \code{SSabund} is its
#' steady state abundance in the absence of any regulation, as specified above. The Normal distribution is truncated to only return positive values.
#'
#' @param nInd Integer. The number of in silico individuals to create.
#' @param insilicosystem An \code{insilicosystem} object. The in silico system based on which individuals are created. See \code{\link{createInSilicoSystem}}.
#' @param genvariants A named list giving the variants segregating in the population for each gene. Each element corresponds to one gene in the system (name of the element = gene ID).
#' Each element is a matrix, in which each column represents a variant of the gene segregating in the population. The rows represent the QTL effect coefficients of each variant
#' (i.e. the impact of each mutation the variant carries). If none provided, will be automatically generated by the function \code{\link{createVariants}}.
#' @param genvariants.freq A named list giving for each gene the allelic frequency of each segregating variant. Each element corresponds to one gene in the system (name of the element = gene ID).
#' Each element is a vector, of length equal to the number of variants of the gene segregating in the population, giving the allele frequency of each of the variants.
#' If none provided, it is assumed that all variants of a given gene have the same allelic frequency.
#' @param InitVar A list of the multiplicative coefficients to be applied to the initial abundance of the different molecules: elements "R" and "P" of the list giving the coefficients for the RNA
#' and protein form of the genes, respectively (coefficient for gene \code{i} at the \code{i}-th position in the vectors). If NULL, all coefficients set to 1.
#' @param initialNoise Logical. Is stochastic noise applied to the initial abundance of the different molecules? Default value is \code{TRUE} (see Details).
#' @param ... Other arguments to be passed to the function \code{\link{insilicoindividualargs}} (i.e. parameters for the generation of the in silico individuals).
#' @return An object of class \code{insilicopopulation}, that is a list composed of:
#' \itemize{
#' \item \code{GenesVariants} A list of variants segregating in the population for each genes (see \code{\link{createVariants}}).
#' \item \code{individualsList} A list of in silico individuals (i.e. objects of class \code{insilicoindividual}, see \code{\link{createIndividual}}).
#' \item \code{indargs} An object of class \code{\link{insilicoindividualargs}}; the parameters used to create the in silico individuals.
#' }
#' @examples
#' \donttest{
#' ## Creating a first population with 3 diploid individuals,
#' ## with 2 variants of each gene segregating in the population
#' mysystem = createInSilicoSystem(G = 6, ploidy = 2)
#' mypop1 = createInSilicoPopulation(nInd = 3, mysystem, ngenevariants = 2)
#'
#' ## Creating a population with 10 tetraploid individuals
#' mysystem = createInSilicoSystem(G = 6, ploidy = 4)
#' mypop2 = createInSilicoPopulation(nInd = 10, mysystem)
#'
#' ## Creating a population with a given list of gene variants
#' mysystem = createInSilicoSystem(G = 3, PC.p = 1, ploidy = 2)
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
#' genvariants$`2`["qtlTCrate", 2] = 0.33
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
#' ## The third gene is not expressed at the beginning of the simulation
#' ## (its initial abundance is 0)
#' InitVar = list("R" = c(1, 1, 0), "P" = c(1, 1, 0))
#'
#' mypop = createInSilicoPopulation(10, mysystem,
#'                                  genvariants = genvariants,
#'                                  genvariants.freq = genvariants.freq,
#'                                  InitVar = InitVar)
#' }
#' @export
createInSilicoPopulation = function(nInd, insilicosystem, genvariants = NULL, genvariants.freq = NULL, InitVar = NULL, initialNoise = TRUE, ...){

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
    individualsList[[i]] = createIndividual(insilicosystem, genvariants, genvariants.freq, indargs, InitVar = InitVar, initialNoise = initialNoise)
  }

  value = list("GenesVariants" = genvariants, "individualsList" = individualsList, "indargs" = indargs)
  attr(value, "class") = "insilicopopulation"

  return(value)
}

#' Plots the QTL effect coefficients of a population.
#'
#' Plots the QTL effect coefficients for all the genes of a in silico system of the in silico individuals.
#'
#' @param insilicopopulation The in silico population to be simulated (see \code{\link{createInSilicoPopulation}}).
#' @param insilicosystem The in silico system (object of class \code{insilicosystem}, see \code{\link{createInSilicoSystem}}).
#' @param scaleLims A vector of length 2 giving the lower and upper limits of the continuous scale (of the QTL effect coefficient
#' values) to be plotted. QTL effect coefficients with values outside these limits are plotted as NA (gray). If NULL (default), the limits are
#' automatically set to the min and max values in the dataset.
#' @param qtlEffectCoeffs A character vector of QTL effect coefficient names to plot. By default, all QTL effect coefficients are
#' represented.
#' @param inds A character vector giving the names of the individuals to plot. By default, all individuals in the population
#' are represented.
#' @param alleles A character vector giving the names of the alleles to plot. By default, all the alleles are represented.
#' @param genes A character or numeric vector of gene IDs to plot. By default, all the genes in the system are represented.
#' @param nGenesPerRow Integer. Number of genes to plot per row.
#' @param ... Any additional parameter to be passed to \code{\link[ggplot2]{theme}} for the plot.
#' @return A plot representing the value (colour) of each QTL effect coefficient (x-axis) of each allele (y-axis) of the different
#' individuals (rows) for each gene (column) in the system. For noncoding genes, some QTL effect coefficients are not relevant (the ones
#' related to protein or translation) and are represented in gray as NA.
#' @examples
#' \donttest{
#' mysystem = createInSilicoSystem(G = 10, ploidy = 2)
#' mypop = createInSilicoPopulation(10, mysystem)
#' plotMutations(mypop, mysystem)
#' ## Only plot the 1st allele of each genes for the genes 1 to 5 and the individuals 1 to 3
#' plotMutations(mypop, mysystem, alleles = c("GCN1"), genes = 1:5,
#'  inds = c("Ind1", "Ind2", "Ind3"))
#' }
#' @export
plotMutations = function(insilicopopulation, insilicosystem, scaleLims = NULL, qtlEffectCoeffs = insilicopopulation$indargs$qtlnames, inds = names(insilicopopulation$individualsList), alleles = insilicosystem$sysargs$gcnList, genes = 1:length(insilicopopulation$GenesVariants), nGenesPerRow = 10, ...){

  genes = as.numeric(genes)

  allqtlEffectCoeffs = insilicopopulation$indargs$qtlnames
  allGenes = 1:length(insilicopopulation$GenesVariants)

  dflist = lapply(inds, function(i){
    mutationsEffect = unname(unlist(insilicopopulation$individualsList[[i]]$QTLeffects[alleles]))

    res = data.frame(Ind = i,
                     Gene = rep(allGenes, times = length(allqtlEffectCoeffs)*length(alleles)),
                     Allele = rep(alleles, each = length(allqtlEffectCoeffs)*length(allGenes)),
                     QTLeffectCoefficient = rep(allqtlEffectCoeffs, each = length(allGenes), times = length(alleles)),
                     MutationsEffect = mutationsEffect,
                     stringsAsFactors = "F")
    res$Allele = factor(res$Allele, levels = sort(alleles, decreasing = T))
    res = dplyr::filter(res, !!sym("Gene") %in% genes & !!sym("QTLeffectCoefficient") %in% qtlEffectCoeffs)
    res$Gene = sapply(res$Gene, function(x){paste0("Gene ", x)})
    res$Gene = factor(res$Gene, levels = sapply(sort(genes), function(x){paste0("Gene ", x)}))
    return(res)
  })

  df = do.call(rbind, dflist)
  df$QTLeffectCoefficient = factor(df$QTLeffectCoefficient, levels = insilicopopulation$indargs$qtlnames)
  indsID = as.numeric(stringr::str_extract(inds, "\\d+"))
  df$Ind = factor(df$Ind, levels = sapply(sort(indsID), function(x){paste0("Ind", x)}))

  isNC = insilicosystem$genes[insilicosystem$genes$coding == "NC", "id"]
  isNC = sapply(isNC, function(x){paste0("Gene ", x)})

  df[df$Gene %in% isNC & df$QTLeffectCoefficient %in% insilicopopulation$indargs$qtlnames[6:10], "MutationsEffect"] = NA

  if(is.null(scaleLims)) scaleLims = c( min(df$MutationsEffect, na.rm = T), max(df$MutationsEffect, na.rm = T))


  # myplot = ggplot2::ggplot(df, aes_string(x = "QTLeffectCoefficient", y = "Allele", fill = "MutationsEffect")) +
  #   ggplot2::geom_tile() +
  #   ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
  #                                 midpoint = 1, na.value = "gray", limits = scaleLims, oob = scales::censor) +
  #   ggplot2::facet_grid(Ind~Gene) +
  #   ggplot2::theme_classic() +
  #   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90), legend.direction = "horizontal", legend.position = "bottom",
  #                  strip.background = ggplot2::element_rect(colour="black", fill="white"))

  nbrows = length(genes)%/%nGenesPerRow + (length(genes)%%nGenesPerRow !=0)
  geneNames = sapply(sort(genes), function(x){paste0("Gene ", x)})

  plots =  lapply(1:nbrows, function(i){
    genesrow = geneNames[(nGenesPerRow*(i-1) + 1):min(nGenesPerRow*i, length(genes))]
    dfGene = dplyr::filter(df, !!sym("Gene") %in% genesrow)
    myplot = ggplot2::ggplot(dfGene, aes_string(x = "QTLeffectCoefficient", y = "Allele", fill = "MutationsEffect")) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                    midpoint = 1, na.value = "gray", limits = scaleLims, oob = scales::censor, name = "QTL effect coefficient values\n(1 = no mutation)") +
      ggplot2::facet_grid(Ind~Gene) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90), legend.direction = "horizontal", legend.position = "bottom", ...) +
      ggplot2::xlab("QTL effect coefficients") + ggplot2::ylab("Alleles")
    return(myplot)
  })

  resplot = ggpubr::ggarrange(plotlist = plots, nrow = nbrows, ncol = 1, common.legend = T, heights = rep(1, length(plots)), legend = "bottom")

  return(resplot)
}
