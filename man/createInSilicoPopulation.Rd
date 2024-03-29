% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/in_silico_population.R
\name{createInSilicoPopulation}
\alias{createInSilicoPopulation}
\title{Creates a population of in silico individuals.}
\usage{
createInSilicoPopulation(
  nInd,
  insilicosystem,
  genvariants = NULL,
  genvariants.freq = NULL,
  InitVar = NULL,
  initialNoise = TRUE,
  ...
)
}
\arguments{
\item{nInd}{Integer. The number of in silico individuals to create.}

\item{insilicosystem}{An \code{insilicosystem} object. The in silico system based on which individuals are created. See \code{\link{createInSilicoSystem}}.}

\item{genvariants}{A named list giving the variants segregating in the population for each gene. Each element corresponds to one gene in the system (name of the element = gene ID).
Each element is a matrix, in which each column represents a variant of the gene segregating in the population. The rows represent the QTL effect coefficients of each variant
(i.e. the impact of each mutation the variant carries). If none provided, will be automatically generated by the function \code{\link{createVariants}}.}

\item{genvariants.freq}{A named list giving for each gene the allelic frequency of each segregating variant. Each element corresponds to one gene in the system (name of the element = gene ID).
Each element is a vector, of length equal to the number of variants of the gene segregating in the population, giving the allele frequency of each of the variants.
If none provided, it is assumed that all variants of a given gene have the same allelic frequency.}

\item{InitVar}{A list of the multiplicative coefficients to be applied to the initial abundance of the different molecules: elements "R" and "P" of the list giving the coefficients for the RNA
and protein form of the genes, respectively (coefficient for gene \code{i} at the \code{i}-th position in the vectors). If NULL, all coefficients set to 1.}

\item{initialNoise}{Logical. Is stochastic noise applied to the initial abundance of the different molecules? Default value is \code{TRUE} (see Details).}

\item{...}{Other arguments to be passed to the function \code{\link{insilicoindividualargs}} (i.e. parameters for the generation of the in silico individuals).}
}
\value{
An object of class \code{insilicopopulation}, that is a list composed of:
\itemize{
\item \code{GenesVariants} A list of variants segregating in the population for each genes (see \code{\link{createVariants}}).
\item \code{individualsList} A list of in silico individuals (i.e. objects of class \code{insilicoindividual}, see \code{\link{createIndividual}}).
\item \code{indargs} An object of class \code{\link{insilicoindividualargs}}; the parameters used to create the in silico individuals.
}
}
\description{
Creates a population of in silico individuals to be simulated.
}
\details{
initialNoise: by default, the initial abundance of a molecule is equal to its steady state abundance in  the absence of any regulation
(e.g. for the RNA abundance of a gene, it is transcription rate / decay rate). If \code{initialNoise = TRUE}, instead the initial abundance of the
molecule will be sampled from a truncated Normal distribution of mean \code{SSabund} and SD \code{sqrt(SSabund)}, where \code{SSabund} is its
steady state abundance in the absence of any regulation, as specified above. The Normal distribution is truncated to only return positive values.
}
\examples{
\dontrun{
## Creating a first population with 3 diploid individuals,
## with 2 variants of each gene segregating in the population
mysystem = createInSilicoSystem(G = 6, ploidy = 2)
mypop1 = createInSilicoPopulation(nInd = 3, mysystem, ngenevariants = 2)

## Creating a population with 10 tetraploid individuals
mysystem = createInSilicoSystem(G = 6, ploidy = 4)
mypop2 = createInSilicoPopulation(nInd = 10, mysystem)

## Creating a population with a given list of gene variants
mysystem = createInSilicoSystem(G = 3, PC.p = 1, ploidy = 2)

## We will create only 1 variant of gene 1, 3 variants of gene 2 and
## 2 variants of gene 3
nbvariants = c(1, 3, 2)

qtlnames = c("qtlTCrate", "qtlRDrate",
             "qtlTCregbind", "qtlRDregrate",
             "qtlactivity", "qtlTLrate",
             "qtlPDrate", "qtlTLregbind",
             "qtlPDregrate", "qtlPTMregrate")

genvariants = lapply(nbvariants, function(x){
  matrix(1, nrow = length(qtlnames), ncol = x,
         dimnames = list(qtlnames, 1:x))
})
names(genvariants) = mysystem$genes$id

## the 2nd variant of gene 2 has a mutation reducing its transcription rate by 3
genvariants$`2`["qtlTCrate", 2] = 0.33
## and the 3rd variant has an increased translation rate
genvariants$`2`["qtlTLrate", 2] = 1.5

## The 2nd variant of gene 3 has a mutation decreasing the activity of
## its active product
genvariants$`3`["qtlactivity", 2] = 0.7

## Allelic frequency of each variant
genvariants.freq = list('1' = c(1),
                        '2' = c(0.6, 0.3, 0.1),
                        '3' = c(0.9, 0.1))

## The third gene is not expressed at the beginning of the simulation
## (its initial abundance is 0)
InitVar = list("R" = c(1, 1, 0), "P" = c(1, 1, 0))

mypop = createInSilicoPopulation(10, mysystem,
                                 genvariants = genvariants,
                                 genvariants.freq = genvariants.freq,
                                 InitVar = InitVar)
}
}
