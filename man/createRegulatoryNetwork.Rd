% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/in_silico_system.R
\name{createRegulatoryNetwork}
\alias{createRegulatoryNetwork}
\title{Creates an in silico regulatory network.}
\usage{
createRegulatoryNetwork(
  regsList,
  tarsList,
  reaction,
  sysargs,
  ev = getJuliaEvaluator()
)
}
\arguments{
\item{regsList}{A named list of length 2. Element "PC" (resp."NC") is a vector of gene IDs of the protein-coding (resp. noncoding) regulators
for the network.}

\item{tarsList}{A named list of length 2. Element "PC" (resp."NC") is a vector of gene IDs of the potential targets of the protein-coding (resp. noncoding)
regulators.}

\item{reaction}{String. The ID of the reaction targeted by the interactions ("TC", "TL", "RD", "PD" or "PTM").}

\item{sysargs}{An object of class \code{\link{insilicosystemargs}} (i.e. a list with parameters for in silico system generation).}

\item{ev}{A Julia evaluator (for the XRJulia package). If none provided select the current evaluator or create one if no evaluator exists.}
}
\value{
A list of two elements:
\itemize{
\item \code{edg}: a data-frame of edges of the network with the following variables:
\itemize{
\item \code{from}: gene ID of the regulator, as a character;
\item \code{to}: gene ID of the target, as an integer;
\item \code{TargetReaction}: the ID of the reaction (as given by \code{reaction});
\item \code{RegSign}: The sign of the reaction ("1" or "-1");
\item \code{RegBy}: Is the regulator a protein-coding gene ("PC"), a noncoding gene ("NC") or a complex ("C")?
}
\item \code{complexes}: a list of complexes composition (each element is named with the complex ID, the components are given as gene IDs).
\item \code{complexesTargetReaction}: a list defining which expression step the different regulatory complexes target (each element is named with the complex ID, the targeted reaction are given with a reaction ID, e.g. "TC" for transcription).
}
}
\description{
Creates an in silico regulatory network given a list of regulators and targets.
}
\examples{
\dontrun{
## We want to create a small transcription regulatory network
## In this example, genes 1 and 2 are protein-coding regulators (say transcription factors),
## gene 3 is a noncoding regulator (say an miRNA), and genes 4-6 are the genes to be regulated
## (all protein-coding, e.g. all encoding enzymes)
createRegulatoryNetwork(regsList = list("PC" = c(1:2), "NC" = c(3)),
     tarsList = list("PC" = c(4:6), "NC" = integer(0)), reaction = "TC",
     sysargs = insilicosystemargs(G = 6))
}
}
