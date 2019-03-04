% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stochastic_simulation.R
\name{plotHeatMap}
\alias{plotHeatMap}
\title{Plots the result of a simulation as a heatmap}
\usage{
plotHeatMap(simdf, inds = unique(simdf$Ind),
  trials = unique(simdf$trial), timeMin = min(simdf$time),
  timeMax = max(simdf$time), mergeAllele = T, mergePTM = T,
  mergeComplexes = F, yLogScale = T, nIndPerRow = 3,
  VirPalOption = "plasma", ...)
}
\arguments{
\item{simdf}{The dataframe with the result of the simulation (see \code{\link{simulateInSilicoSystem}})}

\item{inds}{A vector of in silico individual names for which to plot the expression profiles}

\item{trials}{A vector of trials ID (=number) to use for the plot (see details).}

\item{timeMin}{Numeric. The minimum simulation time to plot. Default value set to the minimum time in the simulation.}

\item{timeMax}{Numeric. The maximum simulation time to plot. Default value set to the maximum time in the simulation.}

\item{mergeAllele}{Are the gene products originating from different alleles merged? Default TRUE. Also see \code{\link{mergeAlleleAbundance}}}

\item{mergePTM}{Are the modified and non-modified versions of the proteins merged? Default TRUE. Also see \code{\link{mergePTMAbundance}}}

\item{mergeComplexes}{Are the free and in complex gene products merged? Default FALSE Also see \code{\link{mergeComplexesAbundance}}}

\item{yLogScale}{Is the y-axis of the plot in log10-scale? If so, the abundance of each species at each time-point is increased by 1 to avoid zero values. Default TRUE.}

\item{nIndPerRow}{Positive integer, the number of individuals to plot per row. Default 3.}

\item{VirPalOption}{String, palette name option to be passed to \code{\link[ggplot2]{scale_fill_viridis_c}}; can be one of "magma", "inferno", "plasma", "viridis" or "cividis". Default value is "plasma".}

\item{...}{Any additional parameter to be passed to \code{\link[ggplot2]{theme}} for the plot of each individual.}
}
\value{
A plot from \code{\link[ggpubr]{ggarrange}}.
}
\description{
Automatically plots the result of a simulation for the selected in silico individuals as a heatmap.
}
\details{
If more than one trial is to be plotted, the mean abundance of each molecule over the different trials is plotted with a solid line,
and the min and max abundances represented as coloured areas around the mean.
}
\examples{
\donttest{
mysystem = createInSilicoSystem(G = 5, regcomplexes = "none")
mypop = createInSilicoPopulation(15, mysystem, ploidy = 2)
sim = simulateParallelInSilicoSystem(mysystem, mypop, 100, ntrials = 5)
plotHeatMap(sim$Simulation, c("Ind1", "Ind2", "Ind3", "Ind4"),
 axis.title = element_text(color = "red"))
}
}