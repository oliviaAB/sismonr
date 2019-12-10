# sismonr 2.0.0

## New features
* Added a function to transform abundance data into RNAseq-like data (read counts)

## Stochastic model
* Added the degradation of molecules even if bound to another molecule or in a complex
* Ploidy is now defined for a system, not for a population
* Initial abundances are now available when generating a population, calculated accounting for the ploidy of the system
* Corrected a mistake in the addEdge function for sampling of kinetic parameters

## Visualisation features
* Added the option to focus on a subset of genes/complexes for plotSimulation and plotHeatMap functions
* Corrected y-axis title for plotSimulation function

# sismonr 1.1.4
* tests fixed for Solaris installation.

# sismonr 1.1.3
* .onAttach function fixed for Solaris installation.

# sismonr 1.1.2
* Added tests (with condition for skipping tests if Julia v1.0 or later is not installed)
* Re-submission to CRAN


# sismonr 1.1.1
* Added a check of the existing version of Julia when attaching the package.
* Submission to CRAN
