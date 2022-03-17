## sismonr 2.3.0
* On loading, check of 1) presence of Julia, 2) correct version of Julia installed (v1.6.5 or higher) and 3) presence of the required Julia modules. If the required Julia modules are not available, a message will be issued, but the modules will not be installed. This can be done by the user through the call to a separate function.
* When creating a new Julia evaluator, if the required Julia modules are not available, calls a function to install them. The user has to provide explicit consent (through a prompt in the command line) before any Julia module is installed.
* Internal Julia syntax changed for consistency with Julia v1.6.5.

## sismonr 2.1.0
* Automatically updates BioSimulator to v0.9.3
* Accounts for changes of syntax in BioSimulator v0.9.3
* Names of possible stochastic simulation algorithms changed

## sismonr 2.0.1
* Forces the installation of BioSimulator v0.5.1

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
