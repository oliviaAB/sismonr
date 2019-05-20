## Resubmission

This is a resubmission. In this version I have:

* Added the reference to a conference paper published on our algorithm. A more complete paper describing the algorithm is in preparation.

* Only allowed functions to create text files if the user has enabled the writing and explicitely given an output directory.

* Created a series of tests for the exported functions, with a condition skipping the test if Julia v1.0 or later is not installed.

## Comments about previous feedback:

The sismonr package relies on opening a connection to a Julia process to execute Julia code from within the functions (I am using for that the CRAN package XRJulia). As such, any example will fail if Julia v1.0 (or later) is not installed on the server running the test. Also, as each  example is run independently, a connection to a Julia process must be started at the beginning of each example which takes a few seconds (and must be closed after the example). Consequently, the examples will take more than 5 seconds to run.

This is why I used \donttest for all examples in the functions documentation and instead created tests with a first condition checking if Julia is installed.

## Test environments
* local Windows 10, R 3.5.1, Julia 1.1.0
* local (VirtualBox) Ubuntu 16.04, R 3.6.0, Julia 1.1.0
* Ubuntu 14.04.5 (on travis-ci), R.6.0
* Winbuilder: x86_64-w64-mingw32 (64-bit) R 3.6.0 and R-devel
* Rhub: Windows Server 2008 R2 SP1, R-devel, 32/64 bit - Ubuntu Linux 16.04 LTS, R-release, GCC - Fedora Linux, R-devel, clang, gfortran

##R CMD check results
There were no ERRORs or WARNINGs.

On Winbuilder and Rhub, there was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE

Maintainer: 'Olivia Angelin-Bonnet <olivia.angelinbonnet@gmail.com>'
New submission

Possibly mis-spelled words in DESCRIPTION:
  Angelin (11:945)
  Omic (3:38)
  Silico (3:25)
  Vignes (11:979)
  Biggs (11:964)
  noncoding (11:184)
  ploidy (11:488)
  silico (11:75, 11:464, 11:754, 11:789)
  translational (11:315)

This is a new (re)submission to CRAN + all words are biological terms or last names and correctly spelled.
