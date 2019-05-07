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

New submission
Possibly mis-spelled words in DESCRIPTION:
  translational (11:315)
  noncoding (11:184)
  Maintainer: 'Olivia Angelin-Bonnet <olivia.angelinbonnet@gmail.com>'
  silico (11:75, 11:464, 11:754, 11:789)
  Omic (3:38)
  Silico (3:25)
  ploidy (11:488)

This is a new submission to CRAN + all words are biological terms and correctly spelled.
