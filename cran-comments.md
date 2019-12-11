## Resubmission

This is a resubmission (sismonr 1.1.4 has been accepted on CRAN). In this version I have:

* Added a function to transform abundance data into RNAseq-like data (read counts)
* Improved the stochastic model generated
* Added some new options for the visualisation functions

## Test environments
- local Windows 10, R 3.5.1, Julia 1.1.0
- local (VirtualBox) Ubuntu 18.04, R 3.6.1, Julia 1.2.0
- Winbuilder: x86_64-w64-mingw32 (64-bit) R 3.6.1
- Rhub:
  - windows-x86_64-devel: Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  - ubuntu-gcc-release: Ubuntu Linux 16.04 LTS, R-release, GCC 
  - fedora-clang-devel: Fedora Linux, R-devel, clang, gfortran
  - solaris-x86-patched: Oracle Solaris 10, x86, 32 bit, R-patched (experimental)

##R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Possibly mis-spelled words in DESCRIPTION:
    GitHub (11:1255)
    sismonr (11:1039, 11:1154)

On Rhub - Fedora Linux, there was 1 NOTE:

* checking CRAN incoming feasibility ...NB: need Internet access to use CRAN incoming checks
   NOTE
   
  Possibly mis-spelled words in DESCRIPTION:
    Angelin (11:945)
    Biggs (11:964)
    Omic (3:38)
    Silico (3:25)
    Vignes (11:979)
    noncoding (11:184)
    ploidy (11:488)
    silico (11:75, 11:464, 11:754, 11:789)
    sismonr (11:1039, 11:1154)
    translational (11:315)
    
All words are biological terms or last names and correctly spelled.
