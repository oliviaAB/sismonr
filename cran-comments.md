## Resubmission

This is a resubmission (sismonr 1.1.3 has been accepted on CRAN). In this version I have:

* Corrected the tests to resolve the Solaris installation failure.

## Comments about previous feedback:

This resubmission is to fix the Solaris installation failure showing on the CRAN check results. Please note that these additional errors did not appear when building the package on the Solaris installation with `devtools::check_rhub`. This is why I am re-submitting a patch.


## Test environments
- local Windows 10, R 3.5.1, Julia 1.1.0
- local (VirtualBox) Ubuntu 16.04, R 3.6.0, Julia 1.1.0
- Winbuilder: x86_64-w64-mingw32 (64-bit) R 3.6.0 and R-devel
- Rhub:
  - windows-x86_64-devel: Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  - ubuntu-gcc-release: Ubuntu Linux 16.04 LTS, R-release, GCC 
  - fedora-clang-devel: Fedora Linux, R-devel, clang, gfortran
  - solaris-x86-patched: Oracle Solaris 10, x86, 32 bit, R-patched (experimental)

##R CMD check results
There were no ERRORs or WARNINGs.

On Rhub, there was 1 NOTE:

* checking CRAN incoming feasibility ...NB: need Internet access to use CRAN incoming checks
 NOTE
Maintainer: ‘Olivia Angelin-Bonnet <olivia.angelinbonnet@gmail.com>’

Possibly mis-spelled words in DESCRIPTION:
  Angelin (11:945)
  Biggs (11:964)
  Omic (3:38)
  Silico (3:25)
  Vignes (11:979)
  noncoding (11:184)
  ploidy (11:488)
  silico (11:75, 11:464, 11:754, 11:789)
  translational (11:315)
  
All words are biological terms or last names and correctly spelled.
