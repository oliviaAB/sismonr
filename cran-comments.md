## Resubmission

This is a resubmission (sismonr 2.1.0 has been accepted on CRAN) to fix some issues from the CRAN checks. In this version I have:

* replaced the .onAttach() function by .onLoad() function
* The.onLoad() function checks for the presence of the required Julia modules and issues a message if they are not available.
* Provided a separate function to install the required Julia modules; the function requires explicit user consent (through command line)
before installing the Julia modules.
* Fixed Julia syntax and corresponding tests to accommodate for Julia 1.6.5.

I also apologize to the CRAN maintainers for the last version of sismonr installing some Julia modules without explicit user consent, which is in violation of the CRAN policies. This has been corrected in the current version as described above. I am truly sorry for this oversight on my part.

## Test environments

- local Windows 10, R 4.1.2, Julia 1.6.4
- Ubuntu 18.04, R 4.1.2, Julia 1.6.5
- Winbuilder: x86_64-w64-mingw32 (64-bit) R 4.3.1
- Winbuilder: x86_64-w64-mingw32 (64-bit) R Under development (unstable) (2022-03-15 r81903 ucrt)
- Rhub:
  - windows-x86_64-devel: Windows Server 2022, R-devel, 64 bit
  - ubuntu-gcc-release: Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - fedora-clang-devel: Fedora Linux, R-devel, clang, gfortran
  - macos-highsierra-release-cran: macOS 10.13.6 High Sierra, R-release, CRAN's setup
  - debian-gcc-release: Debian Linux, R-release, GCC

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

```
* checking R code for possible problems ... [11s] NOTE
File 'sismonr/R/zzz.R':
  .onLoad calls:
    packageStartupMessage("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable.\n")
    packageStartupMessage("Error when checking the existing Julia version. Please check that Julia is correctly installed.")
    packageStartupMessage("The current version of Julia is < to v",     min_version, ". Please install Julia v", min_version, " or later.")
    packageStartupMessage("Required Julia modules not installed. Please run 'installJuliaModules()' to install them.")

See section 'Good practice' in '?.onAttach'.
```

After checking the good practice section of `?.onAttach`, I confirm that the startup messages are necessary to inform the user if the correct version of Julia is not available (as per the system requirements) or the required Julia modules not present. 


## Downstream dependencies

There are currently no downstream dependencies for this package
