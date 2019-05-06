## Test environments
* local Windows 10, R 3.5.1, Julia 1.1.0
* local (VirtualBox) Ubuntu 16.04, R 3.6.0, Julia 1.1.0
* Ubuntu 14.04.5 (on travis-ci), R.6.0

##R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Comments
We used Winbuilder to perform R CMD check. This returned 4 WARNINGs and 4 NOTEs, all of them due to a Julia error:

```
RROR: LoadError: ArgumentError: Module Pkg not found in current path.
Run `Pkg.add("Pkg")` to install the Pkg package.
Stacktrace:
 [1] _require(::Symbol) at .\loading.jl:428
 [2] require(::Symbol) at .\loading.jl:398
 [3] include_from_node1(::String) at .\loading.jl:569
 [4] include(::String) at .\sysimg.jl:14
 [5] process_options(::Base.JLOptions) at .\client.jl:305
 [6] _start() at .\client.jl:371
while loading d:\RCompile\CRANguest\R-devel\lib\sismonr\julia\testModules.jl, in expression starting on line 1
```
This is probably due to the presence of Julia 0.6.4 instead of 1.0, and was not observed when running R CMD check on a local Windows 10 machine with Julia 1.1.0.

We used `check_rhub()`, which returned 4 WARNINGs and 5 NOTEs, all:

```
INFO: Could not find files for the given pattern(s).
```
Which is probably due to julia not being available.
