# sismonr - Simulation of In Silico Multi-Omic Networks R package

## IMPORTANT - Julia requirements
This current version of sismonr requires Julia version 0.6.4 and does not work with later versions of Julia. This will be fixed as soon as possible.

To install Julia v0.6.4, go to [https://julialang.org/downloads/oldreleases.html](https://julialang.org/downloads/oldreleases.html) and follow the instructions.
Please make sure to include the Julia executable in your PATH. For Linux users you can use the following command:
```
sudo ln -s path_to_julia_folder/bin/julia /usr/local/bin/julia
```
to create a symbolic link to julia inside the `/usr/local/bin` folder.

## sismonr

`sismonr` is a package for the simulation of gene expression profiles for *in silico* regulatory systems. Some inovative features of the model include:
- Simulation of protein-coding and **noncoding** genes;
- Simulation of transcriptional and **post-transcriptional** regulation;
- Simulation of the **ploidy** (number of gene copies) of the system.

The algorithm relies on the programming language `Julia` for computational speed. Therefore, you need to install `Julia` on your computer to use the package, but no skills in Julia programming are required!

3 steps are necessary to simulate an *in silico* system:
1. **Generate the "genome" and corresponding regulatory network** constituting your *in silico* system. This is done with the function `createInSilicoSystem(...)`. The user controls the number of genes in the system, the proportion of protein-coding vs noncoding genes, etc., but default values are provided.

2. **Generate the population of *in silico* "individuals"** for which you want to simulate the expression profiles of the abovementioned genes. This is done with the function `createInSilicoPopulation(...)`. The user can decide the ploidy of the individuals, that is the number of copies of the genes that each individual carries. The different individuals can be genetically identical or can carry different variant of the genes. 

3. **Simulate the expression profiles of the different genes** (that is the abundance of the corresponding RNAs and proteins if applicable) for each individual in the *in silico* population. This is done with the function `simulateInSilicoSystem`. Alternatively, for large systems, a parallel version of this function, `simulateParallelInSilicoSystem(...)`, is available (only for Linux).
