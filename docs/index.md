*An R package for generating and simulating in silico biological systems.*

* TOC
{:toc}

# Introduction

*Explain the concepts of in silico system and individual*

The code used throughout this documentation as well as the created variables saved as `.RData` objects are available [here](https://github.com/oliviaAB/sismonr/tree/master/docs/example).

## Abbreviations

## Quickstart

*Write a quick tuto?*

# Creating an *in silico* system

The first step is to generate an *in silico* system. An *in silico* system is composed of a set of genes, and a gene regulatory network or GRN describing the different regulatory interactions occuring between the genes. The system is created by a call to the function `createInSilicoSystem`. The user can control different aspects of the system with the arguments passed to the function. For example:
```r
mysystem = createInSilicoSystem(G = 10, PC.p = 0.7)
```
generates an *in silico* system with 10 genes, and during the generation process each gene has a probability of 0.7 to be a protein-coding gene (as opposed to noncoding gene). The system returned by the function is an object of class `insilicosystem`, i.e. a list whose different attributes are presented below.

## The list of genes

The different genes constituting the system are stored in a data-frame, and can be accessed with:
```r
> mysystem$genes

   id coding TargetReaction PTMform ActiveForm       TCrate     TLrate       RDrate       PDrate
1   1     NC             TL       0         R1 8.887839e-03 0.00000000 1.736133e-04 0.000000e+00
2   2     PC             TL       0         P2 6.711473e-05 0.13226328 8.512016e-05 1.164180e-06
3   3     PC             TC       0         P3 3.193124e-04 0.05770822 1.553645e-04 7.507520e-06
4   4     PC             TL       0         P4 3.748385e-04 0.50706887 6.013606e-04 6.906099e-06
5   5     PC             TL       0         P5 4.415198e-04 0.05302814 4.046843e-04 7.347904e-06
6   6     PC             TC       0         P6 4.691504e-04 0.06439765 6.992371e-04 2.385040e-06
7   7     PC             TL       1        Pm7 7.079025e-04 0.02980232 1.160875e-03 1.338248e-05
8   8     PC            PTM       0         P8 1.611630e-03 0.16976215 1.779581e-04 7.994026e-06
9   9     NC             TC       0         R9 1.325744e-03 0.00000000 1.185721e-03 0.000000e+00
10 10     PC             TC       0        P10 2.972445e-03 0.03415530 3.167154e-04 1.828237e-05
```

Each gene is labeled with an ID (column `id`) and possesses the following parameters:

- `coding`: gives the **coding status** of the gene, i.e. describes if the gene is protein-coding (`coding = "PC"`) or noncoding (`coding = "NC"`);
- `TargetReaction`: describes the biological function of the gene, i.e. the type of regulation that it performs on its targets. This parameter can take one of the following values:
  - `"TC"`: regulator of transcription
  - `"TL"`: regulator of translation
  - `"RD"`: regulator of RNA decay
  - `"PD"`: regulator of protein decay
  - `"PTM"`: regulator of protein post-translational modification
  - `"MR"`: metabolic enzyme (only for protein-coding genes, indicates that the gene cannot regulate the expression of another gene);
- `PTMform`: does the protein of the gene has a (post-translationally) modified form? Will always be `"0"` for noncoding genes;
- `ActiveForm`: what is the active form of the gene? If the gene is noncoding, then `ActiveForm` is set to `"R[gene ID]"` (e.g. `R1` for gene 1), and for a protein-coding gene `ActiveForm = "P[gene ID]"`. If the protein of a protein-coding gene is targeted for post-translational modification (`PTMform = "1"`) then we assume that only the modified form of the protein is active, and thus `ActiveForm = Pm[gene ID]`;
- `TCrate`, `TLrate`, `RDrate` and `PDrate`: give the transcription (in RNA.s$$^{-1}$$), translation (in protein.RNA$$^{-1}$$.s$$^{-1}$$), RNA decay and protein decay (in s$$^{-1}$$) rates of the genes, respectively. `TLrate` and `PDrate` are set to 0 for noncoding genes.

## The Gene Regulatory Network

The regulatory network describing the regulatory interactions between the genes are stored in a data-frame, and accessible with:
```r
> mysystem$edg

   from to TargetReaction RegSign RegBy
1     9  2             TC       1    NC
2     9  3             TC      -1    NC
3    10  4             TC      -1    PC
4     9  5             TC       1    NC
5     3  6             TC      -1    PC
6     9  7             TC      -1    NC
7     9  8             TC       1    NC
8     6  9             TC      -1    PC
9     1 10             TL      -1    NC
10    1  2             TL       1    NC
11    5  3             TL       1    PC
12    4  4             TL       1    PC
13    2  5             TL       1    PC
14    2  7             TL      -1    PC
15 CTL1  8             TL      -1     C
16    8  7            PTM       1    PC
```
Each edge in the network is characterised by the following parameters:

- `from`: ID of the regulator gene. Note that for the edge at row 15, the regulator is not a gene but a regulatory complex (identifiable by its ID starting with `'C'`);
- `RegBy`: type of regulator (`"PC"` for a protein-coding regulator, `"NC"` for a noncoding regulator and `"C"` for a regulatory complex);
- `to`: ID of the target gene;
- `TargetReaction`: type of the regulation, i.e. which expression step of the target is controlled. For example an edge for which `TargetReaction = "TC"` represents a regulation of transcription, etc (see the [Abbreviations](#abbreviations) section);
- `RegSign`: sign of the regulation (`"1"` for an activation and `"-1"` for a repression). Edges corresponding to the regulation of RNA or protein decay always have `RegSign = "1"`, meaning that the regulator increases the decay rate of the target.

The `edge` dataframe shows the global GRN, with all the different types of regulations. The element `mosystem` of the `insilicosystem` object contains the same edges but grouped by type of regulation:
```r
> names(mysystem$mosystem)

[1] "TCRN_edg"  "TLRN_edg"  "RDRN_edg"  "PDRN_edg"  "PTMRN_edg"

> mysystem$mosystem$TCRN_edg

  from to TargetReaction RegSign RegBy TCbindingrate TCunbindingrate TCfoldchange
1    9  2             TC       1    NC   0.006863351     0.009372557     5.173425
2    9  3             TC      -1    NC   0.008071258     0.004499718     0.000000
3   10  4             TC      -1    PC   0.004963256     0.001410222     0.000000
4    9  5             TC       1    NC   0.001369649     0.003318900     7.341233
5    3  6             TC      -1    PC   0.003139301     0.003830472     0.000000
6    9  7             TC      -1    NC   0.001967162     0.002181207     0.000000
7    9  8             TC       1    NC   0.008571345     0.001387365    21.336259
8    6  9             TC      -1    PC   0.005159756     0.009424001     0.000000
```
Here you can only see the edges of the GRN corresponding to the regulation of transcription. Each of the sub-GRNs stored in the `mosystem` element contains additional information about the kinetic parameters associated with each edge. For example, we can see that each edge corresponding to a regulatory interaction targeting the transcription is assigned a binding rate (`TCbindingrate`) and an unbinding rate (`TCunbindingrate`) of the regulator to and from the binding site on the target gene's promoter. The parameter `TCfoldchange` corresponds to the coefficient by which is multiplied the basal transcription rate of the target when the regulator is bound to its binding site (notice that for edges for which `RegSign = "-1"`, i.e. corresponding to a repression, `TCfoldchange = 0`). The kinetic parameters associated with each edge depend on the type of regulation:

```r
> lapply(mysystem$mosystem, function(x){colnames(x)[-(1:5)]})

$`TCRN_edg`
[1] "TCbindingrate"   "TCunbindingrate" "TCfoldchange"   

$TLRN_edg
[1] "TLbindingrate"   "TLunbindingrate" "TLfoldchange"   

$RDRN_edg
[1] "RDregrate"

$PDRN_edg
[1] "PDregrate"

$PTMRN_edg
[1] "PTMregrate"
```

## The regulatory complexes

In the GRN, some regulations can be performed by regulatory complexes. The composition of these complexes is stored in the `complexes` element of the `insilicosystem` object.

```r
> mysystem$complexes

$`CTL1`
[1] "4" "7"
```
This means that the products of genes 4 and 7 assemble in the system to form a regulatory complex labelled "CTL1". The kinetic parameters associated with these complexes (i.e. association and dissociation rates of the components) are available with:

```r
> mysystem$complexeskinetics

$`CTL1`
$`CTL1`$`formationrate`
[1] 0.009995349

$`CTL1`$dissociationrate
[1] 0.001773972
```
The element `complexesTargetReaction` of the `insilicosystem` object simply gives the type of regulation that the complexes accomplish, here regulation of translation:

```r
> mysystem$complexesTargetReaction

$`CTL1`
[1] "TL"
```


## The `sysargs` element

The different parameters used to generate the *in silico* system are stored in the `sysargs` element of the `insilicosystem` object. You can specify a value for each of these parameters during the construction of the system, by passing them to the function `generateInSilicoSystem`.

## Empty *in silico* system

The argument `"empty"` of the `generateInSilicoSystem` function allows you to generate a system without any regulatory interactions:

```r
> emptysystem = createInSilicoSystem(G = 7, empty = T)
> emptysystem$edg

[1] from           to             TargetReaction RegSign        RegBy         
<0 rows> (or 0-length row.names)
```

# Creating an *in silico* population

We will next create a population of *in silico* individuals. Each individual possess different copies of the genes specified in the *in silico* system generated in the previous step. You can decide the ploidy of the individuals, that is the number of copies of each gene that they will carry, the number of different variants of each gene that segregate in this in silico population, etc. For example:
```r
mypop = createInSilicoPopulation(3, mysystem, ngenevariants = 4, ploidy = 2)
```

creates a population of 3 *in silico* diploid individuals (ploidy = 2), assuming that there exist 4 different (genetically speaking) versions of each gene.  The population returned by the function is an object of class `insilicopopulation`, i.e. a list whose different attributes are presented below.

## The gene variants
A gene variant is represented as a vector containing the quantitative effects of its mutations on different kinetic properties of the gene, termed **QTL effect coefficients**. The variants segregating in this population are stored in the element `GenesVariants` of the `insilicopopulation` object returned by the function.
```r
> mypop$GenesVariants[1:2]

$`1`
              1        2         3         4
qtlTCrate     1 1.000000 0.8874961 1.0012797
qtlRDrate     1 1.000000 0.9701102 0.9554292
qtlTCregbind  1 1.000000 0.9334462 1.0101256
qtlRDregrate  1 1.000000 0.9950473 1.1977350
qtlactivity   1 1.080341 1.2668271 0.9128124
qtlTLrate     0 0.000000 0.0000000 0.0000000
qtlPDrate     0 0.000000 0.0000000 0.0000000
qtlTLregbind  0 0.000000 0.0000000 0.0000000
qtlPDregrate  0 0.000000 0.0000000 0.0000000
qtlPTMregrate 0 0.000000 0.0000000 0.0000000

$`2`
              1         2         3        4
qtlTCrate     1 1.0000000 1.0000000 1.000000
qtlRDrate     1 1.0000000 1.2642279 1.000000
qtlTCregbind  1 0.9535243 0.8852958 1.000000
qtlRDregrate  1 1.0000000 0.8183264 1.010594
qtlactivity   1 0.9140379 1.0000000 1.000000
qtlTLrate     1 1.0000000 1.0000000 1.000000
qtlPDrate     1 1.0000000 1.0000000 1.000000
qtlTLregbind  1 1.0000000 1.0000000 1.000000
qtlPDregrate  1 1.0000000 1.0000000 1.000000
qtlPTMregrate 1 1.0000000 1.0000000 1.000000
```
This shows the variants (columns of the dataframe) of genes 1 and 2 that have been generated by the function. The first variant of each gene corresponds to the "original" version of the gene, i.e. all QTL effect coefficients are set to 1. A QTL effect coefficient is a multiplicative coefficient that will be applied to the corresponding kinetic parameter of the gene during the construction of the stochastic model to simulate the expression profiles for the different individuals. As some of the QTL effect coefficients apply to translation- or protein-related steps of the gene expression, they are set to 0 for noncoding genes (see gene 1).

Here, the 2nd variant of gene 1 carries one mutation (one QTL effect coefficient different from 1). The mutation increases the activity of the gene's active product. Gene 1 being a noncoding gene (`coding = "NC"`) encoding a regulatory RNA targeting the translation of its targets (`TargetReaction = "TL"`), it means that the RNAs of this variants have a decreased affinity (i.e. binding rate) for their targets, by around 8%.

QTL effect coefficient name | Effect
--------------------------- | ------
`qtlTCrate` | Affects the basal transcription rate of the gene
`qtlRDrate` | Affects the basal RNA decay rate of the gene
`qtlTCregbind` | Affects the binding rate of the regulators of transcription on the gene's promoter (affects all transcription regulators targeting this gene)
`qtlRDregbind` | Affects the rate at which regulators of RNA decay encountering the RNAs of the gene trigger their degradation (affects all RNA decay regulators targeting this gene)
`qtlactivity` | Affects the activity of the active product of the gene. If the gene is encoding for a regulator of transcription or translation, this affects the binding rate of its active products (i.e. RNAs or proteins) to their binding sites on their targets (affects the binding to all targets of the gene). If the gene encodes a regulator of RNA or protein decay or of protein post-translational modification, this affects the rate at which its active products (i.e. RNAs or proteins) trigger the degradation/transformation of their targets (effect for all targets of the gene).
`qtlTLrate` | Affects the basal translation rate of the gene 
`qtlPDrate` | Affects the basal protein decay rate of the gene
`qtlTLregbind` | Affects the binding rate of the regulators of translation on the gene's RNA binding sites (affects all translation regulators targeting this gene)
`qtlPDregbind` | Affects the rate at which regulators of protein decay encountering the proteins of the gene trigger their degradation (affects all protein decay regulators targeting this gene)
`qtlPTMregbind` | Affects the rate at which regulators of protein post-translational modification encountering the proteins of the gene trigger their modification (affects all protein post-translational modification regulators targeting this gene)

## The *in silico* individuals

The different generated *in silico* individuals are stored in the element `individualsList` of the `insilicopopulation` object. Each individual is represented by a list with the following elements: 

- `haplotype`: gives for each gene (rows) the variants that the individual carries (columns). The different alleles of the genes are denoted "GCNi", (with `i` ranging from 1 to P - the ploidy of the organism). In this example we defined the ploidy of the individuals to be 2, so each individual carries 2 copies of each gene.
```r
> mypop$individualsList$Ind1$haplotype

   GCN1 GCN2
1     1    4
2     1    3
3     3    4
4     2    2
5     4    1
6     4    3
7     2    4
8     4    1
9     1    3
10    1    4
```
Here the first individuals (`Ind1`) carries the variants 1 and 4 of gene 1, and two copies of variant 2 of gene 4.

- `QTLeffects`: gives for each allele (i.e. "GCN1", "GCN2", etc) the value of each QTL effect coefficient for the genes (value for gene `i` at the `i`-st position in the vector of QTL effect coefficients).

```r
> mypop$individualsList$Ind1$QTLeffects$GCN1

$`GCN1`
$`GCN1`$`qtlTCrate`
 [1] 1.0000000 1.0000000 0.9571972 0.9017158 0.8814155 1.0000000 1.0951769 1.0000000 1.0000000 1.0000000

$`GCN1`$qtlRDrate
 [1] 1.0000000 1.0000000 1.1494022 1.0000000 1.0000000 1.0123999 1.0000000 0.9455702 1.0000000 1.0000000

$`GCN1`$qtlTCregbind
 [1] 1.0000000 1.0000000 0.7909081 0.8483935 1.0000000 0.9909619 1.0000000 0.8050945 1.0000000 1.0000000

$`GCN1`$qtlRDregrate
 [1] 1.0000000 1.0000000 1.0000000 0.8034261 1.0000000 1.1755527 1.0000000 0.9967105 1.0000000 1.0000000

$`GCN1`$qtlactivity
 [1] 1.0000000 1.0000000 1.0000000 1.0000000 0.8736387 1.0000000 1.0000000 0.9926014 1.0000000 1.0000000

$`GCN1`$qtlTLrate
 [1] 0.000000 1.000000 1.050869 1.000000 1.000000 1.000000 1.000000 1.000000 0.000000 1.000000

$`GCN1`$qtlPDrate
 [1] 0.0000000 1.0000000 0.8573324 1.0000000 1.0000000 1.0000000 1.0000000 0.9918177 0.0000000 1.0000000

$`GCN1`$qtlTLregbind
 [1] 0.0000000 1.0000000 1.0814660 1.0000000 0.9588442 1.0000000 0.9635237 1.0760147 0.0000000 1.0000000

$`GCN1`$qtlPDregrate
 [1] 0.0000000 1.0000000 1.0404279 0.9307016 1.0000000 0.9828826 1.0000000 1.0000000 0.0000000 1.0000000

$`GCN1`$qtlPTMregrate
 [1] 0.0000000 1.0000000 1.0713693 0.8816208 0.9678285 1.0000000 1.0707949 1.0289422 0.0000000 1.0000000
```

As individual `Ind1`'s first allele ("GCN1") of gene 1 is gene 1's variant 1 (see `mypop$GenesVariants[[1]][,1]`), the first element of each QTL effect coefficient vector for `GCN1` is `1.0` (or `0.0` for QTL effect coefficients that only affect protein-coding genes). Similarly, the 4th element of the different QTL effect coefficients vectors for `GCN1` correspond to the values in `mypop$GenesVariants[[4]][,2]` (as the first gene 4 allele of `Ind1` is gene 4's variant 2).

- `InitVar`: the list of **initial abundance variation coefficients** for the RNAs and proteins of the genes. During the simulation of the individuals gene profiles, a stochastic system is generated and an initial abundance is automatically computed for each product (RNA and protein) of the different genes. When constructing the *in silico* population, the parameter `sameInit` controls whether or not these initial abundances are the same for all individuals in the population. If not (default behaviour of the `createInSilicoPopulation` function), then each gene product is assigned an initial abundance variation coefficient that will be multiplied with the automatically computed initial abundance of the molecule, to give its initial abundance for the corresponding individual.

```r
> mypop$individualsList$Ind1$InitVar
$`GCN1`
$`GCN1`$`R`
 [1] 0.8809851 1.0474897 0.8500242 0.9365774 0.9403471 1.0005447 0.8752230 1.0001330 0.8909920 1.1357460

$`GCN1`$P
 [1] 1.1350323 1.2276211 1.1378157 0.9011492 1.0271452 0.8608892 1.0298241 1.0099649 0.9657831 0.9047659


$GCN2
$GCN2$`R`
 [1] 0.9665166 0.9213717 0.9512336 0.9348600 0.9855229 0.9706291 1.0009749 1.0637908 1.0002845 1.0288500

$GCN2$P
 [1] 1.1190859 1.1086633 1.2548449 1.0356028 1.0033309 1.0531737 0.9510813 1.1071123 1.0277100 1.0018235
```

For example, if \[RNA1<sup>GCN1</sup>\]<sub>0</sub> corresponds to the automatically computed initial abundance of the RNAs produced by the 1st allele of gene 1, then the initial abundance of the gene 1's first allele RNAs is $$\sim$$ 0.88 * \[RNA1<sup>GCN1</sup>\]<sub>0</sub> for individual `Ind1`.

## The `indargs` element

The different parameters used to generate the in silico individuals are stored in the `indargs` element of the `insilicopopulation` object. You can specify a value for each of these parameters during the construction of the system, by passing them to the function `generateInSilicoPopulation`.

# Simulating the system

Once the system and the population have been defined, we can simulate the expression of the genes in the system for each *in silico* individual. We use the function:
```r
sim = simulateInSilicoSystem(mysystem, mypop, simtime = 2000, ntrials = 5)
```

`simtime` allows you to control the simulation end time in seconds (here we simulate the expression of the genes for 500s). `ntrials` correspond to the number of repetitions of the simulation that will be computed for each individual. To speed-up the running time, Linux and MacOS users can use a parallelised version of the simulation function:

```r
sim = simulateParallelInSilicoSystem(mysystem, mypop, simtime = 2000, ntrials = 5)
```

The output of the simulation is a list of 3 elements. The element `runningtime` gives the elapsed time between the beginning and the end of the simulations (all repetitions) for each individual.

```r
> sim$runningtime

[1] 12.58 12.48 12.78
```

The `stochmodel` element is a XRJulia proxy object giving Julia object that stores the stochastic model of the system (do not try to read it, it is not really useful in its current form).

The result of the simulation, that is the abundance of the different species in the system over time for each *in silico* individual, is returned in the `Simulation` element. This is a data-frame, giving for each individual (column `Ind`), for each repetition of the simulation (column `trial`) the abundance of the different species over time (column `time`).

```r
> head(sim$Simulation)

  time trial R5GCN2 P5GCN2 R7GCN2 P7GCN2 Pm7GCN2 R3GCN1 P3GCN1 R1GCN2 R9GCN1 R6GCN2 P6GCN2 R10GCN2 P10GCN2 R1GCN1 R4GCN2
1    0     1      1   7900      1   1292       0      2  17975     49      1      1  19079      10   17565     45      1
2    1     1      1   7898      1      0       0      2  17975     43      0      1  19078      10   17564     39      1
3    2     1      1   7899      1      0       0      2  17976     41      0      1  19077      10   17566     36      1
4    3     1      1   7900      1      0       0      2  17978     40      0      1  19077      10   17566     34      1
5    4     1      1   7900      1      0       0      2  17981     40      0      1  19077      10   17565     33      1
6    5     1      1   7900      1      0       1      2  17985     40      0      1  19077      10   17566     33      1
  P4GCN2 R8GCN1 P8GCN1 R6GCN1 P6GCN1 R10GCN1 P10GCN1 R2GCN2 P2GCN2 R8GCN2 P8GCN2 R5GCN1 P5GCN1 R4GCN1 P4GCN1 R3GCN2 P3GCN2
1  47395      9 194235      1  15596      11   15864      1  99312     10 212919      1   8087      1  41242      2  19824
2  45923      9 194234      1  15595      11   15862      1  99310     10 212914      1   8085      1  40021      2  19827
3  45923      9 194231      1  15595      11   15861      1  99310     10 212913      1   8086      1  40024      2  19831
4  45925      9 194229      1  15595      11   15862      1  99311     10 212913      1   8086      1  40026      2  19833
5  45928      9 194228      1  15594      11   15862      1  99313     10 212912      1   8088      1  40023      2  19833
6  45926      9 194226      1  15594      11   15862      1  99316     10 212911      1   8088      1  40027      2  19836
  R2GCN1 P2GCN1 R9GCN2 R7GCN1 P7GCN1 Pm7GCN1 CTL1_P4GCN1_Pm7GCN2 CTL1_P4GCN2_Pm7GCN1 CTL1_P4GCN1_Pm7GCN1 CTL1_P4GCN2_Pm7GCN2
1      1 109969      1      1   1399       0                   0                   0                   0                   0
2      1 109969      1      1      0       0                 579                 759                 631                 703
3      1 109970      1      1      0       0                 579                 761                 628                 704
4      1 109971      1      1      0       0                 580                 760                 629                 703
5      1 109971      0      1      0       1                 580                 757                 631                 703
6      1 109973      0      1      0       0                 577                 759                 630                 705
   Ind
1 Ind1
2 Ind1
3 Ind1
4 Ind1
5 Ind1
6 Ind1
```

By default, the simulation distinguishes the different gene products (RNAs, proteins and regulatory complexes) according to their allele of origin (e.g. the RNAs arising from the first and second allele of gene 1 will be separately counted in the columns `R1GCN1` and `R1GCN2`, respectively). To obtain results that ignore the allele of origin, you can use:

```r
> simNoAllele = mergeAlleleAbundance(sim$Simulation)
> head(simNoAllele)

  time trial  Ind R5    P5 R7   P7 Pm7 R3    P3 R1 R9 R6    P6 R10   P10 R4    P4 R8     P8 R2     P2 CTL1_P4_Pm7
1    0     1 Ind1  2 15987  2 2691   0  4 37799 94  2  2 34675  21 33429  2 88637 19 407154  2 209281           0
2    1     1 Ind1  2 15983  2    0   0  4 37802 82  1  2 34673  21 33426  2 85944 19 407148  2 209279        2672
3    2     1 Ind1  2 15985  2    0   0  4 37807 77  1  2 34672  21 33427  2 85947 19 407144  2 209280        2672
4    3     1 Ind1  2 15986  2    0   0  4 37811 74  1  2 34672  21 33428  2 85951 19 407142  2 209282        2672
5    4     1 Ind1  2 15988  2    0   1  4 37814 73  0  2 34671  21 33427  2 85951 19 407140  2 209284        2671
6    5     1 Ind1  2 15988  2    0   1  4 37821 73  0  2 34671  21 33428  2 85953 19 407137  2 209289        2671
```

A gene product bound into a regulatory complex is not accounted for when computing the abundance for this species (e.g. if all existing proteins of gene 7 are in a regulatory complex then the abundance for `P7` will be 0). It is possible to ignore the regulatory complexes and compute the abundance of a species by counting each molecule whether it is in a free form or bound into a complex: 

```r
> simNoComplex = mergeComplexesAbundance(sim$Simulation)
> head(simNoComplex)

  time trial R5GCN2 P5GCN2 R7GCN2 P7GCN2 Pm7GCN2 R3GCN1 P3GCN1 R1GCN2 R9GCN1 R6GCN2 P6GCN2 R10GCN2 P10GCN2 R1GCN1 R4GCN2
1    0     1      1   7900      1   1292       0      2  17975     49      1      1  19079      10   17565     45      1
2    1     1      1   7898      1      0    1282      2  17975     43      0      1  19078      10   17564     39      1
3    2     1      1   7899      1      0    1283      2  17976     41      0      1  19077      10   17566     36      1
4    3     1      1   7900      1      0    1283      2  17978     40      0      1  19077      10   17566     34      1
5    4     1      1   7900      1      0    1283      2  17981     40      0      1  19077      10   17565     33      1
6    5     1      1   7900      1      0    1283      2  17985     40      0      1  19077      10   17566     33      1
  P4GCN2 R8GCN1 P8GCN1 R6GCN1 P6GCN1 R10GCN1 P10GCN1 R2GCN2 P2GCN2 R8GCN2 P8GCN2 R5GCN1 P5GCN1 R4GCN1 P4GCN1 R3GCN2 P3GCN2
1  47395      9 194235      1  15596      11   15864      1  99312     10 212919      1   8087      1  41242      2  19824
2  47385      9 194234      1  15595      11   15862      1  99310     10 212914      1   8085      1  41231      2  19827
3  47388      9 194231      1  15595      11   15861      1  99310     10 212913      1   8086      1  41231      2  19831
4  47388      9 194229      1  15595      11   15862      1  99311     10 212913      1   8086      1  41235      2  19833
5  47388      9 194228      1  15594      11   15862      1  99313     10 212912      1   8088      1  41234      2  19833
6  47390      9 194226      1  15594      11   15862      1  99316     10 212911      1   8088      1  41234      2  19836
  R2GCN1 P2GCN1 R9GCN2 R7GCN1 P7GCN1 Pm7GCN1  Ind
1      1 109969      1      1   1399       0 Ind1
2      1 109969      1      1      0    1390 Ind1
3      1 109970      1      1      0    1389 Ind1
4      1 109971      1      1      0    1389 Ind1
5      1 109971      0      1      0    1389 Ind1
6      1 109973      0      1      0    1389 Ind1

```

Lastly, non-modified and modified forms of proteins are counted separately. We merge their abundance with:
```r
> simNoPTM = mergePTMAbundance(simNoAllele)
> head(simNoPTM)

  time trial  Ind R5    P5 R7   P7 R3    P3 R1 R9 R6    P6 R10   P10 R4    P4 R8     P8 R2     P2 CTL1_P4_Pm7
1    0     1 Ind1  2 15987  2 2691  4 37799 94  2  2 34675  21 33429  2 88637 19 407154  2 209281           0
2    1     1 Ind1  2 15983  2    0  4 37802 82  1  2 34673  21 33426  2 85944 19 407148  2 209279        2672
3    2     1 Ind1  2 15985  2    0  4 37807 77  1  2 34672  21 33427  2 85947 19 407144  2 209280        2672
4    3     1 Ind1  2 15986  2    0  4 37811 74  1  2 34672  21 33428  2 85951 19 407142  2 209282        2672
5    4     1 Ind1  2 15988  2    1  4 37814 73  0  2 34671  21 33427  2 85951 19 407140  2 209284        2671
6    5     1 Ind1  2 15988  2    1  4 37821 73  0  2 34671  21 33428  2 85953 19 407137  2 209289        2671
```

All merging functions presented above can be used one after the other or independently, e.g.:

```r
> simNothing = mergeComplexesAbundance(simNoAllele)
> head(simNothing)

  time trial  Ind R5    P5 R7   P7  Pm7 R3    P3 R1 R9 R6    P6 R10   P10 R4    P4 R8     P8 R2     P2
1    0     1 Ind1  2 15987  2 2691    0  4 37799 94  2  2 34675  21 33429  2 88637 19 407154  2 209281
2    1     1 Ind1  2 15983  2    0 2672  4 37802 82  1  2 34673  21 33426  2 88616 19 407148  2 209279
3    2     1 Ind1  2 15985  2    0 2672  4 37807 77  1  2 34672  21 33427  2 88619 19 407144  2 209280
4    3     1 Ind1  2 15986  2    0 2672  4 37811 74  1  2 34672  21 33428  2 88623 19 407142  2 209282
5    4     1 Ind1  2 15988  2    0 2672  4 37814 73  0  2 34671  21 33427  2 88622 19 407140  2 209284
6    5     1 Ind1  2 15988  2    0 2672  4 37821 73  0  2 34671  21 33428  2 88624 19 407137  2 209289
```

## Plotting the simulation

It is possible to visualise the results of the simulation with:

```r
plotSimulation(sim$Simulation)
```
![plotSimulation](images/plotSimulation.png)

This returns a plot of the abundance of the different species (separated by RNAs -bottom-, proteins -middle- and regulatory complexes -top-) over time. As the simulation has been repeated 5 times (`ntrials` = 5), the mean abundance over the different repetitions or trials of the molecules is plotted as a solid lines, and the minimum and maximum values are represented by the coloured areas. By default the abundances are plotted on a log10 scale, but you can change that with the option `yLogScale = F`i in the `plotSimulation`.

By default, the different alleles are merged before plotting (`mergeAllele = T`), and similarly the non-modified and modified versions of the proteins are merged before plotting (`mergePTM = T`). On the contrary, the free and in complex components of the system are not merged (`mergeComplexes = F`).

If you want to focus on one *in silico* individual, and zoom on a smaller time-period, you can use:
```r
plotSimulation(sim$Simulation, inds = c("Ind1"), timeMin = 200, timeMax = 300)
```
![plotSimulation2](images/plotSimulation2.png)

Alernatively, you can plot the abundance of the different components as a heatmap:
```r
plotHeatMap(sim$Simulation)
```
![plotHeatMap](images/plotHeatMap.png)

The settings of this function are the same as those of the `plotSimulation` function presented above. Plotting this specific simulation as a heatmap is not really useful, but such plot can get really interesting for more complex case, as for example the anthocyanin biosynthesis regulation pathway (included as an example):
![plotcolpw](images/heatmap_colpw_tmax100.png)

## The stochastic model

 If you want to see the list of species and reactions (i.e. the stochastic model) of the system, you can use the option `writefile = T`
of the `simulate(Parallel)inSilicoSystem` functions. This generates two text files: one listing the different species in the system, each line giving a species name and its initial abundance, and one listing the biochemical reactions and associated rates. The initial abundances and reaction rates are written in a general form (i.e. giving the QTL effect coefficients/initial abundance variation coefficients to be used to compute numerical values for each individual).

### The species

As the individuals are diploid, there exist two versions of each gene (and gene product): those originating from the first allele (GCN1) and those originating from the second allele (GCN2), e.g.:

```
R1GCN1
R1GCN2	
```


The DNA sequence of genes is not explicitely modelled, except if the gene is regulated at the transcription level. In this case, the gene's DNA form is modelled as the sum of the binding sites of its different regulators. These binding sites can exist in a free or bound state. Morevoer, the binding site of a specific regulator can be occupied by the regulator's product arising from either of the regulator alleles. For example, gene 2 transcription is regulated by gene 9, so the DNA form of gene 2 first allele is:
```
Pr2GCN1reg9F	1 ## free binding site for regulator 9 on gene 2 first allele 
Pr2GCN1reg9GCN1B	0 ## binding site occupied by one of regulator 9's products originating from the first allele of gene 9
Pr2GCN1reg9GCN2B	0 ## binding site occupied by one of regulator 9's products originating from the second allele of gene 9
```
The same scheme is repeated for the second allele of gene 2:

```
Pr2GCN2reg9F	1
Pr2GCN2reg9GCN1B	0
Pr2GCN2reg9GCN2B	0
```

At the beginning of the simulation, all binding sites are in a free state (initial abundance 1 for the free form of the binding sites, 0 for the occupied forms).

The same modelling applies to the RNA form of genes. If the gene is not targeted by regulators of translation (e.g. gene 1), we simply have:
```
R1GCN1	51.19330725498743*InitVar["GCN1"]["R"][1]
R1GCN2	51.19330725498743*InitVar["GCN2"]["R"][1]
```
The initial abundance of gene 1's RNAs is $$\sim 51 \times$$ the initial abundance variation coefficient for the corresponding allele of the considered individual. If on the contrary the gene is targeted by regulator of translation, the RNA form of the gene is modelled as the sum of the RNA binding sites for the different translation regulators. One example is gene 10 whose translation is regulated by gene 1:
```
RBS10GCN1reg1F	9.385223382716024*InitVar["GCN1"]["R"][10]
RBS10GCN1reg1GCN1B	0
RBS10GCN1reg1GCN2B	0
RBS10GCN2reg1F	9.385223382716024*InitVar["GCN2"]["R"][10]
RBS10GCN2reg1GCN1B	0
RBS10GCN2reg1GCN2B	0
```
Again, at the beginning of the simulation, all RNA binding sites are in a free state, and the initial abundance of the RNAs must account for the initial abundance variation coefficient for the corresponding allele of the considered individual.

The proteins are modelled as follow:
```
P2GCN1	(89578.69434518841)*InitVar["GCN1"]["P"][2]
P2GCN2	(89578.69434518841)*InitVar["GCN2"]["P"][2]
```
If a gene is targeted in the GRN by post-translational modification, there also exists a modified form of the protein, e.g. for gene 7:
```
P7GCN1	(1358.0048533970976)*InitVar["GCN1"]["P"][7]
P7GCN2	(1358.0048533970976)*InitVar["GCN2"]["P"][7]
Pm7GCN1	0
Pm7GCN2	0
```

At the beginning of the simulation all proteins are in their original (non-modified) form.

Recall that in our system, the products of genes 4 and 7 for a regulatory complex named `CTL1`. As there exist two versions of the proteins of these genes (arising either from the first or the second allele of the genes), there exist 4 versions of the complex, including all possible combinations of the different protein allelic versions. Initially no regulatory complex is formed in the system.
```
CTL1_P4GCN1_Pm7GCN1	0
CTL1_P4GCN1_Pm7GCN2	0
CTL1_P4GCN2_Pm7GCN1	0
CTL1_P4GCN2_Pm7GCN2	0
```

### The reactions

Each reaction is characterised by a name, a biochemical formula in the form "$$\sum$$ reactants --> $$\sum$$ products", and a rate. For example the formation and dissociations of the regulatory complex CTL1 are (from the section above we know that there are 4 versions of the complex, hence 8 formation/dissociation reactions):
```
formationCTL1_P4GCN1_Pm7GCN1	P4GCN1 + Pm7GCN1 --> CTL1_P4GCN1_Pm7GCN1	0.00999534945189953
dissociationCTL1_P4GCN1_Pm7GCN1	CTL1_P4GCN1_Pm7GCN1 --> P4GCN1 + Pm7GCN1	0.00177397170383483
formationCTL1_P4GCN1_Pm7GCN2	P4GCN1 + Pm7GCN2 --> CTL1_P4GCN1_Pm7GCN2	0.00999534945189953
dissociationCTL1_P4GCN1_Pm7GCN2	CTL1_P4GCN1_Pm7GCN2 --> P4GCN1 + Pm7GCN2	0.00177397170383483
formationCTL1_P4GCN2_Pm7GCN1	P4GCN2 + Pm7GCN1 --> CTL1_P4GCN2_Pm7GCN1	0.00999534945189953
dissociationCTL1_P4GCN2_Pm7GCN1	CTL1_P4GCN2_Pm7GCN1 --> P4GCN2 + Pm7GCN1	0.00177397170383483
formationCTL1_P4GCN2_Pm7GCN2	P4GCN2 + Pm7GCN2 --> CTL1_P4GCN2_Pm7GCN2	0.00999534945189953
dissociationCTL1_P4GCN2_Pm7GCN2	CTL1_P4GCN2_Pm7GCN2 --> P4GCN2 + Pm7GCN2	0.00177397170383483
```

# Appendix
*Create an appendix to list all arguments of insilicosystemargs and insilicoindivargs*


{% include lib/mathjax.html %}
