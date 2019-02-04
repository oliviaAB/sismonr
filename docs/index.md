# R package sismonr

*A R package for generating and simulating in silico biological systems.*

## Introduction

*Explain the concepts of in silico system and individual*

### Abbreviations

## Creating an *in silico* system

The first step is to generate an *in silico* system. An *in silico* system is composed of a set of genes, and a gene regulatory network or GRN describing the different regulatory interactions occuring between the genes. The network is created by a call to the function `createInSilicoSystem`. The user can control different aspects of the system with the arguments passed to the function. For example,
```r
mysystem = createInSilicoSystem(G = 10, PC.p = 0.7)
```
generates an *in silico* system with 10 genes,and during the generation process each of the genes has a probability of 0.7 to be designated protein-coding gene (as opposed to noncoding gene). The system returned by the function is a list of class `createInSilicoSystem`. The different attributes of the system are presented below.

### The list of genes 

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

Each gene is labeled with an ID (column `id`) and possess the following parameters:

- `coding`: gives the **coding status** of the gene, i.e. describes if the gene is protein-coding (`coding = "PC"`) or noncoding (`coding = "NC"`);
- `TargetReaction`: describes the biological function of the gene, i.e. the type of regulation that it performs on its target. This parameter can take one of the following values:
  - `"TC"`: regulator of transcription
  - `"TL"`: regulator of translation
  - `"RD"`: regulator of RNA decay
  - `"PD"`: regulator of protein decay
  - `"PTM"`: regulator of protein post-translational modification
  - `"MR"`: metabolic enzyme (only for protein-coding genes, indicates that the gene cannot regulate the expression of another gene)
- `PTMform`: does the protein of the gene has a modified form? Will always be `"0"` for noncoding genes.
- `ActiveForm`: what is the active form of the gene? If the gene is noncoding, `ActiveForm = R[gene ID]`, and for a protein-coding gene `ActiveForm = "P[gene ID]"`. If the protein of a protein-coding gene is targeted for post-translational modification (`PTMform = "1"`) then we assume that only the modified form of the protein is active, and thus `ActiveForm = Pm[gene ID]`.
- `TCrate`, `TLrate`, `RDrate` and `PDrate`: give the transcription, translation, RNA decay and protein decay rates of the genes, respectively. `TLrate` and `PDrate` are set to 0 for noncoding genes.

### The GRN

The regulatory network describing the regulatory interaction between the genes are stored in a data-frame, and accessible with:
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

- `from`: ID of the regulator gene. Note that for the edge 15, the regulator is not a gene but a regulatory complex (identifiable by its ID starting with `'C'`);
- `RegBy`: type of regulator (`"PC"` for a protein-coding regulator, `"NC"` for a noncoding regulator and `"C"` for a regulatory complex);
- `to`: ID of the target gene;
- `TargetReaction`: type of the regulation, i.e. which expression step of the target is controlled? An edge for which `TargetReaction = "TC"` represents a regulation of transcription, etc (see the [Abbreviations](###Abbreviations) section);
- `RegSign`: sign of the regulation (`"1"` for an activation and `"-1"` for a repression). Edges corresponding to the regulation of RNA or protein decay always have `RegSign = "1"`, meaning that the regulator increases the decay rate of the target.

This shows the global GRN, with the different types of regulations. The element `mosystem` of the `insilicosystem` object contains the same edges but grouped by type of regulation:
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
Here you can only see the edges of the GRN corresponding to the regulation of transcription. In each of the sub-GRNs stored in the `mosystem` element contains the kinetic parameters associated with each edge. For example as we can see here each edge is assigned a binding rate (`TCbindingrate`) and an unbinding rate (`TCunbindingrate`) of the regulator to and from the binding site on the target gene's promoter. The parameter `TCfoldchange` corresponds to the coefficient by which is multiplied the basal transcription rate of the target when the regulator is bound to its binding site (notice that for edges for which `RegSign = "-1"`, i.e. corresponding to a repression, `TCfoldchange = 0`).

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

## Creating an *in silico* population
## Appendix
*Create an appendix to list all arguments of insilicosystemargs and insilicoindivargs*
