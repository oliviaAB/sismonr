# R package sismonr

*A R package for generating and simulating in silico biological systems.*

## Introduction

*Explain the concepts of in silico system and individual*

## Creating an *in silico* system

The first step is to generate an *in silico* system. An *in silico* system is composed of a set of genes, and a gene regulatory network or GRN describing the different regulatory interactions occuring between the genes. The network is created by a call to the function `createInSilicoSystem`. The user can control different aspects of the system with the arguments passed to the function. For example,
```r
mysystem = createInSilicoSystem(G = 10, PC.p = 0.7)
```
generates an *in silico* system with 10 genes,and during the generation process each of the genes has a probability of 0.7 to be designated protein-coding gene (as opposed to noncoding gene). The system returned by the function is a list of class `createInSilicoSystem`. The different attributes of the system are presented below.

### The list of genes 

The different genes consituting the system are stored in a data-frame, and can be accessed with:
```r
> mysystem$genes

   id coding TargetReaction PTMform ActiveForm       TCrate     TLrate       RDrate       PDrate
1   1     PC             PD       0         P1 0.0063075448 0.12473471 5.092704e-05 1.295350e-05
2   2     PC             TC       0         P2 0.0034184515 0.05547504 1.340618e-03 2.507492e-06
3   3     PC             TC       0         P3 0.0003385266 1.44035103 1.847389e-04 2.215386e-06
4   4     PC             TL       0         P4 0.0003175917 9.36918946 3.064079e-04 5.333123e-06
5   5     PC             PD       0         P5 0.0032860392 0.14981765 6.464755e-04 7.416136e-06
6   6     PC             TL       0         P6 0.0003027845 0.09037380 3.936931e-03 2.394177e-06
7   7     PC             RD       0         P7 0.0010848779 0.18605544 1.911332e-04 3.115656e-06
8   8     NC             RD       0         R8 0.0003252939 0.00000000 1.129150e-03 0.000000e+00
9   9     PC             TC       0         P9 0.0002982142 3.02971483 1.483975e-03 3.205987e-06
10 10     NC             RD       0        R10 0.0045610720 0.00000000 1.421933e-03 0.000000e+00
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
- `PTMform`: does the protein of the gene has a modified form? Will always be `"0"` for noncoding genes, and for now is also `"0"` for protein-coding gene. This will be modified later;
- `ActiveForm`: what is the active form of the gene? If the gene is noncoding, `ActiveForm = R[gene ID]`, and for a protein-coding gene `ActiveForm = P[gene ID]`
- `TCrate`, `TLrate`, `RDrate` and `PDrate`: give the transcription, translation, RNA decay and protein decay rates of the genes, respectively. `TLrate` and `PDrate` are set to 0 for noncoding genes.

## Appendix
*Create an appendix to list all arguments of insilicosystemargs and insilicoindivargs*
