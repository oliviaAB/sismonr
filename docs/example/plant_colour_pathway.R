##########################################################################################################################################
##                                                EXAMPLE - PLANT COLOUR PATHWAY                                                        ##
##      This is an example of the use of the sismonr package, in which we model the anthocyanin biosynthesis regulation pathway,        ##
##        as described in Albert, et al. (2014). "A Conserved Network of Transcriptional Activators and Repressors Regulates            ##
##                  Anthocyanin Pigmentation in Eudicots. The Plant Cell. (https://doi.org/10.1105/tpc.113.122069)                      ##
##########################################################################################################################################


library(sismonr)
library(tidyverse)

## Gene ID - name correspondence
genes.name2id = data.frame("ID" = as.character(1:7),
                           "name" = c("MYB", ## 1
                                      "bHLH1", ## 2
                                      "WDR", ## 3
                                      "bHLH2", ## 4
                                      "MYBrep", ## 5
                                      "R3-MYB", ## 6
                                      "DFR"), ## 7
                           stringsAsFactors = F)

## Complex ID - name correspondence
complexes.name2id = data.frame("ID" = sapply(1:5, function(i){paste0("CTC", i)}),
                               "name" = c("MBW1", ## CTC1
                                          "MBW2", ## CTC2
                                          "MBWr", ## CTC3
                                          "R3-bHLH1", ## CTC4
                                          "R3-bHLH2"),
                               stringsAsFactors = F) ## CTC5

id2names = c(genes.name2id$name, complexes.name2id$name)
names(id2names) = c(genes.name2id$ID, complexes.name2id$ID)



## ----------------------------- ##
## Creating the in silico system ##
## ----------------------------- ## ----

## We create a system with 7 genes, and no regulatory interactions (they will be added manually later)
colsystem = createInSilicoSystem(empty = T, G = 7, PC.p = 1, PC.TC.p = 1, PC.TL.p = 0, PC.RD.p = 0, PC.PD.p = 0, PC.PTM.p = 0, PC.MR.p = 0)

## Changing the kinetic parameters of the genes
kineticgenes = data.frame("id" = 1:7,
                          "TCrate" = c(1, 0.1, 0.1, 0.01, 0.01, 0.1, 0.5),
                          "TLrate" = c(0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001),
                          "RDrate" = c(0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                          "PDrate" = c(0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001))

colsystem$genes  = colsystem$genes %>% mutate(TCrate = kineticgenes$TCrate,
                                              TLrate = kineticgenes$TLrate,
                                              RDrate = kineticgenes$RDrate,
                                              PDrate = kineticgenes$PDrate)


## Adding regulatory complexes in the system
compo = list(list("compo" = c(1, 2, 2, 3), "formationrate" = 1, "dissociationrate" = 0.1),
             list("compo" = c(1, 3, 4,4), "formationrate" = 2, "dissociationrate" = 0.1),
             list("compo" = c("CTC2", 5), "formationrate" = 2.5, "dissociationrate" = 0.1),
             list("compo" = c(2, 6), "formationrate" = 1.5, "dissociationrate" = 0.1),
             list("compo" = c(4, 6), "formationrate" = 1.5, "dissociationrate" = 0.1))
for(comp in compo){
  colsystem = addComplex(colsystem, comp$compo, formationrate = comp$formationrate, dissociationrate = comp$dissociationrate)
}

## Adding  regulatory interactions in the system
interactions = list(list("edge" = c("CTC1", 4), "regsign" = "1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2, "TCfoldchange" = 50)),
                    list("edge" = c("CTC2", 4), "regsign" = "1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2, "TCfoldchange" = 4)),
                    list("edge" = c("CTC2", 5), "regsign" = "1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2, "TCfoldchange" = 50)),
                    list("edge" = c("CTC2", 6), "regsign" = "1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2, "TCfoldchange" = 50)),
                    list("edge" = c("CTC2", 7), "regsign" = "1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2, "TCfoldchange" = 10)),
                    list("edge" = c("CTC3", 4), "regsign" = "-1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2)),
                    list("edge" = c("CTC3", 5), "regsign" = "-1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2)),
                    list("edge" = c("CTC3", 6), "regsign" = "-1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2)),
                    list("edge" = c("CTC3", 7), "regsign" = "-1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2)))

for(inter in interactions){
  colsystem = addEdge(colsystem, inter$edge[1], inter$edge[2], regsign = inter$regsign, kinetics = inter$kinetics)
}


## To visualise the GRN
plotGRN(colsystem)

## ---------------------------------- ##
## Creating the in silico individuals ##
## ---------------------------------- ## ----

## We are going to simulate two different individuals or plants
## One is a wild-type plant (no mutation in any of its genes)
## The second is a mutant, in which gene 5 (the MYBrep gene) is overexpressed (here we increase its transcription rate by 5)

plants = createInSilicoPopulation(3, colsystem, sameInit = T, ngenevariants = 1)

## We add the QTL effect coefficient for the second individual such that the transcription rate of gene 5 is increased
plants$individualsList$Ind2$QTLeffects$GCN1$qtlTCrate[5] = 50
plants$individualsList$Ind2$QTLeffects$GCN2$qtlTCrate[5] = 50
plants$individualsList$Ind2$QTLeffects$GCN1$qtlTCregbind[5] = 0
plants$individualsList$Ind2$QTLeffects$GCN2$qtlTCregbind[5] = 0

plants$individualsList$Ind3$QTLeffects$GCN1$qtlRDrate[5] = 6
plants$individualsList$Ind3$QTLeffects$GCN2$qtlRDrate[5] = 6
## Changing the initial conditions
## As specified in Albert et al., 2014, only gene 2 and 3 (bHLH1 and WDR) are constitutively expressed (see Fig. 8).
for(g in names(plants$individualsList$Ind1$InitVar)){
  plants$individualsList$Ind1$InitVar[[g]] = list("R" = c(0, 1, 1, 0, 0, 0, 0),
                                                  "P" = c(0, 1, 1, 0, 0, 0, 0))
  plants$individualsList$Ind2$InitVar[[g]] = list("R" = c(0, 1, 1, 0, 0, 0, 0),
                                                  "P" = c(0, 1, 1, 0, 0, 0, 0))
  plants$individualsList$Ind3$InitVar[[g]] = list("R" = c(0, 1, 1, 0, 0, 0, 0),
                                                  "P" = c(0, 1, 1, 0, 0, 0, 0))
}


## -------------------------------------------------------------- ##
## Simulating the expression profiles of the genes for all plants ##
## -------------------------------------------------------------- ## ----

sim = simulateInSilicoSystem(colsystem, plants, 1000, nepochs = 1000, ntrials = 10)

## We merge the abundance of the different allelic versions of the same genes (i.e. for each gene the abundance of RNAs - and proteins- from the two different alleles are added)
simres = mergeAlleleAbundance(sim$Simulation)



## --------------------------------------------------------- ##
## Plotting the expression profiles over time for all plants ##
## --------------------------------------------------------- ## ----

plotSimulation(sim$Simulation)

## Plot as a heatmap
plotHeatMap(sim$Simulation)

## Plot shown in the sismonr online documentation https://oliviaab.github.io/sismonr/#plotting-the-simulation
plotHeatMap(sim$Simulation, timeMax = 100, text = element.text(size = 8), strip.text.y = element.text(size = 8))