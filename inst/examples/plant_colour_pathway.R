##########################################################################################################################################
##                                                EXAMPLE - PLANT COLOUR PATHWAY                                                        ##
##      This is an example of the use of the sismonr package, in which we model the anthocyanin biosynthesis regulation pathway,        ##
##        as described in Albert, et al. (2014). "A Conserved Network of Transcriptional Activators and Repressors Regulates            ##
##                  Anthocyanin Pigmentation in Eudicots. The Plant Cell. (https://doi.org/10.1105/tpc.113.122069)                      ##
##########################################################################################################################################


library(sismonr)
library(tidyverse)
# library(RColorBrewer)
library(wesanderson)

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
                          "TCrate" = c(2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
                          "TLrate" = c(0.2, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                          "RDrate" = c(0.2, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                          "PDrate" = c(0.02, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001))

colsystem$genes  = colsystem$genes %>% mutate(TCrate = kineticgenes$TCrate,
                                              TLrate = kineticgenes$TLrate,
                                              RDrate = kineticgenes$RDrate,
                                              PDrate = kineticgenes$PDrate)


## Adding regulatory complexes in the system
compo = list(list("compo" = c(1, 2, 3), "formationrate" = 1, "dissociationrate" = 0.1),
             list("compo" = c(1, 3, 4), "formationrate" = 2, "dissociationrate" = 0.1),
             list("compo" = c("CTC2", 5), "formationrate" = 2.5, "dissociationrate" = 0.1),
             list("compo" = c(2, 6), "formationrate" = 1.5, "dissociationrate" = 0.1),
             list("compo" = c(4, 6), "formationrate" = 1.5, "dissociationrate" = 0.1))

for(comp in compo){
  colsystem = addComplex(colsystem, comp$compo, formationrate = comp$formationrate, dissociationrate = comp$dissociationrate)
}

## Adding  regulatory interactions in the system
interactions = list(list("edge" = c("CTC1", 4), "regsign" = "1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2, "TCfoldchange" = 5)),
                    list("edge" = c("CTC2", 4), "regsign" = "1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2, "TCfoldchange" = 4)),
                    list("edge" = c("CTC2", 5), "regsign" = "1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2, "TCfoldchange" = 5)),
                    list("edge" = c("CTC2", 6), "regsign" = "1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2, "TCfoldchange" = 5)),
                    list("edge" = c("CTC2", 7), "regsign" = "1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2, "TCfoldchange" = 5)),
                    list("edge" = c("CTC3", 4), "regsign" = "-1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2)),
                    list("edge" = c("CTC3", 5), "regsign" = "-1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2)),
                    list("edge" = c("CTC3", 6), "regsign" = "-1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2)),
                    list("edge" = c("CTC3", 7), "regsign" = "-1", "kinetics" = list("TCbindingrate" = 0.1, "TCunbindingrate" = 2)))

for(inter in interactions){
  colsystem = addEdge(colsystem, inter$edge[1], inter$edge[2], regsign = inter$regsign, kinetics = inter$kinetics)
}



## ---------------------------------- ##
## Creating the in silico individuals ##
## ---------------------------------- ## ----

## We are going to simulate two different individuals or plants
## One is a wild-type plant (no mutation in any of its genes)
## The second is a mutant, in which gene 5 (the MYBrep gene) is overexpressed (here we increase its transcription rate by 5)

plants = createInSilicoPopulation(2, colsystem, sameInit = T, ngenevariants = 1)

## We add the QTL effect coefficient for the second individual such that the transcription rate of gene 5 is increased
plants$individualsList$Ind2$QTLeffects$GCN1$qtlTCrate[5] = 5
plants$individualsList$Ind2$QTLeffects$GCN2$qtlTCrate[5] = 5

## Changing the initial conditions
## As specified in Albert et al., 2014, only gene 2 and 3 (bHLH1 and WDR) are constitutively expressed (see Fig. 8).
for(g in names(plants$individualsList$Ind1$InitVar)){
  plants$individualsList$Ind1$InitVar[[g]] = list("R" = c(0, 1, 1, 0, 0, 0, 0),
                                                  "P" = c(0, 1, 1, 0, 0, 0, 0))
  plants$individualsList$Ind2$InitVar[[g]] = list("R" = c(0, 1, 1, 0, 0, 0, 0),
                                                  "P" = c(0, 1, 1, 0, 0, 0, 0))
}



## --------------------------------------------------------------- ##
## Simulating the expression profiles of the genes for both plants ##
## --------------------------------------------------------------- ## ----

sim = simulateInSilicoSystem(colsystem, plants, 1000, nepochs = 1000, ntrials = 5)

## We merge the abundance of the different allelic versions of the same genes (i.e. for each gene the abundance of RNAs - and proteins- from the two different alleles are added)
simres = mergeAlleleAbundance(sim$Simulation)



## ---------------------------------------------------------- ##
## Plotting the expression profiles over time for both plants ##
## ---------------------------------------------------------- ## ----

## Transforming the data for ggplot2
toplot = simres %>%
  gather(key = "ID", value = "Abundance", setdiff(names(simres), c("time", "trial", "Ind"))) %>%
  mutate(ID = stringr::str_replace(ID, "_.+", "")) %>%
  mutate(Type = case_when(str_detect(ID, "^R") ~ "RNAs",
                          str_detect(ID, "^P") ~ "Proteins",
                          str_detect(ID, "^C") ~ "Complexes"),
         Components = stringr::str_replace(ID, "^R|^P", "")) %>%
  mutate(Components = id2names[Components]) %>%
  group_by(Ind, time, Components, Type, ID) %>%
  summarise("mean" = mean(Abundance), "LB" = min(Abundance), "UB" = max(Abundance))


plotbreaks = c("WDR", "MYB", "bHLH1", "bHLH2", "MBW1", "MBW2", "MYBrep", "MBWr", "R3-MYB", "R3-bHLH1", "R3-bHLH2", "DFR")

## Choosing the colours
cols = c('#9a6324', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#ff0000')
names(cols) = plotbreaks
# pie(rep(1, length(cols)), col = cols, labels = names(cols))


colourpwplot = ggplot(toplot, aes(x = time)) +
  geom_ribbon(aes(ymin = LB, ymax = UB, fill = Components), alpha = 0.5) +
  geom_line(aes(y = mean, colour = Components)) +
  scale_colour_manual(values = cols, breaks = plotbreaks) +
  scale_fill_manual(values = cols, breaks = plotbreaks) +
  facet_grid(Type~Ind, scales = "free_y", labeller = as_labeller(c(`Ind1` = "Wild type", `Ind2` = "MYBrep overexpressed", `RNAs` = "RNAs", `Proteins` = "Proteins", `Complexes` = "Complexes"))) +
  xlab("Time (s)") + ylab("Components absolute abundance") +
  theme_minimal() +
  theme(legend.text = element_text(size = 7), strip.text = element_text(size = 10))

print(colourpwplot)
