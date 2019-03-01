##########################################################################################################################################
##                                                EXAMPLE - PLANT COLOUR PATHWAY                                                        ##
##      This is an example of the use of the sismonr package, in which we model the anthocyanin biosynthesis regulation pathway,        ##
##        as described in Albert, et al. (2014). "A Conserved Network of Transcriptional Activators and Repressors Regulates            ##
##                  Anthocyanin Pigmentation in Eudicots. The Plant Cell. (https://doi.org/10.1105/tpc.113.122069)                      ##
##########################################################################################################################################


library(sismonr)
library(tidyverse)

rm(list = ls(all = T))

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



## --------------------------------------------------------------- ##
## Simulating the expression profiles of the genes for both plants ##
## --------------------------------------------------------------- ## ----

sim = simulateParallelInSilicoSystem(colsystem, plants, 2000, nepochs = 2000, ntrials = 100)

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
  summarise("mean" = mean(Abundance), "LB" = quantile(Abundance, probs = c(0.025)), "UB" = quantile(Abundance, probs = c(0.975)))


plotbreaks = c("WDR", "MYB", "bHLH1", "bHLH2", "MBW1", "MBW2", "MYBrep", "MBWr", "R3-MYB", "R3-bHLH1", "R3-bHLH2", "DFR")

## Choosing the colours
cols = c('#9a6324', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#ff0000')
names(cols) = plotbreaks
# pie(rep(1, length(cols)), col = cols, labels = names(cols))

colourpwplot = ggplot(toplot, aes(x = time)) +
  geom_ribbon(aes(ymin = LB, ymax = UB, fill = Components), alpha = 0.5, show.legend = F) +
  geom_line(aes(y = mean, colour = Components), show.legend = F) +
  scale_colour_manual(values = cols, breaks = plotbreaks) +
  scale_fill_manual(values = cols, breaks = plotbreaks) +
  scale_y_log10() +
  facet_grid(Type~Ind, scales = "free_y", labeller = as_labeller(c(`Ind1` = "Wild type", `Ind2` = "MYBrep overexpressed", `Ind3` = "MYBrep silenced", `RNAs` = "RNAs", `Proteins` = "Proteins", `Complexes` = "Complexes"))) +
  xlab("Time (s)") + ylab("Components absolute abundance") +
  theme_minimal() +
  theme(legend.text = element_text(size = 7), strip.text = element_text(size = 10), legend.direction = "horizontal", legend.position = "bottom")

## create legend

nbg = length(genes.name2id$name)
nbc = length(complexes.name2id$name)
legend_points = rbind(data.frame("Components" = rep(genes.name2id$name, each = 2), "x" = rep(1:nbg, each = 2), "y" = rep(1:2, nbg)),
                      data.frame("Components" = complexes.name2id$name, "x" = (1:nbc)+nbg, "y" = rep(3, nbc)))
legend_text = rbind(data.frame("Names" = c(genes.name2id$name, complexes.name2id$name), "x" = 1:(nbg+nbc), "y" = 4.1, "angle" = 30, "size" = 3, "hjust" = 0.5),
                    data.frame("Names" = c("RNAs", "Proteins", "Complexes"), "x" = 0.2, "y" = 1:3, "angle" = 0, "size" = 4, "hjust" = 1))
legend_lines = rbind(data.frame("x" = seq(0.5, 12.5, by = 1), "xend" = seq(0.5, 12.5, by = 1), "y" = 0.5, "yend" = 3.9),
                     data.frame("x" = 0.4, "xend" = 12.6, "y" = seq(0.5, 3.5, by = 1), "yend" = seq(0.5, 3.5, by = 1)))

legendplot = ggplot() +
  geom_point(data = legend_points, size = 5, shape = 15, aes(x = x, y = y, group = Components, colour = Components), show.legend = F) +
  geom_text(data = legend_text, aes(label = Names, x = x, y = y, angle = angle, size = size, hjust = hjust), show.legend = F) +
  geom_segment(data = legend_lines, colour = "gray90", aes(x = x, xend = xend, y = y, yend = yend)) +
  scale_colour_manual(values = cols, breaks = plotbreaks) +
  scale_size(range = c(2.5, 3.5), guide = F) +
  xlim(-1.6, 13) + ylim(0.45, 4.5) +
  theme_void()

finalplot = ggdraw() +
  draw_plot(colourpwplot, x = 0.15, y = 0.2 , width = 0.7, height = 0.8) +
  draw_plot(legendplot, x = 0.175, y = 0, width = 0.65, height = 0.2)
print(finalplot)

## --------------------------------------------------------------------------------- ##
#### Computing the experimental and simulated RNA concentration ratios for the genes ##
## --------------------------------------------------------------------------------- ## ----

expdata = data.frame("simID" = c("WDR", "MYB", "bHLH1", "bHLH2", "MYBrep", "R3-MYB", "DFR"),
           "expID" = c("AN11", "PHZ", "JAF13", "AN1", "MYB27", "MYBx", "DFR"),
           "WT" = c(0.016, 0.012, 0.002, 0.0023, 0.075, 0.0015, 0.1),
           "OE" = c(0.025, 0.014, 0.0018, 0.00026, 3, 0.00008, 0.01),
           "RNAi" = c(0.05, 0.012, 0.0014, 0.005, 0.025, 0.05, 0.33))
expdata = mutate(expdata, ratioOE = OE/WT, ratioRNAi = RNAi/WT)

simdata = toplot %>% ungroup %>%
  filter(time == max(toplot$time), Type == "RNAs") %>%
  select(Ind, Components, mean) %>%
  mutate(mean = mean) %>%
  spread(Ind, mean) %>%
  rename(WT = Ind1, OE = Ind2, RNAi = Ind3) %>%
  mutate(ratioOE = OE/WT, ratioRNAi = RNAi/WT)

ratios = merge(expdata, simdata, by.x = "simID", by.y = "Components", suffixes = c("_exp", "_sim")) %>%
  select(simID, expID, ratioOE_exp, ratioOE_sim, ratioRNAi_exp, ratioRNAi_sim)
ratios

