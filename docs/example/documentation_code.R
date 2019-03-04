library(sismonr)

## ---------------- ##
## In silico system ##
## ---------------- ## ----

mysystem = createInSilicoSystem(G = 10, PC.p = 0.7)

## Or, to use the same system as presented as example:
load("mysystem.RData")

mysystem$genes

mysystem$edg

names(mysystem$mosystem)
mysystem$mosystem$TCRN_edg

lapply(mysystem$mosystem, function(x){colnames(x)[-(1:5)]})

mysystem$complexes

mysystem$complexeskinetics

mysystem$complexesTargetReaction

emptysystem = createInSilicoSystem(G = 7, empty = T)
emptysystem$edg

## -------------------- ##
## In silico population ##
## -------------------- ## ----

mypop = createInSilicoPopulation(3, mysystem, ngenevariants = 4, ploidy = 2)

## Or, to use the same population as presented as example:
load("mypop.RData")

mypop$GenesVariants[1:2]

mypop$individualsList$Ind1$haplotype

mypop$individualsList$Ind1$QTLeffects$GCN1

mypop$individualsList$Ind1$InitVar


## ---------- ##
## Simulation ##
## ---------- ## ----

sim = simulateInSilicoSystem(mysystem, mypop, simtime = 2000, ntrials = 5)
#save(sim, file = "sim.RData")
load("sim.RData")

sim$runningtime

head(sim$Simulation)

simNoAllele = mergeAlleleAbundance(sim$Simulation)
head(simNoAllele)

simNoComplex = mergeComplexesAbundance(sim$Simulation)
head(simNoComplex)

simNoPTM = mergePTMAbundance(simNoAllele)
head(simNoPTM)

simNothing = mergeComplexesAbundance(simNoAllele)
head(simNothing)


plotSimulation(sim$Simulation)
ggplot2::ggsave("plotSimulation.png")

plotSimulation(sim$Simulation, inds = c("Ind1"), timeMin = 200, timeMax = 300)
ggplot2::ggsave("plotSimulation2.png")

plotHeatMap(sim$Simulation)
ggplot2::ggsave("plotHeatMap.png")
