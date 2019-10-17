library(sismonr)

## ---------------------------- ##
## Creating an in silico system ##
## ---------------------------- ## ----

mysystem = createInSilicoSystem(G = 10, PC.p = 0.7, ploidy = 2)

# Or, to use the same system as presented in the tutorial,
# download sim.RData from https://github.com/oliviaAB/sismonr/blob/master/docs/example/mysystem.Rdata:
load("mysystem.RData")


## The list of genes

mysystem$genes


## The Gene Regulatory Network

mysystem$edg

plotGRN(mysystem)
plotGRN(mysystem, plotType = "interactive2D")
plotGRN(mysystem, edge.arrow.size = 0.5)

names(mysystem$mosystem)
mysystem$mosystem$TCRN_edg

lapply(mysystem$mosystem, function(x){colnames(x)[-(1:5)]})

plotGRN(mysystem, "TC")


## The regulatory complexes

mysystem$complexes

mysystem$complexeskinetics

mysystem$complexesTargetReaction


## Modifying the in silico system

mysystem2 = addGene(mysystem)
mysystem2$genes

mysystem2 = addGene(mysystem, coding = "PC", TargetReaction = "TL", TCrate = 0.005, PDrate = 0.0007)
mysystem2$genes

mysystem2 = addEdge(mysystem, regID = 3, tarID = 5, regsign = "1",
                    kinetics = list("TLbindingrate" = 0.01,
                                    "TLunbindingrate" = 0.1,
                                    "TLfoldchange" = 10))
tail(mysystem$edg) ## Original system
tail(mysystem2$edg) ## System with the new edge
mysystem$mosystem$TLRN_edg ## Original system
mysystem2$mosystem$TLRN_edg ## System with the new edge

mysystem2 = removeEdge(mysystem2, regID = 3, tarID = 5)
mysystem2$mosystem$TLRN_edg

mysystem2 = addComplex(mysystem, compo = c(6, 9))
mysystem2$complexes

mysystem2 = addComplex(mysystem, compo = c(6, 7))

mysystem2 = removeComplex(mysystem, "CTC1")
mysystem2$complexes
mysystem$mosystem$TCRN_edg ## Original system
mysystem2$mosystem$TCRN_edg ## System without the complex

emptysystem = createInSilicoSystem(G = 7, empty = T)
emptysystem$edg



## ----------------------------- ##
## Creating an silico population ##
## ----------------------------- ## ----

mypop = createInSilicoPopulation(3, mysystem, ngenevariants = 4)

# Or, to use the same population as presented in the tutorial,
# download sim.RData from https://github.com/oliviaAB/sismonr/blob/master/docs/example/mypop.RData
load("mypop.RData")


## The gene variants
mypop$GenesVariants[1:2]


## The in silico individuals

mypop$individualsList$Ind1$haplotype

mypop$individualsList$Ind1$QTLeffects$GCN1

mypop$individualsList$Ind1$InitAbundance

plotMutations(mypop, mysystem, nGenesPerRow = 5)

plotMutations(mypop, mysystem, scaleLims = c(0.85, 1.15),
              qtlEffectCoeffs = c("qtlTCrate", "qtlTLrate", "qtlRDrate", "qtlPDrate"),
              inds = c("Ind1", "Ind2"),
              alleles = "GCN2",
              genes = 1:3)

mypop$individualsList$Ind2$QTLeffects$GCN2$qtlTLrate[2]

## --------------------- ##
## Simulating the system ##
## --------------------- ## ----

sim = simulateInSilicoSystem(mysystem, mypop, simtime = 2000, ntrials = 5)

# Or, to use the same simulation as presented in the tutorial,
# download sim.RData from https://github.com/oliviaAB/sismonr/blob/master/docs/example/sim.RData
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


## Plotting the simulation

plotSimulation(sim$Simulation)
# ggplot2::ggsave("plotSimulation.png", width = 20, height = 20, units = "cm")

plotSimulation(sim$Simulation, mergePTM = F)
# ggplot2::ggsave("plotSimulationPTM.png", width = 20, height = 20, units = "cm")

sumtable = summariseSimulation(sim$Simulation)
head(sumtable)

plotSimulation(sim$Simulation, inds = c("Ind1"), timeMin = 200, timeMax = 300)
# ggplot2::ggsave("plotSimulation2.png", width = 20, height = 20, units = "cm")

sumtable = summariseSimulation(sim$Simulation, inds = c("Ind1"), timeMin = 200, timeMax = 300)
head(sumtable)

plotHeatMap(sim$Simulation)
# ggplot2::ggsave("plotHeatMap.png")

sim = simulateInSilicoSystem(mysystem, mypop, simtime = 100, ntrials = 5, writefile = T, filepath = getwd())
