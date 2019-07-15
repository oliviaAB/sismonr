context("Stochastic simulation")
library(sismonr)
library(XRJulia)

check_julia = function(){
  if (!findJuliaNoError(test = T)) {
    skip("Julia not installed.")
  }else{
    ## Test that the correct version of Julia is installed
    julia_bin = XRJulia::findJulia()
    if (.Platform$OS.type == "windows") cmd_args = "-E \"VERSION < v\\\"1.0.0\\\"\"" else cmd_args = "-E 'VERSION < v\"1.0.0\"'"
    version_check = tryCatch(system2(julia_bin, args = cmd_args, stdout = TRUE, stderr = TRUE), warning = function(e) "", error = function(e) "")
    if(version_check == ""){
      skip("Error when checking the existing Julia version. Please check that Julia is correctly installed.")
    }
    else if(base::grepl("true", version_check, ignore.case = TRUE)){
      skip("Julia version < v1.0 installed, require v1.0 or later")
    }
  }
}

test_that("creation of stochastic system from empty system works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 5, empty = T, PC.p = 1)
  indargs = insilicoindividualargs(ploidy = 3)
  stochsys = createStochSystem(mysystem, indargs)
  dictkeys = juliaGet(juliaEval("collect(keys(%s))", stochsys))
  species = unlist(juliaGet(juliaEval("%s[\"species\"]", stochsys)))
  reactions = unlist(juliaGet(juliaEval("%s[\"reactions\"]", stochsys)))

  mysystemNC = createInSilicoSystem(G = 5, empty = T, PC.p = 0)
  stochsysNC = createStochSystem(mysystemNC, indargs)
  speciesNC = unlist(juliaGet(juliaEval("%s[\"species\"]", stochsysNC)))
  reactionsNC = unlist(juliaGet(juliaEval("%s[\"reactions\"]", stochsysNC)))

  expect_equal(dictkeys, c("reactionsnames", "propensities", "TCproms", "TLproms", "reactions", "initialconditions", "species"))
  expect_equal(juliaEval("length(%s[\"initialconditions\"])", stochsys), length(species))
  expect_equal(juliaEval("length(%s[\"reactionsnames\"])", stochsys), length(reactions))
  expect_equal(juliaEval("length(%s[\"propensities\"])", stochsys), length(reactions))
  expect_equal(juliaGet(juliaEval("[length(v) for v in values(%s[\"TCproms\"])]", stochsys)), rep(0, 3*5))
  expect_equal(juliaGet(juliaEval("[length(v) for v in values(%s[\"TLproms\"])]", stochsys)), rep(0, 3*5))
  expect_equal(length(species), 5*3*2) ## RNA and protein species of each gene in 3 copies (ploidy = 3)
  expect_equal(length(reactions), 5*4*3) ## TC, TL, RD and PD reactions for each of the 3 copies of each gene
  expect_equal(length(speciesNC), 5*3) ## RNA species of each gene in 3 copies (ploidy = 3)
  expect_equal(length(reactionsNC), 5*2*3) ## TC, and RD and reactions for each of the 3 copies of each gene
})

test_that("creation of stochastic system for TC reaction works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 3, empty = T, PC.p = 1, PC.TC.p = 1)
  mysystem = addEdge(mysystem, 1, 3, regsign = "1")
  mysystem = addEdge(mysystem, 2, 3, regsign = "-1")
  indargs = insilicoindividualargs(ploidy = 2)
  stochsys = createStochSystem(mysystem, indargs)
  species = unlist(juliaGet(juliaEval("%s[\"species\"]", stochsys)))
  reactions = unlist(juliaGet(juliaEval("%s[\"reactions\"]", stochsys)))

  expect_equal(length(species), 3*2*2 + ## RNA and protein species of each gene in 2 copies (ploidy = 2)
                 1*2*2*(1+2)) ## 1 binding site per regulator per copy of gene 2 that can be free or bound by either of the 2 copies of the regulator
  expect_equal(length(reactions), 4*2*2 + ## TC, TL, RD and PD rate for each of the 2 copies of each regulator gene (genes 1 and 2)
                 2*2*2*2 + # binding + unbinding reaction of each of the 2 copies of each regulator for each copy of the target
                 3*2 + # TC of gene 3: combination of 2 binding sites each in one of 3 possible states, but only the free state of binding site for reg 2 is active
                 3*2) # TL, RD and PD of gene 3
})

test_that("creation of stochastic system for TL reaction works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 3, empty = T, PC.p = 1, PC.TL.p = 1)
  mysystem = addEdge(mysystem, 1, 3, regsign = "1")
  mysystem = addEdge(mysystem, 2, 3, regsign = "-1")
  indargs = insilicoindividualargs(ploidy = 2)
  stochsys = createStochSystem(mysystem, indargs)
  species = unlist(juliaGet(juliaEval("%s[\"species\"]", stochsys)))
  reactions = unlist(juliaGet(juliaEval("%s[\"reactions\"]", stochsys)))

  expect_equal(length(species), 2*2*2 + ## RNA and protein species of each regulator in 2 copies (ploidy = 2)
                 2 + ## protein species of the 2 copies of the target
                 1*2*2*(1+2)) ## 1 binding site per regulator per copy of gene 2 that can be free or bound by either of the 2 copies of the regulator
  expect_equal(length(reactions), 4*2*2 + ## TC, TL, RD and PD rate for each of the 2 copies of each regulator gene (genes 1 and 2)
                 2*2*2*2 + # binding + unbinding reaction of each of the 2 copies of each regulator for each copy of the target
                 3*2 + # TL of gene 3: combination of 2 binding sites each in one of 3 possible states, but only the free state of binding site for reg 2 is active
                 3*2) # TC, RD and PD of gene 3
})

test_that("creation of stochastic system for RD reaction works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 3, empty = T, PC.p = 1, PC.RD.p = 1)
  mysystem = addEdge(mysystem, 1, 3)
  mysystem = addEdge(mysystem, 2, 3)
  indargs = insilicoindividualargs(ploidy = 2)
  stochsys = createStochSystem(mysystem, indargs)
  species = unlist(juliaGet(juliaEval("%s[\"species\"]", stochsys)))
  reactions = unlist(juliaGet(juliaEval("%s[\"reactions\"]", stochsys)))

  expect_equal(length(species), 3*2*2) ## RNA and protein species of each gene in 2 copies (ploidy = 2)
  expect_equal(length(reactions), 4*2*3 + ## TC, TL, RD and PD rate for each of the 2 copies of each gene
                 2*2*2) ## 2 copies of 2 regulators triggering the decay of each of the 2 copies of the target protein
})


test_that("creation of stochastic system for PD reaction works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 3, empty = T, PC.p = 1, PC.PD.p = 1)
  mysystem = addEdge(mysystem, 1, 3)
  mysystem = addEdge(mysystem, 2, 3)
  indargs = insilicoindividualargs(ploidy = 2)
  stochsys = createStochSystem(mysystem, indargs)
  species = unlist(juliaGet(juliaEval("%s[\"species\"]", stochsys)))
  reactions = unlist(juliaGet(juliaEval("%s[\"reactions\"]", stochsys)))

  expect_equal(length(species), 3*2*2) ## RNA and protein species of each gene in 2 copies (ploidy = 2)
  expect_equal(length(reactions), 4*2*3 + ## TC, TL, RD and PD rate for each of the 2 copies of each gene
                 2*2*2) ## 2 copies of 2 regulators triggering the decay of each of the 2 copies of the target RNA
})

test_that("creation of stochastic system for PTM reaction works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 3, empty = T, PC.p = 1, PC.PTM.p = 1)
  mysystem = addEdge(mysystem, 1, 3)
  mysystem = addEdge(mysystem, 2, 3)
  indargs = insilicoindividualargs(ploidy = 2)
  stochsys = createStochSystem(mysystem, indargs)
  species = unlist(juliaGet(juliaEval("%s[\"species\"]", stochsys)))
  reactions = unlist(juliaGet(juliaEval("%s[\"reactions\"]", stochsys)))

  expect_equal(length(species), 3*2*2 + ## RNA and protein species of each gene in 2 copies (ploidy = 2)
                 2) ## a PTM version of each version of the proteins of gene 3
  expect_equal(length(reactions), 4*2*3 + ## TC, TL, RD and PD rate for each of the 2 copies of each gene
                 2*2*2 + ## 2 copies of 2 regulators triggering the transformation of each of the 2 copies of the target RNA
                 2) ## Decay reaction for the PTM versions of the proteins
})

test_that("creation of stochastic system for regulatory complex works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 3, empty = T, PC.p = 1, PC.TC.p = 1)
  mysystem = addComplex(mysystem, c(1, 2))
  mysystem = addEdge(mysystem, "CTC1", 3, regsign = "1")
  indargs = insilicoindividualargs(ploidy = 2)
  stochsys = createStochSystem(mysystem, indargs)
  species = unlist(juliaGet(juliaEval("%s[\"species\"]", stochsys)))
  reactions = unlist(juliaGet(juliaEval("%s[\"reactions\"]", stochsys)))

  expect_equal(length(species), 3*2*2 + ## RNA and protein species of each gene in 2 copies (ploidy = 2)
                 2*2 + ## the regulatory complex has 4 versions: its 2 components each have 2 versions
                 2*5) ## binding site for the complex on each version of the target gene - each either in free form or bound by one of the
                      ## 4 versions of the complex
  expect_equal(length(reactions), 4*2*2 + ## TC, TL, RD and PD rate for each of the 2 copies of each regulatory gene
                 4*2 + ## association and dissociation reaction for each version of the complex
                 4*2*2 + ## binding and unbinding reaction of each version of the complex to each of the 2 versions of the target
                 5*2  + ## TC reaction with each of the 5 forms of the binding site for the 2 versions of the target gene
                 3*2) ## TL, RD and PD reactions for both versions of the target gene
})

test_that("simulation of in silico system works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 5, regcomplexes = "none", empty = T, PC.p = 1)
  mypop = createInSilicoPopulation(3, mysystem, ploidy = 2)
  sim = simulateInSilicoSystem(mysystem, mypop, simtime = 10, ntrials = 2, nepochs = 5)

  expect_equal(dim(sim$Simulation), c(2*3*6, ## 2 trials, for 3 individuals, with 5+1 time-points recorded
    20+3)) ## 20 species + column for time, ind and trial

})

test_that("merging functions for simulation results works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 3, empty = T, PC.p = 1, PC.PTM.p = 1)
  mysystem = addComplex(mysystem, c(1, 2))
  mysystem = addEdge(mysystem, "CPTM1", 3, regsign = "1")
  mypop = createInSilicoPopulation(1, mysystem, ploidy = 2)
  sim = simulateInSilicoSystem(mysystem, mypop, simtime = 10, ntrials = 1, nepochs = 5)
  simNoAllele = mergeAlleleAbundance(sim$Simulation)
  simNoPTM = mergePTMAbundance(sim$Simulation)
  simNoComplex = mergeComplexesAbundance(sim$Simulation)

  expect_equal(ncol(sim$Simulation), 3*2*2 + 4 + 2 + 3) ##RNA and protein form of each of the 2 versions of each gene + 4 versions of the complex
                                                        ## + PTM version of each of the 2 versions of the protein of gene 3 + column time, ind, trial
  expect_equal(ncol(simNoAllele), 3*2 + 1 + 1 + 3) ##RNA and protein form of each gene + 1 complex + 1 PTM form of protein 3 + columns time, ind, trial
  expect_equal(ncol(simNoPTM), 3*2*2 + 4  + 3) ##RNA and protein form of each of the 2 versions of each gene + 4 versions of the complex
  ## + column time, ind, trial
  expect_equal(ncol(simNoComplex), 3*2*2 + 2 + 3) ##RNA and protein form of each of the 2 versions of each gene
  ## + PTM version of each of the 2 versions of the protein of gene 3 + column time, ind, trial

  expect_equal(simNoAllele[,"R1"], sim$Simulation[,"R1GCN1"]+sim$Simulation[,"R1GCN2"])
  expect_equal(simNoPTM[,"P3GCN1"], sim$Simulation[,"P3GCN1"]+sim$Simulation[,"Pm3GCN1"])
  expect_equal(simNoComplex[,"P1GCN1"], rowSums(sim$Simulation[,grep("^CPTM1_P1GCN1", names(sim$Simulation))]) + sim$Simulation[,"P1GCN1"])
})
