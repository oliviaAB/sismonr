context("Creating in silico populations")
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

test_that("creation of variants works",{
  check_julia()
  indargs = insilicoindividualargs()
  genes = createGenes(insilicosystemargs(G = 5, PC.p = 0))
  variants = createVariants(genes, indargs)

  expect_equal(length(variants), nrow(genes))
  expect_equal(ncol(variants[[1]]), indargs$ngenevariants)
  expect_true(all(variants[[1]][c("qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregrate", "qtlPTMregrate"),] == 0))
})

test_that("creation of in silico individuals works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 3, ploidy = 4)
  indargs = insilicoindividualargs(ngenevariants = 2)
  genvariants = createVariants(mysystem$genes, indargs)
  genvariants.freq = list(`1` = c(0, 1),
                          `2` = rep(0.5, 2),
                          `3` = rep(0.5, 2))
  myind = createIndividual(mysystem, genvariants, genvariants.freq, indargs)
  qtlvalues = sapply(rownames(genvariants[[1]]), function(x){myind$QTLeffects$GCN1[[x]][1]})

  expect_equal(dim(myind$haplotype), c(3, 4))
  expect_false(any(myind$haplotype[1, ] == 1))
  expect_equal(qtlvalues, genvariants[[1]][, 2])
  expect_is(myind, "insilicoindividual")
})

test_that("creation of in silico individuals works - Initial abundance", {
  check_julia()
  mysystem = createInSilicoSystem(G = 3, ploidy = 2)
  indargs = insilicoindividualargs(ngenevariants = 1)
  genvariants = createVariants(mysystem$genes, indargs)
  genvariants.freq = list(`1` = c(1),
                          `2` = rep(1),
                          `3` = rep(1))
  myind1 = createIndividual(mysystem, genvariants, genvariants.freq, indargs, initialNoise = F)
  myind2 = createIndividual(mysystem, genvariants, genvariants.freq, indargs, initialNoise = F)

  myInitVar = list("R" = c(1, 1, 0), "P" = c(1, 2, 1))
  myind3 = createIndividual(mysystem, genvariants, genvariants.freq, indargs, initialNoise = F, InitVar = myInitVar)

  myind4 = createIndividual(mysystem, genvariants, genvariants.freq, indargs, initialNoise = T)

  expect_equal(myind1$InitAbundance$GCN1$R, myind1$InitAbundance$GCN2$R)
  expect_equal(myind1$InitAbundance$GCN1$P, myind1$InitAbundance$GCN2$P)

  expect_equal(myind1$InitAbundance$GCN1$R, myind2$InitAbundance$GCN1$R)
  expect_equal(myind1$InitAbundance$GCN1$P, myind2$InitAbundance$GCN1$P)

  expect_equal(myind3$InitAbundance$GCN1$R[3], 0)
  expect_equal(myind3$InitAbundance$GCN2$R[3], 0)

  expect_equal(abs(myind3$InitAbundance$GCN1$P[2] - 2 * myind1$InitAbundance$GCN1$P[2]) <= 1, TRUE)
  expect_equal(abs(myind3$InitAbundance$GCN2$P[2] - 2 * myind1$InitAbundance$GCN2$P[2]) <= 1, TRUE)

  ## expect_false(all(myind4$InitAbundance$GCN1$R == myind4$InitAbundance$GCN2$R)) removed just in case
  ##expect_false(all(myind4$InitAbundance$GCN1$P == myind4$InitAbundance$GCN2$P))

  ##expect_false(all(myind4$InitAbundance$GCN1$R == myind1$InitAbundance$GCN1$R))
  ##expect_false(all(myind4$InitAbundance$GCN1$P == myind1$InitAbundance$GCN1$P))
})

test_that("creation of in silico population works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 5, ploidy = 4)
  mypop = createInSilicoPopulation(nInd = 3, mysystem, ngenevariants = 2)

  expect_equal(length(mypop$individualsList), 3)
  expect_equal(length(mypop$GenesVariants), 5)
  expect_equal(ncol(mypop$GenesVariants[[1]]), 2)
  expect_equal(dim(mypop$individualsList$Ind1$haplotype), c(5, 4))
})

#removeJuliaEvaluator(getJuliaEvaluator())
