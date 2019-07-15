context("Creating in silico populations")
library(sismonr)

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
  indargs = insilicoindividualargs()
  genes = createGenes(insilicosystemargs(G = 5, PC.p = 0))
  variants = createVariants(genes, indargs)

  expect_equal(length(variants), nrow(genes))
  expect_equal(ncol(variants[[1]]), indargs$ngenevariants)
  expect_true(all(variants[[1]][c("qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregrate", "qtlPTMregrate"),] == 0))
})

test_that("creation of in silico individuals works", {
  indargs = insilicoindividualargs(ploidy = 4, ngenevariants = 2)
  genes = createGenes(insilicosystemargs(G = 3))
  genvariants = createVariants(genes, indargs)
  genvariants.freq = list(`1` = c(0, 1),
                          `2` = rep(0.5, 2),
                          `3` = rep(0.5, 2))
  myind = createIndividual(genvariants, genvariants.freq, indargs)
  qtlvalues = sapply(rownames(genvariants[[1]]), function(x){myind$QTLeffects$GCN1[[x]][1]})

  expect_equal(dim(myind$haplotype), c(3, 4))
  expect_false(any(myind$haplotype[1, ] == 1))
  expect_equal(qtlvalues, genvariants[[1]][, 2])
  expect_is(myind, "insilicoindividual")
})

test_that("creation of in silico population works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 5)
  mypop = createInSilicoPopulation(nInd = 3, mysystem, ploidy = 4, ngenevariants = 2, sameInit = T)

  expect_equal(length(mypop$individualsList), 3)
  expect_equal(length(mypop$GenesVariants), 5)
  expect_equal(ncol(mypop$GenesVariants[[1]]), 2)
  expect_equal(dim(mypop$individualsList$Ind1$haplotype), c(5, 4))
  expect_true(all(unlist(mypop$individualsList$Ind1$InitVar) == 1))
})
