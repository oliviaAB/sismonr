context("Creating system and population arguments lists")
library(sismonr)

test_that("creation of insilicosystemargs works", {
  sysargs = insilicosystemargs(G = 3,
                               PC.p = 0.4,
                               PC.TC.p = 0.6,
                               PC.TL.p = 0,
                               NC.RD.p = 0.8,
                               NC.PD.p = 0.8,
                               basal_translation_rate_samplingfct = function(x){runif(x, 0.1, 0.2)})

  expect_equal(sysargs$G, 3)
  expect_equal(c(sysargs$PC.p, sysargs$NC.p), c(0.4, 0.6))
  expect_equal(c(sysargs$PC.TC.p, sysargs$PC.TL.p, sysargs$PC.RD.p, sysargs$PC.PD.p, sysargs$PC.PTM.p, sysargs$PC.MR.p), c(0.6, 0, 0.1, 0.1, 0.1, 0.1))
  expect_equal(c(sysargs$NC.TC.p, sysargs$NC.TL.p, sysargs$NC.RD.p, sysargs$NC.PD.p, sysargs$NC.PTM.p), c(0, 0, 0.5, 0.5, 0))
  expect_equal(sysargs$RD.pos.p, 1)
  expect_equal(sysargs$PD.pos.p, 1)
  expect_true(all(sysargs[["basal_translation_rate_samplingfct"]](100) <= 0.2 & sysargs[["basal_translation_rate_samplingfct"]](100) >= 0.1))
  expect_equal(length(sysargs[["basal_translation_rate_samplingfct"]](10)), 10)
  expect_is(sysargs, "insilicosystemargs")
})

test_that("creation of insilicoindividualargs works", {
  indargs = insilicoindividualargs(ploidy = 4)

  expect_equal(length(indargs$gcnList), indargs$ploidy)
  expect_equal(length(indargs[["qtleffect_samplingfct"]](10)), 10)
  expect_is(indargs, "insilicoindividualargs")
})
