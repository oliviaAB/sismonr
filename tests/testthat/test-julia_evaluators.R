context("Julia evaluators")
library(sismonr)

test_that("creation of a Julia evaluator works",{
  check_julia()
  ev = newJuliaEvaluator()

  expect_true(XRJulia::juliaEval("@isdefined juliatest"))
})

test_that("using existing Julia evaluators works", {
  check_julia()
  ev = RJulia(.makeNew = T)
  expect_false(XRJulia::juliaEval("@isdefined juliatest", evaluator = ev))
  ev2 = getJuliaEvaluator()
  expect_identical(ev, ev2)
  expect_true(XRJulia::juliaEval("@isdefined juliatest", evaluator = ev))
})

test_that("closing a Julia evaluator works", {
  check_julia()
  ev = getJuliaEvaluator()
  conn = showConnections()
  expect_match(paste0(conn[,"description"], collapse = "; "), paste0(ev$port))
  removeJuliaEvaluator(ev)
  conn = showConnections()
  expect_false(grepl(paste0(ev$port), paste0(conn[,"description"], collapse = "; ")))
})
