context("Julia evaluators")
library(sismonr)

check_julia = function(){
  if (!findJuliaNoError(test = T)) {
    skip("Julia not installed")
  }else{
    ## Test that the correct version of Julia is installed
    julia_bin = XRJulia::findJulia()
    if (.Platform$OS.type == "windows") cmd = paste0('"', julia_bin, '" ', "-E \"VERSION < v\\\"1.0.0\\\"\"") else cmd = paste0(julia_bin, " ", "-E 'VERSION < v\"1.0.0\"'")
    if(base::system(cmd, intern = T, ignore.stdout = F)){
      skip("Julia version < v1.0 installed, require v1.0 or later")
    }
  }
}

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
