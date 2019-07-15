context("Creating in silico systems")
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

test_that("creation of gene data-frame works with different parameters", {
  genes = createGenes(insilicosystemargs(G = 5, PC.p = 1, PC.TC.p = 1, basal_transcription_rate_samplingfct = function(x){rep(1, x)}))

  expect_equal(nrow(genes), 5)
  expect_equal(genes$coding, rep("PC", 5))
  expect_equal(genes$TargetReaction, rep("TC", 5))
  expect_equal(genes$TCrate, rep(1, 5))
})

test_that("sampling of gene kinetic params works for noncoding genes",{
  genes = createGenes(insilicosystemargs(G = 10, PC.p = 0.5))

  expect_true(all(genes$ActiveForm[genes$coding == "PC"] == "P"))
  expect_true(all(genes$ActiveForm[genes$coding == "NC"] == "R"))
  expect_true(all(genes$TLrate[genes$coding == "NC"] == 0))
  expect_true(all(genes$PDrate[genes$coding == "NC"] == 0))
})

test_that("creation of one regulatory network works",{
  check_julia()
  nw = createRegulatoryNetwork(regsList = list("PC" = c(1:2), "NC" = c(3)),
                               tarsList = list("PC" = c(4:6), "NC" = integer(0)), reaction = "TL",
                               sysargs = insilicosystemargs(G = 6, regcomplexes = "none"))

  expect_equal(setdiff(c(1, 2, 3), nw$edg$from), 3)
})

test_that("creation of a multi-omic network works",{
  check_julia()
  mysysargs = insilicosystemargs(G = 5)
  mygenes = createGenes(mysysargs)
  mynetwork = createMultiOmicNetwork(mygenes, mysysargs)

  expect_moedg_equal = function(reaction){
    expect_equal(dplyr::filter(mynetwork$edg, TargetReaction == reaction), mynetwork$mosystem[[paste0(reaction, "RN_edg")]][, 1:5])
  }

  expect_is(mynetwork$edg$from, "character")
  expect_is(mynetwork$edg$to, "character")
  expect_moedg_equal("TC")
  expect_moedg_equal("TL")
  expect_moedg_equal("RD")
  expect_moedg_equal("PD")
  expect_moedg_equal("PTM")
})


test_that("creation of an empty multi-omic network works",{
  check_julia()
  mysysargs = insilicosystemargs(G = 5)
  mygenes = createGenes(mysysargs)
  mynetwork = createEmptyMultiOmicNetwork(mygenes)

  expect_moedg_empty = function(reaction){
    expect_equal(nrow(mynetwork$mosystem[[paste0(reaction, "RN_edg")]]), 0)
  }

  expect_equal(nrow(mynetwork$edg), 0)
  expect_moedg_empty("TC")
  expect_moedg_empty("TL")
  expect_moedg_empty("RD")
  expect_moedg_empty("PD")
  expect_moedg_empty("PTM")
})

test_that("creation of an in silico system works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 5, empty = T)

  expect_equal(nrow(mysystem$genes), 5)
  expect_equal(nrow(mysystem$edg), 0)
  expect_is(mysystem, "insilicosystem")
})

test_that("adding a gene works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 5)
  mysystem2 = addGene(mysystem, "PC", "TC", TCrate = 0.0001, TLrate = 0.001)

  expect_equal(nrow(mysystem2$genes), 6)
  expect_equal(mysystem2$genes[6, "coding"], "PC")
  expect_equal(mysystem2$genes[6, "TargetReaction"], "TC")
  expect_equal(mysystem2$genes[6, "TCrate"], 0.0001)
  expect_equal(mysystem2$genes[6, "TLrate"], 0.001)
  expect_is(mysystem2, "insilicosystem")
})

test_that("adding a regulatory complex works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 10, PC.p = 1, PC.TC.p = 1)
  mysystem2 = addComplex(mysystem, c(1, 2, 3), formationrate = 0.1, dissociationrate = 0.01)
  co = rev(names(mysystem2$complexes))[1]

  expect_equal(mysystem2$complexes[[co]], c("1", "2", "3"))
  expect_equal(mysystem2$complexeskinetics[[co]]$formationrate, 0.1)
  expect_equal(mysystem2$complexeskinetics[[co]]$dissociationrate, 0.01)
  expect_equal(mysystem2$complexesTargetReaction[[co]], "TC")
  expect_is(mysystem2, "insilicosystem")
  expect_error(addComplex(mysystem, c(1, 11)), "The components of the complex do not exist in the system.")
})

test_that("removing a complex works", {
  check_julia()
  cond = T
  ## We want to make sure that there is a regulatory complex in the system
  while(cond){
    mysystem = createInSilicoSystem(G = 10, PC.p = 1, PC.TC.p = 1, regcomplexes.p = 1)
    cond = length(mysystem$complexes) == 0
  }
  mysystem2 = removeComplex(mysystem, "CTC1")

  expect_false("CTC1" %in% names(mysystem2$complexes))
  expect_false("CTC1" %in% names(mysystem2$complexeskinetics))
  expect_false("CTC1" %in% names(mysystem2$complexesTargetReaction))
  expect_equal(nrow(mysystem2$edg), nrow(mysystem$edg) - sum(mysystem$edg$from == "CTC1"))
  expect_equal(nrow(mysystem2$mosystem$TCRN_edg), nrow(mysystem$mosystem$TCRN_edg) - sum(mysystem$mosystem$TCRN_edg$from == "CTC1"))
  expect_equal(sum(mysystem2$edg$from == "CTC1"), 0)
  expect_equal(sum(mysystem2$mosystem$TCRN_edg$from == "CTC1"), 0)
  expect_is(mysystem2, "insilicosystem")
  expect_error(removeComplex(mysystem, "CTL1"), "Complex CTL1 does not exist in the system.")
})

test_that("adding an edge works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 5, PC.p = 1, PC.TC.p = 1, empty = TRUE)
  mysystem2 = addEdge(mysystem, 1, 2, regsign = "1", kinetics = c("TCbindingrate"= 0.01, "TCunbindingrate" = 0.1, "TCfoldchange" = 10))

  expect_equal(nrow(mysystem2$edg), 1)
  expect_equal(mysystem2$edg[1, ], data.frame(from = "1", to = "2", TargetReaction = "TC", RegSign = "1", RegBy = "PC", stringsAsFactors = F))
  expect_equal(mysystem2$mosystem$TCRN_edg[1, ], data.frame(from = "1", to = "2", TargetReaction = "TC", RegSign = "1", RegBy = "PC", TCbindingrate = 0.01, TCunbindingrate = 0.1, TCfoldchange = 10, stringsAsFactors = F))
  expect_error(addEdge(mysystem, 1, 10), "Target 10 does not exist in the system.")
  expect_error(addEdge(mysystem, 10, 1), "Regulator 10 does not exist in the system.")
  expect_error(addEdge(mysystem, 1, "CTC1"), "Complexes cannot be regulated. Please provide a gene ID as target.")
})

test_that("removing an edge works", {
  check_julia()
  mysystem = createInSilicoSystem(G = 10)
  edgtoremove = mysystem$edg[1, c("from", "to", "TargetReaction")]
  mysystem2 = removeEdge(mysystem, edgtoremove$from, edgtoremove$to)

  expect_equal(nrow(dplyr::filter(mysystem2$edg, from == edgtoremove$from & to == edgtoremove$to)), 0)
  expect_equal(nrow(dplyr::filter(mysystem2$mosystem[[paste0(edgtoremove$TargetReaction, "RN_edg")]], from == edgtoremove$from & to == edgtoremove$to)), 0)
  expect_message(removeEdge(mysystem, 11, 1), "No edge exists from gene 11 to gene 1.")
})

# removeJuliaEvaluator(getJuliaEvaluator())
