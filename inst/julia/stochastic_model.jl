# ------------------------------------------------------------------------------------------------ #
##                    FUNCTIONS FOR CREATING THE TEMPLATE STOCHASTIC SYSTEM                       ##
# ------------------------------------------------------------------------------------------------ #


## Gives all the possible combinations of length l of the values of vector vect
function allposscomb(vect, l)

  nb = length(vect)
  comb = Array{eltype(vect)}(nb^l, l)

  for i in 1:l
    comb[:,i] = repeat(vect, inner = [nb^(l-i)], outer = [nb^(i-1)])
  end

  return comb
end

## Gives all possible combinations of the possible states of the different binding sites of a promoter/RNA
## proms = array where each element corresponds to one binding site, and is an array with all possible states of the binding site (the 1st one being the free state)
function combinallpromstates(proms)
  elsize = [length(t) for t in proms]
  combproms = Array{String}(prod(elsize), length(elsize))

  for i in eachindex(elsize)
    combproms[:, i] = repeat(proms[i], inner = prod(elsize[(i+1):end]), outer = prod(elsize[1:(i-1)]))
  end

  return combproms
end



## Gives all the possible combinations of the possible active states of the different binding sites of a promoter/RNA
## promList: an 1D array, each element being a 2D array. Each element corresponds to the possible active states of a 
## given promoter binding site (1st column) and the associated fold change (2nd column)
function combinactivepromstates(promList)
  proms = [t[:, 1] for t in promList]
  fcs = [t[:, 2] for t in promList]
  elsize = [length(t) for t in proms]
  combproms = Array{String}(prod(elsize), length(elsize))
  combfcs = Array{Float64}(prod(elsize), length(elsize))

  for i in eachindex(elsize)
    combproms[:, i] = repeat(proms[i], inner = prod(elsize[(i+1):end]), outer = prod(elsize[1:(i-1)]))
    combfcs[:, i] = repeat(fcs[i], inner = prod(elsize[(i+1):end]), outer = prod(elsize[1:(i-1)]))
  end

  return Dict("proms" => combproms, "fcs" => combfcs)
end



## Creates a chemical reaction for the BioSimulator module 
function reactBioSim(reacList, prodList)
  return join(reacList, " + ") * " --> " * join(prodList, " + ")
end


## Creates the regulatory complexes binding and unbding reactions
function createComplexesReactions(complexes, complexeskinetics, complexsize, complexvariants, activeform)

  spec = []
  initcond = []
  reac = []
  reacnames = []
  prop = []

  for compl in keys(complexes), t in 1:size(complexvariants)[1]
    complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
    push!(spec, complvar) ## adding the complex form to the list of species
    push!(initcond, "0") ## adding the initial abundance of the complex form to the list of initial conditions. For complexes, initial abundance set to 0
    ## Creates the reaction of complex formation
    push!(reac, reactBioSim([activeform[string(complexes[compl][i])]*complexvariants[t, i] for i in 1:complexsize], [complvar])) ## sum of complex components -> compl
    push!(reacnames, "formation"*complvar) 
    push!(prop, """$(complexeskinetics[compl]["formationrate"])""")
    push!(reac, reactBioSim([complvar], [activeform[string(complexes[compl][i])]*complexvariants[t, i] for i in 1:complexsize])) ## compl -> sum of complex components
    push!(reacnames, "dissociation"*complvar) 
    push!(prop, """$(complexeskinetics[compl]["dissociationrate"])""")
  end

  return Dict("species" => spec, "initialconditions" => initcond, "reactions" => reac, "reactionsnames" => reacnames, "propensities" => prop)
end



## Creates the transcription regulators binding and unbinding reactions + promoter binding sites of the genes
function createTCregReactions(edg, genes, activeform, complexes, complexsize, complexvariants, gcnList)
  spec = []
  initcond = []
  reac = []
  reacnames = []
  prop = []

  regsingl = find( x -> !in('C', x), edg["from"]) ## identify single-molecule regulators (i.e. not regulatory complexes)
  regcompl = find( x -> in('C', x), edg["from"]) ## identify regulatory complexes

  promActiveStates = Dict(string(i)*j => [] for i in genes["id"], j in gcnList) ## dictionary, 1 element for each allele version of each gene. Value is an array with m elements, 1 for each regulator. 
                                                                              ## The j-th element of this array is a matrix: rows: all possible active states of the binding site of regulator j (column 1: name of the binding site state, column 2 fold-change associated with this state)
  promAllStates = Dict(string(i)*j => [] for i in genes["id"], j in gcnList) ## dictionary, 1 element for each allele version of each gene. Value is an array with m elements, 1 for each regulator. 
                                                                              ## The j-th element of this array is an array, listing all the possible states of the binding site of regulator j (even if not active)

  ## Loop over all the edges in the transcription regulatory network with single genes as regulators
  for r in regsingl, gcn in gcnList 
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    reg = edg["from"][r]

    ## Create the binding site for regulator
    prom = "Pr" * tar * "reg" * reg

    ## tempactivestates is a matrix where each row correspond to an active state of the promoter: 
    ##    1st column = name of species corresponding to the promoter state, 2nd column = fold change associated with the promoter state
    tempactivestates = [prom*"F" 1]
    boundactive = edg["RegSign"][r] == "1" ## are the bound states of the promoter actives i.e. able to transcribe? No if the regulator is a repressor
    tempallstates = [prom*"F"]

    ## Add the free form of the promoter to the list of species (the bound form is specific to the bound molecule)
    push!(spec, prom*"F")
    ## Add the initial abundance of the free form of the promoter to the list of initial conditions. (1 promoter binding site per gene copy)
    push!(initcond, "1")

    for gcnreg in gcnList ## loop over all possible gene copies of the regulator
      prombound = prom*gcnreg*"B"
      push!(spec, prombound) ## add the bound promoter to the list of species
      push!(initcond, "0") ## add its initial abundance to the list of initial conditions. Initial abundance of bound promoter set to 0
      
      if boundactive
        tempactivestates = vcat(tempactivestates, [prombound edg["TCfoldchange"][r]]) ## if the regulator is an activator add this bound state of the promoter site to the list of active states + the induced fold change
      end
      push!(tempallstates, prombound) ## add this bound state of the promoter site to the list of all sites
      
      ## Add the binding reaction of the regulator to the binding site
      push!(reac, reactBioSim([prom*"F", activeform[reg]*gcnreg], [prombound])) ## promF + reg -> promregB
      push!(reacnames, "binding"*prom*gcnreg)
      push!(prop, """$(edg["TCbindingrate"][r])*QTLeffects["$(gcn)"]["$("qtlTCregbind")"][$(tarid)]*QTLeffects["$(gcnreg)"]["qtlactivity"][$(parse(Int64, reg))]""")
      
      ## Add the unbinding reaction of the regulator from the binding site
      push!(reac, reactBioSim([prombound], [prom*"F", activeform[reg]*gcnreg])) ## promregB -> promF + reg
      push!(reacnames, "unbinding"*prom*gcnreg)
      push!(prop, """$(edg["TCunbindingrate"][r])""")
    end

    push!(promActiveStates[tar], tempactivestates)
    push!(promAllStates[tar], tempallstates)
  end


  ## Loop over all the edges in the transcription regulatory network with complexes as regulators
  for r in regcompl, gcn in gcnList 
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    compl = edg["from"][r]

    ## Create the binding site for regulator
    prom = "Pr" * tar* "reg" * compl

    ## tempactivestates is a matrix where each lign correspond to an active state of the promoter: 
    ##    1st column = name of species corresponding to the promoter state, 2nd column = fold change associated with the promoter state
    tempactivestates = [prom*"F" 1]
    boundactive = edg["RegSign"][r] == "1" ## are the bound states of the promoter actives i.e. able to transcribe? No if the regulator is a repressor
    tempallstates = [prom*"F"]

    ## Add the free and bound forms of the promoter to the list of species
    push!(spec, prom*"F")
    ## Add the initial abundance of the free form of the promoter to the list of initial conditions. (1 promoter binding site per gene copy)
    push!(initcond, "1")

    for t in 1:size(complexvariants)[1]
      complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
        
      prombound = prom*join(complexvariants[t, :])*"B"
      push!(spec, prombound)
      push!(initcond, "0") ## add its initial abundance to the list of initial conditions. Initial abundance of bound promoter set to 0

      if boundactive
        tempactivestates = vcat(tempactivestates, [prombound edg["TCfoldchange"][r]])
      end
      push!(tempallstates, prombound)

      ## Add the binding reaction of the regulator to the binding site
      push!(reac, reactBioSim([prom*"F", complvar], [prombound])) ## promF + complvar -> promB
      push!(reacnames, "binding"*prom*complvar)
      tempprop = """$(edg["TCbindingrate"][r])*QTLeffects["$(gcn)"]["$("qtlTCregbind")"][$(tarid)]"""*join(["""*QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
      push!(prop, tempprop)

      ## add the unbinding of the regulator from the binding site
      push!(reac, reactBioSim([prombound], [prom*"F", complvar])) ## promB -> promF + complvar
      push!(reacnames, "unbinding"*prom*complvar)
      push!(prop, """$(edg["TCunbindingrate"][r])""")
    end
    
    push!(promActiveStates[tar], tempactivestates)
    push!(promAllStates[tar], tempallstates)
  end

  return Dict("species" => spec, "initialconditions" => initcond, "reactions" => reac, "reactionsnames" => reacnames, "propensities" => prop, "promActiveStates" => promActiveStates, "promAllStates" => promAllStates)
end



## Creates the translation regulators binding and unbinding reactions + promoter binding sites of the RNAs
function createTLregReactions(edg, genes, activeform, complexes, complexsize, complexvariants, gcnList)
  spec = []
  initcond = []
  reac = []
  reacnames = []
  prop = []

  regsingl = find( x -> !in('C', x), edg["from"]) ## identify single-molecule regulators (i.e. not regulatory complexes)
  regcompl = find( x -> in('C', x), edg["from"]) ## identify regulatory complexes

  promActiveStates = Dict(string(i)*j => [] for i in genes["id"], j in gcnList) ## dictionary, 1 element for each allele version of each gene. Value is an array with m elements, 1 for each regulator. 
                                                                              ## The j-th element of this array is a matrix: rows: all possible active states of the binding site of regulator j (column 1: name of the binding site state, column 2 fold-change associated with this state)
  promAllStates = Dict(string(i)*j => [] for i in genes["id"], j in gcnList) ## dictionary, 1 element for each allele version of each gene. Value is an array with m elements, 1 for each regulator. 
                                                                              ## The j-th element of this array is an array, listing all the possible states of the binding site of regulator j (even if not active)

  ## Loop over all the edges in the translation regulatory network with single genes as regulators
  for r in regsingl, gcn in gcnList 
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    reg = edg["from"][r]

    ## Create the binding site for regulator
    prom = "RBS" * tar * "reg" * reg

    ## tempactivestates is a matrix where each row correspond to an active state of the binding site: 
    ##    1st column = name of species corresponding to the binding site state, 2nd column = fold change associated with the binding site state
    tempactivestates = [prom*"F" 1]
    boundactive = edg["RegSign"][r] == "1" ## are the bound states of the binding site actives i.e. able to translate? No if the regulator is a repressor
    tempallstates = [prom*"F"]

    ## Add the free form of the binding site to the list of species (the bound form is specific to the bound molecule)
    push!(spec, prom*"F")
    ## Add the initial abundance of the free form of the binding site to the list of initial conditions.
    push!(initcond, """$(genes["TCrate"][tarid]/genes["RDrate"][tarid])*InitVar["$(gcn)"]["R"][$(tarid)]""")

    for gcnreg in gcnList ## loop over all possible gene copies of the regulator
      prombound = prom*gcnreg*"B"
      push!(spec, prombound) ## add the bound binding site to the list of species
      push!(initcond, "0") ## add its initial abundance to the list of initial conditions. Initial abundance of bound binding site set to 0
      
      if boundactive
        tempactivestates = vcat(tempactivestates, [prombound edg["TLfoldchange"][r]]) ## if the regulator is an activator add this bound state of the binding site to the list of active states + the induced fold change
      end
      push!(tempallstates, prombound) ## add this bound state of the binding site to the list of all sites
      
      ## Add the binding reaction of the regulator to the binding site
      push!(reac, reactBioSim([prom*"F", activeform[reg]*gcnreg], [prombound])) ## promF + reg -> promregB
      push!(reacnames, "binding"*prom*gcnreg)
      push!(prop, """$(edg["TLbindingrate"][r])*QTLeffects["$(gcn)"]["$("qtlTLregbind")"][$(tarid)]*QTLeffects["$(gcnreg)"]["qtlactivity"][$(parse(Int64, reg))]""")
      
      ## Add the unbinding reaction of the regulator from the binding site
      push!(reac, reactBioSim([prombound], [prom*"F", activeform[reg]*gcnreg])) ## promregB -> promF + reg
      push!(reacnames, "unbinding"*prom*gcnreg)
      push!(prop, """$(edg["TLunbindingrate"][r])""")
    end

    push!(promActiveStates[tar], tempactivestates)
    push!(promAllStates[tar], tempallstates)
  end


  ## Loop over all the edges in the translation regulatory network with complexes as regulators
  for r in regcompl, gcn in gcnList 
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    compl = edg["from"][r]

    ## Create the binding site for regulator
    prom = "RBS" * tar* "reg" * compl

    ## tempactivestates is a matrix where each lign correspond to an active state of the binding site: 
    ##    1st column = name of species corresponding to the binding site state, 2nd column = fold change associated with the binding site state
    tempactivestates = [prom*"F" 1]
    boundactive = edg["RegSign"][r] == "1" ## are the bound states of the binding site actives i.e. able to transcribe? No if the regulator is a repressor
    tempallstates = [prom*"F"]

    ## Add the free and bound forms of the binding site to the list of species
    push!(spec, prom*"F")
    ## Add the initial abundance of the free form of the binding site to the list of initial conditions.
    push!(initcond, """$(genes["TCrate"][tarid]/genes["RDrate"][tarid])*InitVar["$(gcn)"]["R"][$(tarid)]""")

    for t in 1:size(complexvariants)[1]
      complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
        
      prombound = prom*join(complexvariants[t, :])*"B"
      push!(spec, prombound)
      push!(initcond, "0") ## add its initial abundance to the list of initial conditions. Initial abundance of bound binding site set to 0

      if boundactive
        tempactivestates = vcat(tempactivestates, [prombound edg["TLfoldchange"][r]])
      end
      push!(tempallstates, prombound)

      ## Add the binding reaction of the regulator to the binding site
      push!(reac, reactBioSim([prom*"F", complvar], [prombound])) ## promF + complvar -> promB
      push!(reacnames, "binding"*prom*complvar)
      tempprop = """$(edg["TLbindingrate"][r])*QTLeffects["$(gcn)"]["$("qtlTLregbind")"][$(tarid)]"""*join(["""*QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
      push!(prop, tempprop)

      ## add the unbinding of the regulator from the binding site
      push!(reac, reactBioSim([prombound], [prom*"F", complvar])) ## promB -> promF + complvar
      push!(reacnames, "unbinding"*prom*complvar)
      push!(prop, """$(edg["TLunbindingrate"][r])""")
    end
    
    push!(promActiveStates[tar], tempactivestates)
    push!(promAllStates[tar], tempallstates)
  end

  return Dict("species" => spec, "initialconditions" => initcond, "reactions" => reac, "reactionsnames" => reacnames, "propensities" => prop, "promActiveStates" => promActiveStates, "promAllStates" => promAllStates)
end


## Creates the transcription reactions of the genes
function createTranscriptionReactions(genes, promactiveTC, promallTL, gcnList)
  spec = []
  initcond = []
  reac = []
  reacnames = []
  prop = []

  rnaforms = Dict()

  for g in genes["id"], gcn in gcnList
    gname = string(g) * gcn
    TCactiveproms = promactiveTC[gname]
    TLallproms = promallTL[gname]  

    ## What is the RNA form of the gene? 
    if length(TLallproms) == 0 ## R[gene] if the gene is not controlled at the translation level
      RNAform = ["R"*gname] 
      push!(spec, "R"*gname)
      push!(initcond, """$(genes["TCrate"][g]/genes["RDrate"][g])*InitVar["$(gcn)"]["R"][$(g)]""")
    else ## otherwise transcription produces the free (unbound) form of each binding site present on the RNA (1 site per TL regulator)
      RNAform = [t[1] for t in TLallproms] ## The TLallproms dictionary is constructed such that each free binding site name is at the position 1
    end

    if length(TCactiveproms) == 0
      ## Create transcription of the gene 
      push!(reac, reactBioSim([0], RNAform))
      push!(reacnames, "transcription"*gname)
      push!(prop, """$(genes["TCrate"][g])*QTLeffects["$(gcn)"]["qtlTCrate"][$(g)]""")
    else
      promcomb = combinactivepromstates(TCactiveproms) ## gives all possible active combinations of the different binding site states
      for t in 1:size(promcomb["proms"])[1]
        push!(reac, reactBioSim(promcomb["proms"][t,:], vcat(promcomb["proms"][t,:], RNAform)))
        push!(reacnames, "transcription"*join(promcomb["proms"][t,:]))
        push!(prop, """$(genes["TCrate"][g])*QTLeffects["$(gcn)"]["qtlTCrate"][$(g)]*$(prod(promcomb["fcs"][t,:]))""")
      end
    end
    rnaforms[gname] = RNAform
  end

  return Dict("species" => spec, "initialconditions" => initcond, "reactions" => reac, "reactionsnames" => reacnames, "propensities" => prop, "RNAforms" => rnaforms)
end


## Creates the RNA decay reactions of the genes
function createRNADecayReactions(genes, RNAforms, gcnList)
  reac = []
  reacnames = []
  prop = []

  for g in genes["id"], gcn in gcnList
    gname = string(g) * gcn
    rnaform = RNAforms[gname]
    push!(reac, reactBioSim(rnaform, [0]))
    push!(reacnames, "RNAdecay"*join(rnaform))
    push!(prop, """$(genes["RDrate"][g])*QTLeffects["$(gcn)"]["qtlRDrate"][$(g)]""")
  end

  return Dict("reactions" => reac, "reactionsnames" => reacnames, "propensities" => prop)
end


## Creates the regulator-mediated RNA decay reactions of the genes
function createRDregReactions(edg, genes, RNAforms, activeform, complexes, complexsize, complexvariants, gcnList)
  reac = []
  reacnames = []
  prop = []

  regsingl = find( x -> !in('C', x), edg["from"]) ## identify single-molecule regulators (i.e. not regulatory complexes)
  regcompl = find( x -> in('C', x), edg["from"]) ## identify regulatory complexes

  for r in regsingl, gcn in gcnList 
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    tarRNA = RNAforms[tar]
    reg = edg["from"][r]

    for gcnreg in gcnList
      push!(reac, reactBioSim(vcat(tarRNA, activeform[reg]*gcnreg), [activeform[reg]*gcnreg]))
      push!(reacnames, "RNAdecay"*join(tarRNA)*"reg"*reg*gcnreg)
      push!(prop, """$(edg["RDregrate"][r])*QTLeffects["$(gcn)"]["$("qtlRDregrate")"][$(tarid)]*QTLeffects["$(gcnreg)"]["qtlactivity"][$(parse(Int64, reg))]""")
    end
  end

  for r in regcompl, gcn in gcnList
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    tarRNA = RNAforms[tar]
    compl = edg["from"][r]

    for t in 1:size(complexvariants)[1]
      complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
      push!(reac, reactBioSim(vcat(tarRNA, complvar), [complvar]))
      push!(reacnames, "RNAdecay"*join(tarRNA)*"reg"*complvar)
      tempprop = """$(edg["RDregrate"][r])*QTLeffects["$(gcn)"]["$("qtlRDregrate")"][$(tarid)]"""*join(["*"*"""QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
      push!(prop, tempprop)
    end
  end

  return Dict("reactions" => reac, "reactionsnames" => reacnames, "propensities" => prop)
end


## Creates the translation reactions of the genes
function createTranslationReactions(genes, promactiveTL, gcnList)
  spec = []
  initcond = []
  reac = []
  reacnames = []
  prop = []

  for g in genes["id"][genes["coding"] .== "PC"], gcn in gcnList
    gname = string(g) * gcn
    TLactiveproms = promactiveTL[gname]

    push!(spec, "P"*gname)
    push!(initcond, """($(genes["TCrate"][g]*genes["TLrate"][g]/(genes["RDrate"][g]*genes["PDrate"][g])))*InitVar["$(gcn)"]["P"][$(g)]""")

    if length(TLactiveproms) == 0
      ## create translation of the gene
      push!(reac, reactBioSim(["R"*gname], ["R"*gname, "P"*gname]))
      push!(reacnames, "translation"*gname)
      push!(prop, """$(genes["TLrate"][g])*QTLeffects["$(gcn)"]["qtlTLrate"][$(g)]""")
    else
      promcomb = combinactivepromstates(TLactiveproms) ## gives all possible active combinations of the different binding site states
      for t in 1:size(promcomb["proms"])[1]
        push!(reac, reactBioSim(promcomb["proms"][t,:], vcat(promcomb["proms"][t,:], "P"*gname)))
        push!(reacnames, "translation"*join(promcomb["proms"][t,:]))
        push!(prop, """$(genes["TLrate"][g])*QTLeffects["$(gcn)"]["qtlTLrate"][$(g)]*$(prod(promcomb["fcs"][t,:]))""")
      end
    end
  end

  return Dict("species" => spec, "initialconditions" => initcond, "reactions" => reac, "reactionsnames" => reacnames, "propensities" => prop)
end


## Creates the protein decay reactions of the genes
function createProteinDecayReactions(genes, gcnList)
  spec = []
  initcond = []
  reac = []
  reacnames = []
  prop = []

  for g in genes["id"][genes["coding"] .== "PC"], gcn in gcnList
    gname = string(g) * gcn
    push!(reac, reactBioSim(["P"*gname], [0]))
    push!(reacnames, "proteindecay"*"P"*gname)
    push!(prop, """$(genes["PDrate"][g])*QTLeffects["$(gcn)"]["qtlPDrate"][$(g)]""")

    if genes["PTMform"][g] =="1"
      ## Add to the list of species the PTM form of the protein
      push!(spec, "Pm"*gname)
      push!(initcond, "0") ## initial abundance of modified proteins set to 0
      push!(reac, reactBioSim(["Pm"*gname], [0]))
      push!(reacnames, "proteindecay"*"Pm"*gname)
      push!(prop, """$(genes["PDrate"][g])*QTLeffects["$(gcn)"]["qtlPDrate"][$(g)]""")
    end
  end

  return Dict("species" => spec, "initialconditions" => initcond, "reactions" => reac, "reactionsnames" => reacnames, "propensities" => prop)
end


## Creates the regulator-mediated protein decay reactions of the genes
function createPDregReactions(edg, genes, activeform, complexes, complexsize, complexvariants, gcnList)
  reac = []
  reacnames = []
  prop = []

  regsingl = find( x -> !in('C', x), edg["from"]) ## identify single-molecule regulators (i.e. not regulatory complexes)
  regcompl = find( x -> in('C', x), edg["from"]) ## identify regulatory complexes

  for r in regsingl, gcn in gcnList 
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    reg = edg["from"][r]

    pform = ["P"*tar]
    if genes["PTMform"][tarid] == "1"
      push!(pform, "Pm"*tar)
    end
    
    for p in pform, gcnreg in gcnList
      push!(reac, reactBioSim([p, activeform[reg]*gcnreg], [activeform[reg]*gcnreg]))
      push!(reacnames, "proteindecay"*p*"reg"*reg*gcnreg)
      push!(prop, """$(edg["PDregrate"][r])*QTLeffects["$(gcn)"]["$("qtlPDregrate")"][$(tarid)]*QTLeffects["$(gcnreg)"]["qtlactivity"][$(parse(Int64, reg))]""")
    end
  end

  for r in regcompl, gcn in gcnList
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    compl = edg["from"][r]

    pform = ["P"*tar]
    if genes["PTMform"][tarid] == "1"
      push!(pform, "Pm"*tar)
    end

    for p in pform, t in 1:size(complexvariants)[1]
      complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
      push!(reac, reactBioSim([p, complvar], [complvar]))
      push!(reacnames, "proteindecay"*p*"reg"*complvar)
      tempprop = """$(edg["PDregrate"][r])*QTLeffects["$(gcn)"]["$("qtlPDregrate")"][$(tarid)]"""*join(["*"*"""QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
      push!(prop, tempprop)
    end
  end

  return Dict("reactions" => reac, "reactionsnames" => reacnames, "propensities" => prop)
end

## Creates the protein post-translational modification reactions of the genes
function createPTMregReactions(edg, genes, activeform, complexes, complexsize, complexvariants, gcnList)
  reac = []
  reacnames = []
  prop = []

  regsingl = find( x -> !in('C', x), edg["from"]) ## identify single-molecule regulators (i.e. not regulatory complexes)
  regcompl = find( x -> in('C', x), edg["from"]) ## identify regulatory complexes

  for r in regsingl, gcn in gcnList 
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    reg = edg["from"][r]
    ispos = edg["RegSign"][r] == "1" ## is the regulator transforming the orginal protein into its modified form (RegSign = "1") or the opposite (RegSign = "-1")
    for gcnreg in gcnList
      if ispos
        push!(reac, reactBioSim([ "P"*tar, activeform[reg]*gcnreg], ["Pm"*tar, activeform[reg]*gcnreg]))
        push!(reacnames, "PTM"*tar*"reg"*reg*gcnreg)
      else
        push!(reac, reactBioSim([ "Pm"*tar, activeform[reg]*gcnreg], ["P"*tar, activeform[reg]*gcnreg]))
        push!(reacnames, "de-PTM"*tar*"reg"*reg*gcnreg)
      end
      push!(prop, """$(edg["PTMregrate"][r])*QTLeffects["$(gcn)"]["$("qtlPTMregrate")"][$(tarid)]*QTLeffects["$(gcnreg)"]["qtlactivity"][$(parse(Int64, reg))]""") 
    end
  end

  for r in regcompl, gcn in gcnList
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    compl = edg["from"][r]
    ispos = edg["RegSign"][r] == "1" ## is the regulator transforming the orginal protein into its modified form (RegSign = "1") or the opposite (RegSign = "-1")

    for t in 1:size(complexvariants)[1]
      complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])

      if ispos
        push!(reac, reactBioSim([ "P"*tar, complvar], ["Pm"*tar, complvar]))
        push!(reacnames, "PTM"*tar*"reg"*complvar)
      else
        push!(reac, reactBioSim([ "Pm"*tar, complvar], ["P"*tar, complvar]))
        push!(reacnames, "de-PTM"*tar*"reg"*complvar)
      end
      tempprop = """$(edg["PTMregrate"][r])*QTLeffects["$(gcn)"]["$("qtlPTMregrate")"][$(tarid)]"""*join(["*"*"""QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
      push!(prop, tempprop)
    end
  end
    
  return Dict("reactions" => reac, "reactionsnames" => reacnames, "propensities" => prop)
end

function juliaCreateStochasticSystem(genes, edgTCRN, edgTLRN, edgRDRN, edgPDRN, edgPTMRN, complexes, complexeskinetics, complexsize, gcnList, writefile, filepath, filename)

  ## Ensures that gcnList is a vector (not true if there is only one allele for each gene)
  if typeof(gcnList) == String
    gcnList = [gcnList]
  end

  ## Active form of the genes
  activeform = Dict(zip(map(string, genes["id"]), genes["ActiveForm"]))

  ## Creates the list of all possible allele combinations of the components of a complex
  complexvariants = allposscomb(gcnList, complexsize)

  ## Generates the formation and dissociation reactions of the different regulatory complexes
  complexesReactions = createComplexesReactions(complexes, complexeskinetics, complexsize, complexvariants, activeform)
  # println("complexes done")

  ## Generates the list of all possible binding and unbinding reactions of transcription regulators on promoter binding sites
  regTCreactions = createTCregReactions(edgTCRN, genes, activeform, complexes, complexsize, complexvariants, gcnList)
  promactiveTC = regTCreactions["promActiveStates"]
  promallTC = regTCreactions["promAllStates"]
  # println("regTCreactions done")

  ## Generates the list of all possible binding and unbinding reactions of transcription regulators on promoter binding sites
  regTLreactions = createTLregReactions(edgTLRN, genes, activeform, complexes, complexsize, complexvariants, gcnList)
  promactiveTL = regTLreactions["promActiveStates"]
  promallTL = regTLreactions["promAllStates"]
  # println("regTLreactions done")

  ## Generates the list of all possible transcription reactions for the genes
  transcriptionReactions = createTranscriptionReactions(genes, promactiveTC, promallTL, gcnList)
  RNAforms = transcriptionReactions["RNAforms"]
  # println("transcriptionReactions done")

  ## Generates the list of all possible RNA decay reactions for the genes
  ## In this model, when RNAs are represented by the sum of RBS (translation regulator binding sites), only the free form of the RNA
  ## (= all binding sites are free) decays
  rnaDecayReactions = createRNADecayReactions(genes, RNAforms, gcnList)
  # println("rnaDecayReactions done")

  ## Generates the list of all possible regulator-mediated RNA decay reactions for the genes
  regRDreactions = createRDregReactions(edgRDRN, genes, RNAforms, activeform, complexes, complexsize, complexvariants, gcnList)
  # println("regRDreactions done")

  ## Generates the list of all possible translation reactions for the genes
  translationReactions = createTranslationReactions(genes, promactiveTL, gcnList)
  # println("translationReactions done")

  ## Generates the list of all possible protein decay reactions for the genes
  proteinDecayReactions = createProteinDecayReactions(genes, gcnList)
  # println("proteinDecayReactions done")

  ## Generates the list of all possible regulator-mediated protein decay reactions for the genes
  regPDreactions = createPDregReactions(edgPDRN, genes, activeform, complexes, complexsize, complexvariants, gcnList)
  # println("regPDreactions done")

  ## Generates the list of all possible protein post-translational modification reactions for the genes
  regPTMreactions = createPTMregReactions(edgPTMRN, genes, activeform, complexes, complexsize, complexvariants, gcnList)
  # println("regPTMreactions done")

  ## Output of the function
  species = vcat(complexesReactions["species"],
                 regTCreactions["species"],
                 regTLreactions["species"],
                 transcriptionReactions["species"],
                 translationReactions["species"],
                 proteinDecayReactions["species"])

  initialconditions = vcat(complexesReactions["initialconditions"],
                           regTCreactions["initialconditions"],
                           regTLreactions["initialconditions"],
                           transcriptionReactions["initialconditions"],
                           translationReactions["initialconditions"],
                           proteinDecayReactions["initialconditions"])

  reactions = vcat(complexesReactions["reactions"],
                   regTCreactions["reactions"],
                   regTLreactions["reactions"],
                   transcriptionReactions["reactions"],
                   rnaDecayReactions["reactions"],
                   regRDreactions["reactions"],
                   translationReactions["reactions"],
                   proteinDecayReactions["reactions"],
                   regPDreactions["reactions"],
                   regPTMreactions["reactions"])

  reactionsnames = vcat(complexesReactions["reactionsnames"],
                        regTCreactions["reactionsnames"],
                        regTLreactions["reactionsnames"],
                        transcriptionReactions["reactionsnames"],
                        rnaDecayReactions["reactionsnames"],
                        regRDreactions["reactionsnames"],
                        translationReactions["reactionsnames"],
                        proteinDecayReactions["reactionsnames"],
                        regPDreactions["reactionsnames"],
                        regPTMreactions["reactionsnames"])

  propensities = vcat(complexesReactions["propensities"],
                      regTCreactions["propensities"],
                      regTLreactions["propensities"],
                      transcriptionReactions["propensities"],
                      rnaDecayReactions["propensities"],
                      regRDreactions["propensities"],
                      translationReactions["propensities"],
                      proteinDecayReactions["propensities"],
                      regPDreactions["propensities"],
                      regPTMreactions["propensities"])

  if writefile

    if !ismatch(r"/$", filepath) ## because R function getwd() gives the path to the directory folder without a "/" at the end
      filepath = filepath*"/"
    end

    open(filepath*filename*"_species.txt", "w") do f
      for i in 1:length(species)
        write(f, species[i]*"\t"*initialconditions[i]*"\n")
      end
    end

    open(filepath*filename*"_reactions.txt", "w") do f
      for i in 1:length(reactions)
        write(f, reactionsnames[i]*"\t"*reactions[i]*"\t"*propensities[i]*"\n")
      end
    end
  end

  return Dict("species" => species, "initialconditions" => initialconditions, "reactions" => reactions, "reactionsnames" => reactionsnames, "propensities" => propensities, "TCproms" => promallTC, "TLproms" => promallTL)
end



  #  save("/home/oangelin/Documents/testsismonr/mysystem.jld", "genes", genes, "edgTCRN", edgTCRN, "edgTLRN", edgTLRN, "edgRDRN", edgRDRN,
  #                                                            "edgPDRN", edgPDRN, "edgPTMRN", edgPTMRN, "complexes", complexes, "complexeskinetics", complexeskinetics,
  #                                                            "complexsize", complexsize, "gcnList", gcnList)

#=
mysystem = load("/home/oangelin/Documents/testsismonr/mysystem.jld")
for i in collect(keys(mysystem))
       s = Symbol(i)
       v = mysystem[i]
       @eval (($s) = ($v))
end
=#