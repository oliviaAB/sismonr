using StatsBase
# using JLD ## temporary
using BioSimulator
using DataFrames

# ------------------------------------------------------------------------------------------------ #
##                       FUNCTIONS FOR SIMULATING THE STOCHASTIC SYSTEM                           ##
# ------------------------------------------------------------------------------------------------ #

function createBioSimModel(stochmodel, QTLeffects, InitVar, modelname)
       
  model = BioSimulator.Network(modelname)

  ## Add the species in the model, with their initial abundance
  for i in 1:length(stochmodel["species"])
    i0 = replace(stochmodel["initialconditions"][i], "QTLeffects" => "$QTLeffects")
    i0 = replace(i0, "InitVar" => "$InitVar")
    i0 = eval(Meta.parse(i0))
    #println(stochmodel["species"][i]* "\t"*string(i0))

    if !isa(i0, Number)
      println(stochmodel["initialconditions"][i])
      error("Pb during the construction of the Julia BioSimulator model: initial abundance is not numeric.")
    end
    if typeof(stochmodel["species"][i]) != String
      println(stochmodel["species"][i])
      error("Pb during the construction of the Julia BioSimulator model: species name is not a String.")
    end

    if round(Int, i0)<0
      error("Initial condition of species $(stochmodel["species"][i]) is negative ($(i0)).")
    end

    model <= BioSimulator.Species(stochmodel["species"][i], round(Int, i0))

  end


  ## Add the reactions in the model, with their name and rate
  for i in eachindex(stochmodel["reactions"])
    prop = replace(stochmodel["propensities"][i], "QTLeffects" => "$QTLeffects")
    #println(prop)
    prop = eval(Meta.parse(prop))

    if !isa(prop, Number)
      println(stochmodel["propensities"][i])
      error("Pb during the construction of the Julia BioSimulator model: propensity is not numeric.")
    end
    if typeof(stochmodel["reactionsnames"][i]) != String
      println(stochmodel["reactionsnames"][i])
      error("Pb during the construction of the Julia BioSimulator model: reaction name is not a String.")
    end
    if typeof(stochmodel["reactions"][i]) != String
      println(stochmodel["reactions"][i])
      error("Pb during the construction of the Julia BioSimulator model: reaction is not a String.")
    end

    if prop<0
      error("Reaction rate of reaction $(stochmodel["reactionsnames"][i]) is negative ($(prop)).")
    end

    model <= BioSimulator.Reaction(stochmodel["reactionsnames"][i], prop, stochmodel["reactions"][i])

  end

  return model
end

function allequal(x)
    all(y->y==x[1], x)
end

## Transforms the results of the simulation into the abundance profile of the different molecules over time
## This is beause in the simulation each regulator binding site is a distinct molecule, but in our system 
## they can belong to a same RNA molecule
function transformSimRes2Abundance(resultdf, genes, stochmodel)
  
  abundancedf = resultdf[:, [:time, :trial]]
  #abundancedf = Dict("time" => resultdf[:, :time], "trial" => resultdf[:, :trial])

  for g in collect(keys(stochmodel["TCproms"]))
    gid = parse(Int64, replace(g, r"GCN.+$" => "")) ## gives the gene id

    ## Check that for each binding site on the promoter of each gene, at each time the sum of the abundance
    ## of all species corresponding to all possible binding site states equals 1 (bc only 1 binding site per gene)
    for i in eachindex(stochmodel["TCproms"][g])
      prabundance = resultdf[:, map(Symbol, stochmodel["TCproms"][g][i])]
      sumprom = [sum(permutedims(Vector(row))) for row in eachrow(prabundance)]
      if any(sumprom .!=1)
        error("The sum of promoter states for "*stochmodel["TCproms"][g][i][1]*"not equal to 1.")
      end
    end

    if length(stochmodel["TLproms"][g]) > 0        
      ## Check for each RNA that the different binding sites on the RNA are in equal abundance at each time of the simulation
      rbsabundance = [sum(permutedims(Vector(resultdf[t, map(Symbol, x)]))) for t in 1:size(resultdf)[1], x in stochmodel["TLproms"][g]]
                  
      if !all(mapslices(allequal, rbsabundance, dims = 2))
        error("The abundance of the different binding sites on the RNA "*g*"are not equal.")
      end

      ## Add to abundancedf a column corresponding to the abundance of the RNA associated with g
      abundancedf[Symbol("R"*g)] = rbsabundance[:,1]
      #abundancedf["R"*g] = rbsabundance[:,1]
    else
      abundancedf[Symbol("R"*g)] = resultdf[:, Symbol("R"*g)]
      #abundancedf["R"*g] = resultdf[:, Symbol("R"*g)]
    end

    ## MAYBE to change if we don't make the disctinction between original and modified protein
    if genes["coding"][gid] == "PC"
      abundancedf[Symbol("P"*g)] = resultdf[:, Symbol("P"*g)]
      #abundancedf["P"*g] = resultdf[:, Symbol("P"*g)]

      if genes["PTMform"][gid] == "1"
        abundancedf[Symbol("Pm"*g)] = resultdf[:, Symbol("Pm"*g)]
        #abundancedf["Pm"*g] = resultdf[:, Symbol("Pm"*g)]
      end
    end
  end

  ## Include the abundance of the different regulatory complexes
  for comp in names(resultdf)[findall(x -> occursin(r"^C", x), map(String, names(resultdf)))]
    abundancedf[comp] = resultdf[:, comp]
    #abundancedf[string(comp)] = resultdf[:, comp]
  end

  return abundancedf
end


function juliaStochasticSimulation(stochmodel, QTLeffects, InitVar, genes, simtime;  modelname = "MySimulation", ntrials = 1, nepochs = -1, simalgorithm = "Direct")

  try

    ## If no nepochs provided, record abundance at each time units of the simulation
    if nepochs == -1
      nepochs = simtime
    end

    ## Convert the String name of the simulator to use into a BioSimulator object of type inheriting from Algorithm
    if in(simalgorithm, ["Direct", "FirstReaction", "NextReaction", "OptimizedDirect", "TauLeaping", "StepAnticipation"])
      simalgorithm = eval(Meta.parse("BioSimulator."*simalgorithm*"()"))
    else
      error("Specified algorithm is not implemented in module BioSimulator. Must be \"Direct\", \"FirstReaction\", \"NextReaction\", \"OptimizedDirect\", \"TauLeaping\" or \"StepAnticipation\".")
    end

    model = createBioSimModel(stochmodel, QTLeffects, InitVar, modelname)


    #println("JULIA> Running simulation ...")
    result = simulate(model, simalgorithm, time = convert(Float64, simtime), epochs = round(Int64, nepochs), trials = convert(Int64, ntrials))
    #println("JULIA> Done!")

    resultdf = DataFrame(result)

    abundancedf = transformSimRes2Abundance(resultdf, genes, stochmodel)

    return abundancedf

    catch err
        isa(err, InterruptException) || rethrow(err)
        return nothing
    end

end
