using Pkg
using StatsBase
# using JLD ## temporary
using BioSimulator
using DataFrames

# ------------------------------------------------------------------------------------------------ #
##                       FUNCTIONS FOR SIMULATING THE STOCHASTIC SYSTEM                           ##
# ------------------------------------------------------------------------------------------------ #

function createBioSimModel(stochmodel, QTLeffects, InitAbundance, modelname)

  model = BioSimulator.Network(modelname)

  ## Add the species in the model, with their initial abundance
  for i in 1:length(stochmodel["species"])
    i0 = replace(stochmodel["initialconditions"][i], "QTLeffects" => "$QTLeffects")
    i0 = replace(i0, "InitAbundance" => "$InitAbundance")
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

  julia_version = VERSION >= v"1.2.0"


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
      if julia_version
        abundancedf[!, Symbol("R"*g)] = rbsabundance[:,1] ## Julia 1.2.0 syntax
      else
        abundancedf[Symbol("R"*g)] = rbsabundance[:,1] ## Julia 1.1.0 syntax
      end
    else
      if julia_version
        abundancedf[!, Symbol("R"*g)] = resultdf[:, Symbol("R"*g)] ## Julia 1.2.0 syntax
      else
        abundancedf[Symbol("R"*g)] = resultdf[:, Symbol("R"*g)] ## Julia 1.1.0 syntax
      end
    end

    ## MAYBE to change if we don't make the disctinction between original and modified protein
    if genes["coding"][gid] == "PC"
      if julia_version
        abundancedf[!, Symbol("P"*g)] = resultdf[:, Symbol("P"*g)] ## Julia 1.2.0 syntax
      else
        abundancedf[Symbol("P"*g)] = resultdf[:, Symbol("P"*g)] ## Julia 1.1.0 syntax
      end

      if genes["PTMform"][gid] == "1"
        if julia_version
          abundancedf[!, Symbol("Pm"*g)] = resultdf[:, Symbol("Pm"*g)] ## Julia 1.2.0 syntax
        else
          abundancedf[Symbol("Pm"*g)] = resultdf[:, Symbol("Pm"*g)] ## Julia 1.1.0 syntax
        end
      end
    end
  end

  ## Include the abundance of the different regulatory complexes
  for comp in names(resultdf)[findall(x -> occursin(r"^C", x), map(String, names(resultdf)))]
    if julia_version
      abundancedf[!, comp] = resultdf[:, comp] ## Julia 1.2.0 syntax
    else
      abundancedf[comp] = resultdf[:, comp] ## Julia 1.1.0 syntax
    end
  end

  return abundancedf
end


function juliaStochasticSimulation(stochmodel, QTLeffects, InitAbundance, genes, simtime; modelname = "MySimulation", ntrials = 1, nepochs = -1, simalgorithm = "Direct", seed = nothing)

  try

    Random.seed!(seed)

    ## If no nepochs provided, record abundance at each time units of the simulation
    if nepochs == -1
      nepochs = simtime
    end

#     if Pkg.installed()["BioSimulator"] < v"0.9" ## With the older version of BioSimulator
#
#       ## Convert the String name of the simulator to use into a BioSimulator object of type inheriting from Algorithm
#       if in(simalgorithm, ["Direct", "FirstReaction", "NextReaction", "OptimizedDirect", "TauLeaping", "StepAnticipation"])
#         simalgorithm = eval(Meta.parse("BioSimulator."*simalgorithm*"()"))
#       else
#         error("Specified algorithm is not implemented in module BioSimulator. Must be \"Direct\", \"FirstReaction\", \"NextReaction\", \"OptimizedDirect\", \"TauLeaping\" or \"StepAnticipation\".")
#       end
#
#       model = createBioSimModel(stochmodel, QTLeffects, InitAbundance, modelname)

      #println("JULIA> Running simulation ...")
#       result = simulate(model, simalgorithm, time = convert(Float64, simtime), epochs = round(Int64, nepochs), trials = convert(Int64, ntrials))
      #println("JULIA> Done!")

#       resultdf = DataFrame(result)
#     end

    ## Convert the String name of the simulator to use into a BioSimulator object of type inheriting from Algorithm
    if in(simalgorithm, ["Direct", "EnhancedDirect", "SortingDirect", "FirstReaction", "NextReaction", "TauLeapingDG2001", "TauLeapingDGLP2003", "StepAnticipation", "HybridSAL"])
      simalgorithm = eval(Meta.parse("BioSimulator."*simalgorithm*"()"))
    else
      error("Specified algorithm is not implemented in module BioSimulator. Must be \"Direct\", \"EnhancedDirect\", \"SortingDirect\", \"FirstReaction\", \"NextReaction\", \"TauLeapingDG2001\", \"TauLeapingDGLP2003\", \"StepAnticipation\" or \"HybridSAL\".")
    end

    model = createBioSimModel(stochmodel, QTLeffects, InitAbundance, modelname)

    ## We have to compute the time at which we want to record molecules abundance
    t_step = floor(simtime / nepochs)
    t_to_keep = collect(0:t_step:simtime)
    if t_to_keep[end] != simtime
      t_to_keep[end] = simtime
    end

    #println("JULIA> Running simulation ...")
    result = [simulate(model, simalgorithm, tfinal = convert(Float64, simtime), save_points = t_to_keep) for _ in 1:convert(Int64, ntrials)] #  trials = convert(Int64, ntrials)
    #println("JULIA> Done!")

    resultdf = DataFrame(result)
    new_col_names = Dict(zip(["X$i" for i in 1:length(keys(model.species_list))], collect(keys(model.species_list))))
    new_col_names["t"] = :time
    rename!(resultdf, new_col_names)

    abundancedf = transformSimRes2Abundance(resultdf, genes, stochmodel)

    return abundancedf

    catch err
        isa(err, InterruptException) || rethrow(err)
        return nothing
    end

end
