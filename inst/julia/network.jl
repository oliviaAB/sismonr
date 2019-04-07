using StatsBase

## Test if the file is sourced in a Julia evaluator in R
function juliatest()
  println("Hello world from Julia!")
end

# ------------------------------------------- #
## FUNCTIONS FOR SAMPLING FROM DISTRIBUTIONS ##
# ------------------------------------------- #

## Sample from an discrete exponential distribution
function sampleexpon(n, lambda, max)
  if max > 0
    prob = (1/float(lambda))*exp.(-(1:max)/float(lambda))
    cumprob = cumsum(prob / sum(prob))
    rnb = rand(Int(n))
    res = Int64[]
    for nb in rnb
      push!(res, findfirst(x -> x >= nb, cumprob))
    end 
    return res
  else
    return [0 for i in 1:n]
  end
end


## Sample from an discrete power-law distribution
function samplepowerlaw(n, gamma, max)
  if max > 0
    prob = (1:max).^(-float(gamma))
    cumprob = cumsum(prob / sum(prob))
    rnb = rand(Int(n))
    res = Int64[]
    for nb in rnb
      push!(res, findfirst(x -> x >= nb, cumprob))
    end 
    return res
  else
    return [0 for i in 1:n]
  end
end


# ------------------------------------------------------------------------------------------------ #
## FUNCTIONS FOR COMPUTING PROBA FOR PREFERENTIAL ATTACHMENT OR "INVERSE PREFERENTIAL ATTACHMENT" ##
# ------------------------------------------------------------------------------------------------ #

# Linear preferential attachment from Barabasi-Albert
function probaPA(nodes, edg)
  ki = getInDeg(nodes, edg) .+ 1 # in-degree of each node, add 1 so that nodes with 0 in-degree can still be picked
  res = ki/sum(ki)
  res = res/sum(res)
  return res
end

# "Inverse" preferential attachment from Lachgar
function probaInvPA(nodes, edg)
  ki = getInDeg(nodes, edg) # in-degree of each node 
  if sum(ki) == 0
    res = fill(1, length(nodes))
  else
    res = 1 .- ki/sum(ki)
  end
  res = res / sum(res)
  return res
end


# ------------------------------------------------------- #
## FUNCTIONS FOR GETTING THE IN- AND OUT-DEGREE OF NODES ##
# ------------------------------------------------------- #
## Edges of a graph are stored in a mx2 array (where m is the number of edges)

function getInDeg(nodes, edg)
  to = sort(edg[:,2])
  res = [length(searchsorted(to, n)) for n in nodes]
  return res
end


function getOutDeg(nodes, edg)
  from = sort(edg[:,1])
  res = [length(searchsorted(from, n)) for n in nodes]
  return res
end

# Function for checking if an edge between two nodes exist
# from is an array of nodes, to is a single nodes
function isEdge(from, to, edg)
  edgS = edg[sortperm(edg[:,2]),:]
  toedg = edgS[searchsorted(edgS[:,2], to),:]
  toedg = sort(toedg[:,1])
  [length(searchsorted(toedg, f))>0 for f in from]
end



# ----------------------------------------------- #
## FUNCTIONS FOR GENERATING A REGULATORY NETWORK ##
# ----------------------------------------------- #



## Main function called from R
##Inputs: cf function nwgeneration
function juliaCreateNetwork(reacname, regPC, tarPC, indegPC, outdegPC, outdexexpPC, autoregprobaPC, twonodesloopPC,
  regNC, tarNC, indegNC, outdegNC, outdexexpNC, autoregprobaNC, twonodesloopNC,
  regcomplexes, compsize = nothing, compprob = nothing)

  
  edgPC = nwgeneration("PC", regPC, tarPC, indegPC, outdegPC, outdexexpPC, autoregprobaPC, twonodesloopPC)
  edgNC = nwgeneration("NC", regNC, tarNC, indegNC, outdegNC, outdexexpNC, autoregprobaNC, twonodesloopNC)

  if regcomplexes == "none"
    edg = vcat(edgPC, edgNC)
    edg[:,1] = [string(i) for i in edg[:,1]] ## Transform the integer ID of regulators into String ID 
    complexes = Dict()
  end

  if regcomplexes == "prot"
    resComp = compreg(edgPC, compsize, compprob, reacname)
    edg = vcat(resComp["newedg"], edgNC)
    complexes = resComp["complexes"]
  end

  if regcomplexes == "both"
    resComp = compreg(vcat(edgPC, edgNC), compsize, compprob, reacname)
    edg = resComp["newedg"]
    complexes = resComp["complexes"]
  end

  return Dict("edg" => edg, "complexes" => complexes)

end




## Generate a graph with specified in- and out- degree distribution
## Input:
##    - reg: list of regulator nodes
##    - target: list of target nodes
##    - indeg: string variable (either "exponential" or "powerlaw") specifying the type of preferential attachment used to construct the network
##    - outdeg: string variable (either "exponential" or "powerlaw") specifying the type of distribution from which the out-degree of regulators are sampled
##    - outdegexp: the exponent of the out-degree distribution
##    - autoregproba: probability that a regulatory molecule regulates itself
##    - twonodesloop: do we allow 2-nodes loops? can be true or false
##    - optional argument edg: a 2-D array of edges already existing, to take into account; if not given creates an empty array

## Output:
##    - edg: A 2D array of edges, 1st column: from, 2nd column: to, 3rd column: regBy
function nwgeneration(codingStatus, reg, target, indeg, outdeg, outdegexp, autoregproba, twonodesloop, edg = Array{Any}(undef, 0,2))

  # Ensure that reg and target are arrays
  if typeof(reg) == Int64
    reg = [reg]
  end
  if typeof(target) == Int64
    target = [target]
  end


  ## Get the function for sampling from the desired out- degree distribution
  if outdeg == "exponential"
    foutdeg = getfield(@__MODULE__, Symbol("sampleexpon"))
  elseif outdeg == "powerlaw"
    foutdeg = getfield(@__MODULE__, Symbol("samplepowerlaw"))  
  else
    error("Argument outdeg non-valid: must be \"exponential\" or \"powerlaw\"")
  end

  ## Get the function for computing the probability of each target node to receive a new incoming edge for the desired in-degree distribution
  ##  If in-degree is power law, use the model of preferential attachment from Barabasi-Albert
  ##  If in-degree is exponential, use the model of preferential attachment from Lachgar
  if indeg == "exponential"
    findeg = getfield(@__MODULE__, Symbol("probaInvPA"))
  elseif indeg == "powerlaw"
    findeg = getfield(@__MODULE__, Symbol("probaPA"))
  else
    error("Argument indeg non-valid: must be \"exponential\" or \"powerlaw\"")
  end

  ## Sample the number of target (out-degree) for each regulator
  out = foutdeg(length(reg), outdegexp, length(target))
  sort!(out, rev = true)

  ## Create the mx2 array of edges, 1st column = from, 2nd column = to
  # edg = Array{Int64}(undef, 0,2)

  ## For each regulator, sample from the target list its target according to its number of targets specified in the out variable
  for r in eachindex(reg)

    probTar = findeg(target, edg) # compute for each target the probability of being regulated by r
    # probTar = targetweight .* findeg(target, edg) # compute for each target the probability of being regulated by r, weighted by the weight of each target

    ## we exclude the regulator r from the list of potential targets (autoregulation is addressed later)
    if reg[r] in target
      probTar[findfirst(y -> y == reg[r], target)] = 0
    end

    ## How to deal with loops, i.e. if one or more target(s) already control the regulator r
    ## if twonodesloop == false we don't authorise the regulator to target a node that controls it
    if !twonodesloop
      exEdg = isEdge(target, reg[r], edg)
      probTar[exEdg] .= 0
    end

    ## Make sure that the out-degree of regulator r doesn't exceed the number of targets with non-null proba 
    out[r] = min(out[r], sum(probTar .> 0))

    ## Sample targets of regulator
    sa = StatsBase.sample(target, Weights(probTar), out[r], replace = false)

    ## Add the created edges in edg
    edg = vcat(edg, [fill(reg[r], out[r]) sa])

    ## Add an autoregulatory edge with probability autoregproba
    if rand() <= autoregproba
      edg = vcat(edg, [reg[r] reg[r]])
    end

  end

  edg = hcat(edg, repeat([codingStatus], outer = size(edg)[1]))
  return edg

end


## Generates a new graph with complexes as regulators from a graph given as input
## Input:
##   - edg: the existing graph
##   - compsize: size ( i.e. number of components) of the regulatory complexes
##   - compprob: the probability that regulators targeting a common gene form a regulatory complex
##   - reacname: ID of the type of regulation considered

## Output:
##   - edgcomp: a 2D array of edges, 1st column: from, 2nd column: to, 3rd column: regBy 

function compreg(edg, compsize, compprob, reacname)

  complexes = Dict()
  complexid = 1 ## to give a unique ID to each created complex
    rowstoremove = []
    edgtoadd = Array{Any}(undef, 0,3)
    compsize = Int(compsize)
    edg = edg[sortperm(edg[:,2]),:]

    target = unique(edg[:, 2])

    for tar in target
      temp = searchsorted(edg[:, 2], tar) ## row ID of edges coming to gene tar
      ntry = div(length(temp), compsize)  ## max number of complexes you can simultaneously create from the regulators of tar

      for i in 1:ntry
        if rand()<= compprob ## with probability compprob form a complex
          compo = sample(temp, compsize, replace=true) ## sample the components of the complex (in fact sample rows in the edg matrix)
          # temp = setdiff(temp, compo) ## remove the selected rows (regulators) from the list of possible components (future sampling)
          append!(rowstoremove, compo) ## the edges correpsponding to the selected components of the complex will be removed
          compid = string("C", reacname, complexid) ## create the new complex ID
          complexid +=1
          complexes[compid] = edg[compo, 1] ## in the dictionary of complexes add the composition (ie array of components) of the new complex
          edgtoadd = vcat(edgtoadd, [compid tar "C"]) ## create a new regulatory edge from the complex to the target, with regBy = "C" (complex)
        end
      end
    end

    ## remove the individual edges from regulators forming a complex
    edg = edg[setdiff(collect(1:size(edg)[1]), rowstoremove), :]
    ## add the new edges corresponding to regulation by complexes
    edg = vcat(edg, edgtoadd)

    edg[:,1] = [string(i) for i in edg[:,1]] ## Transform the integer ID of regulators into String ID 

    return Dict("newedg" => edg, "complexes" => complexes)
end
