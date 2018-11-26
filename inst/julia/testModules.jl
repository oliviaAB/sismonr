importExprDF = parse("import DataFrames")
try
	eval(importExprDF)
	print("Module DataFrames installed.\n")
catch err
	print("Installing Julia module DataFrames. This can take a few minutes.\n")
	Pkg.add("DataFrames")
end

importExprBS = parse("import BioSimulator")
try
	eval(importExprBS)
	print("Module BioSimulator installed.\n")
catch err
	print("Installing Julia module Biosimulator. This can take a few minutes.\n")
	Pkg.clone("https://github.com/alanderos91/biosimulator.jl.git", "BioSimulator")
 end

importExprSB = parse("import StatsBase")
try
	eval(importExprSB)
	print("Module StatsBase installed.\n")
catch err
	print("Installing Julia module statsBase. This can take a few minutes.\n")
	Pkg.add("StatsBase")
end