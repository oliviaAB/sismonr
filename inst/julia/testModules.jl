using Pkg

# importExprCR = Meta.parse("import ClobberingReload")
# try
# 	eval(importExprCR)
# 	print("Module ClobberingReload installed.\n")
# catch err
# 	print("Installing Julia module ClobberingReload. This can take a few minutes.\n")
# 	Pkg.add("ClobberingReload")
# end

importExprJS = Meta.parse("import JSON");
try
	eval(importExprJS)
	print("Module JSON installed.\n")
catch err
	print("Installing Julia module JSON. This can take a few minutes.\n")
	Pkg.add("JSON")
end

importExprDF = Meta.parse("import DataFrames");
try
	eval(importExprDF)
	print("Module DataFrames installed.\n")
catch err
	print("Installing Julia module DataFrames. This can take a few minutes.\n")
	Pkg.add("DataFrames")
end

importExprBS = Meta.parse("import BioSimulator");
try
	eval(importExprBS)
	print("Module BioSimulator installed.\n")
catch err
	print("Installing Julia module Biosimulator. This can take a few minutes.\n")
	Pkg.add(PackageSpec(url = "https://github.com/alanderos91/BioSimulator.jl.git"))
 end

importExprSB = Meta.parse("import StatsBase");
try
	eval(importExprSB)
	print("Module StatsBase installed.\n")
catch err
	print("Installing Julia module statsBase. This can take a few minutes.\n")
	Pkg.add("StatsBase")
end