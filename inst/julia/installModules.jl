using Pkg

importExprRD = Meta.parse("import Random");
try
	eval(importExprRD)
	print("Module Random installed.\n")
catch err
	print("Installing Julia module Random. This can take a few minutes.\n")
	Pkg.add("Random")
end

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
	if get(Pkg.dependencies(), Base.UUID("6aabf0a6-ec0e-58c4-b8bf-f645d5ac4ee8"), 3).version < v"0.9.3-beta"
		Pkg.rm("BioSimulator")
		Pkg.add(PackageSpec(url = "https://github.com/alanderos91/BioSimulator.jl.git", rev = "ace20d80a49701f6927e3b3441839d42d685b818"))
	end
	print("Module BioSimulator installed.\n")
catch err
	print("Installing Julia module Biosimulator. This can take a few minutes.\n")
	Pkg.add(PackageSpec(url = "https://github.com/alanderos91/BioSimulator.jl.git", rev = "ace20d80a49701f6927e3b3441839d42d685b818"))
 end

importExprSB = Meta.parse("import StatsBase");
try
	eval(importExprSB)
	print("Module StatsBase installed.\n")
catch err
	print("Installing Julia module statsBase. This can take a few minutes.\n")
	Pkg.add("StatsBase")
end
