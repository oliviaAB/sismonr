using Pkg

missingmdl = 0;

importExprRD = Meta.parse("import Random");
try
	eval(importExprRD);
catch err
	global missingmdl = missingmdl + 1;
end

importExprJS = Meta.parse("import JSON");
try
	eval(importExprJS);
catch err
	global missingmdl = missingmdl + 1;
end

importExprDF = Meta.parse("import DataFrames");
try
	eval(importExprDF);
catch err
	global missingmdl = missingmdl + 1;
end

importExprBS = Meta.parse("import BioSimulator");
try
	eval(importExprBS);
	if get(Pkg.dependencies(), Base.UUID("6aabf0a6-ec0e-58c4-b8bf-f645d5ac4ee8"), 3).version < v"0.9.3-beta"
		global missingmdl = missingmdl + 1;
	end
catch err
	global missingmdl = missingmdl + 1;
 end

importExprSB = Meta.parse("import StatsBase");
try
	eval(importExprSB);
catch err
	global missingmdl = missingmdl + 1;
end

print(missingmdl)
