function [] = analyzeExpData(path)

[s] = loadExpData(path,'os');
if isfield(s,'map'), plotExpMaps(s,path); end
if isfield(s,'tr'), plotExpTransects(s,path); end

end