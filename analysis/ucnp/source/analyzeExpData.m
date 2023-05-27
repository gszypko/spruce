function [] = analyzeExpData(path)

[s] = loadExpData(path,'os-init');
if isfield(s,'map'), plotExpMaps(s,path); end

end