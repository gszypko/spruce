function [out] = readConfigFile(path)

C = readfile(path);
ind = startsWith(C,'equation_set');
out = struct;
out.eq_set = C{ind};
out.eq_set = strtrim(extractAfter(out.eq_set,'='));

end