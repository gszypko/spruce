function [C] = readfile(path)

if ~isfile(path), error('<path> does not correspond to a file.'); end

fid = fopen(path,'r');
iter = 0;
C = cell(1000,1);
while ~feof(fid)
    iter = iter + 1;
    C{iter} = fgetl(fid);
    if length(C) == iter, C{length(C)*2} = []; end
end
C(iter+1:end) = [];

end