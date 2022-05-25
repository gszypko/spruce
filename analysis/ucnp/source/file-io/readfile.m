function [C] = readfile(path)

if ~isfile(path), error('<path> does not correspond to a file.'); end

fid = fopen(path,'r');
line_num = 0;
C = cell(1000,1);
while ~feof(fid)
    line_num = line_num + 1;
    C{line_num} = fgetl(fid);
    if length(C) == line_num, C{length(C)*2} = []; end
end
C(line_num+1:end) = [];

end