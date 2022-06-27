function [s] = spreadsheet2struct(inp,fields)

num_rows = size(inp,1) - 1;
names = inp(1,:);

s = struct;
for i = 1:length(fields)
    col_ind = strcmp(names,fields(i));
    for j = 1:num_rows
        s(j).(fields{i}) = inp{j+1,col_ind};        
    end
end

end