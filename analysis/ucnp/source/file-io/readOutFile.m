function [out] = readOutFile(path,Ng)
% directory (string): full path to directory containing 'plasma.settings'
% opt (bool): (true) trim ghost cells from matrices (false) do not
% out (struct): fields of <out> contain plasma quantities in cgs units, see below

% ensure that plasma.settings file exists
if ~endsWith(path,'.out'), error('File extension must be .out'); end


% read in contents from file
C = readfile(path);
dlm = ',';
temp = split(C{2},dlm);
out.Nx = str2double(temp{1}) - 2*Ng;
out.Ny = str2double(temp{2}) - 2*Ng;
zeromat = zeros(out.Ny,out.Nx);

static_vars = {'pos_x','pos_y','be_x','be_y'};
static_ind = zeros(size(static_vars));
iter = 0;
for i = 1:length(C)
    ind = find(strcmp(C{i},static_vars), 1);
    if ~isempty(ind)
        iter = iter + 1;
        static_ind(iter) = i;
        if iter == length(static_vars)
            break
        end
    end
end

for i = 1:length(static_vars)
    grid = zeromat;
    iter = 0;
    for j = static_ind(i)+1+Ng:static_ind(i)+Ng+out.Nx
        iter = iter + 1;
        temp = split(C{j},dlm);
        line = sscanf(sprintf(' %s',temp{:}),'%f',[1,Inf]);
        grid(:,iter) = line(Ng+1:end-Ng);
    end
    out.(static_vars{i}) = grid;
end

iter_t = 0;
for i = 1:length(C)
    if startsWith(C{i},'t=')
        time_start = i;
        break
    end
end

for i = time_start:length(C)
    if startsWith(C{i},'t=')
        iter_t = iter_t + 1;
        out.vars(iter_t).time = str2double(extractAfter(C{i},'t='));
    else
        if length(split(C{i},dlm)) == 1
            out.vars(iter_t).(C{i}) = getgrid(C,i,Ng,out.Nx,out.Ny);
        end
    end      
end

out.time = [out.vars.time];
out.x_vec = out.pos_x(1,:);
out.y_vec = out.pos_y(:,1)';

end