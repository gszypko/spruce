function [out] = readOutFile(path,Ng,grids_to_load,time_window,time_interval)
% directory (string): full path to directory containing 'plasma.settings'
% Ng (int): number of ghost cells to remove from domain edges

% ensure that plasma.settings file exists
if ~endsWith(path,'.out'), error('File extension must be .out'); end


% read contents from file
C = readfile(path);

% read preamble
dlm = ',';
temp = split(C{2},dlm);
out.Nx = str2double(temp{1}) - 2*Ng;
out.Ny = str2double(temp{2}) - 2*Ng;

% read plasma domain grids
static_vars = {'pos_x','pos_y','be_x','be_y','be_z'};
static_vars_found = 0;
for i = 1:length(C)
    ind = find(strcmp(C{i},static_vars), 1);
    if ~isempty(ind)
        static_vars_found = static_vars_found + 1;
        out.(static_vars{ind}) = getgrid(C,i,Ng,out.Nx,out.Ny);
        if static_vars_found == length(static_vars)
            break
        end
    end
end

% determine indices for each time point
time_ind = zeros(1,1e4);
time_iter = 0;
for i = 1:length(C)
    if startsWith(C{i},'t=')
        time_iter = time_iter + 1;
        time_ind(time_iter) = i;
    end
end
time_ind = time_ind(1:time_iter);

% trim time window and apply interval using inputs
time_start = round(time_window(1)*time_iter/100+1);
time_end = round(time_window(2)*time_iter/100);
time_ind = time_ind(time_start:time_interval:time_end);
time_ind_interval = unique(diff(time_ind))/time_interval;
if length(time_ind_interval)>1, error('Time index spacing is not uniform.'); end
C(1:time_ind(1)-1) = [];
time_ind = time_ind - time_ind(1) + 1;

% load the required grids for the chosen time points
for i = 1:length(time_ind)
    disp(['Loading Grids for Time Point: ' num2str(i) '/' num2str(length(time_ind))])
    out.vars(i).time = str2double(extractAfter(C{time_ind(i)},'t='));
    % begin reading grids for this time point
    vars_found = 0;
    for j = time_ind(i)+1:time_ind(i)+time_ind_interval-1
        % identify if current line corresponds to variable name
        ind = find(strcmp(C{j},grids_to_load));
        if length(ind)==1
            out.vars(i).(grids_to_load{ind}) = getgrid(C,j,Ng,out.Nx,out.Ny);
            vars_found = vars_found + 1;
            % terminate grid reading for this time point if all variables found
            if vars_found == length(grids_to_load)
                % trim contents of mhd.out file 'C' to save on RAM
                if i ~= length(time_ind)
                    C(1:time_ind(i+1)-1) = [];
                    time_ind = time_ind - time_ind(i+1) + 1;
                end
                break
            end
        end
    end
    % ensure that all variables were found for this time point
    if vars_found ~= length(grids_to_load), error('Not all grids were found.'); end
end

out.time = [out.vars.time];
out.x_vec = out.pos_x(1,:);
out.y_vec = out.pos_y(:,1)';

% check if grids are uniform or non-uniform
dx = diff(out.x_vec);
dy = diff(out.y_vec);
uniform_threshold = .001;
ind_x = min(abs(dx - mean(dx)) < uniform_threshold*mean(dx));
ind_y = min(abs(dy - mean(dy)) < uniform_threshold*mean(dy));
out.is_uniform = ind_x && ind_y;

% handle when grids are non-uniform
out.Nx_uni = (max(out.x_vec)-min(out.x_vec))/min(diff(out.x_vec));
out.Ny_uni = (max(out.y_vec)-min(out.y_vec))/min(diff(out.y_vec));
out.x_uni = linspace(min(out.x_vec),max(out.x_vec),out.Nx_uni);
out.y_uni = linspace(min(out.y_vec),max(out.y_vec),out.Ny_uni);
[X_uni,Y_uni] = meshgrid(out.x_uni,out.y_uni);
out.uni_grid = @(var) interp2(out.x_vec,out.y_vec,var,X_uni,Y_uni);

end