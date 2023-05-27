function [] = plotGridEvol(data,flags)
f = filesep;

% get list of grid names from folder
[grid_names,grid_str] = readGridNames();

% check that user specified the variables to be plotted
if isfield(flags,'vars')
    if isempty(flags.vars), error('Variables must be specified for plotting.'); end
else
    error('Variables must be specified for plotting.')
end

% determine the label string for each variable to be plotted
varnames = flags.vars;
varstr = varnames;
for i = 1:length(varnames)
    ind = find(strcmp(grid_names,flags.vars{i}));
    if ~isempty(ind)
        varstr{i} = grid_str{ind};
    end
end

% generate figure
num = length(varnames);
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.117748618784530,0.943062827225131,0.764751381215472,0.066106464113451];

frames = cell(1,length(data));
for k = 1:length([data.grids.vars.time])
    disp(['Plotting: ' num2str(k) '/' num2str(length([data.grids.vars.time]))])
    
    % update axis information
    iter = 0;
    for i = 1:size(ax,1)
        for j = 1:size(ax,2)
            iter = iter + 1;
            if iter > num, break, end
            
            % interpolate data to uniform grid if base grid is non-uniform
            cax = get_axis(fig,ax{i,j});
            if data.grids.is_uniform
                xdata = data.grids.x_vec;
                ydata = data.grids.y_vec;
                zdata = data.grids.vars(k).(varnames{iter});
            else
                xdata = data.grids.x_uni;
                ydata = data.grids.y_uni;
                zdata = data.grids.uni_grid(data.grids.vars(k).(varnames{iter}));
            end
            
            % trim the domain for plotting according to user input
            ind_x = abs(xdata) < flags.plot_window(1);
            ind_y = abs(ydata) < flags.plot_window(2);
            xdata = xdata(ind_x);
            ydata = ydata(ind_y);
            zdata = zdata(ind_y,ind_x);            
            
            % do the plotting and update axis information
            imagesc(xdata,ydata,zdata)
            colorbar
            cax.YDir = 'Normal';
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 10;
            if i == size(ax,1), xlabel('x (cm)'), end
            if j == 1, ylabel('y (cm)'), end
            title(varstr{iter},'FontWeight','normal')
            cmin = min(zdata,[],'all');
            cmax = max(zdata,[],'all')*1.001;
            if cmax ~= 0, cax.CLim = [cmin cmax]; end
        end
    end

    dlm = ' - ';
    str1 = ['Iter = ' num2str(k-1)];
    str2 = ['t = ' num2str(data.grids.vars(k).time*1e6,'%.3g') ' \mus = ' num2str(data.grids.vars(k).time/data.tau,'%.3g') '\tau_e_x_p'];
    an.String = [str2];
    an.Position = [0.159500000000000,0.929069291338582,0.723000000000000,0.080100000000000];

    frames{k} = getframe(fig);
end
close(fig)
filepath = [data.folder f 'grid-evol'];
write_video(filepath,frames);

end

