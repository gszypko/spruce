function [] = plotGridEvol(data,plotFreq)
% os (struct): contains grid info from 'os.grids.out' file

if nargin < 2, plotFreq = 1; end
data.grids.vars = data.grids.vars(1:plotFreq:end);
data.grids.time = data.grids.time(1:plotFreq:end);

% for each time point, plot 2D grids from MHD simulation
gridnames = {'i_n', 'dn', 'i_temp', 'e_temp', 'i_v_x', 'e_v_x', 'j_x', 'E_x'};
gridstr = {'n_i', 'n_i-n_e', 'T_i', 'T_e', 'v_i_x', 'v_e_x', 'j_x', 'E_x'};

% generate figure
f = filesep;
filepath = [data.folder f 'grid-evol'];
row = 2; 
col = 4; 
num = length(gridnames);
[fig,ax,an] = open_subplot(row,col,num,'Visible','on');
fig.Position = [1.303333333333333e+02,1.866000000000000e+02,1.194666666666667e+03,4.817333333333330e+02];
an.Position = [0.1595    0.9084    0.7230    0.0801];

frames = cell(1,length(data));
for k = 1:length([data.grids.vars.time])
    disp(['Plotting: ' num2str(k) '/' num2str(length([data.grids.vars.time]))])
    
    % update axis information
    iter = 0;
    for i = 1:size(ax,1)
        for j = 1:size(ax,2)
            iter = iter + 1;
            if iter > num, break, end

            cax = get_axis(fig,ax{i,j});
            xdata = data.grids.x_vec;
            ydata = data.grids.y_vec;
            zdata = data.grids.vars(k).(gridnames{iter});
            imagesc(xdata,ydata,zdata)
            colorbar
            cax.YDir = 'Normal';
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 10;
            if i == size(ax,1), xlabel('x (cm)'), end
            if j == 1, ylabel('y (cm)'), end
            title(gridstr{iter},'FontWeight','normal')

            if strcmp(gridnames{iter},'E_x'), c_lims = cax.CLim; end
            if strcmp(gridnames{iter},'E_ideal'), cax.CLim = c_lims; end
        end
    end

    dlm = ' - ';
    str1 = ['Iter = ' num2str(k-1)];
    str2 = ['t = ' num2str(data.grids.vars(k).time*1e6,'%.3g') '\mus = ' num2str(data.grids.vars(k).time/data.tau,'%.3g') '\tau_e_x_p'];
    an.String = [str1 dlm str2];
    an.Position = [0.159500000000000,0.929069291338582,0.723000000000000,0.080100000000000];

    frames{k} = getframe(fig);
end
close(fig)
write_video(filepath,frames);

end

