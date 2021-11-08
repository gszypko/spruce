function [] = plotGridEvol(os)
% os (struct): contains grid info from 'os.grids.out' file

%% Plot 2D Grids from MHD Simulation
% set up save directories
f = filesep;
savedir = [os.folder f 'figs-grid-evol'];
mkdir(savedir)

% for each time point, plot 2D grids from MHD simulation
gridNames = {'b_x','b_y','b_magnitude','mag_press','mag_press_x','mag_press_y','divBB_x','divBB_y','mag_diff_x','mag_diff_y','divB','press_x','press_y','lorentz_force_x','lorentz_force_y'};
for k = 1:length([os.grids.time])
    disp(['Plotting: ' num2str(k) '/' num2str(length([os.grids.time]))])
    
    fig = figure('Visible','off');
    fig.Position = [163         124        1466         677];
    fig.Color = [1 1 1];

    row = 4;
    col = 4;
    num = length(gridNames);

    ax = cell(row,col);
    iter = 0;
    for i = 1:row
        for j = 1:col
            if iter > num - 1, break, end
            iter = iter + 1;
            ax{i,j} = subplot(row,col,iter);

            xdata = os.grids.pos_x;
            ydata = os.grids.pos_y;
            zdata = os.grids.vars(k).(gridNames{iter});
            s = pcolor(xdata,ydata,zdata);
            s.EdgeColor = 'flat';
            
            cb = colorbar;
            ax{i,j}.YDir = 'Normal';

            ax{i,j}.PlotBoxAspectRatio = [1 1 1];
            ax{i,j}.FontSize = 12;

            if i == row, xlabel('x (cm)'), end
            if j == 1, ylabel('y (cm)'), end
            title(gridNames{iter})

        end
    end

    an = annotation('textbox');
    an.Position = [0.0058    0.9219    0.9919    0.0686];
    an.HorizontalAlignment = 'center';
    an.VerticalAlignment = 'middle';
    an.LineStyle = 'none';
    an.FontSize = 12;

    dlm = ' - ';
    str1 = ['Iter = ' num2str(k-1)];
    str2 = ['t = ' num2str(os.grids.vars(k).time*1e6) '\mus'];
    an.String = [str1 dlm str2];

    savepath = [savedir f 'tEvol' num2str(k-1) '.png'];
    saveas(fig,savepath)
    close(fig)
end
    
end

