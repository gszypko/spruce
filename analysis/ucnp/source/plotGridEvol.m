function [] = plotGridEvol(os)
% os (struct): contains grid info from 'os.grids.out' file

%% Plot 2D Grids from MHD Simulation
% set up save directories
f = filesep;
savedir = [os.folder f 'figs-grid-evol'];
mkdir(savedir)

% for each time point, plot 2D grids from MHD simulation
gridNames = {'n','temp','press','v_x','v_y','thermal_energy'};
for k = 1:length([os.grids.time])
    disp(['Plotting: ' num2str(k) '/' num2str(length([os.grids.time]))])
    
    fig = figure('Visible','off');
    fig.Position = [432         207        1039         525];
    fig.Color = [1 1 1];

    row = 2;
    col = 3;
    num = row*col;

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
            pcolor(xdata,ydata,zdata)
            shading interp
            
            cb = colorbar;
            cb.FontSize = 12;
            cb.Label.String = gridNames{iter};
            ax{i,j}.YDir = 'Normal';

            ax{i,j}.PlotBoxAspectRatio = [1 1 1];
            ax{i,j}.FontSize = 12;

            if i == row, xlabel('x (cm)'), end
            if j == 1, ylabel('y (cm)'), end

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

