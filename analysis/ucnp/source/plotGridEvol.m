function [] = plotGridEvol(os)
% os (struct): contains grid info from 'os.grids.out' file

% set up save directories
f = filesep;
savedir = [os.folder f 'figs-grid-evol'];
mkdir(savedir)

% for each time point, plot 2D grids from MHD simulation
gridnames = {'n','temp','press','v_x','v_y','thermal_energy'};
gridstr = {'n','T','P','v_x','v_y','E_t_h'};
for k = 1:length([os.grids.vars.time])
    disp(['Plotting: ' num2str(k) '/' num2str(length([os.grids.vars.time]))])
    
    fig = figure('Visible','off');
    fig.Position = [158   329   894   467];
    fig.Color = [1 1 1];

    row = 2;
    col = 3;
    num = length(gridnames);

    ax = cell(row,col);
    iter = 0;
    for i = 1:row
        for j = 1:col
            if iter > num - 1, break, end
            iter = iter + 1;
            ax{i,j} = subplot(row,col,iter);

            xdata = os.grids.pos_x;
            ydata = os.grids.pos_y;
            zdata = os.grids.vars(k).(gridnames{iter});
            surf = pcolor(xdata,ydata,zdata);
            surf.EdgeColor = surf.FaceColor;
            
            cb = colorbar;
            ax{i,j}.YDir = 'Normal';

            ax{i,j}.PlotBoxAspectRatio = [1 1 1];
            ax{i,j}.FontSize = 12;

            if i == row, xlabel('x (cm)'), end
            if j == 1, ylabel('y (cm)'), end
            title(gridstr{iter},'FontWeight','normal')

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

