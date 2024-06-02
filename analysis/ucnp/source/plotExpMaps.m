function [] = plotExpMaps(s,path)
figdir = [path filesep 'map'];
mkdir(figdir)
%% Plot Raw LIF Fits
q = {'n','v','Ti'};
qstr = {'n (cm^-^3)','v (cm/s)','T_i (K)'};

num = length(q);
[fig,ax,an] = open_subplot(num);
fig.Position = [68.600000000000010,2.626000000000000e+02,1.162400000000000e+03,3.452000000000001e+02];
an.Position = [0.233667584308328,0.855739976825029,0.516400000000001,0.102000000000000];

for k = 1:length([s.map.t])
    if ~isempty(s.map(k).t)
        iter = 0;
        for i = 1:size(ax,1)
            for j = 1:size(ax,2)
                iter = iter + 1;
                if iter > num, break, end
                
                cax = get_axis(fig,ax{i,j});
        
                xdata = [s.map(k).x];
                ydata = [s.map(k).y];
                zdata = [s.map(k).(q{j})];
                imagesc(xdata,ydata,zdata)
                colorbar;
                cax.YDir = 'Normal';
                cax.PlotBoxAspectRatio = [1 1 1];
        
                cax.FontSize = 11;
                if i == size(ax,1), xlabel('x (cm)'), end
                if j == 1, ylabel('y (cm)'), end
                title(qstr{j},'FontWeight','normal')
            end
        end
        an.String = ['LIF Fits: t = ' num2str(s.map(k).t*1e6,'%.2g') ' \mus'];
        exportgraphics(fig,[figdir filesep 'maps-' num2str(k) '.png'],'Resolution',300)
    end
end
close(fig)

%% Plot Density Hot Pixel Filtering
q = {'n','n_hpr'};
qstr = {'n_l_i_f (cm^-^3)','n_h_p_r (cm^-^3)'};

num = length(q);
[fig,ax,an] = open_subplot(num);
fig.Position = [2.106000000000000e+02,3.302000000000001e+02,9.512000000000002e+02,374];
an.Position = [0.223000000000000,0.896023529411764,0.516400000000000,0.102000000000000];

for k = 1:length([s.map.t])
    if ~isempty(s.map(k).t)
        iter = 0;
        for i = 1:size(ax,1)
            for j = 1:size(ax,2)
                iter = iter + 1;
                if iter > num, break, end
                
                cax = get_axis(fig,ax{i,j});
        
                xdata = [s.map(k).x];
                ydata = [s.map(k).y];
                zdata = [s.map(k).(q{j})];
                imagesc(xdata,ydata,zdata)
                colorbar;
                cax.YDir = 'Normal';
                cax.PlotBoxAspectRatio = [1 1 1];
        
                cax.FontSize = 11;
                if i == size(ax,1), xlabel('x (cm)'), end
                if j == 1, ylabel('y (cm)'), end
                title(qstr{j},'FontWeight','normal')
            end
        end
        an.String = ['Hot Pixel Filtering of n_l_i_f: t = ' num2str(s.map(k).t*1e6,'%.2g') ' \mus'];
        exportgraphics(fig,[figdir filesep 'maps-n-hpr-' num2str(k) '.png'],'Resolution',300)
    end
end
close(fig)

%% Compare LIF and Integrated Density Transects

num = 1;
[fig,ax,an] = open_subplot(num);
an.Position = [0.2230    0.8372    0.5164    0.1020];

for k = 1:length([s.map.t])
    if ~isempty(s.map(k).t)
        iter = 0;
        for i = 1:size(ax,1)
            for j = 1:size(ax,2)
                iter = iter + 1;
                if iter > num, break, end
                
                cax = get_axis(fig,ax{i,j});
        
                hold off
                q = {'img','img','map','map'};
                q1 = {'x','x','x','x_hpr'};
                q2 = {'n_x','n_x_sg','n_x','n_x_hpr'};
                q3 = {'n_i_n_t','n_i_n_t_,_s_g','n_l_i_f','n_l_i_f_,_h_p_r'};
                l = get_line_specs(length(q));
                for m = 1:length(q)
                    xdata = [s.(q{m})(k).(q1{m})];
                    ydata = [s.(q{m})(k).(q2{m})];
                    plot(xdata,ydata,'LineWidth',2,'MarkerSize',4,'Color',l(m).col,'MarkerFaceColor',l(m).col,'MarkerEdgeColor',l(m).col)
                    hold on
                end
        
                cax.FontSize = 11;
                if i == size(ax,1), xlabel('x (cm)'), end
                if j == 1, ylabel('n(y=0) (cm^-^3)'), end
                title(['t = ' num2str(s.map(k).t*1e6,'%.2g') ' \mus'],'FontWeight','normal')
            end
        end
        grid minor
        lgd = legend(q3);
        exportgraphics(fig,[figdir filesep 'tr-n-cmpr-' num2str(k) '.png'],'Resolution',300)
    end
end
close(fig)
end