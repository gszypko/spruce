function [] = plotExpTransects(s,path)
figdir = [path filesep 'tr'];
mkdir(figdir)

%% Plot Raw LIF Fits
q = {'n','v','Ti'};
qstr = {'n (cm^-^3)','v (cm/s)','T_i (K)'};

num = length(q);
[fig,ax,an] = open_subplot(num);
fig.Position = [68.600000000000010,2.558000000000000e+02,1.162400000000000e+03,3.520000000000002e+02];
an.Position = [0.263261527873367,0.909963125938689,0.516400000000001,0.102000000000000];

for k = 1:length([s.tr.t])
    if ~isempty(s.tr(k).t)
        iter = 0;
        for i = 1:size(ax,1)
            for j = 1:size(ax,2)
                iter = iter + 1;
                if iter > num, break, end
                
                cax = get_axis(fig,ax{i,j});
        
                l = get_line_specs(2);
                xdata = [s.tr(k).x];
                ydata = [s.tr(k).(q{j})];
                plot(xdata,ydata,'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
                grid minor
        
                cax.FontSize = 11;
                cax.PlotBoxAspectRatio = [1 1 1];
                if i == size(ax,1), xlabel('x (cm)'), end
                if j == 1, ylabel('y (cm)'), end
                title(qstr{j},'FontWeight','normal')
            end
        end
        an.String = ['LIF Fits: t = ' num2str(s.tr(k).t*1e6,'%.2g') ' \mus'];
        exportgraphics(fig,[figdir filesep 'tr-' num2str(k) '.png'],'Resolution',300)
    end
end
close(fig)

end