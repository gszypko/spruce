num = 1;
[fig,ax,an] = open_subplot(num);
iter = 0;
for i = 1:size(ax,1)
    for j = 1:size(ax,2)
        iter = iter + 1;
        if iter > num, break, end
        
        cax = get_axis(fig,ax{i,j});

        xdata = [];
        ydata = [];
        zdata = [];
        imagesc(xdata,ydata,zdata)
        colorbar;
        cax.YDir = 'Normal';
        cax.PlotBoxAspectRatio = [1 1 1];

        hold off
        l = get_line_specs(length(vars));
        for k = 1:length(vars)
            xdata = [];
            ydata = [];
            plot(xdata,ydata,'LineWidth',2,'MarkerSize',4,'Color',l(k).col,'MarkerFaceColor',l(k).col,'MarkerEdgeColor',l(k).col)
            hold on
        end

        cax.FontSize = 11;
        if i == size(ax,1), xlabel('label'), end
        if j == 1, ylabel('label'), end
        title('title','FontWeight','normal')
    end
end
an.String = [];

lgd = legend(vars);
