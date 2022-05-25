% opening a subplot
row = 2; 
col = 2; 
num = row*col;
[fig,ax,an] = open_subplot(row,col,num);
fig.Position = [471.6667  349.0000  754.0000  574.6667];
move_ax(ax(1:end-1,:),0,.05)
an.Position = [0.1595    0.9084    0.7230    0.0801];

% updating axis information
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
        cax.FontSize = 12;
        if i == size(ax,1), xlabel('label'), end
        if j == 1, ylabel('label'), end
        title('title','FontWeight','normal')
    end
end
an.String = [];