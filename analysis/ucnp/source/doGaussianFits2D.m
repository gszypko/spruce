function [data] = doGaussianFits2D(data,flags)

% determine which variable corresponds to ion density
n = data.i.n;

% fit 2D Gaussian to ion density distribution for each time point
gauss_fits = struct();
for i = 1:length(data.grids.time)    
    disp(['2D Gaussian Fits: ' num2str(i) '/' num2str(length(data.grids.time))])
    
    x = data.grids.x_vec;
    y = data.grids.y_vec;
    img = data.grids.vars(i).(n);
    fit = fitImgWithGaussian(x,y,img);
    
    gauss_fits(i).time = data.grids.vars(i).time;
    fields = {'x','y','img','imgfit','imgres','fit','amp','x0','y0','sigx','sigy','sigxErr','sigyErr','offset'};
    for j = 1:length(fields)
        gauss_fits(i).(fields{j}) = fit.(fields{j});
    end
end

% plot gaussian fits
% generate figure
num = 5;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];
ax{2}.Position(1) = ax{2}.Position(1) - .0192;
ax{4}.Position(1) = ax{4}.Position(1) - .0192;
an = annotation('textbox');
an.Position = [0.687883435582822,0.114100000000000,0.230216564417178,0.329900000000000];
an.HorizontalAlignment = 'left';
an.VerticalAlignment = 'top';
an.LineStyle = 'none';
an.FontSize = 12;

frames = cell(1,length(gauss_fits));
for i = 1:length(gauss_fits)    
    disp(['Plotting 2D Gaussian Fits: ' num2str(i) '/' num2str(length(gauss_fits))])
    
    hold([ax{:}],'off')
    [~,indx] = min(abs(gauss_fits(i).x));
    [~,indy] = min(abs(gauss_fits(i).y));
    imgx = gauss_fits(i).img(indy,:);
    imgy = gauss_fits(i).img(:,indx)';
    imgfitx = gauss_fits(i).imgfit(indy,:);
    imgfity = gauss_fits(i).imgfit(:,indx)';
    
    lgdstr = {'Sim.','Fit'};
    l = get_line_specs(length(lgdstr));

    cax = get_axis(fig,ax{1});
    if data.grids.is_uniform
        zdata = gauss_fits(i).img;
    else
        zdata = data.grids.uni_grid(gauss_fits(i).img);
    end
    imagesc(x,y,zdata./1e8)
    cax.PlotBoxAspectRatio = [1 1 1];
    cax.YDir = 'Normal';
    shading interp
    cb = colorbar;
    xlabel('x (cm)')
    ylabel('y (cm)')
    title('Simulation','FontWeight','normal')
    
    cax = get_axis(fig,ax{3});
    if data.grids.is_uniform
        zdata = gauss_fits(i).imgfit;
    else
        zdata = data.grids.uni_grid(gauss_fits(i).imgfit);
    end
    imagesc(x,y,zdata./1e8);
    cax.PlotBoxAspectRatio = [1 1 1];
    cax.YDir = 'Normal';
    shading interp
    cb = colorbar;
    xlabel('x (cm)')
    title('Fit','FontWeight','normal')

    cax = get_axis(fig,ax{5});
    if data.grids.is_uniform
        zdata = gauss_fits(i).imgres;
    else
        zdata = data.grids.uni_grid(gauss_fits(i).imgres);
    end
    imagesc(x,y,zdata./1e8);
    cax.PlotBoxAspectRatio = [1 1 1];
    cax.YDir = 'Normal';
    shading interp
    cb = colorbar;
    cb.Label.String = 'n (10^8 cm^-^3)';
    xlabel('x (cm)')
    title('Residual','FontWeight','normal')

    cax = get_axis(fig,ax{2});
    plot(x,imgx./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
    hold on
    plot(x,imgfitx./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
    cax.PlotBoxAspectRatio = [1 1 1];
    cax.YDir = 'Normal';
    title('x axis','FontWeight','normal')
    xlabel('x (cm)')
    ylabel('n (10^8 cm^-^3)')
    
    cax = get_axis(fig,ax{4});
    plot(y,imgy./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
    hold on
    plot(y,imgfity./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
    cax.PlotBoxAspectRatio = [1 1 1];
    cax.YDir = 'Normal';
    title('y axis','FontWeight','normal')
    lgd = legend(lgdstr);
    lgd.Position = [0.735103990127022,0.127535831274251,0.072661043511578,0.060403964650340];
    lgd.FontSize = 12;
    xlabel('y (cm)')

    str1 = ['t = ' num2str(gauss_fits(i).time*1e6,'%.2g') ' \mus = ' num2str(gauss_fits(i).time/data.tau,'%.2g') ' \tau'];
    str2 = ['n_0 = ' num2str(gauss_fits(i).amp/1e8,'%.2g') '\times10^8 cm^{-3}'];
    str3 = ['\sigma_x = ' num2str(gauss_fits(i).sigx*10,'%.2g') ' \pm ' num2str(gauss_fits(i).sigxErr*10,'%.2g') ' mm'];
    str4 = ['\sigma_y = ' num2str(gauss_fits(i).sigy*10,'%.2g') ' \pm ' num2str(gauss_fits(i).sigyErr*10,'%.2g') ' mm'];
    an.String = {str1,str2,str3,str4};
        
    frames{i} = getframe(fig);
end
close(fig)
write_video([data.folder filesep 'gauss-fits'],frames);
save([data.folder filesep 'gauss-fits.mat'],"gauss_fits",'-mat');

end