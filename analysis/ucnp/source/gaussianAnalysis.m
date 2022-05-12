function [] = gaussianAnalysis(data)

% set up save directories
f = filesep;
savedir = [data.folder f 'gauss-fits'];
mkdir(savedir)

s = struct;
s.t = data.grids.time;
s.x = data.grids.x_vec;
s.y = data.grids.y_vec;
s.Te0 = data.settings.Te;
s.sig0 = data.sig0;
s.tau = data.tau;

figiter = 0;
%% Fit Density Distribution with 2D Gaussian to Extract RMS Width

for i = 1:length(data.grids.time)
    disp(['2D Gaussian Fits: ' num2str(i) '/' num2str(length(data.grids.time))])
    x = data.grids.x_vec;
    y = data.grids.y_vec;
    X = data.grids.pos_x;
    Y = data.grids.pos_y;
    img = data.grids.vars(i).n;
    [fit(i)] = fitImgWithGaussian(x,y,img);
    
    [~,indx] = min(abs(x));
    [~,indy] = min(abs(y));
    imgx = img(indy,:);
    imgy = img(:,indx)';
    imgfitx = fit(i).imgfit(indy,:);
    imgfity = fit(i).imgfit(:,indx)';
    
    lgdstr = {'MHD','Fit'};
    l = getLineSpecs(length(lgdstr));
    
    fig = figure('Visible','off');
    figiter = figiter + 1;
    fig.Position = [195         259        1027         491];
    fig.Color = [1 1 1];

    ax = subplot(2,3,1);
    imagesc(x,y,img./1e8)
    ax.PlotBoxAspectRatio = [1 1 1];
    shading interp
    cb = colorbar;
    
    ax = subplot(2,3,2);
    imagesc(x,y,fit(i).imgfit./1e8);
    ax.PlotBoxAspectRatio = [1 1 1];
    shading interp
    cb = colorbar;
    
    ax = subplot(2,3,3);
    imagesc(x,y,fit(i).imgres./1e8);
    ax.PlotBoxAspectRatio = [1 1 1];
    shading interp
    cb = colorbar;
    cb.Label.String = 'n (10^8 cm^-^3)';
    
    ax = subplot(2,3,4);
    ax.Position(1) = ax.Position(1) - .025;
    plot(x,imgx./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
    hold on
    plot(x,imgfitx./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
    ax.PlotBoxAspectRatio = [1 1 1];
    title('x axis','FontWeight','normal')
    ylabel('n (10^8 cm^-^3)')
    
    ax = subplot(2,3,5);
    ax.Position(1) = ax.Position(1) - .025;
    plot(y,imgy./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
    hold on
    plot(y,imgfity./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
    ax.PlotBoxAspectRatio = [1 1 1];
    title('y axis','FontWeight','normal')
    legend(lgdstr)
    
    an = annotation('textbox');
    an.Position = [0.6654    0.1141    0.2527    0.3299];
    an.HorizontalAlignment = 'left';
    an.VerticalAlignment = 'top';
    an.LineStyle = 'none';
    an.FontSize = 12;
    
    str1 = ['t = ' num2str(data.grids.time(i)*1e6,'%.2g') ' \mus = ' num2str(data.grids.time(i)/s.tau,'%.2g') ' \tau_{exp}'];
    str2 = ['n_{max} = ' num2str(fit(i).amp/1e8,'%.2g') '\times10^8 cm^{-3}'];
    str3 = ['\sigma_x = ' num2str(fit(i).sigx*10,'%.2g') ' \pm ' num2str(fit(i).sigxErr*10,'%.2g') ' mm'];
    str4 = ['\sigma_y = ' num2str(fit(i).sigy*10,'%.2g') ' \pm ' num2str(fit(i).sigyErr*10,'%.2g') ' mm'];
    an.String = {str1,str2,str3,str4};
    
    savepath = [savedir f 'fig' num2str(figiter) '-tEvol' num2str(i-1) '.png'];
    saveas(fig,savepath)
    close(fig)
end

%% Compile MHD Data and Vlasov Theory for Comparison

s.data = struct;
for i = 1:length(s.t)
    s.data(i).sigx = fit(i).sigx;
    s.data(i).sigy = fit(i).sigy;
    s.data(i).n = fit(i).amp;
    [~,indx] = min(abs(data.grids.x_vec));
    [~,indy] = min(abs(data.grids.y_vec));
    s.data(i).Te = data.grids.vars(i).temp(indy,indx);
    s.data(i).vx = data.grids.vars(i).v_x(indy,:);
    s.data(i).vy = data.grids.vars(i).v_y(:,indx)';
end

s.theory = struct;
[sig,Te,v] = getUCNPExpansion(s.t,s.sig0,s.Te0,s.x,s.y);
for i = 1:length(s.t)
    s.theory(i).sigx = sig(i);
    s.theory(i).sigy = sig(i);
    s.theory(i).Te = Te(i);
    s.theory(i).vx = v(i).vx;
    s.theory(i).vy = v(i).vy;
end

%% Plot Summary of Expansion for MHD Sims Against Vlasov Theory

fig = figure('Visible','off');
figiter = figiter + 1;
fig.Position = [432   361   996   371];
fig.Color = [1 1 1];

q = {'data','theory'};
qstr = {'MHD','Vlasov'};
colvar = {'sigx','sigy','Te'};
colstr = {'\sigma_x (cm)','\sigma_y (cm)','T_e (K)'};
rowvar = {''};
row = 1;
col = length(colvar);
num = row*col;

ax = cell(row,col);
iter = 0;
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        ax{i,j} = subplot(row,col,iter);
        hold on
        lgdstr = qstr;
        l = getLineSpecs(length(lgdstr));
        
        for k = 1:length(q)
            xdata = [s.t]./s.tau;
            ydata = [s.(q{k}).(colvar{j})];
            lp = l(k);
            plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
        end
        ax{i,j}.PlotBoxAspectRatio = [1 1 1];
        ax{i,j}.FontSize = 12;

        if i == row, xlabel('t / \tau_{exp}'), end
        ylabel(colstr{j})
%         title('title')

    end
end

% add legend
lgd = legend(lgdstr);
lgd.Position = [0.7824    0.6379    0.1021    0.1092];

savename = [savedir f 'fig' num2str(figiter) '-expansion-summary.png'];
saveas(fig,savename)
close(fig)

%% Plot Velocity Transects


for k = 1:length(s.t)
    disp(['Plotting velocity transects: ' num2str(k) '/' num2str(length(s.t))])
    fig = figure('Visible','off');
    figiter = figiter + 1;
    fig.Position = [432   391   665   341];
    fig.Color = [1 1 1];
    
    q = {'data','theory'};
    qstr = {'MHD','Vlasov'};
    colvar1 = {'x','y'};
    colvar1str = {'x (cm)','y (cm)'};
    colvar2 = {'vx','vy'};
    colvar2str = {'v_x (cm/s)','v_y (cm/s)'};
    rowvar = {''};
    row = 1;
    col = 2;
    num = row*col;

    ax = cell(row,col);
    iter = 0;
    for i = 1:row
        for j = 1:col
            if iter > num - 1, break, end
            iter = iter + 1;
            ax{i,j} = subplot(row,col,iter);
            hold on
            lgdstr = qstr;
            l = getLineSpecs(length(lgdstr));
            
            for indq = 1:length(q)
                xdata = [s.(colvar1{j})];
                ydata = [s.(q{indq})(k).(colvar2{j})];
                lp = l(indq);
                plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
            end
            ax{i,j}.PlotBoxAspectRatio = [1 1 1];
            ax{i,j}.FontSize = 12;

            xlabel(colvar1str{j})
            ylabel(colvar2str{j})

        end
    end

    % add legend
    lgd = legend(lgdstr);
    lgd.Position = [0.7473    0.6845    0.1529    0.1188];

    an = annotation('textbox');
    an.Position = [0.0058    0.8762    0.9919    0.1143];
    an.HorizontalAlignment = 'center';
    an.VerticalAlignment = 'middle';
    an.LineStyle = 'none';
    an.FontSize = 12;
    an.String = ['Central Velocity Transects: t = ' num2str(s.t(k),'%.2g') ' \mus = ' num2str(s.t(k)/s.tau,'%.2g') ' \tau_{exp}'];

    
    savename = [savedir f 'fig' num2str(figiter) '-velocity-tEvol' num2str(k) '.png'];
    saveas(fig,savename)
    close(fig)
end
end