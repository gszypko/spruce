function [] = gaussianAnalysis(data)

% set up save directories
f = filesep;

s = struct;
s.t = data.grids.time;
s.x = data.grids.x_vec;
s.y = data.grids.y_vec;
s.Te0 = data.settings.Te;
s.Ti0 = data.settings.Ti;
s.sig0 = data.sig0;
s.tau = data.tau;

%% Fit Density Distribution with 2D Gaussian to Extract RMS Width

file_path = [data.folder f 'gauss-fits'];
vid = VideoWriter(file_path,'MPEG-4');
vid.FrameRate = 10;
open(vid)

fig = figure('Visible','off');
fig.Position = [195         259        1027         491];
fig.Color = [1 1 1];

row = 2;
col = 3;
num = 5;

ax = cell(row,col);
iter = 0;
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        ax{j,i} = subplot(row,col,iter);
        if i == row, ax{j,i}.Position(1) = ax{j,i}.Position(1) - .025; end
    end
end

an = annotation('textbox');
an.Position = [0.6654    0.1141    0.2527    0.3299];
an.HorizontalAlignment = 'left';
an.VerticalAlignment = 'top';
an.LineStyle = 'none';
an.FontSize = 12;
    
for i = 1:length(data.grids.time)
    hold([ax{:}],'off')
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
    l = get_line_specs(length(lgdstr));

    set(fig,'CurrentAxes',ax{1})
    imagesc(x,y,img./1e8)
    ax{1}.PlotBoxAspectRatio = [1 1 1];
    shading interp
    cb = colorbar;
    
    set(fig,'CurrentAxes',ax{2})
    imagesc(x,y,fit(i).imgfit./1e8);
    ax{2}.PlotBoxAspectRatio = [1 1 1];
    shading interp
    cb = colorbar;
    
    set(fig,'CurrentAxes',ax{3})
    imagesc(x,y,fit(i).imgres./1e8);
    ax{3}.PlotBoxAspectRatio = [1 1 1];
    shading interp
    cb = colorbar;
    cb.Label.String = 'n (10^8 cm^-^3)';
    
    set(fig,'CurrentAxes',ax{4})
    plot(x,imgx./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
    hold on
    plot(x,imgfitx./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
    ax{4}.PlotBoxAspectRatio = [1 1 1];
    title('x axis','FontWeight','normal')
    ylabel('n (10^8 cm^-^3)')
    
    set(fig,'CurrentAxes',ax{5})
    plot(y,imgy./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
    hold on
    plot(y,imgfity./1e8,'LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
    ax{5}.PlotBoxAspectRatio = [1 1 1];
    title('y axis','FontWeight','normal')
    legend(lgdstr)
    
    str1 = ['t = ' num2str(data.grids.time(i)*1e6,'%.2g') ' \mus = ' num2str(data.grids.time(i)/s.tau,'%.2g') ' \tau_{exp}'];
    str2 = ['n_{max} = ' num2str(fit(i).amp/1e8,'%.2g') '\times10^8 cm^{-3}'];
    str3 = ['\sigma_x = ' num2str(fit(i).sigx*10,'%.2g') ' \pm ' num2str(fit(i).sigxErr*10,'%.2g') ' mm'];
    str4 = ['\sigma_y = ' num2str(fit(i).sigy*10,'%.2g') ' \pm ' num2str(fit(i).sigyErr*10,'%.2g') ' mm'];
    an.String = {str1,str2,str3,str4};
        
    frame = getframe(fig);
    writeVideo(vid,frame);
end
close(fig)

%% Compile MHD Data and Vlasov Theory for Comparison

s.data = struct;
for i = 1:length(s.t)
    s.data(i).sigx = fit(i).sigx;
    s.data(i).sigy = fit(i).sigy;
    s.data(i).n = fit(i).amp;
    [~,indx] = min(abs(data.grids.x_vec));
    [~,indy] = min(abs(data.grids.y_vec));
    s.data(i).Ti = data.grids.vars(i).temp_i(indy,indx);
    s.data(i).Te = data.grids.vars(i).temp_e(indy,indx);
    s.data(i).vx = data.grids.vars(i).v_x(indy,:);
    s.data(i).vy = data.grids.vars(i).v_y(:,indx)';
end

s.theory = struct;
[sig,Ti,Te,v] = getUCNPExpansion(s.t,s.sig0,s.Ti0,s.Te0,s.x,s.y);
for i = 1:length(s.t)
    s.theory(i).sigx = sig(i);
    s.theory(i).sigy = sig(i);
    s.theory(i).Ti = Ti(i);
    s.theory(i).Te = Te(i);
    s.theory(i).vx = v(i).vx;
    s.theory(i).vy = v(i).vy;
end

%% Plot Summary of Expansion for MHD Sims Against Vlasov Theory

fig = figure('Visible','on');
fig.Position = [195         259        1027         491];
fig.Color = [1 1 1];

q = {'data','theory'};
qstr = {'MHD','Vlasov'};
colvar = {'sigx','sigy','Ti','Te'};
colstr = {'\sigma_x (cm)','\sigma_y (cm)','T_i (K)','T_e (K)'};
rowvar = {''};
row = 2;
col = 2;
num = length(colvar);

ax = cell(row,col);
iter = 0;
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        ax{i,j} = subplot(row,col,iter);
        hold on
        lgdstr = qstr;
        l = get_line_specs(length(lgdstr));
        
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
lgd.Position = [0.777965695280062,0.604462568894955,0.206112384107097,0.086259707997385];

frame = getframe(fig);
writeVideo(vid,frame);
close(fig)

%% Plot Velocity Transects

fig = figure('Visible','off');
fig.Position = [195         259        1027         491];
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
        if i == row, ax{i,j}.Position(1) = ax{i,j}.Position(1) - .025; end
    end
end

an = annotation('textbox');
    an.Position = [0.0058    0.8762    0.9919    0.1143];
    an.HorizontalAlignment = 'center';
    an.VerticalAlignment = 'middle';
    an.LineStyle = 'none';
    an.FontSize = 12;
for k = 1:length(s.t)
    disp(['Plotting velocity transects: ' num2str(k) '/' num2str(length(s.t))])
    hold([ax{:}],'off')
    
    

    iter = 0;
    for i = 1:row
        for j = 1:col
            if iter > num - 1, break, end
            iter = iter + 1;
            set(fig,'CurrentAxes',ax{i,j})
            
            lgdstr = qstr;
            l = get_line_specs(length(lgdstr));
            
            for indq = 1:length(q)
                xdata = [s.(colvar1{j})];
                ydata = [s.(q{indq})(k).(colvar2{j})];
                lp = l(indq);
                plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
                hold on
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

    
    an.String = ['Central Velocity Transects: t = ' num2str(s.t(k)*1e6,'%.2g') ' \mus = ' num2str(s.t(k)/s.tau,'%.2g') ' \tau_{exp}'];

    
    frame = getframe(fig);
    writeVideo(vid,frame);
    
end

close(fig)
end