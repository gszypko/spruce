function [] = gaussianAnalysis(data,flags)
%% handle inputs
f = filesep;

s = struct;
s.t = data.grids.time;
s.x = data.grids.x_vec;
s.y = data.grids.y_vec;
s.Te0 = data.Te;
s.Ti0 = data.Ti;
s.T = data.Te + data.Ti;
s.sig0 = data.sig0;
s.tau = data.tau;

%% Fit Density Distribution with 2D Gaussian to Extract RMS Width
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


% do gaussian fits
for i = 1:length(data.grids.time)    
    x = data.grids.x_vec;
    y = data.grids.y_vec;
    if strcmp(data.config.eq_set,'ideal_2F')
        img = data.grids.vars(i).i_n;
    else
        img = data.grids.vars(i).n;
    end
    [fit(i)] = fitImgWithGaussian(x,y,img,0);
    disp(['2D Gaussian Fits: ' num2str(i) '/' num2str(length(data.grids.time))])
end

% plot gaussian fits
fignum = 0;
frames = cell(1,length(data.grids.time));
for i = 1:length(data.grids.time)    
    disp(['Plotting 2D Gaussian Fits: ' num2str(i) '/' num2str(length(data.grids.time))])
    
    fignum = fignum + 1;
    hold([ax{:}],'off')
    [~,indx] = min(abs(fit(i).x));
    [~,indy] = min(abs(fit(i).y));
    imgx = fit(i).img(indy,:);
    imgy = fit(i).img(:,indx)';
    imgfitx = fit(i).imgfit(indy,:);
    imgfity = fit(i).imgfit(:,indx)';
    
    lgdstr = {'Simulation','Fit'};
    l = get_line_specs(length(lgdstr));

    cax = get_axis(fig,ax{1});
    if data.grids.is_uniform
        zdata = fit(i).img;
    else
        zdata = data.grids.uni_grid(fit(i).img);
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
        zdata = fit(i).imgfit;
    else
        zdata = data.grids.uni_grid(fit(i).imgfit);
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
        zdata = fit(i).imgres;
    else
        zdata = data.grids.uni_grid(fit(i).imgres);
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
    xlabel('y (cm)')

    str1 = ['t = ' num2str(data.grids.time(i)*1e6,'%.2g') ' \mus = ' num2str(data.grids.time(i)/s.tau,'%.2g') ' \tau'];
    str2 = ['n_0 = ' num2str(fit(i).amp/1e8,'%.2g') '\times10^8 cm^{-3}'];
    str3 = ['\sigma_x = ' num2str(fit(i).sigx*10,'%.2g') ' \pm ' num2str(fit(i).sigxErr*10,'%.2g') ' mm'];
    str4 = ['\sigma_y = ' num2str(fit(i).sigy*10,'%.2g') ' \pm ' num2str(fit(i).sigyErr*10,'%.2g') ' mm'];
    an.String = {str1,str2,str3,str4};
        
    frames{i} = getframe(fig);
end
close(fig)
write_video([data.folder f 'gauss-fits'],frames);

%% Compile MHD Data and Vlasov Theory for Comparison

s.data = struct;
for i = 1:length(s.t)
    s.data(i).sigx = fit(i).sigx;
    s.data(i).sigy = fit(i).sigy;
    s.data(i).sig = sqrt(fit(i).sigx*fit(i).sigy);
    s.data(i).n2 = fit(i).amp;
    s.data(i).n3 = fit(i).amp;
    [~,indx] = min(abs(data.grids.x_vec));
    [~,indy] = min(abs(data.grids.y_vec));
    if max(strcmp(data.config.eq_set,{'ideal_mhd','ideal_mhd_cons'}))
        s.data(i).Ti = data.grids.vars(i).temp(indy,indx);
        s.data(i).Te = data.grids.vars(i).temp(indy,indx);
        s.data(i).xt.n = data.grids.vars(i).n(indy,:);
        s.data(i).yt.n = data.grids.vars(i).n(:,indx)';
        s.data(i).xt.v = data.grids.vars(i).v_x(indy,:);
        s.data(i).yt.v = data.grids.vars(i).v_y(:,indx)';
        s.data(i).xt.Te = data.grids.vars(i).temp(indy,:);
        s.data(i).yt.Te = data.grids.vars(i).temp(:,indx)';
    elseif strcmp(data.config.eq_set,{'ideal_2F'})
        s.data(i).Ti = data.grids.vars(i).i_temp(indy,indx);
        s.data(i).Te = data.grids.vars(i).e_temp(indy,indx);
        s.data(i).xt.n = data.grids.vars(i).i_n(indy,:);
        s.data(i).yt.n = data.grids.vars(i).i_n(:,indx)';
        s.data(i).xt.v = data.grids.vars(i).i_v_x(indy,:);
        s.data(i).yt.v = data.grids.vars(i).i_v_x(:,indx)';
        s.data(i).xt.Te = data.grids.vars(i).e_temp(indy,:);
        s.data(i).yt.Te = data.grids.vars(i).e_temp(:,indx)';
    elseif strcmp(data.config.eq_set,{'ideal_mhd_2E'})
        s.data(i).Ti = data.grids.vars(i).i_temp(indy,indx);
        s.data(i).Te = data.grids.vars(i).e_temp(indy,indx);
        s.data(i).xt.n = data.grids.vars(i).n(indy,:);
        s.data(i).yt.n = data.grids.vars(i).n(:,indx)';
        s.data(i).xt.v = data.grids.vars(i).v_x(indy,:);
        s.data(i).yt.v = data.grids.vars(i).v_y(:,indx)';
        s.data(i).xt.Te = data.grids.vars(i).e_temp(indy,:);
        s.data(i).yt.Te = data.grids.vars(i).e_temp(:,indx)';
    else
        error('Error: equation set is not valid.')
    end
    
end

s.theory = struct;
[~,sig,gam,Ti,Te,n2,n3] = kinetic_model(s.t,s.sig0,s.Ti0,s.Te0,s.data(1).n2,data.config.eic_opt);
for i = 1:length(s.t)
    s.theory(i).sigx = sig(i);
    s.theory(i).sigy = sig(i);
    s.theory(i).sig = sig(i);
    s.theory(i).Ti = Ti(i);
    s.theory(i).Te = Te(i);
    s.theory(i).n2 = n2(i);
    s.theory(i).n3 = n3(i);
    s.theory(i).gam = gam(i);
    s.theory(i).vx = gam(i).*s.x;
    s.theory(i).vy = gam(i).*s.y;
end

%% Fit Velocity Transects for Effective Temperature

% define functions for velocity along a given axis
tau = @(sig,T) getTauExp(sig,T);
sig = @(t,sig0,T) sig0.*sqrt(1+t.^2./tau(sig0,T).^2);
gam = @(t,sig,T) t./tau(sig,T).^2/(1+t.^2./tau(sig,T).^2);
v = @(r,t,sig,T) r.*gam(t,sig,T);

% for each time point, fit velocity along x and y axis with model to extract T
for i = 1:length(s.data)
    
    % do x axis fitting
    fun_xt = @(x,data) v(data,s.t(i),s.data(1).sigx,x);
    x0 = s.T;
    ind_x = abs(s.x) < s.data(1).sigx*3;
    xdata = s.x(ind_x);
    ydata = s.data(i).xt.v(ind_x);
    lb = 0; ub = 1e3;
    x = lsqcurvefit(fun_xt,x0,xdata,ydata,lb,ub);
    s.data(i).xt.v_fit = fun_xt(x,s.x);
    s.data(i).T_x = x;
    s.data(i).gam_x = gam(s.t(i),s.data(1).sigx,x);
    s.data(i).gam_x_norm = s.data(i).gam_x/gam(s.t(i),s.data(1).sigx,data.Te+data.Ti);

    % do y axis fitting
    fun_yt = @(x,data) v(data,s.t(i),s.data(1).sigy,x);
    x0 = s.T;
    ind_y = abs(s.y) < s.data(1).sigy*3;
    xdata = s.y(ind_y);
    ydata = s.data(i).yt.v(ind_y);
    lb = 0; ub = 1e3;
    x = lsqcurvefit(fun_yt,x0,xdata,ydata,lb,ub);
    s.data(i).yt.v_fit = fun_yt(x,s.y);
    s.data(i).T_y = x;
    s.data(i).gam_y = gam(s.t(i),s.data(1).sigy,x);
    s.data(i).gam_y_norm = s.data(i).gam_y/gam(s.t(i),s.data(1).sigy,data.Te+data.Ti);

end

%% Fit Size Evolution for Effective Temperature

% fit for temperature on x axis
fit_model_x = @(c,data) sig(data,s.data(1).sigx,c);
xdata = [s.t];
ydata = [s.data.sigx];
T_x = lsqcurvefit(fit_model_x,s.T,xdata,ydata,0,1e3);
s.T_x = T_x;

% fit for temperature on y axis
fit_model_y = @(c,data) sig(data,s.data(1).sigy,c);
xdata = [s.t];
ydata = [s.data.sigy];
T_y = lsqcurvefit(fit_model_y,s.T,xdata,ydata,0,1e3);
s.T_y = T_y;

% fit for rms temperature
fit_model_2D = @(c,data) sig(data,s.data(1).sig,c);
xdata = [s.t];
ydata = [s.data.sig];
T_2D = lsqcurvefit(fit_model_2D,s.T,xdata,ydata,0,1e3);
s.T_2D = T_2D;

% store fits in data structure
sigx_fit = fit_model_x(s.T_x,s.t);
sigy_fit = fit_model_y(s.T_y,s.t);
sig_fit = fit_model_2D(s.T_2D,s.t);
for i = 1:length(s.t)
    s.data(i).sigx_fit = sigx_fit(i);
    s.data(i).sigy_fit = sigy_fit(i);
    s.data(i).sig_fit = sig_fit(i);
end

%% Plot Fits to Plasma Size

p.x.data = [s.data.sigx];
p.x.fit = [s.data.sigx_fit];
p.y.data = [s.data.sigy];
p.y.fit = [s.data.sigy_fit];
p.gm.data = [s.data.sig];
p.gm.fit = [s.data.sig_fit];
fields = fieldnames(p);

colvar = {'x','y','gm'};
colstr = {'\sigma_x (cm)','\sigma_y (cm)','\sigma (cm)'};
vars = {'data','fit'};
title_var = {'T_x','T_y','T_2D'};

num = length(fields);
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
iter = 0;
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        cax = get_axis(fig,ax{i,j});
        hold on

        l = get_line_specs(length(vars));
        for k = 1:length(vars)
            xdata = [s.t].*1e6;
            ydata = [p.(colvar{iter}).(vars{k})];
            lp = l(k);
            plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
        end
        cax.PlotBoxAspectRatio = [1 1 1];
        cax.FontSize = 12;
        grid on
        grid minor

        if i == row, xlabel('t (\mus)'), end
        ylabel(colstr{iter})
        title(['T_e = ' num2str(s.(title_var{iter})) ' K'],'FontWeight','normal')
    end
end
lgd = legend(vars);
lgd.Position = [0.814763192769338,0.212150191817874,0.084684301063469,0.152255642324462];
saveas(fig,[data.folder f 'sig-fits.png']);
close(fig)

%% Plot Velocity Transects with Fits

q1 = {'xt','yt'};
q2 = {'x','y'};
q3 = {'x (cm)','y (cm)'};
vars = {'v','v_fit'};
vars_str = {'data','fit'};
title_var = {'gam_x_norm','gam_y_norm'};
title_str = {'\gamma_x','\gamma_y'};

frames = cell(length([s.t]));
num = length(q1);
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.250927754415475,0.924355890669176,0.516400000000000,0.102000000000000];

for m = 1:length(s.t)
    iter = 0;
    for i = 1:row
        for j = 1:col
            if iter > num - 1, break, end
            iter = iter + 1;
            cax = get_axis(fig,ax{i,j});
            hold off
    
            l = get_line_specs(length(vars));
            for k = 1:length(vars)
                xdata = [s.(q2{iter})].*1e6;
                ydata = [s.data(m).(q1{iter}).(vars{k})];
                lp = l(k);
                plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
                hold on
            end
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 11;
            grid on
            grid minor
    
            if i == row, xlabel(q3{iter}), end
            ylabel('v (cm/s)')
            title([title_str{iter} ' = ' num2str(s.data(m).(title_var{iter})) ' \gamma'],'FontWeight','normal')
        end
    end
    lgd = legend(vars_str);
    lgd.Position = [0.814763192769338,0.212150191817874,0.084684301063469,0.152255642324462];

    an.String = ['t = ' num2str(s.t(m)*1e6) ' \mus'];

    frames{m} = getframe(fig);
end
close(fig)
write_video([data.folder f 'velocity-fits'],frames);

%% Plot Summary of Expansion for MHD Sims Against Vlasov Theory


q = {'data','theory'};
qstr{1} = 'Simulation';
if data.config.eic_opt, qstr{2} = 'Vlasov';
else, qstr{2} = 'Vlasov w/EIC'; end
colvar = {'sigx','sigy','sig','n2','Ti','Te'};
colstr = {'\sigma_x (cm)','\sigma_y (cm)','\sigma (cm)','n_0 (10^8 cm^-^3)','T_i (K)','T_e (K)'};
num = length(colvar);
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];
iter = 0;
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        cax = get_axis(fig,ax{i,j});
        hold on

        lgdstr = qstr;
        l = get_line_specs(length(lgdstr));
        
        for k = 1:length(q)
            xdata = [s.t]./s.tau;
            ydata = [s.(q{k}).(colvar{iter})];
            lp = l(k);
            plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
        end
        cax.PlotBoxAspectRatio = [1 1 1];
        cax.FontSize = 12;
        grid on
        grid minor

        if i == row, xlabel('t / \tau'), end
        ylabel(colstr{iter})
    end
end

% add legend
lgd = legend(lgdstr);
lgd.Position = [0.8954    0.4610    0.0911    0.0787];

saveas(fig,[data.folder f 'vlasov-evol.png']);
close(fig)

%% Plot Central Transects for Density, Velocity, and Temperature

% generate figure
frames = cell(length([s.t]));
num = 6;
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];

col_var = {'n','v','Te'};
row_var_y = {'xt','yt'};
row_var_x = {'x','y'};
col_str = {'n (cm^-3)','v (cm/s)','T_e (K)'};
row_str = {'x (cm)','y (cm)'};

for k = 1:length(s.t)
    disp(['Plotting Transects: ' num2str(k) '/' num2str(length(s.t))])
    hold([ax{:}],'off')

    iter = 0;
    for i = 1:row
        for j = 1:col
            if iter > num - 1, break, end
            iter = iter + 1;
            cax = get_axis(fig,ax{i,j});
            hold off
            
            l = get_line_specs(2);
            for indq = 1:length(q)
                xdata = [s.(row_var_x{i})];
                ydata = [s.data(k).(row_var_y{i}).(col_var{j})];
                lp = l(1);
                plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
            end
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 10;
            grid on
            grid minor

            xlabel(row_str{i})
            ylabel(col_str{j})

        end
    end


    an.String = ['Central Transects: t = ' num2str(s.t(k)*1e6,'%.2g') ' \mus = ' num2str(s.t(k)/s.tau,'%.2g') ' \tau_{exp}'];
    frames{k} = getframe(fig);
end
close(fig)
write_video([data.folder f 'transects'],frames);

end