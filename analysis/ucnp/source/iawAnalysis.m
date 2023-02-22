function [] = iawAnalysis(data,flags)
%% Extract 1D Density Transects
% For each time point, record the density distribution along the x axis
s = struct();
s(length(data.grids.time)).t = [];
for i = 1:length([data.grids.time])
    s(i).t = data.grids.time(i);
    s(i).x = data.grids.x_vec;
    [~,indy] = min(abs(data.grids.y_vec));
    s(i).n = data.grids.vars(i).n(indy,:);
    s(i).gam = sqrt(1+s(i).t^2/data.tau^2);
end

%% Fit 1D Transects with Gaussian
% Fit 1D density transects with a Gaussian to extract size of overall plasma distribution
G = @(x,amp,x0,sig) amp.*exp(-(x-x0).^2./(2*sig^2));
for i = 1:length([s.t])
    % c = [amp x0 sig]
    xdata = s(i).x;
    ydata = s(i).n;
    fun = @(c,d) G(d,c(1),c(2),c(3));
    x0(1) = max(ydata); lb(1) = 0; ub(1) = 2*x0(1);
    x0(2) = 0; lb(2) = min(xdata); ub(2) = max(xdata);
    dx = max(xdata) - min(xdata);
    x0(3) = data.sig0; lb(3) = 0; ub(3) = 2*dx;
    if strcmp(data.settings.n_dist,'uniform')
        x0(3) = 1000; lb(3) = x0(3)/2; ub(3) = x0(3)*2;
    end
    x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
    s(i).ng = fun(x,xdata);
    s(i).n0 = x(1);
    s(i).sig = x(3);
    s(i).dn = s(i).n - s(i).ng;
    s(i).dntild = s(i).dn./s(i).n0;
end

% Plot Gaussian fit results
num = 1;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1595    0.9117    0.7230    0.0801];
frames = cell(1,length([s.t]));
for i = 1:length([s.t])
    cax = get_axis(fig,ax{1});
    q = {'n','ng'};
    qstr = {'Data','Fit'};
    lgdstr = qstr;
    l = get_line_specs(length(lgdstr));
    hold off
    for k = 1:length(q)
        xdata = [s(i).x];
        ydata = [s(i).(q{k})];
        lp = l(k);
        plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
        hold on
    end
    cax.PlotBoxAspectRatio = [1 1 1];
    cax.FontSize = 12;
    grid on
    grid minor
    
    xlabel('x (cm)')
    ylabel('n (cm^-^3)')
    
    lgd = legend(lgdstr);
    lgd.Position = [0.6815    0.8083    0.1797    0.1004];
    
    dlm = ' - ';
    str1 = ['Iter = ' num2str(i-1)];
    str2 = ['t = ' num2str(s(i).t*1e6,'%.3g') '\mus = ' num2str(s(i).t/data.tau,'%.3g') '\tau_e_x_p'];
    an.String = [str1 dlm str2];

    frames{i} = getframe(fig);
end
close(fig)
filepath = [data.folder filesep 'iaw-gauss-fits'];
write_video(filepath,frames);

%% Fit Density Perturbation to Extract Amplitude and Wavenumber
G = @(x,amp,x0,sig) amp.*exp(-(x-x0).^2./(2*sig^2));
for i = 1:length([s.t])
    % c = [A sig k phi]
    xdata = s(i).x;
    ydata = s(i).dntild;
    fun = @(c,d) G(d,1,0,s(i).sig).*(c(1).*cos(c(2).*d+c(3)));
    dy = max(ydata)-min(ydata);
    x0 = []; lb = []; ub = [];
    x0(1) = dy/2; lb(1) = 0; ub(1) = 5*x0(1);
    L = data.settings.n_iaw_sig*s(i).gam;
    x0(2) = 2*pi/L; lb(2) = x0(2)/5; ub(2) = x0(2)*5;
    x0(3) = pi/2; lb(3) = 0; ub(3) = 2*pi;
    x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
    s(i).dntild_fit = fun(x,xdata);
    s(i).dntild_guess = fun(x0,xdata);
    s(i).A = x(1);
    s(i).k = x(2);
    s(i).lam = 2*pi/s(i).k;
    s(i).phi = x(3);
end

% Plot Density Perturbation
num = 1;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1595    0.9117    0.7230    0.0801];
frames = cell(1,length([s.t]));
for i = 1:length([s.t])
    cax = get_axis(fig,ax{1});
    q = {'dntild','dntild_fit'};
    qstr = {'Data','Fit'};
    lgdstr = qstr;
    l = get_line_specs(length(lgdstr));
    hold off
    for k = 1:length(q)
        xdata = [s(i).x];
        ydata = [s(i).(q{k})];
        lp = l(k);
        plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
        hold on
    end
    cax.PlotBoxAspectRatio = [1 1 1];
    cax.FontSize = 12;
    grid on
    grid minor
    
    xlabel('x (cm)')
    ylabel('$\delta\tilde{n}$','Interpreter','latex')
    
    lgd = legend(lgdstr,'Interpreter','latex');
    lgd.Position = [0.6815    0.8083    0.1797    0.1004];
    
    dlm = ' - ';
    str1 = ['Iter = ' num2str(i-1)];
    str2 = ['t = ' num2str(s(i).t*1e6,'%.3g') '\mus = ' num2str(s(i).t/data.tau,'%.3g') '\tau_e_x_p'];
    an.String = [str1 dlm str2];

    frames{i} = getframe(fig);
end
close(fig)
filepath = [data.folder filesep 'iaw-dntild'];
write_video(filepath,frames);

%% Plot Lengthscale Evolution - Fitted Wavelength and Plasma Size

% Plot size evolution
num = 1;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1595    0.9117    0.7230    0.0801];
cax = get_axis(fig,ax{1});
q = {'sig','lam','gam','A'};
qstr = {'$\sigma / \sigma_0$','$\lambda / \lambda_0$','$\sqrt{1+t^2/\tau^2}$','$A / A_0$'};
lgdstr = qstr;
l = get_line_specs(length(lgdstr));
hold off
for k = 1:length(q)
    xdata = [s.t];
    ydata = [s.(q{k})];
    ydata = ydata./ydata(1);
    lp = l(k);
    plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
    hold on
end
cax.PlotBoxAspectRatio = [1 1 1];
cax.FontSize = 12;
grid on
grid minor

xlabel('x (cm)')
ylabel('Length Scale Evolution','Interpreter','latex')

lgd = legend(lgdstr,'Interpreter','latex');
lgd.Position = [0.1821    0.7597    0.2681    0.1444];

dlm = ' - ';
str1 = ['Iter = ' num2str(i-1)];
str2 = ['t = ' num2str(s(i).t*1e6,'%.3g') '\mus = ' num2str(s(i).t/data.tau,'%.3g') '\tau_e_x_p'];
an.String = [str1 dlm str2];

filepath = [data.folder filesep 'iaw-length-evol.png'];
saveas(fig,filepath);
close(fig)

%% Fit for Dispersion Relation
    function [amp] = amp_fun(t,A0,g0,w0,tau)
        temp = @(t) w0.*integral(@(tp) 1./(1+tp.^2./tau^2),0,t);
        phi = zeros(size(t));
        for i = 1:length(t)
            phi(i) = temp(t(i));
        end
        amp = A0.*exp(-g0.*t).*abs(cos(phi));
    end
fun = @(c,d) amp_fun(d,c(1),c(2),c(3),data.tau);
x0 = []; lb = []; ub = [];
x0(1) = s(1).A; lb(1) = x0(1)/2; ub(1) = x0(1)*2;
x0(2) = 0.2/data.tau; lb(2) = 0; ub(2) = 2*x0(2);
cs = getSoundSpeed(data.settings.m_i,data.settings.Ti,data.settings.Te,data.adiabatic_index_e,data.adiabatic_index_i);
x0(3) = s(1).k*cs/2; lb(3) = x0(3)/10; ub(3) = x0(3)*10;
xdata = [s.t];
ydata = [s.A];
num = 1;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1595    0.9117    0.7230    0.0801];
cax = get_axis(fig,ax{1});
hold on
plot(xdata,ydata)
plot(xdata,fun(x0,xdata))
x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
plot(xdata,fun(x,xdata))
legend({'data','guess','fit'})

end