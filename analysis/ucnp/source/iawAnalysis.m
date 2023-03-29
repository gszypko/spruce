function [] = iawAnalysis(data,flags)
%% Extract 1D Density Transects
% For each time point, record the density distribution along the x axis
s = struct();
s(length(data.grids.time)).t = [];
for i = 1:length([data.grids.time])
    s(i).t = data.grids.time(i);
    s(i).x = data.grids.x_vec;
    indy = abs(data.grids.y_vec) < data.grids.gauss_fits(i).sigy/2;
    s(i).n = mean(data.grids.vars(i).n(indy,:),1);
    s(i).time_fac = sqrt(1+s(i).t^2/data.tau^2);
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
    x0(3) = data.sigx; lb(3) = data.sigx/2; ub(3) = data.sigx*10;
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
    L = data.settings.n_iaw_sig*s(i).time_fac;
    k = 2*pi/L;
    xdata = s(i).x;
    ydata = s(i).dntild;
    fun = @(c,d) G(d,1,0,c(3)).*(-c(1).*cos(c(2).*d)+c(4).*cos(2*c(2).*d));
    dy = max(ydata)-min(ydata);
    x0 = []; lb = []; ub = [];
    x0(1) = dy/2; lb(1) = -5*x0(1); ub(1) = 5*x0(1);
    x0(4) = dy/10; lb(4) = -x0(2); ub(4) = x0(2);
    x0(3) = s(i).sig; lb(3) = x0(3).*.85; ub(3) = x0(3)*1.15;
    x0(2) = k; lb(2) = x0(2).*0.75; ub(2) = x0(2)*1.25;
%     x0(5) = pi; lb(5) = 0; ub(5) = 2*pi;
%     x0(3) = k; lb(3) = x0(3)*.1; ub(3) = x0(3)*10;
%     x0(3) = pi/2; lb(3) = 0; ub(3) = 2*pi;
%     x0(4) = s(i).sig; lb(4) = x0(4)/5; ub(4) = x0(4)*5;
    x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
    s(i).dntild_fit = fun(x,xdata);
    s(i).dntild_guess = fun(x0,xdata);
    s(i).A = x(1);
    s(i).k = k;
    s(i).lam = 2*pi/s(i).k;
%     s(i).sig = x(3);
    
%     ydata_fft = fft(ydata);
%     xdata_fft = (0:length(ydata_fft)-1)*(2*pi/mean(diff(xdata)))/length(ydata_fft);
%     plot(xdata_fft,real(ydata_fft).^2);
%     hold on
%     plot(xdata_fft,imag(ydata_fft).^2);
%     hold off
%     xlim([0,100])
%     shg
%     pause(0.1)

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
q = {'sig','lam','time_fac','A'};
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
        for i1 = 1:length(t)
            phi(i1) = temp(t(i1));
        end
        amp = A0.*exp(-g0.*t).*cos(phi);
    end
fun = @(c,d) amp_fun(d,c(1),c(2),c(3),data.tau);
x0 = []; lb = []; ub = [];
x0(1) = s(1).A; lb(1) = -x0(1)*5; ub(1) = x0(1)*5;
x0(2) = 0.2/data.tau; lb(2) = 0; ub(2) = 2*x0(2);
cs = getSoundSpeed(data.settings.m_i,data.settings.Ti,data.settings.Te,data.adiabatic_index_e,0);
x0(3) = s(1).k*cs; lb(3) = x0(3)/10; ub(3) = x0(3)*10;
xdata = [s.t];
ydata = [s.A];

num = 1;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.065026362038664,0.944444444444444,0.929701230228474,0.047355555555555];
cax = get_axis(fig,ax{1});
hold on
l = get_line_specs(2);     
lp = l(1);
plot(xdata.*1e6,ydata,'o','LineWidth',2,'MarkerSize',2,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
lp = l(2);
plot(xdata.*1e6,fun(x,xdata),'-','LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
legend({'data','fit'})
xlabel('t (\mus)')
ylabel('Amplitude, A(t)')
str1 = ['\omega_0 (fit) = ' num2str(x(3)) ' s^-^1'];
str2 = ['\omega_0 (theory) = ' num2str(x0(3)) ' s^-^1'];
an.String = [str1 ' : ' str2];
saveas(fig,[data.folder filesep 'iaw-amp-fit.png'])
close(fig)

end