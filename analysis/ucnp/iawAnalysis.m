function [] = iawAnalysis(data,flags)
%% Extract 1D Density Transects
% For each time point, record the density distribution along the x axis
s = struct();
s(length(data.grids.time)).t = [];
for i = 1:length([data.grids.time])
    s(i).t = data.grids.time(i);
    s(i).x = data.grids.x_vec;
    [~,indx] = min(abs(s(i).x));
    s(i).n = data.grids.vars(i).n(indx,:);
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
    x0(3) = dx/5; lb(3) = 0; ub(3) = 2*dx;
    x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
    s(i).ng = fun(x,xdata);
    s(i).n0 = x(1);
    s(i).sig = x(3);
    s(i).dn = s(i).n - s(i).ng;
    s(i).dntild = s(i).dn./s(i).n0;
end

% Plot fit results
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
    fun = @(c,d) c(1).*G(d,1,0,c(2)).*cos(c(3).*d+c(4));
    dy = max(ydata)-min(ydata);
    x0 = []; lb = []; ub = [];
    x0(1) = dy/2; lb(1) = 0; ub(1) = 5*x0(1);
    x0(2) = s(i).sig; lb(2) = x0(2)/2; ub(2) = x0(2)*2;
    L = s(i).sig/2;
    x0(3) = 2*pi/L; lb(3) = x0(3)/5; ub(3) = x0(3)*5;
    x0(4) = 0; lb(4) = 0; ub(4) = 2*pi;
    x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
    s(i).dntild_fit = fun(x,xdata);
    s(i).dntild_guess = fun(x0,xdata);
    s(i).A = x(1);
    s(i).sig2 = x(2);
    s(i).k = x(3);
    s(i).lam = 2*pi/s(i).k;
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

%% Plot Fitted Wavelength and Plasma Size

end