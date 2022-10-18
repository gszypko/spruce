function [] = compareExpAndSimData(data,flags)

% load experimental data from 'os.mat' file
f = filesep;
base_folder = extractBefore(data.folder,[f 'set_']);
load([base_folder f 'os.mat'],'os');

% find indices of simulation time points that are closest to experimental time points
% ind_t(i) holds the index of data.grids.vars corresponding to experimental time point os(i).delays
ind_t = zeros(size(os)); 
for i = 1:length([os.delays])
    [~,ind_t(i)] = min(abs(os(i).delays*1e-9 - data.grids.time));
end

% create struct to hold info to be plotted
imgs = struct();
imgs(length([os.delays])).t = [];
interpolate_sims = false;
for i = 1:length([os.delays])
    imgs(i).t = os(i).delays*1e-9;
    imgs(i).time_fac = sqrt(1+imgs(i).t^2/data.tau^2);
    if max(data.grids.x_vec) > max(os(i).imgs.xRelInMM/10), interpolate_sims = true; end
    if interpolate_sims
        imgs(i).x = os(i).imgs.xRelInMM/10;
        imgs(i).y = os(i).imgs.yRelInMM/10;
        [imgs(i).X, imgs(i).Y] = meshgrid(imgs(i).x, imgs(i).y);
        imgs(i).n_exp = os(i).imgs.density*1e8;
        imgs(i).n_sim = interp2(data.grids.pos_x,data.grids.pos_y,data.grids.vars(ind_t(i)).n,imgs(i).X,imgs(i).Y);
    elseif ~interpolate_sims
        imgs(i).x = data.grids.x_vec;
        imgs(i).y = data.grids.y_vec;
        [imgs(i).X, imgs(i).Y] = meshgrid(imgs(i).x, imgs(i).y);
        imgs(i).n_exp = interp2(os(i).imgs.xRelInMM/10,os(i).imgs.yRelInMM/10,os(i).imgs.density,imgs(i).X,imgs(i).Y)*1e8;
        imgs(i).n_sim = data.grids.vars(ind_t(i)).n;
    end
    imgs(i).n_res = imgs(i).n_sim - imgs(i).n_exp;
    
    fit_exp = fitImgWithGaussian(imgs(i).x,imgs(i).y,imgs(i).n_exp);
    fit_sim = fitImgWithGaussian(imgs(i).x,imgs(i).y,imgs(i).n_sim);
    imgs(i).sigx_exp = fit_exp.sigx;
    imgs(i).sigx_sim = fit_sim.sigx;
    imgs(i).sigy_exp = fit_exp.sigy;
    imgs(i).sigy_sim = fit_sim.sigy;
    fac_exp = imgs(1).sigx_exp*imgs(1).sigy_exp^2/(imgs(i).sigx_exp*imgs(i).sigy_exp^2);
    fac_sim = imgs(1).sigx_sim*imgs(1).sigy_sim/(imgs(i).sigx_sim*imgs(i).sigy_sim);
    imgs(i).n_exp_norm = imgs(i).n_exp./fac_exp;
    imgs(i).n_sim_norm = imgs(i).n_sim./fac_sim;
    imgs(i).n_res_norm = imgs(i).n_sim_norm - imgs(i).n_exp_norm;
end

% plot density with residuals
colvar = {'n_sim','n_exp','n_res'};
colstr = {'Sim','Exp','Sim - Exp'};
num = length(colvar);
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1595    0.9084    0.7230    0.0801];

for k = 1:length([imgs.t])
    disp(['Plotting: ' num2str(k) '/' num2str(length([imgs.t]))])
    
    % update axis information
    iter = 0;
    for i = 1:size(ax,1)
        for j = 1:size(ax,2)
            iter = iter + 1;
            if iter > num, break, end
            
            cax = get_axis(fig,ax{i,j});
            xdata = imgs(k).x;
            ydata = imgs(k).y;
            zdata = imgs(k).(colvar{j});
            imagesc(xdata,ydata,zdata)
            colorbar
            cax.YDir = 'Normal';
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 10;
            if i == size(ax,1), xlabel('x (cm)'), end
            if j == 1, ylabel('y (cm)'), end
            title(colstr{j},'FontWeight','normal')
            if strcmp(colvar{j},'n_sim'), cmax = max(zdata,[],'all'); end
            if strcmp(colvar{j},'n_exp'), cax.CLim = [0 cmax]; end
            if strcmp(colvar{j},'n_res'), cax.CLim = [-cmax cmax]/4; end

        end
    end
    str1 = ['t = ' num2str(imgs(k).t*1e6,'%.3g') '\mus = ' num2str(imgs(k).t/data.tau,'%.3g') '\tau_e_x_p'];
    an.String = str1;
    an.Position = [0.159500000000000,0.929069291338582,0.723000000000000,0.080100000000000];
    saveas(fig,[data.folder f 'imgs-' num2str(k) '.png']);
end
close(fig)


% do same thing but for lif
% create struct to hold info to be plotted
map = struct();
map(length([os.delays])).t = [];
for i = 1:length([os.delays])
    map(i).t = os(i).delays*1e-9;
    map(i).x = os(i).map.x/10;
    map(i).y = os(i).map.y/10;
    [map(i).X, map(i).Y] = meshgrid(map(i).x, map(i).y);
    map(i).n_exp = os(i).map.nFit*1e8;
    map(i).n_sim = interp2(data.grids.pos_x,data.grids.pos_y,data.grids.vars(ind_t(i)).n,map(i).X,map(i).Y);
    map(i).n_res = map(i).n_sim - map(i).n_exp;
    map(i).v_exp = os(i).map.vExp*100;
    map(i).v_sim = interp2(data.grids.pos_x,data.grids.pos_y,data.grids.vars(ind_t(i)).v_x,map(i).X,map(i).Y);
    map(i).v_res = map(i).v_sim - map(i).v_exp;
    fit_exp = fitImgWithGaussian(map(i).x,map(i).y,map(i).n_exp);
    fit_sim = fitImgWithGaussian(map(i).x,map(i).y,map(i).n_sim);
    map(i).sigx_exp = fit_exp.sigx;
    map(i).sigx_sim = fit_sim.sigx;
    map(i).sigy_exp = fit_exp.sigy;
    map(i).sigy_sim = fit_sim.sigy;
    fac_exp = map(1).sigx_exp*map(1).sigy_exp^2/(map(i).sigx_exp*map(i).sigy_exp^2);
    fac_sim = map(1).sigx_sim*map(1).sigy_sim/(map(i).sigx_sim*map(i).sigy_sim);
    map(i).n_exp_norm = map(i).n_exp./fac_exp;
    map(i).n_sim_norm = map(i).n_sim./fac_sim;
    map(i).n_res_norm = map(i).n_sim_norm - map(i).n_exp_norm;
end

% plot lif with residuals
colvar = {'n_sim','n_exp','n_res','v_sim','v_exp','v_res'};
colstr = {'n_s_i_m','n_e_x_p','n_r_e_s','v_s_i_m','v_e_x_p','v_r_e_s'};
num = length(colvar);
[fig,ax,an] = open_subplot(num,'Visible','on');
an.Position = [0.1595    0.9084    0.7230    0.0801];

for k = 1:length([map.t])
    disp(['Plotting: ' num2str(k) '/' num2str(length([map.t]))])
    
    % update axis information
    iter = 0;
    for i = 1:size(ax,1)
        for j = 1:size(ax,2)
            iter = iter + 1;
            if iter > num, break, end

            cax = get_axis(fig,ax{i,j});
            xdata = map(k).x;
            ydata = map(k).y;
            zdata = map(k).(colvar{iter});

            imagesc(xdata,ydata,zdata)
            colorbar
            cax.YDir = 'Normal';
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 10;
            if i == size(ax,1), xlabel('x (cm)'), end
            if j == 1, ylabel('y (cm)'), end
            title(colstr{iter},'FontWeight','normal')
            if strcmp(colvar{iter},'n_sim'), cmax = max(zdata,[],'all'); end
            if strcmp(colvar{iter},'n_exp'), cax.CLim = [0 cmax]; end
            if strcmp(colvar{iter},'n_res'), cax.CLim = [-cmax cmax]/4; end
            if strcmp(colvar{iter},'v_sim'), cmax = max(zdata,[],'all'); end
            if strcmp(colvar{iter},'v_exp'), cax.CLim = [-cmax cmax+.00001]; end
            if strcmp(colvar{iter},'v_res'), cax.CLim = [-cmax cmax+.00001]/4; end
        end
    end
    str1 = ['t = ' num2str(map(k).t*1e6,'%.3g') '\mus = ' num2str(map(k).t/data.tau,'%.3g') '\tau_e_x_p'];
    an.String = str1;
    an.Position = [0.159500000000000,0.929069291338582,0.723000000000000,0.080100000000000];

    saveas(fig,[data.folder f 'lif-' num2str(k) '.png']);
end
close(fig)

end