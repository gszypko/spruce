function [] = ionHolesAnalysis(data,flags)
%% Check for Required Grids
required_grids = {'n','v_x','v_y','i_temp','e_temp'};
vars = fieldnames(data.grids.vars);
for i = 1:length(required_grids)
    if ~max(strcmp(vars,required_grids{i}))
        error(['Expected variable <' required_grids{i} '> is missing.'])
    end
end

%% Determine Ion Hole Orientation
[~,ind_t] = min(data.grids.time);
x = data.grids.x_vec;
y = data.grids.y_vec;
n = data.grids.vars(ind_t).n;
[gauss_fit,Rx,Ry] = extractHoleOrientation(x,y,n,flags.figvis,flags.hole_orientation,data.folder);

%% Project Data on Grid Aligned with Ion Hole Propagation Axis
% Uses grid transformation from 'extractHoleOrientation'
vars_rot = struct;
vars_rot(length(data.grids.time)).t = [];
for k = 1:length(data.grids.time)
    % record current time
    vars_rot(k).t = data.grids.time(k);

    % original spatial grids
    x = data.grids.x_vec;
    y = data.grids.y_vec;
    [X,Y] = meshgrid(x,y);
    
    % rotated spatial grids
    x_rot = data.grids.x_vec-gauss_fit.x0;
    y_rot = data.grids.y_vec-gauss_fit.y0;
    [X_rot,Y_rot] = meshgrid(x_rot,y_rot);
    
    % rotate grids to obtain values along the ion hole propagation axis - perpendicular quantities are ignored
    vars_rot(k).x = x_rot;
    vars_rot(k).y = y_rot;
    vars_rot(k).n = interp2(X,Y,data.grids.vars(k).n,Rx(X_rot,Y_rot,-gauss_fit.theta),Ry(X_rot,Y_rot,-gauss_fit.theta));
    vars_rot(k).Te = interp2(X,Y,data.grids.vars(k).e_temp,Rx(X_rot,Y_rot,-gauss_fit.theta),Ry(X_rot,Y_rot,-gauss_fit.theta));
    vars_rot(k).Ti = interp2(X,Y,data.grids.vars(k).i_temp,Rx(X_rot,Y_rot,-gauss_fit.theta),Ry(X_rot,Y_rot,-gauss_fit.theta));
    v_x_rot = interp2(X,Y,data.grids.vars(k).v_x,Rx(X_rot,Y_rot,-gauss_fit.theta),Ry(X_rot,Y_rot,-gauss_fit.theta));
    v_y_rot = interp2(X,Y,data.grids.vars(k).v_y,Rx(X_rot,Y_rot,-gauss_fit.theta),Ry(X_rot,Y_rot,-gauss_fit.theta));
    vars_rot(k).v = Rx(v_x_rot,v_y_rot,gauss_fit.theta);
    vars_rot(k).pi = data.settings.m_i.*vars_rot(k).n.*vars_rot(k).v;
    
end

% average data over direction perpendicular to ion hole propagation
vars_1D = struct();
vars_1D(length(data.grids.time)).t = [];
for k = 1:length(data.grids.time)
    vars_1D(k).t = vars_rot(k).t;
    ind_x = abs(vars_rot(k).x) < 3.5*data.sig0;
    ind_y = abs(vars_rot(k).y) < gauss_fit.sigy_g/2;
    fields = fieldnames(vars_rot);
    vars_1D(k).r = vars_rot(k).x(ind_x);
    for i = 4:length(fields)
        vars_1D(k).(fields{i}) = mean(vars_rot(k).(fields{i})(ind_y,ind_x),1); 
    end
end

% create anonymous functions to interpolate 1D variables in time and space
[R,T] = meshgrid(vars_1D(1).r,[vars_1D.t]);
n = reshape([vars_1D.n]',[],length(vars_1D))';
Te = reshape([vars_1D.Te]',[],length(vars_1D))';
Ti = reshape([vars_1D.Ti]',[],length(vars_1D))';
v = reshape([vars_1D.v]',[],length(vars_1D))';
int.n = @(t,r) interp2(R,T,n,r,t);
int.v = @(t,r) interp2(R,T,v,r,t);
int.Te = @(t,r) interp2(R,T,Te,r,t);
int.Ti = @(t,r) interp2(R,T,Ti,r,t);
int.hole = @(t) ionHoleProp(t,int.v,int.Ti,int.Te,data.settings.m_i,data.settings.adiabatic_index);

%% Fit 1D Data to Extract Ion Hole Location

for k = 1:length(vars_1D)
    disp(['Extracting Hole Positions: ' num2str(k) '/' num2str(length(vars_1D))])
    [hole_fit(k)] = extractHoleLocation(vars_1D(k).t,vars_1D(k).r,vars_1D(k).n,vars_1D(k).v,data.sig0,data.tau,int.hole(vars_1D(k).t));
end

%% Plot Transects with Hole Locations
num = 4;
[fig,ax,an] = open_subplot(num,'Visible','on');
an.Position = [0.1567    0.9611    0.7230    0.0339];
frames = cell(size(vars_1D));
q = {'n','v','pi','Te'};
qstr = {'n','v','\pi','T_e'};

% update axis information
for k = 1:length(vars_1D)
    disp(['Plotting transect fits: ' num2str(k) '/' num2str(length([vars_1D.t]))])
    iter = 0;
    for i = 1:size(ax,1)
        for j = 1:size(ax,2)
            iter = iter + 1;
            if iter > num, break, end

            cax = get_axis(fig,ax{i,j});
            cla(cax)
            hold on
            ydata = vars_1D(k).(q{iter});
            plot(vars_1D(k).r,ydata,'.')
            plot([hole_fit(k).p_xL hole_fit(k).p_xL],[min(ydata) max(ydata)])
            plot([hole_fit(k).p_xR hole_fit(k).p_xR],[min(ydata) max(ydata)])
            cax.PlotBoxAspectRatio = [1 1 1];
            if i == size(ax,1), xlabel('x (cm)'), end
            title(qstr{iter},'FontWeight','normal')
            an.String = ['t = ' num2str(vars_1D(k).t*1e6) '\mus = ' num2str(vars_1D(k).t/data.tau,'%.2g') ' \tau'];
        end
    end
    frames{k} = getframe(fig);
end
close(fig)
write_video([data.folder filesep 'ion-hole-transects'],frames);

%% Plot Hole Position vs Analytic Model

num = 1;
[fig,~,an] = open_subplot(num,'Visible','on');
an.Position = [0.1567    0.9611    0.7230    0.0339];
hold on
grid on
grid minor
l = get_line_specs(4);
lc = l(1);
plot([hole_fit.t]./data.tau,[hole_fit.p_xL],'LineWidth',2,'MarkerSize',2,'Color',lc.col,'MarkerFaceColor',lc.col,'MarkerEdgeColor',lc.col)
lc = l(2);
plot([hole_fit.t]./data.tau,[hole_fit.p_xR],'LineWidth',2,'MarkerSize',2,'Color',lc.col,'MarkerFaceColor',lc.col,'MarkerEdgeColor',lc.col)
xan = int.hole([hole_fit.t]);
lc = l(3);
plot([hole_fit.t]./data.tau,-xan,'LineWidth',2,'MarkerSize',2,'Color',lc.col,'MarkerFaceColor',lc.col,'MarkerEdgeColor',lc.col)
lc = l(4);
plot([hole_fit.t]./data.tau,xan,'LineWidth',2,'MarkerSize',2,'Color',lc.col,'MarkerFaceColor',lc.col,'MarkerEdgeColor',lc.col)
lgd = legend({'Fit: x_L','Fit: x_R','Model: x_L','Model: x_R'});
lgd.Position = [0.1748    0.6723    0.2291    0.2252];
xlabel('t/\tau')
ylabel('Hole Pos. (cm)')
saveas(fig,[data.folder filesep 'hole-prop.png'])
close(fig)
end