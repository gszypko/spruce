function [] = ionHolesAnalysis2(data,flags)
f = filesep; % operating system file separator
theta = 20; % initial guess for rotation angle that defines principle axes (degree)

%% Check for Required Grids
% In this section, I compute any grids that I would like to work with that were not saved in the .out file
% I am expecting for grids n, v_x, v_y, i_temp, e_temp, dt, dPdx, dPdy to exist already, so I will run a check for that
expected_grids = {'n','v_x','v_y','i_temp','e_temp','dt','dPdx','dPdy'};
vars = fieldnames(data.grids.vars);
for i = 1:length(expected_grids)
    if ~max(strcmp(vars,expected_grids{i}))
        error(['Expected variable <' expected_grids{i} '> is missing.'])
    end
end

%% Fit Gaussian with Hole to Initial Time Point
[~,ind_t] = min(data.grids.time);
x = data.grids.x_vec;
y = data.grids.y_vec;
n = data.grids.vars(ind_t).n;
[fit,x_rot,y_rot] = extractHoleOrientation(x,y,n);

%% Project Data on Grid Aligned with Perturbation
prop = struct;
prop(length(data.grids.time)).t = [];
num = 1;
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
fig.Position = [257.0000  234.6000  784.0000  641.6000];
for k = 1:length(data.grids.time)
    % record current time
    prop(k).t = data.grids.time(k);

    % project density grid onto rotated coordinate system
    xp = data.grids.x_vec-hole.xp0;
    yp = data.grids.y_vec-hole.yp0;
    [Xp,Yp] = meshgrid(xp,yp);
    np = interp2(X,Y,data.grids.vars(k).n,x_rot(Xp,Yp,-hole.theta),y_rot(Xp,Yp,-hole.theta),'spline',min(data.grids.vars(k).n,[],'all'));

    % fit with 2D Gaussian
    [fit] = fitImgWithGaussian(xp,yp,np,false);
    ind_y = abs(yp) < fit.sigy/2;
    yp = yp(ind_y);
    np = np(ind_y,:);
    np = mean(np,1);
    

    % plot density distribution along rotated axis
    
    plot(xp,np)
    pause(0.5)
end



end