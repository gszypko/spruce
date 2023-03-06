function [] = chargeNeutralityAnalysis(data,flags)

%% Check that Grids are Found
vars = {'E_x','E_1F_x'};
vars_found = zeros(size(vars));
for i = 1:length(vars)
    vars_found(i) = max(strcmp(fieldnames(data.grids.vars),vars{i}));
end
vars_found = min(vars_found);
if ~vars_found, error('Not all variables found for charge neutrality analysis.'); end

%% Get Electric Field Along X Axis

for i = 1:length(data.grids.time)
    p(i).t = data.grids.time(i);
    p(i).x = data.grids.x_vec;
    p(i).n = interp2(data.grids.x_vec,data.grids.y_vec,data.grids.vars(i).i_n,p(i).x,zeros(size(p(i).x)));
    p(i).E_x = p(i).n.*interp2(data.grids.x_vec,data.grids.y_vec,data.grids.vars(i).E_x,p(i).x,zeros(size(p(i).x)));
    p(i).E_1F_x = p(i).n.*interp2(data.grids.x_vec,data.grids.y_vec,data.grids.vars(i).E_1F_x,p(i).x,zeros(size(p(i).x)));
    p(i).E_res = p(i).E_x - p(i).E_1F_x;
end

%% Compute Sum of Squared Residuals

s = struct();
s(length(data.grids.time)).t = [];
for i = 1:length(data.grids.time)
    s(i).t = data.grids.time(i);
    s(i).nE = data.grids.vars(i).i_n*data.grids.vars(i).E_x;
    s(i).nE1F = data.grids.vars(i).i_n*data.grids.vars(i).E_1F_x;
    s(i).nE1F_res = s(i).nE - s(i).nE1F;
    s(i).nE1_res_norm = s(i).nE1F_res;
    s(i).nE1_res_norm_sum = sum(s(i).nE1_res_norm.^2,'all')/sum(s(i).nE1F.^2,"all");


end

%% Plot Electric Field with Residual

vars = {'E_x','E_1F_x','E_res'};
varstr = {'nE_x','nE_1_F_,_x','nE_r_e_s'};
num = 1;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
frames = cell(1,length(p));
for i = 1:length(p)
    hold off
    l = get_line_specs(length(vars));
    for k = 1:length(vars)
        xdata = [p(i).x];
        ydata = [p(i).(vars{k})];
        plot(xdata,ydata,'LineWidth',2,'MarkerSize',4,'Color',l(k).col,'MarkerFaceColor',l(k).col,'MarkerEdgeColor',l(k).col)
        hold on
    end
    xlabel('x (cm)')
    ylabel('Electric Field')
    lgd = legend(varstr);
    frames{i} = getframe(fig);
end

close(fig)
filepath = [data.folder filesep 'e-field'];
write_video(filepath,frames);

end