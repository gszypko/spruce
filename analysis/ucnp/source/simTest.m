function [out] = simTest(data,flags)

g = struct();
g(length(data.grids.time)).t = [];
for i = 1:length(data.grids.time)
    g(i).t = data.grids.time(i);
    x = data.grids.x_vec;
    y = data.grids.y_vec;
    g(i).N = trapz(y,trapz(x,data.grids.vars(i).n,2),1);
    g(i).eps_th = trapz(y,trapz(x,data.grids.vars(i).e_thermal_energy,2),1);
    n = data.grids.vars(i).n;
    Tg = interp2(x,y,data.grids.vars(i).e_temp_g,0,0);
    g(i).Tg = Tg;
    eps_g = n*cts.cgs.kB*Tg/(data.settings.adiabatic_index-1);
    g(i).eps_g = trapz(y,trapz(x,eps_g,2),1);
    g(i).Tc = interp2(x,y,data.grids.vars(i).e_temp,0,0);

end

num=3;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);

cax = get_axis(fig,ax{1,1});
plot([g.t],[g.N]./max([g.N]))
lgd = legend({'N'});

cax = get_axis(fig,ax{1,2});
plot([g.t],[g.eps_th])
hold on
plot([g.t],[g.eps_g])
lgd = legend({'eps_t_h','eps_g'});

cax = get_axis(fig,ax{1,3});
hold on
plot([g.t],[g.Tc])
plot([g.t],[g.Tg])
lgd = legend({'Tc','Tg'});

close(fig)

end