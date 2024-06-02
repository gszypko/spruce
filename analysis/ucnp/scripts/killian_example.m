exp_folder = 'C:\Users\grant\OneDrive\Research\projects\23.12.12-UCNP-Shock-Paper-1\Phase-7.11';
sim_folder = 'set_0';
load([exp_folder filesep 'os.mat'],'os')
load([exp_folder filesep sim_folder filesep 'data.mat'],'data');

Te = zeros(size(data.grids.time));
for i = 1:length(data.grids.time)
    Te(i) = mean(data.grids.vars(i).e_temp,'all');
end