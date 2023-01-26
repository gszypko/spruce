%% Initiate Program
clc, clearvars -except data, close all, f = filesep; setpath;

inp = {'folder';
    'C:\Users\Grant\OneDrive\Research\mhd\data-expsims\an-mhd-09.26.22\HIH-0\set_24';
    'C:\Users\Grant\OneDrive\Research\mhd\data-expsims\an-mhd-09.26.22\HIH-0\set_25';
      };
s = spreadsheet2struct(inp,inp(1,:));

% simulation flags
flags.sim.plot_grids = true;
flags.sim.vlasov_analysis = true;
flags.sim.cmpr_exp_sim = false;
flags.sim.iaw_analysis = false;
flags.sim.ion_holes = false;
flags.sim.hole_orientation = 0;
% flags.sim.vars = {'n', 'v_x', 'v_y', 'temp', 'dt'};
flags.sim.vars = {'i_n', 'dn', 'i_v_x', 'e_v_x', 'i_temp', 'e_temp', 'E_x', 'dt'};
flags.sim.plot_freq = 1;
flags.sim.ghost_cells = 2;
flags.sim.eic_opt = true;
flags.sim.figvis = 'on';

% experimental flags
flags.exp.gen_state = true; % turns on the script to generate init.state from experimental conditions
flags.exp.set = 22; % unique set identifier - modifies directory structure as \...\set_XXX
flags.exp.t_max = 1.25; % units of tau
flags.exp.num_time_pts = 100; % sets time interval for recording grids during simulation
flags.exp.apply_sg_filt = true; % whether or not to apply an SG filter to the density distribution
flags.exp.sg_filt_length = 0.025; % length scale (cm) for SG kernel when filtering density distribution
flags.exp.trim_exp_domain = [.55 .55]; % [x y] experimental domain limits in cm, symmetric about origin
flags.exp.sim_domain = [1 1]; % [x y] domain limits in cm, symmetric about origin
flags.exp.num_grids = 301; % number of points on x axis for simulation domain
flags.exp.n_dist = 'gauss'; % gauss or exp
flags.exp.n_min = 1e6; % minimum in simulation density distribution cm^-^3
flags.exp.is_uniform = false;
flags.exp.grid_growth = 1.06;
flags.exp.grid_spread = 0.99994;
flags.exp.apply_sg_bgd = true; % whether or not to apply an SG filter to the density distribution
flags.exp.sg_bgd_length = 0.075; % length scale (cm) for SG kernel when filtering density distribution
flags.exp.sg_radius = 0.5; % radius outside of which the background sg filter is applied
flags.Te = []; % electron temperature (K) to use for this data set (if empty, use Te in exp data settings)
flags.exp.useLIFFits = false;
flags.exp.figvis = 'off';
flags.exp.config_path = 'C:\Users\grant\Documents\GitHub\mhd\ucnp.config';
flags.exp.settings_path = 'C:\Users\grant\Documents\GitHub\mhd\ucnp.settings';


%% Read in and Analyze Data
disp('Starting Analysis...')
for i = 1:length(s)
    disp(['Data set: ' num2str(i-1) '/' num2str(length(s)-1)])
    
    % determine whether given folder is a simulation or experimental data folder
    files = dir(s(i).folder);
    is_exp_folder = max(strcmp({files.name},'os.mat'));
    is_sim_folder = max(strcmp({files.name},'mhd.out'));
    
    % run flagged options
    if is_sim_folder
        data = loadData(s(i).folder,flags.sim);
        if flags.sim.plot_grids, plotGridEvol(data,flags.sim); end
        if flags.sim.vlasov_analysis, gaussianAnalysis(data,flags.sim); end
        if flags.sim.cmpr_exp_sim, compareExpAndSimData(data,flags.sim); end
        if flags.sim.ion_holes, ionHolesAnalysis(data,flags.sim); end
        if flags.sim.iaw_analysis, iawAnalysis(data,flags.sim); end
    end
    if is_exp_folder
        flags.exp.Te = [];
        if ~isempty(flags.Te)
            if length(flags.Te) == length(s), flags.exp.Te = Te(i);
            else, error('Length of flags.exp.Te must be empty or equal to number of folders.');
            end
        end
        if flags.exp.gen_state, genExpState(s(i).folder,flags.exp); end
    end
end
disp('Analysis Complete.')


