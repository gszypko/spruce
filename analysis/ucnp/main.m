%% Initiate Program
clc, clearvars -except data, close all, f = filesep; setpath;

inp = {'folder';
    'C:\Users\grant\OneDrive\Research\mhd\projects\iaw-density-dist\Killian2012-Fig4a\set_1';
    'C:\Users\grant\OneDrive\Research\mhd\projects\iaw-density-dist\Killian2012-Fig4b\set_1';
    'C:\Users\grant\OneDrive\Research\mhd\projects\iaw-density-dist\Killian2012-Fig4c\set_1';
    'C:\Users\grant\OneDrive\Research\mhd\projects\iaw-density-dist\Killian2012-Fig4d\set_1';
      };
s = spreadsheet2struct(inp,inp(1,:));

% main options
flags.load_mat_file = false;

% loading options
flags.sim.time_window = [0 100]; % [start end] percentage of time window to load
flags.sim.time_interval = 1; % interval for time indices to keep
flags.sim.plot_window = [1 1];
flags.sim.ghost_cells = 2;

% simulation flags
flags.sim.doGaussianFits2D = true;
flags.sim.test = false;
flags.sim.plot_grids = false;
flags.sim.vars = {'n', 'v_x', 'v_y','i_temp','e_temp','dt'};
% flags.sim.vars = {'n','e_temp'};
flags.sim.vlasov_analysis = false;
flags.sim.cmpr_exp_sim = false;
flags.sim.iaw_analysis = true;
flags.sim.ion_holes = false;
flags.sim.charge_neutrality = false;
flags.sim.hole_orientation = 0;

flags.sim.figvis = 'on';

% experimental flags
flags.exp.analyze_experiment = false;
flags.exp.gen_state = true; % turns on the script to generate init.state from experimental conditions
flags.exp.set = 16; % unique set identifier - modifies directory structure as \...\set_XXX
flags.exp.t_max = 2; % units of tau
flags.exp.num_time_pts = 100; % sets time interval for recording grids during simulation
flags.exp.apply_sg_filt = true; % whether or not to apply an SG filter to the density distribution
flags.exp.trim_exp_domain = [.55 .55]; % [x y] experimental domain limits in cm, symmetric about origin
flags.exp.sim_domain = [.93 .93]; % [x y] domain limits in cm, symmetric about origin
flags.exp.num_grids = 301; % number of points on x axis for simulation domain
flags.exp.n_dist = 'gauss'; % gauss or exp
flags.exp.n_min = 1e-4; % minimum in simulation density distribution cm^-^3
flags.exp.is_uniform = true;
flags.exp.grid_growth = 1.04;
flags.exp.grid_spread = 0.9997;
flags.exp.bgd_radius = .45; % radius outside of which the background sg filter is applied
flags.Te = []; % electron temperature (K) to use for this data set (if empty, use Te in exp data settings)
flags.exp.useLIFFits = false;
flags.exp.figvis = 'on';
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
        % load simulation data
        files = dir(s(i).folder);
        data_mat_found = max(strcmp({files.name},'data.mat'));
        if data_mat_found && flags.load_mat_file
            load([s(i).folder f 'data.mat'],'data');
        else
            data = loadData(s(i).folder,flags.sim);
            save([data.folder f 'data.mat'],"data",'-mat');
        end
        
        % run simulation options
        if flags.sim.doGaussianFits2D, data = doGaussianFits2D(data,flags); end
        if flags.sim.test, simTest(data,flags.sim); end
        if flags.sim.plot_grids, plotGridEvol(data,flags.sim); end
        if flags.sim.vlasov_analysis, gaussianAnalysis(data,flags.sim); end
        if flags.sim.cmpr_exp_sim, [out] = compareExpAndSimData(data,flags.sim); end
        if flags.sim.ion_holes, ionHolesAnalysis(data,flags.sim); end
        if flags.sim.iaw_analysis, iawAnalysis(data,flags.sim); end
        if flags.sim.charge_neutrality, chargeNeutralityAnalysis(data,flags.sim); end
    end
    if is_exp_folder
        flags.exp.Te = [];
        if ~isempty(flags.Te)
            if length(flags.Te) == length(s), flags.exp.Te = flags.Te(i);
            else, error('Length of flags.exp.Te must be empty or equal to number of folders.');
            end
        end
        if flags.exp.gen_state, genExpState(s(i).folder,flags.exp); end
        if flags.exp.analyze_experiment,[gauss_fit_int,vel_data,gauss_fit_rot] = analyzeExperiment(s(i).folder,flags.exp); end
    end
end
disp('Analysis Complete.')


