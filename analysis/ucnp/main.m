%% Initiate Program
clc; clearvars -except data; close all; setpath;

inp = {'folder';
      };
s = spreadsheet2struct(inp,inp(1,:));

%% loading options
flags.load_mat_file = true;
flags.sim.time_window = [0 100]; % [start end] percentage of time window to load
flags.sim.time_interval = 1; % interval for time indices to keep
flags.sim.plot_window = [1 1];
flags.sim.ghost_cells = 2;

%% simulation flags
flags.sim.plot_grids = true;
flags.sim.vars = {'n', 'v_x', 'v_y', 'i_temp', 'e_temp', 'dt'};
flags.sim.doGaussianFits2D = false;
flags.sim.vlasov_analysis = false;
flags.sim.cmpr_exp_sim = false;
flags.sim.iaw_analysis = false;
flags.sim.ion_holes = false;
flags.sim.charge_neutrality = false;
flags.sim.test = false;
flags.sim.hole_orientation = 0;
flags.sim.figvis = 'on';

%% experimental flags
flags.exp.analyze_exp_data = true;
flags.exp.gen_state = true; % turns on the script to generate init.state from experimental conditions
flags.exp.set = 0; % unique set identifier - modifies directory structure as \...\set_XXX
flags.exp.n_dist = 'gauss'; % gauss or exp
flags.exp.t_max = 1.5; % units of tau
flags.exp.num_time_pts = 150; % sets time interval for recording grids during simulation
flags.Te = []; % electron temperature (K) to use for this data set (if empty, use Te in exp data settings)
flags.exp.trim_exp_domain = [.6 .6]; % [x y] experimental domain limits in cm, symmetric about origin
flags.exp.sim_domain = [1.08 1.08]; % [x y] domain limits in cm, symmetric about origin
flags.exp.n_min = 1e-4; % minimum in simulation density distribution cm^-^3
flags.exp.num_grids = 301; % number of points on x axis for simulation domain
flags.exp.is_uniform = false;
flags.exp.grid_growth = 1.05;
flags.exp.grid_spread = 0.99985;
flags.exp.bgd_radius = .5; % radius outside of which the background sg filter is applied
flags.exp.figvis = 'off';
flags.exp.config_path = 'C:\Users\grant\Documents\GitHub\mhd\ucnp.config';
flags.exp.settings_path = 'C:\Users\grant\Documents\GitHub\mhd\ucnp.settings';


%% Read and Analyze Data
disp('Starting Analysis...')
for i = 1:length(s)
    %% Determine whether experimental or simulation data set
    disp(['Data set: ' num2str(i-1) '/' num2str(length(s)-1)])
    
    % determine whether given folder is a simulation or experimental data folder
    files = dir(s(i).folder);
    is_exp_folder = max(strcmp({files.name},'os.mat'));
    is_sim_folder = max(strcmp({files.name},'mhd.out'));
    
    %% Handle when simulation folder
    if is_sim_folder
        % load simulation data
        files = dir(s(i).folder);
        data_mat_found = max(strcmp({files.name},'data.mat'));
        if data_mat_found && flags.load_mat_file
            load([s(i).folder filesep 'data.mat'],'data');
        else
            data = loadData(s(i).folder,flags.sim);
            save([data.folder filesep 'data.mat'],"data",'-mat');
        end
        
        % run simulation options
        if flags.sim.plot_grids, plotGridEvol(data,flags.sim); end
        if flags.sim.doGaussianFits2D, data = doGaussianFits2D(data,flags.sim); end
        if flags.sim.test, simTest(data,flags.sim); end
        if flags.sim.vlasov_analysis, vlasovAnalysis(data,flags.sim); end
        if flags.sim.cmpr_exp_sim, [out] = compareExpAndSimData(data,flags.sim); end
        if flags.sim.ion_holes, ionHolesAnalysis(data,flags.sim); end
        if flags.sim.iaw_analysis, iawAnalysis(data,flags.sim); end
        if flags.sim.charge_neutrality, chargeNeutralityAnalysis(data,flags.sim); end
    end
    %% Handle when experimental data folder
    if is_exp_folder
        flags.exp.Te = [];
        if ~isempty(flags.Te)
            if length(flags.Te) == length(s), flags.exp.Te = flags.Te(i);
            else, error('Length of flags.exp.Te must be empty or equal to number of folders.');
            end
        end

        if flags.exp.analyze_exp_data, analyzeExpData(s(i).folder); end
        if flags.exp.gen_state, genExpState(s(i).folder,flags.exp); end
    end
end
disp('Analysis Complete.')


