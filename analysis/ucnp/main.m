%% Initiate Program
clc, clearvars -except data, close all, f = filesep; setpath;

inp = {'folder';
    'C:\Users\Grant\OneDrive\Research\mhd\data-expsims\an-mhd-09.26.22\HIH-0\set_18';
    'C:\Users\Grant\OneDrive\Research\mhd\data-expsims\an-mhd-09.26.22\HIH-2\set_18';
      };
s = spreadsheet2struct(inp,inp(1,:));
Te = [100 100];

% simulation flags
flags.sim.plot_grids = true;
flags.sim.vlasov_analysis = false;
flags.sim.cmpr_exp_sim = false;
flags.sim.iaw_analysis = false;
flags.sim.ion_holes = false;
flags.sim.hole_orientation = 0;

% flags.sim.vars = {'n', 'v_x', 'v_y', 'i_temp', 'e_temp', 'dt'};
flags.sim.vars = {'i_n', 'dn', 'i_v_x', 'j_x', 'i_temp', 'e_temp', 'E_x', 'dt'};
flags.sim.plot_freq = 1;
flags.sim.ghost_cells = 2;
flags.sim.eic_opt = true;
flags.sim.figvis = 'on';

% experimental flags
flags.exp.gen_state = true; % turns on the script to generate init.state from experimental conditions
flags.exp.set = 18; % unique set identifier - modifies directory structure as \...\set_XXX
flags.exp.num_time_pts = 1000; % sets time interval for recording grids during simulation
flags.exp.smooth_density = true; % whether or not to apply an SG filter to the density distribution
flags.exp.sg_imgs_length = 0.02; % length scale (cm) for SG kernel when filtering density distribution
flags.exp.sim_window = [.55 .55]; % [x y] domain limits in cm, symmetric about zero
flags.exp.exp_window = [.55 .55]; % [x y] experimental domain limits in cm, symmetric about zero
flags.exp.sg_back_length = .05;
flags.exp.Nx = 301; % number of points on x axis for simulation domain
flags.exp.Ny = 301; % number of points on y axis for simulation domain
flags.exp.n_min = 1e5; % minimum in simulation density distribution cm^-^3
flags.exp.Te = [];
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
        data = loadData(s(i).folder,flags.sim);
        if flags.sim.plot_grids, plotGridEvol(data,flags.sim); end
        if flags.sim.vlasov_analysis, gaussianAnalysis(data,flags.sim); end
        if flags.sim.cmpr_exp_sim, compareExpAndSimData(data,flags.sim); end
        if flags.sim.ion_holes, ionHolesAnalysis(data,flags.sim); end
        if flags.sim.iaw_analysis, iawAnalysis(data,flags.sim); end
    end
    if is_exp_folder
        if ~isempty(Te), flags.exp.Te = Te(i); end
        if flags.exp.gen_state, genExpState(s(i).folder,flags.exp); end
    end
end
disp('Analysis Complete.')