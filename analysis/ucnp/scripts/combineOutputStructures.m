% initialize script workspace
clc, clearvars -except inp, close all, setpath;
% clear and open <inp> cell
if ~exist('inp','var')
    %%
    inp = {}; open('inp');
end
% LIF analysis subfolder
andir = 'an-mhd-09.26.22';
% combine output structures for each row of <inp>
os_combined = struct();
for i = 1:length(inp)
    % load output structure for current data set
    load([inp{i} filesep andir filesep 'os.mat']);
    % copy fields from loaded output structure into combined output structure
    fields = fieldnames(os);
    for j = 1:length(fields)
        % loop over all time points contained within the loaded output structure
        t = [os.delays];
        for k = 1:length(t)
            os_combined(i).(fields{j}) = os(k).(fields{j});
        end
    end
end
os = os_combined;
clear os_combined;