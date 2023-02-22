function [out] = readConfigFile(path)

% read in contents of config file
C = readfile(path);

% identify the active equation set
eqn_sets = {'ideal_mhd','ideal_mhd_cons','ideal_mhd_2E','ideal_2F'};
eqn_set_found = false;
for i = 1:length(C)
    % process current line of config file, determine whether eqn set and if active
    lhs = strtrim(extractBefore(C{i},'='));
    rhs = strtrim(extractAfter(C{i},'='));
    is_eqn_set = max(strcmp(lhs,eqn_sets));
    is_active = strcmp(rhs,"true");

    % throw error if more than one active eqn set is found
    if is_eqn_set && is_active && eqn_set_found
        error("More than one equation set is active.");
    end

    % store eqn set name if it is active
    if is_eqn_set && is_active
        active_eqn_set = lhs;
        eqn_set_found = true;
    end
end

% identify if EIC option is on or off
ind = startsWith(C,'eic_thermalization');
eic_opt = C(ind);
if length(eic_opt) > 1, error('More than one config exists for <eic_thermalization>.'); end
eic_opt = strtrim(extractAfter(eic_opt,'='));
if strcmp(eic_opt,'true'), eic_opt = true;
elseif strcmp(eic_opt,'false'), eic_opt = false;
else, error('<eic_thermalization> config must be either <true> or <false>.');
end
out.eq_set = active_eqn_set;
out.eic_opt = eic_opt;

% identify global temperature option
out.global_temp_e = false;
out.global_temp_i = false;
ind = startsWith(C,'global_temperature');
global_temp = C(ind);
if length(global_temp) > 1, error('More than one config exists for <global_temperature>.'); end
global_temp = strtrim(extractAfter(global_temp,'='));
if strcmp(global_temp,'true')
    ind = contains(C,'gt_species');
    species = strtrim(extractAfter(C(ind),'='));
    species = strsplit(species{1},',');
    for i = 1:length(species)
        if strcmp(species{i},'e'), out.global_temp_e = true;
        elseif strcmp(species{i},'i'), out.global_temp_i = true;
        end
    end
end

end