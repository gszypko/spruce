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

% output results
out = struct;
out.eq_set = active_eqn_set;

end