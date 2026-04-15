% Load the workspace containing grid sweep variables and simulation constants
% Using LunarCompactionResults.mat instead of SweepResults.mat per project logic
if isfile('LunarCompactionResults.mat')
    load('LunarCompactionResults.mat', 'M_ROV_GRID', 'R_ROLLER_GRID', 'SOIL', 'g_moon', 'c_f');
elseif isfile('SweepResults.mat')
    load('SweepResults.mat', 'M_ROV_GRID', 'R_ROLLER_GRID', 'SOIL', 'g_moon', 'c_f');
else
    error('Could not find results .mat file.');
end

% 1. Filter for all valid designs that pass the 0.45 traction threshold
% Using Max_Traction_results ./ (M_ROV_grid * g_moon) to get mu_req
mu_req_all = Max_Traction_results ./ (M_ROV_GRID * g_moon);
candidates = find(Valid_mask & (mu_req_all < 0.45));

if isempty(candidates)
    error('No designs found below the 0.45 traction limit. Consider increasing M_ROV_grid or R_ROLLER_grid ranges.');
end

% --- Index 1: The "Minimum Mass" Design (The Bare Minimum) ---
% The lightest possible rover that doesn't get stuck.
[~, local_min_mass_idx] = min(M_ROV_GRID(candidates));
idx_min_mass = candidates(local_min_mass_idx);

% --- Index 2: The "Highest Robustness" Design (The Tank) ---
% The design with the absolute lowest required traction (highest margin).
[~, local_low_mu_idx] = min(mu_req_all(candidates));
idx_best_traction = candidates(local_low_mu_idx);

% --- Index 3: The "Middle Ground" (Balanced) ---
% A design roughly halfway between the min and max mass of the successful candidates.
target_mid_mass = (max(M_ROV_GRID(candidates)) + min(M_ROV_GRID(candidates))) / 2;
[~, local_mid_idx] = min(abs(M_ROV_GRID(candidates) - target_mid_mass));
idx_balanced = candidates(local_mid_idx);

% --- Index 4: The "Mass-Efficient Robust" Design ---
% A design that is light but has a better safety margin than the bare minimum.
% We'll look for the design with the lowest mu_req among the bottom 25% of masses.
mass_threshold = quantile(M_ROV_GRID(candidates), 0.25);
light_candidates = candidates(M_ROV_GRID(candidates) <= mass_threshold);
[~, local_light_robust_idx] = min(mu_req_all(light_candidates));
idx_light_robust = light_candidates(local_light_robust_idx);

% Consolidate for the Comparison Script
comparison_indices = [idx_min_mass, idx_light_robust, idx_balanced, idx_best_traction];
labels = {'Minimum Mass Survivor', 'Light-Robust Peak', 'Balanced Design', 'Maximum Traction Margin'};

% Display results
fprintf('\n--- Selection Results ---\n');
for i = 1:length(comparison_indices)
    idx = comparison_indices(i);
    fprintf('%-25s | Index: %d | Mass: %5.1f kg | mu_req: %5.3f\n', ...
        labels{i}, idx, M_ROV_GRID(idx), mu_req_all(idx));
end