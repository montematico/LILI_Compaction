clear; clc
% Load the workspace containing grid sweep variables and simulation constants
if isfile('LunarCompactionResults.mat')
    load('LunarCompactionResults.mat');
else
    error('Could not find results .mat file.');
end

% --- Design Requirements ---
slope_deg = 10;
safety_factor = 1.8;
F_slope = M_ROV_GRID * g_moon * sin(deg2rad(slope_deg));

% 1. Filter for all valid designs
% mu_req includes the force to overcome rollers + slope, scaled by safety factor
mu_req_all = ((Max_Traction_results + F_slope) * safety_factor) ./ (M_ROV_GRID * g_moon);

% Increased threshold to 1.1 to account for 1.8x SF and 10deg slope
mu_threshold = 1.1; 
candidates = find(Valid_mask & (mu_req_all < mu_threshold));

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

% --- Index 5: The "Utopia Compromise" (4D Optimization) ---
% Extract values for safe candidates
m_vals = M_ROV_GRID(candidates);
u_vals = mu_req_all(candidates);
e_vals = Energy_results(candidates);
d_vals = R_ROLLER_GRID(candidates) * n_rollers;

% Normalize to 0-1 scale
norm_m = (m_vals - min(m_vals)) ./ (max(m_vals) - min(m_vals));
norm_u = (u_vals - min(u_vals)) ./ (max(u_vals) - min(u_vals));
norm_e = (e_vals - min(e_vals)) ./ (max(e_vals) - min(e_vals));
norm_d = (d_vals - min(d_vals)) ./ (max(d_vals) - min(d_vals));

% Calculate 4D Euclidean distance to (0,0,0,0)
dist_utopia = sqrt(norm_m.^2 + norm_u.^2 + norm_e.^2 + norm_d.^2);
[~, local_utopia_idx] = min(dist_utopia);
idx_utopia = candidates(local_utopia_idx);

% Consolidate for the Comparison Script
comparison_indices = [idx_min_mass, idx_best_traction, idx_balanced, idx_light_robust, idx_utopia];
labels = {'Minimum Mass Survivor', 'Maximum Traction Margin', 'Balanced Design', 'Light-Robust Peak', 'Utopia Compromise'};

% Display results
fprintf('\n--- Selection Results ---\n');
% Dynamically construct header with wheel count
dbp_header = sprintf('DBP/%dw (N)', DRIVE.Nw);
fprintf('%-25s | %-6s | %-8s | %-6s | %-12s | %-12s | %-10s | %-20s\n', ...
    'Design Label', 'Index', 'Mass', 'mu_req', dbp_header, 'Energy (kWh)', 'Time (hr)', 'Spec. Energy (kWh/kg)');
fprintf('%s\n', repmat('-', 1, 125));
for i = 1:length(comparison_indices)
    idx = comparison_indices(i);
    % DBP per wheel = (F_rollers + F_slope) * SF / DRIVE.Nw
    F_total_req = (Max_Traction_results(idx) + F_slope(idx)) * safety_factor;
    dbp_per_wheel = F_total_req / DRIVE.Nw;
    
    fprintf('%-25s | %6d | %5.1f kg | %6.3f | %12.2f | %12.2f | %10.1f | %18.4f\n', ...
        labels{i}, idx, M_ROV_GRID(idx), mu_req_all(idx), dbp_per_wheel, Energy_results(idx), Time_results(idx), Energy_results(idx) / M_ROV_GRID(idx));
end

% --- Task 2: Tradeoff Scatter Plots ---
figure('Name', 'Design Space Tradeoffs', 'Color', 'w');

hold on; grid on;
% scatter(M_ROV_GRID(Valid_mask), mu_req_all(Valid_mask), 10, [0.8 0.8 0.8], 'filled'); % Background
scatter(M_ROV_GRID(candidates), mu_req_all(candidates), 20, 'b'); % Safe
plot(M_ROV_GRID(idx_utopia), mu_req_all(idx_utopia), 'p', 'MarkerFaceColor', 'r', 'MarkerSize', 15); % Utopia
% yline(0.45, '--r', 'LineWidth', 1.5);
xlabel('Total Rover Mass (kg)'); ylabel('Required Traction Coef. (\mu_{req})');
title('Mass vs. Traction Coef.');

