function plot_compaction_results(dataFile)
close all

if nargin < 1; dataFile = 'LunarCompactionResults.mat'; end
load(dataFile);

% Plotting Configuration
color_percentile = 85; % Upper percentile limit for color axis
show_pareto = false;    % Toggle for overlaying Pareto optimal points

% 1. Calculate Required Traction Coefficient
slope_deg = 10;
safety_factor = 1.8;
g_moon_val = 1.62; % Using local g_moon if available, else 1.62
if exist('g_moon', 'var'); g_moon_val = g_moon; end

% mu_req = ((Max_Traction + (M_ROV * g * sin(theta))) * SF) / (M_ROV * g)
mu_req_all = ((Max_Traction_results(valid_idx) + (M_ROV_GRID(valid_idx) * g_moon_val * sin(deg2rad(slope_deg)))) * safety_factor) ./ (M_ROV_GRID(valid_idx) * g_moon_val);

% 2. Grid for Contours
% Create a regular mesh for interpolation
m_vec = linspace(min(M_ROV_GRID(valid_idx)), max(M_ROV_GRID(valid_idx)), 100);
t_vec = linspace(min(Time_results(valid_idx)), max(Time_results(valid_idx)), 100);
[M_MESH, T_MESH] = meshgrid(m_vec, t_vec);

% 3. Design Robustness Histogram 
unique_masses = unique(M_ROV_GRID);
success_rates = zeros(size(unique_masses));
for i = 1:length(unique_masses)
    mask = (M_ROV_GRID == unique_masses(i));
    success_rates(i) = (sum(Valid_mask(mask), 'all') / sum(mask, 'all')) * 100;
end

subs = sprintf("%d Iterations per mass",(total_iters/length(m_rov_range)));
figure('Name', 'Configuration Success Rate');
bar(unique_masses, success_rates);
xlim([min(m_rov_range), max(m_rov_range)]);
grid on;
xlabel('Total Rover Mass (kg)');
ylabel('Success Rate within Mass Bucket (%)');
title('Design Robustness: Success Rate vs. Rover Mass');
subtitle(subs);

% 4. Calculate Specific Energy (kWh/kg)
Specific_Energy = Energy_results ./ M_ROV_GRID;

% Data for plots
X_data = M_ROV_GRID(valid_idx);
Y_data = Time_results(valid_idx);
P_X = M_ROV_GRID(pareto_idx);
P_Y = Time_results(pareto_idx);

% Figure 1: Mass Ratio
figure('Name', 'Trade Space: Mass Ratio');
plot_contour(X_data, Y_data, M_RATIO_GRID(valid_idx), M_MESH, T_MESH, ...
    'Total Mass vs. Mission Time (Mass Ratio)', 'Mass Ratio', flipud(parula), show_pareto, pareto_idx, P_X, P_Y);

% Figure 2: Traction Requirements (Force)
figure('Name', 'Traction Requirements (Force)');
plot_contour(X_data, Y_data, Max_Traction_results(valid_idx), M_MESH, T_MESH, ...
    'Traction Force Requirement (N)', 'Max Traction Force (N)', 'hot', show_pareto, pareto_idx, P_X, P_Y);

% Figure 3: Required Traction Coefficient (mu_req) - NEW
figure('Name', 'Required Traction Coefficient (\mu)');
plot_contour(X_data, Y_data, mu_req_all, M_MESH, T_MESH, ...
    'Required Traction Coefficient (incl. 10^\circ Slope + SF 1.8)', '\mu_{req}', 'hot', show_pareto, pareto_idx, P_X, P_Y);

% Figure 5: Trade Space - Specific Energy
figure('Name', 'Trade Space: Specific Energy');
plot_contour(X_data, Y_data, Specific_Energy(valid_idx), M_MESH, T_MESH, ...
    'Total Mass vs. Mission Time (Specific Energy)', 'Specific Energy (kWh/kg)', 'turbo', show_pareto, pareto_idx, P_X, P_Y);

% Figure 6: Trade Space - Energy
figure('Name', 'Trade Space: Total Energy');
plot_contour(X_data, Y_data, Energy_results(valid_idx), M_MESH, T_MESH, ...
    'Total Mass vs. Mission Time (Energy)', 'Total Energy (kWh)', 'turbo', show_pareto, pareto_idx, P_X, P_Y);

end

% Helper function for contour plotting
function plot_contour(X, Y, Z, M_MESH, T_MESH, title_str, z_label, cmap, show_pareto, p_idx, p_x, p_y)
    % Switch to 'natural' neighbor interpolation for smoother surfaces
    Z_MESH = griddata(X, Y, Z, M_MESH, T_MESH, 'linear');

    % Filled contours with no edge lines
    contourf(M_MESH, T_MESH, Z_MESH, 20, 'LineColor', 'none');
    hold on;

    % Constant value lines (contour) and clabel removed for clarity per user request

    if show_pareto
        scatter(p_x, p_y, 100, 'p', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
    end
    xlabel('Total Rover Mass (kg)');
    ylabel('Total Mission Time (hr)');
    title(title_str);
    colormap(gca, cmap);
    cb = colorbar;
    ylabel(cb, z_label);
    grid on;
    hold off;
end

