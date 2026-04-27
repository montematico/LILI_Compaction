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

% Figure 1: Unified Fraction (Roller Mass Fraction)
figure('Name', 'Trade Space: Unified Fraction');
plot_contour(X_data, Y_data, UNIFIED_FRAC_GRID(valid_idx), M_MESH, T_MESH, ...
    'Total Mass vs. Mission Time (Unified Fraction)', 'Unified Fraction', flipud(parula), show_pareto, pareto_idx, P_X, P_Y);

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

% Figure 7: Normalized Traction Margin - NEW
% Calculate Dynamic Limit (MU_max) for valid indices
h_g = 0.008; n_g = 18; slip_target = 0.20; 
if ~exist('SOIL','var'); SOIL.K = 0.02; end
MU_max_valid = zeros(size(valid_idx));
for i = 1:length(valid_idx)
    idx = valid_idx(i);
    W_drive_total = M_ROV_GRID(idx) * g_moon_val - (M_ROV_GRID(idx) * g_moon_val * UNIFIED_FRAC_GRID(idx));
    W_wheel = W_drive_total / DRIVE.Nw;
    perf = calculate_wheel_performance(W_wheel, D_WHEEL_GRID(idx), DRIVE.b, h_g, n_g, slip_target, v_sim, SOIL);
    MU_max_valid(i) = perf.MU_max;
end
Traction_Margin = MU_max_valid - mu_req_all;
Norm_Margin = (Traction_Margin - min(Traction_Margin)) ./ (max(Traction_Margin) - min(Traction_Margin));

figure('Name', 'Trade Space: Traction Margin');
plot_contour(X_data, Y_data, Norm_Margin, M_MESH, T_MESH, ...
    'Total Mass vs. Mission Time (Normalized Traction Margin)', 'Normalized Margin (0-1)', 'parula', show_pareto, pareto_idx, P_X, P_Y);

% --- New: 9x9 Performance Correlation Matrix ---
% Variables: Mass, Wheel Diameter, Wheel Width, Frequency, Mass Ratio, Time, Max Traction Force, Total Energy, Margin
M_vals  = M_ROV_GRID(valid_idx);
DW_vals = D_WHEEL_GRID(valid_idx);
BW_vals = B_WHEEL_GRID(valid_idx);
F_vals  = F_GRID(valid_idx);
MR_vals = UNIFIED_FRAC_GRID(valid_idx);
T_vals  = Time_results(valid_idx);
TF_vals = Max_Traction_results(valid_idx);
E_vals  = Energy_results(valid_idx);
N_vals  = Norm_Margin(:);

% Consolidated matrix (9 columns)
vars_matrix = [M_vals(:), DW_vals(:), BW_vals(:), F_vals(:), MR_vals(:), T_vals(:), TF_vals(:), E_vals(:), N_vals];
var_labels = {'Mass (kg)', 'Whl D (m)', 'Whl B (m)', 'Freq (hz)', 'Mass Ratio (%)', 'Time (hr)', 'Trac. (N)', 'Energy (kWh)', 'Traction Margin (%)'};

figure('Name', '9x9 Performance Matrix', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], 'Color', 'w');
[S, AX, BigAx, H, HAx] = plotmatrix(vars_matrix);

% Remove all histograms explicitly
for i = 1:length(HAx)
    set(HAx(i), 'Visible', 'off');
    delete(get(HAx(i), 'Children'));
end

num_vars = length(var_labels);
for i = 1:num_vars
    for j = 1:num_vars
        if i >= j
            % Remove lower triangular subplots and redundant diagonal axes
            set(AX(i,j), 'Visible', 'off');
            delete(get(AX(i,j), 'Children'));
        else
            % Enable grid and numeric scales for all visible scatter plots
            grid(AX(i,j), 'on');
            set(AX(i,j), 'Box', 'on', 'Visible', 'on');
            set(S(i,j), 'MarkerSize', 8); % Slightly increase marker size
            
            % Axis numbers (ticks) only once per row/column (on the staircase)
            if j == i + 1
                set(AX(i,j), 'YTickLabelMode', 'auto');
            else
                set(AX(i,j), 'YTickLabel', '');
            end
            
            if i == j - 1
                set(AX(i,j), 'XTickLabelMode', 'auto');
                xtickangle(AX(i,j), 45); % Rotate x-tick labels to prevent overlap
            else
                set(AX(i,j), 'XTickLabel', '');
            end
        end
    end
end

% Label placement per request: Move axis labels to the far left/bottom edges
for i = 1:num_vars
    % Y-labels on the far left (column 1)
    ylab = get(AX(i, 1), 'YLabel');
    set(ylab, 'String', var_labels{i}, 'FontWeight', 'bold', 'FontSize', 9, 'Visible', 'on');
    
    % X-labels on the far bottom (row num_vars)
    xlab = get(AX(num_vars, i), 'XLabel');
    set(xlab, 'String', var_labels{i}, 'FontWeight', 'bold', 'FontSize', 9, 'Visible', 'on');
end
title(BigAx, 'Performance Sensitivity & Correlation Matrix');

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

