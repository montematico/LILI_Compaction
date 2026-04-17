function plot_compaction_results(dataFile)
close all

if nargin < 1; dataFile = 'LunarCompactionResults.mat'; end
load(dataFile);

% --- Plotting Configuration ---
color_percentile = 85; % Upper percentile limit for color axis
show_pareto = false;    % Toggle for overlaying Pareto optimal points

% 1. Design Robustness Histogram (New Figure)
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

% 2. Calculate Specific Energy (kWh/kg)
Specific_Energy = Energy_results ./ M_ROV_GRID;

% 3. Scatter Plots Styled to Match Reference Graphics
pareto_m_rov = M_ROV_GRID(pareto_idx);
pareto_energy = Energy_results(pareto_idx);
pareto_mass_ratios = M_RATIO_GRID(pareto_idx);
pareto_radius = R_ROLLER_GRID(pareto_idx);

% Figure 1: 2D Pareto Front Scatter (Mass Ratio)
figure('Name', 'Pareto Front: Mass vs. Mission Time');
scatter(M_ROV_GRID(valid_idx), Time_results(valid_idx), 80, M_RATIO_GRID(valid_idx), 'filled', 'MarkerEdgeColor', 'k');
hold on;
if show_pareto
    scatter(M_ROV_GRID(pareto_idx), Time_results(pareto_idx), 120, 'r', 'LineWidth', 2);
end
xlabel('Total Rover Mass (kg)');
ylabel('Total Mission Time (hr)');
grid on;
title('Total Mass vs. Mission Time (Mass Ratio Color)');
colormap(flipud(parula));

% Apply Color Axis Limits
color_data = M_RATIO_GRID(valid_idx);
upper_limit = prctile(color_data, color_percentile);
clim([min(color_data), upper_limit]);

h = colorbar;
ylabel(h, 'Mass Ratio');
if show_pareto
    legend('Feasible Designs', 'Pareto Optimal Front', 'Location', 'northeast');
else
    legend('Feasible Designs', 'Location', 'northeast');
end
hold off;

% Figure 2: Traction Requirements
figure('Name', 'Traction Requirements');
scatter(M_ROV_GRID(valid_idx), Time_results(valid_idx), 80, Max_Traction_results(valid_idx), 'filled', 'MarkerEdgeColor', 'k');
hold on;
if show_pareto
    scatter(M_ROV_GRID(pareto_idx), Time_results(pareto_idx), 120, 'r', 'LineWidth', 2);
end
xlabel('Total Rover Mass (kg)');
ylabel('Total Mission Time (hr)');
title('Traction Requirement: Mass vs. Time (Traction Force Color)');
grid on;

% Apply Color Axis Limits
color_data = Max_Traction_results(valid_idx);
upper_limit = prctile(color_data, color_percentile);
clim([min(color_data), upper_limit]);

h2 = colorbar;
ylabel(h2, 'Max Traction Force (N)');
if show_pareto
    legend('Feasible Designs', 'Pareto Optimal Front', 'Location', 'northeast');
else
    legend('Feasible Designs', 'Location', 'northeast');
end
hold off;

% Figure 5: Trade Space - Mass vs. Time (Specific Energy)
figure('Name', 'Trade Space: Mass vs. Time (Specific Energy)');
scatter(M_ROV_GRID(valid_idx), Time_results(valid_idx), 80, Specific_Energy(valid_idx), 'filled', 'MarkerEdgeColor', 'k');
hold on;
if show_pareto
    scatter(M_ROV_GRID(pareto_idx), Time_results(pareto_idx), 120, 'r', 'LineWidth', 2); 
end
xlabel('Total Rover Mass (kg)');
ylabel('Total Mission Time (hr)');
grid on;
title('Total Mass vs. Mission Time (Specific Energy Color)');
colormap(turbo);

% Apply Color Axis Limits
color_data = Specific_Energy(valid_idx);
upper_limit = prctile(color_data, color_percentile);
clim([min(color_data), upper_limit]);

h5 = colorbar;
ylabel(h5, 'Specific Energy (kWh/kg)');
if show_pareto
    legend('Feasible Designs', 'Pareto Optimal Front', 'Location', 'northeast');
else
    legend('Feasible Designs', 'Location', 'northeast');
end
hold off;

% Figure 6: Trade Space - Mass vs. Time (Energy Color)
figure('Name', 'Trade Space: Mass vs. Time (Energy)');
scatter(M_ROV_GRID(valid_idx), Time_results(valid_idx), 80, Energy_results(valid_idx), 'filled', 'MarkerEdgeColor', 'k');
hold on;
if show_pareto
    scatter(M_ROV_GRID(pareto_idx), Time_results(pareto_idx), 120, 'r', 'LineWidth', 2); 
end
xlabel('Total Rover Mass (kg)');
ylabel('Total Mission Time (hr)');
grid on;
title('Total Mass vs. Mission Time (Energy Color)');
colormap(turbo);

% Apply Color Axis Limits
color_data = Energy_results(valid_idx);
upper_limit = prctile(color_data, color_percentile);
clim([min(color_data), upper_limit]);

h6 = colorbar;
ylabel(h6, 'Total Energy (kWh)');
if show_pareto
    legend('Feasible Designs', 'Pareto Optimal Front', 'Location', 'northeast');
else
    legend('Feasible Designs', 'Location', 'northeast');
end
hold off;

end
