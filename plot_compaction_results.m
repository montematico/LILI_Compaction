function plot_compaction_results(dataFile)
close all

if nargin < 1; dataFile = 'LunarCompactionResults.mat'; end
load(dataFile);

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

% 2. Scatter Plots Styled to Match Reference Graphics
pareto_m_rov = M_ROV_GRID(pareto_idx);
pareto_energy = Energy_results(pareto_idx);
pareto_mass_ratios = M_RATIO_GRID(pareto_idx);
pareto_radius = R_ROLLER_GRID(pareto_idx);

% Figure 1: 2D Pareto Front Scatter (Mass Ratio)
figure('Name', 'Pareto Front: Mass vs. Mission Time');
scatter(M_ROV_GRID(valid_idx), Time_results(valid_idx), 80, M_RATIO_GRID(valid_idx), 'filled', 'MarkerEdgeColor', 'k');
hold on;
scatter(M_ROV_GRID(pareto_idx), Time_results(pareto_idx), 120, 'r', 'LineWidth', 2);
xlabel('Total Rover Mass (kg)');
ylabel('Total Mission Time (hr)');
grid on;
title('Total Mass vs. Mission Time (Mass Ratio Color)');
colormap(flipud(parula));
h = colorbar;
ylabel(h, 'Mass Ratio');
legend('Feasible Designs', 'Pareto Optimal Front', 'Location', 'northeast');
hold off;

% Figure 2: Traction Requirements
figure('Name', 'Traction Requirements');
scatter(M_ROV_GRID(valid_idx), Time_results(valid_idx), 80, Max_Traction_results(valid_idx), 'filled', 'MarkerEdgeColor', 'k');
hold on;
scatter(M_ROV_GRID(pareto_idx), Time_results(pareto_idx), 120, 'r', 'LineWidth', 2);
xlabel('Total Rover Mass (kg)');
ylabel('Total Mission Time (hr)');
title('Traction Requirement: Mass vs. Time (Traction Force Color)');
grid on;
h2 = colorbar;
ylabel(h2, 'Max Traction Force (N)');
legend('Feasible Designs', 'Pareto Optimal Front', 'Location', 'northeast');
hold off;

% Figure 6: Trade Space - Mass vs. Time (Energy Color)
figure('Name', 'Trade Space: Mass vs. Time (Energy)');
scatter(M_ROV_GRID(valid_idx), Time_results(valid_idx), 80, Energy_results(valid_idx), 'filled', 'MarkerEdgeColor', 'k');
hold on;
scatter(M_ROV_GRID(pareto_idx), Time_results(pareto_idx), 120, 'r', 'LineWidth', 2); 
xlabel('Total Rover Mass (kg)');
ylabel('Total Mission Time (hr)');
grid on;
title('Total Mass vs. Mission Time (Energy Color)');
colormap(turbo);
h6 = colorbar;
ylabel(h6, 'Total Energy (kWh)');
legend('Feasible Designs', 'Pareto Optimal Front', 'Location', 'northeast');
hold off;

end
