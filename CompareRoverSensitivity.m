% CompareRoverSensitivity.m
% Evaluates and compares the traction robustness of multiple lunar compaction 
% rover designs under degraded soil conditions.

clear; clc; close all;

%% Context & Inputs
% Load the workspace containing grid sweep variables and simulation constants
% Using LunarCompactionResults.mat instead of SweepResults.mat per project logic
if isfile('LunarCompactionResults.mat')
    load('LunarCompactionResults.mat', 'M_ROV_GRID', 'R_ROLLER_GRID', 'ROL_FRAC_GRID', 'SOIL', 'g_moon', 'c_f');
elseif isfile('SweepResults.mat')
    load('SweepResults.mat', 'M_ROV_GRID', 'R_ROLLER_GRID', 'ROL_FRAC_GRID', 'SOIL', 'g_moon', 'c_f');
else
    error('Could not find results .mat file.');
end

% Define target indices for light, medium, and heavy rover designs
% target_indices = [12137, 53, 65, 716]; 
target_indices = [36137,12065,12716,24266];
labels = {'Minimum mass', ...
          'Robust in bottom 25% of mass', ...
          'Traction robust', ...
          'Utopia Compromise'};

% Threshold for traction failure
mu_threshold = 0.45;

%% Core Task: Soil Degradation Loop
% ... (rest of the calculation loop remains same)
% Soil modulus shifts from -40% weaker to +10% stronger
soil_mod_shifts = linspace(-0.4, 0.1, 10);

% Initialize results matrix (Rows: Rovers, Columns: Shifts)
mu_req_results = zeros(length(target_indices), length(soil_mod_shifts));

% Predefined assumptions based on the prompt
n_rollers = 2;
b_roller = 0.3;

for i = 1:length(target_indices)
    idx = target_indices(i);
    
    % Extract base parameters from grid
    m_opt = M_ROV_GRID(idx);
    roller_fraction = ROL_FRAC_GRID(idx);
    
    % Update label to include mass
    labels{i} = sprintf('%s (%.0f kg)', labels{i}, m_opt);
    
    R_roller_opt = R_ROLLER_GRID(idx);
    
    % Derived parameters
    D_roller = R_roller_opt * 2;
    W_roller = (m_opt * g_moon * roller_fraction) / n_rollers;
    
    for j = 1:length(soil_mod_shifts)
        shift = soil_mod_shifts(j);
        
        % Modify base moduli
        kc_test = SOIL.kc * (1 + shift);
        kphi_test = SOIL.kphi * (1 + shift);
        
        % Calculate Pass 1 sinkage
        z_p1 = calculate_tandem_sinkage(W_roller, b_roller, D_roller, kc_test, kphi_test, n_rollers);
        z_safe_p1 = max(z_p1, 1e-6);
        
        % Calculate locomotion resistance forces (Bekker/Terzaghi)
        % Clamp alpha_arg to avoid imaginary numbers
        alpha_arg_p1 = max(-1, min(1, 1 - 2*z_safe_p1/D_roller));
        alpha_p1 = max(acos(alpha_arg_p1), 1e-6);
        l_o_p1 = z_safe_p1 * tan(pi/4 - SOIL.phi/2)^2;
        
        R_r = (W_roller * n_rollers) * c_f;
        R_c = 0.5 * (kc_test + b_roller * kphi_test) * z_safe_p1^2 * n_rollers;
        
        term1_b = (b_roller * sin(alpha_p1 + SOIL.phi)) / (2 * sin(alpha_p1) * cos(SOIL.phi));
        term2_b = 2 * z_safe_p1 * SOIL.c * SOIL.K_c + SOIL.gamma * z_safe_p1^2 * SOIL.K_gamma;
        term3_b = (l_o_p1^3 * SOIL.gamma / 3) * (pi/2 - SOIL.phi);
        term4_b = SOIL.c * l_o_p1^2 * (1 + tan(pi/4 + SOIL.phi/2));
        R_b = term1_b * term2_b + term3_b + term4_b;
        
        F_loco_p1 = R_r + R_c + R_b;
        
        % Calculate required coefficient of friction
        mu_req = F_loco_p1 / (m_opt * g_moon);
        mu_req_results(i, j) = mu_req;
    end
end

%% Plotting: The Robustness "Money Plot"
figure('Name', 'Rover Robustness', 'NumberTitle', 'off', 'Color', 'w');
hold on; grid on;

% Plot results for each rover
colors = lines(length(target_indices));
markers = {'o', 's', '^', 'd', 'v', 'p', 'h'};
for i = 1:length(target_indices)
    plot(soil_mod_shifts * 100, mu_req_results(i, :), '-', ...
         'Color', colors(i,:), 'Marker', markers{i}, 'MarkerFaceColor', colors(i,:), ...
         'LineWidth', 2, 'DisplayName', labels{i});
end

% Add horizontal dashed red line at y = mu_threshold with dynamic TeX label
ytext = sprintf('Traction Failure Limit (\\mu = %0.2f)', mu_threshold);
hY = yline(mu_threshold, 'r--', 'LineWidth', 2, ...
           'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
set(hY, 'Label', ytext, 'Interpreter', 'tex');

% Optional: Shaded red patch covering the "Danger Zone"
x_lims = [min(soil_mod_shifts * 100), max(soil_mod_shifts * 100)];
y_lims = ylim;
if y_lims(2) < mu_threshold + 0.1
    y_lims(2) = mu_threshold + 0.1;
    ylim(y_lims);
end
patch([x_lims(1) x_lims(2) x_lims(2) x_lims(1)], ...
      [mu_threshold mu_threshold y_lims(2) y_lims(2)], ...
      'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% IMPORTANT: Set X-axis direction to 'reverse'
set(gca, 'XDir', 'reverse');

% Formatting
xlabel('Soil Strength Degradation (%)', 'FontWeight', 'bold');
ylabel('Required Traction Coefficient (\mu_{req})', 'FontWeight', 'bold');
title('Rover Robustness: Traction Risk in Degraded Lunar Soil', 'FontWeight', 'bold', 'FontSize', 12);
legend('Location', 'best');
xlim(x_lims);

hold off;

%% Helper Functions
function z = calculate_tandem_sinkage(W_roller, b, D, kc, kphi, N_total_current)
    % Calculates sinkage for N tandem passes
    z_prev = 0;
    for i = 1:N_total_current
        term_load = (3 * W_roller) / (2 * (kc + b * kphi) * sqrt(D));
        % Max(0, ...) prevents complex numbers
        z_i = (term_load + max(0, z_prev)^(1.5))^(2/3); 
        z_prev = z_i;
    end
    z = z_prev;
end