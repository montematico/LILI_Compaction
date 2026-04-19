% CompareRoverSensitivity.m
% Evaluates and compares the traction robustness of multiple lunar compaction 
% rover designs under degraded soil conditions.

clear; clc; close all;

%% Context & Inputs
% Load the workspace containing grid sweep variables and simulation constants
% Using LunarCompactionResults.mat instead of SweepResults.mat per project logic
if isfile('LunarCompactionResults.mat')
    load('LunarCompactionResults.mat', 'M_ROV_GRID', 'R_ROLLER_GRID', 'ROL_FRAC_GRID', 'SOIL', 'g_moon', 'c_f','n_rollers','b_roller','DRIVE','v_sim');
else
    error('Could not find results .mat file.');
end

% Standard Traction Configuration
slope_deg = 10;      % degrees - assumed design slope
safety_factor = 1.8; % SF for traction

% Define target indices for light, medium, and heavy rover designs
target_indices = [12126,12716,12053,12150,36494];
labels = {'Minimum mass', 'Maximum Traction Margin','Balanced (Mass \mu)', 'Balanced  25%th mass', 'Balanced (Mass,\mu,wheel d,energy)'};

%% Soil Degradation Loop
% Soil modulus shifts from -40% weaker to +40% stronger
soil_mod_shifts = linspace(-0.4, 0.4, 30);

% Initialize results matrix (Rows: Rovers, Columns: Shifts)
mu_req_results = zeros(length(target_indices), length(soil_mod_shifts));
mu_limits = zeros(length(target_indices), 1); % Store dynamic MU_max for each design

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

    % Calculate design-specific MU_max using nominal soil
    % W_wheel for traction check: assume mass is distributed over drive wheels
    W_wheel_nom = (m_opt * g_moon) / DRIVE.Nw;
    perf = calculate_wheel_performance(W_wheel_nom, D_roller, DRIVE.b, DRIVE.h_g, DRIVE.n_g, DRIVE.slip, v_sim, SOIL);
    mu_limits(i) = perf.MU_max;
    
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
        
        % Calculate required coefficient of friction (incl. slope and SF)
        mu_req = ((F_loco_p1 + (m_opt * g_moon * sin(deg2rad(slope_deg)))) * safety_factor) / (m_opt * g_moon);
        mu_req_results(i, j) = mu_req;
    end
end

%% Plotting: The Robustness "Money Plot"
figure('Name', 'Rover Robustness', 'NumberTitle', 'off', 'Color', 'w');
hold on; grid on;

% Plot results for each rover
colors = lines(length(target_indices));
for i = 1:length(target_indices)
    % Required Traction Curve
    plot(soil_mod_shifts * 100, mu_req_results(i, :), '-', ...
       'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', labels{i});
   
    % Individual Traction Limit Line
    % yline(mu_limits(i), '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

% IMPORTANT: Set X-axis direction to 'reverse'
set(gca, 'XDir', 'reverse');

% Formatting
xlabel('Soil Strength Degradation (%)', 'FontWeight', 'bold');
ylabel('Required Traction Coefficient (\mu_{req})', 'FontWeight', 'bold');
title('Rover Robustness: Traction Risk (Dynamic Limits)', 'FontWeight', 'bold', 'FontSize', 12);
subtitle(sprintf('Limits (dashed) vs Req. Traction (solid). Incl. %d^o Slope + SF %.1f', slope_deg, safety_factor));
legend('Location', 'best');
xlim([min(soil_mod_shifts * 100), max(soil_mod_shifts * 100)]);

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