%% SingleMassAnalysis.m
% Deep-dive engineering analysis on a single, optimal lunar compaction rover design.
% Based on grid sweep results and pass-by-pass simulation.

clear; clc; close all;

%% 1. Parameters & Data Loading
target_idx = 42; % Define target index at the top

% Load grid results
if exist('SweepResults.mat', 'file')
    load('SweepResults.mat');
elseif exist('LunarCompactionResults.mat', 'file')
    load('LunarCompactionResults.mat');
else
    error('SweepResults.mat not found. Please ensure grid sweep results are available.');
end

% Extract optimal parameters for this index
try
    % Try uppercase first as seen in LunarCompactionResults.mat
    if exist('M_ROV_GRID', 'var')
        m_opt = M_ROV_GRID(target_idx);
        f_opt = F_GRID(target_idx);
        mass_ratio_opt = M_RATIO_GRID(target_idx);
        R_roller_opt = R_ROLLER_GRID(target_idx);
        m_eccentric_opt = M_ECC_GRID(target_idx);
        roller_fraction_opt = ROL_FRAC_GRID(target_idx);
        D_wheel_opt = D_WHEEL_GRID(target_idx);
    else
        m_opt = M_ROV_grid(target_idx);
        f_opt = F_grid(target_idx);
        mass_ratio_opt = M_RATIO_grid(target_idx);
        R_roller_opt = R_ROLLER_grid(target_idx);
        m_eccentric_opt = M_ECC_grid(target_idx);
        roller_fraction_opt = ROL_FRAC_grid(target_idx);
        D_wheel_opt = D_WHEEL_grid(target_idx);
    end
catch
    error('Target index %d out of bounds or grid variables not found in mat file.', target_idx);
end

% Re-derive dependent parameters for the simulation
D_roller_opt = 2 * R_roller_opt;
W_roller_opt = (m_opt * g_moon * roller_fraction_opt) / n_rollers; 

% Update DRIVE struct for this iteration
DRIVE_opt = DRIVE;
DRIVE_opt.D = D_wheel_opt;

fprintf('--- Single Mass Analysis (Target Index: %d) ---\n', target_idx);
fprintf('Mass: %.2f kg, Freq: %.2f Hz, Mass Ratio: %.2f, Radius: %.2f m\n', m_opt, f_opt, mass_ratio_opt, R_roller_opt);

%% Core Task 1: Pass-by-Pass History
[hist, is_valid] = run_detailed_sim(W_roller_opt, b_roller, D_roller_opt, m_eccentric_opt, f_opt, v_sim, mass_ratio_opt, m_opt, target_relative_density, bounce_margin, h_layer, SOIL, eta_mech, nu, rho_min_lunar, rho_max_lunar, rho_i, Apad, n_rollers, c_f, g_moon, P_hotel, battery_density_Wh_kg, max_battery_fraction, t_work_cycle_h, max_sinkage_ratio, DRIVE_opt);

if ~is_valid
    warning('The simulation for the target design is invalid or failed.');
end

%% Core Task 2: Power & Battery Analysis
% Specific Energy (SE): Total energy / (Volume * Density Change)
total_energy_J = sum(hist.P_loco_hist + hist.P_vib_hist + hist.P_hotel_hist) .* (hist.d_pass_hist ./ v_sim);
total_energy_kWh = sum(total_energy_J) / 3.6e6;

rho_final = hist.rho_hist(end);
rho_initial = hist.rho_hist(1);
V_compacted = Apad * h_layer; % Initial layer volume

SE = sum(total_energy_J) / (V_compacted * (rho_final - rho_initial)); % J / (m^3 * kg/m^3) = J/kg

% Battery Sizing (80% DoD, 110 Wh/kg system-level)
dod = 0.80;
e_density_sys = 110; % Wh/kg
total_energy_Wh = total_energy_kWh * 1000;
energy_per_cycle_Wh = (total_energy_Wh / (hist.t_total / 3600)) * t_work_cycle_h; % Extrapolate to full work cycle
m_battery_req = energy_per_cycle_Wh / (e_density_sys * dod);

% C-Rate
P_total_pass = hist.P_loco_hist + hist.P_vib_hist + hist.P_hotel_hist;
P_peak = max(P_total_pass);
C_rate_peak = P_peak / energy_per_cycle_Wh; % P / E_total_capacity

fprintf('\n--- Power & Battery Analysis ---\n');
fprintf('Specific Energy (SE): %.2f J/kg\n', SE);
fprintf('Required Battery Mass: %.2f kg\n', m_battery_req);
fprintf('Peak C-Rate: %.2f C\n', C_rate_peak);

%% Core Task 3: Soil Uncertainty & Traction Margin
soil_mod_shifts = linspace(-0.3, 0.3, 5);
mu_req_hist = zeros(size(soil_mod_shifts));
traction_fail = false(size(soil_mod_shifts));

for i = 1:length(soil_mod_shifts)
    shift = soil_mod_shifts(i);
    kc_test = SOIL.kc * (1 + shift);
    kphi_test = SOIL.kphi * (1 + shift);
    
    % Test max required traction force (from Pass 1 or similar logic)
    % Re-calculating F_loco for pass 1 with shifted moduli
    z_p1 = calculate_tandem_sinkage(W_roller_opt, b_roller, D_roller_opt, kc_test, kphi_test, n_rollers);
    
    alpha_arg = max(-1, min(1, 1 - 2*z_p1/D_roller_opt));
    alpha = max(acos(alpha_arg), 1e-6);
    l_o = z_p1 * tan(pi/4 - SOIL.phi/2)^2;
    
    R_r = (W_roller_opt * n_rollers) * c_f;
    R_c = 0.5 * (kc_test + b_roller * kphi_test) * z_p1^2 * n_rollers;
    
    term1_b = (b_roller * sin(alpha + SOIL.phi)) / (2 * sin(alpha) * cos(SOIL.phi));
    term2_b = 2 * z_p1 * SOIL.c * SOIL.K_c + SOIL.gamma * z_p1^2 * SOIL.K_gamma;
    term3_b = (l_o^3 * SOIL.gamma / 3) * (pi/2 - SOIL.phi);
    term4_b = SOIL.c * l_o^2 * (1 + tan(pi/4 + SOIL.phi/2));
    R_b = term1_b * term2_b + term3_b + term4_b;
    
    F_loco_p1 = R_r + R_c + R_b;
    mu_req = F_loco_p1 / (m_opt * g_moon);
    
    mu_req_hist(i) = mu_req;
    if mu_req > 0.45
        traction_fail(i) = true;
    end
end

fprintf('\n--- Soil Sensitivity Analysis ---\n');
for i = 1:length(soil_mod_shifts)
    status = "OK";
    if traction_fail(i), status = "TRACTION FAILURE"; end
    fprintf('Shift: %+.0f%%, Required Mu: %.3f [%s]\n', soil_mod_shifts(i)*100, mu_req_hist(i), status);
end

%% Core Task 4: Plotting Generation
% 1. Compaction Evolution (2-Panel Subplot)
figure('Name', 'Compaction Evolution', 'Color', 'w');
subplot(2,1,1);
% Convert absolute density to Relative Density (ID)
ID_hist = (hist.rho_hist - rho_min_lunar) / (rho_max_lunar - rho_min_lunar);
pass_vec = 0:length(hist.rho_hist)-1;

plot(pass_vec, ID_hist, '-o', 'LineWidth', 1.5);
grid on; xlabel('Pass Number'); ylabel('Relative Density (0-1)');
title('Compaction Progress (Relative Density)');
xticks(pass_vec); % Ensure integer ticks for passes

subplot(2,1,2);
% Convert sinkage from meters to millimeters
z_mm_hist = hist.z_hist * 1000;
plot(pass_vec, z_mm_hist, '-s', 'LineWidth', 1.5, 'Color', 'r');
grid on; xlabel('Pass Number'); ylabel('Sinkage (mm)');
title('Sinkage vs. Pass Number');
xticks(pass_vec);

% 2. Power Profile (Stacked Bar)
figure('Name', 'Power Profile', 'Color', 'w');
bar_data = [hist.P_loco_hist', hist.P_vib_hist', hist.P_hotel_hist'];
bar(1:length(hist.P_loco_hist), bar_data, 'stacked');
grid on; xlabel('Pass Number'); ylabel('Power (W)');
legend('Locomotion', 'Vibration', 'Hotel');
title('Power Breakdown per Pass');

% 3. Soil Sensitivity Plot
figure('Name', 'Soil Sensitivity', 'Color', 'w');
plot(soil_mod_shifts*100, mu_req_hist, '-d', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
yline(0.45, 'r--', 'Traction Limit (0.45)', 'LineWidth', 2);
grid on; xlabel('% Change in Soil Moduli (kc, kphi)'); ylabel('Required \mu');
title('Soil Uncertainty & Traction Margin');
ylim([0, max(max(mu_req_hist)*1.1, 0.5)]);

%% Local Function: run_detailed_sim
function [hist, is_valid] = run_detailed_sim(W_roller, b_roller, D_roller, m_eccentric, f, v_sim, mass_ratio, m_rov, target_relative_density, bounce_margin, h_layer, SOIL, eta_mech, nu, rho_min_lunar, rho_max_lunar, rho_i, Apad, n_rollers, c_f, g_moon, P_hotel, battery_density_Wh_kg, max_battery_fraction, t_work_cycle_h, max_sinkage_ratio, DRIVE)
    % Initializations
    R_roller = D_roller / 2;
    omega = 2 * pi * f;
    m_total_wheel = m_rov / n_rollers;
    m_drum = m_total_wheel * mass_ratio;
    F0 = m_eccentric * omega^2;
    rho_f = rho_min_lunar + target_relative_density * (rho_max_lunar - rho_min_lunar);
    rho = rho_i;
    
    % Precomputed constants
    rho_slope = (2.15 - 1.63) / (1.95 - 1.27);
    A_chen_inv = 1 / (pi * 0.151^2);
    E_mult = 6.498e-4;
    nu_term = 1 - nu;
    omega_sq_mdrum_half = (m_drum * omega^2) / 2;
    term2_fixed_part = (F0^2 * m_drum^2) / (m_total_wheel^2);
    xcr_threshold = 0.5 * h_layer;
    is_valid = true;
    
    % Power Model (Vibration)
    KE_peak = (m_eccentric * omega)^2 / (2 * m_drum);
    req_vib_power = (KE_peak * f) / eta_mech;
    
    % History Tracking
    hist.rho_hist = [rho];
    hist.z_hist = [0];
    hist.F_loco_hist = [];
    hist.P_loco_hist = [];
    hist.P_vib_hist = [];
    hist.P_hotel_hist = [];
    hist.d_pass_hist = [];
    hist.t_total = 0;
    
    current_pass = 0;
    z_prev_pass = 0;
    n_total = 0;
    
    while rho < rho_f
        current_pass = current_pass + 1;
        rho_before_pass = rho;
        
        frac = max(0, min(1, (rho - rho_min_lunar) / (rho_max_lunar - rho_min_lunar)));
        kc_dyn = SOIL.kc + frac * (SOIL.kc_dense - SOIL.kc);
        kphi_dyn = SOIL.kphi + frac * (SOIL.kphi_dense - SOIL.kphi);

        N_total_current = current_pass * n_rollers; 
        z_current = calculate_tandem_sinkage(W_roller, b_roller, D_roller, kc_dyn, kphi_dyn, N_total_current);
        h_pass = z_current - z_prev_pass;
        h_eff = max(h_pass, 0.005*R_roller); 
        
        z_safe = max(z_current, 1e-6);
        if z_safe >= (R_roller * max_sinkage_ratio)
            is_valid = false;
            break;
        end
        
        alpha_arg = max(-1, min(1, 1 - 2*z_safe/D_roller));
        alpha = max(acos(alpha_arg), 1e-6);
        l_o = z_safe * tan(pi/4 - SOIL.phi/2)^2;
        
        R_r = (W_roller * n_rollers) * c_f;
        R_c = 0.5 * (kc_dyn + b_roller * kphi_dyn) * z_safe^2 * n_rollers;
        
        term1_b = (b_roller * sin(alpha + SOIL.phi)) / (2 * sin(alpha) * cos(SOIL.phi));
        term2_b = 2 * z_safe * SOIL.c * SOIL.K_c + SOIL.gamma * z_safe^2 * SOIL.K_gamma;
        term3_b = (l_o^3 * SOIL.gamma / 3) * (pi/2 - SOIL.phi);
        term4_b = SOIL.c * l_o^2 * (1 + tan(pi/4 + SOIL.phi/2));
        R_b = term1_b * term2_b + term3_b + term4_b; 
        
        F_loco = R_r + R_c + R_b;
        d_pass = Apad / b_roller;
        t_pass = d_pass / v_sim;
        
        % Power for this pass
        P_loco = F_loco * v_sim; % F * v = Power
        P_vib = req_vib_power;
        P_hotel_val = P_hotel;
        
        [rp_eff, A_col, lc] = rp_fun(R_roller, b_roller, h_eff);
        cycles_in_pass = ceil(f * (lc / v_sim));
        
        for k = 1:cycles_in_pass
            n_total = n_total + 1;
            rho_e = 1.63 + (rho - 1.27) * rho_slope;
            if rho_e < 1.63, rho_e = 1.63; elseif rho_e > 2.15, rho_e = 2.15; end
            Pb = 0.07932 * (100*rho_e - 184)^2 * (A_col * A_chen_inv);
            E_val = E_mult * exp(12.07 * rho_e);
            kcr = 5.686e12 * rho_e^(-41.58) + 0.9079;
            ks = (2 * rp_eff * E_val) / nu_term;
            ksu = ks * kcr;

            if F0 <= Pb
                xcr = 0;
            else
                term1_xcr = omega_sq_mdrum_half / (ks^2 * Pb);
                term2_xcr = term2_fixed_part - Pb^2;
                xcr = term1_xcr * term2_xcr + (Pb/ks - Pb/ksu);
                xcr = max(0, xcr);
            end

            if xcr > xcr_threshold || (F0 <= Pb && rho < rho_f)
                is_valid = false;
                break;
            end

            if xcr > 0
                h_layer_new = h_layer - xcr;
                if h_layer_new <= 1e-6, h_layer_new = 1e-6; end
                rho = rho * (h_layer / h_layer_new);
                h_layer = h_layer_new;
            end
            if rho >= rho_f, break; end
        end
        
        if rho <= rho_before_pass
            is_valid = false;
            break;
        end
        
        % Log pass history
        hist.rho_hist(end+1) = rho;
        hist.z_hist(end+1) = z_current;
        hist.F_loco_hist(end+1) = F_loco;
        hist.P_loco_hist(end+1) = P_loco;
        hist.P_vib_hist(end+1) = P_vib;
        hist.P_hotel_hist(end+1) = P_hotel_val;
        hist.d_pass_hist(end+1) = d_pass;
        hist.t_total = hist.t_total + t_pass;

        if ~is_valid || rho >= rho_f, break; end
        if n_total > 1e6, is_valid = false; break; end 
        z_prev_pass = z_current;
    end
end

%% Helper Functions (Copied from ChenSweep.m)
function z = calculate_tandem_sinkage(W_roller, b, D, kc, kphi, N_total_current)
    z_prev = 0;
    for i = 1:N_total_current
        term_load = (3 * W_roller) / (2 * (kc + b * kphi) * sqrt(D));
        z_i = (term_load + max(0, z_prev)^(1.5))^(2/3); 
        z_prev = z_i;
    end
    z = z_prev;
end

function [rp_eff,A_contact,lc] = rp_fun(r,b,z)
    z = min(z, r * 0.999);
    theta = acos((r-z)/r);
    lc = 2*theta*r;
    A_contact = b * lc;
    rp_eff = sqrt(A_contact/pi); 
end
