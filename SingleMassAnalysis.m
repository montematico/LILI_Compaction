%% SingleMassAnalysis.m
% Deep-dive engineering analysis on a single, optimal lunar compaction rover design.
% Based on grid sweep results and pass-by-pass simulation.

clear; clc; close all;

%% 1. Parameters & Data Loading
target_idx = 12065; % Define target index at the top

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

% Add Variable Fallbacks
if ~exist('n_rollers', 'var'), n_rollers = 2; end
if ~exist('b_roller', 'var'), b_roller = 0.3; end

% Add Mission Parameters & Grousers
h_g = 0.008; % m - Height of grousers on drive wheels
n_g = 18;    % dimensionless - Number of grousers per wheel
slip_target = 0.20; % dimensionless - Target slip for traction check
SOIL.K = 0.02; % m - Shear deformation modulus

W_roller_opt = (m_opt * g_moon * roller_fraction_opt) / n_rollers; 

% Update DRIVE struct for this iteration
DRIVE_opt = DRIVE;
DRIVE_opt.D = D_wheel_opt;
DRIVE_opt.h_g = h_g;
DRIVE_opt.n_g = n_g;
DRIVE_opt.slip = slip_target;

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

%% Core Task 3: Print Specification Table
max_sinkage_mm = max(hist.z_hist) * 1000;
max_F_loco = max(hist.F_loco_hist);
total_passes = length(hist.rho_hist) - 1;
total_time_hr = hist.t_total / 3600;

% Calculate Final Wheel Specs for Table
W_drive_wheel = (m_opt * g_moon - W_roller_opt * n_rollers) / DRIVE_opt.Nw;
final_perf = calculate_wheel_performance(W_drive_wheel, DRIVE_opt.D, DRIVE_opt.b, DRIVE_opt.h_g, DRIVE_opt.n_g, DRIVE_opt.slip, v_sim, SOIL);

fprintf('\n========================================================\n');
fprintf('         LUNAR COMPACTION ROVER SPECIFICATIONS          \n');
fprintf('========================================================\n');

fprintf('\n--- Physical & Geometric ---\n');
fprintf('%-30s : %8.2f kg\n', 'Total Rover Mass', m_opt);
fprintf('%-30s : %8.2f kg\n', 'Static Load per Roller', W_roller_opt/g_moon);
fprintf('%-30s : %8.2f m\n', 'Roller Diameter', D_roller_opt);
fprintf('%-30s : %8.2f m\n', 'Drive Wheel Diameter', D_wheel_opt);
fprintf('%-30s : %8.2f %%\n', 'Roller Mass Fraction', roller_fraction_opt * 100);
fprintf('%-30s : %8.2f mm\n', 'Grouser Height', h_g * 1000);

fprintf('\n--- Vibration System ---\n');
fprintf('%-30s : %8.2f Hz\n', 'Operating Frequency', f_opt);
fprintf('%-30s : %8.2f kg\n', 'Eccentric Mass/Unbalance', m_eccentric_opt);
fprintf('%-30s : %8.2f\n', 'Dynamic Drum Mass Ratio', mass_ratio_opt);

fprintf('\n--- Performance & Mission ---\n');
fprintf('%-30s : %8d\n', 'Total Passes', total_passes);
fprintf('%-30s : %8.2f hr\n', 'Total Mission Time', total_time_hr);
fprintf('%-30s : %8.2f mm\n', 'Maximum Sinkage', max_sinkage_mm);
fprintf('%-30s : %8.2f N\n', 'Max Roller Resistance', max_F_loco);
fprintf('%-30s : %8.3f\n', 'Max Traction Coeff (MU_max)', final_perf.MU_max);
fprintf('%-30s : %8.2f N\n', 'Net Drawbar Pull (DBP)', final_perf.DBP * DRIVE_opt.Nw);

fprintf('\n--- Power & Energy ---\n');
fprintf('%-30s : %8.2f W\n', 'Max Locomotion Power', max(hist.P_loco_hist));
fprintf('%-30s : %8.2f W\n', 'Max Vibration Power', max(hist.P_vib_hist));
fprintf('%-30s : %8.2f W\n', 'Peak Total System Power', P_peak);
fprintf('%-30s : %8.2f kWh\n', 'Total Energy Consumption', total_energy_kWh);
fprintf('%-30s : %8.2f J/kg\n', 'Specific Energy', SE);
fprintf('%-30s : %8.2f kg\n', 'Required Battery Mass', m_battery_req);
fprintf('%-30s : %8.2f C\n', 'Peak C-Rate', C_rate_peak);
fprintf('========================================================\n');

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
        
        % Calculate Wheel Performance & Power
        W_drive_wheel = (m_rov * g_moon - W_roller * n_rollers) / DRIVE.Nw;
        % The wheels must overcome roller resistance PLUS their own resistance
        % We call the model with zero slip to get the "required" effort, 
        % or more accurately, we calculate resistance and then add it to roller load.
        perf = calculate_wheel_performance(W_drive_wheel, DRIVE.D, DRIVE.b, DRIVE.h_g, DRIVE.n_g, DRIVE.slip, v_sim, SOIL);
        
        % Total Power = (Roller Force + Wheel Resistance) * v / eta
        P_loco = (F_loco + perf.R_wheel * DRIVE.Nw) * v_sim / (0.85 * 0.9);
        
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
