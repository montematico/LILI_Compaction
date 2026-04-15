clear; clc; close all;
clear updateProgress;

%% 1. Constants and Setup
g_moon   = 1.62;      % m/s^2
Apad     = 500;       % m^2
Tavail_h = 1000;    % hours, total allowed compaction time for 1 pad - 6mon
energy_max_kWh = 400;
target_relative_density = 0.85;
% roller_fraction = 0.3;
n_rollers   = 2;
bounce_margin = 1.8; %apparently this is fine
h_layer  = 0.12; %depth of layer being compacted
v_sim = 0.05;
max_sinkage_ratio = 1.0; % Max allowable sinkage as a fraction of roller radius
save_figures_to_disk = false; %Save .fig files
dry_run = false; % Set to true to run only 10 combinations for testing
P_hotel = 80; % Watts, baseline power
battery_density_Wh_kg = 150; % Wh/kg
max_battery_fraction = 0.8; % Max mass fraction for battery
t_work_cycle_h = 6; % hours of work before recharge
t_charge_cycle_h = 6; % hours to recharge

% Soil properties (from CompactionChen.m)
SOIL.rho   = 1600;   % kg/m^3
SOIL.gamma = 2470;   % N/m^3
SOIL.n     = 1;
SOIL.kc    = 1400;   % N/m^2
SOIL.kphi  = 830000; % N/m^3
SOIL.phi   = 0.576;  % rad
SOIL.c     = 170;    % N/m^2

% Locomotion and Dense Soil Parameters
c_f = 0.05;                % Coefficient of rolling friction
SOIL.kc_dense   = SOIL.kc*1.5;    % N/m^2 (Assumed 5x "loose" state)
SOIL.kphi_dense = SOIL.kphi*1.5; % N/m^3 (Assumed 5x "loose" state)
SOIL.K_c = 33.37;          % Terzaghi cohesive modulus for lunar soil
SOIL.K_gamma = 72.77;      % Terzaghi frictional modulus for lunar soil

% Other static parameters
eta_mech    = 0.6; %Assumed tractive efficiency
nu            = 0.35; %Regolith Poisson
rho_min_lunar = 1.27; %Used for converting RD's
rho_max_lunar = 1.95; % "
rho_i         = 1.3;

% Pareto Optimization Weights [Total Mass, Energy, Mass Ratio, Freq Diff from 50Hz]
pareto_weights = [1.0, 0.5, 0.0, 0.6]; 

% Drive Wheel Locomotion Parameters
DRIVE.Nw = 4;                     % Number of drive wheels
% DRIVE.D = 1.50;                   % Drive wheel diameter [m]
DRIVE.b = 0.25;                   % Drive wheel width [m]
DRIVE.kwheel = 1e6;               % Radial stiffness [N/m]
DRIVE.slip = 0.40;                % Slip ratio [0 to 1]
DRIVE.bulldozing_factor = 0.20;   % Rb = factor * Rc
DRIVE.K = 0.02;                   % Shear deformation modulus [m]

%% 2. Define the 7D Sweep Grid
f_range = linspace(30, 80, 3);
m_eccentric_range = logspace(-4, log10(0.5e-2), 4);
m_rov_range = linspace(20, 100, 10);
mass_ratio_range = linspace(0.2, 0.8, 3);
R_roller_range = linspace(0.07, 0.20, 3);
roller_fraction_range = linspace(0.2, 0.6, 4);
D_wheel_range = linspace(0.3, 0.8, 3);

[F_GRID, M_ECC_GRID, M_ROV_GRID, M_RATIO_GRID, R_ROLLER_GRID, ROL_FRAC_GRID, D_WHEEL_GRID] = ndgrid(f_range, m_eccentric_range, m_rov_range, mass_ratio_range, R_roller_range, roller_fraction_range, D_wheel_range);

% Initialize results arrays
Rho_results = nan(size(F_GRID));
Power_results = nan(size(F_GRID));
Passes_results = nan(size(F_GRID));
Cycles_results = nan(size(F_GRID));
Energy_results = nan(size(F_GRID));
Time_results = nan(size(F_GRID));
Max_Traction_results = nan(size(F_GRID));
Valid_mask = false(size(F_GRID));

if dry_run
    fprintf('DRY RUN ENABLED: Slicing grid to first 10 combinations.\n');
    F_GRID = F_GRID(1:min(10, end));
    M_ECC_GRID = M_ECC_GRID(1:min(10, end));
    M_ROV_GRID = M_ROV_GRID(1:min(10, end));
    M_RATIO_GRID = M_RATIO_GRID(1:min(10, end));
    R_ROLLER_GRID = R_ROLLER_GRID(1:min(10, end));
    ROL_FRAC_GRID = ROL_FRAC_GRID(1:min(10, end));
    D_WHEEL_GRID = D_WHEEL_GRID(1:min(10, end));
    
    % Re-initialize results arrays for sliced grid
    Rho_results = nan(size(F_GRID));
    Power_results = nan(size(F_GRID));
    Passes_results = nan(size(F_GRID));
    Cycles_results = nan(size(F_GRID));
    Energy_results = nan(size(F_GRID));
    Time_results = nan(size(F_GRID));
    Max_Traction_results = nan(size(F_GRID));
    Valid_mask = false(size(F_GRID));
end

% Parpool setup
% poolobj = gcp('nocreate'); 
% if isempty(poolobj)
%     localCluster = parcluster('local'); 
% 
%     % Request exactly 14 independent process workers.
%     % This maps 1-to-1 with your physical cores and prevents RAM exhaustion.
%     num_workers = min(14, localCluster.NumWorkers); 
% 
%     fprintf('Starting optimized process pool with %d workers...\n', num_workers);
%     parpool(localCluster, num_workers); 
% end
% 1. Set up the DataQueue
D = parallel.pool.DataQueue;
% 2. Define total iterations for percentage calculation
total_iters = numel(F_GRID);
% 3. Use an anonymous function with a persistent counter to track progress
afterEach(D, @(x) updateProgress(total_iters));

%% 3. The Main Sweep Loop
fprintf('Running simulation for %d combinations...\n', numel(F_GRID));
success_count = 0;

tic
parfor i = 1:numel(F_GRID)
    % Extract current iteration's parameters
    f = F_GRID(i);
    m_eccentric = M_ECC_GRID(i);
    m_rov = M_ROV_GRID(i);
    mass_ratio = M_RATIO_GRID(i);
    R_roller = R_ROLLER_GRID(i);
    roller_fraction = ROL_FRAC_GRID(i);
    
    loop_DRIVE = DRIVE;
    loop_DRIVE.D = D_WHEEL_GRID(i);
    
    % Derived parameters
    D_roller = R_roller * 2;
    b_roller = 0.3; % Fixed width
    W_roller = (m_rov * g_moon * roller_fraction) / n_rollers;
    
    % Call the simulation function (Updated with battery and hotel power)
    [rho_final, total_avg_power_W, total_passes, total_cycles, lc_avg_dynamic, total_energy_kWh, is_valid, max_F_loco] = ...
        run_compaction_sim(W_roller, b_roller, D_roller, m_eccentric, f, v_sim, mass_ratio, m_rov, ...
        target_relative_density, bounce_margin, h_layer, SOIL, eta_mech, nu, rho_min_lunar, rho_max_lunar, rho_i, ...
        Apad, n_rollers, c_f, g_moon, P_hotel, battery_density_Wh_kg, max_battery_fraction, t_work_cycle_h, max_sinkage_ratio, loop_DRIVE);
    
    %Time and power sanity check
    if is_valid
        % 1. Calculate final pad mission time including recharge cycles
        distance_to_cover_pad = Apad / b_roller; 
        time_per_pass_seconds = distance_to_cover_pad / v_sim;
        total_active_time_hours = (total_passes * time_per_pass_seconds) / 3600;
        
        charge_cycles = max(0, ceil(total_active_time_hours / t_work_cycle_h) - 1);
        total_time_hours = total_active_time_hours + (charge_cycles * t_charge_cycle_h);
        
        % 2. Apply Sanity Filters
        if total_time_hours > Tavail_h || total_energy_kWh > energy_max_kWh
            is_valid = false; % Discard if it takes too long or uses too much energy
        end
    end

    if is_valid && ~isnan(rho_final)
        % Pad Energy & Time Calculation (Repeat logic for storage)
        distance_to_cover_pad = Apad / b_roller; 
        time_per_pass_seconds = distance_to_cover_pad / v_sim;
        total_active_time_hours = (total_passes * time_per_pass_seconds) / 3600;
        charge_cycles = max(0, ceil(total_active_time_hours / t_work_cycle_h) - 1);
        total_time_hours = total_active_time_hours + (charge_cycles * t_charge_cycle_h);
        
        % Store successful results
        Rho_results(i) = rho_final;
        Power_results(i) = total_avg_power_W;
        Passes_results(i) = total_passes;
        Cycles_results(i) = total_cycles;
        Energy_results(i) = total_energy_kWh;
        Time_results(i) = total_time_hours;
        Max_Traction_results(i) = max_F_loco;
        Valid_mask(i) = true;
        success_count = success_count + 1;
    end
    %Intermittent Status updates
    if mod(i, 50) == 0
        send(D,i);
    end
end
fprintf('Sweep complete. Successes: %d\n\n', success_count);
toc
%% 4. Pareto Optimization
valid_idx = find(Valid_mask);
if isempty(valid_idx)
    fprintf('No successful points found. Cannot perform Pareto analysis.\n');
else
    % Extract valid vectors
    m_rov_vec = M_ROV_GRID(valid_idx);
    energy_vec = Energy_results(valid_idx);
    mass_ratio_vec = M_RATIO_GRID(valid_idx);
    f_vec = F_GRID(valid_idx);
    
    % 4-Objective Matrix: [Minimize Mass, Minimize Energy, Maximize Ratio, Target 50Hz]
    objectives = [m_rov_vec, energy_vec, 1./mass_ratio_vec, abs(f_vec - 50)];
    
    % Find the non-dominated set (Pareto Front)
    is_pareto = true(length(valid_idx), 1);
    for i = 1:length(valid_idx)
        is_dominated_by = any(all(objectives <= objectives(i,:), 2) & any(objectives < objectives(i,:), 2));
        if is_dominated_by
            is_pareto(i) = false;
        end
    end
    
    pareto_idx = valid_idx(is_pareto);
    pareto_objectives = objectives(is_pareto, :);
    
    fprintf('Found %d Pareto-optimal points.\n', length(pareto_idx));

    % 1. Calculate Area Rates
    b_roller = 0.3; % Fixed width
    net_area_rate = v_sim * b_roller;
    gross_area_rates = Apad ./ (Time_results(pareto_idx) * 3600);
    
    % 2. Extract Pareto vectors for stats
    pareto_max_traction = Max_Traction_results(pareto_idx);
    pareto_radius = R_ROLLER_GRID(pareto_idx);
    pareto_m_rov = M_ROV_GRID(pareto_idx);
    pareto_mass_ratios = M_RATIO_GRID(pareto_idx);
    pareto_rol_frac = ROL_FRAC_GRID(pareto_idx);
    pareto_d_wheel = D_WHEEL_GRID(pareto_idx);

    % 3. Print Summary Table
    fprintf('\n================ PARETO OPTIMAL SUMMARY (N=%d) ================\n', length(pareto_idx));
    fprintf('%-25s | %-10s | %-10s | %-10s\n', 'Variable', 'Mean', 'Min', 'Max');
    fprintf('---------------------------------------------------------------\n');
    fprintf('%-25s | %-10.2f | %-10.2f | %-10.2f\n', 'Max Traction Force (N)', mean(pareto_max_traction), min(pareto_max_traction), max(pareto_max_traction));
    fprintf('%-25s | %-10.3f | %-10.3f | %-10.3f\n', 'Roller Radius (m)', mean(pareto_radius), min(pareto_radius), max(pareto_radius));
    fprintf('%-25s | %-10.1f | %-10.1f | %-10.1f\n', 'Total Rover Mass (kg)', mean(pareto_m_rov), min(pareto_m_rov), max(pareto_m_rov));
    fprintf('%-25s | %-10.2f | %-10.2f | %-10.2f\n', 'Mass Ratio', mean(pareto_mass_ratios), min(pareto_mass_ratios), max(pareto_mass_ratios));
    fprintf('%-25s | %-10.4f | %-10.4f | %-10.4f\n', 'Gross Area Rate (m^2/s)', mean(gross_area_rates), min(gross_area_rates), max(gross_area_rates));
    fprintf('%-25s | %-10.2f | %-10.2f | %-10.2f\n', 'Roller Fraction', mean(pareto_rol_frac), min(pareto_rol_frac), max(pareto_rol_frac));
    fprintf('%-25s | %-10.3f | %-10.3f | %-10.3f\n', 'Drive Wheel Diameter (m)', mean(pareto_d_wheel), min(pareto_d_wheel), max(pareto_d_wheel));
    fprintf('---------------------------------------------------------------\n');
    fprintf('Constant Net Area Rate: %.4f m^2/s\n', net_area_rate);
    fprintf('===============================================================\n\n');

    save('LunarCompactionResults.mat');
    %plot_compaction_results('LunarCompactionResults.mat');
end

%% Functions (Copied from CompactionChen.m)
function [rho_final, total_avg_power_W, total_passes, total_cycles, lc_avg_dynamic, total_energy_kWh, is_valid, max_F_loco] = ...
    run_compaction_sim(W_roller, b_roller, D_roller, m_eccentric, f, v_sim, mass_ratio, m_rov, target_relative_density, bounce_margin, h_layer, SOIL, eta_mech, nu, rho_min_lunar, rho_max_lunar, rho_i, Apad, n_rollers, c_f, g_moon, P_hotel, battery_density_Wh_kg, max_battery_fraction, t_work_cycle_h, max_sinkage_ratio, DRIVE)
    
    % Internal Initializations
    max_F_loco = 0;
    R_roller = D_roller / 2;
    omega = 2 * pi * f;
    
    % --- Physical Mass Derivations ---
    m_total_wheel = m_rov / n_rollers;        % Total mass acting on one wheel
    m_drum = m_total_wheel * mass_ratio;      % Mass of the bouncing steel drum (Chen's mp)
    
    F0 = m_eccentric * omega^2;
    
    % % FAST FAIL: Bounce check 
    % if F0 >= (W_roller * bounce_margin)
    %     is_valid = false; 
    %     rho_final = NaN; total_avg_power_W = NaN; total_passes = NaN; total_cycles = NaN; lc_avg_dynamic = NaN; total_energy_kWh = NaN;
    %     return;
    % end
    
    rho_f = rho_min_lunar + target_relative_density * (rho_max_lunar - rho_min_lunar);
    rho = rho_i; 
    
    % --- PRECOMPUTED CONSTANTS FOR OPTIMIZATION ---
    rho_slope = (2.15 - 1.63) / (1.95 - 1.27);
    A_chen_inv = 1 / (pi * 0.151^2);
    E_mult = 6.498e-4; % Combined 6.498e-10 * 1e6
    nu_term = 1 - nu;
    
    % Updated for new mass variables
    omega_sq_mdrum_half = (m_drum * omega^2) / 2;
    term2_fixed_part = (F0^2 * m_drum^2) / (m_total_wheel^2);
    % ----------------------------------------------

    % Loop trackers
    current_pass = 0;
    z_prev_pass = 0; 
    n_total = 0; %number of cycles
    sum_lc = 0; % for running sum
    
    xcr_threshold = 0.5 * h_layer; %plastic deformation thresh.
    is_valid = true;
    total_energy_Joule = 0; % Track cumulative work done

    % Power Model (Vibration)
    KE_peak = (m_eccentric * omega)^2 / (2 * m_drum);
    req_vib_power = (KE_peak * f) / eta_mech;

    % =========================================================================
    % PASS 1 FAST-BREAK TRACTION CHECK
    % Calculate the resistance the roller will experience on the very first pass 
    % using baseline (loose) soil moduli. If the drive wheels cannot overcome 
    % this initial resistance, discard the design immediately.
    % =========================================================================
    
    % 1. Roller Sinkage & Resistance on Pass 1
    z_pass1 = calculate_tandem_sinkage(W_roller, b_roller, D_roller, SOIL.kc, SOIL.kphi, n_rollers);
    z_safe_p1 = max(z_pass1, 1e-6);
    
    alpha_arg_p1 = max(-1, min(1, 1 - 2*z_safe_p1/D_roller));
    alpha_p1 = max(acos(alpha_arg_p1), 1e-6);
    l_o_p1 = z_safe_p1 * tan(pi/4 - SOIL.phi/2)^2;
    
    R_r_p1 = (W_roller * n_rollers) * c_f;
    R_c_p1 = 0.5 * (SOIL.kc + b_roller * SOIL.kphi) * z_safe_p1^2 * n_rollers;
    
    term1_b_p1 = (b_roller * sin(alpha_p1 + SOIL.phi)) / (2 * sin(alpha_p1) * cos(SOIL.phi));
    term2_b_p1 = 2 * z_safe_p1 * SOIL.c * SOIL.K_c + SOIL.gamma * z_safe_p1^2 * SOIL.K_gamma;
    term3_b_p1 = (l_o_p1^3 * SOIL.gamma / 3) * (pi/2 - SOIL.phi);
    term4_b_p1 = SOIL.c * l_o_p1^2 * (1 + tan(pi/4 + SOIL.phi/2));
    R_b_p1 = term1_b_p1 * term2_b_p1 + term3_b_p1 + term4_b_p1;
    
    F_loco_required_p1 = R_r_p1 + R_c_p1 + R_b_p1;
    
    % 2. Calculate Drive Wheel Capabilities
    % Calculate normal load per drive wheel in Newtons (mass not used by rollers)
    W_drive_wheel = (m_rov * g_moon - W_roller * n_rollers) / DRIVE.Nw;
    
    % Call the external wheel model with configuration parameters
    drive_result = lunar_wheel_model(W_drive_wheel, DRIVE.D, DRIVE.b, DRIVE.kwheel, DRIVE.slip, DRIVE.bulldozing_factor, SOIL.kc, SOIL.kphi, SOIL.n, SOIL.c, rad2deg(SOIL.phi), DRIVE.K);
    
    % Drawbar pull is the excess traction after wheel resistance is subtracted.
    total_available_DBP = drive_result.drawbar_pull_DP_N * DRIVE.Nw;
    
    % 3. Fast Break Evaluation
    if total_available_DBP < F_loco_required_p1
        is_valid = false;
        rho_final = NaN; total_avg_power_W = NaN; total_passes = NaN; total_cycles = NaN; lc_avg_dynamic = NaN; total_energy_kWh = NaN;
        return; % Exit the function early to save computation time
    end
    % =========================================================================

    while rho < rho_f
        current_pass = current_pass + 1;
        rho_before_pass = rho; % Track density before the pass starts
        
        % 1. Interpolate Soil Moduli based on Relative Density
        frac = max(0, min(1, (rho - rho_min_lunar) / (rho_max_lunar - rho_min_lunar)));
        kc_dyn = SOIL.kc + frac * (SOIL.kc_dense - SOIL.kc);
        kphi_dyn = SOIL.kphi + frac * (SOIL.kphi_dense - SOIL.kphi);

        % 2. Calculate TOTAL sinkage with dynamic moduli
        N_total_current = current_pass * n_rollers; 
        z_current = calculate_tandem_sinkage(W_roller, b_roller, D_roller, kc_dyn, kphi_dyn, N_total_current);
        h_pass = z_current - z_prev_pass;
        h_eff = max(h_pass, 0.005*R_roller); 
        
        % 3. Calculate Locomotion Resistances for this pass
        z_safe = max(z_current, 1e-6); % Prevent div-by-zero in bulldozing
        if z_safe >= (R_roller * max_sinkage_ratio)
            is_valid = false;
            break;
        end
        
        % Clamp the argument to [-1, 1] to prevent imaginary angles if sinkage exceeds wheel radius
        alpha_arg = max(-1, min(1, 1 - 2*z_safe/D_roller));
        alpha = max(acos(alpha_arg), 1e-6);
        
        l_o = z_safe * tan(pi/4 - SOIL.phi/2)^2;
        
        R_r = (W_roller * n_rollers) * c_f; % Rolling resistance
        R_c = 0.5 * (kc_dyn + b_roller * kphi_dyn) * z_safe^2 * n_rollers; % Compression resistance
        
        % Bulldozing resistance (applied to leading roller only)
        term1_b = (b_roller * sin(alpha + SOIL.phi)) / (2 * sin(alpha) * cos(SOIL.phi));
        term2_b = 2 * z_safe * SOIL.c * SOIL.K_c + SOIL.gamma * z_safe^2 * SOIL.K_gamma;
        term3_b = (l_o^3 * SOIL.gamma / 3) * (pi/2 - SOIL.phi);
        term4_b = SOIL.c * l_o^2 * (1 + tan(pi/4 + SOIL.phi/2));
        R_b = term1_b * term2_b + term3_b + term4_b; 
        
        F_loco = R_r + R_c + R_b; % Total forward force required
        max_F_loco = max(max_F_loco, F_loco);
        
        % Calculate Work for this pass
        d_pass = Apad / b_roller; % Distance traveled per pass
        work_loco_pass = F_loco * d_pass; % Joules
        
        t_pass = d_pass / v_sim; % Seconds driving
        work_vib_pass = req_vib_power * t_pass; % Joules
        work_hotel_pass = P_hotel * t_pass; % Joules
        
        total_energy_Joule = total_energy_Joule + work_loco_pass + work_vib_pass + work_hotel_pass;
        
        % 4. Recalculate rp_eff based on incremental sinkage
        [rp_eff, A_col, lc] = rp_fun(R_roller, b_roller, h_eff);
        cycles_in_pass = ceil(f * (lc / v_sim));
        
       
        for k = 1:cycles_in_pass
            n_total = n_total + 1;
            sum_lc = sum_lc + lc;
            
            % --- INLINED PHYSICS ---
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
            % -----------------------

            if xcr > xcr_threshold
                is_valid = false;
                break;
            elseif F0 <= Pb && rho < rho_f
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
        
        % Fail early if this pass made no progress
        if rho <= rho_before_pass
            is_valid = false;
            break;
        end

        if ~is_valid || rho >= rho_f, break; end
        if n_total > 1e6, is_valid = false; break; end 
        z_prev_pass = z_current;
    end

    % --- BATTERY MASS SANITY CHECK ---
    if is_valid
        total_active_time_sec = current_pass * (Apad / b_roller) / v_sim;
        total_avg_power_W = total_energy_Joule / total_active_time_sec;
        
        % Size battery mass based on work cycle
        energy_per_cycle_Wh = total_avg_power_W * t_work_cycle_h;
        m_battery = energy_per_cycle_Wh / battery_density_Wh_kg; % kg
        
        m_non_roller = m_rov - (W_roller * n_rollers / g_moon);
        max_m_battery = m_non_roller * max_battery_fraction;
        
        if m_battery > max_m_battery
            is_valid = false;
        end
    end

    % NaN Mapping for Failures
    if is_valid
        rho_final = rho;
        total_passes = current_pass;
        total_cycles = n_total;
        lc_avg_dynamic = sum_lc / n_total;
        total_energy_kWh = total_energy_Joule / 3.6e6; % Convert J to kWh
    else
        rho_final = NaN;
        total_avg_power_W = NaN;
        total_passes = NaN;
        total_cycles = NaN;
        lc_avg_dynamic = NaN;
        total_energy_kWh = NaN;
        max_F_loco = NaN;
    end
end

function xcr = chen_residual_deformation(F0, omega, ks, ksu, Pb, m_drum, m_total_wheel)
    if F0 <= Pb
        xcr = 0;
        return;
    end
    term1 = (m_drum * omega^2) / (2 * ks^2 * Pb);
    term2 = (F0^2 * m_drum^2) / (m_total_wheel^2) - Pb^2;
    term3 = Pb/ks - Pb/ksu;
    xcr = term1 * term2 + term3;
    xcr = max(0, xcr);
end

function [rp_eff,A_contact,lc] = rp_fun(r,b,z)
    z = min(z, r * 0.999);
    theta = acos((r-z)/r); %half angle
    lc = 2*theta*r; %full arc length of contact patch
    A_contact = b * lc;               % m^2, contact area
    rp_eff = sqrt(A_contact/pi); 
end

function z = calculate_tandem_sinkage(W_roller, b, D, kc, kphi, N_total_current)
    z_prev = 0;
    for i = 1:N_total_current
        term_load = (3 * W_roller) / (2 * (kc + b * kphi) * sqrt(D));
        % Max(0, ...) prevents complex numbers if floating-point drift pushes z_prev slightly negative
        z_i = (term_load + max(0, z_prev)^(1.5))^(2/3); 
        z_prev = z_i;
    end
    z = z_prev;
end

function rho_earth = translate_density(rho_lunar)
    slope = (2.15 - 1.63) / (1.95 - 1.27);
    rho_earth = 1.63 + (rho_lunar - 1.27) * slope;
    rho_earth = max(1.63, min(2.15, rho_earth));
end

%Soil density function
    % Soil moduli functions (Earth-mapped density rho_e)
function E = E_fun(rho_e); E = 6.498e-10 .* 1e6 .* exp(12.07.*rho_e); end
function k_cr = k_cr_fun(rho_e); k_cr = 5.686e12 * rho_e^(-41.58) + 0.9079; end
function ks = ks_fun(rho_e,rp,nu); ks = 2*rp * E_fun(rho_e) / (1-nu); end
function ksu = ksu_fun(rho_e,rp,kcr,nu); ksu = 2*rp * (kcr * E_fun(rho_e)) / (1-nu); end


function updateProgress(total)
    persistent count
    if isempty(count)
        count = 0;
    end
    count = count + 50;
    pct = (count / total) * 100;
    % Use \r to overwrite the line for a cleaner look
    fprintf('Progress: %.2f%% (%d of %d)\n', pct, count, total);
end