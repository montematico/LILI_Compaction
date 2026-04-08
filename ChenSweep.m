clear; clc;

%% 1. Constants and Setup
g_moon   = 1.62;      % m/s^2
Apad     = 500;       % m^2
Tavail_h = 365*24;    % hours, total allowed compaction time
target_relative_density = 0.85;
roller_fraction = 0.5;
n_rollers   = 1;
bounce_margin = 1.0;
h_layer  = 0.1;
v_sim = 0.05;

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
SOIL.kc_dense   = SOIL.kc*5;    % N/m^2 (Assumed 5x "loose" state)
SOIL.kphi_dense = SOIL.phi*5; % N/m^3 (Assumed 5x "loose" state)
SOIL.K_c = 33.37;          % Terzaghi cohesive modulus for lunar soil
SOIL.K_gamma = 72.77;      % Terzaghi frictional modulus for lunar soil

% Other static parameters
eta_mech    = 0.6; %Assumed tractive efficiency
nu            = 0.35; %Regolith Poisson
rho_min_lunar = 1.27; %Used for converting RD's
rho_max_lunar = 1.95; % "
rho_i         = 1.3;

% Pareto Optimization Weights [Total Mass, Energy, Mass Ratio, Freq Diff from 50Hz]
pareto_weights = [1.0, 0.5, 0.0, 0.0]; 

%% 2. Define the 5D Sweep Grid
f_range = linspace(20, 50, 5);
m_eccentric_range = logspace(-5, -1, 5);
m_rov_range = linspace(20, 100, 5);
mass_ratio_range = linspace(0.2, 0.8, 5);
R_roller_range = linspace(0.15, 0.4, 5);

[F_grid, M_ECC_grid, M_ROV_grid, M_RATIO_grid, R_ROLLER_grid] = ndgrid(f_range, m_eccentric_range, m_rov_range, mass_ratio_range, R_roller_range);

% Initialize results arrays
Rho_results = nan(size(F_grid));
Power_results = nan(size(F_grid));
Passes_results = nan(size(F_grid));
Cycles_results = nan(size(F_grid));
Energy_results = nan(size(F_grid));
Time_results = nan(size(F_grid));
Valid_mask = false(size(F_grid));

%% 3. The Main Sweep Loop
fprintf('Running simulation for %d combinations...\n', numel(F_grid));
success_count = 0;

% Parpool setup
poolobj = gcp('nocreate'); % Check if a pool already exists
if isempty(poolobj)
    parpool("Threads",[6,20]);    % Create a thread-based pool only if none exists
end

tic
parfor i = 1:numel(F_grid)
    % Extract current iteration's parameters
    f = F_grid(i);
    m_eccentric = M_ECC_grid(i);
    m_rov = M_ROV_grid(i);
    mass_ratio = M_RATIO_grid(i);
    m_roller = M_RATIO_grid(i) * M_RATIO_grid(i) / n_rollers; 
    R_roller = R_ROLLER_grid(i);
    
    % Derived parameters
    D_roller = R_roller * 2;
    b_roller = 0.3; % Fixed width
    W_roller = (m_rov * g_moon * roller_fraction) / n_rollers;
    
    % Call the simulation function
    [rho_final, req_power, total_passes, total_cycles, lc_avg_dynamic, total_energy_kWh, is_valid] = ...
        run_compaction_sim(W_roller, b_roller, D_roller, m_eccentric, f, v_sim, mass_ratio, m_rov, target_relative_density, bounce_margin, h_layer, SOIL, eta_mech, m_roller, nu, rho_min_lunar, rho_max_lunar, rho_i, Apad, n_rollers, c_f, g_moon);
    
    if is_valid && ~isnan(rho_final)
        % Pad Energy & Time Calculation
        distance_to_cover_pad = Apad / b_roller; 
        time_per_pass_seconds = distance_to_cover_pad / v_sim;
        total_time_hours = (total_passes * time_per_pass_seconds) / 3600;
        
        % Store successful results
        Rho_results(i) = rho_final;
        Power_results(i) = req_power;
        Passes_results(i) = total_passes;
        Cycles_results(i) = total_cycles;
        Energy_results(i) = total_energy_kWh;
        Time_results(i) = total_time_hours;
        Valid_mask(i) = true;
        success_count = success_count + 1;
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
    m_rov_vec = M_ROV_grid(valid_idx);
    energy_vec = Energy_results(valid_idx);
    mass_ratio_vec = M_RATIO_grid(valid_idx);
    f_vec = F_grid(valid_idx);
    
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

    %% 5. Graphics and Plots
    
    % Figure 1: 2D Pareto Front Scatter
    figure('Name', 'Pareto Front: Mass vs. Pad Energy');
    pareto_m_rov = M_ROV_grid(pareto_idx);
    pareto_energy = Energy_results(pareto_idx);
    pareto_mass_ratios = M_RATIO_grid(pareto_idx);
    
    scatter(pareto_m_rov, pareto_energy, 80, pareto_mass_ratios, 'filled');
    
    xlabel('Total Rover Mass (kg)');
    ylabel('Pad Energy (kWh)');
    grid on;
    title('Pareto Front: Total Mass vs. Pad Energy');
    colormap(flipud(parula));
    h = colorbar;
    ylabel(h, 'Mass Ratio');

    % Figure 2: Feasibility Map (Frequency vs Mass)
    figure('Name', 'Pareto Optimal Configurations');
    plot(F_grid(pareto_idx), M_ROV_grid(pareto_idx), '^b', 'MarkerSize', 8, 'MarkerFaceColor', 'blue');
    xlabel('Frequency (Hz)');
    ylabel('Total Rover Mass (kg)');
    title('Pareto Optimal Configurations in Frequency vs Mass Space');
    grid on;

    % Figure 3: Pareto Table
    figure('Name', 'Pareto Optimal Set');
    tbl_data_points = [
        M_ROV_grid(pareto_idx), ...
        M_RATIO_grid(pareto_idx), ...
        R_ROLLER_grid(pareto_idx), ...
        F_grid(pareto_idx), ...
        M_ECC_grid(pareto_idx), ...
        Time_results(pareto_idx), ...
        Energy_results(pareto_idx), ...
        Power_results(pareto_idx)
    ];
    
    [~, sort_order] = sort(tbl_data_points(:,1), 'ascend'); % Sort by Total Mass
    tbl_data = tbl_data_points(sort_order, :);
    
    cnames = {'Total Mass (kg)', 'Mass Ratio', 'Radius (m)', 'Freq (Hz)', 'Ecc Mom (kg-m)', 'Pad Time (hr)', 'Energy (kWh)', 'Total Power (W)'};
    t = uitable('Data', tbl_data, 'ColumnName', cnames, 'RowName',[], 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);

    % Figure 4: Parallel Coordinates Plot
    figure('Name', 'Parallel Coordinates: Trade Space', 'Position', [100, 100, 1000, 500]);
    pareto_table = table(F_grid(pareto_idx), M_ECC_grid(pareto_idx), M_ROV_grid(pareto_idx), R_ROLLER_grid(pareto_idx), M_RATIO_grid(pareto_idx), Energy_results(pareto_idx), ...
        'VariableNames', {'Freq_Hz', 'Ecc_Moment', 'Rover_Mass', 'Radius_m', 'Mass_Ratio', 'Energy_kWh'});
    
    p = parallelplot(pareto_table);
    try
        p.ColorVariable = 'Energy_kWh';
        p.LineAlpha = 0.4; 
    catch
    end
    title('Parallel Coordinates: Pareto Optimal Configurations');

    % Figure 5: 3D Pareto Scatter (explicitly showing m_eccentric)
    figure('Name', '3D Pareto Front: Mass, Energy, and Eccentric Moment');
    pareto_m_ecc = M_ECC_grid(pareto_idx);
    scatter3(pareto_m_rov, pareto_energy, pareto_m_ecc, 80, pareto_mass_ratios, 'filled');
    xlabel('Total Rover Mass (kg)');
    ylabel('Pad Energy (kWh)');
    zlabel('Eccentric Moment (kg-m)');
    grid on;
    set(gca, 'ZScale', 'log'); % Fixes negative axis and scales eccentric moment properly
    title('3D Pareto Front: Mass, Energy, and Eccentric Moment');
    h3 = colorbar;
    ylabel(h3, 'Mass Ratio');
    view(45, 30);
end

%% Functions (Copied from CompactionChen.m)
function [rho_final, req_power, total_passes, total_cycles, lc_avg_dynamic, total_energy_kWh, is_valid] = ...
    run_compaction_sim(W_roller, b_roller, D_roller, m_eccentric, f, v_sim, mass_ratio, m_rov, target_relative_density, bounce_margin, h_layer, SOIL, eta_mech, m_vibrator, nu, rho_min_lunar, rho_max_lunar, rho_i, Apad, n_rollers, c_f, g_moon)
    
    % Internal Initializations
    R_roller = D_roller / 2;
    omega = 2 * pi * f;
    m_roller_dynamic = (m_rov / 2) * mass_ratio;
    F0 = m_eccentric * omega^2;
    
    rho_f = rho_min_lunar + target_relative_density * (rho_max_lunar - rho_min_lunar);
    rho = rho_i; % Set active tracking variable to the passed initial density
    
    % Soil moduli functions (Earth-mapped density rho_e)
    E_fun = @(rho_e) 6.498e-10 .* 1e6 .* exp(12.07.*rho_e);
    k_cr_fun =  @(rho_e) 5.686e12 * rho_e^(-41.58) + 0.9079;
    ks_fun = @(rho_e,rp) 2*rp * E_fun(rho_e) / (1-nu);
    ksu_fun = @(rho_e,rp,kcr) 2*rp * (kcr * E_fun(rho_e)) / (1-nu);

    % Loop trackers
    current_pass = 0;
    z_prev_pass = 0;
    n_total = 0;
    lc_history = [];
    xcr_threshold = 0.15 * h_layer;
    is_valid = true;
    req_power = 0;
    total_energy_Joule = 0; % Track cumulative work done

    while rho < rho_f
        current_pass = current_pass + 1;
        
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
        
        % Clamp the argument to [-1, 1] to prevent imaginary angles if sinkage exceeds wheel radius
        alpha_arg = max(-1, min(1, 1 - 2*z_safe/D_roller));
        alpha = max(acos(alpha_arg), 1e-6);
        
        l_o = z_safe * tan(pi/4 - SOIL.phi/2)^2;
        
        R_r = m_rov * g_moon * c_f; % Rolling resistance
        R_c = 0.5 * (kc_dyn + b_roller * kphi_dyn) * z_safe^2 * n_rollers; % Compression resistance
        
        % Bulldozing resistance (applied to leading roller only)
        term1_b = (b_roller * sin(alpha + SOIL.phi)) / (2 * sin(alpha) * cos(SOIL.phi));
        term2_b = 2 * z_safe * SOIL.c * SOIL.K_c + SOIL.gamma * z_safe^2 * SOIL.K_gamma;
        term3_b = (l_o^3 * SOIL.gamma / 3) * (pi/2 - SOIL.phi);
        term4_b = SOIL.c * l_o^2 * (1 + tan(pi/4 + SOIL.phi/2));
        R_b = term1_b * term2_b + term3_b + term4_b; 
        
        F_loco = R_r + R_c + R_b; % Total forward force required
        
        % Calculate Work for this pass
        d_pass = Apad / (b_roller * n_rollers); % Distance traveled per pass
        work_loco_pass = F_loco * d_pass; % Joules
        total_energy_Joule = total_energy_Joule + work_loco_pass;
        
        % 4. Recalculate rp_eff based on incremental sinkage
        [rp_eff, A_col, lc] = rp_fun(R_roller, b_roller, h_eff);
        cycles_in_pass = ceil(f * (lc / v_sim));
        
        for k = 1:cycles_in_pass
            n_total = n_total + 1;
            lc_history(n_total) = lc;
            
            rho_e = translate_density(rho);
            
            Pb0 = 0.07932 * (100*rho_e - 184)^2;
            r_chen = 0.151;
            A_chen = pi * r_chen^2;
            Pb = Pb0 * (A_col/A_chen);

            kcr = k_cr_fun(rho_e);
            ks = ks_fun(rho_e, rp_eff);
            ksu = ksu_fun(rho_e, rp_eff, kcr);

            xcr = chen_residual_deformation(F0, omega, ks, ksu, Pb, m_vibrator, m_roller_dynamic);
            
            % Sanity Checks
            if xcr > xcr_threshold
                is_valid = false;
                break;
            elseif F0 >= (W_roller * bounce_margin)
                is_valid = false;
                break;
            elseif F0 <= Pb && rho < rho_f
                is_valid = false;
                break;
            end

            % Power Model
            KE_peak = (m_eccentric * omega)^2 / (2 * m_roller_dynamic);
            req_power = (KE_peak * f) / eta_mech;

            if xcr > 0
                h_layer_new = h_layer - xcr;
                if h_layer_new <= 1e-6, h_layer_new = 1e-6; end
                rho = rho * (h_layer / h_layer_new);
                h_layer = h_layer_new;
            end
            
            if rho >= rho_f, break; end
        end
        
        % Add vibratory work for the pass
        t_pass = d_pass / v_sim; % Seconds driving
        work_vib_pass = req_power * t_pass; % Joules
        total_energy_Joule = total_energy_Joule + work_vib_pass;

        if ~is_valid || rho >= rho_f, break; end
        if n_total > 1e6, is_valid = false; break; end 
        z_prev_pass = z_current;
    end

    % NaN Mapping for Failures
    if is_valid
        rho_final = rho;
        total_passes = current_pass;
        total_cycles = n_total;
        lc_avg_dynamic = mean(lc_history);
        total_energy_kWh = total_energy_Joule / 3.6e6; % Convert J to kWh
    else
        rho_final = NaN;
        req_power = NaN;
        total_passes = NaN;
        total_cycles = NaN;
        lc_avg_dynamic = NaN;
        total_energy_kWh = NaN;
    end
end

function xcr = chen_residual_deformation(F0, omega, ks, ksu, Pb, m_vibrator, m_roller_dynamic)
    if F0 <= Pb
        xcr = 0;
        return;
    end
    term1 = (m_vibrator * omega^2) / (2 * ks^2 * Pb);
    term2 = (F0^2 * m_vibrator^2) / (m_roller_dynamic^2) - Pb^2;
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
