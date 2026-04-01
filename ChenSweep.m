clear; clc
%% CONFIGURATION PARAMETERS (EDIT HERE)
% Pareto Optimization Weights [Total Mass, Power, Mass Ratio, Freq Diff from 50Hz]
pareto_weights = [1.0, 0.5, 1.0, 0.2]; 

% Simulation safety limits
xcr_cap_frac = 0.15; % reject if any single xcr > 15% of h_layer

% System Efficiency
eta_sys = 0.05; % Assumed 5% efficiency from electrical input to plastic soil deformation

% ======================= Environment Params =========================
g_moon   = 1.62;      % m/s^2

% Regolith layer being compacted (representative 1D column)
h_layer  = 0.1;        % m, layer thickness to compact in one lift
rho_i    = 1.3;       % g/cm^3, initial density (1.8 g/cm^3 from Chen)
rho_f    = 1.8;       % g/cm^3, target FSSD ~2.14 g/cm^3

nu = 0.35; %LHS poisson ratio

% Pad mission requirements
Apad     = 500;       % m^2
Tavail_h = 365*24*5;  % hours, total allowed compaction time (2 years but with 6 month MOE)

% --- 1. DEFINE SWEEP VECTORS ---
f_sweep  = linspace(5, 65, 10);      % Hz frequency of eccentric mass
m0e_sweep = logspace(-3, -1, 10);    % kg-m, eccentric moment
m_rov_sweep = linspace(20, 100, 10); % kg, total rover mass
m_roller_frac_sweep = linspace(0.05, 0.5, 10); % fraction of rover mass dedicated to roller
A_vibrator_sweep = logspace(-2, log10(0.5), 10); % m^2, vibrator contact area

% Simulation settings
N_cycles    = 1e4;     % max vibration cycles to simulate

%% MAIN ELASTO-PLASTIC COMPACTION LOOP

% --- 1. SANITY CHECK: SIMULATION TIME BOUNDS ---
A_contact_init = min(A_vibrator_sweep);
max_spot_time_s = N_cycles * (1 / min(f_sweep));
max_pad_time_hr = (max_spot_time_s * (Apad / A_contact_init)) / 3600;

fprintf('--- Sanity Check ---\n');
fprintf('Max simulated real-time per spot: %.2f seconds.\n', max_spot_time_s);
fprintf('Theoretical max time to complete %dm^2 pad: %.1f hours.\n', Apad, max_pad_time_hr);

if max_pad_time_hr > Tavail_h
    fprintf('WARNING: The maximum simulated time could exceed the allowed mission time.\n\n');
else
    fprintf('Check passed: Max simulation bounds are well within mission time.\n\n');
end

% --- 2. BUILD N-DIMENSIONAL GRID ---
[F_grid, M0E_grid, M_ROV_grid, M_ROLLER_FRAC_grid, A_VIB_grid] = ndgrid(f_sweep, m0e_sweep, m_rov_sweep, m_roller_frac_sweep, A_vibrator_sweep);
T_results   = nan(size(F_grid));
Rho_results = nan(size(F_grid));
Ncyc_results = nan(size(F_grid));
Power_results = nan(size(F_grid));
MassRatio_results = nan(size(F_grid));
Proxy_results = nan(size(F_grid));

% Dictionary for labeling and extracting active dimensions
sweep_vars = {
    'Frequency f (Hz)', f_sweep;
    'Eccentric Moment m0e (kg-m)', m0e_sweep;
    'Total Rover Mass m_rov (kg)', m_rov_sweep;
    'Roller Mass Fraction', m_roller_frac_sweep;
    'Vibrator Area (m^2)', A_vibrator_sweep
};

% --- 3. EXECUTE SWEEP ---
fprintf('Running simulation for %d combinations...\n', numel(F_grid));
fail_bounce = 0; fail_yield = 0; fail_time = 0; fail_xcrcap = 0; success = 0;

for i = 1:numel(F_grid)
    % Calculate masses for this iteration
    m_rov_i = M_ROV_grid(i);
    m_roller_frac_i = M_ROLLER_FRAC_grid(i);
    m_roller_i = m_rov_i * m_roller_frac_i;
    m_sprung_i = m_rov_i * (1 - m_roller_frac_i);
    
    % Calculate dynamic forces
    omega_i = 2 * pi * F_grid(i);
    F0_i    = M0E_grid(i) * omega_i^2;
    
    % Bounce Check (compare against total rover weight)
    W_rover_i = m_rov_i * g_moon;
    if F0_i >= W_rover_i
        fail_bounce = fail_bounce + 1;
        continue; 
    end
    
    % Run Simulation
    [rho_final, t_point, n_cycles, P_ss_avg, xcr_max] = run_compaction_sim(...
        F_grid(i), M0E_grid(i), m_roller_i, m_sprung_i, ...
        h_layer, rho_i, rho_f, A_VIB_grid(i), nu, N_cycles, eta_sys);
    
    Rho_results(i) = rho_final;
    
    % XCR Cap Filter
    if xcr_max > xcr_cap_frac * h_layer
        fail_xcrcap = fail_xcrcap + 1;
        continue;
    end
    
    % Time & Constraints Check
    if rho_final >= rho_f
        t_pad_hr = (t_point * (Apad / A_VIB_grid(i))) / 3600;
        
        if t_pad_hr <= Tavail_h
            T_results(i) = t_pad_hr;
            MassRatio_results(i) = m_roller_i / m_sprung_i;
            Power_results(i) = P_ss_avg;
            Ncyc_results(i) = n_cycles;
            Proxy_results(i) = (F0_i^2) * omega_i;
            success = success + 1;
        else
            fail_time = fail_time + 1;
        end
    else
        fail_yield = fail_yield + 1; 
    end
end

fprintf('Sweep complete. Diagnostics:\n');
fprintf('  Successes: %d\n', success);
fprintf('  Failed (Rover Bounced): %d\n', fail_bounce);
fprintf('  Failed (No Soil Yielding): %d\n', fail_yield);
fprintf('  Failed (Exceeded Time): %d\n', fail_time);
fprintf('  Failed (XCR Cap): %d\n\n', fail_xcrcap);

% --- 4. PARETO OPTIMIZATION ---
valid_idx = find(~isnan(T_results));
if isempty(valid_idx)
    fprintf('No successful points found. Cannot perform Pareto analysis.\n');
else
    % Extract valid vectors
    m_rov_vec = M_ROV_grid(valid_idx);
    power_vec = Power_results(valid_idx);
    mass_ratio_vec = MassRatio_results(valid_idx);
    f_vec = F_grid(valid_idx);
    Ncyc_vec = Ncyc_results(valid_idx); % Kept for tracking/plotting
    
    % 4-Objective Matrix: [Minimize Mass, Minimize Power, Maximize Ratio, Target 50Hz]
    objectives = [m_rov_vec, power_vec, 1./mass_ratio_vec, abs(f_vec - 50)];
    
    % Find the non-dominated set
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
    
    % Find "knee" point using weighted normalization
    min_vals = min(pareto_objectives, [], 1);
    max_vals = max(pareto_objectives, [], 1);
    
    range_vals = max_vals - min_vals;
    range_vals(range_vals == 0) = 1; 
    
    norm_objectives = (pareto_objectives - min_vals) ./ range_vals;
    weighted_scores = norm_objectives * pareto_weights';
    [~, knee_local_idx] = min(weighted_scores);
    knee_global_idx = pareto_idx(knee_local_idx);
    
    % Recommended Design Output
    rec_m_rov = M_ROV_grid(knee_global_idx);
    rec_m_roller_frac = M_ROLLER_FRAC_grid(knee_global_idx);
    fprintf('--- Recommended Optimal Design ---\n');
    fprintf('  Total Rover Mass: %.2f kg\n', rec_m_rov);
    fprintf('  Roller Mass:  %.2f kg (%.1f %%)\n', rec_m_rov * rec_m_roller_frac, rec_m_roller_frac*100);
    fprintf('  Sprung Mass:  %.2f kg\n', rec_m_rov * (1-rec_m_roller_frac));
    fprintf('  Frequency:    %.1f Hz\n', F_grid(knee_global_idx));
    fprintf('  Eccentric Moment: %.4f kg-m\n', M0E_grid(knee_global_idx));
    fprintf('  Vibrator Area: %.3f m^2\n', A_VIB_grid(knee_global_idx));
    fprintf('  Cycles to Compact: %d\n', Ncyc_results(knee_global_idx));
    fprintf('  Total Pad Time:    %.1f hours\n', T_results(knee_global_idx));
    fprintf('  Avg. SS Power:     %.2f Watts\n\n', Power_results(knee_global_idx));

    % --- 5. PLOTTING ---
    
    % Figure 1: 2D Pareto Front Scatter
    figure('Name', 'Pareto Front: Mass vs. Power');
    pareto_m_rov = M_ROV_grid(pareto_idx);
    pareto_power = Power_results(pareto_idx);
    pareto_mass_ratios = MassRatio_results(pareto_idx);
    
    scatter(pareto_m_rov, pareto_power, 80, pareto_mass_ratios, 'filled');
    hold on;
    
    % Mark the knee point
    knee_m_rov = M_ROV_grid(knee_global_idx);
    knee_p = Power_results(knee_global_idx);
    plot(knee_m_rov, knee_p, 'p', 'MarkerSize', 20, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'k');
    hold off;
    
    xlabel('Total Rover Mass (kg)');
    ylabel('Steady-State Power (W)');
    set(gca, 'YScale', 'log');
    grid on;
    title('Pareto Front: Total Mass vs. Power');
    colormap(flipud(parula));
    h = colorbar;
    ylabel(h, 'Mass Ratio (m_{roller} / m_{sprung})');

    % Figure 2: Feasibility Map
    figure('Name', 'Feasibility Map');
    max_rho_map = squeeze(max(Rho_results, [], [2 4 5]));
    [X, Y] = ndgrid(f_sweep, m_rov_sweep);
    pcolor(X, Y, max_rho_map'); shading flat;
    colormap(gca, 'jet');
    h_feas = colorbar;
    ylabel(h_feas, 'Max Final Density (g/cm^3)');
    hold on;
    plot(F_grid(pareto_idx), M_ROV_grid(pareto_idx), '^w', 'MarkerSize', 8, 'MarkerFaceColor', 'white');
    plot(F_grid(knee_global_idx), M_ROV_grid(knee_global_idx), '*r', 'MarkerSize', 12, 'LineWidth', 1.5);
    hold off;
    xlabel('Frequency (Hz)');
    ylabel('Total Rover Mass (kg)');
    title('Feasibility Map (\Delta=Pareto, \star=recommended)');

    % Figure 3: Pareto Table
    figure('Name', 'Pareto Optimal Set');
    tbl_data_points = [
        M_ROV_grid(pareto_idx), ...
        M_ROV_grid(pareto_idx) .* M_ROLLER_FRAC_grid(pareto_idx), ...
        M_ROV_grid(pareto_idx) .* (1-M_ROLLER_FRAC_grid(pareto_idx)), ...
        MassRatio_results(pareto_idx), ...
        A_VIB_grid(pareto_idx), ...
        F_grid(pareto_idx), ...
        M0E_grid(pareto_idx), ...
        Ncyc_results(pareto_idx), ...
        T_results(pareto_idx), ...
        Power_results(pareto_idx), ...
        Proxy_results(valid_idx(is_pareto))
    ];
    
    [~, sort_order] = sort(tbl_data_points(:,1), 'ascend'); % Sort by Total Mass
    tbl_data = tbl_data_points(sort_order, :);
    
    cnames = {'m_rov (kg)', 'm_roller (kg)', 'm_sprung (kg)', 'Mass Ratio', 'Area (m^2)', 'f (Hz)', 'm0e (kg-m)', 'N_cycles', 'T_pad (hr)', 'Power (W)', 'Mech Proxy'};
    t = uitable('Data', tbl_data, 'ColumnName', cnames, 'RowName',[], 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
    
    knee_sorted_idx = find(sort_order == knee_local_idx);
    s = uistyle('BackgroundColor', 'yellow');
    addStyle(t, s, 'row', knee_sorted_idx);

    % --- Figure 4: Parallel Coordinates Plot ---
    figure('Name', 'Parallel Coordinates: Trade Space', 'Position', [100, 100, 1000, 500]);
    
    pareto_f = F_grid(pareto_idx);
    pareto_m0e = M0E_grid(pareto_idx);
    pareto_m_rov = M_ROV_grid(pareto_idx);
    pareto_frac = M_ROLLER_FRAC_grid(pareto_idx);
    pareto_A = A_VIB_grid(pareto_idx); 
    pareto_mr = MassRatio_results(pareto_idx);
    pareto_cyc = Ncyc_results(pareto_idx);
    
    % Log scale the power
    pareto_log_pwr = log10(Power_results(pareto_idx)); 

    pareto_table = table(pareto_f, pareto_m0e, pareto_m_rov, pareto_frac, pareto_A, pareto_mr, pareto_log_pwr, pareto_cyc, ...
        'VariableNames', {'Freq_Hz', 'Ecc_Moment', 'Rover_Mass', 'Roller_Frac', 'Area_m2', 'Mass_Ratio', 'Log10_Power_W', 'Cycles'});
    
    p = parallelplot(pareto_table);
    try
        p.ColorVariable = 'Log10_Power_W';
        p.LineAlpha = 0.4; 
    catch
    end
    title('Parallel Coordinates: Pareto Optimal Configurations');
    p.CoordinateVariables = {'Rover_Mass', 'Roller_Frac', 'Area_m2', 'Freq_Hz', 'Ecc_Moment', 'Cycles', 'Log10_Power_W', 'Mass_Ratio'};
end

%% Functions
function [rho, t_point, n, P_ss_avg, xcr_max] = run_compaction_sim(f, m0e, m_roller, m_sprung, h_layer_init, rho_i, rho_f, A_vibrator, nu, N_cycles, eta_sys)
    omega = 2*pi*f;
    F0    = m0e * omega^2;
    t_cycle = 1/f;
    
    rho = rho_i;
    h_layer = h_layer_init;
    P_ss_history = zeros(1,N_cycles);
    xcr_max = 0;
    
    A_col = A_vibrator;
    rp_eff = sqrt(A_vibrator / pi);
    
    for n = 1:N_cycles
        % --- SIMULANT TRANSLATION (LHS-1 to Chen) ---
        rho_eq = 1.63 + (rho - 1.27) * ((2.15 - 1.63) / (1.95 - 1.27));
        rho_eq = min(max(rho_eq, 1.63), 2.15);
        
        % --- CALCULATE MODULI ---
        Pb0    = 0.07932 * (100*rho_eq - 184)^2;
        A_chen = pi * 0.151^2;
        Pb     = Pb0 * (A_col/A_chen);
        E_val = 6.498e-10 * 1e6 * exp(12.07*rho_eq);
        k_cr_val = 5.686e12 * rho_eq^(-41.58) + 0.9079;
        E_su_val = k_cr_val * E_val;
        
        ks  = 2*rp_eff * E_val / (1-nu);
        ksu = 2*rp_eff * E_su_val / (1-nu);
        
        xcr = chen_residual_deformation(F0, omega, ks, ksu, Pb, m_roller, m_sprung);
        xcr_max = max(xcr_max, xcr);
        
        E_plastic_cycle = Pb * xcr; 
        P_ss_cycle = (E_plastic_cycle * f) / eta_sys; 
        P_ss_history(n) = P_ss_cycle;
        
        if xcr > 0
            h_target = h_layer_init * (rho_i / rho_f);
            h_layer_new = h_layer - xcr;
            
            if h_layer_new <= h_target
                h_layer_new = h_target;
            end
            
            V_old = A_col * h_layer;
            V_new = A_col * h_layer_new;
            
            rho = rho * (V_old / V_new);
            h_layer = h_layer_new;
            
            if rho >= rho_f
                break;
            end
        else
            break; 
        end
    end
    t_point = n * t_cycle;
    P_ss_avg = mean(P_ss_history(1:n));
end

function xcr = chen_residual_deformation(F0, omega, ks, ksu, Pb, m_roller, m_sprung)
    if F0 <= Pb
        xcr = 0;
        return;
    end
    term1 = (m_roller * omega^2) / (2 * ks^2 * Pb);
    term2 = (F0^2 * m_roller^2) / (m_sprung^2) - Pb^2;
    term3 = Pb/ks - Pb/ksu;
    xcr = term1 * term2 + term3;
    xcr = max(0, xcr);
end

function [rp_eff,A_contact,lc] = rp_fun(r,b,h)
    %Input:
    % r = roller radius
    % b = roller width
    % h = depth of compactor pass (how far into the soil are you pressng the thing?
    %Output:
    %rp_eff = effective radius ofa a circular puck with same area (chen did it this way but we're rolling so I needed to find equivilants.)
    %A_contact = contact patch
    %lc = arc length for contact area
    
    theta = acos((r-h)/r); %half angle
    lc = 2*theta*r; %full arc length of contact patch
    A_contact = b * lc;               % m^2, contact area
    % pi*r_eff^2 = b*lc => r_eff = sqrt(a_contact/pi)
    rp_eff = sqrt(A_contact/pi); %equivilant circular radius interfacing with chen et al
end