%TODO: Derive scaling laws for me, f, ... etc params to optimize

%% CONFIGURATION PARAMETERS (EDIT HERE)
% Environment
g_moon   = 1.62;      % m/s^2

% Pad mission requirements
Apad     = 500;       % m^2
Tavail_h = 365*24;    % hours, total allowed compaction time
Tavail   = Tavail_h * 3600;  % seconds

% Rover constraints
m_rov       = 80;                  % kg
n_rollers   = 1;                   % total number of rollers on the rover
target_relative_density = 0.85;    % desired final compaction state (0-1)
roller_fraction = 0.5;             % fraction of rover weight on rollers
bounce_margin = 1.0;               % Multiplier for the bounce safety check
W_rov       = m_rov * g_moon;      % N, lunar weight
W_roller    = (m_rov * g_moon * roller_fraction) / n_rollers; % N, static load on single roller

% Roller params (initial guess)
b_roller = 0.3;        % m, roller width
D_roller = 0.4;

% Mass & Vibration
m_eccentric = 0.001;    % kg*m, renamed from m0e
mass_ratio  = 0.5;      % logic for participating mass fraction
m_vibrator  = 5;        % kg, mass of roller directly vibrates
eta_mech    = 0.6;      % Mechanical drivetrain efficiency

% Regolith targets (Lunar)
rho_min_lunar = 1.27; % g/cm^3
rho_max_lunar = 1.95; % g/cm^3
rho_i         = 1.3;  % g/cm^3, initial density
nu            = 0.35; % Poisson ratio

SOIL.rho   = 1600;   % kg/m^3
SOIL.gamma = 2470;   % N/m^3
SOIL.n     = 1;
SOIL.kc    = 1400;   % N/m^2
SOIL.kphi  = 830000; % N/m^3
SOIL.phi   = 0.576;  % rad
SOIL.c     = 170;    % N/m^2

h_layer  = 0.1;        % m, layer thickness to compact in one lift

% Vibration / eccentric mass
f = 30; % Hz

% Simulation settings
v_sim = 0.05; % m/s, speed for pass-based tracking

%% RUN SIMULATION WRAPPER
[rho_final, req_power, total_passes, total_cycles, lc_avg_dynamic, is_valid] = ...
    run_compaction_sim(W_roller, b_roller, D_roller, m_eccentric, f, v_sim, mass_ratio, m_rov, target_relative_density, bounce_margin, h_layer, SOIL, eta_mech, m_vibrator, nu, rho_min_lunar, rho_max_lunar, rho_i);

if is_valid
    fprintf('Target reached! Final density: %.2f g/cm^3, Peak Power: %.1f W, Total passes: %d\n', rho_final, req_power, total_passes);
else
    fprintf('Simulation terminated due to invalidity or stall.\n');
end

%% PAD-LEVEL PRODUCTIVITY TABLE (CLEAN FORMAT)
if is_valid
    v_rover_range = logspace(-3, log10(0.5), 20);
    track_length = sqrt(Apad);

    fprintf('\n');
    fprintf('═%s═%s═%s═%s═%s═\n', repmat('═',1,12), repmat('═',1,10), ...
                               repmat('═',1,8), repmat('═',1,12), repmat('═',1,10));
    fprintf('| %8s | %8s | %7s | %9s | %9s |\n', 'Speed(m/s)', 'Cycles/Pass', ...
           'Passes', 'Time/Pass(min)', 'Total(hr)');
    fprintf('─%s─%s─%s─%s─%s─\n', repmat('─',1,12), repmat('─',1,10), ...
                               repmat('─',1,8), repmat('─',1,12), repmat('─',1,10));

    for i = 1:length(v_rover_range)
        v = v_rover_range(i);
        
        cycles_pass = f * (lc_avg_dynamic / v);
        N_passes = ceil(total_cycles / cycles_pass);
        t_pass = track_length / v / 60;     % min
        t_total_hr = N_passes * t_pass / 60; % hr
        
        fprintf('| %8.4f | %8.0f | %7.0f | %9.1f | %9.1f |\n', ...
                v, cycles_pass, N_passes, t_pass, t_total_hr);
    end

    fprintf('─%s─%s─%s─%s─%s─\n', repmat('─',1,12), repmat('─',1,10), ...
                               repmat('─',1,8), repmat('─',1,12), repmat('─',1,10));
    fprintf('| lc_avg: %.4f m | Total cycles: %d | Budget: %.0f hr |\n', ...
            lc_avg_dynamic, total_cycles, Tavail/3600);
    fprintf('═%s═%s═%s═%s═%s═\n', repmat('═',1,12), repmat('═',1,10), ...
                               repmat('═',1,8), repmat('═',1,12), repmat('═',1,10));
end

%% Functions
function [rho_final, req_power, total_passes, total_cycles, lc_avg_dynamic, is_valid] = ...
    run_compaction_sim(W_roller, b_roller, D_roller, m_eccentric, f, v_sim, mass_ratio, m_rov, target_relative_density, bounce_margin, h_layer, SOIL, eta_mech, m_vibrator, nu, rho_min_lunar, rho_max_lunar, rho_i)
    
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

    while rho < rho_f
        current_pass = current_pass + 1;
        
        % 1. Calculate TOTAL sinkage
        N_total_current = current_pass * (W_roller > 0); % Assume n_rollers effect captured in W_roller scale
        z_current = calculate_tandem_sinkage(W_roller, b_roller, D_roller, SOIL.kc, SOIL.kphi, N_total_current);
        h_pass = z_current - z_prev_pass;
        
        % 2. Recalculate rp_eff based on incremental sinkage
        [rp_eff, A_col, lc] = rp_fun(R_roller, b_roller, h_pass);
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
    else
        rho_final = NaN;
        req_power = NaN;
        total_passes = NaN;
        total_cycles = NaN;
        lc_avg_dynamic = NaN;
    end
end

function xcr = chen_residual_deformation(F0, omega, ks, ksu, Pb, m_vibrator, m_roller_dynamic)
    % Implement Eq. (15) from Chen et al. using ks and ksu.
    %
    % Inputs:
    %   F0   - excitation force amplitude (N)
    %   omega- angular frequency (rad/s)
    %   ks   - compressive stiffness (N/m)
    %   ksu  - resilient stiffness (N/m)
    %   Pb   - plastic limit (N)
    %   m_vibrator - renamed from mp (kg)
    %   m_roller_dynamic - renamed from mt (kg)
    %
    % Output:
    %   xcr  - residual plastic deformation per cycle (m)

    % If dynamic force cannot exceed Pb, no plasticity
    if F0 <= Pb
        xcr = 0;
        return;
    end

    % Eq.(15) structure:
    term1 = (m_vibrator * omega^2) / (2 * ks^2 * Pb);
    term2 = (F0^2 * m_vibrator^2) / (m_roller_dynamic^2) - Pb^2;
    term3 = Pb/ks - Pb/ksu;

    xcr = term1 * term2 + term3;

    % Ensure non-negative residual
    xcr = max(0, xcr);
end

function [rp_eff,A_contact,lc] = rp_fun(r,b,z)
    %Input:
    % r = roller radius
    % b = roller width
    % z = sinkage (depth) into soil
    %Output:
    %rp_eff = effective radius of a circular puck with same area
    %A_contact = contact patch
    %lc = arc length for contact area

    % Clamp z to radius to avoid complex numbers in acos
    z = min(z, r * 0.999);
    theta = acos((r-z)/r); %half angle
    lc = 2*theta*r; %full arc length of contact patch
    A_contact = b * lc;               % m^2, contact area
    rp_eff = sqrt(A_contact/pi); 
end

function z = calculate_tandem_sinkage(W_roller, b, D, kc, kphi, N_total_current)
    % Recursive sinkage equation for lunar regolith (n=1)
    % N_total_current: cumulative number of rollers that have passed
    z_prev = 0;
    for i = 1:N_total_current
        % Equation: z_i = [3*W_roller / (2*(kc + b*kphi)*sqrt(D)) + z_prev^(3/2)]^(2/3)
        term_load = (3 * W_roller) / (2 * (kc + b * kphi) * sqrt(D));
        z_i = (term_load + z_prev^(1.5))^(2/3);
        z_prev = z_i;
    end
    z = z_prev;
end

function rho_earth = translate_density(rho_lunar)
    % Linear mapping from Lunar simulant bounds to Earth bounds for Chen formulas
    % Map 1.27 -> 1.63
    % Map 1.95 -> 2.15
    slope = (2.15 - 1.63) / (1.95 - 1.27);
    rho_earth = 1.63 + (rho_lunar - 1.27) * slope;
    % Clamp output
    rho_earth = max(1.63, min(2.15, rho_earth));
end
