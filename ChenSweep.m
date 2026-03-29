%TODO: Derive scaling laws for me, f, ... etc params to optimize
clear
%% CONFIGURATION PARAMETERS (EDIT HERE)




% ======================= Environment Params =========================
g_moon   = 1.62;      % m/s^2

% Regolith layer being compacted (representative 1D column)
h_layer  = 0.1;        % m, layer thickness to compact in one lift
rho_i    = 1.3;       % g/cm^3, initial density (1.8 g/cm^3 from Chen)
rho_f    = 1.8;       % g/cm^3, target FSSD ~2.14 g/cm^3
rho = rho_i; %current density init (updates in loop)

nu = 0.35; %LHS poisson ratio, technically changes with density but not by much (and I couldn't find a good closed form hueristic)
%https://ascelibrary.org/doi/10.1061/%28ASCE%29AS.1943-5525.0000848 <-poisson ratio source

% Contact area estimate
%[rp_eff,A_col,~] = rp_fun(R_roller,b_roller,h_layer); %replaces the below snippet
% theta = acos((R_roller-h_layer)/R_roller); %half angle
% lc = 2*theta*R_roller; %full arc length of contact patch
% A_contact = b_roller * lc;               % m^2, contact area
% % pi*r_eff^2 = b*lc => r_eff = sqrt(a_contact/pi)
% rp_eff = sqrt(A_contact/pi); %equivilant circular radius interfacing with chen et al

%Soil moduli (fx of density) -- ALL RHOS HERE IN G/CM3 -- all these
E = @(rho) 6.498e-10 .* 1e6 .* exp(12.07.*rho); %Comp modulus (Pa) Chen et al eq4
k_cr =  @(rho) 5.686e12 * rho^(-41.58) + 0.9079; %deformation ratio Chen et al eq6
E_su = @(rho) k_cr(rho)*E(rho); %Resilient modulus Chen et al eq5
ks_fun = @(rho,rp) 2*rp * E(rho) / (1-nu);
ksu_fun = @(rho,rp) 2*rp * E_su(rho) / (1-nu);

%============== Mission / Rover Params =======================
% Pad mission requirements
Apad     = 500;       % m^2
Tavail_h = 365*24*5;    % hours, total allowed compaction time (2 years but with 6 month MOE)
Tavail   = Tavail_h * 3600;  % seconds

% Rover constraints
m_rov       = 400;                  % kg
Roller_frac = 0.5;                 % fraction of rover normal force on roller
mt   = Roller_frac * m_rov; % kg effective vibrating mass. Mass in DOF (roller + attached structure that "bounces")

% Roller params (initial guess)
mp = 150; %kg -- mass of roller directly vibrates
b_roller = 0.3;        % m, roller width
D_roller = 0.4;
R_roller = D_roller/2;       % m, roller diameter
%Vibrator
f = 50;
me = 0.25; %kg eccentric mass
re = 0.10;% m eccentric radius


% --- 1. DEFINE SWEEP VECTORS ---
% Swap these out easily: set a parameter to an array (linspace/logspace) to sweep it,
% or leave it as a scalar to keep it constant. (Array -> Swept)
f_sweep  = f;      % Hz frequency of eccentric mass
mp_sweep = linspace(2,7,5); %kg -- mass of roller directly vibrates
me_sweep = linspace(0.001,0.5,30); % kg eccentric mass' mass (Scalar -> Not swept)
re_sweep = re; %m eccentric mass' radius (Scalar -> Not swept)
mt_sweep = linspace(25,60,20); % kg effective vibrating mass. Mass in DOF (roller + attached structure that "bounces")

% Simulation settings
N_cycles    = 1e5;     % max vibration cycles to simulate
t_cycle     = 1/f;     % s, period of one vibration cycle

%% Feasability check (this should go in the loop) and discard to NaN
%Will the force make the rover bounce? (not good)
%todo
% % EXCITATION FORCE
% F0    = me * re * omega^2; % N, vibratory force amplitude
% F_peak = W_roller + F0;    % N, static + dynamic peak load (used only in logic, optional)

%Time Check
%T to complete A_pad < T_Avail

%% MAIN ELASTO-PLASTIC COMPACTION LOOP
% Refactored for multi-dimensional parameter sweeping.

% --- 1. SANITY CHECK: SIMULATION TIME BOUNDS ---
% Estimate initial contact area to scale single-spot time to full pad
[~, A_contact_init, ~] = rp_fun(R_roller, b_roller, h_layer);

% Calculate max time bounds (worst-case scenario if it hits N_cycles at lowest frequency)
max_spot_time_s = N_cycles * (1 / min(f_sweep));
max_pad_time_hr = (max_spot_time_s * (Apad / A_contact_init)) / 3600;

fprintf('--- Sanity Check ---\n');
fprintf('Max simulated real-time per spot: %.2f seconds.\n', max_spot_time_s);
fprintf('Theoretical max time to complete %dm^2 pad: %.1f hours.\n', Apad, max_pad_time_hr);

% Only throw a warning if it actually exceeds the mission constraints
if max_pad_time_hr > Tavail_h
    fprintf('WARNING: The maximum simulated time could exceed the allowed mission time (%.1f hrs).\n', Tavail_h);
    fprintf('Consider reducing N_cycles or increasing the minimum frequency to prevent stalling the sweep.\n\n');
else
    fprintf('Check passed: Max simulation bounds are well within mission time (%.1f hrs).\n\n', Tavail_h);
end

% --- 2. BUILD N-DIMENSIONAL GRID ---
[F_grid, MP_grid, ME_grid, RE_grid, MT_grid] = ndgrid(f_sweep, mp_sweep, me_sweep, re_sweep, mt_sweep);
T_results   = zeros(size(F_grid));
Rho_results = zeros(size(F_grid));

% Dictionary for labeling and extracting active dimensions
sweep_vars = {
    'Frequency f (Hz)', f_sweep;
    'Vibrating Mass mp (kg)', mp_sweep;
    'Eccentric Mass me (kg)', me_sweep;
    'Eccentric Radius re (m)', re_sweep;
    'Total Mass mt (kg)', mt_sweep
};

% Identify which parameters are being swept (length > 1)
active_idx = [];
for i = 1:size(sweep_vars,1)
    if length(sweep_vars{i,2}) > 1
        active_idx(end+1) = i;
    end
end

% --- 3. EXECUTE SWEEP ---
fprintf('Running simulation for %d combinations...\n', numel(F_grid));

% For the bounce check
W_rover = m_rov * g_moon; 

% Diagnostics counters
fail_bounce = 0; fail_yield = 0; fail_time = 0; success = 0;

for i = 1:numel(F_grid)
    omega_i = 2 * pi * F_grid(i);
    F0_i    = ME_grid(i) * RE_grid(i) * omega_i^2;
    
    % A. Bounce Check
    if F0_i >= W_rover
        T_results(i) = NaN; Rho_results(i) = NaN;
        fail_bounce = fail_bounce + 1;
        continue; 
    end

    % B. Run Simulation
    [rho_final, t_point, n_cycles, lc_avg_dynamic] = run_compaction_sim(...
        F_grid(i), ME_grid(i), RE_grid(i), MP_grid(i), MT_grid(i), ...
        h_layer, rho_i, rho_f, R_roller, b_roller, nu, N_cycles);
    
    Rho_results(i) = rho_final;
    
    % C. Time & Constraints Check
    if rho_final >= rho_f
        A_contact_avg = b_roller * lc_avg_dynamic;
        t_pad_hr = (t_point * (Apad / A_contact_avg)) / 3600;
        
        if t_pad_hr <= Tavail_h
            T_results(i) = t_pad_hr;
            success = success + 1;
        else
            T_results(i) = NaN; 
            fail_time = fail_time + 1;
        end
    else
        T_results(i) = NaN; 
        fail_yield = fail_yield + 1; % Failed to reach density (likely xcr=0)
    end
end

fprintf('Sweep complete. Diagnostics:\n');
fprintf('  Successes: %d\n', success);
fprintf('  Failed (Rover Bounced): %d\n', fail_bounce);
fprintf('  Failed (No Soil Yielding): %d\n', fail_yield);
fprintf('  Failed (Exceeded Time): %d\n\n', fail_time);

% --- 4. PLOTTING (Robust against sparse successes) ---
grids = {F_grid, MP_grid, ME_grid, RE_grid, MT_grid};

if length(active_idx) == 2
    % 2D Sweep Data Extraction
    X = squeeze(grids{active_idx(1)});
    Y = squeeze(grids{active_idx(2)});
    Z = squeeze(T_results);
    R = squeeze(Rho_results);
    
    nameX = sweep_vars{active_idx(1), 1};
    nameY = sweep_vars{active_idx(2), 1};
    
    % Find where we actually have successful data (Not NaN)
    valid_idx = ~isnan(Z);
    
    % Plot 1: Diagnostic Physics Boundary Map
    % This shows exactly WHY points failed
    figure('Name', 'Physics Boundaries');
    pcolor(X, Y, R); shading flat;
    colormap parula; h = colorbar; ylabel(h, 'Final Density (g/cm^3)');
    xlabel(nameX); ylabel(nameY); 
    title('Compaction Regimes (Blank/White = Rover Bounced)');
    
    % Plot 2: Scatter Plot of Valid Times
    % scatter() ignores NaNs, unlike surf() and contourf()
    figure('Name', 'Time to Complete Pad');
    if any(valid_idx, 'all')
        scatter(X(valid_idx), Y(valid_idx), 80, Z(valid_idx), 'filled', 'MarkerEdgeColor', 'k');
        colormap jet; h2 = colorbar; ylabel(h2, 'Total Pad Time (hr)');
        xlabel(nameX); ylabel(nameY); 
        title('Successful Compaction Times (Hours)');
        grid on;
    else
        % Fallback if 0 successes exist
        text(0.5, 0.5, '0 Successes Found in this range', 'HorizontalAlignment', 'center', 'FontSize', 14);
        axis off;
    end
    
elseif length(active_idx) == 3
    % 3D Sweep Data Extraction
    X = squeeze(grids{active_idx(1)});
    Y = squeeze(grids{active_idx(2)});
    Z = squeeze(grids{active_idx(3)});
    V = squeeze(T_results);
    
    nameX = sweep_vars{active_idx(1), 1};
    nameY = sweep_vars{active_idx(2), 1};
    nameZ = sweep_vars{active_idx(3), 1};
    
    valid_idx = ~isnan(V);
    
    figure('Name', '4D Scatter Heatmap');
    if any(valid_idx, 'all')
        scatter3(X(valid_idx), Y(valid_idx), Z(valid_idx), 80, V(valid_idx), 'filled', 'MarkerEdgeColor', 'k');
        colormap jet; h3 = colorbar; ylabel(h3, 'Total Pad Time (hr)');
        xlabel(nameX); ylabel(nameY); zlabel(nameZ);
        title('Successful Pad Compaction Times (hr)');
        grid on;
    else
        text(0.5, 0.5, '0 Successes Found', 'HorizontalAlignment', 'center', 'FontSize', 14);
        axis off;
    end
else
    fprintf('Set exactly 2 or 3 sweep arrays to a length > 1 to view automatic plots.\n');
end

%% Functions
function [rho, t_point, n, lc_avg_dynamic] = run_compaction_sim(f, me, re, mp, mt, h_layer_init, rho_i, rho_f, R_roller, b_roller, nu, N_cycles)
    % Recalculate dynamic properties for this specific iteration
    omega = 2*pi*f;
    F0    = me * re * omega^2;
    t_cycle = 1/f;
    
    rho = rho_i;
    h_layer = h_layer_init;
    lc_history = zeros(1,N_cycles);
    
    for n = 1:N_cycles
        [rp_eff, A_col, lc] = rp_fun(R_roller, b_roller, h_layer);
        lc_history(n) = lc;
        
        % --- SIMULANT TRANSLATION (LHS-1 to Chen) ---
        % LHS-1 bounds: ~1.27 (fluff) to 1.95 (max compaction)
        % Chen bounds:  1.63 (fluff) to 2.15 (max compaction)
        
        % Calculate linear scaling factor between the two simulants using
        % Relative density
        rho_eq = 1.63 + (rho - 1.27) * ((2.15 - 1.63) / (1.95 - 1.27));
        
        % Clamp the equivalent density to prevent exceeding Chen's tested bounds
        rho_eq = min(max(rho_eq, 1.63), 2.15);

        % --- CALCULATE MODULI USING EQUIVALENT DENSITY ---
        Pb0    = 0.07932 * (100*rho_eq - 184)^2;
        A_chen = pi * 0.151^2;
        Pb     = Pb0 * (A_col/A_chen);

        E_val = 6.498e-10 * 1e6 * exp(12.07*rho_eq);
        k_cr_val = 5.686e12 * rho_eq^(-41.58) + 0.9079;
        E_su_val = k_cr_val * E_val;
        
        ks  = 2*rp_eff * E_val / (1-nu);
        ksu = 2*rp_eff * E_su_val / (1-nu);

        xcr = chen_residual_deformation(F0, omega, ks, ksu, Pb, mp, mt);

        if xcr > 0
            % Calculate the exact minimum thickness needed to hit rho_f
            h_target = h_layer_init * (rho_i / rho_f);
            
            h_layer_new = h_layer - xcr;
            
            % CLAMP: Prevent single-cycle overshoots past target density
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
            % Performance Optimization: If no plastic deformation occurs, 
            % it never will at this state. Break to prevent infinite stall.
            break; 
        end
    end
    t_point = n * t_cycle;
    lc_avg_dynamic = mean(lc_history(1:n));
end

function xcr = chen_residual_deformation(F0, omega, ks, ksu, Pb, mp, mt)
    % Implement Eq. (15) from Chen et al. using ks and ksu. [file:1]
    %
    % Inputs:
    %   F0   - excitation force amplitude (N)
    %   omega- angular frequency (rad/s)
    %   ks   - compressive stiffness (N/m)
    %   ksu  - resilient stiffness (N/m)
    %   Pb   - plastic limit (N)
    %   mp   - "mass of mass" (kg)
    %   mt   - table/structure mass (kg)
    %
    % Output:
    %   xcr  - residual plastic deformation per cycle (m)

    % If dynamic force cannot exceed Pb, no plasticity
    if F0 <= Pb
        xcr = 0;
        return;
    end

    % Eq.(15) structure:
    % xcr = (mp * omega^2)/(2 * ks^2 * Pb) * (F0^2 * mp^2 / mt^2 - Pb^2) ...
    %        + Pb/ks - Pb/ksu;
    term1 = (mp * omega^2) / (2 * ks^2 * Pb);
    term2 = (F0^2 * mp^2) / (mt^2) - Pb^2;
    term3 = Pb/ks - Pb/ksu;

    xcr = term1 * term2 + term3;

    % Ensure non-negative residual
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
    % pi*r_eff^2 = b*lc => r_eff = sqrt(a_contact/pi) (does flat vs curved surface matter?)
    rp_eff = sqrt(A_contact/pi); %equivilant circular radius interfacing with chen et al
end
