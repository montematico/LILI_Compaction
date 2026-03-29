%TODO: Derive scaling laws for me, f, ... etc params to optimize

%% CONFIGURATION PARAMETERS (EDIT HERE)
% Environment
g_moon   = 1.62;      % m/s^2

% Pad mission requirements
Apad     = 500;       % m^2
Tavail_h = 365*24;    % hours, total allowed compaction time
Tavail   = Tavail_h * 3600;  % seconds

% Rover constraints
m_rov       = 35;                  % kg
W_rov       = m_rov * g_moon;      % N, lunar weight
Roller_frac = 0.5;                 % fraction of rover normal force on roller
W_roller    = Roller_frac * W_rov; % N

% Roller params (initial guess)
b_roller = 0.3;        % m, roller width
D_roller = 0.4;
R_roller = D_roller/2;       % m, roller diameter

mp = 5; %kg -- mass of roller directly vibrates
mt = W_roller / g_moon;% kg, effective vibrating mass. Mass in DOF (roller + attached structure that "bounces")

% Regolith layer being compacted (representative 1D column)
h_layer  = 0.1;        % m, layer thickness to compact in one lift
rho_i    = 1.3;       % g/cm^3, initial density (1.8 g/cm^3 from Chen)
rho_f    = 1.8;       % g/cm^3, target FSSD ~2.14 g/cm^3
rho = rho_i; %current density init (updates in loop)

nu = 0.35; %LHS poisson ratio, technically changes with density but not by much (and I couldn't find a good closed form hueristic)
%https://ascelibrary.org/doi/10.1061/%28ASCE%29AS.1943-5525.0000848 <-poisson ratio source

% Vibration / eccentric mass (design variables)
f = 30; % Hz
omega = 2*pi*f; % rad/s

me = 0.2; % kg, eccentric mass (guess)
re = 0.05; % m, eccentric radius (guess)

% Simulation settings
N_cycles    = 1e5;     % max vibration cycles to simulate
t_cycle     = 1/f;     % s, period of one vibration cycle

% Contact area estimate
[rp_eff,A_col,~] = rp_fun(R_roller,b_roller,h_layer); %replaces the below snippet
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

% EXCITATION FORCE
F0    = me * re * omega^2; % N, vibratory force amplitude
F_peak = W_roller + F0;    % N, static + dynamic peak load (used only in logic, optional)

%% MAIN ELASTO-PLASTIC COMPACTION LOOP
% Uses Chen-style plastic limit Pb(rho) and residual deformation xcr.

lc_history = zeros(1,N_cycles); %Pre-allocate this?
for n = 1:N_cycles
    %0. Recalculate rp_eff for new h_layer depth.
    [rp_eff,A_col,lc] = rp_fun(R_roller,b_roller,h_layer);
    lc_history(n) = lc;

    %1. Plastic limit for current density (Chen Eq. 18)
    % Eq. (18): Pb = 0.07932 * (100*rho - 184)^2, where rho is density in g/cm^3
    Pb0     = 0.07932 * (100*rho - 184)^2;   % N (model uses this as load threshold) elastic deformation limit
    r_chen = 0.151; %Chen rp (table 3)
    A_chen = pi * r_chen^2;
    Pb = Pb0 * (A_col/A_chen); %Scale by area ratio

    %2. stiffness at current density
    ks = ks_fun(rho,rp_eff);
    ksu = ksu_fun(rho,rp_eff);

    % --- 2. Check if peak load exceeds plastic limit ---
    xcr = chen_residual_deformation(F0, omega, ks, ksu, Pb, mp, mt);

    % Update thickness and density if plastic deformation occurs
    if xcr > 0
        % Assume settlement xcr reduces the layer thickness from h to (h - xcr)
        % while mass stays constant => density increases.
        h_layer_new = h_layer - xcr;
        if h_layer_new <= 0
            h_layer_new = 1e-6; % avoid zero/negative thickness
        end
        
        %Mass conservation in column
        V_old = A_col * h_layer;
        V_new = A_col * h_layer_new;

        rho = rho * (V_old / V_new); %still in g/cm3
        h_layer = h_layer_new;


        % Stop if we reach or exceed target density
        if rho >= rho_f
            fprintf('Target density reached at cycle %d\n', n);
            break;
        end
    end
end

% Total vibration time for this point
t_point = n * t_cycle;   % seconds

fprintf('Final density: %.1f g/cm^3, time at one point: %.1f s\n', rho, t_point);
%% PAD-LEVEL PRODUCTIVITY TABLE (CLEAN FORMAT)
v_rover_range = logspace(-3, log10(0.5), 20);
track_length = sqrt(Apad);

fprintf('\n');
fprintf('═%s═%s═%s═%s═%s═\n', repmat('═',1,12), repmat('═',1,10), ...
                           repmat('═',1,8), repmat('═',1,12), repmat('═',1,10));
fprintf('| %8s | %8s | %7s | %9s | %9s |\n', 'Speed(m/s)', 'Cycles/Pass', ...
       'Passes', 'Time/Pass(min)', 'Total(hr)');
fprintf('─%s─%s─%s─%s─%s─\n', repmat('─',1,12), repmat('─',1,10), ...
                           repmat('─',1,8), repmat('─',1,12), repmat('─',1,10));

lc_avg_dynamic = mean(lc_history);
for i = 1:length(v_rover_range)
    v = v_rover_range(i);
    
    cycles_pass = f * (lc_avg_dynamic / v);
    N_passes = ceil(n / cycles_pass);
    t_pass = track_length / v / 60;     % min
    t_total_hr = N_passes * t_pass / 60; % hr
    
    fprintf('| %8.4f | %8.0f | %7.0f | %9.1f | %9.1f |\n', ...
            v, cycles_pass, N_passes, t_pass, t_total_hr);
end

fprintf('─%s─%s─%s─%s─%s─\n', repmat('─',1,12), repmat('─',1,10), ...
                           repmat('─',1,8), repmat('─',1,12), repmat('─',1,10));
fprintf('| lc_avg: %.4f m | Total cycles: %d | Budget: %.0f hr |\n', ...
        lc_avg_dynamic, n, Tavail/3600);
fprintf('═%s═%s═%s═%s═%s═\n', repmat('═',1,12), repmat('═',1,10), ...
                           repmat('═',1,8), repmat('═',1,12), repmat('═',1,10));

%% Functions
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
    % pi*r_eff^2 = b*lc => r_eff = sqrt(a_contact/pi)
    rp_eff = sqrt(A_contact/pi); %equivilant circular radius interfacing with chen et al
end
