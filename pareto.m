%% ChenPareto.m  –  Pareto optimizer for rover vibratory compaction
%  Sweeps mp, f, me, re, mt. Filters bounce + physical sanity.
%  Outputs: Pareto front (mass vs power, colored by N_passes),
%           feasibility map, recommendation table.
clear; clc; close all;

%% =================== CONFIGURATION (EDIT HERE) ===================
% ---- Environment ----
g_moon   = 1.62;        % m/s^2
Apad     = 500;         % m^2
Tavail_h = 365*24*5;    % hours (2 yr + 6 mo margin)
Tavail   = Tavail_h * 3600;

% ---- Regolith ----
h_layer  = 0.1;         % m, lift thickness
rho_i    = 1.3;         % g/cm^3, initial density
rho_f    = 1.8;         % g/cm^3, target density
nu       = 0.35;        % Poisson ratio

% ---- Rover / Roller geometry ----
m_rov       = 400;      % kg, total rover mass
Roller_frac = 0.5;      % fraction of rover weight on roller
W_roller    = Roller_frac * m_rov * g_moon;  % N
b_roller    = 0.3;      % m, roller width
D_roller    = 0.4;      % m
R_roller    = D_roller/2;

% ---- Simulation ----
N_cycles = 1e5;         % max cycles per spot sim

%% =============== PARETO SWEEP CONFIG (EDIT RANGES) ===============
% Arrays  => swept.   Scalars => held fixed.
% Start with coarse grids, then narrow around the Pareto knee.
mp_sweep = logspace(log10(50),  log10(250), 8);   % kg, roller vibrating mass
f_sweep  = linspace(25, 50, 7);                   % Hz, vibration frequency
me_sweep = logspace(log10(0.05), log10(0.8), 6);  % kg, eccentric mass
re_sweep = [0.05, 0.10];                          % m,  eccentric radius
mt_sweep = Roller_frac * m_rov;                   % kg, fixed (scalar)

% ---- Physical sanity filter ----
% Reject any run where a SINGLE xcr > this fraction of the full layer.
% Prevents 1-cycle "magic compaction" artefacts from polluting the Pareto.
xcr_cap_frac = 0.15;   % e.g. 0.15 = xcr > 15% of h_layer is suspicious

%% =============== BUILD GRID & RUN SWEEP ===============
[F_grid, MP_grid, ME_grid, RE_grid, MT_grid] = ...
    ndgrid(f_sweep, mp_sweep, me_sweep, re_sweep, mt_sweep);

N_pts = numel(F_grid);

% Preallocate result arrays (flat vectors – easier for Pareto later)
res_mp      = zeros(N_pts,1);
res_f       = zeros(N_pts,1);
res_me      = zeros(N_pts,1);
res_re      = zeros(N_pts,1);
res_mt      = zeros(N_pts,1);
res_Ncyc    = nan(N_pts,1);   % cycles to reach rho_f  (= N_passes proxy)
res_t_hr    = nan(N_pts,1);   % total pad time (hr)
res_power   = nan(N_pts,1);   % power proxy:  F0^2 / (2*mt)  [W-ish]
res_rho     = zeros(N_pts,1); % final density reached
res_flag    = zeros(N_pts,1); % 0=ok, 1=bounce, 2=no yield, 3=overtime, 4=xcr_cap

W_rover = m_rov * g_moon;

fprintf('Running %d combinations...\n', N_pts);
t_start = tic;

for i = 1:N_pts
    f_i  = F_grid(i);  mp_i = MP_grid(i);
    me_i = ME_grid(i); re_i = RE_grid(i); mt_i = MT_grid(i);

    omega_i = 2*pi*f_i;
    F0_i    = me_i * re_i * omega_i^2;

    res_mp(i)=mp_i; res_f(i)=f_i; res_me(i)=me_i;
    res_re(i)=re_i; res_mt(i)=mt_i;

    % --- A. Bounce check ---
    if F0_i >= W_rover
        res_flag(i) = 1;
        continue
    end

    % --- B. Sim ---
    [rho_f_sim, t_pt, n_cyc, lc_avg, xcr_max] = run_compaction_sim(...
        f_i, me_i, re_i, mp_i, mt_i, ...
        h_layer, rho_i, rho_f, R_roller, b_roller, nu, N_cycles);

    res_rho(i) = rho_f_sim;

    % --- C. xcr sanity check ---
    if xcr_max > xcr_cap_frac * h_layer
        res_flag(i) = 4;   % suspicious – single cycle too large
        continue
    end

    % --- D. Density check ---
    if rho_f_sim < rho_f
        res_flag(i) = 2;
        continue
    end

    % --- E. Time check ---
    A_avg    = b_roller * lc_avg;
    t_pad_hr = (t_pt * (Apad / A_avg)) / 3600;
    if t_pad_hr > Tavail_h
        res_flag(i) = 3;
        continue
    end

    % --- F. Store results ---
    res_Ncyc(i)  = n_cyc;
    res_t_hr(i)  = t_pad_hr;
    % Power proxy: mechanical power of eccentric drive ~ F0^2 / (2*mt) * (1/omega)
    % Simpler dimensionally consistent proxy: P ~ F0^2 * omega  (force × velocity scale)
    res_power(i) = F0_i^2 * omega_i;   % [N^2 * rad/s]  — relative measure only
end

elapsed = toc(t_start);
fprintf('Sweep done in %.1f s.\n', elapsed);

% Diagnostics
fprintf('\nDiagnostics:\n');
fprintf('  Success:          %d\n', sum(res_flag==0 & ~isnan(res_Ncyc)));
fprintf('  Bounce:           %d\n', sum(res_flag==1));
fprintf('  No yield:         %d\n', sum(res_flag==2));
fprintf('  Overtime:         %d\n', sum(res_flag==3));
fprintf('  xcr cap (artefact filtered): %d\n\n', sum(res_flag==4));

%% =============== PARETO FILTER ===============
% Feasible set: flag==0 and not NaN
ok = (res_flag==0) & ~isnan(res_Ncyc);

if sum(ok) < 2
    warning('Fewer than 2 feasible points – widen sweep ranges.');
    return
end

% Objectives to MINIMISE: [mass, power, N_passes]
obj = [res_mp(ok), res_power(ok), res_Ncyc(ok)];

% Normalise each objective to [0,1] so units don't dominate
obj_n = (obj - min(obj)) ./ (max(obj) - min(obj) + eps);

% Find non-dominated points (true Pareto front)
n_ok = size(obj_n,1);
is_pareto = true(n_ok,1);
for p = 1:n_ok
    for q = 1:n_ok
        if p==q, continue; end
        % q dominates p if it is <= in all objectives AND < in at least one
        if all(obj_n(q,:) <= obj_n(p,:)) && any(obj_n(q,:) < obj_n(p,:))
            is_pareto(p) = false;
            break
        end
    end
end

% Pull out Pareto-front data
ok_idx        = find(ok);
pareto_idx    = ok_idx(is_pareto);
pareto_mp     = res_mp(pareto_idx);
pareto_power  = res_power(pareto_idx);
pareto_Ncyc   = res_Ncyc(pareto_idx);
pareto_f      = res_f(pareto_idx);
pareto_me     = res_me(pareto_idx);
pareto_re     = res_re(pareto_idx);
pareto_mt     = res_mt(pareto_idx);
pareto_t_hr   = res_t_hr(pareto_idx);

fprintf('%d Pareto-optimal designs found.\n\n', sum(is_pareto));

%% =============== KNEE POINT (balanced recommendation) ===============
% Minimise normalised distance to origin (ideal = 0 in all objectives)
obj_pareto_n = (obj(is_pareto,:) - min(obj)) ./ (max(obj) - min(obj) + eps);
% Priority weights matching your stated priorities: mass >> power > passes
w = [1.0, 0.5, 0.1];
score = obj_pareto_n * w(:);
[~, knee_idx] = min(score);

fprintf('===== RECOMMENDATION (Pareto knee) =====\n');
fprintf('  Roller mass  mp  = %.1f kg\n',  pareto_mp(knee_idx));
fprintf('  Frequency    f   = %.1f Hz\n',  pareto_f(knee_idx));
fprintf('  Eccentric    me  = %.3f kg\n',  pareto_me(knee_idx));
fprintf('  Ecc. radius  re  = %.3f m\n',   pareto_re(knee_idx));
fprintf('  Cycles to compact = %.0f\n',    pareto_Ncyc(knee_idx));
fprintf('  Total pad time    = %.1f hr\n', pareto_t_hr(knee_idx));
fprintf('  Power proxy       = %.3e\n',    pareto_power(knee_idx));
fprintf('=========================================\n\n');

%% =============== PLOTS ===============

% --- Figure 1: Pareto front – mass vs power, color = N_passes ---
figure('Name','Pareto Front: Mass vs Power','Position',[50 500 700 500]);
scatter(pareto_mp, pareto_power, 120, pareto_Ncyc, 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.8);
colormap(flipud(parula));   % low passes = bright (good)
cb = colorbar; ylabel(cb,'N_{cycles} to compact (fewer = better)','FontSize',10);
hold on;
% Highlight knee point
scatter(pareto_mp(knee_idx), pareto_power(knee_idx), 250, 'rp', ...
        'filled','MarkerEdgeColor','k','LineWidth',1.5);
text(pareto_mp(knee_idx), pareto_power(knee_idx), ...
     sprintf('  Recommended\n  m_p=%.0f kg, f=%.0f Hz', ...
     pareto_mp(knee_idx), pareto_f(knee_idx)), ...
     'FontSize',9,'Color','r','FontWeight','bold');
xlabel('Roller vibrating mass m_p (kg)','FontSize',12);
ylabel('Power proxy  F_0^2\cdot\omega  (relative)','FontSize',12);
title('Pareto Front: minimize mass & power,  color = cycles to compact','FontSize',11);
set(gca,'YScale','log','FontSize',11); grid on; box on;

% --- Figure 2: Feasibility map – mp vs f, color = final density ---
% Collapse other dims by taking max density achieved for each (mp,f) pair
mp_vals = unique(res_mp);   f_vals = unique(res_f);
rho_map = zeros(length(mp_vals), length(f_vals));
flag_map = ones(length(mp_vals), length(f_vals));  % default: infeasible
for ii = 1:length(mp_vals)
    for jj = 1:length(f_vals)
        sel = abs(res_mp - mp_vals(ii))<1e-6 & abs(res_f - f_vals(jj))<1e-6;
        rho_map(ii,jj)  = max(res_rho(sel));
        if any(res_flag(sel)==0 & ~isnan(res_Ncyc(sel)))
            flag_map(ii,jj) = 0;
        end
    end
end

figure('Name','Feasibility Map: mp vs f','Position',[50 50 700 450]);
imagesc(f_vals, mp_vals, rho_map);
axis xy; colormap(parula);
cb2 = colorbar; ylabel(cb2,'Max final density (g/cm^3)','FontSize',10);
clim([rho_i, rho_f]);
hold on;
% Overlay infeasible cells as hatching via low-alpha red patch
[jj_g, ii_g] = meshgrid(1:length(f_vals), 1:length(mp_vals));
inf_mask = flag_map > 0;
if any(inf_mask(:))
    scatter(f_vals(jj_g(inf_mask)), mp_vals(ii_g(inf_mask)), 60, ...
            'rx','LineWidth',2);
end
% Pareto points on this map
scatter(pareto_f, pareto_mp, 80, 'w^','filled','MarkerEdgeColor','k');
scatter(pareto_f(knee_idx), pareto_mp(knee_idx), 180, 'rp', ...
        'filled','MarkerEdgeColor','k');
xlabel('Frequency f (Hz)','FontSize',12);
ylabel('Roller mass m_p (kg)','FontSize',12);
title('Feasibility map (X = infeasible, triangles = Pareto, star = recommended)','FontSize',10);
set(gca,'FontSize',11); grid off; box on;

% --- Figure 3: Pareto table as formatted text ---
figure('Name','Pareto Design Table','Position',[760 200 560 500]);
axis off;
headers = {'mp (kg)','f (Hz)','me (kg)','re (m)','N_cyc','T_pad (hr)','Power proxy'};





%% Functions
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
