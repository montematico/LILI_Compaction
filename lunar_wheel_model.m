function result = lunar_wheel_model(W, D, b, kwheel, slip, bulldozing_factor, kc, kphi, n, c, phi_deg, K)
% LUNAR_WHEEL_MODEL
% Simple first-pass terramechanics model for a lunar rover wheel.
%
% Outputs a struct with sinkage, traction, resistances, and drawbar pull.
%
% This is a beginner-friendly model:
% - Bekker pressure-sinkage for sinkage
% - simple compaction resistance estimate
% - simple bulldozing approximation
% - grade resistance
% - simplified traction using Janosi-Hanamoto shear law
%
%% -----------------------------
% 1) USER INPUTS
% ------------------------------
% Rover
% m = 20; % rover mass [kg]
% Nw = 4; % number of wheels
beta_deg = 0.0; % slope angle [deg]
% Wheel geometry
% D = 0.30; % wheel diameter [m]
% b = 0.12; % wheel width [m]
% Wheel flexibility
% kwheel = 1e6; % radial stiffness [N/m]
% higher = more rigid
% Slip
% slip = 0.20; % slip ratio [0 to 1]
% Soil parameters
% kc = 1200.0; % cohesive modulus
% kphi = 15000.0; % frictional modulus
% n = 1.0; % sinkage exponent
% c = 150.0; % cohesion [Pa]
% phi_deg = 35.0; % friction angle [deg]
% K = 0.02; % shear deformation modulus [m]
% Simple bulldozing approximation
% bulldozing_factor = 0.20; % Rb = bulldozing_factor * Rc
%% -----------------------------
% 2) CONSTANTS
% ------------------------------
g_moon = 1.62; % lunar gravity [m/s^2]
beta = deg2rad(beta_deg);
phi = deg2rad(phi_deg);
%% -----------------------------
% 3) LOAD PER WHEEL
% ------------------------------
% W = m * g_moon / Nw; % normal load per wheel [N]
%% -----------------------------
% 4) BASIC GEOMETRY
% ------------------------------
r = D / 2.0; % wheel radius [m]
%% -----------------------------
% 5) WHEEL DEFLECTION
% ------------------------------
delta = W / kwheel; % wheel radial deflection [m]
%% -----------------------------
% 6) SOLVE FOR SOIL SINKAGE
% ------------------------------
% We solve:
% W / (b*l(z)) = (kc/b + kphi) * z^n
%
% where contact length:
% l(z) = 2*sqrt(r^2 - (r-z)^2)
sinkage_fun = @(z) sinkage_equation(z, W, r, b, kc, kphi, n);
z_guess = 0.01; % initial guess [m]
options = optimset('Display','off');
z_soil = fzero(sinkage_fun, z_guess, options);
% keep it physically reasonable
z_soil = max(z_soil, 0.0);
z_soil = min(z_soil, 0.99 * 2 * r);
%% -----------------------------
% 7) TOTAL EFFECTIVE SINKAGE
% ------------------------------
z_total = z_soil + delta;
%% -----------------------------
% 8) CONTACT GEOMETRY
% ------------------------------
l = contact_length(z_soil, r); % contact length [m]
A = b * l; % approximate contact area [m^2]
p = W / A; % average normal pressure [Pa]
%% -----------------------------
% 9) COMPACTION RESISTANCE
% ------------------------------
% Beginner-friendly first estimate:
% Rc = b * (kc/b + kphi) * z^(n+1)/(n+1)
Rc = b * (kc / b + kphi) * z_soil^(n + 1) / (n + 1);
%% -----------------------------
% 10) BULLDOZING RESISTANCE
% ------------------------------
Rb = bulldozing_factor * Rc;
%% -----------------------------
% 11) GRADE RESISTANCE
% ------------------------------
Rg = W * sin(beta);
%% -----------------------------
% 12) TRACTION MODEL
% ------------------------------
% Simple beginner approximation:
% j = slip * l
% tau = (c + p*tan(phi)) * (1 - exp(-j/K))
% H = tau * A
j = slip * l; % approximate shear displacement [m]
tau = (c + p * tan(phi)) * (1.0 - exp(-j / K));
H = tau * A;
%% -----------------------------
% 13) DRAWBAR PULL
% ------------------------------
DP = H - (Rc + Rb + Rg);
%% -----------------------------
% 14) STORE OUTPUTS
% ------------------------------
result.W_per_wheel_N = W;
result.radius_m = r;
result.soil_sinkage_m = z_soil;
result.wheel_deflection_m = delta;
result.total_effective_sinkage_m = z_total;
result.contact_length_m = l;
result.contact_area_m2 = A;
result.avg_pressure_Pa = p;
result.shear_displacement_m = j;
result.shear_stress_Pa = tau;
result.gross_traction_H_N = H;
result.compaction_resistance_Rc_N = Rc;
result.bulldozing_resistance_Rb_N = Rb;
result.grade_resistance_Rg_N = Rg;
result.drawbar_pull_DP_N = DP;
%% -----------------------------
% 15) PRINT RESULTS
% ------------------------------
% fprintf('\nLUNAR WHEEL MODEL RESULTS\n');
% fprintf('-------------------------------\n');
% fprintf('Wheel load W = %.6f N\n', result.W_per_wheel_N);
% fprintf('Wheel radius r = %.6f m\n', result.radius_m);
% fprintf('Soil sinkage z_soil = %.6f m\n', result.soil_sinkage_m);
% fprintf('Wheel deflection delta = %.6f m\n', result.wheel_deflection_m);
% fprintf('Total effective sinkage = %.6f m\n', result.total_effective_sinkage_m);
% fprintf('Contact length l = %.6f m\n', result.contact_length_m);
% fprintf('Contact area A = %.6f m^2\n', result.contact_area_m2);
% fprintf('Average pressure p = %.6f Pa\n', result.avg_pressure_Pa);
% fprintf('Gross traction H = %.6f N\n', result.gross_traction_H_N);
% fprintf('Compaction resistance Rc = %.6f N\n', result.compaction_resistance_Rc_N);
% fprintf('Bulldozing resistance Rb = %.6f N\n', result.bulldozing_resistance_Rb_N);
% fprintf('Grade resistance Rg = %.6f N\n', result.grade_resistance_Rg_N);
% fprintf('Drawbar pull DP = %.6f N\n', result.drawbar_pull_DP_N);
end
%% =========================================
% LOCAL FUNCTIONS
% ==========================================
function f = sinkage_equation(z, W, r, b, kc, kphi, n)
z = max(z, 1e-8);
l = contact_length(z, r);
A = b * l;
p_load = W / A;
p_bekker = (kc / b + kphi) * z^n;
f = p_load - p_bekker;
end
function l = contact_length(z, r)
z = max(z, 0.0);
z = min(z, 0.9999 * 2 * r);
inside = r^2 - (r - z)^2;
inside = max(inside, 1e-12);
l = 2.0 * sqrt(inside);
end