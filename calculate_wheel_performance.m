function perf = calculate_wheel_performance(W_wheel, D, b, h_g, n_g, slip, v, SOIL)
% CALCULATE_WHEEL_PERFORMANCE 
% Consolidated terramechanics model for lunar rover wheels with grousers.
%
% Inputs:
%   W_wheel : Vertical load per wheel [N]
%   D       : Wheel diameter [m]
%   b       : Wheel width [m]
%   h_g     : Grouser height [m]
%   n_g     : Number of grousers [dimensionless]
%   slip    : Slip ratio [0 to 1]
%   v       : Velocity [m/s]
%   SOIL    : Struct with kc, kphi, n, c, phi (rad), gamma, K (shear mod)
%
% Outputs:
%   perf    : Struct with MU_max, DBP, Power_req, z_soil, etc.

%% 1. Constants & Setup
r = D / 2;
eta_drive = 0.85 * 0.9; % Motor * Transmission efficiency
c_f = 0.05; % Coefficient of rolling friction (internal resistance)

%% 2. Solve for Soil Sinkage (Bekker)
% W / (b*l(z)) = (kc/b + kphi) * z^n
sinkage_fun = @(z) sinkage_eq(z, W_wheel, r, b, SOIL.kc, SOIL.kphi, SOIL.n);
z_guess = 0.01;
options = optimset('Display','off');
try
    z = fzero(sinkage_fun, z_guess, options);
catch
    z = 0.01; % Fallback
end
z = max(min(z, 0.5 * D), 1e-6); % Physical limits

%% 3. Contact Geometry
l = 2 * sqrt(r^2 - (r - z)^2); % Contact length
A = b * l; % Contact area
p = W_wheel / A; % Average pressure

%% 4. Base Traction (Janosi-Hanamoto)
j = slip * l; % Shear displacement
tau_max = SOIL.c + p * tan(SOIL.phi);
tau = tau_max * (1 - exp(-j / SOIL.K));
F_base = tau * A;

%% 5. Grouser Enhancement (Tim's Model)
% Mobility factor based on sinkage
mob = tanh(3 * h_g / z);
% Effective grouser width/spacing contribution
p_g = (2 * pi * r) / n_g; % Pitch
w_g = 0.25 * p_g; % Effective face width (assumed)
F_g = n_g * (tau * (w_g * b) * mob);

F_gross = F_base + F_g;

%% 6. Motion Resistance
% Compaction resistance
Rc = b * (SOIL.kc / b + SOIL.kphi) * z^(SOIL.n + 1) / (SOIL.n + 1);
% Bulldozing (Simplified factor based on DRIVE struct defaults)
Rb = 0.20 * Rc; 
% Internal/Rolling resistance
Rr = W_wheel * c_f;

R_total = Rc + Rb + Rr;

%% 7. Outputs
perf.z_soil = z;
perf.MU_max = F_gross / W_wheel; % Maximum traction coefficient
perf.DBP = F_gross - R_total; % Net pull per wheel
perf.Power_req = (F_gross * v) / eta_drive; % Power per wheel
perf.R_wheel = R_total; % Resistance per wheel

end

function f = sinkage_eq(z, W, r, b, kc, kphi, n)
    z = max(z, 1e-8);
    l = 2 * sqrt(r^2 - (r - max(0, r-z))^2);
    A = b * l;
    p_load = W / A;
    p_bekker = (kc / b + kphi) * z^n;
    f = p_load - p_bekker;
end
