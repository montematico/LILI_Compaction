clc; clear; close all;

%% ================= SOIL =================
% Standardized SI constants from ChenSweep.m
SOIL.kc    = 1400;   % N/m^2 - Cohesive modulus
SOIL.kphi  = 830000; % N/m^3 - Frictional modulus
SOIL.n     = 1.0;    % Dimensionless - Sinkage exponent
SOIL.phi   = 0.576;  % rad - Internal friction angle
SOIL.c     = 170;    % Pa - Cohesion
SOIL.gamma = 2470;   % N/m^3 - Soil density/gravity
SOIL.K     = 0.02;   % m - Shear deformation modulus

g_moon = 1.62;       % m/s^2 - Lunar gravity
m_rov = 64.4;        % kg - Nominal rover mass
W_total = m_rov * g_moon; % N - Total rover weight
Nw = 6;              % Dimensionless - Number of wheels

%% ================= DESIGN SPACE =================
D_vec = linspace(0.10, 0.40, 30); % m - Wheel Diameter
b_vec = linspace(0.10, 0.40, 30); % m - Wheel Width

[Dg, bg] = meshgrid(D_vec, b_vec);

%% ================= FIXED DESIGN ASSUMPTIONS =================
h_g = 0.008;       % m - Grouser height
n_g = 18;          % Dimensionless - Number of grousers
slip = 0.20;       % Dimensionless - Target slip ratio
v = 0.1;           % m/s - Nominal velocity

%% ================= OUTPUT MATRICES =================
DBP_grouser = zeros(size(Dg));
MU_max = zeros(size(Dg));
Power_req = zeros(size(Dg));
z_soil = zeros(size(Dg));

%% ================= MODEL SWEEP =================
fprintf('Running wheel performance sweep...\n');
for i = 1:numel(Dg)
    D = Dg(i);
    b = bg(i);
    
    W_wheel = W_total / Nw; % N - Load per wheel
    
    % Call unified physics function
    perf = calculate_wheel_performance(W_wheel, D, b, h_g, n_g, slip, v, SOIL);
    
    % Store Results
    DBP_grouser(i) = perf.DBP * Nw; % Total DBP for rover
    MU_max(i) = perf.MU_max;
    Power_req(i) = perf.Power_req * Nw; % Total power for rover
    z_soil(i) = perf.z_soil;
end

%% ================= VISUALIZATION =================

% 1. 3D Surface: Drawbar Pull
figure('Name', 'Drawbar Pull Analysis', 'Color', 'w');
surf(Dg, bg, DBP_grouser, 'EdgeColor','none');
xlabel('Wheel Diameter (m)'); ylabel('Wheel Width (m)'); zlabel('Total DBP (N)');
title('Net Drawbar Pull vs. Wheel Geometry');
colorbar; grid on; view(135, 25);

% 2. Heatmap: Traction Coefficient
figure('Name', 'Traction Coefficient', 'Color', 'w');
imagesc(D_vec, b_vec, MU_max);
set(gca, 'YDir', 'normal');
xlabel('Wheel Diameter (m)'); ylabel('Wheel Width (m)');
title('Max Traction Coefficient (\mu_{max})');
cb = colorbar; ylabel(cb, '\mu_{max}');

% 3. Heatmap: Power Requirement
figure('Name', 'Locomotion Power', 'Color', 'w');
imagesc(D_vec, b_vec, Power_req);
set(gca, 'YDir', 'normal');
xlabel('Wheel Diameter (m)'); ylabel('Wheel Width (m)');
title('Total Locomotion Power (W)');
cb = colorbar; ylabel(cb, 'Watts');

% 4. Heatmap: Soil Sinkage
figure('Name', 'Soil Sinkage', 'Color', 'w');
imagesc(D_vec, b_vec, z_soil * 1000);
set(gca, 'YDir', 'normal');
xlabel('Wheel Diameter (m)'); ylabel('Wheel Width (m)');
title('Wheel Sinkage (mm)');
cb = colorbar; ylabel(cb, 'mm');

fprintf('Analysis complete.\n');
