clc; clear; close all;

%% ================= SOIL =================
kc   = 1.4e6;
kphi = 8.3e6;
n    = 1.1;
phi = deg2rad(35);
c_soil = 170;

g = 1.62;

W = 64.4 * g;
Nw = 6;

%% ================= DESIGN SPACE =================
D_vec = linspace(0.10,0.25,40);
b_vec = linspace(0.08,0.25,40);

[Dg, bg] = meshgrid(D_vec, b_vec);

%% ================= FIXED DESIGN ASSUMPTIONS =================
h = 0.008;
n_g = 18;

eta = 0.85 * 0.9;
v = 0.1;

%% ================= OUTPUT MATRICES =================
DBP_smooth = zeros(size(Dg));
DBP_grouser = zeros(size(Dg));
MU_grouser = zeros(size(Dg));
E_per_m = zeros(size(Dg));

%% ================= MODEL SWEEP =================
for i = 1:numel(Dg)

    D = Dg(i);
    b = bg(i);

    R = D/2;

    k = kc + kphi*b;

    z = (W/(k*D))^(1/(n+1));
    l = sqrt(R*z);

    A = b*l;

    sigma = k*z^n;
    tau = c_soil + sigma*tan(phi);

    %% ---- BASE (NO GROUSERS) ----
    F_base = tau * A;

    %% ---- GROUSER ENHANCEMENT ----
    p_g = (2*pi*R)/n_g;
    w_g = 0.25*p_g;

    mob = tanh(3*h/max(z,1e-6));

    F_g = n_g * (tau*(w_g*b)*mob);

    F_total = F_base + F_g;

    DBP_grouser(i) = Nw * F_total;
    DBP_smooth(i)  = Nw * F_base;

    MU_grouser(i) = F_total / W;

    %% ---- ENERGY PER METER ----
    F_resist = DBP_grouser(i);

    E_per_m(i) = F_resist / eta;   % J/m
end

%% ================= PLOT 1: GROUSER EFFECT =================
figure;
surf(Dg, bg, DBP_smooth, 'EdgeColor','none'); hold on;
surf(Dg, bg, DBP_grouser, 'EdgeColor','none','FaceAlpha',0.7);

xlabel('Wheel Diameter (m)');
ylabel('Wheel Width (m)');
zlabel('DBP (N)');
title('Effect of Grousers on Drawbar Pull');
legend('Smooth Wheel','With Grousers');
grid on; view(135,25);

%% ================= HEATMAP: DBP =================
figure;
imagesc(D_vec, b_vec, DBP_grouser);
set(gca,'YDir','normal');
colorbar;
xlabel('Wheel Diameter (m)');
ylabel('Wheel Width (m)');
title('DBP Heatmap (With Grousers)');

%% ================= HEATMAP: TRACTION COEFFICIENT =================
figure;
imagesc(D_vec, b_vec, MU_grouser);
set(gca,'YDir','normal');
colorbar;
xlabel('Wheel Diameter (m)');
ylabel('Wheel Width (m)');
title('Traction Coefficient (\mu) Heatmap');

%% ================= HEATMAP: ENERGY PER METER =================
figure;
imagesc(D_vec, b_vec, E_per_m);
set(gca,'YDir','normal');
colorbar;
xlabel('Wheel Diameter (m)');
ylabel('Wheel Width (m)');
title('Energy per Meter (J/m) Heatmap');

%% ================= GROUSER EFFECT PLOT (SIMPLE) =================
figure;
plot(D_vec, mean(DBP_grouser,1), 'LineWidth',2); hold on;
plot(D_vec, mean(DBP_smooth,1), '--','LineWidth',2);

xlabel('Wheel Diameter (m)');
ylabel('DBP (N)');
legend('With Grousers','Smooth Wheel');
title('Why Grousers Matter (Averaged Effect)');
grid on;

%% ================= POWER MAP (2D) =================
DBP_range = linspace(min(DBP_smooth(:)), max(DBP_grouser(:)), 50);
v_range   = linspace(0.05, 0.5, 50);

[DBPg, vg] = meshgrid(DBP_range, v_range);

P = (DBPg .* vg) ./ eta;   % Watts

P_kW = P / 1000;

%% ================= HEATMAP =================
figure;
imagesc(DBP_range, v_range, P_kW);

set(gca,'YDir','normal');
colorbar;

xlabel('Drawbar Pull (N)');
ylabel('Velocity (m/s)');
title('Power Requirement Map (DBP vs Velocity)');

%% ================= HEATMAP =================
figure;

imagesc(DBP_range, v_range, P_kW);
set(gca,'YDir','normal');

cb = colorbar;
ylabel(cb,'Power (kW)', ...
    'FontName','Times New Roman', ...
    'FontSize',12);

xlabel('Drawbar Pull (N)', ...
    'FontName','Times New Roman', ...
    'FontSize',12);

ylabel('Velocity (m/s)', ...
    'FontName','Times New Roman', ...
    'FontSize',12);

title('Power Requirement Map (DBP vs Velocity)', ...
    'FontName','Times New Roman', ...
    'FontSize',12);

hold on;

%% ================= CONTOURS =================
levels = [0.1 0.25 0.5 1 2 3]; % kW

[C,h] = contour(DBPg, vg, P_kW, levels, 'k', 'LineWidth',1.2);

clabel(C,h,'FontSize',10,'Color','k');
 