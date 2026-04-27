import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.widgets import Slider

# --- 1. Physics Constants ---
g_m = 1.62  # m/s^2 (Lunar gravity)
m_r_default = 75.0  # kg (Total Rover mass)
W_m_default = m_r_default * g_m  # N 

# Static Weight Distribution (Used ONLY for static ground pressure)
# 0.30 means 30% of the weight rests on the roller, 70% on the drive wheels
roller_fraction = 0.30  

# Decoupled Contact Areas
N_wheels = 6
A_wheel = 0.0025 # m^2 (Wheel contact patch)

# Roller specifications
w_drum = 0.5 # m
d_drum = 0.25 # m
l_contact = 0.05 # m
A_roller = w_drum * l_contact # 0.025 m^2

# --- 2. Dynamic Empirical Models ---
def relative_density(Fc, f=70.0):
    """Compaction depends on absolute force AND frequency (Chen et al. 2021)."""
    D_min, D_max, F_ref = 0.35, 0.95, 200.0
    eta = 1 / (1 + np.exp(-0.5 * (f - 15)))
    F_eff = Fc * eta
    return D_max - (D_max - D_min) * np.exp(-F_eff / F_ref)

def critical_rebound_limit(Dr, m_r=m_r_default):
    """Traction loss limit based on the wheels' static share of the rigid rover's weight."""
    W_total = m_r * g_m
    W_wheels = W_total * (1 - roller_fraction) # Wheels support 70% of the static load
    P_wheel = W_wheels / (N_wheels * A_wheel) 
    
    ku_loose, ku_dense = 12e6, 60e6 # MN/m^3
    ku = ku_loose + (ku_dense - ku_loose) * ((Dr - 0.35) / (0.95 - 0.35))**2
    return (P_wheel / ku) * 1000 # mm

def micro_bounce(Fc, f, m_r=m_r_default, Dr=0.65):
    """1-DOF Rigid Body Model. The entire rover mass vibrates together."""
    w = 2 * np.pi * f
    W_total = m_r * g_m
    
    # RIGID BODY ASSUMPTION: The entire mass of the rover oscillates
    m_dynamic = m_r 
    
    # Soil stiffness under the ROLLER
    ku_loose, ku_dense = 12e6, 60e6
    ku = ku_loose + (ku_dense - ku_loose) * ((Dr - 0.35) / (0.95 - 0.35))**2
    k_soil = ku * A_roller
    
    # Soil Damping (using the full rigid mass)
    c_soil = 2 * 0.3 * np.sqrt(k_soil * m_dynamic)
    
    # Steady-state continuous contact amplitude (Regime A)
    den = np.sqrt((k_soil - m_dynamic * w**2)**2 + (c_soil * w)**2)
    Y = Fc / den
    
    # Dynamic transmitted force
    F_d = Y * np.sqrt(k_soil**2 + (c_soil * w)**2)
    
    # RIGID BODY LIFT-OFF: Dynamic force must exceed the TOTAL rover weight
    X = np.where(F_d <= W_total, Y, Y + (F_d - W_total) / (m_dynamic * w**2))
    return X * 1000  # Convert to mm

def micro_bounce_from_mer(mer, f, m_r=m_r_default, Dr=0.65):
    """Dynamic displacement isolated by Eccentric Moment (me * r)."""
    w = 2 * np.pi * f
    Fc = mer * w**2
    return micro_bounce(Fc, f, m_r, Dr)

# PLOT A: Theoretical Heuristic Upper Bound (m_e * r)
mer_arr = np.linspace(0.001, 0.025, 500)
figA, axA = plt.subplots(figsize=(10, 6))

frequencies_plotA = np.linspace(10, 80, 8)
colors_plotA = plt.cm.RdYlGn_r(np.linspace(0, 1, len(frequencies_plotA)))

for freq, color in zip(frequencies_plotA, colors_plotA):
    axA.plot(mer_arr, micro_bounce_from_mer(mer_arr, freq), color=color, lw=2, label=f'{freq:.0f} Hz', alpha=0.8)

# The Asymptote
asymptote = (mer_arr / m_r_default) * 1000 
axA.plot(mer_arr, asymptote, 'k--', lw=3, label='Theoretical Heuristic Upper Bound ($\\omega \\rightarrow \\infty$)')

axA.set_xlabel('Eccentric Moment, $m_e \\times r$ (kg$\\cdot$m)')
axA.set_ylabel('Micro-Bounce Displacement (mm)')
axA.set_title('Heuristic Bounce Limit vs. Eccentric Moment (Independent of Speed)')
axA.grid(True)
axA.legend(loc='upper left', fontsize=8)
plt.tight_layout()

# PLOT B: Displacement vs Centrifugal Force
Fc_abs = np.linspace(0, 500, 500)
figB, axB = plt.subplots(figsize=(8, 5))

frequencies_plotB = np.linspace(30, 80, 10)
colors_plotB = plt.cm.RdYlGn_r(np.linspace(0, 1, len(frequencies_plotB)))

for freq, color in zip(frequencies_plotB, colors_plotB):
    axB.plot(Fc_abs, micro_bounce(Fc_abs, freq), color=color, lw=2, label=f'{freq:.0f} Hz', alpha=0.8)

axB.axvline(x=W_m_default, color='k', linestyle='--', label=f'Lunar Weight ({W_m_default:.1f} N)')
axB.set_title(f'Max Displacement vs. Absolute Centrifugal Force ({m_r_default}kg Rover)')
axB.set_xlabel('Centrifugal Force, Fc (N)')
axB.set_ylabel('Max Displacement (mm)')
axB.set_ylim(0, 1.0)
axB.grid(True)
axB.legend(fontsize=8)
plt.tight_layout()

# PLOT C: Mass Independence - Displacement vs Force Ratio 
Fc_ratio = np.linspace(1, 15, 500)
figC, axC = plt.subplots(figsize=(8, 5))

masses_plotC = [20, 25, 35, 50, 65, 80, 100]
colors_plotC = plt.cm.viridis(np.linspace(0, 1, len(masses_plotC)))
f_test = 70 

for m, c in zip(masses_plotC, colors_plotC):
    W = m * g_m
    Fc_scaled = Fc_ratio * W
    axC.plot(Fc_ratio, micro_bounce(Fc_scaled, f_test, m), color=c, lw=2.5, alpha=0.8, label=f'Rover: {m}kg')

axC.set_title(f'Displacement vs. Force Ratio @ {f_test} Hz\n(Curves converge as mass-inertia dominates)')
axC.set_xlabel('Force-to-Weight Ratio ($F_c / W_{moon}$)')
axC.set_ylabel('Max Displacement (mm)')
axC.grid(True)
axC.legend(fontsize=8)
plt.tight_layout()


# PLOT D: Density vs Force Ratio for Different Weights 

figD, axD = plt.subplots(figsize=(8, 5))

for m, c in zip(masses_plotC, colors_plotC):
    W = m * g_m
    Fc_array = Fc_ratio * W
    Dr_array = relative_density(Fc_array, f=70)
    axD.plot(Fc_ratio, Dr_array * 100, color=c, lw=2, label=f'Rover Weight: {m} kg', alpha=0.8)

axD.set_xlabel('Force-to-Weight Ratio ($F_c / W_{moon}$)')
axD.set_ylabel('Achieved Relative Density (%)')
axD.set_title('Absolute Force dictates Compaction (Heavier rovers are more efficient)')
axD.grid(True)
axD.legend(fontsize=8)
plt.tight_layout()



Dr_array = np.linspace(0.35, 0.95, 100)

fig11, ax11 = plt.subplots(figsize=(8, 5))

for m, color in zip(masses_plotC, colors_plotC):
    X_lim_array = critical_rebound_limit(Dr_array, m)
    ax11.plot(Dr_array * 100, X_lim_array, color=color, lw=2.5, label=f'Limit ({m} kg)')

ax11.set_xlabel('Relative Density, $D_r$ (%)')
ax11.set_ylabel('Max Allowable Micro-Bounce (mm)')
ax11.set_title('Traction Loss Limit as a Function of Soil Density')
ax11.grid(True, alpha=0.4)
ax11.legend()
fig11.tight_layout()


# PLOT E: Heatmap of Passes vs Centrifugal Force and Relative Density

Fc_vals = np.linspace(10, 500, 200) 
Dr_vals = np.linspace(0.35, 0.95, 200)
F_grid_E, D_grid_E = np.meshgrid(Fc_vals, Dr_vals)

# Calculate asymptotic limit for each Fc at 30 Hz
Dr_inf = relative_density(F_grid_E, f=30.0)

# Compaction rate constant (~90% achieved in 6 passes)
k_compaction = 0.4 

# Mask out unreachable densities (Yield Stress Limit)
reachable_mask = D_grid_E < Dr_inf
N_grid = np.full_like(F_grid_E, np.nan)

# Calculate passes for reachable area
fraction = (D_grid_E[reachable_mask] - 0.35) / (Dr_inf[reachable_mask] - 0.35)
fraction = np.clip(fraction, 0, 0.9999) 
N_grid[reachable_mask] = - (1 / k_compaction) * np.log(1 - fraction)

figE, axE = plt.subplots(figsize=(10, 6))
cmapE = plt.get_cmap('plasma_r').copy()
cmapE.set_bad(color='dimgray') # Unachievable zone

imE = axE.pcolormesh(F_grid_E, D_grid_E * 100, N_grid, cmap=cmapE, shading='auto', vmin=0, vmax=15)
axE.plot(Fc_vals, relative_density(Fc_vals, f=30.0) * 100, 'k--', lw=3, label='Absolute Max Density Limit ($N \\rightarrow \\infty$) @ 30 Hz')

# Add faded dotted lines for max density at 10-25 Hz
linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 5)), (0, (1, 10))]
for idx, f in enumerate([15, 16, 17, 18, 20, 25]):
    Dr_inf_other = relative_density(Fc_vals, f=f) * 100
    axE.plot(Fc_vals, Dr_inf_other, 'gray', linestyle=linestyles[idx], alpha=0.4, linewidth=1.5, label=f'Max Density Limit @ {f} Hz')

axE.set_xlabel('Centrifugal Force, Fc (N)')
axE.set_ylabel('Target Relative Density (%)')
axE.set_title('Compaction Effort (Passes) vs. Force & Density @ 30+ Hz\n(Gray area unreachable; Dotted lines: Max Density at 10-25 Hz)')
cbarE = figE.colorbar(imE, ax=axE, label='Number of Passes Required')
cbarE.extend = 'max'
axE.grid(True, alpha=0.3)
axE.legend(loc='lower right', fontsize=9)
plt.tight_layout()


# PLOT F: Micro-Bounce vs. Mass-Dependent Traction Limits

figF, axF = plt.subplots(figsize=(10, 7))

masses = [15, 25, 35, 50, 65, 80, 100, 130]
mass_colors = plt.cm.Blues(np.linspace(0.45, 1.0, len(masses)))
limit_colors = {m: c for m, c in zip(masses, mass_colors)}

for m in masses:
    Fc = Fc_ratio * m * g_m
    Dr = relative_density(Fc, f=70)
    X_limit = critical_rebound_limit(Dr, m)
    axF.plot(Fc_ratio, X_limit, color=limit_colors[m], linestyle='--', lw=2.5, label=f'Traction Limit ({m}kg)')

frequencies_plotF = np.linspace(10, 80, 8)
freq_colors_plotF = plt.cm.RdYlGn_r(np.linspace(0, 1, len(frequencies_plotF)))

for f, c in zip(frequencies_plotF, freq_colors_plotF):
    w = 2 * np.pi * f
    bounce = (np.maximum(0, Fc_ratio - 1) * g_m) / (w**2) * 1000
    axF.plot(Fc_ratio, bounce, color=c, lw=2.0, alpha=0.85, label=f'Actual Bounce @ {f:.0f} Hz')
    
    for m in masses:
        Fc = Fc_ratio * m * g_m
        Dr = relative_density(Fc, f=70)
        X_limit = critical_rebound_limit(Dr, m)
        diff = X_limit - bounce
        if np.any(diff < 0):
            idx = np.where(diff < 0)[0][0]
            axF.scatter(Fc_ratio[idx], bounce[idx], color=limit_colors[m], edgecolor=c, linewidth=2, s=70, zorder=5)

axF.set_xlabel('Force-to-Weight Ratio ($F_c / W_{moon}$)')
axF.set_ylabel('Displacement (mm)')
axF.set_ylim(0, 0.7)
axF.set_xlim(1, 15)
axF.set_title('Micro-Bounce vs. Mass-Dependent Traction Limits\n(Intersections mark Traction Loss)')
axF.grid(True, alpha=0.4)
axF.legend(loc='upper right', fontsize=8, ncol=3)
plt.tight_layout()


# PLOT G: 3D Boundary Surface (Mass vs Freq vs Max Ratio)

figG = plt.figure(figsize=(10, 7))
axG = figG.add_subplot(111, projection='3d')

M_grid_3d, F_grid_3d = np.meshgrid(np.linspace(10, 130, 40), np.linspace(10, 50, 40))
Max_R_grid = np.zeros_like(M_grid_3d)

for i in range(M_grid_3d.shape[0]):
    for j in range(M_grid_3d.shape[1]):
        m = M_grid_3d[i, j]
        f = F_grid_3d[i, j]
        
        low, high = 1.0, 30.0
        for _ in range(25): 
            mid = (low + high) / 2
            b = micro_bounce(mid * m * g_m, f, m)
            l = critical_rebound_limit(relative_density(mid * m * g_m, f=f), m)
            if b < l:
                low = mid
            else:
                high = mid
        Max_R_grid[i, j] = low

surfG = axG.plot_surface(M_grid_3d, F_grid_3d, Max_R_grid, cmap='viridis', edgecolor='none', alpha=0.9)
axG.set_title('Preliminary Operating Map:\nMax Allowable Force Ratio vs. Mass & Frequency')
axG.set_xlabel('Rover Mass (kg)')
axG.set_ylabel('Frequency (Hz)')
axG.set_zlabel('Max Allowable Ratio ($F_c/W$)')
axG.set_zlim(1, 20)
figG.colorbar(surfG, shrink=0.5, aspect=5, label='Max Ratio Threshold')
plt.tight_layout()

# PLOT H: Interactive Operating Map (Mass Slider)
figH, axH = plt.subplots(figsize=(10, 7))
figH.subplots_adjust(left=0.10, right=0.84, bottom=0.22, top=0.93)
axH_force = axH.twinx()

frequencies_grid = np.linspace(10, 80, 100)
ratios_grid = np.linspace(1, 15, 100)
F_grid, R_grid = np.meshgrid(frequencies_grid, ratios_grid)

cmapH = plt.get_cmap('coolwarm')
normH = mcolors.TwoSlopeNorm(vmin=0, vcenter=100, vmax=200)

def update_force_axis(m):
    axH_force.clear()
    axH_force.set_ylim(axH.get_ylim())
    axH_force.yaxis.set_label_position('right')
    axH_force.yaxis.tick_right()
    yticks = axH.get_yticks()
    axH_force.set_yticks(yticks)
    axH_force.set_yticklabels([f'{tick * m * g_m:.0f}' for tick in yticks])
    axH_force.set_ylabel('Centrifugal Force (N)')
    axH_force.grid(False)

def draw_heatmap(m):
    axH.clear()
    axH_force.clear()
    axH_force.yaxis.set_label_position('right')
    axH_force.yaxis.tick_right()

    Fc_grid = R_grid * m * g_m
    B_grid = micro_bounce(Fc_grid, F_grid, m)
    L_grid = critical_rebound_limit(relative_density(Fc_grid, f=F_grid), m)
    P_grid = np.clip((B_grid / L_grid) * 100, 0, 200)

    contour = axH.contourf(F_grid, R_grid, P_grid, levels=50, cmap=cmapH, norm=normH)
    axH.contour(F_grid, R_grid, P_grid, levels=[100], colors='k', linewidths=2.0, linestyles='--')

    axH.set_xlim(10, 80)
    axH.set_ylim(1, 10)
    axH.set_title(f'Operating Map ({m:.1f} kg Rover)')
    axH.set_xlabel('Frequency (Hz)')
    axH.set_ylabel('Force-to-Weight Ratio ($F_c / W_{moon}$)')
    axH.grid(True, alpha=0.25)

    update_force_axis(m)
    return contour

initial_mass = 65.0
contourH = draw_heatmap(initial_mass)
# Place colorbar in a dedicated axis to keep it fully outside the plot area.
caxH = figH.add_axes([0.91, 0.22, 0.02, 0.71])
cbarH = figH.colorbar(contourH, cax=caxH, label='Percent of Limit (%)')
cbarH.ax.axhline(100, color='k', lw=2, linestyle='--')

ax_mass = plt.axes([0.15, 0.08, 0.65, 0.03])
mass_slider = Slider(
    ax=ax_mass,
    label='Rover Mass (kg)',
    valmin=10.0,
    valmax=130.0,
    valinit=initial_mass,
)

def update_plot(_):
    contour = draw_heatmap(mass_slider.val)
    cbarH.update_normal(contour)
    figH.canvas.draw_idle()

mass_slider.on_changed(update_plot)




# ============================================================
# NEW SUPPORT FUNCTIONS FOR GLOBAL HOP BOUNDARY PLOTS
# ============================================================

def dynamic_response(Fc, f, m_r=m_r_default, Dr=0.65):
    """
    1-DOF rigid-body vertical oscillator.
    Returns:
      Y_mm   = steady-state displacement amplitude [mm]
      Fd_N   = transmitted dynamic force amplitude [N]
      Rhop   = global hop severity = Fd / W_total
    """
    w = 2.0 * np.pi * f
    W_total = m_r * g_m

    # Rigid-body assumption: entire rover oscillates together
    m_dynamic = m_r

    ku_loose, ku_dense = 12e6, 60e6 # MN/m^3
    ku = ku_loose + (ku_dense - ku_loose) * ((Dr - 0.35) / (0.95 - 0.35))**2
    k_soil = ku * A_roller
    c_soil = 2.0 * 0.3 * np.sqrt(k_soil * m_dynamic)

    den = np.sqrt((k_soil - m_dynamic * w**2)**2 + (c_soil * w)**2)
    Y = Fc / den   # [m]

    # transmitted dynamic force amplitude
    Fd = Y * np.sqrt(k_soil**2 + (c_soil * w)**2)

    Rhop = Fd / W_total
    return Y * 1000.0, Fd, Rhop

def hop_metric(Fc, f, m_r):
    """
    Global rigid-body hop severity ratio.
    Rhop = 1 is the rigid-body hop onset boundary.
    """
    Dr = relative_density(Fc, f=f)
    _, _, Rhop = dynamic_response(Fc, f, m_r=m_r, Dr=Dr)
    return Rhop


# ============================================================
# NEW PLOT F-HOP:
# Like Plot F, but markers show GLOBAL HOP ONSET instead of
# wheel traction limit crossings.
# ============================================================

figF_hop, axF_hop = plt.subplots(figsize=(10, 7))

Fc_ratio_hop = np.linspace(1, 15, 600)
frequencies_plotF_hop = np.linspace(10, 80, 8)
freq_colors_plotF_hop = plt.cm.RdYlGn_r(np.linspace(0, 1, len(frequencies_plotF_hop)))

m_ref_hop = 65.0  # choose a reference rover mass for the plot
W_ref_hop = m_ref_hop * g_m

for f, c in zip(frequencies_plotF_hop, freq_colors_plotF_hop):
    Fc_array = Fc_ratio_hop * W_ref_hop

    # bounce amplitude curve
    Dr_array = relative_density(Fc_array, f=f)
    B_array = micro_bounce(Fc_array, f, m_r=m_ref_hop, Dr=Dr_array)
    axF_hop.plot(Fc_ratio_hop, B_array, color=c, lw=2.2, alpha=0.9, label=f'{f:.0f} Hz')

    # hop severity
    Rhop_array = hop_metric(Fc_array, f, m_ref_hop)

    # mark first hop-onset crossing if it exists
    if np.any(Rhop_array >= 1.0):
        idx = np.where(Rhop_array >= 1.0)[0][0]
        axF_hop.scatter(
            Fc_ratio_hop[idx], B_array[idx],
            color=c, edgecolor='k', linewidth=1.5, s=80, zorder=6
        )

axF_hop.axvline(1.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
axF_hop.set_xlabel('Force-to-Weight Ratio ($F_c / W_{moon}$)')
axF_hop.set_ylabel('Micro-Bounce Displacement (mm)')
axF_hop.set_title(f'Global Hop Onset vs. Micro-Bounce ({m_ref_hop:.0f} kg Rover)\n'
                  'Markers indicate first hop-regime crossing ($R_{hop}=1$)')
axF_hop.grid(True, alpha=0.4)
axF_hop.legend(loc='upper right', fontsize=9, ncol=2)
plt.tight_layout()


# ============================================================
# NEW PLOT G-HOP:
# Frequency on x-axis, force-to-weight ratio on left y-axis,
# centrifugal force on right y-axis.
# Multiple rover masses shown as hop-limit contour lines.
# No slider, no heatmap.
# ============================================================

figG_hop, axG_hop = plt.subplots(figsize=(10, 7))
# axG_hop_force = axG_hop.twinx()

masses_hop = [15, 25, 35, 50, 65, 80, 100, 130]
markers = ['o', 's', '^', 'v', 'D', 'p', '*', 'h']

norm = mcolors.Normalize(vmin=0.35, vmax=0.95)
cmap = plt.cm.viridis

freqs_hop = np.linspace(10, 80, 250)
ratio_search = np.linspace(0, 15.0, 1200)

for m, marker in zip(masses_hop, markers):
    hop_limit_ratio = np.full_like(freqs_hop, np.nan, dtype=float)

    for i, f in enumerate(freqs_hop):
        Fc_array = ratio_search * m * g_m
        Rhop_array = hop_metric(Fc_array, f, m)

        # find first ratio where Rhop >= 1
        if np.any(Rhop_array >= 1.0):
            idx = np.where(Rhop_array >= 1.0)[0][0]
            hop_limit_ratio[i] = ratio_search[idx]

    # Create LineCollection with gradient colors
    valid_idx = ~np.isnan(hop_limit_ratio)
    x = freqs_hop[valid_idx]
    y = hop_limit_ratio[valid_idx]
    if len(x) < 2:
        continue

    # Compute Dr for each point
    dr_vals = []
    for i in range(len(x)):
        f = x[i]
        ratio = y[i]
        Fc = ratio * m * g_m
        dr = relative_density(Fc, f)
        dr_vals.append(dr)

    # Segments
    points = np.column_stack((x, y))
    segments = np.array([points[:-1], points[1:]]).transpose(1, 0, 2)

    # Colors for segments, use dr of the first point of each segment
    colors = [cmap(norm(dr)) for dr in dr_vals[:-1]]

    lc = LineCollection(segments, colors=colors, linewidth=2.5)
    axG_hop.add_collection(lc)

    # Add markers at intervals
    marker_indices = np.arange(0, len(x), 10)  # every 10th point
    for idx in marker_indices:
        if idx < len(dr_vals):
            axG_hop.scatter(x[idx], y[idx], color=cmap(norm(dr_vals[idx])), marker=marker, s=50, edgecolor='black', linewidth=0.5, zorder=5)

# For legend
for m, marker in zip(masses_hop, markers):
    axG_hop.plot([], [], marker=marker, color='black', linestyle='none', markersize=8, label=f'{m:.0f} kg')

axG_hop.legend(loc='upper right', fontsize=9, ncol=2)

# Primary axis formatting
axG_hop.set_xlim(10, 80)
axG_hop.set_ylim(0, 8)
axG_hop.set_xlabel('Frequency (Hz)')
axG_hop.set_ylabel('Force-to-Weight Ratio ($F_c/W_{r}$)')
axG_hop.set_title('Bounce Regime Boundaries vs. Frequency for Rover Masses')
axG_hop.grid(True, alpha=0.3)

# Add colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = figG_hop.colorbar(sm, ax=axG_hop, label='Relative Density after 1 Pass')

# Secondary axis: show centrifugal force corresponding to a chosen reference mass
# m_ref_axis = 65.0
# yticks = axG_hop.get_yticks()
# axG_hop_force.set_ylim(axG_hop.get_ylim())
# axG_hop_force.set_yticks(yticks)
# axG_hop_force.set_yticklabels([f'{tick * m_ref_axis * g_m:.0f}' for tick in yticks])
# axG_hop_force.set_ylabel(f'Equivalent Centrifugal Force for {m_ref_axis:.0f} kg Rover (N)')

plt.tight_layout()












# plt.show()









plt.show()