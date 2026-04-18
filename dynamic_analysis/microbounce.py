import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
from matplotlib.widgets import Slider

# --- 1. Physics Constants ---
g_m = 1.62  # m/s^2 (Lunar gravity)
m_r_default = 35.0  # kg (Default Rover mass)
W_m_default = m_r_default * g_m  # N (Static Lunar Weight)

# --- 2. Corrected Empirical Models ---
def relative_density(Fc):
    """Density is a function of absolute force applied to the soil."""
    D_min = 0.35
    D_max = 0.95
    F_ref = 200.0  # N (Empirical force scaling for significant compaction)
    return D_max - (D_max - D_min) * np.exp(-Fc / F_ref)

def critical_rebound_limit(Dr):
    """Traction loss is bounded by the elastic rebound of the soil matrix."""
    X_loose = 0.50 # mm (Springy loose soil)
    X_dense = 0.10 # mm (Rigid locked soil)
    return X_loose - (X_loose - X_dense) * ((Dr - 0.35) / (0.95 - 0.35))

def micro_bounce(Fc, f, m_r=m_r_default):
    """Standard displacement formula based on Net Force and Frequency."""
    w = 2 * np.pi * f
    W_m = m_r * g_m
    F_net = np.maximum(0, Fc - W_m)
    return (F_net / (m_r * w**2)) * 1000  # Convert to mm

def micro_bounce_from_mer(mer, f, m_r=m_r_default):
    """Displacement isolated by Eccentric Moment (me * r)."""
    w = 2 * np.pi * f
    Fc = mer * w**2
    W_m = m_r * g_m
    F_net = np.maximum(0, Fc - W_m)
    return (F_net / (m_r * w**2)) * 1000


# ==========================================
# PLOT 1: Micro-Bounce vs Traction Loss Limit (35kg Rover)
# ==========================================
Fc_ratio = np.linspace(1, 15, 1000)
Fc = Fc_ratio * W_m_default
Dr = relative_density(Fc)
X_limit = critical_rebound_limit(Dr)

fig1, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(Fc_ratio, X_limit, 'b--', lw=3, label='Traction Loss Limit (Elastic Rebound)')
frequencies_plot1 = np.linspace(10, 50, 8)
colors_plot1 = plt.cm.RdYlGn_r(np.linspace(0, 1, len(frequencies_plot1)))
for freq, color in zip(frequencies_plot1, colors_plot1):
    bounce = micro_bounce(Fc, freq)
    ax1.plot(Fc_ratio, bounce, color=color, lw=2.5, label=f'Bounce @ {freq:.0f} Hz', alpha=0.8)
    diff = X_limit - bounce
    if np.any(diff < 0):
        idx = np.where(diff < 0)[0][0]
        ax1.scatter(Fc_ratio[idx], bounce[idx], color=color, s=100, zorder=5)
ax1.set_xlabel('Force-to-Weight Ratio ($F_c / W_{moon}$)')
ax1.set_ylabel('Displacement (mm)')
ax1.set_ylim(0, 0.8)
ax1.set_title('1. Micro-Bounce vs. Elastic Soil Rebound (35kg Rover)')
ax1.grid(True, alpha=0.4)
ax1.legend(loc='upper left', fontsize=8)
fig1.tight_layout()

# ==========================================
# PLOT 2: Safe Traction Zone vs Frequency
# ==========================================
frequencies = np.linspace(10, 50, 200)
opt_ratios = []
for freq in frequencies:
    b_arr = micro_bounce(Fc, freq)
    diff = X_limit - b_arr
    if np.any(diff < 0):
        idx = np.where(diff < 0)[0][0]
        opt_ratios.append(Fc_ratio[idx])
    else:
        opt_ratios.append(15.0)

fig2, ax2 = plt.subplots(figsize=(9, 5))
ax2.plot(frequencies, opt_ratios, 'k-', lw=3)
ax2.fill_between(frequencies, 1, opt_ratios, color='green', alpha=0.2, label='Safe Traction Zone')
ax2.fill_between(frequencies, opt_ratios, 15, color='red', alpha=0.1, label='Traction Loss (Rover Bounces)')
ax2.set_xlabel('Vibration Frequency (Hz)')
ax2.set_ylabel('Maximum Allowed Force Ratio ($F_c / W_{moon}$)')
ax2.set_ylim(1, 15)
ax2.set_xlim(10, 50)
ax2.set_title('2. Traction Limit vs. Frequency (35kg Rover)')
ax2.grid(True)
ax2.legend(loc='upper left')
fig2.tight_layout()

# ==========================================
# PLOT 3: Absolute Physical Limit (m_e * r)
# ==========================================
mer_arr = np.linspace(0.001, 0.025, 500)
fig3, ax3 = plt.subplots(figsize=(10, 6))
frequencies_plot3 = np.linspace(10, 50, 8)
colors_plot3 = plt.cm.RdYlGn_r(np.linspace(0, 1, len(frequencies_plot3)))
for freq, color in zip(frequencies_plot3, colors_plot3):
    ax3.plot(mer_arr, micro_bounce_from_mer(mer_arr, freq), color=color, lw=2, label=f'{freq:.0f} Hz', alpha=0.8)
asymptote = (mer_arr / m_r_default) * 1000 
ax3.plot(mer_arr, asymptote, 'k--', lw=3, label='Absolute Physical Limit ($\\omega \\rightarrow \\infty$)')
ax3.set_xlabel('Eccentric Moment, $m_e \\times r$ (kg$\\cdot$m)')
ax3.set_ylabel('Micro-Bounce Displacement (mm)')
ax3.set_title('3. The Absolute Bounce Limit: Independent of Speed')
ax3.grid(True)
ax3.legend(loc='upper left', fontsize=8)
fig3.tight_layout()

# ==========================================
# PLOT 4: Displacement vs Centrifugal Force
# ==========================================
Fc_abs = np.linspace(0, 500, 500)
fig4, ax4 = plt.subplots(figsize=(8, 5))
frequencies_plot4 = np.linspace(10, 50, 8)
colors_plot4 = plt.cm.RdYlGn_r(np.linspace(0, 1, len(frequencies_plot4)))
for freq, color in zip(frequencies_plot4, colors_plot4):
    ax4.plot(Fc_abs, micro_bounce(Fc_abs, freq), color=color, lw=2, label=f'{freq:.0f} Hz', alpha=0.8)
ax4.axvline(x=W_m_default, color='k', linestyle='--', label=f'Lunar Weight ({W_m_default:.1f} N)')
ax4.set_title('4. Max Displacement vs. Absolute Centrifugal Force (35kg Rover)')
ax4.set_xlabel('Centrifugal Force, Fc (N)')
ax4.set_ylabel('Max Displacement (mm)')
ax4.set_ylim(0, 1.0)
ax4.grid(True)
ax4.legend(fontsize=8)
fig4.tight_layout()

# ==========================================
# PLOT 5: Mass Independence - Displacement vs Force Ratio 
# ==========================================
fig5, ax5 = plt.subplots(figsize=(8, 5))
masses_plot5 = [15, 25, 35, 50, 75, 100]
colors_plot5 = plt.cm.viridis(np.linspace(0, 1, len(masses_plot5)))
f_test = 70 
for m, c in zip(masses_plot5, colors_plot5):
    W = m * g_m
    Fc_scaled = Fc_ratio * W
    ax5.plot(Fc_ratio, micro_bounce(Fc_scaled, f_test, m), color=c, lw=2.5, alpha=0.8, label=f'Rover: {m}kg')
ax5.set_title(f'5. Displacement vs. Force Ratio @ {f_test} Hz\nBounce Height is Mass-Independent')
ax5.set_xlabel('Force-to-Weight Ratio (Fc / W_moon)')
ax5.set_ylabel('Max Displacement (mm)')
ax5.grid(True)
ax5.legend(fontsize=8)
fig5.tight_layout()

# ==========================================
# PLOT 6: Density vs Force Ratio for Different Weights 
# ==========================================
fig6, ax6 = plt.subplots(figsize=(8, 5))
masses_plot6 = [15, 25, 35, 50, 75, 100]
colors_plot6 = plt.cm.viridis(np.linspace(0, 1, len(masses_plot6)))
for m, c in zip(masses_plot6, colors_plot6):
    W = m * g_m
    Fc_array = Fc_ratio * W
    Dr_array = relative_density(Fc_array)
    ax6.plot(Fc_ratio, Dr_array * 100, color=c, lw=2, label=f'Rover Weight: {m} kg', alpha=0.8)
ax6.set_xlabel('Force-to-Weight Ratio ($F_c / W_{moon}$)')
ax6.set_ylabel('Achieved Relative Density (%)')
ax6.set_title('6. Rover Compaction by Weight')
ax6.grid(True)
ax6.legend(fontsize=8)
fig6.tight_layout()

# ==========================================
# PLOT 7: Heatmap of Passes vs Centrifugal Force and Relative Density
# ==========================================
Fc_vals = np.linspace(10, 500, 200) 
Dr_vals = np.linspace(0.35, 0.95, 200)
F_grid_7, D_grid_7 = np.meshgrid(Fc_vals, Dr_vals)

# Calculate asymptotic limit for each Fc
Dr_inf = relative_density(F_grid_7)

# Calculate required passes N
k_compaction = 0.4 # Empirical Compaction rate constant (e.g., ~90% achieved in 6 passes)

# Mask out unreachable densities
reachable_mask = D_grid_7 < Dr_inf

# Initialize N_grid with NaNs to hide unreachable areas
N_grid = np.full_like(F_grid_7, np.nan)

# Safe calculation for reachable area
fraction = (D_grid_7[reachable_mask] - 0.35) / (Dr_inf[reachable_mask] - 0.35)
fraction = np.clip(fraction, 0, 0.9999) 
N_grid[reachable_mask] = - (1 / k_compaction) * np.log(1 - fraction)

fig7, ax7 = plt.subplots(figsize=(10, 6))
cmap7 = plt.get_cmap('plasma_r').copy()
cmap7.set_bad(color='dimgray') # Color for the unachievable zone

im7 = ax7.pcolormesh(F_grid_7, D_grid_7 * 100, N_grid, cmap=cmap7, shading='auto', vmin=0, vmax=15)

# Overlay the asymptotic limit line
ax7.plot(Fc_vals, relative_density(Fc_vals) * 100, 'k--', lw=3, label='Absolute Max Density Limit ($N \\rightarrow \\infty$)')

ax7.set_xlabel('Centrifugal Force, Fc (N)')
ax7.set_ylabel('Target Relative Density (%)')
ax7.set_title('7. Compaction Effort (Passes) vs. Force & Density\n(Gray area is physically unachievable due to yield stress)')
cbar7 = fig7.colorbar(im7, ax=ax7, label='Number of Passes Required')
cbar7.extend = 'max'
ax7.grid(True, alpha=0.3)
ax7.legend(loc='lower right', fontsize=9)
fig7.tight_layout()

# ==========================================
# PLOT 8: 3D Boundary Surface (Mass vs Freq vs Max Ratio)
# ==========================================
fig8 = plt.figure(figsize=(10, 7))
ax8 = fig8.add_subplot(111, projection='3d')

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
            l = critical_rebound_limit(relative_density(mid * m * g_m))
            if b < l:
                low = mid
            else:
                high = mid
        Max_R_grid[i, j] = low

surf8 = ax8.plot_surface(M_grid_3d, F_grid_3d, Max_R_grid, cmap='viridis', edgecolor='none', alpha=0.9)
ax8.set_title('8. The Ultimate Design Envelope:\nMax Allowable Force Ratio vs. Mass & Frequency')
ax8.set_xlabel('Rover Mass (kg)')
ax8.set_ylabel('Frequency (Hz)')
ax8.set_zlabel('Max Allowable Ratio ($F_c/W$)')
ax8.set_zlim(1, 20)
fig8.colorbar(surf8, shrink=0.5, aspect=5, label='Max Ratio Threshold')
fig8.tight_layout()

# ==========================================
# PLOT 9: Interactive Safe Operating Envelope
# ==========================================
fig9, ax9 = plt.subplots(figsize=(10, 7))
plt.subplots_adjust(bottom=0.25) 
ax9_force = ax9.twinx()

frequencies_grid = np.linspace(10, 50, 100)
ratios_grid = np.linspace(1, 15, 100)
F_grid, R_grid = np.meshgrid(frequencies_grid, ratios_grid)

cmap8 = plt.get_cmap('coolwarm')
norm8 = mcolors.TwoSlopeNorm(vmin=0, vcenter=100, vmax=200)

def update_force_axis(m):
    ax9_force.clear()
    ax9_force.set_ylim(ax9.get_ylim())
    ax9_force.yaxis.set_label_position('right')
    ax9_force.yaxis.tick_right()
    yticks = ax9.get_yticks()
    ax9_force.set_yticks(yticks)
    ax9_force.set_yticklabels([f'{tick * m * g_m:.0f}' for tick in yticks])
    ax9_force.set_ylabel('Centrifugal Force (N)')
    ax9_force.grid(False)

def draw_heatmap(m):
    ax9.clear()
    ax9_force.clear()
    ax9_force.yaxis.set_label_position('right')
    ax9_force.yaxis.tick_right()
    Fc_grid = R_grid * m * g_m
    B_grid = micro_bounce(Fc_grid, F_grid, m)
    L_grid = critical_rebound_limit(relative_density(Fc_grid))
    
    P_grid = np.clip((B_grid / L_grid) * 100, 0, 200)
    
    contour = ax9.contourf(F_grid, R_grid, P_grid, 50, cmap=cmap8, norm=norm8)
    ax9.contour(F_grid, R_grid, P_grid, levels=[100], colors='k', linewidths=2.5, linestyles='--')
    ax9.set_ylim(1, 15)
    update_force_axis(m)
    
    ax9.set_title(f'9. Interactive Safe Operating Envelope ({m:.1f} kg Rover)')
    ax9.set_xlabel('Frequency (Hz)')
    ax9.set_ylabel('Force-to-Weight Ratio ($F_c / W_{moon}$)')
    return contour

initial_mass = 35.0
c_plot = draw_heatmap(initial_mass)
cbar9 = fig9.colorbar(c_plot, ax=ax9, shrink=0.8, pad=0.12, label='Percent of Limit (%)')
cbar9.ax.axhline(100, color='k', lw=2.5, linestyle='--')

ax_mass = plt.axes([0.15, 0.1, 0.65, 0.03])
mass_slider = Slider(
    ax=ax_mass,
    label='Rover Mass (kg)',
    valmin=10.0,
    valmax=130.0,
    valinit=initial_mass,
)

def update(val):
    draw_heatmap(mass_slider.val)
    fig9.canvas.draw_idle()

mass_slider.on_changed(update)

# ==========================================
# PLOT 10: Micro-Bounce vs. Mass-Dependent Traction Limits
# ==========================================
Fc_ratio = np.linspace(1, 15, 500)
fig10, ax10 = plt.subplots(figsize=(10, 7))

masses = [15, 25, 35, 50, 75, 100, 130]
mass_colors = plt.cm.Blues(np.linspace(0.45, 1.0, len(masses)))
limit_colors = {m: c for m, c in zip(masses, mass_colors)}

for m in masses:
    Fc = Fc_ratio * m * g_m
    Dr = relative_density(Fc)
    X_limit = critical_rebound_limit(Dr)
    ax10.plot(Fc_ratio, X_limit, color=limit_colors[m], linestyle='--', lw=2.5, label=f'Traction Limit ({m}kg)')

frequencies_plot9 = np.linspace(10, 50, 8)
freq_colors_plot9 = plt.cm.RdYlGn_r(np.linspace(0, 1, len(frequencies_plot9)))
for f, c in zip(frequencies_plot9, freq_colors_plot9):
    w = 2 * np.pi * f
    bounce = (np.maximum(0, Fc_ratio - 1) * g_m) / (w**2) * 1000
    ax10.plot(Fc_ratio, bounce, color=c, lw=2.0, alpha=0.85, label=f'Actual Bounce @ {f:.0f} Hz')
    
    for m in masses:
        Fc = Fc_ratio * m * g_m
        Dr = relative_density(Fc)
        X_limit = critical_rebound_limit(Dr)
        diff = X_limit - bounce
        if np.any(diff < 0):
            idx = np.where(diff < 0)[0][0]
            ax10.scatter(Fc_ratio[idx], bounce[idx], color=limit_colors[m], edgecolor=c, linewidth=2, s=70, zorder=5)

ax10.set_xlabel('Force-to-Weight Ratio ($F_c / W_{moon}$)')
ax10.set_ylabel('Displacement (mm)')
ax10.set_ylim(0, 0.7)
ax10.set_xlim(1, 15)
ax10.set_title('10. Micro-Bounce vs. Mass-Dependent Traction Limits')
ax10.grid(True, alpha=0.4)
ax10.legend(loc='upper right', fontsize=8, ncol=3)
fig10.tight_layout()

# Display all plots at once
plt.show()