import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

# choose your DOE
masses = [50, 75, 100]              # lunar rover masses [kg]
freqs = [30, 50, 75]            # Hz
betas = [0.50, 0.75, 0.95]          # fractions of hop boundary

ratio_search = np.linspace(1.0, 15.0, 2000)
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

rows = []

for m in masses:
    W_moon = m * g_m

    for f in freqs:
        Fc_array = ratio_search * W_moon
        Rhop_array = hop_metric(Fc_array, f, m)   # from your code

        if np.any(Rhop_array >= 1.0):
            idx = np.where(Rhop_array >= 1.0)[0][0]
            R_bnd = ratio_search[idx]
        else:
            R_bnd = np.nan

        for beta in betas:
            R_test = beta * R_bnd
            Fc_test = R_test * W_moon
            earth_equiv_mass = Fc_test / 9.81

            rows.append({
                "mass_kg_lunar": m,
                "freq_Hz": f,
                "boundary_ratio": R_bnd,
                "beta": beta,
                "test_ratio": R_test,
                "centrifugal_force_N": Fc_test,
                "earth_equiv_deadweight_kg": earth_equiv_mass
            })

df = pd.DataFrame(rows)
print(df.round(3))

freqs_2 = np.linspace(30, 75, 200)

colors = {50: 'tab:blue', 75: 'tab:orange', 100: 'tab:green'}
markers = {0.50: 'o', 0.75: 's', 0.95: '^'}

Fc_const = 100.0   # N  <-- adjustable constant centrifugal force

plt.figure(figsize=(10, 7))

for m in masses:
    W = m * g_m
    ratio_search = np.linspace(1.0, 15.0, 1500)
    boundary = []

    for f in freqs_2:
        Fc_array = ratio_search * W
        Rhop_array = hop_metric(Fc_array, f, m)
        if np.any(Rhop_array >= 1.0):
            idx = np.where(Rhop_array >= 1.0)[0][0]
            boundary.append(ratio_search[idx])
        else:
            boundary.append(np.nan)

    boundary = np.array(boundary)
    plt.plot(freqs_2, boundary, color=colors[m], lw=2.5, label=f'{m} kg boundary')

    # --- ADD CONSTANT Fc LINE FOR THIS MASS ---
    R_const = Fc_const / W
    plt.axhline(
        y=R_const,
        color=colors[m],
        linestyle='--',
        linewidth=1.8,
        alpha=0.8,
        label=f'{m} kg @ $F_c$={Fc_const:.0f} N'
    )

    # chosen experimental frequencies
    for f0 in freqs:
        idxf = np.argmin(np.abs(freqs_2 - f0))
        Rb = boundary[idxf]
        for beta in betas:
            plt.scatter(f0, beta * Rb, color=colors[m], marker=markers[beta], s=80)

plt.xlabel('Frequency (Hz)')
plt.ylabel('Force-to-Weight Ratio ($F_c/W_{r}$)')
plt.title(f'Proposed Experimental Points Relative to Bounce Boundary')
plt.grid(True, alpha=0.3)

# Create legend for mass boundaries + constant-Fc lines
# leg1 = plt.legend(loc='upper left', fontsize=9)
# plt.gca().add_artist(leg1)

# Create legend for beta markers
from matplotlib.lines import Line2D
marker_legend_elements = [
    Line2D([0], [0], marker=markers[beta], color='w',
           markerfacecolor='black', markersize=8,
           label=f'β = {beta}')
    for beta in sorted(betas)
]
plt.legend(handles=marker_legend_elements, loc='lower right',
           fontsize=10, title='Beta (Fraction of Boundary)')

plt.show()