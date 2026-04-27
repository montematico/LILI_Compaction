import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 1. USER INPUTS
# ============================================================

# DOE you selected
masses = [50, 75, 100]      # lunar rover masses [kg]
freqs = [30, 50, 75]        # Hz
betas = [0.50, 0.75, 0.95]  # fractions of hop boundary

# Lunar gravity
g_m = 1.62

# Adjustable constant-force design curve
Fc_const = 100.0   # N  <-- change this
f_const  = 50.0    # Hz <-- change this

# Plot bounds for point-mass vs radius
r_min, r_max = 0.005, 0.2   # m
m_min, m_max = 0.01, 0.20    # kg

# ============================================================
# 2. MODEL CONSTANTS (matching your current dynamic model)
# ============================================================

A_roller = 0.5 * 0.05
ku_loose, ku_dense = 12e6, 60e6
D_min, D_max, F_ref = 0.35, 0.95, 200.0

# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

def relative_density(Fc, f=70.0):
    """
    Heuristic compaction law used in your current model.
    """
    eta = 1.0 / (1.0 + np.exp(-0.5 * (np.array(f) - 15.0)))
    F_eff = np.array(Fc) * eta
    return D_max - (D_max - D_min) * np.exp(-F_eff / F_ref)

def dynamic_response(Fc, f, m_r, Dr=0.65):
    """
    1-DOF rigid-body vertical oscillator.
    Returns:
      Y_mm = displacement amplitude [mm]
      Fd   = transmitted dynamic force amplitude [N]
      Rhop = hop severity ratio = Fd / W_total
    """
    w = 2.0 * np.pi * np.array(f)
    W_total = m_r * g_m

    ku = ku_loose + (ku_dense - ku_loose) * ((np.array(Dr) - 0.35) / (0.95 - 0.35))**2
    k_soil = ku * A_roller
    c_soil = 2.0 * 0.3 * np.sqrt(k_soil * m_r)

    den = np.sqrt((k_soil - m_r * w**2)**2 + (c_soil * w)**2)
    Y = np.array(Fc) / den
    Fd = Y * np.sqrt(k_soil**2 + (c_soil * w)**2)
    Rhop = Fd / W_total
    return Y * 1000.0, Fd, Rhop

def hop_metric(Fc, f, m_r):
    """
    Global rigid-body hop severity ratio.
    Boundary occurs at Rhop = 1.
    """
    Dr = relative_density(Fc, f=f)
    return dynamic_response(Fc, f, m_r, Dr)[2]

# ============================================================
# 4. SOLVE THE 27 DOE POINTS
# ============================================================

ratio_search = np.linspace(0.5, 15.0, 5000)
rows = []

for m in masses:
    W_moon = m * g_m

    for f in freqs:
        Fc_array = ratio_search * W_moon
        Rhop_array = hop_metric(Fc_array, f, m)

        idx = np.where(Rhop_array >= 1.0)[0][0]
        R_bnd = ratio_search[idx]

        for beta in betas:
            R_test = beta * R_bnd
            Fc_test = R_test * W_moon
            mer_req = Fc_test / (2 * np.pi * f)**2

            rows.append({
                "mass_kg_lunar": m,
                "freq_Hz": f,
                "beta": beta,
                "hop_boundary_ratio": R_bnd,
                "test_ratio": R_test,
                "centrifugal_force_N": Fc_test,
                "eccentric_moment_kg_m": mer_req
            })

df = pd.DataFrame(rows).sort_values(["mass_kg_lunar", "freq_Hz", "beta"]).reset_index(drop=True)

# Optional: save CSV
# df.to_csv("hop_boundary_doe_points.csv", index=False)

# ============================================================
# 5. HEATMAP GRID: eccentric moment = m_e * r
# ============================================================

r = np.linspace(r_min, r_max, 500)
m_e = np.linspace(m_min, m_max, 500)
R, M = np.meshgrid(r, m_e)
MER = M * R

levels = np.linspace(MER.min(), MER.max(), 60)

# ============================================================
# 6. PLOT 1: BASE HEATMAP
# ============================================================

fig1, ax1 = plt.subplots(figsize=(10, 7))
im1 = ax1.contourf(R, M, MER, levels=levels, cmap="viridis")
cbar1 = fig1.colorbar(im1, ax=ax1, label="Eccentric Moment, $m_e r$ (kg·m)")

ax1.set_xlabel("Radius from Shaft, $r$ (m)")
ax1.set_ylabel("Point Mass, $m_e$ (kg)")
ax1.set_title("Point Mass vs. Radius\nHeatmap = Eccentric Moment $m_e r$")
ax1.grid(True, alpha=0.25)

fig1.tight_layout()
plt.show()

# ============================================================
# 7. PLOT 2: HEATMAP + 27 DOE MOMENT CURVES
# ============================================================

fig2, ax2 = plt.subplots(figsize=(11, 8))
im2 = ax2.contourf(R, M, MER, levels=levels, cmap="viridis")
cbar2 = fig2.colorbar(im2, ax=ax2, label="Eccentric Moment, $m_e r$ (kg·m)")

colors = plt.cm.tab20(np.linspace(0, 1, len(df)))

for i, row in df.iterrows():
    mer_req = row["eccentric_moment_kg_m"]
    mass_curve = mer_req / r
    mask = (mass_curve >= m_min) & (mass_curve <= m_max)

    if np.any(mask):
        label = f'{int(row["mass_kg_lunar"])}kg, {int(row["freq_Hz"])}Hz, β={row["beta"]:.2f}'
        ax2.plot(r[mask], mass_curve[mask], lw=1.6, color=colors[i], label=label)

ax2.set_xlabel("Radius from Shaft, $r$ (m)")
ax2.set_ylabel("Point Mass, $m_e$ (kg)")
ax2.set_title("Required Eccentric Moment Curves for the 27 Design Points")
ax2.grid(True, alpha=0.25)
# ax2.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=True)

fig2.tight_layout()
plt.show()

# ============================================================
# 8. PLOT 3: HEATMAP + DOE CURVES + CONSTANT FORCE CURVE
# ============================================================

fig3, ax3 = plt.subplots(figsize=(11, 8))
im3 = ax3.contourf(R, M, MER, levels=levels, cmap="viridis")
cbar3 = fig3.colorbar(im3, ax=ax3, label="Eccentric Moment, $m_e r$ (kg·m)")

# DOE curves
for i, row in df.iterrows():
    mer_req = row["eccentric_moment_kg_m"]
    mass_curve = mer_req / r
    mask = (mass_curve >= m_min) & (mass_curve <= m_max)

    if np.any(mask):
        label = f'{int(row["mass_kg_lunar"])}kg, {int(row["freq_Hz"])}Hz, β={row["beta"]:.2f}'
        ax3.plot(r[mask], mass_curve[mask], lw=1.3, color=colors[i], alpha=0.8)

# Constant centrifugal force curve
mer_const = Fc_const / (2 * np.pi * f_const)**2
mass_curve_const = mer_const / r
mask_const = (mass_curve_const >= m_min) & (mass_curve_const <= m_max)

ax3.plot(
    r[mask_const],
    mass_curve_const[mask_const],
    color="white",
    linestyle="--",
    linewidth=3,
    label=f'Constant $F_c$ = {Fc_const:.0f} N @ {f_const:.0f} Hz'
)

ax3.set_xlabel("Radius from Shaft, $r$ (m)")
ax3.set_ylabel("Point Mass, $m_e$ (kg)")
ax3.set_title("DOE Eccentric Moment Curves + Constant Centrifugal Force Design Curve")
ax3.grid(True, alpha=0.25)
ax3.legend(loc="upper right", fontsize=9, frameon=True)

fig3.tight_layout()
plt.show()

# ============================================================
# 9. PRINT TABLE OF DOE POINTS
# ============================================================

print(df.round({
    "hop_boundary_ratio": 3,
    "test_ratio": 3,
    "centrifugal_force_N": 2,
    "eccentric_moment_kg_m": 6
}).to_string(index=False))