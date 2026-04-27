import numpy as np
import pandas as pd
from itertools import combinations, product

# ============================================================
# 1. USER SETTINGS
# ============================================================

masses = [50, 75, 100]          # lunar rover masses [kg]
freqs = [30, 50, 75]            # Hz
betas = [0.50, 0.75, 0.95]      # fraction of hop boundary

Fc_const = 100.0                # N, common absolute centrifugal force

g_m = 1.62                      # lunar gravity [m/s^2]

# Hole design search
n_holes = 3                 # start with 3 holes
candidate_radii_mm = np.arange(75, 250, 5)  # possible hole radii from shaft [mm]

# Put all holes on same radial line by default so moments add constructively.
# If you want angularly separated holes, edit this later.
hole_angles_deg = [0, 0, 0]

# Disc / hole assumptions
disc_thickness_mm = 6.35        # 1/4 inch
hole_diameter_mm = 5.5          # M5 clearance hole
pla_density_g_cm3 = 1.24
infill_fraction = 0.30
effective_density_g_cm3 = pla_density_g_cm3 * infill_fraction

# Hardware masses
screw_masses_g = {
    8: 3.2,
    10: 3.7,
    12: 4.2,
    16: 5.2,
    20: 6.1,
    25: 7.3,
    30: 8.5,
    35: 9.7,
}

nut_mass_g = 1.1
washer_mass_g = 0.5

# Approximate washer limits by screw length.
max_washers_by_length = {
    8: 0,
    10: 1,
    12: 2,
    16: 4,
    20: 6,
    25: 6,
    30: 6,
    35: 6,
}

# nut_counts = [1, 2]   # allow one or two nuts per screw
nut_counts = [1]      

# ============================================================
# 2. YOUR EXISTING HOP MODEL FUNCTIONS
# ============================================================

A_roller = 0.5 * 0.05
ku_loose, ku_dense = 12e6, 60e6
D_min, D_max, F_ref = 0.35, 0.95, 200.0

def relative_density(Fc, f=70.0):
    eta = 1 / (1 + np.exp(-0.5 * (np.array(f) - 15)))
    F_eff = np.array(Fc) * eta
    return D_max - (D_max - D_min) * np.exp(-F_eff / F_ref)

def dynamic_response(Fc, f, m_r, Dr=0.65):
    w = 2 * np.pi * np.array(f)
    W_total = m_r * g_m

    ku = ku_loose + (ku_dense - ku_loose) * ((np.array(Dr) - 0.35) / (0.95 - 0.35))**2
    k_soil = ku * A_roller
    c_soil = 2 * 0.3 * np.sqrt(k_soil * m_r)

    den = np.sqrt((k_soil - m_r * w**2)**2 + (c_soil * w)**2)
    Y = np.array(Fc) / den
    Fd = Y * np.sqrt(k_soil**2 + (c_soil * w)**2)
    Rhop = Fd / W_total

    return Y * 1000, Fd, Rhop

def hop_metric(Fc, f, m_r):
    Dr = relative_density(Fc, f=f)
    return dynamic_response(Fc, f, m_r, Dr)[2]

def find_hop_boundary_ratio(m, f, ratio_min=0.5, ratio_max=15.0, n=5000):
    W = m * g_m
    ratio_search = np.linspace(ratio_min, ratio_max, n)
    Fc_array = ratio_search * W
    Rhop = hop_metric(Fc_array, f, m)

    if not np.any(Rhop >= 1.0):
        return np.nan

    idx = np.where(Rhop >= 1.0)[0][0]
    return ratio_search[idx]

# ============================================================
# 3. TARGET ECCENTRIC MOMENTS
# ============================================================

target_rows = []

for m in masses:
    W = m * g_m

    for f in freqs:
        R_bnd = find_hop_boundary_ratio(m, f)

        for beta in betas:
            R_test = beta * R_bnd
            Fc_test = R_test * W
            mer_req = Fc_test / (2 * np.pi * f)**2

            target_rows.append({
                "type": "beta_hop_boundary",
                "mass_kg": m,
                "freq_Hz": f,
                "beta": beta,
                "Fc_N": Fc_test,
                "target_mer_kg_m": mer_req,
            })

# Constant-force cases.
# Note: target eccentric moment depends only on frequency, not rover mass.
for f in freqs:
    mer_req = Fc_const / (2 * np.pi * f)**2
    target_rows.append({
        "type": "constant_Fc",
        "mass_kg": "all",
        "freq_Hz": f,
        "beta": np.nan,
        "Fc_N": Fc_const,
        "target_mer_kg_m": mer_req,
    })

targets_df = pd.DataFrame(target_rows)

# ============================================================
# 4. HARDWARE MASS OPTIONS
# ============================================================

def removed_pla_mass_g(
    hole_diameter_mm=hole_diameter_mm,
    disc_thickness_mm=disc_thickness_mm,
    density_g_cm3=effective_density_g_cm3
):
    d_cm = hole_diameter_mm / 10
    t_cm = disc_thickness_mm / 10
    volume_cm3 = np.pi * (d_cm / 2)**2 * t_cm
    return volume_cm3 * density_g_cm3

m_removed_g = removed_pla_mass_g()

def make_hardware_options():
    """
    Returns list of available hole states.
    Net mass means: hardware mass added minus PLA removed.
    Empty hole means: only removed PLA, so negative mass.
    """
    options = []

    # Empty drilled hole
    options.append({
        "label": "empty",
        "inserted_mass_g": 0.0,
        "net_mass_g": -m_removed_g,
        "screw_len_mm": None,
        "nuts": 0,
        "washers": 0,
    })

    for L, screw_g in screw_masses_g.items():
        max_w = max_washers_by_length[L]

        for nuts in nut_counts:
            for washers in range(max_w + 1):
                inserted = screw_g + nuts * nut_mass_g + washers * washer_mass_g
                net = inserted - m_removed_g

                options.append({
                    "label": f"M5x{L}, {nuts} nut, {washers} washer",
                    "inserted_mass_g": inserted,
                    "net_mass_g": net,
                    "screw_len_mm": L,
                    "nuts": nuts,
                    "washers": washers,
                })

    return options

options = make_hardware_options()
options_df = pd.DataFrame(options)

print(f"Removed PLA per M5 hole: {m_removed_g:.4f} g")
print(f"Number of available hole states: {len(options)}")

# ============================================================
# 5. MOMENT CALCULATION
# ============================================================

def config_moment_vector_kg_m(radii_mm, angles_deg, option_indices):
    """
    Computes vector eccentric moment from selected hardware in each hole.

    Eccentric moment vector:
        M = sum_i m_i r_i [cos(theta_i), sin(theta_i)]

    where m_i is net mass = hardware mass - removed PLA.
    """
    Mx = 0.0
    My = 0.0

    for r_mm, theta_deg, opt_idx in zip(radii_mm, angles_deg, option_indices):
        r_m = r_mm / 1000
        theta = np.deg2rad(theta_deg)
        m_kg = options[opt_idx]["net_mass_g"] / 1000

        Mx += m_kg * r_m * np.cos(theta)
        My += m_kg * r_m * np.sin(theta)

    return np.array([Mx, My])

def config_moment_mag_kg_m(radii_mm, angles_deg, option_indices):
    vec = config_moment_vector_kg_m(radii_mm, angles_deg, option_indices)
    return np.linalg.norm(vec)

# ============================================================
# 6. OPTIMIZE HOLE RADII
# ============================================================

target_mers = targets_df["target_mer_kg_m"].values

best = None
best_assignments = None

# radii combinations are sorted, e.g. [30, 65, 100] mm
for radii_tuple in combinations(candidate_radii_mm, n_holes):
    print(f"Evaluating: radii {radii_tuple} mm")
    radii = np.array(radii_tuple)
    angles = np.array(hole_angles_deg[:n_holes])

    all_combo_indices = list(product(range(len(options)), repeat=n_holes))

    moments = np.array([
        config_moment_mag_kg_m(radii, angles, combo)
        for combo in all_combo_indices
    ])

    assignments = []
    errors = []

    for target in target_mers:
        idx = np.argmin(np.abs(moments - target))
        achieved = moments[idx]
        rel_err = abs(achieved - target) / target

        assignments.append({
            "target": target,
            "achieved": achieved,
            "rel_error": rel_err,
            "combo": all_combo_indices[idx],
        })
        # print(assignments[-1])
        errors.append(rel_err)

    errors = np.array(errors)

    score = np.sqrt(np.mean(errors**2)) + 0.5 * np.max(errors)

    if best is None or score < best["score"]:
        best = {
            "score": score,
            "radii_mm": radii,
            "mean_error": np.mean(errors),
            "rms_error": np.sqrt(np.mean(errors**2)),
            "max_error": np.max(errors),
        }
        best_assignments = assignments

# ============================================================
# 7. REPORT RESULTS
# ============================================================

print("\nBest hole radii:")
print(best)

result_rows = []

for target_row, assign in zip(targets_df.to_dict("records"), best_assignments):
    combo = assign["combo"]

    hole_desc = []
    for h, opt_idx in enumerate(combo):
        opt = options[opt_idx]
        hole_desc.append(
            f"Hole {h+1} @ {best['radii_mm'][h]} mm: {opt['label']} "
            f"(net {opt['net_mass_g']:.2f} g)"
        )

    result_rows.append({
        **target_row,
        "achieved_mer_kg_m": assign["achieved"],
        "abs_error_kg_m": assign["achieved"] - assign["target"],
        "rel_error_percent": 100 * assign["rel_error"],
        "hardware_config": " | ".join(hole_desc),
    })

results_df = pd.DataFrame(result_rows)

pd.set_option("display.max_colwidth", 200)
print("\nBest hardware assignments:")
print(results_df[[
    "type",
    "mass_kg",
    "freq_Hz",
    "beta",
    "Fc_N",
    "target_mer_kg_m",
    "achieved_mer_kg_m",
    "rel_error_percent",
    "hardware_config"
]].round({
    "Fc_N": 2,
    "target_mer_kg_m": 6,
    "achieved_mer_kg_m": 6,
    "rel_error_percent": 2,
}).to_string(index=False))

results_df.to_csv("optimized_eccentric_hole_assignments.csv", index=False)
options_df.to_csv("hardware_mass_options.csv", index=False)