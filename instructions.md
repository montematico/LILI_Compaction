Here are the instructions to update your sweep script and the GeminiCLI prompt to generate your multi-rover sensitivity analysis.

### Part 1: Updating `ChenSweep.m` for Traction Margin

You will need to make two modifications to your main sweep script to include traction in the Pareto optimization.

**1. Update the Initial Setup (Section 1)**
Find the `pareto_weights` array at the end of Section 1 and add a 5th weight for the new metric. Since we want to prioritize traction margin, give it a significant weight.

```matlab
% Pareto Optimization Weights 
% [Total Mass, Energy, Mass Ratio, Freq Diff from 50Hz, Minimize Required Friction (Traction)]
pareto_weights = [1.0, 0.5, 0.0, 0.0, 1.0]; 
```

**2. Update the Pareto Optimization (Section 4)**
Find the beginning of Section 4 (`%% 4. Pareto Optimization`). Right before you define the `objectives` matrix, calculate the required coefficient of friction ($\mu_{req}$) and add it as the 5th objective.

```matlab
    % Extract valid vectors
    m_rov_vec = M_ROV_grid(valid_idx);
    energy_vec = Energy_results(valid_idx);
    mass_ratio_vec = M_RATIO_grid(valid_idx);
    f_vec = F_grid(valid_idx);
    
    % --- NEW: Calculate Required Traction Coefficient ---
    % mu_req = F_loco / W_total
    mu_req_vec = Max_Traction_results(valid_idx) ./ (m_rov_vec .* g_moon);
    
    % 5-Objective Matrix: 
    % [Minimize Mass, Minimize Energy, Maximize Ratio, Target 50Hz, Minimize mu_req]
    objectives = [m_rov_vec, energy_vec, 1./mass_ratio_vec, abs(f_vec - 50), mu_req_vec];
    
    % Apply Weights to Objectives (Optional, but useful if using a weighted sum approach later)
    objectives = objectives .* pareto_weights;
```

---

### Part 2: GeminiCLI Prompt for Multi-Rover Analysis

Copy and paste the following prompt into GeminiCLI to generate the `CompareRoverSensitivity.m` script.

```text
Please write a MATLAB script named `CompareRoverSensitivity.m` to evaluate and compare the traction robustness of multiple lunar compaction rover designs under degraded soil conditions.

**Context & Inputs:**
1. Load the workspace `SweepResults.mat` which contains the grid sweep variables (e.g., `M_ROV_grid`, `R_ROLLER_grid`) and the simulation constants.
2. Define an array of target indices at the top of the script (e.g., `target_indices = [12, 452, 1084]`) representing a light, medium, and heavy rover.
3. Define an array of labels for the legend (e.g., `labels = {'20 kg Rover', '35 kg Rover', '50 kg Rover'}`).
4. Define a threshold for traction failure: `mu_threshold = 0.45`.

**Core Task: Soil Degradation Loop**
Create a vector `soil_mod_shifts` that ranges from -0.4 to +0.1 in 10 steps (representing -40% weaker to +10% stronger soil).
For each rover index:
1. Extract `m_opt` and `R_roller_opt` from the grid. Assume `roller_fraction = 0.5` and `n_rollers = 2`. Derive `D_roller`, `b_roller = 0.3`, and `W_roller`.
2. Loop through the `soil_mod_shifts`. For each shift, modify the base moduli: `kc_test = SOIL.kc * (1 + shift)` and `kphi_test = SOIL.kphi * (1 + shift)`.
3. Calculate Pass 1 sinkage (`z_p1`) using the `calculate_tandem_sinkage` helper function (you will need to include this helper function at the bottom of the script).
4. Calculate the locomotion resistance forces (`R_r`, `R_c`, `R_b`) using the Bekker/Terzaghi equations just like in the main sweep's first pass. Sum these to get `F_loco_p1`.
5. Calculate the required coefficient of friction: `mu_req = F_loco_p1 / (m_opt * g_moon)`. Store this in a results matrix.

**Plotting: The Robustness "Money Plot"**
Create a single figure with the following specifications:
1. Plot the `mu_req` results for each rover against the `soil_mod_shifts * 100` (percentage). Use distinct colors/markers and include a legend using the predefined labels.
2. Add a horizontal dashed red line at `y = mu_threshold` labeled "Traction Failure Limit (\mu = 0.45)".
3. IMPORTANT: Set the X-axis direction to 'reverse' (`set(gca, 'XDir', 'reverse')`) so the plot reads left-to-right as the soil gets softer/weaker.
4. X-axis label: "Soil Strength Degradation (%)". Y-axis label: "Required Traction Coefficient (\mu_{req})". Title: "Rover Robustness: Traction Risk in Degraded Lunar Soil".
5. Optionally, add a shaded red patch covering the area above the `mu_threshold` to visually indicate the "Danger Zone".

Ensure the code is well-commented and handles variable extraction safely.
```