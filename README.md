# Lunar Regolith Compaction Rover Analysis

This repository contains a suite of MATLAB scripts and functions designed to analyze, simulate, and optimize the design of a lunar rover equipped with vibratory compactors. The primary goal is to determine the optimal rover configurations for compressing lunar regolith to build infrastructure, such as landing pads or roads, on the Moon.

## Table of Contents
- [Overview and Purpose](#overview-and-purpose)
- [File Architecture](#file-architecture)
- [Physics and Mathematical Models](#physics-and-mathematical-models)
  - [1. Vibratory Compaction Model](#1-vibratory-compaction-model)
  - [2. Sinkage and Locomotion Resistance Model](#2-sinkage-and-locomotion-resistance-model)
  - [3. Drive Wheel Traction Model](#3-drive-wheel-traction-model)
  - [4. Power and Energy Model](#4-power-and-energy-model)
- [Critical Analysis](#critical-analysis)

---

## Overview and Purpose

The primary purpose of this codebase is to perform a comprehensive systems engineering analysis and trade study for a lunar compaction rover. 

To create stable foundations on the Moon, the loose, powdery lunar regolith must be compacted. The code models a rover driving back and forth over a designated area (a "pad"), utilizing heavy vibratory rollers to increase the soil's relative density. 

The analysis aims to find the "sweet spot" in a massive 7-dimensional design space (rover mass, roller size, vibration frequency, eccentric mass, etc.) that balances competing objectives:
1. **Minimizing Mass**: Launching mass to the Moon is extremely expensive.
2. **Minimizing Energy/Power**: Solar panels and batteries add mass; the rover must be energy-efficient.
3. **Minimizing Time**: The rover must compact the area within a reasonable timeframe (e.g., 6 months).
4. **Maximizing Traction Robustness**: The rover must not get stuck in the loose regolith.

---

## File Architecture

The codebase is organized into several interconnected scripts:

- **`ChenSweep.m`**: The core simulation engine. It performs a massive 7D grid search over design parameters. It uses parallel computing to evaluate thousands of rover configurations, filters out invalid designs (e.g., those that take too long, use too much energy, or get stuck), and performs Pareto optimization to identify the optimal frontier of designs.
- **`CompactionChen.m`**: A simplified, single-configuration test script. It runs the simulation for one specific set of parameters and outputs a clean productivity table. Great for debugging the core physics loop.
- **`lunar_wheel_model.m`**: A standalone terramechanics model that evaluates the tractive capabilities of the rover's drive wheels. It calculates drawbar pull based on Bekker pressure-sinkage and Janosi-Hanamoto shear laws.
- **`IdxSelect.m`**: A helper script that parses the massive results dataset from `ChenSweep.m` and selects specific "archetype" designs of interest, such as the "Minimum Mass Survivor," the "Balanced Design," and the "Maximum Traction Margin" tank.
- **`CompareRoverSensitivity.m`**: Evaluates the selected optimal rovers under degraded, highly uncertain soil conditions. It alters the soil moduli ($k_c, k_\phi$) by $\pm X\%$ to see which designs are most likely to fail (get stuck) if the regolith is looser than expected.
- **`SingleMassAnalysis.m`**: A deep-dive script that takes a single optimal design and runs a high-fidelity, pass-by-pass analysis. It tracks how density, sinkage, power draw, and required traction evolve over time as the soil stiffens.
- **`plot_compaction_results.m`**: Generates high-quality visualizations of the Pareto fronts, success rates, and trade spaces from the sweep results.

---

## Physics and Mathematical Models

The code fuses two major physical domains: **Soil Compaction Physics** and **Vehicle Terramechanics**.

### 1. Vibratory Compaction Model

The compaction physics are based on empirical models developed by Chen et al., mapped to lunar regolith parameters. The vibrator uses an eccentric mass $m_e$ rotating at frequency $f$ (angular velocity $\omega = 2\pi f$). 

The excitation force amplitude generated is:
$$F_0 = m_e \omega^2$$

During each vibration cycle, if the dynamic force exceeds the plastic limit of the soil $P_b$, a residual plastic deformation $x_{cr}$ (sinkage/compaction) occurs. The plastic limit is defined empirically as:
$$P_b = 0.07932(100\rho_e - 184)^2 \left( \frac{A_{col}}{A_{chen}} \right)$$
*(Note: Lunar densities are linearly mapped to Earth equivalent densities $\rho_e$ to utilize these empirical formulas).*

The residual deformation per cycle is calculated as:
$$x_{cr} = \frac{m_d \omega^2}{2 k_s^2 P_b} \left( \frac{F_0^2 m_d^2}{m_t^2} - P_b^2 \right) + \left( \frac{P_b}{k_s} - \frac{P_b}{k_{su}} \right)$$
Where:
- $m_d$: Mass of the bouncing steel drum.
- $m_t$: Total mass acting on the wheel.
- $k_s$: Compressive stiffness of the soil.
- $k_{su}$: Resilient stiffness of the soil.

As the layer compresses by $x_{cr}$, the density of the layer increases according to conservation of mass:
$$\rho_{new} = \rho_{old} \frac{h_{layer}}{h_{layer} - x_{cr}}$$

### 2. Sinkage and Locomotion Resistance Model

To move, the rover must overcome the resistance of the soil. As multiple rollers pass over the soil, cumulative sinkage occurs. 

**Tandem Sinkage:**
The model uses a recursive Bekker-based formulation for $N$ tandem passes. The sinkage $z_i$ after the $i$-th pass is:
$$z_i = \left( \frac{3W}{2(k_c + b k_\phi)\sqrt{D}} + z_{i-1}^{1.5} \right)^{2/3}$$

**Locomotion Resistances:**
The total force required to move the rover forward ($F_{loco}$) is the sum of three components:
$$F_{loco} = R_r + R_c + R_b$$

1. **Rolling Resistance ($R_r$)**: Standard friction in the bearings/mechanisms.
   $$R_r = W_{roller} \cdot n_{rollers} \cdot c_f$$
2. **Compaction Resistance ($R_c$)**: The work done to compress the soil vertically.
   $$R_c = \frac{1}{2} (k_c + b k_\phi) z^2 n_{rollers}$$
3. **Bulldozing Resistance ($R_b$)**: The force required to push the mound of soil building up in front of the leading roller. This involves complex Terzaghi bearing capacity parameters ($K_c, K_\gamma$, cohesion $c$, friction angle $\phi$):
   $$R_b = \frac{b \sin(\alpha + \phi)}{2 \sin\alpha \cos\phi} \left( 2z c K_c + \gamma z^2 K_\gamma \right) + \frac{l_o^3 \gamma}{3} \left( \frac{\pi}{2} - \phi \right) + c l_o^2 \left( 1 + \tan\left(\frac{\pi}{4} + \frac{\phi}{2}\right) \right)$$

### 3. Drive Wheel Traction Model

The `lunar_wheel_model.m` calculates if the drive wheels can generate enough Gross Traction ($H$) to overcome the resistances. It utilizes the Janosi-Hanamoto shear stress equation:
$$\tau = (c + p \tan\phi) (1 - e^{-j/K})$$
Where:
- $p = W/A$ is the average normal pressure.
- $j = slip \cdot l$ is the shear displacement.
- $K$ is the shear deformation modulus.

The Gross Traction is $H = \tau A$. The **Drawbar Pull (DP)**—the actual usable force for moving the vehicle—is:
$$DP = H - (R_c + R_b + R_g)$$
If the required locomotion force $F_{loco}$ exceeds the available Drawbar Pull, the rover gets stuck, and the design is marked invalid.

### 4. Power and Energy Model

The total power requirement dictates the battery and solar panel mass. It is the sum of:
1. **Locomotion Power**: $P_{loco} = F_{loco} \cdot v$
2. **Vibration Power**: The kinetic energy imparted to the drum per cycle, scaled by frequency and drivetrain efficiency.
   $$P_{vib} = \frac{1}{\eta_{mech}} \left( \frac{(m_e \omega)^2}{2 m_d} \right) f$$
3. **Hotel Power**: Static power draw from computers, sensors, and thermal systems ($P_{hotel}$).

---

## Critical Analysis

### Strengths and Purposes
1. **Comprehensive Trade Space**: The code brilliantly explores the complex non-linear trade-offs in rover design. For example, a heavier rover compacts faster, but requires much more traction, risking getting stuck.
2. **Dynamic Soil Response**: The simulation correctly models soil as a dynamic material. As density increases pass-by-pass, the soil moduli ($k_c, k_\phi$) are interpolated to represent "stiffening," which realistically reduces sinkage and resistance on subsequent passes.
3. **Risk Mitigation (Traction)**: The emphasis on extracting the "Required Friction Coefficient ($\mu_{req}$)" and the inclusion of `CompareRoverSensitivity.m` shows a strong engineering focus on the greatest risk to lunar rovers: getting permanently stuck in loose soil.

### Limitations and What it Doesn't Consider
While the systems-level abstraction is excellent for early-stage design, it has several limitations:
1. **Empirical Earth-to-Moon Mapping**: The compaction model relies on mapping lunar density parameters to Earth parameters to use Chen's equations. The interaction of vibration, low gravity, and hard vacuum on regolith particle cohesion is highly uncertain and may diverge from Earth-based empirical models.
2. **Flat Terrain Assumption**: The code assumes a perfectly flat compaction pad ($0^\circ$ grade). It does not account for craters, boulders, or slopes which drastically alter bulldozing resistance and required traction.
3. **Thermal Dynamics**: Vibrating a heavy drum in a vacuum generates immense waste heat. The model accounts for the energy cost but completely ignores the thermal mass (radiators) required to dissipate this heat.
4. **Turning and Kinematics**: The locomotion model simulates a straight-line pass. Turning a multi-roller vibratory rover induces massive shear forces on the regolith (often loosening it) and requires significantly more torque.
5. **Wear and Degradation**: The model assumes constant mechanical efficiency ($\eta_{mech}$) and does not account for the rapid wear caused by abrasive lunar dust getting into vibratory bearings.