# Dynamic Analysis: Micro-Bounce and Traction Limits

This sub-project focuses on the high-fidelity dynamic modeling of the lunar compaction rover's interaction with regolith. While the main simulation (`ChenSweep.m`) performs a broad kinematic screening, this analysis uses a **1-Degree-of-Freedom (1-DOF) Spring-Mass-Damper model** to establish the physical boundaries of stable operation.

## Core Problem: The Traction Paradox
Vibratory compaction requires high centrifugal forces ($F_c$) to crush soil particles. However, in the low-gravity environment of the Moon ($1.62\text{ m/s}^2$), high vibration forces can easily cause the rover to "hop" or lose contact with the ground. If the wheels lift out of their compacted ruts, forward traction is lost, and the mission fails.

## Physical and Mathematical Models

### 1. 1-DOF Dynamic Contact Model
The rover is modeled as a forced harmonic oscillator interacting with the soil matrix:
$$m_r \ddot{y} + c\dot{y} + ky = -W_{moon} + F_c \sin(\omega t)$$

The steady-state displacement amplitude ($Y$) is determined by the system's impedance:
$$Y = \frac{F_c}{\sqrt{(k - m_r \omega^2)^2 + (c\omega)^2}}$$

**Contact Regimes:**
- **Regime A (Continuous Contact):** $F_d \le W_{moon}$. The roller stays on the ground.
- **Regime B/C (Partial Unloading / Lift-off):** $F_d > W_{moon}$. The normal force reaches zero, causing "micro-bounce." The code heuristically approximates this lift-off amplitude by scaling the excess dynamic force.

### 2. Frequency-Dependent Compaction ($D_r$)
Based on Chen et al. (2021), compaction efficiency is not just about force; it is sensitive to frequency. We apply a logistic frequency modifier $\eta(f)$ (centered at 15 Hz) to calculate an *effective* compaction force:
$$F_{eff} = F_c \cdot \eta(f)$$
$$D_r(F_c, f) = D_{max} - (D_{max} - D_{min}) e^{-F_{eff} / F_{ref}}$$

### 3. Elastic Rebound & Traction Limit
Forward traction requires continuous normal force. The allowable "micro-bounce" ($X_{limit}$) before wheels lose traction is defined by the elastic rebound of the soil under the wheel's contact pressure ($P_{wheel}$):
$$X_{limit} = \frac{P_{wheel}}{k_u(D_r)}$$
Where $k_u$ is the unloading modulus, which scales quadratically with soil density ($D_r$).

### 4. Theoretical High-Speed Asymptote
At very high frequencies ($\omega \to \infty$), mass-inertia dominates the system. The bounce amplitude is bounded by the ratio of eccentric moment to rover mass:
$$X_{max} \approx \frac{m_e r}{m_r}$$
This allows engineers to size the eccentric mass to ensure the bounce never exceeds the traction limit, regardless of speed.

---

## File Overview

- **`bounce2.ipynb`**: The primary research notebook. It contains the 1-DOF derivations, soil stiffness interpolations, and generated "Operating Maps."
- **`microbounce2.py`**: A standalone Python implementation of the dynamic models. It features an **Interactive Operating Map** using a mass slider to visualize how the safe vibration envelope changes with rover weight.
- **`bounce_og.ipynb` / `microbounce.py`**: Legacy versions and initial kinematic screening models.

## Critical Insights

- **Mass vs. Frequency**: Heavier rovers are more stable (lower bounce) but also "lock" the soil into a rigid state faster, which reduces the allowable bounce limit.
- **The Sweet Spot**: Stable compaction typically occurs at Force-to-Weight ratios between **1.5 and 4.0**, provided the frequency is high enough ($>30\text{ Hz}$) to minimize displacement.
- **Traction Failure**: The "Danger Zone" occurs when the soil becomes too dense/rigid; the reduced elastic rebound means even tiny bounces can cause total traction loss.

## Limitations
- **Damping Uncertainty**: The model assumes a constant damping ratio ($\zeta = 0.3$). In reality, soil damping is highly non-linear and changes with compaction.
- **Unilateral Approximation**: The current model uses a heuristic to approximate lift-off behavior rather than a full non-smooth dynamic solver (like the Moreau-Jean scheme).
- **Decoupled Assumption**: The model separates roller and wheel contact. In reality, the vibration from the roller travels through the rover frame and affects wheel contact dynamically.
