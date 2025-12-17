# Light-Fold Cosmology: LF-EFT
### Effective Field Theory for Information-Driven Gravity (Negentropy)

**A Python-based framework for testing screened scalar-tensor theories against Solar System and Laboratory constraints.**

---

## 🌌 Abstract

**Light-Fold Cosmology** explores the hypothesis that the thermodynamic tendency of the universe to organize (Negentropy/Information) manifests physically as a scalar field coupled to matter. 

This repository implements a **Chameleon-Screened Effective Field Theory (EFT)** to model this interaction. It solves the "Fifth Force" problem by demonstrating that such a field can be:
1.  **Screened** in high-density environments (Solar System/Earth), passing all General Relativity tests (Cassini, Eöt-Wash).
2.  **Unscreened** in low-density cosmic voids, acting as a massive gravitational booster ($+2 \times 10^{10}\%$ force).

This mechanism offers a candidate explanation for the **Early Galaxy Formation Mystery** (JWST observations) and aligns with the **Second Law of Infodynamics** (Vopson).

## 📂 Repository Contents

### 1. The Physics Core
*   `lf_viability_pipeline.py`: The main scanner. It maps the parameter space ($n, \beta, \Lambda$) against physical constraints:
    *   **Gate A (Cassini):** Solar System PPN bounds ($\gamma - 1$).
    *   **Gate B (Eöt-Wash):** Laboratory torsion balance exclusion curves (Short-range Yukawa).
    *   **Gate D (Ephemerides):** Planetary orbital perturbations.
    
### 2. The Prediction Engine
*   `lf_void_prediction.py`: Takes the surviving models and transports them to a cosmic void ($\rho \approx 0.1 \bar{\rho}_{crit}$). It calculates the effective range ($\lambda_{void}$) and the force boost ($\alpha$) on diffuse gas.
*   `lf_filter_cosmic.py`: Filters for long-range "Architect" candidates capable of influencing Cosmic Web filament formation.

### 3. Data
*   `constraints_*.csv`: Digitized exclusion curves from physics literature (Kapner et al. 2007, Adelberger et al. 2009).
*   `viability_scan.csv`: The output of the viability pipeline (Validating the existence of a survival parameter space).

## 🛠️ Usage

**Prerequisites:** Python 3.9+, `numpy`, `pandas`, `scipy`.

1. **Run the Viability Scan** (Stress tests the theory):
   ```bash
   python lf_viability_pipeline.py
   ```
   *Output: Filters models that survive Eöt-Wash and Cassini.*

2. **Run the Prediction** (Simulate Void Physics):
   ```bash
   python lf_void_prediction.py
   ```
   *Output: Calculates gravity boost in the early universe/voids.*

## 🔬 Key Results

Our analysis identifies a specific parameter region ($n=1$, Strong Coupling $\beta \gg 1$, Small Scale $\Lambda \ll 1$ eV) that:
1.  **Passes** Laboratory tests via the Thin-Shell Mechanism (Screening).
2.  **Predicts** a short-range, high-intensity attractive force in cosmic voids ($+20,000,000,000\%$ gravity boost at $\sim 200$ AU).
3.  **Predicts** a long-range, gentle architectural force in supercluster filaments ($+0.02\%$ at $45$ Mpc).

## 📜 License
MIT License.
Open Science for the acceleration of Truth.

---
*Created by Jean Charbonneau & The Brotherhood (00003).*
*Inspired by the Second Law of Infodynamics.*
