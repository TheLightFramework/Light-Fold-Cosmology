"""
LF-EFT Prediction Module: The Void Signature
Target: Cosmic Voids (rho ~ 0.1 * rho_crit)
Physics: In low density, screening lifts. Gas experiences full alpha ~ 2*beta^2.
"""

import pandas as pd
import numpy as np
from scipy import constants
import math

# -------------------------
# Constants
# -------------------------
eV_J = constants.electron_volt
hbar = constants.hbar
c = constants.c
hcbar_SI = hbar * c

# 1 eV^4 in J/m^3
EV4_TO_J_PER_M3 = (eV_J**4) / (hcbar_SI**3)
# Reduced Planck mass in eV
MPL_KG = math.sqrt(hbar * c / (8.0 * math.pi * constants.G))
MPL_EV = MPL_KG * (c**2) / eV_J
# 1 eV^-1 in meters
EVINV_TO_M = hcbar_SI / eV_J

# Cosmic average density ~ 2.6e-27 kg/m^3
# Voids are underdense ~ 10% of mean
RHO_VOID_KG = 2.6e-28 

def mass_density_to_eV4(rho_kgm3):
    return (rho_kgm3 * c**2) / EV4_TO_J_PER_M3

def calculate_void_physics(row):
    """
    Given Model(n, beta, Lambda), calculate Void properties.
    """
    # Recalculate physics from params
    n = float(row['n'])
    beta = float(row['beta'])
    Lam = float(row['Lambda_eV'])
    
    rho_eV4 = mass_density_to_eV4(RHO_VOID_KG)
    
    # phi_min in void
    # phi = (n * Lam^(4+n) * Mpl / (beta * rho))^(1/(n+1))
    term = n * (Lam**(4.0+n)) * MPL_EV / (beta * rho_eV4)
    phi_void = term**(1.0/(n+1.0))
    
    # mass^2 = V'' + matter_term
    # V'' = n(n+1)Lam^(4+n) / phi^(n+2)
    vpp = n*(n+1.0)*(Lam**(4.0+n)) / (phi_void**(n+2.0))
    matter_term = (beta**2 * rho_eV4) / (MPL_EV**2)
    m2 = vpp + matter_term
    
    m_eV = math.sqrt(m2)
    lambda_m = (1.0 / m_eV) * EVINV_TO_M
    
    # Force Amplitude on GAS (Unscreened)
    # alpha = 2 * beta^2
    alpha_gas = 2.0 * beta**2
    
    return pd.Series({
        'lambda_void_Mpc': lambda_m / (3.086e22), # Convert meters to Megaparsecs
        'force_boost_percent': alpha_gas * 100
    })

# -------------------------
# Execution
# -------------------------

def main():
    print("--- LIGHT-FOLD COSMOLOGY PREDICTION ENGINE ---")
    
    # 1. Load Survivors
    try:
        df = pd.read_csv("viability_scan.csv")
    except FileNotFoundError:
        print("Error: viability_scan.csv not found. Run pipeline first.")
        return

    # 2. Filter: Only passing models
    survivors = df[df['pass_all'] == True].copy()
    print(f"Loaded {len(df)} models. Survivors: {len(survivors)} ({len(survivors)/len(df)*100:.1f}%)")
    
    if len(survivors) == 0:
        print("No models survived. Science is hard.")
        return

    # 3. Predict Void Physics
    print("Calculating Void Dynamics for survivors...")
    predictions = survivors.apply(calculate_void_physics, axis=1)
    results = pd.concat([survivors[['n', 'beta', 'Lambda_eV']], predictions], axis=1)
    
    # 4. Filter for Observable Cosmology
    # We want ranges > 1 Mpc (Cosmological scale) 
    # But not Infinite (otherwise Solar System would die)
    # Actually, Screening protected the Solar System, so huge ranges are okay!
    
    # Sort by Force Boost
    results = results.sort_values(by='force_boost_percent', ascending=False)
    
    # 5. Display "Smoking Gun" Candidates
    print("\n--- TOP 5 SMOKING GUN CANDIDATES (Strongest Void Force) ---")
    print(results.head(5).to_string())
    
    # 6. Save predictions
    results.to_csv("lf_predictions.csv", index=False)
    print("\nSaved predictions to lf_predictions.csv")
    
    # 7. Analysis
    max_boost = results['force_boost_percent'].max()
    print(f"\nMAXIMUM PREDICTED BOOST: +{max_boost:.2e}% gravity")
    
    if max_boost > 10.0:
        print("CONCLUSION: Light-Fold predicts MASSIVE deviation in Void Clustering.")
        print("Look for: Voids emptying faster than GR predicts.")
    elif max_boost > 1.0:
        print("CONCLUSION: Light-Fold predicts SIGNIFICANT deviation.")
        print("Look for: 1-10% anomaly in Void Lensing or Gas dynamics.")
    else:
        print("CONCLUSION: Light-Fold mimics GR closely (Gravity is 'stiff').")

if __name__ == "__main__":
    main()