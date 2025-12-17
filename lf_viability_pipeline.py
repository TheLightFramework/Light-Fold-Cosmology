"""
LF‑EFT viability pipeline (chameleon / screened scalar–tensor) -- AGGRESSIVE SCAN MODE

- Natural units: ħ=c=1 inside the EFT formulas.
- SI → eV^4 conversion is done explicitly.
- Gates:
  (A) Cassini-style PPN gate on Sun effective coupling (always on)
  (B) Optional Yukawa exclusion curves from CSV: alpha_max(lambda)

Outputs:
  viability_scan.csv

This is a research scaffold: it is intended to be *edited* to plug in better environment models,
more realistic density profiles, and more precise observable mappings.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Dict, Tuple, List

import math
import numpy as np
import pandas as pd
from scipy import constants


# -----------------------------
# Unit conversions (SI <-> natural units with eV)
# -----------------------------

eV_J = constants.electron_volt
hbar = constants.hbar
c = constants.c
hcbar_SI = hbar * c  # J*m

# 1 eV^4 in J/m^3 = (eV_J^4)/(ħc)^3
EV4_TO_J_PER_M3 = (eV_J**4) / (hcbar_SI**3)

# reduced Planck mass (in kg): sqrt(ħc / (8πG))
MPL_KG = math.sqrt(hbar * c / (8.0 * math.pi * constants.G))
# convert to eV (energy): M c^2 in Joules then / eV_J
MPL_EV = MPL_KG * (c**2) / eV_J

# 1 eV^-1 in meters
EVINV_TO_M = hcbar_SI / eV_J


def mass_density_kgm3_to_eV4(rho_kgm3: float) -> float:
    """Convert mass density (kg/m^3) to energy density in eV^4 via rho*c^2."""
    rho_Jm3 = rho_kgm3 * c**2
    return rho_Jm3 / EV4_TO_J_PER_M3


def number_density_cm3_to_kgm3(n_cm3: float, particle_mass_kg: float = constants.m_p) -> float:
    """Convert number density in cm^-3 to kg/m^3."""
    n_m3 = n_cm3 * 1e6
    return n_m3 * particle_mass_kg


# -----------------------------
# Physics model core
# -----------------------------

@dataclass(frozen=True)
class ModelParams:
    n: float
    beta: float
    Lambda_eV: float


def phi_min_eV(rho_eV4: float, p: ModelParams) -> float:
    """phi_min(rho) under A(phi)≈1 approximation."""
    n, beta, Lam = p.n, p.beta, p.Lambda_eV
    return (n * (Lam ** (4.0 + n)) * MPL_EV / (beta * rho_eV4)) ** (1.0 / (n + 1.0))


def m_eff_eV(rho_eV4: float, p: ModelParams) -> float:
    """Effective mass m_eff(rho) under A(phi)≈1 approximation."""
    n, beta, Lam = p.n, p.beta, p.Lambda_eV
    phi = phi_min_eV(rho_eV4, p)
    Vpp = n * (n + 1.0) * (Lam ** (4.0 + n)) / (phi ** (n + 2.0))
    matter_term = (beta**2) * rho_eV4 / (MPL_EV**2)
    return math.sqrt(Vpp + matter_term)


@dataclass(frozen=True)
class Body:
    name: str
    mass_kg: float
    radius_m: float
    rho_avg_kgm3: float

    @property
    def Phi_N(self) -> float:
        """Dimensionless Newtonian potential at surface: GM/(Rc^2)."""
        return constants.G * self.mass_kg / (self.radius_m * c**2)

    @property
    def rho_avg_eV4(self) -> float:
        return mass_density_kgm3_to_eV4(self.rho_avg_kgm3)


@dataclass(frozen=True)
class Environment:
    name: str
    rho_kgm3: float  # mass density
    @property
    def rho_eV4(self) -> float:
        return mass_density_kgm3_to_eV4(self.rho_kgm3)


def thin_shell_deltaR_over_R(body: Body, env: Environment, p: ModelParams) -> float:
    """ΔR/R = (φ∞ - φc) / (6 β Mpl Φ_N)."""
    phi_inf = phi_min_eV(env.rho_eV4, p)
    phi_c = phi_min_eV(body.rho_avg_eV4, p)
    return (phi_inf - phi_c) / (6.0 * p.beta * MPL_EV * body.Phi_N)


def beta_eff(body: Body, env: Environment, p: ModelParams) -> float:
    """Effective coupling including thin-shell screening."""
    ts = thin_shell_deltaR_over_R(body, env, p)
    if ts <= 0:
        return 0.0
    if ts < 1.0:
        return 3.0 * p.beta * ts
    return p.beta


def yukawa_alpha(body_i: Body, body_j: Body, env: Environment, p: ModelParams) -> float:
    """α_ij ≈ 2 β_eff,i β_eff,j."""
    return 2.0 * beta_eff(body_i, env, p) * beta_eff(body_j, env, p)


def lambda_ambient_m(env: Environment, p: ModelParams) -> float:
    """Ambient range λ∞ in meters."""
    m_inf_eV = m_eff_eV(env.rho_eV4, p)
    if m_inf_eV == 0:
        return float("inf")
    return (1.0 / m_inf_eV) * EVINV_TO_M


# -----------------------------
# Constraint curves (optional)
# -----------------------------

@dataclass
class AlphaLambdaConstraint:
    """A constraint of the form alpha <= alpha_max(lambda)."""
    name: str
    alpha_max_of_lambda: Callable[[float], float]  # lambda in meters
    r_probe_m: Optional[float] = None  # separation at which to apply exp(-r/lambda), if desired

    def passes(self, alpha: float, lambda_m: float) -> bool:
        if not math.isfinite(lambda_m) or lambda_m <= 0:
            return False
        suppress = 1.0
        if self.r_probe_m is not None:
            suppress = math.exp(-self.r_probe_m / lambda_m)
        return alpha * suppress <= self.alpha_max_of_lambda(lambda_m)


def load_constraint_csv(path: Path, name: str, r_probe_m: Optional[float] = None) -> Optional[AlphaLambdaConstraint]:
    """
    Load a CSV with columns: lambda_m, alpha_max
    Produces a log-log interpolation alpha_max(lambda).
    """
    if not path.exists():
        return None
    df = pd.read_csv(path)
    if "lambda_m" not in df.columns or "alpha_max" not in df.columns:
        raise ValueError(f"{path} must contain columns lambda_m, alpha_max")
    df = df.sort_values("lambda_m")
    lam = df["lambda_m"].values.astype(float)
    amax = df["alpha_max"].values.astype(float)

    # log-log interpolation with safe bounds
    loglam = np.log(lam)
    loga = np.log(amax)

    def alpha_max(lmb: float) -> float:
        if lmb <= lam[0]:
            return float(amax[0])
        if lmb >= lam[-1]:
            return float(amax[-1])
        x = np.log(lmb)
        y = np.interp(x, loglam, loga)
        return float(np.exp(y))

    return AlphaLambdaConstraint(name=name, alpha_max_of_lambda=alpha_max, r_probe_m=r_probe_m)


# -----------------------------
# Cassini-style PPN gate
# -----------------------------

@dataclass
class CassiniGate:
    """
    Conservative gate based on Cassini gamma measurement.
    Uses: 2 * beta_eff_sun^2 <= |gamma-1|_max  (long-range approximation).
    """
    gamma_minus_1_mean: float = 2.1e-5
    gamma_minus_1_sigma: float = 2.3e-5
    nsigma: float = 2.0  # 2σ by default

    def gamma_max_abs(self) -> float:
        return abs(self.gamma_minus_1_mean) + self.nsigma * self.gamma_minus_1_sigma

    def passes(self, beta_eff_sun: float) -> bool:
        return 2.0 * beta_eff_sun**2 <= self.gamma_max_abs()


# -----------------------------
# Default bodies and environments (conservative placeholders)
# -----------------------------

SUN = Body("Sun", mass_kg=1.98847e30, radius_m=6.9634e8, rho_avg_kgm3=1408.0)
EARTH = Body("Earth", mass_kg=5.9722e24, radius_m=6.371e6, rho_avg_kgm3=5515.0)

# Bennu: very approximate bulk parameters; replace with precise values if needed.
BENNU = Body("Bennu", mass_kg=7.329e10, radius_m=246.0, rho_avg_kgm3=1190.0)

# Solar wind at 1 AU: number density ~ 5 cm^-3 (order of magnitude).
SOLAR_WIND = Environment("SolarWind_1AU", rho_kgm3=number_density_cm3_to_kgm3(5.0))

# Lab vacuum placeholder (very rough): equivalent to ~1e-12 kg/m^3 (edit to your experiment).
LAB_VAC = Environment("LabVacuum", rho_kgm3=1e-12)

# Lab test mass (e.g. Tungsten sphere in torsion balance)
# Radius ~ 2cm, Mass ~ 700g, Density ~ 19300 kg/m^3
LAB_TEST_BODY = Body("LabSphere", mass_kg=0.7, radius_m=0.02, rho_avg_kgm3=19300.0)


# -----------------------------
# Scan driver
# -----------------------------

def run_scan(
    n: float = 1.0,
    beta_grid: np.ndarray = None,
    Lambda_grid_eV: np.ndarray = None,
    env_solar: Environment = SOLAR_WIND,
    env_lab: Environment = LAB_VAC,
    out_csv: Path = Path("viability_scan.csv"),
) -> pd.DataFrame:

    if beta_grid is None:
        # AGGRESSIVE SCAN requested by 00003
        # Push beta up to 10,000 (Strong Coupling)
        beta_grid = np.logspace(-2, 4, 25) 
    if Lambda_grid_eV is None:
        # AGGRESSIVE SCAN requested by 00003
        # Range: from 1e-10 eV (very light) to 1 eV (heavy)
        Lambda_grid_eV = np.logspace(-10, 0, 25)

    # Optional constraint curves (CSV files in cwd)
    eotwash = load_constraint_csv(Path("constraints_eotwash.csv"), "EotWash_Yukawa", r_probe_m=1e-3)
    eph = load_constraint_csv(Path("constraints_ephemerides.csv"), "Ephemerides_Yukawa", r_probe_m=1.0*constants.au)
    # bennu = load_constraint_csv(Path("constraints_bennu.csv"), "Bennu_Yukawa", r_probe_m=1.0*constants.au)

    constraints = [c for c in [eotwash, eph] if c is not None]

    cassini = CassiniGate()

    rows: List[Dict] = []
    for beta in beta_grid:
        for Lam in Lambda_grid_eV:
            p = ModelParams(n=n, beta=float(beta), Lambda_eV=float(Lam))

            # Ambient ranges
            lam_solar_m = lambda_ambient_m(env_solar, p)
            lam_lab_m = lambda_ambient_m(env_lab, p)

            # Screening / couplings
            beff_sun = beta_eff(SUN, env_solar, p)
            beff_earth = beta_eff(EARTH, env_solar, p)
            beff_bennu = beta_eff(BENNU, env_solar, p)

            # Core amplitudes
            alpha_sun_earth = 2.0 * beff_sun * beff_earth
            # alpha_sun_bennu = 2.0 * beff_sun * beff_bennu

            # Gate A: Cassini
            pass_cassini = cassini.passes(beff_sun)

            # Optional Yukawa constraints:
            # - For lab, use lab environment and a representative test pair (Earth-like materials placeholder)
            beff_lab_src = beta_eff(EARTH, env_lab, p)  # crude placeholder using Earth properties for shielding calc
# Optional Yukawa constraints:
            # CORRECTED: Use small test body, not Earth, to check for failed screening
            beff_lab_test = beta_eff(LAB_TEST_BODY, env_lab, p)
            alpha_lab = 2.0 * beff_lab_test * beff_lab_test

            pass_all_optional = True
            optional_results = {}
            for con in constraints:
                if con.name.lower().startswith("eotwash"):
                    ok = con.passes(alpha_lab, lam_lab_m)
                # elif con.name.lower().startswith("bennu"):
                #     ok = con.passes(alpha_sun_bennu, lam_solar_m)
                else:
                    ok = con.passes(alpha_sun_earth, lam_solar_m)
                optional_results[con.name] = ok
                pass_all_optional = pass_all_optional and ok

            rows.append({
                "n": n,
                "beta": beta,
                "Lambda_eV": Lam,
                "env_solar": env_solar.name,
                "rho_solar_kgm3": env_solar.rho_kgm3,
                "lambda_solar_m": lam_solar_m,
                "lambda_lab_m": lam_lab_m,
                "beta_eff_sun": beff_sun,
                "beta_eff_earth": beff_earth,
                "alpha_sun_earth": alpha_sun_earth,
                "pass_cassini": pass_cassini,
                **{f"pass_{k}": v for k, v in optional_results.items()},
                "pass_all": pass_cassini and pass_all_optional,
                "skipped_optional_gates": len(constraints) == 0,
            })

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)

    # console summary
    print(f"Saved: {out_csv.resolve()}")
    print(f"Total points: {len(df)}")
    print(f"Pass Cassini: {df['pass_cassini'].mean():.3f}")
    if any(col.startswith("pass_") and col != "pass_all" and col != "pass_cassini" for col in df.columns):
        for col in [c for c in df.columns if c.startswith("pass_") and c not in ("pass_all","pass_cassini")]:
            print(f"{col}: {df[col].mean():.3f}")
    print(f"Pass ALL gates: {df['pass_all'].mean():.3f}")
    if df["skipped_optional_gates"].all():
        print("NOTE: Optional Yukawa gates were skipped (no constraint CSVs found).")

    return df


if __name__ == "__main__":
    run_scan()