"""
Simplified model of a urea–nitrate direct cell (UOR || NO3RR): core calculations (no plotting)
S1 parameters are fitted to the original figure data of Yao 2025 (including the reaction order; HER from the 0 mM curve)

This file does four things, each of which can be checked by hand (see validation.py):
  1. Thermodynamics: standard potentials of the half-reactions from Gibbs energies of formation
  2. Single electrode: potential → current (Tafel kinetics + diffusion limit)
  3. Full cell: current → both electrode potentials (by inversion) → cell voltage
  4. Batch treatment: nitrate concentration falling over time (mass balance + Faraday's law)

Units throughout: potential V vs RHE (pH 14), current density mA cm-2, concentration mM, diffusion layer µm, resistance Ω cm2
Sign convention: V_cell = E_c − E_a − jR_Ω; > 0 means the cell delivers power, < 0 means an external bias is needed
"""
import numpy as np
from scipy.optimize import brentq

F = 96485.0           # Faraday constant, C mol-1
RT_F = 8.314 * 298.15 / F   # 0.02569 V
D_NO3 = 1.902e-5      # diffusion coefficient of NO3- in water, cm2 s-1 (25 °C, CRC Handbook)
N_E_NO3 = 8           # electrons transferred, NO3- -> NH3
N_E_UOR = 6           # electrons transferred, CO(NH2)2 -> N2 + CO2

# ===========================================================================
# 1. Thermodynamics
# ===========================================================================
# Standard Gibbs energies of formation ΔfG° (kJ mol-1, 298.15 K; NBS tables, Wagman et al. 1982)
DG_F = {"H2O": -237.13, "OH-": -157.24, "CO2(g)": -394.36, "CO3^2-": -527.81,
        "NO3-": -111.25, "NH3(aq)": -26.50,
        # urea(aq) = urea(s) -197.33 + Gibbs energy of solution (ΔH_sol 16.91 kJ, ΔS_sol 62.59 J/K)
        # Note: the enthalpy basis in thermo_check() uses the NBS ΔfH°(urea, aq) = -317.65; with
        #       ΔfH°(urea, s) = -333.51 this implies ΔH_sol = 15.86 kJ, a different source from the 16.91
        #       used here. Using 15.86 would move E°(UOR, CO2) from 0.0723 to 0.0741 V (+1.8 mV) and
        #       change no conclusion, so the original value is kept.
        "urea(aq)": -197.33 + (16.91 - 298.15 * 62.59e-3)}

# Standard potentials (V vs RHE, pH 14); thermo_check() recomputes them for comparison with this table
E0 = {
    "HER": 0.000,
    "UOR_lit": 0.37,      # commonly cited value (Boggs 2009: -0.46 V vs SHE)
    "UOR_calc": 0.072,    # recomputed from ΔfG°, product CO2(g) (same value as Protsenko 2023)
    "UOR_CO3": -0.025,    # product CO3^2- (the actual case at pH 14; Protsenko 2024 reports an EMF of 0.024 V)
    "NO3RR": 0.708,       # NO3- -> NH3(aq)
    "OER": 1.229,
}


def thermo_check(verbose=False):
    """Recompute the standard potentials from ΔfG°; returns a dict (V vs RHE).
    Reduction potential E° = -ΔG_red / (nF); vs RHE = vs SHE − E°(HER, pH 14)."""
    G, k = DG_F, 96.485   # k = F in kJ V-1 mol-1
    e_her = -(2 * G["OH-"] - 2 * G["H2O"]) / (2 * k)                    # 2H2O+2e- -> H2+2OH-
    e_uor = -((G["urea(aq)"] + 6 * G["OH-"]) - (5 * G["H2O"] + G["CO2(g)"])) / (6 * k)
    e_uor3 = -((G["urea(aq)"] + 8 * G["OH-"]) - (6 * G["H2O"] + G["CO3^2-"])) / (6 * k)
    e_no3 = -((G["NH3(aq)"] + 9 * G["OH-"]) - (G["NO3-"] + 6 * G["H2O"])) / (8 * k)
    out = dict(HER_SHE=e_her, UOR_calc=e_uor - e_her, UOR_CO3=e_uor3 - e_her,
               NO3RR=e_no3 - e_her)
    # Enthalpy basis: CO(NH2)2(aq) + H2O(l) -> N2 + CO2(g) + 3H2, ΔH°/6F (thermoneutral voltage); urea ΔfH°: see the note at DG_F
    dH = -393.51 - (-317.65) - (-285.83)
    out["UOR_dH"] = dH / (6 * k)
    if verbose:
        for kk, v in out.items():
            print(f"  {kk:9s} {v:+.4f} V")
    return out


def cell_E0(E0_uor):
    """Standard cell potential and ΔG° (per mol NO3-, n = 8)."""
    E = E0["NO3RR"] - E0_uor
    return E, -N_E_NO3 * F * E / 1e3   # V, kJ mol-1


def nernst_shift(factor, n=N_E_NO3):
    """Shift of the equilibrium potential (V) when the NO3- activity changes by `factor`, other activities unchanged."""
    return RT_F / n * np.log(factor)


# ===========================================================================
# 2. Current–potential relation of a single electrode
# ===========================================================================
def tafel_current(eta, eta10, b):
    """Tafel-type kinetic current (mA cm-2, positive).
    j = j0 (10^(η/b) − 1), where j0 follows from j = 10 mA cm-2 at η = η10:
        j0 = 10 / (10^(η10/b) − 1)
    The current is 0 for η ≤ 0 (forward reaction only; the reverse reaction is negligible at the large overpotentials here).
    Note: j0 is only a scaling constant, not a measurable exchange current density."""
    if eta <= 0:
        return 0.0
    return 10.0 * (10 ** (eta / b) - 1) / (10 ** (eta10 / b) - 1)


def j_lim_no3(C_mM, delta_um):
    """Diffusion-limited current density of NO3- (mA cm-2): j_lim = nFDC/δ."""
    C = C_mM * 1e-6            # mol cm-3
    delta = delta_um * 1e-4    # cm
    return N_E_NO3 * F * D_NO3 * C / delta * 1e3


def j_uor(E, p):
    if not p.get("UOR_on", True):        # control case: no urea at the anode
        return 0.0
    return tafel_current(E - p["E0_UOR"], p["E10_UOR"] - p["E0_UOR"], p["b_UOR"])


def j_oer(E, p):
    return tafel_current(E - E0["OER"], p["E10_OER"] - E0["OER"], p["b_OER"])


def j_no3(E, p, C):
    """NO3RR current: kinetics × concentration^γ, in series with the diffusion limit (Koutecký–Levich).
    E10_NO3 refers to C_ref (100 mM); γ is the reaction order."""
    jk = tafel_current(E0["NO3RR"] - E, E0["NO3RR"] - p["E10_NO3"], p["b_NO3"])
    jk *= (C / p["C_ref"]) ** p["gamma"]
    if jk <= 0:
        return 0.0
    jl = j_lim_no3(C, p["delta_um"])
    return jk * jl / (jk + jl)


def j_her(E, p):
    return tafel_current(E0["HER"] - E, E0["HER"] - p["E10_HER"], p["b_HER"])


# ===========================================================================
# 3. Full cell
# ===========================================================================
def anode_potential(j, p):
    """Solve j_UOR(E) + j_OER(E) = j. The current rises monotonically with E, so brentq has a unique root."""
    lo = min(p["E0_UOR"], E0["OER"])
    return brentq(lambda E: j_uor(E, p) + j_oer(E, p) - j, lo, 3.0, xtol=1e-12)


def cathode_potential(j, p, C):
    """Solve j_NO3(E) + j_HER(E) = j. The current rises monotonically as E decreases."""
    return brentq(lambda E: j_no3(E, p, C) + j_her(E, p) - j, -3.0, E0["NO3RR"], xtol=1e-12)


def solve_cell(j, p, C=None):
    """For a current density j (mA cm-2), return dict(V, Ea, Ec, FE_UOR, FE_NO3)."""
    C = p["C_NO3"] if C is None else C
    Ea = anode_potential(j, p)
    Ec = cathode_potential(j, p, C)
    V = Ec - Ea - j * 1e-3 * p["R_ohm"]
    return dict(V=V, Ea=Ea, Ec=Ec, FE_UOR=j_uor(Ea, p) / j, FE_NO3=j_no3(Ec, p, C) / j)


# ===========================================================================
# 4. Parameters (each with its source)
# ===========================================================================
# Source tags:
#   [Y-reported]  value stated in the text of Yao 2025
#   [Y-fitted]    fitted to digitized data from Yao 2025 Fig. 3b / 11a (calibration.py)
#   [lit]         other literature / handbooks
#   [assumed]     assumption of this study, varied in the sensitivity analysis
S1 = dict(
    E0_UOR=E0["UOR_lit"],
    E10_UOR=1.335, b_UOR=0.0442,   # [Y-reported] ANP, 1 M NaOH + 0.5 M urea (Fig. 11b, 11c)
    E10_OER=1.515, b_OER=0.060,    # E10 [Y-reported]; b [Y-fitted] 55–65 mV/dec (Fig. 11a)
    E10_NO3=0.037, b_NO3=0.099,    # [Y-fitted] ANP, 100 mM NaNO3 + 1 M NaOH (Fig. 3a/3b), iR-corrected
    gamma=1.0, C_ref=100.0,        # reaction order [Y-fitted] 0.87 ± 0.14, first order used
    E10_HER=-0.202, b_HER=0.081,   # [Y-fitted] ANP, 0 mM NaNO3 (Fig. 3b, black line)
    delta_um=43.0,                 # diffusion layer [Y-fitted] (Fig. 3b 1 mM and 10 mM curves; stirred)
    C_NO3=100.0,                   # concentration for the polarization comparison = calibration concentration
    R_ohm=1.0,                     # [assumed] zero-gap AEM cell, 1 M KOH
)
SCEN = {
    "S1": S1,                                                           # current catalyst (Yao 2025 ANP)
    "S2": dict(S1, E10_UOR=0.60, b_UOR=0.040, E10_NO3=0.30, b_NO3=0.060),   # [assumed] improved
    "S3": dict(S1, E10_UOR=0.47, b_UOR=0.040, E10_NO3=0.60, b_NO3=0.060),   # [assumed] ideal
}
SCEN_LABEL = {"S1": "S1 current catalyst (Yao 2025)", "S2": "S2 improved (assumed)",
              "S3": "S3 ideal (assumed)"}
# Control: same cathode, no urea at the anode (OER only)
OER_ONLY = dict(S1, UOR_on=False)
# Other measured catalyst points (for the operating-mode map)
NF_POINT = dict(E10_UOR=1.423,    # [Y-reported] Ni foam UOR
                E10_NO3=-0.24)    # [Y-fitted] Ni foam NO3RR, 100 mM (Fig. 3a black line, includes about 0.02 V of iR)


# ===========================================================================
# 5. Sensitivity analysis (one parameter at a time)
# ===========================================================================
J_SENS = 20.0
SENS = [
    # (key, low, high)
    ("E10_UOR", 1.30, 1.42),
    ("E10_NO3", -0.06, 0.14),
    ("b_NO3", 0.070, 0.130),
    ("b_UOR", 0.023, 0.080),
    ("R_ohm", 0.5, 5.0),
    ("C_NO3", 10.0, 200.0),
    ("gamma", 0.5, 1.0),
    ("delta_um", 10.0, 100.0),
    ("E10_HER", -0.30, -0.10),
    ("E0_UOR", E0["UOR_CO3"], E0["UOR_lit"]),
]


def sensitivity(base=None, j=J_SENS):
    base = base or S1
    V0 = solve_cell(j, base)["V"]
    rows = []
    for key, lo, hi in SENS:
        dv = [solve_cell(j, dict(base, **{key: v}))["V"] - V0 for v in (lo, hi)]
        rows.append((key, lo, hi, dv[0], dv[1]))
    rows.sort(key=lambda r: abs(r[4] - r[3]))
    return V0, rows


# ===========================================================================
# 6. Batch treatment (nitrate concentration in the cathode chamber falling over time)
#    dC/dt = − FE_NO3 · I / (8 F V)    (Faraday's law; urea at the anode taken to be in excess)
# ===========================================================================
BATCH = dict(V_L=0.10, A_cm2=10.0, C0=60 / 14.007, dt=10.0, removal=0.90)  # C0 = 60 mg N/L


def simulate_batch(p, strategy, j_set=10.0, frac=0.8, **kw):
    """strategy = 'cc': constant current j_set; 'follow': j = min(j_set, frac × j_lim(C)).
    Integrated with explicit Euler until the removal target is reached. Energy E_J is positive when power is produced."""
    b = dict(BATCH, **kw)
    C, t = b["C0"], 0.0
    T, Cs, FEs, Vs, Js = [0.0], [C], [], [], []
    E_J = Q = Q_no3 = Q_uor = 0.0
    while C > (1 - b["removal"]) * b["C0"] and t < 72 * 3600:
        j = j_set if strategy == "cc" else min(j_set, frac * j_lim_no3(C, p["delta_um"]))
        r = solve_cell(j, p, C)
        I = j * 1e-3 * b["A_cm2"]                                   # A
        dC = r["FE_NO3"] * I / (N_E_NO3 * F) / b["V_L"] * 1e3 * b["dt"]   # mM
        C = max(C - dC, 0.0)
        t += b["dt"]
        E_J += r["V"] * I * b["dt"]
        Q += I * b["dt"]
        Q_no3 += r["FE_NO3"] * I * b["dt"]
        Q_uor += r["FE_UOR"] * I * b["dt"]
        T.append(t / 3600); Cs.append(C); FEs.append(r["FE_NO3"]); Vs.append(r["V"]); Js.append(j)
    n_N_no3 = (b["C0"] - C) * 1e-3 * b["V_L"]          # mol N (nitrate)
    n_N_urea = Q_uor / (N_E_UOR * F) * 2               # mol N (urea: 2 N per molecule)
    kWh = E_J / 3.6e6
    return dict(t=np.array(T), C=np.array(Cs), FE=np.array(FEs), V=np.array(Vs), j=np.array(Js),
                t_end=t / 3600, E_kWh=kWh, Q=Q, Q_no3=Q_no3, n_N_no3=n_N_no3, n_N_urea=n_N_urea,
                FE_avg=Q_no3 / Q,
                kWh_per_kgN_no3=kWh / (n_N_no3 * 14.007e-3),
                kWh_per_kgN_all=kWh / ((n_N_no3 + n_N_urea) * 14.007e-3))


if __name__ == "__main__":
    print("Thermodynamics recomputed (V vs RHE):"); thermo_check(verbose=True)
    for k in ("UOR_lit", "UOR_calc", "UOR_CO3"):
        E, G = cell_E0(E0[k])
        print(f"  E°cell (UOR = {k}) = {E:+.3f} V, ΔG° = {G:.0f} kJ/mol NO3-")
    print(f"  Nernst shift for a 1000x NO3- dilution: {nernst_shift(1e-3)*1e3:+.1f} mV")
    print("\nScenarios at 10 mA cm-2:")
    for k, p in SCEN.items():
        r = solve_cell(10.0, p)
        print(f"  {k}: V={r['V']:+.3f} Ea={r['Ea']:.3f} Ec={r['Ec']:+.3f} "
              f"FE_NO3={r['FE_NO3']*100:.1f}%")
    r = solve_cell(10.0, OER_ONLY)
    print(f"  OER control: V={r['V']:+.3f} Ea={r['Ea']:.3f}")
