"""
Model validation, in two kinds of rows:
  check   model.py against a hand calculation or literature value, with a tolerance (PASS / FAIL)
  metric  comparison with the digitized Yao 2025 LSV data the S1 parameters are fitted to; the number is
          reported as it is, with no pass/fail threshold (a fit to the same data is not a test of the model)
Usage: run_all() returns all rows (list of dict); from the command line use scripts/validate.py.
"""
import re
import numpy as np
from . import model as m
from .yao2025 import REPORTED, FIG3B, HGO_TO_RHE
from . import calibration as cal

ROWS = []
VERBOSE = False


def check(cat, name, value, expect, tol, how, note=""):
    ok = abs(value - expect) <= tol
    ROWS.append(dict(kind="check", cat=cat, name=name, value=value, expect=expect, tol=tol, ok=bool(ok),
                     how=how, note=note))
    if VERBOSE:
        print(f"[{'PASS' if ok else 'FAIL'}] {cat:6s} {name:42s} model {value:>11.5g}  expected {expect:>11.5g}  (±{tol:g})")


def flag(cat, name, ok, how, note=""):
    ROWS.append(dict(kind="check", cat=cat, name=name, value=None, expect=None, tol=None, ok=bool(ok), how=how,
                     note=note))
    if VERBOSE:
        print(f"[{'PASS' if ok else 'FAIL'}] {cat:6s} {name}")


def metric(cat, name, value, ref, how, note=""):
    """Reported number without a threshold; ref is the published or model value it is compared with (or None)."""
    ROWS.append(dict(kind="metric", cat=cat, name=name, value=value, expect=ref, tol=None, ok=None,
                     how=how, note=note))
    if VERBOSE:
        ref_txt = "" if ref is None else f"  reference {ref:>11.5g}"
        print(f"[METRIC] {cat:6s} {name:42s} model {value:>11.5g}{ref_txt}")


def _run():
    # ---------------------------------------------------------------------------
    # A. Stoichiometry: conservation of atoms and charge
    # ---------------------------------------------------------------------------
    def atoms(formula):
        """Minimal formula parser; supports parentheses, e.g. CO(NH2)2."""
        def parse(s):
            cnt, i = {}, 0
            while i < len(s):
                if s[i] == "(":
                    depth, j = 1, i + 1
                    while depth:
                        depth += {"(": 1, ")": -1}.get(s[j], 0); j += 1
                    inner = parse(s[i + 1:j - 1])
                    k = re.match(r"\d*", s[j:]).group(); mult = int(k) if k else 1
                    for e, n in inner.items():
                        cnt[e] = cnt.get(e, 0) + n * mult
                    i = j + len(k)
                else:
                    mm = re.match(r"([A-Z][a-z]?)(\d*)", s[i:])
                    e, k = mm.group(1), mm.group(2)
                    cnt[e] = cnt.get(e, 0) + (int(k) if k else 1)
                    i += len(mm.group(0))
            return cnt
        return parse(formula)


    def balance(lhs, rhs):
        """lhs / rhs: [(coefficient, formula, charge)]; returns the non-zero imbalances {element or 'charge': difference}."""
        tot = {}
        for sgn, side in ((1, lhs), (-1, rhs)):
            for n, f, z in side:
                for e, k in atoms(f).items():
                    tot[e] = tot.get(e, 0) + sgn * n * k
                tot["charge"] = tot.get("charge", 0) + sgn * n * z
        return {k: v for k, v in tot.items() if v != 0}


    E = ("e", "", -1)
    REACTIONS = {
        "UOR: CO(NH2)2 + 6OH- → N2 + 5H2O + CO2 + 6e-":
            ([(1, "CO(NH2)2", 0), (6, "OH", -1)], [(1, "N2", 0), (5, "H2O", 0), (1, "CO2", 0), (6, "", -1)]),
        "NO3RR: NO3- + 6H2O + 8e- → NH3 + 9OH-":
            ([(1, "NO3", -1), (6, "H2O", 0), (8, "", -1)], [(1, "NH3", 0), (9, "OH", -1)]),
        "Overall: 4CO(NH2)2 + 3NO3- → 4N2 + 4CO2 + 3NH3 + 2H2O + 3OH-":
            ([(4, "CO(NH2)2", 0), (3, "NO3", -1)],
             [(4, "N2", 0), (4, "CO2", 0), (3, "NH3", 0), (2, "H2O", 0), (3, "OH", -1)]),
        "UOR (carbonate): CO(NH2)2 + 8OH- → N2 + 6H2O + CO3^2- + 6e-":
            ([(1, "CO(NH2)2", 0), (8, "OH", -1)], [(1, "N2", 0), (6, "H2O", 0), (1, "CO3", -2), (6, "", -1)]),
    }
    for nm, (l, r) in REACTIONS.items():
        d = balance(l, r)
        flag("Stoichiometry", nm, not d, "count C, H, N, O atoms and charge on each side; the difference must be 0", str(d) if d else "")

    # ---------------------------------------------------------------------------
    # B. Thermodynamics
    # ---------------------------------------------------------------------------
    th = m.thermo_check()
    check("Thermodynamics", "E°(NO3-/NH3) vs RHE", th["NO3RR"], m.E0["NO3RR"], 0.001, "recomputed from ΔfG° vs tabulated value")
    check("Thermodynamics", "E°(UOR, CO2) vs RHE", th["UOR_calc"], 0.072, 0.001, "compared with 0.072 V from Protsenko 2023")
    check("Thermodynamics", "E°(UOR, CO3^2-) vs RHE", th["UOR_CO3"], -0.024, 0.002,
          "compared with the 0.024 V EMF (UOR + HER) of Protsenko 2024")
    check("Thermodynamics", "ΔH°/6F (urea electrolysis)", th["UOR_dH"], 0.37, 0.01,
          "tests the guess that 0.37 V is an enthalpy-based value", "0.363 V; numerically close only, literature source not confirmed")
    check("Thermodynamics", "E°cell (UOR = 0.37 V)", m.cell_E0(0.37)[0], 0.338, 0.001, "0.708 − 0.37")
    check("Thermodynamics", "ΔG° (UOR = 0.37 V), kJ/mol NO3-", m.cell_E0(0.37)[1], -8 * 96.485 * 0.338, 0.5,
          "ΔG = −nFE, n = 8")
    check("Thermodynamics", "Nernst: 1000x NO3- dilution (mV)", m.nernst_shift(1e-3) * 1e3,
          -8.314 * 298.15 / 96485 / 8 * np.log(1000) * 1e3, 0.01, "(RT/8F) ln 1000 = 22.2 mV")

    # ---------------------------------------------------------------------------
    # C. Kinetics and mass transport (vs hand calculation)
    # ---------------------------------------------------------------------------
    for k, p in m.SCEN.items():
        check("Kinetics", f"{k}: j_UOR(E10) = 10 mA cm-2", m.j_uor(p["E10_UOR"], p), 10.0, 1e-9,
              "by construction of the Tafel expression the current at E10 is 10")
        check("Kinetics", f"{k}: j_NO3,k(E10) = 10 mA cm-2 (C = C_ref)",
              m.tafel_current(m.E0["NO3RR"] - p["E10_NO3"], m.E0["NO3RR"] - p["E10_NO3"], p["b_NO3"]),
              10.0, 1e-9, "as above")
    # the numerical slope must equal b
    p = m.S1
    E1, E2 = 1.40, 1.45
    slope = (E2 - E1) / np.log10(m.j_uor(E2, p) / m.j_uor(E1, p))
    check("Kinetics", "UOR numerical Tafel slope (mV/dec)", slope * 1e3, p["b_UOR"] * 1e3, 0.01,
          "slope from two points of E vs log j; must equal the input b")
    jl = m.j_lim_no3(100, 43.0)
    hand = 8 * 96485 * 1.902e-5 * (100e-3 / 1000) / (43e-4) * 1e3
    check("Mass transport", "j_lim(100 mM, 43 µm) (mA cm-2)", jl, hand, 1e-6,
          "hand calculation nFDC/δ: 8 × 96485 × 1.902e-5 × 1e-4 / 4.3e-3")
    check("Mass transport", "j_lim proportional to concentration (10 mM / 1 mM)", m.j_lim_no3(10, 43) / m.j_lim_no3(1, 43), 10, 1e-9,
          "linear relation")
    q = dict(m.S1, delta_um=43.0)
    check("Mass transport", "K–L: j_NO3 → j_lim at very negative potential (20 mM)", m.j_no3(-1.5, q, 20) / m.j_lim_no3(20, 43),
          1.0, 1e-3, "with fast kinetics the current is set by diffusion")
    check("Mass transport", "K–L: j_NO3 ≈ j_k when j ≪ j_lim (100 mM, 0.65 V)",
          m.j_no3(0.65, q, 100) / m.tafel_current(m.E0["NO3RR"] - 0.65, m.E0["NO3RR"] - q["E10_NO3"], q["b_NO3"]),
          1.0, 1e-3, "with fast diffusion the current is set by kinetics")

    # ---------------------------------------------------------------------------
    # D. Numerics
    # ---------------------------------------------------------------------------
    js = np.linspace(0.5, 40, 80)
    worst_a = worst_c = 0.0
    mono = bounds = True
    for k, p in list(m.SCEN.items()) + [("OER", m.OER_ONLY)]:
        prev = None
        for j in js:
            r = m.solve_cell(j, p)
            worst_a = max(worst_a, abs(m.j_uor(r["Ea"], p) + m.j_oer(r["Ea"], p) - j))
            worst_c = max(worst_c, abs(m.j_no3(r["Ec"], p, p["C_NO3"]) + m.j_her(r["Ec"], p) - j))
            if prev and not (r["Ea"] > prev["Ea"] and r["Ec"] < prev["Ec"]):
                mono = False
            if not (0 <= r["FE_NO3"] <= 1 + 1e-9 and 0 <= r["FE_UOR"] <= 1 + 1e-9
                    and r["Ec"] <= m.E0["NO3RR"] and r["V"] < m.cell_E0(p["E0_UOR"])[0]):
                bounds = False
            prev = r
    check("Numerics", "max anode current-balance residual (mA cm-2)", worst_a, 0, 1e-8, "substituted back: j_UOR + j_OER − j")
    check("Numerics", "max cathode current-balance residual (mA cm-2)", worst_c, 0, 1e-8, "substituted back: j_NO3 + j_HER − j")
    flag("Numerics", "E_a rises and E_c falls monotonically with current (4 sets × 80 points)", mono, "polarization curves must not fold back")
    flag("Numerics", "0 ≤ FE ≤ 1, E_c ≤ E°(NO3RR), V_cell < E°cell (all points)", bounds, "physical bounds")
    # compare with the analytic inverse (S3, 10 mA, side reactions negligible)
    p = m.SCEN["S3"]
    eta10 = p["E10_UOR"] - p["E0_UOR"]
    Ea_exact = p["E0_UOR"] + p["b_UOR"] * np.log10(10 / 10 * (10 ** (eta10 / p["b_UOR"]) - 1) + 1)
    check("Numerics", "S3 anode potential: brentq vs analytic inverse (V)", m.anode_potential(10, p), Ea_exact, 1e-6,
          "the Tafel expression inverts directly: E = E° + b log10[(j/j0) + 1]")

    # ---------------------------------------------------------------------------
    # E. Comparison with the Yao 2025 data (calibration and held-out curve): metrics, no threshold
    # ---------------------------------------------------------------------------
    an = cal.fit_anode()
    metric("Calibration", "anode LSV: model parameters reproduce Fig. 11a (RMS relative error)", an["used"]["rms"], None,
           "UOR + OER curves, 2–175 mA cm-2, only R_u fitted", f"R_u = {an['used']['Ru']:.2f} Ω cm2")
    metric("Calibration", "anode E10(UOR) after iR correction (V)", an["fixed_b"]["E10_UOR"], REPORTED["UOR_E10_ANP"],
           "fit with b fixed at 44.2 mV/dec; reference: 1.335 V reported by Yao 2025")
    metric("Calibration", "UOR–OER potential gap at 10 mA (mV)",
           (m.anode_potential(10, m.OER_ONLY) - m.anode_potential(10, m.S1)) * 1e3, 180,
           "reference: 180 mV marked in Fig. 11a")
    ca = cal.fit_cathode(gamma=1.0)
    metric("Calibration", "cathode LSV: fit RMS, 0/1/10/100 mM", ca["rms_fit"], None,
           f"4 curves, {ca['n_points']} points, {len(ca['std'])} parameters")
    metric("Calibration", "cathode LSV: 100 mM curve RMS", ca["rms_by_C"][100], None, "the concentration mainly used by the model")
    # where the published 50 mM curve carries more current than the 100 mM curve; no model with a positive
    # reaction order can follow it there, since every current term rises with concentration
    E50, J50 = FIG3B[50]
    used = np.isfinite(J50) & (J50 >= 3.0)                 # the points cathode_resid uses, in the same order
    inverted = (J50 > FIG3B[100][1])[used]
    r50 = cal.cathode_resid(ca["params"], (50,))
    rms_pct = lambda r: 100 * float(np.sqrt(np.mean(r ** 2)))
    e_lo, e_hi = E50[used][inverted].min(), E50[used][inverted].max()
    metric("Held-out test", "cathode LSV: 50 mM curve (not in the fit) RMS", ca["rms_by_C"][50], None,
           "independent prediction",
           f"between {e_lo:+.2f} and {e_hi:+.2f} V vs RHE ({e_lo - HGO_TO_RHE:.2f} to {e_hi - HGO_TO_RHE:.2f} V vs "
           f"Hg/HgO, the scale of Fig. 3b) Yao's 50 mM curve lies above the 100 mM curve, against the concentration "
           f"trend: RMS {rms_pct(r50[inverted]):.0f} % over these {inverted.sum()} points, "
           f"{rms_pct(r50[~inverted]):.0f} % over the other {(~inverted).sum()}")
    metric("Calibration", "fitted E10(NO3RR) vs value used (V)", ca["params"]["E10_NO3"], m.S1["E10_NO3"],
           "reference: value in model.py")
    metric("Calibration", "fitted b(NO3RR) vs value used (V/dec)", ca["params"]["b_NO3"], m.S1["b_NO3"],
           "reference: value in model.py")
    metric("Calibration", "fitted E10(HER) vs value used (V)", ca["params"]["E10_HER"], m.S1["E10_HER"],
           "reference: value in model.py")
    metric("Calibration", "fitted δ vs value used (µm)", ca["params"]["delta_um"], m.S1["delta_um"],
           "reference: value in model.py")
    cf = cal.fit_cathode()
    metric("Calibration", "reaction order γ (free fit)", cf["params"]["gamma"], 1.0,
           "reference: γ = 1 used in the model", f"γ = {cf['params']['gamma']:.2f} ± {cf['std']['gamma']:.2f} (standard error)")

    # ---------------------------------------------------------------------------
    # F. Batch simulation: mass balance, time step, energy sign
    # ---------------------------------------------------------------------------
    r = m.simulate_batch(m.S1, "cc")
    check("Batch", "mass balance: NO3- removed × 8F = charge used for NO3RR", r["n_N_no3"] * 8 * m.F / r["Q_no3"], 1.0,
          1e-9, "both sides accumulated independently")
    flag("Batch", "concentration decreases monotonically with time", bool(np.all(np.diff(r["C"]) <= 0)), "")
    # analytic time for FE = 1: t = 8F ΔC V / I
    p_ideal = dict(m.S1, delta_um=1e-3, E10_HER=-5.0)       # makes FE ≡ 1
    ri = m.simulate_batch(p_ideal, "cc")
    t_hand = 8 * m.F * (0.9 * m.BATCH["C0"] * 1e-3 * m.BATCH["V_L"]) / (10e-3 * m.BATCH["A_cm2"]) / 3600
    check("Batch", "treatment time at FE = 1 (h)", ri["t_end"], t_hand, m.BATCH["dt"] / 3600,
          "Faraday's law by hand: t = nFΔCV / I")
    r2 = m.simulate_batch(m.S1, "cc", dt=2.0)
    check("Batch", "time step 10 s vs 2 s: relative difference in treatment time", abs(r["t_end"] / r2["t_end"] - 1), 0, 0.005, "convergence")
    check("Batch", "time step 10 s vs 2 s: relative difference in energy", abs(r["E_kWh"] / r2["E_kWh"] - 1), 0, 0.005, "convergence")
    rf = m.simulate_batch(m.S1, "follow")
    r2f = m.simulate_batch(m.S1, "follow", dt=2.0)
    check("Batch", "time step (following strategy): relative difference in energy", abs(rf["E_kWh"] / r2f["E_kWh"] - 1), 0, 0.005, "convergence")
    flag("Batch", "S1: V_cell < 0 throughout (consumes power), energy negative", bool(np.all(r["V"] < 0) and r["E_kWh"] < 0),
         "sign convention: power produced is positive")
    check("Batch", "S3 has a positive voltage at 10 mA cm-2 (produces power)", float(m.solve_cell(10, m.SCEN["S3"])["V"] > 0), 1, 0, "")

    # ---------------------------------------------------------------------------
    # G. Straight-line approximation of the operating-mode map
    # ---------------------------------------------------------------------------
    worst, n_in, n_out = 0.0, 0, 0
    for xa in np.linspace(0.45, 1.6, 24):
        for yc in np.linspace(-0.4, 0.68, 24):
            rr = m.solve_cell(10, dict(m.S1, E10_UOR=xa, E10_NO3=yc))
            if rr["FE_NO3"] >= 0.9 and rr["FE_UOR"] >= 0.9:
                n_in += 1
                worst = max(worst, abs(rr["V"] - (yc - xa - 0.01)))
            else:
                n_out += 1
    check("Map", "main-reaction region (FE ≥ 90 % at both electrodes): max difference, line V = y − x − jR vs full model (mV)", worst * 1e3, 0, 5,
          f"{n_in} of 576 grid points lie in the main-reaction region", f"in the other {n_out} points a side reaction "
          "(HER or OER) carries >10 % of the current and the line fails; the map therefore uses the full model for "
          "the iso-voltage lines and marks the side-reaction region")


def run_all(verbose=False):
    """Run all rows; return a list of dict (kind, cat, name, value, expect, tol, ok, how, note).
    kind is "check" (ok is True/False) or "metric" (ok and tol are None)."""
    global VERBOSE
    VERBOSE = verbose
    ROWS.clear()
    _run()
    return list(ROWS)
