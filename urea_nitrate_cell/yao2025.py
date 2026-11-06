"""
Data read from the figures of Yao et al., Sustainability 2025, 17, 3835 (doi:10.3390/su17093835)
------------------------------------------------------------------------
Digitization: the pixels of the journal's original figure files (pub.mdpi-res.com, *-g003.png,
*-g011.png, *-g014.png) were read in a browser, a linear mapping was built from the pixel positions
of the major axis ticks, and the curve centre was read every 10 mV.
Reading error: about ±1 mA cm-2 (line width) and ±1 mV in potential.

Potential conversion (Yao 2025, Sec. 2): E_RHE = E_Hg/HgO + 0.098 V + 0.059 × pH, pH = 14 (1 M NaOH)
"""
import numpy as np

HGO_TO_RHE = 0.098 + 0.059 * 14          # = 0.924 V

# ---------------------------------------------------------------------------
# Fig. 3b: ANP cathode LSV, 1 M NaOH + x mM NaNO3 (potential vs Hg/HgO; cathodic current, taken as positive)
#          potential of point k = -1.195 + 0.01 k  (k = 0..50)
# ---------------------------------------------------------------------------
_E3 = -1.195 + 0.01 * np.arange(51)
_nan = np.nan
_J3 = {
    0: [_nan, 18.8, 16.6, 14.3, 12.4, 10.7, 9.3, 8.1, 6.9, 6.1, 5.3, 4.7, 4.2, 3.9, 3.6, 3.3, 3.0,
        2.8, 2.8, 2.7, 2.7, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.2, 2.2, 2.2, 1.9, 1.9, 1.6,
        1.6, 1.4, 1.3, 1.3, 1.1, 1.1, 0.8, 0.8, 0.8, 0.5, 0.5, 0.3, 0.2, -0.1, -0.1, -0.3],
    1: [26.9, 24.5, 22.5, 20.2, 18.5, 16.8, 15.4, 14.3, 13.2, 12.2, 11.4, 10.7, 10.1, 9.5, 9.1, 8.6,
        8.3, 7.9, 7.6, 7.2, 6.8, 6.6, 6.4, 6.1, 5.8, 5.5, 5.2, 4.9] + [_nan] * 23,
    10: [41.0, 38.5, 36.2, 33.7, 31.5, 29.7, 27.7, 25.9, 24.3, 22.7, 21.0, 19.6, 18.2, 16.8, 15.4,
         14.3, 12.9, 11.7, 10.6, 9.6, 8.7, 8.1, 7.6, 6.9, 6.5, 6.1, 5.8, 5.4, 5.1, 4.8, 4.4, 4.2,
         4.0, 3.8, 3.5, 3.3, 3.2, 3.0, 2.8, 2.7, 2.6, 2.5, 2.3, 2.0, _nan, 2.4, 2.1, 2.0, 2.0, 1.9, 1.9],
    50: [67.4, 65.5, 63.7, 62.0, 60.2, 58.6, 57.0, 55.3, 53.9, 52.3, 50.8, 49.0, 47.6, 46.2, 44.6,
         43.1, 41.5, 39.9, _nan, 36.7, 34.9, 33.2, 31.5, 29.6, 27.6, 25.9, 24.0, 22.3, 20.4, 18.5,
         16.8, 15.0, 13.5, 12.0, 10.5, 9.3, 8.0, 6.8, 5.8, 5.0, 4.3, 3.8, 2.7, 2.4, 2.1, 1.7, 1.4,
         1.2, 0.9, 0.7, 0.5],
    100: [89.9, 87.2, 84.4, 81.6, 78.5, 75.7, 72.9, 70.1, 67.2, 64.4, 61.2, 58.4, 55.5, 52.7, 49.7,
          46.9, 44.0, 41.0, 38.2, 35.4, 32.8, 30.3, 27.7, 25.1, 22.7, 20.3, 18.3, 16.2, 14.4, 12.7,
          11.2, 9.8, 8.6, 7.6, 6.6, 5.9, 5.2, 4.7, 4.4, 3.8, 3.5, 3.3, 3.3, 3.0, 3.0, 2.9, 2.7,
          2.7, 2.7, 2.7, 2.7],
}
FIG3B = {C: (_E3 + HGO_TO_RHE, np.array(v, float)) for C, v in _J3.items()}

# Fig. 3a: Ni foam (NF) cathode, 100 mM NaNO3 + 1 M NaOH (same potential grid as above, k = 1..27)
_JNF = [12.60, 10.84, 9.52, 8.31, 7.43, 6.55, 5.89, 5.23, 5.01, 4.35, 4.02, 3.47, 3.36, 3.03,
        2.81, 2.48, 2.48, 2.15, 2.15, 1.82, 1.82, 1.60, 1.60, 1.49, 1.49, 1.27, 1.27]
FIG3A_NF = (_E3[1:28] + HGO_TO_RHE, np.array(_JNF))

# ---------------------------------------------------------------------------
# Fig. 11a: ANP anode LSV (potential already vs RHE), anodic current
#   UOR: 1 M NaOH + 0.5 M urea; OER: 1 M NaOH
# ---------------------------------------------------------------------------
_EU = 1.31 + 0.01 * np.arange(31)
_JU = [-0.3, 2.3, 5.9, 9.3, 13.4, 17.9, 22.4, 27.4, 32.4, 37.9, 43.4, 48.8, 54.7, 60.4, 66.5,
       72.6, 78.6, 84.9, 91.2, 97.6, 104.5, 110.6, 117.6, 124.1, 130.9, 138.1, 144.8, 152.2,
       159.2, 166.9, 174.7]
FIG11A_UOR = (_EU, np.array(_JU))
_EO = 1.40 + 0.01 * np.arange(22)
_JO = [-0.3, 0.2, 0.3, 0.8, 1.3, 1.7, 2.3, 2.8, 3.6, 4.5, 6.1, 7.6, 9.7, 12.6, 15.4, 18.9, 22.9,
       27.2, 31.8, 36.4, 41.7, 47.2]
FIG11A_OER = (_EO, np.array(_JO))

# ---------------------------------------------------------------------------
# Values reported directly in the text (quoted verbatim, for cross-checking)
# ---------------------------------------------------------------------------
REPORTED = {
    "UOR_E10_ANP": 1.335,       # V vs RHE @10 mA cm-2 (Sec. 3.4; Fig. 11b)
    "UOR_E10_NF": 1.423,
    "UOR_E10_Pt": 1.562,
    "UOR_b_ANP": 0.0442,        # V/dec (Fig. 11c)
    "UOR_b_NF": 0.0954,
    "OER_E10_ANP": 1.515,       # V vs RHE @10 mA cm-2
    "E_cath_ERN": -1.2 + HGO_TO_RHE,   # electrolysis operating potential -1.2 V vs Hg/HgO = -0.276 V vs RHE
    "I_avg_ERN_UOR": 19.0,      # mA (1 cm2 electrode), average current of the ERN + UOR system
    "Vcell_ratio_UOR_OER": 0.9441,
    "EC_saving_UOR": 0.175,     # abstract: 17.5 % lower energy consumption
    "NO3_N0_mg_L": 60.0,        # initial concentration for electrolysis, 60 mg/L as N (Fig. 3c)
}
# Fig. 14d readings (V): ERN+OER 7.14, ERN+UOR 6.76 (ratio 0.947, consistent with the 94.41 % in the text)
FIG14D_VCELL = {"OER": 7.14, "UOR": 6.76, "S(IV)OR": 6.24, "S(IV)/UOR": 6.16}
