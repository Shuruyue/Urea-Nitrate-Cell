"""Write the data block of the interactive page web/cell-explorer.html from the package.

The page computes the cell in the browser with a JavaScript copy of urea_nitrate_cell.model (checked by
scripts/verify_web_page.py). Everything else it shows comes from the package and is written here between the
markers of its data block: the scenario parameters, the standard potentials, the batch settings, the LSV data
digitized from Yao et al. (2025) and the calibration results (uncompensated resistance and background current).

Usage: python scripts/build_web_page.py
"""
import json
from pathlib import Path

import _path  # noqa: F401
import numpy as np

from urea_nitrate_cell import model as m
from urea_nitrate_cell.calibration import fit_anode, fit_cathode
from urea_nitrate_cell.yao2025 import FIG3B, FIG11A_OER, FIG11A_UOR

PAGE = Path(_path.ROOT) / "web" / "cell-explorer.html"
BEGIN = "/* ---- data: begin (written by scripts/build_web_page.py) */"
END = "/* ---- data: end */"


def clean(params):
    return {k: (bool(v) if isinstance(v, bool) else float(v)) for k, v in params.items()}


def pairs(E, J):
    return [[round(float(e), 4), round(float(j), 2)] for e, j in zip(E, J) if np.isfinite(j)]


def page_data():
    anode = fit_anode()
    cathode = fit_cathode(gamma=1.0)
    scenarios = {k: clean(p) for k, p in m.SCEN.items()}
    scenarios["OER"] = clean(m.OER_ONLY)
    labels = dict(m.SCEN_LABEL, OER="Control: S1 with an OER-only anode")
    return {
        "scenarios": scenarios,
        "labels": labels,
        "E0": {k: float(v) for k, v in m.E0.items()},
        "batch": {k: float(v) for k, v in m.BATCH.items()},
        "yao": {
            "cathode": {str(C): pairs(*FIG3B[C]) for C in sorted(FIG3B)},
            "anode_uor": pairs(*FIG11A_UOR),
            "anode_oer": pairs(*FIG11A_OER),
        },
        "fit": {
            "anode_Ru": float(anode["used"]["Ru"]),
            "anode_rms": float(anode["used"]["rms"]),
            "cathode_Ru": float(cathode["params"]["Ru"]),
            "cathode_j_bg": float(cathode["params"]["j_bg"]),
            "cathode_rms_by_C": {str(C): float(v) for C, v in cathode["rms_by_C"].items()},
        },
    }


def main():
    page = PAGE.read_text(encoding="utf-8")
    head, rest = page.split(BEGIN)
    _, tail = rest.split(END)
    block = "\nconst DATA = " + json.dumps(page_data(), ensure_ascii=False, separators=(",", ":")) + ";\n"
    PAGE.write_text(head + BEGIN + block + END + tail, encoding="utf-8", newline="\n")
    print(f"wrote the data block of {PAGE}")


if __name__ == "__main__":
    main()
