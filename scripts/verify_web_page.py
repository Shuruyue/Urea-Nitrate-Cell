"""Check the JavaScript model of the interactive page web/cell-explorer.html against the Python package.

The page holds a JavaScript copy of urea_nitrate_cell/model.py and of the measured-curve functions of
calibration.py between the markers "// ---- physics: begin" and "// ---- physics: end". This script runs that
block in Node.js and compares it with the package for random parameter sets: cell voltage, electrode potentials
and Faradaic efficiencies (solve_cell), limiting currents, the currents of the measured LSV curves (anode_meas,
cathode_meas) and complete batch runs (simulate_batch).

Usage: python scripts/verify_web_page.py [--cases 1500]
Requires Node.js on the PATH.
"""
import argparse
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

import _path  # noqa: F401

from urea_nitrate_cell import model as m
from urea_nitrate_cell.calibration import anode_meas, cathode_meas

PAGE = Path(_path.ROOT) / "web" / "cell-explorer.html"
ABS_TOL, REL_TOL = 1e-9, 1e-8

NODE_HARNESS = r"""
const fs = require("fs"), vm = require("vm");
const page = fs.readFileSync(process.argv[1], "utf8");
const after = page.split("// ---- physics: begin")[1], block = after.slice(after.indexOf("\n")).split("// ---- physics: end")[0];
const ctx = {}; vm.createContext(ctx);
vm.runInContext(block + "\n;this.api = { solveCell, jLimNo3, anodeMeas, cathodeMeas, simulateBatch };", ctx);
const api = ctx.api, input = JSON.parse(fs.readFileSync(0, "utf8")), out = { cells: [], meas: [], batches: [] };
for (const c of input.cells) { const r = api.solveCell(c.j, c.p, c.C); out.cells.push([r.V, r.Ea, r.Ec, r.FE_UOR, r.FE_NO3, api.jLimNo3(c.C, c.p.delta_um)]); }
for (const c of input.meas) out.meas.push(c.kind === "anode" ? api.anodeMeas(c.E, c.p, c.Ru, c.uor) : api.cathodeMeas(c.E, c.C, c.p, c.Ru, c.jbg));
for (const c of input.batches) { const r = api.simulateBatch(c.p, c.strategy, c.j_set, 0.8, c.b); out.batches.push(r && [r.t_end, r.FE_avg, r.E_kWh, r.kWh_per_kgN_no3, r.kWh_per_kgN_all, r.C[r.C.length - 1]]); }
process.stdout.write(JSON.stringify(out, (k, v) => (typeof v === "number" && !isFinite(v)) ? null : v));
"""

RANGES = dict(E10_UOR=(0.45, 1.45), b_UOR=(0.025, 0.12), E10_OER=(1.42, 1.75), b_OER=(0.04, 0.12),
              E10_NO3=(-0.25, 0.6), b_NO3=(0.04, 0.18), gamma=(0.3, 1.4), E10_HER=(-0.45, -0.05), b_HER=(0.04, 0.18),
              delta_um=(8.0, 250.0), R_ohm=(0.0, 8.0))


def random_params(rng):
    p = dict(rng.choice([m.S1, m.SCEN["S2"], m.SCEN["S3"]]))
    for key in rng.sample(sorted(RANGES), 4):
        p[key] = rng.uniform(*RANGES[key])
    p["E0_UOR"] = rng.choice([m.E0["UOR_lit"], m.E0["UOR_calc"], m.E0["UOR_CO3"]])
    if rng.random() < 0.15:
        p["UOR_on"] = False
    return p


def close(a, b):
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= ABS_TOL + REL_TOL * abs(b)


def run(n_cases=1500, seed=0):
    node = shutil.which("node")
    if node is None:
        raise RuntimeError("Node.js is needed to run the page's JavaScript")
    rng = random.Random(seed)
    cells = [dict(p=random_params(rng), j=10 ** rng.uniform(-0.3, 2.3), C=10 ** rng.uniform(-0.3, 3)) for _ in range(n_cases)]
    meas = []
    for _ in range(n_cases // 3):
        p = random_params(rng)
        meas.append(dict(kind="anode", p=p, E=rng.uniform(1.28, 1.6), Ru=rng.uniform(0, 3), uor=rng.random() < 0.7))
        meas.append(dict(kind="cathode", p=p, E=rng.uniform(-0.3, 0.35), C=rng.choice([0, 1, 10, 50, 100]),
                         Ru=rng.uniform(0, 3), jbg=rng.uniform(0, 3)))
    batches = [dict(p=dict(sc), strategy=st, j_set=10.0, b=m.BATCH) for sc in (m.S1, m.SCEN["S2"], m.SCEN["S3"]) for st in ("cc", "follow")]
    batches += [dict(p=random_params(rng), strategy=rng.choice(["cc", "follow"]), j_set=rng.uniform(2, 20), b=m.BATCH) for _ in range(6)]
    payload = json.dumps(dict(cells=cells, meas=meas, batches=batches))
    out = subprocess.run([node, "-e", NODE_HARNESS, str(PAGE)], input=payload, capture_output=True, text=True,
                         check=False)
    if out.returncode != 0:
        raise RuntimeError("Node.js failed:\n" + out.stderr)
    got = json.loads(out.stdout)

    report = {"cells": [0, 0], "measured curves": [0, 0], "batches": [0, 0]}
    unsolved = {"cells": 0, "measured curves": 0, "batches": 0}
    for c, js in zip(cells, got["cells"]):
        try:
            r = m.solve_cell(c["j"], c["p"], c["C"])
            ref = [r["V"], r["Ea"], r["Ec"], r["FE_UOR"], r["FE_NO3"], m.j_lim_no3(c["C"], c["p"]["delta_um"])]
        except ValueError:          # no root in the bracket: the page returns no number
            ref = [None] * 5 + [m.j_lim_no3(c["C"], c["p"]["delta_um"])]
        report["cells"][0] += 1
        unsolved["cells"] += ref[0] is None
        report["cells"][1] += not all(close(a, b) for a, b in zip(js, ref))
    for c, js in zip(meas, got["meas"]):
        q = dict(c["p"], Ru=c["Ru"], j_bg=c.get("jbg", 0.0))
        try:
            ref = anode_meas(c["E"], c["p"], c["Ru"], c["uor"]) if c["kind"] == "anode" else cathode_meas(c["E"], c["C"], q)
        except ValueError:          # no root in the bracket: the page returns no number
            ref = None
        report["measured curves"][0] += 1
        unsolved["measured curves"] += ref is None
        report["measured curves"][1] += not close(js, ref)
    for c, js in zip(batches, got["batches"]):
        try:
            r = m.simulate_batch(c["p"], c["strategy"], j_set=c["j_set"])
            ref = [r["t_end"], r["FE_avg"], r["E_kWh"], r["kWh_per_kgN_no3"], r["kWh_per_kgN_all"], float(r["C"][-1])]
        except ValueError:
            ref = None
        report["batches"][0] += 1
        unsolved["batches"] += ref is None
        report["batches"][1] += not (ref is None and js is None or (js is not None and ref is not None and all(close(a, b) for a, b in zip(js, ref))))
    return report, unsolved


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cases", type=int, default=1500)
    args = parser.parse_args()
    report, unsolved = run(args.cases)
    for key, (n, bad) in report.items():
        print(f"{key:16s} {n:5d} compared ({unsolved[key]} without a solution in both), {bad} differ by more than "
              f"{ABS_TOL:g} + {REL_TOL:g} relative")
    ok = all(bad == 0 for _, bad in report.values())
    print("OK" if ok else "MISMATCH")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
