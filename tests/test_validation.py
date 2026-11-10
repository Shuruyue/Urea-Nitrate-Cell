"""Full validation suite (36 checks + 11 calibration metrics, includes the parameter fits; ~1–2 min)."""
import math

from urea_nitrate_cell.validation import run_all


def test_all_validation_checks_pass():
    rows = run_all()
    checks = [r for r in rows if r["kind"] == "check"]
    metrics = [r for r in rows if r["kind"] == "metric"]
    assert len(checks) == 36 and len(metrics) == 11
    failed = [r["name"] for r in checks if not r["ok"]]
    assert not failed, failed
    # metrics carry no threshold, but each must still be computed
    assert all(r["ok"] is None and r["tol"] is None and math.isfinite(r["value"]) for r in metrics)
