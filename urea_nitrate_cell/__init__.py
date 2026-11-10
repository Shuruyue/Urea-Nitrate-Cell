"""Urea–nitrate direct cell (UOR ‖ NO3RR): thermodynamics, kinetics and mass-transport model.

Modules
-------
model        thermodynamics, single-electrode kinetics, full-cell solver, sensitivity, batch simulation
yao2025      data digitized from the figures of Yao et al., Sustainability 17, 3835 (2025) (CC BY 4.0)
calibration  fit of the model parameters to the yao2025 LSV data
validation   36 automated checks + 11 calibration metrics (run_all)
summary      main results (compute)
plotting     figures
"""
__version__ = "0.1.0"
