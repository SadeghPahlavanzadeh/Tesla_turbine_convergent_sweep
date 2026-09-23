# Tesla_turbine_convergent_sweep
# ANSYS PyFluent Parametric CFD Automation: Tesla Turbine

An advanced, production-grade Python framework built on **ANSYS PyFluent** for executing automated parametric CFD sweeps on a Tesla turbine. This repository features robust boundary layer adaptation and **multi-criteria physical and numerical convergence monitoring**.

---

## Key Features & Improvements (v2.0)

Unlike standard fixed-iteration automation scripts, this framework implements rigorous convergence and quality assurance checks during runtime:

1. **Mass Balance Conservation:** Continuously monitors the relative mass flow rate imbalance between the inlet and outlet boundaries to ensure steady-state continuity.
2. **Torque Generation Stability:** Evaluates the running variance of the aerodynamic torque (integrated shear stress and pressure moments on the co-rotating disks) across a rolling window.
3. **Local Minimum Pressure Fluctuation Tracker:** Following a robust initial solver phase (5,000 iterations), the script isolates the spatial zone of minimum pressure and tracks its localized stability until the pressure delta drops below **100 Pa**.
4. **Adaptive Boundary Layer Refinement:** Automatically monitors $y^+$ on heated/rotating walls and executes dynamic prism-layer splitting when $y^+ > 1.0$ or when critical low-pressure expansion thresholds are crossed.

---

## Parametric Matrix

* **Inlet Pressure Sweep:** 2.0 bar to 3.0 bar absolute (0.1 bar increments, 11 levels).
* **Rotational Speed Sweep:** 12,000 RPM to 25,000 RPM (500 RPM increments, 27 levels).
* **Total Design Points:** 297 automated runs.

---

## Performance Metrics Extracted

* Euler Power Output [W]
* Aerodynamic Torque [N·m]
* Inlet Mass Flow Rate [kg/s]
* Total-to-Static Isentropic Efficiency [%]
* Final Resolved y+ and Minimum Gauge Pressure

---

## Requirements

* ANSYS Fluent (with Python/PyFluent support enabled)
* Python 3.x
* NumPy
* Pandas

---

## Usage

Ensure your baseline case file (`tesla_turbine.cas.h5`) is in the working directory, then run:

```bash
python tesla_turbine_convergent_sweep.py
