# 1D Tidal SWE–Exner Model with Engelund–Hansen Sediment Transport

This repository contains a 1D finite-volume shallow-water morphodynamic solver for tidal overtopping of a coastal barrier.

## Features

- 1D shallow-water equations
- mobile bed via Exner equation
- HLLC approximate Riemann solver
- hydrostatic reconstruction for well-balanced topography treatment
- wetting and drying
- tidal stage forcing at the left boundary
- transmissive or wall boundary at the right boundary
- Shields-stress-based Engelund–Hansen total-load transport
- CSV outputs
- GIF animation

---

## Governing equations

### Hydrodynamics

\[
\partial_t h + \partial_x q = 0
\]

\[
\partial_t q + \partial_x \left(\frac{q^2}{h} + \frac12 g h^2 \right)
=
- g h \partial_x z_b + S_f
\]

where:
- \(h\) = water depth
- \(q = hu\) = unit discharge
- \(u\) = depth-averaged velocity
- \(z_b\) = bed elevation

### Friction

Manning friction is applied semi-implicitly:
\[
\tau_b = \rho_w g n^2 \frac{u |u|}{h^{1/3}}
\]

### Exner bed evolution

\[
(1-\lambda_p)\partial_t z_b + \partial_x q_s = 0
\]

### Sediment transport

A practical signed Engelund–Hansen transport form is used:
\[
q_s = \operatorname{sgn}(u)\,\alpha_{EH}\,\sqrt{(s-1) g d_{50}^3}\,\theta^{5/2}
\]

with optional Shields threshold:
\[
q_s = 0 \quad \text{if } \theta \le \theta_c
\]

and
\[
\theta = \frac{|\tau_b|}{(\rho_s-\rho_w) g d_{50}}
\]

---

## Numerical method

- finite volume
- MUSCL reconstruction of \(\eta\), \(u\), and \(z_b\)
- Audusse-style hydrostatic reconstruction
- HLLC flux
- SSP-RK2 time integration
- semi-implicit Manning friction
- explicit finite-volume Exner update
- upwinded face sediment flux

---

## Installation

```bash
pip install -r requirements.txt
```

---

## Run

```bash
python run_model.py
```

---

## Outputs

### CSV files
- `output/csv/timeseries.csv`
- `output/csv/profile_*.csv`

### GIF
- `output/gif/tidal_cycle.gif`

---

## Suggested workflow

1. Run with the default idealized barrier case.
2. Inspect `timeseries.csv` for flood/ebb reversal and barrier transport.
3. Inspect `profile_*.csv` for free surface, discharge, Shields stress, and bed change.
4. View `tidal_cycle.gif`.

---

## Notes

- Start with **fixed bed validation** if you are calibrating hydraulics.
- Then activate full morphodynamics.
- The Engelund–Hansen coefficient `alpha_eh` may need tuning for your application.
- If the inland boundary is closed, change `right_bc="wall"` in `cases/tidal_barrier_case.py`.
