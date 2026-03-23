# 1D Tidal SWE–Exner Model with Engelund–Hansen Sediment Transport

This repository contains a 1D finite-volume shallow-water morphodynamic solver for tidal overtopping of a coastal barrier. The tidal wave arrives from the left side of the domain, overtops a coastal barrier, and floods the back-barrier basin. The tide also retreats, reversing the flow direction. Both the flood and ebb phases drive sediment transport and evolve the bed via the Exner equation.

## Features

- 1D shallow-water equations (SWE) on a uniform finite-volume grid
- mobile bed via Exner equation coupled to the SWE
- HLLC approximate Riemann solver for robust shock/rarefaction capturing
- MUSCL reconstruction with generalized-minmod limiter (second-order spatial accuracy)
- Audusse-style hydrostatic reconstruction for a well-balanced, positivity-preserving scheme
- SSP-RK2 (Heun's method) time integration (second-order in time)
- wetting and drying with a dry-cell threshold
- tidal stage forcing at the left boundary via a characteristic-based ghost cell
- transmissive or wall boundary at the right boundary
- semi-implicit Manning friction (unconditionally stable friction term)
- Shields-stress-based Engelund–Hansen total-load sediment transport
- CSV outputs for every profile snapshot and a scalar timeseries
- GIF animation of the tidal cycle

---

## Recommended numerical scheme and boundary conditions

### Why HLLC + MUSCL + SSP-RK2 + hydrostatic reconstruction?

For 1D tidal flow over a barrier with wetting/drying:

- **HLLC Riemann solver**: correctly handles shocks (bore-like overtopping fronts), rarefactions (ebb recession), and dry-front propagation without an entropy fix.
- **MUSCL reconstruction**: upgrades to second-order spatial accuracy. The generalized-minmod limiter prevents spurious oscillations near the crest and wet/dry interface.
- **Audusse hydrostatic reconstruction**: ensures the scheme is *well-balanced* (lake-at-rest is preserved exactly over non-flat bathymetry) and *positivity-preserving* (depth never becomes negative during wetting/drying).
- **SSP-RK2 time integration**: strong-stability-preserving two-stage Heun method keeps the TVD property in time and is second-order accurate.
- **Semi-implicit Manning friction**: removes the stiff CFL constraint from the friction term and prevents velocity blow-up in very shallow water.
- **Explicit Exner step**: bed change per time step is small relative to the flow timescale, so an explicit upwind Exner update is stable and accurate.

### Boundary conditions

| Boundary | Condition | Rationale |
|---|---|---|
| Left (offshore) | **Tidal stage** (prescribed η, characteristic velocity) | Forces the prescribed tidal elevation; the ghost-cell velocity is derived from the outgoing Riemann invariant R⁻ = u − 2c, which allows ebb waves to leave freely and prevents spurious reflections |
| Right (inland) | **Transmissive** (zero-gradient extrapolation) | Allows flood water to drain inland without reflection; change to `"wall"` for a closed inland boundary |

---

## Governing equations

### Hydrodynamics

$$
\partial_t h + \partial_x q = 0
$$

$$
\partial_t q + \partial_x \left(\frac{q^2}{h} + \frac12 g h^2 \right)
=
- g h \partial_x z_b + S_f
$$

where:
- $h$ = water depth
- $q = hu$ = unit discharge
- $u$ = depth-averaged velocity
- $z_b$ = bed elevation

### Friction

Manning friction is applied semi-implicitly:
$$
\tau_b = \rho_w g n^2 \frac{u |u|}{h^{1/3}}
$$

### Exner bed evolution

$$
(1-\lambda_p)\partial_t z_b + \partial_x q_s = 0
$$

### Sediment transport

A practical signed Engelund–Hansen transport form is used:
$$
q_s = \operatorname{sgn}(u)\,\alpha_{EH}\,\sqrt{(s-1) g d_{50}^3}\,\theta^{5/2}
$$

with optional Shields threshold:
$$
q_s = 0 \quad \text{if } \theta \le \theta_c
$$

and
$$
\theta = \frac{|\tau_b|}{(\rho_s-\rho_w) g d_{50}}
$$

---

## Numerical method

- finite volume on a uniform 1D grid
- MUSCL reconstruction of $\eta$, $u$, and $z_b$ with generalized-minmod limiter
- Audusse-style hydrostatic reconstruction (well-balanced, positivity-preserving)
- HLLC approximate Riemann flux
- SSP-RK2 (Heun's method) time integration
- characteristic-based ghost cell at the left (tidal) boundary using the outgoing Riemann invariant
- zero-gradient (transmissive) ghost cell at the right boundary
- semi-implicit Manning friction (stable for any dt)
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

Outputs are written to `output/` (created automatically).

### CSV files
- `output/csv/<case_name>/timeseries.csv` — time series of tidal stage, barrier depth, discharge, inundation extent, and bed change extremes
- `output/csv/<case_name>/profile_*.csv` — spatial profiles (x, zb, h, eta, u, q, qs, theta, tau_b) at each output time

### GIF
- `output/gif/<case_name>.gif` — animated tidal cycle showing free surface, discharge, and Shields parameter

---

## Suggested workflow

1. Run with the default idealized barrier case: `python run_model.py`
2. Inspect `timeseries.csv` for flood/ebb reversal and barrier transport.
3. Inspect `profile_*.csv` for free surface, discharge, Shields stress, and bed change.
4. View `tidal_barrier_morphodynamic.gif` in `output/gif/`.

---

## Notes

- Start with **fixed bed validation** if you are calibrating hydraulics.
- Then activate full morphodynamics.
- The Engelund–Hansen coefficient `alpha_eh` may need tuning for your application.
- If the inland boundary is closed, change `right_bc="wall"` in `cases/tidal_barrier_case.py`.


## Validation suite

A fixed-bed validation suite is included to verify the hydrodynamic core before morphodynamic runs.

### Cases

1. **Lake at rest over nonflat topography**
   - verifies well-balancing
   - expected result: no spurious velocity and constant free surface

2. **Dam-break over flat bed**
   - verifies shock capturing
   - expected result: stable propagation of rarefaction and bore

3. **Wetting/drying on a sloping beach**
   - verifies shoreline motion and inundation/recession

4. **Fixed-bed tidal barrier case**
   - verifies tidal forcing, overtopping, and reversal without bed change

### Run validation

```bash
python run_validation.py
```

### Validation outputs

Each case writes outputs to:
- `output/csv/<case_name>/`
- `output/gif/<case_name>.gif`
