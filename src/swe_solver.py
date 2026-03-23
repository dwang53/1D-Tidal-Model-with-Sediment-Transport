import numpy as np

from src.grid import build_grid
from src.bathymetry import initial_bathymetry
from src.state import initial_state, free_surface, velocity
from src.time_integration import compute_dt
from src.swe_operator import hydro_rhs
from src.friction import apply_semiimplicit_manning
from src.exner import exner_step
from src.outputs import ensure_output_dirs, write_profile_csv, write_timeseries_csv
from src.visualization import save_frame, build_gif
from src.boundary_conditions import tidal_stage

def rk2_hydro_step(h, q, zb, t, dt, config, dx):
    rhs_h1, rhs_q1 = hydro_rhs(h, q, zb, t, config, dx)

    h1 = np.maximum(h + dt * rhs_h1, 0.0)
    q1 = q + dt * rhs_q1
    q1 = apply_semiimplicit_manning(h1, q1, config.g, config.manning_n, config.h_dry, dt)
    q1[h1 <= config.h_dry] = 0.0

    rhs_h2, rhs_q2 = hydro_rhs(h1, q1, zb, t + dt, config, dx)

    h2 = np.maximum(0.5 * h + 0.5 * (h1 + dt * rhs_h2), 0.0)
    q2 = 0.5 * q + 0.5 * (q1 + dt * rhs_q2)
    q2 = apply_semiimplicit_manning(h2, q2, config.g, config.manning_n, config.h_dry, dt)
    q2[h2 <= config.h_dry] = 0.0

    return h2, q2

def run_simulation(config):
    ensure_output_dirs()

    x, dx = build_grid(config.x_min, config.x_max, config.nx)
    zb0 = initial_bathymetry(x)
    zb = zb0.copy()

    h, q = initial_state(x, zb, config.eta0, config.h_dry)

    t = 0.0
    step = 0
    next_output = 0.0

    barrier_idx = int(np.argmin(np.abs(x - 500.0)))
    records = []

    while t < config.t_end:
        dt = compute_dt(h, q, config.g, dx, config.cfl, config.dt_max, config.h_dry)
        if t + dt > config.t_end:
            dt = config.t_end - t

        # Hydrodynamic step
        h, q = rk2_hydro_step(h, q, zb, t, dt, config, dx)

        # Morphodynamic step
        zb, diag = exner_step(zb, h, q, dt, dx, config)

        # Wet/dry cleanup
        h[h <= config.h_dry] = 0.0
        q[h <= config.h_dry] = 0.0

        t += dt
        step += 1

        if t >= next_output - 1.0e-12:
            eta = free_surface(h, zb)
            u = velocity(h, q, config.h_dry)

            qs = diag["qs_cell"]
            theta = diag["theta"]
            tau_b = diag["tau_b"]

            write_profile_csv(step, t, x, zb, h, eta, u, q, qs, theta, tau_b)
            save_frame(step, t, x, zb, h, eta, q, theta)

            eta_left = tidal_stage(t, config.eta_mean, config.eta_amp, config.eta_period)
            inundation_extent = np.sum(h > config.h_dry) * dx

            records.append({
                "time": t,
                "eta_left": eta_left,
                "barrier_depth": h[barrier_idx],
                "barrier_discharge": q[barrier_idx],
                "barrier_velocity": u[barrier_idx],
                "barrier_qs": qs[barrier_idx],
                "barrier_theta": theta[barrier_idx],
                "inundation_extent": inundation_extent,
                "max_eta": np.max(eta),
                "min_eta": np.min(eta),
                "max_bed_change": np.max(zb - zb0),
                "min_bed_change": np.min(zb - zb0),
            })

            next_output += config.output_interval

    write_timeseries_csv(records)
    build_gif()

    print("Simulation complete.")
    print("Generated files:")
    print("  output/csv/timeseries.csv")
    print("  output/csv/profile_*.csv")
    print("  output/gif/tidal_cycle.gif")
