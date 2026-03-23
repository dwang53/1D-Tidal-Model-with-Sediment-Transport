import os
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

    case_dir_csv = os.path.join("output", "csv", config.case_name)
    case_dir_frames = os.path.join("output", "frames", config.case_name)
    case_dir_gif = os.path.join("output", "gif")
    os.makedirs(case_dir_csv, exist_ok=True)
    os.makedirs(case_dir_frames, exist_ok=True)
    os.makedirs(case_dir_gif, exist_ok=True)

    while t < config.t_end:
        dt = compute_dt(h, q, config.g, dx, config.cfl, config.dt_max, config.h_dry)
        if t + dt > config.t_end:
            dt = config.t_end - t

        h, q = rk2_hydro_step(h, q, zb, t, dt, config, dx)

        if config.enable_morphodynamics:
            zb, diag = exner_step(zb, h, q, dt, dx, config)
            qs = diag["qs_cell"]
            theta = diag["theta"]
            tau_b = diag["tau_b"]
        else:
            u_tmp = velocity(h, q, config.h_dry)
            qs = np.zeros_like(h)
            theta = np.zeros_like(h)
            tau_b = np.zeros_like(h)

        h[h <= config.h_dry] = 0.0
        q[h <= config.h_dry] = 0.0

        t += dt
        step += 1

        if t >= next_output - 1.0e-12:
            eta = free_surface(h, zb)
            u = velocity(h, q, config.h_dry)

            import pandas as pd
            df = pd.DataFrame({
                "x": x,
                "zb": zb,
                "h": h,
                "eta": eta,
                "u": u,
                "q": q,
                "qs": qs,
                "theta": theta,
                "tau_b": tau_b
            })
            df.to_csv(os.path.join(case_dir_csv, f"profile_{step:05d}.csv"), index=False)

            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

            axes[0].plot(x, zb, color="saddlebrown", linewidth=2, label="bed")
            axes[0].plot(x, eta, color="royalblue", linewidth=2, label="free surface")
            axes[0].fill_between(x, zb, eta, where=(h > 1e-8), color="lightskyblue", alpha=0.5)
            axes[0].set_ylabel("Elevation [m]")
            axes[0].legend(loc="best")
            axes[0].set_title(f"{config.case_name}, t = {t:.2f} s")

            axes[1].plot(x, q, color="darkred", linewidth=1.5, label="discharge q")
            axes[1].axhline(0.0, color="k", linewidth=0.75)
            axes[1].set_ylabel("q")
            axes[1].legend(loc="best")

            axes[2].plot(x, theta, color="darkgreen", linewidth=1.5, label="theta")
            axes[2].axhline(0.0, color="k", linewidth=0.75)
            axes[2].set_ylabel("theta")
            axes[2].set_xlabel("x [m]")
            axes[2].legend(loc="best")

            fig.tight_layout()
            frame_path = os.path.join(case_dir_frames, f"frame_{step:05d}.png")
            fig.savefig(frame_path, dpi=120)
            plt.close(fig)

            eta_left = tidal_stage(t, config.eta_mean, config.eta_amp, config.eta_period) if config.left_bc == "tide" else np.nan
            inundation_extent = np.sum(h > config.h_dry) * dx

            records.append({
                "time": t,
                "eta_left": eta_left,
                "barrier_depth": h[barrier_idx] if barrier_idx < len(h) else np.nan,
                "barrier_discharge": q[barrier_idx] if barrier_idx < len(q) else np.nan,
                "barrier_velocity": u[barrier_idx] if barrier_idx < len(u) else np.nan,
                "inundation_extent": inundation_extent,
                "max_eta": np.max(eta),
                "min_eta": np.min(eta),
                "max_bed_change": np.max(zb - zb0),
                "min_bed_change": np.min(zb - zb0),
            })

            next_output += config.output_interval

    import pandas as pd
    pd.DataFrame(records).to_csv(os.path.join(case_dir_csv, "timeseries.csv"), index=False)

    import imageio.v2 as imageio
    frames = sorted(
        os.path.join(case_dir_frames, f)
        for f in os.listdir(case_dir_frames)
        if f.endswith(".png")
    )
    if frames:
        images = [imageio.imread(f) for f in frames]
        imageio.mimsave(os.path.join(case_dir_gif, f"{config.case_name}.gif"), images, duration=0.15)

    print(f"Simulation complete for case: {config.case_name}")
    print(f"CSV output: output/csv/{config.case_name}/")
    print(f"GIF output: output/gif/{config.case_name}.gif")
