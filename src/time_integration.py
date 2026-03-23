import numpy as np

def compute_dt(h, q, g, dx, cfl, dt_max, h_dry):
    u = np.zeros_like(h)
    wet = h > h_dry
    u[wet] = q[wet] / h[wet]

    max_speed = np.max(np.abs(u) + np.sqrt(g * np.maximum(h, 0.0)))
    if max_speed < 1.0e-14:
        return dt_max

    return min(dt_max, cfl * dx / max_speed)
