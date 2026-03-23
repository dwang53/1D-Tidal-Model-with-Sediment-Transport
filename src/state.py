import numpy as np

def initial_state(x, zb, eta0, h_dry):
    eta = eta0 * np.ones_like(x)
    h = np.maximum(eta - zb, 0.0)
    q = np.zeros_like(x)
    q[h <= h_dry] = 0.0
    return h, q

def free_surface(h, zb):
    return h + zb

def velocity(h, q, h_dry):
    u = np.zeros_like(h)
    wet = h > h_dry
    u[wet] = q[wet] / h[wet]
    return u
