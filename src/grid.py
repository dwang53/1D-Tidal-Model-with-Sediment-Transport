import numpy as np

def build_grid(x_min, x_max, nx):
    dx = (x_max - x_min) / nx
    x = np.linspace(x_min + 0.5 * dx, x_max - 0.5 * dx, nx)
    return x, dx
