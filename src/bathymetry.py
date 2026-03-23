import numpy as np

def barrier_bathymetry(x):
    """
    Simple idealized cross-shore 1D bathymetry:
    - offshore region
    - coastal barrier crest
    - shallow back-barrier basin
    """
    zb = -0.75 + 0.0006 * x

    crest_center = 500.0
    crest_halfwidth = 70.0
    crest_height = 1.35

    shape = np.maximum(0.0, 1.0 - np.abs(x - crest_center) / crest_halfwidth)
    zb += crest_height * shape

    # back-barrier basin
    zb += -0.25 * (x >= 650.0)
    return zb

def initial_bathymetry(x):
    return barrier_bathymetry(x)
