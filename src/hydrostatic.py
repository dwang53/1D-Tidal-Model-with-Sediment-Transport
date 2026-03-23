import numpy as np

def safe_velocity(h, q, h_dry):
    u = np.zeros_like(h)
    wet = h > h_dry
    u[wet] = q[wet] / h[wet]
    return u

def hydrostatic_reconstruction(etaL, uL, zbL, etaR, uR, zbR, h_dry):
    """
    Audusse-style hydrostatic reconstruction.

    Reconstruct eta, u, zb first, then clip depths at interface topography.
    """
    zint = np.maximum(zbL, zbR)

    hL = np.maximum(0.0, etaL - zint)
    hR = np.maximum(0.0, etaR - zint)

    qL = np.zeros_like(hL)
    qR = np.zeros_like(hR)

    wetL = hL > h_dry
    wetR = hR > h_dry

    qL[wetL] = hL[wetL] * uL[wetL]
    qR[wetR] = hR[wetR] * uR[wetR]

    return hL, qL, hR, qR
