import numpy as np

def swe_flux(h, q, g, h_dry):
    f1 = np.zeros_like(h)
    f2 = np.zeros_like(h)

    wet = h > h_dry
    f1[wet] = q[wet]
    f2[wet] = q[wet] ** 2 / h[wet] + 0.5 * g * h[wet] ** 2

    return f1, f2

def hllc_flux(hL, qL, hR, qR, g, h_dry):
    """
    HLLC flux for 1D shallow-water equations.
    """
    n = len(hL)
    f1 = np.zeros(n)
    f2 = np.zeros(n)

    uL = np.zeros(n)
    uR = np.zeros(n)

    wetL = hL > h_dry
    wetR = hR > h_dry

    uL[wetL] = qL[wetL] / hL[wetL]
    uR[wetR] = qR[wetR] / hR[wetR]

    cL = np.sqrt(g * np.maximum(hL, 0.0))
    cR = np.sqrt(g * np.maximum(hR, 0.0))

    FL1, FL2 = swe_flux(hL, qL, g, h_dry)
    FR1, FR2 = swe_flux(hR, qR, g, h_dry)

    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)

    pL = 0.5 * g * hL**2
    pR = 0.5 * g * hR**2

    denom = hL * (SL - uL) - hR * (SR - uR)
    Sstar = np.zeros(n)

    valid = np.abs(denom) > 1.0e-14
    Sstar[valid] = (
        pR[valid] - pL[valid]
        + qL[valid] * (SL[valid] - uL[valid])
        - qR[valid] * (SR[valid] - uR[valid])
    ) / denom[valid]

    hstarL = np.zeros(n)
    hstarR = np.zeros(n)

    denomL = SL - Sstar
    denomR = SR - Sstar

    okL = np.abs(denomL) > 1.0e-14
    okR = np.abs(denomR) > 1.0e-14

    hstarL[okL] = hL[okL] * (SL[okL] - uL[okL]) / denomL[okL]
    hstarR[okR] = hR[okR] * (SR[okR] - uR[okR]) / denomR[okR]

    qstarL = hstarL * Sstar
    qstarR = hstarR * Sstar

    left = 0.0 <= SL
    midL = (SL <= 0.0) & (0.0 <= Sstar)
    midR = (Sstar <= 0.0) & (0.0 <= SR)
    right = SR <= 0.0

    f1[left] = FL1[left]
    f2[left] = FL2[left]

    f1[midL] = FL1[midL] + SL[midL] * (hstarL[midL] - hL[midL])
    f2[midL] = FL2[midL] + SL[midL] * (qstarL[midL] - qL[midL])

    f1[midR] = FR1[midR] + SR[midR] * (hstarR[midR] - hR[midR])
    f2[midR] = FR2[midR] + SR[midR] * (qstarR[midR] - qR[midR])

    f1[right] = FR1[right]
    f2[right] = FR2[right]

    drydry = (hL <= h_dry) & (hR <= h_dry)
    f1[drydry] = 0.0
    f2[drydry] = 0.0

    return f1, f2
