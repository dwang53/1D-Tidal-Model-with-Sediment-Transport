import numpy as np

def minmod3(a, b, c):
    out = np.zeros_like(a)
    same_sign = (np.sign(a) == np.sign(b)) & (np.sign(b) == np.sign(c))
    out[same_sign] = np.sign(a[same_sign]) * np.minimum(
        np.abs(a[same_sign]),
        np.minimum(np.abs(b[same_sign]), np.abs(c[same_sign]))
    )
    return out

def limited_slope(w, theta=1.5):
    """
    Generalized minmod slope.
    Returns limited increment (not divided by dx).
    """
    slope = np.zeros_like(w)
    if len(w) < 3:
        return slope

    dm = w[1:-1] - w[:-2]
    dc = 0.5 * (w[2:] - w[:-2])
    dp = w[2:] - w[1:-1]

    slope[1:-1] = minmod3(theta * dm, dc, theta * dp)
    return slope
