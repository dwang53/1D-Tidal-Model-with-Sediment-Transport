import numpy as np
from src.sediment import face_upwind_transport

def exner_step(zb, h, q, dt, dx, config):
    """
    Explicit finite-volume Exner update with upwind face sediment transport.
    """
    n = len(zb)

    hg = np.zeros(n + 2)
    qg = np.zeros(n + 2)

    hg[1:-1] = h
    qg[1:-1] = q
    hg[0] = hg[1]
    hg[-1] = hg[-2]
    qg[0] = qg[1]
    qg[-1] = qg[-2]

    qs_cell_g, qs_face_g, theta_g, tau_b_g, u_g = face_upwind_transport(
        h=hg,
        q=qg,
        rho_w=config.rho_w,
        rho_s=config.rho_s,
        d50=config.d50,
        g=config.g,
        n_manning=config.manning_n,
        h_dry=config.h_dry,
        alpha_eh=config.alpha_eh,
        theta_cr=config.theta_cr,
        use_threshold=config.use_theta_threshold
    )

    zb_new = zb.copy()
    factor = dt / ((1.0 - config.porosity) * dx)

    for i in range(n):
        qsl = qs_face_g[i]
        qsr = qs_face_g[i + 1]
        zb_new[i] = zb[i] - factor * (qsr - qsl)

    diagnostics = {
        "qs_cell": qs_cell_g[1:-1],
        "theta": theta_g[1:-1],
        "tau_b": tau_b_g[1:-1],
        "u": u_g[1:-1],
        "qs_face_ghost": qs_face_g
    }

    return zb_new, diagnostics
