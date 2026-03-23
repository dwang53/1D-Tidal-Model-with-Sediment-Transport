import numpy as np
from src.friction import bed_shear_stress_manning

def cell_velocity(h, q, h_dry):
    u = np.zeros_like(h)
    wet = h > h_dry
    u[wet] = q[wet] / h[wet]
    return u

def shields_parameter(h, q, rho_w, rho_s, d50, g, n_manning, h_dry):
    """
    Shields parameter based on Manning-derived bed shear.
    """
    u = cell_velocity(h, q, h_dry)
    tau_b = bed_shear_stress_manning(h, u, rho_w, g, n_manning, h_dry)

    denom = (rho_s - rho_w) * g * d50
    theta = np.zeros_like(h)
    if denom > 0.0:
        theta = np.abs(tau_b) / denom
    return theta, tau_b, u

def engelund_hansen_transport(
    h,
    q,
    rho_w,
    rho_s,
    d50,
    g,
    n_manning,
    h_dry,
    alpha_eh,
    theta_cr=0.0,
    use_threshold=False
):
    """
    Signed Engelund-Hansen-style total load transport.

    Practical implementation:
      qs = sign(u) * alpha_eh * sqrt((s-1) g d50^3) * theta^(5/2)

    where:
      s = rho_s / rho_w
      theta = Shields parameter
    """
    theta, tau_b, u = shields_parameter(
        h=h,
        q=q,
        rho_w=rho_w,
        rho_s=rho_s,
        d50=d50,
        g=g,
        n_manning=n_manning,
        h_dry=h_dry
    )

    s_rel = rho_s / rho_w
    scale = np.sqrt((s_rel - 1.0) * g * d50**3)

    theta_eff = theta.copy()
    if use_threshold:
        theta_eff = np.maximum(theta - theta_cr, 0.0)

    qs = np.sign(u) * alpha_eh * scale * theta_eff ** 2.5

    dry = h <= h_dry
    qs[dry] = 0.0
    theta[dry] = 0.0
    tau_b[dry] = 0.0
    u[dry] = 0.0

    return qs, theta, tau_b, u

def face_upwind_transport(
    h,
    q,
    rho_w,
    rho_s,
    d50,
    g,
    n_manning,
    h_dry,
    alpha_eh,
    theta_cr,
    use_threshold
):
    """
    Compute cell-centered qs, theta, tau_b, u and upwind face qs.
    """
    qs_cell, theta, tau_b, u = engelund_hansen_transport(
        h=h,
        q=q,
        rho_w=rho_w,
        rho_s=rho_s,
        d50=d50,
        g=g,
        n_manning=n_manning,
        h_dry=h_dry,
        alpha_eh=alpha_eh,
        theta_cr=theta_cr,
        use_threshold=use_threshold
    )

    uf = 0.5 * (u[:-1] + u[1:])
    qs_face = np.zeros(len(h) - 1)

    pos = uf >= 0.0
    neg = ~pos

    qs_face[pos] = qs_cell[:-1][pos]
    qs_face[neg] = qs_cell[1:][neg]

    return qs_cell, qs_face, theta, tau_b, u
