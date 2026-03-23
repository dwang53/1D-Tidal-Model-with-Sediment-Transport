import numpy as np

def apply_semiimplicit_manning(h, q, g, n_manning, h_dry, dt):
    """
    Semi-implicit Manning damping on momentum.
    """
    qn = q.copy()
    wet = h > h_dry

    if np.any(wet):
        h_eff = np.maximum(h[wet], h_dry)
        alpha = dt * g * n_manning**2 * np.abs(qn[wet]) / (h_eff ** (10.0 / 3.0))
        qn[wet] = qn[wet] / (1.0 + alpha)

    qn[~wet] = 0.0
    return qn

def bed_shear_stress_manning(h, u, rho_w, g, n_manning, h_dry):
    """
    Bed shear stress magnitude signed with flow direction:
      tau_b = rho g n^2 u |u| / h^(1/3)
    """
    tau = np.zeros_like(h)
    wet = h > h_dry
    tau[wet] = rho_w * g * n_manning**2 * u[wet] * np.abs(u[wet]) / (h[wet] ** (1.0 / 3.0))
    return tau
