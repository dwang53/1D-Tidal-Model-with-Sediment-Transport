import numpy as np

def tidal_stage(t, eta_mean, eta_amp, eta_period):
    return eta_mean + eta_amp * np.sin(2.0 * np.pi * t / eta_period)

def build_ghost_cells(h, q, zb, t, config):
    """
    One ghost cell on each side.
    Supported left boundary conditions:
      - tide
      - transmissive
      - wall
    Supported right boundary conditions:
      - transmissive
      - wall
    """
    n = len(h)

    hg = np.zeros(n + 2)
    qg = np.zeros(n + 2)
    zbg = np.zeros(n + 2)

    hg[1:-1] = h
    qg[1:-1] = q
    zbg[1:-1] = zb

    # Left boundary
    zbg[0] = zbg[1]

    if config.left_bc == "tide":
        eta_b = tidal_stage(t, config.eta_mean, config.eta_amp, config.eta_period)
        h_b = max(0.0, eta_b - zbg[0])

        # Characteristic-based velocity: use the outgoing Riemann invariant
        # (u - 2c) from the interior to determine the ghost-cell velocity.
        # R_minus = u_1 - 2*c_1  (outgoing, left-going characteristic)
        # At the boundary: u_b = R_minus + 2*c_b
        if hg[1] > config.h_dry and h_b > config.h_dry:
            u1 = qg[1] / hg[1]
            c1 = np.sqrt(config.g * hg[1])
            c_b = np.sqrt(config.g * h_b)
            u_b = u1 - 2.0 * c1 + 2.0 * c_b
        elif h_b > config.h_dry:
            u_b = 2.0 * np.sqrt(config.g * h_b)
        else:
            u_b = 0.0

        hg[0] = h_b
        qg[0] = h_b * u_b if h_b > config.h_dry else 0.0

    elif config.left_bc == "transmissive":
        hg[0] = hg[1]
        qg[0] = qg[1]

    elif config.left_bc == "wall":
        hg[0] = hg[1]
        qg[0] = -qg[1]

    else:
        raise ValueError(f"Unknown left_bc: {config.left_bc}")

    # Right boundary
    zbg[-1] = zbg[-2]

    if config.right_bc == "transmissive":
        hg[-1] = hg[-2]
        qg[-1] = qg[-2]

    elif config.right_bc == "wall":
        hg[-1] = hg[-2]
        qg[-1] = -qg[-2]

    else:
        raise ValueError(f"Unknown right_bc: {config.right_bc}")

    return hg, qg, zbg
