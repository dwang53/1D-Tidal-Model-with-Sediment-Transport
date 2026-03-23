import numpy as np

from src.boundary_conditions import build_ghost_cells
from src.hydrostatic import safe_velocity, hydrostatic_reconstruction
from src.reconstruction import reconstruct_cell_to_faces
from src.fluxes import hllc_flux

def hydro_rhs(h, q, zb, t, config, dx):
    """
    Semi-discrete operator for 1D SWE with:
      - ghost cells
      - MUSCL reconstruction of eta, u, zb
      - hydrostatic reconstruction
      - HLLC flux
      - Audusse-style flux correction for topography
    """
    hg, qg, zbg = build_ghost_cells(h, q, zb, t, config)

    etag = hg + zbg
    ug = safe_velocity(hg, qg, config.h_dry)

    etaL, etaR = reconstruct_cell_to_faces(etag)
    uL, uR = reconstruct_cell_to_faces(ug)
    zbL, zbR = reconstruct_cell_to_faces(zbg)

    # Original reconstructed depths before hydrostatic truncation
    hL0 = np.maximum(0.0, etaL - zbL)
    hR0 = np.maximum(0.0, etaR - zbR)

    hL, qL, hR, qR = hydrostatic_reconstruction(
        etaL, uL, zbL, etaR, uR, zbR, config.h_dry
    )

    f1, f2 = hllc_flux(hL, qL, hR, qR, config.g, config.h_dry)

    # Audusse-type correction on momentum flux
    corrL = 0.5 * config.g * (hL0**2 - hL**2)
    corrR = 0.5 * config.g * (hR0**2 - hR**2)

    f2L = f2 + corrL
    f2R = f2 + corrR

    n = len(h)
    rhs_h = np.zeros(n)
    rhs_q = np.zeros(n)

    rhs_h[:] = -(f1[1:n+1] - f1[0:n]) / dx
    rhs_q[:] = -(f2L[1:n+1] - f2R[0:n]) / dx

    return rhs_h, rhs_q
