from src.limiters import limited_slope

def reconstruct_cell_to_faces(w, theta=1.5):
    """
    MUSCL reconstruction from cell averages to face states.

    Returns arrays of length N-1:
      wL[i] = state at interface i+1/2 reconstructed from cell i
      wR[i] = state at interface i+1/2 reconstructed from cell i+1
    """
    slope = limited_slope(w, theta=theta)
    wL = w[:-1] + 0.5 * slope[:-1]
    wR = w[1:]  - 0.5 * slope[1:]
    return wL, wR
