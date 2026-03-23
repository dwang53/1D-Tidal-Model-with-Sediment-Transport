from dataclasses import dataclass

@dataclass
class ModelConfig:
    # domain
    x_min: float
    x_max: float
    nx: int

    # hydraulics
    g: float
    manning_n: float
    h_dry: float

    # morphology
    porosity: float
    rho_w: float
    rho_s: float
    d50: float
    theta_cr: float
    alpha_eh: float
    use_theta_threshold: bool
    enable_morphodynamics: bool

    # time
    t_end: float
    cfl: float
    dt_max: float
    output_interval: float

    # tide
    eta_mean: float
    eta_amp: float
    eta_period: float

    # boundary conditions
    left_bc: str     # "tide", "transmissive", "wall"
    right_bc: str    # "transmissive" or "wall"

    # initial condition
    eta0: float

    # optional benchmark metadata
    case_name: str = "unnamed_case"
