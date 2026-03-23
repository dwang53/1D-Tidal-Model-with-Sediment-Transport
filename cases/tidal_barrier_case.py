from src.config import ModelConfig

def build_case():
    return ModelConfig(
        # domain
        x_min=0.0,
        x_max=1000.0,
        nx=600,

        # hydraulics
        g=9.81,
        manning_n=0.025,
        h_dry=1.0e-4,

        # morphology
        porosity=0.4,
        rho_w=1000.0,
        rho_s=2650.0,
        d50=2.0e-4,
        theta_cr=0.047,
        alpha_eh=0.05,
        use_theta_threshold=True,
        enable_morphodynamics=True,

        # time
        t_end=12.42 * 3600.0,
        cfl=0.30,
        dt_max=1.0,
        output_interval=120.0,

        # tide
        eta_mean=0.15,
        eta_amp=1.0,
        eta_period=12.42 * 3600.0,

        # boundary condition
        left_bc="tide",
        right_bc="transmissive",

        # initial condition
        eta0=0.15,

        case_name="tidal_barrier_morphodynamic"
    )from src.config import ModelConfig

def build_case():
    return ModelConfig(
        # domain
        x_min=0.0,
        x_max=1000.0,
        nx=600,

        # hydraulics
        g=9.81,
        manning_n=0.025,
        h_dry=1.0e-4,

        # morphology
        porosity=0.4,
        rho_w=1000.0,
        rho_s=2650.0,
        d50=2.0e-4,
        theta_cr=0.047,
        alpha_eh=0.05,
        use_theta_threshold=True,

        # time
        t_end=12.42 * 3600.0,
        cfl=0.30,
        dt_max=1.0,
        output_interval=120.0,

        # tide
        eta_mean=0.15,
        eta_amp=1.0,
        eta_period=12.42 * 3600.0,

        # boundary condition
        right_bc="transmissive",

        # initial condition
        eta0=0.15,
    )
