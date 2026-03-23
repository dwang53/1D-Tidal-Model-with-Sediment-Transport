import numpy as np
from src.config import ModelConfig

def lake_at_rest_case():
    cfg = ModelConfig(
        x_min=0.0,
        x_max=1000.0,
        nx=400,

        g=9.81,
        manning_n=0.0,
        h_dry=1.0e-6,

        porosity=0.4,
        rho_w=1000.0,
        rho_s=2650.0,
        d50=2.0e-4,
        theta_cr=0.047,
        alpha_eh=0.05,
        use_theta_threshold=True,
        enable_morphodynamics=False,

        t_end=200.0,
        cfl=0.30,
        dt_max=1.0,
        output_interval=20.0,

        eta_mean=1.0,
        eta_amp=0.0,
        eta_period=1.0,

        left_bc="wall",
        right_bc="wall",

        eta0=1.0,
        case_name="validation_lake_at_rest"
    )
    return cfg

def dam_break_case():
    cfg = ModelConfig(
        x_min=0.0,
        x_max=1000.0,
        nx=800,

        g=9.81,
        manning_n=0.0,
        h_dry=1.0e-6,

        porosity=0.4,
        rho_w=1000.0,
        rho_s=2650.0,
        d50=2.0e-4,
        theta_cr=0.047,
        alpha_eh=0.05,
        use_theta_threshold=True,
        enable_morphodynamics=False,

        t_end=60.0,
        cfl=0.25,
        dt_max=0.2,
        output_interval=2.0,

        eta_mean=0.0,
        eta_amp=0.0,
        eta_period=1.0,

        left_bc="wall",
        right_bc="wall",

        eta0=0.0,
        case_name="validation_dam_break"
    )
    return cfg

def wetting_drying_beach_case():
    cfg = ModelConfig(
        x_min=0.0,
        x_max=600.0,
        nx=600,

        g=9.81,
        manning_n=0.01,
        h_dry=1.0e-5,

        porosity=0.4,
        rho_w=1000.0,
        rho_s=2650.0,
        d50=2.0e-4,
        theta_cr=0.047,
        alpha_eh=0.05,
        use_theta_threshold=True,
        enable_morphodynamics=False,

        t_end=600.0,
        cfl=0.25,
        dt_max=0.25,
        output_interval=10.0,

        eta_mean=0.20,
        eta_amp=0.35,
        eta_period=120.0,

        left_bc="tide",
        right_bc="wall",

        eta0=0.20,
        case_name="validation_wetting_drying_beach"
    )
    return cfg

def tidal_barrier_fixed_bed_case():
    cfg = ModelConfig(
        x_min=0.0,
        x_max=1000.0,
        nx=600,

        g=9.81,
        manning_n=0.025,
        h_dry=1.0e-4,

        porosity=0.4,
        rho_w=1000.0,
        rho_s=2650.0,
        d50=2.0e-4,
        theta_cr=0.047,
        alpha_eh=0.05,
        use_theta_threshold=True,
        enable_morphodynamics=False,

        t_end=12.42 * 3600.0,
        cfl=0.30,
        dt_max=1.0,
        output_interval=120.0,

        eta_mean=0.15,
        eta_amp=1.0,
        eta_period=12.42 * 3600.0,

        left_bc="tide",
        right_bc="transmissive",

        eta0=0.15,
        case_name="validation_tidal_barrier_fixed_bed"
    )
    return cfg
