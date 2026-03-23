import numpy as np
from src.swe_solver import run_simulation
from src.grid import build_grid
from src.state import initial_state

def run_case(config, bathymetry_builder=None, initial_condition_builder=None):
    """
    Generic validation helper. Monkey-patches the imported bathymetry/state
    behavior by temporarily wrapping the solver setup pattern externally.

    For simplicity in this repository structure, benchmark-specific setup is
    done by swapping the src.bathymetry.initial_bathymetry function if needed.
    """
    if bathymetry_builder is None and initial_condition_builder is None:
        run_simulation(config)
        return

    from src import bathymetry as bathy_mod
    from src import swe_solver as solver_mod

    old_initial_bathymetry = bathy_mod.initial_bathymetry
    old_solver_initial_bathymetry = solver_mod.initial_bathymetry

    def patched_initial_bathymetry(x):
        if bathymetry_builder is None:
            return old_initial_bathymetry(x)
        return bathymetry_builder(x)

    bathy_mod.initial_bathymetry = patched_initial_bathymetry
    solver_mod.initial_bathymetry = patched_initial_bathymetry

    if initial_condition_builder is None:
        run_simulation(config)
    else:
        old_initial_state = solver_mod.initial_state

        def patched_initial_state(x, zb, eta0, h_dry):
            return initial_condition_builder(x, zb, eta0, h_dry)

        solver_mod.initial_state = patched_initial_state
        try:
            run_simulation(config)
        finally:
            solver_mod.initial_state = old_initial_state

    bathy_mod.initial_bathymetry = old_initial_bathymetry
    solver_mod.initial_bathymetry = old_solver_initial_bathymetry

def lake_at_rest_bathymetry(x):
    return 0.4 * np.exp(-((x - 500.0) / 120.0) ** 2)

def lake_at_rest_initial_condition(x, zb, eta0, h_dry):
    eta = 1.0 * np.ones_like(x)
    h = np.maximum(eta - zb, 0.0)
    q = np.zeros_like(x)
    return h, q

def dam_break_bathymetry(x):
    return np.zeros_like(x)

def dam_break_initial_condition(x, zb, eta0, h_dry):
    h = np.where(x < 500.0, 2.0, 1.0)
    q = np.zeros_like(x)
    return h, q

def wetting_drying_beach_bathymetry(x):
    # linearly rising beach toward the right
    return -0.50 + 0.002 * x

def wetting_drying_beach_initial_condition(x, zb, eta0, h_dry):
    eta = 0.20 * np.ones_like(x)
    h = np.maximum(eta - zb, 0.0)
    q = np.zeros_like(x)
    return h, q

def run_lake_at_rest(config):
    run_case(
        config,
        bathymetry_builder=lake_at_rest_bathymetry,
        initial_condition_builder=lake_at_rest_initial_condition
    )

def run_dam_break(config):
    run_case(
        config,
        bathymetry_builder=dam_break_bathymetry,
        initial_condition_builder=dam_break_initial_condition
    )

def run_wetting_drying_beach(config):
    run_case(
        config,
        bathymetry_builder=wetting_drying_beach_bathymetry,
        initial_condition_builder=wetting_drying_beach_initial_condition
    )

def run_tidal_barrier_fixed_bed(config):
    run_case(config)
