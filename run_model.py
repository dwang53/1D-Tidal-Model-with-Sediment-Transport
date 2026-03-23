from cases.tidal_barrier_case import build_case
from src.swe_solver import run_simulation

if __name__ == "__main__":
    config = build_case()
    run_simulation(config)
