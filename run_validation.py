from validation.benchmarks import (
    lake_at_rest_case,
    dam_break_case,
    wetting_drying_beach_case,
    tidal_barrier_fixed_bed_case,
)
from validation.validation_runner import (
    run_lake_at_rest,
    run_dam_break,
    run_wetting_drying_beach,
    run_tidal_barrier_fixed_bed,
)

if __name__ == "__main__":
    print("Running validation suite...")

    print("\n[1/4] Lake at rest")
    run_lake_at_rest(lake_at_rest_case())

    print("\n[2/4] Dam break")
    run_dam_break(dam_break_case())

    print("\n[3/4] Wetting/drying beach")
    run_wetting_drying_beach(wetting_drying_beach_case())

    print("\n[4/4] Tidal barrier fixed bed")
    run_tidal_barrier_fixed_bed(tidal_barrier_fixed_bed_case())

    print("\nValidation suite complete.")
