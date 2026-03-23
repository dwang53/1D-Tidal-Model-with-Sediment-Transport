import os
import pandas as pd

def ensure_output_dirs():
    os.makedirs("output/csv", exist_ok=True)
    os.makedirs("output/frames", exist_ok=True)
    os.makedirs("output/gif", exist_ok=True)

def write_profile_csv(step, t, x, zb, h, eta, u, q, qs, theta, tau_b):
    df = pd.DataFrame({
        "x": x,
        "zb": zb,
        "h": h,
        "eta": eta,
        "u": u,
        "q": q,
        "qs": qs,
        "theta": theta,
        "tau_b": tau_b
    })
    df.to_csv(f"output/csv/profile_{step:05d}.csv", index=False)

def write_timeseries_csv(records):
    df = pd.DataFrame(records)
    df.to_csv("output/csv/timeseries.csv", index=False)
