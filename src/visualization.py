import os
import imageio.v2 as imageio
import matplotlib.pyplot as plt

def save_frame(step, t, x, zb, h, eta, q, theta):
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    axes[0].plot(x, zb, color="saddlebrown", linewidth=2, label="bed")
    axes[0].plot(x, eta, color="royalblue", linewidth=2, label="free surface")
    axes[0].fill_between(x, zb, eta, where=(h > 1e-8), color="lightskyblue", alpha=0.5)
    axes[0].set_ylabel("Elevation [m]")
    axes[0].legend(loc="best")
    axes[0].set_title(f"1D tidal barrier SWE–Exner model, t = {t:.1f} s")

    axes[1].plot(x, q, color="darkred", linewidth=1.5, label="discharge q")
    axes[1].axhline(0.0, color="k", linewidth=0.75)
    axes[1].set_ylabel("q [m$^2$/s]")
    axes[1].legend(loc="best")

    axes[2].plot(x, theta, color="darkgreen", linewidth=1.5, label="Shields parameter")
    axes[2].axhline(0.0, color="k", linewidth=0.75)
    axes[2].set_ylabel("theta [-]")
    axes[2].set_xlabel("x [m]")
    axes[2].legend(loc="best")

    fig.tight_layout()
    fname = f"output/frames/frame_{step:05d}.png"
    fig.savefig(fname, dpi=120)
    plt.close(fig)
    return fname

def build_gif():
    frames_dir = "output/frames"
    pngs = sorted(
        os.path.join(frames_dir, f)
        for f in os.listdir(frames_dir)
        if f.endswith(".png")
    )

    if not pngs:
        return

    images = [imageio.imread(png) for png in pngs]
    imageio.mimsave("output/gif/tidal_cycle.gif", images, duration=0.15)
