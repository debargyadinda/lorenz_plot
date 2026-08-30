#!/usr/bin/env python3
"""
plot_lorenz.py — visualize the Lorenz attractor, its Lyapunov exponent
convergence, and the butterfly-effect divergence of two nearby trajectories.

Run after:
    ./lorenz_chaos          # writes lorenz_trajectory.csv, lyapunov_convergence.csv
    ./lorenz_divergence     # writes lorenz_divergence.csv
"""
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3D projection)

# ---------- 1. The attractor itself (3D butterfly shape) ----------
traj = pd.read_csv("lorenz_trajectory.csv")

fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection="3d")
ax.plot(traj["x"], traj["y"], traj["z"], color="#2c5f8a", linewidth=0.6)
ax.set_title("Lorenz Attractor — Trajectory in Phase Space")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
fig.tight_layout()
fig.savefig("lorenz_attractor.png", dpi=160)
print("Wrote lorenz_attractor.png")

# ---------- 2. Lyapunov exponent convergence ----------
lyap = pd.read_csv("lyapunov_convergence.csv")

fig2, ax2 = plt.subplots(figsize=(9, 5))
ax2.plot(lyap["renorm_step"], lyap["running_lambda"], color="#1f77b4", linewidth=1.3)
ax2.axhline(0.9056, color="#d62728", linestyle="--", linewidth=1.2,
            label="Known reference value ≈ 0.9056")
ax2.set_xlabel("Renormalization step")
ax2.set_ylabel(r"Running estimate of $\lambda$")
ax2.set_title("Convergence of the Largest Lyapunov Exponent (Benettin's Method)")
ax2.legend()
ax2.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig("lyapunov_convergence.png", dpi=160)
print("Wrote lyapunov_convergence.png")

# ---------- 3. The butterfly effect: two trajectories diverging ----------
div = pd.read_csv("lorenz_divergence.csv")

fig3 = plt.figure(figsize=(9, 7))
ax3 = fig3.add_subplot(111, projection="3d")
ax3.plot(div["x_a"], div["y_a"], div["z_a"], color="#2c5f8a", linewidth=0.8,
         label="Trajectory A")
ax3.plot(div["x_b"], div["y_b"], div["z_b"], color="#d62728", linewidth=0.8,
         label="Trajectory B (started 1e-8 away)")
ax3.set_title("The Butterfly Effect: Two Trajectories, Identical Start (Δ=1e-8)")
ax3.set_xlabel("x")
ax3.set_ylabel("y")
ax3.set_zlabel("z")
ax3.legend()
fig3.tight_layout()
fig3.savefig("lorenz_butterfly.png", dpi=160)
print("Wrote lorenz_butterfly.png")

# ---------- 4. Separation growing exponentially (log scale) ----------
fig4, ax4 = plt.subplots(figsize=(9, 5))
ax4.semilogy(div["t"], div["separation"], color="#2c5f8a", linewidth=1.3)
ax4.set_xlabel("Time")
ax4.set_ylabel("Separation between trajectories (log scale)")
ax4.set_title("Exponential Divergence of Two Initially Identical Trajectories")
ax4.grid(True, which="both", alpha=0.3)
fig4.tight_layout()
fig4.savefig("lorenz_separation.png", dpi=160)
print("Wrote lorenz_separation.png")
