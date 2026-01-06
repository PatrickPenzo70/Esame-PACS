import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
from matplotlib.animation import FuncAnimation, PillowWriter

# --- Load all CSVs ---
files = sorted(glob.glob("static_out_*_grid.csv"))
if not files:
    raise FileNotFoundError("No CSV files found! Run your C++ simulation first.")

# --- Read first file to get grid dimensions ---
data0 = pd.read_csv(files[0])
NX = len(data0['i'].unique())
NY = len(data0['j'].unique())
x = data0['x'].unique()
y = data0['y'].unique()
extent = [x.min()*1e6, x.max()*1e6, y.min()*1e6, y.max()*1e6]  # microns

# --- Preload temperature and phase data ---
T_fields, phi_fields, mask_fields = [], [], []
for f in files:
    data = pd.read_csv(f)
    T_grid = data.pivot(index='j', columns='i', values='T').to_numpy()
    phi_grid = data.pivot(index='j', columns='i', values='phi').to_numpy()
    mask_grid = data.pivot(index='j', columns='i', values='mask').to_numpy()
    T_fields.append(T_grid)
    phi_fields.append(phi_grid)
    mask_fields.append(mask_grid)

# --- Figure setup ---
fig, ax = plt.subplots(figsize=(6,5))
vmin, vmax = np.min(T_fields[0]), np.max(T_fields[-1])
imT = ax.imshow(T_fields[0], origin='lower', extent=extent, cmap='inferno', vmin=vmin, vmax=vmax)
imPhi = ax.imshow(phi_fields[0], origin='lower', extent=extent, cmap='cool', alpha=0.4, vmin=0, vmax=1)
contour = ax.contour(mask_fields[0], levels=[0.5], colors='white',
                     extent=extent, origin='lower', linewidths=1)
ax.set_xlabel('x [μm]')
ax.set_ylabel('y [μm]')
ax.set_title('Temperature + Melting Fraction Evolution')
cbar = fig.colorbar(imT, ax=ax, label='Temperature [K]')
cbar_phi = fig.colorbar(imPhi, ax=ax, label='Melting fraction φ', fraction=0.046, pad=0.04)

# --- Update function ---
def update(frame):
    ax.clear()
    T = T_fields[frame]
    phi = phi_fields[frame]
    mask = mask_fields[frame]

    # Draw temperature
    imT = ax.imshow(T, origin='lower', extent=extent, cmap='inferno', vmin=vmin, vmax=vmax)
    # Overlay φ as semi-transparent
    imPhi = ax.imshow(phi, origin='lower', extent=extent, cmap='cool', alpha=0.45, vmin=0, vmax=1)
    # Particle boundary
    ax.contour(mask, levels=[0.5], colors='white', extent=extent, origin='lower', linewidths=1)
    ax.set_xlabel('x [μm]')
    ax.set_ylabel('y [μm]')
    ax.set_title(f"Step {frame} | File: {files[frame]}")
    return [imT, imPhi]

# --- Create animation ---
anim = FuncAnimation(fig, update, frames=len(T_fields), interval=200, blit=False)

# --- Save as GIF ---
anim.save("particle_heating_melting.gif", writer=PillowWriter(fps=5))
print("✅ Saved animation: particle_heating_melting.gif")

plt.show()

