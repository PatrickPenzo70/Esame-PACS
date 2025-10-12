import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os

# --- PARAMETRI BASE ---
vx_file = "vx_grid.csv"
vy_file = "vy_grid.csv"
particle_pattern = "particle_positions_*.csv"

# --- CARICA CAMPO DEL GAS ---
vx = np.loadtxt(vx_file, delimiter=",")
vy = np.loadtxt(vy_file, delimiter=",")
Ny, Nx = vx.shape

# dominio (coordinate regolari)
x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)
X, Y = np.meshgrid(x, y)

# --- LISTA DEI FILE DI PARTICELLE ---
particle_files = sorted(glob.glob(particle_pattern))
if not particle_files:
    raise FileNotFoundError("Nessun file di particelle trovato (particle_positions_*.csv)")

print(f"Trovati {len(particle_files)} file di particelle")

# --- CREA FIGURA ---
fig, ax = plt.subplots(figsize=(24, 12))
ax.set_title("Simulazione getto di gas con particelle", fontsize=14)
ax.set_xlabel("x")
ax.set_ylabel("y")

# campo colore per vx (velocità gas)
c = ax.pcolormesh(X, Y, vx, shading='auto', cmap='coolwarm', alpha=0.6)

# vettori di velocità gas (ridotti per chiarezza)
skip = 8
ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
          vx[::skip, ::skip], vy[::skip, ::skip],
          color='k', scale=30, alpha=0.4)

# inizializza particelle (scatter)
particles, = ax.plot([], [], 'ro', markersize=3, label='Particelle')

ax.legend(loc='upper right')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

# --- FUNZIONI ANIMAZIONE ---
def init():
    particles.set_data([], [])
    return (particles,)

def update(frame):
    df = pd.read_csv(particle_files[frame])
    x = df['x'].values
    y = df['y'].values
    particles.set_data(x, y)
    ax.set_title(f"Getto di gas - passo {frame} / {len(particle_files)}")
    return (particles,)

# --- CREA ANIMAZIONE ---
ani = animation.FuncAnimation(
    fig, update, frames=len(particle_files),
    init_func=init, blit=True, interval=50, repeat=False
)

plt.tight_layout()
plt.show()

