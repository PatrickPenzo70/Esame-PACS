#!/usr/bin/env python3
# ===============================================================
# Animazione 2D del processo di fusione (da output CSV del C++)
# ===============================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import pandas as pd
from scipy.interpolate import griddata
import os

# ---------------------------------------------------------------
# CONFIGURAZIONE
# ---------------------------------------------------------------
prefix = "heat_out_"        # prefisso file generati dal C++
video_filename = "melting_simulation.gif"  # usa .gif (compatibile senza ffmpeg)
fps = 10
dpi = 120
save_video = True           # salva anche video oltre a mostrare animazione

# ---------------------------------------------------------------
# LETTURA DEI FILE
# ---------------------------------------------------------------
files = sorted(glob.glob(prefix + "*.csv"), key=lambda s: int(s.split("_")[-1].split(".")[0]))
if not files:
    raise FileNotFoundError(f"Nessun file '{prefix}*.csv' trovato nella cartella: {os.getcwd()}")

print(f"Trovati {len(files)} file CSV")

frames = [pd.read_csv(f) for f in files]

# Costruisci una griglia regolare per l’interpolazione
x_unique = np.linspace(min(frames[0]["x"]), max(frames[0]["x"]), 100)
y_unique = np.linspace(min(frames[0]["y"]), max(frames[0]["y"]), 100)
X, Y = np.meshgrid(x_unique, y_unique)

# ---------------------------------------------------------------
# FUNZIONE PER INTERPOLARE SU GRIGLIA
# ---------------------------------------------------------------
def make_field(df, key):
    """Interpola il campo (T o phi) su una griglia regolare."""
    points = np.column_stack([df["x"], df["y"]])
    values = df[key]
    Z = griddata(points, values, (X, Y), method="linear", fill_value=np.nan)
    return Z

# Precalcola tutti i frame
T_fields = [make_field(f, "T") for f in frames]
phi_fields = [make_field(f, "phi") for f in frames]

# ---------------------------------------------------------------
# FIGURA E ANIMAZIONE
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axT, axPhi = axes

# Limiti di colore coerenti nel tempo
Tmin = np.nanmin([np.nanmin(T) for T in T_fields])
Tmax = np.nanmax([np.nanmax(T) for T in T_fields])

imT = axT.imshow(T_fields[0], extent=[X.min(), X.max(), Y.min(), Y.max()],
                 origin='lower', cmap='inferno', vmin=Tmin, vmax=Tmax)
imPhi = axPhi.imshow(phi_fields[0], extent=[X.min(), X.max(), Y.min(), Y.max()],
                     origin='lower', cmap='viridis', vmin=0, vmax=1)

axT.set_title("Temperatura [K]")
axPhi.set_title("Frazione fusa φ")

for ax in axes:
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")

cbT = fig.colorbar(imT, ax=axT, fraction=0.046, pad=0.04)
cbPhi = fig.colorbar(imPhi, ax=axPhi, fraction=0.046, pad=0.04)

# Testo del tempo
txt = fig.text(0.5, 0.92, "t = 0.0 s", ha="center", fontsize=12)

# ---------------------------------------------------------------
# FUNZIONE DI AGGIORNAMENTO FRAME
# ---------------------------------------------------------------
def update(frame_idx):
    imT.set_data(T_fields[frame_idx])
    imPhi.set_data(phi_fields[frame_idx])
    time_est = frame_idx * 0.002  # supponendo output ogni 0.002 s (da C++)
    txt.set_text(f"t ≈ {time_est:.4f} s")
    return imT, imPhi, txt

# ---------------------------------------------------------------
# CREA ANIMAZIONE
# ---------------------------------------------------------------
ani = animation.FuncAnimation(fig, update, frames=len(T_fields),
                              interval=200, blit=False, repeat=True)

plt.tight_layout()

# ---------------------------------------------------------------
# SALVATAGGIO VIDEO (GIF)
# ---------------------------------------------------------------
if save_video:
    print(f"Salvataggio animazione in {video_filename} ...")
    ani.save(video_filename, fps=fps, dpi=dpi)
    print(f"✅ Video salvato: {video_filename}")

plt.show()

