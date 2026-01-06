import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ======================= PARAMETRI ========================
NX, NY = 120, 120      # celle della griglia
Lx, Ly = 1e-3, 1e-3    # dimensioni [m]
dx = Lx / NX
dt = 2e-7              # [s]
steps = 10000           # passi totali
frameskip = 10         # ogni quanti step aggiornare il frame

# Proprietà del materiale (Cr2O3)
rho = 5200.0           # [kg/m³]
cp_s, cp_l = 800.0, 1000.0
k_s, k_l = 5.0, 3.0
Tm = 2435.0            # [K]
Lm = 4e5               # [J/kg]
Kphi = 1e-3            # [1/s/K]

# Condizioni iniziali
T_inf = 15000.0        # plasma [K]
T0 = 300.0             # particella fredda [K]
R = 100e-6             # raggio particella [m]
cx, cy = Lx/2, Ly/2    # centro particella

# ======================= GRIGLIA ==========================
x = np.linspace(0, Lx, NX)
y = np.linspace(0, Ly, NY)
X, Y = np.meshgrid(x, y)
r = np.sqrt((X - cx)**2 + (Y - cy)**2)
mask = (r < R).astype(float)  # 1 = particella, 0 = gas

T = T_inf * np.ones((NY, NX))
T[mask == 1] = T0
phi = np.zeros_like(T)

# ======================= FUNZIONI =========================
def laplacian(T, dx):
    """Laplaciano 2D con differenze finite"""
    return (
        -4*T
        + np.roll(T, 1, axis=0)
        + np.roll(T, -1, axis=0)
        + np.roll(T, 1, axis=1)
        + np.roll(T, -1, axis=1)
    ) / dx**2

# ======================= ANIMAZIONE ========================
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
im_T = axs[0].imshow(T, origin='lower', extent=[0,Lx*1e3,0,Ly*1e3], cmap='inferno', vmin=300, vmax=T_inf)
axs[0].set_title("Temperatura [K]")
im_phi = axs[1].imshow(phi, origin='lower', extent=[0,Lx*1e3,0,Ly*1e3], cmap='viridis', vmin=0, vmax=1)
axs[1].set_title("Frazione fusa φ")

fig.suptitle("Fusione di una particella di Cr₂O₃ nel plasma", fontsize=14)
fig.tight_layout()

def update(frame):
    global T, phi
    # Proprietà locali
    cp = cp_s*(1 - phi) + cp_l*phi
    k = k_s*(1 - phi) + k_l*phi
    alpha = k / (rho * cp)

    # Aggiornamento fusione
    dphi = np.zeros_like(phi)
    idx = (mask > 0) & (T > Tm) & (phi < 1.0)
    dphi[idx] = Kphi * (T[idx] - Tm) * (1 - phi[idx])
    phi[idx] = np.clip(phi[idx] + dphi[idx] * dt, 0, 1)

    # Sorgente latente
    S = -rho * Lm * dphi * mask

    # Diffusione termica esplicita
    lapT = laplacian(T, dx)
    T += dt * (alpha * lapT + S / (rho * cp))

    # Il plasma resta caldo
    T[mask == 0] = T_inf

    # Aggiorna grafica
    im_T.set_data(T)
    im_phi.set_data(phi)
    axs[0].set_title(f"Temperatura [K]  step={frame*frameskip}")
    return [im_T, im_phi]

anim = FuncAnimation(fig, update, frames=steps//frameskip, interval=50, blit=False)
plt.show()

