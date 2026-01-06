import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Folder containing your particle_position_*.csv files
DATA_FOLDER = '.'

# Load all CSV files in order
csv_files = sorted(glob.glob(os.path.join(DATA_FOLDER, 'particle_position_*.csv')))

# Map from label to list of positions
particle_tracks = {}

# Store time steps for synchronization
timesteps = []

for f in csv_files:
    df = pd.read_csv(f)
    time = df['Time'].iloc[0]
    timesteps.append(time)
    for _, row in df.iterrows():
        label = int(row['label'])
        pos = (row['x'], row['y'])
        if label not in particle_tracks:
            particle_tracks[label] = []
        particle_tracks[label].append(pos)

# Get particle count
num_particles = len(particle_tracks)

# Create the animation
fig, ax = plt.subplots(figsize=(8, 4))
sc = ax.scatter([], [], s=10)
trails = [ax.plot([], [], lw=0.5, alpha=0.5)[0] for _ in range(num_particles)]

def init():
    ax.set_xlim(-2, 2)
    ax.set_ylim(-1.5, 1.5)
    return sc, *trails

def update(frame):
    x = []
    y = []
    for i, label in enumerate(sorted(particle_tracks.keys())):
        track = particle_tracks[label]
        if frame < len(track):
            px, py = track[frame]
            x.append(px)
            y.append(py)
            trail_x = [p[0] for p in track[max(0, frame - 50):frame+1]]
            trail_y = [p[1] for p in track[max(0, frame - 50):frame+1]]
            trails[i].set_data(trail_x, trail_y)
    sc.set_offsets(list(zip(x, y)))
    ax.set_title(f"Time: {timesteps[frame]:.4f}s")
    return sc, *trails

ani = animation.FuncAnimation(
    fig, update, frames=len(timesteps),
    init_func=init, blit=True, interval=30
)

plt.tight_layout()
plt.show()

