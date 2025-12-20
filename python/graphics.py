import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import numpy as np
import glob
import os
from datetime import datetime, timedelta
from PIL import Image
import urllib.request

plt.style.use('dark_background')

# ============ CONFIGURATION ============
REFERENCE_DATE = datetime(2024, 1, 1, 0, 0, 0)
EARTH_RADIUS = 6371.0  # km
EARTH_ROTATION_PERIOD = 86400.0  # seconds
SPHERE_RESOLUTION = 50  # Quality of the sphere grid

# ============ DOWNLOAD EARTH TEXTURE ============
def download_earth_texture():
    texture_file = 'earth_texture.jpg'
    if not os.path.exists(texture_file):
        print("Downloading Earth texture map...")
        # NASA Blue Marble: Next Generation texture
        url = 'https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57752/land_shallow_topo_2048.jpg'
        try:
            urllib.request.urlretrieve(url, texture_file)
            print("Download complete!")
        except Exception as e:
            print(f"Could not download texture: {e}. Using solid color.")
            return None
    return texture_file

# ============ DATA LOADING ============
csv_files = glob.glob('output/*.csv')
if not csv_files:
    print("No CSV files found in 'output' directory! Please check your file path.")
    exit()

spacecraft_data = []
colors = plt.cm.tab10(np.linspace(0, 1, len(csv_files)))

for i, csv_file in enumerate(csv_files):
    df = pd.read_csv(csv_file)
    spacecraft_data.append({
        'name': os.path.basename(csv_file).replace('.csv', ''),
        'x': df['ECI_X'].values,
        'y': df['ECI_Y'].values,
        'z': df['ECI_Z'].values,
        'time': df['time'].values,
        'color': colors[i]
    })

# Interpolation Setup
all_times = np.concatenate([sc['time'] for sc in spacecraft_data])
global_time_min, global_time_max = all_times.min(), all_times.max()
num_frames = 100
global_times = np.linspace(global_time_min, global_time_max, num_frames)

for sc in spacecraft_data:
    sc['x_i'] = np.interp(global_times, sc['time'], sc['x'])
    sc['y_i'] = np.interp(global_times, sc['time'], sc['y'])
    sc['z_i'] = np.interp(global_times, sc['time'], sc['z'])

# ============ SPHERE & TEXTURE LOGIC ============
def get_earth_geometry(radius, res):
    """Create sphere and load texture aligned to grid."""
    # Create the meshgrid
    u = np.linspace(0, 2 * np.pi, res)
    v = np.linspace(0, np.pi, res)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

    texture = None
    t_file = download_earth_texture()
    if t_file:
        try:
            img = Image.open(t_file)
            # Resize image so pixels match the u,v grid dimensions exactly
            img = img.resize((res, res))
            texture = np.array(img) / 255.0
        except Exception as e:
            print(f"Texture processing error: {e}")
            
    return x, y, z, texture

sphere_x, sphere_y, sphere_z, earth_texture = get_earth_geometry(EARTH_RADIUS, SPHERE_RESOLUTION)

def rotate_z(x, y, z, angle):
    cos_a, sin_a = np.cos(angle), np.sin(angle)
    return x * cos_a - y * sin_a, x * sin_a + y * cos_a, z

# Initial Greenwich Angle calculation
j2000 = datetime(2000, 1, 1, 12, 0, 0)
days_since_j2000 = (REFERENCE_DATE - j2000).total_seconds() / 86400.0
gmst_at_ref = (280.46 + 360.98564736629 * days_since_j2000) % 360
reference_greenwich_angle = np.radians(gmst_at_ref)

# ============ PLOTTING ============
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')

# Initial Earth Plot
if earth_texture is not None:
    # shade=False prevents the 'dark sphere' effect by using the image's raw colors
    earth_surface = ax.plot_surface(sphere_x, sphere_y, sphere_z, 
                                    facecolors=earth_texture, 
                                    rstride=1, cstride=1, 
                                    antialiased=True, shade=False)
else:
    earth_surface = ax.plot_surface(sphere_x, sphere_y, sphere_z, color='dodgerblue', alpha=0.7)

# Spacecraft elements
for sc in spacecraft_data:
    sc['line'], = ax.plot([], [], [], '-', color=sc['color'], alpha=0.6, label=sc['name'])
    sc['point'], = ax.plot([], [], [], 'o', color=sc['color'], markersize=6)

# Formatting
ax.set_box_aspect([1, 1, 1])
limit = EARTH_RADIUS * 3  # Start with a reasonable zoom, adjustments happen below
all_coords = np.concatenate([sc['x'] for sc in spacecraft_data] + [sc['y'] for sc in spacecraft_data])
dynamic_limit = np.max(np.abs(all_coords)) * 1.1
ax.set_xlim(-dynamic_limit, dynamic_limit)
ax.set_ylim(-dynamic_limit, dynamic_limit)
ax.set_zlim(-dynamic_limit, dynamic_limit)

time_text = ax.text2D(0.05, 0.95, '', transform=ax.transAxes, weight='bold')
ax.legend(loc='upper right', fontsize=8)

# ============ ANIMATION ============
def update(frame):
    global earth_surface
    
    curr_t = global_times[frame]
    curr_date = REFERENCE_DATE + timedelta(seconds=float(curr_t))
    
    # Rotate Earth
    angle = reference_greenwich_angle + (2 * np.pi * curr_t / EARTH_ROTATION_PERIOD)
    rx, ry, rz = rotate_z(sphere_x, sphere_y, sphere_z, angle)
    
    earth_surface.remove()
    if earth_texture is not None:
        earth_surface = ax.plot_surface(rx, ry, rz, facecolors=earth_texture, 
                                        rstride=1, cstride=1, antialiased=True, shade=False)
    else:
        earth_surface = ax.plot_surface(rx, ry, rz, color='dodgerblue', alpha=0.7)
        
    # Update Spacecraft
    artists = [earth_surface]
    for sc in spacecraft_data:
        sc['line'].set_data(sc['x_i'][:frame+1], sc['y_i'][:frame+1])
        sc['line'].set_3d_properties(sc['z_i'][:frame+1])
        sc['point'].set_data([sc['x_i'][frame]], [sc['y_i'][frame]])
        sc['point'].set_3d_properties([sc['z_i'][frame]])
        artists.extend([sc['line'], sc['point']])
        
    time_text.set_text(f'UTC: {curr_date.strftime("%Y-%m-%d %H:%M:%S")}')
    artists.append(time_text)
    return artists

anim = FuncAnimation(fig, update, frames=num_frames, interval=30, blit=False)
plt.show()