import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path


output_directory = Path("output")
names = []

for i, file in enumerate(output_directory.iterdir()):

    print(type(file))
    
    names.append(str(file).split("\\")[1].split("_output.csv")[0])

    sim_df = pd.read_csv(file)

    print(sim_df)

    # Earth sphere
    r_earth = 6378 
    points = np.linspace(0, 2*np.pi, 1000)
    x_earth = r_earth*np.cos(points)
    y_earth = r_earth*np.sin(points)

    plt.style.use('dark_background')

    if i == 0:
        names.insert(0, "Earth")

    # ECI X-Y
    plt.figure(1)
    plt.title("ECI X vs. Y plane")

    if i == 0:
        plt.plot(x_earth, y_earth, color='blue')
    
    plt.plot(sim_df.ECI_X, sim_df.ECI_Y)
    plt.xlabel("X [km]")
    plt.ylabel("Y [km]")
    plt.axis('square')
    plt.legend("Earth")
    plt.legend(names)

    # ECI X-Z
    plt.figure(2)
    plt.title("ECI X vs. Z plane")

    if i == 0:
        plt.plot(x_earth, y_earth, color='blue')
    
    plt.plot(sim_df.ECI_X, sim_df.ECI_Z)
    plt.xlabel("X [km]")
    plt.ylabel("Z [km]")
    plt.axis('square')
    plt.legend(names)

    # ECI y-Z
    plt.figure(3)
    plt.title("ECI Y vs. Z plane")

    if i == 0:
        plt.plot(x_earth, y_earth, color='blue')
    
    plt.plot(sim_df.ECI_Y, sim_df.ECI_Z)
    plt.xlabel("y [km]")
    plt.ylabel("Z [km]")
    plt.axis('square')
    plt.legend(names)

    print(names)
    
    if i == 0:
        names.pop(0)
    

    # Classical Orbital Elements vs. time 
    plt.figure(4)
    plt.title("Semi-major axis vs. time")
    plt.plot(sim_df.time, sim_df.a)
    plt.xlabel("time [s]")
    plt.ylabel("a [km]")
    plt.legend(names)

    plt.figure(5)
    plt.title("Eccentricity vs. time")
    plt.plot(sim_df.time, sim_df.e)
    plt.xlabel("time [s]")
    plt.ylabel("e [-]")
    plt.legend(names)

    plt.figure(6)
    plt.title("Inclination vs. time")
    plt.plot(sim_df.time, sim_df.i)
    plt.xlabel("time [s]")
    plt.ylabel("i [rad]")
    plt.legend(names)

    plt.figure(7)
    plt.title("Longitude of Ascending Node vs. time")
    plt.plot(sim_df.time, sim_df.laan)
    plt.xlabel("time [s]")
    plt.ylabel("Ω [rad]")
    plt.legend(names)

    plt.figure(8)
    plt.title("Argument of Periapsis vs. time")
    plt.plot(sim_df.time, sim_df.gamma)
    plt.xlabel("time [s]")
    plt.ylabel("ω [rad]")
    plt.legend(names)

    plt.figure(9)
    plt.title("True anomaly vs. time")
    plt.plot(sim_df.time, sim_df.f)
    plt.xlabel("time [s]")
    plt.ylabel("f [rad]")
    plt.legend(names)


plt.show()


