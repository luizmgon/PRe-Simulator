import numpy as np
import matplotlib.pyplot as plt
import sys
from commands import *
from trajectories import *
from roblib import *

idx = 0

# Calculate changes in state and velocity based on the current state, velocity, and command.
def model(n_state, v_state, command):
    wind = np.array([[-200],[200],[0]])
    # wind = 0
    dn = J(n_state) @ v_state
    dv = np.linalg.inv(M) @ (command - C(v_state) @ v_state - D(v_state) @ v_state + wind)
    return [dn, dv]

# Find the closest point on the trajectory to the current state.
def get_closest(n_state, traj, path_points):
    x, y, phi = n_state.flatten()

    if traj == "l":
        closest = np.array([[5], [y]])
        return closest
    if traj == "c":
        dir = np.array([[x - 5], [y - 5]]) / np.sqrt((x - 5) ** 2 + (y - 5) ** 2)
        closest = np.array([[5], [5]]) + 2 * dir
        return closest
    else:
        return s_closest(n_state[0,0], n_state[1,0], path_points)

# Get the direction of the trajectory at the closest point.
def get_tangent(closest, traj, path_points):
    if traj == "l":
        tangent = np.array([[0], [1]])
    if traj == "c":
        dir = (closest - np.array([[5], [5]])) / 2
        dir = dir.flatten()
        tangent = np.array([[-dir[1]], [dir[0]]])
    else:
        s = closest
        xs, ys, _, phif, curv, g_c, dx, dy, ddx, ddy = path_interrogation(s, path_points)

        R = np.array([[np.cos(phif), np.sin(phif)], 
                    [-np.sin(phif), np.cos(phif)]],dtype=float,)

        tanx = (np.linalg.inv(R) @ [1, 0]).flatten()[0]
        tany = (np.linalg.inv(R) @ [1, 0]).flatten()[1]
        tangent = np.array([[tanx], [tany]])

    return tangent

# Get the target point to follow on the trajectory.
def get_target(n_state, traj, path_points):
    closest = get_closest(n_state, traj, path_points)
    tangent = get_tangent(closest, traj, path_points)

    if(traj != 'l' and traj != 'c'):
        xs, ys, _, phif, curv, g_c, dx, dy, ddx, ddy = path_interrogation(closest, path_points)
        closest = np.array([[xs], [ys]])


    target = closest + 2 * tangent
    return closest, target

# Draw the chosen trajectory (line or circle).
def draw_traj(traj, path_points):
    if traj == "l":
        return draw_traj_line()
    elif traj == 'c':
        return draw_traj_circle()
    else: draw_path(path_points)

def draw_path(path_points):
    x_values, y_values = path_points[:2]
    plt.plot(x_values, y_values, color="green", linestyle="--")


# Update the state and velocity over time based on the command.
def evolution(n, v, cmd, dt_ctr, positions_x, positions_y, times):
    steps = 100

    for _ in range(steps):
        dn, dv = model(n, v, cmd)
        n += dn * dt_ctr / steps
        v += dv * dt_ctr / steps

        time = times[-1] + dt_ctr / steps
        positions_x.append(n[0][0])
        positions_y.append(n[1][0])
        times.append(time)

    while n[2] > 2 * np.pi:
        n[2] -= 2 * np.pi
    while n[2] < -2 * np.pi:
        n[2] += 2 * np.pi

    return n, v

# Show the current state, target, and trajectory on a plot.
def affichage(target, closest, traj, positions_x, positions_y, n, path_points):


    plt.clf()

    # x_limit = [-1, 10]
    # y_limit = [-1, 10]
    x_limit = [-5, 30]
    y_limit = [-5, 35]
    x_limit = [min(path_points[0]) - 5, max(path_points[0]) + 5]
    y_limit = [min(path_points[1]) - 5, max(path_points[1]) + 5]

    draw_traj(traj, path_points)

    plt.scatter(target[0], target[1], color="red", marker="x")
    plt.scatter(closest[0], closest[1], color="blue", marker="o")
    plt.plot(positions_x, positions_y, color="green")
    plt.xlim(x_limit)
    plt.ylim(y_limit)
    draw_tank(n)

    plt.draw()
    plt.pause(0.001)

# Main function to run the simulation.
def main():

    # traj = sys.argv[1]
    traj = "3"

    # Simulation settings.
    dt_ctr = 0.02
    total_time = 100.0
    steps = int(total_time / dt_ctr)

    # Initial state (position and orientation).
    n = np.array([[0], [20], [0]], dtype=float)  # x, y, Φ
    v = np.array([[0], [0], [0]], dtype=float)  # u, v, r

    # Lists to store the x and y positions and times.
    positions_x = [n[0][0]]
    positions_y = [n[1][0]]
    times = [0]

    u_target = 5

    previous_theta_d = None

    path_points = get_path_points()

    error_integral = np.zeros((3,1))

    # Run the simulation for each time step.
    for step in range(steps):

        closest, target = get_target(n, traj, path_points)
        previous_theta_d, cmd, error_integral = command_los(n, v, target, previous_theta_d, u_target, dt_ctr, error_integral)

        n, v = evolution(n, v, cmd, dt_ctr, positions_x, positions_y, times)

        affichage(target, closest, traj, positions_x, positions_y, n, path_points)

    plt.show()


if __name__ == "__main__":
    main()
