import numpy as np
import matplotlib.pyplot as plt
import sys
from commands import *
from trajectories import *
from roblib import *

idx = 0

# Calculate changes in state and velocity based on the current state, velocity, and command.
def model(n_state, v_state, command):
    dn = J(n_state) @ v_state
    dv = np.linalg.inv(M) @ (command - C(v_state) @ v_state - D(v_state) @ v_state)
    return [dn, dv]

# Find the closest point on the trajectory to the current state.
def get_closest(n_state, traj):
    x, y, phi = n_state.flatten()

    if traj == "l":
        closest = np.array([[5], [y]])
        return closest
    if traj == "c":
        dir = np.array([[x - 5], [y - 5]]) / np.sqrt((x - 5) ** 2 + (y - 5) ** 2)
        closest = np.array([[5], [5]]) + 2 * dir
        return closest

# Get the direction of the trajectory at the closest point.
def get_tangent(closest, traj):
    if traj == "l":
        tangent = np.array([[0], [1]])
    if traj == "c":
        dir = (closest - np.array([[5], [5]])) / 2
        dir = dir.flatten()
        tangent = np.array([[-dir[1]], [dir[0]]])
    return tangent

# Get the target point to follow on the trajectory.
def get_target(n_state, traj):
    closest = get_closest(n_state, traj)
    tangent = get_tangent(closest, traj)
    target = closest + 2 * tangent
    return closest, target

# Draw the chosen trajectory (line or circle).
def draw_traj(traj):
    if traj == "l":
        return draw_traj_line()
    else:
        return draw_traj_circle()

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
def affichage(target, closest, traj, positions_x, positions_y, n):
    plt.clf()

    x_limit = [-1, 10]
    y_limit = [-1, 10]

    draw_traj(traj)

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

    traj = sys.argv[1]

    # Simulation settings.
    dt_ctr = 0.01
    total_time = 100.0
    steps = int(total_time / dt_ctr)

    # Initial state (position and orientation).
    n = np.array([[0], [0], [0]], dtype=float)  # x, y, Φ
    v = np.array([[1], [0], [0]], dtype=float)  # u, v, r

    # Lists to store the x and y positions and times.
    positions_x = [n[0][0]]
    positions_y = [n[1][0]]
    times = [0]

    previous_theta_d = None

    # Run the simulation for each time step.
    for step in range(steps):

        closest, target = get_target(n, traj)
        previous_theta_d, cmd = command_los(n, v, target, previous_theta_d)

        n, v = evolution(n, v, cmd, dt_ctr, positions_x, positions_y, times)

        affichage(target, closest, traj, positions_x, positions_y, n)

    plt.show()


if __name__ == "__main__":
    main()
