import numpy as np
import matplotlib.pyplot as plt
import sys
from commands import *
from trajectories import *
from roblib import *
from scipy.integrate import solve_ivp


# Calculate changes in state and velocity based on the current state, velocity, and command.
def model(n_state, v_state, command):

    dn = J(n_state) @ v_state
    dv = np.linalg.inv(M_model) @ (command - C(v_state) @ v_state - D(v_state) @ v_state)

    return [dn, dv]


# Draw the chosen trajectory.
def draw_traj(x_values, y_values):
    plt.plot(x_values, y_values, color="green", linestyle="--")


# Update the state and velocity over time based on the command.
def evolution(n, v, cmd, dt_ctr):
    steps = 100

    for _ in range(steps):
        dn, dv = model(n, v, cmd)
        v += dv * dt_ctr / steps
        n += dn * dt_ctr / steps

    n[2] = sawtooth(n[2])

    return n, v

# Show the current state, target, and trajectory on a plot.
def affichage(n, error_state, path_points):
    plt.clf()

    s1, y1, phi, s = error_state.flatten()
    # print(s)
    # xs, ys, phif, curv, dcurv_dt, g_c = path_description(s, ds, s_values, mu_values, phi_values)

    x_values, y_values = path_points[:2]
    xs, ys, _, phif = path_interrogation(s, path_points)[:4]

    x_limit = [-20, 10]
    y_limit = [-5, 35]

    draw_traj(x_values, y_values)

    plt.scatter(xs, ys, color="red", marker="x")

    plt.xlim(x_limit)
    plt.ylim(y_limit)

    R = np.array([[np.cos(phif), np.sin(phif), 0], 
                    [-np.sin(phif), np.cos(phif), 0], 
                    [0, 0, 1]],dtype=float,)

    tanx = xs + (np.linalg.inv(R) @ [5, 0, phi]).flatten()[0]
    tany = ys + (np.linalg.inv(R) @ [5, 0, phi]).flatten()[1]
    perpx = xs + (np.linalg.inv(R) @ [0, 5, phi]).flatten()[0]
    perpy = ys + (np.linalg.inv(R) @ [0, 5, phi]).flatten()[1]

    plt.annotate("",xy=(tanx, tany), xytext=(xs, ys), arrowprops=dict(arrowstyle="->", lw=1.5, color="blue"),)

    plt.annotate("",xy=(perpx, perpy), xytext=(xs, ys), arrowprops=dict(arrowstyle="->", lw=1.5, color="blue"),)

    draw_tank(n)

    plt.draw()
    plt.pause(0.001)


# Main function to run the simulation.
def main():

    path_points = get_path_points()
    x_values, y_values, s_values, mu_values, phi_values = draw_p()

    # Simulation settings.
    dt_ctr = 0.01
    total_time = 100.0
    steps = int(total_time / dt_ctr)

    # Initial state (position and orientation).
    n = np.array([[-8], [-2], [-np.pi]], dtype=float)  # x, y, phiU
    v = np.array([[2], [0], [0]], dtype=float)  # u, v, r

    s0 = 1
    u_target = 5
    error_state = initiate_error_state(n, v, s0, path_points )  # s1, y1, phi, s

    # Run the simulation for each time step.
    for step in range(steps):

        affichage(n, error_state, path_points)

        cmd, ds = command_auv(v, error_state, u_target, path_points, s_values, mu_values, phi_values)

        n, v = evolution(n, v, cmd, dt_ctr)

        error_state = update_error_state(error_state, v, n, ds, dt_ctr, path_points)

if __name__ == "__main__":
    main()
