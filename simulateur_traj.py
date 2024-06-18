import numpy as np
import matplotlib.pyplot as plt
import sys
from commands import *
from trajectories import *
from roblib import *
from scipy.integrate import solve_ivp


def model(n_state, v_state, command):
    dn = J(n_state) @ v_state
    dv = np.linalg.inv(M_model) @ (command - C(v_state) @ v_state - D(v_state) @ v_state)
    return [dn, dv]


def model2(t, y, cmd):
    n_state = y[:3].reshape((3, 1))
    v_state = y[3:].reshape((3, 1))

    dn = J(n_state) @ v_state
    dv = np.linalg.inv(M_model) @ (cmd - C(v_state) @ v_state - D(v_state) @ v_state)

    dydt = np.vstack((dn, dv)).flatten()
    return dydt


def command(com, n_state, v_state, target, d_target, dd_target):
    
    target = target[:2]
    d_target = d_target[:2]
    dd_target = dd_target[:2]  
    dn_state = (J(n_state) @ v_state)[:2]

    k = None

    if com == "fblindiag":
        k = dd_target + 2 * (d_target - dn_state) + (target - n_state[:2])
    
    elif com == "slid1diag":
        k = (dd_target + 1*(d_target - dn_state) + 0.5*np.sign(1*(d_target - dn_state) +1*(target - n_state[:2])))

    elif com == "slid2diag":
        k = 1*np.sign(2*(d_target - dn_state) + 2*(target - n_state[:2]))

    if(k is not None): return command_sous_diag(n_state, v_state, k)

    #Fonctionnent pas sans seuil
    if com == "fblin":
        k = dd_target + 2 * (d_target - dn_state) + 1*(target - n_state[:2])
        # k = k/10
    elif com == "slid1":
        k = 1*dd_target +1*(d_target - dn_state) + 10*np.sign((d_target - dn_state) +1*(target - n_state[:2]))
        k = k/10
    elif com == "slid2":
        k = 4*np.sign((d_target - dn_state) + 1*(target - n_state[:2]))

    if(k is not None): return command_sous(n_state, v_state, k)

    else:
        print("Invalid input")
        exit()


def get_new_target(time, traj, path_points):

    if traj == "l":
        return get_target_line(time)
    elif traj == 'c':
        return get_target_circle(time)
    else: return get_target_path(time, path_points)
    


def draw_traj(traj, path_points):
    if traj == "l":
        return draw_traj_line()
    elif traj == 'c':
        return draw_traj_circle()
    else: return draw_traj_path(path_points)


# Show the current state, target, and trajectory on a plot.
def affichage(target, traj, positions_x, positions_y, n, path_points):
    plt.clf()

    x_limit = [-10, 5]
    y_limit = [-1, 32]

    draw_traj(traj, path_points)

    plt.scatter(target[0], target[1], color="red", marker="o")
    plt.plot(positions_x, positions_y, color="blue")
    plt.xlim(x_limit)
    plt.ylim(y_limit)
    draw_tank(n)

    plt.draw()
    plt.pause(0.001)


def evolution(n, v, cmd, dt_ctr, positions_x, positions_y, times):
    steps = 1

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


def evolution2(n0, v0, cmd, dt_ctr, positions_x, positions_y, times):

    t_span = (0, dt_ctr)

    # Define o vetor de estado inicial concatenando n0 e v0
    y0 = np.hstack((n0.flatten(), v0.flatten()))

    # Chama solve_ivp com os parâmetros necessários
    sol = solve_ivp(model2, t_span, y0, args=(cmd,), max_step=dt_ctr / 3)

    # Ajusta o ângulo para estar no intervalo [-2π, 2π]
    while sol.y[2][-1] > 2 * np.pi:
        sol.y[2][-1] -= 2 * np.pi
    while sol.y[2][-1] < -2 * np.pi:
        sol.y[2][-1] += 2 * np.pi

    n = sol.y[:3, -1].reshape((3, 1))
    v = sol.y[3:, -1].reshape((3, 1))

    positions_x += list(sol.y[0][1:])
    positions_y += list(sol.y[1][1:])

    times += list(sol.t[1:] + times[-1])

    return n, v


def main():

    # com = sys.argv[1]
    # traj = sys.argv[2]
    com = "fblindiag"
    traj = "3"

    # Simulation settings
    dt = 0.1
    total_time = 100.0
    steps = int(total_time / dt)
    time = 0

    # Initial state
    n = np.array([[0], [4], [np.pi / 2]], dtype=float)  # x, y, Φ
    v = np.array([[1], [0], [0]], dtype=float)  # u, v, r

    # Lists to store the x and y positions and times.
    positions_x = [n[0][0]]
    positions_y = [n[1][0]]
    times = [0]

    state_error_integral = 0
    speed_error_integral = 0
    speed_error = None

    path_points = get_path_points()

    # Simulation
    for step in range(steps):

        target, d_target, dd_target = get_new_target(time, traj, path_points)

        # cmd = command(com, n, v, target, d_target, dd_target)
        cmd, state_error_integral, speed_error_integral, speed_error = command(n, v, target, d_target, state_error_integral, dt, speed_error_integral, speed_error)

        n, v = evolution2(n, v, cmd, dt, positions_x, positions_y, times)

        affichage(target, traj, positions_x, positions_y, n, path_points)

        time = times[-1]

    plt.show()


if __name__ == "__main__":
    main()
