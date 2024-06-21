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
    


def draw_traj(traj, path_points, axs):
    if traj == "l":
        return draw_traj_line()
    elif traj == 'c':
        return draw_traj_circle()
    else: return draw_traj_path(path_points, axs)


# Show the current state, target, and trajectory on a plot.
def affichage(target, traj, positions_x, positions_y, n, path_points ,ne_state):
        
        plt.clf()

        x_limit = [-10, 5]
        y_limit = [-1, 32]
        axs = plt.subplot(111)
        draw_traj(traj, path_points, axs)

        plt.scatter(target[0], target[1], color="red", marker="o")
        plt.scatter(ne_state[0], ne_state[1], color="yellow", marker="x")
        plt.plot(positions_x, positions_y, color="blue")
        plt.xlim(x_limit)
        plt.ylim(y_limit)
        draw_tank(n)

        plt.draw()
        plt.pause(0.001)

def affichage_complet(target, traj, positions_x, positions_y, n, path_points, times, errors, error_lin_speed, error_ang_speed, taus_lin, taus_ang):
    # Criar uma figura com 4 subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs1 = plt.subplot(221)
    axs2 = plt.subplot(222)
    axs3 = plt.subplot(425)
    axs4 = plt.subplot(427)
    axs5 = plt.subplot(426)
    axs6 = plt.subplot(428)
    x_limit = [-10, 5]
    y_limit = [-1, 32]

    # Primeiro gráfico
    axs1.plot(positions_x, positions_y, color="blue")
    axs1.scatter(target[0], target[1], color="red", marker="o")
    axs1.set_xlim(x_limit)
    axs1.set_ylim(y_limit)
    axs1.set_title("Trajetória")
    axs1.grid(True)
    # Função para desenhar a trajetória
    draw_traj(traj, path_points, axs1)  # Substitua pela sua função de desenho
    # draw_tank(n)  # Substitua pela sua função de desenho do tanque

    # Segundo gráfico
    axs2.plot(times, errors, color="blue")
    axs2.plot([-10, 100], [errors[-1], errors[-1]], color="green", linestyle='--')
    axs2.set_xlim([times[0], times[-1]])
    axs2.set_title("Erros ao Longo do Tempo")
    axs2.grid(True)

    # Terceiro gráfico
    axs3.plot(range(len(taus_lin)), taus_lin)
    axs3.set_title("Linear acceleration")
    axs3.grid(True)

    axs4.plot(range(len(taus_ang)), taus_ang)
    axs4.set_title("Angular acceleration")
    axs4.grid(True)

    # Quarto gráfico
    axs5.plot(range(len(error_lin_speed)), error_lin_speed)
    axs5.set_title("Linear speed error")
    axs5.grid(True)

    axs6.plot(range(len(error_ang_speed)), error_ang_speed)
    axs6.set_title("Angular speed error")
    axs6.grid(True)

    # Ajustar layout para não sobrepor títulos e eixos
    plt.tight_layout()

    # Mostrar a figura
    plt.show()

def evolution2(n0, v0, cmd, dt_ctr, positions_x, positions_y, times, errors, target):

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
    error = mean(list(np.sqrt((sol.y[0][1:] - target[0])**2 + (sol.y[1][1:] - target[1])**2)))
    errors.append(error)

    # times += list(sol.t[1:] + times[-1])
    times.append(sol.t[-1] + times[-1])

    return n, v


def main():

    # com = sys.argv[1]
    # traj = sys.argv[2]
    com = "fblindiag"
    traj = "3"

    # Simulation settings
    dt = 0.02
    total_time = 10.0
    steps = int(total_time / dt)
    time = 0

    # Initial state
    n = np.array([[4], [10], [np.pi / 2]], dtype=float)  # x, y, Φ
    v = np.array([[1], [0], [0]], dtype=float)  # u, v, r

    # Lists to store the x and y positions and times.
    positions_x = [n[0][0]]
    positions_y = [n[1][0]]
    errors = [0]
    error_lin_speed = []
    error_ang_speed = []
    times = [0]
    taus_lin = []
    taus_ang = []

    state_error_integral = 0
    speed_error_integral = 0
    state_error = None

    path_points = get_path_points()

    # Simulation
    for step in range(steps):


        target, d_target, dd_target = get_new_target(time, traj, path_points)



        # cmd = command(com, n, v, target, d_target, dd_target)
        cmd, state_error_integral, speed_error_integral, state_error, ac, er_sp, ne_state = command_h(n, v, target, d_target, dd_target, state_error_integral, dt, speed_error_integral, state_error)

        error_lin_speed.append(er_sp[0])
        error_ang_speed.append(er_sp[2])
        taus_lin.append(ac[0])
        taus_ang.append(ac[2])

        n, v = evolution2(n, v, cmd, dt, positions_x, positions_y, times, errors, target)

        time = times[-1]
        # print(time)

        affichage(target, traj, positions_x, positions_y, n, path_points, ne_state)


    affichage_complet(target, traj, positions_x, positions_y, n, path_points, times, errors, error_lin_speed, error_ang_speed, taus_lin, taus_ang)


if __name__ == "__main__":
    main()
