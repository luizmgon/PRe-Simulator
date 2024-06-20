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

def model2(t, y, cmd):
    n_state = y[:3].reshape((3, 1))
    v_state = y[3:].reshape((3, 1))

    dn = J(n_state) @ v_state
    dv = np.linalg.inv(M_model) @ (cmd - C(v_state) @ v_state - D(v_state) @ v_state)

    dydt = np.vstack((dn, dv)).flatten()
    return dydt


# Draw the chosen trajectory.
def draw_traj(x_values, y_values, axs):
    axs.plot(x_values, y_values, color="green", linestyle="--")


# Update the state and velocity over time based on the command.
def evolution(n, v, cmd, dt_ctr):
    steps = 100

    for _ in range(steps):
        dn, dv = model(n, v, cmd)
        v += dv * dt_ctr / steps
        n += dn * dt_ctr / steps

    n[2] = sawtooth(n[2])

    return n, v

def evolution2(n0, v0, cmd, dt_ctr, positions_x, positions_y, times, errors, xs, ys):

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
    error = mean(list(np.sqrt((sol.y[0][1:] - xs)**2 + (sol.y[1][1:] - ys)**2)))
    errors.append(error)

    # times += list(sol.t[1:] + times[-1])
    times.append(sol.t[-1] + times[-1])

    return n, v


# Show the current state, target, and trajectory on a plot.
def affichage(n, error_state, path_points):
    plt.clf()

    s1, y1, phi, s = error_state.flatten()
    # print(s)
    # xs, ys, phif, curv, dcurv_dt, g_c = path_description(s, ds, s_values, mu_values, phi_values)

    x_values, y_values = path_points[:2]
    xs, ys, _, phif = path_interrogation(s, path_points)[:4]

    # x_limit = [-20, 10]
    # y_limit = [-5, 35]
    x_limit = [-75, 10]
    y_limit = [-5, 60]

    axs = plt.subplot(111)
    draw_traj(x_values, y_values, axs)

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

def affichage_complet(positions_x, positions_y, n, path_points, times, errors, error_lin_speed, error_ang_speed, taus_lin, taus_ang):
    # Criar uma figura com 4 subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs1 = plt.subplot(221)
    axs2 = plt.subplot(222)
    axs3 = plt.subplot(425)
    axs4 = plt.subplot(427)
    axs5 = plt.subplot(426)
    axs6 = plt.subplot(428)
    # x_limit = [-20, 10]
    # y_limit = [-5, 35]
    x_limit = [-75, 10]
    y_limit = [-5, 60]

    x_values, y_values = path_points[:2]

    # Primeiro gráfico
    axs1.plot(positions_x, positions_y, color="blue")
    axs1.set_xlim(x_limit)
    axs1.set_ylim(y_limit)
    axs1.set_title("Trajetória")
    axs1.grid(True)
    # Função para desenhar a trajetória
    draw_traj(x_values, y_values, axs1)  # Substitua pela sua função de desenho
    # draw_tank(n)  # Substitua pela sua função de desenho do tanque

    # Segundo gráfico
    axs2.plot(times[1:], errors[1:], color="blue")
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


# Main function to run the simulation.
def main():

    path_points = get_path_points()
    x_values, y_values, s_values, mu_values, phi_values = draw_p()

    # Simulation settings.
    dt_ctr = 0.01
    total_time = 20
    steps = int(total_time / dt_ctr)

    # Initial state (position and orientation).
    n = np.array([[-30], [50], [-np.pi]], dtype=float)  # x, y, phiU
    v = np.array([[2], [0], [0]], dtype=float)  # u, v, r

    s0 = 1
    u_target = 2
    error_state = initiate_error_state(n, v, s0, path_points )  # s1, y1, phi, s

    positions_x = [n[0][0]]
    positions_y = [n[1][0]]
    errors = [0]
    error_lin_speed = []
    error_ang_speed = []
    times = [0]
    taus_lin = []
    taus_ang = []

    # Run the simulation for each time step.
    for step in range(steps):

        # affichage(n, error_state, path_points)

        cmd, ds, xs, ys, r_error, u_error, dr, du, y1 = command_auv_model(v, error_state, u_target, path_points, s_values, mu_values, phi_values)

        error_lin_speed.append(y1)
        error_ang_speed.append(r_error)
        taus_lin.append(du)
        taus_ang.append(dr)

        n, v = evolution2(n, v, cmd, dt_ctr, positions_x, positions_y, times, errors, xs, ys)

        error_state = update_error_state(error_state, v, n, ds, dt_ctr, path_points)

    affichage_complet(positions_x, positions_y, n, path_points, times, errors, error_lin_speed, error_ang_speed, taus_lin, taus_ang)


if __name__ == "__main__":
    main()
