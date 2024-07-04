import numpy as np
import matplotlib.pyplot as plt
from commands import *

def get_target_line(time):
    # final_position = np.array([[5], [time], [np.pi/2]], dtype=float)  # Final position is a straight line with slope of 2
    # first_derivative = np.array([[0], [1], [0]], dtype=float)  # First derivative (constant velocity)
    # second_derivative = np.array([[0], [0], [0]], dtype=float)  # Second derivative (constant acceleration)
    # return final_position, first_derivative, second_derivative

    final_position = np.array([[5], [time], [np.pi/2]], dtype=float)  # Final position is a straight line with slope of 2
    first_derivative = np.array([[0], [1], [0]], dtype=float)  # First derivative (constant velocity)
    second_derivative = np.array([[0], [0], [0]], dtype=float)  # Second derivative (constant acceleration)
    return final_position, first_derivative, second_derivative

def get_target_path(time, path_points):
    ds = 2
    s = ds * time
    xs, ys, _, phif, curv, g_c, dx, dy, ddx, ddy = path_interrogation(s, path_points)
    final_position = np.array([[xs], [ys], [0]], dtype=float)
    first_derivative = np.array([[ds * dx], [ds * dy], [0]], dtype=float)
    second_derivative = np.array([[ds * ddx], [ds * ddy], [0]], dtype=float)
    return final_position, first_derivative, second_derivative


def get_target_circle(time):
    
    R = 2  # Radius of the circle
    cx, cy = 5, 5  # Center of the circle
    omega = 0.6 # Angular velocity

    x = cx + R * np.cos(omega * time)
    y = cy + R * np.sin(omega * time)
    dx = -R * omega * np.sin(omega * time)
    dy = R * omega * np.cos(omega * time)
    ddx = -R * omega**2 * np.cos(omega * time)
    ddy = -R * omega**2 * np.sin(omega * time)

    position = np.array([[x], [y], [0]], dtype=float)
    first_derivative = np.array([[dx], [dy], [0]], dtype=float)
    second_derivative = np.array([[ddx], [ddy], [0]], dtype=float)

    return position, first_derivative, second_derivative

def get_chemin_target_line(target, n_state):
    x,y = n_state[:2]
    xd, yd = target[:2]

    distance = np.sqrt((xd - x)**2 + (yd - y)**2)

    if(distance < 2):
        target[1] += 5
    
    return target

def get_chemin_target_circle(n, target, n_state):
    x,y = n_state[:2]
    xd, yd = target[:2]

    distance = np.sqrt((xd - x)**2 + (yd - y)**2)

    if(distance < 0.5):
        circle_x = 5 + 2 * np.cos(0.5*n)
        circle_y = 5 + 2 * np.sin(0.5*n)
        print(n)
        target[0] = circle_x
        target[1] = circle_y
        n += 1
    
    return n, target

def draw_traj_line():
    plt.plot([5, 5], [-10, 100], color='green', linestyle='--')  # Trajectory line from (5,-10) to (5,20)

def draw_traj_circle():
    # Calculando a trajetória circular
    circle_theta = np.linspace(0, 2 * np.pi, 100)
    circle_x = 5 + 2 * np.cos(circle_theta)
    circle_y = 5 + 2 * np.sin(circle_theta)
    
    # Desenhando a trajetória circular
    plt.plot(circle_x, circle_y, color='green', linestyle='--')

def draw_traj_path(path_points, axs):
    x_values, y_values = path_points[:2]
    axs.plot(x_values, y_values, color="green", linestyle="--")
