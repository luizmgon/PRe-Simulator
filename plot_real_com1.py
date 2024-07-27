import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import interp1d

# Substitua 'now_str' pelo valor apropriado
# now_str = "2024-07-15_16-52-20-LOSofc"  # Por exemplo
# now_str = "2024-07-15_16-52-20-LOSofc"  # Por exemplo
# now_str = "2024-07-15_16-43-07-Hofc" 
# now_str = "2024-07-15_16-46-49-AUVofc" 
# now_str = "2024-07-15_16-55-43-fblinOFC" 
# now_str = "2024-07-15_16-59-25-SLIDOFC" 

now_str = "2024-07-23_14-57-08-first_H"
# now_str = "2024-07-23_15-25-57-first_AUV_moinsang"
# now_str ="2024-07-23_15-31-57-first_AUV_plusang"
# now_str = "2024-07-23_15-38-55-first_LOS"
# now_str = "2024-07-23_15-49-30-first_FBLIN"
# now_str = "2024-07-23_16-08-30-first_slid"

filename_pose_rob = f"/home/luiz/log_bboat2/{now_str}/pose_rob.txt"
filename_target = f"/home/luiz/log_bboat2/{now_str}/control_target.txt"
filename_path = f"/home/luiz/log_bboat2/{now_str}/path.txt"
filename_command = f"/home/luiz/log_bboat2/{now_str}/command.txt"
filename_speed = f"/home/luiz/log_bboat2/{now_str}/speed.txt"


# Função para carregar e interpolar os dados
def load_and_interpolate(filename, timestamps_command, num_points=10000):
    data = np.loadtxt(filename, delimiter=',')
    timestamps = data[:, 0]
    x_data = data[:, 1]
    y_data = data[:, 2]
    z_data = data[:, 3]

    new_timestamps = np.linspace(timestamps_command.min(), timestamps_command.max(), num=num_points)
    
    interp_x = interp1d(timestamps, x_data, kind='linear', fill_value="extrapolate")
    interp_y = interp1d(timestamps, y_data, kind='linear', fill_value="extrapolate")
    interp_z = interp1d(timestamps, z_data, kind='linear', fill_value="extrapolate")

    x_interp = interp_x(new_timestamps)
    y_interp = interp_y(new_timestamps)
    z_interp = interp_z(new_timestamps)

    return new_timestamps, x_interp, y_interp, z_interp

data_command = np.loadtxt(filename_command, delimiter=',')
timestamps_command = data_command[:, 0]

# Carregar e interpolar dados de ambos os arquivos
timestamps_rob, x_rob, y_rob, z_rob = load_and_interpolate(filename_pose_rob, timestamps_command)
timestamps_target, x_target, y_target, z_target = load_and_interpolate(filename_target, timestamps_command)
timestamps_speed, speed_u, speed_v, speed_r =  load_and_interpolate(filename_speed, timestamps_command)


# Garantir que os novos timestamps sejam os mesmos para ambos os conjuntos de dados
new_timestamps = np.linspace(min(timestamps_rob.min(), timestamps_target.min(), timestamps_speed.min()), 
                             max(timestamps_rob.max(), timestamps_target.max(), timestamps_speed.min()), num=10000)

# Re-interpolar para garantir a mesma linha do tempo
x_rob = interp1d(timestamps_rob, x_rob, kind='linear', fill_value="extrapolate")(new_timestamps)
y_rob = interp1d(timestamps_rob, y_rob, kind='linear', fill_value="extrapolate")(new_timestamps)
z_rob = interp1d(timestamps_rob, z_rob, kind='linear', fill_value="extrapolate")(new_timestamps)

x_target = interp1d(timestamps_target, x_target, kind='linear', fill_value="extrapolate")(new_timestamps)
y_target = interp1d(timestamps_target, y_target, kind='linear', fill_value="extrapolate")(new_timestamps)
z_target = interp1d(timestamps_target, z_target, kind='linear', fill_value="extrapolate")(new_timestamps)

speed_u = interp1d(timestamps_speed, speed_u, kind='linear', fill_value="extrapolate")(new_timestamps)
speed_v = interp1d(timestamps_speed, speed_v, kind='linear', fill_value="extrapolate")(new_timestamps)
speed_r = interp1d(timestamps_speed, speed_r, kind='linear', fill_value="extrapolate")(new_timestamps)


# speed_u = [0]
# speed_r = [0]

errors = np.sqrt((x_rob - x_target)**2 + (y_rob - y_target)**2)

# Carregar os dados do arquivo de comandos

vel_rot = data_command[:, 1]
vel_frente = data_command[:, 2]
vel_rot = interp1d(timestamps_command, vel_rot, kind='previous', fill_value="extrapolate")(new_timestamps)
vel_frente = interp1d(timestamps_command, vel_frente, kind='previous', fill_value="extrapolate")(new_timestamps)
timestamps_command = new_timestamps

path_data = np.loadtxt(filename_path, delimiter=',')
path_x = path_data[0]
path_y = path_data[1]

# Configurar a figura e os eixos
spacex = (max(np.max(x_rob), np.max(path_x)) - min(np.min(x_rob), np.min(path_x))) * 0.1
spacey = (max(np.max(y_rob), np.max(path_y)) - min(np.min(y_rob), np.min(path_y))) * 0.1

fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(10, 12))
ax1.set_ylim(min(np.min(x_rob), np.min(path_x)) - spacex, 
             max(np.max(x_rob), np.max(path_x)) + spacex)
ax1.set_xlim(min(np.min(y_rob), np.min(path_y)) - spacey, 
             max(np.max(y_rob), np.max(y_target), np.max(path_y)) + spacey)

ax1.set_xlabel('X (m)')
ax1.set_ylabel('Y (m)')

# ax2.set_xlim(0, timestamps_command.max() - timestamps_command.min())
# ax2.set_ylim(min(np.min(vel_frente), np.min(speed_u)) - 1, max(np.max(vel_frente), np.max(speed_u)) + 1)
# ax2.set_xlabel('Time (s)')
# ax2.set_ylabel('Velocity input (m/s)')

ax3.set_xlim(0, timestamps_command.max() - timestamps_command.min())
ax3.set_ylim(min(np.min(vel_rot), np.min(speed_u)) - 1, max(np.max(vel_rot), np.min(speed_r)) + 1)
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Velocity input (m/s)')

# ax4.set_xlim(0, timestamps_command.max() - timestamps_command.min())
# ax4.set_ylim(min(errors) - 0.1*(min(errors)), max(errors) + 0.1*max(errors))
# ax4.set_xlabel('Time (s)')
# ax4.set_ylabel('Position Error (m)')

# Inicializar os pontos que serão atualizados
point_rob, = ax1.plot([], [], 'ro')
point_target, = ax1.plot([], [], 'bx', label='Initial position', ms=10, markeredgewidth=3)
trajectory_rob, = ax1.plot([], [], 'r-', label='Followed trajectory')
trajectory_target, = ax1.plot([], [], 'b-')
path_line, = ax1.plot([], [], 'g--', label='Desired Path')

# line_vel_rot, = ax3.plot([], [], 'y-', label='Angular speed')
# line_vel_frente, = ax3.plot([], [], 'g-', label='Forward speed')
line_com_rot, = ax3.plot([], [], 'r-', label='Angular command')
line_com_frente, = ax3.plot([], [], 'b-', label='Forward command')

# error_plot, = ax4.plot([], [], 'b-')

zero, = ax3.plot([-1000, 1000], [0, 0], ':', color='gray')

# Função de inicialização da animação
def init():
    point_rob.set_data([], [])
    point_target.set_data([], [])
    path_line.set_data(path_y, path_x)
    # line_vel_rot.set_data([], [])
    # line_vel_frente.set_data([], [])
    line_com_rot.set_data([], [])
    line_com_frente.set_data([], [])
    # error_plot.set_data([], [])
    return point_rob, point_target, zero, line_com_frente, line_com_rot,

# Função de atualização da animação
def update(frame):
    point_rob.set_data(y_rob[frame], x_rob[frame])
    point_target.set_data(y_rob[0], x_rob[0])
    # point_target.set_data(y_rob[0], x_rob[0])

    # error_plot.set_data(timestamps_command[:frame] - timestamps_command.min(), errors[:frame])
    
    trajectory_rob.set_data(y_rob[:frame], x_rob[:frame])
    # trajectory_target.set_data(y_target[:frame], x_target[:frame])
    
    # Atualizar as linhas de velocidade
    line_com_rot.set_data(timestamps_command[:frame] - timestamps_command.min(), vel_rot[:frame])
    line_com_frente.set_data(timestamps_command[:frame] - timestamps_command.min(), vel_frente[:frame])

    # line_vel_rot.set_data(timestamps_command[:frame] - timestamps_command.min(), speed_r[:frame])
    # line_vel_frente.set_data(timestamps_command[:frame] - timestamps_command.min(), speed_u[:frame])
    
    return point_rob, point_target, trajectory_rob, trajectory_target

# Criar a animação
ani = FuncAnimation(fig, update, frames=len(timestamps_command), init_func=init, blit=False, interval=0, repeat=False)

# Exibir a animação
ax1.legend()
ax3.legend()
# ax2.legend()

plt.show()
