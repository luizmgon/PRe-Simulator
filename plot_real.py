import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import interp1d

# Substitua 'now_str' pelo valor apropriado
# now_str = "2024-07-15_16-52-20-LOSofc"  # Por exemplo
# now_str = "2024-07-15_16-43-07-Hofc" 
now_str = "2024-07-15_16-46-49-AUVofc" 
# now_str = "2024-07-15_16-55-43-fblinOFC" 
# now_str = "2024-07-15_16-59-25-SLIDOFC" 

name = '0'
filename_pose_rob = f"/home/luiz/log_bboat/{now_str}/pose_rob.txt"
filename_target = f"/home/luiz/log_bboat/{now_str}/control_target.txt"
filename_path = f"/home/luiz/log_bboat/{now_str}/path.txt"


# Função para carregar e interpolar os dados
def load_and_interpolate(filename, num_points=10000):
    data = np.loadtxt(filename, delimiter=',')
    timestamps = data[:, 0]
    x_data = data[:, 1]
    y_data = data[:, 2]
    z_data = data[:, 3]

    new_timestamps = np.linspace(timestamps.min(), timestamps.max(), num=num_points)
    
    interp_x = interp1d(timestamps, x_data, kind='linear', fill_value="extrapolate")
    interp_y = interp1d(timestamps, y_data, kind='linear', fill_value="extrapolate")
    interp_z = interp1d(timestamps, z_data, kind='linear', fill_value="extrapolate")

    x_interp = interp_x(new_timestamps)
    y_interp = interp_y(new_timestamps)
    z_interp = interp_z(new_timestamps)

    return new_timestamps, x_interp, y_interp, z_interp

# Carregar e interpolar dados de ambos os arquivos
timestamps_rob, x_rob, y_rob, z_rob = load_and_interpolate(filename_pose_rob)
timestamps_target, x_target, y_target, z_target = load_and_interpolate(filename_target)

# Garantir que os novos timestamps sejam os mesmos para ambos os conjuntos de dados
new_timestamps = np.linspace(min(timestamps_rob.min(), timestamps_target.min()), 
                             max(timestamps_rob.max(), timestamps_target.max()), num=10000)

# Re-interpolar para garantir a mesma linha do tempo
x_rob = interp1d(timestamps_rob, x_rob, kind='linear', fill_value="extrapolate")(new_timestamps)
y_rob = interp1d(timestamps_rob, y_rob, kind='linear', fill_value="extrapolate")(new_timestamps)
z_rob = interp1d(timestamps_rob, z_rob, kind='linear', fill_value="extrapolate")(new_timestamps)

x_target = interp1d(timestamps_target, x_target, kind='linear', fill_value="extrapolate")(new_timestamps)
y_target = interp1d(timestamps_target, y_target, kind='linear', fill_value="extrapolate")(new_timestamps)
z_target = interp1d(timestamps_target, z_target, kind='linear', fill_value="extrapolate")(new_timestamps)


path_data = np.loadtxt(filename_path, delimiter=',')
path_x = path_data[0]
path_y = path_data[1]

# Configurar a figura e os eixos

spacex = (max(np.max(x_rob), np.max(x_target)) - min(np.min(x_rob), np.min(x_target))) * 0.1
spacey = (max(np.max(y_rob), np.max(y_target)) - min(np.min(y_rob), np.min(y_target))) * 0.1


fig, ax = plt.subplots()
ax.set_ylim(min(np.min(x_rob), np.min(x_target), np.min(path_x)) - spacex, 
            max(np.max(x_rob), np.max(x_target), np.max(path_x)) + spacex)
ax.set_xlim(min(np.min(y_rob), np.min(y_target), np.min(path_y)) - spacey, 
            max(np.max(y_rob), np.max(y_target), np.max(path_y)) + spacey)

ax.set_xlabel('X')
ax.set_ylabel('Y')

# background_image_path = '/home/luiz/PRe/sea'
# background_image = plt.imread(background_image_path)
# ax.imshow(background_image,extent=[np.min(ax.get_xlim()), np.max(ax.get_xlim()), np.min(ax.get_ylim()), np.max(ax.get_ylim())], aspect='auto', alpha = 0.7)


# Inicializar os pontos que serão atualizados
point_rob, = ax.plot([], [], 'ro', label='Pose Rob')
point_target, = ax.plot([], [], 'bo', label='Control Target')
trajectory_rob, = ax.plot([], [], 'r-', label='Trajectory Rob')
trajectory_target, = ax.plot([], [], 'b-', label='Trajectory Target')
path_line, = ax.plot([], [], 'g--', label='Desired Path')


# Função de inicialização da animação
def init():
    point_rob.set_data([], [])
    point_target.set_data([], [])
    path_line.set_data(path_y, path_x)

    return point_rob, point_target

# Função de atualização da animação
def update(frame):
    point_rob.set_data(y_rob[frame], x_rob[frame])
    point_target.set_data(y_target[frame], x_target[frame])
    
    # if frame == len(new_timestamps) - 1:
        # Último frame, plotar trajetórias acumuladas
    trajectory_rob.set_data(y_rob[:frame], x_rob[:frame])
    # trajectory_target.set_data(x_target, y_target)
    
    return point_rob, point_target, trajectory_rob

# Criar a animação
ani = FuncAnimation(fig, update, frames=len(new_timestamps), init_func=init, blit=True, interval=0, repeat=False)

# ax.plot(x_rob, y_rob, 'r-', label='Trajectory Rob')
# ax.plot(x_target, y_target, 'b-', label='Trajectory Target')
# plt.draw()

# Exibir a animação
# plt.legend()
plt.show()
