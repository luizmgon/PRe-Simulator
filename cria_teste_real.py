import numpy as np
import os
import time

# Defina o diretório e os nomes dos arquivos
now_str = "2024-07-15"  # Por exemplo
name = '0'
dir_path = f"/home/luiz/log_bboat/{now_str}"
if not os.path.exists(dir_path):
    os.makedirs(dir_path)
filename_pose_rob = os.path.join(dir_path, f"pose_rob{name}.txt")
filename_target = os.path.join(dir_path, f"control_target{name}.txt")
filename_path = os.path.join(dir_path, f"path{name}.txt")

# Função para gerar dados de teste
def generate_test_data(filename, raio, kind, num_points=10000, freq=0.0001):
    timestamps = np.linspace(0, num_points, num=num_points)
    x_data = np.sin(2 * np.pi * freq * timestamps) * raio  # Variação senoidal em x
    y_data = np.cos(2 * np.pi * freq * timestamps) * raio  # Variação cosenoidal em y
    z_data = np.linspace(0, 5, num=num_points)  # Variação linear em z

    if kind == 1:
        data = np.column_stack((timestamps, x_data, z_data, z_data))
    else:
        data = np.column_stack((timestamps, y_data, z_data, z_data))
    
    np.savetxt(filename, data, delimiter=',', fmt='%.6f')

# Gerar dados de teste para os dois arquivos
generate_test_data(filename_pose_rob, 10, 1)
generate_test_data(filename_target, 6, 2)

# Gerar dados para o arquivo path.txt (trajetória desejada)
path_x = np.sin(np.linspace(0, 2*np.pi, num=100)) * 8  # Exemplo de trajetória desejada em x
path_y = np.cos(np.linspace(0, 2*np.pi, num=100)) * 4  # Exemplo de trajetória desejada em y
path_data = np.vstack((path_x, path_y))
np.savetxt(filename_path, path_data, delimiter=',', fmt='%.6f')

print(f"Arquivos de teste criados:\n{filename_pose_rob}\n{filename_target}\n{filename_path}")
