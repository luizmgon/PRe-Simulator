import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import interp1d


# Dados de exemplo (substitua pelos seus dados reais)
timestamps = np.linspace(0, 10, num=5)
new_timestamps = np.linspace(0,10, num = 40)
dados = [1, 0.5, 3, 2, 1]  # Dados de exemplo (substitua pelos seus dados reais)

# Configurar a figura e os eixos
fig, ax = plt.subplots()

vel_rot = interp1d(timestamps, dados, kind='previous', fill_value="extrapolate")(new_timestamps)

# Plotar os dados como onda quadrada
ax.plot(new_timestamps, vel_rot, label='Onda Quadrada')

# Configurações adicionais do gráfico
ax.set_xlabel('Tempo')
ax.set_ylabel('Amplitude')
ax.set_title('Exemplo de Plot de Onda Quadrada')
ax.legend()

# Exibir o gráfico
plt.show()