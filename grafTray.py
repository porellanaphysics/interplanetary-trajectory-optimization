import numpy as np
import matplotlib.pyplot as plt

# Cargar datos
filename = 'coordenadas.txt'

def cargar_datos(filename):
    datos = []
    with open(filename, 'r') as f:
        for linea in f:
            partes = linea.strip().split()
            if len(partes) == 11:
                datos.append([float(p.replace('E', 'e')) for p in partes])
    return np.array(datos)

data = cargar_datos(filename)

# Asignar variables
t = data[:, 0]
x_cohete, y_cohete = data[:, 1], data[:, 2]
x_tierra, y_tierra = data[:, 3], data[:, 4]
x_marte, y_marte = data[:, 5], data[:, 6]
x_lunaSol, y_lunaSol = data[:, 7], data[:, 8]
x_lunaTierra, y_lunaTierra = data[:, 9], data[:, 10]

# Calcular posiciones relativas
x_cohete_tierra = x_cohete - x_tierra
y_cohete_tierra = y_cohete - y_tierra
x_luna_tierra = x_lunaSol - x_tierra
y_luna_tierra = y_lunaSol - y_tierra

x_cohete_marte = x_cohete - x_marte
y_cohete_marte = y_cohete - y_marte

# Crear figura
fig = plt.figure(figsize=(18, 6))
fig.suptitle('Simulación de Trayectoria del Cohete', fontsize=16, y=1.05)

# Configurar subplots
ax1 = fig.add_subplot(131)  # Vista general
ax2 = fig.add_subplot(132)  # Sistema Tierra-Luna
ax3 = fig.add_subplot(133)  # Sistema Marte

# Colores consistentes
COLOR_SOL = 'gold'
COLOR_TIERRA = 'blue'
COLOR_MARTE = 'red'
COLOR_LUNA = 'black'
COLOR_COHETE = 'magenta'
COLOR_TRAYECTORIA = 'green'

# 1. Vista General (Sistema Solar centrado en el Sol)
ax1.set_title('Sistema Solar (centrado en Sol)')
ax1.set_xlim(-1.8, 1.8)
ax1.set_ylim(-1.8, 1.8)
ax1.grid(True, alpha=0.3)
ax1.set_aspect('equal')

# Dibujar elementos estáticos
ax1.plot(0, 0, 'o', color=COLOR_SOL, markersize=12, label='Sol')
ax1.plot(x_tierra, y_tierra, '--', color=COLOR_TIERRA, linewidth=0.7, alpha=0.4, label='Órbita Tierra')
ax1.plot(x_marte, y_marte, '--', color=COLOR_MARTE, linewidth=0.7, alpha=0.4, label='Órbita Marte')
ax1.plot(x_lunaSol, y_lunaSol, '--', color=COLOR_LUNA, linewidth=0.5, alpha=0.3, label='Órbita Luna')
ax1.plot(x_cohete, y_cohete, '-', color=COLOR_TRAYECTORIA, linewidth=1.5, label='Trayectoria Cohete')
ax1.legend(loc='upper right', fontsize=8)

# ---- Sistema Tierra-Luna (ax2) ----
distancia_max_tierra_luna = np.max(np.sqrt(x_luna_tierra**2 + y_luna_tierra**2))
limite_tierra = max(distancia_max_tierra_luna, 6.371e6 / 1.496e11) * 1.2

ax2.set_title('Sistema Tierra-Luna')
ax2.set_xlim(-limite_tierra, limite_tierra)
ax2.set_ylim(-limite_tierra, limite_tierra)
ax2.grid(True, alpha=0.3)
ax2.set_aspect('equal')
ax2.plot(0, 0, 'o', color=COLOR_TIERRA, markersize=10, label='Tierra')
ax2.plot(x_luna_tierra, y_luna_tierra, '--', color=COLOR_LUNA, linewidth=0.5, alpha=0.3, label='Órbita Luna')
ax2.plot(x_cohete_tierra, y_cohete_tierra, '-', color=COLOR_TRAYECTORIA, linewidth=1.5, label='Trayectoria Cohete')
ax2.legend(loc='upper right', fontsize=8)

# ---- Sistema Marte (ax3) ----
limite_marte = np.max(np.sqrt(x_cohete_marte**2 + y_cohete_marte**2)) * 1.2

ax3.set_title('Sistema Marte')
ax3.set_xlim(-limite_marte, limite_marte)
ax3.set_ylim(-limite_marte, limite_marte)
ax3.grid(True, alpha=0.3)
ax3.set_aspect('equal')
ax3.plot(0, 0, 'o', color=COLOR_MARTE, markersize=10, label='Marte')
ax3.plot(x_cohete_marte, y_cohete_marte, '-', color=COLOR_TRAYECTORIA, linewidth=1.5, label='Trayectoria Cohete')
ax3.legend(loc='upper right', fontsize=8)

plt.tight_layout()
plt.show()