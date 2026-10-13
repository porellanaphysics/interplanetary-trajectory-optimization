import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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


# ---- Radios planetarios en UA ----
radioTierra = 6.371e6 / 1.496e11  # Radio de la Tierra en UA
radioMarte = 3.3895e6 / 1.496e11   # Radio de Marte en UA
radioLuna = 1.7374e6 / 1.496e11    # Radio de la Luna en UA

# ---- Sistema Tierra-Luna (ax2) ----
distancia_max_tierra_luna = np.max(np.sqrt(x_luna_tierra**2 + y_luna_tierra**2))
limite_tierra = max(distancia_max_tierra_luna, radioTierra) * 1.2

ax2.set_xlim(-limite_tierra, limite_tierra)
ax2.set_ylim(-limite_tierra, limite_tierra)
circulo_tierra = plt.Circle((0, 0), radioTierra, color=COLOR_TIERRA, alpha=0.7)
ax2.add_patch(circulo_tierra)
ax2.set_aspect('equal')
ax2.grid(True, alpha=0.3)


# ---- Sistema Marte (ax3) ----
distancia_max_tierra_luna = np.max(np.sqrt(x_luna_tierra**2 + y_luna_tierra**2))
limite_marte = distancia_max_tierra_luna * 1.2  # Usamos la distancia Tierra-Luna como referencia

ax3.set_title('Vista centrada en Marte (escala Tierra-Luna)')
ax3.set_xlim(-limite_marte, limite_marte)
ax3.set_ylim(-limite_marte, limite_marte)
ax3.grid(True, alpha=0.3)
ax3.set_aspect('equal')

# Dibujar Marte (radio amplificado para visibilidad)
radioMarte_amplificado = radioMarte * 1  # Factor de amplificación
circulo_marte = plt.Circle((0, 0), radioMarte_amplificado, color=COLOR_MARTE, alpha=0.7)
ax3.add_patch(circulo_marte)

# Inicializar elementos dinámicos
tray_cohete1, = ax1.plot([], [], '-', color=COLOR_TRAYECTORIA, linewidth=1.5, label='Cohete')
punto_tierra1, = ax1.plot([], [], 'o', color=COLOR_TIERRA, markersize=6)
punto_marte1, = ax1.plot([], [], 'o', color=COLOR_MARTE, markersize=6)
punto_luna1, = ax1.plot([], [], 'o', color=COLOR_LUNA, markersize=3)
punto_cohete1, = ax1.plot([], [], 'o', color=COLOR_COHETE, markersize=5)

tray_cohete2, = ax2.plot([], [], '-', color=COLOR_TRAYECTORIA, linewidth=1.5)
punto_tierra2, = ax2.plot([0], [0], 'o', color=COLOR_TIERRA, markersize=10)
punto_luna2, = ax2.plot([], [], 'o', color=COLOR_LUNA, markersize=6)
punto_cohete2, = ax2.plot([], [], 'o', color=COLOR_COHETE, markersize=6)

tray_cohete3, = ax3.plot([], [], '-', color=COLOR_TRAYECTORIA, linewidth=1.5)
punto_marte3, = ax3.plot([0], [0], 'o', color=COLOR_MARTE, markersize=10)
punto_cohete3, = ax3.plot([], [], 'o', color=COLOR_COHETE, markersize=6)

# Inicializar elementos dinámicos (cohete relativo a Marte)
tray_cohete3, = ax3.plot([], [], '-', color=COLOR_TRAYECTORIA, linewidth=1.5)
punto_cohete3, = ax3.plot([], [], 'o', color=COLOR_COHETE, markersize=6)

# Texto del tiempo
time_text = fig.text(0.5, 0.95, '', ha='center', fontsize=12, 
                     bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray'))

ax1.legend(loc='upper right', fontsize=8)

# Función de inicialización
def init():
    for element in [tray_cohete1, tray_cohete2, tray_cohete3,
                   punto_tierra1, punto_marte1, punto_luna1, punto_cohete1,
                   punto_luna2, punto_cohete2, punto_cohete3]:
        element.set_data([], [])
    time_text.set_text('')
    return (tray_cohete1, tray_cohete2, tray_cohete3,
            punto_tierra1, punto_marte1, punto_luna1, punto_cohete1,
            punto_luna2, punto_cohete2, punto_cohete3,
            time_text)

# Función de animación
def animate(i):
    tray_cohete1.set_data(x_cohete[:i+1], y_cohete[:i+1])
    punto_tierra1.set_data([x_tierra[i]], [y_tierra[i]])
    punto_marte1.set_data([x_marte[i]], [y_marte[i]])
    punto_luna1.set_data([x_lunaSol[i]], [y_lunaSol[i]])
    punto_cohete1.set_data([x_cohete[i]], [y_cohete[i]])

    tray_cohete2.set_data(x_cohete_tierra[:i+1], y_cohete_tierra[:i+1])
    punto_luna2.set_data([x_luna_tierra[i]], [y_luna_tierra[i]])
    punto_cohete2.set_data([x_cohete_tierra[i]], [y_cohete_tierra[i]])

    tray_cohete3.set_data(x_cohete_marte[:i+1], y_cohete_marte[:i+1])
    punto_cohete3.set_data([x_cohete_marte[i]], [y_cohete_marte[i]])

    time_text.set_text(f'Tiempo: {t[i]:.2f} minutos')
    return (tray_cohete1, tray_cohete2, tray_cohete3,
            punto_tierra1, punto_marte1, punto_luna1, punto_cohete1,
            punto_luna2, punto_cohete2, punto_cohete3,
            time_text)

# Aplicar tight_layout antes de la animación
plt.tight_layout()

# Crear animación (con blit=False para evitar errores)
ani = FuncAnimation(fig, animate, frames=len(t), init_func=init,
                    interval=50, blit=False, repeat=False)

plt.show()