import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

# Cargar datos
filename = 'orbitaCoordenada/coordenadas_N500_gen100.txt'

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
fig = plt.figure(figsize=(18, 8))
fig.suptitle('Simulación de Trayectoria del Cohete', fontsize=16, y=1.05)

# Configurar subplots y espacio para el slider
gs = fig.add_gridspec(2, 3, height_ratios=[10, 1], hspace=0.6)
ax1 = fig.add_subplot(gs[0, 0])  # Vista general
ax2 = fig.add_subplot(gs[0, 1])  # Sistema Tierra-Luna
ax3 = fig.add_subplot(gs[0, 2])  # Sistema Marte
ax_slider = fig.add_subplot(gs[1, :])  # Slider (ocupa toda la fila inferior)

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

# Dibujar trayectorias completas (en gris claro)
ax1.plot(x_cohete, y_cohete, '-', color='lightgray', linewidth=0.5, alpha=0.5)
ax2.plot(x_cohete_tierra, y_cohete_tierra, '-', color='lightgray', linewidth=0.5, alpha=0.5)
ax3.plot(x_cohete_marte, y_cohete_marte, '-', color='lightgray', linewidth=0.5, alpha=0.5)

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

# Texto del tiempo
time_text = fig.text(0.5, 0.95, '', ha='center', fontsize=12, 
                     bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray'))

ax1.legend(loc='upper right', fontsize=8)

# Crear slider para controlar el tiempo
time_slider = Slider(
    ax=ax_slider,
    label='Tiempo (horas)',
    valmin=t[0],
    valmax=t[-1],
    valinit=t[0],
    valstep=np.diff(t)[0]  # Usar el paso de tiempo del archivo
)

# Función para actualizar la visualización según el tiempo seleccionado
def update(val):
    # Encontrar el índice más cercano al tiempo seleccionado
    idx = np.abs(t - time_slider.val).argmin()
    
    # Actualizar vista general
    tray_cohete1.set_data(x_cohete[:idx+1], y_cohete[:idx+1])
    punto_tierra1.set_data([x_tierra[idx]], [y_tierra[idx]])
    punto_marte1.set_data([x_marte[idx]], [y_marte[idx]])
    punto_luna1.set_data([x_lunaSol[idx]], [y_lunaSol[idx]])
    punto_cohete1.set_data([x_cohete[idx]], [y_cohete[idx]])
    
    # Actualizar vista Tierra-Luna
    tray_cohete2.set_data(x_cohete_tierra[:idx+1], y_cohete_tierra[:idx+1])
    punto_luna2.set_data([x_luna_tierra[idx]], [y_luna_tierra[idx]])
    punto_cohete2.set_data([x_cohete_tierra[idx]], [y_cohete_tierra[idx]])
    
    # Actualizar vista Marte
    tray_cohete3.set_data(x_cohete_marte[:idx+1], y_cohete_marte[:idx+1])
    punto_cohete3.set_data([x_cohete_marte[idx]], [y_cohete_marte[idx]])
    
    time_text.set_text(f'Tiempo: {t[idx]:.2f} horas')
    fig.canvas.draw_idle()

# Registrar la función de actualización con el slider
time_slider.on_changed(update)

# Botones de control
ax_prev = plt.axes([0.15, 0.02, 0.05, 0.03])
ax_next = plt.axes([0.21, 0.02, 0.05, 0.03])
ax_play = plt.axes([0.27, 0.02, 0.05, 0.03])

btn_prev = Button(ax_prev, '<')
btn_next = Button(ax_next, '>')
btn_play = Button(ax_play, 'Play')

# Variables para la animación automática
is_playing = False
current_idx = 0

# Funciones para los botones
def step_prev(event):
    new_val = max(time_slider.val - time_slider.valstep, time_slider.valmin)
    time_slider.set_val(new_val)

def step_next(event):
    new_val = min(time_slider.val + time_slider.valstep, time_slider.valmax)
    time_slider.set_val(new_val)

def toggle_play(event):
    global is_playing
    is_playing = not is_playing
    if is_playing:
        btn_play.label.set_text('Pause')
        play_animation()
    else:
        btn_play.label.set_text('Play')
        fig.canvas.draw_idle()

def play_animation():
    if is_playing:
        step_next(None)
        plt.pause(0.05)  # Pequeña pausa para controlar la velocidad
        fig.canvas.start_event_loop(0.05)  # Para permitir la interacción durante la reproducción
        fig.canvas.flush_events()
        play_animation()  # Llamada recursiva

# Conectar los botones a sus funciones
btn_prev.on_clicked(step_prev)
btn_next.on_clicked(step_next)
btn_play.on_clicked(toggle_play)

# Mostrar la posición inicial
update(t[0])

plt.tight_layout()
plt.show()