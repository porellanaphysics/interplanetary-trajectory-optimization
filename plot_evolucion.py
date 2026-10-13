import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import colormaps
import numpy as np
import matplotlib.colors as mcolors

# Path to the 'evolucion' folder
evolucion_path = './evolucion/'

# Get all CSV files in the folder
csv_files = [f for f in os.listdir(evolucion_path) if f.endswith('.csv')]

# Extract unique N values from filenames
unique_N_values = sorted(set([
    int(f.split('_')[2][1:]) for f in csv_files if 'mejores_resultados_N' in f
]))

# Define base colors for each N value
base_colors = {
    2000: 'blue',
    1000: 'red',
    500: 'green',
    250: 'orange',
    100: 'purple'
}

# Generate a colormap for each base color to vary tones
N_to_color = {}
for N, base_color in base_colors.items():
    colormap = mcolors.LinearSegmentedColormap.from_list(
        f"{base_color}_shades", [mcolors.to_rgba(base_color, alpha) for alpha in [0.5, 1.0]]
    )
    N_to_color[N] = [colormap(i / 10) for i in range(1, 11)]  # Generate 10 shades for each base color

# Initialize the plot
plt.figure(figsize=(10, 6))

# Loop through each CSV file and plot its data
for idx, csv_file in enumerate(csv_files):
    if 'mejores_resultados_N' not in csv_file:
        continue  # Skip files that don't match the expected format

    file_path = os.path.join(evolucion_path, csv_file)
    data = pd.read_csv(file_path)

    # Adjusted to match the correct column names in the CSV
    generations = data['Generacion'] + 1  # Shift generations so 0 becomes 1
    fitness = data['Fitness']

    # Extract N value from the filename for color assignment
    try:
        N_value = int(csv_file.split('_')[2][1:])
    except (IndexError, ValueError):
        print(f"Skipping file with unexpected format: {csv_file}")
        continue

    # Assign a slightly varied color for each dataset with the same N value
    color_variants = N_to_color.get(N_value, ['black'])
    color = color_variants[idx % len(color_variants)]  # Cycle through shades if more datasets than shades

    # Plot the data with a label
    label = csv_file.replace('mejores_resultados_', '').replace('.csv', '').replace('_', ', ')
    plt.plot(generations, fitness, label=label, color=color)

# Add labels, title, and legend
plt.xlabel('Generation (log scale)')
plt.ylabel('Fitness (log scale)')
plt.title('Generación vs Fitness para diferentes N y gen')

# Reapply logarithmic scale
plt.xscale('log')
plt.yscale('log')

# Sort legend entries alphabetically
handles, labels = plt.gca().get_legend_handles_labels()
sorted_legend = sorted(zip(labels, handles), key=lambda x: x[0])
labels, handles = zip(*sorted_legend)
plt.legend(handles, labels, title='Configuraciones', loc='upper left', bbox_to_anchor=(1, 1))

plt.grid(True)

# Format axes to show normal numbers instead of scientific notation
plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter())
plt.gca().xaxis.set_minor_formatter(ticker.ScalarFormatter())
plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter())
plt.gca().yaxis.set_minor_formatter(ticker.ScalarFormatter())
plt.gca().ticklabel_format(style='plain', axis='both')

# Show the plot
plt.tight_layout()
plt.show()

# Cargar los datos del archivo CSV
data = pd.read_csv('energias_orbitas.csv')

# Calcular las velocidades a partir de las energías
data['Velocidad'] = data['Energia'] / 2  # Ejemplo de cálculo, ajustar según la fórmula correcta

# Agregar la transferencia de Hohmann
hohmann_energy = -0.16659057757377105
hohmann_velocity = hohmann_energy / 2  # Ejemplo de cálculo, ajustar según la fórmula correcta

# Crear el gráfico de barras
plt.figure(figsize=(10, 6))
plt.bar(data['Generacion'], data['Velocidad'], label='Velocidades por Generación')
plt.axhline(y=hohmann_velocity, color='r', linestyle='--', label='Transferencia de Hohmann')

# Configurar el gráfico
plt.title('Velocidades por Generación y Transferencia de Hohmann')
plt.xlabel('Generación')
plt.ylabel('Velocidad')
plt.legend()

# Mostrar el gráfico
plt.show()