import numpy as np
import pandas as pd
import simulacion2

# Cargar datos de energías
filename = 'energias_orbitas.csv'
data = pd.read_csv(filename)

# Verificar las columnas disponibles en el archivo CSV
print("Columnas disponibles en el archivo:", data.columns)

# Filtrar las filas con Fitness igual a 1.0 y orbita igual a 1
filtered_data = data[(data['Fitness'] == 1.0) & (data['orbita'] == -1)]

# Identificar el índice de la mejor energía entre las filas filtradas
best_index = filtered_data['Energia'].idxmin()  # Usar el nombre correcto de la columna 'Energia'

# Obtener los parámetros correspondientes
mejor_energia = data.iloc[best_index]
velocidad = mejor_energia['Velocidad']  # Suponiendo que hay una columna 'Velocidad'
tiempo = mejor_energia['Tiempo']        # Suponiendo que hay una columna 'Tiempo'
n = mejor_energia['N']                  # Suponiendo que hay una columna 'N'
g_max = mejor_energia['Generacion']     # Suponiendo que hay una columna 'Generacion'
res = [0]                               # Inicializar res con un valor por defecto

# Llamar a la función representarorbita
#simulacion2.simulacion2.representarorbita(velocidad, tiempo, n, g_max, res[0])
# Imprimir las características de la mejor coordenada
print("Características de la mejor coordenada:")
print(f"Velocidad: {velocidad}")
print(f"Tiempo: {tiempo}")
print(f"N: {n}")
print(f"Generación máxima: {g_max}")
print(f"Energía: {mejor_energia['Energia']}")  # Usar el nombre correcto de la columna 'Energia'