import numpy as np
import simulacion2
import random
import pandas as pd
from multiprocessing import Pool
import os

# Variables generales del programa
N = [50, 100, 5000, 1000, 2000, 5000]
gen = [50, 100, 200, 500, 1000]
#N= [2000]
#gen = [1000]

numPadres = 20
probMutacion = 0.3


def fitness(x, y):
    fitnesDistancia = 1/(1 + x)  # Fitness basado en la distancia al objetivo
    fitnessEnergia = 1/(1 + y)  # Fitness basado en la energía

    return 0.7 * fitnesDistancia + 0.3 * fitnessEnergia  # Combinación ponderada de ambos fitness




# Constantes físicas
G = 6.67430e-11           # Constante gravitacional en Nm^2/kg^2
Ms = 1.989e30             # Masa del Sol en kg
Mt = 5.972e24             # Masa de la Tierra en kg
Ml = 7.34767309e22        # Masa de la Luna en kg
Mm = 6.4171e23            # Masa de Marte en kg
m = 2000.0                # Masa del objeto en kg
L = 1.496e11              # Distancia en metros del Sol a la Tierra
# Radios reescalados
radioTierra = 6.371e6 / L  # Radio reescalado de la Tierra
radioMarte = 3.3895e6 / L  # Radio reescalado de Marte
radioLuna = 1.7374e6 / L   # Radio reescalado de la Luna

Mtotal = Ms + Mt + Ml + Mm  # Masa total del sistema
tCaract = np.sqrt(L**3 / (G * Mtotal))  # Tiempo característico del sistema

# Velocidades angulares
wt = tCaract * 2.0 * np.pi / 31557600.0  # Velocidad angular de la Tierra alrededor del Sol
wm = tCaract * 2.0 * np.pi / 59353000.0  # Velocidad angular de Marte alrededor del Sol
wl = tCaract * 2.0 * np.pi / 2360592.0   # Velocidad angular de la Luna alrededor de la Tierra

# Distancias orbitales reescaladas
dOrbTierra = 1.0               # Distancia orbital de la Tierra al Sol
dOrbMarte = 2.279e11 / L       # Distancia orbital de Marte al Sol
dOrbLuna = 3.844e8 / L         # Distancia orbital de la Luna a la Tierra





# Definición de los límites de las variables
vMin = 3700*tCaract/L
vMax = 4500*tCaract/L

tMax = 2*30*24*3600/tCaract  # 1 mes en segundos reescalado

print(dir(simulacion2))
# Iteramos sobre cada valor de N y gen
for n in N:
    for g_max in gen:
        # Definicion de las variables
        magnitudes = np.zeros((n, 2))  # Matriz para almacenar magnitudes de cada individuo

        # Genero la primera poblacion
        for i in range(n):
            magnitudes[i, 0] = random.uniform(vMin, vMax)  # Velocidad inicial
            magnitudes[i, 1] = random.uniform(0, tMax)     # Tiempo de impulso

        print(f"Simulación para N = {n}, generaciones = {g_max}")
        # Generamos el bucle principal
        def evaluar_individuo(i, magnitudes, dOrbLuna, dOrbTierra, dOrbMarte, fitness):
            res = np.zeros(6)  # Array para almacenar el resultado de la simulación
            simulacion2.simulacion2.ejecutarsimulacion(magnitudes[i, 1], magnitudes[i, 0], res)
        
            # Pongo las excepciones
            if res[4] < dOrbLuna * 1.5:
                return 1.0e-6
            elif res[3] < dOrbTierra * 0.9:
                return 1.0e-6
            elif res[3] > dOrbMarte * 1.2:
                return 1.0e-6
            elif res[5] > 1.0:
                return 1.0
            else:
                return fitness(res[0], magnitudes[i, 0])  # Calcular el fitness del individuo

        #print(f"Número de núcleos disponibles: {os.cpu_count()}")
        for g in range(g_max):
            # Evaluamos la poblacion
            fitnesses = np.zeros(n)

            with Pool(processes=os.cpu_count()) as pool:  # Use all available CPU cores
                args = [(i, magnitudes, dOrbLuna, dOrbTierra, dOrbMarte, fitness) for i in range(n)]
                results = []
                for arg in args:
                    result = pool.apply_async(evaluar_individuo, arg)
                    results.append(result)

                fitnesses = []
                for result in results:
                    try:
                        fitnesses.append(result.get(timeout=10))  # Set a timeout of 10 seconds
                    except Exception as e:
                        print(f"Error al evaluar individuo: {e}, Índice: {arg[0]}, Magnitudes: {magnitudes[arg[0]]}")
                        fitnesses.append(1.0e-6)  # Assign a default fitness value if an error occurs
            # Seleccionamos los mejores padres
            padres_indices = np.argsort(fitnesses)[-numPadres:]  # Indices de los mejores fitness
            padres_magnitudes = magnitudes[padres_indices]       # Magnitudes correspondientes a los mejores padres

            # Guardamos los datos en un archivo CSV
            best_index = np.argmax(fitnesses)  # Índice del mejor individuo
            data = {
                "N": [n],
                "Generacion": [g],
                "Fitness": [fitnesses[best_index]],
                "Magnitud1": [magnitudes[best_index, 0]],
                "Magnitud2": [magnitudes[best_index, 1]]
            }
            df = pd.DataFrame(data)
            print(df)
            filename = f"evolucion/mejores_resultados_N{n}_gen{g_max}.csv"
            if g == 0:
                df.to_csv(filename, index=False, mode='w')  # Crear archivo en la primera generación
            else:
                df.to_csv(filename, index=False, mode='a', header=False)  # Añadir datos en generaciones posteriores
            # Generamos la nueva poblacion
            nueva_poblacion = np.zeros_like(magnitudes)
            # Generamos nuevos individuos a partir de los padres seleccionados
            for b in range(n - numPadres):
                num1 = random.randint(0, numPadres - 1)
                num2 = random.randint(0, numPadres - 1)
                nueva_poblacion[b, 0] = 0.5 * (padres_magnitudes[num1, 0] + padres_magnitudes[num2, 0])
                nueva_poblacion[b, 1] = 0.5 * (padres_magnitudes[num1, 1] + padres_magnitudes[num2, 1])
                # Aplicamos la mutación
                if random.random() < probMutacion:  # Probabilidad de mutación
                    nueva_poblacion[b, 0] += np.random.uniform(-0.3, 0.3) * (vMax - vMin)
                    nueva_poblacion[b, 0] = np.clip(nueva_poblacion[b, 0], vMin, vMax)
                if random.random() < probMutacion:  # Probabilidad de mutación
                    nueva_poblacion[b, 1] += 5 * 24 * 3600 * random.uniform(-1, 1) / tCaract
                    nueva_poblacion[b, 1] = np.clip(nueva_poblacion[b, 1], 0.0, tMax)
            # Rellenamos la población nueva con los padres seleccionados
            nueva_poblacion[n - numPadres:] = padres_magnitudes
            # Verificamos si hay un fitness de 1.0
            if 1.0 in fitnesses:
                print("Fitness óptimo alcanzado. Terminando simulación.")
                break
            # Actualizamos la población para la siguiente generación
            magnitudes = nueva_poblacion


        # Al finalizar, mostramos el mejor resultado
        best_index = np.argmax(fitnesses)
        print(f"Mejor resultado encontrado para N = {n}, generaciones = {g_max}: Velocidad = {magnitudes[best_index, 0]}, Tiempo = {magnitudes[best_index, 1]}, Fitness = {fitnesses[best_index]}")
        res = simulacion2.simulacion2.representarsimulacion(magnitudes[best_index, 1], magnitudes[best_index, 0], n, g_max)
        print(f"Datos de la mejor simulación: {simulacion2.simulacion2.representarorbita(magnitudes[best_index, 1], magnitudes[best_index, 0], n, g_max, res[0])}")
