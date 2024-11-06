# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 20:51:21 2024

@author: pacor
"""
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 20:51:21 2024

@author: pacor
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress

def load_full_data(rutafile):
    """
    Carga un archivo CSV y devuelve un DataFrame con los datos.

    Parámetros
    ----------
    rutafile (str): Ruta al archivo CSV.

    Devuelve
    --------
    pd.DataFrame: Un DataFrame con los datos del archivo CSV, donde las columnas corresponden
                  a los nombres del encabezado y cada fila representa una entrada de datos.
    """
    datos = pd.read_csv(rutafile, encoding="utf8")
    return datos


def graficar_ajustar_escalalog(datos):
    # Obtener el máximo de la segunda columna
    max_val = datos.iloc[:, 1].max()
    
    # Calcular logaritmos en las dos columnas, guardo el log(W) en la tercera
    datos.iloc[:,0] = datos.iloc[:,0].map(lambda x: np.log10(x))
    datos['log_W'] = datos.iloc[:,1].map(lambda x: np.log10(max_val/x))
    
    # Separar en dos DataFrames
    min_val = datos['log_W'].min()
    
    datos_near_max = datos[datos['log_W'] - min_val < 0.1].copy()
    datos_others = datos[datos['log_W'] >= 0.1].copy()
    
    # Ajuste lineal para cada DataFrame
    m1, b1, *_ = linregress(datos_near_max.iloc[:, 0], datos_near_max['log_W'])
    m2, b2, *_ = linregress(datos_others.iloc[:, 0], datos_others['log_W'])
    
    # Calcular el punto de corte
    x = (b2 - b1) / (m1 - m2)
    y = m1 * x + b1
    
    # Calcula la CCC
    CCC = 10**x
    
    # Graficar
    plt.figure(figsize=(12, 6))
    plt.scatter(datos_near_max.iloc[:, 0], datos_near_max['log_W'], color='blue')
    plt.axline(xy1=(0, b1), slope=m1, label=f'$y = {m1:.1f}x {b1:+.1f}$', color='green')
    plt.scatter(datos_others.iloc[:, 0], datos_others['log_W'], color='blue')
    plt.axline(xy1=(0, b2), slope=m2, label=f'$y = {m2:.1f}x {b2:+.1f}$', color='red')
    plt.plot(x, y, color='red', marker='o', markersize=10, label=f'CCC = {CCC:1f} mM')
    plt.xlabel("log(C)")
    plt.ylabel("log(W)")
    plt.legend()
    plt.xlim(datos_others.iloc[:, 0].min()-0.2, datos_near_max.iloc[:, 0].max()+0.1)
    plt.ylim(datos_near_max['log_W'].min()-0.1, datos_others['log_W'].max()+0.2)
    plt.show()
    
    return datos_near_max, datos_others, datos[['Concentracion','log_W']]


def mostrar_tabla(datos):
    """
    Muestra una tabla con los datos originales y logaritmos calculados como imagen.

    Parámetros
    ----------
    datos : pd.DataFrame
        DataFrame con las columnas 'log_C', 'W' y 'log_W' para mostrar en la tabla.
    """
    fig, ax = plt.subplots(figsize=(8, len(datos) * 0.3 + 1))
    ax.axis('tight')
    ax.axis('off')
    tabla = ax.table(cellText=datos.values.round(2), colLabels=['log(C)','log(W)'], cellLoc='center', loc='center')
    tabla.auto_set_font_size(False)
    tabla.set_fontsize(10)
    tabla.scale(0.5,1.2)
    plt.show()


# Ejecución
datos = load_full_data("SistBNaSCN.csv")
resultados1, resultados2, datos_transformados= graficar_ajustar_escalalog(datos)
mostrar_tabla(datos_transformados)
