import os
import matplotlib.pyplot as plt
import math as m
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


def exponential(x, a):
    """
    Función exponencial para el ajuste: exp(-a * x)
    
    Args:
        x (array): Los valores de tau.
        a (float): El parámetro a que se va a ajustar.
        
    Returns:
        array: Los valores de la función exponencial.
    """
    return np.exp(-a * x)

def process_dls_files(directory,output_file):
    """
    Procesa todos los archivos de un directorio, grafica los datos normalizados.

    Args:
        directory (str): Ruta del directorio que contiene los archivos DLS.
    """

    fig, ax =plt.subplots(figsize=(10, 6))

    colores=["red","green","blue","orange","purple","cyan","darkblue","pink","chocolate","crimson"]
    # Lista para almacenar las curvas y las etiquetas
    curves = []
    i=0
    with open(output_file, 'w') as result_file:

        for filename in os.listdir(directory):
            if not filename.endswith(".INT"):
                continue

            filepath = os.path.join(directory, filename)

            # Leer el archivo
            with open(filepath, 'r') as file:
                lines = file.readlines()

            # Ignorar la primera línea
            lines = lines[1:]

            # Extraer datos de las columnas y el valor de M
            x_values = []
            y_values = []
            m_value = None

            for line in lines:
                line = line.strip()
                if line.startswith("M") and len(line.split()) > 1:
                    m_value = float(line.split()[1])
                elif line and not line[0].isalpha():
                    parts = line.split()
                    if len(parts) >= 2:
                        x_values.append(float(parts[0]))
                        y_values.append(float(parts[1]))

            if m_value is None:
                print(f"Advertencia: No se encontró un valor 'M' en el archivo {filename}.")
                continue
            
            #Ajustar eje X
            x_valuesnorm = [x*20 for x in x_values]

            # Normalizar los valores de Y
            y_values_normalized = [y / m_value for y in y_values]
            print(filename)
            cte_siegert = y_values_normalized[0]-1

            g1_values= [m.sqrt(abs(y-1)/cte_siegert) for y in y_values_normalized]

            # Graficar los datos
            label = filename[1:].replace(".INT", "")
            try:
                label_number = int(label)  # Extraer el número para ordenar
            except ValueError:
                label_number = 0  # Si no se puede extraer un número, poner 0
            
            # Añadir la curva y la etiqueta para ordenarlas después
            curves.append((x_valuesnorm, g1_values, label_number))

            # Filtrar los valores de x_valuesnorm y g1_values donde x_valuesnorm < 500
            filtered_x = [x for x in x_valuesnorm if x < 500]
            filtered_g1 = [g for x, g in zip(x_valuesnorm, g1_values) if x < 500]

            # Ajustar la curva exponencial
            popt, pcov = curve_fit(exponential, filtered_x, filtered_g1, p0=[0.01])

            # Obtener el parámetro 'a' del ajuste
            a_value = popt[0]
            perr = m.sqrt(np.diag(pcov))
            result_file.write(f"{label_number}\t{a_value:.7f}\t{perr:.7f}\n")
            if label_number==100:
                ax.plot(filtered_x, exponential(np.array(filtered_x), *popt), '--', color = colores[9])
            else:
                ax.plot(filtered_x, exponential(np.array(filtered_x), *popt), '--', color = colores[i])
                i=i+1


    
    # Ordenar las curvas según el número extraído de las etiquetas
    curves.sort(key=lambda x: x[2])

    # Graficar las curvas ordenadas
    i=0
    for curve in curves:
        x_valuesnorm, g1_values, _ = curve
        label = f'{_}'  # Puedes usar cualquier formato para la etiqueta
        ax.plot(x_valuesnorm, g1_values, label=f'{label}°', color = colores[i])
        i=i+1

    subfigure = False

    if subfigure:
        axins = inset_axes(ax, width= 3.5, height= 1.5, loc='upper center', bbox_to_anchor=(0.5, 0.43), bbox_transform=ax.figure.transFigure)        
        # Agregar líneas que conectan la región del zoom con el inset
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="gray")

            # Graficar las curvas ordenadas
        i=0
        for curve in curves:
            x_valuesnorm, g1_values, _ = curve
            label = f'{_}'  # Puedes usar cualquier formato para la etiqueta
            axins.plot(x_valuesnorm, g1_values, label=f'{label}', color = colores[i])
            i=i+1

        axins.set_xlim(0, 500) # apply the x-limits
        axins.set_ylim(0.7, 1) # apply the y-limits
        axins.tick_params(axis="both", which="major")

        axins.set_yscale('log')

    # Configurar la gráfica
    ax.set_xlabel(r'$\tau (\mu s)$',fontsize=14)
    ax.set_ylabel(r'$g_1(\tau)$',fontsize=14)
    ax.set_yscale('log')
    ax.set_xlim(0,500)
    ax.set_ylim(0.7,1)
    ax.legend(fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize= 14)
    ax.set_title(rf"Muestra 8, Duración $40$ s, 256 canales, $\tau= 20$ $\mu s$", fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.7,which="both")

    #plt.show()
    plt.savefig("ajusteangulo.png", format='png', dpi=500, bbox_inches = 'tight')

# Ejemplo de uso:
process_dls_files("angulo","ajustes.txt")
