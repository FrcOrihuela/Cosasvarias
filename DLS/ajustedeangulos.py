import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def process_file(input_file):
    """
    Lee un archivo con dos columnas, calcula los valores de la primera columna 
    a partir de la fórmula dada, ajusta los valores de la segunda columna 
    mediante regresión lineal (ordenada en el origen nula) y guarda los resultados.
    
    Args:
        input_file (str): Nombre del archivo de entrada con dos columnas.
        output_file (str): Nombre del archivo de salida con los resultados.
    """
    # Parámetros de la fórmula
    factor = (4 * np.pi * 1.33) / 632
    k = 1.3806488E-23
    T=293
    eta=1.002E-3


    # Leer el archivo de entrada
    data = np.loadtxt(input_file)
    
    # Primera columna (valores de x)
    x_data = data[:, 0]
    
    # Segunda columna (valores observados y_data)
    y_data = data[:, 1]
    y_err = data[:, 2]

    # Calcular la primera columna usando la fórmula (4pi1.33/632*sin(x/2))^2
    q = [(factor * np.sin(0.5*np.radians(x)))**2 for x in x_data]
    print([a for a in q])
    # Regresión lineal sin ordenada en el origen utilizando scipy
    # Ajuste de los datos para obtener la pendiente m
    slope, intercept, r_value, p_value, std_err = stats.linregress(q, y_data)

    # La ordenada en el origen debe ser cero, por lo que no la usamos
    D_0 = slope
    a = k*T/(6*np.pi*eta*D_0*1E-12)*1E9
    err_a = a*std_err/D_0
    print(f"Pendiente de la regresión lineal: {D_0:.6f}")
    print(f"R-squared: {r_value:.3f}")

    # Graficar los resultados
    plt.figure(figsize=(8, 6))
    plt.errorbar(q, y_data, yerr=y_err, color='blue', fmt='.', label='Datos experimentales')
    
    # Dibujar la línea de regresión
    y_fit = [D_0 * a + intercept for a in q]
    plt.plot(q, y_fit, color='green', label='Ajuste lineal')
    texto = fr'$D_0$ = ({D_0:.2f} $\pm$ {std_err:.2f}) nm$^2$/$\mu$s'f'\n'fr'$a = ({a:.0f} \pm {err_a: .0f})$ nm'
    plt.text(0.025, 0.7, texto, transform= plt.gca().transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='left', bbox=dict(boxstyle="round,pad=0.3", facecolor="white",edgecolor="grey", alpha=0.7))

    
    plt.xlabel(r'$q^2$ (m$^{-2}$)',fontsize=14)
    plt.ylabel(r'$\Gamma$ ($\mu$s$^{-1}$)',fontsize=14)
    plt.legend(fontsize=14)
    plt.tick_params(axis="both", which="major", labelsize= 14)
    plt.xticks(np.arange(0, 0.0005, step=0.0001))
    plt.grid(True, linestyle='--', alpha=0.7,which="both")
    #plt.show()
    plt.savefig("ajusteradio.png", format='png', dpi=500, bbox_inches = 'tight')

# Ejemplo de uso:
input_file = "ajustes.txt"  # Nombre del archivo de entrada
process_file(input_file)
