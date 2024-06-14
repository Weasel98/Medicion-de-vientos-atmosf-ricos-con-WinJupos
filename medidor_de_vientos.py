#Este programa lee la tabla de mediciones completa exportada en .txt por WinJupos para calcular la velocidad del objeto que se haya multipuntuado
nombre='L'

import math
import numpy as np
from sklearn.metrics import r2_score
from scipy.stats import linregress
import uncertainties as unc
from uncertainties.umath import sin
from uncertainties.umath import cos
from uncertainties import ufloat, umath
import uncertainties.unumpy as unp

#Se define la función que calcula la distancia entre un punto y una recta. Así saber qué punto es el más lejado de la regresión lineal para posteriormente ser eliminado
def distancia_punto_a_recta(punto, recta):
    # Extraer coordenadas del punto
    x, y = punto

    # Coeficientes de la ecuación de la recta
    a, b, c = recta

    # Calcular la proyección del punto sobre la recta
    x_proyectado = (b * (b * x - a * y) - a * c) / (a**2 + b**2)
    y_proyectado = (a * (-b * x + a * y) - b * c) / (a**2 + b**2)

    # Extraer los valores nominales de los puntos proyectados
    x_proyectado_nominal = x_proyectado.nominal_value #Hace falta .nominal_value porque x_proyectado es un valor con error asociado
    y_proyectado_nominal = y_proyectado.nominal_value

    # Calcular la distancia entre el punto y su proyección sobre la recta
    distancia = np.sqrt((x - x_proyectado_nominal)**2 + (y - y_proyectado_nominal)**2)

    return distancia

def calcular_rmsd(medidas):
    # Calcular la media de las medidas
    media = np.mean(medidas)
    
    # Calcular la suma de los cuadrados de las diferencias
    suma_cuadrados_dif = np.sum((medidas - media) ** 2)
    
    # Calcular la desviación cuadrática media (RMSD)
    rmsd = np.sqrt(suma_cuadrados_dif / len(medidas))
    
    return rmsd

def calcular_rmsd_un(medidas):
    # Calcular la media de las medidas
    media = np.mean(medidas.nominal_value)
    
    # Calcular la suma de los cuadrados de las diferencias
    suma_cuadrados_dif = np.sum((medidas.nominal_value - media) ** 2)
    
    # Calcular la desviación cuadrática media (RMSD)
    rmsd = np.sqrt(suma_cuadrados_dif / len(medidas))
    
    return rmsd

grupos = {} # Diccionario para almacenar las filas agrupadas por los dos últimos dígitos
pi=math.pi
Rp=6.6854*10**7 #Radio polar de Júpiter en metros 
Re=7.1492*10**7 #Radio ecuatorial de Júpiter en metros
t_rotacion=9*60*60+55*60+30 #rotación completa del planeta en segundos

# Abre el archivo .txt de mediciones exportado directamente de WinJupos en modo lectura
#¡_¡_¡_¡_¡_¡_¡_¡ CAMBIAR LA RUTA EN FUNCIÓN DEL DOCUMENTO QUE SE NECESITE LEER !_!_!_!_!_!_!_!_!
with open("D:\\Documents\\Universidad\\TFG\\Programas\\Python\\medicion_"+nombre+"_mea.txt", "r") as archivo:
    # Lee todas las líneas del archivo, excluyendo la primera (en esta se encuantra los nombres de las columnas, no los números que necesitamos)
    lineas = archivo.readlines()[1:]

    # Procesa cada línea
    for linea in lineas:
        # Reemplaza las comas por puntos y divide la línea en una lista de elementos
        elementos_linea = linea.replace(',', '.').replace('+', '').split()

        # Obtiene los dos últimos dígitos de la primera columna, que serán los que se utilicen para identificar los mismos puntos entre distintas imágenes
        dos_ultimos_digitos = elementos_linea[3][-2:]

        # Selecciona los elementos específicos 4, 7, 10, 12 y 13: nombre del punto, fecha de captura en día juliano, longitud, fase, latitud
        elementos_seleccionados = [elementos_linea[i] for i in [3, 6, 9, 11, 12]]

        # Convierte cada elemento a números (enteros o decimales) si es posible, de lo contrario, guarda el elemento como está
        fila = [float(elemento) if elemento.replace('.', '', 1).isdigit() else elemento for elemento in elementos_seleccionados]

        # Agrupa las filas en el diccionario según los dos últimos dígitos. De esta manera, cada item del diccionario es el conjunto de medidas de cada punto
        if dos_ultimos_digitos in grupos:
            grupos[dos_ultimos_digitos].append(fila)
        else:
            grupos[dos_ultimos_digitos] = [fila]

v_media=unc.ufloat(0.0, 0.0)
uves_teta=[]
uves_lambda=[]
tetas_medias=[]
error_min_lambda=error_min_teta=0.2 #grados/segundo. Ireducible por la definición 
datos=[]

#En este ciclo for se agrupa punto por punto en listas separadas las variables a utilizar en los cálculos
for key, value in grupos.items():
    t_array=[]          #lista de tiempos en días jacobianos en los que fueron capturadas las imágenes
    teta_array=[]       #lista de latitudes de los puntos 
    lambda_array=[]     #lista de longitudes de los puntos

    for fila in value:              #ciclo en el que se generan las listas de coordenadas: tiempo, latitud y longitud
        tiempo=fila[1]*24*60*60     #pasa de días jacobianos a segundos
        t_array.append(tiempo)
        lambda_array.append(fila[2])
        teta_array.append(fila[4])

    teta_medio=np.mean(teta_array)
    radio=(Rp*Re)/(Re**2*(sin(teta_medio))**2+Rp**2*(cos(teta_medio))**2)**0.5
    error_min_v_teta=(pi/180)*radio*error_min_teta/t_rotacion
    error_min_v_lambda=-(pi/180)*radio*cos(teta_medio)*error_min_lambda/t_rotacion
    error_permitido=np.sqrt(error_min_v_teta**2+error_min_v_lambda**2)                   #error sobredimensionado de partida para ir ajustándose al error mínimo por definición
    error_velocidad=error_permitido+1.0
    while error_velocidad>error_permitido:  #ciclo que se repite hasta llegar al error mínimo por definición

        v_teta, intercept_teta, r_value, p_value, std_err_teta = linregress(t_array, teta_array)            #regresión para calcular la velocidad latitudinal con error
        v_teta=unc.ufloat(v_teta, std_err_teta)             #velocidad latitudinal con error en grados/segundo

        v_lambda, intercept_lambda, r_value, p_value, std_err_lambda = linregress(t_array, lambda_array)    #regresión para calcular la velocidad longitudinal con error
        v_lambda=unc.ufloat(v_lambda, std_err_lambda)       #velocidad longitudinal con error en grados/segundo

        #como no todos los puntos tienen la misma latitud se hace un promedio junto con la desviación cuadrática media como error: teta_medio
        teta_medio=np.mean(teta_array)
        err_teta=calcular_rmsd(teta_array)
        teta_medio=unc.ufloat(math.radians(teta_medio), math.radians(err_teta))

        radio=(Rp*Re)/(Re**2*(sin(teta_medio))**2+Rp**2*(cos(teta_medio))**2)**0.5 #radio aproximado en el que se encuentran los puntos en función del promedio de latitud
        u=-(pi/180)*radio*cos(teta_medio)*v_lambda #velocidad longitudinal en m/s con su error propagado de la regresión
        v=(pi/180)*radio*v_teta                     #velocidad latitudinal en m/s con su error propagado de la regresión
        v_tot=(u**2+v**2)**0.5                      #velocidad total en m/s con su error propagado de las dos velocidades anteriores
        error_velocidad=v_tot.std_dev #v_tot.std_dev es el error propagado para la velocidad media y el que se compara con el máximo aceptable
        
        if error_velocidad==0.0: #esto significaría que la regresión no ha sido buena y hemos terminado con tan solo dos puntos, de ahí que su error sea 0
            error_velocidad=error_permitido #asignamos el valor mínimo permitido por la definición al error de la medida
            print("Una de las regresiones se ha quedado con, tan solo, 2 puntos.") #avisamos de que un regresión ha contado con pocos puntos

        #se crean dos listas vacías donde almanecar las distancias a las que se encuentra cada punto respecto a la regresión, para saber cual es el punto más alejado y eliminarlo para reducir el error
        d_teta_array=[]
        d_lambda_array=[]
        if error_velocidad>error_permitido:
            for i in range(len(t_array)):
                punto_teta=(t_array[i], teta_array[i])
                recta_teta=(-v_teta,1,-intercept_teta)
                punto_lambda=(t_array[i], lambda_array[i])
                recta_lambda=(-v_lambda,1,-intercept_lambda)
                distancia_teta=distancia_punto_a_recta(punto_teta, recta_teta)
                distancia_lambda=distancia_punto_a_recta(punto_lambda, recta_lambda)
                d_teta_array.append(distancia_teta)
                d_lambda_array.append(distancia_lambda)
            if max(d_teta_array)>max(d_lambda_array):
                del t_array[d_teta_array.index(max(d_teta_array))]
                del teta_array[d_teta_array.index(max(d_teta_array))]
                del lambda_array[d_teta_array.index(max(d_teta_array))]
            if max(d_teta_array)<max(d_lambda_array):
                del t_array[d_lambda_array.index(max(d_lambda_array))]
                del teta_array[d_lambda_array.index(max(d_lambda_array))]
                del lambda_array[d_lambda_array.index(max(d_lambda_array))]

    uves_lambda.append(u)
    uves_teta.append(v)
    tetas_medias.append(umath.degrees(teta_medio))
    v_media = (u**2 + v**2)**0.5
    print(v_media, " es la velocidad media del punto ", key)
    añadir=[key,u.nominal_value,u.std_dev,v.nominal_value,v.std_dev, umath.degrees(teta_medio.nominal_value),umath.degrees(teta_medio.std_dev)]
    datos.append(añadir)

velocidad_media_lambda_prop=np.mean(uves_lambda)
velocidad_media_teta_prop=np.mean(uves_teta)

print("La velocidad longitudinal media es de: ", velocidad_media_lambda_prop, "+-", calcular_rmsd(unp.nominal_values(uves_lambda)))
print("La velocidad latitudinal media  es de: ", velocidad_media_teta_prop, "+-", calcular_rmsd(unp.nominal_values(uves_teta)))

# Supongamos que tienes los nombres de las columnas en una lista llamada 'nombres_columnas'
nombres_columnas = ['Punto', 'u', 'err_u', 'v', 'err_v', 'teta_medio', 'err_teta_medio']

# Define el ancho de cada columna manualmente
ancho_columnas = [20, 15, 15, 15, 15, 15, 15]  # Por ejemplo, 7 para 'Punto', 20 para el resto

# Crea el formato de cadena para el encabezado con el ancho de cada columna
formato_encabezado = ''.join(['{:<{}}'.format(nombre, ancho) for nombre, ancho in zip(nombres_columnas, ancho_columnas)])

# Convertir datos a un array de tipo numérico
datos_numerico = np.array(datos).astype(float)

# Guardar los datos en un archivo de texto
np.savetxt('D:\\Documents\\Universidad\\TFG\\Programas\\Python\\datos_'+nombre+'.txt', datos_numerico, delimiter='\t', header=formato_encabezado, comments='', fmt='%f')


datos_finales=['Hot_Spot_'+nombre,velocidad_media_lambda_prop.nominal_value,calcular_rmsd(unp.nominal_values(uves_lambda)),velocidad_media_teta_prop.nominal_value,calcular_rmsd(unp.nominal_values(uves_teta)), np.mean(unp.nominal_values(tetas_medias))]
# Abre el archivo en modo de escritura (append)

# Define los datos finales como una cadena formateada
linea_datos_finales = '\t'.join(map(str, datos_finales)) + '\n'

# Abre el archivo en modo de escritura (append)
with open('D:\\Documents\\Universidad\\TFG\\Programas\\Python\\datos_finales.txt', 'a') as archivo:
    # Escribe la nueva línea en el archivo
    archivo.write(linea_datos_finales)