import os

# Directorio donde se encuentran los archivos
directorio = ''  # Cambia esto al directorio adecuado

directorio_absoluto = os.path.abspath(directorio)
print(directorio_absoluto) 

# Iterar sobre los nombres de archivo y borrar si cumplen con la condición
for numero in range(1, 316):
    nombre_archivo = f"ExportImg_{numero}.DCM"
    ruta_archivo = os.path.join(directorio_absoluto, nombre_archivo)
    
    if os.path.exists(ruta_archivo):
        os.remove(ruta_archivo)
        print(f"Archivo {nombre_archivo} eliminado.")
    else:
        print(f"Archivo {nombre_archivo} no encontrado en la ruta {ruta_archivo}.")