import matplotlib
import numpy as np
from dicompylercore import dicomparser, dvh, dvhcalc
import pylab as pl
import sys

#para calcular un dvh necesito la informacion anatomica y dosimetrica. La información anatomica y geométrica la puedo sacar del archivo
#de estructuras y la dosimetrica de cada estrucutra compaginandola con el archivo de dosimetría

#Primero importamos archivos y sacamos las estructuras en diccionarios usando la libreria dicompoarser. Este objeto es una clase exclusiva de dicomparser, por eso podemos despues usar el metodo .getstructures()

rtssfile = 'rtstruct.dcm'
rtdosefile = 'rtdose.dcm'
RTss = dicomparser.DicomParser(rtssfile)
RTstructures = RTss.GetStructures()

#print(dir(RTstructures))
#print(dir(RTstructures.items()))
#print(RTstructures.items())
#ahora vamos a utilizar las funciones de la libreria para calcular los histogramas, son algunas de estas funciones las que tendremos que modificar para poder hacer el DFH
#creamos un diccionario abierto y vamos insertando items en el
calcdvhs = {}
#calcdvhs[1] =  dvhcalc.get_dvh(rtssfile, rtdosefile, 1)
#print(len(calcdvhs[1].counts))
#sys.exit()

#for key, structure in RTstructures.items(): #iteramos dentro de cada elemento del diccionario de estructuras (key siendo el orden y estructuras #su diccionario (con todas sus caracteristicas en el array))
#    calcdvhs[key] = dvhcalc.get_dvh(rtssfile, rtdosefile, key) #asignamos el histograma calculado al diccionario, la funcion get_dvh() busca las #estructuras y las compagina con las dosis para calcularlo.
#  
#    if (key in calcdvhs) and (len(calcdvhs[key].counts) and calcdvhs[key].counts[0]!=0): #comprobacion de que la estructura no está vacia
#        print ('DVH found for ' + structure['name'])

key_lungR = [clave for clave, valor in RTstructures.items() if valor['name'] == 'Lung_R']
key_lungL = [clave for clave, valor in RTstructures.items() if valor['name'] == 'Lung_L']

print("Claves con name 'LUNG':", key_lungR)
print("Claves con name 'LUNG':", key_lungL)


structure_lungR = RTstructures[key_lungR[0]]
structure_lungL = RTstructures[key_lungL[0]]

print("estructuras 'LUNG':", structure_lungR)
print("estructuras 'LUNG':", structure_lungL)


calcdvhs[key_lungR[0]] = dvhcalc.get_dvh(rtssfile, rtdosefile, key_lungR[0])
calcdvhs[key_lungL[0]] = dvhcalc.get_dvh(rtssfile, rtdosefile, key_lungL[0])

pl.plot(calcdvhs[key_lungR[0]].counts * 100 / calcdvhs[key_lungL[0]].counts[0], 
         color=dvhcalc.np.array(structure_lungR['color'], dtype=float) / 255,
         label=structure_lungR['name'],
         linestyle='dashed')

pl.plot(calcdvhs[key_lungL[0]].counts * 100 / calcdvhs[key_lungL[0]].counts[0], 
         color=dvhcalc.np.array(structure_lungL['color'], dtype=float) / 255,
         label=structure_lungL['name'],
         linestyle='dashed')     
        
               

pl.xlabel('Dose (cGy)')
pl.ylabel('Percentage Volume')
pl.legend(loc=7, borderaxespad=-5)
pl.setp(pl.gca().get_legend().get_texts(), fontsize='x-small')
pl.savefig('histogramas/dvh.png', dpi = 75)