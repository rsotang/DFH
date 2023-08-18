import os
import sys
import inspect

import numpy as np
import pydicom
import matplotlib.pyplot as plt

def load_dicom_files(rtdose_path, rtstruct_path):
    rtdose = pydicom.dcmread(rtdose_path)
    rtstruct = pydicom.dcmread(rtstruct_path)
    return rtdose, rtstruct

def create_dvh(rtdose, rtstruct, roi_name):
    # Extraer información de dosis
    dose_data = np.array(rtdose.pixel_array)
    dose_grid = dose_data * rtdose.DoseGridScaling

    print(length(dose_data))
  
    clase_objeto = dose_grid.__class__
    metodos_clase = inspect.getmembers(clase_objeto, predicate=inspect.ismethod)

    # Obtener todos los atributos de la clase
    atributos_clase = inspect.getmembers(clase_objeto, predicate=inspect.isattribute)

    print(clase_objeto)

  
    # Extraer información de la estructura de interés
    for roi in rtstruct.StructureSetROISequence:
        if roi.ROIName == roi_name:
            roi_number = roi.ROINumber
    
    contour_data = [contour for contour in rtstruct.ROIContourSequence if contour.ReferencedROINumber == roi_number][0].ContourSequence
    
    # Aquí debes implementar el cálculo del DVH basado en dose_grid y contour_data
    # El cálculo exacto puede variar dependiendo de los detalles de tu implementación y del archivo RTStruct.
    
    return dvh_data

def plot_dvh(dvh_data):
    plt.plot(dvh_data['doses'], dvh_data['volume_percentages'])
    plt.xlabel('Dose (Gy)')
    plt.ylabel('Volume (%)')
    plt.title('Dose-Volume Histogram')
    plt.show()

if __name__ == "__main__":
    rtdose_path = "rtdose.dcm"
    rtstruct_path = "rtstruct.dcm"
    
    rtdose, rtstruct = load_dicom_files(rtdose_path, rtstruct_path)
    
    roi_name = "Lung_L"  # Por ejemplo: "Pulmón"
    
    dvh_data = create_dvh(rtdose, rtstruct, roi_name)
    
    #plot_dvh(dvh_data)