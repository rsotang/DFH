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

    #print(dose_data.shape) 114 114 200  
    #print(dose_data[:,:,100])  

  
    # Extraer información de la estructura de interés
    for roi in rtstruct.StructureSetROISequence:
        if roi.ROIName == roi_name:
            roi_number = roi.ROINumber
    
    contour_data = [contour for contour in rtstruct.ROIContourSequence if contour.ReferencedROINumber == roi_number][0].ContourSequence
    
    # Aquí debes implementar el cálculo del DVH basado en dose_grid y contour_data
    # El cálculo exacto puede variar dependiendo de los detalles de tu implementación y del archivo RTStruct.
    
    return dvh_data

#funcion temporal para comprobar como funcionan las rois
def show_roi(roi_name):
   # Extraer información de la estructura de interés
    for roi in rtstruct.StructureSetROISequence:
        if roi.ROIName == roi_name:
            roi_number = roi.ROINumber
    
    contour_data = [contour for contour in rtstruct.ROIContourSequence if contour.ReferencedROINumber == roi_number][0].ContourSequence

    print(contour_data.__class__)
    print(contour_data)
    # Aquí debes implementar el cálculo del DVH basado en dose_grid y contour_data
    # El cálculo exacto puede variar dependiendo de los detalles de tu implementación y del archivo RTStruct.
    
    return roi


def plot_dvh(dvh_data):
    plt.plot(dvh_data['doses'], dvh_data['volume_percentages'])
    plt.xlabel('Dose (Gy)')
    plt.ylabel('Volume (%)')
    plt.title('Dose-Volume Histogram')
    plt.show()

if __name__ == "__main__":
    rtdose_path = "rtdose.dcm"
    rtstruct_path = "rtstruct_lungs.dcm"
    
    rtdose, rtstruct = load_dicom_files(rtdose_path, rtstruct_path)
    
    roi_name = "Lung_L"  # Por ejemplo: "Pulmón"
    roi = show_roi(roi_name)
    #dvh_data = create_dvh(rtdose, rtstruct, roi_name)
    
    #plot_dvh(dvh_data)