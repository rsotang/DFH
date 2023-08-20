import os
import sys
import inspect

import numpy as np
import pydicom
from skimage.draw import polygon
from dicom_mask.convert import struct_to_mask

def get_mask_from_rtstruct(rtstruct_dcm, rtdose_dcm, organ_name):
    # Encuentra el índice de la estructura deseada
    

    return mask

rtstruct_dcm = pydicom.dcmread("rtstruct_lungs.dcm")
rtdose_dcm = pydicom.dcmread("rtdose.dcm")

#i

# Obtener máscara binaria del órgano deseado
organ_name = "Lung_L"
mask = get_mask_from_rtstruct(rtstruct_dcm, rtdose_dcm.pixel_array.shape, organ_name)

# Enmascarar el archivo RTDOSE con la máscara binaria
masked_dose = mask * rtdose_dcm.pixel_array