#!/usr/bin/env python
# -*- coding: utf-8 -*-
# dvhcalc.py
#"""Calculate dose volume histogram (DVH) from DICOM RT Structure/Dose data."""
# Copyright (c) 2016 gluce
# Copyright (c) 2011-2016 Aditya Panchal
# Copyright (c) 2010 Roy Keyes
# This file is part of dicompyler-core, released under a BSD license.
#    See the file license.txt included with this distribution, also
#    available at https://github.com/dicompyler/dicompyler-core/

from __future__ import division

import logging
import sys
import os

import matplotlib.path
import numpy as np
import numpy.ma as ma
import pylab as pl
from dicompylercore import dicomparser, dvh, dvhcalc
from six import iteritems

logger = logging.getLogger('dicompylercore.dvhcalc')


def get_dvh(structure, dose, roi, limit=None, callback=None):
#    """Calculate a cumulative DVH in Gy from a DICOM RT Structure Set & Dose.
#
#    Parameters
#    ----------
#    structure : pydicom Dataset
#        DICOM RT Structure Set used to determine the structure data.
#    dose : pydicom Dataset
#        DICOM RT Dose used to determine the dose grid.
#    roi : int
#        The ROI number used to uniquely identify the structure in the structure
#        set.
#    limit : int, optional
#        Dose limit in cGy as a maximum bin for the histogram.
#    callback : function, optional
#        A function that will be called at every iteration of the calculation.
#     """
  
    rtss = dicomparser.DicomParser(structure)
    rtdose = dicomparser.DicomParser(dose)
    structures = rtss.GetStructures()
    s = structures[roi]
    s['planes'] = rtss.GetStructureCoordinates(roi)
    s['thickness'] = rtss.CalculatePlaneThickness(s['planes'])
    hist = calculate_dvh(s, rtdose, limit, callback)
    return dvh.DVH(counts=hist,
                   bins=(np.arange(0, 2) if (hist.size == 1) else
                         np.arange(0, hist.size + 1) / 100),
                   dvh_type='differential',
                   dose_units='gy',
                   name=s['name']
                   ).cumulative


def calculate_dvh(structure, dose, limit=None, callback=None):
#    """Calculate the differential DVH for the given structure and dose grid.
#
#    Parameters
#    ----------
#    structure : dict
#        A structure (ROI) from an RT Structure Set parsed using DicomParser
#    dose : DicomParser
#        A DicomParser instance of an RT Dose
#    limit : int, optional
#        Dose limit in cGy as a maximum bin for the histogram.
#    callback : function, optional
#        A function that will be called at every iteration of the calculation.
#    """
    planes = structure['planes']
    logger.debug("Calculating DVH of %s %s", structure['id'], structure['name'])

    # Create an empty array of bins to store the histogram in cGy
    # only if the structure has contour data or the dose grid exists
    if ((len(planes)) and ("PixelData" in dose.ds)):

        # Get the dose and image data information
        dd = dose.GetDoseData()
        id = dose.GetImageData()

        # Generate a 2d mesh grid to create a polygon mask in dose coordinates
        # Code taken from Stack Overflow Answer from Joe Kington:
        # https://stackoverflow.com/q/3654289/74123
        # Create vertex coordinates for each grid cell
        x, y = np.meshgrid(np.array(dd['lut'][0]), np.array(dd['lut'][1]))
        x, y = x.flatten(), y.flatten()
        dosegridpoints = np.vstack((x, y)).T

        maxdose = int(dd['dosemax'] * dd['dosegridscaling'] * 100)
        # Remove values above the limit (cGy) if specified
        if isinstance(limit, int):
            if (limit < maxdose):
                maxdose = limit
        hist = np.zeros(maxdose)
    else:
        return np.array([0])

    n = 0
    planedata = {}
    # Iterate over each plane in the structure
    for z, plane in iteritems(planes):
        # Get the dose plane for the current structure plane
        doseplane = dose.GetDoseGrid(z)
        planedata[z] = calculate_plane_histogram(
            plane, doseplane, dosegridpoints,
            maxdose, dd, id, structure, hist)
        n += 1
        if callback:
            callback(n, len(planes))
    # Volume units are given in cm^3
    volume = sum([p[1] for p in planedata.values()]) / 1000
    # Rescale the histogram to reflect the total volume
    hist = sum([p[0] for p in planedata.values()])
    hist = hist * volume / sum(hist)
    # Remove the bins above the max dose for the structure
    hist = np.trim_zeros(hist, trim='b')

    return hist


def calculate_plane_histogram(plane, doseplane, dosegridpoints,
                              maxdose, dd, id, structure, hist):
  #  """Calculate the DVH for the given plane in the structure."""
    contours = [[x[0:2] for x in c['data']] for c in plane]

    # If there is no dose for the current plane, go to the next plane
    if not len(doseplane):
        return (np.arange(0, maxdose), 0)

    # Create a zero valued bool grid
    grid = np.zeros((dd['rows'], dd['columns']), dtype=np.uint8)

    # Calculate the histogram for each contour in the plane
    # and boolean xor to remove holes
    for i, contour in enumerate(contours):
        m = get_contour_mask(dd, id, dosegridpoints, contour)
        grid = np.logical_xor(m.astype(np.uint8), grid).astype(np.bool)

    hist, vol = calculate_contour_dvh(
        grid, doseplane, maxdose, dd, id, structure)
    return (hist, vol)


def get_contour_mask(dd, id, dosegridpoints, contour):
  #  """Get the mask for the contour with respect to the dose plane."""
    doselut = dd['lut']

    c = matplotlib.path.Path(list(contour))
    grid = c.contains_points(dosegridpoints)
    grid = grid.reshape((len(doselut[1]), len(doselut[0])))
    
    return grid


def calculate_contour_dvh(mask, doseplane, maxdose, dd, id, structure):
   # """Calculate the differential DVH for the given contour and dose plane."""
    # Multiply the structure mask by the dose plane to get the dose mask
    mask = ma.array(doseplane * dd['dosegridscaling'] * 100, mask=~mask)
    # Calculate the differential dvh
    hist, edges = np.histogram(mask.compressed(),
                               bins=maxdose,
                               range=(0, maxdose))

    # Calculate the volume for the contour for the given dose plane
    vol = sum(hist) * ((id['pixelspacing'][0]) *
                       (id['pixelspacing'][1]) *
                       (structure['thickness']))
    return hist, vol



######################################################################################################
#AHORA A HACER COSAS DEL SPECT##
def get_dfh(structure, dose, spect, roi, limit=None, callback=None):
    #"""
    #Calculate a cumulative DFH from a DICOM RT Structure Set, Dose, and SPECT.
    
    #Parameters:
    #- structure : pydicom Dataset
    #    DICOM RT Structure Set used to determine the structure data.
    #- dose : pydicom Dataset
    #    DICOM RT Dose used to determine the dose grid.
    #- spect: numpy array
    #    Scaled SPECT data.
    #- roi : int
    #    The ROI number used to uniquely identify the structure in the structure
    #    set.
    #- limit : int, optional
    #    Dose limit in cGy as a maximum bin for the histogram.
    #- callback : function, optional
    #    A function that will be called at every iteration of the calculation.
    #
    #Returns:
    #- dfh : DVH object
    #    The computed DFH.
    #"""
  
    rtss = dicomparser.DicomParser(structure)
    rtdose = dicomparser.DicomParser(dose)
    structures = rtss.GetStructures()
    s = structures[roi]
    s['planes'] = rtss.GetStructureCoordinates(roi)
    s['thickness'] = rtss.CalculatePlaneThickness(s['planes'])
    hist = calculate_dfh(s, rtdose, spect)
    return dvh.DVH(counts=hist,
                   bins=(np.arange(0, 2) if (hist.size == 1) else
                         np.arange(0, hist.size + 1) / 100),
                   dvh_type='differential',
                   dose_units='gy',
                   name=s['name']
                   ).cumulative

def calculate_dfh(structure, dose_parser, spect_data):
    planes = structure['planes']
    logger.debug("Calculating DFH of %s %s", structure['id'], structure['name'])

    dd = dose_parser.GetDoseData()
    maxdose = int(dd['dosemax'] * dd['dosegridscaling'] * 100)
    hist = np.zeros(maxdose)

    n = 0
    planedata = {}
    for z, plane in iteritems(planes):
        doseplane = dose_parser.GetDoseGrid(z)

        # Determine the slice index for spect_data corresponding to the current z plane
        slice_index = int(float(z))  # Modified line here

        current_spect_slice = spect_data[:, :, slice_index]

        planedata[z] = calculate_plane_dfh(plane, doseplane, current_spect_slice, maxdose, dd, dose_parser, hist, structure)
        n += 1

    volume = sum([p[1] for p in planedata.values()]) / 1000
    hist = sum([p[0] for p in planedata.values()])
    hist = hist * volume / sum(hist)
    hist = np.trim_zeros(hist, trim='b')

    return hist

def calculate_plane_dfh(plane, doseplane, current_spect_slice, maxdose, dd, dose_parser, hist, structure):
    contours = [[x[0:2] for x in c['data']] for c in plane]

    # If there is no dose for the current plane, go to the next plane
    if not len(doseplane):
        return (np.arange(0, maxdose), 0)

    # Create a zero-valued bool grid
    rows, cols = current_spect_slice.shape
    grid = np.zeros((rows, cols), dtype=np.bool_)   # Adjusted grid size

    # Generate a 2D mesh grid to create a polygon mask in dose coordinates
    x, y = np.meshgrid(np.array(dd['lut'][0]), np.array(dd['lut'][1]))
    x, y = x.flatten(), y.flatten()
    dosegridpoints = np.vstack((x, y)).T

    # Calculate the histogram for each contour in the plane
    # and boolean xor to remove holes
    for i, contour in enumerate(contours):
        id = dose_parser.GetImageData()
        m = get_contour_mask(dd, id, dosegridpoints, contour)        
        grid = np.logical_xor(m.astype(np.uint8), grid).astype(np.bool_)

    hist, vol = calculate_contour_dfh(grid, doseplane, current_spect_slice, maxdose, dd, structure)
    return (hist, vol)


def calculate_contour_dfh(mask, doseplane, current_spect_slice, maxdose, dd, structure):
    mask_dose = ma.array(doseplane * dd['dosegridscaling'] * 100, mask=~mask)
    
    mask_spect = ma.array(current_spect_slice, mask=~mask)
    
    hist, edges = np.histogram(mask_dose.compressed(),
                               bins=maxdose,
                               range=(0, maxdose),
                               weights=mask_spect.compressed())

    vol = sum(hist) * ((id['pixelspacing'][0]) *
                       (id['pixelspacing'][1]) *
                       (structure['thickness']))
    return hist, vol
  
def load_and_scale_spect(directory):
  
    #Load multiple SPECT DICOM slices from a directory and scale the values using Min-Max normalization.
    
    #Parameters:
    #- directory: str
    #    Path to the directory containing all the SPECT DICOM slices.
    #    
    #Returns:
    #- scaled_spect: 3D numpy array
    #    The scaled SPECT volume.
    
    # List all the DICOM files in the directory
    dicom_files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
    dicom_files.sort()  # Assuming the files are named in sequential order

    # Load the slices and stack them to create a 3D volume
    slices = [dicomparser.DicomParser(os.path.join(directory, f)).GetImage() for f in dicom_files]
    spect_volume = np.stack(slices, axis=-1)
    
    # Min-Max normalization
    min_val, max_val = spect_volume.min(), spect_volume.max()
    scaled_spect = (spect_volume - min_val) / (max_val - min_val)
    
    return scaled_spect
    
def spatial_to_pixel_coordinates(dp, spatial_coords):
    #"""
    #Convert spatial coordinates to pixel indices using DICOM metadata.
    #
    #Parameters:
    #- dp: DicomParser object
    #    The parsed DICOM data from dicompylercore.
    #- spatial_coords: tuple
    #    The (x, y, z) spatial coordinates.
        
    #Returns:
    #- pixel_coords: tuple
    #    The (i, j, k) pixel indices.
    #"""
    # Retrieve DICOM metadata
    # Retrieve and convert DICOM metadata to float
    image_position = np.array([float(i) for i in dp.ds.ImagePositionPatient])  # (x0, y0, z0)
    pixel_spacing = np.array([float(i) for i in dp.ds.PixelSpacing])  # (dx, dy)
    
    # Check if SliceThickness exists, otherwise use a default value
    if 'SliceThickness' in dp.ds and dp.ds.SliceThickness:
        slice_thickness = float(dp.ds.SliceThickness)  # dz
    else:
        print("SliceThickness not found in DICOM header. Using default value of 1.0 mm.")
        slice_thickness = 1.0  # Default slice thickness in mm
    
    orientation = dp.ds.ImageOrientationPatient
    
    # Debug print statements
    print("Spatial Coords:", spatial_coords, type(spatial_coords[0]))
    print("Image Position:", image_position, type(image_position[0]))
    
    # Assuming the image orientation is the default (1,0,0,0,1,0)
    if not np.array_equal(orientation, [1, 0, 0, 0, 1, 0]):
        raise ValueError("Non-standard image orientation is not supported.")
    
    # Calculate pixel indices
    i = int(round((spatial_coords[0] - image_position[0]) / pixel_spacing[0]))
    j = int(round((spatial_coords[1] - image_position[1]) / pixel_spacing[1]))
    k = int(round((spatial_coords[2] - image_position[2]) / slice_thickness))
    
    return (i, j, k)

def map_spect_to_rt(rt_pixel_coords, spect_data, dp):
    spect_pixel_index = dp.GetPixelCoordinates(rt_pixel_coords)
    spect_intensity = spect_data[spect_pixel_index]
    return spect_intensity




def plot_dfh(dfh):
    # Extract data for plotting
    doses = list(dfh.keys())
    intensities = list(dfh.values())
    
    pl.plot(doses, intensities, label="DFH", linestyle='dotted')
    pl.xlabel('Dose (cGy)')
    pl.ylabel('SPECT Intensity Weighted Volume')
    pl.legend()
    pl.savefig('dfh_plot.png', dpi=75)
# ========================== Test DVH Calculation =========================== #


def main():


    # Read the example RT structure and RT dose files
    # The testdata was downloaded from the dicompyler website as testdata.zip

    # Obtain the structures and DVHs from the DICOM data

    rtssfile = 'rtstruct.dcm'
    rtdosefile = 'rtdose.dcm'
    RTss = dicomparser.DicomParser(rtssfile)
    rtdose = dicomparser.DicomParser(rtdosefile)
    #RTdose = dicomparser.DicomParser("testdata/rtdose.dcm") 
    RTstructures = RTss.GetStructures()

    # Generate the calculated DVHs
    calcdvhs = {}
    key_lungR = [clave for clave, valor in RTstructures.items() if valor['name'] == 'Lung_R']
    key_lungL = [clave for clave, valor in RTstructures.items() if valor['name'] == 'Lung_L']
    
    #print("Claves con name 'LUNG':", key_lungR)
    #print("Claves con name 'LUNG':", key_lungL)
    
    
    structure_lungR = RTstructures[key_lungR[0]]
    structure_lungL = RTstructures[key_lungL[0]]
    
    #print("estructuras 'LUNG':", structure_lungR)
    #print("estructuras 'LUNG':", structure_lungL)
    
    
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

#########################################################
    spect_directory = 'SPECT pulmon'
    scaled_spect = load_and_scale_spect(spect_directory)
    

    # Calculate DFH
    dfh = get_dfh(rtssfile, rtdosefile, scaled_spect, 5)

    # Plot DFH
    plot_dfh(dfh)
if __name__ == "__main__":
    main()