# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 13:29:33 2021

@author: mdp18js
"""
#%% Importing python packages
import os
import numpy as np
from time import perf_counter
from vedo import Plotter, Cylinder, settings, Text2D
from vedo import embedWindow
embedWindow(False)
settings.legendSize = .3
settings.legendFont="VTK"

# Verify working dir
def setWorkingDir (root_path, init = False):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    init = not bool(int(input('Do you want to execute the script all at once or run it by cells? \n\t[0]: all at once (recommended if you are already familiar with the script)\n\t[1]: by cells (recommended if you NOT yet familiar with the script). >>>: ')))
    print("Current working directory: {0}".format(os.getcwd()))

    return root_path, init

root_path, init = setWorkingDir(os.getcwd())

c="k"; font= 'VTK';
save = True; plot = True; plotshow = False
azimuth = 0

#%% Start D_TransformThData
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
    from morphoHeart_modules import morphoHeart_funcAnalysis as fcAn
    tic = perf_counter()

    #%% SELECT FILE AND GET METADATA
    #   This section allows the user to select file to process, get its main directories and metadata,
    #   and define some properties
    #   ================================================================================================================

    # Get main directories
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse)
    # Get file to process and directories
    folder, df_file, file_num, blind = fcBasics.selectFile(df_dataset); filename = folder[0:-3]; dORv = filename[9:10]
    #stage = df_file.loc[file_num,'Stage']
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse, 6. dir_LS_Folder selected
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')

    # Import the metadata to know pixel size and distance between slices
    xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
    res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]
    # Import df_results
    df_res = fcBasics.loadDF(filename = filename, file = 'ResultsDF', dir_results = dir_results)
    file_num = df_res[df_res['Folder']==folder].index.values[0]
    if dORv == 'D':
        azimuth = -90
    # Initialise variables
    dict_shapes = dict()
    txt = Text2D(filename, c="k", font= 'CallingCode')
    
    
    #%% LOAD MESHES, CENTRELINES AND OBJECTS
    #   This section will load all the objects and dataframes needed to unloop the heart and create distribution plots
    #   of the selected heart. Specifically it will load the cardiac jelly thickness and myocardium internal ballooning 
    #   meshes, the thickness points classification dataframe, and the dictionaries to create the heart centreline. 
    #   No interactive plots will pop-up after running this code.
    #   ================================================================================================================

    # Get existing cl_dictionaries
    _, mesh_name = fcBasics.code4vmtkCL(filename = filename, mesh_name = ['myoc_int','endo_ext'],
                                                dir_cl = directories[3], printshow = False)
    # Import dictionaries
    dicts = fcBasics.loadDicts(filename = filename, dicts_name = ['dict_obj']+[txt+'_npcl' for txt in mesh_name],
                                                                    directories = [directories[0]]+ [directories[3]]*len([txt+'_npcl' for txt in mesh_name]))
    dict_obj = fcMeshes.splitDicts(dicts[0])
    if len(dict_obj) == 4:
        [dict_planes, dict_pts, dict_kspl, dict_colour] = dict_obj
    else:
        [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes] = dict_obj
    
    
    names4video, meshes4video, rangeThBall = fcPlot.selectMeshes4Video(names = ['myoc','endo','cj', 'cj_thickness',
                                                'myoc_thickness','endo_thickness','myoc_intBall', 'myoc_extBall'],
                                                                       meshes = [], ranges = [])
    
    
    # Import meshes
    [[m_cjTh, m_myocIntBall], [cj_thickness, myoc_intBall]] = fcMeshes.openThicknessMeshes(filename = filename,
                                                                   meshes_names = ['cj_thickness','myoc_intBall'], 
                                                                   extension = 'vtk', dir_stl = directories[2], 
                                                                   dir_txtNnpy = directories[1])
    # Import thickness points classification dataframe
    df_cjThNmyocIntBall = fcBasics.loadDF(filename = filename, file = 'df_cjThNmyocIntBall', dir_results = dir_results)
    
    # Create centreline
    kspl_CL, linLines, _, _, _, _ = fcMeshes.createCLs(dict_cl = dicts[1:], dict_pts = dict_pts,
                                                       dict_kspl = dict_kspl, dict_shapes = dict_shapes,
                                                       colors = ['deepskyblue', 'tomato'])
    # Get centreline ribbon
    cl_ribbonV, kspl_ext, _, _, _ = fcMeshes.createCLRibbon(filename = filename, kspl_CL2use = kspl_CL[0], linLine = linLines[0],
                                                   mesh = m_cjTh, dict_kspl = dict_kspl, dict_shapes = dict_shapes, 
                                                   dict_planes = dict_planes, clRib_type = 'extV', plotshow = False)
    
    # Get ch_cylinder
    cyl_chamber = dict_shapes['cyl2CutChambers_final']
    disk = Cylinder(pos = cyl_chamber['cyl_centre'],r = cyl_chamber['radius_max'], height = 2*0.225, 
                    axis = cyl_chamber['cyl_axis'], c = 'purple', cap = True, res = 300)


init = True
