# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 13:36:21 2021

@author: mdp18js
"""

#%% Importing python packages
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import count
from vedo import *
from vedo import embedWindow#, settings, Plotter, Text2D
embedWindow(False)
settings.legendSize = .3
# settings.legendPos = 1
settings.legendFont="VTK"

init = True
# Verify working dir
def setWorkingDir (root_path, init):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    # init = True
    print("Current working directory: {0}".format(os.getcwd()))

    return root_path, init

root_path, init = setWorkingDir(os.getcwd(),init)

c="k"; font= 'VTK' # 'CallingCode'
azimuth = 0

#%% Start E_Plot
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    # from morphoHeart_modules import morphoHeart_funcContours as fcCont
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes

    #%% Get directories and file
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = 'R')

    #%% Plot things for just one heart
    # Get file to process and directories
    folder, df_file, file_num, blind = fcBasics.selectFile(df_dataset); filename = folder[2:]; dORv = filename[9:10]
    #stage = df_file.loc[file_num,'Stage']
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = 'R')

    # Import df_results
    df_res = fcBasics.loadDF(filename = filename, file = 'ResultsDF', dir_results = dir_results)
    file_num = df_res[df_res['Folder']== filename+'_2A'].index.values[0]
    if dORv == 'D':
        azimuth = -90
    elevation = df_res.loc[file_num,'ang_Heart']
    xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
    res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]
    
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
        
    kspl_CL, linLines, sph_CL, sph_CL_colour, dict_shapes, dict_kspl = fcMeshes.createCLs(dict_cl = dicts[1:], dict_pts = dict_pts,
                                                                                          dict_kspl = dict_kspl, dict_shapes = dict_shapes,
                                                                                          colors = ['deepskyblue', 'tomato'])
    #%%
    from pathlib import Path
    vtk_in_dir = Path(directories[2])
    
    try: 
        cyl_chamber = dict_shapes['cyl2CutChambers_final']
        disk = Cylinder(pos = cyl_chamber['cyl_centre'],r = cyl_chamber['radius_max'], height = 2*0.225, 
                        axis = cyl_chamber['cyl_axis'], c = 'purple', cap = True, res = 300)
        
        num_pt = dict_pts['numPt_CLChamberCut'] 
        
        atr_meshes = []; vent_meshes = []
        atr_meshes, vent_meshes, dict_shapes, s3_cyl  = fcMeshes.getChamberMeshes(filename = filename,
                                    end_name = ['ch0_cut'], names2cut = ['Myoc'],
                                    kspl_CL = kspl_CL[0], num_pt = num_pt, atr_meshes = atr_meshes, vent_meshes = vent_meshes,
                                    dir_txtNnpy = directories[1], dict_shapes = dict_shapes, dict_pts = dict_pts, 
                                    resolution = res, plotshow = True)
    
        atr_meshes, vent_meshes, dict_shapes, _ = fcMeshes.getChamberMeshes(filename = filename,
                                        end_name = ['ch1_cut', 'cj', 'ch0_cut_ext', 'ch1_cut_int', 'ch0_cut_int'],
                                        names2cut = ['Endo', 'CJ', 'Ext.Myoc', 'Int.Endo', 'Ext.CJ'],
                                        kspl_CL = kspl_CL[0], num_pt = num_pt, atr_meshes = [m_atr[0]], vent_meshes = [m_vent[0]],
                                        dir_txtNnpy = directories[1], dict_shapes = dict_shapes, dict_pts = dict_pts, 
                                        resolution = res, s3_cyl = s3_cyl, plotshow = plotshow)
        
        
        
        
        
        
    except: 
        print(filename, '- No cylinder to use to cut')
                    
    n = 0
    for name in vtk_in_dir.glob('*.vtk*'):
        head, tail = os.path.split(name)
        n +=1
        if 'thickness' not in tail: 
            if 'Ball' not in tail: 
                print('Processing: ', tail, '-', str(n))
                [m_mesh] = fcMeshes.openMeshes(filename = filename, meshes_names = [tail[19:-4]],
                                                                          extension = 'vtk', dir_stl = directories[2],
                                                                          alpha = [1], dict_colour = dict_colour)

    
                vp = Plotter(N=1, axes = 13)
                vp.show(m_mesh, at=0, interactive = True)
        

#%% 
init = True