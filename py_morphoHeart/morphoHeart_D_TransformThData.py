# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 12:43:45 2021

@author: mdp18js
"""
#%% Importing python packages
import os
import numpy as np
from time import perf_counter
from vedo import *
#from vedo import Plotter, Cube, settings, Text2D
from vedo import embedWindow
embedWindow(False)
settings.legendSize = .3
# settings.legendPos = 1
settings.legendFont="VTK"

init = False
# Verify working dir
def setWorkingDir (root_path, init):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    init = True
    print("Current working directory: {0}".format(os.getcwd()))

    return root_path, init

root_path, init = setWorkingDir(os.getcwd(),init)

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
    folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]; dORv = filename[9:10]
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
    #   This section will load all the objects and dataframes need to unloop the heart and create distribution plots
    #   of the selected heart. Specifically it will load the cardiac jelly thickness and myocardium internal ballooning 
    #   meshes, the thickness points classification dataframe, and the dictionaries to create the heart centreline. 
    #   No plots will come up after running this code.
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


    #%% UNLOOP HEART
    #   Now we can unloop the heart to create 2D heatmaps of any of the thickness or ballooning measurements already 
    #   taken from this heart. To do this, initially the centreline of the heart is going to be divided into an atrial
    #   and a ventricular centrelineS, each containing about 600 points. Each of these centrelines will be then
    #   used to create 150 cross sectional planes of the chamber thickness mesh, all perpendicular to the centreline 
    #   and equally spaced. Thickness information about the points of the mesh that are cut by each of the cross-
    #   sectional planes will be saved as well as information regarding its position (e.g. angular position with 
    #   respect to the ventral center of the heart, as well as cross-sectional plane in which the points were found). 
    #   If plotshow is True in function unloopChambers, as the heart is cross-sectioned, the user will be able to see 
    #   a selection of cross-sectional planes and points that are being cut by that plane. The user needs to close those 
    #   3D plots for the code to continue processing. 
    #   Once all the information has been acquired for both chambers, a 3D interactive plot will appear showing a 
    #   selection of the planes used to cross-section each chamber. When closed, the code will transform all the 
    #   acquired data and will create 2D heatmaps for each chamber and thickness measurement. When all the heatmap 
    #   plots have been created, a new 3D interactive plot will appear where you will be able to compare the thickness 
    #   measurements on the 2D heatmap and on the color-coded mesh. 
    #   Note: Make sure you select the 'Plots' Pane in Spyder to be able to see the resulting heatmaps and compare 
    #   them to the color coded meshes. 
    #   The heatmap creation can take longer than 5min/each, and there are four to be created. 
    #   Again, be patient, it is worth it! :)
    #   ================================================================================================================
    
    saveHM = True; savePlot = True; plotshow = True
    df_AtrVent = np.asarray(df_cjThNmyocIntBall['AtrVent'])
    
    fcMeshes.plotPtClassif(filename = filename, mesh = m_cjTh, pts_whole = m_cjTh.points(), 
                           pts_class = [df_cjThNmyocIntBall['AtrVent'], df_cjThNmyocIntBall['DorsVent'], df_cjThNmyocIntBall['LeftRight']])
    
    # Unlooping the heart
    _, kspl_vSurf, kspls_HR, arr_all, arr_valve, dict_pts, dict_shapes, dict_kspl = fcMeshes.unloopChambers(filename = filename, 
                                            mesh = m_cjTh.alpha(0.05), kspl_CL = kspl_CL[0], kspl_ext = kspl_ext,
                                            no_planes = 150, pl_CLRibbon =  dict_planes['pl_Parallel2LinLine'], 
                                            param = [cj_thickness, myoc_intBall], param_name = 'CjTh_myocIntBall',
                                            df_AtrVent = df_AtrVent, dict_kspl = dict_kspl, 
                                            dict_shapes = dict_shapes, dict_pts = dict_pts, 
                                            dict_planes = dict_planes, 
                                            dir_results = dir_results, plotshow = plotshow, tol=0.05)
    
    print('Creating heatmaps for each chamber... Note: this can take a while.')
    heatmaps_cj = fcMeshes.heatmapUnlooped(filename = filename, val2unloop = 'cj_thickness', dataname = 'CjTh',
                                           dir_results = dir_results, dirImgs = directories[4], 
                                           names= ['unloopAtrCjTh_myocIntBall', 'unloopVentCjTh_myocIntBall'],
                                           saveHM = saveHM, savePlot = savePlot)
    del heatmaps_cj
    
    heatmaps_ball = fcMeshes.heatmapUnlooped(filename = filename, val2unloop = 'myoc_intBall', dataname = 'myocIntBall',
                                           dir_results = dir_results, dirImgs = directories[4], 
                                           names= ['unloopAtrCjTh_myocIntBall', 'unloopVentCjTh_myocIntBall'],
                                           saveHM = saveHM, savePlot = savePlot)
    
    if 'kspl_vSurf' not in locals():
        kspls_HR = []; kspl_vSurf = []
    m_cjTh.pointColors(cj_thickness, cmap="jet", vmin=0, vmax=25).addScalarBar()
    m_cjTh.mapper().SetScalarRange(0,25)
    m_myocIntBall.pointColors(myoc_intBall, cmap="jet", vmin=0, vmax=100).addScalarBar()
    m_myocIntBall.mapper().SetScalarRange(0,100)
    settings.legendSize = .2
    txt = Text2D(filename+"\n\n >> Plots to validate heatmaps", c="k", font= font)
    vp = Plotter(N=4, axes = 4)
    vp.show(m_cjTh.alpha(1), cl_ribbonV, disk, txt, at = 0)
    vp.show(m_myocIntBall.alpha(1), cl_ribbonV, disk,  at = 1)
    vp.show(m_cjTh.clone().alpha(0.01), kspl_vSurf, kspls_HR, kspl_ext, arr_all, disk, at = 2)
    vp.show(m_myocIntBall.clone().alpha(0.01), kspl_vSurf, kspls_HR, kspl_ext, arr_valve, disk, at = 3, interactive = True)
    
    del heatmaps_ball
    
    #%% GET THICKNESS REGION PLOTS (PROBALITY DENSITY ESTIMATION PLOTS)
    #   Finally, we need to create distribution plots that will allow us to compare the thickness distribution between 
    #   regions of the heart.
    #   ================================================================================================================
    
    file_num = df_file[df_file['Folder']==folder].index.values[0]
    df_cjPDFs = fcAn.kdeThPlots(filename = filename, df_file = df_file, file_num = file_num, variable = 'cj_thickness', 
                                thData = df_cjThNmyocIntBall, dir2save = directories[4], save = True)
    
    if save: 
        # Append all dicts to one object dict
        dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
                                             names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
        
        fcBasics.saveDF(filename = filename, df2save = df_cjPDFs, df_name = 'cjPDFs', dir2save = dir_results)
        fcBasics.saveDF(filename = filename, df2save = df_cjPDFs, df_name = 'cjPDFs', dir2save = os.path.join(dir_data2Analyse, 'R_All','df_cjPDFs'))
    
    toc = perf_counter()
    fcBasics.printTime(tic, toc, 'Transform thickness data')

#%% Init
init = True
