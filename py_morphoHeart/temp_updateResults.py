# -*- coding: utf-8 -*-
"""
Created on Sun May  2 00:26:00 2021

@author: mdp18js
"""
#%% Importing python packages
import os
import numpy as np
from time import perf_counter
from itertools import count
from vedo import Plotter, Cube, settings, Text2D#, vedo2trimesh, trimesh2vedo, Cylinder
from vedo import embedWindow
embedWindow(False)
settings.legendFont="VTK"

init = False
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

c="k"; font= 'VTK'
azimuth = 0
alpha_cube = 0
s_cube = 350

#%% Start F_MultiPlot
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    #from morphoHeart_modules import morphoHeart_funcPlot as fcPlot
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
    from morphoHeart_modules import morphoHeart_funcAnalysis as fcAn
    
     #%% Get directories and file
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = 'R', out_type = 'xlsx')
    
    plot = True
    save = True
    #%% Select hearts to plot
    df_files = fcBasics.selectHearts(df_dataset)
    
    #%%
    repeat = []
    repeat_atrVent = []
    for index, row in df_files.iterrows():
        tic = perf_counter()
        folder = row['Folder']
        filename = folder[2:]; dORv = filename[9:10]
        print(folder, filename)
        df_file = df_files.loc[df_files['Folder'] == folder]
        dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = 'R')
        xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
        res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]
        # Import df_results
        df_res = fcBasics.loadDF(filename = filename, file = 'ResultsDF', dir_results = dir_results)
        file_num = index
        if dORv == 'D' or 'CJ' in filename:
            azimuth = -90
        else: 
            azimuth = 0
        
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
        # Import meshes
        [m_myoc, m_endo, m_cj, m_cjOut, m_cjIn] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc','endo','cj','cj_out','cj_in'],
                                                                      extension = 'vtk', dir_stl = directories[2],
                                                                      alpha = [0.05,0.05,0.05,1,1], dict_colour = dict_colour)
        [m_myocInt, m_myocExt, m_endoInt, m_endoExt] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc_int','myoc_ext','endo_int', 'endo_ext'],
                                                                      extension = 'vtk', dir_stl = directories[2],
                                                                      alpha = [1,1,1,1], dict_colour = dict_colour)
    
        scale_cube = Cube(pos=m_myoc.centerOfMass(), side=350, c='white', alpha=0.01)
        # Plot meshes
        if plot:
            text = str(filename); txt = Text2D(text, c=c, font=font)
            settings.legendSize = .20
            vp = Plotter(N=6, axes=13)
            vp.show(m_myoc, txt, at=0)
            vp.show(m_endo, at=1)
            vp.show(m_cj, at=2)
            vp.show(m_cjIn, at=3)
            vp.show(m_cjOut, at=4)
            vp.show(m_myoc.clone().alpha(0.01), m_endo.clone().alpha(0.01), m_cj.clone().alpha(0.01), scale_cube, at=5, zoom = 2, azimuth = azimuth, interactive=True)
        
        #-----------------
        #Import atrMyoc y VentMyoc
        try:
            [m_atrMyoc, m_ventMyoc] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc_atr','myoc_vent'],
                                                          extension = 'vtk', dir_stl = directories[2],
                                                          alpha = [0.05,0.05], dict_colour = dict_colour)
            settings.legendSize = .20
            vp = Plotter(N=1, axes = 10)
            vp.show(m_atrMyoc, m_ventMyoc, at = 0, interactive=True)
            
            #Get minimum volume cylinder which contains each of the chambers 
            df_res = fcMeshes.addMinCylVol2df(df_res, file_num, [m_atrMyoc, m_ventMyoc], ['myoc_atr','myoc_vent'])
            fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
            fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = os.path.join(dir_data2Analyse, 'R_All', 'df_meas'))
        
        except:
            repeat_atrVent.append(filename)
            
        #-----------------
        #Rerun cut LnR/Rerun classify points and save 
        # Import meshes
        [[m_cjTh, m_myocIntBall], [cj_thickness, myoc_intBall]] = fcMeshes.openThicknessMeshes(filename = filename,
                                                                       meshes_names = ['cj_thickness','myoc_intBall'], 
                                                                       extension = 'vtk', dir_stl = directories[2], 
                                                                       dir_txtNnpy = directories[1])
        
        # [[m_myocTh, m_endoTh], [myoc_thickness, endo_thickness]] = fcMeshes.openThicknessMeshes(filename = filename,
        #                                                                meshes_names = ['myoc_thickness', 'endo_thickness'], 
        #                                                                extension = 'vtk', dir_stl = directories[2], 
        #                                                                dir_txtNnpy = directories[1])
        # [[m_cjTh], [cj_thickness]] = fcMeshes.openThicknessMeshes(filename = filename,
        #                                                                meshes_names = ['cj_thickness'], 
        #                                                                extension = 'vtk', dir_stl = directories[2], 
        #                                                                dir_txtNnpy = directories[1])
        
        # Import thickness points classification dataframe
        df_cjThNmyocIntBall = fcBasics.loadDF(filename = filename, file = 'df_cjThNmyocIntBall', dir_results = dir_results)
        # Get centreline ribbon
        cl_ribbon, kspl_ext, dict_kspl, dict_shapes, dict_planes = fcMeshes.createCLRibbon(filename = filename, kspl_CL2use = kspl_CL[0], linLine = linLines[0],
                                                                 mesh = m_cjTh, dict_kspl = dict_kspl, dict_shapes = dict_shapes, dict_planes = dict_planes)
        dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
                                                 names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
            
        [m_cjThLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_cjTh], cl_ribbon = cl_ribbon)
        try: 
            [m_atrExtCJ, m_cj] = fcMeshes.openMeshes(filename = filename, meshes_names = ['cjExt_atr', 'cj'],
                                                          extension = 'vtk', dir_stl = directories[2],
                                                          alpha = [1, 0.05], dict_colour = dict_colour)
            continue_class = True
        except: 
            try: 
                [m_atrExtCJ, m_cj] = fcMeshes.openMeshes(filename = filename, meshes_names = ['cjExt_atrmyoc_vent', 'cj'],
                                                              extension = 'vtk', dir_stl = directories[2],
                                                              alpha = [1,0.05], dict_colour = dict_colour)
                dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [m_atrExtCJ.color('orange')], names = ['cjExt_atr'],
                                              dict_colour = dict_colour, dir_stl = directories[2], extension = 'vtk')
                continue_class = True
            except: 
                print('- No atrCJ mesh found! - ', filename)
                repeat.append(filename)
                continue_class = False
        
        if continue_class:
            settings.legendSize = .20
            vp = Plotter(N=2, axes = 10)
            vp.show(m_atrExtCJ, m_cjTh, at = 0)
            vp.show(m_atrExtCJ, m_cj, at=1, interactive=True)
            
            df_cjThNmyocIntBall = fcMeshes.classifyHeartPts(filename = filename, dict_planes = dict_planes,
                                                m_whole = m_cjTh, m_left = m_cjThLnR[0], m_atr = m_atrExtCJ, 
                                                data = [cj_thickness, myoc_intBall], 
                                                names_data = ['cj_thickness', 'myoc_intBall'], plot_show = True)
            fcBasics.saveDF(filename = filename, df2save = df_cjThNmyocIntBall, df_name = 'df_cjThNmyocIntBall',
                                dir2save = dir_results)
            dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
                                                 names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
            
            #-----------------
            #% recalculate unloop chambers?
            df_AtrVent = np.asarray(df_cjThNmyocIntBall['AtrVent'])
            # Unlooping the heart chambers
            _, kspl_vSurf, kspls_HR, arr_all, arr_valve, dict_pts, dict_shapes, dict_kspl = fcMeshes.unloopChambers(filename = filename, 
                                                mesh = m_cjTh.alpha(0.05), kspl_CL = kspl_CL[0], kspl_ext = kspl_ext,
                                                no_planes = 150, pl_CLRibbon =  dict_planes['pl_Parallel2LinLine'], 
                                                param = [cj_thickness, myoc_intBall], param_name = 'CjTh_myocIntBall',
                                                df_AtrVent = df_AtrVent, dict_kspl = dict_kspl, 
                                                dict_shapes = dict_shapes, dict_pts = dict_pts, 
                                                dict_planes = dict_planes, 
                                                dir_results = dir_results, plotshow = False, tol=0.05)
            
            print('- Creating heatmaps for each chamber... Note: this can take a while.')
            heatmaps_cj, scale_cj = fcMeshes.heatmapUnlooped(filename = filename, val2unloop = 'cj_thickness', dataname = 'CjTh',
                                               dir_results = dir_results, dirImgs = directories[4], 
                                               names= ['unloopAtrCjTh_myocIntBall', 'unloopVentCjTh_myocIntBall'],
                                               saveHM = True, savePlot = True, default = True)
            del heatmaps_cj
        
            heatmaps_ball, scale_ball = fcMeshes.heatmapUnlooped(filename = filename, val2unloop = 'myoc_intBall', dataname = 'myocIntBall',
                                                   dir_results = dir_results, dirImgs = directories[4], 
                                                   names= ['unloopAtrCjTh_myocIntBall', 'unloopVentCjTh_myocIntBall'],
                                                   saveHM = True, savePlot = True, default = True)
            del heatmaps_ball
            
            dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
                                                 names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
            
            #-----------------
            #% recalculate dist plots?
            df_cjPDFs = fcAn.kdeThPlots(filename = filename, df_file = df_file, file_num = file_num, variable = 'cj_thickness', 
                                    thData = df_cjThNmyocIntBall, dir2save = directories[4], save = True)
            fcBasics.saveDF(filename = filename, df2save = df_cjPDFs, df_name = 'cjPDFs', dir2save = dir_results)
            fcBasics.saveDF(filename = filename, df2save = df_cjPDFs, df_name = 'cjPDFs', dir2save = os.path.join(dir_data2Analyse, 'R_All','df_cjPDFs'))
        
            fcBasics.alert('bubble', 1)
            toc = perf_counter()
            fcBasics.printTime(tic, toc, 'Reprocess - '+filename)
        
        
        #%%
        # # from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
        # df_files = fcBasics.selectHearts(df_dataset)
        # r_folders = df_files['Folder']
        
        # for f, folder in zip(count(), r_folders):
        #     name = folder[2:]
        #     print(f, name)
        #     dir_results = os.path.join(dir_data2Analyse, folder, 'Results_'+name)
        #     heatmapsf = fcMeshes.filterUnloopedDF(filename = name, thickness = 'cj_thickness', 
        #                                           dir_results = dir_results, dir_data2Analyse = dir_data2Analyse,
        #                                           names= ['unloopAtrCjTh_myocIntBall', 'unloopVentCjTh_myocIntBall'], 
        #                                           save_names= ['unloopAtrCjTh', 'unloopVentCjTh'], 
        #                                           saveHM = True)
            
        # for f, folder in zip(count(), r_folders):
        #     name = folder[2:]
        #     print(f, name)
        #     dir_results = os.path.join(dir_data2Analyse, folder, 'Results_'+name)
        #     heatmapsf = fcMeshes.filterUnloopedDF(filename = name, thickness = 'myoc_intBall', 
        #                                           dir_results = dir_results, dir_data2Analyse = dir_data2Analyse,
        #                                           names= ['unloopAtrCjTh_myocIntBall', 'unloopVentCjTh_myocIntBall'], 
        #                                           save_names= ['unloopAtrmyocIntBall', 'unloopVentmyocIntBall'], 
        #                                           saveHM = True)

    
#%% Init
init = True