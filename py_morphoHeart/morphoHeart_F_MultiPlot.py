# -*- coding: utf-8 -*-
"""
morphoHeart - F. MULTIPLOT HEARTS

@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
import numpy as np
from itertools import count
from vedo import Plotter, Cube, settings
from vedo import embedWindow
embedWindow(False)
settings.legendSize = .3
# settings.legendPos = 1
settings.legendFont="VTK"
#settings.useParallelProjection=True

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
alpha_cube = 0.1
s_cube = 350

#%% Start F_MultiPlot
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcPlot as fcPlot

    #%% Get directories and file
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = 'R')

    #%% Select hearts to plot
    df_files = fcBasics.selectHearts(df_dataset)

    #%% Plot A - One plot per heart
    objs = fcPlot.selectObjects()
    dict_paths = fcPlot.createDictPaths(objs, df_files, dir_data2Analyse)
    meshes = fcPlot.loadMultMeshes(objs, dict_paths)
    
    settings.legendSize = .3
    vp = Plotter(shape = (1, len(dict_paths)), axes=13, sharecam = True); m = 0
    for i, n, file in zip(count(), range(len(dict_paths)), df_files['Folder'].tolist()):
        txt = fcPlot.plotTitle(file, df_files)
        scale_cube = Cube(pos=meshes[m].centerOfMass(), side=s_cube, c='gray', alpha=alpha_cube)
        if n != len(dict_paths)-1:
            vp.show(meshes[m:m+len(objs)], scale_cube, txt, at=n)
            m += len(objs)
        else:
            vp.show(meshes[m:m+len(objs)],scale_cube, txt, at=n, zoom = 1.5, interactive=True)
            # screenshot(filename='screenshot.png', scale=None, returnNumpy=False)

    #%% Plot B - One heart per row, objects plotted in different columns
    objs = fcPlot.selectObjects()#['myoc','endo', 'cj']
    dict_paths = fcPlot.createDictPaths(objs, df_files, dir_data2Analyse)
    meshes = fcPlot.loadMultMeshes(objs, dict_paths)

    vp = Plotter(shape = (len(dict_paths), len(objs)), axes=13, sharecam = True); m = 0
    numbs = np.arange(0,len(dict_paths)*len(objs),len(objs)).tolist()
    for i, n in zip(count(), range(len(dict_paths)*len(objs))):
        file = df_files['Folder'].tolist()[n//len(objs)]
        #print(' n: ', n, ' - file: ', file)
        txt = fcPlot.plotTitle(file, df_files)
        if m in numbs:
            scale_cube = Cube(pos=meshes[m].centerOfMass(), side=s_cube, c='white', alpha=alpha_cube)
            # sph = Sphere(pos=meshes[m].centerOfMass(), r=2, c='red', alpha=1)
        if n != len(dict_paths)*len(objs)-1:
            if m in numbs:
                vp.show(meshes[m], scale_cube, txt, at=n)
            else:
                vp.show(meshes[m], scale_cube, at=n)
            m += 1
        else:
            vp.show(meshes[m], scale_cube, at=n, zoom = 1.5, interactive=True)

#%% 
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
    r_folders = df_files['Folder']
    filename = [r_folders[i][2:] for i in range(len(r_folders))]
    
    for f, folder, name in zip(count(), r_folders, filename):
        print(f, name)
        dir_results = os.path.join(dir_data2Analyse, folder, 'Results_'+name)
        heatmapsf = fcMeshes.filterUnloopedDF(filename = name, thickness = 'cj_thickness', 
                                              dir_results = dir_results, dir_data2Analyse = dir_data2Analyse,
                                              names= ['unloopAtrCjTh_myocIntBall', 'unloopVentCjTh_myocIntBall'], 
                                              save_names= ['unloopAtrCjTh', 'unloopVentCjTh'], 
                                              saveHM = True)
        
    for f, folder, name in zip(count(), r_folders, filename):
        print(f, name)
        dir_results = os.path.join(dir_data2Analyse, folder, 'Results_'+name)
        heatmapsf = fcMeshes.filterUnloopedDF(filename = name, thickness = 'myoc_intBall', 
                                              dir_results = dir_results, dir_data2Analyse = dir_data2Analyse,
                                              names= ['unloopAtrCjTh_myocIntBall', 'unloopVentCjTh_myocIntBall'], 
                                              save_names= ['unloopAtrmyocIntBall', 'unloopVentmyocIntBall'], 
                                              saveHM = True)
        

#%% Init
init = True
# from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
# filename = file[2:]
# fcMeshes.saveMesh(filename, meshes[0], 'endoExt_stl', dir_data2Analyse, 'stl')
# fcMeshes.saveMesh(filename, meshes[1], 'endo_stl', dir_data2Analyse, 'stl')