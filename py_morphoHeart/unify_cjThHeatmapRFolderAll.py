# -*- coding: utf-8 -*-
"""
morphoHeart - F. PLOT

@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
from pathlib import Path
from itertools import count
# from vedo import *
from vedo import embedWindow, settings, Plotter, Text2D
embedWindow(False)
settings.legendSize = .3
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

c="k"; font= 'VTK' # 'CallingCode'
azimuth = 0

#%% Start E_Plot
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes

    #%% Get directories and file
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = '2A')

    #%% Recreate plots for selected hearts
    # Get file to process and directories
    df_files = fcBasics.selectHearts(df_dataset)
    for nn, folder in zip(count(), df_files['Folder']):
        df_file = df_files.iloc[[nn]]
        filename = folder[2:]
        dORv = filename[9:10]
        print(folder)
            
        #folder, df_file, blind = fcBasics.selectFile(df_dataset, end_name = 'R'); filename = folder[2:]; dORv = filename[9:10]
        #stage = df_file.loc[file_num,'Stage']
        # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse, 6. dir_LS_Folder selected
        dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = 'R')
        # df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = 'R', out_type = 'xlsx')
        
        # Import the metadata to know pixel size and distance between slices
        xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
        res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]
        # Import df_results
        df_res = fcBasics.loadDF(filename = filename, file = 'ResultsDF', dir_results = dir_results)
        file_num = df_res[df_res['Folder']==filename+'_2A'].index.values[0]
        if dORv == 'D' or 'CJ' in filename:
            azimuth = -90
        else: 
            azimuth = 0
        # elevation = df_res.loc[file_num,'ang_Heart']
        
         # Get main directories
        # _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
        # df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse)
        # Get file to process and directories
        # folder, df_file, file_num, blind = fcBasics.selectFile(df_dataset); filename = folder[0:-3]; dORv = filename[9:10]
        #stage = df_file.loc[file_num,'Stage']
        # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse, 6. dir_LS_Folder selected
        # dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')
        # Import the metadata to know pixel size and distance between slices
        # xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
        # res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]
        # Import df_results
        # df_res = fcBasics.loadDF(filename = filename, file = 'ResultsDF', dir_results = dir_results)
        # file_num = df_res[df_res['Folder']==folder].index.values[0]
        # if dORv == 'D' or 'CJ' in filename:
        #     azimuth = -90
        # else: 
        #     azimuth = 0
            
         # Get existing cl_dictionaries
        [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes, dicts_cl] = fcBasics.import_dicts('mH_C', filename, directories)
    
        # Initialise variables
        dict_shapes = dict()
        txt = Text2D(filename, c="k", font= 'CallingCode')
        rotAngle =  df_res.loc[file_num,'ang_HeartS']
        saveHM = False; savePlot = True
        # Select analysis to run
        tissue_analysis, tissue_opt, m_myoc, dict_unloop = fcMeshes.load_tissues2unloop(filename, directories, dir_results, dict_colour,
                                                                                        default = True)
        for tissue in tissue_analysis:
            print('- Creating heatmaps for each chamber... Note: this can take a while.')
            heatmaps_th, scale_th = fcMeshes.heatmapUnlooped(filename = filename, 
                                                              val2unloop = dict_unloop[tissue]['param_name'], 
                                                              dir_results = dir_results, dirImgs = directories[4], 
                                                              save_names= dict_unloop[tissue]['save_names'],
                                                              hm_names = dict_unloop[tissue]['hm_names'],
                                                              saveHM = saveHM, savePlot = savePlot, 
                                                              default = True, cmap = 'turbo')

#%% Init
init = True
    