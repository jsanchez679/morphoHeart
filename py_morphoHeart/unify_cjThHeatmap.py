# -*- coding: utf-8 -*-
"""
morphoHeart - F. PLOT

@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
from vedo import *
from vedo import embedWindow#, settings, Plotter, Text2D
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
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = '2A')

    #%% Plot things for just one heart
    # Get file to process and directories
    folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]; dORv = filename[9:10]
    #stage = df_file.loc[file_num,'Stage']
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')

    # Import df_results
    df_res = fcBasics.loadDF(filename = filename, file = 'ResultsDF', dir_results = dir_results)
    file_num = df_res[df_res['Folder']== filename+'_2A'].index.values[0]
    if dORv == 'D':
        azimuth = -90
    elevation = df_res.loc[file_num,'ang_Heart']
    
    #%% Get plot of cardiac jelly thickness heatmap with thickness range given as user input 
    [m_thAll, colour_thAll] = fcMeshes.openThicknessMeshes(filename = filename, meshes_names = ['cj_thickness'], extension = 'vtk',
                                  dir_stl = directories[2], dir_txtNnpy = directories[1])
    m_cjTh = m_thAll[0]
    q_min = fcBasics.ask4input('Enter minimum value for -cardiac jelly thickness- heatmap range: ',float)
    q_max = fcBasics.ask4input('Enter maximum value for -cardiac jelly thickness- heatmap range: ',float)
    
    [cj_thickness] = fcBasics.loadNPY(filename = filename, names = ['cj_thickness'], dir_txtNnpy = directories[1])
    m_cjTh.pointColors(cj_thickness, cmap="jet", vmin=q_min, vmax=q_max)
    m_cjTh.addScalarBar()
    m_cjTh.alpha(1)
    m_cjTh.mapper().SetScalarRange(q_min,q_max)

    vp = Plotter(N=1, axes=13)
    vp.show(m_cjTh, at=0, azimuth = azimuth, elevation = elevation, interactive=True)
    
    saveVideo = fcBasics.ask4input('Do you want to save rotating video with the set range? [0]:no/[1]:yes: ',bool)
    if saveVideo:
        fcMeshes.saveVideo (filename = filename, info = 'cj_thicknessUserRange', meshes4video = [m_cjTh],
                            rotAngle  = df_res.loc[file_num,'ang_Heart'], dir2save =  directories[4], alpha_cube = 0, plotshow=True)

#%% Init
init = True
    