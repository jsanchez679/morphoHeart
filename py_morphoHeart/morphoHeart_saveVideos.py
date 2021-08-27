# -*- coding: utf-8 -*-
"""
morphoHeart - saveVideos

@author: Juliana Sanchez-Posada
Version: 20th July, 2021
"""

#%% Importing python packages
import os
from pathlib import Path
from vedo import embedWindow, settings, Plotter, Text2D
embedWindow(False)

# Verify working dir
def setWorkingDir (root_path, init = False):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    init = not bool(int(input('Do you want to execute the script all at once or run it by cells? \n\t[0]: all at once (recommended if you are already familiar with the script)\n\t[1]: by cells (recommended if you are NOT yet familiar with the script). >>>: ')))
    print("Current working directory: {0}".format(os.getcwd()))
    if not init: 
        print('\nIMPORTANT NOTES:\n- Remember to start running from cell %Start saveVideos.\n- NEVER run as an individual cell the cell called %Importing python packages')
    return root_path, init

root_path, init = setWorkingDir(os.getcwd())

c="k"; font= 'VTK' 
azimuth = 0

#%% Start saveVideos
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes

    #%% Get directories and file
    end_name_opt = ['2A', 'R']
    end_name = end_name_opt[fcBasics.ask4input('Select group of files you want to process: \n\t[0]: 2A \n\t[1]: R >>:', bool)]
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = end_name)

    #%% Plot things for just one heart
    # Get file to process and directories
    folder, df_file, file_num, blind = fcBasics.selectFile(df_dataset, end_name = end_name)
    if end_name == '2A': 
        filename = folder[0:-3]
    else: 
        filename = folder[2:]
    dORv = filename[9:10]
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = end_name)

    # Import df_results
    df_res = fcBasics.loadDF(filename = filename, file = 'ResultsDF', dir_results = dir_results)
    file_num = df_res[df_res['Folder']== filename+'_2A'].index.values[0]
    if dORv == 'D' or 'CJ' in filename:
        azimuth = -90
    else: 
        azimuth = 0
    elevation = df_res.loc[file_num,'ang_Heart']
    
    #% Get plot of cardiac jelly thickness heatmap with thickness range given as user input 
    vtk_in_dir = Path(directories[2])
    vtk_meshes = []
    # Import dictionaries
    [dicts] = fcBasics.loadDicts(filename = filename, dicts_name = ['dict_obj'], directories = [directories[0]])
    
    for name in vtk_in_dir.glob('*.vtk*'):
        mesh_name =  os.path.split(name)
        vtk_meshes.append(mesh_name[1][19:-4])
        
    obj_num = fcBasics.selectFromList (vtk_meshes, 'meshes','save rotating videos of')
    vtk_selected = list(vtk_meshes[i] for i in obj_num)
    print('Selected meshes: ',vtk_selected)
    
    # Change zoom or duration of videos
    zoom = 1; duration = 15
    ch_zoom = fcBasics.ask4input('Do you want to change the zoom of the meshes within the rotating videos? \n\t[0]: no, keep the default (zoom = 1)\n\t[1]: yes, please! >>:', bool)
    if ch_zoom: 
        print('- IMPORTANT NOTE: Videos in all versions of morphoHeart have been saved with a zoom of 1. \n  Making this change will make the meshes in the new videos look bigger (>1) or smaller (<1), so keep this in mind ;)')
        zoom = fcBasics.ask4input('Enter new zoom value: ', float)
    ch_duration = fcBasics.ask4input('Do you want to change the duration of the rotating videos? \n\t[0]: no, keep the default (15 sec)\n\t[1]: yes, please! >>:', bool)
    if ch_duration: 
        print('- IMPORTANT NOTE: Videos in all versions of morphoHeart have been saved with a duration of 15 sec. \n  Making this change will make your new videos rotate at a different speed, so keep this in mind ;)')
        duration = fcBasics.ask4input('Enter new duration value (sec): ', int)
    
    meshes = []
    min_max = []
    for i, m_vtk in enumerate(vtk_selected):
        if 'thickness' in m_vtk or 'Ball' in m_vtk:
            [[mesh], [thickness]] = fcMeshes.openThicknessMeshes(filename = filename,
                                                                meshes_names = [m_vtk], 
                                                                extension = 'vtk', dir_stl = directories[2], 
                                                                dir_txtNnpy = directories[1])
            if 'thickness' in m_vtk:
                scale_default = '0-25 um'; q_min = 0; q_max = 25
            elif 'Ball' in m_vtk:
                scale_default = '0-100 um'; q_min = 0; q_max = 100
                
            correct = False
            while not correct:
                min_val = format(min(thickness),'.2f')
                max_val = format(max(thickness),'.2f')
                scale_bar = fcBasics.ask4input('Do you want to set a range to the scale bar of -'+ m_vtk+'-? \n\t[0]: no, set according to mesh ('+min_val+'-'+max_val+'um)\n\t[1]: no, keep the default values ('+scale_default+')\n\t[2]: yes, I would like to set the range! >>>: ', int)
                if scale_bar == 1: 
                    min_max.append((q_min, q_max)); correct = True
                elif scale_bar == 2: 
                    q_min = fcBasics.ask4input('Minimum value for scale bar of -'+m_vtk+'- [um]: ',float)
                    q_max = fcBasics.ask4input('Maximum value for scale bar of -'+m_vtk+'- [um]: ',float)
                    min_max.append((q_min, q_max)); correct = True
                elif scale_bar == 0: 
                    min_max.append('-'); correct = True
                else: 
                    continue
            
            meshes.append(mesh)
        else: 
            [mesh] = fcMeshes.openMeshes(filename = filename, meshes_names = [m_vtk],
                                         extension = 'vtk', dir_stl = directories[2],
                                         alpha = [1], dict_colour = dicts['dict_colour'])
            meshes.append(mesh)
            min_max.append('-')
            
        settings.legendSize = .20
        text = str(filename); txt = Text2D(text, c=c, font=font)
        vp = Plotter(N=1, axes=10)
        vp.show(txt, mesh, at=0, azimuth = azimuth, interactive=True)
    
    # Save rotating videos
    rotAngle =  df_res.loc[file_num,'ang_Heart']
    fcMeshes.saveMultVideos(filename, info = vtk_selected, meshes4video = meshes, rangeThBall = min_max, 
                                rotAngle = rotAngle, dir2save = directories[4], 
                                dir_txtNnpy = directories[1], plotshow = False, alpha_cube = 0, 
                                duration = duration, zoom = zoom)

#%% Init
init = True

    