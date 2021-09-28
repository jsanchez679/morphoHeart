# -*- coding: utf-8 -*-
"""
morphoHeart - C. CUT ALL HEART LAYERS (ATRIUM/VENTRICLE, DORSAL/VENTRAL, LEFT/RIGHT) AND QUANTIFY THEM
Welcome to the fourth code of morphoHeart!
If you are running this code you must have already created all the volume reconstructions of all the tissue layers of 
the heart you are processing and extracted its centreline and are looking forward to get the cardiac jelly and ballooning 
heatmap of this heart! The objective of this code is to divide all the reconstructed heart tissue layers into sections
(atrium/ventricle, dorsal/ventral, left/right) and get all the heatmaps you need for your future analyses. 
At the end of this code you will have created the chamber meshes for each heart tissue layer, heatmap meshes for the 
external myocardium, endocardium and cardiac jelly, representing its thickness, and another heatmap mesh(es) for the 
internal (and external) myocardium representing its ballooning. If you wish, at the end there is also an option to save 
videos of all these meshes rotating, so that you can include them in your presentations! 

Happy sectioning and heatmapping! 

@author: Juliana Sanchez-Posada
Version: 13th April, 2021
"""

#%% Importing python packages
import os
# import numpy as np
from time import perf_counter
from datetime import datetime
# from vedo import *
from vedo import Plotter, Cube, settings, Text2D
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
    init = not bool(int(input('Do you want to execute the script all at once or run it by cells? \n\t[0]: all at once (recommended if you are already familiar with the script)\n\t[1]: by cells (recommended if you are NOT yet familiar with the script). >>>: ')))
    print("- Current working directory: {0}".format(os.getcwd()))
    if not init: 
        print('\nIMPORTANT NOTE:\n- Remember to start running from cell %Start C_CutAndMeasure.\n- NEVER run as an individual cell the cell called %Importing python packages')
    return root_path, init

root_path, init = setWorkingDir(os.getcwd())

c="k"; font= 'VTK';
save = True; plot = False; plotshow = False

#%% Start C_CutAndMeasure
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
    from morphoHeart_modules import morphoHeart_funcPlot as fcPlot
    tic = perf_counter()

    #% SELECT FILE AND GET METADATA
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
    if dORv == 'D' or 'CJ' in filename:
        azimuth = -90
    else: 
        azimuth = 0
    # Initialise variables
    dict_shapes = dict()
    txt = Text2D(filename, c="k", font= 'CallingCode')

    #% LOAD MESHES, CENTRELINES AND OBJECTS
    #   This section will load all the heart tissue layer meshes already created in previous scripts, as well as the
    #   centreline(s) dictionary(ies) and create a plot with all of the meshes.
    #   ================================================================================================================

    # Get existing cl_dictionaries
    [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes, dicts_cl] = fcBasics.import_dicts('mH_C', filename, directories)

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

    df_res = fcBasics.spAnalysis(df_res, file_num)
    
    #% Create ksplines, points, lines and centrelines
    #   This section will create spline(s) of the centreline(s) using the information from the loaded dictionary(ies) and
    #   return first a 3D interactive plot with the myocardium, endocardium and centreline(s) and second a new plot
    #   with two different visualisations of the maximum inscribed spheres with which the centreline(s) were calculated.
    #   Note: The resulting centreline(s) will comprise 300 points (important for cutting the tissue layers into chambers)
    #   ================================================================================================================

    kSplinesCuts = fcMeshes.createKSpls(dict_kspl, kspl_list = ['ksplCut4CL_inflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Ext.Endo(Cut)', 
                                                                'ksplCut4CL_outflow-Int.Myoc(Cut)', 'ksplCut4CL_inflow-Int.Myoc(Cut)'])
    sphCuts = fcMeshes.createSpheres(dict_pts, pts_list = ['sph_Cut4CL_inflow-Int.Myoc(Cut)', 'sph_Cut4CL_inflow-Ext.Endo(Cut)', 
                                                           'sph_Cut4CL_outflow-Int.Myoc(Cut)', 'sph_Cut4CL_outflow-Ext.Endo(Cut)'])
    kspl_CL, linLines, ksplSph_o, dict_kspl = fcMeshes.createCLs(dict_cl = dicts_cl, dict_pts = dict_pts, dict_kspl = dict_kspl,
                                                                     dict_planes = dict_planes, colors = ['deepskyblue', 'tomato'], myoc = m_myoc)
    
    sph_CL, sph_CL_colour, dict_shapes = fcMeshes.createColouredCL(dict_cl= dicts_cl, dict_shapes = dict_shapes)
    
    # Save measurements and shapes
    df_res = fcMeshes.addLinearMeas2df(df_res = df_res, file_num = file_num, lines = linLines, kspl_CL = kspl_CL)
    
    if plot:
        settings.legendSize = .20
        vp = Plotter(N=1, axes=10)
        vp.show(txt, kSplinesCuts, sphCuts, kspl_CL, ksplSph_o, linLines, m_myoc.alpha(0.01), m_endo.alpha(0.01), at=0, azimuth = azimuth, interactive=True)

        settings.legendSize = .20
        if len(dicts_cl) == 2:
            vp = Plotter(N=4, axes=13)
            vp.show(m_myoc.alpha(0.01), sph_CL_colour[0], txt, at=0)
            vp.show(m_endo.alpha(0.01), sph_CL_colour[1], at=1)
            vp.show(m_myoc.alpha(0.01), sph_CL[0], at=2)
            vp.show(m_endo.alpha(0.01), sph_CL[1], at=3, azimuth = azimuth, interactive=True)
        else:
            vp = Plotter(N=2, axes=13)
            vp.show(m_myoc.alpha(0.01), sph_CL_colour[0], txt, at=0)
            vp.show(m_myoc.alpha(0.01), sph_CL[0], at=1, azimuth = azimuth, interactive=True)
    
    #%% GET LAYERS THICKNESS HEATMAP ***
    #   This section of the code will extract the thickness heatmap of each heart tissue layer calculating the shortest
    #   distance from each of the points that make up the external mesh (external myocardium, endocardium or cardiac
    #   jelly) to the corresponding internal mesh. At the end, a 3D plot showing the three different thickness heatmaps
    #   will pop-up and the resulting meshes will be saved.
    #   ================================================================================================================

    # # Get layers thickness heatmap
    # print('\n- Extracting thickness heatmaps for all heart layers...')
    # # Cardiac Jelly
    # cj_thickness, m_cjTh, cj_minmax = fcMeshes.getDistance2Mesh(filename = filename, m_int = m_cjIn, m_ext = m_cjOut,
    #                                                     title = 'Cardiac Jelly Thickness', plotshow = False)
    # # Myocardium
    # myoc_thickness, m_myocTh, myoc_minmax = fcMeshes.getDistance2Mesh(filename = filename, m_int = m_myocInt, m_ext = m_myocExt,
    #                                                     title = 'Myoc.Thickness', plotshow = False)
    # # Endocardium
    # endo_thickness, m_endoTh, endo_minmax = fcMeshes.getDistance2Mesh(filename = filename, m_int = m_endoInt, m_ext = m_endoExt,
    #                                                     title = 'Endo.Thickness', plotshow = False)

    # if save:
    #     fcMeshes.saveThickness(filename = filename, arrays2save = [myoc_thickness, endo_thickness],
    #                             names = ['myoc_thickness','endo_thickness'], dir2save = directories[1])
    #     dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [m_myocTh, m_endoTh],
    #                             names = ['myoc_thickness','endo_thickness'],
    #                             dict_colour = dict_colour, dir_stl = directories[2], extension = 'vtk')
    
    # % NOTE: Uncomment in case you have already run this section and just want to re-load results
    [[m_cjTh, m_myocTh, m_endoTh], [cj_thickness, myoc_thickness, endo_thickness]] = fcMeshes.openThicknessMeshes(filename = filename, 
                                                                                          meshes_names = ['cj_thickness','myoc_thickness','endo_thickness'], 
                                                                                          extension = 'vtk', dir_stl = directories[2], dir_txtNnpy = directories[1]); 
    cj_minmax = (min(cj_thickness), max(cj_thickness))
    myoc_minmax = (min(myoc_thickness), max(myoc_thickness))
    endo_minmax = (min(endo_thickness), max(endo_thickness))
    plot = True
    if plot:
        settings.legendSize = .20
        vp = Plotter(N=3, axes=10)
        vp.show(m_cjTh.alpha(1), scale_cube, txt, at=0)
        vp.show(m_myocTh.alpha(1), scale_cube, at=1)
        vp.show(m_endoTh.alpha(1), scale_cube, at=2, zoom=2, azimuth = azimuth, elevation = 0, interactive=True)
    
    #% GET BALLOONING HEATMAPS 
    #% NOTE: Uncomment in case you have already run this section and just want to re-load results
    sph_ballonning = fcMeshes.sphInSpline(kspl_CL = kspl_CL[0], name = 'sphs_ball', every = 0.6)
    q_selectMeshes = fcBasics.ask4input('Select the meshes from which you would like to extract the ballooning heatmaps \n\t\t[0]: Int.Myocardium/[1]: Ext.Endocardium/[2]: both: ', int)
    meshes_ball = []
    if q_selectMeshes in [0,2]:
        [[m_myocIntBall], [myoc_intBall]] = fcMeshes.openThicknessMeshes(filename = filename, meshes_names = ['myoc_intBall'], 
                                                                      extension = 'vtk', dir_stl = directories[2], dir_txtNnpy = directories[1]); 
        myocInt_minmax = (min(myoc_intBall), max(myoc_intBall)); meshes_ball.append(m_myocIntBall)
        vp = Plotter(N=1, axes = 13)
        vp.show(m_myocIntBall, at=0, interactive=True)
    if q_selectMeshes in [1,2]:
        [[m_myocExtBall], [myoc_extBall]] = fcMeshes.openThicknessMeshes(filename = filename, meshes_names = ['myoc_extBall'], 
                                                                      extension = 'vtk', dir_stl = directories[2], dir_txtNnpy = directories[1]); 
        myocExt_minmax = (min(myoc_extBall), max(myoc_extBall)); meshes_ball.append(m_myocExtBall)
        vp = Plotter(N=1, axes = 13)
        vp.show(m_myocExtBall, at=0, interactive=True)
        
    if plot:
        settings.legendSize = .2
        if q_selectMeshes in [0,1]: 
            vp = Plotter(N=4, axes = 13)
            vp.show(m_myocInt.alpha(0.01), sph_ballonning, at=0)
            vp.show(m_myocInt.alpha(0.01), m_myocExt.alpha(0.01), sph_ballonning, at=1)
            vp.show(meshes_ball[0].alpha(1), sph_ballonning, at=2)
            vp.show(m_cjTh.alpha(1), scale_cube, at=3, azimuth = azimuth, elevation = 0, zoom = 2, interactive = True)
        elif q_selectMeshes == 2: 
            vp = Plotter(N=6, axes = 13)
            vp.show(m_myocInt.alpha(0.01), sph_ballonning, at=0)
            vp.show(m_myocExt.alpha(0.01), sph_ballonning, at=1)
            vp.show(m_myocInt.alpha(0.01), m_myocExt.alpha(0.01), sph_ballonning, at=2)
            vp.show(meshes_ball[0].alpha(1), sph_ballonning, at=3)
            vp.show(meshes_ball[1].alpha(1), sph_ballonning, at=4)
            vp.show(m_cjTh.alpha(1), scale_cube, at=5, azimuth = azimuth, elevation = 0, zoom = 2, interactive = True)
            

    #%% CUT HEART TISSUE LAYERS INTO CHAMBERS
    #% Calculating Myoc just in case
    num_pt = dict_pts['numPt_CLChamberCut']; atr_meshes = []; vent_meshes = []
    atr_meshes, vent_meshes, dict_shapes, s3_cyl  = fcMeshes.getChamberMeshes(filename = filename,
                                    end_name = ['ch0_cut'], names2cut = ['Myoc'],
                                    kspl_CL = kspl_CL[0], num_pt = num_pt, atr_meshes = atr_meshes, vent_meshes = vent_meshes,
                                    dir_txtNnpy = directories[1], dict_shapes = dict_shapes, dict_pts = dict_pts, 
                                    resolution = res, s3_cyl = [], plotshow = True, mesh2cut = m_myoc)
    m_atrMyoc = atr_meshes[0]; m_ventMyoc = vent_meshes[0]
    settings.legendSize = .20
    vp = Plotter(N=1, axes = 10)
    vp.show(m_atrMyoc.color('lightseagreen'), m_ventMyoc.color('darkturquoise'), at = 0,  zoom = 2.2, azimuth = azimuth, interactive=True)
    
    dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [m_atrMyoc, m_ventMyoc],
                                      names = [ 'myoc_atr', 'myoc_vent'],
                                      dict_colour = dict_colour, dir_stl = directories[2], extension = 'vtk')
    dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
                                              names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
    
    m_atrMyocBU = m_atrMyoc.clone()
    m_ventMyocBU = m_ventMyoc.clone()
    
    #%% GET CHAMBERS AND HEART ORIENTATION
    #   This section of the code will measure the chambers and heart orientations and save them in the dataframe with
    #   the rest of the measurements already made. At the end, a 3D plot showing the chambers orientation will pop-up.
    #   ================================================================================================================

    [m_atrMyoc, m_ventMyoc] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc_atr','myoc_vent'],
                                                                  extension = 'vtk', dir_stl = directories[2],
                                                                  alpha = [1,1], dict_colour = dict_colour)
    num_pt = dict_pts['numPt_CLChamberCut']
    
    sph_orient, lines_orient, dict_pts, dict_kspl, df_res = fcMeshes.getChambersOrientation(filename = filename, file_num = file_num, 
                                                                    num_pt = num_pt, kspl_CL2use = kspl_CL[0], distFromCl = 50,
                                                                    myoc_meshes = [m_atrMyoc, m_ventMyoc], linLine = linLines[0],
                                                                    dict_pts = dict_pts, dict_kspl = dict_kspl, df_res = df_res)
    elevation = df_res.loc[file_num,'ang_HeartS']
    orient_atr, orient_vent, orient_atrX, orient_ventX, linLineX = lines_orient
    
    #%% 
    if save:
        # Save measurements dataframe
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
    
    #% Open ext myoc atrium and ventricle
    # [m_atrExtMyo, m_ventExtMyo] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myocExt_atr', 'myocExt_vent'],
    #                                extension = 'vtk', dir_stl = directories[2],
    #                                alpha = [1,1], dict_colour = dict_colour)
    
    df_res, dict_shapes = fcMeshes.getChambersEllipsoid(filename = filename, df_res = df_res, file_num = file_num, 
                                  lines_orient = lines_orient, meshes =(m_atrMyoc, m_ventMyoc), dict_shapes = dict_shapes)

    fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
    fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = os.path.join(dir_data2Analyse, 'R_All','df_all','df_meas'))
    # Append all dicts to one object dict
    dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
                                              names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
        
    #%
    # Calculating Ext.Endo chambers
    num_pt = dict_pts['numPt_CLChamberCut']; atr_meshes = []; vent_meshes = []
    atr_meshes, vent_meshes, dict_shapes, s3_cyl  = fcMeshes.getChamberMeshes(filename = filename,
                                    end_name = ['ch1_cut_ext'], names2cut = ['Ext.Endo'],
                                    kspl_CL = kspl_CL[0], num_pt = num_pt, atr_meshes = atr_meshes, vent_meshes = vent_meshes,
                                    dir_txtNnpy = directories[1], dict_shapes = dict_shapes, dict_pts = dict_pts, 
                                    resolution = res, s3_cyl = [], plotshow = True, mesh2cut = m_endoExt)
    m_atrExtEndo = atr_meshes[0]; m_ventExtEndo = vent_meshes[0]
    settings.legendSize = .20
    vp = Plotter(N=1, axes = 10)
    vp.show(m_atrExtEndo.color('cyan'), m_ventExtEndo.color('springgreen'), at = 0,  zoom = 2.2, azimuth = azimuth, interactive=True)
    
    dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [m_atrExtEndo, m_ventExtEndo],
                                      names = [ 'endoExt_atr', 'endoExt_vent'],
                                      dict_colour = dict_colour, dir_stl = directories[2], extension = 'vtk')
    
    #%%
    # # Calculating Ext.Myoc chambers
    # num_pt = dict_pts['numPt_CLChamberCut']; atr_meshes = []; vent_meshes = []
    # atr_meshes, vent_meshes, dict_shapes, s3_cyl  = fcMeshes.getChamberMeshes(filename = filename,
    #                                 end_name = ['ch0_cut_ext'], names2cut = ['Ext.Myoc'],
    #                                 kspl_CL = kspl_CL[0], num_pt = num_pt, atr_meshes = atr_meshes, vent_meshes = vent_meshes,
    #                                 dir_txtNnpy = directories[1], dict_shapes = dict_shapes, dict_pts = dict_pts, 
    #                                 resolution = res, s3_cyl = [], plotshow = True, mesh2cut = m_myocExt)
    # m_atrExtMyo = atr_meshes[0]; m_ventExtMyo = vent_meshes[0]
    # settings.legendSize = .20
    # vp = Plotter(N=1, axes = 10)
    # vp.show(m_atrExtMyo.color('tomato'), m_ventExtMyo.color('navy'), at = 0,  zoom = 2.2, azimuth = azimuth, interactive=True)
    
    # dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [m_atrExtMyo, m_ventExtMyo],
    #                                   names = [ 'myocExt_atr', 'myocExt_vent'],
    #                                   dict_colour = dict_colour, dir_stl = directories[2], extension = 'vtk')
    
    # dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
    #                                          names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
    
    #%%

    #% CREATE PLANES AND RIBBON TO DIVIDE HEART INTO SECTIONS (Dors/Vent, Left/Right)
    #   This section of the code will create planes to divide each of the chambers into dorso-ventral regions and will
    #   extend dorso-ventrally the centreline to create a ribbon that divides the heart into left and right. 3D
    #   interactive plots will pop-up showing first the dorso-ventral planes created and second the centreline ribbon.
    #   Once the ribbon has been created it will be used to cut the external cardiac jelly mesh/internal
    #   myocardium mesh (and the external myocardium mesh). Resulting meshes will be shown and the user will be asked 
    #   to confirm the resulting left-right classification. 
    #   ================================================================================================================

    # Create Coronal Planes for each chamber
    # [pl_DnV_Atr, pl_DnV_Vent], dict_planes = fcMeshes.createDVPlanes(filename = filename, sph_orient = sph_orient, mesh = m_myoc,
    #                                               kspl_CL = kspl_CL[0], orient_lines = [orient_atr, orient_vent], dict_planes = dict_planes)
    # Get centreline ribbon
    cl_ribbon, kspl_ext, dict_kspl, dict_shapes, dict_planes = fcMeshes.createCLRibbon(filename = filename, file_num = file_num, 
                                                                                       df_res = df_res, kspl_CL2use = kspl_CL[0], 
                                                                                       linLine = linLines[0], mesh = m_myoc, 
                                                                                       dict_kspl = dict_kspl, dict_shapes = dict_shapes, 
                                                                                       dict_planes = dict_planes)
    
    # Divide meshes using ribbon (Left/Right) 
    # First divide the cardiac jelly and save the volumes of it's left and right sides
    [m_cjLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_cj], cl_ribbon = cl_ribbon, 
                                         file_num = file_num, df_res = df_res, colors = ['salmon', 'brown'])
    vp = Plotter(N=2, axes = 10)
    vp.show(m_cjLnR[0], at = 0)
    vp.show(m_cjLnR[1], at=1, interactive = True)
    m_cjThLnR = m_cjLnR
    df_res = fcMeshes.addLayersVolume2df (df_res = df_res, file_num = file_num, meshes = [m_cj]+m_cjLnR, 
                                          names = ['CJ_total','CJ.Left', 'CJ.Right'])
    fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
    
    # Now divide the tissue thickness meshes (myocardium, endocardium and cardiac jelly)
    myoc_and_endo = fcBasics.ask4input('Do you want to run thickness analysis of the myocardium and endocardium as well? \n\t[0]: no, just of the cardiac jelly \n\t[1]: yes, please! >> : ', bool)
    # if q_selectMeshes in [0,2]:
    #     [m_cjThLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_cjTh], cl_ribbon = cl_ribbon, file_num = file_num, 
    #                                                                                    df_res = df_res)
    # if q_selectMeshes in [1,2]:
    #     [m_cjThLnR, m_myocExtLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_cjTh, m_myocExtBall], cl_ribbon = cl_ribbon, 
    #                                                          file_num = file_num, df_res = df_res)
    if myoc_and_endo:
        [m_myocThLnR, m_endoThLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_myocTh, m_endoTh], 
                                                              cl_ribbon = cl_ribbon, file_num = file_num, 
                                                              df_res = df_res)

    #% CLASSIFY THICKNESS AND BALLOONING POINTS
    #   Now that we have defined the planes and ribbon that section the heart into atrial-ventricular, dorsal-ventral,
    #   left-right sections of the heart, this section of the code will classify each of the points that make up the
    #   external cardiac jelly/internal myocardium (meshes in which the cardiac jelly thickness heatmap and the
    #   internal myocardium ballooning heatmap respectively had been mapped) (and the external myocardium mesh), 
    #   resulting in a dataframe containing all this information. 
    #   Once the classification algorithm has run, a 3D interactive plot will pop-up showing the
    #   resulting classification for a subsample of the points color-coded to each region. At the end, the resulting
    #   dataframe will be saved for future analysis.
    #   ================================================================================================================

    dir2save_df = fcBasics.new_dir(dir_results, 'csv_all'); 
    [m_atrExtCJ] = fcMeshes.openMeshes(filename = filename, meshes_names = ['cjExt_atr'],
                                                                  extension = 'vtk', dir_stl = directories[2],
                                                                  alpha = [1], dict_colour = dict_colour)

    if q_selectMeshes in [0,2]:
        df_cjThNmyocIntBall = fcMeshes.classifyHeartPts(filename = filename, df_res = df_res, file_num = file_num,
                                            dict_planes = dict_planes,
                                            m_whole = m_cjTh, m_left = m_cjThLnR[0], m_atr = m_atrExtCJ, 
                                            data = [cj_thickness, myoc_intBall], 
                                            names_data = ['cj_thickness', 'myoc_intBall'], plot_show = True)
        if save:
            fcBasics.saveDF(filename = filename, df2save = df_cjThNmyocIntBall, df_name = 'df_cjThNmyocIntBall',
                            dir2save = dir2save_df)
    # if q_selectMeshes in [1,2]:
    #     df_myocExtBall = fcMeshes.classifyHeartPts(filename = filename, dict_planes = dict_planes,
    #                                         m_whole = m_myocExtBall, m_left = m_myocExtLnR[0], m_atr = m_atrExtMyo, 
    #                                         data = [myoc_extBall], names_data = ['myoc_extBall'], plot_show = True)
    #     if save:
    #         fcBasics.saveDF(filename = filename, df2save = df_myocExtBall, df_name = 'df_myocExtBall',
    #                         dir2save = dir2save_df)
            
    if myoc_and_endo: 
        [m_atrExtMyo, m_atrExtEndo] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myocExt_atr', 'endoExt_atr'],
                                                                  extension = 'vtk', dir_stl = directories[2],
                                                                  alpha = [1,1], dict_colour = dict_colour)
        df_myocTh = fcMeshes.classifyHeartPts(filename = filename, df_res = df_res, file_num = file_num,
                                            dict_planes = dict_planes,
                                            m_whole = m_myocTh, m_left = m_myocThLnR[0], m_atr = m_atrExtMyo, 
                                            data = [myoc_thickness], names_data = ['myoc_thickness'], plot_show = True)
        df_endoTh = fcMeshes.classifyHeartPts(filename = filename, df_res = df_res, file_num = file_num,
                                            dict_planes = dict_planes,
                                            m_whole = m_endoTh, m_left = m_endoThLnR[0], m_atr = m_atrExtEndo, 
                                            data = [endo_thickness], names_data = ['endo_thickness'], plot_show = True)
        if save:
            fcBasics.saveDF(filename = filename, df2save = df_myocTh, df_name = 'df_myocTh', dir2save = dir2save_df)
            fcBasics.saveDF(filename = filename, df2save = df_endoTh, df_name = 'df_endoTh', dir2save = dir2save_df)

    #% SAVE ALL
    #   This section allows the user to save all the meshes, objects and dataframes that have been created. It will also
    #   ask the user if he/she wants to save videos of all/some of the created meshes (including those with heatmaps)
    #   for which the user will be able to define a scale bar range if desired. 
    #   ================================================================================================================

    if save:
        # Append all dicts to one object dict
        dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
                                              names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
        # Save filled dataframe with measured data
        cols2del = ['cyl_radius-myoc_atr', 'cyl_height-myoc_atr','cyl_radius-myoc_vent', 'cyl_height-myoc_vent',
                                    'ang_Heart', 'ang_Atr', 'ang_Vent', 'ang_BtwChambers']
        for col in cols2del:
            if col in list(df_res.columns):
                print('- ',col, ' was deleted!')
                df_res = df_res.drop([col], axis=1)
        
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
        # if os.path.isdir(os.path.join(dir_data2Analyse, 'R_All', 'df_all','df_meas')) == False:
        #     dir_R_all = os.path.join(dir_data2Analyse, 'R_All'); dir_df_meas = os.path.join(dir_R_all, 'df_meas')
        #     os.mkdir(dir_R_all); os.mkdir(dir_df_meas)
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = os.path.join(dir_data2Analyse, 'R_All','df_all','df_meas'))
        # Save rotating videos
    #     if q_selectMeshes == 0:
    #         names4video = ['myoc','endo','cj', 'cj_thickness','myoc_thickness','endo_thickness','myoc_intBall','myoc_endo']
    #         meshes4video = [m_myoc,m_endo,m_cj,m_cjTh,m_myocTh,m_endoTh, m_myocIntBall.alpha(1),[m_myoc.alpha(0.01), m_endo.alpha(0.01), kspl_CL[0]]]
    #         range4video = ['-', '-', '-', cj_minmax, myoc_minmax, endo_minmax, myocInt_minmax, '']
    #     elif q_selectMeshes == 1:
    #         names4video = ['myoc','endo','cj', 'cj_thickness','myoc_thickness','endo_thickness','myoc_extBall','myoc_endo']
    #         meshes4video = [m_myoc,m_endo,m_cj,m_cjTh,m_myocTh,m_endoTh, m_myocExtBall.alpha(1),[m_myoc.alpha(0.01), m_endo.alpha(0.01), kspl_CL[0]]]
    #         range4video = ['-', '-', '-', cj_minmax, myoc_minmax, endo_minmax, myocExt_minmax,'']
    #     elif q_selectMeshes == 2:
    #         names4video = ['myoc','endo','cj', 'cj_thickness','myoc_thickness','endo_thickness','myoc_intBall', 'myoc_extBall','myoc_endo']
    #         meshes4video = [m_myoc,m_endo,m_cj,m_cjTh,m_myocTh,m_endoTh, m_myocIntBall.alpha(1), m_myocExtBall.alpha(1),[m_myoc.alpha(0.01), m_endo.alpha(0.01), kspl_CL[0]]]
    #         range4video = ['-', '-', '-', cj_minmax, myoc_minmax, endo_minmax, myocInt_minmax, myocExt_minmax,'']
            
    #     names4video, meshes4video, rangeThBall = fcPlot.selectMeshes4Video(names = names4video, meshes = meshes4video, ranges = range4video)
    #     sure = fcBasics.ask4input('Are you sure about orientation for videos? >>: ', bool)
        # if sure: 
        #     dir4videos = fcBasics.new_dir(directories[4], 'videos')
        #     rotAngle =  df_res.loc[file_num,'ang_Heart']#0
        #     fcMeshes.saveMultVideos(filename, info = names4video, meshes4video = meshes4video, rangeThBall = rangeThBall, 
        #                             rotAngle = rotAngle, dir2save = dir4videos, 
        #                             dir_txtNnpy = directories[1], plotshow = False, alpha_cube = 0)
        
    #     toc = perf_counter()
    #     fcBasics.printTime(tic, toc, 'Cut and Measure')
        
    
    # if save:
    #     # Append all dicts to one object dict
    #     dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
    #                                          names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
    #     # Save filled dataframe with measured data
    #     fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results, name = 'ResultsDFf')
    #     dir_df_meas = fcBasics.new_dir(fcBasics.new_dir(fcBasics.new_dir(dir_data2Analyse, 'R_All'), 'df_all'), 'df_meas')
    #     # if os.path.isdir(os.path.join(dir_data2Analyse, 'R_All', 'df_all', 'df_meas')) == False:
    #     #     dir_R_all = os.path.join(dir_data2Analyse, 'R_All'); dir_df_all = os.path.join(dir_R_all, 'df_all')
    #     #     dir_df_meas = os.path.join(dir_df_all, 'df_meas'); os.mkdir(dir_R_all); os.mkdir(dir_df_all); os.mkdir(dir_df_meas)
    #     fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_df_meas, name = 'ResultsDFf')

        
#%% Init
init = True

# dict_colour = fcMeshes.saveMeshes(filename = filename, 
#                                           meshes = [m_myoc, m_endo, m_cjTh],
#                                           names = ['myoc_stl', 'endo_stl', 'cj_stl'],
#                                           dict_colour = dict_colour, dir_stl = directories[2], extension = 'stl')

# mTh_names = ['cj_thickness','myoc_thickness','endo_thickness','myoc_intBall']
# [m_thAll, colour_thAll] = fcMeshes.openThicknessMeshes(filename = filename, meshes_names = mTh_names, extension = 'vtk',
#                               dir_stl = directories[2], dir_txtNnpy = directories[1])
# m_cjTh, m_myocTh, m_endoTh, m_myocIntBall = m_thAll1
# cj_thickness, myoc_thickness, endo_thickness, myoc_intBall = colour_thAll
# cj_minmax = (cj_thickness.min(), cj_thickness.max())
# myoc_minmax = (myoc_thickness.min(), myoc_thickness.max())
# endo_minmax = (endo_thickness.min(), endo_thickness.max())
# myocInt_minmax = (myoc_intBall.min(), myoc_intBall.max())


#%%
# import vedo
# from vedo.applications import FreeHandCutPlotter

# vedo.settings.useParallelProjection = True
# plt = FreeHandCutPlotter(m_myocTh).addHoverLegend()
# #plt.init(some_list_of_initial_pts) #optional!
# plt.start(axes=1, bg2='lightblue').close()
#%%
# import trimesh
# import numpy as np
# from vedo import show, Plane, printc, download, embedWindow
# embedWindow(False)

# # load the mesh from filename, file objects are also supported
# f = download('https://github.com/mikedh/trimesh/raw/master/models/featuretype.STL')
# mesh = trimesh.load_mesh(f)

# # get a single cross section of the mesh
# txt = 'cross section of the mesh'
# mslice = mesh.section(plane_origin=mesh.centroid, plane_normal=[0,0,1])

# pl = Plane(mesh.centroid, normal=[0,0,1], sx=6, sy=4, c='green', alpha=0.3)

# slice_2D, to_3D = mslice.to_planar()

# # show objects on N=2 non-synced renderers:
# show([(mesh, pl), (slice_2D, txt)], N=2, sharecam=False, axes=7).close()
