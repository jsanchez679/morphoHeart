# -*- coding: utf-8 -*-
"""
TO GET IMAGES FOR THESIS
Created on Thu Mar 31 13:37:22 2022

@author: mdp18js
"""

#%% Importing python packages
import os
# import numpy as np
from time import perf_counter
# from datetime import datetime
# from vedo import *
from vedo import Plotter, Cube, settings, Text2D, Cylinder
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
save = True; plot = True; plotshow = False

#%% Start C_CutAndMeasure
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
    from morphoHeart_modules import morphoHeart_funcPlot as fcPlot
    tic = perf_counter()

    #%% SELECT FILE AND GET METADATA
    #   This section allows the user to select file to process, get its main directories and metadata,
    #   and define some properties
    #   ================================================================================================================

    # Get main directories
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = 'R', out_type = 'xlsx')
    # Get file to process and directories
    folder, df_file, file_num, blind = fcBasics.selectFile(df_dataset, end_name = 'R'); filename = folder[2:]; dORv = filename[9:10]
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
    # Initialise variables
    dict_shapes = dict()
    txt = Text2D(filename, c="k", font= 'CallingCode')
    rotAngle =  df_res.loc[file_num,'ang_HeartS']
    
    #% LOAD MESHES, CENTRELINES AND OBJECTS
    #   This section will load all the heart tissue layer meshes already created in previous scripts, as well as the
    #   centreline(s) dictionary(ies) and create a plot with all of the meshes.
    #   ================================================================================================================

    # Get existing cl_dictionaries
    [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes, dicts_cl] = fcBasics.import_dicts('mH_C', filename, directories)

    # Import meshes
    in_meshes01 = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc','endo','cj','cj_out','cj_in'],
                                                                  extension = 'vtk', dir_stl = directories[2],
                                                                  alpha = [0.05,0.05,0.05,1,1], dict_colour = dict_colour)
    out_meshes01 = []
    for mesh in in_meshes01:
        mesh = mesh.rotateX(rotAngle)
        out_meshes01.append(mesh)
    m_myoc, m_endo, m_cj, m_cjOut, m_cjIn = out_meshes01
    scale_cube = Cube(pos=m_myoc.centerOfMass(), side=350, c='white', alpha=0.01).legend('Cube')
    
    in_meshes02 = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc_int','myoc_ext','endo_int', 'endo_ext'],
                                                                  extension = 'vtk', dir_stl = directories[2],
                                                                  alpha = [1,1,1,1], dict_colour = dict_colour)
    out_meshes02 = []
    for mesh in in_meshes02:
        mesh = mesh.rotateX(rotAngle)
        out_meshes02.append(mesh)
    m_myocInt, m_myocExt, m_endoInt, m_endoExt = out_meshes02
    m_myocInt.color('crimson')
    m_myocExt.color('gold')
    m_endoInt.color('deep pink')
    m_endoExt.color('light sky blue')
    
    # Plot meshes
    if plot:
        text = str(filename); txt = Text2D(text, c=c, font=font)
        settings.legendSize = .20
        vp = Plotter(N=6, axes=13)
        vp.show(m_myoc.alpha(0.01),scale_cube, txt, at=0)
        vp.show(m_endo.alpha(0.01), scale_cube, at=1)
        vp.show(m_cj.alpha(0.01), scale_cube,  at=2)
        vp.show(m_endo.clone().alpha(0.01), m_cj.clone().alpha(0.05), scale_cube, at=3)
        vp.show(m_cjOut, scale_cube, at=4)
        vp.show(m_myoc.clone().alpha(0.01), m_endo.clone().alpha(0.01), m_cj.clone().alpha(0.01), scale_cube, at=5, 
                zoom = 1.6, azimuth = 0, interactive=True)

        text = str(filename); txt = Text2D(text, c=c, font=font)
        settings.legendSize = .20
        vp = Plotter(N=6, axes=13)
        vp.show(m_myoc.alpha(0.01),scale_cube, txt, at=0)
        vp.show(m_endo.alpha(0.01), scale_cube, at=1)
        vp.show(m_cj.alpha(0.01), scale_cube,  at=2)
        vp.show(m_endo.clone().alpha(0.01), m_cj.clone().alpha(0.05), scale_cube, at=3)
        vp.show(m_cjOut, scale_cube, at=4)
        vp.show(m_myoc.clone().alpha(0.01), m_endo.clone().alpha(0.01), m_cj.clone().alpha(0.01), scale_cube, at=5, zoom = 1.6, azimuth = 90, interactive=True)


        text = str(filename); txt = Text2D(text, c=c, font=font)
        settings.legendSize = .20
        vp = Plotter(N=6, axes=13)
        vp.show(m_myoc.alpha(0.01),scale_cube, txt, at=0)
        vp.show(m_myocInt.alpha(0.01), scale_cube, at=1)
        vp.show(m_myocExt.alpha(0.01), scale_cube, at=2)
        vp.show(m_endo.alpha(0.01), scale_cube, at=3)
        vp.show(m_endoInt.alpha(0.01), scale_cube, at=4)
        vp.show(m_endoExt.alpha(0.01), scale_cube, at=5, zoom = 1.6, azimuth = azimuth, interactive=True)
        
        text = str(filename); txt = Text2D(text, c=c, font=font)
        settings.legendSize = .20
        vp = Plotter(N=1, axes=13)
        vp.show(m_myoc.clone().alpha(0.01), m_endo.clone().alpha(0.01), scale_cube, at=0, zoom = 1.6, azimuth = 0, interactive=True)




    #%% Create ksplines, points, lines and centrelines

    kSplinesCuts = fcMeshes.createKSpls(dict_kspl, kspl_list = ['ksplCut4CL_inflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Ext.Endo(Cut)', 
                                                                'ksplCut4CL_outflow-Int.Myoc(Cut)', 'ksplCut4CL_inflow-Int.Myoc(Cut)'])
    sphCuts = fcMeshes.createSpheres(dict_pts, pts_list = ['sph_Cut4CL_inflow-Int.Myoc(Cut)', 'sph_Cut4CL_inflow-Ext.Endo(Cut)', 
                                                           'sph_Cut4CL_outflow-Int.Myoc(Cut)', 'sph_Cut4CL_outflow-Ext.Endo(Cut)'])
    kspl_CL, linLines, ksplSph_o, dict_kspl = fcMeshes.createCLs(filename = filename, dict_cl = dicts_cl, dict_pts = dict_pts, 
                                                                 dict_kspl = dict_kspl, dict_planes = dict_planes, 
                                                                 colors = ['deepskyblue', 'tomato'], myoc = m_myoc, 
                                                                 dir_stl = directories[3])
    sph_CL, sph_CL_colour, dict_shapes = fcMeshes.createColouredCL(dict_cl= dicts_cl, dict_shapes = dict_shapes)
    
    in_meshes03 = [sphCuts[0], sphCuts[1], kspl_CL[0], linLines[0], ksplSph_o[0], sph_CL_colour[0], sph_CL[0]]
    out_meshes03 = []
    for mesh in in_meshes03:
        mesh = mesh.rotateX(rotAngle)
        out_meshes03.append(mesh)
    sphCuts0, sphCuts1, kspl_CL, linLines, ksplSph_o0, sph_CL_colour, sph_CL = out_meshes03

    lw = 8
    if plot:
        # settings.legendSize = .20
        # vp = Plotter(N=1, axes=13)
        # vp.show(txt, scale_cube, kspl_CL.color('corn flower blue').lw(lw), #ksplSph_o0.color('darkorange').lw(lw), 
        #         linLines.color('lime').lw(lw), m_myoc.alpha(0.01), at=0, azimuth = 0, zoom = 1.6, interactive=True)
        # settings.legendSize = .20
        # vp = Plotter(N=1, axes=13)
        # vp.show(txt, scale_cube, kspl_CL.color('corn flower blue').lw(lw), #ksplSph_o0.color('darkorange').lw(lw), 
        #         linLines.color('lime').lw(lw), m_myoc.alpha(0.01), at=0, azimuth = 90, zoom = 1.6, interactive=True)

        # settings.legendSize = .20
        # vp = Plotter(N=4, axes=13)
        # vp.show(m_myoc.alpha(0.01), sph_CL_colour, scale_cube, txt, at=0)
        # vp.show(m_myoc.alpha(0.01), sph_CL, scale_cube, at=1) 
        # vp.show(m_myocInt.alpha(0.01), sph_CL_colour, scale_cube, txt, at=2)
        # vp.show(m_myocInt.alpha(0.01), sph_CL, scale_cube, at=3, azimuth = azimuth, zoom = 1.6, interactive=True)

    
        # settings.legendSize = .20
        # vp = Plotter(N=4, axes=13)
        # vp.show(m_myocExt.alpha(0.01),m_myocInt.alpha(0.01),scale_cube, txt, at=0)
        # vp.show(m_endoExt.alpha(0.01),m_endoInt.alpha(0.01),scale_cube, at=1)
        # vp.show(m_endoExt.alpha(0.01),m_myocInt.alpha(0.01),scale_cube, at=2)
        # vp.show(scale_cube, kspl_CL.color('corn flower blue').lw(lw), 
        #         m_myocInt.alpha(0.01), at=3, azimuth = azimuth, zoom = 2, interactive=True)
        
        text = str(filename); txt = Text2D(text, c=c, font=font)
        settings.legendSize = .20
        vp = Plotter(N=1, axes=13)
        vp.show(m_myoc.clone().alpha(0.01), m_endo.clone().alpha(0.01), kspl_CL.color('corn flower blue').lw(lw),
                linLines.color('lime').lw(lw), scale_cube, at=0, zoom = 1.6, azimuth = 0, interactive=True)


    #%% GET LAYERS THICKNESS HEATMAP ***
    [[m_cjTh, m_myocTh, m_endoTh], [cj_thickness, myoc_thickness, endo_thickness]] = fcMeshes.openThicknessMeshes(filename = filename, 
                                                                                          meshes_names = ['cj_thickness','myoc_thickness','endo_thickness'], 
                                                                                          extension = 'vtk', dir_stl = directories[2], dir_txtNnpy = directories[1]); 
    [[m_myocIntBall], [myoc_intBall]] = fcMeshes.openThicknessMeshes(filename = filename, meshes_names = ['myoc_intBall'], 
                                                                      extension = 'vtk', dir_stl = directories[2], dir_txtNnpy = directories[1]); 
    
    m_cjTh.pointColors(cj_thickness, cmap="turbo", vmin=0, vmax=20)
    m_cjTh.addScalarBar()
    m_cjTh.alpha(1)
    m_cjTh.mapper().SetScalarRange(0,20)
    
    m_myocTh.pointColors(myoc_thickness, cmap="turbo", vmin=0, vmax=12)
    m_myocTh.addScalarBar()
    m_myocTh.alpha(1)
    m_myocTh.mapper().SetScalarRange(0,12)
    
    m_endoTh.pointColors(endo_thickness, cmap="turbo", vmin=0, vmax=12)
    m_endoTh.addScalarBar()
    m_endoTh.alpha(1)
    m_endoTh.mapper().SetScalarRange(0,12)
    
    m_myocIntBall.pointColors(myoc_intBall, cmap="turbo", vmin=0, vmax=60)
    m_myocIntBall.addScalarBar()
    m_myocIntBall.alpha(1)
    m_myocIntBall.mapper().SetScalarRange(0,60)
    
    in_meshes04 = [m_cjTh, m_myocTh, m_endoTh, m_myocIntBall]
    out_meshes04 = []
    for mesh in in_meshes04:
        mesh = mesh.rotateX(rotAngle)
        out_meshes04.append(mesh)
    m_cjTh, m_myocTh, m_endoTh, m_myocIntBall = out_meshes04
    

    if plot:
        settings.legendSize = .20
        vp = Plotter(N=4, axes=13)
        vp.show(m_cjTh.alpha(1), scale_cube, txt, at=0)
        vp.show(m_myocTh.alpha(1), scale_cube, at=1)
        vp.show(m_endoTh.alpha(1), scale_cube, at=2)
        vp.show(m_myocIntBall.alpha(1), scale_cube, at=3, zoom=1.6, azimuth = azimuth, elevation = 0, interactive=True)
    
    #%% CUT HEART TISSUE LAYERS INTO CHAMBERS
    m_atrVent = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc_atr', 'myoc_vent','cj_atr','cj_vent'],
                                    extension = 'vtk', dir_stl = directories[2],
                                    alpha = [1,1,1,1], dict_colour = dict_colour)
    cyl_chamber = dict_shapes['cyl2CutChambers_final']
    disc = Cylinder(pos = cyl_chamber['cyl_centre'],r = cyl_chamber['radius_max'], height = 2*0.225, 
                    axis = cyl_chamber['cyl_axis'], c = 'purple', cap = True, res = 300)
    
    in_meshes05 = m_atrVent+[disc]
    out_meshes05 = []
    for mesh in in_meshes05:
        mesh = mesh.rotateX(rotAngle)
        out_meshes05.append(mesh)
    m_atrMyoc, m_ventMyoc, m_atrCJ, m_ventCJ, disc = out_meshes05
    
    settings.legendSize = .20
    vp = Plotter(N=6, axes=13)
    vp.show(m_myoc.alpha(0.01), scale_cube, disc, at=0)
    vp.show(m_atrMyoc.alpha(0.01), scale_cube, at=1)
    vp.show(m_ventMyoc.alpha(0.01), scale_cube, at=2)
    vp.show(m_cj.alpha(0.01), disc, scale_cube, at=3)
    vp.show(m_atrCJ.alpha(0.01), scale_cube, at=4)
    vp.show(m_ventCJ.alpha(0.01), scale_cube, at=5, zoom=1.6, azimuth = azimuth, elevation = 0, interactive=True)


    m_ventMyoc.alpha(0.01).color('darkturquoise')
    

    #%% GET CHAMBERS ORIENTATION AND ELLIPSOIDS
    num_pt = dict_pts['numPt_CLChamberCut']
    sph_orient, lines_orient, dict_pts, dict_kspl, df_res = fcMeshes.getChambersOrientation(filename = filename, file_num = file_num, 
                                                                    num_pt = num_pt, kspl_CL2use = kspl_CL, distFromCl = 50,
                                                                    myoc_meshes = [m_atrMyoc, m_ventMyoc], linLine = linLines,
                                                                    dict_pts = dict_pts, dict_kspl = dict_kspl, df_res = df_res,
                                                                    scale_cube = scale_cube)
   
    orient_atr, orient_vent, orient_atrX, orient_ventX, linLineX = lines_orient
    df_res, dict_shapes = fcMeshes.getChambersEllipsoid(filename = filename, df_res = df_res, file_num = file_num, 
                                  lines_orient = lines_orient, meshes =(m_atrMyoc, m_ventMyoc), dict_shapes = dict_shapes)

    
    
    #%% CREATE PLANES AND RIBBON TO DIVIDE HEART INTO SECTIONS (Dors/Vent, Left/Right)
    # Get centreline ribbon
    cl_ribbon, kspl_ext, dict_kspl, dict_shapes, dict_planes = fcMeshes.createCLRibbon(filename = filename, file_num = file_num, 
                                                                                       df_res = df_res, kspl_CL2use = kspl_CL, 
                                                                                       linLine = linLines, mesh = m_myoc, 
                                                                                       dict_kspl = dict_kspl, dict_shapes = dict_shapes, 
                                                                                       dict_planes = dict_planes, scale_cube = scale_cube)

    # Divide meshes using ribbon (Left/Right) 
    # First divide the cardiac jelly and save the volumes of it's left and right sides
    [m_cjLnR], names_LnR = fcMeshes.divideMeshesLnR_new(filename = filename, meshes = [m_cj], cl_ribbon = cl_ribbon, 
                                         file_num = file_num, df_res = df_res, 
                                         scale_cube = scale_cube, colors = ['salmon', 'brown'])
    
    df_res = fcMeshes.addLayersVolume2df (df_res = df_res, file_num = file_num, meshes = [m_cj]+m_cjLnR, 
                                          names = names_LnR)
    
    
    settings.legendSize = .20
    vp = Plotter(N=1, axes=8)
    vp.show(m_cj.alpha(0.01), cl_ribbon, scale_cube, at=0, zoom=1, azimuth = 0, elevation = 0, interactive=True)
    
    settings.legendSize = .20
    vp = Plotter(N=1, axes=8)
    vp.show(m_cj, cl_ribbon, scale_cube, at=0, zoom=1, azimuth = 90, elevation = 0, interactive=True)
    
    # m_cjLnR[0].color('salmon')
    # m_cjLnR[1].color('brown')
    
    settings.legendSize = .20
    vp = Plotter(N=2, axes=13)
    vp.show(m_cjLnR[0], scale_cube, at=0)
    vp.show(m_cjLnR[1], scale_cube, at=1, zoom=1.6, azimuth = 0, elevation = 0, interactive=True)
    settings.legendSize = .20
    vp = Plotter(N=2, axes=13)
    vp.show(m_cjLnR[0], scale_cube, at=0)
    vp.show(m_cjLnR[1],scale_cube, at=1, zoom=1.6, azimuth = 90, elevation = 0, interactive=True)
    
    #%% test atr?vent cj into left and right 
    m_atrVent = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc_atr', 'myoc_vent','cj_atr','cj_vent'],
                                    extension = 'vtk', dir_stl = directories[2],
                                    alpha = [1,1,1,1], dict_colour = dict_colour)
    cyl_chamber = dict_shapes['cyl2CutChambers_final']
    disc = Cylinder(pos = cyl_chamber['cyl_centre'],r = cyl_chamber['radius_max'], height = 2*0.225, 
                    axis = cyl_chamber['cyl_axis'], c = 'purple', cap = True, res = 300)
    
    in_meshes05 = m_atrVent+[disc]
    out_meshes05 = []
    for mesh in in_meshes05:
        mesh = mesh.rotateX(rotAngle)
        out_meshes05.append(mesh)
    m_atrMyoc, m_ventMyoc, m_atrCJ, m_ventCJ, disc = out_meshes05
    
    # Divide meshes using ribbon (Left/Right) 
    # First divide the cardiac jelly and save the volumes of it's left and right sides
    [m_cjLnR], names_LnR = fcMeshes.divideMeshesLnR_new(filename = filename, meshes = [m_atrCJ], cl_ribbon = cl_ribbon, 
                                         file_num = file_num, df_res = df_res, 
                                         scale_cube = scale_cube, colors = ['salmon', 'brown'])
    
    
    
      
    #%%
    cl_ribbonV, kspl_ext, _, _, _ = fcMeshes.createCLRibbon(filename = filename, file_num = file_num, 
                                                df_res = df_res, kspl_CL2use = kspl_CL, linLine = linLines,
                                                mesh = m_myoc, dict_kspl = dict_kspl, dict_shapes = dict_shapes, 
                                                dict_planes = dict_planes, clRib_type = 'extV', plotshow = True)
    
    [m_cjThLnR], _ = fcMeshes.divideMeshesLnR_new(filename = filename, meshes = [m_cjTh], cl_ribbon = cl_ribbon, 
                                               file_num = file_num, df_res = df_res)
    
    #%% CLASSIFY THICKNESS AND BALLOONING POINTS
    [[m_cjTh0, m_myocIntBall0], [cj_thickness, myoc_intBall]] = fcMeshes.openThicknessMeshes(filename = filename, 
                                                                                          meshes_names = ['cj_thickness','myoc_intBall'], 
                                                                                          extension = 'vtk', dir_stl = directories[2], dir_txtNnpy = directories[1]); 
    m_cjTh0.pointColors(cj_thickness, cmap="turbo", vmin=0, vmax=20)
    m_cjTh0.addScalarBar()
    m_cjTh0.alpha(1)
    m_cjTh0.mapper().SetScalarRange(0,20)
    
    [m_atrExtCJ] = fcMeshes.openMeshes(filename = filename, meshes_names = ['cjExt_atr'],
                                                                  extension = 'vtk', dir_stl = directories[2],
                                                                  alpha = [1], dict_colour = dict_colour)
    

    kSplinesCuts = fcMeshes.createKSpls(dict_kspl, kspl_list = ['ksplCut4CL_inflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Ext.Endo(Cut)', 
                                                                'ksplCut4CL_outflow-Int.Myoc(Cut)', 'ksplCut4CL_inflow-Int.Myoc(Cut)'])
    sphCuts = fcMeshes.createSpheres(dict_pts, pts_list = ['sph_Cut4CL_inflow-Int.Myoc(Cut)', 'sph_Cut4CL_inflow-Ext.Endo(Cut)', 
                                                           'sph_Cut4CL_outflow-Int.Myoc(Cut)', 'sph_Cut4CL_outflow-Ext.Endo(Cut)'])
    kspl_CL, linLines, ksplSph_o, dict_kspl = fcMeshes.createCLs(filename = filename, dict_cl = dicts_cl, dict_pts = dict_pts, 
                                                                 dict_kspl = dict_kspl, dict_planes = dict_planes, 
                                                                 colors = ['deepskyblue', 'tomato'], myoc = m_myoc, 
                                                                 dir_stl = directories[3])
    sph_CL, sph_CL_colour, dict_shapes = fcMeshes.createColouredCL(dict_cl= dicts_cl, dict_shapes = dict_shapes)
    
    # Get centreline ribbon
    cl_ribbon, kspl_ext, dict_kspl, dict_shapes, dict_planes = fcMeshes.createCLRibbon(filename = filename, file_num = file_num, 
                                                                                       df_res = df_res, kspl_CL2use = kspl_CL[0], 
                                                                                       linLine = linLines[0], mesh = m_myoc, 
                                                                                       dict_kspl = dict_kspl, dict_shapes = dict_shapes, 
                                                                                       dict_planes = dict_planes)
    
    [m_cjThLnR0], _ = fcMeshes.divideMeshesLnR_new(filename = filename, meshes = [m_cjTh0], cl_ribbon = cl_ribbon, 
                                               file_num = file_num, df_res = df_res)
    
    
    df_cjThNmyocIntBall = fcMeshes.classifyHeartPts(filename = filename, df_res = df_res, file_num = file_num,
                                            dict_planes = dict_planes,
                                            m_whole = m_cjTh0, m_left = m_cjThLnR0[0], m_atr = m_atrExtCJ, 
                                            data = [cj_thickness, myoc_intBall], 
                                            names_data = ['cj_thickness', 'myoc_intBall'], plot_show = True)
    
    #%%
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
    
    