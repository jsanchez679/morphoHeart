# -*- coding: utf-8 -*-
"""
morphoHeart - C. CUT ALL HEART LAYERS (ATRIUM/VENTRICLE, DORSAL/VENTRAL, LEFT/RIGHT) AND QUANTIFY THEM
Welcome to the fourth code of morphoHeart!
If you are running this code you must have already created all the volume reconstructions of all the tissue layers of 
the heart you are processing and extracted its centreline and are looking forward to get the cardiac jelly and ballooning 
heatmap of this heart! The objective of this code is to divide all the reconstructed heart tissue layers into sections
(atrium/ventricle, dorsal/ventral, left/right ) and get all the heatmaps you need for your future analyses. 
At the end of this code you will have created the chamber meshes for each heart tissue layer, heatmap meshes for the 
external myocardium, endocardium and cardiac jelly, representing its thickness, and another heatmap mesh(es) for the 
internal (and external) myocardium representing its balooning. If you wish, at the end there is also an option to save 
videos of all these meshes rotating, so that you can include them in your presentations! 

Happy sectioning and heatmapping! 

@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
from time import perf_counter
from datetime import datetime
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
    # init = True
    print("Current working directory: {0}".format(os.getcwd()))

    return root_path, init

root_path, init = setWorkingDir(os.getcwd(),init)

c="k"; font= 'VTK';
save = True; plot = True; plotshow = False
azimuth = 0

#%% Start C_CutAndMeasure
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    # from morphoHeart_modules import morphoHeart_funcContours as fcCont
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
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
    #   This section will load all the heart tissue layer meshes already created in previous scripts, as well as the
    #   centreline(s) dictionary(ies) and create a plot with all of the meshes.
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
        [dict_planes, dict_pts, dict_kspl, dict_colour,_] = dict_obj

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

    #%% Create ksplines, points, lines and centrelines
    #   This section will create spline(s) of the centreline(s) using the information from the loaded dictionary(ies) and
    #   return first a 3D interactive plot with the myocardium, endocardium and centreline(s) and second a new plot
    #   with two different visualisations of the maximum inscribed spheres with which the centreline(s) were calculated.
    #   Note: The resulting centreline(s) will comprise 300 points (important for next step)
    #   ================================================================================================================

    kSplinesCuts = fcMeshes.createKSpls(dict_kspl, kspl_list = ['ksplCut4CL_inflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Int.Myoc(Cut)', 'ksplCut4CL_inflow-Int.Myoc(Cut)'])
    sphCuts = fcMeshes.createSpheres(dict_pts, pts_list = ['sph_Cut4CL_inflow-Int.Myoc(Cut)', 'sph_Cut4CL_inflow-Ext.Endo(Cut)', 'sph_Cut4CL_outflow-Int.Myoc(Cut)', 'sph_Cut4CL_outflow-Ext.Endo(Cut)'])
    kspl_CL, linLines, sph_CL, sph_CL_colour, dict_shapes, dict_kspl = fcMeshes.createCLs(dict_cl = dicts[1:], dict_pts = dict_pts,
                                                                                          dict_kspl = dict_kspl, dict_shapes = dict_shapes,
                                                                                          colors = ['deepskyblue', 'tomato'])
    # Save measurements and shapes
    df_res = fcMeshes.addLinearMeas2df(df_res = df_res, file_num = file_num, lines = linLines, kspl_CL = kspl_CL)

    if plot:
        settings.legendSize = .20
        vp = Plotter(N=1, axes=10)
        vp.show(txt, kSplinesCuts, sphCuts, kspl_CL, linLines, m_myoc.alpha(0.01), m_endo.alpha(0.01), at=0, azimuth = azimuth, interactive=True)

        settings.legendSize = .20
        if len(dicts[1:]) == 2:
            vp = Plotter(N=4, axes=13)
            vp.show(m_myoc.alpha(0.01), sph_CL_colour[0], txt, at=0)
            vp.show(m_endo.alpha(0.01), sph_CL_colour[1], at=1)
            vp.show(m_myoc.alpha(0.01), sph_CL[0], at=2)
            vp.show(m_endo.alpha(0.01), sph_CL[1], at=3, azimuth = azimuth, interactive=True)
        else:
            vp = Plotter(N=2, axes=13)
            vp.show(m_myoc.alpha(0.01), sph_CL_colour[0], txt, at=0)
            vp.show(m_myoc.alpha(0.01), sph_CL[0], at=1, azimuth = azimuth, interactive=True)

    #%% GET LAYERS THICKNESS HEATMAP
    #   This section of the code will extract the thickness heatmap of each heart tissue layer calculating the shortest
    #   distance from each of the points that make up the external mesh (external myocardium, endocardium or cardiac
    #   jelly) to the corresponding internal mesh. At the end, a 3D plot showing the three different thickness heatmaps
    #   will pop-up and the resulting meshes will be saved.
    #   ================================================================================================================

    # Get layers thickness heatmap
    print('- Extracting thickness heatmaps for all heart layers...')
    # Cardiac Jelly
    cj_thickness, m_cjTh, _ = fcMeshes.getDistance2Mesh(filename = filename, m_int = m_cjIn, m_ext = m_cjOut,
                                                        title = 'Cardiac Jelly Thickness', plotshow = False)
    # Myocardium
    myoc_thickness, m_myocTh, _ = fcMeshes.getDistance2Mesh(filename = filename, m_int = m_myocInt, m_ext = m_myocExt,
                                                        title = 'Myoc.Thickness', plotshow = False)
    # Endocardium
    endo_thickness, m_endoTh, _ = fcMeshes.getDistance2Mesh(filename = filename, m_int = m_endoInt, m_ext = m_endoExt,
                                                        title = 'Endo.Thickness', plotshow = False)

    if plot:
        vp = Plotter(N=3, axes=10)
        vp.show(m_cjTh.alpha(1), txt, at=0)
        vp.show(m_myocTh.alpha(1), at=1)
        vp.show(m_endoTh.alpha(1), scale_cube, at=2, zoom=2, azimuth = azimuth, elevation = 0, interactive=True)

    if save:
        fcMeshes.saveThickness(filename = filename, arrays2save = [cj_thickness, myoc_thickness, endo_thickness],
                                names = ['cj_thickness','myoc_thickness','endo_thickness'], dir2save = directories[1])
        dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [m_cjTh, m_myocTh, m_endoTh],
                                names = ['cj_thickness','myoc_thickness','endo_thickness'],
                                dict_colour = dict_colour, dir_stl = directories[2], extension = 'vtk')

    #%% GET BALLOONING HEATMAPS
    #   This section of the code will extract the balloning heatmaps of the internal and external myocardial meshes
    #   calculating the shortest distance from each of the points that make up these meshes to the centreline.
    #   At the end, a 3D plot showing the two different ballooning heatmaps will pop-up and the resulting meshes
    #   will be saved. Note: Calculation of the ballooning heatmap for each mesh takes about 10-15 min. Be patient! :)
    #   ================================================================================================================
    
    # Get Ballooning
    print('- Extracting ballooning information for Internal (and External) Myocardium... \n\tNOTE: it takes about 10-15 to process each mesh... just be patient :) ')
    print("  > Start time: \t", str(datetime.now())[11:-7])
    sph_ballonning = fcMeshes.sphInSpline(kspl_CL = kspl_CL[0], name = 'sphs_ball', every = 0.6)
    # Add shapes to dict
    dict_shapes = fcMeshes.addShapes2Dict(shapes = [sph_ballonning], dict_shapes = dict_shapes, radius = [[]])

    if plot:
        text = filename+"\n\n >> Int. and Ext. Myocardium and CL"; txt = Text2D(text, c="k", font= 'CallingCode')
        vp = Plotter(N=3, axes = 13)
        vp.show(m_myocInt.alpha(0.01).color(dict_colour['myoc_int']['colour']), sph_ballonning, txt,  at=0)
        vp.show(m_myocExt.alpha(0.01).color('teal'), sph_ballonning, at=1)
        vp.show(m_myocInt, m_myocExt, sph_ballonning, scale_cube, at=2, zoom = 2, elevation = 0, interactive = True)

    myoc_intBall, m_myocIntBall, [myoIntmin, myoIntmax] = fcMeshes.getDistance2Mesh(filename = filename, m_int = sph_ballonning, m_ext = m_myocInt,
                                                                                        title = 'Internal Myocardium Ballooning', alpha = 1)
    if save:
        fcMeshes.saveThickness(filename = filename, arrays2save = [myoc_intBall], names = ['myoc_intBall'], dir2save = directories[1])
        fcMeshes.saveMesh(filename = filename, mesh = m_myocIntBall, mesh_name = 'myoc_intBall', dir_stl = directories[2], extension = 'vtk')

    # myoc_extBall, m_myocExtBall, [myoExtmin, myoExtmax] = fcMeshes.getDistance2Mesh(filename = filename, m_int = sph_ballonning, m_ext = m_myocExt,
    #                                                                                     title = 'External Myocardium Ballooning', alpha = 1)
    # if save:
    #     fcMeshes.saveThickness(filename = filename, arrays2save = [myoc_extBall], names = ['myoc_extBall'], dir2save = directories[1])
    #     fcMeshes.saveMesh(filename = filename, mesh = m_myocExtBall, mesh_name = 'myoc_extBall', dir_stl = directories[2], extension = 'vtk')

    if plot:
        settings.legendSize = .2
        if 'myoc_extBall' not in locals():
            vp = Plotter(N=4, axes = 13)
            vp.show(m_myocInt.alpha(0.01), sph_ballonning, at=0)
            vp.show(m_myocInt.alpha(0.01), m_myocExt.alpha(0.01), sph_ballonning, at=1)
            vp.show(m_myocIntBall.alpha(0.01), sph_ballonning, at=2)
            vp.show(m_cjTh, scale_cube, at=3, azimuth = azimuth, elevation = 0, zoom = 2, interactive = True)
        else: 
            vp = Plotter(N=6, axes = 13)
            vp.show(m_myocInt.alpha(0.01), sph_ballonning, at=0)
            vp.show(m_myocExt.alpha(0.01), sph_ballonning, at=1)
            vp.show(m_myocInt.alpha(0.01), m_myocExt.alpha(0.01), sph_ballonning, at=2)
            vp.show(m_myocIntBall.alpha(0.01), sph_ballonning, at=3)
            vp.show(m_myocExtBall.alpha(0.01), sph_ballonning, at=4)
            vp.show(m_cjTh, scale_cube, at=5, azimuth = azimuth, elevation = 0, zoom = 2, interactive = True)

    #%% CUT HEART TISSUE LAYERS INTO CHAMBERS
    #   This section will allow the user to define a plane through the atrio-ventricular canal to cut all the heart
    #   tissue layers (myocardium, endocardium and cardiac jelly) into atrial and ventricular regions.
    #   Initially, the points that make up the centreline will be used to create a series of spheres (centered in the
    #   centreline) every 10th point starting from the outflow to inflow tract of the heart. A 3D interactive plot will
    #   pop-up showing the myocardium, the centreline and the series of spheres, to allow the user to roughly define
    #   the centreline point number through which the atrio-ventricular canal is located. Using this information a
    #   very think disk will be initialised and a new plot will pop-up in which the user will be able to rotate and 
    #   translate this disk to the position and orientation that will best cut the heart into atrial and ventricular 
    #   regions. When the user is happy with the disk position and orientation, he/she should close the window. A 
    #   new pop-up window will appear showing the selected plane and heart tissue layers. When the window is closed the 
    #   user will be asked if he/she is happy with the defined disk. If so, the radius of the disk will be 
    #   re-calculated based on the myocardial mesh points the disk cuts. The user will then see a new plot of the 
    #   myocardium and the created ring to double check if the radius of the ring is sufficient to split the heart without 
    #   affecting any of its chambers. If the user is happy, a mask of the created ring will be used to divide the masks
    #   of each tissue layer. Once the myocardial mesh has been divided, a 3D interactive plot will appear showing the
    #   resulting meshes. This plot will continue to get filled with the resulting meshes of the two other tissue layers 
    #   (endocardium and cardiac jelly). This process might take between 10-20mins. 
    #   IMPORTANT NOTE: Do not interact with the plot before the atrial and ventricular cardiac jelly meshes have appeared 
    #   in the third plot as this can cause a glitch and might interrupt the whole process! Once these plots appear you 
    #   can move the meshes around and when you are done close the window. The code will continue running cutting the
    #   mesh of the external myocardium, internal endocardium and external cardiac jelly.
    #   When all the volumetric reconstructions of the chambers of all the tissue layers have been generated, volume
    #   measurements will be taken and added to the measurements dataframe.
    #   ================================================================================================================
    
    # Divide heart layers into chambers and save data
    cyl_Chambers, num_pt, dict_shapes, dict_pts = fcMeshes.getRing2CutChambers(filename = filename, kspl_CL = kspl_CL[0],
                                                            mesh2cut = m_myoc, dir_stl = directories[2], dir_txtNnpy = directories[1],
                                                            dict_pts = dict_pts, dict_shapes = dict_shapes)
    
    m_atr, m_vent, dict_shapes = fcMeshes.getChamberMeshes(filename = filename, m_myoc = m_myoc, 
                                    end_name = ['ch0_cut', 'ch1_cut', 'cj', 'ch0_cut_ext', 'ch1_cut_int', 'ch0_cut_int'],
                                    names2cut = ['Myoc','Endo', 'CJ', 'Ext.Myoc', 'Int.Endo', 'Ext.CJ'],
                                    kspl_CL = kspl_CL[0], num_pt = num_pt, dir_txtNnpy = directories[1],
                                    dict_shapes = dict_shapes, resolution = res, plotshow = plotshow)
    
    m_atrMyoc, m_atrEndo, m_atrCJ, m_atrExtMyo, m_atrIntEnd, m_atrExtCJ = m_atr
    m_ventMyoc, m_ventEndo, m_ventCJ, m_ventExtMyo, m_ventIntEnd, m_ventExtCJ = m_vent

    if plot:
        settings.legendSize = .20
        vp = Plotter(N=6, axes = 10)
        vp.show(m_atr[0], m_vent[0], txt, at = 0)
        vp.show(m_atr[1], m_vent[1], at = 1)
        vp.show(m_atr[2], m_vent[2], at = 2)
        vp.show(m_atr[3], m_vent[3], at = 3)
        vp.show(m_atr[4], m_vent[4], at = 4)
        vp.show(m_atr[5], m_vent[5], scale_cube, at = 5, zoom = 2.2, azimuth = azimuth, interactive=True)

    # Add chamber volume info for each layer and save df
    df_res = fcMeshes.addLayersVolume2df (df_res = df_res, file_num = file_num,
                                          meshes = [m_myoc, m_atrMyoc, m_ventMyoc, m_atrExtMyo, m_ventExtMyo,
                                                    m_endo, m_atrEndo, m_ventEndo, m_atrIntEnd, m_ventIntEnd,
                                                    m_cj, m_atrCJ, m_ventCJ])

    if save:
        # Save measurements dataframe
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
        # Update dict_colour
        dict_colour = fcMeshes.saveMeshes(filename = filename, 
                                          meshes = [m_atrMyoc, m_atrEndo, m_atrCJ, m_atrExtMyo, m_atrIntEnd, m_atrExtCJ,
                                                    m_ventMyoc, m_ventEndo, m_ventCJ, m_ventExtMyo, m_ventIntEnd, m_ventExtCJ],
                                          names = ['myoc_atr', 'endo_atr', 'cj_atr', 'myocExt_atr', 'endInt_atr', 'cjExt_atr'
                                                   'myoc_vent', 'endo_vent', 'cj_vent', 'myocExt_vent', 'endInt_vent', 'cjExt_vent'],
                                          dict_colour = dict_colour, dir_stl = directories[2], extension = 'vtk')

    #%% GET CHAMBERS AND HEART ORIENTATION
    #   This section of the code will measure the chambers and heart orientations and save them in the dataframe with
    #   the rest of the measurements already made. At the end, a 3D plot showing the chambers orientation will pop-up.
    #   ================================================================================================================

    sph_orient, lines_orient, dict_pts, dict_kspl, df_res = fcMeshes.getChambersOrientation(filename = filename, file_num = file_num, 
                                                                    num_pt = num_pt, kspl_CL2use = kspl_CL[0],distFromCl = 50,
                                                                    myoc_meshes = [m_myoc, m_atrMyoc, m_ventMyoc], linLine = linLines[0],
                                                                    dict_pts = dict_pts, dict_kspl = dict_kspl, df_res = df_res)
    elevation = df_res.loc[file_num,'ang_Heart']
    orient_atr, orient_vent, orient_atrX, orient_ventX, linLineX = lines_orient

    if save:
        # Save measurements dataframe
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)

    #%% CREATE PLANES AND RIBBON TO DIVIDE HEART INTO SECTIONS (Dors/Vent, Left/Right)
    #   This section of the code will create planes to divide each of the chambers into dorso-ventral regions and will
    #   extend dorso-ventrally the centreline to create a ribbon that divides the heart into left and right. 3D
    #   interactive plots will pop-up showing first the dorso-ventral planes created and second the centreline ribbon.
    #   Once the ribbon has been created it will be used to cut the external cardiac jelly mesh and the internal
    #   myocardium mesh. Resulting meshes will be shown and the user will be asked to confirm the resulting left-right 
    #   classification. 
    #   ================================================================================================================

    # Create Coronal Planes for each chamber
    [pl_DnV_Atr, pl_DnV_Vent], dict_planes = fcMeshes.createDVPlanes(filename = filename, sph_orient = sph_orient, mesh = m_myoc,
                                                  kspl_CL = kspl_CL[0], orient_lines = [orient_atr, orient_vent], dict_planes = dict_planes)
    # Get centreline ribbon
    cl_ribbon, kspl_ext, dict_kspl, dict_shapes, dict_planes = fcMeshes.createCLRibbon(filename = filename, kspl_CL2use = kspl_CL[0], linLine = linLines[0],
                                                                 mesh = m_myoc, dict_kspl = dict_kspl, dict_shapes = dict_shapes, dict_planes = dict_planes)


    # TO DELETEEE!!! Get myocIntBall data
    [[m_cjTh, m_myocIntBall], [cj_thickness, myoc_intBall]] = fcMeshes.openThicknessMeshes(filename = filename, 
                                          meshes_names = ['cj_thickness','myoc_intBall'], extension = 'vtk',
                                          dir_stl = directories[2], dir_txtNnpy = directories[1], print_txt = False)
    
    # Divide meshes using ribbon (Left/Right)
    # [m_cjThLnR, m_myocExtLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_cjTh, m_myocExtBall], cl_ribbon = cl_ribbon)
    [m_cjThLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_cjTh], cl_ribbon = cl_ribbon)
    # [m_myocExtLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_myocExtBall], cl_ribbon = cl_ribbon)

    #%% CLASSIFY THICKNESS AND BALLOONING POINTS
    #   Now that we have defined the planes and ribbon that section the heart into atrial-ventricular, dorsal-ventral,
    #   left-right sections of the heart, this section of the code will classify each of the points that make up the
    #   external cardiac jelly and the internal myocardium (meshes in which the cardiac jelly thickness heatmap and the
    #   internal myocardium ballooning heatmap respectively had been mapped), resulting in a dataframe containing all
    #   this information. Once the classification algorithm has run, a 3D interactive plot will pop-up showing the
    #   resulting classification for a subsample of the points color-coded to each region. At the end, the resulting
    #   dataframe will be saved for future analysis.
    #   ================================================================================================================

    df_cjThNmyocIntBall = fcMeshes.classifyHeartPts(filename = filename, mesh = m_cjTh, dict_planes = dict_planes,
                                            pts_whole = m_cjTh.points(), pts_left = m_cjThLnR[0].points(), 
                                            pts_atr = m_atrExtCJ.points(), data = [cj_thickness, myoc_intBall], 
                                            names_data = ['cj_thickness', 'myoc_intBall'], plot_show = True)
    # df_cjThNmyocIntBall = fcMeshes.classifyHeartPts(filename = filename, mesh = m_cjTh, dict_planes = dict_planes,
    #                                         pts_whole = m_cjTh.points(), pts_left = m_cjThLnR[0].points(), 
    #                                         pts_atr = m_atrExtCJ.points(), data = [cj_thickness], 
    #                                         names_data = ['cj_thickness'], plot_show = True)

    if save:
        fcBasics.saveDF(filename = filename, df2save = df_cjThNmyocIntBall, df_name = 'df_cjThNmyocIntBall',
                        dir2save = dir_results)

    #%% SAVE ALL
    #   This section allows the user to save all the meshes, objects and dataframes that have been created. It will also
    #   ask the user if he/she wants to save videos of all the different created meshes (including those with heatmaps)
    #   ================================================================================================================

    if save:
        # Append all dicts to one object dict
        dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
                                             names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
        # Save filled dataframe with measured data
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = os.path.join(dir_data2Analyse, 'R_All'))

        #Save videos
        saveVideos = fcBasics.ask4input('Do you want to save videos of all the resulting meshes? [0]:no/[1]:yes, please!: ', bool)
        if saveVideos:
            if 'myoc_extBall' not in locals():
                # fcMeshes.saveMultVideos(filename = filename,
                #     info = ['myoc','endo','cj', 'cj_thickness','myoc_thickness','endo_thickness','myoc_intBall'],
                #     meshes4video = [m_myoc,m_endo,m_cj,m_cjTh,m_myocTh,m_endoTh, m_myocIntBall.alpha(1)],
                #     rotAngle= df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow = False)
                fcMeshes.saveMultVideos(filename = filename,
                    info = ['endo_thickness','myoc_intBall'],
                    meshes4video = [m_endoTh, m_myocIntBall.alpha(1)],
                    rotAngle= df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow = False)
            else: 
                fcMeshes.saveMultVideos(filename = filename,
                    info = ['myoc','endo','cj', 'cj_thickness','myoc_thickness','endo_thickness','myoc_intBall','myoc_extBall'],
                    meshes4video = [m_myoc,m_endo,m_cj,m_cjTh,m_myocTh,m_endoTh, m_myocIntBall.alpha(1),m_myocExtBall.alpha(1)],
                    rotAngle= df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow = False)
                

        toc = perf_counter()
        fcBasics.printTime(tic, toc, 'Cut and Measure')

    df_cjThNmyocIntBall.sample(5)

    #%% CODE TO IGNORE!! save more information of cl? cl_ribbon? cl_ext? kspl_vsurf? Chech save shapes function
    #Divide cl for atrium and ventricle, into 300 pts each, use them to create heatmaps of the same length each? 
    #Add myoc_intBall data too to plot as well
    
    # -----------------------------------------------------------------------------------------------------------------
    df_cjThNmyocIntBall = fcBasics.loadDF(filename = filename, file = 'df_cjThNmyocIntBall', dir_results = dir_results)
    df_AtrVent = np.asarray(df_cjThNmyocIntBall['AtrVent'])
    cjTh = np.asarray(df_cjThNmyocIntBall['cj_thickness'])
    
    # Unlooping the heart
    unlooped = fcMeshes.unloopHeart(filename = filename, mesh = m_cjTh, kspl_CL = kspl_CL[0], kspl_ext = kspl_ext,
                                    no_planes = 250, pl_CLRibbon =  dict_planes['pl_Parallel2LinLine'], 
                                    param = cj_thickness, df_AtrVent = df_AtrVent, num_pt = dict_pts['numPt_CLChamberCut'],
                                    plotshow = False, tol=0.05)
    unlooped[100:130, :]
    
    vp = Plotter(N=1)
    vp.show(m_cjTh.alpha(1), at = 0, interactive = True)
    
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    df_unlooped = pd.DataFrame(unlooped, columns=['x','y','z','taken','z_plane','theta','radius','cj_thickness'])
    df_unlooped.sample(10)
    df_unlooped = df_unlooped[df_unlooped['taken']==1]
    
    fig, ax = plt.subplots(figsize=(16, 10))
    heatmap = pd.pivot_table(df_unlooped, values='cj_thickness', columns = 'z_plane', index='theta', aggfunc=np.max)
    ax = sns.heatmap(heatmap, cmap='jet')#, xticklabels=20, yticklabels=550)
    x_pos = np.linspace(0,1,len(ax.get_xticks()))
    y_pos = np.linspace(180,-180,len(ax.get_yticks()), endpoint = True)
    xlabels = [format(x,'.3f') for x in x_pos]
    ylabels = [format(y,'.0f') for y in y_pos]
    plt.show()

init = True

#%%

# df_filt = pd.DataFrame(columns=['x','y','z','taken','z_plane','theta','radius','cj_thickness'])
# z_planes = df_unlooped.z_plane.unique().tolist()#.sort()
# z_planes = sorted(z_planes)
# v_theta = np.linspace(-180,180,180+1)
# #loop for each plane
# for num, z_pl in enumerate(z_planes[0:50]):
#     df_z = df_unlooped[df_unlooped['z_plane'] == z_pl]
#     # print('--PLANE: ',z_pl)
#     for i, theta in enumerate(v_theta[0:-1]):
#         # print('--Theta:', theta)

#         df_t = df_z[(df_z['theta'] >= theta) & (df_z['theta'] <= v_theta[i+1])]
#         x = df_t['x'].mean()
#         y = df_t['y'].mean()
#         z = df_t['z'].mean()
#         taken = df_t['taken'].mean()
#         z_plane = df_t['z_plane'].mean()*100
#         theta = df_t['theta'].mean()
#         radius = df_t['radius'].mean()
#         cj_thickness = df_t['cj_thickness'].max()
#         data = [[x,y,z,taken,z_plane,theta,radius,cj_thickness]]
#         df_r = pd.DataFrame(data, columns=['x','y','z','taken','z_plane','theta','radius','cj_thickness'])

#         df_filt = pd.concat([df_filt, df_r])

# #% REVISARRRR
# fig, ax = plt.subplots(figsize=(8, 5))
# heatmap = pd.pivot_table(df_filt, values='cj_thickness', columns = 'z_plane', index='theta')#, aggfunc=np.max)
# ax = sns.heatmap(heatmap, cmap='jet', xticklabels=20, yticklabels=550)

# x_pos = np.linspace(0,1,len(ax.get_xticks()))
# y_pos = np.linspace(180,-180,len(ax.get_yticks()), endpoint = True)
# xlabels = [format(x,'.3f') for x in x_pos]
# ylabels = [format(y,'.0f') for y in y_pos]

# plt.show()

    
    
    
#     xlabels = [format(x/10,'.1f') for x in ax.get_xticks()]
#     ylabels = [format(y/1800*360-180,'.0f') for y in ax.get_yticks()]
#     ax.set_xticklabels(xlabels)
#     ax.set_yticklabels(ylabels)

  
#     # plt.xlabel('Centreline position\n[Atrium >> Ventricle]', fontsize=10)
#     # plt.ylabel('Angle (\N{DEGREE SIGN})\n[Ventral >> Dorsal >> Ventral]', fontsize=10)
    
#     x_pos = np.linspace(0,1,len(ax.get_xticks()))
#     y_pos = np.linspace(180,-180,len(ax.get_yticks()), endpoint = True)
#     xlabels = [format(x,'.3f') for x in x_pos]
#     ylabels = [format(y,'.0f') for y in y_pos]
    
    
#     flights = sns.load_dataset("flights")

#     flights = flights.pivot("month", "year", "passengers")

#     ax = sns.heatmap(flights)


    # Cut layer with a certain number of planes
    # no_planes = 15
    # sectDict = cjf4.sliceLayerCL(mesh3, no_planes, cl_splineAdd, 'red')
    # cjf4.plotAllLayerSlices(mesh2cut, cl_splineAdd, no_planes, sectDict)

    # Plot a particular slice
    # cjf4.plotLayerSlice(mesh2cut, cl_splineAdd, 5, sectDict)

    # -----------------------------------------------------------------------------------------------------------------
    # getBalloons
    #sph_all, sph_whole = fcMeshes.getBalloonedHeart(myoc_int_npcl)

    # vp = Plotter(N=1, axes=13)
    # vp.show(m_myoc.alpha(0.01), sph_all[0].color('green').alpha(0.05), sph_all[5].alpha(0.05), at=0, interactive=True)

    # sph_whole = booleanOperation(sph_all[0], 'plus',sph_all[2])
    # sph_whole.color('coral')

    # vp = Plotter(N=1, axes=13)
    # vp.show(m_myoc.alpha(0.01), sph_whole, at=0, interactive=True)

    # CODE TO IGNORE!!
    # plane_Ch = Plane(pos=dict_planes['pl2CutMesh_Chamber']['pl_centre'],normal=dict_planes['pl2CutMesh_Chamber']['pl_normal'], sx=500)
    # # Cut cl with plane
    # ksplCL_cut = kspl_CL[0].clone().cutWithMesh(plane_Ch)
    # # Find point of centreline closer to last point of kspline cut
    # ksplCL_cutPt, num_pt = fcMeshes.findClosestPt(ksplCL_cut.points()[0], kspl_CL[0].points())
    # dict_pts['numPt_CLChamberCut'] = num_pt
    # m_names = ['myoc_atr','myoc_vent','endo_atr','endo_vent','cj_atr','cj_vent','cj_in']
    # m_all = fcMeshes.openMeshes(filename = filename, meshes_names = m_names, extension = 'vtk',
    #                             dir_stl = directories[2], alpha = [1]*len(m_names), dict_colour = dict_colour)
    # m_atrMyoc, m_ventMyoc, m_atrEndo, m_ventEndo, m_atrCJ, m_ventCJ, m_cjIn = m_all

    # -----------------------------------------------------------------------------------------------------------------
    #%TO DELETEEE
    # cl_ribbon, dict_kspl, dict_shapes, dict_planes = fcMeshes.createCLRibbon(filename = filename, kspl_CL2use = kspl_CL[0], linLine = linLines[0],
                                                                  # mesh = m_myoc, dict_kspl = dict_kspl, dict_shapes = dict_shapes, dict_planes = dict_planes)
    # mTh_names = ['cj_thickness','myoc_thickness','endo_thickness','myoc_intBall','myoc_extBall']
    # [m_thAll, colour_thAll] = fcMeshes.openThicknessMeshes(filename = filename, meshes_names = mTh_names, extension = 'vtk', dir_stl = directories[2], dir_txtNnpy = directories[1])
    # m_cjTh, m_myocTh, m_endoTh, m_myocIntBall, m_myocExtBall = m_thAll

    # # Dictionary and decode it
    # [dict_cjThNmyocIntBall] = fcBasics.loadDicts(filename = filename, dicts_name = ['cjTh_myocIntBall2class'], directories = [directories[0]])
    # [AnV, DnV_Atr, DnV_Vent, pts_left, pts_whole, meas_param] = fcMeshes.decodeDict (dict2classify = dict_cjThNmyocIntBall,
    #                                                                                  info = ['cj_thickness','myoc_intBall'])
    # # npy arrays
    # cj_thickness, myoc_intBall = meas_param
