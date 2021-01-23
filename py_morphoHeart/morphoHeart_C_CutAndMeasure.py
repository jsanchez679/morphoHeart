# -*- coding: utf-8 -*-
"""
morphoHeart - C. CUT ALL HEART LAYERS (ATRIUM/VENTRICLE, DORSAL/VENTRAL, LEFT/RIGHT) AND QUANTIFY THEM
@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
import psutil
from time import perf_counter
from vtkplotter import *
from vtkplotter import embedWindow
embedWindow(False)

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

c="k"; font= 'CallingCode'
save = False
plot = True
azimuth = 0

#%% Start C_CutAndMeasure
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcContours as fcCont
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
    tic = perf_counter()

    #%% Get main directories (check which ones are actually used)
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

    #%% Load data
    # Import dictionaries
    [dict_obj, myoc_int_npcl, endo_ext_npcl] = fcBasics.loadDicts(filename = filename, dicts_name = ['dict_obj','myoc_int_npcl','endo_ext_npcl'],
                                                                    directories = [directories[0], directories[3], directories[3]])
    dict_obj = fcMeshes.splitDicts(dict_obj)
    if len(dict_obj) == 4:
        [dict_planes, dict_pts, dict_kspl, dict_colour] = dict_obj
    else:
        [dict_planes, dict_pts, dict_kspl, dict_colour,_] = dict_obj

    # Import meshes
    [m_myoc, m_endo, m_cj, m_cjOut, m_cjIn] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc','endo','cj','cj_out','cj_in'],
                                                                  extension = 'vtk', dir_stl = directories[2],
                                                                  alpha = [0.5,0.5,0.5,1,1], dict_colour = dict_colour)
    [m_myocInt, m_myocExt, m_endoInt, m_endoExt] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc_int','myoc_ext','endo_int', 'endo_ext'],
                                                                  extension = 'vtk', dir_stl = directories[2],
                                                                  alpha = [1,1,1,1], dict_colour = dict_colour)

    scale_cube = Cube(pos=m_myoc.centerOfMass(), side=350, c='white', alpha=0.01)
    # Plot meshes
    if plot:
        text = str(filename); txt = Text2D(text, c=c, font=font)
        vp = Plotter(N=6, axes=7)
        vp.show(m_myoc, txt, at=0)
        vp.show(m_endo, at=1)
        vp.show(m_cj, at=2)
        vp.show(m_cjIn, at=3)
        vp.show(m_cjOut, at=4)
        vp.show(m_myoc, m_endo, m_cj, scale_cube, at=5, zoom = 2, azimuth = azimuth, interactive=True)

    #%% Create ksplines, points, lines and centrelines
    kSplinesCuts = fcMeshes.createKSpls(dict_kspl, kspl_list = ['ksplCut4CL_inflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Int.Myoc(Cut)', 'ksplCut4CL_inflow-Int.Myoc(Cut)'])
    sphCuts = fcMeshes.createSpheres(dict_pts, pts_list = ['sph_Cut4CL_inflow-Int.Myoc(Cut)', 'sph_Cut4CL_inflow-Ext.Endo(Cut)', 'sph_Cut4CL_outflow-Int.Myoc(Cut)', 'sph_Cut4CL_outflow-Ext.Endo(Cut)'])
    kspl_CL, linLines, sph_CL, sph_CL_colour, dict_shapes, dict_kspl = fcMeshes.createCLs(dict_cl = [myoc_int_npcl,endo_ext_npcl], dict_pts = dict_pts,
                                                                                          dict_kspl = dict_kspl, dict_shapes = dict_shapes,
                                                                                          colors = ['deepskyblue', 'tomato'])
    # Save measurements and shapes
    df_res = fcMeshes.addLinearMeas2df(df_res = df_res, file_num = file_num, lines = linLines, kspl_CL = kspl_CL)

    if plot:
        vp = Plotter(N=1, axes=10)
        vp.show(txt, kSplinesCuts, sphCuts, kspl_CL, linLines, m_myoc.alpha(0.01), m_endo.alpha(0.01), at=0, azimuth = azimuth, interactive=True)

        vp = Plotter(N=4, axes=7)
        vp.show(m_myoc.alpha(0.01), sph_CL_colour[0], txt, at=0)
        vp.show(m_endo.alpha(0.01), sph_CL_colour[1], at=1)
        vp.show(m_myoc.alpha(0.01), sph_CL[0], at=2)
        vp.show(m_endo.alpha(0.01), sph_CL[1], at=3, azimuth = azimuth, interactive=True)

    #%% Cut layers into chambers
    # CHANGE!!! - Ellipsoid instead?? - https://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib
    # Divide heart layers into chambers and save data
    pl_Chambers, dict_planes, dict_pts, num_pt = fcMeshes.getInfo2CutChambers(filename = filename, kspl_CL = kspl_CL[0],
                                                                              mesh2cut = m_myoc, dict_planes = dict_planes, dict_pts = dict_pts)

    m_atr, m_vent = fcMeshes.getChamberMeshes(filename = filename, dir_txtNnpy = directories[1],
                                                end_name = ['ch0_cut', 'ch1_cut', 'cj', 'ch0_cut_ext', 'ch1_cut_int'],
                                                names2cut = ['Myoc','Endo', 'CJ', 'Ext.Myoc', 'Int.Endo'],
                                                kspl_CL = kspl_CL[0], num_pt = num_pt, dict_planes = dict_planes, resolution = res)

    m_atrMyoc, m_atrEndo, m_atrCJ, m_atrExtMyo, m_atrIntEnd = m_atr
    m_ventMyoc, m_ventEndo, m_ventCJ, m_ventExtMyo, m_ventIntEnd = m_vent

    if plot:
        vp = Plotter(N=6, axes = 10)
        vp.show(m_atr[0:3], pl_Chambers.alpha(0.05), txt, at = 0)
        vp.show(m_vent[0:3], pl_Chambers.alpha(0.05), at = 3)
        vp.show(m_atr[3], m_vent[3], pl_Chambers.alpha(0.05), at = 1)
        vp.show(m_atr[4], m_vent[4], pl_Chambers.alpha(0.05), at = 4)
        vp.show(m_atr[3], m_vent[3], pl_Chambers.alpha(0.05), at = 2)
        vp.show(m_atr[4], m_vent[4], pl_Chambers.alpha(0.05), scale_cube, at = 5, zoom = 2.2, azimuth = azimuth, interactive=True)

    # Add chamber volume info for each layer and save df
    df_res = fcMeshes.addLayersVolume2df (df_res = df_res, file_num = file_num,
                                          meshes = [m_myoc, m_atrMyoc, m_ventMyoc, m_atrExtMyo, m_ventExtMyo,
                                                    m_endo, m_atrEndo, m_ventEndo, m_atrIntEnd, m_ventIntEnd,
                                                    m_cj, m_atrCJ, m_ventCJ])

    if save:
        # Save measurements data frame
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
        # Update dict_colour
        dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [m_atrMyoc, m_atrEndo, m_atrCJ, m_ventMyoc, m_ventEndo, m_ventCJ],
                            names = ['myoc_atr', 'endo_atr', 'cj_atr', 'myoc_vent', 'endo_vent', 'cj_vent'],
                            dict_colour = dict_colour, dir_stl = directories[2], extension = 'vtk')

    #%% Get chambers and heart orientation (sph_valve has the position of the cutting point in cl)
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

    sph_orient, lines_orient, dict_pts, dict_kspl, df_res = fcMeshes.getChambersOrientation(filename = filename, file_num = file_num, num_pt = num_pt, kspl_CL2use = kspl_CL[0],
                                                                                    myoc_meshes = [m_myoc, m_atrMyoc, m_ventMyoc], linLine = linLines[0],
                                                                                    dict_pts = dict_pts, dict_kspl = dict_kspl, df_res = df_res)
    elevation = df_res.loc[file_num,'ang_Heart']
    orient_atr, orient_vent, orient_atrX, orient_ventX, linLineX = lines_orient

    if save:
        # Save measurements data frame
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)

    #%% Get layers thickness heatmap
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
        vp.show(m_endoTh.alpha(1), scale_cube, at=2, zoom=2, azimuth = azimuth, elevation = elevation, interactive=True)

    if save:
        fcMeshes.saveThickness(filename = filename, arrays2save = [cj_thickness, myoc_thickness, endo_thickness],
                                names = ['cj_thickness','myoc_thickness','endo_thickness'], dir2save = directories[1])
        dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [m_cjTh, m_myocTh, m_endoTh],
                                names = ['cj_thickness','myoc_thickness','endo_thickness'],
                                dict_colour = dict_colour, dir_stl = directories[2], extension = 'vtk')

    #%% Get Ballooning
    sph_ballonning = fcMeshes.sphInSpline(kspl_CL = kspl_CL[0], name = 'sphs_ball', every = 0.6)
    # Add shapes to dict
    dict_shapes = fcMeshes.addShapes2Dict(shapes = [sph_ballonning], info = [''], dict_shapes = dict_shapes, radius = [[]])

    if plot:
        text = filename+"\n\n >> Int. and Ext. Myocardium and CL"; txt = Text2D(text, c="k", font= 'CallingCode')
        vp = Plotter(N=3, axes = 7)
        vp.show(m_myocInt.alpha(0.01).color(dict_colour['myoc_int']['colour']), sph_ballonning, txt,  at=0)
        vp.show(m_myocExt.alpha(0.01).color('teal'), sph_ballonning, at=1)
        vp.show(m_myocInt, m_myocExt, sph_ballonning, scale_cube, at=2, zoom = 2, elevation = elevation, interactive = True)

    print('- Extracting ballooning information for Internal and External Myocardium... \n\tit takes about 10-15 to process each mesh...')
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
        vp = Plotter(N=6, axes = 7)
        vp.show(m_myocInt.alpha(0.01), sph_ballonning, at=0)
        vp.show(m_myocExt.alpha(0.01), sph_ballonning, at=1)
        vp.show(m_myocInt.alpha(0.01), m_myocExt.alpha(0.01), sph_ballonning, at=2)
        vp.show(m_myocIntBall.alpha(0.01), sph_ballonning, at=3)
        vp.show(m_myocExtBall.alpha(0.01), sph_ballonning, at=4)
        vp.show(m_cjTh, scale_cube, at=5, azimuth = azimuth, elevation = elevation, zoom = 2, interactive = True)

    #%% Define planes and ribbons to divide heart into sections (Left/Right; Atr/Vent; Dors/Vent) (earlier stages???)
    # Create Coronal Planes for each chamber
    [pl_DnV_Atr, pl_DnV_Vent], dict_planes = fcMeshes.createDVPlanes(filename = filename, sph_orient = sph_orient, mesh = m_myoc,
                                                  kspl_CL = kspl_CL[0], orient_lines = [orient_atr, orient_vent], dict_planes = dict_planes)
    # Get centreline ribbon
    cl_ribbon, dict_kspl, dict_shapes, dict_planes = fcMeshes.createCLRibbon(filename = filename, kspl_CL2use = kspl_CL[0], linLine = linLines[0],
                                                                 mesh = m_myoc, dict_kspl = dict_kspl, dict_shapes = dict_shapes, dict_planes = dict_planes)

    #%% TO DELETEEE
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

    #%% Divide meshes using ribbon (Left/Right)
    [m_cjThLnR, m_myocExtLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_cjTh, m_myocExtBall], cl_ribbon = cl_ribbon)
    [m_cjThLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_cjTh], cl_ribbon = cl_ribbon)
    # [m_myocExtLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_myocExtBall], cl_ribbon = cl_ribbon)

    #%% Classify thickness and ballooning points
    df_cjThNmyocIntBall = fcMeshes.classifyHeartPts(filename = filename, mesh = m_cjTh, dict_planes = dict_planes,
                                            pts_whole = m_cjTh.points(), pts_left = m_cjThLnR[0].points(),
                                            data = [cj_thickness, myoc_intBall], names_data = ['cj_thickness', 'myoc_intBall'], plot_show = True)
    # df_cjThNmyocIntBall = fcMeshes.classifyHeartPts(filename = filename, mesh = m_cjTh, dict_planes = dict_planes,
    #                                         pts_whole = m_cjTh.points(), pts_left = m_cjThLnR[0].points(),
    #                                         data = [cj_thickness], names_data = ['cj_thickness'], plot_show = True)

    if save:
        fcBasics.saveDF(filename = filename, df2save = df_cjThNmyocIntBall, df_name = 'df_cjThNmyocIntBall', dir2save = dir_results)

    #%% Save all results
    if save:
        # Append all dicts to one object dict
        dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
                                             names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
        # Save filled dataframe with measured data
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
        fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = os.path.join(dir_data2Analyse, 'R_All'))

        #VIDEOS!!! # NEW VERSION CHECK!!!
        fcMeshes.saveMultVideos(filename = filename,
                    info = ['myoc','endo','cj', 'cj_thickness','myoc_thickness','endo_thickness','myoc_intBall','myoc_extBall'],
                    meshes4video = [m_myoc,m_endo,m_cj,m_cjTh,m_myocTh,m_endoTh, m_myocIntBall.alpha(1),m_myocExtBall.alpha(1)],
                    rotAngle= df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow = False)

        toc = perf_counter()
        fcBasics.printTime(tic, toc, 'Cut and Measure')

init = True

    #%% Unlooping the heart
    # unlooped = fcMeshes.unloopHeart(filename = filename, mesh = m_cjTh, kspl_CL2use = kspl_CL[0], cl_ribbon = cl_ribbon, no_planes = 100,
    #                                 pl_CLRibbon =  dict_planes['pl_Parallel2LinLine'], param = cj_thickness, tol=0.05)

    # pl_CLRibbon = dict_planes['pl_Parallel2LinLine']

    #%% #%% Cut layer with n number of planes perpendicular to spline and plot
    # mesh2cut = mesh_ch1
    # #Transform mesh from vtk structure to trimesh
    # print("Transforming mesh into trimesh object...")
    # mesh3= vtk2trimesh(mesh2cut)
    # cjf4.alert("jump",1)

    # Cut layer with a certain number of planes
    # no_planes = 15
    # sectDict = cjf4.sliceLayerCL(mesh3, no_planes, cl_splineAdd, 'red')
    # cjf4.plotAllLayerSlices(mesh2cut, cl_splineAdd, no_planes, sectDict)

    # Plot a particular slice
    # cjf4.plotLayerSlice(mesh2cut, cl_splineAdd, 5, sectDict)

    #%% getBalloons
    #sph_all, sph_whole = fcMeshes.getBalloonedHeart(myoc_int_npcl)

    # vp = Plotter(N=1, axes=7)
    # vp.show(m_myoc.alpha(0.01), sph_all[0].color('green').alpha(0.05), sph_all[5].alpha(0.05), at=0, interactive=True)

    # sph_whole = booleanOperation(sph_all[0], 'plus',sph_all[2])
    # sph_whole.color('coral')

    # vp = Plotter(N=1, axes=7)
    # vp.show(m_myoc.alpha(0.01), sph_whole, at=0, interactive=True)
