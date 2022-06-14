# -*- coding: utf-8 -*-
"""
TO GET IMAGES FOR THESIS
Created on Thu Mar 31 13:37:22 2022

@author: mdp18js
"""

#%% Importing python packages
import os
import numpy as np
from time import perf_counter
from itertools import count
# from datetime import datetime
# from vedo import *
from vedo import Plotter, Cube, settings, Text2D, Cylinder, Plane
from vedo import embedWindow
import pandas as pd

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
to_reprocess = []

#%% Start C_CutAndMeasure
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
    from morphoHeart_modules import morphoHeart_funcContours as fcCont
    tic = perf_counter()

#%% func - getCubeRibbon
    def getCubeRibbonMask(filename, file_num, df_res, cl_ribbon, res, dir_txtNnpy):
        
        spaw_analysis = False
        if 'spaw_ct' in df_res.loc[file_num,'spAnalysis']:
            spaw_analysis = True
        
        cl_ribbonR = cl_ribbon.clone().x(res[0])
        # cl_ribbonL = cl_ribbon.clone().x(-res[0])
        cl_ribbonF = cl_ribbon.clone().y(res[1])
        # cl_ribbonB = cl_ribbon.clone().y(-res[1])
        cl_ribbonT = cl_ribbon.clone().z(res[1])
        # cl_ribbonD = cl_ribbon.clone().z(-res[1])
        
        vp = Plotter(N=1, axes = 4)
        vp.show(cl_ribbon, cl_ribbonR,cl_ribbonF, cl_ribbonT, at=0, interactive = True)
        
        dir_txtNnpy = directories[1]
        #Load stack shape 
        [[xdim, ydim, zdim]] = fcBasics.loadNPY(filename, ['stackShape'], dir_txtNnpy, print_txt = False)
        # stack_shape = [xdim, ydim, zdim+2]
        
        rib_pts = cl_ribbon.points()
        rib_ptsR = cl_ribbonR.points()
        rib_ptsF = cl_ribbonF.points()
        rib_ptsT = cl_ribbonT.points()
        
        rib_points_rot = np.zeros_like(rib_pts)
        rib_points_rotR = np.zeros_like(rib_ptsR)
        rib_points_rotF = np.zeros_like(rib_ptsF)
        rib_points_rotT = np.zeros_like(rib_ptsT)
        
        # Rotate the points that make up the HR disc, to convert them to a stack
        resolution = res
        if 'CJ' not in filename: 
            axis = [0,0,1]
        else: 
            axis = [1,0,0]
            
        for i, pt in enumerate(rib_pts):
            rib_points_rot[i] = (np.dot(fcMeshes.rotation_matrix(axis = axis, theta = np.radians(90)),pt))
        for i, pt in enumerate(rib_ptsR):
            rib_points_rotR[i] = (np.dot(fcMeshes.rotation_matrix(axis = axis, theta = np.radians(90)),pt))
        for i, pt in enumerate(rib_ptsF):
            rib_points_rotF[i] = (np.dot(fcMeshes.rotation_matrix(axis = axis, theta = np.radians(90)),pt))
        for i, pt in enumerate(rib_ptsT):
            rib_points_rotT[i] = (np.dot(fcMeshes.rotation_matrix(axis = axis, theta = np.radians(90)),pt))
            
        rib_pix = np.transpose(np.asarray([rib_points_rot[:,i]//resolution[i] for i in range(len(resolution))]))
        rib_pix = rib_pix.astype(int)
        rib_pix = np.unique(rib_pix, axis =0)
        
        rib_pixR = np.transpose(np.asarray([rib_points_rotR[:,i]//resolution[i] for i in range(len(resolution))]))
        rib_pixR = rib_pixR.astype(int)
        rib_pixR = np.unique(rib_pixR, axis =0)
        
        rib_pixF = np.transpose(np.asarray([rib_points_rotF[:,i]//resolution[i] for i in range(len(resolution))]))
        rib_pixF = rib_pixF.astype(int)
        rib_pixF = np.unique(rib_pixF, axis =0)
        
        rib_pixT = np.transpose(np.asarray([rib_points_rotT[:,i]//resolution[i] for i in range(len(resolution))]))
        rib_pixT = rib_pixT.astype(int)
        rib_pixT = np.unique(rib_pixT, axis =0)
        # print(rib_pix.shape)
        
        rib_pix_out = rib_pix.copy()
        rib_pix_outR = rib_pixR.copy()
        rib_pix_outF = rib_pixF.copy()
        rib_pix_outT = rib_pixT.copy()
        
        index_out = []
        # Clean rib_pix if out of stack shape
        for index, pt in enumerate(rib_pix):
            if pt[0] > xdim-2 or pt[0] < 0:
                delete = True
            elif pt[1] > ydim-2 or pt[1] < 0:
                delete = True
            elif pt[2] > zdim+2-1 or pt[2] < 0:
                delete = True
            else: 
                delete = False
            
            if delete:
                # print(pt)
                index_out.append(index)
                
        rib_pix_out = np.delete(rib_pix_out, index_out, axis = 0)
        
        index_out = []
        # Clean rib_pix if out of stack shape
        for index, pt in enumerate(rib_pixR):
            if pt[0] > xdim-2 or pt[0] < 0:
                delete = True
            elif pt[1] > ydim-2 or pt[1] < 0:
                delete = True
            elif pt[2] > zdim+2-1 or pt[2] < 0:
                delete = True
            else: 
                delete = False
            
            if delete:
                # print(pt)
                index_out.append(index)
                
        rib_pix_outR = np.delete(rib_pix_outR, index_out, axis = 0)
        
        index_out = []
        # Clean rib_pix if out of stack shape
        for index, pt in enumerate(rib_pixF):
            if pt[0] > xdim-2 or pt[0] < 0:
                delete = True
            elif pt[1] > ydim-2 or pt[1] < 0:
                delete = True
            elif pt[2] > zdim+2-1 or pt[2] < 0:
                delete = True
            else: 
                delete = False
            
            if delete:
                # print(pt)
                index_out.append(index)
                
        rib_pix_outF = np.delete(rib_pix_outF, index_out, axis = 0)
    
        index_out = []
        # Clean rib_pix if out of stack shape
        for index, pt in enumerate(rib_pixT):
            if pt[0] > xdim-2 or pt[0] < 0:
                delete = True
            elif pt[1] > ydim-2 or pt[1] < 0:
                delete = True
            elif pt[2] > zdim+2-1 or pt[2] < 0:
                delete = True
            else: 
                delete = False
            
            if delete:
                # print(pt)
                index_out.append(index)
                
        rib_pix_outT = np.delete(rib_pix_outT, index_out, axis = 0)
        
        # Create mask of ring
        s3_rib = np.zeros((xdim, ydim, zdim+2))
        s3_rib[rib_pix_out[:,0],rib_pix_out[:,1],rib_pix_out[:,2]] = 1
        s3_rib[rib_pix_outR[:,0],rib_pix_outR[:,1],rib_pix_outR[:,2]] = 1
        s3_rib[rib_pix_outF[:,0],rib_pix_outF[:,1],rib_pix_outF[:,2]] = 1
        s3_rib[rib_pix_outT[:,0],rib_pix_outT[:,1],rib_pix_outT[:,2]] = 1
        
        # xdim = 10; ydim = 10; zdim = 5
        s3_filledCube = np.zeros((xdim, ydim, zdim+2))
        s3_filledCube[1:xdim-1,1:ydim-1,1:zdim+1] = 1
        s3_filledCube[rib_pix_out[:,0],rib_pix_out[:,1],rib_pix_out[:,2]] = 0
        s3_filledCube[rib_pix_outR[:,0],rib_pix_outR[:,1],rib_pix_outR[:,2]] = 0
        s3_filledCube[rib_pix_outF[:,0],rib_pix_outF[:,1],rib_pix_outF[:,2]] = 0
        s3_filledCube[rib_pix_outT[:,0],rib_pix_outT[:,1],rib_pix_outT[:,2]] = 0
        
        # test_rib = fcMeshes.createExtLayerMesh(filename, s3_rib, resolution, 'ribbon', 'ribbon', extractLargest = False, plotshow = False)
        # test_cube = fcMeshes.createExtLayerMesh(filename, s3_filledCube, resolution, 'cube', 'cube', extractLargest = True, plotshow = False)
        # test_cube.alpha(0.05)
        
        # settings.legendSize = .20
        # vp = Plotter(N=3, axes = 1)
        # vp.show(test_cube, test_rib, at=0)
        # vp.show(test_cube, at = 1)
        # vp.show(test_rib, at = 2, interactive = True)
        
        # print('rib:', test_rib.bounds())
        # print('cube:', test_cube.bounds())
            
        if not spaw_analysis: 
            for xpos in range(0,xdim):
                for zpos in range(0,zdim+2): 
                    yline = s3_filledCube[xpos,0:ydim,zpos]
                    index_y = np.where(yline == 0)[0]
                    index_y = list(index_y)
                    index_y.pop(0);index_y.pop(-1)
                    if len(index_y) > 0:
                        # print(index_x, ypos, zpos)
                        s3_filledCube[xpos,index_y[0]:ydim,zpos] = 0
        else: 
            
            for ypos in range(0,ydim):
                for zpos in range(0,zdim+2): 
                    xline = s3_filledCube[0:xdim,ypos,zpos]
                    index_x = np.where(xline == 0)[0]
                    index_x = list(index_x)
                    index_x.pop(0);index_x.pop(-1)
                    if len(index_x) > 0:
                        # print(index_x, ypos, zpos)
                        s3_filledCube[index_x[0]:xdim,ypos,zpos] = 0
                    
    
        test_rib = fcMeshes.createExtLayerMesh(filename, s3_rib, resolution, 'ribbon', 'ribbon', extractLargest = False, plotshow = False)
        test_cube = fcMeshes.createExtLayerMesh(filename, s3_filledCube, resolution, 'cube', 'cube', extractLargest = True, plotshow = False)
        test_cube.alpha(0.05)
        
        settings.legendSize = .20
        vp = Plotter(N=3, axes = 1)
        vp.show(test_cube, test_rib, at=0)
        vp.show(test_cube, at = 1)
        vp.show(test_rib, at = 2, interactive = True)
        
        s3_filledCubeBoolA = s3_filledCube.astype(bool)
        
        return s3_filledCubeBoolA

#%% def function - cutLR then atr/vent
    def divideMeshesLnR_newUp(filename, df_res, directories, dir_data2Analyse, mesh, s3_filledCubeBoolA, res, 
                              scale_cube = [], colors =  ['skyblue','darkblue'], saveStacks = True, chambers = True):
        
        spaw_analysis = False
        if 'spaw_ct' in df_res.loc[file_num,'spAnalysis']:
            spaw_analysis = True
        
        mesh_legend = mesh._legend
        if not spaw_analysis: 
            print('\n- Dividing '+mesh_legend+' into left and right sides')
            names_LnR = ['CJ_total','CJ.Left', 'CJ.Right']
            names_LnRm = ['CJ_total','-Left', '-Right']
        else: 
            print('\n- Dividing '+mesh_legend+' into dorsal and ventral sides')
            names_LnR = ['CJ_total','CJ.Dorsal', 'CJ.Ventral']
            names_LnRm = ['CJ_total','-Dorsal', '-Ventral']
            
         # s3_filledCubeBoolA = s3_filledCube.astype(bool)
        [s3], _ = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1], end_name = ['cj'], print_txt = False)
        
        masked_s3A = s3.copy()
        masked_s3A[s3_filledCubeBoolA] = 0
        test_cutA = fcMeshes.createExtLayerMesh(filename, masked_s3A, res, 'cjcutA', 'cjcutA', extractLargest = False, plotshow = False)
        test_cutA.legend(mesh_legend+names_LnRm[1]).color(colors[0])
        
        s3_filledCubeBoolB = np.invert(s3_filledCubeBoolA)
        masked_s3B = s3.copy()
        masked_s3B[s3_filledCubeBoolB] = 0
        test_cutB = fcMeshes.createExtLayerMesh(filename, masked_s3B, res, 'cjcutB', 'cjcutB', extractLargest = False, plotshow = False)
        test_cutB.legend(mesh_legend+names_LnRm[2]).color(colors[1])
        
        settings.legendSize = .20
        vp = Plotter(N=3, axes = 1)
        vp.show(mesh, at=0)
        vp.show(test_cutA, scale_cube, at = 1)
        vp.show(test_cutB, scale_cube, at = 2, interactive = True)
        
        if saveStacks:
            fcCont.save_s3(filename = filename, s3 = masked_s3A, dir_txtNnpy = directories[1], layer = 'cj_AOCVIC')
            fcCont.save_s3(filename = filename, s3 = masked_s3B, dir_txtNnpy = directories[1], layer = 'cj_AICVOC')
    
        print(m_cj.volume())
        print(test_cutA.volume())
        print(test_cutB.volume())
        print(test_cutA.volume()+test_cutB.volume())
        print(((test_cutA.volume()+test_cutB.volume())-m_cj.volume())/m_cj.volume()*100)
        
        meshes_LnR = [test_cutA, test_cutB]
        
        return meshes_LnR, names_LnR

    #%% SELECT FILE AND GET METADATA
    #   This section allows the user to select file to process, get its main directories and metadata,
    #   and define some properties
    #   ================================================================================================================

    # Get main directories
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = 'R', out_type = 'xlsx')
    # Get file to process and directories
    df_files = fcBasics.selectHearts(df_dataset)
    for nn, folder in zip(count(), df_files['Folder']):
        df_file = df_files.iloc[[nn]]
        filename = folder[2:]
        dORv = filename[9:10]
            
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
            
        # Initialise variables
        dict_shapes = dict()
        txt = Text2D(filename, c="k", font= 'CallingCode')
        rotAngle =  df_res.loc[file_num,'ang_HeartS']
        plot = True
        
        #%%
        # Get existing cl_dictionaries
        [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes, dicts_cl] = fcBasics.import_dicts('mH_C', filename, directories)
        
         # Import meshes
        [m_myoc, m_cj] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc','cj'],
                                                                      extension = 'vtk', dir_stl = directories[2],
                                                                      alpha = [1,1], dict_colour = dict_colour)
        
        
        
        scale_cube = Cube(pos=m_myoc.centerOfMass(), side=350, c='white', alpha=0.01).legend('Cube')
        
        kspl_CL, linLines, ksplSph_o, dict_kspl = fcMeshes.createCLs(filename = filename, dict_cl = dicts_cl, dict_pts = dict_pts, 
                                                                      dict_kspl = dict_kspl, dict_planes = dict_planes, 
                                                                      colors = ['deepskyblue', 'tomato'], myoc = m_myoc, 
                                                                      dir_stl = directories[3])
        sph_CL, sph_CL_colour, dict_shapes = fcMeshes.createColouredCL(dict_cl= dicts_cl, dict_shapes = dict_shapes)
        
        # Create a disc with better resolution to transform into pixels to mask stack
        cl_ribbon, kspl_ext, dict_kspl, dict_shapes, dict_planes = fcMeshes.createCLRibbon(filename = filename, file_num = file_num, 
                                                                                            df_res = df_res, kspl_CL2use = kspl_CL[0], 
                                                                                            linLine = linLines[0], mesh = m_myoc, 
                                                                                            dict_kspl = dict_kspl, dict_shapes = dict_shapes, 
                                                                                            dict_planes = dict_planes, scale_cube = scale_cube, 
                                                                                            plotshow = True)
        
        s3_filledCubeBoolA = getCubeRibbonMask(filename, file_num, df_res, cl_ribbon, res, directories[1])
        
        # Just Cardiac Jelly (Total)
        mesh = m_cj
        m_cjLnR , names_LnR = divideMeshesLnR_newUp(filename, df_res, directories, dir_data2Analyse, mesh, s3_filledCubeBoolA, res, 
                                                      scale_cube = [], colors =  ['skyblue','darkblue'], chambers = True)
        
        print(names_LnR)
                   
        df_res = fcMeshes.addLayersVolume2df (df_res = df_res, file_num = file_num, meshes = [m_cj]+m_cjLnR, 
                                              names = names_LnR)
        txt0 = Text2D("\n >> "+names_LnR[1]+": "+format(m_cjLnR[0].volume(), '.1f'), c=c, font=font)
        txt1 = Text2D("\n >> "+names_LnR[2]+": "+format(m_cjLnR[1].volume(), '.1f'), c=c, font=font)
        vp = Plotter(N=2, axes=13)
        vp.show(m_cjLnR[0], scale_cube, txt0, at=0)
        vp.show(m_cjLnR[1], scale_cube, txt1, at=1, zoom=1.6, azimuth = 0, elevation = 0, interactive=True)
        
        save = fcBasics.ask4input('Save? [0]:no, [1]:yes! >>:', bool)
        if save: 
            fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
            dir_df_meas = fcBasics.new_dir(fcBasics.new_dir(fcBasics.new_dir(dir_data2Analyse, 'R_All'), 'df_all'), 'df_meas')
            fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_df_meas)#, name = 'ResultsDFf')
            
        
        # Get Atrium and Ventricle LnR
        cyl_chamber = dict_shapes['cyl2CutChambers_final']
        disc = Cylinder(pos = cyl_chamber['cyl_centre'],r = cyl_chamber['radius_max'], height = 2*0.225, 
                    axis = cyl_chamber['cyl_axis'], c = 'purple', cap = True, res = 300)
        num_pt = dict_pts['numPt_CLChamberCut']
        
        s32cut = ['cj_AOCVIC','cj_AICVOC']
        s32cutn = ['CJ.AOCVIC','CJ.AICVOC']
        m_atr, m_vent, dict_shapes, _ = fcMeshes.getChamberMeshes(filename = filename,
                                    end_name = s32cut,
                                    names2cut = s32cutn,
                                    kspl_CL = kspl_CL[0], num_pt = num_pt, atr_meshes = [], vent_meshes = [],
                                    dir_txtNnpy = directories[1], dict_shapes = dict_shapes, dict_pts = dict_pts, 
                                    resolution = res, s3_cyl = [], plotshow = plotshow,
                                    colour_anv = ['lightseagreen','darkturquoise'])
        m_AOC, m_AIC = m_atr
        m_VIC, m_VOC = m_vent
        
        m_AOC.legend('CJ.AOC')
        m_AIC.legend('CJ.AIC')
        m_VOC.legend('CJ.VOC')
        m_VIC.legend('CJ.VIC')
        
        settings.legendSize = .20
        vp = Plotter(N=4, axes = 1)
        vp.show(m_AOC, scale_cube, at=0)
        vp.show(m_VIC, scale_cube, at = 1)
        vp.show(m_AIC, scale_cube, at = 2)
        vp.show(m_VOC, scale_cube, at = 3, interactive = True)
        
        df_res = fcMeshes.addLayersVolume2df (df_res = df_res, file_num = file_num, 
                                              meshes = [m_cjLnR[0],m_AOC, m_VIC], 
                                              names = ['CJ.AOC-VIC','CJ.AOC', 'CJ.VIC'])
        df_res = fcMeshes.addLayersVolume2df (df_res = df_res, file_num = file_num, 
                                              meshes = [m_cjLnR[1],m_AIC, m_VOC], 
                                              names = ['CJ.AIC-VOC','CJ.AIC', 'CJ.VOC'])
        
        save = True#fcBasics.ask4input('Save? [0]:no, [1]:yes! >>:', bool)
        if save: 
            fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
            dir_df_meas = fcBasics.new_dir(fcBasics.new_dir(fcBasics.new_dir(dir_data2Analyse, 'R_All'), 'df_all'), 'df_meas')
            fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_df_meas)#, name = 'ResultsDFf')
            
        # q_nextHeart = fcBasics.ask4input('Next Heart [0]:no, [1]:yes! >>:', bool)
        # if q_nextHeart:
        #     continue
        # else: 
        #     print('Last heart number ', nn)
        #     break
    
    fcBasics.alert('frog', 1)
    fcBasics.alert('jump', 2)
    
#%% Init
init = True    
    
    
    
    
    
    
    
    
    
    