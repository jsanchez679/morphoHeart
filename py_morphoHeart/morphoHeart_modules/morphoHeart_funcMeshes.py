# -*- coding: utf-8 -*-
"""
morphoHeart_funcMeshes

Version: September 30, 2020
@author: Juliana Sanchez-Posada

"""

#%% Importing python packages
import numpy as np
import random
import math
import os
import seaborn as sns
import matplotlib.pyplot as plt

from itertools import count

from skimage import measure#, io
from scipy.interpolate import splprep, splev, interpn
from scipy.spatial.distance import cdist
import pandas as pd

# from datetime import datetime
from progress.bar import Bar
suffix = '%(index)d/%(max)d - %(elapsed)ds'

c="k"
font= 'VTK'

# from vedo import *
from vedo import embedWindow, Plotter, settings, load, Text2D, Mesh, KSpline, Sphere, Plane, Ribbon, Points, merge, shapes
from vedo import fitPlane, Spheres, Line, Cylinder, colorMap, Cube, Arrow, Video, vedo2trimesh, trimesh2vedo, Ellipsoid
import trimesh
from time import perf_counter
embedWindow(False)

import json
#from json import JSONEncoder

#%% Importing morphoHeart packages
from .morphoHeart_funcBasics import alert, ask4input, loadNPY, saveDF, loadDF, getInputNumbers, new_dir #, saveDict
from .morphoHeart_funcContours import save_s3, loadStacks, drawLine #, plt_s3, save_s3s

#%% class - NumpyArrayEncoder
# Definition of class to save dictionary
class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NumpyArrayEncoder, self).default(obj)

#%% - LOADING
#%% func - openMeshes
def openMeshes(filename, meshes_names, extension, dir_stl, alpha, dict_colour):
    """
    Function to load a list of meshes

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    meshes_names : list of strings
        List with the names of the meshes to load.
    extension : str
        Extension to saved mesh 'vtk'/'stl'.
    dir_stl : path
        Path to the folder where the meshes are saved.
    alpha : list of floats
        List of opacity values to assign to each mesh
    dict_colour : dictionary
        dictionary with information about the meshes colour

    Returns
    -------
    meshes_out : list of meshes
        output List of meshes in same order as meshes_names.

    """

    names_all = ['myoc','myoc_ext','myoc_int','myoc_atr','myoc_vent',
                  'endo','endo_ext','endo_int','endo_atr','endo_vent',
                  'cj','cj_out','cj_in','cj_atr','cj_vent',
                  'cjExt_atr','cjExt_vent', 'cjExt_atrmyoc_vent', 
                  'myocExt_atr', 'myocExt_vent', 'endoExt_atr', 'endoExt_vent']

    legend_all = ['Myocardium','Ext.Myoc', 'Int.Myoc','Atrium(Myoc)','Ventricle(Myoc)',
                  'Endocardium', 'Ext.Endo', 'Int.Endo','Atrium(Endo)','Ventricle(Endo)',
                  'CardiacJelly','Ext.CJ','Int.CJ','Atrium(CJ)','Ventricle(CJ)',
                  'CJ.ExtAtr', 'CJ.ExtVent','CJ.ExtAtr', 
                  'Myoc.ExtAtr', 'Endo.ExtAtr', 'Myoc.ExtVent', 'Endo.ExtVent']

    print('- Loading meshes...')
    meshes_out = []

    for i, name in enumerate(meshes_names):

        index = names_all.index(name)
        mesh_title = filename+"_"+name+"."+extension
        mesh_dir = os.path.join(dir_stl, mesh_title)
        mesh_out = load(mesh_dir)
        print("\t>> Mesh loaded - "+mesh_title+"!")
        try:    
            mesh_colour = dict_colour[name]['colour']
        except: 
            mesh_colour = 'chocolate'
        mesh_out.alpha(alpha[i]).legend(legend_all[index]).wireframe().color(mesh_colour)

        meshes_out.append(mesh_out)

    alert("wohoo",1)
    # print('\n')

    return meshes_out

#%% func - openThicknessMeshes
def openThicknessMeshes(filename, meshes_names, extension, dir_stl, dir_txtNnpy, print_txt = True):
    """
    Function to load a list of meshes that are coloured by a particular property (e.g. thickness, ballooning)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    meshes_names : list of strings
        List with the names of the meshes to load.
    extension : str
        Extension to saved mesh 'vtk'/'stl'.
    dir_stl :  path
        Path to the folder where the meshes are saved.
    dir_txtNnpy :  path
        Path to the folder where the np arrays are saved.
    print_txt : Bool
        True if text is to be printed, else False. The default is True.

    Returns
    -------
    meshes_out : list of meshes
        output List of meshes in same order as meshes_names.
    colour_arrays : list of numpy arrays
        list of numpy arrays containing information of thickness/ballooning to color meshes

    """

    names_all = ['myoc_intBall', 'myoc_extBall', 'myoc_thickness','endo_thickness','cj_thickness']
    legend_all = ['Int.Myoc Ball.','Ext.Myoc Ball.','Myoc.Thickness','Endo.Thickness','CJ.Thickness']
    bar_names_all = ['Int.Myoc\nBalloning\n[um]','Ext.Myoc\nBalloning\n[um]','Myoc.Thickness\n[um]','Endo.Thickness\n[um]','CJ.Thickness\n[um]']
    alpha_all = [1,1,0.1,0.1,0.1]

    if print_txt:
        print('- Loading thickness meshes...')
    meshes_out = []
    colour_arrays = []

    for i, name in zip(count(), meshes_names):
        # Get index
        index = names_all.index(name)
        # Load mesh
        mesh_title = filename+"_"+name+"."+extension
        mesh_dir = os.path.join(dir_stl, mesh_title)
        mesh_out = load(mesh_dir)
        if print_txt:
           print("\t>> Mesh loaded - "+mesh_title+"!")
        # Load colour array
        colour_title = filename+"_"+name+".npy"
        colour_dir = os.path.join(dir_txtNnpy, colour_title)
        sp_colour = np.load(colour_dir)


        mesh_out.pointColors(sp_colour, cmap="jet")
        mesh_out.addScalarBar(title=bar_names_all[index])

        mesh_out.alpha(alpha_all[index]).legend(legend_all[index]).wireframe()

        meshes_out.append(mesh_out)
        colour_arrays.append(sp_colour)

    if print_txt:
        alert("wohoo",1)
        print('\n')

    return meshes_out, colour_arrays

#%% func - load_tissues2unloop 
def load_tissues2unloop(filename, directories, dir_results, dict_colour):
    
    q_tissue_analysis = ask4input('Select the heart tissues you want to unloop and get 2D heatmaps out: \n\t[0]: cardiac jelly thickness \n\t[1]: balloning of internal myocardium \n\t[2]: myocardial thickness \n\t[3]: endocardial thickness \n\t[all]: All of them!>> :', str)
    tissue_opt = ['CjTh','myocIntBall','MyocTh','EndoTh']
    tissue_analysis = [tissue_opt[i] for i in getInputNumbers(q_tissue_analysis, tissue_opt)]
    if 'CjTh' in tissue_analysis and 'myocIntBall' in tissue_analysis: 
        tissue_analysis.remove('CjTh'); tissue_analysis.remove('myocIntBall')
        tissue_analysis.append('CjThNmyocIntBall')
    tissue_analysis = sorted(tissue_analysis, key=str.casefold)
    tissue_opt.append('CjThNmyocIntBall')
    
    # Import myocardium
    [m_myoc] = openMeshes(filename = filename, meshes_names = ['myoc'], extension = 'vtk', dir_stl = directories[2],
                                                                  alpha = [0.05], dict_colour = dict_colour)
    # Meshes names to load
    meshes_names = [['cj_thickness'],['myoc_intBall'], ['myoc_thickness'], ['endo_thickness'],['cj_thickness','myoc_intBall']]
    # Dataframe names of points classification
    df_names = ['df_cjThNmyocIntBall', 'df_cjThNmyocIntBall', 'df_myocTh' ,'df_endoTh', 'df_cjThNmyocIntBall']
    # Select chambers to unloop
    unloopnames = ['unloopAtr', 'unloopVent']
    chambers = ['Atrium', 'Ventricle']
    q_both_chs = ask4input('Select the chambers you would like to unloop for the selected tissues: \n\t[0]: atrium\n\t[1]: ventricle\n\t[all/0,1/0-1]: both! >> : ', str)
    index_selected_chambers = getInputNumbers(q_both_chs, chambers)

    dict_unloop = dict()
    for n, tissue, m_name, df_name in zip(count(), tissue_opt, meshes_names, df_names): 
        if tissue in tissue_analysis:
            dict_unloop[tissue] = dict_tissue = dict()
            print(m_name, df_name)
            [m_Th, df_thickness] = openThicknessMeshes(filename = filename, meshes_names = m_name, 
                                            extension = 'vtk', dir_stl = directories[2], dir_txtNnpy = directories[1])
            
            df_Th = loadDF(filename = filename, file = df_name, dir_results = os.path.join(dir_results, 'csv_all'))
            
            dict_tissue['dataname'] = tissue#[tissue]
            dict_tissue['mesh'] = m_Th[0]
            dict_tissue['pts_whole'] = m_Th[0].points()
            dict_tissue['pts_class'] =  [df_Th['AtrVent'], df_Th['DorsVent'], df_Th['LeftRight']]
            # if tissue == 'cjThNmyocIntBall':
            #     dict_tissue['param'] = df_thickness
            # else: 
            dict_tissue['param'] = df_thickness
            dict_tissue['param_name'] = m_name
            dict_tissue['df_AtrVent'] = df_Th['AtrVent']
            selected_chambers = [chambers[cc] for cc in index_selected_chambers]
            dict_tissue['selected_chambers'] = selected_chambers
            save_names = [unloopnames[cc]+tissue  for cc in index_selected_chambers]
            dict_tissue['save_names'] = save_names
            
            if tissue == 'CjThNmyocIntBall':
                var_hm_names = []
                for sp_tissue in ['CjTh','myocIntBall']:
                    hm_names = [unloopnames[cc]+'_'+sp_tissue  for cc in index_selected_chambers]
                    var_hm_names.append(hm_names)
            else: 
                var_hm_names = [[unloopnames[cc]+'_'+tissue  for cc in index_selected_chambers]]
            dict_tissue['hm_names'] = var_hm_names

    return tissue_analysis, tissue_opt, m_myoc, dict_unloop, 
    
#%% - DICTIONARIES
#%% func - fillNsaveObjDict
def fillNsaveObjDict(filename, dicts, names, dir2save):
    """
    Function to create an object dictionary of all the objects created for each of the analysed hearts

    Parameters
    ----------
    filename :  str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    dicts : list of dictionaries
        [dict_planes, dict_pts, dict_kspl, dict_colour(, dict_shapes)].
    names : list of names for dicts
        List of the names to use as keys to create a dict of dicts  ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour'(, 'dict_shapes')].
    dir2save : path
        Path to the folder where the dictionaries are saved.

    Returns
    -------
    dict_obj : dictionary
        Dictionary of dictionaries [dict_planes, dict_pts, dict_kspl, dict_colour(, dict_shapes)].

    """
    print('- Merging dictionaries...')
    dict_obj = dict()
    for i, sp_dict, name in zip(count(), dicts, names):
        print('\t >> '+name)
        dict_obj[name] = sp_dict

    jsonDict_name = filename+"_dict_obj.json"
    json2save_dir = os.path.join(dir2save,jsonDict_name)

    with open(json2save_dir, "w") as write_file:
        json.dump(dict_obj, write_file, cls=NumpyArrayEncoder)
    print("\t>> dict_obj saved correctly!\n\t> File: "+jsonDict_name);
    alert("countdown",1)

    return dict_obj

#%% func - splitDicts - moved to Basics 
# def splitDicts(dict_obj):
#     """
#     Function that splits the dictionary of dictionaries into a list of dictionaries

#     Parameters
#     ----------
#     dict_obj : dictionary
#         Dictionary of dictionaries [dict_planes, dict_pts, dict_kspl, dict_colour(, dict_shapes)].

#     Returns
#     -------
#     dicts_all : list of dictionaries
#         [dict_planes, dict_pts, dict_kspl, dict_colour(, dict_shapes)].

#     """

#     dicts_all = []
#     for i, name in enumerate(list(dict_obj.keys())):
#         ext_dict = dict_obj[name]
#         dicts_all.append(ext_dict)

#     return dicts_all

#%% func - addPlanes2Dict
def addPlanes2Dict (planes, pls_centre, pls_normal, info, dict_planes, print_txt = True):
    """
     Function that adds a list of planes to dict_planes

    Parameters
    ----------
    planes : List of Planes
        Planes to add to dict - (vedo Planes)
    pls_centre : list of list of floats
        list of List with the x,y,z coordinates of each of the planes' centre
    pls_normal : list of list of floats
        list of List with the x,y,z coordinates of each of the planes' normal
    info :  list of strs
        Additional text to include in the planes names.
    dict_planes : dictionary
        Initialised dictionary with planes information
    print_txt :  bool, optional
        True if confirmation of action is needed, else False. The default is True.

    Returns
    -------
    dict_planes : dictionary
        Resulting dictionary with planes information updated

    """

    for i, plane, inf in zip(count(), planes, info):
        if inf =='':
            specific_plane = dict_planes[plane._legend]= dict()
        else:
            specific_plane = dict_planes[plane._legend+'-'+inf]= dict()

        specific_plane['pl_normal'] = pls_normal[i]
        specific_plane['pl_centre'] = pls_centre[i]
        specific_plane['color'] = plane.color()

    if print_txt:
        print('\n>> Planes have been added to dict_planes!')

    return dict_planes

#%% func - addPoints2Dict
def addPoints2Dict (spheres, info, dict_pts, print_txt = True):
    """
    Function that adds a list of spheres to dict_pts

    Parameters
    ----------
    spheres : List of Spheres
        List of Spheres to add to dict - (vedo Spheres)
    info : list of strs
        Additional text to include in the planes names.
    dict_pts : dictionary
        Initialised dictionary with points information
    print_txt :  bool, optional
        True if confirmation of action is needed, else False. The default is True.

    Returns
    -------
    dict_pts : dictionary
        Resulting dictionary where the list of 'spheres' was added

    """

    for i, sph, inf in zip(count(), spheres, info):
        if inf =='':
            specific_pt = dict_pts[sph._legend]= dict()
        else:
            specific_pt = dict_pts[sph._legend+'-'+inf]= dict()

        specific_pt['sph_position'] = sph.pos()
        specific_pt['color'] = sph.color()

    if print_txt:
        print('\n>> Points have been added to dict_pts!')

    return dict_pts

#%% func - addKSplines2Dict
def addKSplines2Dict (kspls, info, dict_kspl, print_txt = True):
    """
    Function that adds a list of splines to dict_kspl

    Parameters
    ----------
    kspls : List of Ksplines
        List of Kspline to add to dict (vedo KSplines)
    info : list of strs
        Additional text to include in the planes names.
    dict_kspl : dictionary
        Initialised dictionary with kspline information
    print_txt : bool, optional
        True if confirmation of action is needed, else False. The default is True.

    Returns
    -------
    dict_kspl : dictionary
        Resulting dictionary where the list of 'ksplines' was added

    """

    for i, kspl, inf in zip(count(), kspls, info):
        #print(i)
        if inf =='':
            specific_kspl = dict_kspl[kspl._legend]= dict()
        else:
            specific_kspl = dict_kspl[kspl._legend+'-'+inf]= dict()

        specific_kspl['kspl_pts'] = kspl.points()
        specific_kspl['color'] = kspl.color()

    if print_txt:
        print('\n>> KSplines have been added to dict_kspl!')

    return dict_kspl

#%% func - addShapes2Dict
def addShapes2Dict(shapes, dict_shapes, radius, print_txt = True):
    """
    Function that adds a list of shapes to dict_kspl

    Parameters
    ----------
    shapes : list of shapes
        List of shapes (ribbons, spheres, etc) to add to dict- (vedo Objects)
    dict_shapes : dictionary
        Initialised dictionary with shapes information
    radius : list of floats
        List of floats with the radii of the spheres (if applicable)
    print_txt : bool, optional
        True if confirmation of action is needed, else False. The default is True.

    Returns
    -------
    dict_shapes : dictionary
        Resulting dictionary where the list of 'shapes' was added

    """
    
    for i, sp_shape, radii in zip(count(), shapes, radius):
        name = sp_shape._legend
        specific_shape = dict_shapes[sp_shape._legend]= dict()
        
        if 'sphs_maxInsSphRadC' in name:
            specific_shape['maxInscSphRad'] = radii
        if 'sphs_maxInsSphRad' in name:
            specific_shape['kspl_pts'] = sp_shape.points()
        elif 'cyl2CutChambers_o' in name: 
            specific_shape['radius_max'] = radii[0]
            specific_shape['radius_min'] = radii[1]
            specific_shape['cyl_axis'] = radii[2]
            specific_shape['cyl_centre'] = radii[3]
        elif 'cyl2CutChambers_final' in name: 
            specific_shape['radius_max'] = radii[0]
            specific_shape['radius_min'] = radii[1]
            specific_shape['radius_num'] = radii[2]
            
            specific_shape['height_max'] = radii[3]
            specific_shape['height_min'] = radii[4]
            specific_shape['height_num'] = radii[5]
            
            specific_shape['cyl_axis'] = radii[6]
            specific_shape['cyl_centre'] = radii[7]
            specific_shape['resolution'] = radii[8]
        elif 'Ellip' in name: 
            specific_shape['center'] = radii[0]
            specific_shape['x_dim'] = radii[1]
            specific_shape['y_dim'] = radii[2]
            specific_shape['z_dim'] = radii[3]
        else:
            specific_shape['Points'] = sp_shape.points()
            specific_shape['NoPoints'] = len(sp_shape.points())
            
        specific_shape['color'] = sp_shape.color()

    if print_txt:
        print('\n>> Shapes have been added to dict_shapes!')

    return dict_shapes


    
#%% - DATAFRAMES
#%% func - addSurfArea2df
def addSurfArea2df (df_res, file_num,  meshes):
    """
    Function that measures the surface area of a list of meshes given as input and adds it to an existing dataframe

    Parameters
    ----------
    df_res : dataframe
        Dataframe with the measured information of the heart being processed.
    file_num : int
        Index number of the selected heart being processed.
    meshes : list of meshes (vedo Meshes)
        List of meshes from which the surface area wants to be measured.

    Returns
    -------
    df_resFilled : dataframe
        Dataframe with the updated measured information of the heart being processed.

    """

    names = ['Myoc', 'Int.Myoc', 'Ext.Myoc', 
             'Endo', 'Int.Endo', 'Ext.Endo', 
             'Myoc', 'Int.Myoc', 'Ext.Myoc', 
             'Endo', 'Int.Endo', 'Ext.Endo', 
             'CJ', 'Int.CJ', 'Ext.CJ', 'CJ',
             'Atr.ExtMyoc','Atr.IntEndo','Atr.ExtCJ',
             'Vent.ExtMyoc','Vent.IntEndo','Vent.ExtCJ']
    legends = ['Myocardium', 'Int.Myoc', 'Ext.Myoc', 
               'Endocardium', 'Int.Endo', 'Ext.Endo', 
               'Cut (Myoc)', 'Cut (Int.Myoc)', 'Cut (Ext.Myoc)', 
               'Cut (Endo)', 'Cut (Int.Endo)', 'Cut (Ext.Endo)', 
               'CardiacJelly', 'Int.CJ', 'Ext.CJ', 'CJ',
               'Ext.Myoc_Atr','Int.Endo_Atr','Ext.CJ_Atr', 
               'Ext.Myoc_Vent','Int.Endo_Vent','Ext.CJ_Vent']
    
    df_resFilled = df_res.copy()

    for n, mesh in enumerate(meshes):
        leg_mesh = mesh._legend
        name = names[legends.index(leg_mesh)]
        area = mesh.area()
        
        df_resFilled.loc[file_num,'SurfArea_'+name] = area

        area_print = format(area, '.1f')
        print('-', name, '- SurfArea: ', area_print, 'um^2')

    return df_resFilled

#%% func - addLayersVolume2df
def addLayersVolume2df (df_res, file_num,  meshes, names = ['Myoc', 'Atr.Myoc', 'Vent.Myoc','Atr.ExtMyoc', 'Vent.ExtMyoc',
             'Endo', 'Atr.Endo', 'Vent.Endo', 'Atr.IntEndo', 'Vent.IntEndo','CJ', 'Atr.CJ', 'Vent.CJ']):
    """
    Function that measures the volume of a list of meshes given as input and adds it to an existing dataframe

    Parameters
    ----------
    df_res : dataframe
        Dataframe with the measured information of the heart being processed.
    file_num : int
        Index number of the selected heart being processed.
    meshes :  list of meshes (vedo Meshes)
        List of meshes from which the volume wants to be measured.

    Returns
    -------
    df_resFilled : dataframe
        Resulting dataframe in which the volume information was added.

    """

    # names = ['Myoc', 'Atr.Myoc', 'Vent.Myoc','Atr.ExtMyoc', 'Vent.ExtMyoc',
    #          'Endo', 'Atr.Endo', 'Vent.Endo', 'Atr.IntEndo', 'Vent.IntEndo',
    #          'CJ', 'Atr.CJ', 'Vent.CJ']

    df_resFilled = df_res

    print('- Getting volumes...')
    for n, name, mesh in zip(count(), names, meshes):
        volume = mesh.volume()
        df_resFilled.loc[file_num,'Vol_'+name] = volume

        vol_print = format(volume, '.1f')
        #print('Name:', name, '\t\t- Volume: ', vol_print, 'um^3')

        val = (n+1)%5
        #print(val)
        if val == 1:
            vol_all = volume
            print('\n-> ', name, '\t\t- Volume: ', vol_print, 'um^3')
        elif val == 2:
            vol_atr = volume
            print('-> ', name, '\t- Volume: ', vol_print, 'um^3')
        elif val == 3:
            vol_vent = volume
            print('-> ', name, '\t- Volume: ', vol_print, 'um^3')
            diff_print = format(vol_all - (vol_atr+vol_vent), '.1f')
            diff_perc = format(100*(1-(vol_atr+vol_vent)/vol_all), '.1f')
            print('\t - Vol difference: ', diff_print,  'um^3 - (', diff_perc,'%)')
        else:
            print('-> ', name, '\t- Volume: ', vol_print, 'um^3')

    return df_resFilled

#%% func - addIntExtVol2df
def addIntExtVol2df(df_res, file_num,  meshes):
    """
    Function that measures the volume of a list of external meshes given as input and adds it to an existing dataframe

    Parameters
    ----------
    df_res : dataframe
        Dataframe with the measured information of the heart being processed.
    file_num : int
        Index number of the selected heart being processed.
    meshes : list of meshes (vedo Meshes)
        List of meshes from which the volume wants to be measured.

    Returns
    -------
    df_resFilled : dataframe
        Resulting dataframe in which the volume information was added.

    """
    names = ['Int.Myoc', 'Ext.Myoc', 'Int.Endo', 'Ext.Endo']
    df_resFilled = df_res

    print('- Getting internal and external volumes of heart layers...')
    for n, name, mesh in zip(count(), names, meshes):
        volume = mesh.volume()
        df_resFilled.loc[file_num,'Vol_'+name] = volume

    return df_resFilled

#%% func - addLinearMeas2df
def addLinearMeas2df(df_res, file_num, lines, kspl_CL):
    """
    Function that measures the length of a list of lines/ksplines given as input and adds it to an existing dataframe

    Parameters
    ----------
    df_res : dataframe
        Dataframe with the measured information of the heart being processed.
    file_num : int
        Index number of the selected heart being processed.
    lines : list of lines (vedo Lines)
        List of lines from which the length wants to be measured.
    kspl_CL : list of ksplines (vedo KSplines)
        List of ksplines from which the length wants to be measured.

    Returns
    -------
    df_resFilled :  dataframe
        Resulting dataframe in which the length measurements were added.

    """
    df_resFilled = df_res
    for nl, line in enumerate(lines):
        length_l = line.length()
        name_l = line._legend
        df_resFilled.loc[file_num, name_l] = length_l

    for ncl, cl in enumerate(kspl_CL):
        length_cl = cl.length()
        name_cl = cl._legend
        df_resFilled['Length_'+name_cl] = length_cl

    return df_resFilled

#%% func - addOrientationAngles2df
def addOrientationAngles2df(df_res, file_num, angles, names):
    """
    Function that adds a list of angles (orientations) to an existing dataframe

    Parameters
    ----------
    df_res : dataframe
        Dataframe with the measured information of the heart being processed.
    file_num : int
        Index number of the selected heart being processed.
    angles : list of floats
        List with the angle values to be added to dataframe.
    names : list of str
        List with the names that will be given to the angle values in the dataframe.

    Returns
    -------
    df_resFilled : dataframe
        Resulting dataframe in which the angle measurements were added.

    """

    df_resFilled = df_res
    for i, angle, name in zip(count(), angles, names):
        df_resFilled.loc[file_num, name] = angle

    return df_resFilled

#%% - S3s MANIPULATION
#%% func - maskChamberS3s
def maskChamberS3s (s3_mask, pl_normal, pl_centre, resolution):
    """
    Function used to cut the heart into chambers (atrium and ventricle).
    For this, the s3_mask given as input is cut with a plane whose normal and centre are also given, creating two
    new arrays containing each of the chambers.

    Parameters
    ----------
    s3_mask : numpy array of booleans
        Array with information about the heart layer (int/ext/all). The size of this array corresponds to the size of the original stack.
    pl_normal : list of floats
        List with the x,y,z coordinatesof the plane's normal
    pl_centre : list of floats
        List with the x,y,z coordinates of the plane's centre
    resolution : list of floats
        List with the x,y, z scaling values of the images taken. This information is taken from the metadata of the original file.

    Returns
    -------
    s3_atr : numpy array of booleans
        Resulting array with information of the atrium of the heart. Final size is same size as the input array.
    s3_vent : numpy array of booleans
        Resulting array with information of the ventricle of the heart. Final size is same size as the input array.

    """

    # Get dimensions of stack
    xdim, ydim, zdim = s3_mask.shape
    # Reshape stack as a vector
    s3_mask_v = s3_mask.reshape(-1)
    # Get vectors of x,y and z positions in array
    pix_coord_pos = np.where(s3_mask >= 0)
    del s3_mask

    # Trasform coordinate positions to um using resolution
    pix_um = np.transpose(np.asarray([pix_coord_pos[i]*resolution[i] for i in range(len(resolution))]))
    del pix_coord_pos

    # Make normal of the input plane a unit vector
    normal_unit = unit_vector(pl_normal)
    # Find all the d values of pix_um
    d_pix_um = np.dot(np.subtract(pix_um,np.array(pl_centre)),np.array(normal_unit))
    del pix_um

    # Clear vector d_pix_um using only those that are 1 in stack
    d_pve_pix_um = s3_mask_v*d_pix_um
    del d_pix_um

    # Duplicate s3_mask_v to initialise atrium and ventricle
    s3_vent_v = np.copy(s3_mask_v)
    s3_vent_v = s3_vent_v.astype('uint8')
    s3_atr_v = np.copy(s3_mask_v)
    del s3_mask_v

    s3_atr_v = s3_atr_v.astype('uint8')
    # Find all positions in d_pve_pix_um that are at either side of the plane
    pos_atr = np.where(d_pve_pix_um < 0)[0]
    pos_vent = np.where(d_pve_pix_um > 0)[0]
    del d_pve_pix_um

    # Remove the points that belong to the other chamber
    s3_vent_v[pos_atr] = 0
    s3_atr_v[pos_vent] = 0
    del pos_atr, pos_vent

    # Reshape vector into matrix/stack
    s3_vent = s3_vent_v.reshape((xdim, ydim, zdim))
    s3_atr = s3_atr_v.reshape((xdim, ydim, zdim))
    del s3_vent_v, s3_atr_v

    # bar.finish()
    alert('wohoo',1)

    return s3_atr, s3_vent

#%% #%% func - maskRing2CutChamberS3s
def maskRing2CutChamberS3s (s3_mask, s3_cyl, stack_shape):
    """
    Function used to cut the heart layers into chambers (atrium and ventricle).
    For this, the s3_mask given as input is cut with a ring whose mask is alalso given, creating one array 
    with a ring dividing the input mask into two distinctive chambers

    Parameters
    ----------
    s3_mask : numpy array of booleans
        Array with information about the heart layer (int/ext/all). The size of this array corresponds to the size of the original stack.
    s3_cyl : numpy array of booleans
        Array with ring to cut tissue layer into chambers. The size of this array corresponds to the size of the original stack.
    stack_shape : list
        List containing the stack shape [x,y,z]

    Returns
    -------
    s3_mask : numpy array of booleans
        Resulting array masked with ring to cut tissue layer into chambers. Final size is same size as the input array.

    """

    xdim, ydim, zdim = stack_shape
    
    bar = Bar('- Masking cut', max=zdim, suffix = suffix, check_tty=False, hide_cursor=False)
    for slc in range(zdim):
        im_cyl =s3_cyl[:,:,slc]
        pos_pts = np.where(im_cyl == 1)
        im = s3_mask[:,:,slc]
    
        clicks = [(pos_pts[0][i], pos_pts[1][i]) for i in range(pos_pts[0].shape[0])]
        if len(clicks+clicks) > 200:
            clicks_random = random.sample(clicks+clicks, 200)#2*len(clicks))
        else: 
            clicks_random = random.sample(clicks+clicks, 2*len(clicks))
        myIm = drawLine(clicks_random, im, '0')
        s3_mask[:,:,slc] = myIm
        bar.next()
    bar.finish()

    return s3_mask

#%% func - selectCutS3sOptMxLoad- underDev
def selectCutS3sOptMxLoad_UD(filename, meshes_cut, dict_planes, resolution, dir_txtNnpy, save):

    if 'CJ' in filename:
        azimuth = -90
        cj = True
    else: 
        azimuth = 0
        cj = False
        
    s3_myoc_names_in = ['ch0_all','ch0_int', 'ch0_ext']
    s3_myoc_names_out = ['ch0_cut', 'ch0_cut_int', 'ch0_cut_ext']
    myoc_names = ['Myoc', 'Int.Myoc', 'Ext.Myoc']
    s3_endo_names_in = ['ch1_cut', 'ch1_int', 'ch1_cut_ext']
    s3_endo_names_out = ['ch1_cut', 'ch1_cut_int', 'ch1_cut_ext']
    endo_names = ['Endo', 'Int.Endo', 'Ext.Endo']

    cut_type = ['inflow', 'outflow']
    # cuts = []
    pls_normal = []
    pls_centre = []
    cuts_selected = []

    # Define what to cut from each layer
    myoc_cuts = []; endo_cuts = []
    # Define lists to save final meshes
    meshes_myoc = []; meshes_endo = []

    for cut in cut_type:
        labels = [mesh._legend for mesh in meshes_cut]
        text = filename+"\n\n >> Take a closer look at the -" +cut + "- of the meshes \n\tto decide which layer to cut\n >> Close the window when done"
        txt = Text2D(text, c="k", font= font)
        settings.legendSize = .15
        vp = Plotter(N=len(meshes_cut)+1, axes=4)
        for n, mesh in enumerate(meshes_cut):
            if n == 0:
                vp.show(mesh, txt, at=n)
            else: 
                vp.show(mesh, at=n)
        vp.show(meshes_cut, at=n+1, zoom=1, azimuth = azimuth, interactive=True)
        
        text_q_cuts = 'Select the layer from which you want to cut the -'+ cut + '- tract' 
        for i, m_name in enumerate(labels):
            text_q_cuts = text_q_cuts+'\n\t['+str(i)+']: '+ m_name
        text_q_cuts = text_q_cuts+'\n\t['+str(i+1)+']: both'+'\n\t['+str(i+2)+']: none? >>: '
        q_cuts = fcBasics.ask4input(text_q_cuts,int) 
        if q_cuts == len(labels):
            q_cuts = 'all'
        elif q_cuts == len(labels)+1:
            q_cuts = 'none'
        else: 
            q_cuts = str(q_cuts)
            
        if q_cuts != 'none':
            # q_cutsf will contain the index number(s) of the mesh(es) that need(s) to be cut
            q_cutsf = fcBasics.getInputNumbers(q_cuts, labels)
            meshesSel2cut = [meshes_cut[i] for i in q_cutsf]
            labelsSel2cut = [labels[i] for i in q_cutsf]

        if q_cuts != 'none':
            cuts_selected.append(cut)
            # Get plane to cut
            if len(meshesSel2cut) >= 2: 
                mesh_out = meshes_cut[labels.index('Myocardium')]
                mesh_in = meshes_cut[labels.index('Endocardium')]
            else: 
                mesh_out = ''
                mesh_in = meshes_cut[0]
                
            plane_cut, pl_cut_centre, pl_cut_normal = getPlane(filename = filename, type_cut = cut, info = '', mesh_in = mesh_in,
                                                                        mesh_out = mesh_out)
            # Reorient plane to images (s3)
            plane_im, pl_im_centre, pl_im_normal = rotatePlane2Images(pl_cut_centre, pl_cut_normal, type_cut = cut, cj = cj)
            pls_normal.append(pl_im_normal); pls_centre.append(pl_im_centre)
            #Save planes to dict
            dict_planes = addPlanes2Dict(planes = [plane_cut, plane_im], pls_centre = [pl_cut_centre ,pl_im_centre],
                                                    pls_normal = [pl_cut_normal, pl_im_normal], info = ['',''], dict_planes = dict_planes, print_txt = False)
            if 'Myocardium' in labelsSel2cut: #q_cuts == 0 or q_cuts == 2:
                myoc_cuts.append(cut)
                m_myoc = meshes_cut[labels.index('Myocardium')]
            if 'Endocardium' in labelsSel2cut: #q_cuts == 1 or q_cuts == 2:
                endo_cuts.append(cut)
                m_endo = meshes_cut[labels.index('Endocardium')]

    # Cut Myocardial layers
    if len(myoc_cuts) == 2:
        # Cut myocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - inf&outf (Myoc)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        for n, s3_name, myoc_name, s3_name_out in zip(count(), s3_myoc_names_in, myoc_names, s3_myoc_names_out):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            s3 = cutInfAndOutfOptMx(s3, pls_normal, pls_centre, resolution, '(Myoc)')
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_myoc, layer = myoc_name, plotshow = False)
            meshes_myoc.append(mesh_out)
            save_s3(filename = filename, s3 = s3, dir_txtNnpy = dir_txtNnpy, layer = s3_name_out)
            bar.next()
        bar.finish()

    elif len(myoc_cuts) == 1:
        index_myoc = cuts_selected.index(myoc_cuts[0])
        # Cut myocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - ' + myoc_cuts[0]+' (Myoc)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        for n, s3_name, myoc_name, s3_name_out in zip(count(), s3_myoc_names_in, myoc_names, s3_myoc_names_out):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            s3 = cutInfOrOutfOptMx(s3, pls_normal[index_myoc], pls_centre[index_myoc], resolution = resolution,
                                                                      option = myoc_cuts[0], mesh_name = '(Myoc)')
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_myoc, layer = myoc_name, plotshow = False)
            meshes_myoc.append(mesh_out)
            save_s3(filename = filename, s3 = s3, dir_txtNnpy = dir_txtNnpy, layer = s3_name_out)
            bar.next()
        bar.finish()
    else:
        print('- No cuts made to Myocardium!')
        for n, s3_name, myoc_name in zip(count(), s3_myoc_names_in, myoc_names):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_myoc, layer = myoc_name, plotshow = False)
            meshes_myoc.append(mesh_out)

    alert('whistle', 1)

    # Cut Endocardial layers
    if len(endo_cuts) == 2:
        # Cut endocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - inf&outf (Endo)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        for n, s3_name, endo_name, s3_name_out in zip(count(), s3_endo_names_in, endo_names, s3_endo_names_out):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            s3 = cutInfAndOutfOptMx(s3, pls_normal, pls_centre, resolution, '(Endo)')
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_endo, layer = endo_name, plotshow = False)
            meshes_endo.append(mesh_out)
            save_s3(filename = filename, s3 = s3, dir_txtNnpy = dir_txtNnpy, layer = s3_name_out)
            bar.next()
        bar.finish()

    elif len(endo_cuts) == 1:
        index_endo = cuts_selected.index(endo_cuts[0])
        # Cut myocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - ' + endo_cuts[0]+' (Endo)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        for n, s3_name, endo_name, s3_name_out in zip(count(), s3_endo_names_in, endo_names, s3_endo_names_out):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            s3 = cutInfOrOutfOptMx(s3, pls_normal[index_endo], pls_centre[index_endo], resolution = resolution,
                                                                      option = endo_cuts[0], mesh_name = '(Endo)')
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_endo, layer = endo_name, plotshow = False)
            meshes_endo.append(mesh_out)
            save_s3(filename = filename, s3 = s3, dir_txtNnpy = dir_txtNnpy, layer = s3_name_out)
            bar.next()
        bar.finish()
    else:
        print('- No cuts made to Endocardium!')
        for n, s3_name, endo_name in zip(count(), s3_endo_names_in, endo_names):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_endo, layer = endo_name, plotshow = False)
            meshes_endo.append(mesh_out)
    alert('jump', 1)

    meshes_cut = meshes_myoc+meshes_endo
    text= filename+"\n\n >> Resulting meshes"
    txt = Text2D(text, c="k", font= font)
    settings.legendSize = .3
    vp = Plotter(N=6, axes=10)
    for i, mesh in enumerate(meshes_cut):
        if i == 0:
            vp.show(mesh, txt, at = i, zoom = 1.2)
        elif i > 0 and i < 5:
            vp.show(mesh, at = i, zoom = 1.2)
        else:
            vp.show(mesh, at = i, zoom = 1.2, interactive = True)

    return meshes_cut, dict_planes

#%% func - selectCutS3sOptMxLoad
def selectCutS3sOptMxLoad(filename, m_endo, m_myoc, dict_planes, resolution, dir_txtNnpy, save):
    """
    Function used to cut inflow and/or outflow tract of the s3 masks (s3s2cut) given as input

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    m_endo : mesh
        Endocardial mesh. (vedo Mesh)
    m_myoc : mesh
        Myocardial mesh. (vedo Mesh)
    dict_planes : dictionary
        Initialised dictionary with planes information
    resolution : list of floats
        List with the x,y, z scaling values of the images taken. This information is taken from the metadata of the original file.
    dir_txtNnpy : path
        Path to the folder where the np arrays are saved.
    save : boolean
        True if you want to save the final masks, else False.

    Returns
    -------
    meshes_cut : list of meshes (vedo Meshes)
        List of all the surface reconstructions created with the cut masks
    dict_planes :  dictionary
        Resulting dictionary with planes information updated

    """

    if 'CJ' in filename:
        azimuth = -90
        cj = True
    else: 
        azimuth = 0
        cj = False
        
    s3_myoc_names_in = ['ch0_all','ch0_int', 'ch0_ext']
    s3_myoc_names_out = ['ch0_cut', 'ch0_cut_int', 'ch0_cut_ext']
    myoc_names = ['Myoc', 'Int.Myoc', 'Ext.Myoc']
    s3_endo_names_in = ['ch1_cut', 'ch1_int', 'ch1_cut_ext']
    s3_endo_names_out = ['ch1_cut', 'ch1_cut_int', 'ch1_cut_ext']
    endo_names = ['Endo', 'Int.Endo', 'Ext.Endo']

    cut_type = ['inflow', 'outflow']
    cuts = []
    pls_normal = []
    pls_centre = []
    cuts_selected = []

    # Define what to cut from each layer
    myoc_cuts = []; endo_cuts = []
    # Define lists to save final meshes
    meshes_myoc = []; meshes_endo = []

    for cut in cut_type:
        text = filename+"\n\n >> Take a closer look at the -" +cut + "- of both meshes \n\tto decide which layer to cut\n >> [0]:myoc/[1]:endo/[2]:both/[3]:none\n >> Close the window when done"
        txt = Text2D(text, c="k", font= font)
        settings.legendSize = .15
        vp = Plotter(N=3, axes=4)
        vp.show(m_myoc, txt, at=0, zoom=1)
        vp.show(m_endo, at=1, zoom=1)
        vp.show(m_myoc, m_endo, at=2, zoom=1, azimuth = azimuth, interactive=True)

        q_cuts = ask4input('Select the layer from which you want to cut the -'+ cut + '- tract \n  [0]:myoc/[1]:endo/[2]:both/[3]:none?: ',int)
        cuts.append(q_cuts)

        if q_cuts == 0 or q_cuts == 1 or q_cuts == 2:
            cuts_selected.append(cut)
            # Get plane to cut
            plane_cut, pl_cut_centre, pl_cut_normal = getPlane(filename = filename, type_cut = cut, info = '', mesh_in = m_endo,
                                                                        mesh_out = m_myoc)
            # Reorient plane to images (s3)
            plane_im, pl_im_centre, pl_im_normal = rotatePlane2Images(pl_cut_centre, pl_cut_normal, type_cut = cut, cj = cj)
            pls_normal.append(pl_im_normal); pls_centre.append(pl_im_centre)
            #Save planes to dict
            dict_planes = addPlanes2Dict(planes = [plane_cut, plane_im], pls_centre = [pl_cut_centre ,pl_im_centre],
                                                    pls_normal = [pl_cut_normal, pl_im_normal], info = ['',''], dict_planes = dict_planes, print_txt = False)
            if q_cuts == 0 or q_cuts == 2:
                myoc_cuts.append(cut)
            if q_cuts == 1 or q_cuts == 2:
                endo_cuts.append(cut)

    # Cut Myocardial layers
    if len(myoc_cuts) == 2:
        # Cut myocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - inf&outf (Myoc)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        for n, s3_name, myoc_name, s3_name_out in zip(count(), s3_myoc_names_in, myoc_names, s3_myoc_names_out):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            s3 = cutInfAndOutfOptMx(s3, pls_normal, pls_centre, resolution, '(Myoc)')
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_myoc, layer = myoc_name, plotshow = False)
            meshes_myoc.append(mesh_out)
            save_s3(filename = filename, s3 = s3, dir_txtNnpy = dir_txtNnpy, layer = s3_name_out)
            bar.next()
        bar.finish()

    elif len(myoc_cuts) == 1:
        index_myoc = cuts_selected.index(myoc_cuts[0])
        # Cut myocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - ' + myoc_cuts[0]+' (Myoc)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        for n, s3_name, myoc_name, s3_name_out in zip(count(), s3_myoc_names_in, myoc_names, s3_myoc_names_out):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            s3 = cutInfOrOutfOptMx(s3, pls_normal[index_myoc], pls_centre[index_myoc], resolution = resolution,
                                                                      option = myoc_cuts[0], mesh_name = '(Myoc)')
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_myoc, layer = myoc_name, plotshow = False)
            meshes_myoc.append(mesh_out)
            save_s3(filename = filename, s3 = s3, dir_txtNnpy = dir_txtNnpy, layer = s3_name_out)
            bar.next()
        bar.finish()
    else:
        print('- No cuts made to Myocardium!')
        for n, s3_name, myoc_name in zip(count(), s3_myoc_names_in, myoc_names):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_myoc, layer = myoc_name, plotshow = False)
            meshes_myoc.append(mesh_out)

    alert('whistle', 1)

    # Cut Endocardial layers
    if len(endo_cuts) == 2:
        # Cut endocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - inf&outf (Endo)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        for n, s3_name, endo_name, s3_name_out in zip(count(), s3_endo_names_in, endo_names, s3_endo_names_out):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            s3 = cutInfAndOutfOptMx(s3, pls_normal, pls_centre, resolution, '(Endo)')
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_endo, layer = endo_name, plotshow = False)
            meshes_endo.append(mesh_out)
            save_s3(filename = filename, s3 = s3, dir_txtNnpy = dir_txtNnpy, layer = s3_name_out)
            bar.next()
        bar.finish()

    elif len(endo_cuts) == 1:
        index_endo = cuts_selected.index(endo_cuts[0])
        # Cut myocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - ' + endo_cuts[0]+' (Endo)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        for n, s3_name, endo_name, s3_name_out in zip(count(), s3_endo_names_in, endo_names, s3_endo_names_out):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            s3 = cutInfOrOutfOptMx(s3, pls_normal[index_endo], pls_centre[index_endo], resolution = resolution,
                                                                      option = endo_cuts[0], mesh_name = '(Endo)')
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_endo, layer = endo_name, plotshow = False)
            meshes_endo.append(mesh_out)
            save_s3(filename = filename, s3 = s3, dir_txtNnpy = dir_txtNnpy, layer = s3_name_out)
            bar.next()
        bar.finish()
    else:
        print('- No cuts made to Endocardium!')
        for n, s3_name, endo_name in zip(count(), s3_endo_names_in, endo_names):
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
            mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                        mesh_original = m_endo, layer = endo_name, plotshow = False)
            meshes_endo.append(mesh_out)
    alert('jump', 1)

    meshes_cut = meshes_myoc+meshes_endo
    text= filename+"\n\n >> Resulting meshes"
    txt = Text2D(text, c="k", font= font)
    settings.legendSize = .3
    vp = Plotter(N=6, axes=10)
    for i, mesh in enumerate(meshes_cut):
        if i == 0:
            vp.show(mesh, txt, at = i, zoom = 1.2)
        elif i > 0 and i < 5:
            vp.show(mesh, at = i, zoom = 1.2)
        else:
            vp.show(mesh, at = i, zoom = 1.2, interactive = True)

    return meshes_cut, dict_planes

#%% func - cutInfAndOutfOptMx
def cutInfAndOutfOptMx(s3_cut, pls_normal, pls_centre, resolution, mesh_name):
    """
    Function used to cut inflow AND outflow tract of the s3 mask (s3_cut) given as input

    Parameters
    ----------
    s3_cut : numpy array
        Mask of the tissue layer that wants to be cut
    pls_normal : list of list of floats
        list of List with the x,y,z coordinatesof each of the planes' normal
    pls_centre : list of list of floats
        list of List with the x,y,z coordinates of each of the planes' centre
    resolution : list of floats
        List with the x,y, z scaling values of the images taken. This information is taken from the metadata of the original file.
    mesh_name : str
        Names of the mesh being cut

    Returns
    -------
    s3f_cut : numpy array
        Final mask of the cut tissue layer

    """

    # print('- Cutting s3 - inf&outf '+mesh_name)

    # Get dimensions of external stack
    xdim, ydim, zdim = s3_cut.shape
    # Reshape stacks as a vector
    s3_cut_v = s3_cut.reshape(-1)

    # Get vectors of x,y and z positions
    pix_coord_pos = np.where(s3_cut >= 0)
    del s3_cut
    # Trasform coordinate positions to um using resolution
    pix_um = np.transpose(np.asarray([pix_coord_pos[i]*resolution[i] for i in range(len(resolution))]))
    del pix_coord_pos

    normal_inf = unit_vector(pls_normal[0])
    normal_outf = unit_vector(pls_normal[1])

    # Find all the d values of pix_um
    d_pix_um_Inf = np.dot(np.subtract(pix_um,np.array(pls_centre[0])),np.array(normal_inf))
    d_pix_um_Outf = np.dot(np.subtract(pix_um,np.array(pls_centre[1])),np.array(normal_outf))
    del pix_um

    # Clear vector d_pix_um using only those that are 1 in stack
    d_pve_pix_um_Inf = s3_cut_v*d_pix_um_Inf
    d_pve_pix_um_Outf = s3_cut_v*d_pix_um_Outf
    del d_pix_um_Inf, d_pix_um_Outf

    # Duplicate s3f_v to initialise stacks without inflow
    s3f_all_v = np.copy(s3_cut_v)
    s3f_all_v.astype('uint8')
    del s3_cut_v

    # Find all positions in d_pve_pix_um that are at either side of the planes (outside of mesh)
    pos_outside_inf = np.where(d_pve_pix_um_Inf < 0)[0]
    pos_outside_outf = np.where(d_pve_pix_um_Outf > 0)[0]
    del d_pve_pix_um_Inf, d_pve_pix_um_Outf

    # Remove the points that are outside of the mesh (inflow)
    s3f_all_v[pos_outside_inf] = 0
    del pos_outside_inf

    # Remove the points that are outside of the mesh (ouflow)
    s3f_all_v[pos_outside_outf] = 0
    del pos_outside_outf

    # Reshape vector into matrix/stack
    s3f_cut = s3f_all_v.reshape((xdim, ydim, zdim))

    # alert('wohoo',1)

    return s3f_cut

#%% func - cutInfOrOutfOptMx
def cutInfOrOutfOptMx (s3_cut, pls_normal, pls_centre, resolution, option, mesh_name):
    """
    Function used to cut inflow OR outflow tract of the s3 mask (s3_cut) given as input

    Parameters
    ----------
    s3_cut : numpy array
        Mask of the tissue layer that wants to be cut
    pls_normal : list of floats
        List with the x,y,z coordinatesof each the plane's normal
    pls_centre : list of floats
        List with the x,y,z coordinatesof each the plane's centre
    resolution : list of floats
        List with the x,y, z scaling values of the images taken.
    option : str
        'Inflow'/'Outflow', depending on the type fo cut being made to the mask/mesh.
    mesh_name : str
        Names of the mesh being cut

    Returns
    -------
    s3f_cut : numpy array
        Final mask of the cut tissue layer

    """

    # print('- Cutting s3 - ' + option+' '+mesh_name)

    # Get dimensions of external stack
    xdim, ydim, zdim = s3_cut.shape
    # Reshape stacks as a vector
    s3_cut_v = s3_cut.reshape(-1)

    # Get vectors of x,y and z positions
    pix_coord_pos = np.where(s3_cut >= 0)
    del s3_cut
    # Trasform coordinate positions to um using resolution
    pix_um = np.transpose(np.asarray([pix_coord_pos[i]*resolution[i] for i in range(len(resolution))]))
    del pix_coord_pos

    normal  = unit_vector(pls_normal)
    # Find all the d values of pix_um
    d_pix_um = np.dot(np.subtract(pix_um,np.array(pls_centre)),np.array(normal))

    # Clear vector d_pix_um using only those that are 1 in stack
    d_pve_pix_um = s3_cut_v*d_pix_um
    del pix_um

    # Duplicate s3f_v to initialise stacks without inflow/outflow
    s3f_all_v = np.copy(s3_cut_v)
    s3f_all_v.astype('uint8')
    del s3_cut_v

    # Find all positions in d_pve_pix_um that are at either side of the planes (outside of mesh)
    if option == 'inflow':
        pos_outside = np.where(d_pve_pix_um < 0)[0]
    elif option == 'outflow':
        pos_outside = np.where(d_pve_pix_um > 0)[0]
    del d_pve_pix_um

    # Remove the points that are outside of the mesh (inflow/outflow)
    s3f_all_v[pos_outside] = 0
    del pos_outside

    # Reshape vector into matrix/stack
    s3f_cut = s3f_all_v.reshape((xdim, ydim, zdim))
    del s3f_all_v

    # alert('wohoo',1)

    return s3f_cut

#%% - CREATE/MODIFY OBJECTS
#%% >>> MESHES
#%% func - createAll3LayerMeshes
def createAll3LayerMeshes(filename, s3_all, s3_in, s3_out, resolution, layer):
    """
    Function that creates the 3 meshes of a heart layer (int, ext, all), using the masks given as input

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    s3_all : numpy array
        Mask of the tissue layer
    s3_in : numpy array
        Mask of the internal contours of the tissue layer
    s3_out : numpy array
        Mask of the external contours of the tissue layer
    resolution : list of floats
        List with the x,y, z scaling values of the images taken. This information is taken from the metadata of the original file.
    layer : str
        Name of the heart layer being reconstructed

    Returns
    -------
    mesh_all : mesh
        Mesh of the heart layer. (vedo Mesh)
    mesh_in : mesh
        Mesh of the filled internal contours of the heart layer. (vedo Mesh)
    mesh_out : mesh
        Mesh of the filled external contours of the heart layer. (vedo Mesh)

    """
    if 'CJ' in filename:
        azimuth = -90
    else: 
        azimuth = 0

    print('- Creating surface reconstruction of '+ layer +' ')

    # Extract vertices, faces, normals and values of each mesh
    verts_all, faces_all, _, _ = measure.marching_cubes_lewiner(s3_all, spacing=resolution)
    verts_in, faces_in, _, _ = measure.marching_cubes_lewiner(s3_in, spacing=resolution)
    verts_out, faces_out, _, _ = measure.marching_cubes_lewiner(s3_out, spacing=resolution)
    alert('frog',1)

    # Create meshes
    mesh_all = Mesh([verts_all, faces_all])
    mesh_in = Mesh([verts_in, faces_in])
    mesh_out = Mesh([verts_out, faces_out])
    if not 'CJ' in filename: 
        mesh_all.rotateZ(-90).wireframe(True)
        mesh_in.rotateZ(-90).wireframe(True)
        mesh_out.rotateZ(-90).wireframe(True)
    else: 
        mesh_all.rotateX(-90).wireframe(True)
        mesh_in.rotateX(-90).wireframe(True)
        mesh_out.rotateX(-90).wireframe(True)
        
    alert('clown',1)

    if layer == 'Myoc' or layer == 'Endoc':
        mesh_all = mesh_all.extractLargestRegion()
        mesh_in = mesh_in.extractLargestRegion()
        mesh_out = mesh_out.extractLargestRegion()

    if layer == 'CJ':
        mesh_all.color("darkorange").alpha(1).wireframe().legend(layer)
        mesh_in.color("gold").alpha(0.05).wireframe().legend('Int.'+layer)
        mesh_out.color("limegreen").alpha(1).wireframe().legend('Ext.'+layer)
    elif layer == 'Myoc':
        mesh_all.color("darkcyan").alpha(0.05).wireframe().legend(layer)
        mesh_in.color("gold").alpha(0.05).wireframe().legend('Int.'+layer)
        mesh_out.color("limegreen").alpha(1).wireframe().legend('Ext.'+layer)
    elif layer == 'Endo':  #Endocardium
        mesh_all.color("darkmagenta").alpha(0.05).wireframe().legend(layer)
        mesh_in.color("gold").alpha(0.05).wireframe().legend('Int.'+layer)
        mesh_out.color("limegreen").alpha(1).wireframe().legend('Ext.'+layer)
    else:  #Endocardium
        # print('AJA')
        mesh_all.color("cornflowerblue").alpha(1).wireframe().legend(layer)
        mesh_in.color("indigo").alpha(1).wireframe().legend('Int.'+layer)
        mesh_out.color("turquoise").alpha(1).wireframe().legend('Ext.'+layer)

    text = filename+"\n\n >> "+layer
    txt = Text2D(text, c=c, font=font)
    settings.legendSize = .3
    if layer in ['CJ', 'Myoc', 'Endo']:
        vp = Plotter(N=3, axes=13)
        vp.show(mesh_all, txt, at=0, zoom=1.2)
        vp.show(mesh_all, mesh_out, at=1, zoom=1.2)
        vp.show(mesh_all, mesh_in, at=2, zoom=1.2, azimuth = azimuth, interactive=True)
    else: 
        # print('AJA')
        vp = Plotter(N=3, axes=13)
        vp.show(mesh_in, txt, at=0, zoom=1.2)
        vp.show(mesh_out, at=1, zoom=1.2)
        vp.show(mesh_all, at=2, zoom=1.2, azimuth = azimuth, interactive=True)

    return mesh_all, mesh_in, mesh_out

#%% func - createExtLayerMesh
def createExtLayerMesh(filename, s3_ext, resolution, layer, info, extractLargest = True, plotshow = True):
    """
    Function that creates the external mesh of a heart layer, using the mask given as input

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    s3_ext : numpy array
        Mask of the external contours of the tissue layer
    resolution : list of floats
        List with the x,y, z scaling values of the images taken. This information is taken from the metadata of the original file.
    layer : str
        Name of the heart layer being reconstructed
    info : str
        Text with additional information.
    extractLargest : boolean, optional
        True if you only want to keep the largest piece of the reconstructed mesh, else False. The default is True.
    plotshow : boolean, optional
        True if you want to see the resulting mesh in a plot, else False. The default is True.

    Returns
    -------
    mesh_ext : mesh
        Mesh of the filled external contours of the heart layer. (vedo Mesh)

    """

    print('- Creating surface reconstruction of '+ layer +' ('+info+')')

    # Extract vertices, faces, normals and values of each mesh
    verts_ext, faces_ext, _, _ = measure.marching_cubes_lewiner(s3_ext, spacing=resolution)
    alert('frog',1)

    mesh_ext = Mesh([verts_ext, faces_ext])
    if 'CJ' not in filename: 
        mesh_ext.rotateZ(-90).wireframe(True)
    else: 
        mesh_ext.rotateX(-90).wireframe(True)
    alert('clown',1)

    layers_cut = ['Myoc', 'Endo']

    if layer in layers_cut and extractLargest:
        mesh_ext = mesh_ext.extractLargestRegion()

    if layer == 'CJ':
        mesh_ext.color("darkorange").alpha(1).wireframe().legend('Ext.'+layer)
    elif layer == 'Myoc':
        mesh_ext.color("darkcyan").alpha(1).wireframe().legend('Ext.'+layer)
    else:  #'Endo'
        mesh_ext.color("darkmagenta").alpha(1).wireframe().legend('Ext.'+layer)

    if plotshow:
        text = filename+"\n\n >> External "+layer+' ('+info+')'
        txt = Text2D(text, c=c, font=font)
        settings.legendSize = .3
        vp = Plotter(N=1, axes=13)
        vp.show(mesh_ext, txt, at=0, zoom=1.2, interactive=True)

    return mesh_ext

#%% func - createLayerMesh
def createLayerMesh(filename, s3, resolution, layer, name, colour, alpha, plotshow = False):
    """
    Function that creates the mesh of a heart layer, using the masks given as input

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    s3 : numpy array
        Mask of the tissue layer
    resolution : list of floats
        List with the x,y, z scaling values of the images taken. This information is taken from the metadata of the original file.
    layer : str
        Name of the heart layer being reconstructed
    name : str
        Name of the heart layer being reconstructed
    colour : str
        Colour
    alpha : float
        Opacity value.
    plotshow : boolean, optional
        True if you want to see the resulting mesh in a plot, else False. The default is False.

    Returns
    -------
    mesh : mesh
        Mesh of the heart layer. (vedo Mesh)

    """

    print('- Creating surface reconstructions of '+ name )

    # Extract vertices, faces, normals and values of each mesh
    verts, faces, _, _ = measure.marching_cubes_lewiner(s3, spacing=resolution)
    alert('frog',1)

    # Create meshes
    mesh = Mesh([verts, faces])
    if layer == 'Myoc' or layer == 'Endo':
        mesh = mesh.extractLargestRegion()

    if 'CJ' not in filename: 
        mesh.rotateZ(-90).color(colour).alpha(alpha).legend(name).wireframe()
    else: 
        mesh.rotateX(-90).color(colour).alpha(alpha).legend(name).wireframe()
    alert('clown',1)

    if plotshow:
        text = filename+"\n\n >> "+layer
        txt = Text2D(text, c=c, font=font)
        settings.legendSize = .3
        vp = Plotter(N=1, axes=13)
        vp.show(mesh, txt, at=0, zoom=1.2, interactive=True)

    return mesh

#%% func - getCutMesh
def getCutMesh(filename, s3_cut, resolution, mesh_original, layer, plotshow):
    """
    Function to create surface reconstruction of the cut mask (s3_cut) given as input.

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    s3_cut : numpy array
        Mask of the tissue layer.
    resolution : list of floats
        List with the x,y, z scaling values of the images taken. This information is taken from the metadata of the original file.
    mesh_original : mesh
        Original mesh of the heart layer. (vedo Mesh).
    layer :  str
        Name of the heart layer being reconstructed.
    plotshow : boolean
        True if you want to see the resulting mesh in a plot, else False.

    Returns
    -------
    mesh_cut : mesh
        Resulting mesh. (vedo Mesh).

    """
    if plotshow:
        print('- Creating surface reconstructions of cut '+ layer)

    verts_cut, faces_cut, _, _ = measure.marching_cubes_lewiner(s3_cut, spacing=resolution)
    mesh_cut = Mesh([verts_cut, faces_cut])
    if layer not in ['Int.Endo', 'CJ','CardiacJelly']: #!= 'Int.Endo' or layer != 'CJ' or layer != 'CardiacJelly':
        # if 'CJ' not in filename: 
        mesh_cut = mesh_cut.extractLargestRegion()
        # print('aja')
    # else: 
    #     print('aja2')
    if plotshow:
        alert('wohoo',1)

    if layer == 'Cardiac Jelly' or layer == 'CJ' or layer == 'CardiacJelly':
        color = 'darkorange'
    elif layer == 'Myocardium' or layer == 'Myoc':
        color = 'darkcyan'
    elif layer == 'Endocardium' or layer == 'Endo':
        color = 'darkmagenta'
    else:
        color = 'indigo'

    if 'CJ' not in filename: 
        mesh_cut.rotateZ(-90).color(color).legend('Cut ('+layer+')')
    else: 
        mesh_cut.rotateX(-90).color(color).legend('Cut ('+layer+')')

    text = filename
    txt = Text2D(text, c=c, font=font)

    if plotshow:
        mesh_original.legend(layer).alpha(0.1).color('tomato')
        settings.legendSize = .3
        vp = Plotter(N=1, axes=4)
        vp.show(mesh_original, mesh_cut, txt, at=0, zoom=1.2, interactive=True)

    return mesh_cut

#%% func - recreateCutMesh
def recreateCutMesh(filename, name, resolution, dir_txtNnpy, dict_colour):
    
    names_out_if_cut = ['ch0_cut', 'ch0_cut_int', 'ch0_cut_ext', 'ch1_cut', 'ch1_cut_int', 'ch1_cut_ext', 
                        'cj', 'cj_int', 'ch0_cut_int']
    names_out_no_cut = ['ch0_all','ch0_int', 'ch0_ext','ch1_cut', 'ch1_int', 'ch1_cut_ext', 
                        'cj', 'cj_int', 'ch0_cut_int']
    names = ['ch0_all','ch0_int', 'ch0_ext','ch1_all', 'ch1_int', 'ch1_ext', 
                 'cj', 'cj_in', 'cj_out']
    m_names = ['myoc','myoc_int','myoc_ext','endo','endo_int','endo_ext',
               'cj', 'cj_in', 'cj_out']
    leg_names = ['Myocardium','Int.Myoc','Ext.Myoc','Endocardium','Int.Endo','Ext.Endo',
               'CardiacJelly', 'Int.CJ', 'Ext.CJ']
    
    index = names.index(name)
    try:
        s3_name = names_out_if_cut[index]
        [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
    except: 
        s3_name = names_out_no_cut[index]
        [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
    # print(leg_names[index])
    mesh_out = getCutMesh(filename = filename, s3_cut = s3, resolution = resolution,
                                mesh_original = '', layer = leg_names[index], plotshow = False)
    alert('wohoo', 1)
    mesh_out.color(dict_colour[m_names[index]]['colour']).legend(leg_names[index]).alpha(0.1)
    
    vp = Plotter(N=1, axes = 10)
    vp.show(mesh_out, at=0, interactive = True)
    
    return mesh_out
            
#%% func - createMeshes4CL
def createMeshes4CL(filename, meshes, plotshow):
    """
    Function that cleans and smooths meshes given as input to get centreline using VMTK

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    meshes : list of meshes
        List of meshes to .
    plotshow : boolean
        True if you want to see the resulting mesh in a plot, else False.

    Returns
    -------
    meshes4cl_out : list of meshes
        List of cleaned and smoothed meshes. (vedo Mesh)
    meshes4cl_names : list of str
        List of names with which to save exported meshes to obtain centreline

    """

    if 'CJ' in filename:
        azimuth = -90
    else: 
        azimuth = 0
        
    meshes_selected = []
    meshes4cl = []
    meshes4cl_names = []
    mesh_colors = ['springgreen', 'violet']
    names = ['Int.Myoc(Cut)', 'Ext.Endo(Cut)']
    names_export = ['myoc_int_cut4cl','endo_ext_cut4cl']

    q_selectMeshes = ask4input('Select the meshes you would like to export to extract centreline \n\t\t[0]: Int.Myocardium/[1]: Ext.Endocardium/[2]: both: ', int)
    if q_selectMeshes in [0,2]:
        meshes_selected.append(meshes[0])
        meshes4cl_names.append(names_export[0])
    if q_selectMeshes in [1,2]:
        meshes_selected.append(meshes[1])
        meshes4cl_names.append(names_export[1])

    for i, mesh_cl in enumerate(meshes_selected):
        print("- File being processed for centreline: ", filename +' - '+names[i])
        mesh_cl.legend(names[i])

        if plotshow:
            settings.legendSize = .3
            vp = Plotter(N=2, axes=4)
            text1 = filename+"\n- Original mesh - "+names[i]
            txt1 = Text2D(text1, c=c, font=font)
            vp.show(mesh_cl.clone(),txt1, at=0)

        print('- Cleaning mesh '+names[i])
        print("\t- Original number of points making up mesh: ",mesh_cl.NPoints())
        # mesh_cl = meshes_cl[i]
        # Reduce the number of points that make up the mesh
        mesh_cl.clean(tol=0.005)
        print("\t- Number of points after cleaning surface: ",mesh_cl.NPoints(),'\n- Smoothing mesh...', names[i])
        # Smooth mesh
        meshCL_cut = mesh_cl.clone().smoothMLS2D(f=0.5)
        meshCL_cut.color(mesh_colors[i]).legend(names[i]+"-C&S")
        print('- Mesh smoothed!');alert('wohoo',1)

        if plotshow:
            text2 = " \n- Cleaned and smoothed mesh"
            txt2 = Text2D(text2, c=c, font=font)
            vp.show(meshCL_cut, txt2, at=1, azimuth = azimuth, interactive=True)

        meshes4cl.append(meshCL_cut)

    meshes4cl_out = meshes4cl[-len(meshes_selected):]

    settings.legendSize = .3
    vp = Plotter(N=len(meshes4cl_out), axes=13)
    for i in range(len(meshes4cl_out)):
        if i != len(meshes4cl_out)-1:
            vp.show(meshes4cl[i], at=i, zoom=1)
        else:
            vp.show(meshes4cl[i], at=i, zoom=1, azimuth = azimuth, interactive = True)

    return meshes4cl_out, meshes4cl_names

#%% func - cutMeshes4CL
def cutMeshes4CL(filename, meshes, cuts, cut_direction, dicts, plotshow):
    """
    Funtion that cuts the inflow and outflow tract of meshes from which the centreline will be obtained.

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    meshes : list of meshes
        List of cleaned and smoothed meshes. (vedo Mesh)
    cuts : list of str
        List with type of cuts ['inflow', 'outflow']
    cut_direction : list of booleans
        list of booleans indicating the cut direction depending on the 'cuts'
    dicts : list of dicts
        [dict_planes, dict_pts, dict_kspl].
    plotshow : boolean
        True if you want to see the resulting mesh in a plot, else False.

    Returns
    -------
    meshes[-len(meshes):] : list of meshes
        List of cut meshes with same length as input meshes (vedo Mesh)
    dicts_f :  list of dicts
        List of updated dictionaries

    """
    if 'CJ' in filename:
        azimuth = -90
    else: 
        azimuth = 0
    
    mark_colors = ['deepskyblue', 'tomato']*3
    mesh_colors = ['springgreen', 'violet']*3

    dict_planes, dict_pts, dict_kspl = dicts
    ksplines = []; spheres = []
    num_meshes_out = len(meshes)

    for n, cut in enumerate(cuts):
        print('- Cutting: '+cut)
        #Get plane to cut
        if isinstance(meshes, list):
            if len(meshes) == 2:
                planeCL_cut, plCL_cut_centre, plCL_cut_normal = getPlane(filename = filename, type_cut = cut,
                                                                     info = '4CL', mesh_in = meshes[1],
                                                                     mesh_out = meshes[0], dict_planes = dict_planes)
            else:
                planeCL_cut, plCL_cut_centre, plCL_cut_normal = getPlane(filename = filename, type_cut = cut,
                                                                     info = '4CL', mesh_in = meshes[0],
                                                                     mesh_out = '', dict_planes = dict_planes)
        else:
            planeCL_cut, plCL_cut_centre, plCL_cut_normal = getPlane(filename = filename, type_cut = cut,
                                                                 info = '4CL', mesh_in = meshes,
                                                                 mesh_out = '', dict_planes = dict_planes)
        dict_planes = addPlanes2Dict(planes = [planeCL_cut], pls_centre = [plCL_cut_centre],
                                            pls_normal = [plCL_cut_normal], info = [''], dict_planes = dict_planes)

        for i, mesh4cl in enumerate(meshes[-num_meshes_out:]):
            # mesh_name = mesh4cl._legend
            #print(mesh_name)
            pts2cut, _ = getPointsAtPlane(points = mesh4cl.points(), pl_normal = plCL_cut_normal,
                                       pl_centre = plCL_cut_centre)
            ordpts, _ = order_pts(points = pts2cut)
            #print('ordpts', ordpts)
            # Create spline around cut
            kspl = KSpline(ordpts, continuity=0, tension=0, bias=0, closed=True)
            kspl.color(mark_colors[i]).legend('ksplCut4CL_'+cut).lw(2)
            ksplines.append(kspl)
            # Get centroid of kspline to add to the centreline
            kspl_bounds = kspl.bounds()
            pt_centroid = np.mean(np.asarray(kspl_bounds).reshape((3, 2)),axis=1)
            sph_centroid = Sphere(pos=pt_centroid, r=2, c=mark_colors[i]).legend('sph_Cut4CL_'+cut)
            spheres.append(sph_centroid)

            dict_pts = addPoints2Dict(spheres = [sph_centroid], info = [mesh4cl._legend[:-4]], dict_pts = dict_pts)
            dict_kspl = addKSplines2Dict(kspls = [kspl], info = mesh4cl._legend[:-4], dict_kspl = dict_kspl)

            # Cutmesh using created plane
            mesh4cl_new = mesh4cl.clone().cutWithMesh(planeCL_cut, invert=cut_direction[n])
            mesh4cl_new = mesh4cl_new.extractLargestRegion()
            mesh4cl_new.color(mesh_colors[i]).alpha(0.05).wireframe(True).legend(mesh4cl._legend)

            settings.legendSize = .3
            vp = Plotter(N=1)
            text = "- Resulting mesh after cutting"
            txt = Text2D(text, c=c, font=font)
            vp.show(mesh4cl_new, kspl, sph_centroid, txt, at=0, viewup="y", azimuth = azimuth, interactive=True)

            meshes.append(mesh4cl_new)

    if plotshow:
        settings.legendSize = .3
        vp = Plotter(N=1, axes=13)
        text = filename+"\n\n >> Resulting mesh after cutting inflow & outflow tract"
        txt = Text2D(text, c=c, font=font)
        vp.show(meshes[-num_meshes_out:], ksplines, spheres, txt, at=0, azimuth = azimuth, interactive=True)

    dicts_f = dict_planes, dict_pts, dict_kspl

    return meshes[-num_meshes_out:], dicts_f

#%% func - divideMeshesLnR
def divideMeshesLnR(filename, meshes, cl_ribbon, file_num, df_res, colors =  ['skyblue','darkblue']):
    """
    Function that divides meshes into Left and Right using the extended centreline (ribbon)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    meshes : list of meshes
        List of meshes to divide. (vedo Mesh)
    cl_ribbon : ribbon
        Dorso-ventral extended centreline (vedo Ribbon)

    Returns
    -------
    meshes_LnR : list of list of meshes
        List of list of the divided meshes [left, right] (vedo Mesh)

    """
    
    spaw_analysis = False
    if 'spaw' in df_res.loc[file_num,'spAnalysis']:
        spaw_analysis = True
    
    meshes_LnR = []
    for i, mesh in enumerate(meshes):

        meshes2cutLR = mesh
        mesh_legend = mesh._legend
        if not spaw_analysis: 
            print('\n- Dividing '+mesh_legend+' into left and right sides')
        else: 
            print('\n- Dividing '+mesh_legend+' into dorsal and ventral sides')
            
        mesh_1 = meshes2cutLR.clone().cutWithMesh(cl_ribbon, invert=True)
        mesh_1.alpha(0.05).wireframe(True).color(colors[0])

        mesh_2 = meshes2cutLR.clone().cutWithMesh(cl_ribbon, invert=False)
        mesh_2.alpha(0.05).wireframe(True).color(colors[1])
        alert("wohoo",1)

        text = filename+"\n\n >> Resulting mesh after cutting with centreline \n >> Before closing the window make sure you confirm \n\t"
        if not spaw_analysis: 
            mesh_1.legend(mesh_legend+'-Left')
            mesh_2.legend(mesh_legend+'-Right')
            txt2 = colors[0]+": Left, AND "+colors[1]+": Right"
        else: 
            mesh_1.legend(mesh_legend+'-Dorsal')
            mesh_2.legend(mesh_legend+'-Ventral')
            txt2 = colors[0]+": Dorsal, AND "+colors[1]+": Ventral"

        txt = Text2D(text, c=c, font=font)
        settings.legendSize = .3
        vp = Plotter(N=3, axes=13)
        vp.show(meshes2cutLR, cl_ribbon, txt, at=0)
        vp.show(mesh_1, at=1)
        vp.show(mesh_2, at=2, zoom = 1.2, interactive=True)
        
        # if not spaw_analysis: 
        q_happy = ask4input('Are meshes classified correctly -'+ txt2 +'- [0]:no/[1]:yes: ', bool)
        # else: 
        #     q_happy = ask4input('Are meshes classified correctly -A(light blue):dorsal, B(dark blue):ventral? [0]:no/[1]:yes: ', bool)
            
        if not q_happy:
            leftMesh = ask4input('Select the mesh number that corresponds to the left side. \n\t[1]: the one in the middle plot \n\t[2]: the one in the right plot', int)
            if leftMesh == 1:
                if not spaw_analysis: 
                    mesh_1.legend(mesh_legend+'-Left')
                    mesh_2.legend(mesh_legend+'-Right')
                else: 
                    mesh_1.legend(mesh_legend+'-Dorsal')
                    mesh_2.legend(mesh_legend+'-Ventral')
                mesh_LnR = [mesh_1, mesh_2]
            elif leftMesh == 2:
                if not spaw_analysis:
                    mesh_2.legend(mesh_legend+'-Left')
                    mesh_1.legend(mesh_legend+'-Right')
                else: 
                    mesh_2.legend(mesh_legend+'-Dorsal')
                    mesh_1.legend(mesh_legend+'-Ventral')
                mesh_LnR = [mesh_2, mesh_1]

        else:
            mesh_LnR = [mesh_1, mesh_2]

        meshes_LnR.append(mesh_LnR)

    return meshes_LnR

#%% func - getRing2CutChambers
def getRing2CutChambers(filename, kspl_CL, mesh2cut, resolution, dir_stl, dir_txtNnpy, dict_pts, dict_shapes):
    """
    Function to define the ring needed to cut the meshes into atrium and ventricle

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    kspl_CL : Kspline
        Centreline (vedo KSpline)
    mesh2cut : Mesh
        Mesh (vedo Mesh)
    resolution : list of floats
        List with the x,y, z scaling values of the images taken. This information is taken from the metadata of the original file.
    dir_stl : path
        Path to the folder where the meshes are saved.
    dir_txtNnpy : path
        Path to the folder where the np arrays are saved.
    dict_pts : dictionary
        Initialised dictionary with points information
    dict_shapes : dictionary
        Initialised dictionary with shapes information

    Returns
    -------
    cyl_final : Cylinder
        Cylinder/ring defined to cut meshes into chambers
    num_pt : int
        Index of the centreline point closer to the plane that cuts the meshes into the two chambers.
    atr_meshes: list of meshes
        list of atrial meshes (vedo Meshes)
    vent_meshes: list of meshes
        list of ventricular meshes (vedo Meshes)
    dict_shapes : dictionary
        Resulting dictionary with shapes information updated
    dict_pts : dictionary
        Resulting dictionary with points information updated
    s3_cyl : numpy array of booleans
        Array with disc mask to cut chambers. The size of this array corresponds to the size of the original stack. 

    """
    dORv = filename[9:10]
    if dORv == 'D':
        azimuth = -90
    else:
        azimuth = 0
    
    # Plot spheres through centreline inside myocardium
    spheres_spl = sphInSpline(kspl_CL = kspl_CL, colour = True)
    settings.legendSize = .3
    vp = Plotter(N=1, axes = 7)
    text = filename+"\n\n >> Define the centreline point number to use to initialise \n  disc to divide heart into chambers \n  [NOTE: Spheres appear in centreline every 10 points, starting from \n  outflow (blue) to inflow (red) tract]"
    txt = Text2D(text, c="k", font= font)
    vp.show(mesh2cut.alpha(0.01), kspl_CL, spheres_spl, txt, at=0, azimuth = azimuth, interactive=True)
    
    # Get myocIntBall data
    [[m_myocIntBall], [myoc_intBall]] = openThicknessMeshes(filename = filename, meshes_names = ['myoc_intBall'], extension = 'vtk',
                                      dir_stl = dir_stl, dir_txtNnpy = dir_txtNnpy, print_txt = False)
        
    # Get disc position and orientation to cut heart layers
    mesh2cut.alpha(0.05)
    happyWithMyocCut = False
    happyWithDisc = False
    while not happyWithMyocCut:
        while not happyWithDisc:
            # Create plane
            num_pt = ask4input('Enter the centreline point number you want to use to initialise the disc to divide the heart into chambers: ', int)
            # Use flat disc or centreline orientation at point to define disc?
            centrelineORFlat = ask4input('Do you want to use the centreline orientation at the selected point to define disc orientation or \n  initialise the disc in a plane normal to the y-z plane? \n\t[0]: Use centreline orientation at the selected point\n\t[1]: Use plane perpendicular to the y-z plane >>>: ', bool)
            if not centrelineORFlat: 
                pl_Ch_normal, pl_Ch_centre = getPlaneNormal2Pt (pt_num = num_pt, spline_pts = kspl_CL.points())
            else: 
                pl_Ch_centre = kspl_CL.points()[num_pt]
                pl_Ch_normal = [1,0,0]
            # print('- Modifying disc position to cut chambers. Initially defined disc radius is 60um. \n  Once you have selected the position of the disc a new radius will be calculated based on the mesh points the disc cuts. \n  If you are not happy with the disc radius, you will be able to modify it just before proceeding to the cut.')
            # Modify (rotate and move cylinder/disc)
            cyl_test, sph_test, rotX, rotY, rotZ = modifyDisc (filename = filename,
                                                    pl_normal = pl_Ch_normal, pl_centre = pl_Ch_centre, 
                                                    radius = 60, type_cut = 'Chamber',
                                                    mesh1 = mesh2cut, xyz_bounds = mesh2cut.bounds())
            # print('pl_Ch_centre', pl_Ch_centre)
            # print('sph_test.pos()',sph_test.pos())
            # print('cyl_test.pos()',cyl_test.pos())
            # print('num_pt:', num_pt)
            
            # Get new normal of rotated disc
            pl_Ch_normal_corrected = newNormal3DRot(normal = pl_Ch_normal, rotX = rotX, rotY = rotY, rotZ = rotZ)
            normal_unit = unit_vector(pl_Ch_normal_corrected)*10
            # Get central point of newly defined disc
            pl_Ch_centre = sph_test.pos()
            
            # Get points at plane
            # Myocardium
            pts2cut, _ = getPointsAtPlane(points = mesh2cut.points(), pl_normal = pl_Ch_normal_corrected,
                                                pl_centre = pl_Ch_centre, tol = 1)
            # Internal myocardium with ballooning data
            _, data2cut = getPointsAtPlane(points = m_myocIntBall.points(), pl_normal = pl_Ch_normal_corrected,
                                                pl_centre = pl_Ch_centre, tol = 1, addData = myoc_intBall)
            plane_Ch = Plane(pos = pl_Ch_centre, normal = pl_Ch_normal, sx = 300)
            # Cut cl with plane
            ksplCL_cut = kspl_CL.clone().cutWithMesh(plane_Ch, invert=True)
            # ksplCL_cut.lw(5).color('tomato')
            # print(ksplCL_cut.points()[0], ksplCL_cut.points()[-1])
            # Find point of centreline closer to last point of kspline cut
            ksplCL_cutPt, num_pt = findClosestPtGuess(ksplCL_cut.points(), kspl_CL.points(), index_guess = num_pt)#findClosestPt(ksplCL_cut.points()[-1], kspl_CL.points())
            # print('num_pt:', num_pt)
            
            # Newly defined centreline point cut by plane/disc
            cl_point = pl_Ch_centre# kspl_CL.points()[num_pt]
            sph_cut = Sphere(pos = cl_point, r=4, c='gold').legend('sph_ChamberCut')
        
            # Order the points
            ordpts, _ = order_pts(points = pts2cut)
            # Find distance between points at plane and cl_point to define radius
            dist = []
            for pt in ordpts:
                dist.append(findDist(cl_point, pt))
            
            r_circle_max = np.mean(dist)*1.2
            r_circle_min = min(dist)*0.8
            if r_circle_max < max(data2cut):
                r_circle_max = max(data2cut)*1.5
            
            # Build new disc to confirm
            cyl_final = Cylinder(pos = pl_Ch_centre,r = r_circle_max, height = 2*0.225, axis = normal_unit, c = 'purple', cap = True, res = 300)
            r_circle_max_str = r_circle_max
            
            text = filename+"\n\n >> Check the position and the radius of the disc to cut the heart into chambers.\n  Make sure it is cutting the heart through the AVC and hopefully not cutting any other chamber regions. \n >> Close the window when done"
            txt = Text2D(text, c=c, font=font)
            settings.legendSize = .3
            vp = Plotter(N=1, axes=4)
            vp.show(mesh2cut, cyl_final, ksplCL_cut, txt, at=0, viewup="y", azimuth=0, elevation=0, interactive=True)
            happy = ask4input('Are you happy with the position of the disc [radius: '+format(r_circle_max_str,'.2f')+"um] to cut heart into chambers? \n  [0]: no, I would like to define a new position for the disc\n  [1]: yes, but I would like to redefine the disc radius \n  [2]: yes, I am happy with both, disc position and radius :", int)
            if happy == 1:
                happy_rad = False
                while not happy_rad:
                    r_circle_max = ask4input('Input disc radius [um]: ', float)
                    text = filename+"\n\n >> New radius \n  Check the radius of the disc to cut the myocardial tissue into \n   chambers. Make sure it is cutting through the AVC and not catching any other chamber regions. \n >> Close the window when done"
                    cyl_final = Cylinder(pos = pl_Ch_centre,r = r_circle_max, height = 2*0.225, axis = normal_unit, c = 'purple', cap = True, res = 300)
                    r_circle_max_str = r_circle_max
                    txt = Text2D(text, c="k", font= font)
                    settings.legendSize = .15
                    vp = Plotter(N=1, axes = 10)
                    vp.show(mesh2cut.alpha(1), kspl_CL, cyl_final, sph_cut, txt, at = 0, interactive=True)
                    r_circle_max_str = r_circle_max
                    happy_rad = ask4input('Is the selected radius ['+format(r_circle_max_str,'.2f')+"um] sufficient to cut heart into chambers? \n  [0]: no, I would like to change its value \n  [1]: yes, it cuts the heart without disrupting too much the chambers!: ", bool)
                happyWithDisc = True
            elif happy == 2:
                happyWithDisc = True
            
        cyl_final.legend('cyl2CutChambers_o')
        cyl_data = [r_circle_max, r_circle_min, normal_unit, pl_Ch_centre]
        dict_shapes = addShapes2Dict (shapes = [cyl_final], dict_shapes = dict_shapes, radius = [cyl_data], print_txt = False)
        
        #Define atrium and ventricle
        try: 
            avc_minus50y = kspl_CL.points()[num_pt-50]
        except: 
            avc_minus50y = kspl_CL.points()[0]
            print('-> First centreline point got selected as AVC-50')
        try: 
            avc_plus50y = kspl_CL.points()[num_pt+50]
        except: 
            avc_plus50y = kspl_CL.points()[-1]
            print('-> Last centreline point got selected as  AVC+50')
            
        sph_atr = Sphere(pos = avc_plus50y, r=4, c='darkorange').legend('sph_CentreOfAtrium')
        sph_vent = Sphere(pos = avc_minus50y, r=4, c='deeppink').legend('sph_CentreOfVentricle')
        sph_cut = Sphere(pos = kspl_CL.points()[num_pt], r=4, c='gold').legend('sph_ChamberCut')
        
        text = filename+"\n\n >> Have a look at the spheres inside the heart \n   and confirm they are correctly named according to the chamber \n   in which they are positioned."
        txt = Text2D(text, c="k", font= font)
        settings.legendSize = .15
        vp = Plotter(N=1, axes = 10)
        vp.show(mesh2cut.alpha(0.01), kspl_CL, sph_cut, sph_atr, sph_vent, txt, at = 0, interactive=True)
        
        atrVentSph_corr = ask4input('Were the atrium and ventricle centre spheres named correctly? \n\t[0]: no, orange is ventricle, and pink is atrium \n\t[1]: yes, orange is atrium, and pink is ventricle! >>>:', bool)
        if not atrVentSph_corr: 
            sph_atr = Sphere(pos = avc_minus50y, r=4, c='deeppink').legend('sph_CentreOfAtrium')
            sph_vent = Sphere(pos = avc_plus50y, r=4, c='darkorange').legend('sph_CentreOfVentricle')
            
        # Add pt to dict
        dict_pts = addPoints2Dict(spheres = [sph_cut, sph_atr, sph_vent], info = ['', 'ChCut','ChCut'], dict_pts = dict_pts)
        dict_pts['numPt_CLChamberCut'] = num_pt
        
        # Now cut myocardium
        print('\n')
        atr_meshes = []; vent_meshes = []
        atr_meshes, vent_meshes, dict_shapes, s3_cyl  = getChamberMeshes(filename = filename,
                                    end_name = ['ch0_cut'], names2cut = ['Myoc'],
                                    kspl_CL = kspl_CL, num_pt = num_pt, atr_meshes = atr_meshes, vent_meshes = vent_meshes,
                                    dir_txtNnpy = dir_txtNnpy, dict_shapes = dict_shapes, dict_pts = dict_pts, 
                                    resolution = resolution, plotshow = True)
        
        happyWithMyocCut = ask4input('Has the myocardium been divided correctly into chambers with the defined disc? \n\t [0]:no, I would like to redefine the disc \n\t [1]:yes, I am happy with the cut made! >>>: ', bool)
        if not happyWithMyocCut:
            happyWithDisc = False
            
    return cyl_final, num_pt, atr_meshes, vent_meshes, dict_shapes, dict_pts, s3_cyl

#%% func - getChamberMeshes
def getChamberMeshes(filename, end_name, names2cut, kspl_CL, num_pt, atr_meshes, vent_meshes, dir_txtNnpy, 
                     dict_shapes, dict_pts, resolution, s3_cyl = [], plotshow = False, mesh2cut = []):
    """
    Function to cut meshes and get its chambers (atrium/ventricle) using the cylinder/disc information given as 
    input (dict_shapes)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX)
    end_name : list of str
        List of names given to the s3 masks saved in dir_txtNnpy
    names2cut : list of str
        List of mesh names being cut
    kspl_CL : Kspline
        Centreline (vedo KSpline)
    num_pt : int
        Index of the centreline point closer to the plane that cuts the meshes into the two chambers.
    atr_meshes: list of meshes
        list of atrial meshes (vedo Meshes)
    vent_meshes: list of meshes
        list of ventricular meshes (vedo Meshes)
    dir_txtNnpy : path
        Path to the folder where the np arrays are saved.
    dict_shapes :  dictionary
        Initialised dictionary with shapes information
    resolution : list of floats
        List with the x,y, z scaling values of the images taken. This information is taken from the metadata of the original file.
    s3_cyl : numpy array of booleans, optional
        Array with disc mask to cut chambers. The size of this array corresponds to the size of the original stack. The default is [].
    plotshow : boolean, optional
        True if you want to see the resulting mesh in a plot, else False. The default is False.
        

    Returns
    -------
    atr_meshes : list of meshes
        List of atrial meshes (vedo Meshes)
    vent_meshes : list of meshes
        List of ventricular meshes (vedo Meshes)
    dict_shapes : dictionary
        Resulting dictionary with shapes information updated
    s3_cyl : numpy array of booleans
        Array with disc mask to cut chambers. The size of this array corresponds to the size of the original stack. 
        
    """

    dORv = filename[9:10]
    if dORv == 'D':
        azimuth = -90
    else:
        azimuth = 0
        
    poss_names = ['ch0_cut','ch1_cut', 'cj', 'ch0_cut_ext', 'ch1_cut_int', 'ch0_cut_int', 'ch1_cut_ext']
    noCut_names = ['ch0_all','ch1_cut', 'cj', 'ch0_ext', 'ch1_int', 'ch0_int', 'ch1_ext']
    
    #Load stack shape 
    [[xdim, ydim, zdim]] = loadNPY(filename, ['stackShape'], dir_txtNnpy, print_txt = False)
    stack_shape = [xdim, ydim, zdim+2]
    
    # If it is the myocardium the one that is being cut, as this is the first mesh that gets cut, it should go into
    # this if. Hence, here define the atrial and ventricular centre spheres and save them or give them as output
    if isinstance(s3_cyl, list) and len(s3_cyl) == 0:
        
        print('- Creating disc mask to cut chambers... (this takes about 2mins)')
        r_circle_max = dict_shapes['cyl2CutChambers_o']['radius_max']
        r_circle_min = dict_shapes['cyl2CutChambers_o']['radius_min']
        normal_unit = dict_shapes['cyl2CutChambers_o']['cyl_axis']
        cl_point = dict_shapes['cyl2CutChambers_o']['cyl_centre']
        
        res_cyl = 2000
        num_rad = int(3*int(r_circle_max)) #int(((r_circle_max-r_circle_min)/0.225)+1)
        num_h = 9; h_min = 0.225/2; h_max = 0.225*2
        for j, rad in enumerate(np.linspace(r_circle_min/2, r_circle_max, num_rad)):
            for i,h in enumerate(np.linspace(h_min,h_max, num_h)):
                cyl = Cylinder(pos = cl_point,r = rad, height = h, axis = normal_unit, c = 'lime', cap = True, res = res_cyl)#.wireframe(True)
                if i == 0 and j == 0:
                    cyl_pts = cyl.points()
                else: 
                    cyl_pts = np.concatenate((cyl_pts, cyl.points()))
        
        cyl.legend('cyl2CutChambers_final')
        cyl_data = [r_circle_max, r_circle_min/2, num_rad, h_max, h_min, num_h, normal_unit, cl_point, res_cyl]
        dict_shapes = addShapes2Dict (shapes = [cyl], dict_shapes = dict_shapes, radius = [cyl_data], print_txt = False)
        
        cyl_points_rot = np.zeros_like(cyl_pts)
        if 'CJ' not in filename: 
            axis = [0,0,1]
        else: 
            axis = [1,0,0]
            
        for i, pt in enumerate(cyl_pts):
            cyl_points_rot[i] = (np.dot(rotation_matrix(axis = axis, theta = np.radians(90)),pt))
            
        cyl_pix = np.transpose(np.asarray([cyl_points_rot[:,i]//resolution[i] for i in range(len(resolution))]))
        cyl_pix = cyl_pix.astype(int)
        cyl_pix = np.unique(cyl_pix, axis =0)
        # print(cyl_pix.shape)
        
        cyl_pix_out = cyl_pix.copy()
        index_out = []
        # Clean cyl_pix if out of stack shape
        for index, pt in enumerate(cyl_pix):
            # print(index, pt)
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
                
        cyl_pix_out = np.delete(cyl_pix_out, index_out, axis = 0)
    
        # Create mask of ring
        s3_cyl = np.zeros((xdim, ydim, zdim+2))
        s3_cyl[cyl_pix_out[:,0],cyl_pix_out[:,1],cyl_pix_out[:,2]] = 1
    
    cl_point = dict_shapes['cyl2CutChambers_o']['cyl_centre']
    
    # Create empty lists to save atrium and ventricles
    if len(end_name) == 1: 
        atr_color = ['lightseagreen']
        vent_color = ['darkturquoise']
    else: 
        atr_color = ['purple', 'orange', 'darkblue', 'maroon', 'tomato', 'cyan']
        vent_color = ['mediumvioletred', 'chocolate', 'indigo', 'crimson', 'skyblue', 'springgreen']

    try: 
        atr_point = dict_pts['sph_CentreOfAtrium-ChCut']['sph_position']
        vent_point = dict_pts['sph_CentreOfVentricle-ChCut']['sph_position']
    except: 
        #Define atrium and ventricle
        try: 
            avc_minus50y = kspl_CL.points()[num_pt-50]
        except: 
            avc_minus50y = kspl_CL.points()[0]
            print('-> First centreline point got selected as AVC-50')
        try: 
            avc_plus50y = kspl_CL.points()[num_pt+50]
        except: 
            avc_plus50y = kspl_CL.points()[-1]
            print('-> Last centreline point got selected as  AVC+50')
            
        sph_atr = Sphere(pos = avc_plus50y, r=4, c='darkorange').legend('sph_CentreOfAtrium')
        sph_vent = Sphere(pos = avc_minus50y, r=4, c='deeppink').legend('sph_CentreOfVentricle')
        sph_cut_o = Sphere(pos = kspl_CL.points()[num_pt], r=4, c='gold').legend('sph_ChamberCut')
        
        text = filename+"\n\n >> Have a look at the spheres inside the heart \n   and confirm they are correctly named according to the chamber \n   in which they are positioned."
        txt = Text2D(text, c="k", font= font)
        settings.legendSize = .15
        vp = Plotter(N=1, axes = 10)
        vp.show(mesh2cut.alpha(0.01), kspl_CL, sph_cut_o, sph_atr, sph_vent, txt, at = 0, interactive=True)
        
        atrVentSph_corr = ask4input('Were the atrium and ventricle centre spheres named correctly? \n\t[0]: no, orange is ventricle, and pink is atrium \n\t[1]: yes, orange is atrium, and pink is ventricle! >>>:', bool)
        if not atrVentSph_corr: 
            sph_atr = Sphere(pos = avc_minus50y, r=4, c='deeppink').legend('sph_CentreOfAtrium')
            sph_vent = Sphere(pos = avc_plus50y, r=4, c='darkorange').legend('sph_CentreOfVentricle')
            
        # Add pt to dict
        dict_pts = addPoints2Dict(spheres = [sph_cut_o, sph_atr, sph_vent], info = ['', 'ChCut','ChCut'], dict_pts = dict_pts)
        atr_point = dict_pts['sph_CentreOfAtrium-ChCut']['sph_position']
        vent_point = dict_pts['sph_CentreOfVentricle-ChCut']['sph_position']
        
    sph_cut = Sphere(pos = cl_point, r=4, c='gold').legend('sph_ChamberCut')
    
    tic = perf_counter()

    for n, s3_name, name in zip(count(), end_name, names2cut):
        # Mask s3s vent and atrium
        settings.legendSize = .3
        if len(end_name) == 1: 
            vp = Plotter (N=1, axes=10)
            text2 = filename+"\n\n >>  Result of dividing heart layers into chambers. "+name
            txt2 = Text2D(text2, c="k", font= font)
        print('\n- Cutting s3 (', name,')')
        try: 
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
        except: 
            end_name = noCut_names[poss_names.index(s3_name)]
            [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [end_name], print_txt = False)
        s3_ring = maskRing2CutChamberS3s(s3_mask = s3, s3_cyl = s3_cyl, stack_shape = stack_shape)
        # plt_s3(start_slc = 0, end_slc = zdim-1, im_every = 20,
        #                 s3_int = s3, s3_ext = s3_cyl, plotshow = plotshow, option = "chamber ring")
        # Create chamber meshes
        m_cut = createLayerMesh(filename = filename, s3 = s3_ring, resolution = resolution, layer = name+'_ChCut', name = name+'_ChCut',
                                        colour = atr_color[n], alpha = 1, plotshow = False)
        
        m_chambers = m_cut.splitByConnectivity(maxdepth=1000)
        alert('clown', 1)
        
        meshes_atr = []
        meshes_vent = []
        bar = Bar('- Asigning meshes to chambers', max=len(m_chambers), suffix = suffix, check_tty=False, hide_cursor=False)
        for mesh in m_chambers:
            centOfMass = mesh.centerOfMass()
            dist_AtrOrVent = [findDist(atr_point, centOfMass),findDist(vent_point, centOfMass)]
            index_AtrOrVent = np.where(dist_AtrOrVent == min(dist_AtrOrVent))[0][0]
            if index_AtrOrVent == 0:
                meshes_atr.append(mesh)
                # print('Atrium')
            else: #1
                meshes_vent.append(mesh)
                # print('Ventr')
            bar.next()
        bar.finish()
        
        print('- Number of individual meshes making part of the atrium: ', len(meshes_atr))
        print('- Number of individual meshes making part of the ventricle: ', len(meshes_vent))
        
        atr_whole = merge(meshes_atr)
        atr_whole.legend(name+'_Atr').color(atr_color[n])
        atr_meshes.append(atr_whole)
        
        vent_whole = merge(meshes_vent)
        try: 
            vent_whole.legend(name+'_Vent').color(vent_color[n])
            vent_meshes.append(vent_whole)
        except: 
            print('- Cut resulted in just one mesh. Redifine disc and try again...')
        
        alert('whistle', 1)
        if len(end_name) == 1: 
            vp.show(atr_whole, vent_whole, kspl_CL, sph_cut, txt2,  at=0, azimuth = azimuth, interactive = True)
            
    toc = perf_counter()
    time = toc-tic
    if len(end_name) != 1: 
        print("- All layers have been cut!  > Total time taken to cut meshes = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")
        alert('jump', 1)

    return atr_meshes, vent_meshes, dict_shapes, s3_cyl

#%% func - getLRMeshes - test in case cutting with ribbon doesn't work!
# def getLRMeshes(filename, end_name, names2cut, kspl_CL, num_pt, atr_meshes, vent_meshes, dir_txtNnpy, 
#                      dict_shapes, dict_pts, resolution, s3_cyl = [], plotshow = False, mesh2cut = []):
    
#     dORv = filename[9:10]
#     if dORv == 'D':
#         azimuth = -90
#     else:
#         azimuth = 0
        
#     #Load stack shape 
#     [[xdim, ydim, zdim]] = loadNPY(filename, ['stackShape'], dir_txtNnpy, print_txt = False)
#     stack_shape = [xdim, ydim, zdim+2]
    
#     rib_pts = createCLRibbon(filename = filename, file_num = file_num, 
#                             df_res = df_res, kspl_CL2use = kspl_CL[0], linLine = linLines[0],
#                             mesh = m_myoc, dict_kspl = dict_kspl, dict_shapes = dict_shapes, 
#                             dict_planes = dict_planes, clRib_type = 'HDStack', plotshow = True)
    
#     rib_pts_rot = np.zeros_like(rib_pts)
#     if 'CJ' not in filename: 
#         axis = [0,0,1]
#     else: 
#         axis = [1,0,0]
    
#     for i, pt in enumerate(rib_pts):
#         rib_pts_rot[i] = (np.dot(rotation_matrix(axis = axis, theta = np.radians(90)),pt))
        
#     rib_pix = np.transpose(np.asarray([rib_pts_rot[:,i]//resolution[i] for i in range(len(resolution))]))
#     rib_pix = rib_pix.astype(int)
#     rib_pix = np.unique(rib_pix, axis =0)
#     # print(cyl_pix.shape)
    
#     rib_pix_out = rib_pix.copy()
#     index_out = []
#     # Clean cyl_pix if out of stack shape
#     for index, pt in enumerate(rib_pix):
#         # print(index, pt)
#         if pt[0] > xdim-2 or pt[0] < 0:
#             delete = True
#         elif pt[1] > ydim-2 or pt[1] < 0:
#             delete = True
#         elif pt[2] > zdim+2-1 or pt[2] < 0:
#             delete = True
#         else: 
#             delete = False
        
#         if delete:
#             # print(pt)
#             index_out.append(index)
            
#     rib_pix_out = np.delete(rib_pix_out, index_out, axis = 0)
#     # Create mask of ring
#     s3_rib = np.zeros((xdim, ydim, zdim+2))
#     s3_rib[rib_pix_out[:,0],rib_pix_out[:,1],rib_pix_out[:,2]] = 1
            
    
#     for n, s3_name, name in zip(count(), end_name, names2cut):
#         # Mask s3s vent and atrium
#         settings.legendSize = .3
#         if len(end_name) == 1: 
#             vp = Plotter (N=1, axes=10)
#             text2 = filename+"\n\n >>  Result of dividing heart layers into chambers. "+name
#             txt2 = Text2D(text2, c="k", font= font)
#         print('\n- Cutting s3 (', name,')')
#         try: 
#             [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [s3_name], print_txt = False)
#         except: 
#             end_name = noCut_names[poss_names.index(s3_name)]
#             [s3], _ = loadStacks(filename = filename, dir_txtNnpy = dir_txtNnpy, end_name = [end_name], print_txt = False)
#         s3_ring = maskRing2CutChamberS3s(s3_mask = s3, s3_cyl = s3_cyl, stack_shape = stack_shape)
#         # plt_s3(start_slc = 0, end_slc = zdim-1, im_every = 20,
#         #                 s3_int = s3, s3_ext = s3_cyl, plotshow = plotshow, option = "chamber ring")
#         # Create chamber meshes
#         m_cut = createLayerMesh(filename = filename, s3 = s3_ring, resolution = resolution, layer = name+'_ChCut', name = name+'_ChCut',
#                                         colour = atr_color[n], alpha = 1, plotshow = False)
        
#%% >>> PLANES
#%% func - createPlane
def createPlane(dict_planes, name):
    """
    Function that creates the 'name' plane from the dict_plane

    Parameters
    ----------
    dict_planes : dictionary
        Initialised dictionary with planes information
    name : str
        Plane name.

    Returns
    -------
    plane_out : Plane
        vedo Plane.

    """

    normal = dict_planes[name]['pl_normal']
    centre = dict_planes[name]['pl_centre']
    color = dict_planes[name]['color']

    plane_out = Plane(pos = centre, normal = normal, sx = 300).color(color).alpha(1)
    plane_out.legend(name)

    return plane_out

#%% func - getPlane
def getPlane(filename, type_cut, info, mesh_in, mesh_out = '', option = [True,True,True,True,True,True], dict_planes = []):
    """
    Function that creates a plane defined by the user

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    type_cut : str
        Text defining the type of cut that is going to be made with the defined plane
    info : str
        Text with additional information to name plane.
    mesh_in : mesh
        Internal mesh (vedo Mesh)
    mesh_out : mesh
        External mesh (vedo Mesh)
    option : list of booleans
        List of booleans indicating the sliders to use in getPlanePos function.
        [sliderX, sliderY, sliderZ, sliderRotX, sliderRotY, sliderRotZ]
        The default is [True,True,True,True,True,True].

    Returns
    -------
    plane_new : Plane
        Final plane defined by the user
    pl_centre : list of floats
        List with the x,y,z coordinates of the plane's centre
    normal_corrected : list of floats
        List with the x,y,z coordinatesof the plane's normal rotated in 3D

    """

    print('- Getting plane to cut '+ type_cut)
    #mesh1.alpha(0.05); mesh2.alpha(0.05)
    while True:
        # Create plane
        if info == '4CL' or info == 'redefCL':
            try: 
                centre_o = dict_planes['pl2CutMesh_'+type_cut]['pl_centre']
                normal_o = dict_planes['pl2CutMesh_'+type_cut]['pl_normal']
                # print('predef plane found')
                plane, normal, rotX, rotY, rotZ = getPlanePos(filename, type_cut, mesh_in.bounds(), option,  mesh_in, mesh_out, centre_o, normal_o)
            except: 
                # print('NO predef plane found')
                plane, normal, rotX, rotY, rotZ = getPlanePos(filename, type_cut, mesh_in.bounds(), option,  mesh_in, mesh_out)
        else: 
            plane, normal, rotX, rotY, rotZ = getPlanePos(filename, type_cut, mesh_in.bounds(), option,  mesh_in, mesh_out)
        # Get new normal of rotated plane
        normal_corrected = newNormal3DRot(normal, rotX, rotY, rotZ)
        #normal_corrected = np.asarray(newNormal(normal, rotX))
        # Get central point of new plane and create sphere
        pl_centre = plane.pos()
        sph_centre = Sphere(pos=pl_centre,r=2,c="purple")
        # Build new plane to confirm
        plane_new = Plane(pos=pl_centre,normal=normal_corrected, sx=500).color("green").alpha(1).legend('pl2CutMesh'+info+'_'+type_cut)

        text = filename+"\n\n >> Confirm plane position to proceed with the cut ("+type_cut+")\n >> Close the window when done"
        txt = Text2D(text, c=c, font=font)

        settings.legendSize = .2
        vp = Plotter(N=1, axes=4)
        if mesh_out != '':
            vp.show(mesh_in, mesh_out, plane, plane_new, sph_centre, txt, at=0, viewup="y", azimuth=0, elevation=0, interactive=True)
        else:
            vp.show(mesh_in, plane, plane_new, sph_centre, txt, at=0, viewup="y", azimuth=0, elevation=0, interactive=True)
        if info == 'redefCL':
            add_text = 'smoothed mesh with the defined plane to define new centreline end point?'
        else: 
            add_text = type_cut +' with the defined plane?'
        happy = ask4input('Do you want to cut the '+add_text+' \n  [0]:no, I would like to define a new plane/[1]:yes, continue!: ', bool)
        if happy:
            break

    return plane_new, pl_centre, normal_corrected

#%% func - getPlanePos
def getPlanePos (filename, type_cut, xyz_bounds, option, mesh_in, mesh_out = '', centre = [], normal = (0,1,0)):
    """
    Function that shows a plot so that the user can define a plane (mesh opacity can be changed)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    type_cut : str
        Text defining the type of cut that is going to be made with the defined plane
    xyz_bounds : list of floats
        x,y,z boundaries of mesh_out
    option : list of booleans
        List of booleans indicating the sliders to use in getPlanePos function.
        [sliderX, sliderY, sliderZ, sliderRotX, sliderRotY, sliderRotZ]
    mesh_in : mesh
        Internal mesh (vedo Mesh)
    mesh_out : mesh, optional
        If given as input, External mesh (vedo Mesh). The default is ''.
    centre : list of floats, optional
        List with the x,y,z coordinates of the initial plane's centre. The default is [].
    normal : list of floats, optional
        List with the x,y,z coordinates of the initial plane's normal. The default is (0,1,0).
        
    Returns
    -------
    plane : Plane
        Final plane defined by the user
    normal : list of floats
        List with the x,y,z coordinates of the plane's normal
    rotX : list of floats
        List of angles (deg) of the resulting rotation around the x-axis.
    rotY : list of floats
        List of angles (deg) of the resulting rotation around the y-axis.
    rotZ : list of floats
        List of angles (deg) of the resulting rotation around the z-axis.

    """
    

    xmin, xmax, ymin, ymax, zmin, zmax = xyz_bounds
    x_size = xmax - xmin
    y_size = ymax - ymin
    z_size = zmax - zmin
    
    xval = sorted([xmin-0.3*x_size,xmax+0.3*x_size])
    yval = sorted([ymin-0.3*y_size,ymax+0.3*y_size])
    zval = sorted([zmin-0.3*z_size,zmax+0.3*z_size])
    
    box_size = max(x_size, y_size, z_size)*1.2
    
    if centre == []: 
        centre = (x_size/2+xmin, ymin, z_size/2+zmin)
    # normal = (0,1,0)
    #print("centre:", centre)

    rotX = [0]
    rotY = [0]
    rotZ = [0]

    # Functions to move and rotate plane
    def sliderX(widget, event):
        valueX = widget.GetRepresentation().GetValue()
        plane.x(valueX)

    def sliderY(widget, event):
        valueY = widget.GetRepresentation().GetValue()
        plane.y(valueY)

    def sliderZ(widget, event):
        valueZ = widget.GetRepresentation().GetValue()
        plane.z(valueZ)

    def sliderRotX(widget, event):
        valueRX = widget.GetRepresentation().GetValue()
        rotX.append(valueRX)
        plane.rotateX(valueRX, rad=False)

    def sliderRotY(widget, event):
        valueRY = widget.GetRepresentation().GetValue()
        rotY.append(valueRY)
        plane.rotateY(valueRY, rad=False)

    def sliderRotZ(widget, event):
        valueRZ = widget.GetRepresentation().GetValue()
        rotZ.append(valueRZ)
        plane.rotateZ(valueRZ, rad=False)

    def sliderAlphaMeshOut(widget, event):
        valueAlpha = widget.GetRepresentation().GetValue()
        mesh_out.alpha(valueAlpha)

    settings.legendSize = .2
    vp = Plotter(N=1, axes=8)
    plane = Plane(pos=centre, normal=normal, sx=box_size*1.2).color("gainsboro").alpha(1)

    if option[0]: #sliderX
        vp.addSlider2D(sliderX, xval[0], xval[1], value=centre[0],
                    pos=[(0.1,0.15), (0.3,0.15)], title="- > x position > +", c="crimson" )
    if option[1]: #sliderY
        vp.addSlider2D(sliderY, yval[0], yval[1], value=centre[1],
                    pos=[(0.4,0.15), (0.6,0.15)], title="- > y position > +", c="dodgerblue" )
    if option[2]: #sliderZ
        vp.addSlider2D(sliderZ, zval[0], zval[1], value=centre[2],
                    pos=[(0.7,0.15), (0.9,0.15)], title="- > z position > +", c="limegreen")
    if option[3]: #sliderRotX
        vp.addSlider2D(sliderRotX, -1, +1, value=0,
                    pos=[(0.1,0.05), (0.3,0.05)], title="- > x rotation > +", c="deeppink")
    if option[4]: #sliderRotY
        vp.addSlider2D(sliderRotY, -1, +1, value=0,
                    pos=[(0.4,0.05), (0.6,0.05)], title="- > y rotation > +", c="gold")
    if option[5]: #sliderRotZ
        vp.addSlider2D(sliderRotZ, -1, +1, value=0,
                    pos=[(0.7,0.05), (0.9,0.05)], title="- > z rotation > +", c="teal")

    if mesh_out != '':
        vp.addSlider2D(sliderAlphaMeshOut, xmin=0.01, xmax=0.99, value=0.01,
                   pos=[(0.95,0.25), (0.95,0.45)], c="blue", title="Ext.Mesh Opacity)")

    text = filename+"\n\n >> Define plane position to make cut ("+type_cut+")\n >> Close the window when done"
    txt = Text2D(text, c="k", font= font)
    if mesh_out != '':
        vp.show(mesh_in, mesh_out, plane, txt, viewup="y", zoom=1, interactive=True)
    else:
        vp.show(mesh_in, plane, txt, viewup="y", zoom=1, interactive=True)
    #azimuth=-90, elevation=0,

    return plane, normal, rotX, rotY, rotZ

#%% func - modifyPlane
def modifyPlane(filename, pl_normal, pl_centre, type_cut, mesh1, xyz_bounds, option):
    """
    Function that shows a plot so that the user can define a plane

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    pl_normal : list of floats
        List with the x,y,z coordinatesof the plane's normal
    pl_centre : list of floats
        List with the x,y,z coordinates of the plane's centre
    type_cut : str
        Text defining the type of cut that is going to be made with the modified plane
    mesh1 : mesh
        Myocardial mesh to use as background to define plane
    xyz_bounds : list of floats
        x,y,z boundaries of mesh_out
    option : list of booleans
        List of booleans indicating the sliders to use in getPlanePos function.
        [sliderX, sliderY, sliderZ, sliderRotX, sliderRotY, sliderRotZ]

    Returns
    -------
    plane : Plane
        Final plane defined by the user
    rotX : list of floats
        List of angles (deg) of the resulting rotation around the x-axis.
    rotY : list of floats
        List of angles (deg) of the resulting rotation around the y-axis.
    rotZ : list of floats
        List of angles (deg) of the resulting rotation around the z-axis.

    """

    xmin, xmax, ymin, ymax, zmin, zmax = xyz_bounds
    x_size = xmax - xmin
    y_size = ymax - ymin
    z_size = zmax - zmin

    box_size = max(x_size, y_size, z_size)*1.2
    # centre = (x_size/2+xmin, ymin, z_size/2+zmin)
    # normal = (0,1,0)
    #print("centre:", centre)

    rotX = [0]
    rotY = [0]
    rotZ = [0]

    # Functions to move and rotate plane
    def sliderX(widget, event):
        valueX = widget.GetRepresentation().GetValue()
        plane.x(valueX)

    def sliderY(widget, event):
        valueY = widget.GetRepresentation().GetValue()
        plane.y(valueY)

    def sliderZ(widget, event):
        valueZ = widget.GetRepresentation().GetValue()
        plane.z(valueZ)

    def sliderRotX(widget, event):
        valueRX = widget.GetRepresentation().GetValue()
        rotX.append(valueRX)
        plane.rotateX(valueRX, rad=False)

    def sliderRotY(widget, event):
        valueRY = widget.GetRepresentation().GetValue()
        rotY.append(valueRY)
        plane.rotateY(valueRY, rad=False)

    def sliderRotZ(widget, event):
        valueRZ = widget.GetRepresentation().GetValue()
        rotZ.append(valueRZ)
        plane.rotateZ(valueRZ, rad=False)

    settings.legendSize = .3
    vp = Plotter(N=1, axes=8)
    plane = Plane(pos=pl_centre, normal=pl_normal, sx=box_size*1.2).color("gainsboro").alpha(1)

    if option[0]: #sliderX
        vp.addSlider2D(sliderX, xmin*0.8, xmax*1.2, value=pl_centre[0],
                    pos=[(0.1,0.15), (0.3,0.15)], title="- > x position > +", c="crimson" )
    if option[1]: #sliderY
        vp.addSlider2D(sliderY, ymin*1.2, ymax*0.8, value=pl_centre[1],
                    pos=[(0.4,0.15), (0.6,0.15)], title="- > y position > +", c="dodgerblue" )
    if option[2]: #sliderZ
        vp.addSlider2D(sliderZ, zmin*1.2, zmax*0.8, value=pl_centre[2],
                    pos=[(0.7,0.15), (0.9,0.15)], title="- > z position > +", c="limegreen")
    if option[3]: #sliderRotX
        vp.addSlider2D(sliderRotX, -1, +1, value=0,
                    pos=[(0.1,0.05), (0.3,0.05)], title="- > x rotation > +", c="deeppink")
    if option[4]: #sliderRotY
        vp.addSlider2D(sliderRotY, -1, +1, value=0,
                    pos=[(0.4,0.05), (0.6,0.05)], title="- > y rotation > +", c="gold")
    if option[5]: #sliderRotZ
        vp.addSlider2D(sliderRotZ, -1, +1, value=0,
                    pos=[(0.7,0.05), (0.9,0.05)], title="- > z rotation > +", c="teal")

    text = filename+"\n\n >> Define plane position to make cut ("+type_cut+")\n >> Close the window when done"
    txt = Text2D(text, c="k", font= font)
    vp.show(mesh1, plane, txt, viewup="y", zoom=1, interactive=True)
    #azimuth=-90, elevation=0,

    return plane, rotX, rotY, rotZ

#%% func - modifyDisc
def modifyDisc (filename, pl_normal, pl_centre, radius, type_cut, mesh1, xyz_bounds, option = [True,True,True,True,True,True]):
    """
    Function that shows a plot so that the user can define a cylinder (disc)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    pl_normal : list of floats
        List with the x,y,z coordinatesof the plane's normal
    pl_centre : list of floats
        List with the x,y,z coordinates of the plane's centre
    radius: float
        Initial radius value to create disc to cut mesh
    type_cut : str
        Text defining the type of cut that is going to be made with the modified plane
    mesh1 : mesh
        Myocardial mesh to use as background to define plane
    xyz_bounds : list of floats
        x,y,z boundaries of mesh_out
    option : list of booleans
        List of booleans indicating the sliders to use in getPlanePos function.
        [sliderX, sliderY, sliderZ, sliderRotX, sliderRotY, sliderRotZ]

    Returns
    -------
    cyl_test : Cylinder
        Final cylinder position defined by the user
    sph_cyl : Sphere
        Final disc centre position defined by the user
    rotX : list of floats
        List of angles (deg) of the resulting rotation around the x-axis.
    rotY : list of floats
        List of angles (deg) of the resulting rotation around the y-axis.
    rotZ : list of floats
        List of angles (deg) of the resulting rotation around the z-axis.

    """

    xmin, xmax, ymin, ymax, zmin, zmax = xyz_bounds
    # x_size = xmax - xmin
    # y_size = ymax - ymin
    # z_size = zmax - zmin

    rotX = [0]
    rotY = [0]
    rotZ = [0]

    # Functions to move and rotate cyl_test
    def sliderX(widget, event):
        valueX = widget.GetRepresentation().GetValue()
        cyl_test.x(valueX)
        sph_cyl.x(valueX)

    def sliderY(widget, event):
        valueY = widget.GetRepresentation().GetValue()
        cyl_test.y(valueY)
        sph_cyl.y(valueY)

    def sliderZ(widget, event):
        valueZ = widget.GetRepresentation().GetValue()
        cyl_test.z(valueZ)
        sph_cyl.z(valueZ)

    def sliderRotX(widget, event):
        valueRX = widget.GetRepresentation().GetValue()
        rotX.append(valueRX)
        cyl_test.rotateX(valueRX, rad=False)
        sph_cyl.rotateX(valueRX, rad=False)

    def sliderRotY(widget, event):
        valueRY = widget.GetRepresentation().GetValue()
        rotY.append(valueRY)
        cyl_test.rotateY(valueRY, rad=False)
        sph_cyl.rotateY(valueRY, rad=False)

    def sliderRotZ(widget, event):
        valueRZ = widget.GetRepresentation().GetValue()
        rotZ.append(valueRZ)
        cyl_test.rotateZ(valueRZ, rad=False)
        sph_cyl.rotateZ(valueRZ, rad=False)
        
    def sliderAlphaMeshOut(widget, event):
        valueAlpha = widget.GetRepresentation().GetValue()
        mesh1.alpha(valueAlpha)

    settings.legendSize = .3
    vp = Plotter(N=1, axes=8)
    cyl_test = Cylinder(pos = pl_centre,r = radius, height = 2*0.225, axis = pl_normal, c = 'purple', cap = True, res = 300)
    sph_cyl = Sphere(pos = pl_centre, r=4, c='gold')
    
    if option[0]: #sliderX
        vp.addSlider2D(sliderX, xmin*0.8, xmax*1.2, value=pl_centre[0],
                    pos=[(0.1,0.15), (0.3,0.15)], title="- > x position > +", c="crimson" )
    if option[1]: #sliderY
        vp.addSlider2D(sliderY, ymin*1.2, ymax*0.8, value=pl_centre[1],
                    pos=[(0.4,0.15), (0.6,0.15)], title="- > y position > +", c="dodgerblue" )
    if option[2]: #sliderZ
        vp.addSlider2D(sliderZ, zmin*1.2, zmax*0.8, value=pl_centre[2],
                    pos=[(0.7,0.15), (0.9,0.15)], title="- > z position > +", c="limegreen")
    if option[3]: #sliderRotX
        vp.addSlider2D(sliderRotX, -1, +1, value=0,
                    pos=[(0.1,0.05), (0.3,0.05)], title="- > x rotation > +", c="deeppink")
    if option[4]: #sliderRotY
        vp.addSlider2D(sliderRotY, -1, +1, value=0,
                    pos=[(0.4,0.05), (0.6,0.05)], title="- > y rotation > +", c="gold")
    if option[5]: #sliderRotZ
        vp.addSlider2D(sliderRotZ, -1, +1, value=0,
                    pos=[(0.7,0.05), (0.9,0.05)], title="- > z rotation > +", c="teal")
        
    vp.addSlider2D(sliderAlphaMeshOut, xmin=0.01, xmax=0.99, value=0.01,
               pos=[(0.95,0.25), (0.95,0.45)], c="blue", title="Myocardial Opacity)")

    text = filename+"\n\n>> Define disc position to make cut ("+type_cut+")\n   Make sure it cuts through the AVC and effectively separates the chambers. \n>> Close the window when done. \n>> Note: Initially defined disc radius is 60um. Once you have selected the\n   position of the disc a new radius will be calculated based \n   on the mesh points the disc cuts. \n   If you are not happy with the disc radius, \n   you will be able to modify it just before proceeding to the cut."
    txt = Text2D(text, c="k", font= font)
    vp.show(mesh1, cyl_test, sph_cyl, txt, viewup="y", zoom=1, interactive=True)
    #azimuth=-90, elevation=0,

    return cyl_test, sph_cyl, rotX, rotY, rotZ

#%% func - getPlaneNormal2Pt
def getPlaneNormal2Pt (pt_num, spline_pts):
    """
    Funtion that gets a plane normal to a point in a spline

    Parameters
    ----------
    pt_num : int
        Index of the centreline point closer to the plane that cuts the meshes into the two chambers.
    spline_pts : list of x,y,z points coordinates
        x,y,z points coordinates of spline

    Returns
    -------
    normal : list of floats
        List with the x,y,z coordinatesof the final plane's normal
    pt_centre : list of floats
        List with the x,y,z coordinates of the final plane's centre

    """

    pt_centre = spline_pts[pt_num]
    normal = spline_pts[pt_num-1]-spline_pts[pt_num+1]

    return normal, pt_centre

#%% func - rotatePlane2Images
def rotatePlane2Images (pl_centre, pl_normal, type_cut, cj = False):
    """
    Function that rotates the planes defined in the surface reconstructions to the images mask

    Parameters
    ----------
    pl_centre : list of floats
        List with the x,y,z coordinates of the plane's centre defined in the mesh
    pl_normal : list of floats
        List with the x,y,z coordinates of the plane's normal defined in the mesh
    type_cut : str
        Text defining the type of cut that is going to be made with the defined plane
    cj : bool, optional
        True if analysing Anjalie's data, else False. Default is False. 

    Returns
    -------
    plane_im : Plane
        Final plane to cut masks
    pl_im_centre : list of floats
        List with the x,y,z coordinates of the plane's centre to cut the masks
    pl_im_normal : list of floats
        List with the x,y,z coordinates of the plane's normal to cut the masks

    """
    if not cj: 
        pl_im_centre = (np.dot(rotation_matrix(axis = [0,0,1], theta = np.radians(90)), pl_centre))
        pl_im_normal = (np.dot(rotation_matrix(axis = [0,0,1], theta = np.radians(90)), pl_normal))
    else: 
        pl_im_centre = (np.dot(rotation_matrix(axis = [1,0,0], theta = np.radians(90)), pl_centre))
        pl_im_normal = (np.dot(rotation_matrix(axis = [1,0,0], theta = np.radians(90)), pl_normal))

    plane_im = Plane(pos=pl_im_centre,normal=pl_im_normal, sx=500).color("blue").alpha(0.5).legend('pl2CutIm'+'_'+type_cut)

    return plane_im, pl_im_centre, pl_im_normal

#%% func - createDVPlanes
def createDVPlanes(filename, sph_orient, mesh, kspl_CL, orient_lines, dict_planes):
    """
    Function that creates the dorso-ventral planes dividing each of the heart chambers

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    sph_orient : list of spheres
        list of spheres positioned in particular regions of the heart [sph_atr, sph_vent, sph_outf, sph_inf, sph_valve]
    mesh : mesh
        Mesh (vedo Mesh)
    kspl_CL : Kspline
        Centreline (vedo KSpline)
    orient_lines : list of lines
        List of lines defining the orientation of each chamber
    dict_planes : dictionary
        Initialised dictionary with planes information

    Returns
    -------
    pl_DnV : list of planes
        [dorso-ventral plane for the atrium, dorso-ventral plane for the ventricle]
    dict_planes : dictionary
        Resulting dictionary with planes information updated

    """

    sph_atr, sph_vent, sph_outf, sph_inf, sph_valve = sph_orient

    # Atrium
    ptsDnV_Atr = Points([sph_inf.pos(), sph_atr.pos(), sph_valve.pos()])
    pl_DnV_Atr = fitPlane(ptsDnV_Atr.points()).scale(4).c("palegreen").alpha(1).legend("pl_AtrCoronal")
    pl_DnVAtr_centre = pl_DnV_Atr.center
    pl_DnVAtr_normal = pl_DnV_Atr.normal

    # Ventricle
    ptsDnV_Vent = Points([sph_vent.pos(), sph_outf.pos(), sph_valve.pos()])
    pl_DnV_Vent = fitPlane(ptsDnV_Vent.points()).scale(4).c("lightskyblue").alpha(1).legend("pl_VentCoronal")
    pl_DnVVent_centre = pl_DnV_Vent.center
    pl_DnVVent_normal = pl_DnV_Vent.normal

    # Save on dict_planes
    dict_planes = addPlanes2Dict (planes = [pl_DnV_Atr, pl_DnV_Vent], pls_centre = [pl_DnVAtr_centre, pl_DnVVent_centre],
                                  pls_normal = [pl_DnVAtr_normal, pl_DnVVent_normal], info =['',''], dict_planes = dict_planes)

    text = filename+"\n\n >> Creating coronal planes to divide each chamber"
    txt = Text2D(text, c="k", font= font)
    settings.legendSize = .3
    vp = Plotter(N=1, axes = 13)
    vp.show(mesh.alpha(0.01), sph_orient, kspl_CL, orient_lines, pl_DnV_Atr, pl_DnV_Vent, txt, at=0, interactive=True)

    pl_DnV = [pl_DnV_Atr, pl_DnV_Vent]

    return pl_DnV, dict_planes

#%% >>> SHAPES
#%% func - sphInSpline
def sphInSpline(kspl_CL, colour = False, name = '', every = 10):
    """
    Function that creates a group of spheres through a spline given as input.

    Parameters
    ----------
    kspl_CL : Kspline
        Centreline (vedo KSpline)
    colour : bool, optional
        If True color spheres according to position within centreline, else False. The default if False
    name : str, optional
        Name given to the group of spheres. The default is ''.
    every : int (1) or float, optional
        Value that defines how close together the created spheres will be. The default is 10.

    Returns
    -------
    spheres_spline : list of spheres/Spheres
        list of spheres (vedo Sphere) / Spheres (vedo Sphere)

    """
    if colour: 
        # ypos = [kspl_CL.points()[i][1] for i in range(len(kspl_CL.points()))]
        # ymin, ymax = np.min(ypos), np.max(ypos)
        # print("min and max of variances:", ymin, ymax)
        vcols = [colorMap(v, "jet", 0, len(kspl_CL.points())) for v in list(range(len(kspl_CL.points())))]  # scalars->colors
    
    if every > 1:
        spheres_spline = []
        for num, point in enumerate(kspl_CL.points()):
            if num % every == 0 or num == kspl_CL.NPoints()-1:
                if colour:
                    sphere_pt = Sphere(pos=point, r=2, c=vcols[num]).addScalarBar(title='Centreline\nPoint Number')
                else:
                    sphere_pt = Sphere(pos=point, r=2, c='coral')
                spheres_spline.append(sphere_pt)
    else:
        kspl_new = KSpline(kspl_CL.points(), res = round(kspl_CL.NPoints()/every))
        spheres_spline = Spheres(kspl_new.points(), c='coral', r=2)
        if name != '':
            spheres_spline.legend(name)

    return spheres_spline

#%% func - createKSpls
def createKSpls(dict_kspl, kspl_list):
    """
    Function that creates a KSpline using the points given as input in the dict_kspl

    Parameters
    ----------
    dict_kspl :  dictionary
        Initialised dictionary with kspline information
    kspl_list : list of str
        List of names of the ksplines to create

    Returns
    -------
    kSplines_out : List of Ksplines
        Resulting list of ksplines created

    """

    print('- Creating ksplines...')
    kSplines_out = []
    for num, kspl_name in enumerate(list(dict_kspl.keys())):
        if kspl_name in kspl_list:
            pts_kspl = dict_kspl[kspl_name]['kspl_pts']
            color_set = dict_kspl[kspl_name]['color'] # color_set = colors[num] #
            if 'Ext' in kspl_name or 'CL_' in kspl_name or 'ext' in kspl_name:
                kspl = KSpline(pts_kspl, res = 200)
            elif 'lin' in kspl_name:
                kspl = KSpline(pts_kspl, continuity=0, tension=0, bias=0, closed=False)
            else:
                kspl = KSpline(pts_kspl, continuity=0, tension=0, bias=0, closed=True)

            kspl.color(color_set).lw(2).legend(kspl_name)

            kSplines_out.append(kspl)

    print('\t>> KSplines created')

    return kSplines_out

#%% func - createSpheres
def createSpheres(dict_pts, pts_list):
    """
    Function that creates spheres using the points given as input in the dict_pts

    Parameters
    ----------
    dict_pts : dictionary
        Initialised dictionary with points information
    pts_list : list of str
        List of names of the spheres to create

    Returns
    -------
    sph_out : list of spheres
        list of resulting spheres (vedo Sphere)

    """

    print('- Creating spheres...')
    sph_out = []
    for num, pt_name in enumerate(list(dict_pts.keys())):
        if pt_name in pts_list:
            pt_pos = dict_pts[pt_name]['sph_position']
            #print('original pos:', pt_pos)
            color_set = dict_pts[pt_name]['color'] # color_set = colors[num] #
            sph = Sphere(pos=pt_pos, r=2, c=color_set)
            sph.legend(pt_name)

            sph_out.append(sph)

    print('\t>> Spheres created')

    return sph_out

#%% func - createCLs
def createCLs(filename, dict_cl, dict_pts, dict_kspl, dict_planes, colors, myoc, dir_stl):
    """
    Function that creates the centrelines using the points given as input in the dict_cl

    Parameters
    ----------
    dict_cl : dictionary
        Initialised dictionary with centreline information.
    dict_pts : dictionary
        Initialised dictionary with points information.
    dict_kspl :  dictionary
        Initialised dictionary with kspline information.
    dict_shapes : dictionary
        Initialised dictionary with shapes information.
    colors : list of str
        List of colours.

    Returns
    -------
    kspl_CL : Kspline
        Centreline (vedo KSpline)
    linLines : Linear lines
        Lines (vedo Line)
    spheres_CL : Spheres
        Group of spheres with maximum inscribed centreline radius and positions given by VMTK and color coded by actual sphere radius.
    spheres_CL_col : Spheres
        Group of spheres with actual sphere radius in the positions given by VMTK.
    dict_shapes : dictionary
        Resulting dictionary with shapes information updated.
    dict_kspl : dictionary
        Resulting dictionary with ksplines information updated.

    """

    print('- Creating centreline(s)...')
    kspl_CL = []
    linLines = []

    layerNames = ['Int.Myoc(Cut)', 'Ext.Endo(Cut)']
    sph_genName  = ['sph_Cut4CL_inflow-', 'sph_Cut4CL_outflow-']
    color_linLines = ['lawngreen', 'blueviolet']

    # Iterate through layer
    for i, dictLayer in enumerate(dict_cl):
        #Get cl points from vmtk
        pts_cl = np.asarray(dictLayer['Points'])
        # Interpolate points of original centreline
        pts_int_o = getInterpolatedPts(points=pts_cl, nPoints = 300)
        # Create kspline with original points
        kspl_o = KSpline(pts_int_o, res = 300).color('darkorange').legend('CLo_'+layerNames[i]).lw(5)
        # Create IFT and OFT spheres to show original points position 
        sph_inf_o = Sphere(pts_int_o[-1], r = 3, c='tomato')
        sph_outf_o = Sphere(pts_int_o[0], r = 3, c='navy')
        
        # Get outflow point for that layer
        sph_outf = sph_genName[1]+layerNames[i]
        pt2add_outf = dict_pts[sph_outf]['sph_position']
        pts_withOutf = np.insert(pts_cl, 0, np.transpose(pt2add_outf), axis=0)

        #Get inflow point for that layer (three options)
        # - Option 1 (add point obtained when creating mesh for centreline - vmtk)
        sph_inf = sph_genName[0]+layerNames[i]
        pt2add_inf = dict_pts[sph_inf]['sph_position']
        pts_all_opt1 = np.insert(pts_withOutf, len(pts_withOutf), np.transpose(pt2add_inf), axis=0)

        # Interpolate points
        pts_int_opt1 = getInterpolatedPts(points=pts_all_opt1, nPoints = 300)
        # Create kspline with points
        kspl_opt1 = KSpline(pts_int_opt1, res = 300)
        kspl_opt1.color(colors[i]).legend('CL_'+layerNames[i]).lw(5)
        
        # - Option 2 (add point of extended original centreline)
        pt_m10 = pts_int_o[-10]
        sph_m10 = Sphere(pos = pt_m10, r=4, c='lime')
        pt_m1 = pts_int_o[-1]
        sph_m1 = Sphere(pos = pt_m1, r=4, c='tomato')
        dir_v = unit_vector(pt_m1-pt_m10)*3
        pts_ext = np.array([pt_m10, pt_m1, pt_m1+dir_v*2,pt_m1+dir_v*4,pt_m1+dir_v*6,pt_m1+dir_v*8,pt_m1+dir_v*10])
        # print(pts_ext)
        xd = np.diff(pts_ext[:,0])
        yd = np.diff(pts_ext[:,1])
        zd = np.diff(pts_ext[:,2])
        dist = np.sqrt(xd**2+yd**2+zd**2)
        u = np.cumsum(dist)
        u = np.hstack([[0],u])
        t = np.linspace(0, u[-1],100)
        resamp_pts = interpn((u,), pts_ext, t)
        # print(len(resamp_pts))
        kspl_resamp = KSpline(resamp_pts, res = 100).lw(5).color('deeppink')
        
        pl_centre = dict_planes['pl2CutMesh4CL_inflow']['pl_centre']
        pl_normal = dict_planes['pl2CutMesh4CL_inflow']['pl_normal']
        pl_IFT = Plane(pos = pl_centre, normal = pl_normal)
        kspl_testB = kspl_resamp.clone().cutWithMesh(pl_IFT, invert = True).color('gold')
        pt2add = kspl_testB.points()[-1]

        pts_all_opt2 = np.insert(pts_withOutf, len(pts_withOutf), np.transpose(pt2add), axis=0)

        # Interpolate points
        pts_int_opt2 = getInterpolatedPts(points=pts_all_opt2, nPoints = 300)
        # Create kspline with points
        kspl_opt2 = KSpline(pts_int_opt2, res = 300)
        kspl_opt2.color(colors[i]).legend('CL_'+layerNames[i]).lw(5)
        
        # - Option 3 (add point of extended original centreline midline between chamber centre and in/outf tract)
        pt_m30 = pts_int_o[-50]
        sph_m30 = Sphere(pos = pt_m30, r=4, c='purple')
        pt_m1 = pts_int_o[-1]
        dir_v30 = unit_vector(pt_m1-pt_m30)*3
        pts_ext = np.array([pt_m30, pt_m1, pt_m1+dir_v30*2,pt_m1+dir_v30*4,pt_m1+dir_v30*6,pt_m1+dir_v30*8,pt_m1+dir_v30*10])
        # print(pts_ext)
        xd = np.diff(pts_ext[:,0])
        yd = np.diff(pts_ext[:,1])
        zd = np.diff(pts_ext[:,2])
        dist = np.sqrt(xd**2+yd**2+zd**2)
        u = np.cumsum(dist)
        u = np.hstack([[0],u])
        t = np.linspace(0, u[-1],100)
        resamp_pts = interpn((u,), pts_ext, t)
        # print(len(resamp_pts))
        kspl_resamp = KSpline(resamp_pts, res = 100).lw(5).color('deeppink')
    
        kspl_testC = kspl_resamp.clone().cutWithMesh(pl_IFT, invert = True).color('gold')
        pt2add30 = kspl_testC.points()[-1]

        pts_all_opt3 = np.insert(pts_withOutf, len(pts_withOutf), np.transpose(pt2add30), axis=0)

        # Interpolate points
        pts_int_opt3 = getInterpolatedPts(points=pts_all_opt3, nPoints = 300)
        # Create kspline with points
        kspl_opt3 = KSpline(pts_int_opt3, res = 300)
        kspl_opt3.color(colors[i]).legend('CL_'+layerNames[i]).lw(5)
        
        ksp_sel = [kspl_opt1, kspl_opt2, kspl_opt3]
        pt2add_sel = [pt2add_inf, pt2add, pt2add30]
        
        #Select the CL to use
        # First look if one has already been selected 
        try: 
            q_select = dict_kspl['kspl_opt_sel']
        # If not, plot all three options and ask for selection
        except: 
            settings.legendSize = .20
            txt1 = Text2D('Option 0', c=c, font = font)
            txt2 = Text2D('Option 1', c=c, font = font)
            txt3 = Text2D('Option 2', c=c, font = font)
            vp = Plotter(N=3, axes=10)
            vp.show(kspl_opt1, sph_inf_o, sph_outf_o, myoc.alpha(0.01), txt1, at=0)
            vp.show(kspl_opt2, sph_inf_o, sph_outf_o, sph_m1, sph_m10, myoc.alpha(0.01), txt2, at=1)
            vp.show(kspl_opt3, sph_inf_o, sph_outf_o, sph_m1, sph_m30, myoc.alpha(0.01), txt3,  at=2, interactive=True)
            q_select = ask4input('Select the preferred centreline for processing this heart \n\t-[0]: Option 0 - left \n\t-[1]: Option 1 - centre \n\t-[2]: Option 2 - right \n\t-[3]: Define IFT last point manually >>:', int)
            
            dict_kspl['kspl_opt_sel'] = q_select
        
        # If selection is diff to 3 then select that CL
        if q_select !=3:
            kspl_CL.append(ksp_sel[q_select])
            pt2add_final = pt2add_sel[q_select]
            
        # If selection is 3
        else: 
            # See first if the process of cutting the mesh has already been done
            try: 
                # - Option 4 (add point obtained when creating mesh for centreline - vmtk)
                pt2add_inf4 = dict_kspl['linLine_'+layerNames[i]]['kspl_pts'][0]
                pt2add_inf5 = pts_withOutf[-1]+unit_vector(pt2add_inf4-pts_withOutf[-1])*(findDist(pt2add_inf4,pts_withOutf[-1])/2)
                sph_inf5 = Sphere(pos=pt2add_inf5, r=2, c='lightcoral').legend('new_interm_pt')
                sph_centroid = Sphere(pos=pt2add_inf4, r=2, c='deeppink').legend('new_pt')
                pts_all_opt4 = np.insert(pts_withOutf, len(pts_withOutf), np.transpose(pt2add_inf5), axis=0)
                pts_all_opt5 = np.insert(pts_all_opt4, len(pts_all_opt4), np.transpose(pt2add_inf4), axis=0)
                
                # Interpolate points
                pts_int_opt5 = getInterpolatedPts(points=pts_all_opt5, nPoints = 300)
                # Create kspline with points
                kspl_opt4 = KSpline(pts_int_opt5, res = 300)
                kspl_opt4.color(colors[i]).legend('CL_'+layerNames[i]).lw(5)
                
                vp = Plotter(N=1)
                text = "- Selected Centreline"
                txt = Text2D(text, c=c, font=font)
                vp.show(myoc, kspl_opt4, sph_centroid, sph_inf5, txt, at=0, viewup="y", interactive=True)
                
            # If not, then cut the mesh and define final point for centreline
            except: 
            # if True: 
                pt_happy = False
                
                mesh_title = filename+'_myoc_int_cut4cl.stl'
                mesh_dir = os.path.join(dir_stl, mesh_title)
                meshCL_cut = load(mesh_dir)
                meshCL_cut.color(' mediumvioletred').alpha(0.5)
                
                while not pt_happy:
                    planeCL_cut, plCL_cut_centre, plCL_cut_normal = getPlane(filename = filename +' - Redefine CL', type_cut = 'inflow',
                                                                        info = 'redefCL', mesh_in = meshCL_cut,
                                                                        mesh_out = kspl_o, dict_planes = dict_planes)
                    pts2cut, _ = getPointsAtPlane(points = meshCL_cut.points(), pl_normal = plCL_cut_normal,
                                              pl_centre = plCL_cut_centre)
                    ordpts, _ = order_pts(points = pts2cut)
    
                    # Create spline around cut
                    kspl = KSpline(ordpts, continuity=0, tension=0, bias=0, closed=True)
                    kspl.color('purple').lw(2)
    
                    # Get centroid of kspline to add to the centreline
                    kspl_bounds = kspl.bounds()
                    pt_centroid = np.mean(np.asarray(kspl_bounds).reshape((3, 2)),axis=1)
                    sph_centroid = Sphere(pos=pt_centroid, r=2, c='deeppink').legend('new_pt')
                    
                      # - Option 4 (add point obtained when re-cutting exported mesh4cl)
                    pt2add_inf4 = pt_centroid
                    pt2add_inf5 = pts_withOutf[-1]+unit_vector(pt_centroid-pts_withOutf[-1])*(findDist(pt_centroid,pts_withOutf[-1])/2)
                    # print(pts_withOutf[-1], pt2add_inf5, pt2add_inf4)
                    sph_inf5 = Sphere(pos=pt2add_inf5, r=2, c='lightcoral').legend('new_interm_pt')
                    pts_all_opt4 = np.insert(pts_withOutf, len(pts_withOutf), np.transpose(pt2add_inf5), axis=0)
                    pts_all_opt5 = np.insert(pts_all_opt4, len(pts_all_opt4), np.transpose(pt2add_inf4), axis=0)
            
                    # Interpolate points
                    pts_int_opt5 = getInterpolatedPts(points=pts_all_opt5, nPoints = 300)
                    # Create kspline with points
                    kspl_opt4 = KSpline(pts_int_opt5, res = 300)
                    kspl_opt4.color(colors[i]).legend('CL_'+layerNames[i]).lw(5)
                    
                    vp = Plotter(N=1)
                    text = "- Resulting mesh after cutting"
                    txt = Text2D(text, c=c, font=font)
                    vp.show(meshCL_cut, kspl, kspl_opt4, sph_centroid, sph_inf5, txt, at=0, viewup="y", interactive=True)
                    
                    pt_happy = ask4input('Are you happy with the new defined centreline? \n\t-[0]:no, I want to define a new cuting plane to redefine final point\n\t-[1]:yes! >>:', bool)
                    
            kspl_CL.append(kspl_opt4)
            pt2add_final = pt2add_inf4
            
        # Get linear lines
        linearLine = Line(pt2add_final, pt2add_outf, c=color_linLines[i], lw=3)
        linearLine.legend('linLine_'+layerNames[i])
        linLines.append(linearLine)

        dict_kspl = addKSplines2Dict(kspls = [linearLine, kspl_CL[-1]], info = ['',''], dict_kspl = dict_kspl)

    return kspl_CL, linLines, [kspl_o, sph_inf_o, sph_outf_o], dict_kspl

#%% func -createColouredCL
def createColouredCL(dict_cl, dict_shapes):
    """
    Function that creates the centrelines using the points given as input in the dict_cl

    Parameters
    ----------
    dict_cl : dictionary
        Initialised dictionary with centreline information.
    dict_shapes : dictionary
        Initialised dictionary with shapes information.

    Returns
    -------
    spheres_CL : Spheres
        Group of spheres with maximum inscribed centreline radius and positions given by VMTK and color coded by actual sphere radius.
    spheres_CL_col : Spheres
        Group of spheres with actual sphere radius in the positions given by VMTK.
    dict_shapes : dictionary
        Resulting dictionary with shapes information updated.

    """

    print('- Creating centreline associated objects...')

    spheres_CL = []
    spheres_CL_col = []

    layerNames = ['Int.Myoc(Cut)', 'Ext.Endo(Cut)']

    # Iterate through layer
    for i, dictLayer in enumerate(dict_cl):
        #Get cl points from vmtk
        pts_cl = np.asarray(dictLayer['Points'])

        #Get cl points from vmtk
        pts_cl = np.asarray(dictLayer['Points'])
        
        sphData_cl = np.asarray(dictLayer['PointData']['MaximumInscribedSphereRadius'])
        #Plot spheres through centreline inside mesh.
        vcols = [colorMap(v, "jet", sphData_cl.max(), sphData_cl.min()) for v in sphData_cl]
        sph_cl = Spheres(pts_cl, c='red', r=sphData_cl).legend('sphs_maxInsSphRad_'+layerNames[i])
        sph_cl_colour = Spheres(pts_cl, c=vcols, r=sphData_cl.min()).addScalarBar(title='Spheres Radius\n[um]').legend('sphs_maxInsSphRadC_'+layerNames[i])
        sph_cl_colour.mapper().SetScalarRange(sphData_cl.min(),sphData_cl.max())
        spheres_CL.append(sph_cl)
        spheres_CL_col.append(sph_cl_colour)

        dict_shapes = addShapes2Dict(shapes = [sph_cl, sph_cl_colour], dict_shapes = dict_shapes, radius = [[],sphData_cl])

    return spheres_CL, spheres_CL_col, dict_shapes

#%% func - createCLRibbon
def createCLRibbon(filename, file_num, df_res, kspl_CL2use, linLine, mesh, dict_kspl, dict_shapes, dict_planes, clRib_type = 'extDV', plotshow = True):
    """
    Function that creates dorso-ventral extended centreline ribbon

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    file_num : int
        Index number of the selected heart being processed.
    df_res : dataframe
        Dataframe with the measured information of the heart being processed.
    kspl_CL2use : Kspline
        Centreline (vedo KSpline)
    linLine : Line
        Line defining the linear length of the heart (inflow-outflow)
    mesh : mesh
        Mesh (vedo Mesh)
    dict_kspl :  dictionary
        Initialised dictionary with kspline information
    dict_shapes :  dictionary
        Initialised dictionary with shapes information
    dict_planes : dictionary
        Initialised dictionary with planes information
    clRib_type : str ['extDV', 'extV']
        String to define the centreline ribbon the function should return. If 'extDV': then a Dorso-Ventrally extended 
        centreline ribbon will be returned, else if 'extV', only a ventrally extended centreline ribbon with much 
        higher resolution will be returned.The default is 'extDV'.
    plotshow : boolean, optional
        True if you want to see the resulting mesh in a plot, else False. The default is True.

    Returns
    -------
    cl_ribbon : ribbon
        Dorso-ventral extended centreline (vedo Ribbon)
    kspl_ext : Kspline
        Kspline of extended centreline
    dict_kspl : dictionary
        Resulting dictionary with ksplines information updated
    dict_shapes : dictionary
        Resulting dictionary with shapes information updated
    dict_planes : dictionary
        Resulting dictionary with planes information updated

    """

    pts_cl = kspl_CL2use.points()
    # Extended centreline
    inf_ext_normal = (pts_cl[-1]+(pts_cl[-1]-pts_cl[-3])*70)
    outf_ext_normal = (pts_cl[0]+(pts_cl[0]-pts_cl[1])*70)
    inf_ext_sphere = Sphere(pos=inf_ext_normal, r=3, c='purple').legend("sph_infCLExtended")
    outf_ext_sphere = Sphere(pos=outf_ext_normal, r=3, c='purple').legend("sph_outfCLExtended")

    pts_cl_ext = np.insert(pts_cl,0,np.transpose(outf_ext_normal), axis=0)
    pts_cl_ext = np.insert(pts_cl_ext,len(pts_cl_ext),np.transpose(inf_ext_normal), axis=0)
    
    # Increase the resolution of the extended centreline and interpolate to unify sampling
    xd = np.diff(pts_cl_ext[:,0])
    yd = np.diff(pts_cl_ext[:,1])
    zd = np.diff(pts_cl_ext[:,2])
    dist = np.sqrt(xd**2+yd**2+zd**2)
    u = np.cumsum(dist)
    u = np.hstack([[0],u])
    t = np.linspace(0, u[-1], 601)
    resamp_pts = interpn((u,), pts_cl_ext, t)
    kspl_ext = KSpline(resamp_pts, res=601).color('purple').legend('kspl_extended')

    spaw_analysis = False
    if 'spaw' in df_res.loc[file_num,'spAnalysis']:
        spaw_analysis = True #ask4input('You are processing a heart that came from an incross of spaw heterozygous.\n  Please, select the way this heart is looping to continue processing: \n\t[0]: right-left\n\t[1]: dorso-ventral >>>: ', bool)
        
    # Create plane to project centreline
    dORv = filename[9:10]

    if dORv == 'D' or 'CJ' in filename: 
        azimuth = -90
        linLineX = linLine.clone().projectOnPlane('z').c(linLine.color()).z(0)
    elif spaw_analysis: 
        print('spaw')
        azimuth = 0
        linLineX = linLine.clone().projectOnPlane('z').c(linLine.color()).z(150)
    elif dORv == 'V':
        azimuth = 0
        linLineX = linLine.clone().projectOnPlane('x').c(linLine.color()).x(0)

    ptsPl_linLine = Points([linLineX.points()[0], linLine.points()[0], linLine.points()[1]])
    pl_linLine = fitPlane(ptsPl_linLine.points()).scale(4).c('mediumaquamarine').alpha(1).legend('pl_Parallel2LinLine')
    pl_linLine_normal = pl_linLine.normal
    pl_linLine_centre = pl_linLine.center
    dict_planes = addPlanes2Dict (planes = [pl_linLine], pls_centre = [pl_linLine_centre],
                                            pls_normal = [pl_linLine_normal], info = [''], dict_planes = dict_planes)

    pl_linLine_unitNormal = unit_vector(pl_linLine_normal)
    x_ul, y_ul, z_ul = pl_linLine_unitNormal*2
    x_ucl, y_ucl, z_ucl = pl_linLine_unitNormal*15
    pl_linLine_unitNormal120 = pl_linLine_unitNormal*120
    x_cl, y_cl, z_cl = pl_linLine_unitNormal120
    
    if clRib_type == 'extDV': # Names are switched but it works
        kspl_ext_D = kspl_ext.clone().x(x_cl).y(y_cl).z(z_cl).legend('kspl_CLExtD')
        kspl_ext_V = kspl_ext.clone().x(-x_cl).y(-y_cl).z(-z_cl).legend('kspl_CLExtV')
        cl_ribbon = Ribbon(kspl_ext_D, kspl_ext_V, alpha=0.2, res=(220, 5))
        cl_ribbon = cl_ribbon.wireframe(True).legend("rib_ExtCL(D-V)")
    
    elif clRib_type == 'extV':
        cl_ribbon = []
        for i in range(10):
            kspl_ext_DA = kspl_ext.clone().x(i*x_ucl).y(i*y_ucl).z(i*z_ucl)
            kspl_ext_DB = kspl_ext.clone().x((i+1)*x_ucl).y((i+1)*y_ucl).z((i+1)*z_ucl)
            cl_ribbon2un = Ribbon(kspl_ext_DA, kspl_ext_DB, alpha=0.2, res=(220, 5))
            cl_ribbon.append(cl_ribbon2un)
        cl_ribbon = merge(cl_ribbon)
        cl_ribbon.legend('rib_ExtCL(V)').wireframe(True)
    
    elif clRib_type == 'HDStack':
        cl_ribbon = []
        for i in range(100):
            kspl_ext_D = kspl_ext.clone().x(x_cl-i*x_ul).y(y_cl-i*y_ul).z(z_cl-i*z_ul)
            kspl_ext_V = kspl_ext.clone().x(-x_cl+i*x_ul).y(-y_cl+i*y_ul).z(-z_cl+i*z_ul)
            cl_ribbon2un = Ribbon(kspl_ext_D, kspl_ext_V, alpha=0.2, res=(220, 20))
            if i == 0: 
                rib_pts = cl_ribbon2un.points()
            else: 
                rib_pts = np.concatenate((rib_pts,cl_ribbon2un.points()))
            cl_ribbon.append(cl_ribbon2un)
        cl_ribbon = merge(cl_ribbon)
        cl_ribbon.legend('HDStack').wireframe(True)
        
    dict_shapes = addShapes2Dict (shapes = [cl_ribbon], dict_shapes = dict_shapes, radius = [[]])
    if clRib_type == 'extDV':
        dict_kspl = addKSplines2Dict(kspls = [kspl_ext_D, kspl_ext, kspl_ext_V], info = ['','', ''], dict_kspl = dict_kspl)

    if plotshow:
        text = filename+"\n\n >> Creating Extended Centreline ("+clRib_type+")"
        txt = Text2D(text, c="k", font= font)
    
        elevation = 0#df_res.loc[file_num,'ang_Heart']
        settings.legendSize = .3
        vp = Plotter(N=1, axes=1)
        vp.show(mesh, linLine, linLineX, kspl_CL2use, kspl_ext, inf_ext_sphere, outf_ext_sphere, cl_ribbon, txt, at=0, azimuth = azimuth, elevation = elevation, interactive=1)

    if clRib_type == 'HDStack':
        return rib_pts, cl_ribbon
    else: 
        return cl_ribbon, kspl_ext, dict_kspl, dict_shapes, dict_planes

#%% - MEASURE
#%% func - getChambersOrientation
def getChambersOrientation(filename, file_num, num_pt, kspl_CL2use, distFromCl, myoc_meshes, linLine, dict_pts, dict_kspl, df_res):

    print('- Measuring chamber orientations')
    m_atrMyoc, m_ventMyoc = myoc_meshes
    atr_pt = num_pt + distFromCl
    vent_pt = num_pt - distFromCl

    # Create spheres
    sph_atr = Sphere(pos = kspl_CL2use.points()[atr_pt], r=4, c='navy').legend('sph_AtrCentre')
    sph_vent = Sphere(pos = kspl_CL2use.points()[vent_pt], r=4, c='red').legend('sph_VentCentre')
    sph_inf = Sphere(pos = kspl_CL2use.points()[-1], r=4, c='dodgerblue').legend('sph_OutflowCentre')
    sph_outf = Sphere(pos = kspl_CL2use.points()[0], r=4, c='tomato').legend('sph_InflowCentre')
    sph_valve = Sphere(pos = kspl_CL2use.points()[num_pt], r=4, c='darkorange').legend('sph_AVCCentre')

    sph_orient = [sph_atr, sph_vent, sph_outf, sph_inf, sph_valve]
    dict_pts = addPoints2Dict(spheres = sph_orient, info = ['Orient','Orient','Orient','Orient','Orient'], dict_pts = dict_pts)

    # >> Heart Orientation
    # Find orientation in which images where taken
    dORv = filename[9:10]
    if dORv == 'D' or 'CJ' in filename:
        side_plane = 'z'
        ventral_plane = 'x'
        azimuth = -135; elevation = 0

    elif dORv == 'V':# or spaw_analysis == False: 
        side_plane = 'x'
        ventral_plane = 'z'
        azimuth = -45; elevation = 0
        
    # FROM THE SIDE
    linLineX = linLine.clone().projectOnPlane(side_plane).c(linLine.color()).legend('linLine(ProjX)')
    # [0: v_head = OFT, 1: v_tail = IFT]
    pts_heartS =  orientVectors(linLineX) 
    # print('pts_heartS', pts_heartS)
    ang_heartS = findAngleBtwVectorsZ(pts_heartS, np.array([[0,1,0],[0,0,0]])) 
    # rotate the heart orientation line to confirm direction 
    linLineXRot = linLineX.clone().rotateX(ang_heartS).color('gold')
    ang_linLineXRot =  findAngleBtwVectorsZ(linLineXRot.points(), np.array([[0,0,0],[0,1,0]]))
    
    if ang_linLineXRot != 0:
        ang_heartS = -ang_heartS
        linLineXRot = linLineX.clone().rotateX(ang_heartS).color('gold')
        ang_linLineXRot =  findAngleBtwVectorsZ(linLineXRot.points(), np.array([[0,0,0],[0,1,0]]))
        print('In to if')
        
    if ang_linLineXRot == 0:
        hrt_txt = 'Heart orientation from the side (deg): '+ format(ang_heartS,'.1f')
        print('>> Orientation from the side\n\t>> '+hrt_txt)
    else: 
        print('Check!')
    
    # FROM THE FRONT/VENTRALLY
    if dORv == 'D' or 'CJ' in filename:
        linLineRot = linLine.clone().rotateZ(ang_heartS).color('forestgreen')
    elif dORv == 'V':# or spaw_analysis == False: 
        linLineRot = linLine.clone().rotateX(ang_heartS).color('forestgreen')
        
    linLineRotProj = linLineRot.clone().projectOnPlane(ventral_plane).color('springgreen')
    # [0: v_head = OFT, 1: v_tail = IFT]
    pts_heartV =  orientVectors(linLineRotProj) 
        
    # >> Atrial orientation
    orient_atr = Line(sph_inf.pos(), sph_atr.pos(), c="steelblue", lw=3).legend('lin_OrientAtr')
    orient_atrX = orient_atr.clone().projectOnPlane(side_plane).c('steelblue').legend('lin_OrientAtr(ProjX)')
    if dORv == 'D' or 'CJ' in filename:
        orient_atrRot = orient_atr.clone().rotateZ(ang_heartS).color('navy')
    elif dORv == 'V':# or spaw_analysis == False:
        orient_atrRot = orient_atr.clone().rotateX(ang_heartS).color('navy')
        
    orient_atrRotProj = orient_atrRot.clone().projectOnPlane(ventral_plane).color('turquoise')
    
    # FROM THE SIDE
    # point closer to v_head (OFT) goes first, to make sure it is in the same orientation as pts_heartS
    pts_atrS = orientVectors(orient_atrX, ref_pt = pts_heartS[0])
    ang_atrS = findAngleBtwVectorsZ(pts_atrS, pts_heartS)
    atr_txtS = 'Atrial orientation with respect to heart from the side (deg): '+ format(ang_atrS,'.1f')
    print('\t>> '+atr_txtS)
    # FROM THE FRONT/VENTRALLY
    # point closer to v_head (OFT) goes first, to make sure it is in the same orientation as pts_heartV
    pts_atrV = orientVectors(orient_atrRotProj, ref_pt = pts_heartV[0]) 
    ang_atrV = findAngleBtwVectorsZ(pts_atrV, pts_heartV)
    atr_txtV = 'Atrial orientation with respect to heart ventrally (deg): '+ format(ang_atrV,'.1f')
    
    # >> Ventricular orientation
    orient_vent = Line(sph_vent.pos(), sph_outf.pos(), c="hotpink", lw=3).legend('lin_OrientVent')
    orient_ventX = orient_vent.clone().projectOnPlane(side_plane).c('hotpink').legend('lin_OrientVent(ProjX)')

    if dORv == 'D' or 'CJ' in filename:
        orient_ventRot = orient_vent.clone().rotateZ(ang_heartS).color('crimson')
    elif dORv == 'V':# or spaw_analysis == False:
        orient_ventRot = orient_vent.clone().rotateX(ang_heartS).color('crimson')
    
    orient_ventRotProj = orient_ventRot.clone().projectOnPlane(ventral_plane).color('tomato')
    
    # FROM THE SIDE
    # point closer to v_tail (IFT) goes first, to make sure it is in the inverse orientation as pts_heartS
    pts_ventS = orientVectors(orient_ventX, ref_pt = pts_heartS[1]) 
    ang_ventS = findAngleBtwVectorsZ(pts_ventS, pts_heartS)
    vent_txtS = 'Ventricular orientation with respect to heart from the side (deg): '+ format(ang_ventS,'.1f')
    print('\t>> '+vent_txtS)
    # FROM THE FRONT/VENTRALLY
    # point closer to v_tail (IFT) goes first, to make sure it is in the inverse orientation as pts_heartS
    pts_ventV = orientVectors(orient_ventRotProj, ref_pt = pts_heartV[1])
    ang_ventV = findAngleBtwVectorsZ(pts_ventV, pts_heartV)
    vent_txtV = 'Ventricular orientation with respect to heart ventrally (deg): '+ format(ang_ventV,'.1f')

    # >> Angle between chambers
    # FROM THE SIDE
    ang_chsS = findAngleBtwVectorsZ(pts_atrS, pts_ventS)
    chs_txtS = 'Angle between chambers from the side (deg): '+ format(ang_chsS,'.1f')
    print('\t>> '+chs_txtS)
    # FROM THE FRONT/VENTRALLY
    ang_chsV =  findAngleBtwVectorsZ(pts_atrV, pts_ventV)
    chs_txtV = 'Angle between chambers ventrally (deg): '+ format(ang_chsV,'.1f')
    
    print('>> Ventral orientation \n\t>> '+atr_txtV)
    print('\t>> '+vent_txtV)
    print('\t>> '+chs_txtV)
    
    lines_orient = [orient_atr, orient_vent, orient_atrX, orient_ventX, linLineX]
    dict_kspl = addKSplines2Dict(kspls = lines_orient, info = '', dict_kspl = dict_kspl)
    
    df_res = addOrientationAngles2df(df_res = df_res, file_num = file_num, 
                                     angles = [ang_heartS, ang_atrS, ang_ventS, ang_chsS]+[ang_atrV, ang_ventV, ang_chsV], 
                                     names = ['ang_HeartS', 'ang_AtrS', 'ang_VentS', 'ang_BtwChambersS']+['ang_AtrV', 'ang_VentV', 'ang_BtwChambersV'])

    z_maxpos = pts_heartS[pts_heartS[:,2].argsort()][-1,-1]+50
    
    text0 = filename+"\n\n >> Measuring chamber orientation from the side\n - "+atr_txtS+"\n - "+vent_txtS+"\n - "+chs_txtS+"\n - "+hrt_txt
    text1 = filename+"\n\n >> Measuring chamber orientation ventrally \n - "+atr_txtV+"\n - "+vent_txtV+"\n - "+chs_txtV
    txt0 = Text2D(text0, c="k", font= font)
    txt1 = Text2D(text1, c="k", font= font)
    settings.legendSize = .15
    vp = Plotter(N=2, axes = 8)
    vp.show(m_atrMyoc.alpha(0.01), m_ventMyoc.alpha(0.01), sph_orient, kspl_CL2use, orient_atr, orient_vent, linLine,
            orient_atrX.alpha(1), orient_ventX.alpha(1), linLineX.alpha(1), txt0, 
            # azimuth = azimuth, elevation = elevation, 
            zoom = 0.8, at=0)
    vp.show(m_atrMyoc.alpha(0.01), m_ventMyoc.alpha(0.01), sph_orient, kspl_CL2use, orient_atr, orient_vent, linLine,
            orient_atrRot.alpha(1), orient_ventRot.alpha(1), linLineRot.alpha(1),
            orient_atrRotProj.alpha(1).z(z_maxpos), orient_ventRotProj.alpha(1).z(z_maxpos), linLineRotProj.alpha(1).z(z_maxpos), txt1,
            zoom = 0.8, at=1, azimuth = azimuth, elevation = elevation, interactive=True)

    return sph_orient, lines_orient, dict_pts, dict_kspl, df_res

#%% func - getChambersEllipsoid
def getChambersEllipsoid(filename, df_res, file_num, lines_orient, meshes, dict_shapes, plot_show = True):
    
    df_resFilled = df_res.copy()
    orient_atr, orient_vent, _, _, _ = lines_orient
    m_atr, m_vent = meshes
    fill_df = False
    
    rotAngle =  df_res.loc[file_num,'ang_HeartS']
    atrAngle = df_res.loc[file_num,'ang_AtrS']
    ventAngle = df_res.loc[file_num,'ang_VentS']
    
    text = filename+'\n\n >> Fitting best ellipse into '
    
    # Find orientation in which images where taken or if it is a spaw mutant
    spaw_analysis = False
    if 'spaw' in df_res.loc[file_num,'spAnalysis']:
        spaw_analysis = True #ask4input('You are processing a heart that came from an incross of spaw heterozygous.\n  Please, select the way this heart is looping to continue processing: \n\t[0]: right-left\n\t[1]: dorso-ventral >>>: ', bool)
    
    # Heart Orientation
    dORv = filename[9:10]
    for m_ch in meshes: 
        ch_legend = m_ch._legend
        if  'Atr' in ch_legend: 
            chamber = 'Atrium'
            orient_ch = orient_atr.clone()
            s_name = 'Atr'
            
        elif 'Vent' in ch_legend:
            chamber = 'Ventricle'
            orient_ch = orient_vent.clone()
            s_name = 'Vent'
        print('- Reorienting ',chamber,'...')
        
        if dORv == 'D' or 'CJ' in filename:
            print('\t-> Dorsal or CJ')
            pts_o = orient_ch.clone().points()
            # print(pts_o,'\n\n', orient_ch.points())
            pts_o[0,2] = 0; pts_o[1,2] = 0
            ang_test = findAngleBtwVectorsZ(pts_o, np.array([[0,0,0],[0,1,0]]))
            rotX = ang_test
            print('\t- D rotX:', rotX, ' -', chamber)
            
            orient_chRot = orient_ch.clone().rotateZ(rotX).color('gold')
            ang_chX = findAngleBtwVectorsZ(orient_chRot.points(), np.array([[0,0,0],[0,1,0]]))
            orient_chRotRot = orient_chRot.clone().rotateX(ang_chX).color('navy')
            ang_chF = findAngleBtwVectorsZ(orient_chRotRot.points(), np.array([[0,0,0],[0,1,0]]))
            pl_rotX = 'z'
            pl_rotZ = 'x'
            
            if ang_chF != 0:
                 ang_chX = -ang_chX
                 orient_chRotRot = orient_chRot.clone().rotateX(ang_chX).color('navy')
                 ang_chF = findAngleBtwVectorsZ(orient_chRotRot.points(), np.array([[0,0,0],[0,1,0]]))
                 print('\t- Note: in D if')
            
            if ang_chF == 0:
                #Rotate mesh in X and Z
                m_chRot = m_ch.clone().rotateZ(rotX).rotateX(ang_chX)
                
                m_chRot_projX = m_chRot.clone().projectOnPlane(plane = pl_rotX).c('palegreen').alpha(0.05)
                m_chRot_projZ = m_chRot.clone().projectOnPlane(plane = pl_rotZ).c('skyblue').alpha(0.05)
                xa_min, xa_max, ya_min, ya_max, za_min, za_max = m_chRot_projX.bounds()
                _, _, _, _,  za_min, za_max = m_chRot_projZ.bounds()
                xa_pos, ya_pos, _ = m_chRot_projX.centerOfMass()
                _, _, za_pos = m_chRot_projZ.centerOfMass()
                
                ellip_ch = Ellipsoid(pos=(xa_pos, ya_pos, za_pos), axis1=(xa_max-xa_min, 0,0), axis2=(0, ya_max-ya_min, 0), axis3=(0, 0, za_max-za_min), c='salmon', alpha=0.5, res=24)
                dict_shapes['Ellipsoid_'+s_name] = {'pos': (xa_pos, ya_pos, za_pos), 'axis1': (xa_max-xa_min, 0,0), 
                                                    'axis2': (0, ya_max-ya_min, 0), 'axis3': (0, 0, za_max-za_min)}
                
                width_ch = za_max-za_min
                length_ch = ya_max-ya_min
                depth_ch = xa_max-xa_min
                fill_df = True
                
                txt_ch = Text2D(text+chamber+' -0 angF: {:.1f}'.format(ang_chF), c="k", font= font)
                txt_ch2 = Text2D(chamber+': \n\t Width: {:.1f} um, Length: {:.1f} um, Depth: {:.1f} um'.format(width_ch, length_ch, depth_ch), c="k", font= font)
                
                vp = Plotter(N=2, axes = 1)
                vp.show(m_chRot, m_chRot_projX, m_chRot_projZ, txt_ch, at = 0)
                vp.show(ellip_ch, m_chRot, txt_ch2, at=1, interactive = True)
                
        elif dORv == 'V':# or spaw_analysis == False: 
            print('\t-> Ventral')
            # if not spaw_analysis:
            if chamber == 'Atrium':
                rotX = rotAngle-atrAngle
            elif chamber == 'Ventricle':
                rotX = rotAngle+(180-ventAngle)
                print(' \t- rotX:', rotX)
            # else: 
            #     pts_o = orient_ch.clone().points()
            #     # print(pts_o,'\n\n', orient_ch.points())
            #     pts_o[0,0] = 0; pts_o[1,0] = 0
            #     ang_test = findAngleBtwVectorsZ(pts_o, np.array([[0,0,0],[0,1,0]]))
            #     rotX = ang_test  
            #     print('\t- Note: in V if spaw')
                # print('\t- rotX:', ang_test)
                
            print('\t- V rotX:', rotX, '-', chamber)
            orient_chRot = orient_ch.clone().rotateX(rotX).color('gold')
            ang_chX =  findAngleBtwVectorsZ(orient_chRot.points(), np.array([[0,0,0],[0,1,0]]))
            orient_chRotRot = orient_chRot.clone().rotateZ(ang_chX).color('navy')
            ang_chF =  findAngleBtwVectorsZ(orient_chRotRot.points(), np.array([[0,0,0],[0,1,0]]))
            pl_rotX = 'x'
            pl_rotZ = 'z'
            print('\t- ang_chX:', ang_chX)
            
            if ang_chF != 0:
                print('\t- Note: in V if ang_chF != 0 (opt1)')
                ang_chX = -ang_chX
                print('\t- ang_chX:', ang_chX)
                orient_chRotRot = orient_chRot.clone().rotateZ(ang_chX).color('navy')
                ang_chF = findAngleBtwVectorsZ(orient_chRotRot.points(), np.array([[0,0,0],[0,1,0]]))
            
            # if still angle is not 0
            if ang_chF !=0 and not spaw_analysis:
                pts_o = orient_ch.clone().points()
                pts_o[0,0] = 0; pts_o[1,0] = 0
                ang_test = findAngleBtwVectorsZ(pts_o, np.array([[0,0,0],[0,1,0]]))
                rotX = ang_test  
                print('\t- Note: in V if not spaw but ang_chF !=0')
                
                print('\t- V rotX:', rotX, '-', chamber)
                orient_chRot = orient_ch.clone().rotateX(rotX).color('gold')
                ang_chX =  findAngleBtwVectorsZ(orient_chRot.points(), np.array([[0,0,0],[0,1,0]]))
                orient_chRotRot = orient_chRot.clone().rotateZ(ang_chX).color('navy')
                ang_chF =  findAngleBtwVectorsZ(orient_chRotRot.points(), np.array([[0,0,0],[0,1,0]]))
                pl_rotX = 'x'
                pl_rotZ = 'z'
                print('\t- ang_chX:', ang_chX)
                
                if ang_chF != 0:
                    print('\t- Note: in V if ang_chF != 0 (opt2)')
                    ang_chX = -ang_chX
                    print('\t- ang_chX:', ang_chX)
                    orient_chRotRot = orient_chRot.clone().rotateZ(ang_chX).color('navy')
                    ang_chF = findAngleBtwVectorsZ(orient_chRotRot.points(), np.array([[0,0,0],[0,1,0]]))
            
            if ang_chF == 0:
                #Rotate mesh in X and Z
                m_chRot = m_ch.clone().rotateX(rotX).rotateZ(ang_chX)
                
                m_chRot_projX = m_chRot.clone().projectOnPlane(plane = pl_rotX).c('palegreen').alpha(0.05)
                m_chRot_projZ = m_chRot.clone().projectOnPlane(plane = pl_rotZ).c('skyblue').alpha(0.05)
                _, _, ya_min, ya_max, za_min, za_max = m_chRot_projX.bounds()
                xa_min, xa_max, _, _, _, _ = m_chRot_projZ.bounds()
                _, ya_pos, za_pos = m_chRot_projX.centerOfMass()
                xa_pos, _, _ = m_chRot_projZ.centerOfMass()
                
                ellip_ch = Ellipsoid(pos=(xa_pos, ya_pos, za_pos), axis1=(xa_max-xa_min, 0,0), axis2=(0, ya_max-ya_min, 0), axis3=(0, 0, za_max-za_min), c='salmon', alpha=0.5, res=24)
                dict_shapes['Ellipsoid_'+s_name] = {'pos': (xa_pos, ya_pos, za_pos), 'axis1': (xa_max-xa_min, 0,0), 
                                                    'axis2': (0, ya_max-ya_min, 0), 'axis3': (0, 0, za_max-za_min)}
                
                width_ch = xa_max-xa_min
                length_ch = ya_max-ya_min
                depth_ch = za_max-za_min
                fill_df = True

                txt_ch = Text2D(text+chamber+' - angF: {:.1f}'.format(ang_chF), c="k", font= font)
                txt_ch2 = Text2D(chamber+': \n\t Width: {:.1f} um, Length: {:.1f} um, Depth: {:.1f} um'.format(width_ch, length_ch, depth_ch), c="k", font= font)
                
                vp = Plotter(N=2, axes = 1)
                vp.show(m_chRot, m_chRot_projX, m_chRot_projZ, txt_ch, at = 0)
                vp.show(ellip_ch, m_chRot,  orient_ch, orient_chRot, orient_chRotRot, txt_ch2, at=1, interactive = True)
            
            else: 
                print('-ERROR: Something went wrong reorienting the '+chamber+'!')
                txt_ch = Text2D(text+chamber+' - angF: {:.1f}'.format(ang_chF), c="k", font= font)
                vp = Plotter(N=1, axes = 1)
                vp.show(m_ch, orient_ch, orient_chRot, orient_chRotRot, txt_ch, at = 0, interactive = True)
            
        if fill_df: 
            df_resFilled.loc[file_num, 'Ellip'+s_name+'_Width'] = width_ch
            df_resFilled.loc[file_num, 'Ellip'+s_name+'_Length'] = length_ch
            df_resFilled.loc[file_num, 'Ellip'+s_name+'_Depth'] = depth_ch
            df_resFilled.loc[file_num, 'Ellip'+s_name+'_Asphericity'] = ellip_ch.asphericity()
            
    if fill_df: 
        alert('wohoo', 1)
        print('- All ellipsoid parameters of the imported chamber(s) have been saved in df_res and dict_shapes!')
            
    return df_resFilled, dict_shapes
   
#%% func - getDistance2Mesh
def getDistance2Mesh(filename, m_int, m_ext, title, alpha = 1, plotshow = True):
    """
    Function that gets the distance between m_ext and m_int and color codes m_ext accordingly

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    m_int : Mesh/Object
        Internal mesh/object (vedo Mesh/Object)
    m_ext : Mesh
        External mesh (vedo Mesh)
    title : str
        Name given to the distance being calculated
    alpha : float, optional
        Opacity value. The default is 0.1.
    plotshow : boolean, optional
        True if you want to see the resulting mesh in a plot, else False. The default is True.

    Returns
    -------
    thickness : array of floats
        Numpy array with distance being measured
    m_ext_out : Mesh
        Final external mesh color coded (vedo Mesh)
    min_max : list of floats
        List with [min,max] values of the measured distance

    """

    thickness_meshes = ['Cardiac Jelly Thickness', 'Myoc.Thickness', 'Endo.Thickness']
    tic = perf_counter()
    print('- Extracting '+title)
    # print("  > Start time: \t", str(datetime.now())[11:-7])
    m_ext_out = m_ext.clone()
    m_ext_out.distanceToMesh(m_int, signed=True, negate=False)
    if title in thickness_meshes:
        thickness = - m_ext_out.getPointArray("Distance")
    else:
        thickness = m_ext_out.getPointArray("Distance")

    vmin, vmax = np.min(thickness),np.max(thickness)
    toc = perf_counter()
    time = toc-tic
    # print("  > End time: \t\t", str(datetime.now())[11:-7])
    print("\t>> Total time taken to get heatmap = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")

    alert('jump',1)

    val_min = vmin; val_max = vmax

    title_split = title.split(' ')
    title_bar = title_split[0]
    for txt in title_split[1:]:
        title_bar = title_bar +'\n'+ txt

    m_ext_out.pointColors(thickness, cmap="jet", vmin=val_min, vmax=val_max)
    m_ext_out.mapper().SetScalarRange(val_min,val_max)
    m_ext_out.addScalarBar(title=title_bar+' [um]', pos=(0.8, 0.05))
    m_ext_out.mapper().SetScalarRange(val_min,val_max)

    m_ext_out.alpha(alpha).legend(m_ext._legend)
    m_int.color("white").alpha(1).wireframe()

    if plotshow:
        text = filename+"\n\n >>"+title+" Heat map\n\t   - min and max: ("+format(vmin, '.2f')+" , "+format(vmax, '.2f')+")"
        txt = Text2D(text, c="k", font= font)

        cube = Cube(pos=m_ext_out.centerOfMass(), side=300, c='white', alpha=0.01)
        settings.legendSize = .3
        vp = Plotter(N=1, axes=10)
        vp.show(m_ext_out, cube, txt, at=0, zoom=1.2)#vp.show(m_ext_out, m_int, cube, txt, at=0, zoom=1.2)

    min_max = (vmin, vmax)

    return thickness, m_ext_out, min_max

#%% func - getExtCLonSurf
def getExtCLonSurf(filename, mesh, kspl_ext, pl_CLRibbon, process_plotshow = True):
    
    # - Get unitary normal of plane to create CL_ribbon
    pl_normCLRibbon = unit_vector(pl_CLRibbon['pl_normal'])
    
    # Increase the resolution of the extended centreline and interpolate to unify sampling
    xd = np.diff(kspl_ext.points()[:,0])
    yd = np.diff(kspl_ext.points()[:,1])
    zd = np.diff(kspl_ext.points()[:,2])
    dist = np.sqrt(xd**2+yd**2+zd**2)
    u = np.cumsum(dist)
    u = np.hstack([[0],u])
    t = np.linspace(0, u[-1],1000)
    resamp_pts = interpn((u,), kspl_ext.points(), t)
    kspl_ext = KSpline(resamp_pts, res = 1000).lw(5).color('deeppink').legend('kspl_extHR')
    
    #Find the points that intersect with the ribbon
    pts_int = []
    for num in range(len(kspl_ext.points())):
        try: 
            cl_pt_test = kspl_ext.points()[num]
            pt_int = mesh.intersectWithLine(cl_pt_test, cl_pt_test+150*pl_normCLRibbon)
            rad_pts = [np.linalg.norm(x- cl_pt_test) for x in pt_int]
            # if len(rad_pts)>1:
                # print(rad_pts)
            ind_pt = np.where(rad_pts == max(rad_pts))[0][0]
            # print(ind_pt)
            pts_int.append(pt_int[ind_pt])

        except: 
            # print('exc')
            if num > 750:
                dist_pts = kspl_ext.points()[num] - kspl_ext.points()[num-1]
                try: 
                    pt_int = pts_int[-1]+dist_pts
                    pts_int.append(pt_int)
                except:
                    # print('pass')
                    pass
            else: 
                pass
    
    # KSpline on surface
    kspl_vSurf = KSpline(pts_int).color('black').lw(4).legend('kspl_VSurfaceIntMyoc')
    
    if process_plotshow: 
        vp = Plotter(N=1, axes=4)
        vp.show(kspl_vSurf, kspl_ext, mesh, at = 0, interactive = True)
    
    return kspl_vSurf

#%% func - getExtCLHighRes
def getExtCLHighRes(filename, mesh, kspl_ext, kspl_CL, dict_planes, process_plotshow = True):
    
    inv = [True, False]
    # Create kspline_extended with higher resolution, but cut it using the inflow and outflow  planes defined to 
    # cut inf anf outf tracts of hearts to extract centreline
    
    add_pts = 50
    try: 
        # print('-CLHR Opt A')
        kspl_CLnew = kspl_ext.clone()
        # Cut with inflow plane
        n_points_In = -10
        for invert in inv:
            kspl_test = kspl_ext.clone().cutWithMesh(Plane(pos=dict_planes['pl2CutMesh_inflow']['pl_centre'], normal=dict_planes['pl2CutMesh_inflow']['pl_normal'], sx=300), invert=invert)
            # print('kspl_test In',kspl_test.NPoints())
            if kspl_test.NPoints() > n_points_In: 
                # print('in-In')
                inv_fin_In = invert
                kspl_CLnew_cutIn = kspl_ext.clone().cutWithMesh(Plane(pos=dict_planes['pl2CutMesh_inflow']['pl_centre'], normal=dict_planes['pl2CutMesh_inflow']['pl_normal'], sx=300), invert=inv_fin_In).color('darkorange')
                n_points_In = kspl_CLnew_cutIn.NPoints() 
        _, num_pt_inf = findClosestPtGuess(kspl_CLnew_cutIn.points(), kspl_ext.points(), index_guess = -1)
    except: 
        num_pt_inf = kspl_ext.NPoints()-1
    
    try: 
        # print('-CLHR Opt B')
        n_points_Out = -10
        for invert in inv:
            kspl_test = kspl_ext.clone().cutWithMesh(Plane(pos=dict_planes['pl2CutMesh_outflow']['pl_centre'], normal=dict_planes['pl2CutMesh_outflow']['pl_normal'], sx=300), invert=invert)
            # print('kspl_test Out',kspl_test.NPoints())
            if kspl_test.NPoints() > n_points_Out: 
                # print('in-Out')
                inv_fin_Out = invert
                kspl_CLnew_cutOut = kspl_ext.clone().cutWithMesh(Plane(pos=dict_planes['pl2CutMesh_outflow']['pl_centre'], normal=dict_planes['pl2CutMesh_outflow']['pl_normal'], sx=300), invert=inv_fin_Out).color('lime')
                n_points_Out = kspl_CLnew_cutOut.NPoints() 
        # _, num_pt_outf =  findClosestPt(kspl_CLnew_cutOut.points()[-1], kspl_ext.points())#, index_guess = 0)
        _, num_pt_outf =  findClosestPtGuess(kspl_CLnew_cutOut.points(), kspl_ext.points(), index_guess = 0)
    except: 
        num_pt_outf = 0
    # print('Initial def:',num_pt_outf, num_pt_inf)
    
    # Define starting point of new kspline
    if (num_pt_outf-add_pts) < 0:
        ind_outf = 0
    else: 
        ind_outf = num_pt_outf - add_pts
        
    # Define ending point of new kspline
    if (num_pt_inf+add_pts) > len(kspl_ext.points()):
        ind_inf = len(kspl_ext.points())
    else: 
        ind_inf = num_pt_inf+add_pts
    # print('Final def:',ind_outf, ind_inf)

    # Create new extended and cut kspline with higher resolution
    kspl_CLnew = KSpline(kspl_ext.points()[ind_outf:ind_inf], res = 600).lw(5).color('deeppink').legend(kspl_CL._legend+'_HighRes')
    
    if process_plotshow: 
        vp = Plotter(N=3, axes=4)
        vp.show(kspl_CLnew, kspl_ext, mesh, at = 0)
        vp.show(kspl_CLnew_cutIn, mesh, at = 1)
        vp.show(kspl_CLnew_cutOut, mesh, at = 2, interactive = True)
    
    return kspl_CLnew

#%% func - ksplChamberCut
def ksplChamberCut(mesh, kspl_CLnew, dict_shapes, dict_pts, process_plotshow = True): 
    
    cyl_chamber = dict_shapes['cyl2CutChambers_final']
    disc = Cylinder(pos = cyl_chamber['cyl_centre'],r = cyl_chamber['radius_max'], height = 2*0.225, 
                    axis = cyl_chamber['cyl_axis'], c = 'purple', cap = True, res = 300)
    
    # Find position within new kspline cut by disc used to cut chambers (num_pt new)
    ksplCL_cut = kspl_CLnew.clone().cutWithMesh(Plane(pos=dict_shapes['cyl2CutChambers_final']['cyl_centre'], normal=dict_shapes['cyl2CutChambers_final']['cyl_axis'], sx=300), invert=True)
    # Find point of new kspline closer to last point of kspline cut
    # print('-ChCut')
    _, num_pt = findClosestPtGuess(ksplCL_cut.points(), kspl_CLnew.points(),350)
    # print('num_pt',num_pt)
    
    # Add pt to dict
    sph_cut = Sphere(pos = kspl_CLnew.points()[num_pt], r=4, c='gold').legend('sph_ChamberCutCLHighRes')
    dict_pts = addPoints2Dict(spheres = [sph_cut], info = [''], dict_pts = dict_pts)
    dict_pts['numPt_CLChamberCutHighRes'] = num_pt
    
    # Create kspline for each chamber
    n_vent_pts = 0
    res_v = 595
    while n_vent_pts < 610:
        res_v +=1
        kspl_vent_pts = kspl_CLnew.points()[0:num_pt]
        kspl_vent = KSpline(points = kspl_vent_pts[::-1], res = res_v).color('tomato').lw(15).legend('kspl_vent')
        # kspl_vent = KSpline(points = kspl_CLnew.points()[0:num_pt], res = res_v).color('tomato').lw(8).legend('kspl_vent')
        n_vent_pts = kspl_vent.NPoints()
        
    n_atr_pts = 0
    res_a = 595
    while n_atr_pts < 610:
        res_a +=1
        kspl_atr = KSpline(points = kspl_CLnew.points()[num_pt:], res = res_a).color('goldenrod').lw(15).legend('kspl_atr')
        n_atr_pts = kspl_atr.NPoints()
    

        
    if process_plotshow: 
        vp = Plotter(N=2, axes=4)
        vp.show(kspl_CLnew, kspl_vent, sph_cut, disc, mesh, at = 0)
        vp.show(kspl_CLnew, kspl_atr, mesh, sph_cut, disc, at = 1, interactive = True)
        
    return kspl_vent, kspl_atr, dict_pts, sph_cut

#%% func - unloopChambers
def unloopChambers(filename, mesh, kspl_CL, kspl_ext, no_planes, pl_CLRibbon, param, param_name, df_AtrVent, selected_chambers, 
                   dict_kspl, dict_shapes, dict_pts, dict_planes, dir_results, save_names, plotshow = True, tol = 0.05, print_txt = False):
   
    print('\n\n- Unlooping the heart chambers...')
    dfs_unlooped = []
    spheres_zeroDeg = []
    arr_vectZeroDeg = []
    
    if plotshow:
        plotevery = no_planes // 5
        print('- Plotting every X number of planes:', plotevery)
    
    # KSpline on surface
    kspl_vSurf = getExtCLonSurf(filename, mesh, kspl_ext, pl_CLRibbon, False)
    # print(len(kspl_vSurf.points()))

    # Create new extended and cut kspline with higher resolution
    kspl_CLnew = getExtCLHighRes(filename, mesh, kspl_ext, kspl_CL, dict_planes, False) 

    # # Create kspline for each chamber
    kspl_vent, kspl_atr, dict_pts, sph_cut =  ksplChamberCut(mesh, kspl_CLnew, dict_shapes, dict_pts, False)
    # Find position within new kspline cut by disc used to cut chambers (num_pt new)
    num_pt_AVC = dict_pts['numPt_CLChamberCutHighRes']
    
    #Add all ksplines created to dictionary
    dict_kspl = addKSplines2Dict(kspls = [kspl_ext, kspl_vSurf, kspl_atr, kspl_vent], info = ['','','',''], dict_kspl = dict_kspl)
    
    if plotshow: 
        vp = Plotter(N=3, axes=4)
        vp.show(kspl_CLnew, kspl_vSurf, kspl_ext, mesh, sph_cut, at = 0)
        # vp.show(ksplCL_cut, arr_vectPlCut, kspl_vSurf, mesh, sph_cut, at = 0)
        vp.show(mesh, kspl_vent, sph_cut, at=1)
        vp.show(mesh, kspl_atr, sph_cut, at=2, interactive = True)
        
    # colours_sphs = range(no_planes+2)
    # NOW UNLOOP EACH CHAMBER
    kspls_HR = [kspl_atr, kspl_vent]
    texts = ['- Unlooping atrium', '- Unlooping ventricle']
    # names= ['unloopAtr', 'unloopVent']
    # chambers = ['Atrium', 'Ventricle',['Atrium', 'Ventricle']]
    # q_both_chs = ask4input('Select the chambers you would like to unloop: \n\t[0]: atrium\n\t[1]: ventricle\n\t[2]: both! >> : ', int)
    # selected_chambers = chambers[q_both_chs]
    
    if len(param_name) == 1:
        matrix_len = 8
        add_column = False
    else: 
        matrix_len = 9
        add_column = True
    
    num_ch_unlooped = -1
    for ksp_num, kspl in enumerate(kspls_HR):
    # for ksp_num in [0,1]:
        if ksp_num == 0 and 'Atrium' in selected_chambers or ksp_num == 1 and 'Ventricle' in selected_chambers:
            num_ch_unlooped+=1
            pass
        else: 
            continue
        
        # print(ksp_num)
        # Create matrix with all data
        #0:x, 1:y, 2:z, 3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: parameters
        matrix_unlooped = np.zeros((len(mesh.points()),matrix_len))
        matrix_unlooped[:,0:3] = mesh.points()
        matrix_unlooped[:,7] = param[0]
        if add_column: 
            matrix_unlooped[:,8] = param[1]
    
        # Get normals and centres of planes to cut heart (number of planes given as input+2)
        pl_normals, pl_centres = getPlaneNormals(no_planes = no_planes, spline_pts = kspl.points())
        # Give a number between 1-2 to each atrial plane, and between 0-1 to each ventricular plane
        plane_num = np.linspace(2-ksp_num,1-ksp_num,len(pl_normals))
        # Initialise index for HR centreline to identify chamber in which planes are cutting and flip plane values for ventricle
        # - Ventricle
        if ksp_num == 1:
            plane_num = plane_num[::-1]
        # - Atrium
        else: # ksp_num == 0
            index_guess = len(kspl_CLnew.points())-1
        
        bar = Bar(texts[ksp_num], max=len(pl_normals), suffix = suffix, check_tty=False, hide_cursor=False)
        # Iterate through planeS
        for i, normal, centre in zip(count(), pl_normals, pl_centres):
            # print('\n>-Plane num:', i, 'normal:', normal)
            # Initialise variables for each chamber
            if ksp_num == 1:
                normal = -normal
            # print('new_normal:', normal)
            # A. Get cut plane info
            # - Info Plane (Plane orientation/normal representation)
            arr_vectPlCut = Arrow(centre, centre+normal*10, s = 0.1, c='orange') #*
            #arr_vectPlCut_all.append(arr_vectPlCut)  #*
            
            # B. Cut high resolution centreline with plane and define chamber
            ksplCL_cut = kspl_CLnew.clone().cutWithMesh(Plane(pos=centre, normal=normal, sx=300), invert=True).lw(5).color('tomato')
            # Find point of centreline closer to last point of kspline that was cut by plane (ksplCL_cut)
            # If kspl_CLnew wasn't cut, then initialise it with last index of cl
            if len(ksplCL_cut.points()) == 0:
                pts_o = [kspl_CLnew.points()[-1]] #???
                # print('IF: length of cut centreline: ',len(pts_o))
            # Else, initialise it with all the points of the cut cl
            else: 
                pts_o = ksplCL_cut.points(); 
                # print('ELSE: length of cut centreline: ',len(pts_o))
            
            # If the ventricle is the one being unlooped, then initialise index as length of centreline
            if ksp_num == 1:
                index_guess= len(ksplCL_cut.points())
            # print('index_guess:',index_guess)
            
            #Find closest point between the pts_o and the high resolution centreline, initialising it using the index_guess prev defined
            pt_out, pt_num = findClosestPtGuess(pts_o, kspl_CLnew.points(), index_guess)
            index_guess = pt_num
            # print(index_guess)
            
            if pt_num < num_pt_AVC:
                chamber = 'ventricle'
            else: 
                chamber = 'atrium'
            
            # C. Cut surface centreline (kspl_vSurf) with plane and identify 0 deg angle point
            # Find point of surf_centreline cut by plane (ksplCL_vSurf_cut)
            kspl_vSurf_split_TF = []
            try: 
                kspl_vSurf_cut_T = kspl_vSurf.clone().cutWithMesh(Plane(pos=centre, normal=normal, sx=300), invert=True).lw(5).color('magenta')
                kspl_vSurf_split_T = kspl_vSurf_cut_T.splitByConnectivity()
                kspl_vSurf_split_TF.append(kspl_vSurf_split_T)
                # print('lenT',len(kspl_vSurf_split_T))
            except: 
                pass
            try: 
                kspl_vSurf_cut_F = kspl_vSurf.clone().cutWithMesh(Plane(pos=centre, normal=normal, sx=300), invert=False).lw(5).color('deepskyblue')
                kspl_vSurf_split_F = kspl_vSurf_cut_F.splitByConnectivity()
                kspl_vSurf_split_TF.append(kspl_vSurf_split_F)
                # print('lenF',len(kspl_vSurf_split_F))
            except: 
                pass
            # print('lenTF',len(kspl_vSurf_split_TF))

            # Initilise index_vSurf_cut when i == 0 for each chamber
            # print('>>>>>>>> i:', i)
            if i == 0:
                if ksp_num == 0:
                    print('\n\t-IN_ATRIUM')
                    index_vSurf_cut = len(kspl_vSurf.points())-1
                    sph_pt_uniq = Sphere(pos=kspl_vSurf.points()[index_vSurf_cut], r=3,c='hotpink')
                    plane = Plane(pos=centre, normal=normal, sx=300, alpha=0.5).color('medium orchid')
                    tol2use = 0#0.01
                elif ksp_num == 1: 
                    print('\n\t-IN_VENTRICLE')
                    index_vSurf_cut = 0
                    sph_pt_uniq = Sphere(pos=kspl_vSurf.points()[index_vSurf_cut], r=3,c='hotpink')
                    tol2use = tol
            else: 
                tol2use = tol
                    
            # Get all points of kspl_vSurf_cut in plane
            colors = ['orangered','gold', 'olive','lime', 'teal', 'aqua', 'dodgerblue', 'navy', 'indigo', 'purple', 'hotpink','chocolate']*10
            oo = 0
            kspl_split_plane = []
            pts_in_plane = []
            for jj, split_TF in enumerate(kspl_vSurf_split_TF):
                for kk, kspl_split in  enumerate(split_TF):
                    kspl_split_plane.append(kspl_split.color(colors[oo]).lw(2))
                    oo+=1
                    pts_split = [kspl_split.points()[0], kspl_split.points()[-1]]
                    for pt in pts_split: 
                        if isPointinInPlane(normal, centre, pt, tol = 0.001):
                            pts_in_plane.append(pt)

            pts_in_plane = np.unique(np.array(pts_in_plane), axis = 0)
            # print('\npts_in_plane:', pts_in_plane, '- shape:', pts_in_plane.shape)
            # indexes_min_dist = []
            if print_txt: 
                print('index_vSurf_cut_in:', index_vSurf_cut)
            
            # One or more points in plane
            if pts_in_plane.shape[0] > 0:
                # Initialised previous index
                if isinstance(index_vSurf_cut, int) or isinstance(index_vSurf_cut, np.int64):
                    sph_pt_uniq, index_vSurf_cut = find_sph_pt_uniq(pts_in_plane, kspl_vSurf, index_vSurf_cut = index_vSurf_cut)

                # Not initilised previous index
                else: 
                    if print_txt: 
                        print('aja! - Not initilised previous index!')
                    sphs_in_plane = []
                    for pp, pt in enumerate(pts_in_plane):
                        sph_pt = Sphere(pos = pt, c=colors[pp], r=3)
                        sph_pt.name = f"No.{pp}"
                        sphs_in_plane.append(sph_pt)
                    
                    ind_sph_sel = []
                    plane_new = Plane(pos=centre, normal=normal, sx=300, alpha=0.5).color('light steel blue')
                    def select_sph_in_plane(evt):
                        if not evt.actor: return
                        if isinstance(evt.actor, shapes.Sphere): 
                            sil = evt.actor.silhouette().lineWidth(6).c('red5')
                            print("\n\tYou clicked: Sphere "+evt.actor.name)
                            plt.remove(silcont.pop()).add(sil)
                            silcont.append(sil)
                            ind_sph_sel.append(evt.actor.name.split('.')[-1])
                
                    silcont = [None]
                    
                    txt_sel = Text2D('> Select best point cutting plane...', c=c, font=font)
                    plt = Plotter(axes=1, bg='gainsboro')
                    plt.addCallback('mouse click', select_sph_in_plane)
                    plt.show(mesh, kspl, plane_new, plane, kspl_vSurf, sphs_in_plane, txt_sel, zoom=1.2).close()
                    # print('ind_sph_sel:', ind_sph_sel)
                    if len(ind_sph_sel) > 0:
                        pts_in_plane_us = np.array([pts_in_plane[int(ind_sph_sel[-1])]])
                        sph_pt_uniq, index_vSurf_cut = find_sph_pt_uniq(pts_in_plane_us, kspl_vSurf, index_vSurf_cut = 0)
                    else: 
                        index_vSurf_cut = 'empty_notSelected'
            
            # No points in plane
            else: 
                if print_txt: 
                    print('aja! - No points in plane')
                index_vSurf_cut = 'empty2'
            
            plane = Plane(pos=centre, normal=normal, sx=300, alpha=0.5).color('medium orchid')
            
            if print_txt: 
                print('FINAL: index_vSurf_cut_out:', index_vSurf_cut, '- type:', type(index_vSurf_cut))
            
            if not isinstance(index_vSurf_cut, str):
                spheres_zeroDeg.append(sph_pt_uniq)
                # Vector from centre to cl_surface point being cut by plane
                v_zero = unit_vector(kspl_vSurf.points()[index_vSurf_cut] - centre)
                arr_vectZeroDeg.append(Arrow(centre, kspl_vSurf.points()[index_vSurf_cut], s = 0.1, c='dodgerblue'))
                # print('tol2use', tol2use)
                # D. Get points of mesh at plane
                d_points = np.absolute(np.dot(np.subtract(matrix_unlooped[:,0:3],np.asarray(centre)),np.asarray(normal)))
                # Find the indexes of the points that have not been yet taken, are at the plane and are in the 
                # chamber being analysed
                index_ptsAtPlane = np.where((d_points <= tol2use) & (matrix_unlooped[:,3] == 0) & (df_AtrVent == chamber))
                # print(d_points.min(), d_points.max(),'-lenptsatplane:',len(index_ptsAtPlane[0]))
                # Define new matrix just with the points on plane
                new_matrix = matrix_unlooped[index_ptsAtPlane,:][0]
                # - Get points of mesh that are on plane, centered on centreline point
                ptsC = np.subtract(new_matrix[:,0:3],np.asarray(centre))
                # - Get the radius of those points
                radius = [np.linalg.norm(x) for x in ptsC]
                
                # E. Find direction of point with respect to plane that includes central point, vC and the normal of the cutting plane
                # Vector normal to plane normal and v_zero (vector defined from centre of plane to cut-pt in centreline surface)
                normal_divLR = np.cross(normal, v_zero)
                # Define using these vectors if the points are all lying in the same side or not (1, -1)
                lORr = np.sign(np.dot(ptsC, np.asarray(normal_divLR)))
    
                # F. Get angle of points in that plane using v_zero
                av = np.dot(ptsC,v_zero)
                cosTheta = np.divide(av, radius) # vectors of magnitude one
                theta = np.arccos(cosTheta)*180/np.pi
                theta_corr = np.multiply(lORr, theta)
                
                # - Save all obtained values in matrix_unlooped
                for num, index in enumerate(index_ptsAtPlane[0]):
                    #3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: param
                    matrix_unlooped[index,3] = 1
                    matrix_unlooped[index,4] = plane_num[i]
                    matrix_unlooped[index,5] = theta_corr[num]
                    matrix_unlooped[index,6] = radius[num]
            
                # Plot stuff every X planes
                if plotshow:
                    if i % plotevery == 0:
                        sphL = []; sphR = []
                        for num, pt in enumerate(ptsC):
                            if num % 20 == 0:
                                if lORr[num] == 1:
                                    sphL.append(Sphere(pt+centre, r=2, c='blueviolet'))
                                else:
                                    sphR.append(Sphere(pt+centre, r=2, c='gold'))

                        text = filename+'\n\n >> Unlooping the heart (chamber: '+chamber+') - Plane No: '+str(i)+'/'+str(no_planes+2)
                        txt = Text2D(text, c=c, font=font)
                        plane = Plane(pos=centre, normal=normal, sx=300, alpha=0.5).color('medium orchid')
                        sph_centre = Sphere(centre, r=2, c='red')
                        arr_centre2vzero = Arrow(centre, kspl_vSurf.points()[index_vSurf_cut], s = 0.1, c='light green')
                        # kspl_vSurf_cut_all.append(kspl_vSurf_cut)
                        vp= Plotter(N=1, axes=13)
                        vp.show(mesh, sphL, sphR, kspl, plane, arr_vectPlCut, kspl_vSurf, kspl_split_plane, sph_pt_uniq, sph_centre, arr_centre2vzero, txt, at=0, interactive=True)
            
            else: 
                if print_txt:
                    print('-Plane No.', i, ' - PASS!')

            bar.next()
        
        # Save shapes 
        dict_shapes = addShapes2Dict (shapes = [kspl_vSurf, kspl_CLnew, kspl_atr, kspl_vent], dict_shapes = dict_shapes, radius = [[],[],[],[]], print_txt = False)
        
        df_unlooped = pd.DataFrame(matrix_unlooped, columns=['x','y','z','taken','z_plane','theta','radius']+param_name)
        df_unlooped = df_unlooped[df_unlooped['taken']==1]
        df_unlooped_f = df_unlooped.drop(['x', 'y','z'], axis=1)
        if add_column:
            df_unlooped_f.astype({'taken': 'bool','z_plane':'float16','theta':'float16','radius':'float16','cj_thickness':'float16','myoc_intBall':'float16' }).dtypes
        else: 
            try: 
                df_unlooped_f.astype({'taken': 'bool','z_plane':'float16','theta':'float16','radius':'float16',param_name[0]:'float16'}).dtypes
            except: 
                print('test did not work!')
        print('\n')
        dir2savef = os.path.join(dir_results, 'csv_all')
        saveDF(filename = filename, df2save = df_unlooped_f, df_name = 'df_'+save_names[num_ch_unlooped],#+dataname,
                        dir2save = dir2savef)
        dfs_unlooped.append(df_unlooped_f)
        
        bar.finish()
        return_list = [spheres_zeroDeg, arr_vectZeroDeg]
        
    return return_list,  kspl_vSurf, dfs_unlooped
            
#%% func - heatmapUnlooped
def heatmapUnlooped (filename, val2unloop, dir_results, dirImgs, save_names, hm_names, saveHM = True, 
                     savePlot = True, default = False, cmap = 'turbo'):
    """
    Function to create heatmap of unlooped data, specifically of the columns given as input by the val2unloop variable.  

    Parameters
    ----------
    filename :  str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    val2unloop: str
        String with the name of the column from the input dataframe that wants to be unlooped and represented as a 2D
        heatmap. 
    dataname : str
        String referencing the type of information that will be saved in the resulting dataframe. 
    dir_results : path
        Path to the folder where the results are saved.
    dirImgs : path
        Path to the folder where the images and videos are saved.
    names : list of str
        List of strings with the names of the dataframe files to use.
    saveHM : boolean, optional
        True if you want to save the resulting heatmap dataframe, else False. The default is True
    savePlot : boolean, optional
        True if you want to save the resulting heatmap plots, else False. The default is True
        
    Returns
    -------
    heatmaps : List of dataframes
        List of dataframes containing the heatmap representation of the val2unloop ('cj_thickness', 'myoc_intBall')
    scale_set : tuple
        Tuple of minimum and maximum value set for heatmap range. 

    """

    dir2savef = os.path.join(dir_results, 'csv_all')
    dirImgsf = new_dir(dirImgs, 'hm')
    heatmaps = []; scale_all = []
    for nn, val in zip(count(), val2unloop): 
        print('\n- Creating heatmaps for '+val+'...')
        
        for i, s_name in zip(count(), save_names): 
            print('\n> LOADING: ','df_'+s_name +'-i:'+str(i))
            print('- hm_name:', hm_names[nn][i])
            df = loadDF(filename, 'df_'+s_name, dir2savef)
            df = df.drop(['taken'], axis=1)
            df.astype('float16').dtypes
            
            heatmap = pd.pivot_table(df, values= val, columns = 'theta', index='z_plane', aggfunc=np.max)
            heatmap.astype('float16').dtypes
            
            if val == 'cj_thickness':
                title = filename +' - Cardiac jelly thickness [um] '
            elif val == 'myoc_intBall': 
                title = filename +' - Myocardium ballooning [um] '
            elif val == 'myoc_thickness':
                title = filename +' - Myocardial thickness [um] '
            elif val == 'endo_thickness':
                title = filename +' - Endocardial thickness [um] '
                
            if 'Atr' in s_name:
                # dir4heatmap = os.path.join(dirImgs,filename+'_hm_'+name[0:9]+'_'+dataname+'.png')
                title = title +'- Atrium'
            else: 
                # dir4heatmap = os.path.join(dirImgs,filename+'_hm_'+name[0:10]+'_'+dataname+'.png')
                title = title +'- Ventricle'
            print('\t- title:', title)
            
            if i == 0: 
                if val == 'cj_thickness':
                    if not default: 
                        scale_bar = ask4input('Do you want to set a range to the scale bar of -'+ val+'-? \n\t[0]: no, use the default values 0-25 um \n\t[1]: yes, I would like to set the range! >>>: ', bool)
                        if scale_bar: 
                            vmin = ask4input('Minimum value for scale bar of -'+val+'-[um]: ',float)
                            vmax = ask4input('Maximum value for scale bar of -'+val+'-[um]: ',float)
                        else: 
                            vmin = 0
                            vmax = 25
                    else: 
                        vmin = 0
                        vmax = 25
                elif val in ['myoc_thickness', 'endo_thickness']:
                    if not default: 
                        scale_bar = ask4input('Do you want to set a range to the scale bar of -'+ val+'-? \n\t[0]: no, use the default values 0-15 um \n\t[1]: yes, I would like to set the range! >>>: ', bool)
                        if scale_bar: 
                            vmin = ask4input('Minimum value for scale bar of -'+val+'-[um]: ',float)
                            vmax = ask4input('Maximum value for scale bar of -'+val+'-[um]: ',float)
                        else: 
                            vmin = 0
                            vmax = 15
                    else: 
                        vmin = 0
                        vmax = 15
                elif val == 'myoc_intBall': 
                    if not default: 
                        scale_bar = ask4input('Do you want to set a range to the scale bar of -'+ val+'-? \n\t[0]: no, use the default values 0-100 um \n\t[1]: yes, I would like to set the range! >>>: ', bool)
                        if scale_bar: 
                            vmin = ask4input('Minimum value for scale bar of -'+val+'-[um]: ',float)
                            vmax = ask4input('Maximum value for scale bar of -'+val+'-[um]: ',float)
                        else: 
                            vmin = 0
                            vmax = 100
                    else: 
                        vmin = 0
                        vmax = 100
                        
            # Make figure
            fig, ax = plt.subplots(figsize=(16, 10))
            ax = sns.heatmap(heatmap, cmap=cmap, vmin = vmin, vmax = vmax)#, xticklabels=20, yticklabels=550)
            
            max_val = df.z_plane.max()
            if max_val > 1.1: 
                y_labels = [1,2]
                y_text = '[Inflow tract >> Valve]'
            else: # Ventricle
                y_labels = [0,1]
                y_text = '[Valve >> Outflow tract]'
                
            x_pos = ax.get_xticks()
            # x_lab = ax.get_xticklabels()
            x_pos_new = np.linspace(x_pos[0], x_pos[-1], 19)
            x_lab_new = np.arange(-180,200,20)
            ax.set_xticks(x_pos_new) 
            ax.set_xticklabels(x_lab_new, rotation=30)
            
            y_pos = ax.get_yticks()
            # y_lab = ax.get_yticklabels()
            y_pos_new = np.linspace(y_pos[0], y_pos[-1], 11)
            y_lab_new = np.linspace(y_labels[0],y_labels[1],11)
            y_lab_new = [format(y,'.2f') for y in y_lab_new]
            ax.set_yticks(y_pos_new) 
            ax.set_yticklabels(y_lab_new, rotation=0)
            
            plt.ylabel('Centreline position '+y_text+'\n', fontsize=10)
            plt.xlabel('Angle (\N{DEGREE SIGN}) [Dorsal >> Right >> Ventral >> Left >> Dorsal]', fontsize=10)
        
            plt.title(title, fontsize = 15)
            
            # if isinstance(hm_name, list):
            #     sp_hm_name = hm_name[i]
            # elif isinstance(hm_name, str):
            #     sp_hm_name = hm_name
            sp_hm_name = hm_names[nn][i]
            print('\t- sp_hm_name:',sp_hm_name)
            
            if savePlot: 
                dir4heatmap = os.path.join(dirImgsf,filename+'_hm_'+sp_hm_name+'.png')
                plt.savefig(dir4heatmap, dpi=300, bbox_inches='tight', transparent=True)
                # print(dir4heatmap)
                
            plt.show()
            if saveHM: 
                saveDF(filename = filename, df2save = heatmap, df_name = 'hm_'+sp_hm_name, dir2save = dir2savef)
                print('\t- Saved hm_'+sp_hm_name)
            
            heatmaps.append(heatmap)
            alert('frog', 1)
        alert('wohoo', 1)
        scale_set = (vmin,vmax)
        scale_all.append(scale_set)
    
    return heatmaps, scale_all
    
#%% func - filterUnloopedDF
def filterUnloopedDF(filename, val2unloop, dir_results, dir_data2Analyse, save_names, hm_names, saveHM = True, cmap= 'turbo'):
    """
    #https://stackoverflow.com/questions/38940946/average-of-multiple-dataframes-with-the-same-columns-and-indices
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    thickness : TYPE
        DESCRIPTION.
    dir_results : TYPE
        DESCRIPTION.
    dir_data2Analyse : TYPE
        DESCRIPTION.
    names : TYPE
        DESCRIPTION.
    save_names : TYPE
        DESCRIPTION.
    saveHM : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    heatmaps : TYPE
        DESCRIPTION.

    """
    
    print('\n>>> Filtering dataframes...')
    dir2savef = os.path.join(dir_results, 'csv_all')
    heatmaps = []
    for nn, val in zip(count(), val2unloop): 
        print('\t- Filtering ', val)
        for ii, name in zip(count(), save_names): 
            print('- LOADING: ','df_'+name )
            print('- hm_name:', hm_names[nn][ii])
            df = loadDF(filename, 'df_'+name, dir2savef)
            df = df.drop(['taken'], axis=1)
            df.astype('float16').dtypes
            df_cols = ['z_plane', 'theta', 'radius']+[val],# list(df.columns)
            
            angles_f = np.linspace(-180,180,num = 360*6+1, endpoint = True)
            step = round((angles_f[1]-angles_f[0])/2, 4)
            
            th_list = []
            # myocIntBall_list = []
            rad_list = []
            zplane_list = []
            theta_list = []
            
            for j, z_plane in enumerate(sorted(df.z_plane.unique())):
                df_z = df[df['z_plane'] == z_plane]
                #max? mean?
                for i, ang in enumerate(angles_f):
                    theta_list.append(round(ang,3))
                    filt = df_z[(df_z['theta'] >= ang-step) & (df_z['theta'] < ang+step)]
                    th_val = filt[val].max()
                    th_list.append(th_val)
                    # if thickness == 'cj_thickness': 
                    #     myocIntBall_val = filt['myoc_intBall'].max()
                    #     myocIntBall_list.append(myocIntBall_val)
                    
                    rad_val = filt['radius'].max()
                    rad_list.append(rad_val)
                    
                    zplane_list.append(round(z_plane, 3))
            # if val == 'cj_thickness': 
            #     df_filt = pd.DataFrame(list(zip(zplane_list, theta_list, rad_list, th_list, myocIntBall_list)), 
            #                            columns =df_cols) 
            # else: 
            df_filt = pd.DataFrame(list(zip(zplane_list, theta_list, rad_list, th_list)), columns =df_cols[0]) 
            
            alert('wohoo', 1)
            df_filt.astype('float16').dtypes
            heatmap = pd.pivot_table(df_filt, values= val, columns = 'theta', index='z_plane', aggfunc=np.max)
            heatmap.astype('float16').dtypes
            fig, ax = plt.subplots(figsize=(16, 10))
            # ax = sns.heatmap(heatmap, cmap='jet')
            sns.heatmap(heatmap, cmap=cmap)#, vmin = vmin, vmax = vmax)#, xticklabels=20, yticklabels=550)
            plt.show()
            
            # if isinstance(hm_name, list):
            #     sp_hm_name = hm_name[n]
            # elif isinstance(hm_name, str):
            #     sp_hm_name = hm_name
            sp_hm_name = hm_names[nn][ii]
            print('\t- sp_hm_name:',sp_hm_name)
            
            if saveHM: 
                saveDF(filename = filename, df2save = heatmap, df_name = 'hmf_'+sp_hm_name, dir2save = dir2savef)
                saveDF(filename = filename, df2save = heatmap, df_name = 'hmf_'+sp_hm_name, dir2save = os.path.join(dir_data2Analyse, 'R_All','df_all','df_hmf'))
                print('\t- Saved hmf_'+sp_hm_name)
                
            alert('frog', 1)
            heatmaps.append(heatmap)
        
    return heatmaps

#%% func - normUnloopedDF
def normUnloopedDF(filename, thickness, dir_results, dir_data2Analyse, names, save_names, saveHM = True):
    """
    #https://stackoverflow.com/questions/38940946/average-of-multiple-dataframes-with-the-same-columns-and-indices
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    thickness : TYPE
        DESCRIPTION.
    dir_results : TYPE
        DESCRIPTION.
    dir_data2Analyse : TYPE
        DESCRIPTION.
    names : TYPE
        DESCRIPTION.
    save_names : TYPE
        DESCRIPTION.
    saveHM : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    heatmaps : TYPE
        DESCRIPTION.

    """
    df_res = loadDF(filename, 'ResultsDF', dir_results)
    file_num = df_res[df_res['Folder']==filename+'_2A'].index.values[0]
    
    dict_norm = {'cj_thickness': {'unloopAtrCjTh_myocIntBall': df_res.loc[file_num,'Vol_Atr.CJ'], 'unloopVentCjTh_myocIntBall': df_res.loc[file_num,'Vol_Vent.CJ']},
                 'myoc_intBall': {'unloopAtrCjTh_myocIntBall': df_res.loc[file_num,'Vol_Atr.ExtMyoc'], 'unloopVentCjTh_myocIntBall': df_res.loc[file_num,'Vol_Vent.ExtMyoc']}}
    print(dict_norm)
    
    heatmaps = []
    for n, name in enumerate(names): 
        df = loadDF(filename, 'df_'+name, dir_results)
        df = df.drop(['taken'], axis=1)
        df.astype('float16').dtypes
        df_cols = list(df.columns)
        df_cols.append('cj_thickness_norm')
        df_cols.append('myoc_intBall_norm')
        
        angles_f = np.linspace(-180,180,num = 360*6+1, endpoint = True)
        step = round((angles_f[1]-angles_f[0])/2, 4)
        
        cjTh_list = []
        myocIntBall_list = []
        cjTh_norm = []
        myocIntBall_norm = []
        rad_list = []
        zplane_list = []
        theta_list = []
        
        for j, z_plane in enumerate(sorted(df.z_plane.unique())):
            df_z = df[df['z_plane'] == z_plane]
            #max? mean?
            for i, ang in enumerate(angles_f):
                theta_list.append(round(ang,3))
                filt = df_z[(df_z['theta'] >= ang-step) & (df_z['theta'] < ang+step)]
                cjTh_val = filt['cj_thickness'].max()
                cjTh_list.append(cjTh_val)
                cjTh_norm.append(cjTh_val*1000000/dict_norm['cj_thickness'][name])
                
                myocIntBall_val = filt['myoc_intBall'].max()
                myocIntBall_list.append(myocIntBall_val)
                myocIntBall_norm.append(myocIntBall_val*1000000/dict_norm['myoc_intBall'][name])
                
                rad_val = filt['radius'].max()
                rad_list.append(rad_val)
                
                zplane_list.append(round(z_plane, 3))
        df_filt = pd.DataFrame(list(zip(zplane_list, theta_list, rad_list, cjTh_list, myocIntBall_list, cjTh_norm, myocIntBall_norm)), 
                               columns =df_cols) 
        
        alert('wohoo', 1)
        df_filt.astype('float16').dtypes
        heatmap = pd.pivot_table(df_filt, values= thickness+'_norm', columns = 'theta', index='z_plane', aggfunc=np.max)
        heatmap.astype('float16').dtypes
        fig, ax = plt.subplots(figsize=(16, 10))
        # ax = sns.heatmap(heatmap, cmap='jet')
        sns.heatmap(heatmap, cmap='jet')#, vmin = vmin, vmax = vmax)#, xticklabels=20, yticklabels=550)
        plt.show()

        if saveHM: 
            saveDF(filename = filename, df2save = heatmap, df_name = 'hmN_'+save_names[n], dir2save = dir_results)
            saveDF(filename = filename, df2save = heatmap, df_name = 'hmN_'+save_names[n], dir2save = os.path.join(dir_data2Analyse, 'R_All','df_hmN'))
            
            alert('frog', 1)
        heatmaps.append(heatmap)
        
    return heatmaps
        

#%% - PROCESS DATA
#%% func - rotation_matrix
def rotation_matrix(axis, theta):
    """
    Returns the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    https://stackoverflow.com/questions/6802577/rotation-of-3d-vector

    Parameters
    ----------
    axis : numpy array
        np array with the x,y,z coordinates of the original's plane axis
    theta : float
        Angle of rotation around the axis.

    Returns
    -------
    numpy array
        Rotation matrix.

    """

    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

#%% func - unit_vector
def unit_vector(v):
    """
    Function that returns the unit vector of the vector given as input

    Parameters
    ----------
    v : numpy array
        np array with the x,y,z coordinates of a plane axis

    Returns
    -------
    numpy array
        np array with the x,y,z coordinates of a unitary plane axis

    """

    #mag = math.sqrt(sum(i**2 for i in x))
    sqrs = [i**2 for i in v]
    mag = math.sqrt(sum(sqrs))
    v_unit = [j/mag for j in v]

    return np.asarray(v_unit)

#%% func - newNormal3DRot
def newNormal3DRot (normal, rotX, rotY, rotZ):
    """
    Function that returns a vector rotated around X, Y and Z axis

    Parameters
    ----------
    normal : numpy array
        np array with the x,y,z coordinates of the original's plane axis
    rotX : list of floats
        List of angles (deg) of the resulting rotation around the x-axis.
    rotY : list of floats
        List of angles (deg) of the resulting rotation around the Y-axis.
    rotZ : list of floats
        List of angles (deg) of the resulting rotation around the Z-axis.

    Returns
    -------
    normal_rotZ : numpy array
        np array with the x,y,z coordinates of the rotated axis

    """

    ang_X = np.radians(sum(rotX))
    ang_Y = np.radians(sum(rotY))
    ang_Z = np.radians(sum(rotZ))
    # print(sum(rotX), sum(rotY), sum(rotZ))

    normal_rotX = (np.dot(rotation_matrix(axis = [1,0,0], theta = ang_X), normal))
    normal_rotY = (np.dot(rotation_matrix(axis = [0,1,0], theta = ang_Y), normal_rotX))
    normal_rotZ = (np.dot(rotation_matrix(axis = [0,0,1], theta = ang_Z), normal_rotY))

    return normal_rotZ

#%% func - findClosestPt
def findClosestPt(pt_o, pts):
    """
    Function that finds the closest point of the centreline where a plane cuts it

    Parameters
    ----------
    pt_o : numpy array
        np array with the x,y,z coordinates of the point that cuts kspline
    pts : array of coordinates
        Array with x,y,z coordinates of the centreline points

    Returns
    -------
    pt_out : numpy array
        np array with the x,y,z coordinates of the centreline point closest to the kspline cut
    num_pt : int
        Index of the centreline point closer to the plane that cuts the kspline.

    """

    minima_distancia = 9999999
    # for n in range(len(pts)):
    for index, value in enumerate(pts):
        # pt = pts[n]
        squared_dist = np.sum((value-pt_o)**2, axis=0)
        dist = np.sqrt(squared_dist)

        if dist < minima_distancia:
            minima_distancia = dist
            pt_out = value
            num_pt = index
        print(f"num:{minima_distancia}")
        if dist==0:
            print(pt_o)
            print(value)
            print(squared_dist)
            print(dist)
            input()
            
    print(f"Distancia:{dist}")
   
    if not num_pt:
        raise Exception("no num")
    return pt_out, num_pt

#%% func - findClosestPtGuess
def findClosestPtGuess(pts_cut, pts, index_guess):
    """
    Function that finds the closest point of the centreline where a plane cuts given a point index as an initial guess

    Parameters
    ----------
    pts_cut : array of coordinates
        np array with the x,y,z coordinates of all the points that make the cut section of the kspline
    pts : array of coordinates
        Array with x,y,z coordinates of ALL the centreline points
    index_guess : int
        Index of the point guessed to be closest to new closest point

    Returns
    -------
    pt_out : numpy array
        np array with the x,y,z coordinates of the centreline point closest to the kspline cut
    num_pt : int
        Index of the centreline point closer to the plane that cuts the kspline.

    """
    # First find the closes point of the input pts_cut to pt_guess = pts[index]
    pt_guess = pts[index_guess]
    min_dist_guess = 999999999999999999
    for pt_1 in pts_cut:
        squared_dist_guess = np.sum((pt_1-pt_guess)**2, axis=0)
        dist_guess = np.sqrt(squared_dist_guess)

        if dist_guess < min_dist_guess:
            min_dist_guess = dist_guess
            pt_o = pt_1
            
    # print('>> pt_o:', pt_o,'- pt_guess:', pt_guess, '- min_dist_guess:', min_dist_guess)
    
    # Now that we know the coordinates (pt_o) of the cut centreline that is cut by the plane and closest to the pt_guess then
    # find the closest point of all the centreline points closest to pt_o and get also its index
    min_dist = 999999999999999999
    for n, pt in enumerate(pts):#range(len(pts)):
        # pt = pts[n]
        squared_dist = np.sum((pt-pt_o)**2, axis=0)
        dist = np.sqrt(squared_dist)

        if dist < min_dist:
            min_dist = dist
            pt_out = pt
            num_pt = n
    
    # print('>> pt_out:', pt_out,'- num_pt:', num_pt, '- min_dist:', min_dist)
    # input()

    return pt_out, num_pt

#%% func - isPointinInPlane
def isPointinInPlane(normal, centre, point, tol=0.001):
    # print('tol:',tol)
    diff_points = point - centre
    # print('dot:',np.dot(diff_points,normal))
    if abs(np.dot(diff_points,normal)) < tol:
        return True
    else: 
        return False
    
#%% func - findDist
def findDist(pt1, pt2):
    """
    Function that returns the distance between two points given as input

    Parameters
    ----------
    pt1 : numpy array
        np array with the x,y,z coordinates of a point
    pt2 : numpy array
        np array with the x,y,z coordinates of a point

    Returns
    -------
    dist : int
        Resulting distance between pt1 and pt2

    """
    squared_dist = np.sum((pt1-pt2)**2, axis=0)
    dist = np.sqrt(squared_dist)

    return dist

#%% func - find_sph_pt_uniq
def find_sph_pt_uniq(pts_in_plane, kspl_vSurf, index_vSurf_cut = 0):
    
    indexes_min_dist = []
    dists = cdist(pts_in_plane,kspl_vSurf.points())
    # print('shape:', dists.shape)
    min_dists = np.amin(dists, axis=1)
    for x in range(pts_in_plane.shape[0]):
        where = np.where(dists[x,:] == min_dists[x])
        indexes_min_dist.append(where[0][0])
    # print('indexes_min_dist:', indexes_min_dist)
    diff_indexes = abs(np.array(indexes_min_dist) - index_vSurf_cut)
    # print('diff_indexes:', diff_indexes)
    ind_min_diff_indexes = np.where(diff_indexes == np.min(diff_indexes))[0][0]
    index_vSurf_cut = indexes_min_dist[ind_min_diff_indexes]
    # print('index_vSurf_cut_out:', index_vSurf_cut, '- type:', type(index_vSurf_cut))
    pt_uniq = kspl_vSurf.points()[index_vSurf_cut]
    sph_pt_uniq = Sphere(pos=pt_uniq, r=3,c='hotpink')
    
    return sph_pt_uniq, index_vSurf_cut
    
#%% func - getPointsAtPlane
def getPointsAtPlane (points, pl_normal, pl_centre, tol=2, addData = []):
    """
    Function to get points within mesh at certain heights (y positions) to create kspline

    Parameters
    ----------
    points : array of coordinates
        Array with x,y,z coordinates of the mesh4cl points
    pl_normal : list of floats
        List with the x,y,z coordinatesof the plane's normal
    pl_centre : list of floats
        List with the x,y,z coordinatesof the plane's centre
    tol : float, optional
        Tolerance defined to get points in plane. The default is 2.
    addData : list, optional
        List of additional parameter associated to input points that wants to be extracted from points at plane. 
        Default is [] (empty list).

    Returns
    -------
    pts_cut : array of coordinates
        Array with x,y,z coordinates of the mesh4cl points that are cut by plane
    data_cut : array of cut data - same length as pts_cut
        List of additional data of only points found at plane. If addData = [], then data_cut = []
    
    """

    pts_cut = []
    data_cut = []
    
    d = pl_normal.dot(pl_centre)
    #print('d for tol:', d)
    d_range = [d-tol, d+tol]
    #print('d range:', d_range)
    d_range.sort()

    for i, pt in enumerate(points):
        d_pt = pl_normal.dot(pt)
        if d_pt>d_range[0] and d_pt<d_range[1]:
            pts_cut.append(pt)
            if addData != []:
                data_cut.append(addData[i])
            
            #print(pt)

    pts_cut = np.asarray(pts_cut)
    data_cut = np.asarray(data_cut)

    return pts_cut, data_cut

#%% func - order_pts
def order_pts (points):
    """
    Function that returns an ordered array of points

    Parameters
    ----------
    points :  array of coordinates
        Array with x,y,z coordinates to order

    Returns
    -------
    ordered_pts : array of coordinates
        Array with the ordered x,y,z coordinates
    angle_deg : list of floats
        List with the associated angle of each of the ordered points

    """

    center_pt = np.mean(points, axis=0)
    cent_pts = np.zeros_like(points)
    for num in range(len(points)):
        cent_pts[num][0]=points[num][0]-center_pt[0]
        cent_pts[num][2]=points[num][2]-center_pt[2]

    angle_deg = np.zeros((len(points),1))
    for num, pt in enumerate(cent_pts):
        angle_deg[num] = np.arctan2(pt[2],pt[0])*(180/np.pi)

    index_sort = np.argsort(angle_deg, axis=0)
    ordered_pts = np.zeros_like(points)
    for i, pos in enumerate(index_sort):
        ordered_pts[i]=points[pos]

    return ordered_pts, angle_deg

#%% func - getInterpolatedPts
def getInterpolatedPts(points, nPoints):
    """
    Function that interpolates input points

    Parameters
    ----------
    points : array of coordinates
        Array with original x,y,z coordinates of points
    nPoints : int
        Number of final points in spline

    Returns
    -------
    pts_interp : array of coordinates
        Final array with interpolated x,y,z coordinates of points

    """

    minx, miny, minz = np.min(points, axis=0)
    maxx, maxy, maxz = np.max(points, axis=0)
    maxb = max(maxx - minx, maxy - miny, maxz - minz)
    smooth = 0.5*maxb/2  # must be in absolute units

    x = points[:,0]
    y = points[:,1]
    z = points[:,2]

    #https://stackoverflow.com/questions/47948453/scipy-interpolate-splprep-error-invalid-inputs
    okay = np.where(np.abs(np.diff(x)) + np.abs(np.diff(y)) + np.abs(np.diff(z)) > 0)
    xp = np.r_[x[okay], x[-1]]
    yp = np.r_[y[okay], y[-1]]
    zp = np.r_[z[okay], z[-1]]

    tck, u = splprep([xp, yp, zp], s=smooth)
    new_array = np.linspace(0, 1, nPoints)
    xnew, ynew, znew = splev(new_array, tck)

    pts_interp = np.c_[xnew, ynew, znew]

    return pts_interp

#%% func - findAngleBtwVectorsZ
def findAngleBtwVectorsZ(pts1, pts2):
    """
    Function that returns the angle between two vectors on the XY-plane 

    Parameters
    ----------
    pts1 : array of coordinates defining a vector
        Coordinates defining the head and tail of vector
    pts2 : array of coordinates defining a vector
        Coordinates defining the head and tail of vector

    Returns
    -------
    angle : float
        Angle between vectors

    """

    mag_v1 = findDist(pts1[0],pts1[1])
    mag_v2 = findDist(pts2[0],pts2[1])
    
    vect1 = pts1[1]-pts1[0]
    vect2 = pts2[1]-pts2[0]

    dotProd = np.dot(vect1,vect2)

    angle = math.degrees(math.acos(dotProd/(mag_v1*mag_v2)))

    return angle

#%% func - orientVectors
def orientVectors(line, ref_pt = []):
    """
    Function that orients the input line in a particular direction

    Parameters
    ----------
    line : line
        Line defining orientation (vedo line)
    ref_pt : array of coordinates defining a point to use as reference, optional
        Coordinates defining the reference point to use to orient the vector. If line is atr_or, the reference point 
        should be the point of linLine in the outflow tract. if line is vent_or, the reference point 
        should be the point of linLine in the inflow tract. The default is [].
    
    Returns
    -------
    pts_linef : array of coordinates defining a vector
        Coordinates defining the head and tail of the reoriented vector

    """
    pts_line = line.points()
    if ref_pt !=[]:#line._legend  == 'lin_OrientAtr(ProjX)' or line._legend == 'lin_OrientVent(ProjX)':
        # print('OpA-Atrium')
        pts_linef = np.ones_like(pts_line)
        dist_0 = findDist(ref_pt, pts_line[0]) 
        dist_1 = findDist(ref_pt, pts_line[1])
        # print(dist_0, dist_1)
        if dist_0 < dist_1:
            pts_linef[0] = pts_line[0]
            pts_linef[1] = pts_line[1]
        else: 
            pts_linef[0] = pts_line[1]
            pts_linef[1] = pts_line[0]

    else: # Linear line "linLine(ProjX)" # head first, tail second
        # print('OpC-Heart')
        pts_linef = pts_line[pts_line[:,1].argsort()[::-1]]
        # print('aja')

    #print(np.diff(pts_line, axis = 0))

    return  pts_linef

#%% func - classifyPtsMx
def classifyPtsMx(dict_planes, pl_name, pts_whole, spaw_analysis):
    """
    Function that classifies the input points (pts_whole) as atrium/ventricle or dorsal/ventral depending on the rest of the inputs given

    Parameters
    ----------
    dict_planes : dictionary
        Initialised dictionary with planes information
    pl_name : str
        Name of the plane to use - this defines the type of classification that will be made
    pts_whole : array of coordinates
        Array with x,y,z coordinates of whole mesh

    Returns
    -------
    pts_classFinal : list of str
        List with the final classification for each of the points given as input

    """

    if pl_name == 'pl2CutMesh_Chamber':
        ptA = 'atrium'
        ptB = 'ventricle'
        name = 'Atrium-Ventricle'
        
    elif pl_name == 'pl_AtrCoronal':
        if not spaw_analysis: 
            ptA = 'dorsal'
            ptB = 'ventral'
            name = 'Dorsal-Ventral (Atr)'
        else: 
            ptA = 'right'
            ptB = 'left'
            name = 'Left-Right (Atr)'
            
    elif pl_name == 'pl_VentCoronal':
        if not spaw_analysis: 
            ptA = 'ventral'
            ptB = 'dorsal'
            name = 'Dorsal-Ventral (Vent)'
        else: 
            ptA = 'left'
            ptB = 'right'
            name = 'Left-Right (Vent)'
            
    print('--> '+name)

    pts_classFinal = np.empty(len(pts_whole), dtype='object')
    # - AnV
    pl_normal = dict_planes[pl_name]['pl_normal']
    # Make normal a unit vector
    normal_unit = unit_vector(pl_normal)
    pl_centre = dict_planes[pl_name]['pl_centre']

    # Find all the d values of pix_um
    d_pts = np.dot(np.subtract(pts_whole,np.array(pl_centre)),np.array(normal_unit))

    # Find all positions in d_pts that are at each side of the plane
    pos_A = np.where(d_pts < 0)[0]
    pos_B = np.where(d_pts > 0)[0]

    # Remove the points that belong to the other side
    pts_classFinal[pos_A] = ptA
    pts_classFinal[pos_B] = ptB

    return pts_classFinal

#%% func - classifyHeartPts
def classifyHeartPts(filename, df_res, file_num, dict_planes, m_whole, m_left, m_atr, data, names_data, plot_show = True):
    """
    Function that classifies the points that make up a mesh as atrium/ventricle, dorsal/ventral and left/right

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    dict_planes : dictionary
        Initialised dictionary with planes information
    m_whole : mesh
        Mesh whose points are going to be classified
    m_left : mesh
        Left side of mesh to be classified
    m_atr : mesh
        Atrial mesh to be classified
    data : list of arrays
        List of arrays with distance data to be classified
    names_data : list of str
        List with the names of the data corresponding to each point
    plot_show : boolean, optional
        True if you want to see the resulting mesh in a plot, else False. The default is True.

    Returns
    -------
    df_classPts : dataframe
        Dataframe including the points classification and corresponding distance data

    """

    print('\n- Classifying points for ', names_data)
    tic = perf_counter()

    cols = ['AtrVent','DorsVent','LeftRight'] + names_data
    df_classPts = pd.DataFrame(columns=cols)

    spaw_analysis = False
    if 'spaw' in df_res.loc[file_num,'spAnalysis']:
        spaw_analysis = True #ask4input('You are processing a heart that came from an incross of spaw heterozygous.\n  Please, select the way this heart is looping to continue processing: \n\t[0]: right-left\n\t[1]: dorso-ventral >>>: ', bool)
        
    # Classify AnV
    print('--> Atrium-Ventricle')
    pts_atr = m_atr.points()
    pts_whole = m_whole.points()
    av = pts_whole.view([('', pts_whole.dtype)] * pts_whole.shape[1]).ravel()
    cv = pts_atr.view([('', pts_atr.dtype)] * pts_atr.shape[1]).ravel()
    d_isin = np.isin(av,cv)
    index_atr = np.where(d_isin == True)[0]
    pts_classAnV = np.empty(len(pts_whole), dtype='object')
    pts_classAnV[:] = 'ventricle'
    pts_classAnV[index_atr] = 'atrium'
    del av, cv, pts_whole
    
    # pts_classAnV = classifyPtsMx(dict_planes = dict_planes, pl_name = 'pl2CutMesh_Chamber', pts_whole = pts_whole)
    df_classPts['AtrVent'] = pts_classAnV
    
    if not spaw_analysis:
        # Classify DnV
        pts_whole = m_whole.points()
        # - Classify DnV_Atr
        pts_classDnV_Atr = classifyPtsMx(dict_planes = dict_planes, pl_name = 'pl_AtrCoronal', pts_whole =  pts_whole, spaw_analysis = spaw_analysis)
        # - Classify DnV_Vent
        pts_classDnV_Vent = classifyPtsMx(dict_planes = dict_planes, pl_name = 'pl_VentCoronal', pts_whole =  pts_whole, spaw_analysis = spaw_analysis)
    
        atr = np.where(pts_classAnV == 'atrium')[0]
        vent = np.where(pts_classAnV == 'ventricle')[0]
    
        pts_classDnV_Atr[vent] = ''
        pts_classDnV_Vent[atr] = ''
    
        pts_classDnV = pts_classDnV_Atr+pts_classDnV_Vent
        df_classPts['DorsVent'] = pts_classDnV
        del pts_whole
    
        # Classify LnR
        print('--> Left-Right')
        pts_left = m_left.points()
        pts_whole = m_whole.points()
        av = pts_whole.view([('', pts_whole.dtype)] * pts_whole.shape[1]).ravel()
        bv = pts_left.view([('', pts_left.dtype)] * pts_left.shape[1]).ravel()
        # cint = np.intersect1d(av, bv).view(a.dtype).reshape(-1, a.shape[1])
        c_isin = np.isin(av,bv)
        index_left = np.where(c_isin == True)[0]
    
        pts_classLnR = np.empty(len(pts_whole), dtype='object')
        pts_classLnR[:] = 'right'
        pts_classLnR[index_left] = 'left'
    
        df_classPts['LeftRight'] = pts_classLnR
    
        for i, name, dat in zip(count(), names_data, data):
            df_classPts[name] = dat
    
    # spaw analysis
    else: 
        # Classify LnR
        pts_whole = m_whole.points()
        # - Classify LnR_Atr
        pts_classLnR_Atr = classifyPtsMx(dict_planes = dict_planes, pl_name = 'pl_AtrCoronal', pts_whole =  pts_whole, spaw_analysis = spaw_analysis)
        # - Classify DnV_Vent
        pts_classLnR_Vent = classifyPtsMx(dict_planes = dict_planes, pl_name = 'pl_VentCoronal', pts_whole =  pts_whole, spaw_analysis = spaw_analysis)
    
        atr = np.where(pts_classAnV == 'atrium')[0]
        vent = np.where(pts_classAnV == 'ventricle')[0]
    
        pts_classLnR_Atr[vent] = ''
        pts_classLnR_Vent[atr] = ''
    
        pts_classLnR = pts_classLnR_Atr+pts_classLnR_Vent
        df_classPts['LeftRight'] = pts_classLnR
        del pts_whole
    
        # Classify DnV
        print('--> Dorsal-Ventral')
        pts_left = m_left.points()
        pts_whole = m_whole.points()
        av = pts_whole.view([('', pts_whole.dtype)] * pts_whole.shape[1]).ravel()
        bv = pts_left.view([('', pts_left.dtype)] * pts_left.shape[1]).ravel()
        # cint = np.intersect1d(av, bv).view(a.dtype).reshape(-1, a.shape[1])
        c_isin = np.isin(av,bv)
        index_left = np.where(c_isin == True)[0]
    
        pts_classDnV = np.empty(len(pts_whole), dtype='object')
        pts_classDnV[:] = 'ventral'
        pts_classDnV[index_left] = 'dorsal'
    
        df_classPts['DorsVent'] = pts_classDnV
    
        for i, name, dat in zip(count(), names_data, data):
            df_classPts[name] = dat

    toc = perf_counter()
    time = toc-tic
    print("- All Done - points have been classified!\n- Sample of classified points")
    print(df_classPts.sample(10))
    print("- Time taken to classify = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")
    alert('whistle',1)

    if plot_show:
        plotPtClassif(filename, m_whole, pts_whole, [pts_classAnV, pts_classDnV, pts_classLnR])

    return df_classPts

#%% func - classifyHeartPts_BU
def classifyHeartPts_BU(filename, dict_planes, m_whole, m_left, m_atr, data, names_data, plot_show = True):
    """
    Function that classifies the points that make up a mesh as atrium/ventricle, dorsal/ventral and left/right

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    dict_planes : dictionary
        Initialised dictionary with planes information
    m_whole : mesh
        Mesh whose points are going to be classified
    m_left : mesh
        Left side of mesh to be classified
    m_atr : mesh
        Atrial mesh to be classified
    data : list of arrays
        List of arrays with distance data to be classified
    names_data : list of str
        List with the names of the data corresponding to each point
    plot_show : boolean, optional
        True if you want to see the resulting mesh in a plot, else False. The default is True.

    Returns
    -------
    df_classPts : dataframe
        Dataframe including the points classification and corresponding distance data

    """

    print('- Classifying points for ', names_data)
    tic = perf_counter()

    cols = ['AtrVent','DorsVent','LeftRight'] + names_data
    df_classPts = pd.DataFrame(columns=cols)

    # Classify AnV
    print('--> Atrium-Ventricle')
    pts_atr = m_atr.points()
    pts_whole = m_whole.points()
    av = pts_whole.view([('', pts_whole.dtype)] * pts_whole.shape[1]).ravel()
    cv = pts_atr.view([('', pts_atr.dtype)] * pts_atr.shape[1]).ravel()
    d_isin = np.isin(av,cv)
    index_atr = np.where(d_isin == True)[0]
    
    pts_classAnV = np.empty(len(pts_whole), dtype='object')
    pts_classAnV[:] = 'ventricle'
    pts_classAnV[index_atr] = 'atrium'
    del av, cv, pts_whole
    
    # pts_classAnV = classifyPtsMx(dict_planes = dict_planes, pl_name = 'pl2CutMesh_Chamber', pts_whole = pts_whole)
    df_classPts['AtrVent'] = pts_classAnV
    
    # Classify DnV
    pts_whole = m_whole.points()
    # - Classify DnV_Atr
    pts_classDnV_Atr = classifyPtsMx(dict_planes = dict_planes, pl_name = 'pl_AtrCoronal', pts_whole =  pts_whole)
    # - Classify DnV_Vent
    pts_classDnV_Vent = classifyPtsMx(dict_planes = dict_planes, pl_name = 'pl_VentCoronal', pts_whole =  pts_whole)

    atr = np.where(pts_classAnV == 'atrium')[0]
    vent = np.where(pts_classAnV == 'ventricle')[0]

    pts_classDnV_Atr[vent] = ''
    pts_classDnV_Vent[atr] = ''

    pts_classDnV = pts_classDnV_Atr+pts_classDnV_Vent
    df_classPts['DorsVent'] = pts_classDnV
    del pts_whole

    # Classify LnR
    print('--> Left-Right')
    pts_left = m_left.points()
    pts_whole = m_whole.points()
    av = pts_whole.view([('', pts_whole.dtype)] * pts_whole.shape[1]).ravel()
    bv = pts_left.view([('', pts_left.dtype)] * pts_left.shape[1]).ravel()
    # cint = np.intersect1d(av, bv).view(a.dtype).reshape(-1, a.shape[1])
    c_isin = np.isin(av,bv)
    index_left = np.where(c_isin == True)[0]

    pts_classLnR = np.empty(len(pts_whole), dtype='object')
    pts_classLnR[:] = 'right'
    pts_classLnR[index_left] = 'left'

    df_classPts['LeftRight'] = pts_classLnR

    for i, name, dat in zip(count(), names_data, data):
        df_classPts[name] = dat

    toc = perf_counter()
    time = toc-tic
    print("- All Done - points have been classified!\n- Sample of classified points")
    print(df_classPts.sample(10))
    print("- Time taken to classify = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")
    alert('whistle',1)

    if plot_show:
        plotPtClassif(filename, m_whole, pts_whole, [pts_classAnV, pts_classDnV, pts_classLnR])

    return df_classPts

#%% func - getPlaneNormals
def getPlaneNormals (no_planes, spline_pts):
    """
    Function that returns a list with normal vectors to create cutting planes. 
    Note: the input spline is checked in the inverse order, from last to first point

    Parameters
    ----------
    no_planes : int
        Number of planes that will be used to get transverse sections of heart
    spline_pts : list of coordinates
        List of centreline coordinates

    Returns
    -------
    normals : list of list of floats
        list of List with the x,y,z coordinates of each of the planes' normal
    pt_centre : list of list of floats
        list of List with the x,y,z coordinates of each of the planes' centre

    """

    normals = []
    pt_centre = []
    every = len(spline_pts)//no_planes
    #print(len(spline_pts), no_planes, every)
    list_index = list(range(len(spline_pts)-2,1,-every))
    # print(len(list_index))
    for i in list_index:
        pt_centre.append(spline_pts[i])
        normal = spline_pts[i-1]-spline_pts[i]
        normals.append(normal)

    return normals, pt_centre

#%% - SAVING
#%% func - saveMesh
def saveMesh(filename, mesh, mesh_name, dir_stl, extension):
    """
    Function to save mesh

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    mesh : mesh
        Mesh to save (vedo Mesh)
    mesh_name : str
        Mesh name
    dir_stl : path
        Path to the folder where the meshes are saved.
    extension : str
        Extension to saved mesh 'vtk'/'stl'.

    Returns
    -------
    None.

    """

    mesh_title = filename+"_"+mesh_name+"."+extension
    mesh_dir = os.path.join(dir_stl, mesh_title)
    mesh.write(mesh_dir)
    print("- Saved mesh - "+mesh_title+"!")
    alert("countdown",1)

#%% func - saveMeshes
def saveMeshes(filename, meshes, names, dict_colour, dir_stl, extension):
    """
    Function to save meshes

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    meshes : list of meshes
        List of meshes to save (vedo Meshes)
    names : list of str
        List of meshes names
    dict_colour : dictionary
        dictionary with information about the meshes colour
    dir_stl : path
        Path to the folder where the meshes are saved.
    extension : str
        Extension to saved mesh 'vtk'/'stl'.

    Returns
    -------
    dict_colour : dictionary
        Updated dictionary with information about the meshes colour

    """

    for i, mesh, name in zip(count(), meshes, names):
        mesh_title = filename+"_"+name+"."+extension
        mesh_dir = os.path.join(dir_stl, mesh_title)
        mesh.write(mesh_dir)

        sp_dict = dict_colour[name] = dict()
        sp_dict['colour'] = mesh.color()

    print("- All meshes have been saved!")
    alert("countdown",1)

    return dict_colour

#%% func - saveLayertMeshes
def saveLayerMeshes (filename, mesh_all, mesh_in, mesh_out, layer_name, dir_mesh, extension):
    """
    Function to save internal, external and actual heart layer mesh

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    mesh_all : mesh
        Mesh of the heart layer. (vedo Mesh)
    mesh_in : mesh
        Mesh of the filled internal contours of the heart layer. (vedo Mesh)
    mesh_out : mesh
        Mesh of the filled external contours of the heart layer. (vedo Mesh)
    layer_name : list of str
        List of meshes names
    dir_mesh : path
        Path to the folder where the meshes are saved.
    extension : str
        Extension to saved mesh 'vtk'/'stl'.

    Returns
    -------
    None.

    """

    print('- Saving '+ layer_name + ' meshes (All/In/Out) ')
    mesh_all_title = filename+"_"+layer_name+"."+extension
    mesh_all_dir = os.path.join(dir_mesh, mesh_all_title)
    mesh_in_title = filename+"_"+layer_name+"In."+extension
    mesh_in_dir = os.path.join(dir_mesh, mesh_in_title)
    mesh_out_title = filename+"_"+layer_name+"Out."+extension
    mesh_out_dir = os.path.join(dir_mesh, mesh_out_title)

    mesh_all.write(mesh_all_dir)
    mesh_in.write(mesh_in_dir)
    mesh_out.write(mesh_out_dir)

    print("- All meshes for "+layer_name+" were saved! -ext:"+extension)
    alert("countdown",1)


#%% func - saveThickness
def saveThickness(filename, arrays2save, names, dir2save):
    """
    Function to save distance information about mesh

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    arrays2save : array
        Array with distance data to save
    names : str
        Name of data to save
    dir2save : path
        Path to the folder where the dictionaries are saved.

    Returns
    -------
    None.

    """

    for i, arrayTh, name in zip(count(), arrays2save, names):
        thickness_title = filename+"_"+name
        thick_cj_dir = os.path.join(dir2save,thickness_title)
        np.save(thick_cj_dir, arrayTh)

    print("- All arrays (thickness/ballooning) have been saved!")
    alert("countdown",1)

#%% func - saveVideo
def saveVideo (filename, info, meshes4video, rotAngle, dir2save, alpha_cube, plotshow=True, zoom = 1, duration = 15):
    """
    Function that saves video of rotating mesh

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    info : str
        Text with additional information to print.
    meshes4video : list of meshes
        List of meshes to include in video
    rotAngle : float
        Heart inclination angle
    dir2save : path
        Path to the folder where the videos will be saved.
    alpha_cube : float
        Float number indicating the opacity of the scaling_cube
    plotshow : boolean, optional
        True if you want to see the resulting mesh in a plot, else False. The default is True.
    zoom : float
        Float number indicating zoom
    duration : int, optional
        Video duration in seconds. The default is 15.

    Returns
    -------
    None.

    """
    # alpha_cube = 0.1
    s_cube = 350

    if rotAngle > 90:
        rotAngle = 180 - rotAngle

    dORv = filename[9:10]
    print('- Saving video of ', info, '...')
    rotMeshes4video = []
    for j, mesh in enumerate(meshes4video):
        if dORv == 'V':
            rot_mesh = mesh.rotateX(rotAngle)
        elif dORv == 'D':
            rot_mesh = mesh.rotateY(90)
            rot_mesh = rot_mesh.rotateZ(rotAngle)
        rotMeshes4video.append(rot_mesh)
    
    scale_cube = Cube(pos=rotMeshes4video[0].centerOfMass(), side=s_cube, c='white', alpha=0.1).legend('Cube')
    rotMeshes4video.append(scale_cube)

    if plotshow:
        text = filename+"\n\n >> Rotated meshes"; txt = Text2D(text, c="k", font= font)
        vp1 = Plotter(N=1, axes=8)
        vp1.show(rotMeshes4video, txt, at=0, axes=8, zoom=1.2)

    rotMeshes4video[-1].alpha(alpha_cube)
    
    text2 = filename; txt2 = Text2D(text2, c="k", font= font)
    
    name = info.split(' ')
    settings.legendSize = .3
    vp = Plotter(bg='white', axes=10, offscreen=True)
    vp.show(rotMeshes4video, txt2, zoom=1.4)
    video_name = os.path.join(dir2save, filename+"_"+name[0]+".mp4")
    video = Video(video_name, duration=duration, backend='opencv')
    for i in range(180):
        vp.show(elevation=0, azimuth=2, zoom = zoom)  # render the scene
        video.addFrame()
    video.close()

    for k, mesh in enumerate(rotMeshes4video):
        if dORv == 'V':
            rot_mesh = mesh.rotateX(-rotAngle)
        elif dORv == 'D':
            rot_mesh = mesh.rotateZ(-rotAngle)
            rot_mesh = rot_mesh.rotateY(-90)

    alert('wohoo',1)

#%% func - getLocals  - REMOVE???
# def getLocals():
#     """
#     Function that returns a list of the meshes that are locals (open and saved in temporary memory)

#     Returns
#     -------
#     mesh_locals : list
#         List with the mesh names that are in locals.

#     """
#     meshes = ['myoc','endo','cj', 'cj_thickness','myoc_thickness','endo_thickness','myoc_intBall','myoc_extBall']
#     mesh_locals = []
#     for mesh in meshes: 
#         if mesh in locals():
#             mesh_locals.append(mesh)
            
#     return mesh_locals

#%% func- saveMultVideos
def saveMultVideos(filename, info, meshes4video, rangeThBall, rotAngle, dir2save, dir_txtNnpy, plotshow, alpha_cube = 0, zoom = 1, duration = 15):
    """
    Function that saves videos of rotating mesh given as input (meshes4video)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    info : str
        Text with additional information to print.
    meshes4video : list of list of meshes
        lisf of List of meshes to include in video
    rangeThBall : list of strs/tuples
        List of tuples indicating the min and max value for scalar range of thickness/ballooning meshes. 
    rotAngle : floats
        Heart inclination angle
    dir2save : path
        Path to the folder where the videos will be saved.
    dir_txtNnpy: 
        Path to the folder where the npy arrays with thickness are saved.
    plotshow : boolean
        True if you want to see the resulting mesh in a plot, else False.
    alpha_cube : float, optional
        Float number indicating the opacity of the scaling_cube. The default is 0.

    Returns
    -------
    None.

    """

    print('\n- Saving videos... this might take a while (about 2-3 min/video)')
    # bar = Bar('- Saving' , max = len(info), suffix = suffix, check_tty=False, hide_cursor=False)
    for i, name, mesh, scale in zip(count(), info, meshes4video, rangeThBall):
        # if 'thickness' in name or 'Ball' in name:
        if isinstance(scale, tuple): 
                [thickness] = loadNPY(filename = filename, names = [name], dir_txtNnpy = dir_txtNnpy, print_txt = False)
                if isinstance(mesh, list):
                    mesh[0].pointColors(thickness, cmap="jet", vmin=scale[0], vmax=scale[1])
                    mesh[0].addScalarBar()
                    mesh[0].alpha(1)
                    mesh[0].mapper().SetScalarRange(scale[0],scale[1])
                    mesh4video = mesh
                else: 
                    mesh.pointColors(thickness, cmap="jet", vmin=scale[0], vmax=scale[1])
                    mesh.addScalarBar()
                    mesh.alpha(1)
                    mesh.mapper().SetScalarRange(scale[0],scale[1])
                    mesh4video = [mesh]
        elif isinstance(mesh,list) and len(mesh)>1: 
            mesh4video = mesh
        else: 
            mesh4video = [mesh]
            
        saveVideo(filename = filename, info = name+' ['+str(i+1)+'/'+str(len(info))+']', meshes4video = mesh4video,
                        rotAngle = rotAngle, dir2save = dir2save, alpha_cube = alpha_cube, plotshow = plotshow, zoom = zoom, duration = duration)
        # bar.next()
    # bar.finish()
    alert('whistle',1)
    print('- All the videos have been saved!')

#%% - PLOTTERS
#%% func - plotPtClassif
def plotPtClassif(filename, mesh, pts_whole, pts_class):
    """
    Function that plots a subset of the classified points

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    mesh : mesh
        Mesh whose points were classified
    pts_whole : array of coordinates
        Array with x,y,z coordinates of whole mesh
    pts_class : list of str
        List with the final classification for points

    Returns
    -------
    None.

    """

    [pts_classAnV, pts_classDnV, pts_classLnR] = pts_class
    names_class = [['atrium', 'ventricle'],['dorsal', 'ventral'],['left', 'right']]
    color = [['tomato','gold'], ['greenyellow', 'darkgreen'],['deepskyblue', 'darkblue']]
    text = filename+'\n\n >> Point classification'
    txt = Text2D(text, c="k", font= font)

    if 'CJ' not in filename:
        every = 50
    else: 
        every = 20
        
    settings.legendSize = .2
    vp = Plotter(shape = (1, 4), axes = 1)#13
    vp.show(mesh, txt, at = 0)
    for i, sp_class in enumerate(pts_class):

        ind_first = np.where(sp_class == names_class[i][0])[0]
        pts_first = pts_whole[ind_first]
        sph_first = Spheres(pts_first[0::every,:], c = color[i][0], r=1).legend(names_class[i][0])


        ind_second = np.where(sp_class == names_class[i][1])[0]
        pts_second = pts_whole[ind_second]
        sph_second = Spheres(pts_second[0::every,:], c = color[i][1], r=1).legend(names_class[i][1])

        if i != 2:
            vp.show(mesh, sph_first, sph_second, at = i+1)
        else:
            vp.show(mesh, sph_first, sph_second, at = i+1, interactive = True)

#%% - PRINT
#%% func - decodeDict
def decodeDict (dict2classify, info):
    """
    Function to decode dictionary used to classify points

    Parameters
    ----------
    dict2classify : dict
        Dictionary with data to be classified
    info : list of str
        List of names of the data added to the dict (cj_thickness, myoc_intBall)

    Returns
    -------
    list
        List with information of the planes being added to the dictionary
        [AnV, DnV_Atr, DnV_Vent, pts_left, pts_whole, meas_param].

    """

    # Decode dictionary
    # - AnV
    d_AnV = dict2classify['AnV']['d']
    normal_AnV = dict2classify['AnV']['normal']
    classSorted_AnV = dict2classify['AnV']['classSorted']
    AnV = [d_AnV, normal_AnV, classSorted_AnV]
    # - DnV_Atr
    d_DnV_Atr = dict2classify['DnV_Atr']['d']
    normal_DnV_Atr = dict2classify['DnV_Atr']['normal']
    classSorted_DnV_Atr = dict2classify['DnV_Atr']['classSorted']
    DnV_Atr = [d_DnV_Atr, normal_DnV_Atr, classSorted_DnV_Atr]
    # - DnV_Vent
    d_DnV_Vent = dict2classify['DnV_Vent']['d']
    normal_DnV_Vent = dict2classify['DnV_Vent']['normal']
    classSorted_DnV_Vent = dict2classify['DnV_Vent']['classSorted']
    DnV_Vent = [d_DnV_Vent, normal_DnV_Vent, classSorted_DnV_Vent]

    # Pts to classify
    pts_left = np.asarray(dict2classify['pts_Left'])
    pts_whole = dict2classify['pts_Whole']

    meas_param = []
    for i, inf in enumerate(info):
        param = np.asarray(dict2classify['param_'+inf])
        meas_param.append(param)

    return [AnV, DnV_Atr, DnV_Vent, pts_left, pts_whole, meas_param]

#%% - ALERT WHEN IMPORTED
print ("IMPORTED: morphoHeart_funcMeshes")
alert('jump',1)

#%% Others (back-up / developing)
# A --------
# import numpy as np
# def random_arr(low, high, size, rand_type):
#     if rand_type == float:
#         arr = [np.random.uniform(low, high) for _ in range(size)]
#     elif rand_type == int:
#         arr = [np.random.randint(low, high) for _ in range(size)]
#     return  arr

# B --------
# Cut layer with n number of planes perpendicular to spline and plot
# mesh2cut = mesh_ch1
#Transform mesh from vtk structure to trimesh
# print("Transforming mesh into trimesh object...")
# mesh3= vtk2trimesh(mesh2cut)
# cjf4.alert("jump",1)

# Cut layer with a certain number of planes
# no_planes = 15
# sectDict = cjf4.sliceLayerCL(mesh3, no_planes, cl_splineAdd, 'red')
# cjf4.plotAllLayerSlices(mesh2cut, cl_splineAdd, no_planes, sectDict)

# Plot a particular slice
# cjf4.plotLayerSlice(mesh2cut, cl_splineAdd, 5, sectDict)

# #%% func - sliceLayerCL
# def sliceLayerCL(mesh, no_planes, spline, color):
#     sections = dict()
#     allPlanes = []
#     allSlice_2D = []
#     allSlice_3D = []
#     allCenters = []
    
#     # Get normals and centres of planes
#     normals, pt_centre = getPlaneNormals (no_planes, spline.points())
    
#     for slice_num in range(no_planes):
#         print("Processing slice no.", slice_num+1, "/", no_planes)
#         # Get mesh cuts
#         plane, slice_2D, slice_3D, sphere, mslice = getSlices(mesh, normals[slice_num], pt_centre[slice_num], color)
#         allPlanes.append(plane)
#         allSlice_2D.append(slice_2D)
#         allSlice_3D.append(slice_3D)
#         allCenters.append(sphere)
#         print("Mesh sliced...")
        
#     print("Mesh has been sliced!")
        
#     sections['Planes'] = allPlanes
#     sections['Spheres'] = allCenters
#     sections['slice2D'] = allSlice_2D
#     sections['slice3D'] = allSlice_3D
    
#     return sections

# #%% func - plotAllLayerSlices
# def plotAllLayerSlices (mesh, spline, no_planes, sections):
#     N = no_planes*2
#     n_plot1 = list(range(0,no_planes,1))
#     n_plot2 = list(range(no_planes,no_planes*2,1))
#     vpA = Plotter(N=N, sharecam=False) 
#     for n in range(no_planes):
#         plane = sections['Planes'][n]
#         sphere = sections['Spheres'][n]
#         slice3D = sections['slice3D'][n]
#         slice2D = sections['slice2D'][n]
#         vpA.show(mesh, plane, slice3D, spline, sphere, at=n_plot1[n], zoom=2)
#         if n == no_planes-1:
#             vpA.show(slice2D, at=n_plot2[n], interactive=1)
#         else:
#             vpA.show(slice2D, at=n_plot2[n])

# #%% func - plotSlice 
# def plotLayerSlice (mesh, spline, slice_num, sections):
#     plane = sections['Planes'][slice_num]
#     slice3D = sections['slice3D'][slice_num]
#     slice2D = sections['slice2D'][slice_num]
#     sphere = sections['Spheres'][slice_num]
    
#     vpA = Plotter(N=2, sharecam=False) 
#     vpA.show(mesh, plane, slice3D, spline, sphere, at=0, zoom=2)
#     vpA.show(slice2D, at=1, interactive=1)

# C --------
# getBalloons
#sph_all, sph_whole = fcMeshes.getBalloonedHeart(myoc_int_npcl)

# vp = Plotter(N=1, axes=13)
# vp.show(m_myoc.alpha(0.01), sph_all[0].color('green').alpha(0.05), sph_all[5].alpha(0.05), at=0, interactive=True)

# sph_whole = booleanOperation(sph_all[0], 'plus',sph_all[2])
# sph_whole.color('coral')

# vp = Plotter(N=1, axes=13)
# vp.show(m_myoc.alpha(0.01), sph_whole, at=0, interactive=True)

#%% func - unloopChambers - Old but previous version to current
# def unloopChambersF(filename, mesh, kspl_CL, kspl_ext, no_planes, pl_CLRibbon, param, param_name, df_AtrVent, dict_kspl, 
#                    dict_shapes, dict_pts, dict_planes, dir_results, plotshow = False, tol = 0.05):
#     """
#     Function to unloop the heart and get two dataframes with the data for the unlooped hearts

#     Parameters
#     ----------
#     filename : str
#         Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
#     mesh : mesh
#         Color coded mesh that will be unlooped
#     kspl_CL : Kspline
#         Centreline (vedo KSpline)
#     kspl_ext : Kspline
#         Extended centreline (vedo KSpline).
#     no_planes : int
#         Number of planes that will be used to get transverse sections of each chamber (heart = chamber*2)
#     pl_CLRibbon :  Plane
#         Plane used to extend dorso-ventrally the centreline
#     param : list of arrays of floats
#         List of numpy arrays with distance measurements that want to be unlooped
#     param_name : str
#         Name of distance parameter(s) to add to saving file
#     df_AtrVent : numpy array of objects 
#         Array with points classification of atrium and ventricle
#     dict_kspl : dictionary
#         Initialised dictionary with ksplines information
#     dict_shapes : dictionary
#         Initialised dictionary with shapes information
#     dict_pts : dictionary
#         Initialised dictionary with points information
#     dict_planes : dictionary
#         Initialised dictionary with planes information
#     dir_results : path
#         Path to the folder where the results are saved.
#     plotshow : boolean, optional
#         True if you want to see the resulting mesh in a plot, else False. The default is False.
#     tol : float, optional
#         Tolerance defined to get points in plane. The default is 0.05.

#     Returns
#     -------
#     df_unlooped : dataframe
#         Dataframe containing unlooped heart data.
#     kspl_vSurf : Kspline
#         KSpline in the surface of the heart indicating the 0 degree direction.
#     kspls_HR : List of KSplines
#         List of ksplines for atrium and ventricle with high resolution.
#     arr_all : list of arrows 
#         List of arrows from the centreline to the ventral centre point of plane
#     arr_valve : list of arrows 
#         List of arrows from the centreline to the valve
#     dict_pts : dictionary
#         Resulting dictionary with points information updated
#     dict_shapes : dictionary
#         Resulting dictionary with shapes information updated
#     dict_kspl : dictionary
#         Resulting dictionary with ksplines information updated

#     """

#     print('\n\n- Unlooping the heart chambers...')
#     dfs_unlooped = []
#     planes_all = []
#     spheres_all = []
#     if plotshow:
#         plotevery = no_planes // 3
#         print('- Plotting every X number of planes:', plotevery)
    
#     # KSpline on surface
#     kspl_vSurf = getExtCLonSurf(filename, mesh, kspl_ext, pl_CLRibbon)

#     # Create new extended and cut kspline with higher resolution
#     kspl_CLnew = getExtCLHighRes(filename, mesh, kspl_ext, kspl_CL, dict_planes) 

#     # # Create kspline for each chamber
#     kspl_vent, kspl_atr, dict_pts, sph_cut =  ksplChamberCut(mesh, kspl_CLnew, dict_shapes, dict_pts)
#     # Find position within new kspline cut by disc used to cut chambers (num_pt new)
#     num_pt = dict_pts['numPt_CLChamberCutHighRes']
    
#     #Add all ksplines created to dictionary
#     dict_kspl = addKSplines2Dict(kspls = [kspl_ext, kspl_vSurf, kspl_atr, kspl_vent], info = ['','','',''], dict_kspl = dict_kspl)
    
#     vp = Plotter(N=3, axes=4)
#     vp.show(kspl_CLnew, kspl_vSurf, kspl_ext, mesh, sph_cut, at = 0)
#     # vp.show(ksplCL_cut, arr_vectPlCut, kspl_vSurf, mesh, sph_cut, at = 0)
#     vp.show(mesh, kspl_vent, sph_cut, at=1)
#     vp.show(mesh, kspl_atr, sph_cut, at=2, interactive = True)
    
#     # NOW UNLOOP EACH CHAMBER
#     kspls_HR = [kspl_atr, kspl_vent]
#     texts = ['- Unlooping atrium', '- Unlooping ventricle']
#     names= ['unloopAtr', 'unloopVent']
#     arr_all = []
#     arr_valve = []
#     for ksp_num, kspl in enumerate(kspls_HR):
#         # Create matrix with all data
#         #0:x, 1:y, 2:z, 3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: parameters
#         matrix_unlooped = np.zeros((len(mesh.points()),9))
#         matrix_unlooped[:,0:3] = mesh.points()
#         matrix_unlooped[:,7] = param[0]
#         matrix_unlooped[:,8] = param[1]
    
#         # Get normals and centres of planes to cut heart (number of planes given as input+2)
#         pl_normals, pl_centres = getPlaneNormals(no_planes = no_planes, spline_pts = kspl.points())
#         # Give a number between 1-2 to each atrial plane, and between 0-1 to each ventricular plane
#         plane_num = np.linspace(2-ksp_num,1-ksp_num,len(pl_normals))
#         if ksp_num == 1:
#             plane_num = plane_num[::-1]
#         else: 
#             index_guess = len(kspl_CLnew.points())-1
#         planes = []
#         spheres = []
        
#         bar = Bar(texts[ksp_num], max=len(pl_normals), suffix = suffix, check_tty=False, hide_cursor=False)
#         # Iterate through each plane
#         for i, normal, centre in zip(count(), pl_normals, pl_centres):
#             if ksp_num == 1:
#                 normal = -normal
#             # A. Get cut plane info
#             # - Info Plane (Plane orientation/normal representation)
#             arr_vectPlCut = Arrow(centre, centre+normal*20, s = 0.1, c='orange')
#             # Cut cl with plane
#             ksplCL_cut = kspl_CLnew.clone().cutWithMesh(Plane(pos=centre, normal=normal, sx=300), invert=True).lw(5).color('tomato')
            
#             # Find point of centreline closer to last point of kspline that was cut by plane (ksplCL_cut)
#             # If kspl_CLnew wasn't cut, then initialise it with last index of cl
#             if len(ksplCL_cut.points()) == 0:
#                 pts_o = [kspl_CLnew.points()[-1]] #???
#             # Else, then initialise it with all the points of the cut cl
#             else: 
#                 pts_o = ksplCL_cut.points(); 
#                 # print('len: ',len(pts_o))
            
#             # If the ventricle is the one being unlooped, then initialise index as length of centreline
#             if ksp_num == 1:
#                 index_guess= len(ksplCL_cut.points())
#             # print('index_guess:',index_guess)
            
#             #Find closest point between the pts_o and the high resolution centreline, initialising it using th index_guess prev defined
#             pt_out, pt_num = findClosestPtGuess(pts_o, kspl_CLnew.points(), index_guess)
#             index_guess = pt_num
            
#             #Find closest point between the pts_o and the vSurface centreline
#             if ksp_num == 0 and i ==0:
#                 ind_vSurf = -1
#             elif ksp_num ==1 and i ==0:
#                 ind_vSurf = 0
#             pt_vSurf, ind_vSurf = findClosestPtGuess(pts_o, kspl_vSurf.points(), ind_vSurf)
            
#             print(pt_vSurf)
            
#             if pt_num < num_pt:
#                 chamber = 'ventricle'
#             else: 
#                 chamber = 'atrium'
#             # print('pt_num:',pt_num, 'num_pt:',num_pt, chamber, pt_out)
            
#             try: 
#                 if ksp_num == 0 and i == 0:
#                     tol2use = 0.01
#                 else:
#                     tol2use = tol
#                 # print('tol: ', tol2use)
                
#                 # B. Get vector that defines 0 deg angle
#                 pts_zero, _ = getPointsAtPlane(points = kspl_vSurf.points(), pl_normal = normal,
#                                                 pl_centre = centre, tol=0.25)
#                 if ksp_num == 0:# and i == 0:
#                     pt4dist = centre
                    
#                 # pt4dist = pt_vSurf
#                 min_dist = 100000
#                 for ik, pt in enumerate(pts_zero):
#                     dist = findDist(pt,pt4dist)
#                     if dist < min_dist:
#                         ind = ik; min_dist = dist
#                 pt_zero = pts_zero[ind]
#                 pt4dist = pt_zero
                
#                 # Vector from centre to cl_surface point being cut by plane
#                 v_zero = unit_vector(pt_zero - centre)
        
#                 # C. Get points of mesh at plane
#                 d_points = np.absolute(np.dot(np.subtract(matrix_unlooped[:,0:3],np.asarray(centre)),np.asarray(normal)))
#                 # Find the indexes of the points that have not been yet taken, are at the plane and are in the 
#                 # chamber being analysed
#                 index_ptsAtPlane = np.where((d_points <= tol2use) & (matrix_unlooped[:,3] == 0) & (df_AtrVent == chamber))
#                 # print(d_points.min(), d_points.max(),'-lenptsatplane:',len(index_ptsAtPlane[0]))
                
#                 # Define new matrix just with the points on plane
#                 new_matrix = matrix_unlooped[index_ptsAtPlane,:][0]
#                 # - Get points of mesh that are on plane, centered on centreline point
#                 ptsC = np.subtract(new_matrix[:,0:3],np.asarray(centre))
#                 # - Get the radius of those points
#                 radius = [np.linalg.norm(x) for x in ptsC]
    
#                 # D. Find direction of point with respect to plane that includes central point, vC and the normal of the cutting plane
#                 # Vector normal to plane normal and v_zero (vector defined from centre of plane to cut-pt in centreline surface)
#                 normal_divLR = np.cross(normal, v_zero)
#                 # Define using these vectors if the points are all lying in the same side or not (1, -1)
#                 lORr = np.sign(np.dot(ptsC, np.asarray(normal_divLR)))
    
#                 # E. Get angle of points in that plane using v_zero
#                 av = np.dot(ptsC,v_zero)
#                 cosTheta = np.divide(av, radius) # vectors of magnitude one
#                 theta = np.arccos(cosTheta)*180/np.pi
#                 theta_corr = np.multiply(lORr, theta)
    
#                 pl_cut = Plane(pos = centre, normal = normal, sx = 300).color('light steel blue').alpha(0.5)
#                 sph_C = Sphere(centre, r=2, c='red')
#                 arr_vectXYZ = Arrow(centre, pt_zero, s = 0.1)
#                 arr_valv = Arrow(centre, normal_divLR-centre, s = 0.1)
#                 arr_all.append(arr_vectXYZ)
#                 arr_valve.append(arr_valv)
#                 spheres_cut = []
#                 if plotshow:
#                     # Create stuff to plot
#                     sph_C = Sphere(centre, r=2, c='red')
#                     sph_VXYZ = Sphere(pt_zero, r=2, c='green')
#                     arr_vectXYZ = Arrow(centre, pt_zero, s = 0.1)
#                     if i % plotevery == 5:
#                         sphL = []; sphR = []
#                         for num, pt in enumerate(ptsC):
#                             if num % 50 == 0:
#                                 if lORr[num] == 1:
#                                     sphL.append(Sphere(pt+centre, r=2, c='blueviolet'))
#                                 else:
#                                     sphR.append(Sphere(pt+centre, r=2, c='gold'))
#                         settings.legendSize = .3
#                         text = filename+"\n\n >> Unlooping the heart"
#                         txt = Text2D(text, c=c, font=font)
#                         vp = Plotter(N=1, axes=4)
#                         # vp.show(mesh,sphL, sphR, kspl_vSurf, kspl, kspl_ext, arr_vectXYZ, arr_vectPlCut, sph_C, sph_VXYZ, pl_cut, txt, at=0, azimuth = 0, interactive=1)
#                         vp.show(mesh,sphL, sphR, kspl_vSurf, kspl, kspl_ext, ksplCL_cut, arr_vectXYZ, arr_vectPlCut, sph_C, sph_VXYZ, pl_cut, txt, at=0, azimuth = 0, interactive=1)
#                         spheres_cut.append([sphL,sphR])
#                 if i % 30 == 0:
#                     planes.append(pl_cut)
#                     spheres.append(sph_C)
#                 # - Save all obtained values in matrix_unlooped
#                 for num, index in enumerate(index_ptsAtPlane[0]):
#                     #3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: param
#                     matrix_unlooped[index,3] = 1
#                     matrix_unlooped[index,4] = plane_num[i]
#                     matrix_unlooped[index,5] = theta_corr[num]
#                     matrix_unlooped[index,6] = radius[num]
                
#             except: 
#                 pass
#             bar.next()

#         dict_shapes = addShapes2Dict (shapes = [kspl_vSurf, kspl_CLnew, kspl_atr, kspl_vent], dict_shapes = dict_shapes, radius = [[],[],[],[]], print_txt = False)
        
#         df_unlooped = pd.DataFrame(matrix_unlooped, columns=['x','y','z','taken','z_plane','theta','radius','cj_thickness', 'myoc_intBall'])
#         df_unlooped = df_unlooped[df_unlooped['taken']==1]
#         df_unlooped_f = df_unlooped.drop(['x', 'y','z'], axis=1)
#         df_unlooped_f.astype({'taken': 'bool','z_plane':'float16','theta':'float16','radius':'float16','cj_thickness':'float16','myoc_intBall':'float16' }).dtypes
#         print('\n')
#         saveDF(filename = filename, df2save = df_unlooped_f, df_name = 'df_'+names[ksp_num]+param_name,
#                         dir2save = dir_results)
#         dfs_unlooped.append(df_unlooped_f)
#         planes_all.append(planes)
#         spheres_all.append(spheres)
        
#     bar.finish() 
        
#     if plotshow:
#         vp = Plotter(N=2, axes = 4)
#         vp.show(mesh, kspl_vSurf, kspls_HR[0], planes_all[0], spheres_all[0], Text2D(texts[0], c="k", font= font), at=0)
#         vp.show(mesh, kspl_vSurf, kspls_HR[1], planes_all[1], spheres_all[1], Text2D(texts[1], c="k", font= font), at=1, interactive=True)
        
#     #return df_unlooped, kspl_vSurf, kspls_HR, arr_all, arr_valve, dict_pts, dict_shapes, 
#     return df_unlooped, kspl_vSurf, kspls_HR, arr_all, spheres_cut, dict_pts, dict_shapes, dict_kspl

#%% func - unloopChambers_OLD
# def unloopChambers_OLD(filename, mesh, kspl_CL, kspl_ext, no_planes, pl_CLRibbon, param, param_name, df_AtrVent, dict_kspl, 
#                    dict_shapes, dict_pts, dict_planes, dir_results, plotshow = False, tol = 0.05):
#     """
#     Function to unloop the heart and get two dataframes with the data for the unlooped hearts

#     Parameters
#     ----------
#     filename : str
#         Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
#     mesh : mesh
#         Color coded mesh that will be unlooped
#     kspl_CL : Kspline
#         Centreline (vedo KSpline)
#     kspl_ext : Kspline
#         Extended centreline (vedo KSpline).
#     no_planes : int
#         Number of planes that will be used to get transverse sections of each chamber (heart = chamber*2)
#     pl_CLRibbon :  Plane
#         Plane used to extend dorso-ventrally the centreline
#     param : list of arrays of floats
#         List of numpy arrays with distance measurements that want to be unlooped
#     param_name : str
#         Name of distance parameter(s) to add to saving file
#     df_AtrVent : numpy array of objects 
#         Array with points classification of atrium and ventricle
#     dict_kspl : dictionary
#         Initialised dictionary with ksplines information
#     dict_shapes : dictionary
#         Initialised dictionary with shapes information
#     dict_pts : dictionary
#         Initialised dictionary with points information
#     dict_planes : dictionary
#         Initialised dictionary with planes information
#     dir_results : path
#         Path to the folder where the results are saved.
#     plotshow : boolean, optional
#         True if you want to see the resulting mesh in a plot, else False. The default is False.
#     tol : float, optional
#         Tolerance defined to get points in plane. The default is 0.05.

#     Returns
#     -------
#     df_unlooped : dataframe
#         Dataframe containing unlooped heart data.
#     kspl_vSurf : Kspline
#         KSpline in the surface of the heart indicating the 0 degree direction.
#     kspls_HR : List of KSplines
#         List of ksplines for atrium and ventricle with high resolution.
#     arr_all : list of arrows 
#         List of arrows from the centreline to the ventral centre point of plane
#     arr_valve : list of arrows 
#         List of arrows from the centreline to the valve
#     dict_pts : dictionary
#         Resulting dictionary with points information updated
#     dict_shapes : dictionary
#         Resulting dictionary with shapes information updated
#     dict_kspl : dictionary
#         Resulting dictionary with ksplines information updated

#     """

#     print('\n\n')
#     dfs_unlooped = []
#     planes_all = []
#     spheres_all = []
#     if plotshow:
#         plotevery = no_planes // 5
#         print('- Plotting every X number of planes:', plotevery)

#     # - Get unitary normal of plane to create CL_ribbon
#     pl_normCLRibbon = unit_vector(pl_CLRibbon['pl_normal'])
    
#     # Increase the resolution of the extended centreline and interpolate to unify sampling
#     xd = np.diff(kspl_ext.points()[:,0])
#     yd = np.diff(kspl_ext.points()[:,1])
#     zd = np.diff(kspl_ext.points()[:,2])
#     dist = np.sqrt(xd**2+yd**2+zd**2)
#     u = np.cumsum(dist)
#     u = np.hstack([[0],u])
#     t = np.linspace(0, u[-1], 1000)
#     resamp_pts = interpn((u,), kspl_ext.points(), t)
#     kspl_ext = KSpline(resamp_pts, res = 1000).lw(5).color('deeppink').legend('kspl_extHR')
    
#     #Find the points that intersect with the ribbon
#     pts_int = []
#     for num in range(len(kspl_ext.points())):
#         try: 
#             cl_pt_test = kspl_ext.points()[num]
#             pt_int = mesh.intersectWithLine(cl_pt_test, cl_pt_test+60*pl_normCLRibbon)
#             rad_pts = [np.linalg.norm(x- cl_pt_test) for x in pt_int]
#             # if len(rad_pts)>1:
#                 # print(rad_pts)
#             ind_pt = np.where(rad_pts == max(rad_pts))[0][0]
#             # print(ind_pt)
#             pts_int.append(pt_int[ind_pt])
#             # cl_pt_test = kspl_ext.points()[num]
#             # pt_int = mesh.intersectWithLine(cl_pt_test, cl_pt_test+60*pl_normCLRibbon)
#             # if len(pts_int) > 0:
#             #     dist2last_pts = [np.linalg.norm(x- pst_int[-1]) for x in pt_int]
#             #     ind_pt = np.where(dist2last_pts == min(dist2last_pts))[0][0]
#             # else: 
#             #     rad_pts = [np.linalg.norm(x- cl_pt_test) for x in pt_int]
#             #     ind_pt = np.where(rad_pts == max(rad_pts))[0][0]
#             # pts_int.append(pt_int[ind_pt])
#         except: 
#             # print('exc')
#             if num > 750:
#                 dist_pts = kspl_ext.points()[num] - kspl_ext.points()[num-1]
#                 try: 
#                     pt_int = pts_int[-1]+dist_pts
#                     pts_int.append(pt_int)
#                 except:
#                     # print('pass')
#                     pass
#             else: 
#                 pass
    
#     # KSpline on surface
#     kspl_vSurf = KSpline(pts_int).color('black').lw(4).legend('kspl_VSurfaceIntMyoc')
    
#     #Create kspline_ext cut with inf and outflow with higher resolution
#     try: 
#         add_pts = 50
#         kspl_CLnew = kspl_ext.clone()
#         # Cut with inflow plane
#         kspl_CLnew_cutIn = kspl_ext.clone().cutWithMesh(Plane(pos=dict_planes['pl2CutMesh_inflow']['pl_centre'], normal=dict_planes['pl2CutMesh_inflow']['pl_normal'], sx=300), invert=True).color('navy')
#         _, num_pt_inf = findClosestPtGuess(kspl_CLnew_cutIn.points(), kspl_ext.points(), index_guess = -1)
#         kspl_CLnew_cutOut = kspl_ext.clone().cutWithMesh(Plane(pos=dict_planes['pl2CutMesh_outflow']['pl_centre'], normal=dict_planes['pl2CutMesh_outflow']['pl_normal'], sx=300), invert=True).color('lime')
#         _, num_pt_outf =  findClosestPt(kspl_CLnew_cutOut.points()[-1], kspl_ext.points())#, index_guess = 0)
#         # print(num_pt_outf, num_pt_inf)
        
#         if (num_pt_outf-add_pts) < 0:
#             ind_outf = 0
#         else: 
#             ind_outf = num_pt_outf - add_pts
    
#         if (num_pt_inf+add_pts) > len(kspl_ext.points()):
#             ind_inf = len(kspl_ext.points())
#         else: 
#             ind_inf = num_pt_inf+add_pts
#         print(ind_outf, ind_inf)
    
#         kspl_CLnew = KSpline(kspl_ext.points()[ind_outf:ind_inf], res = 600).lw(5).color('deeppink').legend(kspl_CL._legend+'_HighRes')
#         # print('try1')
#     except: 
#         kspl_CLnew = KSpline(points = kspl_CL.points(), res= 600).color('deepskyblue').lw(5).legend(kspl_CL._legend+'_HighRes')
#         # print('try2')
    
#     # print('A')
#     #Get num_pt new
#     ksplCL_cut = kspl_CLnew.clone().cutWithMesh(Plane(pos=dict_shapes['cyl2CutChambers_o']['cyl_centre'], normal=dict_shapes['cyl2CutChambers_o']['cyl_axis'], sx=300), invert=True)
#     # ksplCL_cut = kspl_ext.clone().cutWithMesh(Plane(pos=cyl_centre, normal=cyl_normal, sx=300), invert=True)
#     # Find point of centreline closer to last point of kspline cut
#     # _, num_pt = findClosestPt(ksplCL_cut.points()[-1], kspl_CLnew.points())
#     _, num_pt = findClosestPtGuess(ksplCL_cut.points(), kspl_CLnew.points(),300)
#     print(num_pt)
#     # print(num_pt)
#     # Add pt to dict
#     sph_cut = Sphere(pos = kspl_CLnew.points()[num_pt], r=4, c='gold').legend('sph_ChamberCutCLHighRes')
#     dict_pts = addPoints2Dict(spheres = [sph_cut], info = [''], dict_pts = dict_pts)
#     dict_pts['numPt_CLChamberCutHighRes'] = num_pt
    
#     print('B')
#     # Create kspline for atrium and ventricle
#     n_vent_pts = 0
#     res_v = 595
#     while n_vent_pts < 610:
#         res_v +=1
#         kspl_vent = KSpline(points = kspl_CLnew.points()[0:num_pt], res = res_v).color('tomato').lw(8).legend('kspl_vent')
#         n_vent_pts = kspl_vent.NPoints()
#         print('vent')
        
#     n_atr_pts = 0
#     res_a = 595
#     while n_atr_pts < 610:
#         res_a +=1
#         kspl_atr = KSpline(points = kspl_CLnew.points()[num_pt:], res = res_a).color('goldenrod').lw(8).legend('kspl_atr')
#         n_atr_pts = kspl_atr.NPoints()
#         # print('n_atr_pts:', n_atr_pts)
        
#     # sph_1st_vent = Sphere(pos = kspl_vent.points()[0], r = 3, c = 'lime')
#     # sph_1st_atr = Sphere(pos = kspl_atr.points()[0], r = 3, c = 'navy')
    
#     dict_kspl = addKSplines2Dict(kspls = [kspl_ext, kspl_vSurf, kspl_atr, kspl_vent], info = ['','','',''], dict_kspl = dict_kspl)
#     print('C')
#     vp = Plotter(N=3, axes=4)
#     vp.show(kspl_CLnew, kspl_atr, kspl_vent, kspl_vSurf, kspl_ext, mesh, sph_cut, at = 0)
#     vp.show(mesh, kspl_CLnew_cutIn, kspl_CLnew,  at=1)
#     vp.show(mesh, kspl_CLnew_cutOut, kspl_ext, at=2, interactive = True)
    
#     kspls_HR = [kspl_atr, kspl_vent]
#     texts = ['- Unlooping atrium', '- Unlooping ventricle']
#     names= ['unloopAtr', 'unloopVent']
#     arr_all = []
#     arr_valve = []
#     for ksp_num, kspl in enumerate(kspls_HR):
#         # Create matrix with all data
#         #0:x, 1:y, 2:z, 3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: parameters
#         matrix_unlooped = np.zeros((len(mesh.points()),9))
#         matrix_unlooped[:,0:3] = mesh.points()
#         matrix_unlooped[:,7] = param[0]
#         matrix_unlooped[:,8] = param[1]
    
#         # Get normals and centres of planes
#         pl_normals, pl_centres = getPlaneNormals(no_planes = no_planes, spline_pts = kspl.points())
#         # pl_normals = pl_normals[1:-2]
#         # pl_centres = pl_centres[1:-2]
#         plane_num = np.linspace(2-ksp_num,1-ksp_num,len(pl_normals))
#         print('D')
#         planes = []
#         spheres = []
#         bar = Bar(texts[ksp_num], max=len(pl_normals), suffix = suffix, check_tty=False, hide_cursor=False)
#         index_guess = len(kspl_CLnew.points())-1
#         # Iterate through each plane
#         for i, normal, centre in zip(count(), pl_normals, pl_centres):
#             # A. Get cut plane info
#             # - Info Plane
#             arr_vectPlCut = Arrow(centre, centre+normal*20, s = 0.1, c='orange')
#             # Cut cl with plane
#             ksplCL_cut = kspl_CLnew.clone().cutWithMesh(Plane(pos=centre, normal=normal, sx=300), invert=True).lw(5).color('tomato')
#             # Find point of centreline closer to last point of kspline cut
#             # print(len(ksplCL_cut.points()), len(kspl_CLnew.points()), index_guess)
#             # print('val i:', i)
#             if len(ksplCL_cut.points()) == 0:
#                 pts_o = [kspl_CLnew.points()[-1]] #???
#             else: 
#                 pts_o = ksplCL_cut.points()
#             _, pt_num = findClosestPtGuess(pts_o, kspl_CLnew.points(), index_guess)
#             index_guess = pt_num
            
#             if pt_num < num_pt:
#                 chamber = 'ventricle'
#             else: 
#                 chamber = 'atrium'
#             print(pt_num, num_pt, chamber)
            
#             try: 
#                 print('E')
#                 if ksp_num == 0 and i == 0:
#                     tol2use = 0.01
#                 else:
#                     tol2use = tol
#                 # print('tol: ', tol2use)
                
#                 # B. Get vector that defines 0 deg angle
#                 pts_zero, _ = getPointsAtPlane(points = kspl_vSurf.points(), pl_normal = normal,
#                                                 pl_centre = centre, tol=0.25)
#                 min_dist = 100000
#                 for ik, pt in enumerate(pts_zero):
#                     dist = findDist(pt, centre)
#                     if dist < min_dist:
#                         ind = ik; min_dist = dist
#                 pt_zero = pts_zero[ind]
                
#                 # Vector from centre to cl_surface point being cut by plane
#                 v_zero = unit_vector(pt_zero - centre)
        
#                 # C. Get points of mesh at plane
#                 d_points = np.absolute(np.dot(np.subtract(matrix_unlooped[:,0:3],np.asarray(centre)),np.asarray(normal)))
#                 # Find the indexes of the points that have not been yet taken, are at the plane and are in the 
#                 # chamber being analysed
#                 index_ptsAtPlane = np.where((d_points <= tol2use) & (matrix_unlooped[:,3] == 0) & (df_AtrVent == chamber))
#                 # print(d_points.min(), d_points.max(),'-lenptsatplane:',len(index_ptsAtPlane[0]))
                
#                 # Define new matrix just with the points on plane
#                 new_matrix = matrix_unlooped[index_ptsAtPlane,:][0]
#                 # - Get points of mesh that are on plane, centered on centreline point
#                 ptsC = np.subtract(new_matrix[:,0:3],np.asarray(centre))
#                 # - Get the radius of those points
#                 radius = [np.linalg.norm(x) for x in ptsC]
    
#                 # D. Find direction of point with respect to plane that includes central point, vC and the normal of the cutting plane
#                 # Vector normal to plane normal and v_zero (vector defined from centre of plane to cut-pt in centreline surface)
#                 normal_divLR = np.cross(normal, v_zero)
#                 # Define using these vectors if the points are all lying in the same side or not (1, -1)
#                 lORr = np.sign(np.dot(ptsC, np.asarray(normal_divLR)))
    
#                 # E. Get angle of points in that plane using v_zero
#                 av = np.dot(ptsC,v_zero)
#                 cosTheta = np.divide(av, radius) # vectors of magnitude one
#                 theta = np.arccos(cosTheta)*180/np.pi
#                 theta_corr = np.multiply(lORr, theta)
    
#                 pl_cut = Plane(pos = centre, normal = normal, sx = 300).color('light steel blue').alpha(0.5)
#                 sph_C = Sphere(centre, r=2, c='red')
#                 arr_vectXYZ = Arrow(centre, pt_zero, s = 0.1)
#                 arr_valv = Arrow(centre, normal_divLR-centre, s = 0.1)
#                 arr_all.append(arr_vectXYZ)
#                 arr_valve.append(arr_valv)
#                 if plotshow:
#                     # Create stuff to plot
#                     sph_C = Sphere(centre, r=2, c='red')
#                     sph_VXYZ = Sphere(pt_zero, r=2, c='green')
#                     arr_vectXYZ = Arrow(centre, pt_zero, s = 0.1)
#                     if i % plotevery == 5:
#                         sphL = []; sphR = []
#                         for num, pt in enumerate(ptsC):
#                             if num % 50 == 0:
#                                 if lORr[num] == 1:
#                                     sphL.append(Sphere(pt+centre, r=2, c='blueviolet'))
#                                 else:
#                                     sphR.append(Sphere(pt+centre, r=2, c='gold'))
#                         settings.legendSize = .3
#                         text = filename+"\n\n >> Unlooping the heart"
#                         txt = Text2D(text, c=c, font=font)
#                         vp = Plotter(N=1, axes=4)
#                         # vp.show(mesh,sphL, sphR, kspl_vSurf, kspl, kspl_ext, arr_vectXYZ, arr_vectPlCut, sph_C, sph_VXYZ, pl_cut, txt, at=0, azimuth = 0, interactive=1)
#                         vp.show(mesh,sphL, sphR, kspl_vSurf, kspl, kspl_ext, ksplCL_cut, arr_vectXYZ, arr_vectPlCut, sph_C, sph_VXYZ, pl_cut, txt, at=0, azimuth = 0, interactive=1)
                
#                 if i % 40 == 0:
#                     planes.append(pl_cut)
#                     spheres.append(sph_C)
#                 print('F')
#                 # - Save all obtained values in matrix_unlooped
#                 for num, index in enumerate(index_ptsAtPlane[0]):
#                     #3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: param
#                     matrix_unlooped[index,3] = 1
#                     matrix_unlooped[index,4] = plane_num[i]
#                     matrix_unlooped[index,5] = theta_corr[num]
#                     matrix_unlooped[index,6] = radius[num]
                
#             except: 
#                 pass
#             bar.next()
        
#         print('G-out')
#         dict_shapes = addShapes2Dict (shapes = [kspl_vSurf, kspl_CLnew, kspl_atr, kspl_vent], dict_shapes = dict_shapes, radius = [[],[],[],[]], print_txt = False)
        
#         df_unlooped = pd.DataFrame(matrix_unlooped, columns=['x','y','z','taken','z_plane','theta','radius','cj_thickness', 'myoc_intBall'])
#         df_unlooped = df_unlooped[df_unlooped['taken']==1]
#         df_unlooped_f = df_unlooped.drop(['x', 'y','z'], axis=1)
#         df_unlooped_f.astype({'taken': 'bool','z_plane':'float16','theta':'float16','radius':'float16','cj_thickness':'float16','myoc_intBall':'float16' }).dtypes
#         print('\n')
#         saveDF(filename = filename, df2save = df_unlooped_f, df_name = 'df_'+names[ksp_num]+param_name,
#                         dir2save = dir_results)
#         dfs_unlooped.append(df_unlooped_f)
#         planes_all.append(planes)
#         spheres_all.append(spheres)
        
#     bar.finish() 
        
#     if plotshow:
#         vp = Plotter(N=2, axes = 4)
#         vp.show(mesh, kspl_vSurf, kspls_HR[0], planes_all[0], spheres_all[0], Text2D(texts[0], c="k", font= font), at=0)
#         vp.show(mesh, kspl_vSurf, kspls_HR[1], planes_all[1], spheres_all[1], Text2D(texts[1], c="k", font= font), at=1, interactive = True)
        
#     return df_unlooped, kspl_vSurf, kspls_HR, arr_all, arr_valve, dict_pts, dict_shapes, dict_kspl

#%% func - unloopChambers_OLD2
# def unloopChambers_OLD2(filename, mesh, kspl_CL, kspl_ext, no_planes, pl_CLRibbon, param, param_name, df_AtrVent, dict_kspl, 
#                     dict_shapes, dict_pts, dict_planes, dir_results, plotshow = False, tol = 0.05):
#     """
#     Function to unloop the heart and get two dataframes with the data for the unlooped hearts

#     Parameters
#     ----------
#     filename : str
#         Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
#     mesh : mesh
#         Color coded mesh that will be unlooped
#     kspl_CL : Kspline
#         Centreline (vedo KSpline)
#     kspl_ext : Kspline
#         Extended centreline (vedo KSpline).
#     no_planes : int
#         Number of planes that will be used to get transverse sections of each chamber (heart = chamber*2)
#     pl_CLRibbon :  Plane
#         Plane used to extend dorso-ventrally the centreline
#     param : list of arrays of floats
#         List of numpy arrays with distance measurements that want to be unlooped
#     param_name : str
#         Name of distance parameter(s) to add to saving file
#     df_AtrVent : numpy array of objects 
#         Array with points classification of atrium and ventricle
#     dict_kspl : dictionary
#         Initialised dictionary with ksplines information
#     dict_shapes : dictionary
#         Initialised dictionary with shapes information
#     dict_pts : dictionary
#         Initialised dictionary with points information
#     dict_planes : dictionary
#         Initialised dictionary with planes information
#     dir_results : path
#         Path to the folder where the results are saved.
#     plotshow : boolean, optional
#         True if you want to see the resulting mesh in a plot, else False. The default is False.
#     tol : float, optional
#         Tolerance defined to get points in plane. The default is 0.05.

#     Returns
#     -------
#     df_unlooped : dataframe
#         Dataframe containing unlooped heart data.
#     kspl_vSurf : Kspline
#         KSpline in the surface of the heart indicating the 0 degree direction.
#     kspls_HR : List of KSplines
#         List of ksplines for atrium and ventricle with high resolution.
#     arr_all : list of arrows 
#         List of arrows from the centreline to the ventral centre point of plane
#     arr_valve : list of arrows 
#         List of arrows from the centreline to the valve
#     dict_pts : dictionary
#         Resulting dictionary with points information updated
#     dict_shapes : dictionary
#         Resulting dictionary with shapes information updated
#     dict_kspl : dictionary
#         Resulting dictionary with ksplines information updated

#     """

#     print('\n\n- Unlooping the heart chambers...')
#     dfs_unlooped = []
#     planes_all = []
#     spheres_all = []
#     if plotshow:
#         plotevery = no_planes // 3
#         print('- Plotting every X number of planes:', plotevery)

#     # # - Get unitary normal of plane to create CL_ribbon
#     # pl_normCLRibbon = unit_vector(pl_CLRibbon['pl_normal'])
    
#     # # Increase the resolution of the extended centreline and interpolate to unify sampling
#     # xd = np.diff(kspl_ext.points()[:,0])
#     # yd = np.diff(kspl_ext.points()[:,1])
#     # zd = np.diff(kspl_ext.points()[:,2])
#     # dist = np.sqrt(xd**2+yd**2+zd**2)
#     # u = np.cumsum(dist)
#     # u = np.hstack([[0],u])
#     # t = np.linspace(0, u[-1], 1000)
#     # resamp_pts = interpn((u,), kspl_ext.points(), t)
#     # kspl_ext = KSpline(resamp_pts, res = 1000).lw(5).color('deeppink').legend('kspl_extHR')
    
#     # #Find the points that intersect with the ribbon
#     # pts_int = []
#     # for num in range(len(kspl_ext.points())):
#     #     try: 
#     #         cl_pt_test = kspl_ext.points()[num]
#     #         pt_int = mesh.intersectWithLine(cl_pt_test, cl_pt_test+60*pl_normCLRibbon)
#     #         rad_pts = [np.linalg.norm(x- cl_pt_test) for x in pt_int]
#     #         # if len(rad_pts)>1:
#     #             # print(rad_pts)
#     #         ind_pt = np.where(rad_pts == max(rad_pts))[0][0]
#     #         # print(ind_pt)
#     #         pts_int.append(pt_int[ind_pt])

#     #     except: 
#     #         # print('exc')
#     #         if num > 750:
#     #             dist_pts = kspl_ext.points()[num] - kspl_ext.points()[num-1]
#     #             try: 
#     #                 pt_int = pts_int[-1]+dist_pts
#     #                 pts_int.append(pt_int)
#     #             except:
#     #                 # print('pass')
#     #                 pass
#     #         else: 
#     #             pass
    
#     # # KSpline on surface
#     # kspl_vSurf = KSpline(pts_int).color('black').lw(4).legend('kspl_VSurfaceIntMyoc')
    
#     # KSpline on surface
#     kspl_vSurf = getExtCLonSurf(filename, mesh, kspl_ext, pl_CLRibbon)
    
#     # inv = [True, False]
#     # # Create kspline_extended with higher resolution, but cut it using the inflow and outflow  planes defined to 
#     # # cut inf anf outf tracts of hearts to extract centreline
#     # try: 
#     #     add_pts = 50
#     #     kspl_CLnew = kspl_ext.clone()
#     #     # Cut with inflow plane
#     #     n_points_In = -10
#     #     for invert in inv:
#     #         kspl_test = kspl_ext.clone().cutWithMesh(Plane(pos=dict_planes['pl2CutMesh_inflow']['pl_centre'], normal=dict_planes['pl2CutMesh_inflow']['pl_normal'], sx=300), invert=invert)
#     #         # print('kspl_test In',kspl_test.NPoints())
#     #         if kspl_test.NPoints() > n_points_In: 
#     #             # print('in-In')
#     #             inv_fin_In = invert
#     #             kspl_CLnew_cutIn = kspl_ext.clone().cutWithMesh(Plane(pos=dict_planes['pl2CutMesh_inflow']['pl_centre'], normal=dict_planes['pl2CutMesh_inflow']['pl_normal'], sx=300), invert=inv_fin_In).color('darkorange')
#     #             n_points_In = kspl_CLnew_cutIn.NPoints() 
#     #     _, num_pt_inf = findClosestPtGuess(kspl_CLnew_cutIn.points(), kspl_ext.points(), index_guess = -1)
#     # except: 
#     #     num_pt_inf = kspl_ext.NPoints()-1
    
#     # try: 
#     #     n_points_Out = -10
#     #     for invert in inv:
#     #         kspl_test = kspl_ext.clone().cutWithMesh(Plane(pos=dict_planes['pl2CutMesh_outflow']['pl_centre'], normal=dict_planes['pl2CutMesh_outflow']['pl_normal'], sx=300), invert=invert)
#     #         # print('kspl_test Out',kspl_test.NPoints())
#     #         if kspl_test.NPoints() > n_points_Out: 
#     #             # print('in-Out')
#     #             inv_fin_Out = invert
#     #             kspl_CLnew_cutOut = kspl_ext.clone().cutWithMesh(Plane(pos=dict_planes['pl2CutMesh_outflow']['pl_centre'], normal=dict_planes['pl2CutMesh_outflow']['pl_normal'], sx=300), invert=inv_fin_Out).color('lime')
#     #             n_points_Out = kspl_CLnew_cutOut.NPoints() 
#     #     # _, num_pt_outf =  findClosestPt(kspl_CLnew_cutOut.points()[-1], kspl_ext.points())#, index_guess = 0)
#     #     _, num_pt_outf =  findClosestPtGuess(kspl_CLnew_cutOut.points(), kspl_ext.points(), index_guess = 0)
#     # except: 
#     #     num_pt_outf = 0
#     # # print('Initial def:',num_pt_outf, num_pt_inf)
    
#     # # Define starting point of new kspline
#     # if (num_pt_outf-add_pts) < 0:
#     #     ind_outf = 0
#     # else: 
#     #     ind_outf = num_pt_outf - add_pts
        
#     # # Define ending point of new kspline
#     # if (num_pt_inf+add_pts) > len(kspl_ext.points()):
#     #     ind_inf = len(kspl_ext.points())
#     # else: 
#     #     ind_inf = num_pt_inf+add_pts
#     # # print('Final def:',ind_outf, ind_inf)

#     # Create new extended and cut kspline with higher resolution
#     # kspl_CLnew = KSpline(kspl_ext.points()[ind_outf:ind_inf], res = 600).lw(5).color('deeppink').legend(kspl_CL._legend+'_HighRes')

#     # Create new extended and cut kspline with higher resolution
#     kspl_CLnew = getExtCLHighRes(filename, mesh, kspl_ext, kspl_CL, dict_planes) 
    
#     # # Find position within new kspline cut by disc used to cut chambers (num_pt new)
#     # ksplCL_cut = kspl_CLnew.clone().cutWithMesh(Plane(pos=dict_shapes['cyl2CutChambers_o']['cyl_centre'], normal=dict_shapes['cyl2CutChambers_o']['cyl_axis'], sx=300), invert=True)
#     # # Find point of new kspline closer to last point of kspline cut
#     # _, num_pt = findClosestPtGuess(ksplCL_cut.points(), kspl_CLnew.points(),300)
#     # # print('num_pt',num_pt)
    
#     # # Add pt to dict
#     # sph_cut = Sphere(pos = kspl_CLnew.points()[num_pt], r=4, c='gold').legend('sph_ChamberCutCLHighRes')
#     # dict_pts = addPoints2Dict(spheres = [sph_cut], info = [''], dict_pts = dict_pts)
#     # dict_pts['numPt_CLChamberCutHighRes'] = num_pt
    
#     # # Create kspline for each chamber
#     # n_vent_pts = 0
#     # res_v = 595
#     # while n_vent_pts < 610:
#     #     res_v +=1
#     #     kspl_vent_pts = kspl_CLnew.points()[0:num_pt]
#     #     kspl_vent = KSpline(points = kspl_vent_pts[::-1], res = res_v).color('tomato').lw(8).legend('kspl_vent')
#     #     # kspl_vent = KSpline(points = kspl_CLnew.points()[0:num_pt], res = res_v).color('tomato').lw(8).legend('kspl_vent')
#     #     n_vent_pts = kspl_vent.NPoints()
        
#     # n_atr_pts = 0
#     # res_a = 595
#     # while n_atr_pts < 610:
#     #     res_a +=1
#     #     kspl_atr = KSpline(points = kspl_CLnew.points()[num_pt:], res = res_a).color('goldenrod').lw(8).legend('kspl_atr')
#     #     n_atr_pts = kspl_atr.NPoints()

#     kspl_vent, kspl_atr, dict_pts, sph_cut =  ksplChamberCut(mesh, kspl_CLnew, dict_shapes, dict_pts)
#     num_pt = dict_pts['numPt_CLChamberCutHighRes']
    
#     #Add all ksplines created to dictionary
#     dict_kspl = addKSplines2Dict(kspls = [kspl_ext, kspl_vSurf, kspl_atr, kspl_vent], info = ['','','',''], dict_kspl = dict_kspl)
    
#     vp = Plotter(N=3, axes=4)
#     vp.show(kspl_CLnew, kspl_vSurf, kspl_ext, mesh, sph_cut, at = 0)
#     # vp.show(ksplCL_cut, arr_vectPlCut, kspl_vSurf, mesh, sph_cut, at = 0)
#     vp.show(mesh, kspl_vent, sph_cut, at=1)
#     vp.show(mesh, kspl_atr, sph_cut, at=2, interactive = True)
    
#     # NOW UNLOOP EACH CHAMBER
#     kspls_HR = [kspl_atr, kspl_vent]
#     texts = ['- Unlooping atrium', '- Unlooping ventricle']
#     names= ['unloopAtr', 'unloopVent']
#     arr_all = []
#     arr_valve = []
#     for ksp_num, kspl in enumerate(kspls_HR):
#         # Create matrix with all data
#         #0:x, 1:y, 2:z, 3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: parameters
#         matrix_unlooped = np.zeros((len(mesh.points()),9))
#         matrix_unlooped[:,0:3] = mesh.points()
#         matrix_unlooped[:,7] = param[0]
#         matrix_unlooped[:,8] = param[1]
    
#         # Get normals and centres of planes to cut heart (number of planes given as input+2)
#         pl_normals, pl_centres = getPlaneNormals(no_planes = no_planes, spline_pts = kspl.points())
#         # Give a number between 1-2 to each atrial plane, and between 0-1 to each ventricular plane
#         plane_num = np.linspace(2-ksp_num,1-ksp_num,len(pl_normals))
#         if ksp_num == 1:
#             plane_num = plane_num[::-1]
#         else: 
#             index_guess = len(kspl_CLnew.points())-1
#         planes = []
#         spheres = []
#         bar = Bar(texts[ksp_num], max=len(pl_normals), suffix = suffix, check_tty=False, hide_cursor=False)
#         # Iterate through each plane
#         for i, normal, centre in zip(count(), pl_normals, pl_centres):
#             if ksp_num == 1:
#                 normal = -normal
#             # A. Get cut plane info
#             # - Info Plane (Plane orientation/normal representation)
#             arr_vectPlCut = Arrow(centre, centre+normal*20, s = 0.1, c='orange')
#             # Cut cl with plane
#             ksplCL_cut = kspl_CLnew.clone().cutWithMesh(Plane(pos=centre, normal=normal, sx=300), invert=True).lw(5).color('tomato')
#             # Find point of centreline closer to last point of kspline cut
#             # print(len(ksplCL_cut.points()), len(kspl_CLnew.points()), index_guess)
#             if len(ksplCL_cut.points()) == 0:
#                 pts_o = [kspl_CLnew.points()[-1]] #???
#             else: 
#                 pts_o = ksplCL_cut.points(); 
#                 # print('len: ',len(pts_o))
            
#             if ksp_num == 1: 
#                 index_guess= len(ksplCL_cut.points())
#             # print('index_guess:',index_guess)
#             pt_out, pt_num = findClosestPtGuess(pts_o, kspl_CLnew.points(), index_guess)
#             index_guess = pt_num
            
#             if pt_num < num_pt:
#                 chamber = 'ventricle'
#             else: 
#                 chamber = 'atrium'
#             # print('pt_num:',pt_num, 'num_pt:',num_pt, chamber, pt_out)
            
#             try: 
#                 if ksp_num == 0 and i == 0:
#                     tol2use = 0.01
#                 else:
#                     tol2use = tol
#                 # print('tol: ', tol2use)
                
#                 # B. Get vector that defines 0 deg angle
#                 pts_zero, _ = getPointsAtPlane(points = kspl_vSurf.points(), pl_normal = normal,
#                                                 pl_centre = centre, tol=0.25)
#                 if ksp_num == 0 and i == 0:
#                     pt4dist = centre
                    
#                 min_dist = 100000
#                 for ik, pt in enumerate(pts_zero):
#                     dist = findDist(pt,pt4dist)
#                     if dist < min_dist:
#                         ind = ik; min_dist = dist
#                 pt_zero = pts_zero[ind]
#                 pt4dist = pt_zero
                
#                 # Vector from centre to cl_surface point being cut by plane
#                 v_zero = unit_vector(pt_zero - centre)
        
#                 # C. Get points of mesh at plane
#                 d_points = np.absolute(np.dot(np.subtract(matrix_unlooped[:,0:3],np.asarray(centre)),np.asarray(normal)))
#                 # Find the indexes of the points that have not been yet taken, are at the plane and are in the 
#                 # chamber being analysed
#                 index_ptsAtPlane = np.where((d_points <= tol2use) & (matrix_unlooped[:,3] == 0) & (df_AtrVent == chamber))
#                 # print(d_points.min(), d_points.max(),'-lenptsatplane:',len(index_ptsAtPlane[0]))
                
#                 # Define new matrix just with the points on plane
#                 new_matrix = matrix_unlooped[index_ptsAtPlane,:][0]
#                 # - Get points of mesh that are on plane, centered on centreline point
#                 ptsC = np.subtract(new_matrix[:,0:3],np.asarray(centre))
#                 # - Get the radius of those points
#                 radius = [np.linalg.norm(x) for x in ptsC]
    
#                 # D. Find direction of point with respect to plane that includes central point, vC and the normal of the cutting plane
#                 # Vector normal to plane normal and v_zero (vector defined from centre of plane to cut-pt in centreline surface)
#                 normal_divLR = np.cross(normal, v_zero)
#                 # Define using these vectors if the points are all lying in the same side or not (1, -1)
#                 lORr = np.sign(np.dot(ptsC, np.asarray(normal_divLR)))
    
#                 # E. Get angle of points in that plane using v_zero
#                 av = np.dot(ptsC,v_zero)
#                 cosTheta = np.divide(av, radius) # vectors of magnitude one
#                 theta = np.arccos(cosTheta)*180/np.pi
#                 theta_corr = np.multiply(lORr, theta)
    
#                 pl_cut = Plane(pos = centre, normal = normal, sx = 300).color('light steel blue').alpha(0.5)
#                 sph_C = Sphere(centre, r=2, c='red')
#                 arr_vectXYZ = Arrow(centre, pt_zero, s = 0.1)
#                 arr_valv = Arrow(centre, normal_divLR-centre, s = 0.1)
#                 arr_all.append(arr_vectXYZ)
#                 arr_valve.append(arr_valv)
#                 spheres_cut = []
#                 if plotshow:
#                     # Create stuff to plot
#                     sph_C = Sphere(centre, r=2, c='red')
#                     sph_VXYZ = Sphere(pt_zero, r=2, c='green')
#                     arr_vectXYZ = Arrow(centre, pt_zero, s = 0.1)
#                     if i % plotevery == 5:
#                         sphL = []; sphR = []
#                         for num, pt in enumerate(ptsC):
#                             if num % 50 == 0:
#                                 if lORr[num] == 1:
#                                     sphL.append(Sphere(pt+centre, r=2, c='blueviolet'))
#                                 else:
#                                     sphR.append(Sphere(pt+centre, r=2, c='gold'))
#                         settings.legendSize = .3
#                         text = filename+"\n\n >> Unlooping the heart"
#                         txt = Text2D(text, c=c, font=font)
#                         vp = Plotter(N=1, axes=4)
#                         # vp.show(mesh,sphL, sphR, kspl_vSurf, kspl, kspl_ext, arr_vectXYZ, arr_vectPlCut, sph_C, sph_VXYZ, pl_cut, txt, at=0, azimuth = 0, interactive=1)
#                         vp.show(mesh,sphL, sphR, kspl_vSurf, kspl, kspl_ext, ksplCL_cut, arr_vectXYZ, arr_vectPlCut, sph_C, sph_VXYZ, pl_cut, txt, at=0, azimuth = 0, interactive=1)
#                         spheres_cut.append([sphL,sphR])
#                 if i % 30 == 0:
#                     planes.append(pl_cut)
#                     spheres.append(sph_C)
#                 # - Save all obtained values in matrix_unlooped
#                 for num, index in enumerate(index_ptsAtPlane[0]):
#                     #3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: param
#                     matrix_unlooped[index,3] = 1
#                     matrix_unlooped[index,4] = plane_num[i]
#                     matrix_unlooped[index,5] = theta_corr[num]
#                     matrix_unlooped[index,6] = radius[num]
                
#             except: 
#                 pass
#             bar.next()

#         dict_shapes = addShapes2Dict (shapes = [kspl_vSurf, kspl_CLnew, kspl_atr, kspl_vent], dict_shapes = dict_shapes, radius = [[],[],[],[]], print_txt = False)
        
#         df_unlooped = pd.DataFrame(matrix_unlooped, columns=['x','y','z','taken','z_plane','theta','radius','cj_thickness', 'myoc_intBall'])
#         df_unlooped = df_unlooped[df_unlooped['taken']==1]
#         df_unlooped_f = df_unlooped.drop(['x', 'y','z'], axis=1)
#         df_unlooped_f.astype({'taken': 'bool','z_plane':'float16','theta':'float16','radius':'float16','cj_thickness':'float16','myoc_intBall':'float16' }).dtypes
#         print('\n')
#         saveDF(filename = filename, df2save = df_unlooped_f, df_name = 'df_'+names[ksp_num]+param_name,
#                         dir2save = dir_results)
#         dfs_unlooped.append(df_unlooped_f)
#         planes_all.append(planes)
#         spheres_all.append(spheres)
        
#     bar.finish() 
        
#     if plotshow:
#         vp = Plotter(N=2, axes = 4)
#         vp.show(mesh, kspl_vSurf, kspls_HR[0], planes_all[0], spheres_all[0], Text2D(texts[0], c="k", font= font), at=0)
#         vp.show(mesh, kspl_vSurf, kspls_HR[1], planes_all[1], spheres_all[1], Text2D(texts[1], c="k", font= font), at=1, interactive=True)
        
#     #return df_unlooped, kspl_vSurf, kspls_HR, arr_all, arr_valve, dict_pts, dict_shapes, 
#     return df_unlooped, kspl_vSurf, kspls_HR, arr_all, spheres_cut, dict_pts, dict_shapes, dict_kspl

#%% func - getChambersOrientationBU
# def getChambersOrientationBU(filename, file_num, num_pt, kspl_CL2use, distFromCl, myoc_meshes, linLine, dict_pts, dict_kspl, df_res):
#     """
#     Function to get chambers orientation

#     Parameters
#     ----------
#     filename : str
#         Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
#     file_num : int
#         Index number of the selected heart being processed.
#     num_pt : int
#         Index of the centreline point closer to the plane that cuts the meshes into the two chambers.
#     kspl_CL2use : Kspline
#         Centreline (vedo KSpline)
#     distFromCl : int
#         Number of points apart from the AV Canal (num_pt) in which to define centre of chambe.
#     myoc_meshes : list of meshes
#         Myocardial meshes [m_myoc, m_atrMyoc, m_ventMyoc]
#     linLine : Line
#         Line defining the linear length of the heart (inflow-outflow)
#     dict_pts : dictionary
#         Initialised dictionary with points information
#     dict_kspl :  dictionary
#         Initialised dictionary with kspline information
#     df_res : dataframe
#         Dataframe with the measured information of the heart being processed.

#     Returns
#     -------
#     sph_orient : TYPE
#         DESCRIPTION.
#     lines_orient : TYPE
#         DESCRIPTION.
#     dict_pts : dictionary
#         Resulting dictionary with points information updated
#     dict_kspl : dictionary
#         Resulting dictionary with ksplines information updated
#     df_res : dataframe
#         Dataframe with the updated measured information of the heart being processed.

#     """

#     print('- Measuring chamber orientations')
#     m_myoc, m_atrMyoc, m_ventMyoc = myoc_meshes
#     atr_pt = num_pt + distFromCl
#     vent_pt = num_pt - distFromCl

#     # Create spheres
#     sph_atr = Sphere(pos = kspl_CL2use.points()[atr_pt], r=4, c='navy').legend('sph_AtrCentre')
#     sph_vent = Sphere(pos = kspl_CL2use.points()[vent_pt], r=4, c='red').legend('sph_VentCentre')
#     sph_inf = Sphere(pos = kspl_CL2use.points()[-1], r=4, c='dodgerblue').legend('sph_OutflowCentre')
#     sph_outf = Sphere(pos = kspl_CL2use.points()[0], r=4, c='tomato').legend('sph_InflowCentre')
#     sph_valve = Sphere(pos = kspl_CL2use.points()[num_pt], r=4, c='darkorange').legend('sph_AVCCentre')

#     sph_orient = [sph_atr, sph_vent, sph_outf, sph_inf, sph_valve]
#     dict_pts = addPoints2Dict(spheres = sph_orient, info = ['Orient','Orient','Orient','Orient','Orient'], dict_pts = dict_pts)

#     # Find orientation in which images where taken or if it is a spaw mutant
#     spaw_analysis = False
#     if 'spaw' in df_res.loc[file_num,'spAnalysis']:
#         spaw_analysis = True #ask4input('You are processing a heart that came from an incross of spaw heterozygous.\n  Please, select the way this heart is looping to continue processing: \n\t[0]: right-left\n\t[1]: dorso-ventral >>>: ', bool)
        
#     dORv = filename[9:10]
#     # Heart Orientation
#     if dORv == 'D' or 'CJ' in filename:
#         linLineX = linLine.clone().projectOnPlane('z').c(linLine.color()).x(0).legend('linLine(ProjX)')
#         pts_heart = orientVectors(linLineX)
#         # print('pts_heart', pts_heart)
#         ang_heart = findAngleBtwVectorsZ(pts_heart, np.array([[0,0,0],[0,1,0]]))
#         azimuth = -135; elevation = 0
#     elif spaw_analysis == True:
#         linLineX = linLine.clone().projectOnPlane('y').c(linLine.color()).x(0).legend('linLine(ProjX)')
#         pts_heart = orientVectors(linLineX)
#         # print('pts_heart', pts_heart)
#         ang_heart = findAngleBtwVectorsZ(pts_heart, np.array([[0,0,0],[0,1,0]]))
#         azimuth = -135; elevation = 0
#     elif dORv == 'V' or spaw_analysis == False: 
#         linLineX = linLine.clone().projectOnPlane('x').c(linLine.color()).x(0).legend('linLine(ProjX)')
#         pts_heart = orientVectors(linLineX)
#         # print('pts_heart', pts_heart)
#         ang_heart = findAngleBtwVectorsZ(pts_heart, np.array([[0,1,0],[0,0,0]]))
#         azimuth = -45; elevation = 0

#     # Atrial orientation
#     orient_atr = Line(sph_inf.pos(), sph_atr.pos(), c="steelblue", lw=3).legend('lin_OrientAtr')
#     # _, y_atrInf, z_atrInf = sph_inf.pos()
#     # _, y_atrCent, z_atrCent = sph_atr.pos()

#     if dORv == 'D' or 'CJ' in filename:
#         orient_atrX = orient_atr.clone().projectOnPlane('z').c('steelblue').x(0).legend('lin_OrientAtr(ProjX)')
#     elif spaw_analysis == True:
#         orient_atrX = orient_atr.clone().projectOnPlane('y').c('steelblue').x(0).legend('lin_OrientAtr(ProjX)')
#     elif dORv == 'V' or spaw_analysis == False:
#         orient_atrX = orient_atr.clone().projectOnPlane('x').c('steelblue').x(0).legend('lin_OrientAtr(ProjX)')

#     pts_atr = orientVectors(orient_atrX, ref_pt = pts_heart[0])
#     ang_atr = findAngleBtwVectorsZ(pts_atr, pts_heart)
#     atr_txt = 'Atrial orientation with respect to heart (deg): '+ format(ang_atr,'.1f')
#     print('\t>> '+atr_txt)

#     # Ventricular orientation
#     orient_vent = Line(sph_vent.pos(), sph_outf.pos(), c="hotpink", lw=3).legend('lin_OrientVent')
#     # _, y_ventCent, z_ventCent = sph_vent.pos()
#     # _, y_ventOutf, z_ventOutf = sph_outf.pos()
    
#     if dORv == 'D' or 'CJ' in filename:
#         orient_ventX = orient_vent.clone().projectOnPlane('z').c('hotpink').x(0).legend('lin_OrientVent(ProjX)')
#     elif spaw_analysis == True:
#         orient_ventX = orient_vent.clone().projectOnPlane('y').c('hotpink').x(0).legend('lin_OrientVent(ProjX)')
#     elif dORv == 'V' or spaw_analysis == False:
#         orient_ventX = orient_vent.clone().projectOnPlane('x').c('hotpink').x(0).legend('lin_OrientVent(ProjX)')
        
#     pts_vent = orientVectors(orient_ventX, ref_pt = pts_heart[1])
#     ang_vent = findAngleBtwVectorsZ(pts_vent, pts_heart)
#     vent_txt = 'Ventricular orientation with respect to heart (deg): '+ format(ang_vent,'.1f')
#     print('\t>> '+vent_txt)

#     # Angle between chambers
#     ang_chs = findAngleBtwVectorsZ(pts_atr, pts_vent)
#     chs_txt = 'Angle between chambers (deg): '+ format(ang_chs,'.1f')
#     print('\t>> '+chs_txt)

#     lines_orient = [orient_atr, orient_vent, orient_atrX, orient_ventX, linLineX]
#     dict_kspl = addKSplines2Dict(kspls = lines_orient, info = '', dict_kspl = dict_kspl)
#     df_res = addOrientationAngles2df(df_res = df_res, file_num = file_num, angles = [ang_heart, ang_atr, ang_vent, ang_chs], names = ['ang_Heart', 'ang_Atr', 'ang_Vent', 'ang_BtwChambers'])

#     text = filename+"\n\n >> Measuring chamber orientation\n - "+atr_txt+"\n - "+vent_txt+"\n - "+chs_txt
#     txt = Text2D(text, c="k", font= font)
#     settings.legendSize = .15
#     vp = Plotter(N=1, axes = 8)
#     vp.show(m_myoc.alpha(0.01), m_atrMyoc.alpha(0.01), m_ventMyoc.alpha(0.01), sph_orient, kspl_CL2use, orient_atr, orient_vent, linLine,
#             orient_atrX.alpha(1), orient_ventX.alpha(1), linLineX.alpha(1), txt, azimuth = azimuth, elevation = elevation, zoom = 0.8, at=0, interactive=True) #
#     # vp.show()

#     return sph_orient, lines_orient, dict_pts, dict_kspl, df_res

#%% THE END