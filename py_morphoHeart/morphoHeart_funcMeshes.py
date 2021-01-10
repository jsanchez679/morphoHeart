# -*- coding: utf-8 -*-
"""
morphoHeart_funcMeshes
Functions included: 
    - alert(sound, times)
    - 
    
Version: September 30, 2020    
@author: Juliana Sanchez-Posada

"""

#%% Importing python packages
import numpy as np
import math
import os

from itertools import count

from skimage import measure, io
import numpy as np
from scipy.interpolate import splprep, splev
import pandas as pd

from datetime import datetime
from time import perf_counter
from progress.bar import Bar
suffix = '%(index)d/%(max)d - %(elapsed)ds'

c="k"
font= 'CallingCode'

from vtkplotter import *
from vtkplotter import embedWindow
from time import perf_counter
embedWindow(False)

import json
from json import JSONEncoder

#%% Importing morphoHeart packages
from morphoHeart_funcBasics import alert, ask4input
from morphoHeart_funcContours import save_s3s

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
    
    names_all = ['myoc','myoc_ext','myoc_int','myoc_atr','myoc_vent',
                  'endo','endo_ext','endo_int','endo_atr','endo_vent', 
                  'cj','cj_out','cj_in','cj_atr','cj_vent']
    
    legend_all = ['Myocardium','Ext.Myoc', 'Int.Myoc','Atrium(Myoc)','Ventricle(Myoc)',
                  'Endocardium', 'Ext.Endo', 'Int.Endo','Atrium(Endo)','Ventricle(Endo)',
                  'CardiacJelly','Ext.CJ','Int.CJ','Atrium(CJ)','Ventricle(CJ)']
    
    print('- Loading meshes...')
    meshes_out = []
    
    for i, name in enumerate(meshes_names):    
        
        index = names_all.index(name)
        mesh_title = filename+"_"+name+"."+extension
        mesh_dir = os.path.join(dir_stl, mesh_title)
        mesh_out = load(mesh_dir)
        print("\t>> Mesh loaded - "+mesh_title+"!")
        
        mesh_colour = dict_colour[name]['colour']
        mesh_out.alpha(alpha[i]).legend(legend_all[index]).wireframe().color(mesh_colour)
    
        meshes_out.append(mesh_out)
        
    alert("wohoo",1)
    print('\n')
    
    return meshes_out

#%% func - openThicknessMeshes 
def openThicknessMeshes(filename, meshes_names, extension, dir_stl, dir_txtNnpy):
    
    names_all = ['myoc_intBall', 'myoc_extBall', 'myoc_thickness','endo_thickness','cj_thickness']
    legend_all = ['Int.Myoc Ball.','Ext.Myoc Ball.','Myoc.Thickness','Endo.Thickness','CJ.Thickness']
    bar_names_all = ['Int.Myoc\nBalloning\n[um]','Ext.Myoc\nBalloning\n[um]','Myoc.Thickness\n[um]','Endo.Thickness\n[um]','CJ.Thickness\n[um]']
    alpha_all = [1,1,0.1,0.1,0.1]
    
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
        
        
    alert("wohoo",1)
    print('\n')
    
    return meshes_out, colour_arrays

#%% - DICTIONARIES 
#%% func - fillNsaveObjDict
def fillNsaveObjDict(filename, dicts, names, dir2save):
    
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

#%% func - splitDicts
def splitDicts(dict_obj):
    
    dicts_all = []
    for i, name in enumerate(list(dict_obj.keys())):
        ext_dict = dict_obj[name]
        dicts_all.append(ext_dict)
        
    return dicts_all

#%% func - saveDict
def saveDict(filename, dict2save, name, dir2save):
    
    jsonDict_name = filename+"_"+name+".json"
    json2save_dir = os.path.join(dir2save,jsonDict_name)
    
    with open(json2save_dir, "w") as write_file:
        json.dump(dict2save, write_file, cls=NumpyArrayEncoder)
    print("\t>> Dictionary saved correctly!\n\t> File: "+jsonDict_name); 
    alert("countdown",1)

#%% func - decodeDict
def decodeDict (dict2classify, info):
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

#%% func - addPlane2Dict
def addPlane2Dict (plane, pl_centre, pl_normal, info, dict_planes, print_txt = True):
    
    if info =='':
        specific_plane = dict_planes[plane.legend()]= dict()
    else: 
        specific_plane = dict_planes[plane.legend()+'-'+info]= dict()
        
    specific_plane['pl_normal'] = pl_normal
    specific_plane['pl_centre'] = pl_centre
    specific_plane['color'] = plane.color()
    
    if print_txt:
        print('\t>> Plane has been added to dict_planes!')
    
    return dict_planes

#%% func - addPlanes2Dict
def addPlanes2Dict (planes, pls_centre, pls_normal, info, dict_planes, print_txt = True):
    
    for i, plane, inf in zip(count(), planes, info):
        if inf =='':
            specific_plane = dict_planes[plane.legend()]= dict()
        else: 
            specific_plane = dict_planes[plane.legend()+'-'+inf]= dict()
            
        specific_plane['pl_normal'] = pls_normal[i]
        specific_plane['pl_centre'] = pls_centre[i]
        specific_plane['color'] = plane.color()
    
    if print_txt:
        print('\t>> Planes have been added to dict_planes!')
    
    return dict_planes

#%% func - addPoint2Dict
def addPoint2Dict (sphere, info, dict_pts, print_txt = True):
    
    if info =='':
        specific_pt = dict_pts[sphere.legend()]= dict()
    else: 
        specific_pt = dict_pts[sphere.legend()+'-'+info]= dict()
        
    specific_pt['sph_position'] = sphere.pos()
    specific_pt['color'] = sphere.color()
    
    if print_txt:
        print('\t>> Point has been added to dict_pts!')
    
    return dict_pts

#%% func - addPoints2Dict
def addPoints2Dict (spheres, info, dict_pts, print_txt = True):
    
    for i, sph, inf in zip(count(), spheres, info):
        if inf =='':
            specific_pt = dict_pts[sph.legend()]= dict()
        else: 
            specific_pt = dict_pts[sph.legend()+'-'+inf]= dict()
            
        specific_pt['sph_position'] = sph.pos()
        specific_pt['color'] = sph.color()
    
    if print_txt:
        print('\t>> Points have been added to dict_pts!')
    
    return dict_pts

#%% func - addKSpline2Dict
def addKSpline2Dict (kspl, info, dict_kspl, print_txt = True):
    
    if info =='':
        specific_kspl = dict_kspl[kspl.legend()]= dict()
    else: 
        specific_kspl = dict_kspl[kspl.legend()+'-'+info]= dict()
        
    specific_kspl['kspl_pts'] = kspl.points()
    specific_kspl['color'] = kspl.color()
    
    if print_txt:
        print('\t>> KSpline has been added to dict_kspl!')
    
    return dict_kspl

#%% func - addKSplines2Dict
def addKSplines2Dict (kspls, info, dict_kspl, print_txt = True):
    
    for i, kspl, inf in zip(count(), kspls, info):
        #print(i)
        if inf =='':
            specific_kspl = dict_kspl[kspl.legend()]= dict()
        else: 
            specific_kspl = dict_kspl[kspl.legend()+'-'+inf]= dict()
            
        specific_kspl['kspl_pts'] = kspl.points()
        specific_kspl['color'] = kspl.color()
    
    if print_txt:
        print('\t>> KSplines have been added to dict_kspl!')
    
    return dict_kspl

#%% func - addShapes2Dict 
def addShapes2Dict(shapes, info, dict_shapes, radius, print_txt = True):
    
    for i, sp_shape, inf, radii in zip(count(), shapes, info, radius):
        if inf =='':
            specific_shape = dict_shapes[sp_shape.legend()]= dict()
        else: 
            specific_shape = dict_shapes[sp_shape.legend()+'-'+inf]= dict()
    
        specific_shape['kspl_pts'] = sp_shape.points()
        specific_shape['color'] = sp_shape.color()
        
        if len(radii) != 0:
            specific_shape['maxInscSphRad'] = radii
    
    if print_txt:
        print('\t>> Shapes have been added to dict_shapes!')

    return dict_shapes

#%% - DATAFRAMES
#%% func - addSurfArea2df
def addSurfArea2df (df_res, file_num,  meshes):
    names = ['Myoc', 'Int.Myoc', 'Ext.Myoc', 'Endo', 'Int.Endo', 'Ext.Endo', 'CJ', 'Int.CJ', 'Ext.CJ']
    df_resFilled = df_res
    
    for n, name, mesh in zip(count(), names, meshes):
        area = mesh.area()
        df_resFilled.loc[file_num,'SurfArea_'+name] = area
        
        area_print = format(area, '.1f')
        print('Name:', name, '- SurfArea: ', area_print, 'um^2')
        
    return df_resFilled

#%% func - addLayersVolume2df
def addLayersVolume2df (df_res, file_num,  meshes):
    names = ['Myoc', 'Atr.Myoc', 'Vent.Myoc','Atr.ExtMyoc', 'Vent.ExtMyoc', 
             'Endo', 'Atr.Endo', 'Vent.Endo', 'Atr.IntEndo', 'Vent.IntEndo', 
             'CJ', 'Atr.CJ', 'Vent.CJ']
    
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
            print('\t - Vol difference: ', diff_print,  'um^3')
        else:
            print('-> ', name, '\t- Volume: ', vol_print, 'um^3')
        
    return df_resFilled

#%% func - addIntExtVol2df
def addIntExtVol2df(df_res, file_num,  meshes):
    names = ['Int.Myoc', 'Ext.Myoc', 'Int.Endo', 'Ext.Endo']
    df_resFilled = df_res
    
    print('- Getting internal and external volumes of heart layers...')
    for n, name, mesh in zip(count(), names, meshes):
        volume = mesh.volume()
        df_resFilled.loc[file_num,'Vol_'+name] = volume
        
    return df_resFilled

#%% func - addLinearMeas2df
def addLinearMeas2df(df_res, file_num, lines, kspl_CL):
    
    df_resFilled = df_res
    for nl, line in enumerate(lines):
        length_l = line.length()
        name_l = line.legend()
        df_resFilled.loc[file_num, name_l] = length_l
    
    for ncl, cl in enumerate(kspl_CL):
        length_cl = cl.length()
        name_cl = cl.legend()
        df_resFilled['Length_'+name_cl] = length_cl
    
    return df_resFilled

#%% func - addOrientationAngles2df
def addOrientationAngles2df(df_res, file_num, angles, names):
    
    df_resFilled = df_res
    for i, angle, name in zip(count(), angles, names):
        df_resFilled.loc[file_num, name] = angle

    return df_resFilled
    
#%% - S3s MANIPULATION
#%% func - maskChamberS3s 
def maskChamberS3s (s3_mask, pl_normal, pl_centre, resolution):
    
    # Get dimensions of stack
    xdim, ydim, zdim = s3_mask.shape
    # Reshape stack as a vector
    s3_mask_v = s3_mask.reshape(-1)
    # Get vectors of x,y and z positions
    pix_coord_pos = np.where(s3_mask >= 0)
    # Trasform coordinate positions to um using resolution
    pix_um = np.transpose(np.asarray([pix_coord_pos[i]*resolution[i] for i in range(len(resolution))]))
    # Make normal a unit vector
    normal_unit = unit_vector(pl_normal)
    # Find all the d values of pix_um
    d_pix_um = np.dot(np.subtract(pix_um,np.array(pl_centre)),np.array(normal_unit))
    # Clear vector d_pix_um using only those that are 1 in stack
    d_pve_pix_um = s3_mask_v*d_pix_um
    # Duplicate s3_mask_v to initialise atrium and ventricle
    s3_vent_v = np.copy(s3_mask_v)
    s3_atr_v = np.copy(s3_mask_v)
    # Find all positions in d_pve_pix_um that are at either side of the plane
    pos_atr = np.where(d_pve_pix_um < 0)[0]
    pos_vent = np.where(d_pve_pix_um > 0)[0]
    # Remove the points that belong to the other chamber
    s3_vent_v[pos_atr] = 0
    s3_atr_v[pos_vent] = 0
    # Reshape vector into matrix/stack
    s3_vent = s3_vent_v.reshape((xdim, ydim, zdim))
    s3_atr = s3_atr_v.reshape((xdim, ydim, zdim))
    
    # OLD VERSION
    # s3_atr =  np.zeros_like(s3_mask, dtype='bool')
    # s3_vent =  np.zeros_like(s3_mask, dtype='bool')
    
    # #print('Masking Chambers s3 - ' + option)
    # # Get coordinates of pixels that have a value = 1
    # pix_coord_pos = np.where(s3_mask == 1)
    # # Make normal of plane unitary 
    # normal_unit = unit_vector(pl_normal)
    # normalx10 = [k*100 for k in normal_unit]
    # # Create points at either side of the plane
    # pt_atr = pl_centre+normalx10
    # pt_vent = pl_centre-normalx10
    
    # # Get plane values
    # d = normal_unit.dot(pl_centre)
    # a, b, c = normal_unit
    # #print("Normal_unit", normal_unit)
    # #print("Plane d", d)
    
    # x_atr,y_atr,z_atr = pt_atr
    # dotProd_pt_atr = dot((a,b,c), (x_atr,y_atr,z_atr))
    # #print("dotProd_pt_atr:", dotProd_pt_atr)
    # x_vent,y_vent,z_vent = pt_vent
    # dotProd_pt_vent = dot((a,b,c), (x_vent,y_vent,z_vent))
    # #print("dotProd_pt_vent:", dotProd_pt_vent)
    
    # bar = Bar('Masking Chambers s3', max = 100, suffix = suffix, check_tty=False, hide_cursor=False)
    # num_pts = round(len(pix_coord_pos[0])/100)
    # for i, x, y, z in zip(count(), pix_coord_pos[0], pix_coord_pos[1], pix_coord_pos[2]):
    #     # Pix coordinates = [x,y,z]
    #     # [x,y,z] * resolution = Coordinates in um [x_um, y_um, z_um]
    #     pt_um = [x,y,z]*np.array(resolution)
    #     x_um, y_um, z_um = pt_um
    #     # Find dot product to know on which side of the plane each pix is
    #     dotProd_pt = dot((a,b,c), (x_um, y_um, z_um))
    #     if dotProd_pt > d and dotProd_pt_atr > d or dotProd_pt < d and dotProd_pt_atr < d:
    #         s3_atr[x,y,z] = 1
    #     elif dotProd_pt > d and dotProd_pt_vent > d or dotProd_pt < d and dotProd_pt_vent < d:
    #         s3_vent[x,y,z] = 1
            
    #     if i % num_pts == 0:
    #         bar.next()
        
    # bar.finish()
    alert('wohoo',1)
            
    return s3_atr, s3_vent

#%% func - selectCutS3sOptMx
def selectCutS3sOptMx(filename, s3s2cut, m_endo, m_myoc, dict_planes, resolution, dir_txtNnpy, save):
    
    s3_ch0_int, s3_ch0_ext, s3_ch0, s3_ch1_int, s3_ch1_ext_cl, s3_ch1_cl = s3s2cut
    
    # Create new s3s to save cuts
    s3_ch0_cut = s3_ch0; s3_ch0_int_cut = s3_ch0_int; s3_ch0_ext_cut = s3_ch0_ext
    s3_ch1_cut = s3_ch1_cl; s3_ch1_int_cut = s3_ch1_int; s3_ch1_ext_cut = s3_ch1_ext_cl
    
    cut_type = ['inflow', 'outflow']
    cuts = []
    pls_normal = []
    pls_centre = []
    cuts_selected = []
    
    # Define what to cut from each layer
    myoc_cuts = []; endo_cuts = []
    
    for cut in cut_type:
        text = filename+"\n\n >> Take a closer look at the -" +cut + "- of both meshes \n\tto decide which layer to cut\n >> [0]:myoc/[1]:endo/[2]:both/[3]:none\n >> Close the window when done"
        txt = Text2D(text, c="k", font= 'CallingCode')
        vp = Plotter(N=3, axes=4)
        vp.show(m_myoc, txt, at=0, zoom=1)
        vp.show(m_endo, at=1, zoom=1)
        vp.show(m_myoc, m_endo, at=2, zoom=1, interactive=True)
        
        q_cuts = ask4input('Select the layer from which you want to cut the -'+ cut + '- tract \n  [0]:myoc/[1]:endo/[2]:both/[3]:none?: ',int)
        cuts.append(q_cuts)
        
        if q_cuts == 0 or q_cuts == 1 or q_cuts == 2:
            cuts_selected.append(cut)
            # Get plane to cut
            plane_cut, pl_cut_centre, pl_cut_normal = getPlane(filename = filename, type_cut = cut, info = '', mesh_in = m_endo, 
                                                                        mesh_out = m_myoc)
            # Reorient plane to images (s3)
            plane_im, pl_im_centre, pl_im_normal = rotatePlane2Images(pl_cut_centre, pl_cut_normal, type_cut = cut)
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
        # print('- Cutting s3 - inf&outf (Myoc)')
        s3_ch0_cut = cutInfAndOutfOptMx(s3_ch0_cut, pls_normal, pls_centre, resolution, '(Myoc)'); bar.next()
        s3_ch0_int_cut = cutInfAndOutfOptMx(s3_ch0_int_cut, pls_normal, pls_centre, resolution, '(Myoc)'); bar.next()
        s3_ch0_ext_cut = cutInfAndOutfOptMx(s3_ch0_ext_cut, pls_normal, pls_centre, resolution, '(Myoc)'); bar.next()
        bar.finish()
        
    elif len(myoc_cuts) == 1:
        index_myoc = cuts_selected.index(myoc_cuts[0])
        # Cut myocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - ' + myoc_cuts[0]+' (Myoc)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        # print('- Cutting s3 - ' + myoc_cuts[0]+' (Myoc)')
        s3_ch0_cut = cutInfOrOutfOptMx(s3_ch0_cut, pls_normal[index_myoc], pls_centre[index_myoc], resolution = resolution, 
                                                                      option = myoc_cuts[0], mesh_name = '(Myoc)'); bar.next()
        s3_ch0_int_cut = cutInfOrOutfOptMx(s3_ch0_int_cut, pls_normal[index_myoc], pls_centre[index_myoc], resolution = resolution, 
                                                                      option = myoc_cuts[0], mesh_name = '(Myoc)'); bar.next()
        s3_ch0_ext_cut = cutInfOrOutfOptMx(s3_ch0_ext_cut, pls_normal[index_myoc], pls_centre[index_myoc], resolution = resolution, 
                                                                      option = myoc_cuts[0], mesh_name = '(Myoc)'); bar.next()
        bar.finish()
    else: 
        print('- No cuts made to Myocardium!')
    
    # Cut Endocardial layers
    if len(endo_cuts) == 2:
        # Cut endocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - inf&outf (Endo)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        # print('- Cutting s3 - inf&outf (Endo)')
        s3_ch1_cut = cutInfAndOutfOptMx(s3_ch1_cut, pls_normal, pls_centre, resolution, '(Endo)'); bar.next()
        s3_ch1_int_cut = cutInfAndOutfOptMx(s3_ch1_int_cut, pls_normal, pls_centre, resolution, '(Endo)'); bar.next()
        s3_ch1_ext_cut  = cutInfAndOutfOptMx(s3_ch1_ext_cut, pls_normal, pls_centre, resolution, '(Endo)'); bar.next()
        bar.finish()
        
    elif len(endo_cuts) == 1:
        index_endo = cuts_selected.index(endo_cuts[0])
        # Cut myocardial s3_all, s3_int, s3_ext
        bar = Bar('- Cutting s3 - ' + endo_cuts[0]+' (Endo)', max = 3, suffix = suffix, check_tty=False, hide_cursor=False)
        # print('- Cutting s3 - ' + endo_cuts[0]+' (Endo)')
        s3_ch1_cut = cutInfOrOutfOptMx(s3_ch1_cut, pls_normal[index_endo], pls_centre[index_endo], resolution = resolution, 
                                                                      option = endo_cuts[0], mesh_name = '(Endo)'); bar.next()
        s3_ch1_int_cut = cutInfOrOutfOptMx(s3_ch1_int_cut, pls_normal[index_endo], pls_centre[index_endo], resolution = resolution, 
                                                                      option = endo_cuts[0], mesh_name = '(Endo)'); bar.next()
        s3_ch1_ext_cut = cutInfOrOutfOptMx(s3_ch1_ext_cut, pls_normal[index_endo], pls_centre[index_endo], resolution = resolution, 
                                                                      option = endo_cuts[0], mesh_name = '(Endo)'); bar.next()
        bar.finish()
    else: 
        print('- No cuts made to Endocardium!')
        
    alert('whistle', 1)
    # Create meshes
    # Myocardial layers
    if any(x in cuts for x in [0, 2]):
        #Create new mesh myocardial s3_all
        myoc_cut = getCutMesh(filename = filename, s3_cut = s3_ch0_cut, resolution = resolution, 
                                        mesh_original = m_myoc, layer = 'Myoc', plotshow = False)
        #Create new mesh myocardial s3_int
        myoc_cut_int = getCutMesh(filename = filename, s3_cut = s3_ch0_int_cut, resolution = resolution, 
                                        mesh_original = '', layer = 'Int.Myoc', plotshow = False)
        #Create new mesh myocardial s3_ext
        myoc_cut_ext = getCutMesh(filename = filename, s3_cut = s3_ch0_ext_cut, resolution = resolution, 
                                        mesh_original = '', layer = 'Ext.Myoc', plotshow = False)
        if save:
            save_s3s(filename = filename, s3_all = s3_ch0_cut, s3_int = s3_ch0_int_cut, s3_ext = s3_ch0_ext_cut, 
                    dir_txtNnpy = dir_txtNnpy, layer = 'ch0_cut')
    else:
        myoc_cut, myoc_cut_int, myoc_cut_ext = createAll3LayerMeshes(filename = filename, s3_all = s3_ch0, s3_in = s3_ch0_int,
                                    s3_out = s3_ch0_ext, resolution = resolution, layer = 'Myoc')
    # Endocardial layers
    if any(x in cuts for x in [1, 2]):
        #Create new mesh endocardial s3_all
        endo_cut = getCutMesh(filename = filename, s3_cut = s3_ch1_cut, resolution = resolution, 
                                                    mesh_original = m_endo, layer = 'Endo', plotshow = False)
        #Create new mesh endocardial s3_int
        endo_cut_int = getCutMesh(filename = filename, s3_cut = s3_ch1_int_cut, resolution = resolution, 
                                        mesh_original = '', layer = 'Int.Endo', plotshow = False)
        #Create new mesh endocardial s3_ext
        endo_cut_ext = getCutMesh(filename = filename, s3_cut = s3_ch1_ext_cut, resolution = resolution, 
                                        mesh_original = '', layer = 'Ext.Endo', plotshow = False)
        if save:
            save_s3s(filename = filename, s3_all = s3_ch1_cut, s3_int = s3_ch1_int_cut, s3_ext = s3_ch1_ext_cut, 
                    dir_txtNnpy = dir_txtNnpy, layer = 'ch1_cut')
    else:
        endo_cut, endo_cut_int, endo_cut_ext = createAll3LayerMeshes(filename = filename, s3_all = s3_ch1_cl, s3_in = s3_ch1_int,
                                    s3_out = s3_ch1_ext_cl, resolution = resolution, layer = 'Endo')
        
    s3s_cut = [s3_ch0_int_cut, s3_ch0_ext_cut, s3_ch0_cut, s3_ch1_int_cut, s3_ch1_ext_cut, s3_ch1_cut]
    meshes_cut = [myoc_cut, myoc_cut_int, myoc_cut_ext, endo_cut, endo_cut_int, endo_cut_ext]
    alert('jump', 1)
    
    text= filename+"\n\n >> Resulting meshes"
    txt = Text2D(text, c="k", font= 'CallingCode')
    vp = Plotter(N=6, axes=10)
    for i, mesh in enumerate(meshes_cut):
        if i == 0:
            vp.show(mesh, txt, at = i, zoom = 1.2)
        elif i > 0 and i < 5:
            vp.show(mesh, at = i, zoom = 1.2)
        else:
            vp.show(mesh, at = i, zoom = 1.2, interactive = True)
    
    return s3s_cut, meshes_cut, dict_planes

# OLD VERSION
# #%% func - selectCutS3sOpt
# def selectCutS3sOpt(filename, s3s2cut, m_endo, m_myoc, dict_planes, resolution, dir_txtNnpy, save):
    
#     s3_ch0_int, s3_ch0_ext, s3_ch0, s3_ch1_int, s3_ch1_ext_cl, s3_ch1_cl = s3s2cut
    
#     # Create new s3s to save cuts
#     s3_ch0_cut = s3_ch0; s3_ch0_int_cut = s3_ch0_int; s3_ch0_ext_cut = s3_ch0_ext
#     s3_ch1_cut = s3_ch1_cl; s3_ch1_int_cut = s3_ch1_int; s3_ch1_ext_cut = s3_ch1_ext_cl
    
#     cut_type = ['inflow', 'outflow']
#     cuts = []
#     pls_normal = []
#     pls_centre = []
#     cuts_selected = []
    
#     # Define what to cut from each layer
#     myoc_cuts = []; endo_cuts = []
    
#     for cut in cut_type:
#         text = filename+"\n\n >> Take a closer look at the -" +cut + "- of both meshes \n\tto decide which layer to cut\n >> [0]:myoc/[1]:endo/[2]:both/[3]:none\n >> Close the window when done"
#         txt = Text2D(text, c="k", font= 'CallingCode')
#         vp = Plotter(N=3, axes=4)
#         vp.show(m_myoc, txt, at=0, zoom=1)
#         vp.show(m_endo, at=1, zoom=1)
#         vp.show(m_myoc, m_endo, at=2, zoom=1, interactive=True)
        
#         q_cuts = ask4input('Select the layer from which you want to cut the -'+ cut + '- tract \n  [0]:myoc/[1]:endo/[2]:both/[3]:none?: ',int)
#         cuts.append(q_cuts)
        
#         if q_cuts == 0 or q_cuts == 1 or q_cuts == 2:
#             cuts_selected.append(cut)
#             # Get plane to cut
#             plane_cut, pl_cut_centre, pl_cut_normal = getPlane(filename = filename, type_cut = cut, info = '', mesh_in = m_endo, 
#                                                                         mesh_out = m_myoc)
#             # Reorient plane to images (s3)
#             plane_im, pl_im_centre, pl_im_normal = rotatePlane2Images(pl_cut_centre, pl_cut_normal, type_cut = cut)
#             pls_normal.append(pl_im_normal); pls_centre.append(pl_im_centre)
#             #Save planes to dict
#             dict_planes = addPlanes2Dict(planes = [plane_cut, plane_im], pls_centre = [pl_cut_centre ,pl_im_centre], 
#                                                     pls_normal = [pl_cut_normal, pl_im_normal], info = ['',''], dict_planes = dict_planes)
#             if q_cuts == 0 or q_cuts == 2:
#                 myoc_cuts.append(cut)
#             if q_cuts == 1 or q_cuts == 2:
#                 endo_cuts.append(cut)
                
#     # Get data from cutting planes
#     # d_cuts, normal_cuts, d_red, d_green = getCuttingPlanesInfoMx(pls_normal, pls_centre)
    
#     # Cut Myocardial layers
#     if len(myoc_cuts) == 2:
#         # Cut myocardial s3_all, s3_int, s3_ext
#         # s3_ch0_cut, s3_ch0_int_cut, s3_ch0_ext_cut = cutInfAndOutfOpt(s3_ch0_cut, s3_ch0_int_cut, s3_ch0_ext_cut, 
#         #                                                             d_cuts, normal_cuts, d_red, d_green, resolution, '(Myoc)')
#         s3_ch0_cut, s3_ch0_int_cut, s3_ch0_ext_cut = cutInfAndOutfOptMx(s3_ch0_cut, s3_ch0_int_cut, s3_ch0_ext_cut, 
#                                                                     pls_normal, pls_centre, resolution, '(Myoc)')
#     elif len(myoc_cuts) == 1:
#         index_myoc = cuts_selected.index(myoc_cuts[0])
#         # Cut myocardial s3_all, s3_int, s3_ext
#         # s3_ch0_cut, s3_ch0_int_cut, s3_ch0_ext_cut = cutInfOrOutfOpt (s3_ch0_cut, s3_ch0_int_cut, s3_ch0_ext_cut, 
#         #                                                               d_cut = d_cuts[index_myoc], normal_cut = normal_cuts[index_myoc], 
#         #                                                               d_red = d_red[index_myoc], d_green = d_green[index_myoc], resolution = resolution, 
#         #                                                               option = myoc_cuts[0], mesh_name = '(Myoc)')
#         s3_ch0_cut, s3_ch0_int_cut, s3_ch0_ext_cut = cutInfOrOutfOptMx(s3_ch0_cut, s3_ch0_int_cut, s3_ch0_ext_cut, 
#                                                                       pls_normal[index_myoc], pls_centre[index_myoc], resolution = resolution, 
#                                                                       option = myoc_cuts[0], mesh_name = '(Myoc)')
#     else: 
#         print('- No cuts made to Myocardium!')
    
#     # Cut Endocardial layers
#     if len(endo_cuts) == 2:
#         # Cut endocardial s3_all
#         # s3_ch1_cut, s3_ch1_int_cut, s3_ch1_ext_cut  = cutInfAndOutfOpt(s3_ch1_cut, s3_ch1_int_cut, s3_ch1_ext_cut, 
#         #                                                               d_cuts, normal_cuts, d_red, d_green, resolution, '(Endo)')
#         s3_ch1_cut, s3_ch1_int_cut, s3_ch1_ext_cut  = cutInfAndOutfOptMx(s3_ch1_cut, s3_ch1_int_cut, s3_ch1_ext_cut, 
#                                                                       pls_normal, pls_centre, resolution, '(Endo)')
        
#     elif len(endo_cuts) == 1:
#         index_endo = cuts_selected.index(endo_cuts[0])
#         # Cut myocardial s3_all
#         # s3_ch1_cut, s3_ch1_int_cut, s3_ch1_ext_cut = cutInfOrOutfOpt (s3_ch1_cut, s3_ch1_int_cut, s3_ch1_ext_cut, 
#         #                                                               d_cut = d_cuts[index_endo], normal_cut = normal_cuts[index_endo], 
#         #                                                               d_red = d_red[index_endo], d_green = d_green[index_endo], resolution = resolution, 
#         #                                                               option = endo_cuts[0], mesh_name = '(Endo)')
#         s3_ch1_cut, s3_ch1_int_cut, s3_ch1_ext_cut = cutInfOrOutfOptMx(s3_ch1_cut, s3_ch1_int_cut, s3_ch1_ext_cut, 
#                                                                       pls_normal[index_endo], pls_centre[index_endo], resolution = resolution, 
#                                                                       option = endo_cuts[0], mesh_name = '(Endo)')
#     else: 
#         print('- No cuts made to Endocardium!')
        
#     alert('whistle', 1)
#     # Create meshes
#     # Myocardial layers
#     if any(x in cuts for x in [0, 2]):
#         #Create new mesh myocardial s3_all
#         myoc_cut = getCutMesh(filename = filename, s3_cut = s3_ch0_cut, resolution = resolution, 
#                                         mesh_original = m_myoc, layer = 'Myoc', plotshow = True)
#         #Create new mesh myocardial s3_int
#         myoc_cut_int = getCutMesh(filename = filename, s3_cut = s3_ch0_int_cut, resolution = resolution, 
#                                         mesh_original = '', layer = 'Int.Myoc', plotshow = False)
#         #Create new mesh myocardial s3_ext
#         myoc_cut_ext = getCutMesh(filename = filename, s3_cut = s3_ch0_ext_cut, resolution = resolution, 
#                                         mesh_original = '', layer = 'Ext.Myoc', plotshow = False)
#         if save:
#             save_s3s(filename = filename, s3_all = s3_ch0_cut, s3_int = s3_ch0_int_cut, s3_ext = s3_ch0_ext_cut, 
#                     dir_txtNnpy = dir_txtNnpy, layer = 'ch0_cut')
#     else:
#         myoc_cut, myoc_cut_int, myoc_cut_ext = createAll3LayerMeshes(filename = filename, s3_all = s3_ch0, s3_in = s3_ch0_int,
#                                     s3_out = s3_ch0_ext, resolution = resolution, layer = 'Myoc')
#     # Endocardial layers
#     if any(x in cuts for x in [1, 2]):
#         #Create new mesh endocardial s3_all
#         endo_cut = getCutMesh(filename = filename, s3_cut = s3_ch1_cut, resolution = resolution, 
#                                                     mesh_original = m_endo, layer = 'Endo', plotshow = False)
#         #Create new mesh endocardial s3_int
#         endo_cut_int = getCutMesh(filename = filename, s3_cut = s3_ch1_int_cut, resolution = resolution, 
#                                         mesh_original = '', layer = 'Int.Endo', plotshow = False)
#         #Create new mesh endocardial s3_ext
#         endo_cut_ext = getCutMesh(filename = filename, s3_cut = s3_ch1_ext_cut, resolution = resolution, 
#                                         mesh_original = '', layer = 'Ext.Endo', plotshow = False)
#         if save:
#             save_s3s(filename = filename, s3_all = s3_ch1_cut, s3_int = s3_ch1_int_cut, s3_ext = s3_ch1_ext_cut, 
#                     dir_txtNnpy = dir_txtNnpy, layer = 'ch1_cut')
#     else:
#         endo_cut, endo_cut_int, endo_cut_ext = createAll3LayerMeshes(filename = filename, s3_all = s3_ch1_cl, s3_in = s3_ch1_int,
#                                     s3_out = s3_ch1_ext_cl, resolution = resolution, layer = 'Endo')
        
#     s3s_cut = [s3_ch0_int_cut, s3_ch0_ext_cut, s3_ch0_cut, s3_ch1_int_cut, s3_ch1_ext_cut, s3_ch1_cut]
#     meshes_cut = [myoc_cut, myoc_cut_int, myoc_cut_ext, endo_cut, endo_cut_int, endo_cut_ext]
#     alert('jump', 1)
    
#     text= filename+"\n\n >> Resulting meshes"
#     txt = Text2D(text, c="k", font= 'CallingCode')
#     vp = Plotter(N=6, axes=10)
#     for i, mesh in enumerate(meshes_cut):
#         if i == 0:
#             vp.show(mesh, txt, at = i, zoom = 1.2)
#         elif i > 0 and i < 5:
#             vp.show(mesh, at = i, zoom = 1.2)
#         else:
#             vp.show(mesh, at = i, zoom = 1.2, interactive = True)
    
#     return s3s_cut, meshes_cut, dict_planes

#%% func - getCuttingPlanesInfo
def getCuttingPlanesInfo(pls_normal, pls_centre):
    
    d_cuts = []
    normal_cuts = []
    d_red = []
    d_green = []
    
    for i in range(len(pls_normal)):
        pl_normal = pls_normal[i]
        pl_centre = pls_centre[i]
        
        # Make normal of plane unitary 
        normal_unit = unit_vector(pl_normal)
        normal_cuts.append(normal_unit)
        normalx10 = [k*100 for k in normal_unit]
        # Create points at either side of the plane
        pt_red = pl_centre+normalx10
        pt_green = pl_centre-normalx10
        
        # Get plane values
        d = normal_unit.dot(pl_centre)
        d_cuts.append(d)
        a, b, c = normal_unit
        
        x_red,y_red,z_red = pt_red
        dotProd_pt_red = dot((a,b,c), (x_red,y_red,z_red))
        d_red.append(dotProd_pt_red)
        x_green,y_green,z_green = pt_green
        dotProd_pt_green = dot((a,b,c), (x_green,y_green,z_green))
        d_green.append(dotProd_pt_green)
        
    return d_cuts, normal_cuts, d_red, d_green

#%% func - cutInfAndOutfOptMx
def cutInfAndOutfOptMx(s3_cut, pls_normal, pls_centre, resolution, mesh_name):
    
    # print('- Cutting s3 - inf&outf '+mesh_name)
    
    # Get dimensions of external stack
    xdim, ydim, zdim = s3_cut.shape
    # Reshape stacks as a vector
    s3_cut_v = s3_cut.reshape(-1)
    
    # Get vectors of x,y and z positions
    pix_coord_pos = np.where(s3_cut >= 0)
    # Trasform coordinate positions to um using resolution
    pix_um = np.transpose(np.asarray([pix_coord_pos[i]*resolution[i] for i in range(len(resolution))]))
    
    normal_inf = unit_vector(pls_normal[0])
    normal_outf = unit_vector(pls_normal[1])
    
    # Find all the d values of pix_um
    d_pix_um_Inf = np.dot(np.subtract(pix_um,np.array(pls_centre[0])),np.array(normal_inf))
    d_pix_um_Outf = np.dot(np.subtract(pix_um,np.array(pls_centre[1])),np.array(normal_outf))
    
    # Clear vector d_pix_um using only those that are 1 in stack
    d_pve_pix_um_Inf = s3_cut_v*d_pix_um_Inf
    d_pve_pix_um_Outf = s3_cut_v*d_pix_um_Outf
    
    # Duplicate s3f_v to initialise stacks without inflow
    s3f_all_v = np.copy(s3_cut_v)
    s3f_all_v.astype('uint8')
    
    # Find all positions in d_pve_pix_um that are at either side of the planes (outside of mesh)
    pos_outside_inf = np.where(d_pve_pix_um_Inf < 0)[0]
    pos_outside_outf = np.where(d_pve_pix_um_Outf > 0)[0]
    
    # Remove the points that are outside of the mesh (inflow)
    s3f_all_v[pos_outside_inf] = 0
    
    # # Duplicate s3f2_v to initialise stacks without inflow and outflow
    # s3f2_all_v = np.copy(s3f_all_v)
    
    # Remove the points that are outside of the mesh (ouflow)
    s3f_all_v[pos_outside_outf] = 0
     
    # Reshape vector into matrix/stack
    s3f_cut = s3f_all_v.reshape((xdim, ydim, zdim))
    
    alert('wohoo',1)
    
    return s3f_cut
    
# OLD VERSION
# #%% func - cutInfAndOutfOpt
# def cutInfAndOutfOpt(s3_cut, s3_cut_int, s3_cut_ext, d_cuts, normal_cuts, d_red, d_green, resolution, mesh_name):
    
#     # Get coordinates of pixels that have a value = 1
#     pix_coord_pos = np.where(s3_cut_ext == 1)
    
#     # Assign classes to points according to option
#     #   Inflow removal: red remove, green keep
#     #   Outflow removal: green remove, red keep
    
#     a_inf, b_inf, c_inf = normal_cuts[0]
#     a_outf, b_outf, c_outf = normal_cuts[1]
#     d_inf, d_outf = d_cuts
#     red_inf, red_outf = d_red
#     green_inf, green_outf = d_green 
    
#     # Cut 
#     bar = Bar('Cutting s3 - inf&outf '+mesh_name, max = 100, suffix = suffix, check_tty=False, hide_cursor=False)
#     num_pts = round(len(pix_coord_pos[0])/100)
#     for i, x, y, z in zip(count(), pix_coord_pos[0], pix_coord_pos[1], pix_coord_pos[2]):
#         # Pix coordinates = [x,y,z]
#         # [x,y,z] * resolution = Coordinates in um [x_um, y_um, z_um]
#         pt_um = [x,y,z]*np.array(resolution)
#         x_um, y_um, z_um = pt_um
#         # Find dot product to know on which side of the plane each pix is (INFLOW)
#         dotProd_pt_inf = dot((a_inf,b_inf,c_inf), (x_um, y_um, z_um))
#         # - Classify point at side of inflow
#         # --> inside mesh, looking from inflow cutting plane
#         if dotProd_pt_inf > d_inf and red_inf > d_inf or dotProd_pt_inf < d_inf and red_inf < d_inf: 
#             do_inf = "inside"
#         # --> outside mesh, looking from inflow cutting plane
#         elif dotProd_pt_inf > d_inf and green_inf > d_inf or dotProd_pt_inf < d_inf and green_inf < d_inf: 
#             do_inf = "outside"
        
#         # Find dot product to know on which side of the plane each pix is (INFLOW)
#         dotProd_pt_outf = dot((a_outf, b_outf, c_outf), (x_um, y_um, z_um))
#         # - Classify point at side of outflow
#         # --> inside mesh, looking from outflow cutting plane
#         if dotProd_pt_outf > d_outf and red_outf > d_outf or dotProd_pt_outf < d_outf and red_outf < d_outf: 
#             do_outf = "outside"
#         # --> outside mesh, looking from outflow cutting plane
#         elif dotProd_pt_outf > d_outf and green_outf > d_outf or dotProd_pt_outf < d_outf and green_outf < d_outf: 
#             do_outf = "inside"
        
#         if do_inf == 'outside' or do_outf == 'outside':
#             s3_cut[x,y,z] = 0
#             s3_cut_int[x,y,z] = 0
#             s3_cut_ext[x,y,z] = 0
            
#         if i % num_pts == 0:
#             bar.next()
        
#     # bar.finish()
#     alert('wohoo',1)
    
#     return s3_cut, s3_cut_int, s3_cut_ext

#%% func - cutInfOrOutfOptMx
def cutInfOrOutfOptMx (s3_cut, pls_normal, pls_centre, resolution, option, mesh_name):
    
    # print('- Cutting s3 - ' + option+' '+mesh_name)
    
    # Get dimensions of external stack
    xdim, ydim, zdim = s3_cut.shape
    # Reshape stacks as a vector
    s3_cut_v = s3_cut.reshape(-1)
    
    # Get vectors of x,y and z positions
    pix_coord_pos = np.where(s3_cut >= 0)
    # Trasform coordinate positions to um using resolution
    pix_um = np.transpose(np.asarray([pix_coord_pos[i]*resolution[i] for i in range(len(resolution))]))
    
    normal  = unit_vector(pls_normal)
    # Find all the d values of pix_um
    d_pix_um = np.dot(np.subtract(pix_um,np.array(pls_centre)),np.array(normal))
    
    # Clear vector d_pix_um using only those that are 1 in stack
    d_pve_pix_um = s3_cut_v*d_pix_um
    
    # Duplicate s3f_v to initialise stacks without inflow/outflow
    s3f_all_v = np.copy(s3_cut_v)
    s3f_all_v.astype('uint8')
    
    # Find all positions in d_pve_pix_um that are at either side of the planes (outside of mesh)
    if option == 'inflow':
        pos_outside = np.where(d_pve_pix_um < 0)[0]
    elif option == 'outflow':
        pos_outside = np.where(d_pve_pix_um > 0)[0]
    
    # Remove the points that are outside of the mesh (inflow/outflow)
    s3f_all_v[pos_outside] = 0
    
    # Reshape vector into matrix/stack
    s3f_cut = s3f_all_v.reshape((xdim, ydim, zdim))
    
    alert('wohoo',1)
            
    return s3f_cut

# OLD VERSION
# #%% func - cutInfOrOutfOpt
# def cutInfOrOutfOpt (s3_cut, s3_cut_int, s3_cut_ext, d_cut, normal_cut, d_red, d_green, resolution, option, mesh_name):
    
#     # Get coordinates of pixels that have a value = 1
#     pix_coord_pos = np.where(s3_cut_ext == 1)
    
#     # Assign classes to points according to option
#     #   Inflow removal: red remove, green keep
#     #   Outflow removal: green remove, red keep
#     if option == "outflow":
#         class_red = "remove"
#         class_green = "keep"
#     elif option == "inflow":
#         class_red = "keep"
#         class_green = "remove"
    
#     a, b, c = normal_cut
#     bar = Bar('Cutting s3 - ' + option+' '+mesh_name, max = 100, suffix = suffix, check_tty=False, hide_cursor=False)
#     num_pts = round(len(pix_coord_pos[0])/100)
#     for i, x, y, z in zip(count(), pix_coord_pos[0], pix_coord_pos[1], pix_coord_pos[2]):
#         # Pix coordinates = [x,y,z]
#         # [x,y,z] * resolution = Coordinates in um [x_um, y_um, z_um]
#         pt_um = [x,y,z]*np.array(resolution)
#         x_um, y_um, z_um = pt_um
#         # Find dot product to know on which side of the plane each pix is
#         dotProd_pt = dot((a,b,c), (x_um, y_um, z_um))
#         if dotProd_pt > d_cut and d_red > d_cut or dotProd_pt < d_cut and d_red < d_cut:
#             do = class_red
#         elif dotProd_pt > d_cut and d_green > d_cut or dotProd_pt < d_cut and d_green < d_cut:
#             do = class_green
        
#         if do == 'remove':
#             s3_cut[x,y,z] = 0
#             s3_cut_int[x,y,z] = 0
#             s3_cut_ext[x,y,z] = 0
            
#         if i % num_pts == 0:
#             bar.next()
        
#     bar.finish()
#     alert('wohoo',1)
            
#     return s3_cut, s3_cut_int, s3_cut_ext

#%% - CREATE/MODIFY OBJECTS 
#%% >>> MESHES
#%% func - createAll3LayerMeshes
def createAll3LayerMeshes(filename, s3_all, s3_in, s3_out, resolution, layer):
    
    print('- Creating surface reconstruction of '+ layer +' ')
    
    # Extract vertices, faces, normals and values of each mesh
    verts_all, faces_all, _, _ = measure.marching_cubes_lewiner(s3_all, spacing=resolution)
    verts_in, faces_in, _, _ = measure.marching_cubes_lewiner(s3_in, spacing=resolution)
    verts_out, faces_out, _, _ = measure.marching_cubes_lewiner(s3_out, spacing=resolution)
    alert('frog',1)
    
    # Create meshes
    mesh_all = Mesh([verts_all, faces_all])
    mesh_all.rotateZ(-90).wireframe(True)
    
    mesh_in = Mesh([verts_in, faces_in])
    mesh_in.rotateZ(-90).wireframe(True)
    
    mesh_out = Mesh([verts_out, faces_out])
    mesh_out.rotateZ(-90).wireframe(True)
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
        mesh_all.color("indigo").alpha(0.05).wireframe().legend(layer)
        mesh_in.color("indigo").alpha(0.05).wireframe().legend('Int.'+layer)
        mesh_out.color("indigo").alpha(1).wireframe().legend('Ext.'+layer)
        
    text = filename+"\n\n >> "+layer
    txt = Text2D(text, c=c, font=font)
    vp = Plotter(N=3, axes=7)
    vp.show(mesh_all, txt, at=0, zoom=1.2)
    vp.show(mesh_all, mesh_out, at=1, zoom=1.2)
    vp.show(mesh_all, mesh_in, at=2, zoom=1.2, interactive=True)
    
    return mesh_all, mesh_in, mesh_out

#%% func - createExtLayerMesh
def createExtLayerMesh(filename, s3_ext, resolution, layer, info, plotshow = True):
    
    print('- Creating surface reconstruction of '+ layer +' ('+info+')')
    
    # Extract vertices, faces, normals and values of each mesh
    verts_ext, faces_ext, _, _ = measure.marching_cubes_lewiner(s3_ext, spacing=resolution)
    alert('frog',1)
    
    mesh_ext = Mesh([verts_ext, faces_ext])
    mesh_ext.rotateZ(-90).wireframe(True)        
    alert('clown',1)
    
    if layer == 'Myoc' or layer == 'Endo':
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
        vp = Plotter(N=1, axes=7)
        vp.show(mesh_ext, txt, at=0, zoom=1.2, interactive=True)
    
    return mesh_ext

#%% func - createLayerMesh
def createLayerMesh(filename, s3, resolution, layer, name, colour, alpha, plotshow = False):
    
    print('- Creating surface reconstructions of '+ name )
    
    # Extract vertices, faces, normals and values of each mesh
    verts, faces, _, _ = measure.marching_cubes_lewiner(s3, spacing=resolution)
    alert('frog',1)
    
    # Create meshes
    mesh = Mesh([verts, faces])
    if layer == 'Myoc' or layer == 'Endo':
        mesh = mesh.extractLargestRegion()
        
    mesh.rotateZ(-90).color(colour).alpha(alpha).legend(name).wireframe()
    alert('clown',1)
    
    if plotshow:
        text = filename+"\n\n >> "+layer
        txt = Text2D(text, c=c, font=font)
        
        vp = Plotter(N=1, axes=7)
        vp.show(mesh, txt, at=0, zoom=1.2, interactive=True)
    
    return mesh

#%% func - colorMesh
def colorMesh(mesh, m_name, alpha = [1, 1, 1, 1, 1]):

    meshes_names = np.array(['myoc', 'endo', 'cj', 'cjOut', 'cjIn'])
    colors = ['darkcyan','darkmagenta','darkorange','limegreen','gold']
    legend = ['Myocardium','Endocardium','Cardiac Jelly','CJ (Out)','CJ (In)']
    
    num = np.where(meshes_names == m_name)[0][0]
    
    mesh.alpha(alpha[num]).legend(legend[num]).wireframe().color(colors[num])
    
    return mesh    

#%% func - getCutMesh
def getCutMesh(filename, s3_cut, resolution, mesh_original, layer, plotshow):
    
    print('- Creating surface reconstructions of cut '+ layer)
    
    verts_cut, faces_cut, _, _ = measure.marching_cubes_lewiner(s3_cut, spacing=resolution)
    mesh_cut = Mesh([verts_cut, faces_cut])
    if layer != 'Int.Endo':
        mesh_cut = mesh_cut.extractLargestRegion()
    alert('wohoo',1)
    
    if layer == 'Cardiac Jelly' or layer == 'CJ': 
        color = 'darkorange'
    elif layer == 'Myocardium' or layer == 'Myoc':
        color = 'darkcyan'
    elif layer == 'Endocardium' or layer == 'Endo':
        color = 'darkmagenta'
    else: 
        color = 'indigo'
    
    mesh_cut.rotateZ(-90).color(color).legend('Cut ('+layer+')')
    
    text = filename 
    txt = Text2D(text, c=c, font=font)
    
    if plotshow:
        mesh_original.legend(layer).alpha(0.1).color('tomato')
        vp = Plotter(N=1, axes=4)
        vp.show(mesh_original, mesh_cut, txt, at=0, zoom=1.2, interactive=True)
    
    return mesh_cut

#%% func - createMeshes4CL
def createMeshes4CL(filename, meshes, names, mesh_colors, plotshow):
    
    meshes4cl = []

    for i, mesh_cl in enumerate(meshes):
        print("- File being processed for centreline: ", filename +' - '+names[i])
        mesh_cl.legend(names[i])
        
        if plotshow:        
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
            vp.show(meshCL_cut, txt2, at=1, interactive=1, axes=4)
        
        meshes4cl.append(meshCL_cut)
    
    vp = Plotter(N=2, axes=7)
    vp.show(meshes4cl[-1], at=0, zoom=1)
    vp.show(meshes4cl[-2], at=1, zoom=1, interactive = True)
    
    meshes4cl_out = meshes4cl[-2:]
        
    return meshes4cl_out

#%% func - cutMeshes4CL
def cutMeshes4CL(filename, meshes, names, cuts, cut_direction, mark_colors, mesh_colors, dicts, plotshow):
    
    dict_planes, dict_pts, dict_kspl = dicts 
    ksplines = []; spheres = []
    
    for n, cut in enumerate(cuts):
        print('- Cutting: '+cut)
        #Get plane to cut
        planeCL_cut, plCL_cut_centre, plCL_cut_normal = getPlane(filename = filename, type_cut = cut, 
                                                                 info = '4CL', mesh_in = meshes[1],
                                                                 mesh_out = meshes[0])
        dict_planes = addPlane2Dict (plane = planeCL_cut, pl_centre = plCL_cut_centre, 
                                            pl_normal = plCL_cut_normal, info = '', dict_planes = dict_planes)
        
        for i, mesh4cl in enumerate(meshes[-2:]):
            pts2cut = getPointsAtPlane(points = mesh4cl.points(), pl_normal = plCL_cut_normal, 
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
            
            # # Plot result 
            # vp = Plotter(N=1)
            # text = filename+"\n\n - Mesh - KSpline and centre of cut "+ cut +"-"+names[i]
            # txt = Text2D(text, c=c, font=font)
            # vp.show(mesh4cl, kspl, sph_centroid, txt, at=0, axes=1, interactive=True)
            
            dict_pts = addPoint2Dict(sphere = sph_centroid, info = names[i], dict_pts = dict_pts)
            dict_kspl = addKSpline2Dict(kspl = kspl, info = names[i], dict_kspl = dict_kspl)
                
            # Cutmesh using created plane
            mesh4cl = mesh4cl.clone().cutWithMesh(planeCL_cut, invert=cut_direction[n])
            mesh4cl = mesh4cl.extractLargestRegion()
            mesh4cl.color(mesh_colors[i]).alpha(0.05).wireframe(True).legend('CL_cutMesh')
            
            vp = Plotter(N=1)
            text = "- Resulting mesh after cutting"
            txt = Text2D(text, c=c, font=font)
            vp.show(mesh4cl, kspl, sph_centroid, txt, at=0, viewup="y",  interactive=True)
            
            meshes.append(mesh4cl)
        
    if plotshow:
        vp = Plotter(N=1, axes=7)
        text = filename+"\n\n >> Resulting mesh after cutting inflow & outflow tract"
        txt = Text2D(text, c=c, font=font)
        vp.show(meshes[-2:], ksplines, spheres, txt, at=0, interactive=True)
    
    dicts_f = dict_planes, dict_pts, dict_kspl
    
    return meshes[-2:], dicts_f

#%% func - divideMeshesLnR
def divideMeshesLnR(filename, meshes, cl_ribbon):
    
    meshes_LnR = []
    for i, mesh in enumerate(meshes):
        
        meshes2cutLR = mesh
        mesh_legend = mesh.legend()
        print('- Dividing '+mesh_legend+' into left and right sides')
        
        mesh_1 = meshes2cutLR.clone().cutWithMesh(cl_ribbon, invert=True)
        mesh_1.alpha(0.05).wireframe(True).color('deepskyblue')
        
        mesh_2 = meshes2cutLR.clone().cutWithMesh(cl_ribbon, invert=False)
        mesh_2.alpha(0.05).wireframe(True).color('darkblue')
        alert("wohoo",1)
        
        mesh_1.legend('Left(A)')
        mesh_2.legend('Right(B)')
            
        text = filename+"\n\n >> Resulting mesh after cutting with centreline \n >> Before closing the window make sure you confirm A: left, B: Right"
        txt = Text2D(text, c=c, font=font)
        vp = Plotter(N=2, axes=7)
        vp.show(meshes2cutLR, cl_ribbon, txt, at=0)
        vp.show(mesh_1, mesh_2, at=1, zoom = 1.2, interactive=True)
        
        q_happy = ask4input('Are meshes classified correctly [A: left, B: Right]? [0]:No / [1]:Yes: ', int)
        if q_happy == 0:
            leftMesh = ask4input('Select the mesh number that corresponds to the left side [1]:'+mesh_legend+'-(1) / [2]:'+mesh_legend+'-(2): ', int)
            if leftMesh == 1:
                mesh_1.legend(mesh_legend+'-Left')
                mesh_2.legend(mesh_legend+'-Right')
                mesh_LnR = [mesh_1, mesh_2] 
            elif leftMesh == 2:
                mesh_2.legend(mesh_legend+'-Left')
                mesh_1.legend(mesh_legend+'-Right')
                mesh_LnR = [mesh_2, mesh_1] 
        
        else: 
            mesh_LnR = [mesh_1, mesh_2] 
            
        meshes_LnR.append(mesh_LnR)
        
    return meshes_LnR

#%% func - getChamberMeshes
def getChamberMeshes(filename, s3s2cut, names2cut, kspl_CL, mesh2cut, resolution, dict_planes, dict_pts):
    
    dORv = filename[9:10]
    if dORv == 'D':
        azimuth = -90
    else: 
        azimuth = 0
    
    # Plot spheres through centreline inside myocardium
    spheres_spl = sphInSpline(kspl_CL = kspl_CL)
    vp = Plotter(N=1, axes = 4)
    text = filename+"\n\n >> Decide the centreline point number to use to initialise plane to divide chambers \n"
    txt = Text2D(text, c="k", font= 'CallingCode')
    vp.show(mesh2cut.alpha(0.01), kspl_CL, spheres_spl, txt, at=0, azimuth = azimuth, interactive=True)
    
    # Get plane to cut heart layers
    plane_Ch, pl_Ch_centre, pl_Ch_normal = getPlane4ChDivision(filename = filename, type_cut = 'Chamber', 
                                                                    mesh1 = mesh2cut, kspl_CL = kspl_CL)
    # Reorient plane to images (s3)
    plane_imCh, pl_imCh_centre, pl_imCh_normal = rotatePlane2Images(pl_Ch_centre, pl_Ch_normal, type_cut = 'Chamber')
    
    # Save planes to dict
    dict_planes = addPlanes2Dict(planes = [plane_Ch, plane_imCh], pls_centre = [pl_Ch_centre ,pl_imCh_centre], 
                             pls_normal = [pl_Ch_normal, pl_imCh_normal], info = ['',''], dict_planes = dict_planes)
    
    # Cut cl with plane
    ksplCL_cut = kspl_CL.clone().cutWithMesh(plane_Ch)
    # Find point of centreline closer to last point of kspline cut
    ksplCL_cutPt, num_pt = findClosestPt(ksplCL_cut.points()[0], kspl_CL.points())
    # Create sphere for that point
    sph_cut = Sphere(pos = ksplCL_cutPt, r=4, c='gold').legend('sph_ChamberCut')
    # Add pt to dict
    dict_pts = addPoint2Dict (sphere = sph_cut, info = '', dict_pts = dict_pts)
    dict_pts['numPt_CLChamberCut'] = num_pt
    
    # Create empty lists to save atrium and ventricles
    atr_meshes = []
    atr_color = ['lightseagreen', 'purple', 'orange', 'darkblue', 'indigo']
    vent_meshes = []
    vent_color = ['darkturquoise', 'mediumvioletred', 'chocolate', 'indigo', 'darkblue']
    
    tic = perf_counter()
    vp = Plotter (N=3, axes=10)
    text2 = filename+"\n\n >>  Result of dividing heart layers into chambers"
    txt2 = Text2D(text2, c="k", font= 'CallingCode')
    for n, s3, name in zip(count(), s3s2cut, names2cut):
        # Mask s3s vent and atrium
        print('- Cutting s3 (', name,')')
        s3_vent, s3_atr = maskChamberS3s(s3_mask = s3, pl_normal = pl_imCh_normal, pl_centre = pl_imCh_centre, resolution = resolution)
        # Create chamber meshes 
        atr = createLayerMesh(filename = filename, s3 = s3_atr, resolution = resolution, layer = name, name = name+'_Atr',
                                       colour = atr_color[n], alpha = 0.01, plotshow = False)
        atr_meshes.append(atr)
        vent = createLayerMesh(filename = filename, s3 = s3_vent, resolution = resolution, layer = name, name = name+'_Vent',
                                        colour = vent_color[n], alpha = 0.01, plotshow = False)
        vent_meshes.append(vent)
        if n == 0:
            vp.show(atr, vent, kspl_CL, spheres_spl, sph_cut, txt2, at=n)
        elif n == 1:
            vp.show(atr, vent, kspl_CL, spheres_spl, sph_cut, at=n)
        elif n == 2:
            vp.show(atr, vent, kspl_CL, spheres_spl, sph_cut, at=n, azimuth = azimuth, interactive = True)
    toc = perf_counter()
    time = toc-tic
    print("- All layers have been cut!  > Total time taken to cut meshes = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")
    alert('jump', 1)
    
    return atr_meshes, vent_meshes, plane_Ch, dict_planes, dict_pts, num_pt

#%% >>> PLANES
#%% func - createPlane
def createPlane(dict_planes, name):
    normal = dict_planes[name]['pl_normal'] 
    centre = dict_planes[name]['pl_centre'] 
    color = dict_planes[name]['color']
    
    plane_out = Plane(pos = centre, normal = normal, sx = 300).color(color).alpha(1)
    plane_out.legend(name)
    
    return plane_out

#%% func - getPlane
def getPlane (filename, type_cut, info, mesh_in, mesh_out, option = [True,True,True,True,True,True]):
    print('- Getting plane to cut '+ type_cut)
    #mesh1.alpha(0.05); mesh2.alpha(0.05)
    while True:
        # Create plane
        plane, normal, rotX, rotY, rotZ = getPlanePos(filename, type_cut, mesh_in, mesh_out, mesh_in.bounds(), option)
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
        
        vp = Plotter(N=1, axes=4)
        vp.show(mesh_in, mesh_out, plane, plane_new, sph_centre, txt, at=0, viewup="y", azimuth=0, elevation=0, interactive=True)
        happy = ask4input('Do you want to cut the '+type_cut+' with the defined plane? \n  [0]:no, I would like to define a new plane / [1]:yes, continue:', int)
        if happy == 1:
            break

    return plane_new, pl_centre, normal_corrected

#%% func - getPlanePos
def getPlanePos (filename, type_cut, mesh_in, mesh_out, xyz_bounds, option):
    
    xmin, xmax, ymin, ymax, zmin, zmax = xyz_bounds
    x_size = xmax - xmin
    y_size = ymax - ymin
    z_size = zmax - zmin

    box_size = max(x_size, y_size, z_size)*1.2
    centre = (x_size/2+xmin, ymin, z_size/2+zmin)
    normal = (0,1,0)
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
    
    vp = Plotter(N=1, axes=8)
    plane = Plane(pos=centre, normal=normal, sx=box_size*1.2).color("gainsboro").alpha(1)
    
    if option[0]: #sliderX
        vp.addSlider2D(sliderX, xmin*0.8, xmax*1.2, value=centre[0], 
                    pos=[(0.1,0.05), (0.3,0.05)], title="- > x position > +", c="crimson" )
    if option[1]: #sliderY
        vp.addSlider2D(sliderY, ymin*1.2, ymax*0.8, value=centre[1], 
                    pos=[(0.4,0.05), (0.6,0.05)], title="- > y position > +", c="dodgerblue" )
    if option[2]: #sliderZ
        vp.addSlider2D(sliderZ, zmin*1.2, zmax*0.8, value=centre[2], 
                    pos=[(0.7,0.05), (0.9,0.05)], title="- > z position > +", c="limegreen")
    if option[3]: #sliderRotX
        vp.addSlider2D(sliderRotX, -1, +1, value=0, 
                    pos=[(0.1,0.15), (0.3,0.15)], title="- > x rotation > +", c="deeppink")
    if option[4]: #sliderRotY
        vp.addSlider2D(sliderRotY, -1, +1, value=0, 
                    pos=[(0.4,0.15), (0.6,0.15)], title="- > y rotation > +", c="gold")
    if option[5]: #sliderRotZ
        vp.addSlider2D(sliderRotZ, -1, +1, value=0, 
                    pos=[(0.7,0.15), (0.9,0.15)], title="- > z rotation > +", c="teal")
    
    vp.addSlider2D(sliderAlphaMeshOut, xmin=0.01, xmax=0.99, value=0.01,
               pos=[(0.95,0.25), (0.95,0.45)], c="blue", title="Ext.Mesh Opacity)")
    
    text = filename+"\n\n >> Define plane position to make cut ("+type_cut+")\n >> Close the window when done"
    txt = Text2D(text, c="k", font= 'CallingCode')
    vp.show(mesh_in, mesh_out, plane, txt, viewup="y", zoom=1, interactive=True)
    #azimuth=-90, elevation=0,
    
    return plane, normal, rotX, rotY, rotZ

#%% func - getPlane4ChDivision
def getPlane4ChDivision (filename, type_cut, mesh1, kspl_CL, option = [True,True,True,True,True,True]):
    print('- Getting plane to divide '+ type_cut)
    mesh1.alpha(0.05)
    while True:
        # Create plane
        pt_num = ask4input('Enter the centreline point number you want to use to initialise plane to divide '+type_cut+': ', int)
        pl_Ch_normal, pt_Ch_centre = getPlaneNormal2Pt (pt_num = pt_num, spline_pts = kspl_CL.points())
        # Modify (rotate and move plane)
        plane_Ch, rotX, rotY, rotZ = modifyPlane (filename = filename, 
                                                pl_normal = pl_Ch_normal,
                                                pl_centre = pt_Ch_centre, 
                                                type_cut = type_cut, 
                                                mesh1 = mesh1,
                                                xyz_bounds = mesh1.bounds(), 
                                                option = option)
        # Get new normal of rotated plane
        pl_Ch_normal_corrected = newNormal3DRot(normal = pl_Ch_normal, rotX = rotX, rotY = rotY, rotZ = rotZ)
        # Get central point of new plane and create sphere
        pl_Ch_centre = plane_Ch.pos()
        sph_centre = Sphere(pos=pl_Ch_centre, r=2, c="purple")
        # Build new plane to confirm
        plane_Ch_final = Plane(pos=pl_Ch_centre,normal=pl_Ch_normal_corrected, sx=400).color("green").alpha(0.5).legend('pl2CutMesh'+'_'+type_cut)
        
        text = filename+"\n\n >> Confirm plane position to proceed with division of "+type_cut+"\n >> Close the window when done and define whether you want to use this plane or define a new one"
        txt = Text2D(text, c=c, font=font)
        vp = Plotter(N=1, axes=4)
        vp.show(mesh1, plane_Ch, plane_Ch_final, sph_centre, txt, at=0, viewup="y", azimuth=0, elevation=0, interactive=True)
        happy = ask4input('Do you want to cut the '+type_cut+' with the defined plane? \n   [0]:no, I would like to define a new plane / [1]:yes, continue with the cut >> :', int)
        if happy == 1:
            break

    return plane_Ch_final, pl_Ch_centre, pl_Ch_normal_corrected

#%% func - modifyPlane
def modifyPlane (filename, pl_normal, pl_centre, type_cut, mesh1, xyz_bounds, option):
    
    xmin, xmax, ymin, ymax, zmin, zmax = xyz_bounds
    x_size = xmax - xmin
    y_size = ymax - ymin
    z_size = zmax - zmin

    box_size = max(x_size, y_size, z_size)*1.2
    #centre = (x_size/2+xmin, ymin, z_size/2+zmin)
    #normal = (0,1,0)
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
    
    vp = Plotter(N=1, axes=8)
    plane = Plane(pos=pl_centre, normal=pl_normal, sx=box_size*1.2).color("gainsboro").alpha(1)
    
    if option[0]: #sliderX
        vp.addSlider2D(sliderX, xmin*0.8, xmax*1.2, value=pl_centre[0], 
                    pos=[(0.1,0.05), (0.3,0.05)], title="- > x position > +", c="crimson" )
    if option[1]: #sliderY
        vp.addSlider2D(sliderY, ymin*1.2, ymax*0.8, value=pl_centre[1], 
                    pos=[(0.4,0.05), (0.6,0.05)], title="- > y position > +", c="dodgerblue" )
    if option[2]: #sliderZ
        vp.addSlider2D(sliderZ, zmin*1.2, zmax*0.8, value=pl_centre[2], 
                    pos=[(0.7,0.05), (0.9,0.05)], title="- > z position > +", c="limegreen")
    if option[3]: #sliderRotX
        vp.addSlider2D(sliderRotX, -1, +1, value=0, 
                    pos=[(0.1,0.15), (0.3,0.15)], title="- > x rotation > +", c="deeppink")
    if option[4]: #sliderRotY
        vp.addSlider2D(sliderRotY, -1, +1, value=0, 
                    pos=[(0.4,0.15), (0.6,0.15)], title="- > y rotation > +", c="gold")
    if option[5]: #sliderRotZ
        vp.addSlider2D(sliderRotZ, -1, +1, value=0, 
                    pos=[(0.7,0.15), (0.9,0.15)], title="- > z rotation > +", c="teal")
    
    text = filename+"\n\n >> Define plane position to make cut ("+type_cut+")\n >> Close the window when done"
    txt = Text2D(text, c="k", font= 'CallingCode')
    vp.show(mesh1, plane, txt, viewup="y", zoom=1, interactive=True)
    #azimuth=-90, elevation=0,
    
    return plane, rotX, rotY, rotZ

#%% func - getPlaneNormal2Pt
def getPlaneNormal2Pt (pt_num, spline_pts):
    pt_centre = spline_pts[pt_num]
    normal = spline_pts[pt_num-1]-spline_pts[pt_num+1]
    
    return normal, pt_centre

#%% func - rotatePlane2Images
def rotatePlane2Images (pl_centre, pl_normal, type_cut):
    pl_im_centre = (np.dot(rotation_matrix(axis = [0,0,1], theta = np.radians(90)), pl_centre)) 
    pl_im_normal = (np.dot(rotation_matrix(axis = [0,0,1], theta = np.radians(90)), pl_normal)) 
    
    plane_im = Plane(pos=pl_im_centre,normal=pl_im_normal, sx=500).color("blue").alpha(0.5).legend('pl2CutIm'+'_'+type_cut)

    return plane_im, pl_im_centre, pl_im_normal

#%% func - createDVPlanes
def createDVPlanes(filename, sph_orient, mesh, kspl_CL, orient_lines, dict_planes):
    
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
    txt = Text2D(text, c="k", font= 'CallingCode')
    vp = Plotter(N=1, axes = 7)
    vp.show(mesh.alpha(0.01), sph_orient, kspl_CL, orient_lines, pl_DnV_Atr, pl_DnV_Vent, txt, at=0, interactive=True)
    
    pl_DnV = [pl_DnV_Atr, pl_DnV_Vent]
    
    return pl_DnV, dict_planes
    
#%% >>> SHAPES
#%% func - sphInSpline
def sphInSpline(kspl_CL, name = '', every = 10):
    
    if every > 1:
        spheres_spline = []
        for num, point in enumerate(kspl_CL.points()):
            if num % every == 0 or num == kspl_CL.NPoints()-1:
                sphere_pt = Sphere(pos=point, r=2, c='coral')
                spheres_spline.append(sphere_pt)
    else: 
        kspl_new = KSpline(kspl_CL.points(), res = round(kspl_CL.NPoints()/every))
        spheres_spline = Spheres(kspl_new.points(), c='coral', r=2)
        if name != '':
            spheres_spline.legend(name)
    
    return spheres_spline

#%% func - createKSpls
def createKSpls(dict_kspl, kspl_list):#, colors):
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
def createSpheres(dict_pts, pts_list):#, colors):
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
def createCLs(dict_cl, dict_pts, dict_kspl, dict_shapes, colors):
    
    print('- Creating centreline and associated objects...')
    kspl_CL = []
    linLines = []
    spheres_CL = []
    spheres_CL_col = []
    
    layerNames = ['Int.Myoc(Cut)', 'Ext.Endo(Cut)']
    sph_genName  = ['sph_Cut4CL_inflow-', 'sph_Cut4CL_outflow-']
    color_linLines = ['lawngreen', 'blueviolet']
    
    # Iterate through layer
    for i, dictLayer in enumerate(dict_cl):
        #Get cl points from vmtk
        pts_cl = np.asarray(dictLayer['Points'])
        sphData_cl = np.asarray(dictLayer['PointData']['MaximumInscribedSphereRadius'])
        
        #Get inflow and outflow point for that layer
        sph_inf = sph_genName[0]+layerNames[i]
        pt2add_inf = dict_pts[sph_inf]['sph_position']
        sph_outf = sph_genName[1]+layerNames[i]
        pt2add_outf = dict_pts[sph_outf]['sph_position']
        
        pts_withOutf = np.insert(pts_cl, 0, np.transpose(pt2add_outf), axis=0)
        pts_all = np.insert(pts_withOutf, len(pts_withOutf), np.transpose(pt2add_inf), axis=0)
        
        # Interpolate points
        pts_int = getInterpolatedPts(points=pts_all, nPoints = 300)
        # Create kspline with points
        kspl = KSpline(pts_int, res = 300)

        kspl.color(colors[i]).legend('CL_'+layerNames[i]).lw(3)
        kspl_CL.append(kspl)
        
        linearLine = Line(pt2add_inf, pt2add_outf, c=color_linLines[i], lw=3)
        linearLine.legend('linLine_'+layerNames[i])
        linLines.append(linearLine)
        
        #Plot spheres through centreline inside mesh.
        vcols = [colorMap(v, "jet", sphData_cl.max(), sphData_cl.min()) for v in sphData_cl]
        sph_cl = Spheres(pts_cl, c='red', r=sphData_cl).legend('sphs_maxInsSphRad_'+layerNames[i])
        sph_cl_colour = Spheres(pts_cl, c=vcols, r=sphData_cl.min()).addScalarBar(title='Sphere Radius\n[um]').legend('sphs_maxInsSphRadC_'+layerNames[i])
        sph_cl_colour.mapper().SetScalarRange(sphData_cl.min(),sphData_cl.max())
        spheres_CL.append(sph_cl)
        spheres_CL_col.append(sph_cl_colour)
        
        dict_shapes = addShapes2Dict(shapes = [sph_cl, sph_cl_colour], info = ['',''], dict_shapes = dict_shapes, radius = [[],sphData_cl])
        dict_kspl = addKSplines2Dict(kspls = [linearLine, kspl], info = ['',''], dict_kspl = dict_kspl)
        
    return kspl_CL, linLines, spheres_CL, spheres_CL_col, dict_shapes, dict_kspl

#%% func - createCLRibbon
def createCLRibbon(filename, kspl_CL2use, linLine, mesh, dict_kspl, dict_shapes, dict_planes):
        
    pts_cl = kspl_CL2use.points()
    # Extended centreline
    inf_ext_normal = (pts_cl[-1]+(pts_cl[-1]-pts_cl[-3])*70)
    outf_ext_normal = (pts_cl[0]+(pts_cl[0]-pts_cl[1])*70)
    inf_ext_sphere = Sphere(pos=inf_ext_normal, r=3, c='purple').legend("sph_infCLExtended")
    outf_ext_sphere = Sphere(pos=outf_ext_normal, r=3, c='purple').legend("sph_outfCLExtended")
    
    pts_cl_ext = np.insert(pts_cl,0,np.transpose(outf_ext_normal), axis=0)
    pts_cl_ext = np.insert(pts_cl_ext,len(pts_cl_ext),np.transpose(inf_ext_normal), axis=0)
    kspl_ext = KSpline(pts_cl_ext, res=220).color('purple')
    
    # Create plane to project centreline
    dORv = filename[9:10]
    if dORv == 'V':
        azimuth = 0
        linLineX = linLine.clone().projectOnPlane('x').c(linLine.color()).x(0)

    elif dORv == 'D': #Acaaa
        azimuth = -90
        linLineX = linLine.clone().projectOnPlane('z').c(linLine.color()).z(0)
    
    ptsPl_linLine = Points([linLineX.points()[0], linLine.points()[0], linLine.points()[1]])
    pl_linLine = fitPlane(ptsPl_linLine.points()).scale(4).c('mediumaquamarine').alpha(1).legend('pl_Parallel2LinLine')
    pl_linLine_normal = pl_linLine.normal
    pl_linLine_centre = pl_linLine.center
    dict_planes = addPlane2Dict (plane = pl_linLine, pl_centre = pl_linLine_centre, 
                                            pl_normal = pl_linLine_normal, info = '', dict_planes = dict_planes)
    
    pl_linLine_unitNormal = unit_vector(pl_linLine_normal)*120
    x_cl, y_cl, z_cl = pl_linLine_unitNormal
    
    kspl_ext_D = kspl_ext.clone().x(x_cl).y(y_cl).z(z_cl).legend('kspl_CLExtD')
    kspl_ext_V = kspl_ext.clone().x(-x_cl).y(-y_cl).z(-z_cl).legend('kspl_CLExtV')
    cl_ribbon = Ribbon(kspl_ext_D, kspl_ext_V, alpha=0.2, res=(220, 5))
    cl_ribbon = cl_ribbon.wireframe(True).legend("rib_ExtCL(D-V)")
    
    dict_shapes = addShapes2Dict (shapes = [cl_ribbon], info = [''], dict_shapes = dict_shapes, radius = [[]])
    dict_kspl = addKSplines2Dict(kspls = [kspl_ext_D, kspl_ext_V], info = ['',''], dict_kspl = dict_kspl)
    
    text = filename+"\n\n >> Creating Extended Centreline to Divide Right/Left"    
    txt = Text2D(text, c="k", font= 'CallingCode')
    
    vp = Plotter(N=1, axes=7)
    vp.show(mesh, kspl_CL2use, kspl_ext, inf_ext_sphere, outf_ext_sphere, cl_ribbon, txt, at=0, azimuth = azimuth, interactive=1)
    
    return cl_ribbon, dict_kspl, dict_shapes, dict_planes

#%% func - createBalloonedHeart # NOT WORKING
# def createBalloonedHeart(dictLayer):
    
#     pts_cl = np.asarray(dictLayer['Points'])
#     sphData_cl = np.asarray(dictLayer['PointData']['MaximumInscribedSphereRadius'])
    
#     for n, pt_pos, radius in zip(count(), pts_cl, sphData_cl):
#         print('n:',n)
#         sph2add = Sphere(pos=pt_pos, r=radius, c='red')
#         if n == 0:
#             sph_whole = sph2add
#         elif n > 1:
#             sph_whole = booleanOperation(sph2add, "plus", sph_whole.clone())
     
#     sph_1 = Sphere(pos=pts_cl[1], r=sphData_cl[1], c='red')
    
#     sph_whole = booleanOperation(sph_1, "plus", sph_whole.clone())
#     sph_whole = sph_whole.color('coral').legend('FilledBalloon')
    
#     return sph_all, sph_whole

#%% - OBJECTS
#%% func - rotateObjects
def rotateObjects(objects, angles, first):
    angX, angY, angZ = angles
    objects_out = []
    
    if first:
        for i, obj in enumerate(objects):
            obj.rotateX(angX).rotateY(angY).rotate(angZ)
            
            objects_out.append(obj)
    
        first = False
    #print('- Objects have been rotated (x,y,z): ' + str(angles)+'!')
    
    return objects_out, first

#%% func - rotateAll
def rotateAll (objects_all, angles, first):
    
    objects_out = []
    for j, obj_list in enumerate(objects_all):
        obj_list, _ = rotateObjects(objects = obj_list, angles = angles, first = first)
        objects_out.append(obj_list)
        
        if j == len(objects_all)-1:
            first_out = False
    print('- All objects have been rotated (x,y,z): ' + str(angles)+'!')
        
    return objects_out, first_out

#%% - MEASURE
#%% func - getMeshInclinationAngle 
def getMeshInclinationAngle(filename, dict_pts, m_myoc):
    
    pts = ['Sph_CL_inflow-Int.Myoc','Sph_CL_outflow-Int.Myoc']
    pts_pos = []
    
    for num, pt_name in enumerate(pts):
        pt_pos = dict_pts[pt_name]['sph_position']
        pts_pos.append(pt_pos)
    
    linearLine = Line(pts_pos[0], pts_pos[1], c='magenta', lw=3).legend('linLine (LL)') 
    _, y_inf, z_inf = pts_pos[0]
    _, y_outf, z_outf = pts_pos[1]
            
    linLineX = linearLine.clone().projectOnPlane('x').c('deeppink').x(-20).legend('LL (ProjX)') 

    angle = math.degrees(math.atan2((y_outf-y_inf),(z_outf-z_inf)))-90
    angle_print = format(angle, '.1f')
    print('- Angle to rotate meshes (deg):', angle_print)
    
    m_myoc.alpha(0.01)
    text = filename+"\n\n >> Angle to rotate meshes (deg):"+angle_print
    txt = Text2D(text, c="k", font= 'CallingCode')
    text2 = ">> Rotated mesh"
    txt2 = Text2D(text2, c="k", font= 'CallingCode')
    vp = Plotter(N=2, axes=1, sharecam = False)
    vp.show(m_myoc, linearLine, linLineX, txt, at=0)
    vp.show(m_myoc.clone().rotateX(angle), linearLine.clone().rotateX(angle), linLineX.clone().rotateX(angle), txt2, viewup="y", at=1, interactive = True)
    
    return angle
    
#%% func - getChambersOrientation
def getChambersOrientation(filename, file_num, num_pt, kspl_CL2use, myoc_meshes, linLine, dict_pts, dict_kspl, df_res):
    
    print('- Measuring chamber orientations')
    m_myoc, m_atrMyoc, m_ventMyoc = myoc_meshes
    atr_pt = num_pt+25
    vent_pt = num_pt-25
    
    # Create spheres
    sph_atr = Sphere(pos = kspl_CL2use.points()[atr_pt], r=4, c='navy').legend('sph_AtrCentre')
    sph_vent = Sphere(pos = kspl_CL2use.points()[vent_pt], r=4, c='red').legend('sph_VentCentre')
    sph_inf = Sphere(pos = kspl_CL2use.points()[-1], r=4, c='dodgerblue').legend('sph_OutflowCentre')
    sph_outf = Sphere(pos = kspl_CL2use.points()[0], r=4, c='tomato').legend('sph_InflowCentre')
    sph_valve = Sphere(pos = kspl_CL2use.points()[num_pt], r=4, c='darkorange').legend('sph_AVCCentre')
    
    sph_orient = [sph_atr, sph_vent, sph_outf, sph_inf, sph_valve]
    dict_pts = addPoints2Dict(spheres = sph_orient, info = ['Orient','Orient','Orient','Orient','Orient'], dict_pts = dict_pts)
    
    # Find orientation in which images where taken
    dORv = filename[9:10]
    # Heart Orientation     
    if dORv == 'V':
        linLineX = linLine.clone().projectOnPlane('x').c(linLine.color()).x(0).legend('linLine(ProjX)') 
        pts_heart = orientVectors(linLineX)
        # print('pts_heart', pts_heart)
        ang_heart = findAngleBtwVectors(pts_heart, np.array([[0,1,0],[0,0,0]]))
        azimuth = -45; elevation = 0
    elif dORv == 'D':
        linLineX = linLine.clone().projectOnPlane('z').c(linLine.color()).x(0).legend('linLine(ProjX)') 
        pts_heart = orientVectors(linLineX)
        # print('pts_heart', pts_heart)
        ang_heart = findAngleBtwVectors(pts_heart, np.array([[0,0,0],[0,1,0]]))
        azimuth = -135; elevation = 0
    
    # Atrial orientation
    orient_atr = Line(sph_inf.pos(), sph_atr.pos(), c="steelblue", lw=3).legend('lin_OrientAtr')
    _, y_atrInf, z_atrInf = sph_inf.pos()
    _, y_atrCent, z_atrCent = sph_atr.pos()
    if dORv == 'V':
        orient_atrX = orient_atr.clone().projectOnPlane('x').c('steelblue').x(0).legend('lin_OrientAtr(ProjX)')
    elif dORv == 'D':
        orient_atrX = orient_atr.clone().projectOnPlane('z').c('steelblue').x(0).legend('lin_OrientAtr(ProjX)')
    
    pts_atr = orientVectors(orient_atrX)
    ang_atr = findAngleBtwVectors(pts_atr, pts_heart)
    atr_txt = 'Atrial orientation with respect to heart (deg): '+ format(ang_atr,'.1f')
    print('\t>> '+atr_txt)
    
    # Ventricular orientation
    orient_vent = Line(sph_vent.pos(), sph_outf.pos(), c="hotpink", lw=3).legend('lin_OrientVent')
    _, y_ventCent, z_ventCent = sph_vent.pos()
    _, y_ventOutf, z_ventOutf = sph_outf.pos()
    if dORv == 'V':
        # azimuth = 90
        orient_ventX = orient_vent.clone().projectOnPlane('x').c('hotpink').x(0).legend('lin_OrientVent(ProjX)')
    elif dORv == 'D': 
        # azimuth = 0
        orient_ventX = orient_vent.clone().projectOnPlane('z').c('hotpink').x(0).legend('lin_OrientVent(ProjX)')
    pts_vent = orientVectors(orient_ventX)
    ang_vent = findAngleBtwVectors(pts_vent, pts_heart)
    vent_txt = 'Ventricular orientation with respect to heart (deg): '+ format(ang_vent,'.1f')
    print('\t>> '+vent_txt)
    
    # Angle between chambers
    ang_chs = findAngleBtwVectors(pts_atr, pts_vent)
    chs_txt = 'Angle between chambers (deg): '+ format(ang_chs,'.1f')
    print('\t>> '+chs_txt)

    lines_orient = [orient_atr, orient_vent, orient_atrX, orient_ventX, linLineX]
    dict_kspl = addKSplines2Dict(kspls = lines_orient, info = '', dict_kspl = dict_kspl)
    df_res = addOrientationAngles2df(df_res = df_res, file_num = file_num, angles = [ang_heart, ang_atr, ang_vent, ang_chs], names = ['ang_Heart', 'ang_Atr', 'ang_Vent', 'ang_BtwChambers'])
    
    text = filename+"\n\n >> Measuring chamber orientation\n - "+atr_txt+"\n - "+vent_txt+"\n - "+chs_txt
    txt = Text2D(text, c="k", font= 'CallingCode')
    vp = Plotter(N=1, axes = 8)
    vp.show(m_myoc.alpha(0.01), m_atrMyoc.alpha(0.01), m_ventMyoc.alpha(0.01), sph_orient, kspl_CL2use, orient_atr, orient_vent, linLine, 
            orient_atrX.alpha(1), orient_ventX.alpha(1), linLineX.alpha(1), txt, azimuth = azimuth, elevation = elevation, zoom = 0.8, at=0, interactive=True) #
    # vp.show()

    return sph_orient, lines_orient, dict_pts, dict_kspl, df_res

#%% func - getDistance2Mesh
def getDistance2Mesh(filename, m_int, m_ext, title, alpha = 0.1, plotshow = True):
    
    thickness_meshes = ['Cardiac Jelly Thickness', 'Myoc.Thickness', 'Endo.Thickness']
    tic = perf_counter()
    print('- Extracting '+title)
    print("  > Start time: \t", str(datetime.now())[11:-7])
    m_ext_out = m_ext.clone()
    m_ext_out.distanceToMesh(m_int, signed=True, negate=False)
    if title in thickness_meshes:
        thickness = - m_ext_out.getPointArray("Distance")
    else:
        thickness = m_ext_out.getPointArray("Distance")
    
    vmin, vmax = np.min(thickness),np.max(thickness)
    toc = perf_counter()
    time = toc-tic
    print("  > End time: \t\t", str(datetime.now())[11:-7])
    print("\t>> Total time taken to get heatmap = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")     

    alert('jump',1)
    
    val_min = vmin; val_max = vmax
        
    title_split = title.split(' ')
    title_bar = title_split[0]
    for txt in title_split[1:]:
        title_bar = title_bar +'\n'+ txt 
        
    m_ext_out.pointColors(thickness, cmap="jet", vmin=val_min, vmax=val_max)
    m_ext_out.mapper().SetScalarRange(val_min,val_max)
    m_ext_out.addScalarBar(title=title_bar+' [um]')
    m_ext_out.mapper().SetScalarRange(val_min,val_max)
    
    m_ext_out.alpha(alpha).legend(m_ext.legend()) 
    m_int.color("white").alpha(1).wireframe()
    
    if plotshow:
        text = filename+"\n\n >>"+title+" Heat map\n\t   - min and max: ("+format(vmin, '.2f')+" , "+format(vmax, '.2f')+")"
        txt = Text2D(text, c="k", font= 'CallingCode')
        
        cube = Cube(pos=m_ext_out.centerOfMass(), side=300, c='white', alpha=0.01)
        vp = Plotter(N=1, axes=10)  
        vp.show(m_ext_out, m_int, cube, txt, at=0, zoom=1.2)
    
    min_max = [vmin, vmax]

    return thickness, m_ext_out, min_max

#%% func - unloopHeart
def unloopHeart(filename, mesh, kspl_CL2use, cl_ribbon, no_planes, pl_CLRibbon, param, plotshow = False, tol=0.05):
    
    # Create matrix with all data
    #0:x, 1:y, 2:z, 3:taken, 4:z_plane, 5:theta, 6: radius, 7: param
    matrix_unlooped = np.zeros((len(mesh.points()),8))
    matrix_unlooped[:,0:3] = mesh.points()
    matrix_unlooped[:,8-1] = param 
    
    # - Get unitary normal of plane to create CL_ribbon
    pl_normCLRibbon = unit_vector(pl_CLRibbon['pl_normal'])
    pl_centCLRibbon = np.asarray(pl_CLRibbon['pl_centre'])
    plane_CLRibbon = Plane(pos = pl_centCLRibbon, normal = pl_normCLRibbon, sx = 300).color('skyblue').alpha(0.5)
    arr_vectNormRib = Arrow(pl_centCLRibbon, pl_centCLRibbon+np.asarray(pl_CLRibbon['pl_normal'])*120, s = 0.1, c = 'cyan')
    
    text = filename+"\n\n >> Is the cyan arrow pointing towards the ventral side of the heart? \n    Check, close the window and answer."    
    txt = Text2D(text, c="k", font= 'CallingCode')
    vp = Plotter(N=1, axes=4)
    vp.show(mesh, kspl_CL2use, cl_ribbon, arr_vectNormRib,txt, at=0, azimuth = 0, interactive=1)

    q_ventralDir = ask4input('Is the cyan arrow pointing towards the ventral side of the heart? [0]:no/[1]:yes : ', int)
    if q_ventralDir == 0:
        pl_normCLRibbon = -120*pl_normCLRibbon
    else: 
        pl_normCLRibbon = 120*pl_normCLRibbon
    
    # Get normals and centres of planes 
    pl_normals, pl_centres = getPlaneNormals (no_planes = no_planes+2, spline_pts = kspl_CL2use.points())
    pl_normals = pl_normals[1:-2]
    pl_centres = pl_centres[1:-2]
    # print(pl_normals)
    
    plane_num = np.linspace(0,1,len(pl_normals))
    
    # Iterate through each plane
    for i, normal, centre in zip(count(), pl_normals, pl_centres):
        # print('normal', normal)
        # A. Get cut plane info 
        # - Info Plane
        d = normal.dot(centre)
        arr_vectPlCut = Arrow(centre, centre+normal*20, s = 0.1, c='orange')
        
        # B. Get vector that defines 0 deg angle
        arr_normCRib = Arrow(centre, centre+pl_normCLRibbon, s = 0.1, c='black')
        # - Get projection of pl_normalCLRibbon on cutting plane
        proj = pl_normCLRibbon - (np.dot(pl_normCLRibbon, normal)/np.linalg.norm(normal)**2)*normal 
        #print(proj)
        
        v_vectC = proj
        # - Unit vector
        v_unitvectC = unit_vector(v_vectC)
        
        # Create stuff to plot 
        pl_cut = Plane(pos = centre, normal = normal, sx = 300).color('coral').alpha(0.5)
        sph_C = Sphere(centre, r=2, c='red')
        sph_VXYZ = Sphere(centre+v_vectC, r=2, c='green')
        arr_vectXYZ = Arrow(centre, centre+v_vectC, s = 0.1)
        
        # C. Get points of mesh at plane
        d_points = np.absolute(np.dot(np.subtract(matrix_unlooped[:,0:3],np.asarray(centre)),np.asarray(normal)))
        # Find the indexes where the d values fit the plane and are not yet taken
        index_ptsAtPlane = np.where((d_points <= tol) & (matrix_unlooped[:,3] == 0))
        
        # Define new matrix just with the points on plane
        new_matrix = matrix_unlooped[index_ptsAtPlane,:][0]
        # - Get points of mesh that are on plane, centered on centreline point
        ptsC = np.subtract(new_matrix[:,0:3],np.asarray(centre))
        # print('ptsc:', ptsC, 'len:', len(ptsC))
        # - Get the radius of those points
        radius = [np.linalg.norm(x) for x in ptsC]
        # print('radius:', radius)
        
        # D. Find direction of point with respect to plane that includes central point, vC and the normal of the cutting plane
        normal_divLR = np.cross(normal, v_unitvectC)
        lORr = np.sign(np.dot(ptsC, np.asarray(normal_divLR)))
        # print(normal_divLR)
        # print('lORr:', lORr,'min:', min(lORr), ' - max:', max(lORr))
        
        # E. Get angle of points in that plane
        av = np.dot(ptsC,v_unitvectC)
        # print('av:', av)
        cosTheta = np.divide(av, radius) 
        # print('cosTheta:', cosTheta)
        theta = np.arccos(cosTheta)*180/np.pi
        
        # print('theta:', theta)
        # print(min(theta), max(theta))
        theta_corr = np.multiply(lORr, theta)
        # print('theta_corr:', theta_corr)
        # print(min(theta_corr), max(theta_corr))
        
        if i % 10 == 0:
            # sphL = Spheres(ptC[::50,:]+centre, r=2, c='blueviolet')
            # sphR = Spheres(ptC[::50,:]+centre, r=2, c='blueviolet')
            
            sphL = []
            sphR = []
            for num, pt in enumerate(ptsC):
                if num % 50 == 0:
                    if lORr[num] == 1:
                        sphL.append(Sphere(pt+centre, r=2, c='blueviolet'))
                    else: 
                        sphR.append(Sphere(pt+centre, r=2, c='gold'))
                    
            vp = Plotter(N=1, axes=4)
            vp.show(mesh,sphL, sphR,kspl_CL2use, cl_ribbon, arr_normCRib, arr_vectXYZ,arr_vectPlCut,arr_vectNormRib,sph_C, sph_VXYZ, pl_cut, at=0, azimuth = 0, interactive=1)

        # - Save all obtained values in matrix_unlooped
        for num, index in enumerate(index_ptsAtPlane[0]):
            #3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: param
            matrix_unlooped[index,3] = 1
            matrix_unlooped[index,4] = plane_num[i]
            matrix_unlooped[index,5] = theta_corr[num]
            matrix_unlooped[index,6] = radius[num]
            matrix_unlooped[index,7] = param[num]
    
    return matrix_unlooped

#%% func - unloopHeart-BU2
# def unloopHeart(mesh, kspl_CL2use, cl_ribbon, no_planes, pl_CLRibbon, param, tol=0.05):
    
#     # Create matrix with all data
#     n_param = 1#param.shape[1]
#     dims = 7+n_param
#     #0:x, 1:y, 2:z, 3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: param
#     matrix_unlooped = np.zeros((len(mesh.points()),dims))
#     matrix_unlooped[:,0:3] = mesh.points()
#     matrix_unlooped[:,dims-1] = param #-n_param:dims] = param
    
#     # Get info cl_ribbon and plane with which it was created
#     matrix_clRibbon = np.asarray(cl_ribbon.points())
#     # - Get unitary normal of plane
#     pl_normCLRibbon = unit_vector(pl_CLRibbon['pl_normal'])
    
#     # print('pl_normCLRibbon ', pl_normCLRibbon)
#     pl_centCLRibbon = np.asarray(pl_CLRibbon['pl_centre'])
#     plane_CLRibbon = Plane(pos = pl_centCLRibbon, normal = pl_normCLRibbon, sx = 300).color('skyblue').alpha(0.5)
#     arr_vectNormRib = Arrow(pl_centCLRibbon, pl_centCLRibbon+np.asarray(pl_CLRibbon['pl_normal'])*20, s = 0.1, c = 'cyan')
#     pos_pt = pl_normCLRibbon*50+pl_centCLRibbon#pl_centCLRibbon+200
#     sph_pos = Sphere(pos_pt, r=2, c='orangered')
#     # - d value of clr_plane
#     # d_clr = np.dot(pl_normCLRibbon, pl_centCLRibbon)
#     d_clr_pve = np.dot(pos_pt,pl_normCLRibbon)
#     print('d_clr_pve:', d_clr_pve)
#     if d_clr_pve > 0:
#         criteria = 'positive'
#     else:
#         criteria = 'negative'
    
    
#     # Get normals and centres of planes 
#     pl_normals, pl_centres = getPlaneNormals (no_planes = no_planes+2, spline_pts = kspl_CL2use.points())
#     pl_normals = pl_normals[1:-2]
#     pl_centres = pl_centres[1:-2]
#     # print(pl_normals)
    
#     plane_num = np.linspace(0,1,len(pl_normals))
    
#     # Iterate through each plane
#     for i, normal, centre in zip(count(), pl_normals, pl_centres):
#         # print('normal', normal)
#         # A. Get cut plane info 
#         # - Info Plane
#         d = normal.dot(centre)
#         arr_vectPlCut = Arrow(centre, centre+normal*20, s = 0.1, c='orange')
        
#         if i % 10 == 0:
#             pl_cut = Plane(pos = centre, normal = normal, sx = 300).color('coral').alpha(0.5)
#             sph_C = Sphere(centre, r=2, c='red')
#             vp = Plotter(N=1, axes=4)
#             vp.show(mesh,sph_pos,kspl_CL2use, cl_ribbon,arr_vectPlCut,arr_vectNormRib,pl_cut, sph_C, at=0, azimuth = 0, interactive=1)

        
#         # B. Get vector that defines 0 deg angle
#         # - Get points of cl_ribbon that intersect plane to cut
#         d_clRpts = np.dot(np.subtract(matrix_clRibbon,np.asarray(centre)),np.asarray(normal))
#         # - Get d values of cl_ribbon with respect to plane of cl extension
#         d_clRpts2 = np.dot(np.subtract(matrix_clRibbon,np.asarray(pl_centCLRibbon)),pl_normCLRibbon)
#         # - Find the indexes of the cl_ribbon points where the d values fit the plane to cut and have positive d_clr values
#         if criteria == 'positive':
#             index_clRAtPlanes = np.where((np.absolute(d_clRpts) <= 0.5)  & (d_clRpts2 > 0))
#         else: 
#             index_clRAtPlanes = np.where((np.absolute(d_clRpts) <= 0.5)  & (d_clRpts2 < 0))
#         # print(index_clRAtPlanes)   
#         cut_matrix = matrix_clRibbon[index_clRAtPlanes]
#         d_clRptsCut = np.dot(np.subtract(cut_matrix,np.asarray(pl_centCLRibbon)),pl_normCLRibbon)
#         # print(d_clRptsCut)
#         if criteria == 'positive':
#             index_v = np.where(d_clRpts2[index_clRAtPlanes] == min(d_clRpts2[index_clRAtPlanes]))
#         else: 
#             index_v = np.where(d_clRpts2[index_clRAtPlanes] == max(d_clRpts2[index_clRAtPlanes]))
#         # - Point that defines the vector with which 0 deg is defined
#         v_vectXYZ = cut_matrix[index_v]
#         # print('v_vectXYZ:', v_vectXYZ)
#         # print('d_vXYZ: ', d_clRptsCut[index_v])
#         # - Coordinates of vector with respect to centre point in cl
#         v_vectC = v_vectXYZ[0]-centre
#         v_unitvectC = unit_vector(v_vectC)
        
#         # Create stuff to plot 
#         pl_cut = Plane(pos = centre, normal = normal, sx = 300).color('coral').alpha(0.5)
#         sph_C = Sphere(centre, r=2, c='red')
#         sph_VXYZ = Sphere(v_vectXYZ[0], r=2, c='green')
#         arr_vectXYZ = Arrow(centre, centre+v_vectC, s = 0.1)
        
#         # if i % 10 == 0:
#         #     vp = Plotter(N=1, axes=4)
#         #     vp.show(mesh, kspl_CL2use, cl_ribbon, arr_vectXYZ,arr_vectPlCut,arr_vectNormRib,sph_C, sph_VXYZ, pl_cut, at=0, azimuth = 0, interactive=1)

#         # C. Get points of mesh at plane
#         d_points = np.absolute(np.dot(np.subtract(matrix_unlooped[:,0:3],np.asarray(centre)),np.asarray(normal)))
#         # Find the indexes where the d values fit the plane and are not yet taken
#         index_ptsAtPlane = np.where((d_points <= tol) & (matrix_unlooped[:,3] == 0))
        
#         # Define new matrix just with the points on plane
#         new_matrix = matrix_unlooped[index_ptsAtPlane,:][0]
#         # - Get points of mesh that are on plane, centered on centreline point
#         ptsC = np.subtract(new_matrix[:,0:3],np.asarray(centre))
#         # print('ptsc:', ptsC, 'len:', len(ptsC))
#         # - Get the radius of those points
#         radius = [np.linalg.norm(x) for x in ptsC]
#         # print('radius:', radius)
        
#         # D. Find direction of point with respect to plane that includes central point, vC and the normal of the cutting plane
#         normal_divLR = np.cross(normal, v_unitvectC)
#         lORr = np.sign(np.dot(ptsC, np.asarray(normal_divLR)))
#         # print(normal_divLR)
#         # print('lORr:', lORr,'min:', min(lORr), ' - max:', max(lORr))
        
#         # E. Get angle of points in that plane
#         av = np.dot(ptsC,v_unitvectC)
#         # print('av:', av)
#         cosTheta = np.divide(av, radius) 
#         # print('cosTheta:', cosTheta)
#         theta = np.arccos(cosTheta)*180/np.pi
        
#         # print('theta:', theta)
#         # print(min(theta), max(theta))
#         theta_corr = np.multiply(lORr, theta)
#         # print('theta_corr:', theta_corr)
#         # print(min(theta_corr), max(theta_corr))
        
#         if i % 10 == 0:
#             sphL = []
#             sphR = []
#             for num, pt in enumerate(ptsC):
#                 if num % 200 == 0:
#                     if lORr[num] == 1:
#                         sphL.append(Sphere(pt+centre, r=2, c='red'))
#                     else: 
#                         sphR.append(Sphere(pt+centre, r=2, c='navy'))
                    
#             vp = Plotter(N=1, axes=4)
#             vp.show(mesh,sphL, sphR, sph_pos,kspl_CL2use, cl_ribbon, arr_vectXYZ,arr_vectPlCut,arr_vectNormRib,sph_C, sph_VXYZ, pl_cut, at=0, azimuth = 0, interactive=1)

#         # - Save all obtained values in matrix_unlooped
#         for num, index in enumerate(index_ptsAtPlane[0]):
#             #3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: param
#             matrix_unlooped[index,3] = 1
#             matrix_unlooped[index,4] = plane_num[i]
#             matrix_unlooped[index,5] = theta_corr[num]
#             matrix_unlooped[index,6] = radius[num]
#             matrix_unlooped[index,dims-1] = param[num]
    
#     return matrix_unlooped

#%% func - unloopHeart-BU
def unloopHeart2(mesh, kspl_CL2use, cl_ribbon, no_planes, pl_CLRibbon, param, tol=0.05):
    
    # Create matrix with all data
    n_param = 1#param.shape[1]
    dims = 7+n_param
    #0:x, 1:y, 2:z, 3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: param
    matrix_unlooped = np.zeros((len(mesh.points()),dims))
    matrix_unlooped[:,0:3] = mesh.points()
    matrix_unlooped[:,dims-1] = param #-n_param:dims] = param
    
    # Get info cl_ribbon and plane with which it was created
    matrix_clRibbon = np.asarray(cl_ribbon.points())
    # - Get unitary normal of plane
    pl_normCLRibbon = unit_vector(pl_CLRibbon['pl_normal'])
    
    print('pl_normCLRibbon ', pl_normCLRibbon)
    pl_centCLRibbon = np.asarray(pl_CLRibbon['pl_centre'])
    plane_CLRibbon = Plane(pos = pl_centCLRibbon, normal = pl_normCLRibbon, sx = 300).color('skyblue').alpha(0.5)
    arr_vectNormRib = Arrow(pl_centCLRibbon, pl_centCLRibbon+np.asarray(pl_CLRibbon['pl_normal'])*20, s = 0.1, c = 'cyan')
    pos_pt = pl_centCLRibbon+200
    sph_pos = Sphere(pos_pt, r=2, c='orangered')
    # - d value of clr_plane
    d_clr = np.dot(pl_normCLRibbon, pl_centCLRibbon)
    d_clr_pve = np.dot(pl_normCLRibbon,pos_pt)-d_clr
    print('d_clr_pve:', d_clr_pve)
    if d_clr_pve > 0:
        criteria = 'positive'
    else:
        criteria = 'negative'
    
    # Get normals and centres of planes 
    pl_normals, pl_centres = getPlaneNormals (no_planes = no_planes+2, spline_pts = kspl_CL2use.points())
    pl_normals = pl_normals[1:-2]
    pl_centres = pl_centres[1:-2]
    # print(pl_normals)
    
    plane_num = np.linspace(0,1,no_planes)
    
    # Iterate through each plane
    for i, normal, centre in zip(count(), pl_normals, pl_centres):
        print('normal', normal)
        # A. Get cut plane info 
        # - Info Plane
        d = normal.dot(centre)
        arr_vectPlCut = Arrow(centre, centre+normal*20, s = 0.1, c='orange')
        
        # B. Get vector that defines 0 deg angle
        # - Get points of cl_ribbon that intersect plane to cut
        d_clRpts = np.dot(matrix_clRibbon,normal)-d
        # - Get d values of cl_ribbon with respect to plane of cl extension
        d_clRpts2 = np.dot(matrix_clRibbon,pl_centCLRibbon)-d_clr
        # - Find the indexes of the cl_ribbon points where the d values fit the plane to cut and have positive d_clr values
        if criteria == 'positive':
            index_clRAtPlanes = np.where((np.absolute(d_clRpts) <= 0.8)  & (d_clRpts2 > 0))
        else: 
            index_clRAtPlanes = np.where((np.absolute(d_clRpts) <= 0.8)  & (d_clRpts2 < 0))
        print('len:', len(index_clRAtPlanes))   
        cut_matrix = matrix_clRibbon[index_clRAtPlanes]
        index_v = np.where(d_clRpts2[index_clRAtPlanes] == min(d_clRpts2[index_clRAtPlanes]))
        # - Point that defines the vector with which 0 deg is defined
        v_vectXYZ = cut_matrix[index_v]
        print('v_vectXYZ:', v_vectXYZ)
        print('d_vXYZ: ', d_clRpts2[index_v])
        # - Coordinates of vector with respect to centre point in cl
        v_vectC = v_vectXYZ[0] - centre
        v_unitvectC = unit_vector(v_vectC)
        
        # Create stuff to plot 
        pl_cut = Plane(pos = centre, normal = normal, sx = 300).color('coral').alpha(0.5)
        sph_C = Sphere(centre, r=2, c='red')
        sph_VXYZ = Sphere(v_vectXYZ[0], r=2, c='green')
        arr_vectXYZ = Arrow(centre, centre+v_vectC, s = 0.1)
        # if i == 0:
        vp = Plotter(N=1, axes=1)#8)
        vp.show(mesh, sph_pos, kspl_CL2use, cl_ribbon, arr_vectXYZ,arr_vectPlCut,arr_vectNormRib,sph_C, sph_VXYZ, pl_cut, plane_CLRibbon, at=0, azimuth = 0, interactive=1)

        # C. Get points of mesh at plane
        d_points = np.absolute(np.dot(matrix_unlooped[:,0:3],normal)-d)
        # Find the indexes where the d values fit the plane and are not yet taken
        index_ptsAtPlane = np.where((d_points <= tol) & (matrix_unlooped[:,3] == 0))
        
        # Define new matrix just with the points on plane
        new_matrix = matrix_unlooped[index_ptsAtPlane,:][0]
        # - Get points of mesh that are on plane, centered on centreline point
        ptsC = np.subtract(new_matrix[:,0:3],np.asarray(centre))
        # - Get the radius of those points
        radius = [np.linalg.norm(x) for x in ptsC]
        # - Get angle
        av = np.dot(ptsC,v_unitvectC)
        cosTheta = np.divide(av, radius) 
        theta = np.arccos(cosTheta)*180/np.pi
        
        # - Save all obtained values in matrix_unlooped
        for num, index in enumerate(index_ptsAtPlane[0]):
            #3:taken, 4:z_plane, 5:theta, 6: radius, 7-8: param
            matrix_unlooped[index,3] = 1
            matrix_unlooped[index,4] = plane_num[i]
            matrix_unlooped[index,5] = theta[num]
            matrix_unlooped[index,6] = radius[num]
            matrix_unlooped[index,dims-1] = param[num]
    
    return matrix_unlooped

#%% - PROCESS DATA
#%% func - newNormal
def newNormal (normal, rotX):
    alpha = np.radians(sum(rotX))
  
    new_normal = normal
    if alpha != 0:
        #Define rotation matrix
        r = np.array(((np.cos(alpha), -np.sin(alpha)),
               (np.sin(alpha),  np.cos(alpha))))
        #Define vector yz
        v = np.array((new_normal[1],new_normal[2]))
        #Apply rotation matrix
        new_v = r.dot(v)
        
        new_normal =  np.asarray((normal[0], new_v[0], new_v[1]))
        
    return new_normal

#%% func - rotation_matrix
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    
    e.g.I have two vectors as Python lists and an angle:
        v = [3,5,0]
        axis = [4,4,1]
        theta = 1.2 #radian
    What is the best/easiest way to get the resulting vector when rotating 
    the v vector around the axis?
    
    print(np.dot(rotation_matrix(axis, theta), v)) 
    # [ 2.74911638  4.77180932  1.91629719]
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
    #mag = math.sqrt(sum(i**2 for i in x))
    sqrs = [i**2 for i in v]
    mag = math.sqrt(sum(sqrs))
    v_unit = [j/mag for j in v]
        
    return np.asarray(v_unit)

#%% func - newNormal3DRot
def newNormal3DRot (normal, rotX, rotY, rotZ):
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
    min_dist = 9999999
    for n in range(len(pts)):
        pt = pts[n]
        squared_dist = np.sum((pt-pt_o)**2, axis=0)
        dist = np.sqrt(squared_dist)
        
        if dist < min_dist:
            min_dist = dist
            pt_out = pt
            num_pt = n
        
    return pt_out, num_pt

#%% func - findDist
def findDist(pt1, pt2):
    squared_dist = np.sum((pt1-pt2)**2, axis=0)
    dist = np.sqrt(squared_dist)
    
    return dist

#%% func - getPointsAtPlane 
# Function to get points within mesh at certain heights (y positions) to create kspline
def getPointsAtPlane (points, pl_normal, pl_centre, tol=2):
    
    pts_cut = []
        
    d = pl_normal.dot(pl_centre)
    #print('d for tol:', d)
    d_range = [d-tol, d+tol]
    #print('d range:', d_range)
    d_range.sort()
    
    for pt in points:
        d_pt = pl_normal.dot(pt)
        if d_pt>d_range[0] and d_pt<d_range[1]:
            pts_cut.append(pt)
            #print(pt)
    
    pts_cut = np.asarray(pts_cut)
    
    return pts_cut

#%% func - getPointsAtPlaneAddParam 
# Function to get points within mesh at certain heights (y positions) to create kspline
# def getPointsAtPlaneAddParam (points, pl_normal, pl_centre, taken, assigned_plane, tol=0.05):
    
#     d = pl_normal.dot(pl_centre)
#     #print('d for tol:', d)
#     d_range = [d-tol, d+tol]
#     #print('d range:', d_range)
#     d_range.sort()
    
#     d_points = np.dot(points,pl_normal)
#     index_ptsAtPlane = np.where((d_points >= -tol) & (d_points <= tol))[0]
    
#     for index in index_ptsAtPlane:
#         if taken[index] != 1:
            
            
#     for num, pt in enumerate(points):
#         d_pt = pl_normal.dot(pt)
#         if d_pt>d_range[0] and d_pt<d_range[1]:
#             isTaken = taken[num]
#             if isTaken == 0:
#                 param_cut.append(param[num])
#                 pts_cut.append(pt)
#                 taken[num] = 1
    
#     pts_cut = np.asarray(pts_cut)
#     param_cut = np.asarray(param_cut)
#     taken = np.asarray(taken)
    
#     return pts_cut, param_cut, taken

#%% func - order_pts
def order_pts (points):
    
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

#%% func - order_ptsAddParam
def order_ptsAddParam (points, center_pt, param, taken):
    
    cent_pts = np.zeros_like(points)
    for num in range(len(points)):
        cent_pts[num][0]=points[num][0]-center_pt[0]
        cent_pts[num][2]=points[num][2]-center_pt[2]
    
    angle_deg = np.zeros((len(points),1))
    for num, pt in enumerate(cent_pts):
        angle_deg[num] = np.arctan2(pt[2],pt[0])*(180/np.pi)
    
    index_sort = np.argsort(angle_deg, axis=0)
    ordered_pts = np.zeros_like(points)
    ordered_param = np.zeros_like(param)
    ordered_angle = np.zeros_like(angle_deg)
    #ordered_taken = np.zeros_like(taken)
    for i, pos in enumerate(index_sort):
        ordered_pts[i]=points[pos]
        ordered_angle[i] = angle_deg[pos]
        ordered_param[i] = param[pos]
        #ordered_taken[i] = taken[pos] # keep the info of taken this far?
    
    return ordered_pts, ordered_angle, ordered_param#, ordered_taken

#%% func - getInterpolatedPts
def getInterpolatedPts(points, nPoints):
   
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

#%% func - findAngleBtwVectors
def findAngleBtwVectors(pts1, pts2):
    
    mag_v1 = findDist(pts1[0],pts1[1])
    mag_v2 = findDist(pts2[0],pts2[1])
    
    vect1 = pts1[1]-pts1[0]
    vect2 = pts2[1]-pts2[0]
    
    dotProd = np.dot(vect1,vect2)
    
    angle = math.degrees(math.acos(dotProd/(mag_v1*mag_v2)))
    
    return angle

#%% func - orientVectors
def orientVectors(line):
    if line.legend()  == 'lin_OrientAtr(ProjX)':
        pts_line = line.points()
        pts_line = pts_line[pts_line[:,2].argsort()]
    elif line.legend() == 'lin_OrientVent(ProjX)':
        pts_line = line.points()
        pts_line = pts_line[pts_line[:,2].argsort()[::-1]]
    else: # Linear line
        pts_line = line.points()
        pts_line = pts_line[pts_line[:,2].argsort()]
    
    #print(np.diff(pts_line, axis = 0))
    
    return  pts_line

#%% func - classifyPtsMx
def classifyPtsMx(dict_planes, pl_name, pts_whole):
    """
    

    Parameters
    ----------
    dict_planes : TYPE
        DESCRIPTION.
    pl_name : TYPE
        DESCRIPTION.
    pts_whole : TYPE
        DESCRIPTION.

    Returns
    -------
    pts_classFinal : TYPE
        DESCRIPTION.

    """
    
    if pl_name == 'pl2CutMesh_Chamber':
        ptA = 'atrium'
        ptB = 'ventricle'
        name = 'Atrium-Ventricle'
    elif pl_name == 'pl_AtrCoronal':
        ptA = 'dorsal'
        ptB = 'ventral'
        name = 'Dorsal-Ventral (Atr)'
    elif pl_name == 'pl_VentCoronal':
        ptA = 'ventral'
        ptB = 'dorsal'
        name = 'Dorsal-Ventral (Vent)'
        
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
def classifyHeartPts(filename, mesh, dict_planes, pts_whole, pts_left, data, names_data, plot_show = True):
    
    print('- Classifying points for ', names_data)
    tic = perf_counter()
    
    cols = ['AtrVent','DorsVent','LeftRight'] + names_data
    df_classPts = pd.DataFrame(columns=cols)
    
    # Classify AnV
    pts_classAnV = classifyPtsMx(dict_planes = dict_planes, pl_name = 'pl2CutMesh_Chamber', pts_whole = pts_whole)
    df_classPts['AtrVent'] = pts_classAnV
    # Classify DnV
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
    
    # Classify LnR
    print('--> Left-Right')
    av = pts_whole.view([('', pts_whole.dtype)] * pts_whole.shape[1]).ravel()
    bv = pts_left.view([('', pts_left.dtype)] * pts_left.shape[1]).ravel()
    # cint = np.intersect1d(av, bv).view(a.dtype).reshape(-1, a.shape[1])
    c_isin = np.isin(av,bv)
    index_isin = np.where(c_isin == True)[0]
    
    pts_classLnR = np.empty(len(pts_whole), dtype='object')
    pts_classLnR[:] = 'right'
    pts_classLnR[index_isin] = 'left'
    
    df_classPts['LeftRight'] = pts_classLnR
    
    for i, name, dat in zip(count(), names_data, data):
        df_classPts[name] = dat
        
    toc = perf_counter()
    time = toc-tic
    print("- All Done - points have been classified!\n")
    print(df_classPts.sample(20))
    print("- Time taken to classify = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")     
    alert('whistle',1)
    
    if plot_show:
        plotPtClassif(filename, mesh, pts_whole, [pts_classAnV, pts_classDnV, pts_classLnR])
    
    return df_classPts
    
#%% func - intersecOfSets
def intersecOfSets(arr1, arr2, arr3): 
    # https://www.geeksforgeeks.org/python-program-find-common-elements-three-lists-using-sets/
    # Converting the arrays into sets 
    s1 = set(arr1) 
    s2 = set(arr2) 
    s3 = set(arr3) 
      
    # Calculates intersection of sets on s1 and s2 
    set1 = s1.intersection(s2) 
    # Calculates intersection of sets on set1 and s3 
    result_set = set1.intersection(s3)
    # Converts resulting set to list 
    final_list = list(result_set) 
    #print(final_list) 
    
    return final_list
        
#%% func - classifyPtAtSidesofPlane
def classifyPtAtSidesofPlane(pt, d, normal, classes_sorted):
    x,y,z = pt
    a,b,c = normal
    
    dotProd_pt = dot((a,b,c), (x,y,z))
    #print("dotProd_pt:", dotProd_pt)
    if dotProd_pt < d:
        pts_class = classes_sorted[0]
    else:
        pts_class = classes_sorted[-1]
    #print("Point classified as:", pts_class)
    
    return pts_class

#%% func - getPlaneNormals
def getPlaneNormals (no_planes, spline_pts):
        
    normals = []
    pt_centre = []
    every = len(spline_pts)//no_planes
    for i in range(len(spline_pts)-1,0,-every):
        pt_centre.append(spline_pts[i])
        normal = spline_pts[i-1]-spline_pts[i]
        normals.append(normal)

    return normals, pt_centre

#%% - SAVING  
#%% func - saveMesh
def saveMesh(filename, mesh, mesh_name, dir_stl, extension):
        
    mesh_title = filename+"_"+mesh_name+"."+extension
    mesh_dir = os.path.join(dir_stl, mesh_title)
    mesh.write(mesh_dir)
    print("- Saved mesh - "+mesh_title+"!")
    alert("countdown",1)
    
#%% func - saveMeshes 
def saveMeshes(filename, meshes, names, dict_colour, dir_stl, extension):
    
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
    
    for i, arrayTh, name in zip(count(), arrays2save, names):
        thickness_title = filename+"_"+name
        thick_cj_dir = os.path.join(dir2save,thickness_title)
        np.save(thick_cj_dir, arrayTh)
        
    print("- All arrays (thickness/ballooning) have been saved!")
    alert("countdown",1)
    
#%% func - saveVideo 
def saveVideo (filename, info, meshes4video, rotAngle, dir2save, plotshow=True, duration = 15):
    
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
    
    if plotshow:
        text = filename+"\n\n >> Rotated meshes"; txt = Text2D(text, c="k", font= 'CallingCode')
        vp1 = Plotter(N=1, axes=8)      
        vp1.show(rotMeshes4video, txt, at=0, axes=8, zoom=1.2)

    text2 = filename; txt2 = Text2D(text2, c="k", font= 'CallingCode')
    
    vp = Plotter(bg='white', axes=10, offscreen=True)
    vp.show(rotMeshes4video, txt2, zoom=1.4)
    video_name = os.path.join(dir2save, filename+"_"+info+".mp4")
    video = Video(video_name, duration=duration, backend='opencv')
    for i in range(180):
        vp.show(elevation=0, azimuth=2)  # render the scene
        video.addFrame()   
    video.close()
            
    for k, mesh in enumerate(rotMeshes4video):
        if dORv == 'V':
            rot_mesh = mesh.rotateX(-rotAngle)
        elif dORv == 'D':
            rot_mesh = mesh.rotateZ(-rotAngle)
            rot_mesh = rot_mesh.rotateY(-90)

    alert('wohoo',1) 

#%% func- saveMultVideos
def saveMultVideos(filename, info, meshes4video, rotAngle, dir2save, plotshow):
    
    print('Saving videos... this might take a while... (about 2min/video')
    bar = Bar('Saving' , max = len(info), suffix = suffix, check_tty=False, hide_cursor=False)
    for i, name, mesh in zip(count(), info, meshes4video):
        saveVideo(filename = filename, info = name, meshes4video = [mesh], 
                        rotAngle = rotAngle, dir2save = dir2save, plotshow = plotshow)
        bar.next()
    bar.finish()
    alert('whistle',1)
        
#%% - PLOTTERS
#%% func - plotPtClassif
def plotPtClassif(filename, mesh, pts_whole, pts_class):
    
    [pts_classAnV, pts_classDnV, pts_classLnR] = pts_class
    names_class = [['atrium', 'ventricle'],['dorsal', 'ventral'],['left', 'right']]
    color = [['tomato','gold'], ['greenyellow', 'darkgreen'],['deepskyblue', 'darkblue']]
    text = filename+'\n\n >> Point classification'
    txt = Text2D(text, c="k", font= 'CallingCode')
    
    vp = Plotter(shape = (1, 4), axes = 7)
    vp.show(mesh, txt, at = 0)
    for i, sp_class in enumerate(pts_class):
    
        ind_first = np.where(sp_class == names_class[i][0])[0]
        pts_first = pts_whole[ind_first]
        sph_first = Spheres(pts_first[0::50,:], c = color[i][0], r=1).legend(names_class[i][0])
        
        
        ind_second = np.where(sp_class == names_class[i][1])[0]
        pts_second = pts_whole[ind_second]
        sph_second = Spheres(pts_second[0::50,:], c = color[i][1], r=1).legend(names_class[i][1])
        
        if i != 2: 
            vp.show(mesh, sph_first, sph_second, at = i+1)
        else: 
            vp.show(mesh, sph_first, sph_second, at = i+1, interactive = True)
            
#%% func - plotterIndivAlpha
def plotterIndivAlpha(filename, N, meshes, alpha):
    
    vp = Plotter(N=N, axes=7)
    
    text = str(filename)
    txt = Text2D(text, c=c, font=font)
    for num, mesh in enumerate(meshes):
        #If object is a mesh
        if not isinstance(mesh, list):
            mesh.alpha(alpha[num])
        #If object is a list of meshes
        else: 
            for num_i, mesh_i in enumerate(mesh):
                mesh_i.alpha(alpha[num][num_i])
        if num == 0:
            vp.show(mesh, txt, at=num)
        elif num != N-1 or num !=0:
            vp.show(mesh, at=num)
        else:
            vp.show(mesh, at=num, interactive=True)

#%% func - plotterIndiv
def plotterIndiv(filename, N, meshes):
    
    vp = Plotter(N=N, axes=7)
    text = str(filename)
    txt = Text2D(text, c=c, font=font)
    for num, mesh in enumerate(meshes):
        if num == 0:
            vp.show(mesh, txt, at=num)
        elif num != N-1 or num !=0:
            vp.show(mesh, at=num)
        else:
            vp.show(mesh, at=num, interactive=True)
            
#%% func - plotterAll
def plotterAll(filename, meshes, alpha):
    
    text = str(filename)
    txt = Text2D(text, c=c, font=font)
    
    meshes_alpha = []
    vp = Plotter(N=1, axes=7)
    for num, mesh in enumerate(meshes):
        mesh.alpha(alpha[num])
        meshes_alpha.append(mesh)
        
    vp.show(meshes_alpha, txt, at=0, interactive=True)

#%% - PRINT
#%% func - code4vmtkCL
def code4vmtkCL(filename, mesh_name, dir_cl, printshow):

    mesh_titleA = filename+"_"+mesh_name[0]+".stl"
    mesh_titleB = filename+"_"+mesh_name[1]+".stl"
    # mesh_dirA = '"'+os.path.join(dir_cl, mesh_titleA)+'"'
    # mesh_dirB = '"'+os.path.join(dir_cl, mesh_titleB)+'"'
    
    meshML_titleA = filename+"_"+mesh_name[0]+"_cut4clML.stl"
    meshML_titleB = filename+"_"+mesh_name[1]+"_cut4clML.stl"
    meshML_dirA = '"'+os.path.join(dir_cl, meshML_titleA)+'"'
    meshML_dirB = '"'+os.path.join(dir_cl, meshML_titleB)+'"'
    
    cl_titleA = filename+"_"+mesh_name[0]+"_cl.vtp"
    cl_titleB = filename+"_"+mesh_name[1]+"_cl.vtp"
    cl_dirA = '"'+os.path.join(dir_cl, cl_titleA)+'"'
    cl_dirB = '"'+os.path.join(dir_cl, cl_titleB)+'"'
    
    if printshow: 
        print("You are done in python now... to get the centreline with each of the meshes follow the next steps:")
        print(">>> 1. Open the files:  -", mesh_titleA,", ", mesh_titleB+" - in Meshlab")
        print(">>> 2. Run Filters > Remeshing, Simpl.. > Screened Poisson Surf Reco (check Pre-clean)")
        print(">>> 3. Cut outflow tract if needed and export the resulting surface adding '_ML' after _cut4cl in the same folder")
        print(">>> 4. Open VMTK, copy the next text and run it...")
        vmtktxtA = "vmtksurfacereader -ifile "+ meshML_dirA +" --pipe vmtksurfacesmoothing -passband 0.1 -iterations 30 --pipe vmtkcenterlines -seedselector openprofiles -ofile"+ cl_dirA+ " --pipe vmtkrenderer --pipe vmtksurfaceviewer -opacity 0.25 --pipe vmtksurfaceviewer -i @vmtkcenterlines.o -array MaximumInscribedSphereRadius"
        vmtktxtB = "vmtksurfacereader -ifile "+ meshML_dirB +" --pipe vmtksurfacesmoothing -passband 0.1 -iterations 30 --pipe vmtkcenterlines -seedselector openprofiles -ofile"+ cl_dirB+ " --pipe vmtkrenderer --pipe vmtksurfaceviewer -opacity 0.25 --pipe vmtksurfaceviewer -i @vmtkcenterlines.o -array MaximumInscribedSphereRadius"
        print(str(vmtktxtA), '\n\n\n', str(vmtktxtB))
    
    alert("wohoo",1)
    
    return cl_dirA, cl_dirB

#%% Others (back-up)
# import numpy as np
# def random_arr(low, high, size, rand_type):
#     if rand_type == float:
#         arr = [np.random.uniform(low, high) for _ in range(size)]
#     elif rand_type == int:
#         arr = [np.random.randint(low, high) for _ in range(size)]
#     return  arr
    
#%% - ALERT WHEN IMPORTED
print ("IMPORTED: morphoHeart_funcMeshes")
alert('jump',1)