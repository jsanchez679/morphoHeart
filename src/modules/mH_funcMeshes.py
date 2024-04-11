'''
morphoHeart_funcMeshes

@author: Juliana Sanchez-Posada
'''
#%% ##### - Imports - ########################################################
import os
from datetime import datetime
from pathlib import Path, WindowsPath, PurePath
import vedo
#print('vedo:', vedo.__version__)
import numpy as np
import math
# from textwrap import wrap
import flatdict
from itertools import count
import json
import vmtk
from vmtk import pypes, vmtkscripts
from scipy.interpolate import splprep, splev, interpn
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation as R
from time import perf_counter
import copy
from typing import Union
from skimage import measure, io
import seaborn as sns
import random
import pandas as pd
import pickle as pl
import vtk
#Following: https://stackoverflow.com/questions/23573707/disable-or-catch-vtk-warnings-in-vtkoutputwindow-when-embedding-mayavi
vtk.vtkObject.GlobalWarningDisplayOff()
import matplotlib.pyplot as plt

#%% morphoHeart Imports - ##################################################
from ..gui.config import mH_config
from .mH_funcBasics import get_by_path, alert, df_reset_index, df_add_value, palette_rbg

# path_fcMeshes = os.path.abspath(__file__)
path_mHImages = mH_config.path_mHImages
# Load logo
path_logo = path_mHImages / 'logo-07.jpg'
# path_bkgd = path_mHImages / 'white_background400.png'
logo = vedo.Picture(str(path_logo))

#%% Set default fonts and sizes for plots
txt_font = mH_config.txt_font
leg_font = mH_config.leg_font
leg_width = mH_config.leg_width
leg_height = mH_config.leg_height
txt_size = mH_config.txt_size
txt_color = mH_config.txt_color
txt_slider_size = mH_config.txt_slider_size

#%% - morphoHeart Classes and Functions to Work with Meshes
class NumpyArrayEncoder(json.JSONEncoder):
    """
    Definition of class to save dictionary
    Args:
        json (_type_): _description_
    """
    def default(self, obj):
        if isinstance(obj, int):
            return int(obj)
        elif isinstance(obj, float):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, Union(WindowsPath, PurePath)):
            return str(obj)
        else:
            return super(NumpyArrayEncoder, self).default(obj)
        
#Functions linked to GUI
def s32Meshes(organ, gui_keep_largest:dict, win, rotateZ_90=True):#

    proc_ms_all = ['MeshesProc','A-Create3DMesh', 'Status']
    run = False
    for ch in organ.obj_imChannels.keys(): 
        win.win_msg('Creating Channel '+ch[-1]+' meshes!')
        #Check if all the meshes for each channel have been created and set new_set accordingly
        im_ch = organ.obj_imChannels[ch]
        proc_im_all = ['ImProc',im_ch.channel_no,'D-S3Create', 'Status']
        new_set = True

        win.prog_bar_range(0,3)
        aa = 1
        for cont in ['int', 'tiss', 'ext']:
            win.win_msg('Creating meshes of Channel '+im_ch.channel_no[-1]+'! ('+str(aa)+'/3)')
            proc_cont_im = ['ImProc',im_ch.channel_no,'D-S3Create','Info', cont, 'Status']
            proc_cont_ms = ['MeshesProc','A-Create3DMesh', im_ch.channel_no, cont, 'Status']
            mesh_created = im_ch.s32Meshes(cont_type=cont,
                                            keep_largest=gui_keep_largest[im_ch.channel_no][cont],
                                            rotateZ_90 = rotateZ_90, new_set = new_set)
            if mesh_created: 
                # Update organ workflow
                organ.update_mHworkflow(proc_cont_im, 'DONE')
                organ.update_mHworkflow(proc_cont_ms, 'DONE')
                #Update progress in main_win
                win.update_workflow_progress()
                aa+=1
                win.prog_bar_update(aa)
                run = True
            else: 
                if win.already_closed_s3s.isChecked():
                    win.win_msg('*Remember to update the meshes within the s3_numpy folder to the ones in which contours have been selected!', win.keeplargest_play)
                else: 
                    win.win_msg('*Something went wrong creating the meshes for '+ch+'!. Make sure all the Channel s3s have been correctly saved to be able to create the meshes.', win.keeplargest_play)
                return
            
        if run: 
            #Update organ workflow
            organ.update_mHworkflow(proc_im_all, 'DONE')      
            #Message User
            win.win_msg(' Channel '+ch[-1]+' meshes were successfully created!')
            #Enable button for plot
            getattr(win, 'keeplargest_plot_'+ch).setEnabled(True)
            getattr(win, 'summary_whole_plot_'+ch).setEnabled(True)
        else: 
            win.win_msg('*fcM.s32Meshes!')
            alert('error_beep')
    if run: 
        # Update organ workflow
        organ.update_mHworkflow(proc_ms_all, 'Initialised')
    else: 
        win.win_msg('*fcM.s32Meshes!')
        alert('error_beep')

def clean_ch(organ, gui_clean, win, plot_settings=(False,None)):#

    workflow = organ.workflow['morphoHeart']
    
    keys = [key for key in gui_clean.keys() if 'ch' in key]
    for ch in keys:
        win.win_msg('Cleaning Channel '+ch[-1]+' meshes!')
        #Get channel and contours to clean
        ch_to_clean = organ.obj_imChannels[ch]
        ch_to_clean.load_chS3s(cont_types = gui_clean[ch]['cont'])

        #Get channel to use as mask for cleaning
        with_ch = gui_clean[ch]['with_ch']
        ch_with = organ.obj_imChannels[with_ch]
        with_cont = gui_clean[ch]['with_cont']
        ch_with.load_chS3s(cont_types = [with_cont])
        s3_mask = getattr(ch_with, 's3_'+with_cont)
        inverted = gui_clean[ch]['inverted']
        process = ['ImProc', ch, 'E-CleanCh', 'Status']
        
        win.prog_bar_range(0,len(gui_clean[ch]['cont']))
        aa = 0
        for cont in gui_clean[ch]['cont']:
            win.win_msg('Cleaning '+ch+'-'+cont+' with '+with_ch+'-'+with_cont+' (inverted: '+str(inverted)+').')
            #Get the contour to clean
            s3 = getattr(ch_to_clean, 's3_'+cont)
            ch_to_clean.ch_clean(s3_mask=s3_mask, s3=s3, inverted=inverted, plot_settings=plot_settings)

            #Update workflow 
            proc_up = ['ImProc',ch,'E-CleanCh','Info',s3.cont_type, 'Status']
            organ.update_mHworkflow(proc_up, 'DONE')
        
            _ = ch_to_clean.s32Meshes(cont_type=cont, 
                                    #   keep_largest=organ.mH_settings['wf_info']['keep_largest'][ch][cont], 
                                    #   rotateZ_90=organ.mH_settings['setup']['rotateZ_90'],
                                      new_set=True)
            aa+=1
            win.prog_bar_update(aa)

        # Update organ workflow
        organ.update_mHworkflow(process, 'DONE')
        organ.check_status(process = 'ImProc')
        print('organ.workflow:', workflow)

        #Message User
        win.win_msg('Contours of channel '+ch[-1]+' were successfully cleaned!')

        #Enable button for plot
        getattr(win, 'cleanup_plot_'+ch).setEnabled(True)
        getattr(win, 'summary_whole_plot_'+ch).setEnabled(True)

    #Update Status in GUI
    win.update_status(workflow, process, win.cleanup_status)
    
    print('\nEND Cleanup')
    print('organ.mH_settings:', organ.mH_settings)
    print('organ.workflow:', workflow)
   
def trim_top_bottom_S3s(organ, meshes, no_cut, cuts_out, win):#
    
    workflow = organ.workflow['morphoHeart']
    #Update workflow for channels that won't be cut
    for ch_cont in no_cut:
        ch, cont = ch_cont.split('_')
        proc_im = ['ImProc', ch, 'E-TrimS3']
        proc_ms = ['MeshesProc', 'B-TrimMesh', ch, cont, 'Status']
        organ.update_mHworkflow(proc_im+['Info',cont,'Status'],'DONE-NoCut')
        organ.update_mHworkflow(proc_im+['Status'],'DONE-NoCut')
        organ.update_mHworkflow(proc_ms,'DONE-NoCut')
    
    #Cut channel s3s and recreate meshes
    proc_ms_all = ['MeshesProc', 'B-TrimMesh', 'Status']
    for mesh in meshes:
        ch_name = mesh.imChannel.channel_no
        im_ch = mesh.imChannel
        proc_im = ['ImProc', ch_name, 'E-TrimS3']
        
        win.prog_bar_range(0,3)
        aa = 0
        for cont in ['int', 'tiss', 'ext']: 
            win.win_msg('Trimming meshes of Channel '+im_ch.channel_no[-1]+'! ('+str(aa+1)+'/3)')
            proc_ms = ['MeshesProc', 'B-TrimMesh', ch_name, cont, 'Status']
            cuts = []
            for side in ['top', 'bottom']: 
                if cuts_out[side]['chs'][ch_name][cont]:
                    cuts.append(side)
            if len(cuts) > 0: 
                im_ch.trimS3(cuts=cuts, cont=cont, cuts_out=cuts_out)
                _ = im_ch.s32Meshes(cont_type=cont, 
                                    new_set=True)

                # Update organ workflow
                process = proc_im+['Info',cont,'Status']
                organ.update_mHworkflow(process,'DONE')
                organ.update_mHworkflow(proc_ms,'DONE')
                #Update progress in main_win
                win.update_workflow_progress()

            aa+=1
            win.prog_bar_update(aa)

        # Update organ workflow
        process_up = proc_im+['Status']
        organ.update_mHworkflow(process_up,'DONE')

        #Message User
        win.win_msg(' Channel '+im_ch.channel_no[-1]+' meshes were successfully trimmed!')

        #Enable button for plot
        getattr(win, 'trimming_plot_'+im_ch.channel_no).setEnabled(True)
        getattr(win, 'summary_whole_plot_'+im_ch.channel_no).setEnabled(True)

    # Update organ workflow
    organ.update_mHworkflow(proc_ms_all, 'DONE')

    # Update mH_settings
    proc_set = ['wf_info']
    update = cuts_out
    organ.update_settings(proc_set, update, 'mH', add='trimming')

    #Update Status in GUI
    win.update_status(workflow, proc_ms_all, win.trimming_status)

    print('\nEND Trimming')
    print('organ.mH_settings:', organ.mH_settings)
    print('organ.workflow:', workflow)

def reset_meshes(controller, option): 
    win = controller.main_win
    temp_s3s = getattr(win, 'temp_s3s_'+option)
    organ = win.organ

    workflow = organ.workflow['morphoHeart']
    if option == 'trim': 
        proc_ms_all = ['MeshesProc', 'B-TrimMesh', 'Status']

    for ch_name in temp_s3s.keys():
        im_ch = organ.obj_imChannels[ch_name]
        if option == 'trim': 
            proc_im = ['ImProc', ch_name, 'E-TrimS3']
            key_word = ' trimmed '
        else: 
            proc_im = ['ImProc', ch_name, 'E-CleanCh']
            key_word = ' cleaned '
        
        win.prog_bar_range(0,3)
        aa = 0
        for cont in ['int', 'tiss', 'ext']: 
            win.win_msg('Resetting'+key_word+'meshes of Channel '+im_ch.channel_no[-1]+'! ('+str(aa+1)+'/3)')
            if option == 'trim':
                proc_ms = ['MeshesProc', 'B-TrimMesh', ch_name, cont, 'Status']

            im_ch.load_chS3s(cont_types=[cont])
            s3 = getattr(im_ch, 's3_'+cont)
            temp_s3 = temp_s3s[ch_name][cont]
            s3.s3_save(temp_s3)
            _ = im_ch.s32Meshes(cont_type=cont, 
                        new_set=True)
            
            # Update organ workflow
            process = proc_im+['Info',cont,'Status']
            organ.update_mHworkflow(process,'NI')
            if option == 'trim': 
                organ.update_mHworkflow(proc_ms,'NI')

            #Update progress in main_win
            win.update_workflow_progress()

            aa+=1
            win.prog_bar_update(aa)
        
    # Update organ workflow 
    process_up = proc_im+['Status']
    organ.update_mHworkflow(process_up,'NI')
    if option == 'trim': 
        organ.update_mHworkflow(proc_ms_all,'NI')
        #Update Status in GUI
        win.update_status(workflow, proc_ms_all, win.trimming_status)
        #Reset buttons too
        win.trimming_play.setChecked(False)
        #Disable button
        win.reset_trimmed_s3s.setEnabled(False)
    else: 
        #Update Status in GUI
        win.update_status(workflow, process_up, win.cleanup_status)
        #Reset buttons too
        win.cleanup_play.setChecked(False)
        #Disable button
        win.reset_cleaned_s3s.setEnabled(False)

    #Message User
    win.win_msg(' Channel '+im_ch.channel_no[-1]+' meshes were successfully reset!')
    #Update progress in main_win
    win.update_workflow_progress()

def get_stack_orientation(organ, gui_orientation, win):# 

    workflow = organ.workflow['morphoHeart']
    colors = [[255,215,0,200],[0,0,205,200],[255,0,0,200]]
    views = organ.mH_settings['setup']['orientation']['stack'].split(', ')
    planar_views, stack_cube = organ.get_orientation(views, colors, mtype='STACK')
    planar_views = get_ref_vectors(planar_views)

    #Update organ workflow
    process = ['MeshesProc', 'A-Set_Orientation']
    organ.update_mHworkflow(process+['Stack'],'DONE')
    roi_status = get_by_path(workflow, process+['ROI'])
    if  roi_status == 'NI':
        organ.update_mHworkflow(process+['Status'],'Initialised')
    elif roi_status == 'DONE':
        organ.update_mHworkflow(process+['Status'],'DONE')

    print('planar_views: ',planar_views)
    # Update mH_settings
    gui_orientation['stack']['planar_views'] = planar_views
    gui_orientation['stack']['stack_cube'] = stack_cube
    proc_set = ['wf_info']
    update = gui_orientation
    organ.update_settings(proc_set, update, 'mH', add='orientation')

    #Enable button for plot
    getattr(win, 'stack_orient_plot').setEnabled(True)

    print('\nEND Stack Orientation')
    print('organ.mH_settings:', organ.mH_settings)
    print('organ.workflow:', workflow)

def get_roi_orientation(organ, gui_orientation:dict, win):#

    workflow = organ.workflow['morphoHeart']
    on_hold = False
    colors = [[255,215,0,200],[0,0,205,200],[255,0,0,200]]
    planar_views = None; settings = None; roi_cube = None
    if gui_orientation['roi']['reorient']: 
        if gui_orientation['roi']['method'] == 'Centreline': 
            centreline = gui_orientation['roi']['centreline']
            #Check if centrelines have been obtained
            if workflow['MeshesProc']['C-Centreline']['Status'] == 'DONE': 
                win.win_msg("Setting Organ/ROI Orientation with "+centreline+"'s centreline...")
                planar_views, settings, roi_cube = organ.get_ROI_orientation(gui_orientation, colors)
                print('Method: Centreline -', planar_views, settings)
            else: 
                on_hold = True
                win.win_msg('!Organ/ROI Orientation will be set once the centreline of the '+centreline+' has been obtained.')
        else: #Manual
            planar_views, settings, roi_cube = organ.get_ROI_orientation(gui_orientation, colors)
            win.rotX.setText("%.1f" % roi_cube['rotate_user']['rotX']+'°')
            win.rotY.setText("%.1f" % roi_cube['rotate_user']['rotY']+'°')
            win.rotZ.setText("%.1f" % roi_cube['rotate_user']['rotZ']+'°')
            print('Method: Manual -', planar_views)
    else: 
        planar_views, settings, roi_cube = organ.get_ROI_orientation(gui_orientation, colors)
        print('Method: No rotation -', planar_views)

    if not on_hold: 
        #Update organ workflow
        process = ['MeshesProc', 'A-Set_Orientation']
        organ.update_mHworkflow(process+['ROI'],'DONE')
        stack_status = get_by_path(workflow, process+['Stack'])
        if  stack_status == 'NI':
            organ.update_mHworkflow(process+['Status'],'Initialised')
        elif stack_status == 'DONE':
            organ.update_mHworkflow(process+['Status'],'DONE')

        # Update mH_settings
        gui_orientation['roi']['planar_views'] = planar_views
        gui_orientation['roi']['roi_cube'] = roi_cube
        if settings != None: 
            gui_orientation['roi']['settings'] = settings
        
        if 'orientation' in organ.mH_settings['wf_info'].keys():
            proc_set = ['wf_info', 'orientation', 'roi']
            update = gui_orientation['roi']
        else: 
            proc_set = ['wf_info']
            update = gui_orientation
        
        #Enable button for plot
        getattr(win, 'roi_orient_plot').setEnabled(True)
        
    else: 
        proc_set = ['wf_info']
        update = gui_orientation
    organ.update_settings(proc_set, update, 'mH', add='orientation')

    print('\nEND ROI Orientation')
    print('organ.mH_settings:', organ.mH_settings)
    print('organ.workflow:', workflow)

    return on_hold

def get_ref_vectors(planar_views): 

    for view in planar_views:
        pl_normal = unit_vector(planar_views[view]['pl_normal'])
        #Find the direction of the normal 
        index = np.argmax(pl_normal)
        if index == 0: 
            ref_vector = np.array([0.0, 1.0, 0.0])
            proj_plane = 'x'
        elif index == 1: 
            ref_vector = np.array([0.0, 0.0, 1.0])
            proj_plane = 'y'
        else: 
            ref_vector = np.array([1.0, 0.0, 0.0])
            proj_plane = 'z'
        print(view, pl_normal, proj_plane)
        planar_views[view]['ref_vector'] = ref_vector
        planar_views[view]['proj_plane'] = proj_plane

    return planar_views

def extract_chNS(organ, rotateZ_90, win, plot_settings):#
    from .mH_classes_new import ImChannelNS

    workflow = organ.workflow['morphoHeart']
    im_ns = ImChannelNS(organ=organ, ch_name='chNS')
    im_ns.create_chNSS3s(win, plot_settings=plot_settings)
    proc_im = ['ImProc',im_ns.channel_no,'D-S3Create','Status']
    proc_mesh = ['MeshesProc','A-Create3DMesh', im_ns.channel_no,'Status']

    gui_keep_largest = {'int': True, 'ext': True, 'tiss': False}
    win.prog_bar_range(0,3)
    aa = 0
    for cont in ['int', 'tiss', 'ext']: 
        win.win_msg('Creating meshes of Channel NS! ('+str(aa+1)+'/3)')
        im_ns.s32Meshes(cont_type=cont, 
                        keep_largest=gui_keep_largest[cont],
                        rotateZ_90 = rotateZ_90,
                        new_set=True)
        aa+=1
        win.prog_bar_update(aa)

    #Message User
    win.win_msg(' ChannelNS meshes were successfully created!')

    # Update organ workflow
    organ.update_mHworkflow(proc_im, 'DONE')
    organ.update_mHworkflow(proc_mesh, 'DONE')

    #Enable button for plot
    plot_btn = getattr(win, 'chNS_plot')
    plot_btn.setEnabled(True)

    #Update Status in GUI
    win.update_status(workflow, proc_mesh, win.chNS_status)

    print('\nEND chNS')
    print('organ.mH_settings:', organ.mH_settings)
    print('organ.workflow:', workflow)
    
#Centreline functions
def proc_meshes4cl(organ, win):#
    """
    Funtion that cuts the inflow and outflow tract of meshes from which 
    the centreline will be obtained.

    """
    same_planes = win.gui_centreline['SimplifyMesh']['same_planes']
    tol = win.gui_centreline['SimplifyMesh']['tol']
    workflow = organ.workflow['morphoHeart']
    
    # Set cut names 
    cuts_names = {'top': {'heart_def': 'outflow tract','other': 'top'},
                'bottom': {'heart_def': 'inflow tract','other': 'bottom'}}
    
    if mH_config.heart_default: 
        name_dict =  'heart_def'     
    else: 
        name_dict = 'other'
        
    # Assign colors to top and bottom splines and spheres
    color_cuts = {'bottom': 'salmon','top':'mediumseagreen' }
    
    # Get organ settings
    filename = organ.user_organName
    ext = 'stl'; directory = organ.dir_res(dir ='centreline')

    # Get the meshes from which the centreline needs to be extracted
    cl_names = list(organ.mH_settings['measure']['CL'].keys())
    # Create lists with mesh_mH, channel number and mesh names
    mH_mesh4cl = []; chs4cl = []
    for n, name in enumerate(cl_names): 
        ch, cont, _ = name.split('_')
        mH_mesh = organ.obj_meshes[ch+'_'+cont]
        mH_mesh4cl.append(mH_mesh)
        chs4cl.append(organ.obj_meshes[ch+'_'+cont].channel_no)
    
    print(cl_names, mH_mesh4cl, chs4cl)

    # Get dictionary with initialised planes to cut top/bottom
    planes_info = {'top':None, 'bottom':None}
    if 'trimming' in organ.mH_settings['wf_info'].keys():
        trimming_info = organ.mH_settings['wf_info']['trimming']
        for side in ['bottom', 'top']: 
            try: 
                planes_info[side] = trimming_info[side]['plane_info_mesh']
            except:
                print('No plane has been created for this side:', side)
                planes_info[side] = None
    else: 
        planes_info = {'bottom': None, 'top': None}
    
    ksplines = []; spheres = []; m4clf={}
    plane_cuts = {'bottom': {'dir': True, 'plane': None, 'pl_dict': None}, 
                    'top': {'dir': False,'plane': None, 'pl_dict': None}}
    plane_cuts_all = {}
    ch_cuts = {}
    proc_all = ['MeshesProc','C-Centreline','SimplifyMesh', 'Status']
    for n, mH_msh in enumerate(mH_mesh4cl): # iterate through meshes2cut
        kspl_ch = {}
        sp_ch = {}
        mesh_ch = mH_msh.channel_no
        mesh_cont = mH_msh.mesh_type
        proc_mesh = ['MeshesProc','C-Centreline','SimplifyMesh', mesh_ch, mesh_cont, 'Status']
        for nn, pl_cut in zip(count(), plane_cuts.keys()): # iterate through cuts (top and bottom)
            if nn == 0 and same_planes: 
                if n == 0: 
                    settings = {'color': {}, 'name':{}}
                    meshes_in = []
                    for aa, msh in enumerate(mH_mesh4cl): 
                        print('>> Smoothing mesh - ',pl_cut, 'direction:', plane_cuts[pl_cut]['dir'])
                        sm_msh_o = msh.mesh4CL()
                        meshes_in.append(sm_msh_o)
                        settings['color'][aa] = msh.color
                        settings['name'][aa] = msh.legend
                        
                sm_msh = meshes_in[n]
            elif nn == 0 and not same_planes:
                print('>> Smoothing mesh - ',pl_cut, 'direction:', plane_cuts[pl_cut]['dir'])
                sm_msh = mH_msh.mesh4CL()
                settings = {'color': {0: mH_msh.color}, 'name':{0: mH_msh.legend}}
            else: 
                pass
            print('settings:', settings)

            # Get planes for first mesh
            if n == 0 and same_planes: 
                if planes_info[pl_cut] == None:
                    #Planes have not been initialised
                    print('-Planes have not been initialised for ', pl_cut)
                    plane, pl_dict = get_plane(filename=filename, 
                                                txt = 'cut '+cuts_names[pl_cut][name_dict],
                                                meshes = meshes_in, settings = settings) 
                else:
                    print('-Planes had been initialised for ', pl_cut)
                    # Planes have been initialised
                    plane, pl_dict = get_plane(filename=filename, 
                                                txt = 'cut '+cuts_names[pl_cut][name_dict],
                                                meshes = meshes_in, settings = settings, 
                                                def_pl = planes_info[pl_cut]) 
                
                plane_cuts[pl_cut]['plane'] = plane
                plane_cuts[pl_cut]['pl_dict'] = pl_dict
            
            elif not same_planes: #This planes are specific for the mesh being cut
                if planes_info[pl_cut] == None:
                    #Planes have not been initialised
                    print('-Planes have not been initialised for ', pl_cut)
                    plane, pl_dict = get_plane(filename=filename, 
                                                txt = 'cut '+cuts_names[pl_cut][name_dict],
                                                meshes = [sm_msh], settings = settings) 
                else:
                    print('-Planes had been initialised for ', pl_cut)
                    # Planes have been initialised
                    plane, pl_dict = get_plane(filename=filename, 
                                                txt = 'cut '+cuts_names[pl_cut][name_dict],
                                                meshes = [sm_msh], settings = settings,
                                                def_pl = planes_info[pl_cut]) 
                
                plane_cuts[pl_cut]['plane'] = plane
                plane_cuts[pl_cut]['pl_dict'] = pl_dict

            print('> Cutting mesh: ', mH_msh.legend, '-', pl_cut)
            pts2cut, _ = get_pts_at_plane(points = sm_msh.points(), 
                                            pl_normal = plane_cuts[pl_cut]['pl_dict']['pl_normal'],
                                            pl_centre = plane_cuts[pl_cut]['pl_dict']['pl_centre'], tol=tol)
            ordpts, _ = order_pts(points = pts2cut)
        
            # Create spline around cut
            kspl = vedo.KSpline(ordpts, continuity=0, tension=0, bias=0, closed=True)
            kspl.color(color_cuts[pl_cut]).legend('cut4CL_'+pl_cut).lw(2)
            organ.add_object(kspl, proc='cut4cl', class_name=[pl_cut,mH_msh.name], name='KSpline')
            ksplines.append(kspl)
            
            # Get centroid of kspline to add to the centreline
            kspl_bounds = kspl.bounds()
            pt_centroid = np.mean(np.asarray(kspl_bounds).reshape((3, 2)),axis=1)
            sph_centroid = vedo.Sphere(pos=pt_centroid, r=2).legend('cut4CL_'+pl_cut).color(color_cuts[pl_cut])
            organ.add_object(sph_centroid, proc='cut4cl', class_name=[pl_cut,mH_msh.name], name='Sphere')
            spheres.append(sph_centroid)
            
            # Cutmesh using created plane
            msh_new = sm_msh.clone().cut_with_mesh(plane_cuts[pl_cut]['plane'], invert=plane_cuts[pl_cut]['dir'])
            msh_new = msh_new.extract_largest_region()
            msh_new.alpha(0.05).wireframe(True).legend(mH_msh.legend)
            
            kspl_ch[pl_cut] = kspl
            sp_ch[pl_cut] = sph_centroid
            sm_msh = msh_new.clone()

        m4clf[mesh_ch+'_'+mesh_cont] = {'mesh': sm_msh, 
                                        'kspl' : kspl_ch,
                                        'centroid': sp_ch}
        
        print('> Saving cut meshes')
        mesh_title = filename+"_"+mH_msh.name+"_cut4cl."+ext
        mesh_titleML = filename+"_"+mH_msh.name+"_cut4clML."+ext
        mesh_dir = directory / mesh_title
        msh_new.write(str(mesh_dir))

        if mesh_ch not in ch_cuts.keys():
            ch_cuts[mesh_ch] = {mesh_cont: {'dir_cleanMesh': mesh_title, 
                                            'dir_meshLabMesh': mesh_titleML}}
        else: 
            ch_cuts[mesh_ch][mesh_cont] = {'dir_cleanMesh': mesh_title, 
                                            'dir_meshLabMesh': mesh_titleML}

        if same_planes: 
            plane_cuts_all = plane_cuts
        else: 
            plane_cuts_all[mesh_ch+'_'+mesh_cont] = plane_cuts

        # Update organ workflow
        organ.update_mHworkflow(proc_mesh, 'DONE')

        #Find the position of the cl in the cl_list of the smoothed mesh that was obtained
        index_cl = cl_names.index(mesh_ch+'_'+mesh_cont+'_whole')+1
        #Enable button for plot
        plot_btn = getattr(win, 'cl_clean_plot'+str(index_cl))
        plot_btn.setEnabled(True)

        #Update Status in GUI
        status_sq = getattr(win, 'clClean_status'+str(index_cl))
        win.update_status(workflow, proc_mesh, status_sq)

    #Message User
    win.win_msg('All meshes for centreline have been successfully smoothed!')

    # Update organ workflow
    organ.update_mHworkflow(proc_all, 'DONE')
    organ.update_mHworkflow(process =  ['MeshesProc','C-Centreline','Status'], update = 'Initialised')

    #Remove plane key from plane cuts to add to mH_settings['wf_info']
    if not same_planes: 
        for mH_msh in mH_mesh4cl: 
            for side in ['bottom', 'top']: 
                plane_cuts_all[mH_msh.name][side].pop('plane', None)
    else: 
        for side in ['bottom', 'top']: 
            plane_cuts_all[side].pop('plane', None)

    #Update mH_settings
    win.gui_centreline['SimplifyMesh']['plane_cuts'] = plane_cuts_all
    win.gui_centreline['dirs'] = ch_cuts
    proc_set = ['wf_info']
    update =  win.gui_centreline
    organ.update_settings(proc_set, update, 'mH', add='centreline')
    print('FFF:',organ.mH_settings['wf_info']['centreline'])
    
    print('\nEND Simplifying Mesh')
    print('organ.mH_settings:', organ.mH_settings)
    print('organ.workflow:', workflow)
    print(' win.gui_centreline:',  win.gui_centreline)

    print('m4clf:', m4clf)
    return m4clf

def extract_cl(organ, win, voronoi=False):#
    
    workflow = organ.workflow['morphoHeart']
    proc_all = ['MeshesProc','C-Centreline','vmtk_CL', 'Status']
    cl_names = list(organ.mH_settings['measure']['CL'].keys())
    nn = 0
    done = []
    for name in cl_names: 
        ch, cont, _ = name.split('_')
        proc_sp = ['MeshesProc','C-Centreline','vmtk_CL', ch, cont, 'Status']
        #Get the vmtk txt
        vmtktxt, dir_npcl = code_vmtk(organ, ch, cont, voronoi)
        
        myArguments = vmtktxt
        try: 
            myPype = pypes.PypeRun(myArguments)
            
            for m in range(len(myPype.ScriptObjectList)):
                if isinstance(myPype.ScriptObjectList[m], vmtk.vmtkcenterlines.vmtkCenterlines):
                    centerlineReader = myPype.ScriptObjectList[m]
            clNumpyAdaptor = vmtkscripts.vmtkCenterlinesToNumpy()
            clNumpyAdaptor.Centerlines = centerlineReader.Centerlines # Surface # cl
            clNumpyAdaptor.Execute()
            numpyCenterlines = clNumpyAdaptor.ArrayDict

            # Copy dictionary into new dictionary to save
            centrelines_dict = dict()
            centrelines_dict['Points'] =  numpyCenterlines['Points']

            cellData = centrelines_dict['CellData'] = dict()
            cellData['CellPointIds'] = numpyCenterlines['CellData']['CellPointIds']
            pointData = centrelines_dict['PointData'] = dict()
            pointData['EdgeArray'] = numpyCenterlines['PointData']['EdgeArray']
            pointData['EdgePCoordArray'] = numpyCenterlines['PointData']['EdgePCoordArray']
            pointData['MaximumInscribedSphereRadius'] = numpyCenterlines['PointData']['MaximumInscribedSphereRadius']

            with open(dir_npcl, "w") as write_file:
                json.dump(centrelines_dict, write_file, cls=NumpyArrayEncoder)
                print('>> Dictionary saved correctly!\n> File: '+ dir_npcl.name);
                alert('countdown')

            cl_data = load_vmtkCL(organ, ch, cont)
            #Get cl points from vmtk
            pts_cl = np.asarray(cl_data['Points'])

        except: 
            win.win_msg('*Something went wrong when trying to obtain the centreline for '+name+'!')
            return False
        
        try: 
            if len(pts_cl) > 0: 
                # Update organ workflow
                organ.update_mHworkflow(proc_sp, 'DONE')
                status_sq = getattr(win, 'vmtk_status'+str(nn+1))
                win.update_status(workflow, proc_sp, status_sq)
                nn+=1
                done.append(True)
        except: 
            done.append(False)
            win.win_msg('*Something went wrong when trying to obtain the centreline for '+name+'!')

    if all(done): 
        #Message User
        win.win_msg('The centreline for all meshes have been successfully obtained!')
        # Update organ workflow
        organ.update_mHworkflow(proc_all, 'DONE')
        return True
    else: 
        win.win_msg('*Something went wrong when trying to obtain the centreline for the selected meshes!')
        return False
    
        # --------
        # ArrayDict
        #     ['Points']                   <-- required, is Nx3 array of N vertexes and x, y, z locations
        #     ['PointData']                <-- required, even if subarrays are empty
        #         ['PointDataArray1']      <-- optional, (ex. MaximumInscribedSphereRadius)
        #         ['PointDataArray2']      <-- optional
        #         ...
        #     ['CellData']                 <-- required
        #         ['CellPointIds']         <-- required, list of Mx1 arrays defining cell connectivity to ['Points']
        #         ['CellDataArray1']       <-- optional, (ex: CenterlineTractId)
        #         ['CellDataArray2']       <-- optional
        #            ...

def code_vmtk(organ, ch, cont, voronoi=False):#
    """
    Function that gets directories information and prints a series of instructions to process the meshes to obtain centreline and
    run the vmtk code.

    """
    dir_cl = organ.dir_res(dir ='centreline')
    dir_meshML = dir_cl / organ.mH_settings['wf_info']['centreline']['dirs'][ch][cont]['dir_meshLabMesh']
    dir_npcl = dir_meshML.name.replace('_cut4clML.stl', '_npcl.json')
    dir_npcl = dir_cl / dir_npcl
    
    dir_meshML_str = str(dir_meshML).replace('\\','/')
    dir_meshML_strF = '"'+dir_meshML_str+'"'
    
    vtp_name = dir_meshML.name.replace('_cut4clML.stl', '_cl.vtp')
    dir_vtp = dir_cl / vtp_name
    dir_vtp_str = str(dir_vtp).replace('\\','/')
    dir_vtp_strF = '"'+dir_vtp_str+'"'
    
    # General code
    # vmtktxt = "vmtksurfacereader -ifile "+ dir_meshML_strF +" --pipe vmtksurfacesmoothing -passband 0.1 -iterations 30 --pipe vmtkcenterlines -seedselector openprofiles -ofile " + dir_vtp_strF + " --pipe vmtkrenderer --pipe vmtksurfaceviewer -opacity 0.25 --pipe vmtksurfaceviewer -i @vmtkcenterlines.o -array MaximumInscribedSphereRadius"
    # code for voronoi diagrams
    # vmtktxt = "vmtksurfacereader -ifile "+ dir_meshML_strF +" --pipe vmtksurfacesmoothing -passband 0.1 -iterations 30 --pipe vmtkcenterlines -seedselector openprofiles -ofile " + dir_vtp_strF + " --pipe vmtkrenderer --pipe vmtksurfaceviewer -opacity 0.25 --pipe vmtksurfaceviewer -i @vmtkcenterlines.voronoidiagram -array MaximumInscribedSphereRadius -i @vmtkcenterlines.o -array MaximumInscribedSphereRadius"

    vmtk_txt1_reader = "vmtksurfacereader -ifile "+ dir_meshML_strF+" "
    vmtk_txt2_smooth = "--pipe vmtksurfacesmoothing -passband 0.1 -iterations 30 "
    vmtk_txt3_cl = "--pipe vmtkcenterlines -seedselector openprofiles -ofile "+ dir_vtp_strF + " "
    if voronoi: 
        vmtk_txt4_viewer = "--pipe vmtkrenderer --pipe vmtksurfaceviewer -opacity 0.25 --pipe vmtksurfaceviewer -i @vmtkcenterlines.voronoidiagram -array MaximumInscribedSphereRadius -i @vmtkcenterlines.o -array MaximumInscribedSphereRadius"
    else: 
        vmtk_txt4_viewer = "--pipe vmtkrenderer --pipe vmtksurfaceviewer -opacity 0.25 --pipe vmtksurfaceviewer -i @vmtkcenterlines.o -array MaximumInscribedSphereRadius"
    vmtk_txtF = vmtk_txt1_reader + vmtk_txt2_smooth + vmtk_txt3_cl + vmtk_txt4_viewer

    return vmtk_txtF, dir_npcl

def load_vmtkCL(organ, ch, cont):#
    
    dir_cl = organ.dir_res(dir ='centreline')#organ.info['dirs']['centreline']
    dir_meshML = dir_cl / organ.mH_settings['wf_info']['centreline']['dirs'][ch][cont]['dir_meshLabMesh']
    dir_npclo = dir_meshML.name.replace('_cut4clML.stl', '_npcl.json')
    dir_npcl = dir_cl / dir_npclo
    
    json2open_dir = dir_npcl
    #print("Started Reading JSON file")
    with open(json2open_dir, "r") as read_file:
        print("\t>> "+dir_npclo+": Converting JSON encoded data into Numpy Array")
        decodedArray = json.load(read_file)

    return decodedArray

def define_CL(organ, name, nPoints): 

    ch, cont, _ = name.split('_')
    cl_data = load_vmtkCL(organ, ch, cont)
    #Get cl points from vmtk
    pts_cl = np.asarray(cl_data['Points'])
    # Interpolate points of original centreline
    pts_int_o = get_interpolated_pts(points=pts_cl, nPoints = nPoints)
    # Create IFT and OFT spheres to show original points position
    sph_inf_o = vedo.Sphere(pts_int_o[-1], r = 3, c='tomato')
    sph_outf_o = vedo.Sphere(pts_int_o[0], r = 3, c='navy')
    
    # Get outflow point for that layer
    pt2add_outf = organ.objects['Spheres']['cut4cl']['top'][ch+'_'+cont]['center']
    pts_withOutf = np.insert(pts_cl, 0, np.transpose(pt2add_outf), axis=0)
    
    pt2add_inf = organ.objects['Spheres']['cut4cl']['bottom'][ch+'_'+cont]['center']
    pts_all_opt1 = np.insert(pts_withOutf, len(pts_withOutf), np.transpose(pt2add_inf), axis=0)

    # Interpolate points
    pts_int_opt1 = get_interpolated_pts(points=pts_all_opt1, nPoints = nPoints)
    # Create kspline with points
    kspl_opt1 = vedo.KSpline(pts_int_opt1, res = nPoints)
    kspl_opt1.color('blueviolet').legend('CL_'+ch+'_'+cont).lw(5)
    dict_clOpt = {'Option 3': {'kspl': kspl_opt1, 'sph_bot': sph_outf_o, 
                                'sph_top': sph_inf_o, 'pt2add': pt2add_inf, 
                                'description': 'Point in meshesCut4Cl'}}
    
    mesh = organ.obj_meshes[ch+'_'+cont]
    kspl_f = mesh.modify_centreline(kspl_o=kspl_opt1, mesh=mesh.mesh, ext=False)
    pts_cl_extf = kspl_f.points()
       
    # Increase the resolution of the extended centreline and interpolate to unify sampling
    xd = np.diff(pts_cl_extf[:,0])
    yd = np.diff(pts_cl_extf[:,1])
    zd = np.diff(pts_cl_extf[:,2])
    dist = np.sqrt(xd**2+yd**2+zd**2)
    u = np.cumsum(dist)
    u = np.hstack([[0],u])
    t = np.linspace(0, u[-1], nPoints)#601
    try: 
        resamp_pts = interpn((u,), pts_cl_extf, t)
        print('Created resampled KSpline using interpn')
    except: 
        # Christian K Answer
        # https://stackoverflow.com/questions/19117660/how-to-generate-equispaced-interpolating-values/19122075#19122075
        xn = np.interp(t, u, pts_cl_extf[:,0])
        yn = np.interp(t, u, pts_cl_extf[:,1])
        zn = np.interp(t, u, pts_cl_extf[:,2])
        resamp_pts = np.vstack((xn,yn,zn))
        resamp_pts = resamp_pts.T
        print('Created resampled KSpline using np.interp')

    kspl_ext = vedo.KSpline(resamp_pts, res=nPoints).color('pink').legend('ExtendedCL').lw(2)#601

    return kspl_ext, dict_clOpt
    
def create_CLs(organ, name, nPoints):#
    """
    Function that creates the centrelines using the points given as input in the dict_cl

    """
    cl_colors = {'Op1':'navy','Op2':'blueviolet'}
    
    ch, cont, _ = name.split('_')
    cl_data = load_vmtkCL(organ, ch, cont)
    #Get cl points from vmtk
    pts_cl = np.asarray(cl_data['Points'])
    # Interpolate points of original centreline
    pts_int_o = get_interpolated_pts(points=pts_cl, nPoints = nPoints)
    # Last CL point
    pt_m1 = pts_int_o[-1]
    sph_m1 = vedo.Sphere(pos = pt_m1, r=4, c='black')
    # Create kspline with original points
    kspl_o = vedo.KSpline(pts_int_o, res = nPoints).color('gold').legend('CLo_'+ch+'_'+cont).lw(5)
    # Create IFT and OFT spheres to show original points position
    sph_inf_o = vedo.Sphere(pts_int_o[-1], r = 3, c='tomato')
    sph_outf_o = vedo.Sphere(pts_int_o[0], r = 3, c='navy')
    
    # Get outflow point for that layer
    pt2add_outf = organ.objects['Spheres']['cut4cl']['top'][ch+'_'+cont]['center']
    pts_withOutf = np.insert(pts_cl, 0, np.transpose(pt2add_outf), axis=0)
    
    #Dictionary with kspl options
    dict_clOpt = {}

    # Get inflow point for that layer (four options)
    # - Option 1, use centroids of kspline cuts
    pt2add_inf = organ.objects['Spheres']['cut4cl']['bottom'][ch+'_'+cont]['center']
    pts_all_opt1 = np.insert(pts_withOutf, len(pts_withOutf), np.transpose(pt2add_inf), axis=0)

    # Interpolate points
    pts_int_opt1 = get_interpolated_pts(points=pts_all_opt1, nPoints = nPoints)
    # Create kspline with points
    kspl_opt1 = vedo.KSpline(pts_int_opt1, res = nPoints)
    kspl_opt1.color(cl_colors['Op1']).legend('(Op1) CL_'+ch+'_'+cont).lw(5)
    dict_clOpt['Option 1'] = {'kspl': kspl_opt1, 'sph_bot': sph_outf_o, 
                                'sph_top': sph_inf_o, 'pt2add': pt2add_inf, 
                                'description': 'Point in meshesCut4Cl'}
    
    # - Option 6 (add mid point between final point of original centreline and inf tract)
    pt_m6 = pt_m1+(pt2add_inf-pt_m1)/2
    sph_m6 = vedo.Sphere(pos = pt_m6, r=4, c='gold')

    pts_all_opt6 = np.insert(pts_withOutf, len(pts_withOutf), np.transpose(pt_m6), axis=0)
    pts_all_opt6f = np.insert(pts_all_opt6, len(pts_all_opt6), np.transpose(pt2add_inf), axis=0)

    # Interpolate points
    pts_int_opt6 = get_interpolated_pts(points=pts_all_opt6f, nPoints = nPoints)
    # Create kspline with points
    kspl_opt6 = vedo.KSpline(pts_int_opt6, res = nPoints)
    kspl_opt6.color(cl_colors['Op2']).legend('(Op2) CL_'+ch+'_'+cont).lw(5)
    dict_clOpt['Option 2'] = {'kspl': kspl_opt6, 'sph_bot': sph_m6, 
                                'sph_top': sph_inf_o, 'pt2add': pt2add_inf,
                                'description': 'Center point between last and meshesCut4Cl'}
    
    obj = []; txt = []
    mesh = organ.obj_meshes[ch+'_'+cont].mesh.clone()
    for num, opt in enumerate(dict_clOpt.keys()):
        obj2add = (mesh.alpha(0.01), sph_m1, kspl_o, dict_clOpt[opt]['kspl'], dict_clOpt[opt]['sph_bot'], dict_clOpt[opt]['sph_top'])
        obj.append(obj2add)
        if num == 0: 
            txt.append((num, organ.user_organName +' - '+opt))
        else: 
            txt.append((num, opt))
    
    plot_grid(obj=obj, txt=txt, axes=5, sc_side=max(organ.get_maj_bounds())) 

    return dict_clOpt     

def extend_CL(pts_int_o, pts_withOutf, num, nPoints, plane_info):
    
    pt_m10 = pts_int_o[-num]
    sph_m10 = vedo.Sphere(pos = pt_m10, r=4, c='lime')
    pt_m1 = pts_int_o[-1]
    # sph_m1 = vedo.Sphere(pos = pt_m1, r=4, c='tomato')
    dir_v = unit_vector(pt_m1-pt_m10)*3
    pts_ext = np.array([pt_m10, pt_m1, pt_m1+dir_v*2,pt_m1+dir_v*4,pt_m1+dir_v*6,pt_m1+dir_v*8,pt_m1+dir_v*10])
    xd = np.diff(pts_ext[:,0])
    yd = np.diff(pts_ext[:,1])
    zd = np.diff(pts_ext[:,2])
    dist = np.sqrt(xd**2+yd**2+zd**2)
    u = np.cumsum(dist)
    u = np.hstack([[0],u])
    t = np.linspace(0, u[-1],nPoints)
    resamp_pts = interpn((u,), pts_ext, t)
    kspl_resamp = vedo.KSpline(resamp_pts, res = nPoints).lw(5)
    
    #Get plane
    pl_centre_bot = plane_info['pl_centre']
    pl_normal_bot = plane_info['pl_normal']
    pl_IFT = vedo.Plane(pos = pl_centre_bot, normal = pl_normal_bot)
    
    #Cut KSpline
    kspl_test = kspl_resamp.clone().cutWithMesh(pl_IFT, invert = True)
    pt2add = kspl_test.points()[-1]
    
    #Insert last point into the original CL
    pts_all_opt = np.insert(pts_withOutf, len(pts_withOutf), np.transpose(pt2add), axis=0)

    # Interpolate points
    pts_int_opt = get_interpolated_pts(points=pts_all_opt, nPoints = nPoints)
    
    return pts_int_opt, pt2add, sph_m10

#Thickness / Ballooning
def extract_ballooning(organ, name, name_cl, setup):#
   
    ch, cont = name
    mesh2ball = organ.obj_meshes[ch+'_'+cont]
    from_cl, from_cl_type = name_cl
    cl4ball = organ.obj_meshes[from_cl+'_'+from_cl_type].get_centreline()
    sph4ball = sphs_in_spline(kspl=cl4ball,every=0.6)
    sph4ball.legend('sphs_ball').alpha(0.1)
    
    mesh_ball, distance, min_max = get_distance_to(mesh_to = mesh2ball, 
                                                   mesh_from = sph4ball, 
                                                   from_name ='CL('+from_cl+'_'+from_cl_type+')',
                                                   range = (setup['min_val'], setup['max_val']),
                                                   color_map = setup['colormap'])
    mesh_ball.alpha(1)
    #Add min-max values to mH_settings
    proc_range = ['measure', 'ball', ch+'_'+cont+'_('+from_cl+'_'+from_cl_type+')']
    upd_range = {'range_o':{'min_val': min_max[0], 'max_val': min_max[1]}, 
                 'range_user': {'min_val': setup['min_val'], 'max_val': setup['max_val']}, 
                 'colormap': setup['colormap']}
    organ.update_settings(proc_range, update = upd_range, mH = 'mH')

    #Add mesh_ball to the mesh_meas attribute
    m_type = 'ballCL('+from_cl+'_'+from_cl_type+')'
    
    if not hasattr(mesh2ball, 'mesh_meas'):
        mesh2ball.mesh_meas = {}
    mesh2ball.mesh_meas[m_type] = mesh_ball
    mesh2ball.save_mesh(m_type=m_type)
    mesh2ball.save_array(array=distance, m_type=m_type)
    
    # Update organ workflow
    cont_proc = cont+'_('+from_cl+'_'+from_cl_type+')'
    proc_wft = ['MeshesProc','D-Ballooning',ch, cont_proc, 'Status']
    organ.update_mHworkflow(process = proc_wft, update = 'DONE')

def get_thickness(organ, name, thck_dict, setup):#
    
    ch, cont = name
    mesh_tiss = organ.obj_meshes[ch+'_tiss']
    if thck_dict['n_type'] == 'int>ext': 
        mesh_to = organ.obj_meshes[ch+'_ext']
        mesh_from = organ.obj_meshes[ch+'_int'].mesh
    else:# n_type == 'ext>int': 
        mesh_to = organ.obj_meshes[ch+'_int']
        mesh_from = organ.obj_meshes[ch+'_ext'].mesh

    mesh_thck, distance, min_max = get_distance_to(mesh_to = mesh_to, 
                                                   mesh_from = mesh_from, 
                                                   from_name = thck_dict['n_type'],  
                                                   range = (setup['min_val'], setup['max_val']),
                                                   color_map = setup['colormap'])
    mesh_thck.alpha(1)
    # Add mesh_ball to the mesh_meas attribute
    n_type_new = thck_dict['n_type'] .replace('>','TO')
    proc_range = ['measure', thck_dict['short'], ch+'_tiss_whole']
    upd_range = {'range_o':{'min_val': min_max[0], 'max_val': min_max[1]}, 
                 'range_user': {'min_val': setup['min_val'], 'max_val': setup['max_val']},
                 'colormap': setup['colormap']}
    organ.update_settings(proc_range, update = upd_range, mH = 'mH')
    
    #Add mesh_thck to the mesh_meas attribute
    m_type = 'thck('+n_type_new+')' #thck(intTOext)

    if not hasattr(mesh_tiss, 'mesh_meas'):
        mesh_tiss.mesh_meas = {}
    mesh_tiss.mesh_meas[m_type] = mesh_thck
    mesh_tiss.save_mesh(m_type=m_type)
    mesh_tiss.save_array(array=distance, m_type=m_type)
    
    # Update organ workflow
    proc_wft = ['MeshesProc', thck_dict['method'], ch, cont, 'Status']
    organ.update_mHworkflow(process = proc_wft, update = 'DONE')

def get_distance_to(mesh_to, mesh_from, from_name, range, color_map='turbo'):#
    """
    Function that gets the distance between m_ext and m_int and color codes m_ext accordingly
    mesh_from = cl/int
    mesh_to = ext - the one that has the colors
    """
    if isinstance(mesh_from,vedo.shapes.Spheres): 
        name = 'CL to Tissue'
        multip = 1
        mesh_name = mesh_to.legend
        suffix = '_CL_to_Tissue'
        title = mesh_name+'\n'+name+' [um]\n'+from_name
    else: 
        name = 'Thickness'
        if from_name == 'int>ext':
            multip = -1
        else: 
            multip = 1
        mesh_name = mesh_to.legend.split('_')
        if len(mesh_name)>2:
            mesh_name = mesh_to.legend.replace('_'+mesh_name[-1],'')
        else: 
            mesh_name = mesh_name[0]
        suffix = '_Thickness'
        title = mesh_name+'\n'+name+' [um]\n('+from_name+')'
    
    tic = perf_counter()
    print('> Start time: \t', str(datetime.now())[11:-7])
    mesh_toB = mesh_to.mesh.clone()
    mesh_toB.distance_to(mesh_from, signed=True)
    toc = perf_counter()
    time = toc-tic
    print('> End time: \t', str(datetime.now())[11:-7])
    print("\t>> Total time taken to get "+name+" = ", format(time/60,'.2f'), "min")
    distance = multip*mesh_toB.pointdata['Distance']
    mesh_toB.pointdata['Distance'] = distance
    vmin, vmax = np.min(distance),np.max(distance)
    min_max = (vmin, vmax)
    print('> vmax:', format(vmax,'.2f'), 'vmin:', format(vmin,'.2f'))
    alert('jump')
    
    mesh_toB.cmap(color_map)
    mesh_toB.add_scalarbar(title=title, pos=(0.8, 0.05))
    if range == (None, None): 
        mesh_toB.mapper().SetScalarRange(vmin,vmax)
    else: 
        mesh_toB.mapper().SetScalarRange(range[0],range[1])
    mesh_toB.legend(mesh_name+suffix)
    
    return mesh_toB, distance, min_max

#Meshes functions
def s3_to_mesh(s3, res, name:str, color='cyan', rotateZ_90=True):

    verts, faces, _, _ = measure.marching_cubes(s3, spacing=res, method='lewiner')
    
    mesh = vedo.Mesh([verts, faces])
    if rotateZ_90: 
        mesh.rotateZ(-90).wireframe(True)
    mesh = mesh.extract_largest_region()
    alert('frog')
    mesh.color(color).alpha(1).wireframe().legend(name)

    return mesh

#Segments functions
def sphs_in_spline(kspl, colour=False, color_map='turbo', every=10):
    """
    Function that creates a group of spheres through a spline given as input.
    """
    if colour:
        # scalars->colors
        vcols = [vedo.color_map(v, color_map, 0, len(kspl.points())) for v in list(range(len(kspl.points())))]  
    
    if every > 1:
        spheres_spline = []
        for num, point in enumerate(kspl.points()):
            if num % every == 0 or num == kspl.npoints-1:
                if num < 100:
                    size = (0.02,0.02)
                else: 
                    size = (0.03,0.02)
                if colour:
                    sphere_pt = vedo.Sphere(pos=point, r=2, c=vcols[num]).addScalarBar(title='Centreline\nPoint Number')
                    sphere_pt.legend('sph.'+str(num))
                    sphere_pt.caption(str(num), point=point,
                              size=size, padding = 3, justify='center', font=txt_font, alpha=0.8)
                else:
                    sphere_pt = vedo.Sphere(pos=point, r=2, c='coral')
                spheres_spline.append(sphere_pt)
    else:
        kspl_new = vedo.KSpline(kspl.points(), res = round(kspl.npoints/every))
        spheres_spline = vedo.Spheres(kspl_new.points(), c='coral', r=2)

    return spheres_spline

def create_disc_mask(organ, cut, s3_shape, h_min = 0.1125): 
    
    #Load stack shape
    xdim, ydim, zdim = s3_shape
    print('s3_shape:',s3_shape, '- xdim:', xdim, '- ydim:',ydim, '- zdim:', zdim)
        
    disc_info = organ.mH_settings['wf_info']['segments']['setup'][cut]['cut_info']
    for disc in disc_info:
        print('Creating mask for '+ disc)
        disc_name = disc.replace(' ', '').replace('.','')
        
        disc_radius = disc_info[disc]['radius']
        normal_unit = disc_info[disc]['normal_unit']
        pl_centre = disc_info[disc]['pl_centre']
        height = disc_info[disc]['height']
        res = disc_info[disc]['res']
        
        # Create a disc with better resolution to transform into pixels to mask stack
        res_cyl = 2000
        num_rad = int(3*int(disc_radius)) #int(((r_circle_max-r_circle_min)/0.225)+1)
        num_h = 9
        for j, rad in enumerate(np.linspace((disc_radius*0.2)/2, disc_radius, num_rad)):
            for i, h in enumerate(np.linspace(h_min,height, num_h)):
                cyl = vedo.Cylinder(pos = pl_centre, r = rad, height = h, axis = normal_unit, c = 'lime', cap = True, res = res_cyl)#.wireframe(True)
                if i == 0 and j == 0:
                    cyl_pts = cyl.points()
                else:
                    cyl_pts = np.concatenate((cyl_pts, cyl.points()))

        # Rotate the points that make up the HR disc, to convert them to a stack
        cyl_points_rot = np.zeros_like(cyl_pts)
        if organ.mH_settings['setup']['rotateZ_90']: #'CJ' not in filename: 
            axis = [0,0,1]
        else: 
            axis = [0,0,0]
    
        for i, pt in enumerate(cyl_pts):
            cyl_points_rot[i] = (np.dot(rotation_matrix(axis = axis, theta = np.radians(90)),pt))
    
        cyl_pix = np.transpose(np.asarray([cyl_points_rot[:,i]//res[i] for i in range(len(res))]))
        cyl_pix = cyl_pix.astype(int)
        cyl_pix = np.unique(cyl_pix, axis=0)
        
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
        s3_cyl = s3_cyl.astype('uint8')
        
        name2save = organ.user_organName + '_mask_'+cut+'_'+disc_name+'.npy'
        dir2save = organ.dir_res(dir ='s3_numpy') / name2save
        np.save(dir2save, s3_cyl)
        
        if not dir2save.is_file():
            print('>> Error: s3 mask of '+disc+' was not saved correctly!\n>> File: mask_'+disc_name+'.npy')
            alert('error_beep')
        else: 
            print('>> s3 mask of '+disc+' saved correctly!')
            alert('countdown')

def segm_ext_ext(organ, mesh, cut, segm_names, palette, win):

    #Mask the s3_stack to create segments
    cut_masked = mesh.mask_segments(cut = cut, check=True)
    if cut_masked != None: 
        #Create a dictionary containing the information of the classified segments 
        dict_segm, colors = organ.dict_segments(cut, palette)
        #Ask user to classify the segments using the interactive plot
        dict_segm = classify_segments(meshes=cut_masked, dict_segm=dict_segm, 
                                    colors_dict=colors)
        # print('dict_segm after classification:', dict_segm)

        meshes_segm = {}; final_subsgm = {}; ext_subsgm_names = {}
        #Create submeshes for the input mesh
        for n, segm, color in zip(count(), segm_names, palette):
            print('\t- Creating segment No.',n, ' for ', cut, '.', segm, color)
            sp_dict_segm = dict_segm[segm]
            # subsgm: instance of class Submesh
            # final_segm_mesh: instance of mesh (vedo)
            subsgm = mesh.create_segment(name = segm, cut = cut, color = color)
            final_segm_mesh = create_subsegment(organ, subsgm, cut, cut_masked,
                                                    'segm', sp_dict_segm, color)
            save_submesh(organ, subsgm, final_segm_mesh, win)
            final_subsgm[segm]=subsgm
            meshes_segm[segm] = final_segm_mesh
            ext_subsgm_names[segm] = subsgm.sub_name_all

        #Add ext_subsm names to load them in the future
        proc_set = ['wf_info', 'segments', 'setup', cut, 'names']
        organ.update_settings(proc_set, ext_subsgm_names, 'mH')
        print('wf_info (segm)-after: ', get_by_path(organ.mH_settings, ['wf_info', 'segments']))
        print('organ.submeshes:', organ.submeshes)

        #final_subsgm:dict of class submesh 
        #meshes_segm: dict of final meshes
        return final_subsgm, meshes_segm
    else: 
        return None, None

def get_segments(organ, mesh, cut, segm_names, palette, ext_subsgm, win=None): 
        
    #Mask the s3_stack to create segments
    cut_masked = mesh.mask_segments(cut = cut)
    #Create a dictionary containing the information of the classified segments 
    dict_segm = organ.dict_segments(cut, other=False)

    #Create submeshes of the input mesh 
    meshes_segm = {}
    for n, segm, color in zip(count(), segm_names, palette):
        sp_dict_segm = classify_segments_from_ext(meshes = cut_masked, 
                                                    dict_segm = dict_segm[segm],
                                                    ext_sub = ext_subsgm[segm])
        # print('dict_segm after classif: ', sp_dict_segm)
        #Create submesh - segment
        subsgm = mesh.create_segment(name = segm, cut = cut, color = color)
        final_segm_mesh = create_subsegment(organ, subsgm, cut, cut_masked, 
                                                    'segm', sp_dict_segm, color)
        if organ.mH_settings['wf_info']['segments']['save_all']:
            save_submesh(organ, subsgm, final_segm_mesh, win)
        meshes_segm[segm] = final_segm_mesh
        
    # print('\n\n FINAL! mH_settings[measure]: ', organ.mH_settings['measure'], '\n\n')
    # print('wf_info (segm)-after: ', get_by_path(organ.mH_settings, ['wf_info', 'segments']))
    # print('organ.submeshes:', organ.submeshes)
    # print('organ.obj_temp:', organ.obj_temp)

    return meshes_segm

def classify_segments(meshes, dict_segm, colors_dict):
    
    flat_segm = flatdict.FlatDict(dict_segm)
    colors = [colors_dict[key] for key in colors_dict]
    mesh_classif = [flat_segm[key] for key in flat_segm if 'meshes_number' in key]
    names = [flat_segm[key] for key in flat_segm if 'user_name' in key]
    
    #https://seaborn.pydata.org/tutorial/color_palettes.html
    palette = sns.color_palette("Set2", len(meshes))
    for n, mesh in enumerate(meshes):
        mesh.color(palette[n])
    
    def func(evt):
        if not evt.actor:
            return
        mesh_no = evt.actor.info['legend']
        mesh_color = evt.actor.color()
        # print('mesh_color:',mesh_color)
        mesh_color = evt.actor.color()*255
        is_in_list = np.any(np.all(mesh_color == colors, axis=1))
        if is_in_list:
            bool_list = np.all(mesh_color == colors, axis=1)
            ind = np.where(bool_list == True)[0][0]
            ind_plus1 = ind+1
            new_ind = (ind_plus1)%len(colors)
            new_color = colors[new_ind]
            evt.actor.color(new_color)
            mesh_classif[ind].remove(mesh_no)
            mesh_classif[new_ind].append(mesh_no)
            msg.text("You clicked Mesh "+mesh_no+' to classify it as '+names[new_ind])
        else: 
            evt.actor.color(colors[1])
            mesh_classif[1].append(mesh_no)
            msg.text("You clicked Mesh "+mesh_no+' to classify it as '+names[1])
        print(mesh_classif)
        
    mks = []; sym = ['o']*len(dict_segm)
    for segm in dict_segm:
        mks.append(vedo.Marker('*').c(colors_dict[segm]).legend(dict_segm[segm]['user_name']))
    
    txA = 'Instructions: Click each segment mesh until it is coloured according to the segment it belongs to.'
    txB = '\n[Note: Colours will loop for each mesh as you click it ]'
    txt0 = vedo.Text2D(txA+txB, c=txt_color, font=txt_font, s=txt_size)
    lb = vedo.LegendBox(mks, markers=sym, font=txt_font, 
                        width=leg_width/1.5, height=leg_height/1.5)
    
    msg = vedo.Text2D("", pos="bottom-center", c=txt_color, font=txt_font, s=txt_size, bg='red', alpha=0.2)
    vpt = vedo.Plotter(axes=5, bg='white')
    vpt.add_icon(logo, pos=(0.1,1), size=0.25)
    vpt.add_callback('mouse click', func)
    vpt.show(meshes, lb, msg, txt0, zoom=1.2)
    
    return dict_segm
    
def classify_segments_from_ext(meshes, dict_segm, ext_sub):
    
    ext_sub_mesh = ext_sub.get_segm_mesh()
    list_meshes = []
    for mesh in meshes: 
        if isinstance(mesh, vedo.Mesh):
            com = mesh.center_of_mass()
            if ext_sub_mesh.is_inside(com):
                list_meshes.append(mesh.info['legend'])
    dict_segm['meshes_number'] = list_meshes
            
    return dict_segm

#Sections/Regions Functions
def get_stack_clRibbon(organ, s3_shape, mesh_cl, cl_ribbon, win):
    
    res = mesh_cl.resolution
    cl_ribbonR = cl_ribbon.clone().x(res[0])
    cl_ribbonF = cl_ribbon.clone().y(res[1])
    cl_ribbonT = cl_ribbon.clone().z(res[1])
    cl_ribbonS = [cl_ribbon, cl_ribbonR, cl_ribbonF, cl_ribbonT]

    #Load stack shape
    xdim, ydim, zdims = s3_shape
    zdim = zdims-2

    # Rotate the points that make up the cl_ribbon, to convert them to a stack
    cust_angle = organ.info['custom_angle']
    print('cust_angle:', cust_angle)
    
    if organ.mH_settings['setup']['rotateZ_90']: #'CJ' not in filename: 
        axis = [0,0,1]
        theta = np.radians(90)
    else: 
        axis = [0,0,0]
        theta = np.radians(0)
    print('axis:',axis, '- theta:',theta)
    
    if int(cust_angle) != 0: 
        print('cust_angle != 0, develop this!')

    s3_rib = np.zeros((xdim, ydim, zdim+2))
    s3_filledCube = np.zeros((xdim, ydim, zdim+2))
    s3_filledCube[1:xdim-1,1:ydim-1,1:zdim+1] = 1

    win.prog_bar_range(1,16)
    n = 0
    # Rotate the points in all ribbons and fit into stack size 
    for cl_rib in cl_ribbonS:
        rib_pts = cl_rib.points()
        rib_points_rot = np.zeros_like(rib_pts)
        for i, pt in enumerate(rib_pts):
            rib_points_rot[i] = (np.dot(rotation_matrix(axis = axis, theta = theta),pt))
        n+=1; win.prog_bar_update(n)
        rib_pix = np.transpose(np.asarray([rib_points_rot[:,i]//res[i] for i in range(len(res))]))
        rib_pix = rib_pix.astype(int)
        rib_pix = np.unique(rib_pix, axis=0)
        n+=1; win.prog_bar_update(n)
        rib_pix_out = rib_pix.copy()
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
                index_out.append(index)
        n+=1; win.prog_bar_update(n)
        rib_pix_out = np.delete(rib_pix_out, index_out, axis = 0)
        # Create mask of cl_ribbon
        s3_rib[rib_pix_out[:,0],rib_pix_out[:,1],rib_pix_out[:,2]] = 1
        # Create filled cube just to one side of cl
        s3_filledCube[rib_pix_out[:,0],rib_pix_out[:,1],rib_pix_out[:,2]] = 0
        n+=1; win.prog_bar_update(n)
    win.prog_bar_update(16)
    alert('clown')
        
    # Create volume of extended centreline mask
    test_rib = s3_to_mesh(s3_rib, res=res, name='Extended CL', color='darkmagenta')
    
    return s3_filledCube, test_rib

def get_cube_clRibbon(organ, s3_shape, cut, s3_filledCube, res, pl_normal, filling_method):

    # #Load stack shape
    xdim, ydim, zdims = s3_shape
    zdim = zdims-2

    #Identify the direction in which the cubes need to be built
    if organ.mH_settings['wf_info']['sections'][cut.title()]['axis_lab'].lower() == 'roi':
        ref_vect = [0, 1, 0]# [[0, 1, 0], [0, 0, 0]]
    else: 
        ref_vect = organ.mH_settings['wf_info']['orientation']['roi']['settings']['ref_vectF'][0]
    
    ext_dir = list(np.cross(ref_vect, pl_normal))
    ext_dir_abs = np.absolute(ext_dir)
    ext_dir_absf = list(ext_dir_abs)
    max_ext_dir = max(ext_dir_abs)
    coord_dir = ext_dir_absf.index(max_ext_dir)
    methodA = False; methodB = False; methodC = False
    if coord_dir == 0: 
        method_def = 'methodA'
        if filling_method == 'No.1':
            methodA = True
        elif filling_method == 'No.2': 
            methodB = True
        else: #'No.3'
            methodC = True
    else:# coord_dir == 1 or coord_dir == 2: 
        method_def = 'methodB'
        if filling_method == 'No.1':
            methodB = True
        elif filling_method == 'No.2': 
            methodC = True
        else: #'No.3'
            methodA = True

    print('pl_normal:', pl_normal, '\nref_vect:', ref_vect, '\next_dir:', ext_dir, '\ncoord_dir:', coord_dir)
    print('method_def: ', method_def, '>> Selected: A: ', methodA, '- B: ', methodB, '- C: ', methodC)
    
    if methodA: #coord_dir == 0: 
        print('Extending cube in the x? direction')
        for xpos in range(0,xdim):
            for zpos in range(0,zdim+2): 
                yline = s3_filledCube[xpos,0:ydim,zpos]
                index_y = np.where(yline == 0)[0]
                index_y = list(index_y)
                index_y.pop(0);index_y.pop(-1)
                if len(index_y) > 0:
                    s3_filledCube[xpos,index_y[0]:ydim,zpos] = 0
        
    elif methodB: #coord_dir == 1 or coord_dir == 2:
        print('Extending cube in the y direction - check!!!')
        for ypos in range(0,ydim):
            for zpos in range(0,zdim+2): 
                xline = s3_filledCube[0:xdim,ypos,zpos]
                index_x = np.where(xline == 0)[0]
                index_x = list(index_x)
                index_x.pop(0);index_x.pop(-1)
                if len(index_x) > 0:
                    s3_filledCube[index_x[-1]:xdim,ypos,zpos] = 0

    elif methodC:# coord_dir == 2: 
        print('Extending cube in the z direction')
        for xpos in range(0,xdim):
            for ypos in range(0,ydim): 
                zline = s3_filledCube[xpos,ypos,0:zdim+2]
                index_z = np.where(zline == 0)[0]
                index_z = list(index_z)
                index_z.pop(0);index_z.pop(-1)
                if len(index_z) > 0:
                    s3_filledCube[xpos,ypos,index_z[0]:zdim+2] = 0

    alert('woohoo')
    
    #Create volume of filled side of extended centreline mask
    mask_cube = s3_to_mesh(s3_filledCube, res=res, name='Filled CLRibbon SideA', color='darkblue')
    mask_cube.alpha(0.05)

    s3_filledCube = s3_filledCube.astype('uint8')

    #Create the inverse section and make the user select the section that corresponds to section 1
    s3_filledCubeBoolA = copy.deepcopy(s3_filledCube).astype(bool)
    s3_filledCubeBoolB = np.full_like(s3_filledCubeBoolA, False)
    s3_filledCubeBoolB[1:-1,1:-1,1:-1] = True
    s3_filledCubeBoolB[s3_filledCubeBoolA] = False

    mask_cubeB = s3_to_mesh(s3_filledCubeBoolB, res=res, name='Filled CLRibbon SideB', color='skyblue')
    mask_cubeB.alpha(0.05).linewidth(1)
    mask_cube.linewidth(1)
    
    mask_cube.color('darkblue'); mask_cubeB.color('skyblue')
    mask_cube_split = [mask_cube, mask_cubeB]

    s3_filledCubes = {'SideA': s3_filledCube, 'SideB': s3_filledCubeBoolB}

    return mask_cube_split, s3_filledCubes
    
def select_ribMask(organ, cut, mask_cube_split, mesh_cl, kspl_ext): 
    
    dict_sides = {}
    for mesh in mask_cube_split:
        dict_sides[mesh.info['legend']] = {'color': mesh.color()*255}

    # dict_sides = copy.deepcopy(dict_sides)
    flat_segm = flatdict.FlatDict(dict_sides)
    colors = [flat_segm[key] for key in flat_segm if 'color' in key]
    color_sel = [210,105,30] # yellow [255,215,0] # chocolate (210,105,30)
    
    def func(evt):
        if not evt.actor:
            return
        mesh_no = evt.actor.info['legend']
        if 'Ribbon' in mesh_no: 
            mesh_color = evt.actor.color()*255
            msg.text(mesh_no +' has been selected as '+name_sect1.upper())
            is_in_list = np.any(np.all(mesh_color == colors, axis=1))
            if is_in_list: 
                new_color = color_sel
                evt.actor.color(new_color)
                selected_mesh['name'] = mesh_no
                for mesh in mask_cube_split:
                    if mesh.info['legend'] != mesh_no:
                        # print('running for mesh: ', mesh.info['legend'])
                        mesh.color(dict_sides[mesh.info['legend']]['color'])
            else: 
                new_color = color_sel
                evt.actor.color(new_color)
                selected_mesh['name'] = mesh_no
                for mesh in mask_cube_split:
                    if mesh.info['legend'] != mesh_no:
                        mesh.color(dict_sides[mesh.info['legend']]['color'])
        else: 
            print('Other mesh was selected')
    
    name_sect1 = organ.mH_settings['setup']['sect'][cut.title()]['name_sections']['sect1']
    mks = [vedo.Marker('*').c(color_sel).legend('Section No.1 ('+name_sect1+')')]
    sym = ['o']
    lb = vedo.LegendBox(mks, markers=sym, font=txt_font, 
                        width=leg_width, height=leg_height/3)

    msg = vedo.Text2D("", pos="bottom-center", c=txt_color, font=txt_font, s=txt_size, bg='red', alpha=0.2)
    selected_mesh = {'name': 'NS'}
    
    txA = 'Instructions '+cut.title()+': Select the mesh that corresponds to Section No.1 ('+name_sect1.upper()+').\nClose the window when you are done.\nNote: If the selection is not working you will be prompted with a message to select the mesh after closing the window,\n just remember the color of the mesh corresponding to the '+name_sect1.upper()+' side.'
    txt0 = vedo.Text2D(txA, c=txt_color, font=txt_font, s=txt_size)
    
    vpt = vedo.Plotter(axes=1)
    vpt.add_icon(logo, pos=(0.1,1), size=0.25)
    vpt.add_callback('mouse click', func)
    vpt.show(mask_cube_split, mesh_cl, kspl_ext, txt0, lb, msg, azimuth=45, elevation=20, zoom=0.8, interactive=True)
    print('selected_mesh:', selected_mesh)

    return selected_mesh

def save_ribMask_side(organ, cut, selected_side, s3_filledCubes):

    if 'SideA' in selected_side['name']:
        mask = s3_filledCubes['SideA']
    else: #'SideB' in sel_side['name']:
        mask = s3_filledCubes['SideB']

    mask_selected = mask.astype('uint8')
    name2save = organ.user_organName +'_mask_'+cut.title()+'_sect1.npy'
    dir2save = organ.dir_res(dir ='s3_numpy') / name2save
    np.save(dir2save, mask_selected)
    
    proc = ['wf_info', 'sections', cut.title(), 'mask_name']
    organ.update_settings(proc, update = name2save, mH = 'mH')
    print(organ.mH_settings['wf_info']['sections'])
    setattr(organ, 'mask_sect_'+cut.lower(), name2save)

    if not dir2save.is_file():
        print('>> Error: s3 mask of sections was not saved correctly!\n>> File: mask_sect.npy')
        alert('error_beep')
    else: 
        print('>> s3 mask of sections saved correctly!')
        alert('countdown')
    
def get_sections(organ, mesh, cut, sect_names, palette, win):
    
    #Create submeshes of the input mesh 
    meshes_sect = {}
    for n, sect, color in zip(count(), sect_names, palette):
        print(mesh, cut, sect, color)
        subsct, final_sect_mesh = create_subsection(organ, mesh, cut, 
                                                    sect, color)
        meshes_sect[sect] = final_sect_mesh
    print('\n\n FINAL! mH_settings[measure]: ', organ.mH_settings['measure'], '\n\n')

    print('wf_info (sect)-after: ', get_by_path(organ.mH_settings, ['wf_info', 'sections']))
    print('organ.submeshes:', organ.submeshes)

    return meshes_sect

def create_subsection(organ, mesh, cut, sect, color):

    #Create submesh - section
    subsct = mesh.create_section(name = sect, cut = cut, color = color)
    print(subsct.__dict__)
    final_sect_mesh = subsct.get_sect_mesh()
    #Add submesh to organ
    organ.add_submesh(subsct)
    #Get segm measurements
    measure_submesh(organ, subsct, final_sect_mesh)

    # Update organ workflow
    cut, ch, cont, sect = subsct.sub_name_all.split('_')
    proc_wft = ['MeshesProc', 'E-Sections', cut, ch, cont, 'Status']
    organ.update_mHworkflow(process = proc_wft, update = 'DONE')

    return subsct, final_sect_mesh

def get_segm_sects(organ, ch_cont, cuts, names, palette, win):
    
    ch, cont = ch_cont.split('_')
    seg_cut, _ = cuts
    # cutsf = 's'+seg_cut+'_'+reg_cut
    segm_names, sect_names = names

    #Get the segment's method 
    method = organ.mH_settings['wf_info']['segments']['setup'][seg_cut]['ch_info'][ch][cont]
    if method in ['ext-ext', 'cut_with_ext-ext', 'cut_with_other_ext-ext']: 
        try:
            ext_subsgm = organ.ext_subsgm
            print('try ext_subsgm')
        except: 
            if seg_cut[0] == 's': 
                seg_cut = seg_cut[1:]
            ext_subsgm = organ.get_ext_subsgm(seg_cut)
            print('except ext_subsgm')
        print('ext_subsgm: ',ext_subsgm)
    else: 
        print('This needs to be coded!')

    #Create submeshes of the input mesh 
    meshes_segm_sect = {}
    meshes_final = []
    for sect in sect_names:
        palette_dict = {}
        for key in palette.keys():
            if sect in key:
                new_key = key.split('_')[0]
                palette_dict[new_key] = palette[key]
        print('palette_dict', palette_dict)
        meshes_segm_sect[sect] = {}
        # palette_sect = [palette_dict[key] for key in palette_dict.keys()]
        for segm in segm_names: 
            segm_sect = segm+'_'+sect
            color = palette[segm+'_'+sect]
            print('Dividing into Segment-Section Intersections '+segm_sect)
            final_segm_sect_mesh = create_sub_segm_sect(organ, cuts, ch_cont, 
                                                        segm, sect, 
                                                        color, palette_dict, 
                                                        ext_subsgm)
            meshes_segm_sect[sect][segm] = final_segm_sect_mesh

            print('organ.submeshes:',organ.submeshes.keys())
            meshes_final.append(final_segm_sect_mesh)

    # vp = vedo.Plotter(N=1, axes=4)
    # vp.show(meshes_final, at = 0, interactive = True)
        
    print('\n\n FINAL! mH_settings[measure]: ', organ.mH_settings['measure'], '\n\n')
    print('organ.submeshes:', organ.submeshes)
    alert('woohoo')

    return meshes_segm_sect

def create_sub_segm_sect(organ, cuts, ch_cont, segm, sect, color, palette_dict, ext_subsgm): 

    seg_cut, reg_cut = cuts
    cutsf = 's'+seg_cut+'_'+reg_cut
    #Extract segment
    segm_sect = segm+'_'+sect
    obj_segm = organ.obj_subm[seg_cut+'_'+ch_cont+'_'+segm]
    #Create subsegm_subsect
    submesh = obj_segm.create_segm_sect(segm_sect = segm_sect, cuts = cutsf, color = color)

    #SECTION/REGION
    #Get the section mask containing only this section
    sect_mask = submesh.get_sect_mesh(output='mask')
    #SEGMENT
    #Mask the segments providing as input the masked section
    cut_masked = submesh.mask_segments(cut = seg_cut, s3 = sect_mask)
    #Create a dictionary containing the information of the classified segments 
    dict_segm = organ.dict_segments(seg_cut, other=False)
    #Classify the resulting segments using ext mesh
    sp_dict_segm = classify_segments_from_ext(meshes = cut_masked, 
                                                dict_segm = dict_segm[segm],
                                                ext_sub = ext_subsgm[segm])
    print('dict_segm after classif: ', sp_dict_segm)
    final_segm_sect_mesh = create_subsegment(organ, submesh, seg_cut, cut_masked, 
                                            'segm-sect', sp_dict_segm, color)
    
    # Update organ workflow
    ch, cont = ch_cont.split('_')
    proc_wft = ['MeshesProc', 'E-Segments_Sections', 's'+seg_cut, reg_cut, ch, cont, 'Status']
    organ.update_mHworkflow(process = proc_wft, update = 'DONE')

    return final_segm_sect_mesh
    
#General SubMeshes
def save_submesh(organ, submesh, mesh, win, ext='.vtk'):
    
    mesh_name = organ.user_organName+'_'+submesh.sub_name_all+ext
    mesh_dir = organ.dir_res(dir ='meshes') / mesh_name
    mesh.write(str(mesh_dir))
    
    #Update mH_settings with directory
    cut, ch, cont, segm = submesh.sub_name_all.split('_')
    proc_set = ['wf_info', 'segments', 'setup', cut, 'dirs', ch, cont]
    organ.update_settings(proc_set+[segm], mesh_name, 'mH')

    print('wf_info (segm)-after: ', get_by_path(organ.mH_settings, ['wf_info', 'segments']))
    print('>> Mesh '+mesh_name+' has been saved!')
    if win != None: 
        win.win_msg('Mesh '+mesh_name+' has been saved!')
    alert('countdown')        

def create_subsegment(organ, subsgm, cut, cut_masked, stype, sp_dict_segm, color):#
    
    #Assign meshes and measure them
    list_meshes = sp_dict_segm['meshes_number']
    #Find meshes and add them to the mesh
    mesh_segm = []; 
    print('len(cut_masked) -'+subsgm.sub_name_all+': ', len(cut_masked))
    for cut_mesh in cut_masked:
        if isinstance(cut_mesh, vedo.Mesh):
            if cut_mesh.info['legend'] in list_meshes: 
                mesh_segm.append(cut_mesh)

    if len(mesh_segm) == 0: 
        print('Something went wrong and no meshes were recognised as part of this segment: ', subsgm.sub_name_all)
        final_segm_mesh = []
    else: 
        #Merge all meshes and apply final settings
        final_segm_mesh = vedo.merge(mesh_segm)
        final_segm_mesh.color(color).alpha(0.1).legend(subsgm.sub_legend)
        subsgm.color = color
        #Add submesh to organ
        organ.add_submesh(subsgm)
        #Get measurements
        measure_submesh(organ, subsgm, final_segm_mesh)
    
        # Update organ workflow
        if stype == 'segm': 
            cut, ch, cont, segm = subsgm.sub_name_all.split('_')
            proc_wft = ['MeshesProc', 'E-Segments', cut, ch, cont, 'Status']
            organ.update_mHworkflow(process = proc_wft, update = 'DONE')
        else: 
            pass

    # final_segm_mesh: instance of mesh (vedo)
    return final_segm_mesh

def create_iso_volume(organ, ch, plot=True):

    #Open image for channel btn
    img_dir = organ.img_dirs[ch]['image']['dir']
    color = organ.mC_settings['setup']['color_chs'][ch]
    res = organ.get_resolution()
    threshold = organ.mC_settings['wf_info']['isosurface'][ch]['threshold']
    alpha = organ.mC_settings['wf_info']['isosurface'][ch]['alpha']

    # Invert stack
    if organ.mC_settings['wf_info']['isosurface']['invert']: 
        res_copy = copy.deepcopy(res)
        if res_copy[1] > 0:
            res_copy[1] = -res_copy[1]
        print('res_copy:', res_copy)

    try: 
        vol = vedo.Volume(str(img_dir), spacing = res_copy).isosurface(int(threshold)).alpha(float(alpha))
    except:
        return False
    
    vol.color(color)

    if hasattr(organ, 'vol_iso'): 
        organ.vol_iso[ch] = vol
    else: 
        organ.vol_iso = {ch: vol}

    if plot: 
        txtf = organ.user_organName+' \n - Isosurface volume reconstruction for Ch'+ch[-1]+': '+organ.mC_settings['setup']['name_chs'][ch]+' [threshold: '+str(threshold)+']'
        txt = vedo.Text2D(txtf,  c=txt_color, font=txt_font, s=txt_size)
        vp = vedo.Plotter(N=1, axes = 1)
        vp.add_icon(logo, pos=(0.1,0.0), size=0.25)
        vp.show(vol, txt, at=0, interactive=True)

    return True

def remove_cells(organ, vol_iso): 

    vols = []; names = []; colors = []
    for ch in vol_iso: 
        vols.append(vol_iso[ch])
        names.append(organ.mC_settings['setup']['name_chs'][ch])
        colors.append(organ.mC_settings['setup']['color_chs'][ch])

    sphs = organ.cellsMC['chA'].cells
    cells2remove = []
    silcont = [None]
    def remove_cells(evt):
        if not evt.actor: return
        if isinstance(evt.actor, vedo.shapes.Sphere): 
            sil = evt.actor.silhouette().lineWidth(6).c('red')
            print("You clicked: "+evt.actor.name)
            msg.text("You clicked: Sphere "+evt.actor.name)
            cell_no = evt.actor.name
            cell_no = int(cell_no.split('.')[-1])
            if cell_no < 1000000: # remove
                cells2remove.append(cell_no)
                sphs[cell_no].color('black')
                new_no = 1000000+cell_no
                sphs[cell_no].name = f"Cell Nr.{new_no}"
            else: # add back
                ind2rem = cells2remove.index(cell_no-1000000)
                cells2remove.pop(ind2rem)
                old_no = cell_no-1000000
                sphs[old_no].color('gold')
                sphs[old_no].name = f"Cell Nr.{old_no}"
                
            plt.remove(silcont.pop()).add(sil)
            silcont.append(sil)

    if len(vol_iso) >= 1: 
        def sliderAlphaMeshOut(widget, event):
            valueAlpha = widget.GetRepresentation().GetValue()
            vols[0].alpha(valueAlpha)
    if len(vol_iso) >= 2: 
        def sliderAlphaMeshOut2(widget, event):
            valueAlpha = widget.GetRepresentation().GetValue()
            vols[1].alpha(valueAlpha)
    if len(vol_iso) >= 3: 
        def sliderAlphaMeshOut3(widget, event):
            valueAlpha = widget.GetRepresentation().GetValue()
            vols[2].alpha(valueAlpha)
    
    txt_slider_size2 = 0.7
    
    msg = vedo.Text2D("", pos="bottom-center", c=txt_color, font=txt_font, s=txt_size, bg='red', alpha=0.2)
    txtf = organ.user_organName+' - Removing cells \n -Click in the cells you would like to remove from the analysis, and close the window when done.' 
    txt = vedo.Text2D(txtf,  c=txt_color, font=txt_font, s=txt_size)
    plt = vedo.Plotter(axes=1)
    plt.add_icon(logo, pos=(0.1,0.0), size=0.25)
    if len(vol_iso) >= 1: 
        plt.addSlider2D(sliderAlphaMeshOut, xmin=0, xmax=0.99, value=0.05,
               pos=[(0.92,0.25), (0.92,0.35)], c= colors[0], 
               title='Opacity\n'+ names[0].title(), title_size=txt_slider_size2)
    if len(vol_iso) >=2:
        plt.addSlider2D(sliderAlphaMeshOut2, xmin=0, xmax=0.99, value=0.05,
               pos=[(0.92,0.40), (0.92,0.50)], c=colors[1], 
               title='Opacity\n'+ names[1].title(), title_size=txt_slider_size2)
    if len(vol_iso) >=3:
        plt.addSlider2D(sliderAlphaMeshOut3, xmin=0, xmax=0.99, value=0.05,
               pos=[(0.72,0.25), (0.72,0.35)], c=colors[2],
               title='Opacity\n'+ names[2].title(), title_size=txt_slider_size2)
        
    plt.addCallback('mouse click', remove_cells)
    plt.show(sphs, vols, msg, txt, zoom=0.8)
    
    # Read file
    cells_position = organ.cellsMC['chA'].df_cells()
    deleted = cells_position['deleted']

    for cell in cells2remove: 
        deleted[cell] = 'YES'
    
    cells_position['deleted'] = deleted
    organ.cellsMC['chA'].save_cells(cells_position)
    organ.cellsMC['chA'].set_cells(cells_position)

#%% - Measuring functions
def measure_centreline(organ, nPoints):
    
    cl_measurements = organ.mH_settings['setup']['params']['2']['measure']
    cl_names = list(organ.mH_settings['measure']['CL'].keys())
    for name in cl_names: 
        df_res = df_reset_index(df=organ.mH_settings['df_res'], 
                                     mult_index= ['Parameter', 'Tissue-Contour'])
        cl_meas = {}
        ch, cont, segm = name.split('_')
        #Create linear line
        if cl_measurements['looped_length']: 
            cl_final = organ.obj_meshes[ch+'_'+cont].get_centreline(nPoints=nPoints)
            cl_length = cl_final.length()
            df_res = df_add_value(df=df_res, index=('Centreline: Looped Length', ch+'_'+cont+'_whole'), value=cl_length)
            cl_meas['looped_length'] = cl_length

        #Create linear line
        if cl_measurements['linear_length']: 
            linLine = organ.obj_meshes[ch+'_'+cont].get_linLine()
            lin_length = linLine.length()
            df_res = df_add_value(df=df_res, index=('Centreline: Linear Length', ch+'_'+cont+'_whole'), value=lin_length)
            cl_meas['lin_length'] = lin_length
            
            process = ['measure', 'CL', ch+'_'+cont+'_'+segm]
            organ.update_settings(process = process, update = cl_meas, mH='mH')
        
        df_res = df_reset_index(df=df_res, mult_index= ['Parameter', 'Tissue-Contour', 'User (Tissue-Contour)'])
        organ.mH_settings['df_res'] = df_res

        # Update organ workflow
        processes = [['MeshesProc', 'F-Measure', 'CL', ch, cont, segm],
                        ['MeshesProc', 'F-Measure','CL', ch, cont, segm],
                        ['MeshesProc', 'F-Measure','CL', ch, cont, segm]]
        for process in processes: 
            organ.update_mHworkflow(process=process, update='DONE')

    print('organ.mH_settings', organ.mH_settings)
    print('organ.workflow', organ.workflow)

def measure_submesh(organ, submesh, mesh): 

    df_res = df_reset_index(df=organ.mH_settings['df_res'], 
                                     mult_index= ['Parameter', 'Tissue-Contour'])
    
    if submesh.sub_mesh_type == 'Segment': 
        name = 'segm'
        param_name = 'Segment'
    elif submesh.sub_mesh_type == 'Section': 
        name = 'sect'
        param_name = 'Region'
    else: 
        name = 'segm-sect'
        param_name = 'Segm-Reg'

    #Get measurements to acquire
    measurements = organ.mH_settings['setup'][name]['measure']
    if measurements['Vol']: 
        vol = mesh.volume()
        df_res = df_add_value(df=df_res, index=('Volume: '+param_name, submesh.sub_name_all), value=vol)
        # organ.mH_settings['measure']['Vol('+name+')'][submesh.sub_name_all] = vol
    if measurements['SA']: 
        area = mesh.area()
        df_res = df_add_value(df=df_res, index=('Surface Area: '+param_name, submesh.sub_name_all), value=area)
        # organ.mH_settings['measure']['SA('+name+')'][submesh.sub_name_all] = area
    if name == 'segm': 
        if measurements['Ellip']: 
            #Do the ellipsoid!
            pass

     #Fill-up results table
    df_res = df_reset_index(df=df_res, mult_index= ['Parameter', 'Tissue-Contour', 'User (Tissue-Contour)'])
    organ.mH_settings['df_res'] = df_res

def measure_ellipse(organ, mesh_segm, segm_no, planar_views, divs, cut, ref_vect, gui_angles, name, plot): 

    #Get mesh and rotate it if rotation has been set
    mesh2rot = mesh_segm.clone()

    #Find the corresponding div
    for div in divs.keys(): 
        if divs[div]['segm'] == segm_no: 
            print('div:', div)
            break
    
    #Get original line
    p0 = np.array(gui_angles[cut]['original_line'][div]['p0'])
    p1 = np.array(gui_angles[cut]['original_line'][div]['p1'])
    
    ang_x, ang_y, ang_z = find_euler_angles_rotation(p0, p1, ref_vect)
    new_mesh = mesh2rot.clone().rotate_x(ang_x).rotate_y(ang_y).rotate_z(ang_z).color('magenta')
    
    if plot: 
        arrow_or = vedo.Arrow(start_pt=p1, end_pt=p0, c='darkred', s=0.1)
        new_arrow = arrow_or.clone().rotate_x(ang_x).rotate_y(ang_y).rotate_z(ang_z).color('cyan')

    x_min = []; x_max=[]; y_min=[]; y_max=[]; z_min=[]; z_max=[]
    segm_projts = {}
    for view in planar_views:
        proj_plane = planar_views[view]['proj_plane']
        color = planar_views[view]['color'][0:3]
        
        segm_proj = new_mesh.clone().project_on_plane(plane = proj_plane).c(color).alpha(0.05)
        segm_proj.legend('Proj '+proj_plane.upper())
        xa_min, xa_max, ya_min, ya_max, za_min, za_max = segm_proj.bounds()
        x_min.append(xa_min); x_max.append(xa_max)
        y_min.append(ya_min); y_max.append(ya_max)
        z_min.append(za_min); z_max.append(za_max)

        segm_projts[proj_plane] = segm_proj

    x_min = list(set([i for i in x_min if i != 0.0]))[0]
    y_min = list(set([i for i in y_min if i != 0.0]))[0]
    z_min = list(set([i for i in z_min if i != 0.0]))[0]
    x_max = list(set([i for i in x_max if i != 0.0]))[0]
    y_max = list(set([i for i in y_max if i != 0.0]))[0]
    z_max = list(set([i for i in z_max if i != 0.0]))[0]

    x_dim = x_max-x_min; 
    y_dim = y_max-y_min; 
    z_dim = z_max-z_min; 

    x_pos = x_min+x_dim/2
    y_pos = y_min+y_dim/2
    z_pos = z_min+z_dim/2

    ellip_ch = vedo.Ellipsoid(pos=(x_pos, y_pos, z_pos), 
                            axis1=(x_dim, 0,0), 
                            axis2=(0, y_dim, 0), 
                            axis3=(0, 0, z_dim), 
                            c='limegreen', alpha=0.5, res=24)
    ellip_ch.legend('Ellipsoid '+name)
    
    dims = {'X-dim': x_dim, 'Y-dim': y_dim, 'Z-dim': z_dim, 'Asphericity': ellip_ch.asphericity()}
    
    val = 80; segm_projts_f = []
    for ax in segm_projts.keys():
        if ax == 'x': 
            segm_projts[ax].x(x_min-val)
        elif ax == 'y': 
            segm_projts[ax].y(y_min-val)
        else: #proj_plane == 'z': 
            segm_projts[ax].z(z_min-val)
        segm_projts_f.append(segm_projts[ax])

    if plot: 
        m2 = tuple([ellip_ch, new_mesh]+segm_projts_f)
        # obj = [(mesh_segm, arrow_or, new_arrow, new_mesh), m2]
        # txt = [(0, organ.user_organName), (1, 'Ellipsoid created for '+name)]
        # plot_grid(obj=obj, txt=txt, axes=1, sc_side=max(organ.get_maj_bounds()))

        obj = [m2]
        txt = [(0, organ.user_organName+' > Ellipsoid created for '+name)]
        plot_grid(obj=obj, txt=txt, axes=1, sc_side=max(organ.get_maj_bounds()))

    else: 
        return dims


#%% - Points classification
def classify_heart_pts(m_whole, obj_segm, data):
    """
    Function that classifies the points that make up a mesh as atrium/ventricle, dorsal/ventral and left/right

    """
    pts_whole = m_whole.points()
    pts_classAnV = np.empty(len(pts_whole), dtype='object')

    if obj_segm != []: 
        pts_classAnV[:] = 'other'
        index_segm_all = []
        names = []
        for segm in obj_segm.keys(): 
            # Classify AnV
            print('Classifying points for '+segm)
            subm = obj_segm[segm].get_segm_mesh()
            no_segm, name = segm.split(':')
            names.append(name)
            pts_segm = subm.points()

            #Classify
            av = pts_whole.view([('', pts_whole.dtype)] * pts_whole.shape[1]).ravel()
            cv = pts_segm.view([('', pts_segm.dtype)] * pts_segm.shape[1]).ravel()
            d_isin = np.isin(av,cv)
            index_segm = np.where(d_isin == True)[0]
            print('length original index_segm: ',len(index_segm))
            if len(index_segm_all) >0: 
                for indx in index_segm_all: 
                    print(indx)
                    intersect = np.intersect1d(index_segm, indx)
                    if len(intersect)!= 0:
                        index_segm = np.setdiff1d(index_segm, intersect)
                        print('length final index_segm: ',len(index_segm))
            else: 
                pass
        
            index_segm_all.append(index_segm)
            
            pts_classAnV[index_segm] = segm
            unique = np.unique(pts_classAnV)
            print('unique:', unique)
        
        class_name = '-'.join(names)
        cols = [class_name]+list(data.keys())
        df_classPts = pd.DataFrame(columns=cols)
        df_classPts[class_name] = pts_classAnV
        print('len(df_classPts):', len(df_classPts))

        if mH_config.dev_plots: 
            plot_classif_pts(m_whole, pts_whole, list(obj_segm.keys()), pts_classAnV)
    
    else: 
        class_name = 'whole'
        cols = [class_name]+list(data.keys())
        df_classPts = pd.DataFrame(columns=cols)
        pts_classAnV[:] = 'whole'
        df_classPts[class_name] = pts_classAnV
        print('len(df_classPts):', len(df_classPts))

    for key, item in data.items():
        print(len(item))
        df_classPts[key] = item

    print("- All Done - points have been classified!\n- Sample of classified points")
    print(df_classPts.sample(10))
    alert('whistle')

    return df_classPts, class_name

#Heatmaps 3D-2D
def get_extCL_on_surf(mesh, kspl_ext, direction):
    
    # - Get unitary normal of plane to create CL_ribbon
    pl_normCLRibbon = unit_vector(direction)
    
    # Increase the resolution of the extended centreline and interpolate to unify sampling
    xd = np.diff(kspl_ext.points()[:,0])
    yd = np.diff(kspl_ext.points()[:,1])
    zd = np.diff(kspl_ext.points()[:,2])
    dist = np.sqrt(xd**2+yd**2+zd**2)
    u = np.cumsum(dist)
    u = np.hstack([[0],u])
    t = np.linspace(0, u[-1],1000)
    try: 
        resamp_pts = interpn((u,), kspl_ext.points(), t)
        print('Created resampled KSpline using interpn')
    except: 
        # Christian K Answer
        # https://stackoverflow.com/questions/19117660/how-to-generate-equispaced-interpolating-values/19122075#19122075
        xn = np.interp(t, u, kspl_ext.points()[:,0])
        yn = np.interp(t, u, kspl_ext.points()[:,1])
        zn = np.interp(t, u, kspl_ext.points()[:,2])
        resamp_pts = np.vstack((xn,yn,zn))
        resamp_pts = resamp_pts.T
        print('Created resampled KSpline using np.interp')

    kspl_ext = vedo.KSpline(resamp_pts, res = 1000).lw(5).color('deeppink').legend('kspl_extHR')
    
    #Find the points that intersect with the ribbon
    pts_int = []
    for num in range(len(kspl_ext.points())):
        cl_pt_test = kspl_ext.points()[num]
        pt_int = mesh.intersect_with_line(cl_pt_test, cl_pt_test+150*pl_normCLRibbon)
        if len(pt_int) != 0: 
            rad_pts = [np.linalg.norm(x- cl_pt_test) for x in pt_int]
            ind_pt = np.where(rad_pts == max(rad_pts))[0][0]
            pts_int.append(pt_int[ind_pt])
        else: 
            pass
    
    # KSpline on surface
    kspl_vSurf = vedo.KSpline(pts_int).color('black').lw(4).legend('kspl_VSurfaceIntMyoc')
    print('Start kspl: ', kspl_vSurf.points()[0],'- End kspl: ', kspl_vSurf.points()[-1])

    if mH_config.dev_plots: 
        vp = vedo.Plotter(N=1, axes=4)
        vp.show(kspl_vSurf, kspl_ext, mesh, at = 0, interactive = True)
    
    return kspl_vSurf

def get_extCL_highRes(organ, mesh, kspl_ext):

    trimming_set = organ.mH_settings['wf_info']['trimming']
    dict_planes = {}
    for side in ['top', 'bottom']: 
        if 'plane_info_mesh' in trimming_set[side]:
            dict_planes[side] = trimming_set[side]['plane_info_mesh']
        else: 
            dict_planes[side] = 'NI'
    print('dict_planes:', dict_planes)
    
    inv = [True, False]
    # Create kspline_extended with higher resolution, but cut it using the inflow and outflow planes defined to 
    # cut inf anf outf tracts of hearts to extract centreline
    kspl_CLnew = kspl_ext.clone()

    # Cut inflow (bottom) first
    # Cut with inflow plane if it exists
    if isinstance(dict_planes['bottom'], dict): 
        n_points_In = -10
        num_pt_inf = False
        for invert in inv:
            plane_in = vedo.Plane(pos=dict_planes['bottom']['pl_centre'], normal=dict_planes['bottom']['pl_normal'], s=(300,300))
            kspl_test = kspl_ext.clone().cutWithMesh(plane_in, invert=invert).color('gold')

            if kspl_test.NPoints() > n_points_In: 
                kspl_CLnew_cutIn = kspl_test.clone().color('darkorange')
                n_points_In = kspl_CLnew_cutIn.NPoints() 
                _, num_pt_inf = find_closest_pt_guess(kspl_CLnew_cutIn.points(), kspl_ext.points(), index_guess = -1)
                print(num_pt_inf)
                break
            else: 
                print('aaa')
                pass

        if num_pt_inf == False: 
            num_pt_inf = kspl_ext.NPoints()-1
    else: 
        num_pt_inf = kspl_ext.NPoints()-1

    # Cut outflow (top) second
    # Cut with outflow plane if it exists
    if isinstance(dict_planes['top'], dict): 
        n_points_Out = -10
        num_pt_outf = False
        for invert in inv:
            plane_out = vedo.Plane(pos=dict_planes['top']['pl_centre'], normal=dict_planes['top']['pl_normal'], s=(300,300))
            kspl_test = kspl_ext.clone().cutWithMesh(plane_out, invert=invert).color('blue')

            if kspl_test.NPoints() > n_points_Out: 
                kspl_CLnew_cutOut = kspl_test.clone().color('lime')
                n_points_Out = kspl_CLnew_cutOut.NPoints() 
                _, num_pt_outf =  find_closest_pt_guess(kspl_CLnew_cutOut.points(), kspl_ext.points(), index_guess = 0)
                print(num_pt_outf)
                break
            else: 
                print('bbb')
                pass

        if num_pt_outf == False: 
            num_pt_outf = 0
    else: 
        num_pt_outf = 0
        
    add_pts = 50
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
    print('Initial def:',num_pt_outf, num_pt_inf)
    print('Final def:',ind_outf, ind_inf)

    # Create new extended and cut kspline with higher resolution
    kspl_CLnew = vedo.KSpline(kspl_ext.points()[ind_outf:ind_inf], res = 600).lw(5).color('deeppink').legend('HighResCL')
    
    if mH_config.dev_plots: 
        vp = vedo.Plotter(N=1, axes=4)
        vp.show(kspl_CLnew, kspl_ext, mesh, at = 0, interactive = True)
    
    return kspl_CLnew

def kspl_chamber_cut(organ, mesh, kspl_CLnew, segm_cuts_info, cut, ordered_segm={}, init=False): 
    
    num_pts = {}; spheres = []; kspls_segm = {}; list_num_pts = []
    for disc in segm_cuts_info.keys(): 
        
        cut_info = segm_cuts_info[disc]
        disc_o = vedo.Cylinder(pos = cut_info['pl_centre'],r = cut_info['radius'], height = cut_info['height'], 
                            axis = cut_info['normal_unit'], c = 'purple', cap = True, res = 300)
        
        # Find position within new kspline cut by disc used to cut segments (num_pt new)
        plane_disc = vedo.Plane(pos=cut_info['pl_centre'], normal=cut_info['normal_unit'], s=(300,300))
        ksplCL_cut = kspl_CLnew.clone().cutWithMesh(plane_disc, invert=True)
        # Find point of new kspline closer to last point of kspline cut
        _, num_pt = find_closest_pt_guess(ksplCL_cut.points(), kspl_CLnew.points(),int(len(kspl_CLnew.points())*0.6))
        list_num_pts.append(num_pt)
        num_pts[disc] = num_pt
        
        # Add pt to dict
        sph_cut = vedo.Sphere(pos = kspl_CLnew.points()[num_pt], r=4, c='gold')
        sph_o = vedo.Sphere(pos = kspl_CLnew.points()[0], r=4, c='lime')
        sph_f = vedo.Sphere(pos = kspl_CLnew.points()[-1], r=4, c='cyan')
        spheres.append(sph_cut)

    if init: 
        ordered_segm = order_segms(organ, kspl_CLnew, list_num_pts, cut)
        # print('ordered_segm:', ordered_segm)
    else: 
        colors = ['tomato','darkblue','yellow','chocolate','purple','gray']
        kspl_list = []
        ordered_kspl = order_segms(organ, kspl_CLnew, list_num_pts, cut)
        print('ordered_kspl (new):', ordered_kspl)

        #Create a kspline for each segm
        nn = 0
        for div in ordered_segm.keys():
            ordered_segm[div].pop('num_pts_range', None)
            ordered_segm[div]['num_pts_range'] = ordered_kspl[div]['num_pts_range']
            numa, numb = ordered_segm[div]['num_pts_range']

            if ordered_segm[div]['invert_plane_num']:
                invert = True
            else: 
                invert = False

            n_pts = 0
            res = 595
            while n_pts < 610:
                res +=1
                kspl_pts = kspl_CLnew.points()[numa:numb]
                if invert: #atrium
                    kspl_final = vedo.KSpline(points = kspl_pts, res = res).color(colors[nn]).lw(15)
                else: #ventricle
                    kspl_final = vedo.KSpline(points = kspl_pts[::-1], res = res).color(colors[nn]).lw(15)
                n_pts = kspl_final.NPoints()
            print('n_pts:',n_pts)
            ordered_segm[div]['kspl'] = kspl_final
            kspl_list.append(kspl_final)
            nn+=1
    
        if mH_config.dev_plots: 
            vp = vedo.Plotter(N=1, axes=4)
            vp.show(kspl_CLnew, kspl_list, sph_cut, sph_o, sph_f, disc_o, mesh, at = 0, interactive = True)
        print('ordered_segm (not init):', ordered_segm)

    return ordered_segm

def order_segms(organ, kspl_CLnew, num_pts, cut): 

    names = organ.mH_settings['setup']['segm'][cut]['name_segments']
    ext_subsgm = organ.get_ext_subsgm(cut)
    ext_meshes = {}
    for sub in ext_subsgm.keys(): 
        ext_m = ext_subsgm[sub].get_segm_mesh()
        ext_meshes[sub] = ext_m

    num_pts.append(0); num_pts.append(kspl_CLnew.NPoints())
    num_pts = sorted(num_pts)

    ordered_segm = {}
    for n, numa in enumerate(num_pts[:-1]):
        numb = num_pts[n+1]
        num_btw = (numb-numa)//2+numa
        kspl_pt = kspl_CLnew.points()[num_btw]
        
        for m_ext in ext_meshes.keys(): 
            if ext_meshes[m_ext].is_inside(kspl_pt):
                ordered_segm['div'+str(n+1)] = {'num_pts_range': (numa, numb),
                                                'segm': m_ext, 
                                                'name': names[m_ext],
                                                'y_axis': (n+1, n),
                                                'kspl': None, 
                                                'invert_plane_num': None}
                break
        
        if 'div'+str(n+1) not in ordered_segm.keys(): 
            numb = num_pts[n+1]
            num_btw = numb-30
            kspl_pt = kspl_CLnew.points()[num_btw]

            for m_ext in ext_meshes.keys(): 
                if ext_meshes[m_ext].is_inside(kspl_pt):
                    ordered_segm['div'+str(n+1)] = {'num_pts_range': (numa, numb),
                                                    'segm': m_ext, 
                                                    'name': names[m_ext],
                                                    'y_axis': (n+1, n),
                                                    'kspl': None, 
                                                    'invert_plane_num': None}
                    break

        if 'div'+str(n+1) not in ordered_segm.keys(): 
            numb = num_pts[n+1]
            num_btw = numa+30
            kspl_pt = kspl_CLnew.points()[num_btw]

            for m_ext in ext_meshes.keys(): 
                if ext_meshes[m_ext].is_inside(kspl_pt):
                    ordered_segm['div'+str(n+1)] = {'num_pts_range': (numa, numb),
                                                    'segm': m_ext, 
                                                    'name': names[m_ext],
                                                    'y_axis': (n+1, n),
                                                    'kspl': None, 
                                                    'invert_plane_num': None}
                    break

    return ordered_segm

def unloop_chamber(organ, mesh, kspl_CLnew, kspl_vSurf,
                                df_classPts, labels, gui_heatmaps2d, kspl_data):

    if gui_heatmaps2d['plot']['plot_planes']: 
        plotevery = gui_heatmaps2d['plot']['every_planes']
        print('- Plotting every X number of planes:', plotevery)

    hmitem, class_name = labels
    param = df_classPts[hmitem]
    classes = df_classPts[class_name]
    
    # Create matrix with all data
    #0:x, 1:y, 2:z, 3:taken, 4:z_plane, 5:theta, 6: radius, 7: parameter
    matrix_unlooped = np.zeros((len(mesh.points()),8))
    matrix_unlooped[:,0:3] = mesh.points()
    matrix_unlooped[:,7] = param
    
    #Segment being unlooped
    if kspl_data['segm'] != 'NA': 
        chamber = kspl_data['segm']+':'+kspl_data['name']
        multip = 1
    else: 
        chamber = kspl_data['name']
        multip = 2

    print('Segment being unlooped:', chamber)

    # Get normals and centres of planes to cut heart
    no_planes = gui_heatmaps2d['nPlanes']*multip
    print('no_planes:', no_planes)
    tol2use = gui_heatmaps2d['tol']
    kspl = kspl_data['kspl']
    select_start = True
    # pl_normals, pl_centres = get_plane_normals(no_planes = no_planes, spline_pts = kspl.points())
    pl_normals, pl_centres = get_plane_normals_to_proj_kspl(organ = organ, no_planes = no_planes, 
                                                            kspl = kspl, gui_heatmaps2d = gui_heatmaps2d)

    # Give a number between 1-2 to each atrial plane, and between 0-1 to each ventricular plane
    top, bottom = kspl_data['y_axis']
    plane_num = np.linspace(top, bottom, len(pl_normals))
    # Initialise index for HR centreline to identify chamber in which planes are cutting and flip plane values for initial segment
    if not kspl_data['invert_plane_num']:
        plane_num = plane_num[::-1]
        t2b = True
        print('aja ventricle 0')
        invert = True
    else:      
        print('aja atrium 0')
        t2b = False
        invert = True

    # Iterate through planes
    not_init = True; started = False
    sph_pt_surfinal = vedo.Sphere(kspl_vSurf.points()[0], r=1, c='violet')

    #Create empty lists to add objects every 20 planes
    colors = ['orangered','gold', 'olive','lime', 'teal', 'aqua', 'dodgerblue', 'navy', 'indigo', 'purple', 'hotpink','chocolate']*10
    list_planes = []; list_CL_sph = []; list_vSurf_sph = []; list_prev_sph = []; 
    arr_zero_deg = []; sph_left = []; sph_right = []; ii = 0
    for i, normal, centre in zip(count(), pl_normals, pl_centres):
        print('\n>-Plane num:', i, 'normal:', normal)
        # Initialise variables for each chamber
        if t2b: 
            normal = -normal

        # A. Get cut plane info
        #   Info Plane (Plane orientation/normal representation)
        arr_vectPlCut = vedo.Arrow(centre, centre+normal*10, s = 0.1, c='orange')
        
        # B. Cut high resolution centreline with plane and define chamber
        plane_cut = vedo.Plane(pos=centre, normal=normal, s=(300,300)).color('light slate gray').alpha(0.4)
        
        # If the ventricle is the one being unlooped, then initialise index as length of the cut centreline
        if i == 0:
            ksplCL_cut_all = kspl_CLnew.clone().cutWithMesh(plane_cut, invert=invert).lw(5).color('tomato')
            ksplCL_cut_split = ksplCL_cut_all.split()
            print('ksplCL_cut_split:', ksplCL_cut_split)
            if len(ksplCL_cut_split) == 1: 
                ksplCL_cut = ksplCL_cut_split[0].lw(8).color('darkblue')
            else: 
                try:
                    for kspl in ksplCL_cut_split:
                        #Find the CL piece that has the last point 
                        #Revisar para la auricula si el invert voltearlo? 
                        if any((kspl.points()[:]==pt_out).all(axis=1)):
                            print('kspl found!')
                            ksplCL_cut = kspl.lw(8).color('gold')
                            break
                except:
                    ksplCL_cut = ksplCL_cut_split[0].lw(8).color('darkblue')

            index_guess= len(ksplCL_cut.points())

        # Find closest point of the high resolution centreline that has been cut by plane and is close to the
        #  previous pt_o (index_guess point)
        pt_out, index_guess = find_closest_pt2pl(normal, centre, kspl_CLnew, index_guess)
        print('pt_out:', pt_out, '- index_guess:', index_guess)
        sph_pt_out = vedo.Sphere(pt_out, r=2, c='turquoise')

        # C. Cut surface centreline (kspl_vSurf) with plane and identify 0 deg angle point
        #   Find point of surf_centreline cut by plane (ksplCL_vSurf_cut)
        kspl_vSurf_cut_all = kspl_vSurf.clone().cutWithMesh(plane_cut, invert=invert).lw(5).color('magenta')
        print(type(kspl_vSurf_cut_all))
        try: 
            kspl_vSurf_cut_split = kspl_vSurf_cut_all.split()
            print('kspl_vSurf_cut_split:', kspl_vSurf_cut_split, ' - len:',len(kspl_vSurf_cut_split))
        except: 
            print('except trying to cut kspl_vSurf')
            kspl_vSurf_cut_split = []

        if len(kspl_vSurf_cut_split) > 0: 
            if not_init: 
                print('not init input')
                if len(kspl_vSurf_cut_split) == 1: 
                    kspl_vSurf_cut = kspl_vSurf_cut_split[0].lw(8).color('white')
                else: 
                    for kspls in kspl_vSurf_cut_split: 
                        if any((kspls.points()[:]==pt_surf).all(axis=1)): 
                            print('ksplsurf found!')
                            kspl_vSurf_cut = kspls.lw(8).color('cyan')
                            break
                # last_pt_cut = kspl_vSurf_cut.points()[-1]
                # print('last_pt_cut:', last_pt_cut)
                # sph_pt_surf2 = vedo.Sphere(last_pt_cut, r=2, c='black')
        
                idx_surf = len(kspl_vSurf_cut.points())
                not_init = False

            pt_surf, idx_surf = find_closest_pt2pl(normal, centre, kspl_vSurf, idx_surf)
            print('pt_surf:', pt_surf, '- idx_surf:', idx_surf)
            sph_pt_surf = vedo.Sphere(pt_surf, r=2, c='gold')

            # if findDist(pt_surf, last_pt_cut) > 1: 
            #     print('Dist: ', findDist(pt_surf,last_pt_cut))
            #     alert('bubble')
            #     pts_in_plane = [pt_surf, last_pt_cut]
            #     idx_in_plane = [idx_surf, len(kspl_vSurf_cut.points())]

            #     objects = [mesh.alpha(0.05), plane_cut, kspl_CLnew, ksplCL_cut, sph_pt_out, kspl_vSurf, kspl_vSurf_cut, sph_pt_surf, sph_pt_surf2, sph_pt_surfinal]
            #     pt_selected, index_selected = select_sph_vSurf(pts_in_plane, idx_in_plane, objects)
                
            #     if index_selected == 'empty': 
            #         print('Skipping this one plane:', i)
            #     else: 
            #         pt_surf = pt_selected
            #         idx_surf = index_selected
            #         sph_pt_surf = vedo.Sphere(pt_surf, r=3, c='darkviolet')

            #         #Find the real idx_surf on the whole kspl_vSurf
            #         idx_surf = find_pt_idx_in_kspline(kspl_vSurf, pt_surf)
            # else: 
            # vp= vedo.Plotter(N=1, axes=13)
            # vp.show(mesh.alpha(0.05), ksplCL_cut, kspl_vSurf_cut, kspl_CLnew, arr_vectPlCut, plane_cut, kspl_vSurf, sph_pt_out, sph_pt_surf, sph_pt_surf2, sph_pt_surfinal, at=0, interactive=True)
            # # pass

            if select_start and not started: 
                pts_in_plane = [pt_surf]; idx_in_plane = [idx_surf]
                objects = [mesh.alpha(0.05), plane_cut, kspl_CLnew, sph_pt_out, kspl_vSurf, sph_pt_surf, sph_pt_surfinal]
                _, index_selected = select_sph_vSurf(pts_in_plane, idx_in_plane, objects)
                if index_selected != 'empty': 
                    print('Skipping this one plane:', i)
                    started = True

            if i % 20 == 0 and i != 0 and (mH_config.dev_hm3d2d or gui_heatmaps2d['plot']['plot_planes']): 
                list_planes.append(plane_cut.alpha(0.2).color(colors[ii])); list_CL_sph.append(sph_pt_out.color(colors[ii]))
                list_vSurf_sph.append(sph_pt_surf.color(colors[ii])); list_prev_sph.append(sph_pt_surfinal)
                ii+=1

            print('Final plane:', i, '- pt_surf:', pt_surf, '- idx_surf:', idx_surf)
            sph_pt_surfinal = vedo.Sphere(pt_surf, r=1, c='violet')

            # if mH_config.dev_hm3d2d: 
            #     vp= vedo.Plotter(N=1, axes=13)
            #     vp.show(mesh.alpha(0.05), plane_cut, kspl_CLnew, kspl_vSurf, arr_vectPlCut, sph_pt_out, sph_pt_surf, at=0, interactive=True)
    
        else: 
            idx_surf = None
            print('No-cut plane: ', i)

        #D. Having all the points, not cut mesh 
        if idx_surf != None and isinstance(idx_surf, int) and started:
            # Vector from centre to cl_surface point being cut by plane
            v_zero = unit_vector(kspl_vSurf.points()[idx_surf] - centre)
            if i % 20 == 0 and i != 0 and mH_config.dev_hm3d2d: 
                arr_zero_deg.append(vedo.Arrow(centre, kspl_vSurf.points()[idx_surf], s = 0.1, c='dodgerblue'))
            # D. Get points of mesh at plane
            d_points = np.absolute(np.dot(np.subtract(matrix_unlooped[:,0:3],np.asarray(centre)),np.asarray(normal)))
            # Find the indexes of the points that have not been yet taken, are at the plane and are in the 
            # chamber being analysed
            index_ptsAtPlane = np.where((d_points <= tol2use) & (matrix_unlooped[:,3] == 0) & (classes == chamber))
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
            if mH_config.dev_plots or gui_heatmaps2d['plot']['plot_planes']:
                if i % plotevery == 0 and i != 0:
                    sphL = []; sphR = []
                    for num, pt in enumerate(ptsC):
                        if num % 20 == 0:
                            if lORr[num] == 1:
                                sphL.append(vedo.Sphere(pt+centre, r=2, c='blueviolet'))
                            else:
                                sphR.append(vedo.Sphere(pt+centre, r=2, c='gold'))

                    text = '>> Unlooping the heart (segment: '+chamber+') - Plane No: '+str(i)+'/'+str(no_planes+2)
                    txt = vedo.Text2D(text, c=txt_color, font=txt_font, s=txt_size)
                    sph_centre = vedo.Sphere(centre, r=2, c='red')
                    arr_centre2vzero = vedo.Arrow(centre, kspl_vSurf.points()[idx_surf], s = 0.1, c='light green')
                    vp= vedo.Plotter(N=1, axes=13)
                    vp.add_icon(logo, pos=(0.1,1), size=0.25)
                    vp.show(mesh, sphL, sphR, plane_cut, arr_vectPlCut, kspl_CLnew, sph_pt_out, kspl_vSurf, sph_pt_surf, sph_centre, arr_centre2vzero, txt, at=0, interactive=True)
                
            if i % 20 == 0 and i != 0 and (mH_config.dev_hm3d2d or gui_heatmaps2d['plot']['plot_planes']):
                for num, pt in enumerate(ptsC):
                    if num % 20 == 0:
                        if lORr[num] == 1:
                            sph_left.append(vedo.Sphere(pt+centre, r=2, c=colors[ii-1]))
                        else:
                            sph_right.append(vedo.Sphere(pt+centre, r=2, c=colors[ii-1]).alpha(0.5))
        else: 
            pass
            # print('-Plane No.', i, ' - PASS!')

    if mH_config.dev_plots or gui_heatmaps2d['plot']['plot_planes']:
        text = '>> Unlooping the heart (segment: '+chamber+')'
        txt = vedo.Text2D(text, c=txt_color, font=txt_font, s=txt_size)
        vp= vedo.Plotter(N=1, axes=13)
        vp.add_icon(logo, pos=(0.1,1), size=0.25)
        vp.show(mesh.alpha(0.05),  kspl_CLnew, kspl_vSurf, list_planes, list_CL_sph, list_vSurf_sph, list_prev_sph, sph_left, sph_right, at=0, interactive=True)

    df_unlooped = pd.DataFrame(matrix_unlooped, columns=['x','y','z','taken','z_plane','theta','radius',hmitem])
    df_unlooped = df_unlooped[df_unlooped['taken']==1]
    df_unloopedf = df_unlooped.drop(['x', 'y','z'], axis=1)
    df_unloopedf.astype({'taken': 'bool','z_plane':'float16','theta':'float16','radius':'float16', hmitem:'float16'}).dtypes
    print(df_unloopedf.sample(10))

    organ_name = organ.user_organName
    title_df = organ_name+'_dfUnloop_'+hmitem+'_'+kspl_data['name']+'.csv'
    dir_df = organ.dir_res(dir='csv_all') / title_df
    df_unloopedf.to_csv(dir_df)

    return df_unloopedf, title_df
      
def heatmap_unlooped(organ, kspl_data, df_unloopedf, hmitem, ch, gui_thball):
    """
    Function to create heatmap of unlooped data, specifically of the columns given as input by the val2unloop variable.  
    """
    print('hmitem/val:', hmitem)
    print(df_unloopedf.sample(10))

    #Get all saving settings
    organ_name = organ.user_organName
    title_df = organ_name+'_hmUnloop_'+hmitem+'_'+kspl_data['name']+'.csv'
    extension = gui_thball['heatmaps2D']['extension']
    title_hm = organ_name+'_hmUnloop_'+hmitem+'_'+kspl_data['name']+extension
    dir_df = organ.dir_res(dir='csv_all') / title_df
    dir_hm = organ.dir_res(dir='imgs_videos') / title_hm
    print(dir_df, '\n', dir_hm)

    print('\n- Creating heatmaps for '+hmitem+'_'+kspl_data['name'].title())
    df_unloopedf = df_unloopedf.drop(['taken'], axis=1)
    df_unloopedf.astype('float16').dtypes
    
    heatmap = pd.pivot_table(df_unloopedf, values= hmitem, columns = 'theta', index='z_plane', aggfunc=np.max)
    heatmap.astype('float16').dtypes
    print(heatmap.sample(10))

    #Create the figure
    tissue_name = organ.mH_settings['setup']['name_chs'][ch]
    if 'th' in hmitem: 
        if 'i2e' in hmitem: 
            title = organ_name +' - '+tissue_name.title()+' Thickness (int2ext) [um] - '+kspl_data['name'].title()
        else: 
            title = organ_name +' - '+tissue_name.title()+' Thickness (ext2int) [um] - '+kspl_data['name'].title()
    else: 
        #Get cl info
        _, infoo = hmitem.split('[')
        ch_cont, cl_info = infoo.split('(')
        _, cont = ch_cont.split('-')
        from_cl, from_cl_type = cl_info[3:-2].split('-')
        title = organ_name +' - CL > Tissue ('+tissue_name.title()+'_'+cont+') - CL: '+from_cl+'_'+from_cl_type+' [um] - '+kspl_data['name'].title()
    print('- title:', title)

    #Get all construction settings
    cmap = gui_thball[hmitem]['colormap']
    vmin = gui_thball[hmitem]['min_val']
    vmax = gui_thball[hmitem]['max_val']

    # Make figure
    fig, ax = plt.subplots(figsize=(16, 10))
    c = ax.pcolor(heatmap, cmap=cmap, vmin = vmin, vmax = vmax)
    cb = fig.colorbar(c, ax=ax)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=10)
    ax.invert_yaxis()
    # b = sns.heatmap(heatmap, cmap=cmap, vmin = vmin, vmax = vmax, ax=ax)

    # set the xticks
    x_pos = ax.get_xticks()
    x_pos_new = np.linspace(x_pos[0], x_pos[-1], 19)
    x_lab_new = np.arange(-180,200,20)
    ax.set_xticks(x_pos_new) 
    ax.set_xticklabels(x_lab_new, rotation=30, fontsize=10)#, fontname='Arial')
    
    xlabels=np.linspace(heatmap.columns.min(), heatmap.columns.max(), len(x_pos)).round(3)
    print('xlabels:', xlabels)

    # set the yticks
    y_pos = ax.get_yticks()
    # y_pos_new = np.linspace(y_pos[0], y_pos[-1], 11)
    # y_labels = sorted(list(kspl_data['y_axis']))
    # y_lab_new = np.linspace(y_labels[0],y_labels[1],11)
    # y_lab_new = [format(y,'.2f') for y in y_lab_new]
    # ax.set_yticks(y_pos_new) 
    # ax.set_yticklabels(y_lab_new, rotation=0, fontsize=10)#, fontname='Arial')

    ylabels=np.linspace(heatmap.index.min(), heatmap.index.max(), len(y_pos)).round(2)
    ax.set_yticks(ticks=y_pos, labels=ylabels)
    ax.set_yticklabels(ylabels, rotation=0, fontsize=10)#, fontname='Arial')
    print('ylabels:', ylabels)
    
    y_text = 'Centreline Position ['+kspl_data['name'].title()+']'
    plt.ylabel(y_text, fontsize=10)
    plt.xlabel('Angle (\N{DEGREE SIGN})', fontsize=10)
    plt.title(title, fontsize = 12)

    for pos in ['top', 'right', 'bottom', 'left']:
        ax.spines[pos].set_visible(False)

    #Save figure and heatmap dataframe
    plt.savefig(dir_hm, dpi=300, bbox_inches='tight', transparent=True)
    alert('clown')

def get_unlooped_heatmap(hmitem, dir_df): 

    df_unloopedf= pd.read_csv(str(dir_df))
    # print(df_unloopedf.sample(10))
    df_unloopedf = df_unloopedf.drop(['taken'], axis=1)
    # print(df_unloopedf.sample(10))
    df_unloopedf.astype('float16').dtypes

    heatmap = pd.pivot_table(df_unloopedf, values= hmitem, columns = 'theta', index='z_plane', aggfunc='max')
    # heatmap.astype('float16').dtypes

    return heatmap

def select_sph_vSurf(pts_in_plane, idx_in_plane, objects): 
    
    colors = ['gold', 'black']#['orangered','gold', 'olive','lime', 'teal', 'aqua', 'dodgerblue', 'navy', 'indigo', 'purple', 'hotpink','chocolate']*10
    # print('pts_in_plane:', pts_in_plane)
    sphs_in_plane = []
    for pp, pt in enumerate(pts_in_plane):
        sph_pt = vedo.Sphere(pos = pt, c=colors[pp], r=3)
        sph_pt.name = f"No.{pp}"
        sphs_in_plane.append(sph_pt)
    
    ind_sph_sel = []
    def select_sph_in_plane(evt):
        if not evt.actor: return
        if isinstance(evt.actor, vedo.Sphere): 
            sil = evt.actor.silhouette().lineWidth(6).c('lawn green')
            print("You clicked: Sphere "+evt.actor.name)
            plt.remove(silcont.pop()).add(sil)
            silcont.append(sil)
            ind_sph_sel.append(evt.actor.name.split('.')[-1])
            msg.text("You clicked: Sphere "+evt.actor.name)
    silcont = [None]
    
    msg = vedo.Text2D("", pos="bottom-center", c=txt_color, font=txt_font, s=txt_size, bg='red', alpha=0.2)
    txt_sel = vedo.Text2D('> Instructions: Click on the yellow sphere on the plane from which you would like to\nstart unlooping the 3D heatmap. To advance in plane number close this window without \nselecting the sphere and a new plot with the next plane will appear.', c=txt_color, font=txt_font, s=txt_size)
    plt = vedo.Plotter(axes=1)
    plt.addCallback('mouse click', select_sph_in_plane)
    plt.show(objects, sphs_in_plane, txt_sel, msg, zoom=1.2)
    
    print('ind_sph_sel:', ind_sph_sel)
    if len(ind_sph_sel) > 0:
        index = int(ind_sph_sel[-1])
        pt_out = pts_in_plane[int(index)]
        idx_out = idx_in_plane[int(index)]
    else: 
        idx_out = 'empty'
        pt_out = None
    
    return pt_out, idx_out

#%% - Vectorial calculations
def find_angle_btw_pts(pts1, pts2):
    """
    Function that returns the angle between two vectors on the input-plane
    """
    
    #Get the magnitude of the vector
    mag_v1 = findDist(pts1[0],pts1[1])
    mag_v2 = findDist(pts2[0],pts2[1])

    vect1 = pts1[1]-pts1[0]
    vect2 = pts2[1]-pts2[0]

    dotProd = np.dot(vect1,vect2)
    angle = math.degrees(math.acos(dotProd/(mag_v1*mag_v2)))

    return angle

def find_angle_btw_pts_atan(vref, vect, normal):
    #https://stackoverflow.com/a/33920320
    #atan2((Va x Vb) . Vn, Va . Vb)
    va = vref
    vb = vect

    vaxvb = np.cross(va, vb)
    dotvn = np.dot(vaxvb, normal)
    vadotvb = np.dot(va, vb)
    angle =  math.degrees(math.atan2(dotvn, vadotvb))
    
    return angle

def findDist(pt1, pt2):
    """
    Function that returns the distance between two points given as input
    """
    squared_dist = np.sum((pt1-pt2)**2, axis=0)
    dist = np.sqrt(squared_dist)

    return dist

def find_pt_idx_in_kspline(kspl, pt_in):

    value = 99999999
    for nn, pt in enumerate(kspl.points()):
        dist = findDist(pt, pt_in)
        if dist < value: 
            value = dist
            idx = nn

    print('Ffinal distance: ',findDist(pt_in, kspl.points()[idx]))

    return idx
     
def find_euler_angles_rotation(p0, p1, ref_vect): 

    vect_zero = p0-p1
    vect_unit = unit_vector(vect_zero)
    mat = rotation_matrix_from_vectors(vect_unit, ref_vect)
    vec1_rot = mat.dot(vect_unit)
    print('Assert rotation:', np.allclose(vec1_rot / np.linalg.norm(vec1_rot), ref_vect / np.linalg.norm(ref_vect)))
    r = R.from_matrix(mat)
    angles = r.as_euler("zyx",degrees=True)
    ang_z = angles[0]
    ang_y = angles[1]
    ang_x = angles[2]

    return ang_x, ang_y, ang_z

def rotation_matrix_from_vectors(vec1, vec2):
    #https://stackoverflow.com/a/59204638
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

    return rotation_matrix

def get_plane_normal2pt (pt_num, points):
    """
    Funtion that gets a plane normal to a point in a spline

    """

    pt_centre = points[pt_num]
    normal = points[pt_num-1]-points[pt_num+1]

    return normal, pt_centre

def find_closest_pt_guess(pts_cut, pts, index_guess):
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
            
    # Now that we know the coordinates (pt_o) of the cut centreline that is cut by the plane and closest to the pt_guess then
    # find the closest point of all the centreline points closest to pt_o and get also its index
    min_dist = 999999999999999999
    for n, pt in enumerate(pts):
        squared_dist = np.sum((pt-pt_o)**2, axis=0)
        dist = np.sqrt(squared_dist)

        if dist < min_dist:
            min_dist = dist
            pt_out = pt
            num_pt = n

    return pt_out, num_pt

def find_closest_pt2pl(normal, centre, kspl_cut, index_guess): 
    tol = len(kspl_cut.points())//25
    # print('tol:', tol, ' - index_guess:',index_guess, type(tol), type(index_guess))
    values = []
    for point in kspl_cut.points():
        diff_points = point - centre
        values.append(abs(np.dot(diff_points,normal)))
    
    # Set the region of the kspl in which to look for the point
    if index_guess-tol < 0: 
        cuta = 0
    else: 
        cuta = index_guess-tol
    if index_guess+tol > len(kspl_cut.points())-1: 
        cutb = len(kspl_cut.points())-1
    else: 
        cutb = index_guess+tol

    # print(cuta, index_guess, cutb)
    values_cut = values[cuta:cutb]
    # print('min(values):', min(values_cut))
    num_pt = values.index(min(values_cut))
    pt_out = kspl_cut.points()[num_pt]

    return pt_out, num_pt

def get_plane_normals(no_planes, spline_pts):
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
    list_index = list(range(len(spline_pts)-2,1,-every))
    for i in list_index:
        pt_centre.append(spline_pts[i])
        normal = spline_pts[i-1]-spline_pts[i]
        normals.append(normal)

    return normals, pt_centre

def get_plane_normals_other_centres(no_planes, proj_pts, kspl_pts):
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
    every = (len(proj_pts)+2)//no_planes
    print(every)
    list_index = list(range(len(proj_pts)-2,1,-every))
    for i in list_index:
        pt_centre.append(kspl_pts[i])
        normal = proj_pts[i-1]-proj_pts[i]
        normals.append(normal)

    return normals, pt_centre

def get_plane_normals_to_proj_kspl(organ, no_planes, kspl, gui_heatmaps2d):

    #Get cube's face centre
    cube_name = gui_heatmaps2d['axis_lab']
    cube2use = getattr(organ, cube_name.lower()+'_cube')
    orient_cube = cube2use['cube']
    cube_face = int(gui_heatmaps2d['direction']['plane_no'])
    ext_plane = gui_heatmaps2d['direction']['plane_normal']

    cell_centre = orient_cube.cell_centers()[int(cube_face)]

    plane_ext = vedo.Plane(pos = cell_centre, normal = ext_plane)
    proj_kspl = kspl.clone().project_on_plane(plane=plane_ext).lw(3).color('black')

    # print('len(kspl.points()):', len(kspl.points()))
    # print('len(proj_kspl.points()):', len(proj_kspl.points()))

    pl_normals, pl_centres = get_plane_normals_other_centres(no_planes = no_planes, 
                                                             proj_pts = proj_kspl.points(), 
                                                             kspl_pts = kspl.points())

    return pl_normals, pl_centres
   
#%% - Plotting functions
def plot_grid(obj:list, txt=[], axes=1, zoom=1, lg_pos='top-left',
              sc_side=350, azimuth = 0, elevation = 0, add_scale_cube=True):
    
    # Create ScaleCube
    if add_scale_cube: 
        if isinstance(obj[0], tuple):
            scale_cube = vedo.Cube(pos=obj[0][0].center_of_mass(), side=sc_side, c='white', alpha=0.01).legend('ScaleCube')
        else: 
            scale_cube = vedo.Cube(pos=obj[0].center_of_mass(), side=sc_side, c='white', alpha=0.01).legend('ScaleCube')
    else: 
        scale_cube = []
    
    # Set logo position
    if lg_pos =='top-left':
        if len(obj)>5:
            pos = (0.1,3)
        elif len(obj)>3:
            pos = (0.1,2)
        else:
            pos = (0.1,1)
    else: 
        if len(obj)>3:
            pos = (0.1,2) # to correct 
        else:
            pos = (0.8,0.05)
    
    n_obj = len(obj)
    # Create tuples for text
    post = [tup[0] for tup in txt]; txt_out = []; n = 0
    for num in range(n_obj):
        if num in post:
            txt_out.append(vedo.Text2D(txt[n][1], c=txt_color, font=txt_font, s=txt_size))
            n += 1
        else: 
            txt_out.append(vedo.Text2D('', c=txt_color, font=txt_font, s=txt_size))
    
    # Add scale_cube to last plot
    if isinstance(obj[n_obj-1], tuple):
        obj_list = list(obj[n_obj-1])
        obj_list.append(scale_cube)
        obj_f = tuple(obj_list)
    else: #vedo.mesh.Mesh
        obj_f = (obj[n_obj-1], scale_cube)
    obj[n_obj-1] = obj_f
    
    # Now plot
    lbox = []
    vp = vedo.Plotter(N=len(obj), axes=axes)
    vp.add_icon(logo, pos=pos, size=0.25)
    for num in range(len(obj)):
        if isinstance(obj[num], tuple):
            try: 
                lbox.append(vedo.LegendBox(list(obj[num]), font=leg_font, width=leg_width))
            except: 
                lbox.append('')
        else:
            lbox.append(vedo.LegendBox([obj[num]], font=leg_font, width=leg_width))
        if num != len(obj)-1:
            vp.show(obj[num], lbox[num], txt_out[num], at=num)
        else: # num == len(obj)-1
            vp.show(obj[num], lbox[num], txt_out[num], at=num, zoom=zoom, azimuth = azimuth, elevation = elevation, interactive=True)

def plot_all_organ(organ):
    
    obj = []
    txt = [(0, organ.user_organName  + ' - All organ meshes')]
    for mesh in organ.meshes.keys():
        # print(organ.obj_meshes[mesh].mesh)
        obj.append((organ.obj_meshes[mesh].mesh))
    
    plot_grid(obj=obj, txt=txt, axes=5, sc_side=max(organ.get_maj_bounds()))

def plot_meas_meshes(organ, meas:list, color_map = 'turbo', alpha=1):
    obj = []
    for param in meas: 
        names = [item for item in organ.parent_project.mH_param2meas if param in item]
        for name in names: 
            ch = name[0]; cont = name[1]
            mesh_mHA = organ.obj_meshes[ch+'_'+cont]
            if 'ball' in param:
                from_cl = organ.mH_settings['wf_info']['MeshesProc']['D-Ballooning'][ch][cont]['from_cl']
                from_cl_type = organ.mH_settings['wf_info']['MeshesProc']['D-Ballooning'][ch][cont]['from_cl_type']
                mesh_out = mesh_mHA.balloon_mesh(n_type=from_cl+'_'+from_cl_type, color_map=color_map, alpha=alpha)
            elif 'thickness' in param:
                n_type = name[3][-7:].replace('>','TO')
                mesh_out = mesh_mHA.thickness_mesh(n_type=n_type, color_map=color_map, alpha=alpha)
            obj.append((mesh_out))
    
    if len(obj) == 0:
        print('>> No meshes were recognised for param', meas)
        alert('error_beep')
    else: 
        txt = [(0, organ.user_organName)]
        plot_grid(obj, txt, axes=5, sc_side=max(organ.get_maj_bounds()))
        
def plot_organCLs(organ, axes=5, plotshow=True):

    organ_info = organ.mH_settings['setup']
    dict_cl = {}; obj = []; txt = []; n = 0
    for item in organ.mH_settings['measure']['CL']: 
        n += 1
        ch, cont, segm = item.split('_')
        name = organ_info['name_chs'][ch]+'-'+cont+' ('+ch+'-'+cont+'-'+segm+')'
        dict_cl[n] = name
        mesh_o = organ.obj_meshes[ch+'_'+cont]
        cl_o = 'cl'# mesh_o.get_centreline(color = 'indigo')
        obj.append((mesh_o.mesh))#, cl_o))
        if n-1==0:
            txt.append((n-1, organ.user_organName+'\n-> '+name))
        else:  
            txt.append((n-1, '\n-> '+name))
    if plotshow: 
        plot_grid(obj=obj, txt=txt, axes=axes, sc_side=max(organ.get_maj_bounds()))

    return dict_cl

def plot_classif_pts(mesh, pts_whole, names_class, pts_class, every=50):
    """
    Function that plots a subset of the classified points
    """
    color = ['tomato','gold','magenta','deepblue','pink','crimson']
    
    text = '>> Point classification'
    txt = vedo.Text2D(text, c=txt_color, font=txt_font, s=txt_size)

    vp = vedo.Plotter(shape = (1, 2), axes = 13)#13
    vp.show(mesh, txt, at = 0)
    spheres = []; text_class = []
    for i, sp_class in enumerate(names_class):
        ind_first = np.where(pts_class == sp_class)[0]
        print(len(ind_first))
        pts_first = pts_whole[ind_first]
        sph_first = vedo.Spheres(pts_first[0::every,:], c = color[i], r=1).legend(names_class[i])
        spheres.append(sph_first)
        text_class.append(names_class[i]+':'+color[i])

    text_classf = ', '.join(text_class)
    txt2 = vedo.Text2D(text_classf,  c=txt_color, font=txt_font, s=txt_size)
    vp.show(mesh, spheres, txt2, at = 1, interactive = True)

#%% - Plane handling functions 
def get_plane(filename, txt:str, meshes:list, settings: dict, def_pl = None, 
                             option = [True,True,True,True,True,True]):#
    '''
    Function that creates a plane defined by the user
    '''
    
    if all([isinstance(mesh, vedo.Mesh) for mesh in meshes]):
        meshes_mesh = meshes
    else: 
        meshes_mesh = [mesh.mesh for mesh in meshes]
        
    # Create plane
    if def_pl != None:
        plane, normal, rotX, rotY, rotZ = get_plane_pos(filename, txt, meshes_mesh, settings, option, def_pl)
    else:
        plane, normal, rotX, rotY, rotZ = get_plane_pos(filename, txt, meshes_mesh, settings, option)
        
    # Get new normal of rotated plane
    normal_corrected = new_normal_3DRot(normal, rotX, rotY, rotZ)
    # Get central point of new plane and create sphere
    pl_centre = plane.pos()
    sph_centre = vedo.Sphere(pos=pl_centre,r=2,c='black')
    # Build new plane to confirm
    plane_new = vedo.Plane(pos=pl_centre,normal=normal_corrected).color('green').alpha(1).legend('New Plane')

    normal_txt = str([' {:.2f}'.format(i) for i in normal_corrected]).replace("'","")
    centre_txt = str([' {:.2f}'.format(i) for i in pl_centre]).replace("'","")
    text = filename+'\n\nUser defined plane to '+ txt.upper() +'.\nPlane normal: '+normal_txt+' - Plane centre: '+centre_txt+'.\nClose the window when done.'
    txt2D = vedo.Text2D(text, c=txt_color, font=txt_font, s=txt_size)

    vp = vedo.Plotter(N=1, axes=4)
    vp.add_icon(logo, pos=(0.8,0.05), size=0.25)
    vp.show(meshes_mesh, plane, plane_new, sph_centre, txt2D, at=0, viewup='y', azimuth=0, elevation=0, interactive=True)

    pl_dict = {'pl_normal': normal_corrected,
                    'pl_centre': pl_centre}
    
    return plane_new, pl_dict

def get_plane_pos(filename, txt, meshes, settings, option, 
                     def_pl= {'pl_normal': (0,1,0), 'pl_centre': []}, zoom = 0.5):
    '''
    Function that shows a plot so that the user can define a plane (mesh opacity can be changed)
    meshes: list (outer mesh in position 0, inner mesh in position 1 or 2)

    '''
    
    xmin, xmax, ymin, ymax, zmin, zmax = meshes[0].bounds()
    x_size = xmax - xmin; y_size = ymax - ymin; z_size = zmax - zmin

    xval = sorted([xmin-0.3*x_size,xmax+0.3*x_size])
    yval = sorted([ymin-0.3*y_size,ymax+0.3*y_size])
    zval = sorted([zmin-0.3*z_size,zmax+0.3*z_size])

    box_size = max(x_size, y_size, z_size)

    if len(def_pl['pl_centre']) != 3:
        centre = (x_size/2+xmin, y_size/2+ymin, z_size/2+zmin)
    else: 
        centre = def_pl['pl_centre']
    normal = def_pl['pl_normal']

    rotX = [0]; rotY = [0]; rotZ = [0]

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
        meshes[0].alpha(valueAlpha)
    
    if len(meshes) >= 2: 
        def sliderAlphaMeshOut2(widget, event):
            valueAlpha = widget.GetRepresentation().GetValue()
            meshes[1].alpha(valueAlpha)
    if len(meshes) >= 3: 
        def sliderAlphaMeshOut3(widget, event):
            valueAlpha = widget.GetRepresentation().GetValue()
            meshes[2].alpha(valueAlpha)
    if len(meshes) >= 4: 
        def sliderAlphaMeshOut4(widget, event):
            valueAlpha = widget.GetRepresentation().GetValue()
            meshes[3].alpha(valueAlpha)
    if len(meshes) >= 5: 
        def sliderAlphaMeshOut5(widget, event):
            valueAlpha = widget.GetRepresentation().GetValue()
            meshes[4].alpha(valueAlpha)

    lbox = vedo.LegendBox(meshes, font=leg_font, width=leg_width, padding=1)
    #vedo.settings.legendSize = .2
    vp = vedo.Plotter(N=1, axes=8)
    vp.add_icon(logo, pos=(0.85,0.75), size=0.10)
    plane = vedo.Plane(pos=centre, normal=normal, 
                       s=(box_size*1.5, box_size*1.5)).color('dimgray').alpha(1)
    if option[0]: #sliderX
        vp.addSlider2D(sliderX, xval[0], xval[1], value=centre[0],
                    pos=[(0.1,0.15), (0.3,0.15)], title='- > x position > +', 
                    c='crimson', title_size=txt_slider_size)
    if option[1]: #sliderY
        vp.addSlider2D(sliderY, yval[0], yval[1], value=centre[1],
                    pos=[(0.4,0.15), (0.6,0.15)], title='- > y position > +', 
                    c='dodgerblue', title_size=txt_slider_size)
    if option[2]: #sliderZ
        vp.addSlider2D(sliderZ, zval[0], zval[1], value=centre[2],
                    pos=[(0.7,0.15), (0.9,0.15)], title='- > z position > +', 
                    c='limegreen', title_size=txt_slider_size)
    if option[3]: #sliderRotX
        vp.addSlider2D(sliderRotX, -1, +1, value=0,
                    pos=[(0.1,0.05), (0.3,0.05)], title='- > x rotation > +', 
                    c='deeppink', title_size=txt_slider_size)
    if option[4]: #sliderRotY
        vp.addSlider2D(sliderRotY, -1, +1, value=0,
                    pos=[(0.4,0.05), (0.6,0.05)], title='- > y rotation > +', 
                    c='gold', title_size=txt_slider_size)
    if option[5]: #sliderRotZ
        vp.addSlider2D(sliderRotZ, -1, +1, value=0,
                    pos=[(0.7,0.05), (0.9,0.05)], title='- > z rotation > +', 
                    c='teal', title_size=txt_slider_size)
        
    txt_slider_size2 = 0.7
    vp.addSlider2D(sliderAlphaMeshOut, xmin=0, xmax=0.99, value=0.01,
               pos=[(0.92,0.25), (0.92,0.35)], c= settings['color'][0], 
               title='Opacity\n'+ settings['name'][0].title(), title_size=txt_slider_size2)
    if len(meshes) >=2:
        vp.addSlider2D(sliderAlphaMeshOut2, xmin=0, xmax=0.99, value=0.01,
               pos=[(0.92,0.40), (0.92,0.50)], c=settings['color'][1], 
               title='Opacity\n'+ settings['name'][1].title(), title_size=txt_slider_size2)
    if len(meshes) >=3:
        vp.addSlider2D(sliderAlphaMeshOut3, xmin=0, xmax=0.99, value=0.01,
               pos=[(0.72,0.25), (0.72,0.35)], c=settings['color'][2],
               title='Opacity\n'+ settings['name'][2].title(), title_size=txt_slider_size2)
    if len(meshes) >=4:
        vp.addSlider2D(sliderAlphaMeshOut4, xmin=0, xmax=0.99, value=0.01,
               pos=[(0.72,0.40), (0.72,0.50)], c=settings['color'][3], 
               title='Opacity\n'+ settings['name'][3].title(), title_size=txt_slider_size2)
        
    text = filename+'\n\nDefine plane position to '+txt.upper()+'. \nClose the window when done'
    txt = vedo.Text2D(text, c=txt_color, font=txt_font, s=txt_size)
    vp.show(meshes, plane, lbox, txt, viewup='y', zoom=zoom, interactive=True)

    return plane, normal, rotX, rotY, rotZ

def modify_disc(filename, txt, mesh, option,  
                    def_pl= {'pl_normal': (0,1,0), 'pl_centre': []}, 
                    radius = 60, height = 0.45, 
                    color = 'purple', res = 300, zoom = 0.5):
    """
    Function that shows a plot so that the user can define a cylinder (disc)

    """
    
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds()
    x_size = xmax - xmin; y_size = ymax - ymin; z_size = zmax - zmin
    
    xval = sorted([xmin-0.3*x_size,xmax+0.3*x_size])
    yval = sorted([ymin-0.3*y_size,ymax+0.3*y_size])
    zval = sorted([zmin-0.3*z_size,zmax+0.3*z_size])

    if len(def_pl['pl_centre']) == 0:
        centre = (x_size/2+xmin, ymin, z_size/2+zmin)
    else: 
        centre = def_pl['pl_centre']
    normal = def_pl['pl_normal']

    rotX = [0];    rotY = [0];    rotZ = [0]

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
        mesh.alpha(valueAlpha)

    lbox = vedo.LegendBox([mesh], font=leg_font, width=leg_width, padding=1)
    vp = vedo.Plotter(N=1, axes=8)
    vp.add_icon(logo, pos=(0.85,0.75), size=0.10)
    
    cyl_test = vedo.Cylinder(pos = centre, r = radius, height = height, 
                             axis = normal, c = color, cap = True, res = res)
    sph_cyl = vedo.Sphere(pos = centre, r=4, c='gold')

    if option[0]: #sliderX
        vp.addSlider2D(sliderX, xval[0], xval[1], value=centre[0],
                    pos=[(0.1,0.15), (0.3,0.15)], title='- > x position > +', 
                    c='crimson', title_size=txt_slider_size)
    if option[1]: #sliderY
        vp.addSlider2D(sliderY, yval[0], yval[1], value=centre[1],
                    pos=[(0.4,0.15), (0.6,0.15)], title='- > y position > +', 
                    c='dodgerblue', title_size=txt_slider_size)
    if option[2]: #sliderZ
        vp.addSlider2D(sliderZ, zval[0], zval[1], value=centre[2],
                    pos=[(0.7,0.15), (0.9,0.15)], title='- > z position > +', 
                    c='limegreen', title_size=txt_slider_size)
    if option[3]: #sliderRotX
        vp.addSlider2D(sliderRotX, -1, +1, value=0,
                    pos=[(0.1,0.05), (0.3,0.05)], title='- > x rotation > +', 
                    c='deeppink', title_size=txt_slider_size)
    if option[4]: #sliderRotY
        vp.addSlider2D(sliderRotY, -1, +1, value=0,
                    pos=[(0.4,0.05), (0.6,0.05)], title='- > y rotation > +', 
                    c='gold', title_size=txt_slider_size)
    if option[5]: #sliderRotZ
        vp.addSlider2D(sliderRotZ, -1, +1, value=0,
                    pos=[(0.7,0.05), (0.9,0.05)], title='- > z rotation > +', 
                    c='teal', title_size=txt_slider_size)

    vp.addSlider2D(sliderAlphaMeshOut, xmin=0.01, xmax=0.99, value=0.01,
               pos=[(0.95,0.25), (0.95,0.45)], c="blue", title="Mesh Opacity")

    text = filename+'\nDefine disc position to '+txt+'.\nMake sure it cuts all the tissue effectively separating it into individual segments.\nClose the window when done.\n[Note: The initially defined disc radius is '+str(radius)+'um.\nIf you are not happy with the disc radius, you will be able to modify it just before proceeding with the cut].'
    txt = vedo.Text2D(text, c=txt_color, font=txt_font, s=txt_size)
    vp.show(mesh, cyl_test, sph_cyl, lbox, txt, azimuth=45, elevation=45, zoom=zoom, interactive=True)

    return cyl_test, sph_cyl, rotX, rotY, rotZ

def modify_cube(filename, txt, mesh, orient_cube, centre, option, zoom = 0.5):
    """
    Function that shows a plot so that the user can define a cylinder (disc)

    """
    xmin, xmax, ymin, ymax, zmin, zmax = orient_cube.bounds()
    x_size = xmax - xmin; y_size = ymax - ymin; z_size = zmax - zmin

    xval = sorted([xmin-0.3*x_size,xmax+0.3*x_size])
    yval = sorted([ymin-0.3*y_size,ymax+0.3*y_size])
    zval = sorted([zmin-0.3*z_size,zmax+0.3*z_size])

    rotX = [0];    rotY = [0];    rotZ = [0]

    # Functions to move and rotate cube
    def sliderX(widget, event):
        valueX = widget.GetRepresentation().GetValue()
        orient_cube.x(valueX)

    def sliderY(widget, event):
        valueY = widget.GetRepresentation().GetValue()
        orient_cube.y(valueY)

    def sliderZ(widget, event):
        valueZ = widget.GetRepresentation().GetValue()
        orient_cube.z(valueZ)

    def sliderRotX(widget, event):
        valueRX = widget.GetRepresentation().GetValue()
        rotX.append(valueRX)
        orient_cube.rotateX(valueRX, rad=False)

    def sliderRotY(widget, event):
        valueRY = widget.GetRepresentation().GetValue()
        rotY.append(valueRY)
        orient_cube.rotateY(valueRY, rad=False)

    def sliderRotZ(widget, event):
        valueRZ = widget.GetRepresentation().GetValue()
        rotZ.append(valueRZ)
        orient_cube.rotateZ(valueRZ, rad=False)
    
    def sliderAlphaMeshOut(widget, event):
        valueAlpha = widget.GetRepresentation().GetValue()
        orient_cube.alpha(valueAlpha)

    lbox = vedo.LegendBox([mesh], font=leg_font, width=leg_width, padding=1)
    vp = vedo.Plotter(N=1, axes=8)
    vp.add_icon(logo, pos=(0.85,0.75), size=0.10)

    
    if option[0]: #sliderX
        vp.addSlider2D(sliderX, xval[0], xval[1], value=centre[0],
                    pos=[(0.1,0.15), (0.3,0.15)], title='- > x position > +', 
                    c='crimson', title_size=txt_slider_size)
    if option[1]: #sliderY
        vp.addSlider2D(sliderY, yval[0], yval[1], value=centre[1],
                    pos=[(0.4,0.15), (0.6,0.15)], title='- > y position > +', 
                    c='dodgerblue', title_size=txt_slider_size)
    if option[2]: #sliderZ
        vp.addSlider2D(sliderZ, zval[0], zval[1], value=centre[2],
                    pos=[(0.7,0.15), (0.9,0.15)], title='- > z position > +', 
                    c='limegreen', title_size=txt_slider_size)
    if option[3]: #sliderRotX
        vp.addSlider2D(sliderRotX, -1, +1, value=0,
                    pos=[(0.1,0.05), (0.3,0.05)], title='- > x rotation > +', 
                    c='deeppink', title_size=txt_slider_size)
    if option[4]: #sliderRotY
        vp.addSlider2D(sliderRotY, -1, +1, value=0,
                    pos=[(0.4,0.05), (0.6,0.05)], title='- > y rotation > +', 
                    c='gold', title_size=txt_slider_size)
    if option[5]: #sliderRotZ
        vp.addSlider2D(sliderRotZ, -1, +1, value=0,
                    pos=[(0.7,0.05), (0.9,0.05)], title='- > z rotation > +', 
                    c='teal', title_size=txt_slider_size)

    vp.addSlider2D(sliderAlphaMeshOut, xmin=0.01, xmax=0.99, value=0.01,
               pos=[(0.95,0.25), (0.95,0.45)], c="blue", title="Mesh Opacity")

    text = filename+'\nDefine cube position to '+txt+'.\nClose the window when done.\n[Note: The initially defined cube position and orientation is the same set for the STACK.'
    txt = vedo.Text2D(text, c=txt_color, font=txt_font, s=txt_size)
    vp.show(mesh, orient_cube, lbox, txt, azimuth=45, elevation=45, zoom=zoom, interactive=True)

    return orient_cube, rotX, rotY, rotZ

#%% - Mesh Operations
def get_pts_at_plane (points, pl_normal, pl_centre, tol=2, addData = []):
    """
    Function to get points within mesh at a certain plane position to create kspline

    """
    pts_cut = []
    data_cut = []

    d = pl_normal.dot(pl_centre)
    d_range = [d-tol, d+tol]
    d_range.sort()

    for i, pt in enumerate(points):
        d_pt = pl_normal.dot(pt)
        if d_pt>d_range[0] and d_pt<d_range[1]:
            pts_cut.append(pt)
            if addData != []:
                data_cut.append(addData[i])

    pts_cut = np.asarray(pts_cut)
    data_cut = np.asarray(data_cut)

    return pts_cut, data_cut

def get_cells_within_planes(controller, organ, cells_position, cut): 

    sphs_pos = cells_position[["Position X", "Position Y", "Position Z"]]
    sphs_pos_tuple = list(sphs_pos.itertuples(index=False, name=None))
    color_segm = organ.mC_settings['setup']['segm_mC'][cut]['colors']

    #Classify cells based on defined plane(s)
    sp_pos_classes = []
    for n, plane in enumerate(organ.mC_settings['wf_info']['segments_cells'][cut]['cut_info']): 
        sp_cells_class = np.empty(len(sphs_pos_tuple), dtype='object')
        pl_normal = organ.mC_settings['wf_info']['segments_cells'][cut]['cut_info'][plane]['pl_normal']
        pl_centre = organ.mC_settings['wf_info']['segments_cells'][cut]['cut_info'][plane]['pl_centre']
        
        # Find all the d values of pix_um
        d_pts = np.dot(np.subtract(sphs_pos_tuple,np.array(pl_centre)),np.array(pl_normal))
        pos_A = np.where(d_pts >= 0)[0]
        pos_B = np.where(d_pts < 0)[0]
    
        # Remove the points that belong to the other side
        sp_cells_class[pos_A] = 1+n
        sp_cells_class[pos_B] = 2+n

        sp_pos_classes.append(sp_cells_class)
        
        if len(sp_pos_classes)>1:
            #Here decide where each cell goes 
            pass
        else: 
            sp_pos_classes = sp_pos_classes[0]

    cells_position = organ.cellsMC['chA'].assign_class(cells_position, sp_pos_classes, col_name = 'Segment-'+cut)
    #Get mean spheres so user can define which segment is which
    # - Create mean spheres
    sphs_mean = create_mean_sphs(cells_position, sp_pos_classes, color_segm, 
                                 mtype = 'segm', cut=cut)    
    #Classify mean spheres
    sphs_classf = get_segm_mean_sph(controller, organ, cut, sphs_mean, color_segm)

    final_class = {}
    for num_class in list(set(sp_pos_classes)): 
        #Find group name 
        for segm in sphs_classf: 
            if num_class == int(sphs_classf[segm]['class'].split('Group No.')[1]): 
                final_class[num_class] = {'segm': segm, 
                                        'name': sphs_classf[segm]['name'], 
                                        'color': color_segm[segm]}

    segm_class = np.empty(len(sphs_pos_tuple), dtype='object')
    color_class = np.empty(len(sphs_pos_tuple), dtype='object')
    for n, val in enumerate(sp_pos_classes): 
        segm_class[n] = final_class[val]['segm']+':'+final_class[val]['name']
        color_class[n] = color_segm[final_class[val]['segm']]

    col_name = 'Segment-'+cut
    cells_position = organ.cellsMC['chA'].assign_class(cells_position, segm_class, col_name = col_name)
    cells_out = organ.cellsMC['chA'].colour_cells(sphs_pos = cells_position, 
                                                    color_class = color_class)
    segm_class = cells_position[col_name]

    organ.cellsMC['chA'].cells = cells_out

    organ.cellsMC['chA'].save_cells(cells_position)
    
    return cells_out, segm_class
        
def create_mean_sphs(cells_position, sp_pos_classes, color_segm, mtype, cut):

    if mtype == 'segm': 
        col_name = 'Segment-'+cut

    sphs_mean = []
    colors_init = palette_rbg('Set1', len(color_segm))
    all_given_num = list(set(sp_pos_classes))
    all_given_num.sort()
    for given_num in all_given_num:
        filt_segm = cells_position[cells_position[col_name] == given_num] 
        filt_segm = filt_segm[["Position X", "Position Y", "Position Z"]]
        mean = list(filt_segm.mean())

        sph_mean = vedo.Sphere(pos = mean, r=10, c=colors_init[given_num-1])
        sph_mean.legend('Group No.'+str(given_num))
        sphs_mean.append(sph_mean)

    return sphs_mean

def get_segm_mean_sph(controller, organ, cut, sphs_mean, color_segm): 

    segm_names = organ.mC_settings['setup']['segm_mC'][cut]['name_segments']
    colors = [color_segm[key] for key in color_segm]
    names = [segm_names[key] for key in segm_names]
    sph_classif = [[] for key in segm_names]
    
    done = False
    while not done:     
        # vols, names_iso, colors_vols
        iso_vols, vol_settings = organ.organ_vol_iso()

        def func(evt):
            if not evt.actor:
                return
            if isinstance(evt.actor, vedo.shapes.Sphere): 
                group_no = evt.actor.info['legend']
                group_color = evt.actor.color()
                group_color = np.round(evt.actor.color()*255,0)
                print(group_no, '- group_color:',group_color, type(group_color))
                is_in_list = np.any(np.all(group_color == colors, axis=1))
                if is_in_list:
                    bool_list = np.all(group_color == colors, axis=1)
                    ind = np.where(bool_list == True)[0][0]#n_out#  
                    ind_plus1 = ind+1
                    new_ind = (ind_plus1)%len(colors)
                    new_color = colors[new_ind]
                    evt.actor.color(new_color)
                    sph_classif[ind].remove(group_no)
                    sph_classif[new_ind].append(group_no)
                    msg.text("You clicked Sphere "+group_no+' to classify it as '+names[new_ind].upper())
                else: 
                    evt.actor.color(colors[0])
                    sph_classif[0].append(group_no)
                    msg.text("You clicked Sphere "+group_no+' to classify it as '+names[0].upper())
                print(sph_classif)
            
        mks = []; sym = ['o']*len(segm_names)
        for segm in segm_names:
            mks.append(vedo.Marker('*').c(color_segm[segm]).legend(segm_names[segm]))
        
        if len(iso_vols) >= 1: 
            def sliderAlphaMeshOut(widget, event):
                valueAlpha = widget.GetRepresentation().GetValue()
                iso_vols[0].alpha(valueAlpha)
        if len(iso_vols) >= 2: 
            def sliderAlphaMeshOut2(widget, event):
                valueAlpha = widget.GetRepresentation().GetValue()
                iso_vols[1].alpha(valueAlpha)
        if len(iso_vols) >= 3: 
            def sliderAlphaMeshOut3(widget, event):
                valueAlpha = widget.GetRepresentation().GetValue()
                iso_vols[2].alpha(valueAlpha)

        txA = 'Instructions: Click each coloured sphere until it is coloured according to the segment it belongs to.'
        txB = '\n[Note: Colours will loop for each sphere as you click it ]'
        txt0 = vedo.Text2D(txA+txB, c=txt_color, font=txt_font, s=txt_size)
        lb = vedo.LegendBox(mks, markers=sym, font=txt_font, 
                            width=leg_width/1.5, height=leg_height/1.5)
        
        txt_slider_size2 = 0.7
        msg = vedo.Text2D("", pos="bottom-center", c=txt_color, 
                          font=txt_font, s=txt_size, bg='red', alpha=0.2)
        
        plt = vedo.Plotter(axes=5, bg='white')
        plt.add_icon(logo, pos=(0.1,1), size=0.25)
        if len(iso_vols) >= 1: 
            plt.addSlider2D(sliderAlphaMeshOut, xmin=0, xmax=0.99, value=0.05,
                pos=[(0.92,0.25), (0.92,0.35)], c= vol_settings['color'][0], 
                title='Opacity\n'+ vol_settings['name'][0].title(), title_size=txt_slider_size2)
        if len(iso_vols) >=2:
            plt.addSlider2D(sliderAlphaMeshOut2, xmin=0, xmax=0.99, value=0.05,
                pos=[(0.92,0.40), (0.92,0.50)], c=vol_settings['color'][1], 
                title='Opacity\n'+ vol_settings['name'][1].title(), title_size=txt_slider_size2)
        if len(iso_vols) >=3:
            plt.addSlider2D(sliderAlphaMeshOut3, xmin=0, xmax=0.99, value=0.05,
                pos=[(0.72,0.25), (0.72,0.35)], c=vol_settings['color'][2],
                title='Opacity\n'+ vol_settings['name'][2].title(), title_size=txt_slider_size2)
            
        plt.add_callback('mouse click', func)
        plt.show(iso_vols, sphs_mean, msg, lb, txt0, zoom=1.2)

        done_sp = [len(sp_class) == 1 for sp_class in sph_classif]
          
        if all(done_sp): 
            done = True
        else:
            controller.main_win.win_msg('*Please classify one sphere per segment.')

    sphs_classf = {}
    for n, segm in enumerate(segm_names.keys()):
        sphs_classf[segm] = {'name': segm_names[segm], 
                             'class': sph_classif[n][0]}

    return sphs_classf

def modify_cell_class(organ, cells, cut, cells_class): 

    segm_names = organ.mC_settings['setup']['segm_mC'][cut]['name_segments']
    segm_colors = organ.mC_settings['setup']['segm_mC'][cut]['colors']

    segm_no = []; colors = []
    for segm in segm_names: 
        segm_no.append(segm+':'+segm_names[segm])
        colors.append(segm_colors[segm])
    
    iso_vols, vol_settings = organ.organ_vol_iso()

    def classify_cells_into_chambers(evt):
        if not evt.actor: 
            return
        if isinstance(evt.actor, vedo.shapes.Sphere): 
            sil = evt.actor.silhouette().lineWidth(20).c('red')
            #Get cell number (cell id)
            cell_no = evt.actor.name
            print(cell_no)
            cell_no = int(cell_no.split('.')[-1])
            current_class = cells_class.at[cell_no]
            for n, ss in enumerate(segm_no): 
                if current_class == ss: 
                    ind = n
                    break
            # ind = np.where(segm_no == current_class)[0]
            ind_plus1 = ind+1
            new_ind = (ind_plus1)%len(colors)
            new_color = colors[new_ind]
            evt.actor.color(new_color)
            cells_class.at[cell_no] = segm_no[new_ind]
            msg.text("Cell "+ evt.actor.name +' is now classified as part of '+segm_no[new_ind])
            plt.remove(silcont.pop()).add(sil)
            silcont.append(sil)

    silcont = [None]

    if len(iso_vols) >= 1: 
        def sliderAlphaMeshOut(widget, event):
            valueAlpha = widget.GetRepresentation().GetValue()
            iso_vols[0].alpha(valueAlpha)
    if len(iso_vols) >= 2: 
        def sliderAlphaMeshOut2(widget, event):
            valueAlpha = widget.GetRepresentation().GetValue()
            iso_vols[1].alpha(valueAlpha)
    if len(iso_vols) >= 3: 
        def sliderAlphaMeshOut3(widget, event):
            valueAlpha = widget.GetRepresentation().GetValue()
            iso_vols[2].alpha(valueAlpha)
    
    txt_slider_size2 = 0.7
    msg = vedo.Text2D("", pos="bottom-center", c=txt_color, font=txt_font, s=txt_size, bg='red', alpha=0.2)

    plt = vedo.Plotter(N=1, axes=1)
    plt.add_icon(logo, pos=(0.1,1), size=0.25)
    if len(iso_vols) >= 1: 
        plt.addSlider2D(sliderAlphaMeshOut, xmin=0, xmax=0.99, value=0.05,
               pos=[(0.92,0.25), (0.92,0.35)], c= vol_settings['color'][0], 
               title='Opacity\n'+ vol_settings['name'][0].title(), title_size=txt_slider_size2)
    if len(iso_vols) >=2:
        plt.addSlider2D(sliderAlphaMeshOut2, xmin=0, xmax=0.99, value=0.05,
               pos=[(0.92,0.40), (0.92,0.50)], c=vol_settings['color'][1], 
               title='Opacity\n'+ vol_settings['name'][1].title(), title_size=txt_slider_size2)
    if len(iso_vols) >=3:
        plt.addSlider2D(sliderAlphaMeshOut3, xmin=0, xmax=0.99, value=0.05,
               pos=[(0.72,0.25), (0.72,0.35)], c=vol_settings['color'][2],
               title='Opacity\n'+ vol_settings['name'][2].title(), title_size=txt_slider_size2)
        
    plt.addCallback('mouse click', classify_cells_into_chambers)
    plt.show(cells, iso_vols, msg, zoom=1.2)

    return cells, cells_class

def count_cells(organ, cells_position, cut, mtype, col_name): 

    if mtype == 'segm': 
        settings_name = 'name_segments'
        mC_name = 'mC_segm'

    names = organ.mC_settings['setup'][mtype+'_mC'][cut][settings_name]

    segm_no = []; 
    for name in names: 
        segm_no.append(name+':'+names[name])

    counts_cells = cells_position[col_name].value_counts()
    if mC_name not in organ.mC_settings['measure']: 
        organ.mC_settings['measure'][mC_name] = {}
    if cut not in organ.mC_settings['measure'][mC_name]:
        organ.mC_settings['measure'][mC_name][cut] = {}
    
    for sp in segm_no: 
        count_sp = counts_cells[sp]
        print(sp, ':', count_sp)
        if sp not in organ.mC_settings['measure'][mC_name][cut]:
            organ.mC_settings['measure'][mC_name][cut][sp] = {}

        organ.mC_settings['measure'][mC_name][cut][sp]['total'] = count_sp
   
#%% - Math operations 
def new_normal_3DRot (normal, rotX, rotY, rotZ):
    '''
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

    '''

    ang_X = np.radians(sum(rotX))
    ang_Y = np.radians(sum(rotY))
    ang_Z = np.radians(sum(rotZ))
    # print(sum(rotX), sum(rotY), sum(rotZ))

    normal_rotX = (np.dot(rotation_matrix(axis = [1,0,0], theta = ang_X), normal))
    normal_rotY = (np.dot(rotation_matrix(axis = [0,1,0], theta = ang_Y), normal_rotX))
    normal_rotZ = (np.dot(rotation_matrix(axis = [0,0,1], theta = ang_Z), normal_rotY))

    return normal_rotZ

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

def rotate_plane2im (pl_centre_o, pl_normal_o, chNS = False):
    """
    Function that rotates the planes defined in the surface reconstructions to the images mask

    """
    pl_centre = copy.deepcopy(pl_centre_o)
    pl_normal = copy.deepcopy(pl_normal_o)
    if not chNS:
        pl_im_centre = (np.dot(rotation_matrix(axis = [0,0,1], theta = np.radians(90)), pl_centre))
        pl_im_normal = (np.dot(rotation_matrix(axis = [0,0,1], theta = np.radians(90)), pl_normal))
    else:
        pl_im_centre = (np.dot(rotation_matrix(axis = [1,0,0], theta = np.radians(90)), pl_centre))
        pl_im_normal = (np.dot(rotation_matrix(axis = [1,0,0], theta = np.radians(90)), pl_normal))

    plane_im = vedo.Plane(pos=pl_im_centre,normal=pl_im_normal, s=(500,500)).color("blue").alpha(0.5).legend('Rotated plane to images')
    
    pl_dict = {'pl_normal': pl_im_normal,
               'pl_centre': pl_im_centre}
    
    return plane_im, pl_dict

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

def order_pts(points):
    """
    Function that returns an ordered array of points
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

def get_interpolated_pts(points, nPoints):
    """
    Function that interpolates input points

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

#%% - morphoCell
def find_dist_btw_all_cells(df_segm): 

    sphs_all = df_segm[df_segm['deleted'] == 'NO'] 
    sphs_pos_val = sphs_all[['Position X', 'Position Y', 'Position Z']].to_numpy()
    dists = cdist(sphs_pos_val, sphs_pos_val)

    return dists

def get_N_closest_cells(n_closestCells, count, dists, df_cellsChamber):
    
    cell_ids = list(df_cellsChamber.index)
    
    # Get the n closest cells per cell and record its distance
    index_closest_cells = []
    min_dist_to_cells = []
    for cell in range(count):
        dist_cell = dists[:,cell]
        index = np.argsort(dist_cell)
        # print(index[0:n_closestCells])
        list_index = []
        for ind in index[0:n_closestCells]:
            list_index.append(cell_ids[ind])
        # print(list_index)
        index_closest_cells.append(list_index)
        min_dist_to_cells.append(np.sort(dist_cell)[1:n_closestCells])
        
    return index_closest_cells, min_dist_to_cells

def option_cluster_number(filename, n_closestCells, avg_min_dist, index_closest_cells):

    num_cells = len(avg_min_dist)
    # Sort the average distance to sort clustering 
    index_sortedCells = list(np.argsort(avg_min_dist))*2

    num_clusters = []
    all_selected_index = []
    all_firstCell_sel_index = []
    all_group_sel_index = []
    all_indexes_org_cells = []
    
    for starting_cell in list(range(num_cells)):
        # print(starting_cell)
        selected_index = []
        firstCell_sel_index = []
        group_sel_index = []
        #Cluster cells in groups with n neighbours 
        num_group = 0
        indexes_org_cells = index_sortedCells[starting_cell:starting_cell+num_cells]
        for j, org_ind in enumerate(indexes_org_cells):
            # print(j, org_ind)
            # Cells indexes of revising group
            indg = index_closest_cells[org_ind]
            first_cell_index = indg[0]
            # print('Cells indexes of revising group (indg):', indg)
            not_intersect = set(firstCell_sel_index)^set(indg)
            # print('not_intersect',not_intersect)
            if len(not_intersect) == len(firstCell_sel_index)+ n_closestCells:
                firstCell_sel_index.append(first_cell_index)
                group_sel_index.append(indg)
                # print('-> First cell number ', first_cell_index, ' added!')
                # print('-  Cells indexes of revising group (indg):', indg)
                for k, indx in enumerate(indg):
                    selected_index.append(indx)
                    # print(indx)
                num_group+=1
                
        num_clusters.append(num_group)
        all_selected_index.append(selected_index)
        all_firstCell_sel_index.append(firstCell_sel_index)
        all_group_sel_index.append(group_sel_index)
        all_indexes_org_cells.append(indexes_org_cells)
        
    #Get the indexes of the groups that have more clusters in them
    index_cells_max_num_clusters = [i for i, j in enumerate(num_clusters) if j == max(num_clusters)]
    selected_index_max_cluster = [all_selected_index[i] for i in index_cells_max_num_clusters]
    firstCell_max_clusters = [all_firstCell_sel_index[i] for i in index_cells_max_num_clusters]
    group_max_clusters = [all_group_sel_index[i] for i in index_cells_max_num_clusters]
    indexes_org_cells_max_clusters = [all_indexes_org_cells[i] for i in index_cells_max_num_clusters]
    
    # #Compare now all the possible groups/clusters to make sure they are different
    # num_possible_clusters = 0
    # central_cells_of_all_diff_selected_clusters = []
    # selected_cells_of_all_diff_selected_clusters = []
    # indexes_org_cells_of_all_diff_selected_clusters = []
    # for m, sel_ind_max_cluster, ind_org_max_cluster in zip(count(), selected_index_max_cluster, indexes_org_cells_max_clusters):
    #     central_max = sel_ind_max_cluster[::n_closestCells]
    #     diff_sets = []
    #     for f, central in enumerate(central_cells_of_all_diff_selected_clusters):
    #         cells_diff = set(central_max)^set(central)
    #         diff_sets.append(len(cells_diff))
    #     # print('diff_sets:',diff_sets)
    #     if not 0 in diff_sets or m == 0:
    #         # print('added:', m, '-', central_max[0])
    #         num_possible_clusters +=1
    #         central_cells_of_all_diff_selected_clusters.append(central_max)
    #         indexes_org_cells_of_all_diff_selected_clusters.append(ind_org_max_cluster)
    #         selected_cells_of_all_diff_selected_clusters.append(sel_ind_max_cluster)
    
    # print('- Number of possible group options (maximizing cluster number): ', len(indexes_org_cells_of_all_diff_selected_clusters))
    
    # if len(indexes_org_cells_of_all_diff_selected_clusters) > 3:
    #     indexes_org_cells_of_all_diff_selected_clusters = indexes_org_cells_of_all_diff_selected_clusters[0:3]
    #     selected_cells_of_all_diff_selected_clusters = selected_cells_of_all_diff_selected_clusters[0:3]
    #     central_cells_of_all_diff_selected_clusters = central_cells_of_all_diff_selected_clusters[0:3]
    
    #return indexes_org_cells_of_all_diff_selected_clusters, selected_cells_of_all_diff_selected_clusters#, central_cells_of_all_diff_selected_clusters
    return firstCell_max_clusters, group_max_clusters, indexes_org_cells_max_clusters, selected_index_max_cluster

def cluster_segm_cells(organ, n_closestCells, selected_cells, mindist2Cells, avg_min_dist, 
                        df_cellsChamber, chamber, plot):
    
    #Add names of vols to plot
    vols, settings = organ.organ_vol_iso()

    # Define cluster colours 
    centre_cell = ['maroon','darkorange','goldenrod','yellowgreen','limegreen','aqua','skyblue','mediumblue',
                   'blueviolet','magenta', 'deeppink', 'sandybrown','dimgray','crimson','coral','gold','lawngreen',
                   'palegreen','teal','dodgerblue','navy','purple', 'mediumvioletred', 'pink','sienna', 'silver']*150
    directions = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1),(1,0,1),(1,1,1)]*150
    
    df_clusters_ALL = []
    sphs_distColour_ALL = []
    sil_distColour_ALL = []
    for sel_cells in selected_cells:
        list_ind_cellsCentre = []
        list_ind_cellsNeighbours = []
        list_averageDist = []
        list_allDistances = []
        num_group = 0
        sphs_distColour = create_segm_cells_from_df(df_cellsChamber, chamber)
        sil_group = [] 
        cell_groups = [sel_cells[x:x+n_closestCells] for x in range(0, len(sel_cells), n_closestCells)]
        for cell_group in cell_groups:
            # print(cell_group)
            for i, cell_no in enumerate(cell_group):
                #Find the sphere whose name is cell_no
                for n, cell in enumerate(sphs_distColour):
                    if int(cell.name.split('Cell Nr.')[1]) == cell_no:
                        ind_cell = n
                        break
                if i == 0: 
                    # print(cell_no)
                    cell.color(centre_cell[num_group]).scale(1.5)
                    list_ind_cellsCentre.append(cell_no)
                    list_averageDist.append(avg_min_dist[ind_cell])
                    list_allDistances.append(mindist2Cells[ind_cell])
                else: 
                    # sphs_distColour[cell_no].color(neighbour_cell[num_group])
                    sil = cell.silhouette(directions[num_group]).lineWidth(6).c(centre_cell[num_group])
                    sil_group.append(sil)
            list_ind_cellsNeighbours.append(cell_group)
            num_group+=1
                
        df_clusters = pd.DataFrame()
        df_clusters['Cluster_No'] = list(range(len(list_ind_cellsCentre)))
        df_clusters['Index_CentralCell'] = list_ind_cellsCentre
        df_clusters['AverageDist'] = list_averageDist
        df_clusters['Index_NeighbourCells'] = list_ind_cellsNeighbours
        df_clusters['Distances']= list_allDistances
        df_clusters['No_NeighbourCells'] = n_closestCells-1
        df_clusters = df_clusters.set_index('Cluster_No')
        
        df_clusters_ALL.append(df_clusters)
        sphs_distColour_ALL.append(sphs_distColour)
        sil_distColour_ALL.append(sil_group)
        # print(num_group)
    
    if plot:
        n_options = len(selected_cells)
        if n_options > 6: 
            sphs_distColour_ALL = sphs_distColour_ALL[0:6]
            sil_distColour_ALL = sil_distColour_ALL[0:6]

        vp = vedo.Plotter(N=len(sphs_distColour_ALL), axes = 4)
        for i in range(len(sphs_distColour_ALL)):
            text = 'Clustering option No. '+str(i+1)
            txt = vedo.Text2D(text, c=txt_color, font=txt_font)
            if i != n_options-1:
                vp.show(vols[0].alpha(0.1), sphs_distColour_ALL[i], sil_distColour_ALL[i], txt, at=i, zoom = 1)
            else: 
                vp.show(vols[0].alpha(0.1), sphs_distColour_ALL[i], sil_distColour_ALL[i],txt, at=i, zoom = 1, interactive = True) 
            
    return df_clusters_ALL, sphs_distColour_ALL, sil_distColour_ALL

def create_segm_cells_from_df(df_cellsChamber, chamber):

    cell_ids = list(df_cellsChamber.index)
    deleted = df_cellsChamber['deleted']
    pos_x = df_cellsChamber['Position X']
    pos_y = df_cellsChamber['Position Y']
    pos_z = df_cellsChamber['Position Z']

    sphs = []
    for i, x, y, z, id, bool_del in zip(count(), pos_x, pos_y, pos_z, cell_ids, deleted):
        if bool_del == 'NO': 
            pos = (x, y, z)
            s = vedo.Sphere(r=2).pos(pos).color('black')
            s.name = f"Cell Nr.{id}"
            sphs.append(s)
        else: 
            print('sph deleted no:', id)
        
    return sphs

def extract_segm_IND(organ, cut, segm, n_closest_cells, segm_count, plot, process = 'IND'):

    cell_distances = organ.cellsMC['chA'].df_cells()
    df_segm = cell_distances[cell_distances['Segment-'+cut] == segm]

    # Find euclidian distance beteen all cells in chamber
    dists = find_dist_btw_all_cells(df_segm)

    # # Get the n closest cells per cell and record its distance
    index_closestCells, min_dist_to_cells = get_N_closest_cells(n_closest_cells, segm_count, dists, df_segm)

    # Get average minimum distance to the n closest cells 
    avg_min_dist = [np.mean(cell_dist) for cell_dist in min_dist_to_cells]
    
    #Cluster cells in groups with n neighbours 
    centre_cells, neighb_cells, ind_orgCells, selected_cells = option_cluster_number(organ.user_organName, n_closest_cells, 
                                                                                        avg_min_dist, index_closestCells)

    #Select clusters
    df_clusters_ALL, sphs_distColour_ALL, sil_distColour_ALL = cluster_segm_cells(organ, n_closest_cells, 
                                                                                  selected_cells, min_dist_to_cells, 
                                                                                  avg_min_dist, df_segm, segm, plot)
    if process == 'IND': 
        return df_clusters_ALL, sphs_distColour_ALL, sil_distColour_ALL
    else: #process == 'zones':
        return index_closestCells, min_dist_to_cells, avg_min_dist

def create_zone(organ, zone, df_zone): 

    zones = organ.mC_settings['setup']['zone_mC'][zone]['name_zones']
    sil_color = organ.mC_settings['setup']['zone_mC'][zone]['colors']
    directions = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1),(1,0,1),(1,1,1)]*3

    #First get the cells from the zone the user wants to use
    if organ.mC_settings['wf_info']['zone_cells'][zone]['use_segm']: 
        segm_cut_info = organ.mC_settings['wf_info']['zone_cells'][zone]['ind_segm']
        cut, names = segm_cut_info.split(':')
        #Get all names
        dict_names = organ.mC_settings['setup']['segm_mC'][cut]['name_segments']
        segm_names = {}
        for key, item in dict_names.items():
            segm_names[key] = key+':'+item
        #See if the segms selected are more than one
        all_names = names.split(', ')
        #Create an empty dictionary where values of all segm will be saved
        all_data_to_run = {}

        #Option 1: Use all segments for cut
        if len(all_names) > 1:
            #Using all segments
            for segm in segm_names: 
                all_data_to_run[segm] = {}
        else: 
            #Get the segm that will be used
            for ss in dict_names: 
                if all_names[0] == dict_names[ss]: 
                    segm = ss
                    break
            all_data_to_run[segm] = {} 

        sphs_all = []
        for ss in all_data_to_run: 
            sphs_all.append(organ.ind['mC_segm'][cut.title()][ss][0])
        sphs_allf = sum(sphs_all, [])
        sphs_zones = sphs_allf.copy()
        sil_zones = []
        num_cells_zone = {}
        for zz in zones: 
            num_cells_zone[zones[zz].strip()] = []

        #Modify cells according to info in df
        for ind, row in df_zone.iterrows():
            cell_no = row['Index_CentralCell']
            zone = row['Zone'].strip()
            print(ind, cell_no, zone)
            #Find cell within sphs
            found = False
            for sph in sphs_zones: 
                sph_name = sph.name.split('Cell Nr.')[1]
                if int(sph_name) == cell_no: 
                    found = True
                    break
            if found: 
                #Find zone number
                for zz in zones:
                    if zones[zz].strip() == zone: 
                        break 
                sil = sph.silhouette((1,0,0)).lineWidth(20).c(sil_color[zz]).alpha(1)
                sil_zones.append(sil)
                num_cells_zone[zone].append(cell_no)
    
    return sphs_zones, sil_zones

def select_cell_for_zones(organ, cut, zone, segm_names, all_data):
    
    iso_vols, vol_settings = organ.organ_vol_iso()
    zones = organ.mC_settings['setup']['zone_mC'][zone]['name_zones']
    sil_color = organ.mC_settings['setup']['zone_mC'][zone]['colors']
    bkgd = organ.mC_settings['wf_info']['zone_cells']['background']

    mks = []; sym = ['o']*len(zones)
    for zn in zones:
        mks.append(vedo.Marker('*').c(sil_color[zn]).legend(zones[zn]))
    
    sphs = 'sphs_distColour'
    index = 'index_closest_cells'
    avg_dist = 'avg_min_dist'
    min_dist = 'min_dist_to_cells'

    all_sphs = [all_data[segm][sphs] for segm in all_data]
    all_sphsf = sum(all_sphs,[])

    #Find the numbers of cells per segment
    df_cells = organ.cellsMC['chA'].df_cells()
    col_name = 'Segment-'+cut
    
    segm_cells = {}; segm_ids = {}
    for segm, uname in segm_names.items(): 
        segm_cells[segm] = df_cells[df_cells[col_name] == uname]
        segm_ids[segm] = list(segm_cells[segm].index)

    list_ind_cellsCentreB = []
    list_ind_cellsNeighboursB = []
    list_averageDistB = []
    list_allDistancesB = []
    list_region = []
    selected_cells = []
    sil_final =[]
    msg_zone = []
    
    for i, zz in zip(count(), zones): 
        print('>>>', i, zz)
        sphs_distColour_copy = all_sphsf.copy()
        cells_region = []
        silcont = []
        def select_cells_per_chamber_region(evt):
            if not evt.actor: return
            if isinstance(evt.actor, vedo.shapes.Sphere): 
                cell_no = evt.actor.name
                cell_no = int(cell_no.split('.')[-1])
                if cell_no not in selected_cells:
                    #If the cell wants fo be deselected
                    if cell_no in cells_region: 
                        ind_cell = cells_region.index(cell_no)
                        sil2rem = silcont[ind_cell]
                        plt.remove(sil2rem)
                        cells_region.pop(ind_cell)
                        msg.text('You removed '+evt.actor.name+' from '+zones[zz].upper()+'\n(Number of cells selected so far: '+str(len(cells_region))+')')
                        silcont.pop(ind_cell)
                    else: 
                        sil = evt.actor.silhouette().lineWidth(20).c(sil_color[zz]).alpha(1)
                        sil.legend('No.'+str(cell_no))
                        cells_region.append(cell_no)
                        msg.text('You seleced '+evt.actor.name+' as part of '+zones[zz].upper()+'\n(Number of cells selected so far: '+str(len(cells_region))+')')
                        plt.add(sil) #plt.remove(silcont.pop()).add(sil)
                        silcont.append(sil)
                    print('cells_region:',cells_region)
                    # print('silcont:', silcont)
                else: 
                    msg.text(evt.actor.name+' has already been selected for another zone. \n\t Please select another cell for -'+zones[zz])

        text = 'Select cells in the Zone: '+zones[zz]+'\n Instructions: \n-Click cells that are part of the zone of interest and a silhouette\nwith the corresponding zone color will appear around the selected cell\n-Play with the mesh opacity to be able to select cells that are within the isovolume \n-If you wish to deselect a cell, click the cell again\n-Close the window when you are done to continue selecting cells for the next zone.'
        txt = vedo.Text2D(text, c=txt_color, font= txt_font)
        lb = vedo.LegendBox(mks, markers=sym, font=txt_font, 
                    width=leg_width/1.5, height=leg_height)
        
        if len(iso_vols) >= 1: 
            def sliderAlphaMeshOut(widget, event):
                valueAlpha = widget.GetRepresentation().GetValue()
                iso_vols[0].alpha(valueAlpha)
        if len(iso_vols) >= 2: 
            def sliderAlphaMeshOut2(widget, event):
                valueAlpha = widget.GetRepresentation().GetValue()
                iso_vols[1].alpha(valueAlpha)
        if len(iso_vols) >= 3: 
            def sliderAlphaMeshOut3(widget, event):
                valueAlpha = widget.GetRepresentation().GetValue()
                iso_vols[2].alpha(valueAlpha)
        
        txt_slider_size2 = 0.7
        msg = vedo.Text2D("", pos="bottom-center", c=txt_color, font=txt_font, s=txt_size, bg='red', alpha=0.2)

        plt = vedo.Plotter(N=1, axes=1, bg=bkgd)
        plt.add_icon(logo, pos=(0.1,0.1), size=0.25)
        if len(iso_vols) >= 1: 
            plt.addSlider2D(sliderAlphaMeshOut, xmin=0, xmax=0.99, value=0.05,
                pos=[(0.92,0.25), (0.92,0.35)], c= vol_settings['color'][0], 
                title='Opacity\n'+ vol_settings['name'][0].title(), title_size=txt_slider_size2)
        if len(iso_vols) >=2:
            plt.addSlider2D(sliderAlphaMeshOut2, xmin=0, xmax=0.99, value=0.05,
                pos=[(0.92,0.40), (0.92,0.50)], c=vol_settings['color'][1], 
                title='Opacity\n'+ vol_settings['name'][1].title(), title_size=txt_slider_size2)
        if len(iso_vols) >=3:
            plt.addSlider2D(sliderAlphaMeshOut3, xmin=0, xmax=0.99, value=0.05,
                pos=[(0.72,0.25), (0.72,0.35)], c=vol_settings['color'][2],
                title='Opacity\n'+ vol_settings['name'][2].title(), title_size=txt_slider_size2)
            
        plt.addCallback('mouse click', select_cells_per_chamber_region)
        plt.show(all_sphsf, iso_vols, sil_final, txt, msg, lb, zoom=1.2)

        for i, cell_no in enumerate(cells_region): 
            #Find the segm the cell is part of
            for ss in segm_ids: 
                if cell_no in segm_ids[ss]: 
                    cell_is_part_of = ss
                    print('cell_is_part_of:', cell_is_part_of)
                    break
            #Find the index where the info of that cell is stored within lists of data
            for ii, aa in enumerate(all_data[cell_is_part_of][index]):
                if cell_no in aa: 
                    iif = ii
                    print(iif)
                    break

            list_ind_cellsCentreB.append(cell_no)
            list_ind_cellsNeighboursB.append(all_data[cell_is_part_of][index][iif])#index_closestCells[cell_no])
            list_averageDistB.append(all_data[cell_is_part_of][avg_dist][iif])#avg_min_dist[cell_no])
            list_allDistancesB.append(all_data[cell_is_part_of][min_dist][iif])#mindist2Cells[cell_no])
            list_region.append(zones[zz].strip())
            selected_cells.append(cell_no)
            sil_final.append(silcont[i])
        
        msg_zone.append(zones[zz]+':'+str(len(cells_region)))

    msg_zonef = ', '.join(msg_zone)
    txt_msg = ' -- Total number of cells selected per zone --\n' + msg_zonef
    msg = vedo.Text2D(txt_msg, pos="bottom-center", c=txt_color, font=txt_font, s=txt_size, bg='red', alpha=0.2)
    text3 = 'CHECKING: Final selection of cells for ALL zones'
    txt3= vedo.Text2D(text3, font=txt_font)
    vp = vedo.Plotter(axes = 1, bg=bkgd)
    vp.add_icon(logo, pos=(0.1,0.1), size=0.25)
    if len(iso_vols) >= 1: 
        vp.addSlider2D(sliderAlphaMeshOut, xmin=0, xmax=0.99, value=0.05,
            pos=[(0.92,0.25), (0.92,0.35)], c= vol_settings['color'][0], 
            title='Opacity\n'+ vol_settings['name'][0].title(), title_size=txt_slider_size2)
    if len(iso_vols) >=2:
        vp.addSlider2D(sliderAlphaMeshOut2, xmin=0, xmax=0.99, value=0.05,
            pos=[(0.92,0.40), (0.92,0.50)], c=vol_settings['color'][1], 
            title='Opacity\n'+ vol_settings['name'][1].title(), title_size=txt_slider_size2)
    if len(iso_vols) >=3:
        vp.addSlider2D(sliderAlphaMeshOut3, xmin=0, xmax=0.99, value=0.05,
            pos=[(0.72,0.25), (0.72,0.35)], c=vol_settings['color'][2],
            title='Opacity\n'+ vol_settings['name'][2].title(), title_size=txt_slider_size2)
    vp.show(sphs_distColour_copy, sil_final, iso_vols, lb, txt3, msg, zoom=1.2, interactive = True)#silcont

    df_zones = pd.DataFrame()
    df_zones['Index_CentralCell'] = list_ind_cellsCentreB
    df_zones['AverageDist'] = list_averageDistB
    df_zones['Index_NeighbourCells'] = list_ind_cellsNeighboursB
    df_zones['Distances']= list_allDistancesB
    df_zones['Zone']= list_region
        
    return df_zones, sphs_distColour_copy, sil_final

#%% Module loaded
print('morphoHeart! - Loaded funcMeshes')