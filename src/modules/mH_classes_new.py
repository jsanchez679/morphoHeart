'''
morphoHeart_classes

@author: Juliana Sanchez-Posada
'''
#%% ##### - Imports - ########################################################
import os
from datetime import datetime
import pathlib
from pathlib import Path
import numpy as np
from skimage import measure, io
import copy
import json
import collections
import pprint
import matplotlib.pyplot as plt
import flatdict
from typing import Union
from time import perf_counter
import vedo as vedo
from vedo import write
from scipy.interpolate import interpn#, splprep, splev
from itertools import count
import random
from skimage.draw import line_aa
import seaborn as sns

#%% morphoHeart Imports - ##################################################
from ..gui.gui_classes import *
from .mH_funcBasics import alert, make_Paths, make_tuples, get_by_path, set_by_path, rename_directory
from .mH_funcMeshes import (unit_vector, plot_organCLs, find_angle_btw_pts, new_normal_3DRot,
                            classify_segments_from_ext, create_subsegment, create_subsection, plot_grid, 
                            modify_cube)
from ..gui.config import mH_config

path_mHImages = mH_config.path_mHImages
# Load logo
path_logo = path_mHImages / 'logo-07.jpg'
logo = vedo.Picture(str(path_logo))

#%% Set default fonts and sizes for plots
txt_font = mH_config.txt_font
leg_font = mH_config.leg_font
leg_width = mH_config.leg_width
leg_height = mH_config.leg_height
txt_size = mH_config.txt_size
txt_color = mH_config.txt_color
txt_slider_size = mH_config.txt_slider_size

#%% ##### - Class definition - ###############################################
# Definition of class to save dictionary
class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, int): 
            return int(obj)
        elif isinstance(obj, np.int64): 
            return int(obj)
        elif isinstance(obj, float):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.float32):
            return obj.tolist()
        elif isinstance(obj, np.float64):
            return obj.tolist()
        elif isinstance(obj, pathlib.WindowsPath):
            return str(obj)
        elif isinstance(obj, pd.DataFrame): 
            return obj.to_json(orient='table')
        else:
            print('type object: ', type(obj), obj)
            return super(NumpyArrayEncoder, self).default(obj)

class Project(): 
    '''
    Project class: 
    A new project would include a group or groups of Organs
    that are part of the same experiment and will then be averaged and/or
    compared. From this a project can contain just one organ or 
    multiple organs with different genotypes and stages. The settings to
    process all the organs contained in a project are set up when starting a 
    new project and can be amended if needed as the organs are processed.
    '''
    def __init__(self, proj_dict:dict, new:bool):#
            
        def create_mHName(self):#
            '''
            func - create name for a morphoHeart project
            This function will assign the newly created project a name using a
            timestamp
            '''
            now_str = datetime.now().strftime('%Y%m%d%H%M')
            return 'mH_Proj-'+now_str

        if new:
            self.mH_projName = create_mHName(self)
            self.user_projName = proj_dict['name']#.replace(' ', '_')
            self.info = {'mH_projName': self.mH_projName,
                            'user_projName': self.user_projName,
                            'user_projNotes': proj_dict['notes'], 
                            'date_created' : proj_dict['date'],
                            'dirs': {}, 
                            'heart_default': proj_dict['heart_default']}
            
            self.analysis = proj_dict['analysis']
            self.dir_proj = Path(proj_dict['dir_proj'])
            self.organs = {}
            self.cellGroups = {}
            self.gui_custom_data = {'strain': [], 
                                    'stage': [],
                                    'genotype': [],
                                    'manipulation': [],
                                    'im_orientation': ['ventral', 'dorsal'], 
                                    'im_res_units': ['um', 'mm', 'cm', 'm']}
            try:
                if proj_dict['heart_default']: 
                    mH_config.heart_default = True
            except: 
                pass
        else: 
            load_dict = {'name': proj_dict['name'], 'dir': proj_dict['dir']}
            self.load_project(load_dict=load_dict)
    
    def load_project(self, load_dict:dict):#
        print('Loading project:', load_dict)
        proj_name = load_dict['name']
        proj_name_us = proj_name.replace(' ', '_')
        dir_res = Path(load_dict['dir'])
        jsonDict_name = 'mH_'+proj_name_us+'_project.json'
        json2open_dir = dir_res / 'settings' / jsonDict_name
        if json2open_dir.is_file():
            with open(json2open_dir, "r") as read_file:
                print(">> "+jsonDict_name+": Opening JSON encoded data")
                load_dict = json.load(read_file)
            
            load_dict = make_Paths(load_dict)
            
            tuple_keys = [['mH_settings','setup','chNS','ch_ext'], 
                          ['mH_settings','setup','chNS','ch_int'],]
            load_dict = make_tuples(load_dict, tuple_keys)
        
            self.info = load_dict['info']
            self.user_projName = load_dict['info']['user_projName']
            self.mH_projName = load_dict['info']['mH_projName']
            self.analysis = load_dict['analysis']
            
            self.mH_settings = load_dict['mH_settings']
            json_str = self.mH_settings['df_res'] 
            self.mH_settings['df_res'] = pd.read_json(json_str, orient='table')
            self.mH_channels = load_dict['mH_channels']
            self.mH_segments = load_dict['mH_segments']
            self.mH_sections = load_dict['mH_sections']
            # self.mH_param2meas = load_dict['mH_param2meas'] #[tuple(item) for item in load_dict['mH_param2meas']]
            
            self.mC_settings = load_dict['mC_settings']
            self.mC_channels = load_dict['mC_channels']
            self.mC_segments = load_dict['mC_segments']
            self.mC_sections = load_dict['mC_sections']
            # self.mC_param2meas = load_dict['mC_param2meas']
            
            self.workflow = load_dict['workflow']
            self.mH_methods = load_dict['mH_methods']
            self.mC_methods = load_dict['mC_methods']
            self.organs = load_dict['organs']

            self.gui_custom_data = load_dict['gui_custom_data']
          
            self.dir_proj = dir_res# load_dict['info']['dir_proj']
            self.dir_info = dir_res / 'settings' / jsonDict_name #load_dict['info']['dir_info']
            self.info['dir_proj'] = self.dir_proj
            self.info['dir_info'] = self.dir_info
            # print(self.dir_info)

        else: 
            print('>> Error: No project with name ',load_dict['name'],' was found!\n Directory: ',str(json2open_dir))
            alert('error_beep')

    def set_settings(self, settings:dict):#
        '''
        func - Create general project settings
        This function will be called when the user creates a new project and 
        fills information regarding the workflow for such project which will get into the 
        function packed as the mH_settings dictionary. 
        The output of this function will create an attribute to the project containing 
        most of the user settings except for the selected parameters. 
        '''
        # print('Settings to create new project:', settings)
        self.set_mH_settings(mH_settings = settings['mH']['settings'], 
                             mH_params = settings['mH']['params'])
        self.set_mC_settings(mC_settings = settings['mC']['settings'], 
                             mC_params = settings['mC']['params'])
        
    def set_mH_settings(self, mH_settings:dict, mH_params:dict):#
        
        if self.analysis['morphoHeart']:
            self.mH_settings = {}
            #Add setup dict containing all mH_settings for the new project
            self.mH_settings['setup'] = mH_settings
            self.mH_settings['setup']['keep_largest'] = {}
            self.mH_settings['setup']['alpha'] = {}

            #Add empty dict to save info related to the user-selected processes
            self.mH_settings['wf_info'] = {}

            #Add attribute with information of channels
            self.mH_channels = mH_settings['chs_all']

            #Add attribute with info regarding segments
            if isinstance(mH_settings['segm'], dict):
                self.mH_segments = {}
                for key in mH_settings['segm']: 
                    if 'Cut' in key:
                        self.mH_segments[key]= mH_settings['segm'][key]['name_segments']
            else: 
                 self.mH_segments = False

            #Add attribute with info regarding sections
            if isinstance(mH_settings['sect'], dict): 
                self.mH_sections = {}
                for key in mH_settings['sect']: 
                    if 'Cut' in key:
                        self.mH_sections[key]= mH_settings['sect'][key]['name_sections']
            else: 
                self.mH_sections = False
        else: 
            self.mH_settings = None
            self.mH_channels = None
            self.mH_segments = None
            self.mH_sections = None
            self.mH_methods = None

        self.clean_False(user_param = mH_params)
        self.set_mH_methods()

    def set_mC_settings(self, mC_settings:dict, mC_params:dict):# 
        
        if self.analysis['morphoCell']:
            self.mC_settings = {}
            #Add setup dict containing all mC_settings for the new project
            self.mC_settings['setup'] = mC_settings

            #Add empty dict to save info related to the user-selected processes
            self.mC_settings['wf_info'] = {}

            #Add attribute with information of channels
            self.mC_channels = mC_settings['name_chs']

            #Add attribute with info regarding segments
            if isinstance(mC_settings['segm_mC'], dict):
                self.mC_segments = {}
                for key in mC_settings['segm_mC']: 
                    if 'Cut' in key:
                        self.mC_segments[key]= mC_settings['segm_mC'][key]['name_segments']
            else: 
                 self.mC_segments = False

            #Add attribute with info regarding sections
            if isinstance(mC_settings['segm-sect_mC'], dict): 
                self.mC_sections = {}
                for key in mC_settings['sect_mC']: 
                    if 'Cut' in key:
                        self.mC_sections[key]= mC_settings['sect_mC'][key]['name_sections']
            else: 
                self.mC_sections = False

            #Add attribute with info regarding zones
            if isinstance(mC_settings['zone_mC'], dict): 
                self.mC_zones = {}
                for key in mC_settings['zone_mC']: 
                    if 'Cut' in key:
                        self.mC_zones[key]= mC_settings['zone_mC'][key]['name_zones']
            else: 
                self.mC_zones = False
            
            self.mC_settings['measure'] = {}

            self.set_mC_methods()

        else: 
            self.mC_settings = None
            self.mC_channels = None
            self.mC_segments = None
            self.mC_sections = None
            self.mC_zones = None
            self.mC_methods = None
        
        
    def clean_False(self, user_param:dict):#

        user_param_new = copy.deepcopy(user_param)
        for param in user_param: 
            for key in user_param[param]: 
                if not user_param[param][key]: 
                    user_param_new[param].pop(key, None)
        
        self.mH_settings['measure'] = user_param_new
    
    def set_mH_methods(self):#
        mH_param2meas = self.mH_settings['measure'] #self.mH_param2meas

        if len(mH_param2meas)>0: 
            methods = ['A-Set_Orientation', 'A-Create3DMesh','B-TrimMesh']
        
            if 'ball' in mH_param2meas.keys():
                methods.append('C-Centreline')
                methods.append('D-Ballooning')
            elif 'CL' in mH_param2meas.keys():
                methods.append('C-Centreline')

            if 'th_i2e' in mH_param2meas.keys():
                methods.append('D-Thickness_int>ext')
            if 'th_e2i' in mH_param2meas.keys():
                methods.append('D-Thickness_ext>int')

            if 'Vol(segm)' in mH_param2meas.keys() or 'SA(segm)' in mH_param2meas.keys() or 'Ellip(segm)' in mH_param2meas.keys():
                methods.append('E-Segments')
            if 'Vol(sect)' in mH_param2meas.keys() or 'SA(sect)' in mH_param2meas.keys(): 
                methods.append('E-Sections')
            if 'Vol(segm-sect)' in mH_param2meas.keys() or 'SA(segm-sect)' in mH_param2meas.keys(): 
                methods.append('E-Segments_Sections')
        else: 
            methods = []
        
        self.mH_methods = methods

    def set_mC_methods(self): 
        methods = ['A-Set_Orientation', 'A-SetExtraChs', 'A-CleanCells']
        if self.mC_segments != None and self.mC_segments != False: 
            methods.append('B-Segments')
        if self.mC_sections != None and self.mC_sections != False: 
            methods.append('B-Regions')
        if isinstance(self.mC_settings['setup']['segm-sect_mC'], dict):
            methods.append('B-Segments_Regions')
        if self.mC_zones != None: 
            methods.append('B-Zones')
        methods.append('C-Measure')
        
        self.mC_methods = methods

    def create_proj_dir(self):#
        self.dir_proj.mkdir(parents=True, exist_ok=True)
        if self.dir_proj.is_dir():
            self.info['dir_proj'] = self.dir_proj
        else: 
            print('>> Error: Project directory could not be created!\n>> Dir: '+str(self.dir_proj))
            alert('error_beep')
            
    def set_workflow(self):#
        '''
        This function will initialise the dictionary that will contain the workflow of the
        project. This workflow will be assigned to each organ that is part of the created project
        and will be updated in each organ as the user advances in the processing. 
        '''
        workflow = dict()
        dict_mH = dict()
        dict_mC = dict()
        mH_param2meas = self.mH_settings['measure']
       
        if self.analysis['morphoHeart']: 
            mH_channels = sorted(self.mH_channels)
    
            dict_ImProc = dict()
            dict_ImProc['Status'] = 'NI'
            dict_MeshesProc = dict()
    
            # Find the meas_param that include the extraction of a centreline
            item_centreline = [tuple(item.split('_')) for item in mH_param2meas['CL'].keys()]

            # Find the meas_param that include the extraction of mH_segments
            if 'Vol(segm)' in mH_param2meas:
                segm_vol = [item for item in mH_param2meas['Vol(segm)'].keys()]
            else: 
                segm_vol = []
            if 'SA(segm)' in mH_param2meas:
                segm_sa = [item for item in mH_param2meas['SA(segm)'].keys()]
            else: 
                segm_sa = []
            if 'Ellip(segm)' in mH_param2meas:
                segm_ellip = [item for item in mH_param2meas['Ellip(segm)'].keys()]
            else:
                segm_ellip = []

            segm_list = list(set(segm_vol) | set(segm_sa) | set(segm_ellip))
            if len(segm_list) > 0: 
                segm_list = [tuple(tup.split('_')) for tup in segm_list]
                cut_segm = sorted(list(set([tup for (tup,_,_,_) in segm_list])))
                ch_segm = sorted(list(set([tup for (_,tup,_,_) in segm_list])))
            else: 
                segm_list = []; cut_segm = []; ch_segm = []

            # Find the meas_param that include the extraction of mH_sections
            if 'Vol(sect)' in mH_param2meas:
                sect_vol = [item for item in mH_param2meas['Vol(sect)'].keys()]
            else: 
                sect_vol = []
            if 'SA(sect)' in mH_param2meas:
                sect_sa = [item for item in mH_param2meas['SA(sect)'].keys()]
            else: 
                sect_sa = []

            sect_list = list(set(sect_vol) | set(sect_sa))
            if len(sect_list) > 0: 
                sect_list = [tuple(item.split('_')) for item in sect_list]
                cut_sect = sorted(list(set([tup for (tup,_,_,_) in sect_list])))
                ch_sect=  sorted(list(set([tup for (_,tup,_,_) in sect_list])))
            else: 
                sect_list = []; cut_sect = []; ch_sect = []

            # Find the meas_param that include the extraction of mH_segm-sect
            if 'Vol(segm-sect)' in mH_param2meas:
                segm_sect_vol = [item for item in mH_param2meas['Vol(segm-sect)'].keys()]
            else: 
                segm_sect_vol = []
            if 'SA(segm-sect)' in mH_param2meas:
                segm_sect_sa = [item for item in mH_param2meas['SA(segm-sect)'].keys()]
            else: 
                segm_sect_sa = []

            segm_sect_list = list(set(segm_sect_vol) | set(segm_sect_sa))
            if len(segm_sect_list) > 0: 
                segm_sect_list = [tuple(item.split('_')) for item in segm_sect_list]
                scut_segm_sect = sorted(list(set([tup for (tup,_,_,_,_,_) in segm_sect_list])))
                rcut_segm_sect = sorted(list(set([tup for (_,tup,_,_,_,_) in segm_sect_list])))
                ch_segm_sect=  sorted(list(set([tup for (_,_,tup,_,_,_) in segm_sect_list])))
            else: 
                segm_sect_list = []; scut_segm_sect = []; rcut_segm_sect = []; ch_segm_sect = []

            # Find the meas_param that include the extraction of ballooning
            ball_opts = self.mH_settings['setup']['params'][5]['measure']
            # Find the meas_param that include the extraction of thickness
            item_thickness_intext = [tuple(item.split('_')) for item in mH_param2meas['th_i2e'].keys()]
            item_thickness_extint = [tuple(item.split('_')) for item in mH_param2meas['th_e2i'].keys()]
    
            dict_MeshesProc = {'Status' : 'NI'}
            for met in self.mH_methods:
                dict_MeshesProc[met] =  {'Status': 'NI'}
                               
            # Project status
            for ch in mH_channels:
                if 'A-Create3DMesh' in dict_MeshesProc.keys():
                    dict_MeshesProc['A-Set_Orientation'] = {'Status': 'NI', 
                                                            'Stack': 'NI', 
                                                            'ROI': 'NI'}
                    if 'NS' not in ch:
                        if self.mH_settings['setup']['mask_ch'][ch]:
                            mask_status = 'NI'
                        else: 
                            mask_status = 'N/A'
                        dict_ImProc[ch] = {'Status': 'NI',
                                            'A-MaskChannel': {'Status': mask_status},
                                            'B-CloseCont':{'Status': 'NI',
                                                            'Steps':{'A-Autom': {'Status': 'NI'},
                                                                    'B-Manual': {'Status': 'NI'},
                                                                    'C-CloseInOut': {'Status': 'NI'}}},
                                            'C-SelectCont':{'Status': 'NI'},
                                            'D-S3Create':{'Status': 'NI',
                                                        'Info': {'tiss':{'Status': 'NI'}, 
                                                                'int':{'Status': 'NI'}, 
                                                                'ext':{'Status': 'NI'}}}}
                        dict_ImProc[ch]['E-CleanCh'] = {'Status': 'NI',
                                                            'Info': {'tiss':{'Status': 'NI'}, 
                                                                    'int':{'Status': 'NI'}, 
                                                                    'ext':{'Status': 'NI'}}}
                        dict_ImProc[ch]['E-TrimS3'] = {'Status': 'NI',
                                                            'Info':{'tiss':{'Status': 'NI'}, 
                                                                    'int':{'Status': 'NI'},
                                                                    'ext':{'Status': 'NI'}}}
                    else: 
                        dict_ImProc[ch] = {'Status': 'NI',
                                            'D-S3Create':{'Status': 'NI'}} 
                        
            for nn, ch in enumerate(mH_channels):
                for process in ['A-Create3DMesh','B-TrimMesh','C-Centreline']:
                    if 'NS' not in ch:
                        if process != 'C-Centreline':
                            dict_MeshesProc[process][ch] = {}
                        for nnn, cont in enumerate(['tiss', 'int', 'ext']):
                            if process == 'A-Create3DMesh' or process == 'B-TrimMesh':
                                dict_MeshesProc[process][ch][cont] = {'Status': 'NI'}
                         
                            if process == 'C-Centreline' and 'C-Centreline' in dict_MeshesProc.keys():
                                if nn == 0 and nnn == 0: 
                                    dict_MeshesProc[process]['Status'] = 'NI'
                                    dict_MeshesProc[process]['SimplifyMesh'] = {'Status':'NI'}
                                    dict_MeshesProc[process]['vmtk_CL'] = {'Status':'NI'}
                                    dict_MeshesProc[process]['buildCL'] = {'Status':'NI'}
                                if (ch, cont, 'whole') in item_centreline:
                                    if ch not in dict_MeshesProc[process]['SimplifyMesh'].keys(): 
                                        dict_MeshesProc[process]['SimplifyMesh'][ch] = {}
                                        dict_MeshesProc[process]['vmtk_CL'][ch] = {}
                                        dict_MeshesProc[process]['buildCL'][ch] = {}
                                    dict_MeshesProc[process]['SimplifyMesh'][ch][cont] = {'Status': 'NI'}
                                    dict_MeshesProc[process]['vmtk_CL'][ch][cont] = {'Status': 'NI'}
                                    dict_MeshesProc[process]['buildCL'][ch][cont] = {'Status': 'NI'}
                    else: 
                        if process == 'A-Create3DMesh':
                            dict_MeshesProc[process][ch] = {'Status': 'NI'}

                for cont in ['tiss', 'int', 'ext']:
                    if (ch, cont, 'whole') in item_thickness_intext:
                        dict_MeshesProc['D-Thickness_int>ext'][ch] = {}
                        dict_MeshesProc['D-Thickness_int>ext'][ch][cont] = {'Status': 'NI'}
                        
                    if (ch, cont, 'whole') in item_thickness_extint:
                         dict_MeshesProc['D-Thickness_ext>int'][ch] = {}
                         dict_MeshesProc['D-Thickness_ext>int'][ch][cont] = {'Status': 'NI'} 

            #Ballooning
            for opt in ball_opts: 
                ch = ball_opts[opt]['to_mesh']
                cont = ball_opts[opt]['to_mesh_type']
                cl_ch = ball_opts[opt]['from_cl']
                cl_cont = ball_opts[opt]['from_cl_type']
            
                if ch in dict_MeshesProc['D-Ballooning'].keys():
                    dict_MeshesProc['D-Ballooning'][ch][cont+'_('+cl_ch+'_'+cl_cont+')'] =  {'Status': 'NI'}
                else: 
                    dict_MeshesProc['D-Ballooning'][ch] = {}
                    dict_MeshesProc['D-Ballooning'][ch][cont+'_('+cl_ch+'_'+cl_cont+')'] =  {'Status': 'NI'}

            # Segments
            for cutg in cut_segm: 
                dict_MeshesProc['E-Segments'][cutg] = {}
                for ch in ch_segm:
                    dict_MeshesProc['E-Segments'][cutg][ch] = {}
                    for cont in ['tiss', 'int', 'ext']:
                        if (cutg, ch, cont, 'segm1') in segm_list:
                            dict_MeshesProc['E-Segments'][cutg][ch][cont] = {'Status': 'NI'}

            # Sections
            for cutc in cut_sect: 
                dict_MeshesProc['E-Sections'][cutc] = {}
                for ch in ch_sect:
                    dict_MeshesProc['E-Sections'][cutc][ch] = {}
                    for cont in ['tiss', 'int', 'ext']:
                        if (cutc, ch, cont, 'sect1') in sect_list:
                            dict_MeshesProc['E-Sections'][cutc][ch][cont] = {'Status': 'NI'}

            # Segments-Sections
            dict_segm_sect = {'Status': 'NI'}
            for cutk in scut_segm_sect: 
                dict_segm_sect[cutk] = {}
                for cutr in rcut_segm_sect: 
                    dict_segm_sect[cutk][cutr] = {}
                    for ch in ch_segm_sect:
                        dict_segm_sect[cutk][cutr][ch] = {}
                        for cont in ['tiss', 'int', 'ext']:
                            if (cutk, cutr, ch, cont, 'segm1', 'sect1') in segm_sect_list:
                                dict_segm_sect[cutk][cutr][ch][cont] = {'Status': 'NI'}

            dict_segm_sectf = flatdict.FlatDict(dict_segm_sect)
            for key in dict_segm_sectf.keys(): 
                if dict_segm_sectf[key] == {}:
                    dict_segm_sectf.pop(key, None)

            dict_MeshesProc['E-Segments_Sections'] = dict_segm_sectf.as_dict()

            # Measure Dictionary
            dict_meas = flatdict.FlatDict(mH_param2meas).as_dict()
            for dicti in dict_meas: 
                if isinstance(dict_meas[dicti], flatdict.FlatDict):
                    dict_meas[dicti] = dict_meas[dicti].as_dict()

            ball_dict = {}
            for keyb in dict_meas['ball']:
                value = dict_meas['ball'][keyb]
                split_key = keyb.split('_')
                ch = split_key[0]
                cont = split_key[1]
                cl_ch = split_key[2].split('(')[1]
                cl_cont = split_key[3].split(')')[0]
                if ch+'_'+cont in ball_dict.keys():
                    ball_dict[ch+'_'+cont][cl_ch+'_'+cl_cont] =  value
                else: 
                    ball_dict[ch+'_'+cont] = {}
                    ball_dict[ch+'_'+cont][cl_ch+'_'+cl_cont] = value
            #Remove ballooning from dict_meas
            dict_meas.pop('ball', None)

            for param in dict_meas: 
                flat_param = flatdict.FlatDict(dict_meas[param], delimiter='_').as_dict()
                dict_meas[param] = flat_param

            #Add balloning dict 
            dict_meas['ball'] = ball_dict
            dict_MeshesProc['F-Measure'] = dict_meas
            dict_MeshesProc['F-Measure']['Status'] = 'NI'
                                                                                       
            dict_mH['ImProc'] = dict_ImProc
            dict_mH['MeshesProc'] = dict_MeshesProc

            #Create new mH_settings['measure'] where to add data
            self.create_df_res()
        
        if self.analysis['morphoCell']:
            mC_channels = sorted(self.mC_channels)

            dict_mC['Status'] = 'NI'
            for met in self.mC_methods: 
                dict_mC[met] = {'Status': 'NI'}

            #Project Status
            for ch in mC_channels: 
                if ch != 'chA':
                    dict_mC['A-SetExtraChs'][ch] = {'Status': 'NI'}

            #Segments
            if 'B-Segments' in self.mC_methods: 
                for cuts in self.mC_settings['setup']['segm_mC'].keys():
                    if 'Cut' in cuts and not 'In2' in cuts:
                        dict_mC['B-Segments'][cuts] = {'Status': 'NI'}

            #Regions
            if 'B-Regions' in self.mC_methods: 
                for cuts in self.mC_settings['setup']['sect_mC'].keys():
                    if 'Cut' in cuts and not 'In2' in cuts:
                        dict_mC['B-Regions'][cuts] = {'Status': 'NI'}
            
            #Zones #check if status: NI are being created for all zones
            if 'B-Zones' in self.mC_methods: 
                for cuts in self.mC_settings['setup']['zone_mC'].keys():
                    if 'Cut' in cuts and not 'In2' in cuts:
                        dict_mC['B-Zones'][cuts] = {'Status': 'NI'}
        
        workflow['morphoHeart'] = dict_mH
        workflow['morphoCell'] = dict_mC
        self.workflow = workflow
    
    def create_df_res(self): 
        mH_params = self.mH_settings['setup']['params']
        for num in mH_params: 
            if mH_params[num]['s'] == 'CL': 
                break
        cl_dict = mH_params[num]['measure']

        measurements = self.mH_settings['measure']
        for param in measurements: 
            if 'CL' in param: 
                for key in measurements[param]: 
                    measurements[param][key] = cl_dict
            if 'Ellip' in param: 
                for key in measurements[param]: 
                    measurements[param][key] = {'ell_width': True, 'ell_length': True, 'ell_depth': True, 'ell_asphericity': True}

        dict_names = {}
        for pp in mH_params:
            var = mH_params[pp]
            dict_names[var['s']] = var['l']
        dict_names['Ellip'] = 'Ellipsoid'
        dict_names['Angles'] = 'Angles'

        df_index = pd.DataFrame.from_dict(measurements, orient='index')
        #Drop variables that don't result in a single measurement
        vars2drop = ['th_e2i', 'th_i2e', 'ball', 'hm3Dto2D']
        vars = list(df_index.index)
        for var in vars: 
            if var in vars2drop: 
                df_index = df_index.drop(var)
        cols = list(df_index.columns)

        #Add column with actual names of variables
        var_names = []
        for index, row in df_index.iterrows(): 
            try: 
                var_names.append(dict_names[index])
            except: 
                var, typpe = index.split('(')
                if typpe == 'segm)': 
                    name = 'Segment'
                elif typpe == 'segm-sect)': 
                    name = 'Segm-Reg'
                else: 
                    name = 'Region'
                var_names.append(dict_names[var]+': '+name)

        df_index['Parameter'] = var_names
        df_index = df_index.reset_index()
        df_index = df_index.drop(['index'], axis=1)
        df_melt = pd.melt(df_index, id_vars = ['Parameter'],  value_vars=cols, value_name='Value')
        df_melt = df_melt.rename(columns={"variable": "Tissue-Contour"})
        df_melt = df_melt.dropna()
        mult_index= ['Parameter', 'Tissue-Contour']
        df_melt = df_melt.set_index(mult_index)
        #Create a copy to modify
        df_final = df_melt.copy(deep=True)

        if 'CL' in vars:
            key_cl = {'linear_length': 'Linear Length', 'looped_length': 'Looped Length'}
            dict_CL = {}
            df_CL = df_melt.loc[[dict_names['CL']]]
            for index, row in df_CL.iterrows():
                row_cl = row['Value']
                df_final.drop(index, axis=0, inplace=True)
                for key, item in row_cl.items():
                    new_index = 'Centreline: '+key_cl[key]
                    new_variable = index[1]
                    dict_CL[(new_index, new_variable)] = item

            if len(dict_CL) != 0: 
                df_CL = pd.DataFrame(dict_CL, index =[0])
                df_CL_melt = pd.melt(df_CL, var_name=mult_index,value_name='Value')
                df_CL_melt = df_CL_melt.set_index(mult_index)
                df_final = pd.concat([df_final, df_CL_melt])
                df_final = df_final.sort_values(by=['Parameter'])
            else: 
                df_final = df_final.sort_values(by=['Parameter'])

        #Add values from Ellipsoids
        if 'Ellip(segm)' in vars: 
            # key_ellip = {'ell_width': 'Width', 'ell_length': 'Length', 'ell_depth': 'Depth', 'ell_asphericity': 'Asphericity'}
            key_ellip = {'ell_width': 'X-dim', 'ell_length': 'Y-dim', 'ell_depth': 'Z-dim', 'ell_asphericity': 'Asphericity'}
            dict_ellip = {}
            df_ellip = df_melt.loc[['Ellipsoid: Segment']]
            for index, row in df_ellip.iterrows():
                row_ell = row['Value']
                df_final.drop(index, axis=0, inplace=True)
                for key, item in row_ell.items():
                    new_index = 'Ellipsoid: '+key_ellip[key]
                    new_variable = index[1]
                    dict_ellip[(new_index, new_variable)] = item

            if len(dict_ellip) != 0: 
                df_ellip = pd.DataFrame(dict_ellip, index =[0])
                df_ellip_melt = pd.melt(df_ellip, var_name=mult_index,value_name='Value')
                df_ellip_melt = df_ellip_melt.set_index(mult_index)
                df_final = pd.concat([df_final, df_ellip_melt])
                df_final = df_final.sort_values(by=['Parameter'])
            else: 
                df_final = df_final.sort_values(by=['Parameter'])

        #Change True Values to TBO
        values_updated = []
        for index, row in df_final.iterrows(): 
            if isinstance(row['Value'], bool): 
                values_updated.append('TBO')
            else: 
                values_updated.append(row['Value'])
        df_final['Value'] = values_updated

        #Add column with better names
        name_chs = self.mH_settings['setup']['name_chs']
        if isinstance(self.mH_settings['setup']['segm'], dict):
            name_segm = {}
            for cut in [key for key in self.mH_settings['setup']['segm'] if 'Cut' in key]:
                name_segm[cut] = self.mH_settings['setup']['segm'][cut]['name_segments']
        if isinstance(self.mH_settings['setup']['sect'], dict):
            name_sect = {}
            for cut in [key for key in self.mH_settings['setup']['sect'] if 'Cut' in key]:
                name_sect[cut] = self.mH_settings['setup']['sect'][cut]['name_sections']
        name_cont = {'int': 'internal', 'tiss': 'tissue', 'ext': 'external'}

        user_tiss_cont = []
        for index, _ in df_final.iterrows(): 
            # print(index)
            param, tiss_cont = index
            split_name = tiss_cont.split('_')
            if len(split_name) == 1 and tiss_cont == 'roi': 
                namef = 'Organ/ROI'
            elif len(split_name) == 3: 
                ch, cont, _ = split_name
                namef = name_chs[ch]+ ' ('+name_cont[cont]+')'
            elif len(split_name) == 4: 
                cut, ch, cont, subm = split_name
                if 'segm' in subm: 
                    namef = cut+': '+name_chs[ch]+ '-'+name_cont[cont]+' ('+name_segm[cut][subm]+')'
                else: 
                    namef = cut+': '+name_chs[ch]+ '-'+name_cont[cont]+' ('+name_sect[cut][subm]+')'
            elif len(split_name) == 6: #Intersections
                scut, rcut, ch, cont, segm, sect = split_name
                namef = scut[1:]+'-'+rcut+': '+name_chs[ch]+ '-'+name_cont[cont]+' ('+name_segm[scut[1:]][segm]+'-'+name_sect[rcut][sect]+')'
            else: 
                namef = 'Check: '+tiss_cont

            user_tiss_cont.append(namef)
        
        df_final['User (Tissue-Contour)'] = user_tiss_cont
        df_finalf = df_final.reset_index()
        df_finalf = df_finalf.set_index(mult_index+['User (Tissue-Contour)'])

        df_finalf = df_finalf.sort_values(['Parameter','Tissue-Contour'],
                                                ascending = [True, True])
        self.mH_settings['df_res'] = df_finalf

    def save_project(self, temp_dir=None, alert_on=True):#
        #Create a new dictionary that contains all the settings
        proj_name = self.user_projName.replace(' ', '_')
        jsonDict_name = 'mH_'+proj_name+'_project.json'
        json2save_par = Path(self.dir_proj) / 'settings'
        json2save_par.mkdir(parents=True, exist_ok=True)
        
        self.dir_info = Path(self.dir_proj) / 'settings' / jsonDict_name
        self.info['dir_info'] = self.dir_info
        
        all_info = {}
        all_info['info'] = self.info#
        all_info['analysis'] = self.analysis#
        all_info['mH_methods'] = self.mH_methods#
        all_info['mH_settings'] = self.mH_settings#
        all_info['mH_channels'] = self.mH_channels#
        all_info['mH_segments'] = self.mH_segments#
        all_info['mH_sections'] = self.mH_sections#
        # all_info['mH_param2meas'] = self.mH_param2meas#
        
        all_info['mC_methods'] = self.mC_methods#
        all_info['mC_settings'] = self.mC_settings#
        all_info['mC_channels'] = self.mC_channels#
        all_info['mC_segments'] = self.mC_segments#
        all_info['mC_sections'] = self.mC_sections#
        # all_info['mC_param2meas'] = self.mC_param2meas#
        
        all_info['workflow'] = self.workflow#
        all_info['organs'] = self.organs#
        all_info['gui_custom_data'] = self.gui_custom_data#
        
        if not json2save_par.is_dir():
            print('>> Error: Settings directory could not be created!\n>> Directory: '+jsonDict_name)
            if alert_on: 
                alert('error_beep')
        else: 
            json2save_dir = json2save_par / jsonDict_name
            with open(str(json2save_dir), "w") as write_file:
                json.dump(all_info, write_file, cls=NumpyArrayEncoder)

            if not json2save_dir.is_file():
                print('>> Error: Project settings file was not saved correctly!\n>> File: '+jsonDict_name)
                if alert_on: 
                    alert('error_beep')
            else: 
                print('>> Project settings file saved correctly! - File: '+jsonDict_name)
                if alert_on: 
                    alert('countdown')

            if temp_dir != None and isinstance(temp_dir, Path): 
                temp_name = temp_dir.stem
                proj_temp = copy.deepcopy(all_info)
                proj_temp['info']['user_projName'] = temp_name
                with open(str(temp_dir), "w") as write_file:
                    json.dump(proj_temp, write_file, cls=NumpyArrayEncoder)
                print('>> Project template file saved correctly!\n>> Path: '+str(temp_dir))
    
    def add_organ(self, organ):#
        organ_name = organ.user_organName
        if organ_name not in self.organs.keys(): 
            #New Organ!
            dict_organ = copy.deepcopy(organ.info)
            dict_organ.pop('project', None)
            dict_organ['dir_res'] = organ.dir_res()
            #Get current workflow
            wf_so_far = self.get_current_wf(organ) #ABC check what exactly this creates
            dict_organ['workflow'] = wf_so_far
            print(' organ wf_so_far:', wf_so_far)
            #Add organ to project's organs
            self.organs[organ_name] = dict_organ
            #Add current organ data to gui data
            strain_it = self.gui_custom_data['strain']
            strain_it.append(organ.info['strain'])
            self.gui_custom_data['strain'] = list(set(strain_it))
            stage_it = self.gui_custom_data['stage']
            stage_it.append(organ.info['stage'])
            self.gui_custom_data['stage'] = list(set(stage_it))
            genot_it = self.gui_custom_data['genotype']
            genot_it.append(organ.info['genotype'])
            self.gui_custom_data['genotype'] =list(set(genot_it))
            manip_it = self.gui_custom_data['manipulation']
            manip_it.append(organ.info['manipulation'])
            self.gui_custom_data['manipulation'] = list(set(manip_it))
            imOr_it = self.gui_custom_data['im_orientation']
            imOr_it.append(organ.info['im_orientation'])
            self.gui_custom_data['im_orientation'] = list(set(imOr_it))
            units_it = self.gui_custom_data['im_res_units']
            units_it.append(organ.info['im_res_units'][0])
            self.gui_custom_data['im_res_units'] = list(set(units_it))
            self.save_project()
        else: 
            #Update organ info
            #Get current workflow
            wf_so_far = self.get_current_wf(organ)
            print('organ wf_so_far:', wf_so_far)
            self.organs[organ_name]['workflow'] = wf_so_far
        
    def delete_organs(self, organs): 
        print('Deleting organs: ', organs)
        for organ in organs: 
            organ_folder = self.organs[organ]['user_organName']
            dir_o = Path(self.dir_proj) / Path(organ_folder)
            dir_f = Path(self.dir_proj) / Path(organ_folder+'_deleted')
            try: 
                rename_directory(dir_o, dir_f)
            except: 
                print('The original directory was not found:', str(dir_o))

            self.organs.pop(organ, None)
        print('self.organs:',self.organs)

    def get_current_wf(self, organ): #

        flat_wf = flatdict.FlatDict(copy.deepcopy(organ.workflow))
        keep_keys = [key for key in flat_wf.keys() if len(key.split(':'))== 4 and 'Status' in key and 'morphoHeart' in key]
        keep_keys_mC = [key for key in flat_wf.keys() if len(key.split(':'))== 3 and 'morphoCell' in key]
        for key in flat_wf.keys(): 
            if key not in keep_keys+keep_keys_mC: 
                flat_wf.pop(key, None)
        out_dict = flat_wf.as_dict()
        for keyi in out_dict: 
            if isinstance(out_dict[keyi], flatdict.FlatDict):
                out_dict[keyi] = out_dict[keyi].as_dict()
        
        out_dict['morphoHeart']['MeshesProc']['F-Measure']['Status'] = organ.measure_status
        print('out_dict: ', out_dict)
        return out_dict

    def remove_organ(self, organ):
        organs = copy.deepcopy(self.organs)
        organ_name = organ.user_organName
        organs.pop(organ_name, None)
        self.organs = organs
        self.save_project()

    def load_organ(self, organ_to_load:str):#

        organ_folder = self.organs[organ_to_load]['user_organName']
        dir_res = Path(self.dir_proj) / organ_folder
        jsonDict_name = 'mH_'+organ_folder+'_organ.json'
        json2open_dir = Path(dir_res) / 'settings' / jsonDict_name
        if json2open_dir.is_file():
            with open(json2open_dir, "r") as read_file:
                print(">> "+jsonDict_name+": Opening JSON encoded data")
                dict_out = json.load(read_file)
            organ_dict = {'load_dict': dict_out}
            organ = Organ(project=self, organ_dict=organ_dict, new=False)
        else: 
            organ = None
            print('>> Error: No organ name with name ',self.organs[organ_to_load]['user_organName'],' was found!\n Directory: ',str(json2open_dir))
            alert('error_beep')
            
        return organ

class Organ():
    'Organ Class'
    
    def __init__(self, project:Project, organ_dict:dict, new:bool):# 
        
        self.parent_project = project
        self.on_hold = False
        
        if new:
            print('\nCreating new organ!')
            self.measure_status = 'NI'
            user_settings = organ_dict['settings']
            img_dirs = organ_dict['img_dirs']
            self.user_organName = user_settings['user_organName']
            self.info = user_settings
            self.info['dirs'] = project.info['dirs']    
            self.img_dirs = img_dirs
            self.create_mHName()
            self.analysis = copy.deepcopy(project.analysis)
            if self.analysis['morphoHeart']:
                self.mH_settings = copy.deepcopy(project.mH_settings)
                self.imChannels = {}
                self.obj_imChannels = {}
                self.imChannelNS = {}
                self.obj_imChannelNS = {}
                self.meshes = {}
                self.submeshes = {}
                self.obj_meshes = {}
                self.obj_subm = {}
                self.objects = {'KSplines': {}, 'Spheres': {}}
                self.obj_temp = {}
                if 'C-Centreline' in self.parent_project.mH_methods:
                    self.objects['KSplines']['cut4cl'] = {'bottom': {}, 'top':{}}
                    self.objects['Spheres']['cut4cl'] = {'bottom': {}, 'top':{}}
                    self.objects['Centreline'] = {}
            if self.analysis['morphoCell']: #ABC 
                self.mC_settings = copy.deepcopy(project.mC_settings)
                self.imChannelsMC = {}
                self.obj_imChannelsMC = {}
                self.cellsMC = {}

            self.workflow = copy.deepcopy(project.workflow) 
            self.create_folders(project.analysis) 
            self.create_channels(project.analysis) #ABC create cells and channels others than those shared wth mH
        else: 
            print('\nLoading organ!')
            load_dict = organ_dict['load_dict']
            self.load_organ(load_dict=load_dict)

    def create_mHName(self): #
        now_str = datetime.now().strftime('%Y%m%d%H%M')
        self.mH_organName = 'mH_Organ-'+now_str

    def create_channels(self, analysis): 
        if analysis['morphoHeart']: 
            channels = [key for key in self.parent_project.mH_channels.keys() if 'NS' not in key]
            for ch_name in channels:
                im_ch = self.create_ch(ch_name=ch_name)

        if analysis['morphoCell']: #ABC
            channelsMC = [key for key in self.parent_project.mC_channels.keys()]
            for ch_nameMC in channelsMC:
                if ch_nameMC != 'chA':
                    im_chMC = self.create_mCch(ch_name=ch_nameMC)
                else:
                    cells_ch = self.create_cells(ch_name=ch_nameMC)

    def load_organ(self, load_dict:dict):#
        load_dict = make_Paths(load_dict)

        self.info = load_dict['Organ']
        self.user_organName = self.info['user_organName']
        self.mH_organName = load_dict['mH_organName']
        
        self.img_dirs = load_dict['img_dirs']
        self.folder = load_dict['folder']
        self.analysis = load_dict['analysis']
        if 'on_hold' in load_dict.keys():
            self.on_hold = load_dict['on_hold']
        else: 
            self.on_hold = False
        
        if self.analysis['morphoHeart']:
            tuple_keys = [['mH_settings','setup','chNS','ch_ext'], 
                            ['mH_settings','setup','chNS','ch_int'],]
        else: 
            tuple_keys = []
        
        for ch in load_dict['imChannels'].keys():
            tuple_keys.append(['imChannels', ch, 'shape'])
            tuple_keys.append(['imChannels', ch, 'shape_s3'])
            for cont in load_dict['imChannels'][ch]['contStack'].keys():
                tuple_keys.append(['imChannels', ch, 'contStack',cont,'shape_s3'])
        
        load_dict = make_tuples(load_dict, tuple_keys)
        
        # Workflow
        self.workflow = load_dict['workflow']

        self.objects = load_dict['objects']
        if self.analysis['morphoHeart']:
            # mH_Settings
            self.mH_settings = load_dict['mH_settings']
            json_str = self.mH_settings['df_res'] 
            self.mH_settings['df_res'] = pd.read_json(json_str, orient='table')
            # imChannels
            self.imChannels = load_dict['imChannels']
            self.load_objImChannels()
            # imChannelNS
            self.imChannelNS = load_dict['imChannelNS']
            self.load_objImChannelNS()
            # meshes
            self.meshes = load_dict['meshes']
            self.load_objMeshes()
            # submeshes
            if 'submeshes' in load_dict.keys():
                submeshes_dict = load_dict['submeshes']
                flat_subm_dict = flatdict.FlatDict(submeshes_dict)
                list_colors = [key.split(':') for key in flat_subm_dict if 'color' in key]
                submeshes_dict_new = make_tuples(submeshes_dict, list_colors)
                self.submeshes = submeshes_dict_new
                self.load_objSubmeshes(submeshes_dict)
            else: 
                self.submeshes = {}
            print('>>>> Loaded submeshes: ', self.submeshes)
            #obj_temp
            if 'obj_temp' in load_dict.keys():
                self.obj_temp = load_dict['obj_temp']
                print('organ.obj_temp:', self.obj_temp)
            
        if self.analysis['morphoCell']:
            # mC_Settings
            self.mC_settings = load_dict['mC_settings']
            try: 
                zz_names = {}
                for zii in self.mC_settings['setup']['zone_mC']['Zone1']['name_zones']: 
                    zz_names[zii] = self.mC_settings['setup']['zone_mC']['Zone1']['name_zones'][zii].strip()
                self.mC_settings['setup']['zone_mC']['Zone1']['name_zones'] = zz_names
                print('Modified mC_settings - name_zones')
            except:
                pass

            if 'measure' not in self.mC_settings: 
                self.mC_settings['measure'] = {} 

            if 'imChannelMC' in load_dict.keys():
                self.imChannelsMC = load_dict['imChannelMC']
                self.load_objImChannelMC()
            if 'cells_MC' in load_dict.keys(): 
                self.cellsMC = load_dict['cells_MC']
                self.load_objCells()

            if 'B-Zones' in self.workflow['morphoCell']: 
                if len(self.workflow['morphoCell']['B-Zones']) == 1: 
                    print('Adding wf to Zones')
                    #Get all zones and add Status NI to all
                    zone_all = [zone for zone in self.mC_settings['setup']['zone_mC'] if '2Zones' not in zone]
                    for zone in zone_all: 
                        self.workflow['morphoCell']['B-Zones'][zone] = {'Status': 'NI'}
        
        if 'orientation' in self.mH_settings['wf_info'].keys():
            self.load_orient_cubes()
        
    def load_orient_cubes(self):
        #Stack
        stack_dict = self.mH_settings['wf_info']['orientation']['stack']
        if 'planar_views' in stack_dict.keys() and 'stack_cube' in stack_dict.keys():
            #Create cube
            pos = stack_dict['stack_cube']['pos']
            side = stack_dict['stack_cube']['side']
            color = stack_dict['stack_cube']['color']
            rotateY = stack_dict['stack_cube']['rotateY']
            orient_cube = vedo.Cube(pos=pos, side=side, c=color)
            orient_cube.linewidth(1).force_opaque()
            if rotateY: 
                cust_angle = stack_dict['stack_cube']['rotateY']
                orient_cube.rotate_y(cust_angle)
            orient_cube.pos(pos)
            orient_cube_clear = orient_cube.clone().alpha(0.5)

            orient_cubef = self.colour_cube(orient_cube = orient_cube, 
                                            axis_dict = stack_dict)
            self.stack_cube = {'cube': orient_cubef, 
                               'clear': orient_cube_clear}

        #ROI
        roi_dict = self.mH_settings['wf_info']['orientation']['roi']
        if 'planar_views' in roi_dict.keys() and 'roi_cube' in roi_dict.keys():
            #Create cube
            pos_r = roi_dict['roi_cube']['pos']
            side_r = roi_dict['roi_cube']['side']
            color_r = roi_dict['roi_cube']['color']

            orient_cube_r = vedo.Cube(pos=pos_r, side=side_r, c=color_r)
            orient_cube_r.linewidth(1).force_opaque()

            rot_x = 0; rot_y = 0; rot_z = 0
            if 'rotate_x' in roi_dict['roi_cube'].keys():
                rot_x = roi_dict['roi_cube']['rotate_x']
                
            if 'rotate_y' in roi_dict['roi_cube'].keys():
                rot_y = roi_dict['roi_cube']['rotate_y']
                
            if 'rotate_z' in roi_dict['roi_cube'].keys():
                rot_z = roi_dict['roi_cube']['rotate_z']
               
            if 'rotate_user' in roi_dict['roi_cube'].keys(): 
                rot_x = rot_x + roi_dict['roi_cube']['rotate_user']['rotX']
                rot_y = rot_y + roi_dict['roi_cube']['rotate_user']['rotY']
                rot_z = rot_z + roi_dict['roi_cube']['rotate_user']['rotZ']
            
            orient_cube_r.rotate_x(rot_x)
            orient_cube_r.rotate_y(rot_y)
            orient_cube_r.rotate_z(rot_z)

            orient_cube_r.pos(pos_r)
            orient_cube_clear_r = orient_cube_r.clone().alpha(0.5)

            orient_cube_rf = self.colour_cube(orient_cube = orient_cube_r, 
                                                axis_dict = roi_dict)
            self.roi_cube = {'cube': orient_cube_rf, 
                               'clear': orient_cube_clear_r}

    def colour_cube(self, orient_cube, axis_dict): 
        normal_dict = {}
        for pt in orient_cube.cell_centers(): 
            idcell = orient_cube.closest_point(pt, return_cell_id=True)
            cells = orient_cube.cells()[idcell]
            points = [orient_cube.points()[cell] for cell in cells]
            plane_fit = vedo.fit_plane(points, signed=True)
            normal_dict[idcell] = plane_fit.normal

        for view in axis_dict['planar_views'].keys(): 
            color = axis_dict['planar_views'][view]['color']
            idcell_o = axis_dict['planar_views'][view]['idcell']
            pl_normal = axis_dict['planar_views'][view]['pl_normal']

            if np.array_equal(normal_dict[int(idcell_o)], pl_normal):
                orient_cube.cellcolors[int(idcell_o)] = color
            else: 
                print("idcells don't match!")
                alert('bubble')
        
        return orient_cube
        
    def load_objImChannels(self):
        self.obj_imChannels = {}
        if len(self.imChannels) > 0:
            for imCh in self.imChannels:
                im_ch = ImChannel(organ=self, ch_name=imCh)
                self.obj_imChannels[imCh] = im_ch
        
    def load_objImChannelNS(self):
        self.obj_imChannelNS = {}
        if len(self.imChannelNS) > 0: 
            for imCh in self.imChannelNS:
                im_ch = ImChannelNS(organ=self, ch_name=imCh)#, new=False)
                self.obj_imChannelNS[imCh] = im_ch

            #Modify XOR if organ was created before 26/10/23
            date_update = '2023-10-26'
            date_format = '%Y-%m-%d'
            date_update_f = datetime.strptime(date_update, date_format)
            #Get date of organ creation 
            date_created = self.parent_project.info['date_created']
            date_created_f = datetime.strptime(date_created, date_format)
            if date_created_f < date_update_f:
                self.mH_settings['setup']['chNS']['operation'] = 'AND-XOR'
                self.parent_project.mH_settings['setup']['chNS']['operation'] = 'AND-XOR'
                print('Updated XOR to AND-XOR in organ and project')
        
    def load_objMeshes(self):#
        self.obj_meshes = {}
        if len(self.meshes) > 0: 
            for mesh in self.meshes:
                ch_no = self.meshes[mesh]['channel_no']
                if 'NS' in mesh: 
                    imCh = self.obj_imChannelNS[ch_no]
                else: 
                    imCh = self.obj_imChannels[ch_no]
                mesh_type = self.meshes[mesh]['mesh_type']
                keep_largest = self.meshes[mesh]['keep_largest']
                rotateZ_90 = self.meshes[mesh]['rotateZ_90'] 
                mesh_prop = {'keep_largest': keep_largest, 'rotateZ_90': rotateZ_90}
                msh = Mesh_mH(imChannel = imCh, mesh_type = mesh_type, 
                              mesh_prop = mesh_prop)
                self.obj_meshes[mesh] = msh
    
    def load_objSubmeshes(self, submeshes_dict): 
        flat_subm_dict = flatdict.FlatDict(submeshes_dict)
        list_colors = [key.split(':') for key in flat_subm_dict if 'color' in key]
        submeshes_dict = make_tuples(submeshes_dict, list_colors)

        self.obj_subm = {}
        for subm in submeshes_dict.keys():
            if 'sCut' not in subm: 
                cut, ch, cont, sub = subm.split('_')
                mesh = self.obj_meshes[ch+'_'+cont]
                if 'segm' in sub: 
                    #Get mesh to create submesh
                    color = submeshes_dict[subm]['color']
                    submesh = mesh.create_segment(name = sub, cut = cut, color =color)
                    self.obj_subm[subm] = submesh
                else: #'sect' in sub: 
                    color = submeshes_dict[subm]['color']
                    submesh = mesh.create_section(name = sub, cut = cut, color =color)
                    self.obj_subm[subm] = submesh
            else: 
                seg_cut, reg_cut, ch, cont, segm, sect = subm.split('_')
                obj_segm = self.obj_subm[seg_cut[1:]+'_'+ch+'_'+cont+'_'+segm]
                #Create subsegm_subsect
                segm_sect = segm+'_'+sect
                cutsf = seg_cut+'_'+reg_cut
                color = submeshes_dict[subm]['color']
                submesh = obj_segm.create_segm_sect(segm_sect = segm_sect, cuts = cutsf, color = color)
                self.obj_subm[subm] = submesh

    def load_objImChannelMC(self): 
        self.obj_imChannelsMC = {}
        if len(self.imChannelsMC) > 0:
            for imCh in self.imChannelsMC:
                im_ch = ImChannelMC(organ=self, ch_name=imCh)
                self.obj_imChannelsMC[imCh] = im_ch

    def load_objCells(self): 
        
        cells = Cells(organ=self)
        self.cellsMC = {'chA': cells}

    def create_folders(self, analysis):#
        if analysis['morphoHeart']: #ABC
            dirResults = ['meshes', 'csv_all', 'imgs_videos', 's3_numpy', 'centreline', 'settings']
        else:
            if analysis['morphoCell']:
                dirResults = ['meshes', 'csv_all', 'imgs_videos', 's3_numpy', 'settings']
            else:
                dirResults = []

        organ_folder = self.user_organName
        for direc in dirResults:
            dir2create = Path(self.parent_project.dir_proj) / organ_folder / direc
            dir2create.mkdir(parents=True, exist_ok=True)
            if dir2create.is_dir():
                self.info['dirs'][direc] = True #dir2create
            else: 
                print('>> Error: Directory ', organ_folder, '/', direc, ' could not be created!')
                alert('error_beep')
        self.folder = organ_folder
        # self.dir_res = self.parent_project.dir_proj / organ_folder

    def dir_res(self, dir=None):# 
        dir_proj = Path(self.parent_project.dir_proj)
        dir_res = dir_proj / self.folder
        if dir == None: 
            return dir_res
        else: 
            sp_dir = self.info['dirs'][dir]
            if sp_dir: 
                dir_sp = dir_res / dir
                return Path(dir_sp)
            else: 
                print("The specified directory doesn't exist: (dir_res: "+ str(dir_res)+", dir: "+dir)
                return None

    def add_channel(self, imChannel):#
        # Check first if the channel has been already added to the organ
        new = False
        if imChannel.channel_no not in self.imChannels.keys():
            new = True
            
        if new: 
            print('>> Adding Channel - ', imChannel.channel_no)
            channel_dict = {}
            channel_dict['parent_organ_name'] = imChannel.parent_organ_name
            channel_dict['channel_no'] = imChannel.channel_no
            channel_dict['user_chName'] = imChannel.user_chName
            channel_dict['ch_relation'] = imChannel.ch_relation
            channel_dict['to_mask'] = imChannel.to_mask
            channel_dict['resolution'] = imChannel.resolution
            channel_dict['dir_cho'] = imChannel.dir_cho
            if imChannel.to_mask:
                channel_dict['dir_mk'] = imChannel.dir_mk
            channel_dict['masked'] = imChannel.masked
            channel_dict['shape'] = imChannel.shape
            channel_dict['process'] = imChannel.process
            channel_dict['contStack'] = imChannel.contStack
            channel_dict['dir_stckproc'] = imChannel.dir_stckproc
            
            self.imChannels[imChannel.channel_no] = channel_dict
            
        else: # just update im_proc 
            self.imChannels[imChannel.channel_no]['process'] = imChannel.process
            self.imChannels[imChannel.channel_no]['contStack'] = imChannel.contStack
            dir_stck = self.dir_res(dir='s3_numpy') / imChannel.dir_stckproc            
            if dir_stck.is_file():
                self.imChannels[imChannel.channel_no]['dir_stckproc'] = imChannel.dir_stckproc
            if hasattr(imChannel, 'shape_s3'):
                self.imChannels[imChannel.channel_no]['shape_s3'] = imChannel.shape_s3
                self.info['shape_s3'] = imChannel.shape_s3

        self.obj_imChannels[imChannel.channel_no] = imChannel
        
    def add_channelNS(self, imChannelNS):
        # Check first if the channel has been already added to the organ
        new = False
        if imChannelNS.channel_no not in self.imChannelNS.keys():
            new = True
            
        if new: 
            print('>> Adding ChannelNS - ', imChannelNS.channel_no)
            channel_dict = {}
            channel_dict['parent_organ_name'] = imChannelNS.parent_organ_name
            channel_dict['channel_no'] = imChannelNS.channel_no
            channel_dict['user_chName'] = imChannelNS.user_chName
            channel_dict['ch_relation'] = imChannelNS.ch_relation
            channel_dict['resolution'] = imChannelNS.resolution
            channel_dict['process'] = imChannelNS.process
            channel_dict['contStack'] = imChannelNS.contStack
            channel_dict['setup_NS'] = imChannelNS.setup_NS
            self.imChannelNS[imChannelNS.channel_no] = channel_dict
            
        else: # just update im_proc 
            self.imChannelNS[imChannelNS.channel_no]['process'] = imChannelNS.process
            self.imChannelNS[imChannelNS.channel_no]['contStack'] = imChannelNS.contStack
           
        self.obj_imChannelNS[imChannelNS.channel_no] = imChannelNS

    def add_mesh(self, mesh):#
    
        new = False
        if mesh.name not in self.meshes.keys():
            new = True
        if new: 
            print('>> Adding Mesh - ', mesh.name)
            self.meshes[mesh.name] = {}
            self.meshes[mesh.name]['parent_organ'] = mesh.parent_organ.user_organName
            self.meshes[mesh.name]['channel_no'] = mesh.imChannel.channel_no
            self.meshes[mesh.name]['user_meshName'] = mesh.user_meshName
            self.meshes[mesh.name]['mesh_type'] = mesh.mesh_type
            self.meshes[mesh.name]['legend'] = mesh.legend
            self.meshes[mesh.name]['name'] = mesh.name
            self.meshes[mesh.name]['resolution'] = mesh.resolution
            self.meshes[mesh.name]['color'] = mesh.color
            self.meshes[mesh.name]['alpha'] = mesh.alpha
            self.meshes[mesh.name]['keep_largest'] = mesh.keep_largest
            self.meshes[mesh.name]['rotateZ_90'] = mesh.rotateZ_90
            self.meshes[mesh.name]['s3_dir'] = mesh.s3_dir
            if hasattr(mesh,'dir_out'):
                self.meshes[mesh.name]['dir_out'] = mesh.dir_out
            self.obj_meshes[mesh.name] = mesh
            if hasattr(mesh, 'dirs'):
                self.meshes[mesh.name]['dirs'] = mesh.dirs
        else: #Just updating things that could change
            self.meshes[mesh.name]['color'] = mesh.color
            self.meshes[mesh.name]['alpha'] = mesh.alpha
            self.meshes[mesh.name]['keep_largest'] = mesh.keep_largest
            self.obj_meshes[mesh.name] = mesh
            if hasattr(mesh, 'dirs'):
                self.meshes[mesh.name]['dirs'] = mesh.dirs
            print('>> Mesh data updated!')
    
    def add_submesh(self, submesh):
        new = False
        if submesh.sub_name_all not in self.submeshes.keys():
            new = True
        if new:
            print('>> Adding SubMesh - ', submesh.sub_name_all)
            self.submeshes[submesh.sub_name_all] = {}
            #New to submesh
            self.submeshes[submesh.sub_name_all]['sub_name'] = submesh.sub_name
            self.submeshes[submesh.sub_name_all]['sub_name_all'] = submesh.sub_name_all
            if submesh.sub_mesh_type != 'Segment-Section': 
                self.submeshes[submesh.sub_name_all]['parent_mesh'] = submesh.parent_mesh.name
            else: 
                self.submeshes[submesh.sub_name_all]['parent_mesh'] = submesh.parent_mesh.parent_mesh.name

            self.submeshes[submesh.sub_name_all]['sub_mesh_type'] = submesh.sub_mesh_type
            self.submeshes[submesh.sub_name_all]['sub_legend'] = submesh.sub_legend
            self.submeshes[submesh.sub_name_all]['sub_user_name'] = submesh.sub_user_name
            self.submeshes[submesh.sub_name_all]['cut'] = submesh.cut
            #Mesh related
            self.submeshes[submesh.sub_name_all]['resolution'] = submesh.resolution
            self.submeshes[submesh.sub_name_all]['color'] = submesh.color
            self.submeshes[submesh.sub_name_all]['alpha'] = submesh.alpha
            self.submeshes[submesh.sub_name_all]['keep_largest'] = submesh.keep_largest
            self.submeshes[submesh.sub_name_all]['rotateZ_90'] = submesh.rotateZ_90
            #Inherited from parent_mesh
            if submesh.sub_mesh_type != 'Segment-Section': 
                self.submeshes[submesh.sub_name_all]['parent_mesh'] = {}
                self.submeshes[submesh.sub_name_all]['parent_mesh']['legend'] = submesh.parent_mesh.legend
                self.submeshes[submesh.sub_name_all]['parent_mesh']['name'] = submesh.parent_mesh.name
                self.submeshes[submesh.sub_name_all]['parent_mesh']['imChannel'] = submesh.parent_mesh.imChannel.channel_no
            else: 
                self.submeshes[submesh.sub_name_all]['parent_mesh'] = {}
                self.submeshes[submesh.sub_name_all]['parent_mesh']['legend'] = submesh.parent_mesh.parent_mesh.legend
                self.submeshes[submesh.sub_name_all]['parent_mesh']['name'] = submesh.parent_mesh.parent_mesh.name
                self.submeshes[submesh.sub_name_all]['parent_mesh']['imChannel'] = submesh.parent_mesh.parent_mesh.imChannel.channel_no
            #Add to obj_submesh
            self.obj_subm[submesh.sub_name_all] = submesh
        else: #Just updating things that could change
            self.submeshes[submesh.sub_name_all]['color'] = submesh.color
            self.submeshes[submesh.sub_name_all]['alpha'] = submesh.alpha
            self.submeshes[submesh.sub_name_all]['keep_largest'] = submesh.keep_largest
            print('>> SubMesh data updated!')
            #Add to obj_submesh
            self.obj_subm[submesh.sub_name_all] = submesh

        if submesh.sub_mesh_type == 'Section' or submesh.sub_mesh_type == 'Segment-Section':
            self.submeshes[submesh.sub_name_all]['s3_invert'] = submesh.s3_invert
        elif submesh.sub_mesh_type == 'Segment':
            if hasattr(submesh, 'dict_segm'):
                self.submeshes[submesh.sub_name_all]['dict_segm'] = submesh.dict_segm
       
    def add_object(self, obj, proc:str, class_name:Union[list,str], name):
        
        if isinstance(obj, vedo.shapes.KSpline):# or name == 'KSpline':
            if proc != 'Centreline': 
                if isinstance(class_name, list):
                    classif, mesh_name = class_name
                    self.objects['KSplines'][proc][classif][mesh_name] = {'points': obj.points(), 
                                                                       'color': obj.color()}
                else: 
                    self.objects['KSplines'][proc][class_name] = {'points': obj.points(), 
                                                                       'color': obj.color()}
            else: 
                self.objects[proc][class_name] = {'points': obj.points(), 
                                                  'color': obj.color()}
            
        if isinstance(obj, vedo.shapes.Sphere):# or name == 'Sphere':
            if isinstance(class_name, list):
                classif, mesh_name = class_name
                self.objects['Spheres'][proc][classif][mesh_name] = {'center': obj.center, 
                                                                   'color': obj.color()}
            else: 
                self.objects['Spheres'][proc][class_name] = {'center': obj.center, 
                                                                   'color': obj.color()}

    def add_channel_mC(self, imChannelMC): 
        # Check first if the channel has been already added to the organ
        new = False
        if imChannelMC.channel_no not in self.imChannelsMC.keys():
            new = True

        if new: 
            print('>> Adding morphoCell extra Channel - ', imChannelMC.channel_no)
            channel_dict = {}
            channel_dict['parent_organ_name'] = imChannelMC.parent_organ_name
            channel_dict['channel_no'] = imChannelMC.channel_no
            channel_dict['user_chName'] = imChannelMC.user_chName
            channel_dict['mH_channel'] = imChannelMC.mH_channel
            channel_dict['resolution'] = imChannelMC.resolution
            channel_dict['dir_cho'] = imChannelMC.dir_cho
            channel_dict['shape'] = imChannelMC.shape
            channel_dict['process'] = imChannelMC.process

            self.imChannelsMC[imChannelMC.channel_no] = channel_dict

        self.obj_imChannelsMC[imChannelMC.channel_no] = imChannelMC

    def add_cells_mC(self, cells): 
        new = False
        if 'chA' not in self.cellsMC.keys(): 
            new = True
        
        if new: 
            print('>> Adding morphoCell Cells - ', cells.channel_no)
            cells_dict = {}
            cells_dict['parent_organ_name'] = cells.parent_organ_name
            cells_dict['channel_no'] = cells.channel_no
            cells_dict['user_chName'] = cells.user_chName
            cells_dict['resolution'] = cells.resolution
            cells_dict['dir_cho'] = cells.dir_cho
            cells_dict['dir_img'] = cells.dir_img
            cells_dict['dir_cells'] = cells.dir_cells
            cells_dict['shape'] = cells.shape

            self.cellsMC[cells.channel_no] = cells_dict

    def load_objTemp(self, proc, key, ch_cont, obj_temp): 

        if proc == 'centreline':
            if 'SimplifyMesh' == key:
                ch, cont = ch_cont.split('_')
                #Mesh
                directory = self.dir_res(dir ='centreline')
                mesh_dir = directory / self.mH_settings['wf_info']['centreline']['dirs'][ch][cont]['dir_cleanMesh']
                mesh_out = vedo.load(str(mesh_dir))
                mesh_out.color(self.mH_settings['setup']['color_chs'][ch][cont])
                obj_temp[proc][key][ch_cont]['mesh'] = mesh_out
                for side in ['bottom', 'top']: 
                    #Kspl
                    points = self.objects['KSplines']['cut4cl'][side][ch_cont]['points']
                    color = self.objects['KSplines']['cut4cl'][side][ch_cont]['color']
                    kspl = vedo.KSpline(points, continuity=0, tension=0, bias=0, closed=True)
                    kspl.color(color)
                    obj_temp[proc][key][ch_cont]['kspl'][side] = kspl
                     #Centroid
                    center = self.objects['Spheres']['cut4cl'][side][ch_cont]['center']
                    color = self.objects['Spheres']['cut4cl'][side][ch_cont]['color']
                    sph_centroid = vedo.Sphere(pos=center, r=2).legend('cut4CL_'+side)
                    sph_centroid.color(color)
                    obj_temp[proc][key][ch_cont]['centroid'][side] = sph_centroid
        elif proc == 'segments': 
            #key = cut, ch_cont = submesh.sub_name_all
            cut, ch, cont, segm = ch_cont.split('_')
            directory = self.dir_res(dir = 'meshes')
            mesh_dir = directory / self.mH_settings['wf_info']['segments']['setup'][key]['dirs'][ch][cont][segm]
            mesh_out = vedo.load(str(mesh_dir))
            mesh_out.color(self.mH_settings['setup']['segm'][cut]['colors'][segm])
            obj_temp[proc][key][ch_cont] = mesh_out
              
        return obj_temp

    def create_ch(self, ch_name:str):#
        print('---- Creating channel ('+ch_name+')! ----')
        image = ImChannel(organ=self, ch_name=ch_name)#,new=True
        return image
    
    def create_mCch(self, ch_name:str):
        print('---- Creating morphoCell channel ('+ch_name+')! ----')
        imageMC = ImChannelMC(organ=self, ch_name=ch_name)#,new=True
        return imageMC

    def create_cells(self, ch_name:str):
        print('---- Creating Cells for morphoCell analysis ('+ch_name+')! ----')
        cells = Cells(organ=self)#,new=True
        return cells

    def save_organ(self, alert_on=True):#
        print('Saving organ')
        print('self.obj_temp:', self.obj_temp)
        all_info = {}
        all_info['Organ'] = self.info
        all_info['img_dirs'] = self.img_dirs
        all_info['folder'] = self.folder
        all_info['analysis'] = self.analysis
        all_info['on_hold'] = self.on_hold
        
        organ_name = self.user_organName#.replace(' ', '_')
        jsonDict_name = 'mH_'+organ_name+'_organ.json'
        json2save_dir = self.dir_res(dir='settings') / jsonDict_name

        if self.analysis['morphoHeart']:
            all_info['mH_settings'] = self.mH_settings
            
            image_dict = copy.deepcopy(self.imChannels)
            for ch in image_dict.keys():
                image_dict[ch].pop('parent_organ', None)
            all_info['imChannels'] = image_dict
            
            imageNS_dict = copy.deepcopy(self.imChannelNS)
            for chNS in imageNS_dict.keys():
                imageNS_dict[chNS].pop('parent_organ', None)
            all_info['imChannelNS'] = imageNS_dict

            all_info['meshes'] = self.meshes
            all_info['submeshes'] = self.submeshes
            all_info['objects'] = self.objects

            #Flatten obj_temp
            flat_obj_temp = flatdict.FlatDict(self.obj_temp)
            for key in flat_obj_temp: 
                flat_obj_temp[key] = None
            obj_temp = flat_obj_temp.as_dict()
            all_info['obj_temp'] = obj_temp
            
        if self.analysis['morphoCell']:
            #ABC what else is there to save from MC 
            all_info['mC_settings'] = self.mC_settings

            image_dictMC = copy.deepcopy(self.imChannelsMC)
            for chMC in image_dictMC.keys():
                image_dictMC[chMC].pop('parent_organ', None)
            all_info['imChannelMC'] = image_dictMC

            try: 
                cells_MC = {'chA': {'parent_organ_name': self.cellsMC['chA']['parent_organ_name'], 
                                    'channel_no': self.cellsMC['chA']['channel_no'], 
                                    'user_chName': self.cellsMC['chA']['user_chName'],
                                    'dir_cells': self.cellsMC['chA']['dir_cells'],
                                    'dir_cho': self.cellsMC['chA']['dir_cho'],
                                    'dir_img': self.cellsMC['chA']['dir_img'],
                                    'resolution': self.cellsMC['chA']['resolution'], 
                                    'shape': self.cellsMC['chA']['shape']}}
                print('try saving cells')
            except:
                print('except saving cells')
                cells_MC = {'chA': {'parent_organ_name': self.cellsMC['chA'].parent_organ_name,
                                'channel_no': self.cellsMC['chA'].channel_no, 
                                'user_chName': self.cellsMC['chA'].user_chName,
                                'dir_cells': self.cellsMC['chA'].dir_cells,
                                'dir_cho': self.cellsMC['chA'].dir_cho,
                                'dir_img': self.cellsMC['chA'].dir_img,
                                'resolution': self.cellsMC['chA'].resolution, 
                                'shape': self.cellsMC['chA'].shape}}

            all_info['cells_MC'] = cells_MC

        all_info['workflow'] = self.workflow

        self.dir_info = Path('settings') / jsonDict_name
        all_info['dir_info'] = self.dir_info
        all_info['mH_organName'] = self.mH_organName

        with open(str(json2save_dir), "w") as write_file:
            json.dump(all_info, write_file, cls=NumpyArrayEncoder)

        if not json2save_dir.is_file():
            print('>> Error: Organ settings file was not saved correctly!\n>> File: '+jsonDict_name)
            if alert_on: 
                alert('error_beep')
        else: 
            print('\n>> Organ settings file saved correctly! - File: '+jsonDict_name)
            if alert_on:
                alert('countdown')
            
    def check_status(self, process:str):
        wf = self.workflow['morphoHeart']
        if process=='ImProc':
            ch_done = []
            for ch in self.imChannels.keys():
                # First check close contours
                close_done = []
                for key_a in ['A-Autom', 'B-Manual', 'C-CloseInOut']:
                    val = get_by_path(wf, [process,ch,'B-CloseCont','Steps',key_a,'Status'])
                    close_done.append(val)
                if all(flag == 'DONE' for flag in close_done):
                    self.update_mHworkflow([process,ch,'B-CloseCont','Status'], 'DONE')

                # Now update all the workflow
                proc_done = []
                for key_b in wf[process][ch].keys():
                    if key_b != 'Status':
                        val_b = get_by_path(wf, [process,ch,key_b,'Status'])
                        proc_done.append(val_b)
                if all('DONE' in flag for flag in proc_done):
                    self.update_mHworkflow([process,ch,'Status'], 'DONE')
                val_c = get_by_path(wf, [process,ch,'Status'])
                ch_done.append(val_c)
            
            if all(flag == 'DONE' for flag in ch_done):
                self.update_mHworkflow([process,'Status'], 'DONE')
                
        if process == 'MeshesProc':
            proc_done = []
            flat_dict = flatdict.FlatDict(copy.deepcopy(wf[process]))
            for key_f in [proc for proc in wf[process].keys() if proc != 'Status']:
                proc_done.append(wf[process][key_f]['Status'])
                if key_f != 'E-Segments':
                    dict_proc = [item for item in flat_dict.keys() if key_f in item]
                    ch_cont_done = [flat_dict[item] for item in dict_proc[1:]]
                    if all(flag == 'DONE' for flag in ch_cont_done):
                        self.update_mHworkflow([process,key_f,'Status'], 'DONE')
            
            if all(flag == 'DONE' for flag in proc_done):
                self.update_mHworkflow([process,'Status'], 'DONE')
            elif any(flag == 'DONE' for flag in proc_done):
                self.update_mHworkflow([process,'Status'], 'Initialised')
            else: 
                pass

    def update_mHworkflow(self, process, update):#
        workflow = self.workflow['morphoHeart']
        set_by_path(workflow, process, update)
        print('> Update:', process, get_by_path(workflow, process))   

    def update_mCworkflow(self, process, update):#
        workflow = self.workflow['morphoCell']
        set_by_path(workflow, process, update)
        print('> Update:', process, get_by_path(workflow, process))       

    def update_settings(self, process, update, mH='mH', add=None):#
        if mH =='mH':
            settings = self.mH_settings
        else:  #=='mC'
            settings = self.mC_settings
        
        if process == ['wf_info'] and add != None: 
            if len(settings[process[0]])>0: 
                set_keys = list(settings[process[0]].keys())
                if add not in set_keys: 
                    settings_proc = {}
                    for key in set_keys: 
                        settings_proc[key] = settings[process[0]][key]
                    settings_proc[add] = {}
                    settingsf = copy.deepcopy(settings)
                    settingsf[process[0]] = settings_proc
                else: 
                    settingsf = settings
            else: 
                settings_proc = {}
                settings_proc[add] = {}
                settingsf = copy.deepcopy(settings)
                settingsf[process[0]] = settings_proc
            process = process+[add]
        else: 
            settingsf = settings

        if mH =='mH':
            self.mH_settings = settingsf
        else:  #=='mC'
            self.mC_settings = settingsf

        set_by_path(settingsf, process, update)
    
    def get_ext_int_chs(self, return_int=False): 
        chs = list(self.imChannels.keys())
        ch_ext = []; ch_int = []
        if self.mH_settings['setup']['all_contained'] or self.mH_settings['setup']['one_contained']:
            if len(chs)>1:# and len(chs)<3:
                for ch in chs:
                    if self.mH_settings['setup']['chs_relation'][ch] == 'external':
                        ch_ext = self.obj_imChannels[ch]
                    if self.mH_settings['setup']['chs_relation'][ch] == 'internal':
                        ch_int = self.obj_imChannels[ch]
            else: #if len(chs) == 1: 
                print('There is only one channel and it is:', chs)
                for ch in chs:
                    ch_ext = self.obj_imChannels[ch]
                ch_int = None
        else: 
            ch_ext = 'independent'
            return_int = False

        if return_int: 
            return ch_ext, ch_int
        else: 
            return ch_ext
    
    def get_ext_subsgm(self, cut): 

        ext_subsgm = {}
        for name in self.mH_settings['wf_info']['segments']['setup'][cut]['names'].items():
            ext_subsgm[name[0]] = self.obj_subm[name[1]]
        self.ext_subsgm = ext_subsgm

        return self.ext_subsgm
    
    def dict_segments(self, cut, palette = None, other=True):
    
        segments_info = self.mH_settings['setup']['segm'][cut]
        # no_segm = segments_info['no_segments']
        #https://seaborn.pydata.org/tutorial/color_palettes.html
        # palette_all = sns.color_palette("cool_r", no_segm*4)
        # palette = random.sample(palette_all, no_segm)
        
        if other: 
            dict_segm = {'other': {'user_name': 'other',
                                'meshes_number': []}}
        else:
            dict_segm = {}
        
        colors = {'other': [128,128,128]}
        for n, segm in enumerate(segments_info['name_segments']): 
            dict_segm[segm] = {}
            dict_segm[segm]['user_name'] = segments_info['name_segments'][segm]
            dict_segm[segm]['meshes_number'] = []
            if palette != None: 
                try: 
                    colors[segm] = palette[segm]
                    print('try colors - dict_segments')
                except:
                    print('except colors - dict_segments')
                    colors[segm] = palette[n]
        
        if palette != None: 
            return dict_segm, colors
        else: 
            return dict_segm

    def check_method(self, method:str):
        if method in self.parent_project.mH_methods:
            return True
        else:
            return False  
    
    def get_orientation(self, views, colors, mtype:str):#CHECkk
       
        im_orient = self.info['im_orientation']
        print('self.info[im_orientation]',im_orient)
        rotateY = False
        if str(self.info['custom_angle']) != '0': 
            cust_angle = self.info['custom_angle']
            rotateY = True

        ext_ch = self.get_ext_int_chs()
        if isinstance(ext_ch, str) and ext_ch == 'independent':
            ext_ch = self.obj_imChannels[list(self.obj_imChannels.keys())[0]]

        mesh_ext = self.obj_meshes[ext_ch.channel_no+'_tiss']
        pos = mesh_ext.mesh.center_of_mass()
        # print('pos:', pos, type(pos))
        side = max(self.get_maj_bounds())
        color_o = [152,251,152,255]
        orient_cube = vedo.Cube(pos=pos, side=side, c=color_o[:-1])
        orient_cube.linewidth(1).force_opaque()

        stack_cube = {'pos': pos, 
                        'side': side, 
                        'color': color_o[:-1], 
                        'rotateY': rotateY}
        
        if rotateY: 
            orient_cube.rotate_y(-cust_angle)
            stack_cube['custom_angle'] = -cust_angle
        else: 
            stack_cube['custom_angle'] = 0
        
        orient_cube.pos(pos)
        orient_cube_clear = orient_cube.clone().alpha(0.5)
        txt0 = vedo.Text2D(self.user_organName+' - Reference cube and mesh to select planar views in '+mtype+'...', c=txt_color, font=txt_font, s=txt_size)
        
        mks = []; sym = ['o']*len(views)
        for n, view, col in zip(count(), views, colors):
            mks.append(vedo.Marker('*').c(col[0:-1]).legend(view))
        lb = vedo.LegendBox(mks, markers=sym, font=txt_font, 
                            width=leg_width/1.5, height=leg_height/1.5)
        
        vpt = MyFaceSelectingPlotter(N=2, axes=1,colors=colors, color_o=color_o, 
                                        views=views)
        vpt.add_icon(logo, pos=(0.1,1), size=0.25)
        vpt.add_callback("key press", vpt.on_key_press)
        vpt.add_callback("mouse click", vpt.select_cube_face)
        vpt.show(mesh_ext.mesh, orient_cube_clear,txt0, at=0)
        vpt.show(orient_cube, lb, vpt.msg, vpt.msg_face, at=1, azimuth=45, elevation=30, zoom=0.8, interactive=True)
            
        self.stack_cube = {'cube': orient_cube,
                           'clear': orient_cube_clear}

        return vpt.planar_views, stack_cube
            
    def get_ROI_orientation(self, gui_orientation, colors):

        views = self.mH_settings['setup']['orientation']['roi'].split(', ')
        if gui_orientation['roi']['reorient']: 
            if gui_orientation['roi']['method'] == 'Centreline': 
                centreline = gui_orientation['roi']['centreline'].split('(')[1].split(')')[0]
                plane = gui_orientation['roi']['plane_orient']
                ref_vect = gui_orientation['roi']['vector_orient']
                planar_views, settings, roi_cube = self.orient_by_cl(views, centreline, plane, ref_vect, colors)
            elif gui_orientation['roi']['method'] == 'Manual': 
                planar_views, settings, roi_cube = self.orient_manual(views, colors)
        else: #No rotation
            settings = None
            planar_views, roi_cube = self.get_orientation(views, colors, mtype='ROI')
        planar_views = self.get_ref_vectors(planar_views, roi_cube)

        return planar_views, settings, roi_cube

    def orient_by_cl(self, views, centreline:str, plane:str, ref_vect:str, colors):
        
        ch, cont = centreline.split('_')
        cl_mesh = self.obj_meshes[ch+'_'+cont]
        linLine = cl_mesh.get_linLine()
        pts = linLine.points()
        
        plane_coord = {'XY': 2, 'YZ': 0, 'XZ': 1}
        if isinstance(ref_vect, str):
            ref_vectAll = {'X+': np.array([[0,0,0],[1,0,0]]),
                           'Y+': np.array([[0,0,0],[0,1,0]]),
                           'Z+': np.array([[0,0,0],[0,0,1]])}
            ref_vectF = ref_vectAll[ref_vect]
            
        coord = plane_coord[plane]
        for pt in pts: 
            pt[coord] = 0
        
        angle = find_angle_btw_pts(pts, ref_vectF)
        if angle > 90: 
            angle = angle - 90
        
        pos = cl_mesh.mesh.center_of_mass()
        side = max(self.get_maj_bounds())  
        color_o = [152,251,152,255]
        orient_cube = vedo.Cube(pos=pos, side=side, c=color_o[:-1])
        orient_cube.linewidth(1).force_opaque()
        
        roi_cube = {'pos': pos, 
                    'side': side, 
                    'color': color_o[:-1]}
        
        if angle != 0: 
            if coord == 0: 
                orient_cube.rotate_x(angle)
                roi_cube['rotate_x'] = angle
            elif coord == 1: 
                orient_cube.rotate_y(angle)
                roi_cube['rotate_y'] = angle
            elif coord == 2: 
                orient_cube.rotate_z(angle)
                roi_cube['rotate_z'] = angle
        
        orient_cube.pos(pos)
        orient_cube_clear = orient_cube.clone().alpha(0.5)
        
        txt0 = vedo.Text2D(self.user_organName+' - Reference cube and mesh to select planar views in ROI (organ)...', c=txt_color, font=txt_font, s=txt_size)
        
        mks = []; sym = ['o']*len(views)
        for n, view, col in zip(count(), views, colors):
            mks.append(vedo.Marker('*').c(col[0:-1]).legend(view))
        lb = vedo.LegendBox(mks, markers=sym, font=txt_font, 
                            width=leg_width/1.5, height=leg_height/1.5)
        
        vpt = MyFaceSelectingPlotter(N=2, axes=1,colors=colors, color_o=color_o, 
                                     views=views)
        vpt.add_icon(logo, pos=(0.1,1), size=0.25)
        vpt.add_callback("key press", vpt.on_key_press)
        vpt.add_callback("mouse click", vpt.select_cube_face)
        vpt.show(cl_mesh.mesh, linLine, orient_cube_clear, txt0, at=0)
        vpt.show(orient_cube, lb, vpt.msg, vpt.msg_face, at=1, azimuth=45, elevation=30, zoom=0.8, interactive=True)
        
        settings = {'proj_plane': plane, 
                    'ref_vect': ref_vect,
                    'ref_vectF': ref_vectF,
                    'orient_vect': pts,
                    'angle_deg': angle}
        
        self.roi_cube = {'cube': orient_cube,
                           'clear': orient_cube_clear}
        
        return  vpt.planar_views, settings, roi_cube

    def orient_manual(self, views, colors):

        im_orient = self.info['im_orientation']
        print('self.info[im_orientation]',im_orient)
        rotateY = False
        if im_orient == 'custom': 
            cust_angle = self.info['custom_angle']
            rotateY = True

        ext_ch = self.get_ext_int_chs()
        if isinstance(ext_ch, str) and ext_ch == 'independent':
            ext_ch = self.obj_imChannels[list(self.obj_imChannels.keys())[0]]

        mesh_ext = self.obj_meshes[ext_ch.channel_no+'_tiss']
        pos = mesh_ext.mesh.center_of_mass()

        side = max(self.get_maj_bounds())
        color_o = [152,251,152,255]
        orient_cube = vedo.Cube(pos=pos, side=side, c=color_o[:-1])
        orient_cube.linewidth(1).force_opaque()
        orient_cube_clear = orient_cube.clone().alpha(0.5)

        roi_cube = {'pos': pos, 
                        'side': side, 
                        'color': color_o[:-1], 
                        'rotateY': rotateY}
        
        if rotateY: 
            orient_cube.rotate_y(-cust_angle)
            roi_cube['custom_angle'] = -cust_angle
        else: 
            roi_cube['custom_angle'] = 0

        orient_cube_clear, rotX, rotY, rotZ = modify_cube(filename = self.user_organName,
                                                            txt = 'set the ROI Orientation', 
                                                            mesh = mesh_ext.mesh,
                                                            orient_cube = orient_cube_clear,
                                                            centre = pos, 
                                                            option = [True,True,True,True,True,True],
                                                            zoom=0.5)
        
        roi_cube['rotate_user'] = {'rotX': sum(rotX), 
                                     'rotY': sum(rotY), 
                                     'rotZ': sum(rotZ)}
        
        orient_cube_clear.pos(pos)

        #Rotate the cube using the x y and z values 
        orient_cube.rotate_x(sum(rotX)).rotate_y(sum(rotY)).rotate_z(sum(rotZ)).pos(pos)
        
        txt0 = vedo.Text2D(self.user_organName+' - Reference cube and mesh to select planar views in ROI (organ)...', c=txt_color, font=txt_font, s=txt_size)
        
        mks = []; sym = ['o']*len(views)
        for n, view, col in zip(count(), views, colors):
            mks.append(vedo.Marker('*').c(col[0:-1]).legend(view))
        lb = vedo.LegendBox(mks, markers=sym, font=txt_font, 
                            width=leg_width/1.5, height=leg_height/1.5)
        
        vpt = MyFaceSelectingPlotter(N=2, axes=1,colors=colors, color_o=color_o, 
                                     views=views)
        vpt.add_icon(logo, pos=(0.1,1), size=0.25)
        vpt.add_callback("key press", vpt.on_key_press)
        vpt.add_callback("mouse click", vpt.select_cube_face)
        vpt.show(mesh_ext.mesh, orient_cube_clear, txt0, at=0)
        vpt.show(orient_cube, lb, vpt.msg, vpt.msg_face, at=1, azimuth=45, elevation=30, zoom=0.8, interactive=True)
        
        settings = {}
        
        self.roi_cube = {'cube': orient_cube,
                           'clear': orient_cube_clear}

        return vpt.planar_views, settings, roi_cube

    def get_ref_vectors(self, planar_views, roi_cube):

        rotX = 0; rotY = 0; rotZ = 0
        #Check the rotations that were made 
        if 'rotate_x' in roi_cube.keys(): 
            rotX = roi_cube['rotate_x']
        if 'rotate_y' in roi_cube.keys(): 
            rotY = roi_cube['rotate_y']
        if 'rotate_z' in roi_cube.keys(): 
            rotZ = roi_cube['rotate_z']
        print(rotX, rotY, rotZ)
        if 'rotate_user' in roi_cube.keys(): 
            rotX = rotX + roi_cube['rotate_user']['rotX']
            rotY = rotY + roi_cube['rotate_user']['rotY']
            rotZ = rotZ + roi_cube['rotate_user']['rotZ']

        for view in planar_views: 
            pl_normal_rot = planar_views[view]['pl_normal']
            #Get the unrotated normal
            pl_normal_o = unit_vector(new_normal_3DRot(pl_normal_rot, [-rotX], [-rotY], [-rotZ]))
            planar_views[view]['normal_unrotated'] = pl_normal_o
            
            #Find the axis of the unrotated normal and set a reference vector
            index = np.argmax(np.absolute(pl_normal_o))
            if index == 0: 
                ref_vector = np.array([0.0, 1.0, 0.0])
                proj_plane = 'x'
            elif index == 1: 
                ref_vector = np.array([0.0, 0.0, 1.0])
                proj_plane = 'y'
            else: 
                ref_vector = np.array([1.0, 0.0, 0.0])
                proj_plane = 'z'

            #Rotate the reference vector
            #Get the unrotated normal
            ref_vect_rot = new_normal_3DRot(ref_vector, [rotX], [rotY], [rotZ])
            planar_views[view]['ref_vector'] = ref_vect_rot
            planar_views[view]['proj_plane'] = proj_plane

        return planar_views

    def get_maj_bounds(self):
        x_b = 0; y_b = 0; z_b = 0
        for mesh_o in self.obj_meshes:
            m_mesh = self.obj_meshes[mesh_o].mesh
            x1,x2,y1,y2,z1,z2 = m_mesh.bounds()
            if x2-x1 > x_b:
                x_b = x2-x1
            if y2-y1 > y_b:
                y_b = y2-y1
            if z2-z1 > z_b:
                z_b = z2-z1   
                
        return [x_b, y_b, z_b]
        
    def organ_vol_iso(self): 

        iso_vols = []; vol_settings = {}
        if hasattr(self, 'vol_iso'): 
            for n, vol in enumerate(self.vol_iso.keys()):
                iso_vols.append(self.vol_iso[vol])
                name = vol+': '+self.mC_settings['setup']['name_chs'][vol]
                color = self.mC_settings['setup']['color_chs'][vol]
                vol_settings = {'color': {n: color}, 'name':{n: name}}
        
        return iso_vols, vol_settings

    #Get all the set mH variables in __init__
    def get_notes(self):
        return self.info['user_organNotes']

    def get_custom_angle(self):
        return self.info['custom_angle']
    
    def get_resolution(self):
        return self.info['resolution']

    def get_units_resolution(self):
        return self.info['units_resolution']

    def get_stage(self):
        return self.info['stage']

    def get_strain(self):
        return self.info['strain']

    def get_genotype(self):
        return self.info['genotype']

    def get_manipulation(self): 
        return self.info['manipulation']
    
    def get_dir_res(self):
        return self.dir_res()

    def get_direc(self, name:str):
        return self.dir_res(dir=name)

class ImChannel(): #channel
    'morphoHeart Image Channel Class (this class will be used to contain the images as tiffs that have been'
    'closed and the resulting s3s that come up from each channel'
    
    def __init__(self, organ:Organ, ch_name:str):#
        
        self.parent_organ = organ
        self.parent_organ_name = organ.user_organName
        self.channel_no = ch_name
        self.user_chName = organ.mH_settings['setup']['name_chs'][ch_name]
        self.ch_relation = organ.mH_settings['setup']['chs_relation'][ch_name]
        
        if self.channel_no not in organ.imChannels.keys():   
            print('>> New Channel-', self.channel_no)
            self.new_ImChannel()
        else: 
            print('>> Loading Channel-', self.channel_no)
            self.load_channel()

    def new_ImChannel(self):
        organ = self.parent_organ
        ch_name = self.channel_no

        self.to_mask = organ.mH_settings['setup']['mask_ch'][ch_name]
        self.resolution = organ.info['resolution']
        self.dir_cho = organ.img_dirs[ch_name]['image']['dir']       
        if self.to_mask:
            self.dir_mk = organ.img_dirs[ch_name]['mask']['dir'] 
        self.masked = False
        self.shape = self.im().shape
        self.process = ['Init']
        self.contStack = {}

        organ.mH_settings['setup']['keep_largest'][ch_name] = {}
        organ.mH_settings['setup']['alpha'][ch_name] = {}

        self.save_channel(im_proc=self.im_proc(new=True))
        organ.add_channel(imChannel=self)
        # organ.save_organ()
        
    def im(self):#
        im = io.imread(str(self.dir_cho))
        if not isinstance(im, np.ndarray):
            print('>> Error: morphoHeart was unable to load tiff.\n>> Directory: ',str(self.dir_cho))
            alert('error_beep')
        else: 
            return im
    
    def im_proc(self, new=False):#
        if new: 
            im_proc =  np.copy(self.im())  
        else: 
            if hasattr(self, 'dir_stckproc'):
                dir_stck = self.parent_organ.dir_res(dir='s3_numpy') / self.dir_stckproc
                im_proc = np.load(str(dir_stck))#self.dir_stckproc))
                if not isinstance(im_proc, np.ndarray):
                    print('>> Error: morphoHeart was unable to load processed tiff.\n>> Directory: ',str(dir_stck))
                    alert('error_beep')
            else: 
                im_proc =  np.copy(self.im())      
        return im_proc
        
    def load_channel(self):
        organ = self.parent_organ
        ch_name = self.channel_no
        
        self.to_mask = organ.imChannels[ch_name]['to_mask']
        self.resolution = organ.imChannels[ch_name]['resolution']
        self.dir_cho = Path(organ.imChannels[ch_name]['dir_cho'])
        if self.to_mask:
                self.dir_mk = Path(organ.imChannels[ch_name]['dir_mk'])
        self.masked = organ.imChannels[ch_name]['masked']
        self.shape = tuple(organ.imChannels[ch_name]['shape'])
        if 'shape_s3' in organ.imChannels[ch_name].keys():
            self.shape_s3 = tuple(organ.imChannels[ch_name]['shape_s3'])
        self.process = organ.imChannels[ch_name]['process']
        contStack_dict = organ.imChannels[ch_name]['contStack']
        for cont in contStack_dict.keys():
            contStack_dict[cont]['s3_file'] = Path(contStack_dict[cont]['s3_file'])
        self.contStack = contStack_dict
        self.dir_stckproc = organ.imChannels[ch_name]['dir_stckproc']

    def get_channel_no(self):
        return self.channel_no

    def get_resolution(self):
        return self.resolution

    def get_shape(self):
        return self.shape
    
    def add_contStack(self, contStack):
        # Check first if the contStack has been already added to the channel
        if contStack.cont_type not in self.contStack.keys():
            new = True
        else: 
            new = False
            
        if new: 
            contStack_dict = copy.deepcopy(contStack.__dict__)
            contStack_dict.pop('im_channel', None)
            self.contStack[contStack.cont_type] = contStack_dict
        else: # just update im_proc 
            self.contStack[contStack.cont_type]['process'] = contStack.process
     
    def maskIm(self):#
        # Workflow process
        workflow = self.parent_organ.workflow['morphoHeart']
        process = ['ImProc', self.channel_no, 'A-MaskChannel','Status']

        #Load images
        im_o = np.copy(self.im_proc(new=True))
        im_mask = io.imread(str(self.dir_mk))
        #Process
        print('\n---- Masking! ----')
        if self.shape == im_mask.shape:
            #Check the dimensions of the mask with those of the image
            im_o[im_mask == False] = 0
            self.masked = True

            #Update organ workflow
            self.parent_organ.update_mHworkflow(process, update = 'DONE')
            
            process_up = ['ImProc', self.channel_no,'Status']
            if get_by_path(workflow, process_up) == 'NI':
                self.parent_organ.update_mHworkflow(process_up, update = 'Initialised')

            #Update channel process
            self.process.append('Masked')
            
            #Update organ imChannels
            self.parent_organ.imChannels[self.channel_no]['masked'] = True
            self.parent_organ.add_channel(self)
                
            process_up2 = ['ImProc','Status']
            if get_by_path(workflow, process_up2) == 'NI':
                self.parent_organ.update_mHworkflow(process_up2, update = 'Initialised')
            
            #Save channel
            self.save_channel(im_proc=im_o)
            
        else: 
            print('For some reason self.shape != im_mask.shape!')
            alert('bubble')

    def closeContours_auto(self, gui_param, gui_plot, win):
        from .mH_funcContours import autom_close_contours
        # Workflow process
        workflow = self.parent_organ.workflow['morphoHeart']
        process = ['ImProc', self.channel_no, 'B-CloseCont','Steps','A-Autom','Status']
        # Load image
        im_proc = self.im_proc()

        #Process
        print('\n---- Closing Contours Auto! ----')
        #> Close contours Automatically
        im_proc = autom_close_contours(stack = im_proc, ch = self.channel_no,
                                        gui_param = gui_param, gui_plot = gui_plot, win = win)

        #Update organ imChannels
        self.parent_organ.add_channel(self)
        #Update channel process
        slc_first_py = gui_param['start_slc']
        slc_last_py = gui_param['end_slc']
        self.process.append('ClosedCont-Auto - Slc'+str(slc_first_py)+'-'+str(slc_last_py))

        #Update organ workflow
        self.parent_organ.update_mHworkflow(process, update = 'Initialised')
        getattr(win, 'autom_close_'+self.channel_no+'_done').setEnabled(True)
        process_up = ['ImProc',self.channel_no,'B-CloseCont','Status']
        if get_by_path(workflow, process_up) == 'NI':
            self.parent_organ.update_mHworkflow(process_up, update = 'Initialised')
        process_up = ['ImProc', self.channel_no,'Status']
        if get_by_path(workflow, process_up) == 'NI':
            self.parent_organ.update_mHworkflow(process_up, update = 'Initialised')
        process_up2 = ['ImProc','Status']
        if get_by_path(workflow, process_up2) == 'NI':
            self.parent_organ.update_mHworkflow(process_up2, update = 'Initialised')

        #Save channel
        self.save_channel(im_proc=im_proc)

        #Plot slice range
        win.plot_all_slices(ch = self.channel_no, slice_range = (slc_first_py, slc_last_py+1))

    def create_chS3s(self, layerDict:dict, win, cont_list=['int', 'ext', 'tiss']):
        # Workflow process
        workflow = self.parent_organ.workflow['morphoHeart']
        process = ['ImProc',self.channel_no,'D-S3Create','Status']
        
        try: 
            win.prog_bar_range(0,3)
        except: 
            pass

        dirs_cont = []; shapes_s3 = []
        aa = 0
        for cont in cont_list:
            win.win_msg('Creating masked stacks for each contour of Channel '+self.channel_no[-1]+' ('+str(aa)+'/'+str(len(cont_list))+').')
            s3 = ContStack(im_channel=self, cont_type=cont, layerDict=layerDict)#new=True,
            self.add_contStack(s3)
            path2file = self.parent_organ.dir_res(dir='s3_numpy') / s3.s3_file
            dirs_cont.append(path2file.is_file())
            shapes_s3.append(s3.shape_s3)
            #Update organ workflow
            process_cont = ['ImProc',self.channel_no,'D-S3Create','Info',cont,'Status']
            self.parent_organ.update_mHworkflow(process_cont, update = 'DONE')
            aa+=1
            try: 
                win.prog_bar_update(aa)
            except:
                pass

        win.win_msg('Creating masked stacks for each contour of Channel '+self.channel_no[-1]+' ('+str(aa)+'/'+str(len(cont_list))+').')
        #Update organ workflow
        if all(flag for flag in dirs_cont):
            if shapes_s3.count(shapes_s3[0]) == len(shapes_s3):
                self.shape_s3 = s3.shape_s3
                self.parent_organ.update_mHworkflow(process, update = 'DONE')
            else: 
                print('>> Error: self.shape_s3 = s3.shape')
        
        #Update channel process
        self.process.append('CreateS3')
        
        #Update organ imChannel
        self.parent_organ.add_channel(self)
        # self.parent_organ.save_organ()
        
        process_up2 = ['ImProc','Status']
        if get_by_path(workflow, process_up2) == 'NI':
            self.parent_organ.update_mHworkflow(process_up2, update = 'Initialised')
            
    def load_chS3s (self, cont_types:list):
        for cont in cont_types:
            s3 = ContStack(im_channel=self, cont_type=cont)#, new=False)
            setattr(self, 's3_'+cont, s3)
            self.add_contStack(s3)
        
        #Update channel process
        self.process.append('LoadS3')
        
        #Update organ imChannel
        self.parent_organ.add_channel(self)
        # self.parent_organ.save_organ()

    def trimS3(self, cuts, cont, cuts_out): 
        
        workflow = self.parent_organ.workflow['morphoHeart']
        process = ['ImProc', self.channel_no, 'E-TrimS3','Status']

        #Load s3s
        self.load_chS3s(cont_types=[cont])
        s3 = getattr(self, 's3_'+cont)

        #Process
        print('\n---- Trimming S3s! ----')                             
        if len(cuts) == 1:
            pl = cuts_out[cuts[0]]['plane_info_image']
            s3.cutW1Plane(pl, cuts[0])
           
        if len(cuts) == 2:
            pl1 = cuts_out['bottom']['plane_info_image']
            pl2 = cuts_out['top']['plane_info_image']
            s3.cutW2Planes(pl1, pl2)

        #Update channel process
        self.process.append('TrimS3')
        
        #Update organ imChannels
        self.parent_organ.add_channel(self)
        
    def s32Meshes(self, cont_type:str, keep_largest=None, rotateZ_90=None, new_set=False):

        mesh_prop = {'keep_largest': keep_largest, 'rotateZ_90': rotateZ_90}
        try: 
            try: 
                mesh = Mesh_mH(imChannel = self, mesh_type = cont_type, 
                            mesh_prop = mesh_prop, new_set = new_set)
                print('!>> ',mesh.__dict__)
                return True
            except RuntimeError: 
                return False
        except: 
            return False

    def save_channel(self, im_proc):
        organ_name = self.parent_organ.user_organName
        im_name = organ_name + '_StckProc_' + self.channel_no + '.npy'
        im_dir = self.parent_organ.dir_res(dir='s3_numpy') / im_name
        np.save(im_dir, im_proc)
        if not im_dir.is_file():
            print('>> Error: Processed channel was not saved correctly!\n>> File: '+im_name)
            alert('error_beep')
        else: 
            print('>> Processed channel saved correctly! - ', im_name)
            alert('countdown')
            self.dir_stckproc = im_name
    
    def ch_clean(self, s3_mask, s3, inverted, plot_settings):#
        """
        Function to clean channel contour using other channel as a mask
        """
        plot, im_every = plot_settings 
        # What happens if the s3() are None? 
        s3_s = s3.s3()
        if not isinstance(s3_s, np.ndarray): 
            print('>> Error: Not isinstance(s3_s, np.array)')
            alert('clown')
            return
        s3_mask_s = s3_mask.s3()
        if not isinstance(s3_mask_s, np.ndarray): 
            print('>> Error: not isinstance(s3_mask_s, np.array)')
            alert('clown')
            return
        s3_bits = np.zeros_like(s3_s, dtype='uint8')
        s3_new =  np.zeros_like(s3_s, dtype='uint8')

        index = list(s3.shape_s3).index(min(s3.shape_s3))
        if index == 2:
            for slc in range(s3.shape_s3[2]):
                mask_slc = s3_mask_s[:,:,slc]
                toClean_slc = s3_s[:,:,slc]

                if inverted:
                    # Invert ch to use as mask 
                    inv_slc = np.where((mask_slc==0)|(mask_slc==1), mask_slc^1, mask_slc)
                else: 
                    # Keep ch to use as mask as it is
                    inv_slc = np.copy(mask_slc)

                # inverted_mask or mask AND ch1_2clean
                toRemove_slc = np.logical_and(toClean_slc, inv_slc)
                # Keep only the clean bit
                cleaned_slc = np.logical_xor(toClean_slc, toRemove_slc)

                if plot and slc in list(range(0,s3.shape_s3[0],im_every)):
                    print('Plotting! slc:', slc)

                s3_bits[:,:,slc] = toRemove_slc
                s3_new[:,:,slc] = cleaned_slc
                
            s3_new = s3_new.astype('uint8')
            s3.s3_save(s3_new)
            alert('whistle')   
        else:
            print('>> (ch_clean) Index different to 2, check!')
            alert('error_beep')
                             
    def slc_plot (self, slc, mask_slc, toClean_slc, toRemove_slc, cleaned_slc, inverted):
        """
        Function to plot mask, original image and result
        """
        
        if inverted: 
            txt = ['ch0_inv','ch1','ch0_inv AND ch1','ch0_inv AND ch1\nxOR ch1']
        else: 
            txt = ['ch0','ch1','ch0 AND ch1','ch0 AND ch1\nxOR ch1']
       
        #Plot
        fig, ax = plt.subplots(1, 4, figsize = (10,2.5))
        fig.suptitle("Slice:"+str(slc), y=1.05, weight="semibold")
        ax[0].imshow(mask_slc)
        ax[1].imshow(toClean_slc)
        ax[2].imshow(toRemove_slc)
        ax[3].imshow(cleaned_slc)
        for num in range(0,4,1):
            ax[num].set_title(txt[num])
            ax[num].set_xticks([])
            ax[num].set_yticks([])

        plt.show()
        
class ImChannelNS(): #channel

    'morphoHeart Image Channel Negative Space'
    
    def __init__(self, organ:Organ, ch_name:str):#, new=True):

        self.parent_organ = organ
        self.parent_organ_name = organ.user_organName
        self.channel_no = ch_name
        self.user_chName = organ.mH_settings['setup'][ch_name]['user_nsChName']
        self.ch_relation = 'negative-space'

        if self.channel_no not in organ.imChannelNS.keys():
            print('>> New ChannelNS')
            self.new_ImChannelNS()
        else: 
            print('>> Loading ChannelNS')
            self.load_channel()
    
    def new_ImChannelNS (self):
        organ = self.parent_organ
        ch_name = self.channel_no
        
        self.resolution = organ.info['resolution']
        self.process = ['Init']
        self.contStack = {}
        
        # external contour
        ext_s3_name = organ.mH_settings['setup'][ch_name]['ch_ext'][0]
        ext_s3_type = organ.mH_settings['setup'][ch_name]['ch_ext'][1]
        # internal contour
        int_s3_name = organ.mH_settings['setup'][ch_name]['ch_int'][0]
        int_s3_type = organ.mH_settings['setup'][ch_name]['ch_int'][1]

        organ.mH_settings['setup'][ch_name]['keep_largest'] = {}
        organ.mH_settings['setup'][ch_name]['alpha'] = {}
        
        self.setup_NS = {'ext':{'name': ext_s3_name, 'type': ext_s3_type}, 
                            'int':{'name': int_s3_name, 'type': int_s3_type}}
        
        organ.add_channelNS(imChannelNS=self)
        organ.check_status(process='ImProc')
        # organ.save_organ()
    
    def load_channel(self):
        
        organ = self.parent_organ
        ch_name = self.channel_no
        
        self.resolution = organ.imChannelNS[ch_name]['resolution']
        self.process = organ.imChannelNS[ch_name]['process']
        contStack_dict = organ.imChannelNS[ch_name]['contStack']
        for cont in contStack_dict.keys():
            contStack_dict[cont]['s3_file'] = Path(contStack_dict[cont]['s3_file'])
            contStack_dict[cont]['shape_s3'] = tuple(contStack_dict[cont]['shape_s3'])
        self.contStack = contStack_dict
        self.setup_NS = organ.imChannelNS[ch_name]['setup_NS']
        
    def create_chNSS3s(self, win, plot_settings=(False, None)):

        organ = self.parent_organ
        win.win_msg('Creating masked stacks for each contour of channel '+self.channel_no+' (0/3).')
        ext_s3_name = self.setup_NS['ext']['name']
        ext_s3_type = self.setup_NS['ext']['type']
        ext_s3 = ContStack(im_channel=organ.obj_imChannels[ext_s3_name], 
                           cont_type=ext_s3_type)#, new=False)
        self.s3_ext = ext_s3
        self.add_contStack(ext_s3, cont_type = 'ext')
        
        win.win_msg('Creating masked stacks for each contour of channel '+self.channel_no+' (1/3).')
        int_s3_name = self.setup_NS['int']['name']
        int_s3_type = self.setup_NS['int']['type']
        int_s3 = ContStack(im_channel=organ.obj_imChannels[int_s3_name], 
                           cont_type=int_s3_type)#, new=False)
        self.s3_int = int_s3
        self.add_contStack(int_s3, cont_type = 'int')
        
        win.win_msg('Creating masked stacks for each contour of channel '+self.channel_no+' (2/3).')
        layerDict = organ.mH_settings['setup']['chNS']
        layerDict['plot_settings'] = plot_settings
        tiss_s3 = ContStack(im_channel = self, cont_type = 'tiss',
                            layerDict=layerDict)#,  new = True)
        self.s3_tiss = tiss_s3
        self.add_contStack(tiss_s3, cont_type = 'tiss')

        process = ['ImProc','chNS','Status']
        organ.update_mHworkflow(process, 'DONE')
        organ.check_status('ImProc')

    def load_chS3s (self, cont_types:list):
        for cont in cont_types:
            s3 = ContStack(im_channel=self, cont_type=cont)#, new=False)
            setattr(self, 's3_'+cont, s3)
            self.add_contStack(s3, cont)
        
        #Update channel process
        self.process.append('LoadS3')
        
        #Update organ imChannel
        self.parent_organ.add_channelNS(imChannelNS=self)
        # self.parent_organ.save_organ()
     
    def get_channel_no(self):
        return self.channel_no

    def get_resolution(self):
        return self.resolution
    
    def add_contStack(self, contStack, cont_type):
        # Check first if the contStack has been already added to the channel
        new = False
        if cont_type not in self.contStack.keys():
            new = True
            
        if new: 
            contStack_dict = {}
            contStack_dict['cont_type'] = cont_type
            contStack_dict['imfilled_name'] = contStack.imfilled_name
            contStack_dict['cont_name'] = contStack.cont_name
            contStack_dict['s3_file'] = contStack.s3_file
            # contStack_dict['s3_dir'] = contStack.s3_dir
            contStack_dict['shape_s3'] = contStack.shape_s3
            contStack_dict['process'] = contStack.process
            self.contStack[cont_type] = contStack_dict
        else: # just update process 
            self.contStack[cont_type]['process'] = contStack.process

    def create_s3_tiss(self, layerDict, plot_settings): 
        """
        Function to extract the negative space channel
        """        
        # Workflow process
        process = ['ImProc', self.channel_no,'D-S3Create','Status']
      
        plot, im_every = plot_settings 
        print('>> Extracting '+self.user_chName+'!')
        operation = layerDict['operation']

        extNS = self.s3_ext.s3()
        intNS = self.s3_int.s3()
        
        s3_bits = np.zeros_like(extNS, dtype='uint8')
        s3_new =  np.zeros_like(extNS, dtype='uint8')

        index = list(extNS.shape).index(min(extNS.shape))
        if index == 2:
            for slc in range(extNS.shape[2]):
                s3_intNS = intNS[:,:,slc]
                s3_extNS = extNS[:,:,slc]
                if operation == 'AND-XOR':# or operation == 'XOR': 
                    # inverted_mask or mask AND ch1_2clean
                    s3_AND = np.logical_and(s3_extNS, s3_intNS)
                    # Keep only the clean bit
                    s3_XOR = np.logical_xor(s3_extNS, s3_AND)

                if plot and slc in list(range(0,extNS.shape[0],im_every)):
                    print('Plotting! slc: ', slc)
                    # self.slc_plot(slc, s3_extNS, s3_intNS, s3_AND, s3_XOR)

                s3_bits[:,:,slc] = s3_AND
                s3_new[:,:,slc] = s3_XOR
                
            s3_new = s3_new.astype('uint8')
            alert('whistle')   
            
        else:
            print('>> Index different to 2, check!')
            alert('error_beep')

        #Update organ workflow
        self.parent_organ.update_mHworkflow(process, update = 'DONE')
        
        return s3_new
    
    def slc_plot (self, slc, s3_intNS, s3_extNS, s3_AND, s3_XOR): #Add this plot 
        """
        Function to plot mask, original image and result
        """
       
        txt = ['ext','int','extANDint','(extANDint)XORext']

        #Plot
        fig, ax = plt.subplots(1, 4, figsize = (10,2.5))
        fig.suptitle("Slice:"+str(slc), y=1.05, weight="semibold")
        ax[0].imshow(s3_extNS)
        ax[1].imshow(s3_intNS)
        ax[2].imshow(s3_AND)
        ax[3].imshow(s3_XOR)
        for num in range(0,4,1):
            ax[num].set_title(txt[num])
            ax[num].set_xticks([])
            ax[num].set_yticks([])

        plt.show()
        
    # func - s32Meshes
    def s32Meshes(self, cont_type:str, keep_largest=None, rotateZ_90=None, new_set=False):

        mesh_prop = {'keep_largest': keep_largest, 'rotateZ_90': rotateZ_90}

        mesh = Mesh_mH(imChannel = self, mesh_type = cont_type, 
                        mesh_prop = mesh_prop, new_set = new_set)
        
        self.parent_organ.check_status('MeshesProc')
        
class ContStack(): 
    'morphoHeart Contour Stack Class'
    
    def __init__(self, im_channel:Union[ImChannel,ImChannelNS], 
                             cont_type:str, layerDict={}):
        
        cont_types = ['int', 'ext', 'tiss']
        names = ['imIntFilledCont', 'imExtFilledCont', 'imAllFilledCont']

        index = cont_types.index(cont_type)
        self.cont_type = cont_type
        self.imfilled_name = names[index]
        self.im_channel = im_channel
        self.cont_name = im_channel.channel_no+'_'+self.cont_type
        
        parent_organ = im_channel.parent_organ

        if im_channel.channel_no == 'chNS':
            if self.cont_type == 'ext' or self.cont_type == 'int': 
                ch_name = im_channel.setup_NS[self.cont_type]['name']
                cont_name = im_channel.setup_NS[self.cont_type]['type']
                self.s3_file = parent_organ.user_organName + '_s3_' + ch_name + '_' + cont_name + '.npy'
            else: 
                self.s3_file = parent_organ.user_organName + '_s3_' + im_channel.channel_no + '_' + self.cont_type + '.npy'
        else: 
            self.s3_file = parent_organ.user_organName + '_s3_' + im_channel.channel_no + '_' + self.cont_type + '.npy'

        s3_file_dir =  parent_organ.dir_res(dir='s3_numpy') / self.s3_file 
        if self.cont_type not in self.im_channel.contStack.keys() or not s3_file_dir.is_file():
            if im_channel.channel_no == 'chNS':
                # TO DO! Add information about the masking process as attributes here!
                plot_settings = layerDict['plot_settings']
                s3 = im_channel.create_s3_tiss(layerDict=layerDict, plot_settings=plot_settings)
                parent_organ.mH_settings['setup'][im_channel.channel_no]['keep_largest'][self.cont_type] = {}
                parent_organ.mH_settings['setup'][im_channel.channel_no]['alpha'][self.cont_type] = {}
            else: 
                s3 = self.s3_create(layerDict = layerDict)
                parent_organ.mH_settings['setup']['keep_largest'][im_channel.channel_no][self.cont_type] = {}
                parent_organ.mH_settings['setup']['alpha'][im_channel.channel_no][self.cont_type] = {}
            self.s3_save(s3)
            self.shape_s3 = s3.shape
            self.process = ['Init']
        else: 
            # s3 = self.s3_create(layerDict = layerDict)
            s3 = self.s3()
            self.shape_s3 = s3.shape
            self.process = im_channel.contStack[cont_type]['process']
            self.process.append('Loaded')
    
    def s3_create(self, layerDict:dict):
        if isinstance(layerDict, dict): 
            x_dim = self.im_channel.shape[0]
            y_dim = self.im_channel.shape[1]
            z_dim = self.im_channel.shape[2]
            # print(x_dim, y_dim, z_dim)
            s3 = np.empty((y_dim,z_dim,x_dim+2))
            s3 = s3.astype('uint8')

        elif isinstance(layerDict, np.ndarray): 
            s3 = layerDict

        return s3
    
    def s3(self):
        dir_s3 = self.im_channel.parent_organ.dir_res(dir='s3_numpy') / self.s3_file

        if dir_s3.is_file():
            s3 = np.load(dir_s3)
        else: 
            print('>> Error: s3 file does not exist!\n>> File: '+self.s3_file)
            alert('error_beep')
            s3 = None
            
        return s3

    def s3_save(self, s3):
        organ = self.im_channel.parent_organ
        dir2save = organ.dir_res(dir='s3_numpy') / self.s3_file
        np.save(dir2save, s3)
        if not dir2save.is_file():
            print('>> Error: s3 file was not saved correctly!\n>> File: '+self.s3_file)
            alert('error_beep')
        else: 
            print('>> s3 file saved correctly! - ', self.im_channel.channel_no, '-', self.cont_type)
            alert('countdown')
                
    def cutW2Planes(self, pl1, pl2):
        """
        Function used to cut inflow AND outflow tract of the s3 mask (s3_cut) given as input
    
        """
        #Load s3 and resolution
        s32cut = self.s3()
        resolution = self.im_channel.resolution
        
        # Get dimensions of external stack
        xdim, ydim, zdim = s32cut.shape
        # Reshape stacks as a vector
        s3_cut_v = s32cut.reshape(-1)
    
        # Get vectors of x,y and z positions
        pix_coord_pos = np.where(s32cut >= 0)
        del s32cut
        # Trasform coordinate positions to um using resolution
        pix_um = np.transpose(np.asarray([pix_coord_pos[i]*resolution[i] for i in range(len(resolution))]))
        del pix_coord_pos
    
        normal_inf = unit_vector(pl1['pl_normal'])#pls_normal[0])
        normal_outf = unit_vector(pl2['pl_normal'])#pls_normal[1])
    
        # Find all the d values of pix_um
        d_pix_um_Inf = np.dot(np.subtract(pix_um,np.array(pl1['pl_centre'])),np.array(normal_inf))
        d_pix_um_Outf = np.dot(np.subtract(pix_um,np.array(pl2['pl_centre'])),np.array(normal_outf))
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
        
        # Save new s3
        self.s3_save(s3f_cut)
        alert('woohoo')
    
        # return s3f_cut
    
    # func - cutInfOrOutfOptMx
    def cutW1Plane (self, pl, cut):
        """
        Function used to cut inflow OR outflow tract of the s3 mask (s3_cut) given as input
    
        """
    
        #Load s3 and resolution
        s32cut = self.s3()
        resolution = self.im_channel.resolution
    
        # Get dimensions of external stack
        xdim, ydim, zdim = s32cut.shape
        # Reshape stacks as a vector
        s3_cut_v = s32cut.reshape(-1)
    
        # Get vectors of x,y and z positions
        pix_coord_pos = np.where(s32cut >= 0)
        del s32cut
        # Trasform coordinate positions to um using resolution
        pix_um = np.transpose(np.asarray([pix_coord_pos[i]*resolution[i] for i in range(len(resolution))]))
        del pix_coord_pos
    
        normal  = unit_vector(pl['pl_normal'])
        # Find all the d values of pix_um
        d_pix_um = np.dot(np.subtract(pix_um,np.array(pl['pl_centre'])),np.array(normal))
    
        # Clear vector d_pix_um using only those that are 1 in stack
        d_pve_pix_um = s3_cut_v*d_pix_um
        del pix_um
    
        # Duplicate s3f_v to initialise stacks without inflow/outflow
        s3f_all_v = np.copy(s3_cut_v)
        s3f_all_v.astype('uint8')
        del s3_cut_v
    
        # Find all positions in d_pve_pix_um that are at either side of the planes (outside of mesh)
        if cut == 'inflow tract' or cut == 'bottom':
            pos_outside = np.where(d_pve_pix_um < 0)[0]
        elif cut == 'outflow tract' or cut == 'top':
            pos_outside = np.where(d_pve_pix_um > 0)[0]
        del d_pve_pix_um
    
        # Remove the points that are outside of the mesh (inflow/outflow)
        s3f_all_v[pos_outside] = 0
        del pos_outside
    
        # Reshape vector into matrix/stack
        s3f_cut = s3f_all_v.reshape((xdim, ydim, zdim))
        del s3f_all_v
        
        # Save new s3
        self.s3_save(s3f_cut)
        alert('woohoo')

class Mesh_mH():
    'morphoHeart Mesh Class'
    
    def __init__(self, imChannel:ImChannel, mesh_type:str, mesh_prop:dict, 
                 new_set=False):#, new=True):
        
        self.parent_organ = imChannel.parent_organ
        self.imChannel = imChannel
        self.channel_no = imChannel.channel_no
        self.user_meshName = self.parent_organ.mH_settings['setup']['name_chs'][self.channel_no]
        self.mesh_type = mesh_type
        self.legend = self.user_meshName+'_'+self.mesh_type
        self.name = self.channel_no +'_'+self.mesh_type
        self.resolution = imChannel.get_resolution()

        if self.name not in self.parent_organ.meshes.keys() or new_set:
            print('>> New mesh - ', self.name)
            new = True
            self.new_mesh(mesh_prop)
        else: 
            print('>> Reload mesh - ', self.name)
            new = False
            self.reload_mesh(mesh_prop, new_set)
            
        self.mesh.color(self.color)
        self.mesh.alpha(self.alpha)
        if new or new_set: 
            self.save_mesh()

    def new_mesh(self, mesh_prop):#

        if mesh_prop['keep_largest'] != None: 
            keep_largest = mesh_prop['keep_largest']
        else: 
            try: 
                keep_largest = self.parent_organ.mH_settings['wf_info']['keep_largest'][self.channel_no][self.mesh_type]
            except: 
                keep_largest = False

        if mesh_prop['rotateZ_90'] != None: 
            rotateZ_90 = mesh_prop['rotateZ_90']
        else:
            rotateZ_90 = self.parent_organ.mH_settings['setup']['rotateZ_90']

        self.keep_largest = keep_largest
        self.rotateZ_90 = rotateZ_90
        print('self.keep_largest:', self.keep_largest, ' - self.rotateZ_90:', self.rotateZ_90)
        
        self.create_mesh(keep_largest = self.keep_largest, rotateZ_90 = self.rotateZ_90)
        if self.channel_no != 'chNS': 
            self.color = self.parent_organ.mH_settings['setup']['color_chs'][self.channel_no][self.mesh_type]
        else: 
            self.color = self.parent_organ.mH_settings['setup'][self.channel_no]['color_chns'][self.mesh_type]
        self.alpha = 0.05
        
        #Update settings
        if self.channel_no != 'chNS': 
            set_proc = [self.channel_no, self.mesh_type]
            self.parent_organ.update_settings(['setup','keep_largest']+set_proc, self.keep_largest, 'mH')
            self.parent_organ.update_settings(['setup','alpha']+set_proc, self.alpha, 'mH')
        else: 
            self.parent_organ.update_settings(['setup',self.channel_no,'keep_largest',self.mesh_type], self.keep_largest, 'mH')
            self.parent_organ.update_settings(['setup',self.channel_no,'alpha',self.mesh_type], self.alpha, 'mH')

        self.dirs = {'mesh': None, 'arrays': None}
        self.parent_organ.check_status(process = 'MeshesProc')
        self.mesh_meas = {}

    def reload_mesh(self, mesh_prop, new_set):#

        if self.channel_no != 'chNS': 
            self.color = self.parent_organ.mH_settings['setup']['color_chs'][self.channel_no][self.mesh_type]
            self.alpha = self.parent_organ.mH_settings['setup']['alpha'][self.channel_no][self.mesh_type]
        else:
            self.color = self.parent_organ.mH_settings['setup'][self.channel_no]['color_chns'][self.mesh_type]
            self.alpha = self.parent_organ.mH_settings['setup'][self.channel_no]['alpha'][self.mesh_type]

        self.s3_dir = self.imChannel.contStack[self.mesh_type]['s3_file']

        if new_set: 
            if mesh_prop['keep_largest'] != None: 
                keep_largest = mesh_prop['keep_largest']
            else: 
                try: 
                    keep_largest = self.parent_organ.mH_settings['wf_info']['keep_largest'][self.channel_no][self.mesh_type]
                except: 
                    keep_largest = False

            if mesh_prop['rotateZ_90'] != None: 
                rotateZ_90 = mesh_prop['rotateZ_90']
            else:
                rotateZ_90 = self.parent_organ.mH_settings['setup']['rotateZ_90']

            self.keep_largest = keep_largest
            self.rotateZ_90 = rotateZ_90
            print('self.keep_largest:', self.keep_largest, ' - self.rotateZ_90:', self.rotateZ_90)
            print('>> Re-creating mesh -', self.name)
            self.create_mesh(keep_largest = keep_largest, rotateZ_90 = rotateZ_90)
        else: 
            if self.channel_no != 'chNS':
                try: 
                    self.keep_largest = self.parent_organ.mH_settings['wf_info']['keep_largest'][self.channel_no][self.mesh_type]
                except: 
                    self.keep_largest = False
            else: 
                try: 
                    self.keep_largest = self.parent_organ.mH_settings['setup'][self.channel_no]['keep_largest'][self.mesh_type]
                except:
                    self.keep_largest = False

            self.rotateZ_90 = self.parent_organ.mH_settings['setup']['rotateZ_90']
            print('self.keep_largest:', self.keep_largest, ' - self.rotateZ_90:', self.rotateZ_90)
            print('>> Loading mesh-', self.name)
            self.load_mesh()
            if self.name in self.parent_organ.objects['Centreline'].keys():
                if self.parent_organ.workflow['morphoHeart']['MeshesProc']['C-Centreline']['buildCL']['Status'] == 'DONE':
                    self.set_centreline()

        self.dirs = self.parent_organ.meshes[self.name]['dirs']
        if self.dirs['mesh'] != None: 
            for n, meas_mesh in enumerate(self.dirs['mesh'].keys()):
                if n == 0:
                    self.mesh_meas = {}
                    n_type = meas_mesh.split('(')[1][:-1]
                if 'ball' in meas_mesh:
                    m_ball = self.balloon_mesh(n_type = n_type)
                    self.mesh_meas[meas_mesh] = m_ball
                elif 'thck' in meas_mesh: 
                    m_thck = self.thickness_mesh(n_type = n_type)
                    self.mesh_meas[meas_mesh] = m_thck
        else: 
            self.mesh_meas = {}

    def create_mesh(self, keep_largest:bool, rotateZ_90:bool):
        # Extract vertices, faces, normals and values of each mesh
        s3_type = 's3_'+self.mesh_type
        try: 
            s3 = getattr(self.imChannel, s3_type)
        except: 
            self.imChannel.load_chS3s(cont_types=[self.mesh_type])
            s3 = getattr(self.imChannel, s3_type)
            
        self.s3_dir = s3.s3_file
        s3s3 = s3.s3()
        verts, faces, _, _ = measure.marching_cubes(s3s3, spacing=self.resolution, method='lewiner')
    
        # Create meshes
        mesh = vedo.Mesh([verts, faces])
        if keep_largest:
            mesh = mesh.extractLargestRegion()
        if rotateZ_90:
            mesh.rotateZ(-90)
        mesh.legend(self.legend).wireframe()
        self.mesh = mesh
    
    def load_mesh(self):
        parent_organ = self.parent_organ
        mesh_name = parent_organ.user_organName+'_'+self.name+'.vtk'
        mesh_dir = parent_organ.dir_res(dir='meshes')  / mesh_name
        mesh_out = vedo.load(str(mesh_dir))
        mesh_out.legend(self.legend).wireframe()
        self.dir_out = mesh_name
        self.mesh = mesh_out

    def save_mesh(self, m_type='self', ext='.vtk'):
        parent_organ = self.parent_organ
        organ_name = parent_organ.user_organName
        if m_type == 'self':
            if ext != '.vtk':
                mesh_name = organ_name+'_'+self.legend+ext
            else: #== .vtk
                mesh_name = organ_name+'_'+self.name+ext
            mesh_dir = parent_organ.dir_res(dir='meshes') / mesh_name
            self.dir_out = mesh_name
            mesh_out = self.mesh
            
        elif 'ball' in m_type or 'thck' in m_type: 
            print('Saving thickness mesh: ', m_type, '('+self.name+')')
            if ext != '.vtk':
                mesh_name = parent_organ.user_organName+'_'+self.legend+'_CM'+m_type+ext
            else: #== .vtk
                mesh_name = parent_organ.user_organName+'_'+self.name+'_CM'+m_type+ext
            mesh_dir = parent_organ.dir_res(dir='meshes') / mesh_name
            mesh_out = self.mesh_meas[m_type]
            if self.dirs['mesh'] == None: 
                self.dirs['mesh'] = {m_type: mesh_name}
            else:
                self.dirs['mesh'][m_type] = mesh_name
            
        mesh_out.write(str(mesh_dir))
        print('>> Mesh '+mesh_name+' has been saved!')
        alert('countdown')        
        self.parent_organ.add_mesh(self)
    
    def save_array(self, array, m_type):
        
        parent_organ = self.parent_organ        
        # title = parent_organ.user_organName+'_'+self.legend+'--'+m_type
        title = parent_organ.user_organName+'_'+self.name+'_CM'+m_type
        np2save_dir = self.parent_organ.dir_res(dir='csv_all') / title
        np.save(np2save_dir, array)
        if self.dirs['arrays'] == None: 
            self.dirs['arrays'] = {m_type: title}
        else:
            self.dirs['arrays'][m_type] = title
            
        np2save_dirf = Path(str(np2save_dir)+'.npy')
        if not np2save_dirf.is_file():
            print('>> Error: Array was not saved correctly!\n>> File: '+title)
            alert('error_beep')
        else: 
            print('>> Project settings file saved correctly!\n>> File: '+title)
            alert('countdown')
            
    def mesh4CL(self):
        """
        Function that cleans and smooths meshes given as input to get centreline using VMTK
    
        """
        mesh4cl = self.mesh.clone()
        print('>> Cleaning mesh '+self.legend)
        print("\t- Original number of points making up mesh: ", mesh4cl.NPoints())
        # Reduce the number of points that make up the mesh
        mesh4cl.subsample(fraction = 0.005)#tol=0.005)
        print("\t- Number of points after cleaning surface: ",mesh4cl.NPoints(),'\n- Smoothing mesh...', self.legend)
                
        # Smooth mesh
        mesh4cl_cut = mesh4cl.clone().smooth_mls_2d(f=0.2)
        mesh4cl_cut.legend(self.legend+"-C&S").color(self.color)
        print('>> Mesh smoothed!')
        alert('woohoo')

        return mesh4cl_cut
    
    def set_centreline(self):
        try: 
            cl_info = self.parent_organ.objects['Centreline'][self.name]
            self.centreline_info = cl_info
            print('>> Centerline has been set for ', self.name)
        except: 
            print('>> No centreline has been created for this mesh - ', self.name)

    def get_channel_no(self):
        return self.channel_no
    
    def get_user_meshName(self):
        return self.user_meshName
    
    def get_legend(self):
        return self.mesh.legend
    
    def change_user_meshName(self, new_user_meshName):
        self.user_meshName = new_user_meshName
    
    def get_mesh_type(self):
        return self.mesh_type
    
    def get_imChannel(self):
        return self.imChannel
    
    def get_organ(self):
        return self.parent_organ
    
    def set_alpha(self, mesh_alpha):      
        self.mesh.alpha(mesh_alpha)
        self.alpha = mesh_alpha
        #Update settings
        if self.channel_no != 'chNS': 
            self.parent_organ.update_settings(['setup','alpha',self.channel_no, self.mesh_type], self.alpha, 'mH')
        else: 
            self.parent_organ.update_settings(['setup',self.channel_no, 'alpha', self.mesh_type], self.alpha, 'mH')
        self.parent_organ.meshes[self.name]['alpha'] = self.alpha
    
    def get_alpha(self):
        return self.alpha
        
    def set_color(self, mesh_color):
        self.mesh.color(mesh_color)
        self.color = mesh_color
        #Update settings
        if self.channel_no != 'chNS': 
            self.parent_organ.update_settings(['setup','color_chs', self.channel_no, self.mesh_type], self.color, 'mH')
        else: 
            self.parent_organ.update_settings(['setup',self.channel_no,'alpha',self.mesh_type], self.color, 'mH')
        
        self.parent_organ.meshes[self.name]['color'] = self.color
        
    def get_color(self):
        return self.mesh_color   
        
    def get_mesh(self):
        try: 
            return self.mesh 
        except:
            self.load_mesh()
            return self.mesh
    
    def thickness_mesh(self, n_type, color_map = 'turbo', alpha=1):
        print('n_type: ', n_type)
        #n_type = thck(intTOext) 
        try: 
            mesh_out = self.mesh_meas['thck('+n_type+')']
            print('>> Extracting mesh from mesh_meas attribute')
            return mesh_out.alpha(alpha)
        except: 
            mesh_name = self.dirs['mesh']['thck('+n_type+')']
            dir_mesh = self.parent_organ.dir_res(dir='meshes') / mesh_name
            title = str(self.dirs['arrays']['thck('+n_type+')'])+'.npy'
            dir_npy = self.parent_organ.dir_res(dir='csv_all') / title
            print('dir_mesh:', dir_mesh, '\n','-dir_npy:', dir_npy)
            if dir_mesh.is_file() and dir_npy.is_file():
                name = 'Thickness'
                title = self.legend+'\n'+name+' [um]\n('+n_type.replace('TO','>')+')'
                print('fix this!!! ')
                if n_type == 'intTOext': 
                    short = 'th_i2e['+self.name.replace('_', '-')+']'
                else:
                    short = 'th_e2i['+self.name.replace('_', '-')+']'
                setup = self.parent_organ.mH_settings['wf_info']['heatmaps'][short]
                mesh_out = self.load_meas_mesh(dir_mesh, dir_npy, title, 
                                               setup=setup, alpha=alpha)
                return mesh_out.alpha(alpha)
            else: 
                print('>> Error: Unable to load mesh', self.name,'-',n_type)
                alert('error_beep')
                return None
        
    def balloon_mesh(self, n_type, color_map = 'turbo', alpha=1):
        #n_type = 'ballCL('+from_cl+'_'+from_cl_type+')'
    
        print('Mesh:', self.name, '-n_type: ', n_type)
        try: 
            mesh_out = self.mesh_meas['ballCL('+n_type+')']
            print('>> Extracting mesh from mesh_meas attribute')
            return mesh_out.alpha(alpha)
        except: 
            mesh_name = self.dirs['mesh']['ballCL('+n_type+')']
            dir_mesh = self.parent_organ.dir_res(dir='meshes') / mesh_name
            title = str(self.dirs['arrays']['ballCL('+n_type+')'])+'.npy'
            dir_npy = self.parent_organ.dir_res(dir='csv_all') / title
            print('dir_mesh:', dir_mesh, '\n','-dir_npy:', dir_npy)
            if dir_mesh.is_file() and dir_npy.is_file():
                name = 'Ballooning'
                title = self.legend+'\n'+name+' [um]\n('+n_type+')'
                cl_info = n_type.replace('_', '-')#split('(')[1].split(')')[0]
                short = 'ball['+self.name.replace('_', '-')+'(CL.'+cl_info+')]'
                #'ball[ch1-int(CL:ch1-int)]'
                setup = self.parent_organ.mH_settings['wf_info']['heatmaps'][short]
                mesh_out = self.load_meas_mesh(dir_mesh, dir_npy, title, 
                                               setup=setup, alpha=alpha)
                return mesh_out.alpha(alpha)
            else: 
                print('>> Error: Unable to load mesh', self.name,'-',n_type)
                alert('error_beep')
                return None
        
    def load_meas_mesh(self, dir_mesh, dir_npy, title, setup, alpha):
        title_print = title.replace('\n', ' ')
        print('>> Loading mesh '+title_print)
        mesh_out = vedo.load(str(dir_mesh))
        npy_colour = np.load(dir_npy)
        color_map = setup['colormap']
        
        # Assign colour
        mesh_out.pointdata['Distance'] = npy_colour
        if setup['default']: 
            vmin, vmax = np.min(npy_colour),np.max(npy_colour)
        else: 
            vmin = setup['min_val']
            vmax = setup['max_val']
        mesh_out.cmap(color_map)
        mesh_out.alpha(alpha)
        mesh_out.add_scalarbar(title=title, pos=(0.7, 0.05))
        mesh_out.mapper().SetScalarRange(vmin,vmax)
        mesh_out.legend(title)
        
        return mesh_out
        
    def get_centreline(self, nPoints=300, color='deepskyblue'): 
        try: 
            points = self.centreline_info['points']
            kspl = vedo.KSpline(points, res = nPoints).color(color).lw(5).legend('CL_'+self.name)
            return kspl
        except: 
            print('>> No centreline has been created for this mesh - ', self.name)
            return None
        
    def get_linLine(self, color='aqua'):
        cl_final =  self.get_centreline()
        cl_points = cl_final.points()
        cl_pt0 = cl_points[0]
        cl_ptm1 = cl_points[-1]
         
        #Create linear line
        linLine = vedo.Line(cl_pt0, cl_ptm1, c=color, lw=5)
        
        return linLine

    def get_clRibbon(self, nPoints, nRes, pl_normal, clRib_type, use_prev=False, 
                     ext_points = None, plot=True, oldV=False):
        """
        Function that creates dorso-ventral extended centreline ribbon
        """

        if oldV: 
            cl = self.get_centreline(nPoints)
            pts_cl = cl.points()
            
            # Extended centreline
            nn = -20
            inf_ext_normal = (pts_cl[nn]+(pts_cl[-1]-pts_cl[nn])*5)#*70
            outf_ext_normal = (pts_cl[0]+(pts_cl[0]-pts_cl[1])*100)#*70 (test for LnR cut Jun14.22)
        
            pts_cl_ext = np.insert(pts_cl,0,np.transpose(outf_ext_normal), axis=0)
            pts_cl_ext = np.insert(pts_cl_ext,len(pts_cl_ext),np.transpose(inf_ext_normal), axis=0)

            kspl_f = vedo.KSpline(pts_cl_ext, res=nRes).color('green').legend('Ksplo').lw(2)#601

        else: 
            if not use_prev: 
                cl = self.get_centreline(nPoints)
                pts_cl = cl.points()
                
                # Extended centreline
                nn = -2
                unit_inf = unit_vector(pts_cl[-1]-pts_cl[nn])
                inf_ext_normal = (pts_cl[-1]+(unit_inf)*100)#*70
                unit_outf = unit_vector(pts_cl[0]-pts_cl[1])
                outf_ext_normal = (pts_cl[0]+(unit_outf)*100)#*70 (test for LnR cut Jun14.22)
            
                pts_cl_ext = np.insert(pts_cl,0,np.transpose(outf_ext_normal), axis=0)
                pts_cl_ext = np.insert(pts_cl_ext,len(pts_cl_ext),np.transpose(inf_ext_normal), axis=0)

                kspl_o = vedo.KSpline(pts_cl_ext, res=nRes).color('green').legend('Ksplo').lw(2)#601
                kspl_f = self.modify_centreline(kspl_o=kspl_o, mesh=self.mesh)
            else: 
                kspl_o = []
                kspl_f = vedo.KSpline(ext_points, res=nRes).color('green').legend('Ksplo').lw(2)#601
        
        pts_cl_extf = kspl_f.points()
       
        # Increase the resolution of the extended centreline and interpolate to unify sampling
        xd = np.diff(pts_cl_extf[:,0])
        yd = np.diff(pts_cl_extf[:,1])
        zd = np.diff(pts_cl_extf[:,2])
        dist = np.sqrt(xd**2+yd**2+zd**2)
        u = np.cumsum(dist)
        u = np.hstack([[0],u])
        t = np.linspace(0, u[-1], nRes)#601
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

        kspl_ext = vedo.KSpline(resamp_pts, res=nRes).color('pink').legend('ExtendedCL').lw(2)#601
        
        pl_linLine_unitNormal = unit_vector(pl_normal)
        maj_bound =(max(self.parent_organ.get_maj_bounds())/2)*2
        pl_linLine_unitNormal120 = pl_linLine_unitNormal*maj_bound

        if clRib_type == 'ext2sides': # Names are switched but it works
            x_cl, y_cl, z_cl = pl_linLine_unitNormal120
            kspl_ext_D = kspl_ext.clone().x(x_cl).y(y_cl).z(z_cl).legend('kspl_CLExt1')
            kspl_ext_V = kspl_ext.clone().x(-x_cl).y(-y_cl).z(-z_cl).legend('kspl_CLExt2')
            cl_ribbon = vedo.Ribbon(kspl_ext_D, kspl_ext_V, alpha=0.2, res=(1500, 1500))
            cl_ribbon = cl_ribbon.wireframe(True).legend("rib_ExtCL(2-sides)")
    
        elif clRib_type == 'ext1side':
            x_ucl, y_ucl, z_ucl = pl_linLine_unitNormal*15
            cl_ribbon = []
            for i in range(10):
                kspl_ext_DA = kspl_ext.clone().x(i*x_ucl).y(i*y_ucl).z(i*z_ucl)
                kspl_ext_DB = kspl_ext.clone().x((i+1)*x_ucl).y((i+1)*y_ucl).z((i+1)*z_ucl)
                cl_ribbon2un = vedo.Ribbon(kspl_ext_DA, kspl_ext_DB, alpha=0.2, res=(220, 5))
                cl_ribbon.append(cl_ribbon2un)
            cl_ribbon = vedo.merge(cl_ribbon)
            cl_ribbon.legend('rib_ExtCL(1-side)').wireframe(True)

        else: 
            print('What? Which function is calling?')
            alert('bubble')

        if plot: 
            text = '> Final Extended Centreline and Ribbon.\n  Are you happy with the extended centreline/ribbon created?\n  Close window to continue'
            txt = vedo.Text2D(text, c=txt_color, font=txt_font, s=txt_size)

            vp = vedo.Plotter(N=1, axes=1)
            vp.add_icon(logo, pos=(0.9,1), size=0.25)
            vp.show(self.mesh, kspl_o, kspl_ext, cl_ribbon, txt, at=0, interactive=True)
 
        return cl_ribbon, kspl_ext
    
    def modify_centreline(self, kspl_o, mesh, ext=True):
        
        #Text
        if ext: 
            text = '>> Modify Extended Centreline Instructions: \n  -Drag extreme centreline points with mouse\n  -Remove them by selecting and pressing -Delete-\n  -Press q (lower q) when ready to proceed.\n  -Note: Only modify the orientation of the \n          extreme points but not the length of the centreline,\n          to make sure the created ribbon cuts the \n          whole stack into two different regions.'
        else: 
            text = '>> Modify Centreline Instructions: \n  -Drag extreme centreline points with mouse\n  -Remove them by selecting and pressing -Delete-\n  -Press q (lower q) when ready to proceed.'
        txt = vedo.Text2D(text, c=txt_color, font=txt_font, s=txt_size)

        #Make the user define the final points for the centreline
        vp = vedo.Plotter(N=1, axes=1)
        vp.add_icon(logo, pos=(0.9,1), size=0.25)
        vp.show(mesh, kspl_o, txt, interactive=False, at=0)
        # Add the spline tool using the same points and interact with it
        sptool = vp.add_spline_tool(kspl_o, closed=False)
        vp.interactive()
        # Switch off the tool
        sptool.off()
        # Extract and visualize the resulting spline
        sp = sptool.spline().lw(4)
        vp.close()

        text = '>> Extended Centreline is ready! \nClose window to continue'
        txt = vedo.Text2D(text, c=txt_color, font=txt_font, s=txt_size)

        vp2 = vedo.Plotter(N=1, axes=1)
        vp2.add_icon(logo, pos=(0.9,1), size=0.25)
        vp2.show(mesh, sp, txt, interactive=True, at=0)

        return sp

    def get_volume(self): 
        mesh_vol = self.mesh.volume()
        return mesh_vol
        
    def get_area(self): 
        mesh_area = self.mesh.area()
        return mesh_area
    
    def create_section(self, name, cut, color, alpha=0.05):
           
        sect_info = self.parent_organ.mH_settings['setup']['sect'][cut]
        if name == 'sect1':
            invert = True
        else: 
            invert = False
        print('name:', name, '- invert:', invert)
        
        submesh = SubMesh(parent_mesh = self, sub_mesh_type='Section', 
                          name = name, cut = cut, user_name = sect_info['name_sections'][name],
                          color = color, alpha = alpha)#,
        
        submesh.s3_invert = invert
        submesh.sub_user_name = sect_info['name_sections'][submesh.sub_name]
        self.parent_organ.add_submesh(submesh)
        
        return submesh
        
    def mask_segments(self, cut, check=False):
          
        # Get segments info
        no_discs = self.parent_organ.mH_settings['setup']['segm'][cut]['no_cuts_4segments']
        no_segm =  self.parent_organ.mH_settings['setup']['segm'][cut]['no_segments']

        # Mask im_channel
        im_ch = self.imChannel
        im_ch.load_chS3s([self.mesh_type])
        cont_tiss = getattr(im_ch, 's3_'+self.mesh_type)
        s3 = cont_tiss.s3()
        masked_s3 = s3.copy()
        
        for nn in range(no_discs):
            name_s3 = self.parent_organ.user_organName+'_mask_'+cut+'_DiscNo'+str(nn)+'.npy'
            s3_dir = self.parent_organ.dir_res(dir='s3_numpy') / name_s3
            s3_mask = np.load(str(s3_dir))
            s3_mask = s3_mask.astype('bool')
            masked_s3 = mask_disc(masked_s3, s3_mask)

        if 'NS' not in im_ch.channel_no: 
            max_depth = 1000
        else: 
            max_depth = 2000
        
        masked_mesh = create_submesh(masked_s3, self.resolution, keep_largest=False, 
                                     rotateZ_90=self.rotateZ_90)
        cut_masked = masked_mesh.split(maxdepth=max_depth)
        print('> Meshes making up tissue: ', len(cut_masked))
        alert('frog')
        if len(cut_masked) < no_segm: 
            if check: 
                obj = cut_masked
                plot_grid(obj=obj, txt=[], axes=5, sc_side=max(self.parent_organ.get_maj_bounds()))
                return None
            else: 
                rotate = True
        else: 
            rotate = True

        if rotate: 
            cut_masked_rot = []
            for n, mesh in enumerate(cut_masked):
                if self.rotateZ_90:
                    cut_masked_rot.append(mesh.rotate_z(-90).alpha(0.1).legend('No.'+str(n)))
                else:
                    cut_masked_rot.append(mesh.alpha(0.1).legend('No.'+str(n)))

            print('Rotated final segments (vedo.Mesh):', len(cut_masked_rot))
            
        return cut_masked_rot
        
    def create_segment(self, name, cut, color):
        
        segm_info = self.parent_organ.mH_settings['setup']['segm'][cut]
        alpha = self.alpha
        submesh = SubMesh(parent_mesh = self, sub_mesh_type='Segment', 
                          name = name, cut = cut, user_name = segm_info['name_segments'][name], 
                          color=color, alpha = alpha)
        
        submesh.sub_user_name = segm_info['name_segments'][submesh.sub_name]
        self.parent_organ.add_submesh(submesh)
            
        return submesh
    
class SubMesh():
    
    def __init__(self, parent_mesh, sub_mesh_type:str, name: str,
                 cut:str, user_name='', color='gold', alpha=0.05):
        
        self.parent_mesh = parent_mesh
        self.sub_name = name # ch_cont_segm/sect
        self.sub_mesh_type = sub_mesh_type # Section, Segment or Segment-Section
        self.keep_largest = False#keep_largest
        self.cut = cut

        if sub_mesh_type == 'Segment-Section': 
            parent_mesh = self.parent_mesh.parent_mesh #grandparent_mesh
            
        self.sub_name_all = cut+'_'+parent_mesh.name+'_'+name
        parent_organ = parent_mesh.parent_organ 
        
        if self.sub_name_all not in parent_organ.submeshes.keys():
            print('>> New submesh - ', self.sub_name_all)
            # new = True
            self.sub_legend = parent_mesh.legend + '_' + user_name # e.g. myoc_ext_atrium
            self.color = color
            self.alpha = alpha
            self.rotateZ_90 = parent_mesh.rotateZ_90
        else: 
            # new = False
            print('>> Recreating submesh - ', self.sub_name_all)
            #Get data from submesh dict
            submesh_dict = parent_organ.submeshes[self.sub_name_all]
            self.sub_legend = submesh_dict['sub_legend']
            self.color = submesh_dict['color']
            self.alpha = submesh_dict['alpha']
            self.rotateZ_90 = submesh_dict['rotateZ_90']
            for attr in ['s3_invert', 's3_mask_dir', 'dict_segm', 'sub_user_name']:
                if attr in submesh_dict.keys():
                    value = submesh_dict[attr]
                    setattr(self, attr, value)
            
        self.imChannel = parent_mesh.imChannel
        self.mesh_type = parent_mesh.mesh_type
        self.resolution = parent_mesh.resolution
                    
    def get_sect_mesh(self, output='mesh'):
        
        print('>>>> get_sect_mesh: ',self.sub_name_all)
        if self.sub_mesh_type != 'Segment-Section': 
            cut = self.cut
        else: 
            cut = self.cut.split('_')[1:][0]

        mask_name = self.imChannel.parent_organ.mH_settings['wf_info']['sections'][cut.title()]['mask_name']#getattr(self.parent_mesh.parent_organ, 'mask_sect_'+cut.lower())
        mask_dir = self.imChannel.parent_organ.dir_res(dir ='s3_numpy') / mask_name

        s3_mask = np.load(str(mask_dir))
        s3_mask = s3_mask.astype('bool')
        if self.s3_invert: 
            maskF = np.invert(s3_mask)
        else: 
            maskF = s3_mask
            
        im_ch = self.imChannel      
        cont = self.mesh_type
        im_ch.load_chS3s([cont])
        cont_tiss = getattr(im_ch, 's3_'+cont)
        s3 = cont_tiss.s3()
        masked_s3 = s3.copy()
    
        masked_s3[maskF] = 0
        if output == 'mask':
            return masked_s3
        elif output == 'mesh':
            mesh = create_submesh(masked_s3, self.resolution, self.keep_largest, self.rotateZ_90)
            mesh.legend(self.sub_legend).wireframe()
            mesh.alpha(self.alpha)
            mesh.color(self.color)
            return mesh
        else: 
            print('get_sect_mesh what?')
            alert('error_beep')
            return None

    def get_segm_mesh(self):
        
        print('>>>> get_segm_mesh: ',self.sub_name_all)
        cut, ch, cont, segm = self.sub_name_all.split('_')
        directory = self.parent_mesh.parent_organ.dir_res(dir = 'meshes')
        mesh_dir = directory / self.parent_mesh.parent_organ.mH_settings['wf_info']['segments']['setup'][cut]['dirs'][ch][cont][segm]
        segm_mesh = vedo.load(str(mesh_dir))
        print('segm_mesh:', type(segm_mesh))
        segm_mesh.color(self.parent_mesh.parent_organ.mH_settings['setup']['segm'][cut]['colors'][segm])
        segm_mesh.legend(self.sub_legend).wireframe().alpha(self.alpha)
            
        return segm_mesh
    
    def get_sect_segm_mesh(self, seg_cut): 
        
        #SECTION/REGION
        segm, sect = self.sub_name.split('_')
        #Get the section mask containing only this section
        sect_mask = self.get_sect_mesh(output='mask')
        #SEGMENT
        #Mask the segments providing as input the masked section
        cut_masked = self.mask_segments(cut = seg_cut, s3 = sect_mask)
        #Create a dictionary containing the information of the classified segments 
        #Heree!!!
        dict_segm = self.imChannel.parent_organ.dict_segments(seg_cut, other=False)
        try: 
            ch, cont = self.parent_mesh.name.split('_')
        except: 
            ch, cont = self.parent_mesh.parent_mesh.name.split('_')
        method = self.imChannel.parent_organ.mH_settings['wf_info']['segments']['setup'][seg_cut]['ch_info'][ch][cont]
        if method in ['ext-ext', 'cut_with_ext-ext', 'cut_with_other_ext-ext']: 
            try:
                ext_subsgm = self.imChannel.parent_organ.ext_subsgm
                print('try ext_subsgm')
            except: 
                if 's' in seg_cut: 
                    seg_cut = seg_cut[1:]
                ext_subsgm = self.imChannel.parent_organ.get_ext_subsgm(seg_cut)
            print('except ext_subsgm')
            print('ext_subsgm: ',ext_subsgm)

            #Classify the resulting segments using ext mesh
            sp_dict_segm = classify_segments_from_ext(meshes = cut_masked, 
                                                        dict_segm = dict_segm[segm],
                                                        ext_sub = ext_subsgm[segm])
            print('dict_segm after classif: ', sp_dict_segm)
            segm_sect_mesh = create_subsegment(self.imChannel.parent_organ, self, seg_cut, cut_masked, 
                                                    'segm-sect', sp_dict_segm, self.color)
        else: 
            segm_sect_mesh = None
            print('This needs to be coded!')

        return segm_sect_mesh
    
    def create_segm_sect(self, segm_sect, cuts, color, alpha=0.05): 
        seg_cut, reg_cut = cuts.split('_')
        seg_name, reg_name = segm_sect.split('_')
        segm_info = self.parent_mesh.parent_organ.mH_settings['setup']['segm'][seg_cut[1:]]
        sect_info = self.parent_mesh.parent_organ.mH_settings['setup']['sect'][reg_cut]
        if reg_name == 'sect1':
            invert = True
        else: 
            invert = False
        print('seg_name:', seg_name, 'reg_name:', reg_name, '- invert:', invert)

        user_name = segm_info['name_segments'][seg_name]+'_'+sect_info['name_sections'][reg_name]
        print('cuts:', cuts, 'user_name:', user_name)
        submesh = SubMesh(parent_mesh = self, sub_mesh_type='Segment-Section', 
                          name = segm_sect, cut = cuts, user_name = user_name,
                          color = color, alpha = alpha)#,

        submesh.s3_invert = invert
        submesh.sub_user_name = segm_info['name_segments'][seg_name]+'_'+sect_info['name_sections'][reg_name]
        submesh.imChannel.parent_organ.add_submesh(submesh)

        return submesh
    
    def mask_segments(self, cut, s3):
          
        # Get segments info
        no_discs = self.imChannel.parent_organ.mH_settings['setup']['segm'][cut]['no_cuts_4segments']
        no_segm =  self.imChannel.parent_organ.mH_settings['setup']['segm'][cut]['no_segments']

        # Mask im_channel
        im_ch = self.imChannel
        masked_s3 = s3.copy()
        
        for nn in range(no_discs):
            name_s3 = self.imChannel.parent_organ.user_organName+'_mask_'+cut+'_DiscNo'+str(nn)+'.npy'
            s3_dir = self.imChannel.parent_organ.dir_res(dir='s3_numpy') / name_s3
            s3_mask = np.load(str(s3_dir))
            s3_mask = s3_mask.astype('bool')
            masked_s3 = mask_disc(masked_s3, s3_mask)

        if 'NS' not in im_ch.channel_no: 
            max_depth = 1000
        else: 
            max_depth = 2000
        
        masked_mesh = create_submesh(masked_s3, self.resolution, keep_largest=False, 
                                     rotateZ_90=self.rotateZ_90)
        cut_masked = masked_mesh.split(maxdepth=max_depth)
        print('> Meshes making up tissue: ', len(cut_masked))
        alert('frog')
        if len(cut_masked) < no_segm: 
            obj = cut_masked
            plot_grid(obj=obj, txt=[], axes=5, sc_side=max(self.imChannel.parent_organ.get_maj_bounds()))
            cut_masked_rot = cut_masked
        else: 
            cut_masked_rot = []
            for n, mesh in enumerate(cut_masked):
                if self.rotateZ_90:
                    cut_masked_rot.append(mesh.rotate_z(-90).alpha(0.1).legend('No.'+str(n)))
                else:
                    cut_masked_rot.append(mesh.alpha(0.1).legend('No.'+str(n)))

        print('Rotated final segments (vedo.Mesh):', len(cut_masked_rot))
        return cut_masked_rot
    
    def get_segm_sect_mesh(self): 
        pass

    def set_alpha(self, mesh_alpha):      
        self.alpha = mesh_alpha
        #Update settings
        self.parent_mesh.parent_organ.submeshes[self.sub_name_all]['alpha'] = self.alpha
        # self.parent_mesh.parent_organ.save_organ()
    
    def get_alpha(self):
        return self.parent_mesh.parent_organ.submeshes[self.sub_name_all]['alpha']
        
    def set_color(self, mesh_color):
        self.color = mesh_color
        #Update settings
        self.parent_mesh.parent_organ.submeshes[self.sub_name_all]['color'] = self.color
        # self.parent_mesh.parent_organ.save_organ()
        
    def get_color(self):
        return self.parent_mesh.parent_organ.submeshes[self.sub_name_all]['color'] 

class ImChannelMC():
    def __init__(self, organ:Organ, ch_name:str):#
        
        self.parent_organ = organ
        self.parent_organ_name = organ.user_organName
        self.channel_no = ch_name
        self.user_chName = organ.mC_settings['setup']['name_chs'][ch_name]
        if organ.mC_settings['setup']['mH_channel'][ch_name] != False: 
            self.mH_channel = organ.mC_settings['setup']['mH_channel'][ch_name].split(' (')[0]
        else: 
            self.mH_channel = False

        print('Create/Load a mC Channel')
        if self.channel_no not in organ.imChannelsMC.keys():   
            print('>> New Channel-', self.channel_no)
            self.new_imChannelMC()
        else: 
            print('>> Loading Channel-', self.channel_no)
            self.load_imChannelMC()
    
    def new_imChannelMC(self):
        organ = self.parent_organ
        ch_name = self.channel_no

        self.resolution = organ.info['resolution']
        self.dir_cho = organ.img_dirs[ch_name]['image']['dir']       
        self.shape = self.im().shape
        self.process = ['Init']
    
        #organ.mC_settings['setup']['alpha'][ch_name] = {}
        if self.mH_channel == False: 
            self.save_channel(im_proc=self.im())
        else: 
            #Connect to ImChannel
            self.imChannel = self.parent_organ.obj_imChannels[self.mH_channel]
        organ.add_channel_mC(imChannelMC=self)

    def load_imChannelMC(self): 
        organ = self.parent_organ
        ch_name = self.channel_no
        
        self.resolution = organ.imChannelsMC[ch_name]['resolution']
        self.dir_cho = Path(organ.imChannelsMC[ch_name]['dir_cho'])
        self.shape = tuple(organ.imChannelsMC[ch_name]['shape'])
        self.process = organ.imChannelsMC[ch_name]['process']

    def im(self):#
        if hasattr(self, 'dir_stckproc'):
            dir_stck = self.parent_organ.dir_res(dir='s3_numpy') / self.dir_stckproc
            im = np.load(str(dir_stck))
        else: 
            im = io.imread(str(self.dir_cho))

        if not isinstance(im, np.ndarray):
            print('>> Error: morphoCell was unable to load tiff.\n>> Directory: ',str(self.dir_cho))
            alert('error_beep')
        else: 
            return im
        
    def save_channel(self, im_proc):
        organ_name = self.parent_organ.user_organName
        im_name = organ_name + '_StckProc_' + self.channel_no + '.npy'
        im_dir = self.parent_organ.dir_res(dir='s3_numpy') / im_name
        np.save(im_dir, im_proc)
        if not im_dir.is_file():
            print('>> Error: Processed channel was not saved correctly!\n>> File: '+im_name)
            alert('error_beep')
        else: 
            print('>> Processed channel saved correctly! - ', im_name)
            alert('countdown')
            self.dir_stckproc = im_name
        
class Cells():
    def __init__(self, organ:Organ):

        self.parent_organ = organ
        self.parent_organ_name = organ.user_organName
        self.channel_no = 'chA'
        self.user_chName = organ.mC_settings['setup']['name_chs']['chA']

        if self.channel_no not in organ.cellsMC.keys():  
            print('>> New Cells Channel-', self.channel_no)
            self.new_Cells()
        else: 
            print('>> Loading Channel-', self.channel_no)
            self.load_Cells()
            
    def new_Cells(self): 
        organ = self.parent_organ
        ch_name = self.channel_no

        self.resolution = organ.info['resolution']
        self.dir_cho = organ.img_dirs[ch_name]['data']['dir'] 
        self.dir_img = organ.img_dirs[ch_name]['image']['dir'] 
        self.shape =  io.imread(str(self.dir_img)).shape
        ext = Path(self.dir_cho).suffix

        #Check the extension of the file
        if ext == '.xlsx' or ext == '.xls':
            cells_position = pd.read_excel(self.dir_cho, header = 3, 
                            usecols = ['Position X','Position Y','Position Z', 'ID'], index_col=3) 
        else: #'.csv'
            cells_position = pd.read_csv(self.dir_cho, header = 2, 
                                         usecols = ['Position X','Position Y','Position Z', 'ID'], index_col=3)
        
        print('cells_position:', cells_position)
        sphs_pos = cells_position.copy()
        sphs_pos['Position Y'] = -cells_position['Position Y']

        self.set_cells(sphs_pos, init=True)
        self.save_cells(cells = sphs_pos, init=True)
        organ.add_cells_mC(cells=self)
    
    def load_Cells(self): 

        organ = self.parent_organ
        ch_name = self.channel_no

        self.resolution = organ.cellsMC[ch_name]['resolution']
        self.dir_cho = organ.cellsMC[ch_name]['dir_cho']
        try: 
            self.dir_img = organ.cellsMC[ch_name]['dir_img']
        except: 
            self.dir_img = organ.img_dirs['chA']['image']
        self.dir_cells = organ.cellsMC[ch_name]['dir_cells']
        self.shape = organ.cellsMC[ch_name]['shape']
        
        cells = self.df_cells()
        self.set_cells(cells)
    
    def df_cells(self): 

        cells_dir = self.parent_organ.dir_res(dir='s3_numpy') / self.dir_cells
        self.cells_position = pd.read_csv(cells_dir, index_col=0)

        return self.cells_position

    def set_cells(self, sphs_pos, init=False): 

        cell_ids = list(sphs_pos.index)
        try:
            deleted = sphs_pos['deleted']
        except: 
            deleted = ['NO']*len(sphs_pos)
        pos_x = sphs_pos['Position X']
        pos_y = sphs_pos['Position Y']
        pos_z = sphs_pos['Position Z']

        cols = range(len(sphs_pos))
        color = self.parent_organ.mC_settings['setup']['color_chs']['chA']
        
        sphs = []
        for i, x, y, z, id, bool_del in zip(count(), pos_x, pos_y, pos_z, cell_ids, deleted):
            if bool_del == 'NO': 
                pos = (x, y, z)
                s = vedo.Sphere(r=2).pos(pos).color(color)#cols[i])
                s.name = f"Cell Nr.{id}"
                sphs.append(s)
            else: 
                print('sph deleted no:', id)
        
        self.cells = sphs

    def remove_cells(self): 

        vol_iso, vol_settings = self.parent_organ.organ_vol_iso()

        sphs = self.cells
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
                    evt.actor.color('black')
                    new_no = 1000000+cell_no
                    evt.actor.name = f"Cell Nr.{new_no}"
                else: # add back
                    ind2rem = cells2remove.index(cell_no-1000000)
                    cells2remove.pop(ind2rem)
                    old_no = cell_no-1000000
                    evt.actor.color('gold')
                    evt.actor.name = f"Cell Nr.{old_no}"
                    
                plt.remove(silcont.pop()).add(sil)
                silcont.append(sil)

        if len(vol_iso) >= 1: 
            def sliderAlphaMeshOut(widget, event):
                valueAlpha = widget.GetRepresentation().GetValue()
                vol_iso[0].alpha(valueAlpha)
        if len(vol_iso) >= 2: 
            def sliderAlphaMeshOut2(widget, event):
                valueAlpha = widget.GetRepresentation().GetValue()
                vol_iso[1].alpha(valueAlpha)
        if len(vol_iso) >= 3: 
            def sliderAlphaMeshOut3(widget, event):
                valueAlpha = widget.GetRepresentation().GetValue()
                vol_iso[2].alpha(valueAlpha)
        
        txt_slider_size2 = 0.7
        
        msg = vedo.Text2D("", pos="bottom-center", c=txt_color, font=txt_font, s=txt_size, bg='red', alpha=0.2)
        txtf = self.parent_organ.user_organName+' - Removing cells \n -Click in the cells you would like to remove from the analysis, and close the window when done.' 
        txt = vedo.Text2D(txtf,  c=txt_color, font=txt_font, s=txt_size)
        plt = vedo.Plotter(axes=1)
        plt.add_icon(logo, pos=(0.1,0.0), size=0.25)
        if len(vol_iso) >= 1: 
            plt.addSlider2D(sliderAlphaMeshOut, xmin=0, xmax=0.99, value=0.05,
                pos=[(0.92,0.25), (0.92,0.35)], c= vol_settings['color'][0], 
                title='Opacity\n'+ vol_settings['name'][0].title(), title_size=txt_slider_size2)
        if len(vol_iso) >=2:
            plt.addSlider2D(sliderAlphaMeshOut2, xmin=0, xmax=0.99, value=0.05,
                pos=[(0.92,0.40), (0.92,0.50)], c=vol_settings['color'][1], 
                title='Opacity\n'+ vol_settings['name'][1].title(), title_size=txt_slider_size2)
        if len(vol_iso) >=3:
            plt.addSlider2D(sliderAlphaMeshOut3, xmin=0, xmax=0.99, value=0.05,
                pos=[(0.72,0.25), (0.72,0.35)], c=vol_settings['color'][2],
                title='Opacity\n'+ vol_settings['name'][2].title(), title_size=txt_slider_size2)
            
        plt.addCallback('mouse click', remove_cells)
        plt.show(sphs, vol_iso, msg, txt, zoom=0.8)
        
        # Read file
        cells_position = self.df_cells()
        deleted = cells_position['deleted']

        for cell in cells2remove: 
            deleted.at[cell] = 'YES'
        
        cells_position['deleted'] = deleted
        self.save_cells(cells_position)
        self.set_cells(cells_position)

    def get_colour_class(self, cut, mtype): 

        if mtype == 'segm': 
            col_name = 'Segment-'+cut
            type_name = 'segm_mC'

        cells_position = self.parent_organ.cellsMC['chA'].df_cells()
        segm_colors = self.parent_organ.mC_settings['setup'][type_name][cut]['colors']
        segm_class = cells_position[col_name]
        color_class = np.empty(len(cells_position), dtype='object')

        for n, val in enumerate(segm_class): 
            if isinstance(val, str):
                segm, _  = val.split(':')
                color_class[n] = segm_colors[segm]

        return color_class

    def colour_cells(self, sphs_pos, color_class):

        cell_ids = list(sphs_pos.index)
        deleted = sphs_pos['deleted']
        sphs_pos = sphs_pos[["Position X", "Position Y", "Position Z"]]
        sphs_pos_tuple = list(sphs_pos.itertuples(index=False, name=None))

        sphs = []
        for i, pos, id, bool_del, color in zip(count(), sphs_pos_tuple, cell_ids, deleted, color_class):
            if bool_del == 'NO': 
                s = vedo.Sphere(r=2).pos(pos).color(color)
                s.name = f"Cell Nr.{id}"
                sphs.append(s)
        
        self.cells = sphs

        return sphs

    def save_cells(self, cells, init=False): 
        #Add a deleted row
        if init:
            n_cells = len(cells)
            deleted = ['NO']*n_cells
            cells['deleted'] = deleted

        organ_name = self.parent_organ.user_organName
        cells_name = organ_name + '_CellsProc_' + self.channel_no + '.csv'
        cells_dir = self.parent_organ.dir_res(dir='s3_numpy') / cells_name
        cells.to_csv(cells_dir, index=True) 

        if not cells_dir.is_file():
            print('>> Error: Cell positions file was not saved correctly!\n>> File: '+cells_name)
            alert('error_beep')
        else: 
            self.dir_cells = cells_name
            print('>> Processed channel saved correctly! - ', cells_name)
            alert('countdown')
    
    def assign_class(self, cells_position, sp_class, col_name,):

        deleted = cells_position['deleted']
        sp_classf = []
        for i, cl, bool_del in zip(count(), sp_class, deleted): 
            if bool_del == 'YES':
                sp_classf.append('NA')
            else: 
                sp_classf.append(cl)
        
        cells_position[col_name] = sp_classf

        return cells_position
        
class TempProj(): 
    def __init__(self, proj_dict:dict):

        self.info = proj_dict['info']
        self.analysis = proj_dict['analysis']
        self.mH_methods = proj_dict['mH_methods']
        self.mH_settings = proj_dict['mH_settings']
        self.mH_channels = proj_dict['mH_channels']
        self.mH_segments = proj_dict['mH_segments']
        self.mH_sections = proj_dict['mH_sections']

        self.mC_methods = proj_dict['mC_methods']
        self.mC_settings = proj_dict['mC_settings']
        self.mC_channels = proj_dict['mC_channels']
        self.mC_segments = proj_dict['mC_segments']
        self.mC_sections = proj_dict['mC_sections']

        self.workflow = proj_dict['workflow']
        self.organs = proj_dict['organs']
        self.gui_custom_data = proj_dict['gui_custom_data']

class MyFaceSelectingPlotter(vedo.Plotter):
    def __init__(self, colors, color_o, views, **kwargs):
        
        # Create planar_views dictionary
        planar_views = {}
        for n, view, color in zip(count(), views, colors): 
            planar_views[view] = {'color': color}
            
        self.planar_views = planar_views
        self.views = views
        self.color_o = color_o
        self.done = False
        
        # Create message that displays instructions
        self.msg = vedo.Text2D("", pos="top-center", c=txt_color, bg='white', font=txt_font, alpha=0.8, s=0.7)
        self.msg_face = vedo.Text2D("", pos="bottom-center", c=txt_color, bg='red', font=txt_font, alpha=0.2, s=0.7)
        
        # Initialise plotter with current planar view
        self.active_n = 0
        self.selected_faces = []
        self.active_color = self.planar_views[self.current_view()]['color']
        self.get_msg()
        vedo.printc('Selecting '+self.current_view().upper()+'...', c="g", invert=True)

        #Initialise Plotter 
        super().__init__(**kwargs)
    
    def current_view(self):
        return self.views[self.active_n]
        
    def get_msg(self):
        if self.active_n < len(self.views)-1:
            msg1 = 'Instructions: Select (click) the cube face that represents the '+self.current_view().upper()+' face'
            msg2 = '\n press -c- and then rotate the cube to continue.'
            self.msg.text(msg1+msg2)
        else:
            if self.check_full():
                msg_close = 'Instructions: You are done selecting planar views. Close the window to continue.'
                self.msg.text(msg_close)
            else: 
                msg1 = 'Instructions: Select (click) the cube face that represents the '+self.current_view().upper()+' face'
                msg2 = '\n  press -c- and then rotate the cube to continue.'
                self.msg.text(msg1+msg2)
    
    def check_full(self):
        check = []
        for pv in self.planar_views:
            if 'pl_normal' in self.planar_views[pv].keys():
                check.append(isinstance(self.planar_views[pv]['pl_normal'], np.ndarray))
            else: 
                check.append(False)

        return all(check)
        
    def on_key_press(self, evt):
        if not self.done: 
            if evt.keypress == "c":
                planar_view = self.current_view()
                if 'idcell' in self.planar_views[planar_view].keys():
                    self.selected_faces.append(self.planar_views[planar_view]['idcell'])
                    self.get_msg()
                    vedo.printc('>> n:'+str(self.active_n)+'-'+str(len(self.views)), c='orange', invert=True)
                    if self.active_n < len(self.views)-1:
                        self.active_n += 1
                        planar_view = self.current_view()
                        self.active_color = self.planar_views[planar_view]['color']
                        vedo.printc('Now selecting '+planar_view.upper()+'...', c="g", invert=True)
                        self.get_msg()
                    else: 
                        self.get_msg()
                        if self.check_full(): 
                            self.done = True
                            vedo.printc('You are done, now close the window!', c='orange', invert=True)
                    vedo.printc('n:'+str(self.active_n), c='orange', invert=True)
                else: 
                    msg_warning = "You need to select a cube's face for the "+planar_view.upper()+" face to continue."
                    self.msg_face.text(msg_warning)
                    vedo.printc('No cell has been selected',c='r', invert=True)
        else: 
            vedo.printc('You are done, now close the window!', c='orange', invert=True)
        
    def select_cube_face(self, evt):
        if not self.done: 
            if isinstance(evt.actor, vedo.shapes.Cube):
                orient_cube = evt.actor
                if not orient_cube:
                    return
                pt = evt.picked3d
                idcell = orient_cube.closest_point(pt, return_cell_id=True)
                vedo.printc('You clicked (idcell):', idcell, c='y', invert=True)
                if set(orient_cube.cellcolors[idcell]) == set(self.color_o):
                    orient_cube.cellcolors[idcell] = self.active_color #RGBA 
                    for cell_no in range(len(orient_cube.cells())):
                        if cell_no != idcell and cell_no not in self.selected_faces: 
                            orient_cube.cellcolors[cell_no] = self.color_o #RGBA 
                planar_view =  self.current_view()
                self.msg_face.text("You selected cube's face number "+str(idcell)+" as the "+planar_view.upper()+" face")
                self.planar_views[planar_view]['idcell'] = idcell
                cells = orient_cube.cells()[idcell]
                points = [orient_cube.points()[cell] for cell in cells]
                plane_fit = vedo.fit_plane(points, signed=True)
                # print(plane_fit.normal, type(plane_fit.normal))
                self.planar_views[planar_view]['pl_normal'] = plane_fit.normal
                self.planar_views[planar_view]['pl_centre'] = orient_cube.cell_centers()[idcell]
        else: 
            vedo.printc('You are done, now close the window!', c='orange', invert=True)
        
#%% - Drawing Functions
#%% func - draw_line
def draw_line (clicks, myIm, color_draw):
    """
    Function that draws white or black line connecting all the clicks received as input
    """
    for num, click in enumerate(clicks):
        if num < len(clicks)-1:
            pt1x, pt1y = click
            pt2x, pt2y = clicks[num+1]
            rr, cc, val = line_aa(int(pt1x), int(pt1y),
                                  int(pt2x), int(pt2y))
            rr1, cc1, val1 = line_aa(int(pt1x)+1, int(pt1y),
                                     int(pt2x)+1, int(pt2y))
            rr2, cc2, val2 = line_aa(int(pt1x)-1, int(pt1y),
                                     int(pt2x)-1, int(pt2y))
            if color_draw == "white" or color_draw == "":
                myIm[rr, cc] = val * 50000
            elif color_draw == "1":
                myIm[rr, cc] = 1
                myIm[rr1, cc1] = 1
                myIm[rr2, cc2] = 1
            elif color_draw == "0":
                myIm[rr, cc] = 0
                myIm[rr1, cc1] = 0
                myIm[rr2, cc2] = 0
            else: #"black"
                myIm[rr, cc] = val * 0
                myIm[rr1, cc1] = val1 * 0
                
    return myIm

#%% - Masking
#%% func - mask_disc
def mask_disc(s3, s3_cyl):
    
    #Load stack shape
    s3_mask = copy.deepcopy(s3)
    zdim = s3_mask.shape[2]
    
    for slc in range(zdim):
        im_cyl =s3_cyl[:,:,slc]
        pos_pts = np.where(im_cyl == 1)
        clicks = [(pos_pts[0][i], pos_pts[1][i]) for i in range(pos_pts[0].shape[0])]
        if len(clicks+clicks) > 200:
            clicks_random = random.sample(clicks+clicks, 200)#2*len(clicks))
        else:
            clicks_random = random.sample(clicks+clicks, 2*len(clicks))
            
        im = s3_mask[:,:,slc]
        myIm = draw_line(clicks_random, im, '0')
        s3_mask[:,:,slc] = myIm

    return s3_mask

#%% - Mesh functions
#%% func - create_mesh
def create_submesh(masked_s3, resolution, keep_largest:bool, rotateZ_90:bool):
    
    verts, faces, _, _ = measure.marching_cubes(masked_s3, spacing=resolution, method='lewiner')
    # Create meshes
    mesh = vedo.Mesh([verts, faces])
    if keep_largest:
        mesh = mesh.extract_largest_region()
    if rotateZ_90:
        mesh.rotate_z(-90)
    alert('woohoo')
    
    return mesh
 
#%% Module loaded
print('morphoHeart! - Loaded mH_classes')