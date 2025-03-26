'''
morphoHeart_funcBasics

@author: Juliana Sanchez-Posada
'''
#%% Imports - ########################################################
import os
from pathlib import Path, PurePath
from playsound import playsound
import numpy as np 
import flatdict
from functools import reduce  
import operator
from itertools import count
import vedo as vedo
import copy
import pandas as pd
import seaborn as sns

path_fcBasics = os.path.abspath(__file__)

#%% morphoHeart Imports - ##################################################
from ..gui.config import mH_config

alert_all=True
heart_default=False
dict_gui = {'alert_all': alert_all,
            'heart_default': heart_default}

#%% - morphoHeart Basic General Functions 
def alert(sound:str):
    '''
    bubble:
    clown:
    connection:
    countdown:
    error:
    error_beep:
    frog: 
    jump:
    whistle:
    woohoo:
    '''
    # print('alert: ', mH_config.gui_sound)
    mp3_basic = ['connection','woohoo', 'error_beep', 'error']
    
    if mH_config.gui_sound[0]: 
        if mH_config.gui_sound[1] == 'Minimal' and sound in mp3_basic:
            sound_on = True
        elif mH_config.gui_sound[1] == 'Minimal' and sound not in mp3_basic:
            sound_on = False
        else: 
            sound_on = True
    else: 
        sound_on = False
    try: 
        if sound_on: 
            path_parentSounds = Path(path_fcBasics).parent.parent.parent
            sound_mp3= sound +'.mp3'
            path = path_parentSounds / 'sounds' / sound_mp3
            playsound(str(path))
            if sound == 'error_beep': 
                print('ERRORRR!')
    except:# playsound.PlaysoundException: 
        print('Sound not played: ', sound)
        pass

def input_range(response):
    try: 
        obj_num = []
        comma_split = response.split(',')
        print(comma_split)
        for string in comma_split:
            if '-' in string:
                minus_split = string.split('-')
                print(minus_split)
                for n in list(range(int(minus_split[0]),int(minus_split[1])+1,1)):
                    obj_num.append(n)
                    print('Appended: ', str(n))
            else:
                obj_num.append(int(string))
                print('Appended2: ', string)
        return obj_num
    except: 
        obj_num = 'error'

# Dataframes
def df_reset_index(df:pd.DataFrame, mult_index:list): 

    df_new = df.reset_index()
    # print(df_new.head(10))
    df_new = df_new.set_index(mult_index)
    # print(df_new.head(10), '\n', df_new.index)

    return df_new

def df_add_value(df:pd.DataFrame, index:tuple, value):

    df.loc[index, 'Value'] = value
    print('New value: ',df.loc[index, 'Value'])

    return df

# Dictionaries
def compare_dicts(dict_1, dict_2, dict_1_name, dict_2_name, path="", ignore_dir=False):
    """Compare two dictionaries recursively to find non mathcing elements
    https://stackoverflow.com/questions/27265939/comparing-python-dictionaries-and-nested-dictionaries

    Args:
        dict_1: dictionary 1
        dict_2: dictionary 2

    Returns:

    """
    from .mH_classes_new import Project, Organ, ImChannel, ImChannelNS, ContStack, Mesh_mH
    
    err = ''
    key_err = ''
    value_err = ''
    old_path = path
    for k in dict_1.keys():
        path = old_path + "[%s]" % k
        if not k in dict_2:
            key_err += "\tKey %s%s not in %s\n" % (dict_1_name, path, dict_2_name)
        else:
            if isinstance(dict_1[k], dict) and isinstance(dict_2[k], dict):
                err += compare_dicts(dict_1[k],dict_2[k], dict_1_name, dict_2_name, path)#'d1','d2', path)
            else:
                if 'dir' in k: 
                    if str(dict_1[k]) != str(dict_2[k]):
                        value_err += "Value of %s%s (%s) not same as \n\t %s%s (%s)\n"\
                            % (dict_1_name, path, dict_1[k], dict_2_name, path, dict_2[k])
                            
                elif isinstance(dict_1[k], np.ndarray) or isinstance(dict_2[k], np.ndarray):
                    # print(k)
                    if not isinstance(dict_1[k], np.ndarray):
                        dd1 = np.array(dict_1[k])
                    else: 
                        dd1 = dict_1[k]
                        
                    if not isinstance(dict_2[k], np.ndarray):
                        dd2 = np.array(dict_2[k])
                    else: 
                        dd2 = dict_2[k]

                    comparison = dd1 == dd2
                    equal_arrays = comparison.all()
                    if not equal_arrays: 
                        value_err += "Value of %s%s (%s) not same as \n\t %s%s (%s)\n"\
                            % (dict_1_name, path, dict_1[k], dict_2_name, path, dict_2[k])
                            
                else: 
                    if not isinstance(dict_1[k], (Project, Organ, ImChannel, ImChannelNS, ContStack, Mesh_mH, vedo.Mesh)):
                        # print(k)
                        if dict_1[k] != dict_2[k]:
                            value_err += "Value of %s%s (%s) not same as \n\t %s%s (%s)\n"\
                                % (dict_1_name, path, dict_1[k], dict_2_name, path, dict_2[k])
                    # else: 
                    #     print(type(dict_1[k]))


    for k in dict_2.keys():
        path = old_path + "[%s]" % k
        if not k in dict_1:
            key_err += "\tKey %s%s not in %s\n" % (dict_2_name, path, dict_1_name)

    res = key_err + value_err + err
    return res

def compare_nested_dicts(dict_1, dict_2, dict_1_name, dict_2_name, path=""):
    
    res = compare_dicts(dict_1, dict_2, dict_1_name, dict_2_name, path="")
    if len(res) ==0:
        return '\tNo differences!'
    else: 
        return res
    
def update_gui_set(loaded:dict, current:dict): 
    flat_loaded = flatdict.FlatDict(loaded)
    flat_current = flatdict.FlatDict(current)
    current_keys = flat_current.keys()
    loaded_keys = flat_loaded.keys()

    changed = False
    final_dict = copy.deepcopy(loaded)
    for key in current_keys: 
        # print('>>', key)
        value_current = get_by_path(current, key.split(':'))
        print('current:', value_current)
        try: #The value already existed in the loaded dictionary
            value_loaded = get_by_path(loaded, key.split(':'))
            print('value_loaded:',value_loaded)
            if value_current != value_loaded:
                print('value_current != value_loaded:', value_current, value_loaded, '- key: ', key)
                if value_current == {} and isinstance(value_loaded, dict): 
                    print('>> NOT Changed')
                    pass
                else: 
                    set_by_path(final_dict, key.split(':'),value_current)
                    changed = True
        except: #The value doesn't exist in the loaded dictionary
            print('add key to loaded!')
            changed = True
            try: 
                set_by_path(final_dict, key.split(':'), value_current)
            except KeyError as e: 
                key_error = e.args[0]
                #Get the position of that key in the flat key
                split_key = key.split(':')
                len_all_keys = len(split_key)
                index = split_key.index(key_error)
                #Get the length of the keys that need to be added
                len_key = len(split_key[index:])
                for num in range(len_key):
                    # print(split_key, index+num)
                    key2add = split_key[:index+num+1]
                    # print('key2add:', key2add)
                    if num != len_key-1:
                        set_by_path(final_dict, key2add, {})
                    else:
                        set_by_path(final_dict, key2add, value_current)
    
    print('GUI:', current)
    print('Loaded_o: ',loaded)
    print('Loaded_f: ',final_dict)

    return final_dict, changed
        
def get_by_path(root_dict, items):
    """Access a nested object in root_dict by item sequence.
    by Martijn Pieters (https://stackoverflow.com/questions/14692690/access-nested-dictionary-items-via-a-list-of-keys)
    """
    return reduce(operator.getitem, items, root_dict)

def set_by_path(root_dict, items, value, add=False):
    """Set a value in a nested object in root_dict by item sequence.
    by Martijn Pieters (https://stackoverflow.com/questions/14692690/access-nested-dictionary-items-via-a-list-of-keys)
    """    
    get_by_path(root_dict, items[:-1])[items[-1]] = value

def del_by_path(root_dict, items):
    """Delete a key-value in a nested object in root_dict by item sequence.
    by Martijn Pieters (https://stackoverflow.com/questions/14692690/access-nested-dictionary-items-via-a-list-of-keys)
    """
    del get_by_path(root_dict, items[:-1])[items[-1]]

def make_Paths(load_dict):
    
    flat_dict = flatdict.FlatDict(copy.deepcopy(load_dict))
    # Make all paths into Path
    dir_keys = [key.split(':') for key in flat_dict.keys() if 'dir' in key and 'direction' not in key and 'extended_dir' not in key]
    # print(dir_keys)
    for key in dir_keys:
        value = get_by_path(load_dict, key)
        # print('key:', key)
        # # print('value:', value)
        # if isinstance(value, dict) and len(value) == 0: 
        #     pass
        # elif value != None and value != 'NotAssigned' and not isinstance(value, bool):
        #     set_by_path(load_dict, key, Path(value))
        if isinstance(value, str):
            if 'R_' in value or '\\' in value or '.tif' in value:
                print('key:', key)
                print('value:', value, '-', type(value))
                set_by_path(load_dict, key, Path(value))
                print('set_by_path: DONE')
            else: 
                pass
        else: 
            pass
    
    return load_dict

def rename_directory(old_name, new_name):
    path = Path(old_name)
    path.rename(new_name)

def make_tuples(load_dict, tuple_keys): 
    flat_dict = flatdict.FlatDict(copy.deepcopy(load_dict))
    #Make all keys from input list into tuples
    separator = ':'
    for tup in tuple_keys:
        str_tup = separator.join(tup)
        if str_tup in flat_dict.keys(): 
            value = get_by_path(load_dict, tup)
            if value != None:
                set_by_path(load_dict, tup, tuple(value))
        
    return load_dict

# Color palette as RGB
def palette_rbg(name:str, num:int, rgb=True):
    #https://projects.susielu.com/viz-palette?colors=[%22#ffa500%22,%22#6495ed%22,%22#dc133b%22,%22#ade64f%22,%22#c36bea%22,%22#36cbd3%22,%22#a23e27%22,%22#FF1493%22,%22#8b008b%22,%22#3ff44c%22,%22#FF6347%22,%22#C0C0C0%22,%22#FFD700%22,%22#006400%22,%22#00FFFF%22,%22#DA70D6%22,%22#D2691E%22,%22#7FFFD4%22,%22#F0E68C%22,%22#DC143C%22]&backgroundColor=%22white%22&fontColor=%22black%22&mode=%22normal%22
    if name != 'mH_default': 
        rgb_colors = []
        palette =  sns.color_palette(name, num)
        if rgb: 
            for color in palette:
                tup = []
                for value in color:
                    tup.append(round(value*255))
                rgb_colors.append(tuple(tup))
        else: 
            rgb_colors = palette
    else: 
        rgb_colors = ["#ffa500","#6495ed","#dc133b","#ade64f","#c36bea","#36cbd3","#a23e27","#ff1493","#8b008b","#3ff44c",
                        "#ff6347","#c0c0c0","#ffd700","#006400","#00ffff","#da70d6","#d2691e","#7fffd4","#f0e68c","#dc143c"]

    return rgb_colors

#Unique value in tuple
def find_unique_index(tup):
    if len(set(tup))==2:
        if tup[0] == tup[1]:
            return 2
        elif tup[0] == tup[2]:
            return 1
        else:
            return 0
    else: 
        print('Error: len(set(tup)) != 2')
        alert('error_beep')
        return list(tup).index(min(tup))
    
#%% Module loaded
print('morphoHeart! - Loaded funcBasics')