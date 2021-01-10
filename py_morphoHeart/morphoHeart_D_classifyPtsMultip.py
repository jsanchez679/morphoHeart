# -*- coding: utf-8 -*-
"""
morphoHeart - D. CLASSIFY POINTS (ATRIUM/VENTRICLE, DORSAL/VENTRAL, LEFT/RIGHT) 
@author: Juliana Sanchez-Posada
"""

#%% Importing packages
# Python
import os
import platform
import numpy as np
import json
import math
import pandas as pd

from datetime import datetime
from time import perf_counter
from progress.bar import Bar
import multiprocessing

#% morphoHeart 
import morphoHeart_funcBasics as fcBasics 
from morphoHeart_funcBasics import alert

def main():
    global AnV; global DnV_Atr; global DnV_Vent; global cj_thickness
    global x_pts_left; global y_pts_left; global z_pts_left
    global pts_classLnR; global pts_classAtnVent; global pts_classDnV
    
    suffix = '%(index)d/%(max)d - %(elapsed)ds'
    root_path = root_path = os.getcwd()
    # Get main directories (check which ones are actually used)
    _, _, dir_lsOngoing, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_lsOngoing, dir_data2Analyse)
    # Get file to process and directories 
    folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse)
    
    # FUNCTIONS
    # func - importDicts
    def importDicts(filename, dict_name, dict_directory):
        
        print('- Importing dictionary '+ dict_name)
    
        jsonDict_name = filename+"_"+dict_name+".json"
        json2open_dir = os.path.join(dict_directory,jsonDict_name)
        #print("Started Reading JSON file")
        with open(json2open_dir, "r") as read_file:
            print("\t>> "+dict_name+": Converting JSON encoded data into Numpy Array")
            decodedArray = json.load(read_file)
        
        alert('jump',1)
        print('\n')
        
        return decodedArray
    
    # func - intersecOfSets
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
    
    # func - classifyPtAtSidesofPlane
    def classifyPtAtSidesofPlane(pt, d, normal, classes_sorted):
        # x,y,z = pt
        # a,b,c = normal
        
        dotProd_pt = np.dot(pt, normal)#dot((a,b,c), (x,y,z))
        #print("dotProd_pt:", dotProd_pt)
        if dotProd_pt < d:
            pts_class = classes_sorted[0]
        else:
            pts_class = classes_sorted[-1]
        #print("Point classified as:", pts_class)
        
        return pts_class
    
    # func - decodeDict
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
        meas_param = np.asarray(dict2classify['param_'+info])
        
        return [AnV, DnV_Atr, DnV_Vent, pts_left, pts_whole, meas_param]
    
    # func - classifyPts
    def classifyPts(pt):
                
        # Classify mesh points at either side of extended centreline (left/right)
        # for index, pt in enumerate(pts_whole): 
        #print('index:', index)
        # look first for each individual coordinate
        x = pt[0]
        where_x = np.where(x_pts_left == x)[0]
        y = pt[1]
        where_y = np.where(y_pts_left == y)[0]
        z = pt[2]
        where_z = np.where(z_pts_left == z)[0]
        # Unify information from the three coordinates 
        where_final = intersecOfSets(where_x, where_y, where_z)
       
        # Classify point as left/right -> If there is an intersection, point belongs to the left side, if not it belongs to the right
        if len(where_final) == 0:
            pts_classLnR.append('right')
        else: 
            pts_classLnR.append('left')
            
        # Classify point as atrium/ventricle
        pt_classAnV = classifyPtAtSidesofPlane(pt = [x,y,z], d = d_AnV, normal = normal_AnV, 
                                               classes_sorted = classSorted_AnV)
        pts_classAtnVent.append(pt_classAnV)
        
        # Classify point as dorsal/ventral using atrial or ventricular coronal plane
        if pt_classAnV == 'atrium':
            pt_classDnV = classifyPtAtSidesofPlane(pt = [x,y,z], d = d_DnV_Atr, normal = normal_DnV_Atr, 
                                                   classes_sorted = classSorted_DnV_Atr)
        elif pt_classAnV == 'ventricle':
            pt_classDnV = classifyPtAtSidesofPlane(pt = [x,y,z], d = d_DnV_Vent, normal = normal_DnV_Vent, 
                                                   classes_sorted = classSorted_DnV_Vent)
        pts_classDnV.append(pt_classDnV)
        
        # if num % num_pts == 0:
        #     bar.next()
        # if index == 100000:
        #     break
        
        pts_classFinal = [pts_classLnR, pts_classAtnVent, pts_classDnV]
        
        return pts_classFinal
    
    def ptsClass2df(param_val, name, pts_classFinal):
    
        pts_classLnR, pts_classAtnVent, pts_classDnV = pts_classFinal
        df_ptsClass = {name:param_val,'AtrVent':pts_classAtnVent,
                    'DorsVent': pts_classDnV, 'LeftRight': pts_classLnR}
        df_ptsClass = pd.DataFrame(df_ptsClass, columns = [name,'AtrVent',
                                                       'DorsVent','LeftRight'])
        
        print(df_ptsClass.sample(10))
    
        return df_ptsClass

    # func - saveDF
    def saveDF(filename, df2save, df_name, dir2save):
        
        name_csv = filename+'_'+df_name+'.csv'
        dir_DF = os.path.join(dir2save, name_csv)
        df2save.to_csv(dir_DF)
        
        print('-', df_name, ' has been saved!')
        alert('countdown',1)
    
    # ------------------------------------------------------------------------------------
    print('Entrooooo!!')

    # Import dictionaries
    dict_cjTh2classify = importDicts(filename = filename, dict_name = 'classCJTh', dict_directory = directories[0])
    # cj_thickness = np.asarray(dict_cjTh2classify['param_'+'cjThickness'])
    [AnV, DnV_Atr, DnV_Vent, pts_left, pts_whole, cj_thickness] = decodeDict (dict_cjTh2classify, info = 'cjThickness')
    # global AnV; global DnV_Atr; global DnV_Vent; global cj_thickness
    
    # - Get x, y, and z arrays of pts that make up the left side of the mesh
    x_pts_left = pts_left[:,0]; y_pts_left = pts_left[:,1]; z_pts_left = pts_left[:,2]
    
    # Create empty lists
    pts_classLnR = []
    pts_classAtnVent = []
    pts_classDnV = []
    
    chunks = [pts_whole[i::8] for i in range(9)]
    
    # Start
    tic = perf_counter()
    pool = multiprocessing.Pool(processes = 8)
    pts_classCJTh = pool.map_async(classifyPts, chunks)
    
    #End
    toc = perf_counter()
    time = toc-tic
    alert('wohoo',1)
    print("- All Done - points have been classified!")
    print("- Time taken to classify = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")     

    #pts_classFinal = [pts_classLnR, pts_classAtnVent, pts_classDnV]
    print('len:',len(pts_classLnR),len(pts_classAtnVent),len(pts_classDnV))
    
    # pts_classCJTh = classifyPts(dict2classify = dict_cjTh2classify, info = 'cjThickness')
    
    df_cjTh = ptsClass2df(param_val = cj_thickness, name = 'cjThickness', pts_classFinal = [pts_classLnR, pts_classAtnVent, pts_classDnV])
    saveDF(filename = filename, df2save = df_cjTh, df_name = 'df_cjThAuto', dir2save = dir_results)
    
if __name__ == '__main__':
    main()