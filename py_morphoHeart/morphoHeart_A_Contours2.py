# -*- coding: utf-8 -*-
"""
morphoHeart - A. CONTOURS

@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
import platform
import numpy as np

# Verify working dir
def setWorkingDir (root_path):
    if platform.system() == 'Windows':
        wd = r'D:\Documents JSP\Dropbox\Dropbox_Juliana\PhD_Thesis\Data_ongoing\LS_ongoing\A_LS_Analysis\py_LSAnalysis\py_morphoHeart'
    else: 
        wd = r'/Users/juliana/Dropbox/Dropbox_Juliana/PhD_Thesis/Data_ongoing/LS_ongoing/A_LS_Analysis/py_LSAnalysis/py_morphoHeart'
    if root_path != wd:
        os.chdir(wd)
        root_path = os.getcwd()
    return root_path

root_path = setWorkingDir(os.getcwd())

from skimage import measure
from vtkplotter import *
from vtkplotter import embedWindow
embedWindow(False)

#%% Importing morphoHeart packages
import morphoHeart_funcBasics as fcBasics
import morphoHeart_funcContours as fcCont
import morphoHeart_funcMeshes as fcMeshes

#%% Get directories and file 
_, _, dir_lsOngoing, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
df_dataset = fcBasics.exportDatasetCSV(dir_lsOngoing, dir_data2Analyse)
# Get file to process and directories 
folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]
# directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse
dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')
# Import the metadata to know pixel size and distance between slices
xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]

initial_runA = True
initial_runB = True
exit_txt = False

#%% i. CLOSE CONTOURS
q_closeCont = fcBasics.ask4input('Do you want to close the contours of '+filename+' [0]:no/[1]:yes?: ','int')
if q_closeCont == 1:
    # Get channel(s)
    channels = fcCont.selectChannel()
    for ch in channels:
        if initial_runA == True: 
            # Import stack, mask and mask stack
            stack, stack_o, _, filetype = fcCont.openStack(filename = filename, chStr = ch, 
                                                    dir_ims2Analyse = directories[5], 
                                                    dir_txtNnpy = directories[1])
            maskSt = fcCont.openMask(filename = filename, chStr = ch, 
                              dir_ims2Analyse = directories[5], filetype = filetype)
            stack, stack_o = fcCont.maskStack(stack, stack_o, maskSt, filetype = filetype)
        # Show contours 
        # fcCont.showContours(myStack = stack, slices = (0,len(stack)), minLenContour = 250, 
        #                               plotEvery = 25, figureSize = (5,5), plotshow= True)
        fcCont.showGridContours(myStack = stack, slices = (0,len(stack)), n_rows = 6)
        
        # >> Save contours v0
        q_saveCont = fcBasics.ask4input('Do you want to save the images with overlay of contours [0]:no/[1]:yes?: ','int')
        if q_saveCont == 1: 
            fcCont.savePltContours(dir_ims2Analyse = directories[5], filename = filename, 
                                         myStack = stack, chStr = ch, 
                                         slices = (0,len(stack)), contVersion = '0')
    
        # Save stack shape
        stack_shape = np.array([stack.shape[1], stack.shape[2], stack.shape[0]])
        fcCont.saveStShape(filename = filename, dir_txtNnpy = directories[1], stack_shape = stack_shape)
        
        # Define the first and last slice in the channel where there are is heart tissue
        print('\nChannel: ',ch,'- Number of slices: ', stack.shape[0])
        slc_first =fcBasics.ask4input('- Heart tissue layer starts appearing in slice (inclusive): ', 'int'); slc_first_o = slc_first
        slc_last = fcBasics.ask4input('- .... and is no longer present from slice (inclusive) + 1 onwards: ','int')
        
        # >> Automatic contour closure 
        q_autom = fcBasics.ask4input('NEXT: Automatic contour closure \n\t-> Continue [0]:no/[1]:yes?: ','int')
        if q_autom == 1:     
            # Automatic close contours of stack
            stack_closed = fcCont.createInitialClosedStack(myStack = stack, slices = (slc_first, slc_last))
            slc_first_close = slc_first 
            slc_last_close = slc_last
            stack_closed = fcCont.automCloseStackContours(myStack = stack_closed, ch = ch, slices = (slc_first_close,slc_last_close), 
                                                                new_stack = stack_closed, plotEvery = 100)
            fcCont.showGridContours(myStack = stack_closed, slices = (0,len(stack)), n_rows = 6)
            # fcCont.showContours(myStack = stack_closed, slices = (0,len(stack_closed)), minLenContour = 250, 
            #                           plotEvery = 1, figureSize = (5,5), plotshow= True)
            # Save contours v1
            fcCont.savePltContours(dir_ims2Analyse = directories[5], filename = filename, 
                                         myStack = stack_closed, chStr = ch, 
                                         slices = (0,len(stack_closed)), contVersion = '1')
            fcCont.saveStackAsNPY(myStack = stack_closed, filename = filename, 
                                        chStr = ch, stage = 'closedCJ', dir2save = directories[1])
            initial_runA = False
        else: 
            stack_closed = stack
            
        # >> Manually close remaining contours
        q_manual = fcBasics.ask4input('NEXT: Finish closing contours manually \n\t-> Continue [0]:no/[1]:yes?: ','int')
        if q_manual == 1: 
            # Process manually each slice and save number of contours per slice 
            stack_closed, slc_first_o = fcCont.manuallyCloseContours (stack_closed = stack_closed, 
                                                stack_o = stack_o, slices =(187, slc_last), 
                                                n_rows = 6, chStr = ch)
            fcCont.showGridContours(myStack = stack_closed, slices = (0,len(stack)), n_rows = 6)
            # fcCont.showContours(myStack = stack_closed, slices = (slc_first,slc_last), minLenContour = 250, 
            #                           plotEvery = 1, figureSize = (5,5), plotshow= True)
            initial_runA = False
            fcCont.savePltContours(directories[5], filename, stack_closed, ch, 
                                              slices = (0,len(stack_closed)), contVersion = "mC")
            fcCont.saveStackAsNPY(stack_closed, filename, ch, 'closedCJ', directories[1])
        else: 
            stack_closed = stack
            
        # >> Close inflow and outflow tracts
        q_infOutf = fcBasics.ask4input('NEXT: Close inflow and outflow tracts \n\t-> Continue [0]:no/[1]:yes?: ','int')
        if q_infOutf == 1: 
            for region in ['Inflow', 'Outflow']:
                q_region = fcBasics.ask4input('>>> Do you want to close -'+region+'- tract [0]:no/[1]:yes?: ','int')
                if q_region == 1:
                    # Close Inflow and outflow tract
                    slc_reg_first =fcBasics.ask4input('- Close '+region+' from slice (inclusive): ','int')
                    slc_reg_last = fcBasics.ask4input('- .... till slice (inclusive): ','int')
                    stack_closed, slices_reg = fcCont.closeInfOutfStack(stack_closed = stack_closed, 
                                              slices = (slc_reg_first,slc_reg_last), chStr = ch)
                    slc_reg_first, slc_reg_last = slices_reg
            
                    fcCont.showContours(myStack = stack_closed, slices = (slc_reg_first,slc_reg_last), minLenContour = 250, 
                                              plotEvery = 1, figureSize = (5,5), plotshow= True)
            initial_runA = False
        if q_manual == 1 or q_infOutf == 1:
            # Save contours vCJ
            fcCont.savePltContours(directories[5], filename, stack_closed, ch, 
                                              slices = (0,len(stack_closed)), contVersion = "CJ")
            fcCont.saveStackAsNPY(stack_closed, filename, ch, 'closedCJ', directories[1])
            
#%% Check
first = 425
last = 460
fcCont.showContours(myStack = stack_closed, slices = (first, last), minLenContour = 250, 
                    plotEvery = 1, figureSize = (5,5), plotshow= True)
#%% ii. SELECT CONTOURS
if q_closeCont == 0:
    # Get channel(s)
    channels = fcCont.selectChannel()
    
for ch in channels:
    # Import stack, mask and mask stack
    stack, stack_o, file2Analyse, filetype = fcCont.openStack(filename = filename, chStr = ch, 
                                                dir_ims2Analyse = directories[5], 
                                                dir_txtNnpy = directories[1])
    maskSt = fcCont.openMask(filename = filename, chStr = ch, 
                          dir_ims2Analyse = directories[5], filetype = filetype)
    stack, _ = fcCont.maskStack(stack, stack_o, maskSt, filetype = filetype)
    stack_shape = np.array([stack.shape[1], stack.shape[2], stack.shape[0]])
    
    # Show contours
    fcCont.showContours(myStack = stack, slices = (0,len(stack)), minLenContour = 250, 
                                  plotEvery = 25, figureSize = (5,5), plotshow= True)
    # Define number of contours per slice
    slcCont, numCont = fcCont.getSlicesContNum(stack)
    # Show tuples according to slcContours
    tuple_slc, numCont, slcCont = fcCont.tuple_pairs (numCont, slcCont, False)
    if initial_runB: 
        # Creation of heartLayer dictionary with basic info
        heartLayer, version_N, update_dict = fcCont.dictCreation(filename = filename, file = file2Analyse, minLenContour = 250, 
                                                chStr = ch, slcCont = slcCont, numCont = numCont, xy_Scaling_um = xy_Scaling_um, 
                                                z_Scaling_um = z_Scaling_um, heartLayer = dict(), update = False, curr_v= 0)
        initial_runB = False
    # Fill dictionary
    while len(tuple_slc) != 0:
        heartLayer, tuple_slc, exit_txt = fcCont.dictFill (tuple_slc = tuple_slc, numContours = numCont, 
                                            stack = stack, chStr = ch, heartLayer = heartLayer, minLenContour = 250)
        if exit_txt: break
    
    if not exit_txt: 
        #% Check slices
        fcCont.plotSelectedContours(imageEvery = 1, stack = stack, heartLayer = heartLayer)
        #% Correct dictionary if needed
        q_modifyDict = fcBasics.ask4input("Do you want to modify the created dictionary [0]:no/[1]:yes?: ",'int')
        if q_modifyDict == 1:
            q_newFile2Analyse = fcBasics.ask4input(">>Do you want to change the image file (file2Analyse)? [0]:no/[1]:yes?: ",'int')
            if q_newFile2Analyse == 1:
                file2Analyse = filename + input("New file2Analyse filename? (e.g. _comp_EDC_close.tif): ")   
            q_newMinNumPts = fcBasics.ask4input(">>Do you want to change the minimum number of points that make up a contour? [0]:no/[1]:yes?: ",'int')
            if q_newMinNumPts == 1:
                minNumPts = fcBasics.ask4input("Enter new minimum number of points: ",'int')
            if q_newFile2Analyse == 1 or q_newMinNumPts == 1:
                heartLayer, version_N, update_dict = fcCont.dictCreation(filename = filename, file = file2Analyse, minLenContour = minNumPts, 
                                                        channel = ch, xy_Scaling_um = xy_Scaling_um, z_Scaling_um = z_Scaling_um, 
                                                        heartLayer = heartLayer, update = update_dict, curr_v= version_N)
            while q_modifyDict == 1:    
                heartLayer = fcCont.modifyDict(stack = stack, chStr = ch, heartLayer = heartLayer, 
                                                     minLenContour = 300)  
                q_modifyDict = fcBasics.ask4input("Do you want to modify any other part of the dictionary? [0]:no/[1]:yes?: ",'int')
                if q_modifyDict == 0:
                    break
        #% Check slices
        fcCont.plotSelectedContours(imageEvery = 10, stack = stack, heartLayer = heartLayer)
        
        q_saveDict = fcBasics.ask4input("Do you want to save heartLayer2S dictionary? [0]:no/[1]:yes?: ",'int')
        if q_saveDict == 1:
            #% Create dictionary to save
            heartLayer2S = fcCont.smallDict2Save(heartLayer)
            
            #% Saving a small version of heartLayer dictionary
            import json
            jsonDict_name = filename+"_heartlayer2S_"+ch+".json"
            json2save_dir = os.path.join(directories[0],jsonDict_name)
            
            with open(json2save_dir, "w") as write_file:
                json.dump(heartLayer2S, write_file, cls=fcCont.NumpyArrayEncoder)
            print("- heartLayer2S saved correctly!\n- File: "+jsonDict_name); fcBasics.alert("countdown",1)
        
        # Create 3D numpy arrays and save them
        q_save_s3s = fcBasics.ask4input("Do you want to save the stacks created? [0]:no/[1]:yes?: ",'int')
        _, _, s3_all = fcCont.save_s3s_fromDict(filename = filename, chStr = ch, stack_shape = stack_shape, 
                  heartLayer = heartLayer, dir_txtNnpy = directories[1], save = q_save_s3s)

        #Find surfaces in 3D for channel analysed
        mesh = fcMeshes.createLayerMesh(filename = filename, s3 = s3_all, resolution = res, layer = 'Channel '+str(ch), 
                                        name = 'Channel '+str(ch), colour = 'cornflowerblue', alpha = 1, plotshow=True)
        
        
        
        