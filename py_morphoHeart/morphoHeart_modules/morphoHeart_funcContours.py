# -*- coding: utf-8 -*-
"""
morphoHeart_funcContours

@author: Juliana Sanchez-Posada
Version: April 13, 2021
"""
#%% Importing python packages
import os
import numpy as np
from scipy.spatial import distance
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import cv2

from pylab import *
# from pylab import rc #*
rc('axes', linewidth=1, edgecolor = 'w')

from progress.bar import Bar
suffix = '%(index)d/%(max)d - %(elapsed)ds - avg: %(avg)ds'

from datetime import datetime
from time import perf_counter
import json
from skimage import measure, io
from skimage.measure import label, regionprops
from skimage.draw import line_aa
import scipy.ndimage as ndimage
from scipy.spatial import ConvexHull

#%% Importing morphoHeart packages
from .morphoHeart_funcBasics import alert, ask4input, saveDict, loadDicts

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

#%% Color definitions
colors =  ['C0', 'C1', 'C2','C3','C4','C5','C6','C8','C9','#B0C4DE',
           '#DAA520','#FF7F50','#C71585','#7CFC00','#20B2AA',
           '#00FFFF','#8A2BE2','#FF1493','#D2691E','#DC143C']*20

#%% -MAIN: CLOSE CONTOURS
#%% func - main_importImages
def main_importImages (filename, channel, directories, n_rows):
    """
    Function that imports stack and stack_shape

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    channel : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    directories : list of paths
        List of paths to the folders where the files are saved.
        [0. dir_dict, 1. dir_txtNnpy, 2. dir_stl, 3. dir_cl, 4. dir_imsNvideos, 5. dir_ims2Analyse, 6. dir_folder2A.]
    n_rows : int
        Number of rows in subplot

    Returns
    -------
    stack : numpy array
        Stack to process.
    stack_o : numpy array
        Copy of stack to process.
    stack_m : numpy array
        Copy of original masked stack to process.
    file : str
        Loaded file name (.tif/.npy)
    stack_shape : numpy array
        Stack shape [x,y,z]
    processDict : dictionary
        Updated or newly created dictionary with information about the processed carried out so far in A_Contours script
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    """

    # Import stack, mask and mask stack
    stack, stack_o, stack_m, file, filetype = openStack(filename = filename, chStr = channel,
                                                        dir_ims2Analyse = directories[5],
                                                        dir_txtNnpy = directories[1])
    maskSt = openMask(filename = filename, chStr = channel,
                      dir_ims2Analyse = directories[5], filetype = filetype)
    stack, stack_o, stack_m = maskStack(stack, stack_o, stack_m, maskSt, filetype = filetype)
    
    # Setting level 
    try:
        [processDict] = loadDicts(filename, ['processDict'], [directories[1]], print_txt = False)
        level = processDict[channel]['level']
    except: 
        if 'LSCJ' in filename: 
            chLevel = ask4input('Do you want to change the level to select contours from 0.5 to 200? [0]: no / [1]: yes!', bool)
            if chLevel:
                level = 200
            else: 
                level = 0.5
        else: 
            level = 0.5
    print('-> Level : '+str(level))
    
    # Show contours
    print('\n- Plotting contours of imported stack... (Plots tab)')
    showGridContours(myStack = stack, slices = (0,len(stack)), n_rows = n_rows, level = level)

    # Save stack shape
    stack_shape = np.array([stack.shape[1], stack.shape[2], stack.shape[0]])
    saveStShape(filename = filename, dir_txtNnpy = directories[1], stack_shape = stack_shape)
    
    try:
        [processDict] = loadDicts(filename, ['processDict'], [directories[1]], print_txt = False)
        try:
            ch_processDict = processDict[channel]
        except: 
            print('\n- Channel: ',channel,'- Number of slices: ', stack.shape[0])
            slc_first =ask4input('- Heart tissue layer ('+channel+') starts appearing in slice (inclusive): ', int)
            slc_last = ask4input('- .... and is no longer present from slice onwards: ',int)#+1???
            ch_processDict = processDict[channel] = dict()
            ch_processDict['G-Slc_tissueLayerFirst'] = slc_first
            ch_processDict['G-Slc_tissueLayerLast'] = slc_last
            ch_processDict['A-AutomCloseContours'] = 'NO'
            ch_processDict['B-ManualCloseContours'] = 'NO'
            ch_processDict['C-ClosedInflowTract'] = 'NO'
            ch_processDict['C-ClosedOutflowTract'] = 'NO'
            ch_processDict['G-LastProcessingStep'] = '0'
            ch_processDict['level'] = level
            
    except: 
        # Define the first and last slice in the channel where there is heart tissue
        print('\n- Channel: ',channel,'- Number of slices: ', stack.shape[0])
        slc_first =ask4input('- Heart tissue layer ('+channel+') starts appearing in slice (inclusive): ', int)
        slc_last = ask4input('- .... and is no longer present from slice onwards: ',int)#+1???
        
        processDict = dict()
        ch_processDict = processDict[channel] = dict()
        ch_processDict['G-Slc_tissueLayerFirst'] = slc_first
        ch_processDict['G-Slc_tissueLayerLast'] = slc_last
        ch_processDict['A-AutomCloseContours'] = 'NO'
        ch_processDict['B-ManualCloseContours'] = 'NO'
        ch_processDict['C-ClosedInflowTract'] = 'NO'
        ch_processDict['C-ClosedOutflowTract'] = 'NO'
        ch_processDict['G-LastProcessingStep'] = '0'
        ch_processDict['level'] = level
    
    
    # # >> Save contours v0
    # q_saveCont = ask4input('Do you want to save the images as PNGs with overlay of contours? [0]:no/[1]:yes: ',bool)
    # if q_saveCont:
    #     savePltContours(dir_ims2Analyse = directories[5], filename = filename,
    #                                   myStack = stack, chStr = channel,
    #                                   slices = (0,len(stack)), contVersion = '0')

    return stack, stack_o, stack_m, file, stack_shape, processDict, level

#%% A. func - main_automCloseCont
def main_automCloseCont(filename, channel, directories, stack, plotEvery, n_rows, processDict, level):
    """
    Function that automatically closes the contours of the stack being processed

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    channel : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    directories : list of paths
        List of paths to the folders where the files are saved.
        [0. dir_dict, 1. dir_txtNnpy, 2. dir_stl, 3. dir_cl, 4. dir_imsNvideos, 5. dir_ims2Analyse, 6. dir_folder2A.]
    stack : numpy array
        Stack to process.
    plotEvery : int
        Number defining the interval in which slices will be plotted.
    n_rows : int
        Number of rows in subplot
    processDict : dictionary
        Dictionary with information about the processed carried out so far in A_Contours script
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    stack_closed : numpy array
        Resulting closed stack.
    processDict : dictionary
        Updated dictionary with information about the processed carried out so far in A_Contours script
    done_autom : boolean
        True if automatic closure of contours has already been performed, else False.

    """
    stack_closed = stack
    done_autom = False
    already_done = False
    print('\n- FIRST: Automatically close contours')
    
    # >> Automatic contour closure
    ch_processDict = processDict[channel]
    if ch_processDict['A-AutomCloseContours'] == 'DONE':
        q_manual = ask4input('You already run the automatic closure of the contours in this file/channel. Do you want to run it again? [0]:no/[1]:yes: ', bool)
        if q_manual:
            slc_firstAutom = ask4input('- Start closing contours automatically from slice (inclusive): ', int)
            slc_lastAutom = ask4input('-  .... until slice (inclusive): ',int)#+1???
            # Automatic close contours of stack
            stack_closed = createInitialClosedStack(myStack = stack, slices = (slc_firstAutom, slc_lastAutom))
            stack_closed, done_autom = automCloseStackContours(myStack = stack, ch = channel, slices = (slc_firstAutom, slc_lastAutom+1),
                                                    new_stack = stack_closed, plotEvery = plotEvery, level = level)
            showGridContours(myStack = stack_closed, slices = (0,len(stack)), n_rows = n_rows, level = level)
            ch_processDict['A-AutomCloseContours'] = 'DONE'
            ch_processDict['G-LastProcessingStep'] = 'A'
        else:
            done_autom = True
            already_done = True
            stack_closed = stack
    
    else: #First time processing:
        q_autom = ask4input('Do you want to AUTOMATICALLY CLOSE THE CONTOURS of '+filename+'/'+channel+'? [0]:no/[1]:yes: ',bool)
        if q_autom:
            # Define the first and last slice in the channel where there are is heart tissue
            print('\n- Channel: ',channel,'- Number of slices: ', stack.shape[0])
            slc_firstAutom = ask4input('- Start closing contours automatically from slice (inclusive): ', int)
            slc_lastAutom = ask4input('-  .... until slice (inclusive): ',int)#+1???

            # Automatic close contours of stack
            stack_closed = createInitialClosedStack(myStack = stack, slices = (slc_firstAutom, slc_lastAutom))
            stack_closed, done_autom = automCloseStackContours(myStack = stack, ch = channel, slices = (slc_firstAutom, slc_lastAutom+1),
                                                    new_stack = stack_closed, plotEvery = plotEvery, level = level)
            showGridContours(myStack = stack_closed, slices = (0,len(stack)), n_rows = n_rows, level = level)
            ch_processDict['A-AutomCloseContours'] = 'DONE'
            ch_processDict['G-LastProcessingStep'] = 'A'
    
    saveDict(filename, processDict , 'processDict', directories[1], False)
    
    # # >> Automatic contour closure
    # try:
    #     [processDict] = loadDicts(filename, ['processDict'], [directories[1]], print_txt = False)
    #     try:
    #         ch_processDict = processDict[channel]
    #         if ch_processDict['A-AutomCloseContours'] == 'DONE':
    #             q_manual = ask4input('You already run the automatic closure of the contours in this file/channel. Do you want to run it again? [0]:no/[1]:yes: ', bool)
    #             if q_manual:
    #                 slc_first = ch_processDict['G-Slc_tissueLayerFirst']
    #                 slc_last = ch_processDict['G-Slc_tissueLayerLast']
    #                 # Automatic close contours of stack
    #                 stack_closed = createInitialClosedStack(myStack = stack, slices = (slc_first, slc_last))
    #                 stack_closed, done_autom = automCloseStackContours(myStack = stack, ch = channel, slices = (slc_first,slc_last),
    #                                                         new_stack = stack_closed, plotEvery = plotEvery)#100
    #                 showGridContours(myStack = stack_closed, slices = (0,len(stack)), n_rows = n_rows)
    #                 ch_processDict['A-AutomCloseContours'] = 'DONE'
    #                 ch_processDict['G-LastProcessingStep'] = 'A'
    #             else:
    #                 done_autom = True
    #                 already_done = True
    #                 stack_closed = stack
    #     except:
    #         q_autom = ask4input('Do you want to AUTOMATICALLY CLOSE THE CONTOURS of '+filename+'/'+channel+'? [0]:no/[1]:yes: ',bool)
    #         if q_autom:
    #             # Define the first and last slice in the channel where there are is heart tissue
    #             print('\n- Channel: ',channel,'- Number of slices: ', stack.shape[0])
    #             slc_first =ask4input('- Heart tissue layer ('+channel+') starts appearing in slice (inclusive): ', int)
    #             slc_last = ask4input('- .... and is no longer present from slice onwards: ',int)#+1???

    #             # Automatic close contours of stack
    #             stack_closed = createInitialClosedStack(myStack = stack, slices = (slc_first, slc_last))
    #             stack_closed, done_autom = automCloseStackContours(myStack = stack, ch = channel, slices = (slc_first,slc_last),
    #                                                                 new_stack = stack_closed, plotEvery = plotEvery)#100
    #             showGridContours(myStack = stack_closed, slices = (0,len(stack)), n_rows = n_rows)

    #             ch_processDict = processDict[channel] = dict()
    #             ch_processDict['G-Slc_tissueLayerFirst'] = slc_first
    #             ch_processDict['G-Slc_tissueLayerLast'] = slc_last
    #             ch_processDict['A-AutomCloseContours'] = 'DONE'
    #             ch_processDict['B-ManualCloseContours'] = 'NO'
    #             ch_processDict['C-ClosedInflowTract'] = 'NO'
    #             ch_processDict['C-ClosedOutflowTract'] = 'NO'
    #             ch_processDict['G-LastProcessingStep'] = 'A'
                
    #     saveDict(filename, processDict , 'processDict', directories[1], False)
        
    # except: #First time processing:
    #     q_autom = ask4input('Do you want to AUTOMATICALLY CLOSE THE CONTOURS of '+filename+'/'+channel+'? [0]:no/[1]:yes: ',bool)
    #     if q_autom:
    #         # Define the first and last slice in the channel where there are is heart tissue
    #         print('\n- Channel: ',channel,'- Number of slices: ', stack.shape[0])
    #         slc_first =ask4input('- Heart tissue layer ('+channel+') starts appearing in slice (inclusive): ', int)
    #         slc_last = ask4input('- .... and is no longer present from slice onwards: ',int)#+1???

    #         # Automatic close contours of stack
    #         stack_closed = createInitialClosedStack(myStack = stack, slices = (slc_first, slc_last))
    #         stack_closed, done_autom = automCloseStackContours(myStack = stack, ch = channel, slices = (slc_first,slc_last),
    #                                                             new_stack = stack_closed, plotEvery = plotEvery)#100
    #         showGridContours(myStack = stack_closed, slices = (0,len(stack)), n_rows = n_rows)

    #         processDict = dict()
    #         ch_processDict = processDict[channel] = dict()
    #         ch_processDict['G-Slc_tissueLayerFirst'] = slc_first
    #         ch_processDict['G-Slc_tissueLayerLast'] = slc_last
    #         ch_processDict['A-AutomCloseContours'] = 'DONE'
    #         ch_processDict['B-ManualCloseContours'] = 'NO'
    #         ch_processDict['C-ClosedInflowTract'] = 'NO'
    #         ch_processDict['C-ClosedOutflowTract'] = 'NO'
    #         ch_processDict['G-LastProcessingStep'] = 'A'
        
    #     else: 
    #         processDict = dict()
    #         ch_processDict = processDict[channel] = dict()
    #         ch_processDict['G-Slc_tissueLayerFirst'] = 0
    #         ch_processDict['G-Slc_tissueLayerLast'] = len(stack)
    #         ch_processDict['A-AutomCloseContours'] = 'DONE'
    #         ch_processDict['B-ManualCloseContours'] = 'NO'
    #         ch_processDict['C-ClosedInflowTract'] = 'NO'
    #         ch_processDict['C-ClosedOutflowTract'] = 'NO'
    #         ch_processDict['G-LastProcessingStep'] = 'A'
                    
    # saveDict(filename, processDict , 'processDict', directories[1], False)

    if done_autom:
        if not already_done:
            # # Save contours v1
            # q_saveCont = ask4input('Do you want to save the images as PNGs with overlay of contours? [0]:no/[1]:yes: ',bool)
            # if q_saveCont:
            #     savePltContours(dir_ims2Analyse = directories[5], filename = filename,
            #                                  myStack = stack_closed, chStr = channel,
            #                                  slices = (0,len(stack_closed)), contVersion = '1')
            # Save contours v1
            # q_saveNPY = ask4input("Do you want to save the stack automatically closed? \n\t[0]:no, I am going to continue processing, so I'll save it on a later stage \n\t[1]:yes, I might stop processing now / yes, I want to save it just in case\n >>> : ",bool)
            # if q_saveNPY:
            saveStackAsNPY(myStack = stack_closed, filename = filename,
                           chStr = channel, stage = 'closedCJ', dir2save = directories[1])

    return stack_closed, processDict, done_autom

#%% B. func - main_manuallyCloseContours
def main_manuallyCloseContours(filename, channel, directories, stack_closed, stack_o, stack_m, processDict, n_rows, level, checking = False):
    """
    Function that allows the user to manually close/clean the contours of the stack being processed

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    channel : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    directories : list of paths
        List of paths to the folders where the files are saved.
        [0. dir_dict, 1. dir_txtNnpy, 2. dir_stl, 3. dir_cl, 4. dir_imsNvideos, 5. dir_ims2Analyse, 6. dir_folder2A.]
    stack_closed : numpy array
        Input closed stack.
    stack_o : numpy array
        Copy of closed stack.
    stack_m : numpy array
        Copy of original masked stack to process.
    processDict : dictionary
        Dictionary with information about the processed carried out so far in A_Contours script
    n_rows : int
        Number of rows in subplot
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours
    checking : boolean
        True if you are running this function to just close some slices after closing inf/out tracts, else False.
        The default is False

    Returns
    -------
    stack_closed : numpy array
        Resulting closed stack.
    processDict : dictionary
        Updated dictionary with information about the processed carried out so far in A_Contours script
    done_manual : boolean
        True if manual closure of contours has already been performed, else False.

    """
    
    exit_code = False
    q_saveNPY = False
    print('\n- NEXT: Manually close contours and clean slices')
    # >> Manually close remaining contours
    done_manual = False
    if processDict[channel]['B-ManualCloseContours'] == 'DONE':
        q_manual = ask4input('You already closed all the contours of this file/channel. Do you want to make any additional changes? [0]:no/[1]:yes: ', bool)
        if q_manual:
            slc_first =ask4input('Enter the slice number from which you want to close or modify the stack: ', int)
            slc_last = ask4input('- .... until slice (inclusive): ',int)+1
        else:
            done_manual = True
    elif processDict[channel]['B-ManualCloseContours'] == 'Ongoing':
        try:
            last_slc_closed = processDict[channel]['B-Slc_manuallyClosedLast']
            print('- You already started manually closing the contours of this file/channel (last slice closed: '+str(last_slc_closed)+')')
            q_startFromLast = ask4input('Do you want to manually close the contours of the previously selected slice range ('+str(last_slc_closed)+'-'+str(processDict[channel]['G-Slc_tissueLayerLast'])+')? [0]:no/[1]:yes: ',bool)
            if q_startFromLast:
                slc_first = last_slc_closed
                slc_last = processDict[channel]['G-Slc_tissueLayerLast']
            else:
                slc_first =ask4input('Enter the slice number from which you want to start manually closing the contours: ', int)
                slc_last = ask4input('- .... until slice (inclusive): ',int)+1
        except:
            slc_first =ask4input('Enter the slice number from which you want to start manually closing the contours: ', int)
            slc_last = ask4input('- .... until slice (inclusive): ',int)+1

    else: #processDict[channel]['B-ManualCloseContours'] == 'NO':
        slc_first =ask4input('Enter the slice number from which you want to start manually closing the contours: ', int)
        slc_last = ask4input('- .... until slice (inclusive): ',int)+1

    if not done_manual:
        # Process manually each slice
        stack_closed, last_slc, exit_code = manuallyCloseContours(stack_closed = stack_closed,
                                            stack_o = stack_o, stack_m = stack_m, slices =(slc_first, slc_last),
                                            n_rows = n_rows, chStr = channel, exit_code = False, level = level)
        if exit_code:
            done_manual = False
            processDict[channel]['B-Slc_manuallyClosedLast'] = last_slc
            processDict[channel]['G-LastProcessingStep'] = processDict[channel]['G-LastProcessingStep']+'-B'
            processDict[channel]['B-ManualCloseContours'] = 'Ongoing'
            q_save = ask4input("EXIT: Do you want to save the stack you have closed so far? \n\t[0]:no, I am going to continue processing, so I'll save it on a later stage \n\t[1]:yes, I might stop processing now / yes, I want to save it just in case\n >>> : ",bool)
            if q_save:
                saveStackAsNPY(stack_closed, filename, channel, 'closedCJ', directories[1])
            saveDict(filename, processDict , 'processDict', directories[1], False)
        else:
            # alert('wohoo', 1)
            # print('- Manually closing contours done!')
            if not checking:
                q_plotAll = ask4input('Do you want to plot all the slices to double check everything has been closed? [0]:no/[1]:yes: ', bool)
                if q_plotAll:
                    showGridContours(myStack = stack_closed, 
                                     slices = (processDict[channel]['G-Slc_tissueLayerFirst'],processDict[channel]['G-Slc_tissueLayerLast']), 
                                     n_rows = n_rows, level = level)
                # # Save contours vMan
                # q_saveCont = ask4input('Do you want to save the images as PNGs with overlay of contours? [0]:no/[1]:yes: ',bool)
                # if q_saveCont:
                #     savePltContours(directories[5], filename, stack_closed, channel,
                #                                   slices = (0,len(stack_closed)), contVersion = "mC")
            # Save closed contours
            q_saveClosedCont = ask4input('Are you done manually closing the contours? [0]:no/[1]:yes: ',bool)
            if q_saveClosedCont:
                done_manual = True
                processDict[channel]['G-LastProcessingStep'] = processDict[channel]['G-LastProcessingStep']+'-B'
                processDict[channel]['B-ManualCloseContours'] = 'DONE'
                processDict[channel]['B-Slc_manuallyClosedLast'] = processDict[channel]['G-Slc_tissueLayerLast']
                # Save contours v1
                if not checking:
                    q_saveNPY = ask4input("Do you want to save the final stack you manually closed? \n\t[0]:no, I am going to continue processing, so I'll save it on a later stage \n\t[1]:yes, I might stop processing now / yes, I want to save it just in case\n >>> : ",bool)
            else:
                continueManClosing = ask4input('Would you like to close some additional slices? [0]:no/[1]:yes: ', bool)
                while continueManClosing: 
                    slc_list, slc_end = getSlices((0, len(stack_closed)), 'you would like to close')
                    exit_code, slc_end, last_slc, stack_closed = manuallyCloseContoursTuple (slc_list, stack_closed, 
                                                                                             stack_o, stack_m, channel, exit_code, level)
                    if exit_code:
                        alert('error',1)
                        print("\n- Exit script - last slice ", last_slc)
                        done_manual = False
                        processDict[channel]['G-LastProcessingStep'] = processDict[channel]['G-LastProcessingStep']+'-B'
                        processDict[channel]['B-ManualCloseContours'] = 'Ongoing'
                        processDict[channel]['B-Slc_manuallyClosedLast'] = last_slc
                        # Save contours v1
                        if not checking:
                            q_saveNPY = ask4input("Do you want to save the stack you manually closed so far? \n\t[0]:no, I am going to continue processing, so I'll save it on a later stage \n\t[1]:yes, I might stop processing now / yes, I want to save it just in case\n >>> : ",bool)
                    continueManClosing = ask4input('Would you like to close some additional slices? [0]:no/[1]:yes: ', bool)
                    
                else:         
                    q_saveClosedCont = ask4input('Are you done manually closing the contours? [0]:no/[1]:yes: ',bool)
                    if q_saveClosedCont: 
                        done_manual = True
                        processDict[channel]['G-LastProcessingStep'] = processDict[channel]['G-LastProcessingStep']+'-B'
                        processDict[channel]['B-ManualCloseContours'] = 'DONE'
                        processDict[channel]['B-Slc_manuallyClosedLast'] = processDict[channel]['G-Slc_tissueLayerLast']
                    else: 
                        done_manual = False
                        processDict[channel]['G-LastProcessingStep'] = processDict[channel]['G-LastProcessingStep']+'-B'
                        processDict[channel]['B-ManualCloseContours'] = 'Ongoing'
                        processDict[channel]['B-Slc_manuallyClosedLast'] = last_slc
                    # Save contours v1
                    if not checking:
                        q_saveNPY = ask4input("Do you want to save the stack you manually closed so far? \n\t[0]:no, I am going to continue processing, so I'll save it on a later stage \n\t[1]:yes, I might stop processing now / yes, I want to save it just in case\n >>> : ",bool)

            if q_saveNPY or checking:
                saveStackAsNPY(stack_closed, filename, channel, 'closedCJ', directories[1])
            saveDict(filename, processDict , 'processDict', directories[1], False)

    return stack_closed, processDict, done_manual

#%% C. func - main_closeInfAndOutfTract
def main_closeInfAndOutfTract(filename, channel, directories, stack_closed, processDict, n_rows, level):
    """
    Function that allows the user to close de inflow and outflow tract of the stack being processed with different methods

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    channel : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    directories : list of paths
        List of paths to the folders where the files are saved.
        [0. dir_dict, 1. dir_txtNnpy, 2. dir_stl, 3. dir_cl, 4. dir_imsNvideos, 5. dir_ims2Analyse, 6. dir_folder2A.]
    stack_closed : numpy array
        Input closed stack.
    processDict : dictionary
        Dictionary with information about the processed carried out so far in A_Contours script
    n_rows : int
        Number of rows in subplot
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    stack_closed : numpy array
        Resulting closed stack.
    processDict : dictionary
        Updated dictionary with information about the processed carried out so far in A_Contours script
    done_infOutf : boolean
        True if inflow and outflow tract closure of contours has already been performed, else False.

    """

    print('\n- NEXT: Close inflow and outflow tracts')
    # >> Close inflow and outflow tract
    done_infOutf = False
    done_InOT = [False, False]
    inf_progress = processDict[channel]['C-ClosedInflowTract']
    outf_progress = processDict[channel]['C-ClosedOutflowTract']
    if inf_progress == 'NO' and outf_progress == 'NO':
        q_infOutf = ask4input('Do you want to CLOSE the INFLOW/OUTFLOW tracts of '+filename+'? [0]:no/[1]:yes: ',bool)
    else:
        q_infOutf = ask4input('Do you want to continue CLOSING the INFLOW/OUTFLOW tracts of '+filename+'? [0]:no/[1]:yes: ',bool)
        if not q_infOutf:
            done_infOutf = ask4input('Are you done CLOSING both the INFLOW/OUTFLOW tracts of '+filename+'? [0]:no/[1]:yes: ',bool)
            
    if q_infOutf:
        for i, region in enumerate(['Inflow', 'Outflow']):
            # print(region, i)
            q_region = ask4input('>>> Do you want to close -'+region+'- tract? [0]:no/[1]:yes: ',bool)
            while q_region:
                try:
                    slc_reg_first = int(processDict[channel]['C-Slc_closing'+region+'TractLast'])
                    slc_reg_last = processDict[channel]['C-Slc_closed'+region+'TractLast']
                    
                    q_startFromLast = ask4input('Do you want to close the last slice range where the '+region+' tract was closed (Slcs:'+str(slc_reg_first)+'-'+str(slc_reg_last)+')? [0]:no/[1]:yes: ',bool)
                    if q_startFromLast:
                        print('- Closing '+ region + ' tract between slices '+str(slc_reg_first)+'-'+str(slc_reg_last))
                    else:
                        slc_reg_first =ask4input('- Close '+region+' tract starting from slice (inclusive): ',int)
                        slc_reg_last =ask4input('- .... until slice (inclusive): ',int)+1

                except:
                    slc_reg_first =ask4input('- Close '+region+' tract starting from slice (inclusive): ',int)
                    slc_reg_last = ask4input('- .... until slice (inclusive): ',int)+1
                    processDict[channel]['C-Slc_closed'+region+'TractFirst'] = slc_reg_first
                    processDict[channel]['C-Slc_closed'+region+'TractLast'] = slc_reg_last

                stack_closed, slices_reg, exit_code = closeInfOutfStack(stack_closed = stack_closed,
                                          slices = (slc_reg_first,slc_reg_last), chStr = channel, exit_code = False, 
                                          region = region, level = level)
                slc_firstOut, _ = slices_reg

                if exit_code:
                    done_InOT[i] = False
                    processDict[channel]['G-LastProcessingStep'] = processDict[channel]['G-LastProcessingStep']+'-C.'+region
                    processDict[channel]['C-Closed'+region+'Tract'] = 'Ongoing'
                    processDict[channel]['C-Slc_closing'+region+'TractLast'] = slc_firstOut
                    print('- EXIT: Closing '+region+' - last slice: '+str(slc_firstOut))
                    if slc_reg_first != slc_firstOut:
                        print('- Printing slices where '+ region + ' tract was closed')
                        showGridContours(myStack = stack_closed, slices = (slc_reg_first,slc_firstOut), n_rows = n_rows, level = level)
                    q_saveNPY = ask4input("Do you want to save the stack you have closed so far before exiting? \n\t[0]:no, I am going to continue processing, so I'll save it on a later stage \n\t[1]:yes, I might stop processing now / yes, I want to save it just in case\n >>> : ",bool)
                    saveDict(filename, processDict , 'processDict', directories[1], False)
                    print('- EXITING!')
                    break
                else:
                    print('- Plotting slices where '+ region + ' tract was closed')
                    showGridContours(myStack = stack_closed, slices = (slc_reg_first,slc_reg_last), n_rows = n_rows, level = level)
                    q_region_done = ask4input('Are you done closing the '+region+' tract? \n\t[0]:no, there is another group of slices that needs to be closed\n\t[1]:yes, continue with next process!: ',bool)
                    if q_region_done:
                        q_region = False
                        processDict[channel]['G-LastProcessingStep'] = processDict[channel]['G-LastProcessingStep']+'-C.'+region
                        processDict[channel]['C-Closed'+region+'Tract'] = 'DONE'
                        processDict[channel]['C-Slc_closing'+region+'TractLast'] = slc_firstOut
                        q_saveNPY = ask4input("Do you want to save the final stack with the "+region+" tract closed? \n\t[0]:no, I am going to continue processing, so I'll save it on a later stage \n\t[1]:yes, I might stop processing now / yes, I want to save it just in case\n >>> : ",bool)
                        saveDict(filename, processDict , 'processDict', directories[1], False)
                    else:
                        q_saveNPY = False
                        processDict[channel]['C-Slc_closing'+region+'TractLast'] = 'New'

                if q_saveNPY:
                    saveStackAsNPY(stack_closed, filename, channel, 'closedCJ', directories[1])

            if not q_region:
                done_InOT[i] = True

        # Save contours inf/Outf
        # q_saveCont = ask4input('Do you want to save the images as PNGs with overlay of contours? [0]:no/[1]:yes: ',bool)
        # if q_saveCont and not exit_code:
        #     savePltContours(directories[5], filename, stack_closed, channel,
        #                                   slices = (0,len(stack_closed)), contVersion = "CJ")

        done_infOutf = done_InOT[0] and done_InOT[1]
        if region == 'Outflow' and done_infOutf:
            alert('wohoo',1)
            print('- All done! - You are done closing inflow and outflow tracts for '+channel+'!')

    return stack_closed, processDict, done_infOutf

#%% - OPEN/LOAD
#%% func - openStack
def openStack (filename, chStr, dir_ims2Analyse, dir_txtNnpy):
    """
    Function that loads the stack from a .tif or a npy file

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    dir_ims2Analyse : path
        Path to folder of images to analyse
    dir_txtNnpy : path
        Path to the folder where the np arrays are saved.

    Returns
    -------
    stack : numpy array
        Stack to process. .
    stack_o : numpy array
        Copy of stack to process. .
    stack_m : numpy array
        Copy of stack to process to use with mask.
    file : str
        Loaded file name (.tif/.npy)
    filetype : str
        Text indicating the extension of the file to process ('tif'/'npy')

    """

    file = filename+"_St_"+chStr+"_closedCJ.npy"
    dir_stack = os.path.join(dir_txtNnpy, file)
    if os.path.exists(dir_stack) == True:
        q_file2Open = ask4input('There is already a processed stack for '+chStr+' - '+ file + ' - \n\t >> Do you want to continue processing this file? [0]:no/[1]:yes: ', bool)
        if q_file2Open:
            filetype = 'npy'
            stack = np.load(dir_stack)
            stack_o = np.load(dir_stack)
            try: 
                file_m = filename + "_"+chStr+"_EDC.tif"
                dir_m = os.path.join(dir_ims2Analyse,str(file_m))
                stack_m = io.imread(dir_m)
                # print('loaded')
            except: 
                # print('not loaded')
                stack_m = np.load(dir_stack)
        else:
            filetype = 'tif'
            file = filename + "_"+chStr+"_EDC.tif"
            q_tif = ask4input('Open '+file+' instead? [0]:no/[1]:yes: ', bool)
            if q_tif:
                dir_stack = os.path.join(dir_ims2Analyse,str(file))
                stack = io.imread(dir_stack)
                stack_o = io.imread(dir_stack)
                stack_m = io.imread(dir_stack)
    else:
        filetype = 'tif'
        file = filename + "_"+chStr+"_EDC.tif"
        dir_stack = os.path.join(dir_ims2Analyse,str(file))
        stack = io.imread(dir_stack)
        stack_o = io.imread(dir_stack)
        stack_m = io.imread(dir_stack)

    print("- Running analysis of \n\t- file: ", file, "\n\t- channel: ", chStr)
    alert('wohoo',1)

    return stack, stack_o, stack_m, file, filetype

#%% func - openMask
def openMask (filename, chStr, dir_ims2Analyse, filetype):
    """
    Function that loads the mask for the channel being processed

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    dir_ims2Analyse : path
        Path to folder of images to analyse
    filetype : str
        Text indicating the extension of the file to process ('tif'/'npy')

    Returns
    -------
    maskSt : numpy array
        Resulting masked stack.

    """
    # Import mask
    if filetype == 'tif':
        maskfile = filename + "_"+chStr+"_mask.tif"
        dir_mask = os.path.join(dir_ims2Analyse,str(maskfile))
        maskSt = io.imread(dir_mask)
    
        # Change mask to bool
        maskSt[maskSt<255] = False
        maskSt[maskSt==255] = True
        maskSt = np.array(maskSt, dtype=bool)

        print("- Using: \n\t- mask: ", maskfile)
        
    else:#'npy'
        try: 
            maskfile = filename + "_"+chStr+"_mask.tif"
            dir_mask = os.path.join(dir_ims2Analyse,str(maskfile))
            maskSt = io.imread(dir_mask)
        
            # Change mask to bool
            maskSt[maskSt<255] = False
            maskSt[maskSt==255] = True
            maskSt = np.array(maskSt, dtype=bool)
    
            print("- Using: \n\t- mask: ", maskfile)
        except: 
            maskSt = []

    return maskSt

#%% func - maskStack
def maskStack (stack, stack_o, stack_m, maskSt, filetype):
    """
    Function that applies the mask to the stack to clean the slices

    Parameters
    ----------
    stack : numpy array
        Stack to process.
    stack_o : numpy array
        Copy of stack to process.
    stack_m : numpy array
        Copy of original stack in which to apply mask.
    maskSt : numpy array
        Resulting masked stack.
    filetype : str
        Text indicating the extension of the file to process ('tif'/'npy')

    Returns
    -------
    stack : numpy array
        Resulting masked stacked.
    stack_o : numpy array
        Copy of resulting masked stacked.

    """

    if filetype == 'tif':
        stack[maskSt == False] = 0
        stack_o[maskSt == False] = 0
        stack_m[maskSt == False] = 0
    else: #'npy'
        if type(maskSt) != list:
            # print('masked...')
            stack_m[maskSt == False] = 0

    return stack, stack_o, stack_m

#%% func - loadStacks
def loadStacks(filename, dir_txtNnpy, end_name, print_txt = True):
    """
    Function that loads a list of stacks to process

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    dir_txtNnpy : path
        Path to the folder where the np arrays are saved.
    end_name : list of str
        List of end_names given to npy stacks to load
    print_txt : bool, optional
        True if confirmation of action is needed, else False. The default is True.

    Returns
    -------
    s3s : list of numpy arrays
        List of loaded stacks
    stackShape : numpy array
        Stack shape [x,y,z]

    """

    #print('- Loading stacks (s3s)')
    s3s = []
    if print_txt:
        bar = Bar('- Loading stacks (s3s)', max=len(end_name), suffix = suffix, check_tty=False, hide_cursor=False)
    for i, name in enumerate(end_name):

        s3_title = filename+"_s3_"+name+".npy"
        s3_dir = os.path.join(dir_txtNnpy, s3_title)
        s3 = np.load(s3_dir)
        s3s.append(s3)
        if print_txt:
            bar.next()

    if print_txt:
        bar.finish()
    shape_txt = filename+"_stackShape.npy"
    shape_dir = os.path.join(dir_txtNnpy, shape_txt)
    stackShape = np.load(shape_dir)

    if print_txt:
        alert('jump',1)

    return s3s, stackShape

#%% func - saveStShape
def saveStShape (filename, dir_txtNnpy, stack_shape):
    """
    Function that saves a npy file with the information about the shape of the stack in dir_txtNnpy given as input

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    dir_txtNnpy : path
        Path to the folder where the np arrays are saved.
    stack_shape : numpy array
        Stack shape [x,y,z]

    Returns
    -------
    None.

    """

    shape_txt = filename+"_stackShape"
    shape_dir = os.path.join(dir_txtNnpy, shape_txt)

    if os.path.isdir(shape_dir) == False:
        np.save(shape_dir, stack_shape)
        print('- '+shape_txt, " file was created!")

#%% func - selectChannel
def selectChannel():
    """
    Function that prompts the user to select the channel that wants to be processed and returns a string

    Returns
    -------
    ch2load : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).

    """

    while True:
        q_ch2load = ask4input('Select channel(s) to load and process [0]:ch0 / [1]:ch1?: ', int)# / [2]:both?: ', int)

        if q_ch2load == 0:
            ch2load = 'ch0'
            break
        elif q_ch2load == 1:
            ch2load = 'ch1'
            break
        else:
            print('\nERROR: Wrong input!')

    return ch2load

#%% - PLOTS
#%% func - showContours
def showContours (myStack, slices, minLenContour, plotEvery, figureSize, plotshow, level):
    """
    Function to plots the contours every X number of slices between the range defined by 'slices'

    Parameters
    ----------
    myStack : numpy array
        Stack being processed.
    slices : tuple
        start_slice, end_slice
    minLenContour : int
        Minimum number of points making up contoura and used to filter contours found in each slice
    plotEvery : int
        Number defining the interval in which slices will be plotted.
    figureSize : tuple
        Figure size (x,y)
    plotshow : boolean
        True to show plot, else False
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    None.

    """

    slc_range = range(slices[0], slices[1])
    print("\n- Plotting contours every "+ str(plotEvery)+" slice(s)\n\t- Minimum number of points per contour: ", minLenContour)
    for index, slc in enumerate(slc_range):
        mult = index % plotEvery
        if mult == 0:
            # Load the slice as myIm
            myIm = myStack[slc,:,:]

            # Find all the contours of the image
            #find_contours(array, level, fully_connected='low', positive_orientation='low', *, mask=None)
            contours = measure.find_contours(myIm, level, 'high', 'high')

            # Display the image and plot an overlay of all contours found that have
            # more than minLenContour points
            fig1, ax1 = plt.subplots(figsize=figureSize)
            ax1.imshow(myIm, cmap=plt.cm.gray)

            # Variable to save the number of contours found
            num_contour = 0

            # Go through all the contours
            for i, contour in enumerate(contours):
                # Get only the contours made up of more than the designated number of points
                if len(contour)>minLenContour:
                    # plot the contour
                    ax1.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color = colors[num_contour])
                    num_contour += 1

            ax1.axis('image')
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_title('Slice num: '+str(slc))

            if plotshow == True:
                plt.show()
            plt.close()

    alert('wohoo',1)

#%% func - showGridContours
def showGridContours(myStack, slices, n_rows, level):
    """
    Function that generates a series of plots with subplots of slices and its contours between the range defined by 'slices'

    Parameters
    ----------
    myStack : numpy array
        Stack being processed.
    slices : tuple
        start_slice, end_slice
    n_rows : int
        Number of rows in subplot
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    None.

    """

    slcs_per_im = n_rows*4

    slices_first = list(range(slices[0],slices[1]+1,slcs_per_im))
    #print(slices_first)
    slices_last =  list(range(slices[0]+slcs_per_im,slices[1]+1,slcs_per_im))
    #print(slices_last)

    if slices_last != slices[-1]:
        slices_last.append(slices[-1])

    for i in range(len(slices_first)):
        slc_tuple = (slices_first[i], slices_last[i])
        plotSlcsRange(myStack, slc_tuple, 'Contours', slcs_per_im, level)

    alert('wohoo',1)

#%% func - savePltContours
def savePltContours (dir_ims2Analyse, filename, myStack, chStr, slices, contVersion, level):
    """
    Function that saves plots of slices with overlay of contours as pngs. Plots are not shown.

    Parameters
    ----------
    dir_ims2Analyse : path
        Path to folder of images to analyse
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    myStack : numpy array
        Stack being processed.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    slices :  tuple
        start_slice, end_slice
    contVersion : int
        version number of contours.
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    None.

    """

    minLenContour = 250#int(input('Enter minimum number of points per contour: '))
    # contVersion = input("Contours folder version v-X-: ")

    # Define folder and tifs name
    dir_TIFs = "_" + chStr + "_Cont" + str(minLenContour) + "v"+str(contVersion)
    tifName = chStr + "_Cont" + str(minLenContour) + "v"+str(contVersion)
    name4dirTIFs = filename + dir_TIFs

    # Create directory to save tifs
    dir4tifs = os.path.join(dir_ims2Analyse,str(name4dirTIFs))
    if os.path.isdir(dir4tifs) == False:
        os.mkdir(dir4tifs)
    tifPath = os.path.join(dir4tifs,tifName)

    print('level:',str(level))
    bar = Bar('Saving slices with contours', max=slices[1]-slices[0], suffix = suffix, check_tty=False, hide_cursor=False)
    for slc in range(slices[0], slices[1]):

        # Load the slice as myIm
        myIm = myStack[slc,:,:]

        # Find all the contours of the image
        contours = measure.find_contours(myIm, level, 'high', 'high')

        # Display the image and plot an overlay of all contours found that have
        # more than 500 points
        fig1, ax1 = plt.subplots(figsize=(10, 10), constrained_layout=False)
        ax1.imshow(myIm, cmap=plt.cm.gray)

        # Variable to save the number of contours found
        num_contour = 0

        # Go through all the contours
        for i, contour in enumerate(contours):
            # Get only the contours made up of more than the designated number of points
            if len(contour)>minLenContour:
                # plot the contour
                ax1.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color = colors[num_contour])
                num_contour += 1

        # if slc % 50 == 0:
            # print("Total number of contours found in slice",slc, "/", slices[1], " is:", num_contour)
        ax1.axis('image')
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_title('Slice num: '+str(slc))

        # Filename for the saved plot using the current slice num
        if slc < 10:
            fileName_slc = tifPath + "_slc00"+str(slc)+".tif"
        elif slc < 100:
            fileName_slc = tifPath + "_slc0"+str(slc)+".tif"
        else:
            fileName_slc = tifPath + "_slc"+str(slc)+".tif"

        # Save the plot into a jpg
        plt.savefig(fileName_slc, bbox_inches='tight', pad_inches=0.2, dpi = 72)
        plt.close()

        bar.next()

    plt.close("all")
    bar.finish()

    alert('wohoo',1)
    print("- All subplots have been saved!")

#%% func - plotSlcsRange
def plotSlcsRange(stack_closed, slices_plot, text, slcs_per_im, level):
    """
    Function that generates a plot with subplots of slices and its contours between the range defined by 'slices_plot'

    Parameters
    ----------
    stack_closed : numpy array
        Input closed stack.
    slices_plot : tuple
        start_slice, end_slice
    text : str
        Text that describes the process being performed at the moment the plot gets displayed
    slcs_per_im : int
        Number of slices as subplots per plot.
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    None.

    """

    slc_plot_list = list(range(slices_plot[0], slices_plot[1]))
    n_im = len(slc_plot_list)

    #Plot
    Im_size = 5
    n_cols = 4
    n_rows = int(slcs_per_im/n_cols)
    h_fig11 = n_rows*Im_size
    w_fig11 = n_cols*Im_size
    #n_im = n_rows*n_cols

    fig11 = plt.figure(figsize=(w_fig11, h_fig11), constrained_layout=False)

    # gridspec inside gridspec
    grid = fig11.add_gridspec(n_rows, n_cols, wspace=0.0, hspace=0.1)

    for im in range(n_im):
        #Get Image and Label
        slc = slc_plot_list[im]
        # print(slc)
        try: 
            myIm = stack_closed[slc][:][:]
            contours, numCont = getContExpCont (myIm, minLenContour = 250, level = level)

            # Plot
            ax = fig11.add_subplot(grid[im])
            ax.imshow(myIm, cmap=plt.cm.gray)
            ax.set_xticks([])
            ax.set_yticks([])
            for n, contour in enumerate(contours):
                ax.plot(contour[:, 1], contour[:, 0], linewidth=1)
            ax.set_title("Slc "+str(slc))
        except:
            continue

    plt.suptitle(text +": Contours for slices ("+ str(slices_plot[0])+'-'+str(slices_plot[1]-1)+')', y=0.925, fontsize=14)
    plt.show()

#%% func - getContExpCont_plt
def getContExpCont_plt (myIm, slcNum, chStr, minLenContour, figSize, level, plot_show = True):
    """
    Function that gets and returns the contours of a particular slice (slcNum) plotting the image with an overlay
    of the contours

    Parameters
    ----------
    myIm : numpy array
        Imaged being processed
    slcNum : int
        Slice number being processed.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    minLenContour : int
        Minimum number of points making up contoura and used to filter contours found in each slice
    figSize : tuple
        Figure size (x,y)
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours
    plot_show : boolean, optional
        True to show plot, else False. The default is True.

    Returns
    -------
    arr_contours : list of arrays
        list of numpy array with coordinates making up the contours of the input slice.

    """

    chNum = int(chStr[-1])
    # Create an empty array to save all the contours of each slice individually
    arr_contours = []

    #Find all the contours of the image
    contours = measure.find_contours(myIm, level, 'high', 'high')

    if plot_show:
        # Display the image and plot an overlay of all contours found that have
        # more than minLenContour points
        fig1, ax1 = plt.subplots(figsize=(figSize, figSize))
        ax1.imshow(myIm, cmap=plt.cm.gray)

    # Go through all the contours
    for num_n, contour in enumerate(contours):
        # Get only the contours made up of more than the designated number of points
        if len(contour)>minLenContour:
            # Append contour to the array
            arr_contours.append(contour)
            if plot_show:
                # plot the contour
                ax1.plot(contour[:, 1], contour[:, 0], linewidth=2)

    if plot_show:
        ax1.set_title("Channel "+str(chNum)+" / Slice "+str(slcNum), fontsize=12)
        ax1.set_xticks([])
        ax1.set_yticks([])
        plt.show()

    return arr_contours

#%% func - plotFilledCont
def plotFilledCont (myIm, allContours, imIntFilledCont, imExtFilledCont, imAllFilledCont, plotshow, slcNum):
    """
    Funtion that plots a subplot of the filled contours for the particular image (myIm) being processed.
    [Internal, External, Layer]

    Parameters
    ----------
    myIm : numpy array
        Imaged being processed
    allContours : list of arrays
        list of numpy array with coordinates making up each contour.
    imIntFilledCont : numpy array of booleans
        Array with filled internal contours of myIm
    imExtFilledCont : numpy array of booleans
        Array with filled external contours of myIm
    imAllFilledCont : numpy array of booleans
        Array with filled heart layer contours of myIm
    plotshow : boolean
        True to show plot, else False
    slcNum : int
        Slice number being processed.

    Returns
    -------
    None

    """

    if plotshow:
        condition = isinstance(allContours, float)
        if condition == False:
            figR, axR = plt.subplots(1,4, figsize=(14, 3.5))
            axR[0].imshow(imIntFilledCont)
            axR[0].set_title("Filled Internal Contours", fontsize=10)
            axR[0].set_xticks([])
            axR[0].set_yticks([])

            axR[1].imshow(imExtFilledCont)
            axR[1].set_title("Filled External Contours", fontsize=10)
            axR[1].set_xticks([])
            axR[1].set_yticks([])

            axR[2].imshow(imAllFilledCont)
            axR[2].set_title("Filled All Contours", fontsize=10)
            axR[2].set_xticks([])
            axR[2].set_yticks([])

            axR[3].imshow(myIm, cmap=plt.cm.gray)
            titleAll = "Image with contours [slc"+str(slcNum)+"]"
            axR[3].set_title(titleAll, fontsize=10)
            for num, contour in enumerate(allContours):
                axR[3].plot(contour[:, 1], contour[:, 0], linewidth=1.5)
            axR[3].set_xticks([])
            axR[3].set_yticks([])

        else:
            figR, axR = plt.subplots(1,4, figsize=(14, 3.5))
            axR[0].imshow(imIntFilledCont)
            axR[0].set_title("Filled Internal Contours", fontsize=10)
            axR[0].set_xticks([])
            axR[0].set_yticks([])

            axR[1].imshow(imExtFilledCont)
            axR[1].set_title("Filled External Contours", fontsize=10)
            axR[1].set_xticks([])
            axR[1].set_yticks([])

            axR[2].imshow(imAllFilledCont)
            axR[2].set_title("Filled All Contours", fontsize=10)
            axR[2].set_xticks([])
            axR[2].set_yticks([])

            axR[3].imshow(myIm, cmap=plt.cm.gray)
            titleAll = "Image (No contours) [slc"+str(slcNum)+"]"
            axR[3].set_title(titleAll, fontsize=10)
            axR[3].set_xticks([])
            axR[3].set_yticks([])


        plt.show()

#%% func - plotSelectedContours
def plotSelectedContours (imageEvery, stack, heartLayer, plot_all = False):
    """
    Function that plots all contours of a particular slice as an overlay on the image

    Parameters
    ----------
    imageEvery : int
        Number defining the interval in which slices will be plotted.
    stack : numpy array
        Stack being processed.
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates
    plot_all : boolean, optional
        True to show plot, else False. The default is False.

    Returns
    -------
    None.

    """
    if not plot_all:
        slc_tuple = (0,len(stack))
        slc_list, _ = getSlices(slc_tuple, 'you would like to check selected contours')
        # print(slc_list)
    else:
        slc_list = list(range(0,len(stack),1))

    index = 0
    for pos, keySlc in enumerate(heartLayer.keys()):
        if keySlc[0:3] == "slc":
            slcNum = int(keySlc[3:])
            if slcNum in slc_list:
                mult = index % imageEvery
                if mult == 0:
                    # print("Slc ", slcNum, keySlc)
                    myIm = stack[slcNum][:][:]
                    allContours = heartLayer[keySlc]["allContours"]
                    imIntFilledCont = heartLayer[keySlc]["imIntFilledCont"]
                    imExtFilledCont = heartLayer[keySlc]["imExtFilledCont"]
                    imAllFilledCont = heartLayer[keySlc]["imAllFilledCont"]

                    plotFilledCont(myIm, allContours, imIntFilledCont, imExtFilledCont, imAllFilledCont, True, slcNum)

                    index +=1

#%% func - plt_s3
def plt_s3 (start_slc, end_slc, im_every, s3_int, s3_ext, plotshow, option):
    """
    Function that plots the s3 array containing the images with filled contours

    Parameters
    ----------
    start_slc : int
        First slice
    end_slc : int
        Final slice
    im_every : int
        Number defining the interval in which slices will be plotted.
    s3_int : numpy array
        s3 array contianing internal contours
    s3_ext : numpy array
        s3 array contianing external contours.
    plotshow : boolean
        True to show plot, else False
    option : str
        "both ch"/"cardiacjelly"/"chamber ring" that changes the title to use on top of the subplots

    Returns
    -------
    None.

    """

    if plotshow:
        if option == "cardiacjelly":
            print('- Plotting ext myoc and int endo to get CJ (Plots tab)')
        elif option == "both ch":
            print('- Plotting myoc and ext endo to clean endo  (Plots tab)')
        elif option == "chamber ring":
            print('- Plotting heart layer and chamber ring to cut chambers  (Plots tab)')

        for slc in range(start_slc, end_slc-1, im_every):
            fig, ax = plt.subplots(1,2, figsize = (8,4))
            fig.suptitle("Filled contours for Slice "+ str(slc))
            ax[0].imshow(s3_int[:,:,slc+1])
            ax[0].set_xticks([])
            ax[0].set_yticks([])

            ax[1].imshow(s3_ext[:,:,slc+1])
            ax[1].set_xticks([])
            ax[1].set_yticks([])

            if option == "cardiacjelly":
                ax[0].set_title("Int.Myoc - ch1", fontsize=10)
                ax[1].set_title("Ext.Endoc - ch2", fontsize=10)

            elif option == "both ch":
                ax[0].set_title("Myocardium - ch1", fontsize=10)
                ax[1].set_title("Endocardium - ch2", fontsize=10)
                
            elif option == "chamber ring":
                ax[0].set_title("Heart tissue layer", fontsize=10)
                ax[1].set_title("Chamber ring", fontsize=10)


            plt.show()

#%% func - ch_clean_plt
def ch_clean_plt (mask_s3, toClean_s3, toRemove_s3, cleaned_s3, plotshow, im_every, option):
    """
    Function to plot mask, original image and result

    Parameters
    ----------
    mask_s3 : numpy array
        Stack of images to use as mask.
    toClean_s3 : numpy array
        Stack of images to clean.
    toRemove_s3 : numpy array
         Stack of images with regions to be removed.
    cleaned_s3 : numpy array
         Stack of cleaned images.
    plotshow : boolean
        True to show plot, else False
    im_every : int
        Number defining the interval in which slices will be plotted.
    option : str
        "clean"/"cardiacjelly" that changes the title to use on top of the subplots.

    Returns
    -------
    None.

    """
    if plotshow:
        if option == "cardiacjelly":
            print('- Extracting cardiac jelly ')
        elif option == "clean":
            print('- Cleaning endocardium ')

        for slc in range(0,toClean_s3.shape[2],im_every):

            #Get images
            mask_slc = mask_s3[:,:,slc]
            toClean_slc = toClean_s3[:,:,slc]
            toRemove_slc = toRemove_s3[:,:,slc]
            cleaned_slc = cleaned_s3[:,:,slc]

            #Plot
            fig, ax = plt.subplots(1, 4, figsize = (10,2.5))
            fig.suptitle("Slice:"+str(slc), y=1.05, weight="semibold")
            ax[0].imshow(mask_slc)
            ax[1].imshow(toClean_slc)
            ax[2].imshow(toRemove_slc)
            ax[3].imshow(cleaned_slc)
            if option == "clean":
                ax[0].set_title("ch0_inv")
                ax[1].set_title("ch1")
                ax[2].set_title("ch0_inv AND ch1")
                ax[3].set_title("(ch0_inv AND ch1)\nxOR ch1")
            elif option == "cardiacjelly":
                ax[0].set_title("ch0_int")
                ax[1].set_title("ch1_ext")
                ax[2].set_title("ch0_int AND ch1_ext")
                ax[3].set_title("Cardiac jelly")
            ax[0].set_xticks([])
            ax[0].set_yticks([])
            ax[1].set_xticks([])
            ax[1].set_yticks([])
            ax[2].set_xticks([])
            ax[2].set_yticks([])
            ax[3].set_xticks([])
            ax[3].set_yticks([])

            plt.show()

        alert('wohoo',1)

#%% - STACKS, SLICES, TUPLES, CLICKS
#%% func - createInitialClosedStack
def createInitialClosedStack (myStack, slices):
    """
    Funtion that initialises the closed stack between the slice range given as input by the user

    Parameters
    ----------
    myStack : numpy array
        Stack being processed.
    slices : tuple
        start_slice, end_slice

    Returns
    -------
    stack_closed : numpy array
        Resulting closed stack.

    """

    stack_closed = np.zeros_like(myStack, dtype='uint16')

    if slices[0] > 0:
        for slc in range(0,slices[0]):
            stack_closed[slc,:,:] = myStack[slc,:,:]
    if slices[1] < len(myStack):
        for slc in range(slices[1],len(myStack)):
            stack_closed[slc,:,:] = myStack[slc,:,:]

    return stack_closed

#%% func - getSlices
def getSlices (slc_tuple, text):
    """
    Funtion that returns a list with the slice numbers to be processed/checked.

    Parameters
    ----------
    slc_tuple : tuple of int
        Tuple defining the initial and final slice.
    text : str
        Additional text given as input to print corresponding to the process being performed.

    Returns
    -------
    slc_list : list of int
        List of slices to be processed/checked.
    slc_end : int
        Last slice of the tuple.

    """

    slc_list = []
    input_slc = ask4input('Select the slices '+text+' from -('+str(slc_tuple[0])+','+str(slc_tuple[1]-1)+')- \n\t\t(e.g. to close slices 5, 9, 10, and 11 type: 5,9-11)/[all/ ]:all/[N/n]:none): ', str)
    if input_slc == 'all' or input_slc == '':
        slc_list = list(range(slc_tuple[0],slc_tuple[1],1))
        # print(slc_list)
    elif input_slc in ['n', 'N']:
        slc_list = []
    else:
        slc_list = []
        comma_split = input_slc.split(',')

        for string in comma_split:
            if '-' in string:
                minus_split = string.split('-')
                #print(minus_split)
                for n in list(range(int(minus_split[0]),int(minus_split[1])+1,1)):
                    #print(n)
                    slc_list.append(n)
            elif string == '':
                pass
            else:
                slc_list.append(int(string))

    slc_end = slc_tuple[1]

    return slc_list, slc_end

#%% func - getSlicesContNum
def getSlicesContNum(stack):
    """
    Funtion that prompts the user the contour groups and the number of contours found in each contour group

    Parameters
    ----------
    stack : numpy array
        Stack being processed.

    Returns
    -------
    slcCont : list of int
        List of slice numbers where contour groups start.
    numCont : list of int
        List of contour numbers in each of the contour groups.

    """

    q_happy = False
    while not q_happy:
        slcCont = []
        numCont = []

        # q_slcCont = input('> Write the number of the first slice of each contour group separated by a comma\n(including 0 and len(stack)) >>>>> ')
        q_slcCont = ask4input('Write the number of the first slice of each contour group separated by a comma\n(include 0 and '+str(len(stack))+') >>>>> ', str)

        q_slcCont = q_slcCont.split(',')
        if int(q_slcCont[-1]) != len(stack):
            q_slcCont.append(str(len(stack)))

        for i in range(len(q_slcCont)):
            slcCont.append(int(q_slcCont[i]))
            if i >= 1:
                # q_numCont = input('> Number of contours found between slices '+q_slcCont[i-1]+ '-'+ q_slcCont[i]+': ')
                q_numCont = ask4input('Number of contours found between slices '+q_slcCont[i-1]+ '-'+ str(int(q_slcCont[i])-1)+': ', int)
                numCont.append(q_numCont)

        print('- Input revision:')
        for n, numC in enumerate(numCont):
            print('\t- Slices '+str(slcCont[n])+'-'+str(slcCont[n+1]-1)+':\t\t'+str(numC))

        q_happy = ask4input('Is this correct? [0]:no/[1]:yes: ', bool)

    return slcCont, numCont

#%% func - tuple_pairs
def tuple_pairs (numCont, slcCont, printshow, max_slc_diff):
    """
    Funtion that subdivides the contour groups into a maximum number of slices (max_slc_diff) and returns a new list
    of numCont and slcCont

    Parameters
    ----------
    numCont : list of int
        List of contour numbers in each of the contour groups.
    slcCont : list of int
        List of slice numbers where contour groups start.
    printshow : boolean
        True to print resulting tuples, else False
    max_slc_diff : int
        Maximum number of slices making up each contour group.

    Returns
    -------
    tuple_slc : list of tuples of int
        List with tuples including all slices.
    numContours_final :list of int
        Updated List of contour numbers in each of the contour groups.
    slcContours_final : list of int
        Updated list of slice numbers where contour groups start.

    """


    # len_numContours = len(numCont)
    numContours_final = []
    slcContours_final = []

    print('\n- slcContours:', slcCont)
    print('- numContours:', numCont)

    for num, slc in enumerate(slcCont):
        if num < len(slcCont)-1:
            first_tuple = slc
            end_tuple = slcCont[num+1]
            #print(num, slc)
            diff = end_tuple - first_tuple
            #print(diff)
            if diff <= max_slc_diff:
                numContours_final.append(numCont[num])
                slcContours_final.append(first_tuple)
            else:
                slcs2add = list(range(first_tuple,end_tuple, max_slc_diff))
                for slc2add in slcs2add:
                    numContours_final.append(numCont[num])
                    slcContours_final.append(slc2add)
        else:
            end_tuple = slc
            slcContours_final.append(end_tuple)

    print('- slcContours_final: len*',len(slcContours_final),'* -', slcContours_final)
    print('- numContours_final: len*',len(numContours_final),'* -', numContours_final)

    tuple_slc = [None]*len(numContours_final)

    for num in range(len(tuple_slc)):
        # if num == 0:
        #     tuple_slc[num]=(0, slcContours_py_final[num+1])
        # else:
        tuple_slc[num]=(slcContours_final[num], slcContours_final[num+1])

        if printshow == True:
            print("- tuple: ",tuple_slc[num])
            print("- range for that tuple: ", list(range(tuple_slc[num][0], tuple_slc[num][1])))

    return tuple_slc, numContours_final, slcContours_final

#%% func - update_tuple_pair
def update_tuple_pair(tuple_slc, numCont, heartLayer):
    """
    Funtion that promts the user to update (or not) the tuple pairs to fill dictionary

    Parameters
    ----------
    tuple_slc : list of tuples of int
        List with tuples including all slices.
    numCont : list of int
        List of contour numbers in each of the contour groups.
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates

    Returns
    -------
    tuple_slc_up : list of tuples of int
        Updated List with tuples including remaining slices.
    numCont_upf : list of int
        Updated List of contour numbers in the remaining contour groups.

    """

    try:
        q_startFromLast = ask4input('Do you want to continue selecting contours between slices -'+str(tuple_slc[0][0])+'-'+str(tuple_slc[-1][1])+' onwards? [0]:no/[1]:yes: ', bool)
        if not q_startFromLast:
            tuple_slc_up, numCont_upf = new_tuples_from_originals(heartLayer)
        else:
            tuple_slc_up = tuple_slc
            numCont_upf = numCont
    except:
        not_new = True
        print('- There are no tuples in the list.')
        while not_new:
            q_newTuples = ask4input('Do you want to select the contours for a new group of slices OR for individual slices? \n\t[0]:no, I am done selecting contours\n\t[1]:yes, group of slices\n\t[2]:yes, individual slices: ', int)
            if q_newTuples == 1:
                tuple_slc_up, numCont_upf = new_tuples_from_originals(heartLayer)
                not_new = False
            elif q_newTuples == 2:
                tuple_slc_up, numCont_upf = new_indiv_tuples(heartLayer)
                not_new = False
            elif q_newTuples == 0:
                tuple_slc_up = []
                numCont_upf = []
                not_new = False
            else:
                print('ERROR: -'+q_newTuples+'- The number entered needs to be a [0],[1] or [2]!')

    return tuple_slc_up, numCont_upf

#%% func - new_tuples_from_originals
def new_tuples_from_originals(heartLayer):
    """
    Function that returns a new group of tuples if user wants to fill the dictionary starting from a different slice

    Parameters
    ----------
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates

    Returns
    -------
    tuple_slc_up : list of tuples of int
        Updated List with tuples including remaining slices.
    numCont_upf : list of int
        Updated List of contour numbers in the remaining contour groups.

    """

    numCont_o = heartLayer['info']['numCont_group']
    numCont_o = [int(num) for num in numCont_o[1:-1].split(',')]
    slcCont_o = heartLayer['info']['slcCont_group']
    slcCont_o = [int(num) for num in slcCont_o[1:-1].split(',')]
    q_startFrom = ask4input('> Start group from slice number: ', int)
    q_EndAt = ask4input('> End group at slice number (inclusive): ', int)+1

    try:
        index_start = slcCont_o.index(q_startFrom)
        slcCont_up = slcCont_o[index_start:]
        numCont_up = numCont_o[index_start:]
    except:
        slcCont_add = slcCont_o.copy()
        slcCont_add.append(q_startFrom)
        slcCont_add.sort()
        index_start = slcCont_add.index(q_startFrom)
        slcCont_up = slcCont_add[index_start:]
        numCont_up = numCont_o[index_start-1:]
    # Cut end
    try:
        index_end = slcCont_up.index(q_EndAt)
        slcCont_up = slcCont_up[:index_end+1]
        numCont_up = numCont_up[:index_end]
    except:
        slcCont_add = slcCont_up.copy()
        slcCont_add.append(q_EndAt)
        slcCont_add.sort()
        index_end = slcCont_add.index(q_EndAt)
        slcCont_up = slcCont_add[:index_end+1]
        numCont_up = numCont_up[:index_end]

    tuple_slc_up, numCont_upf, _ = tuple_pairs(numCont_up, slcCont_up, False, 40) #HEREEE

    return tuple_slc_up, numCont_upf

#%% func - new_indiv_tuples_OLD
def new_indiv_tuples_old(heartLayer):
    """
    Function that returns a new group of tuples if user wants to fill the dictionary of individual slices

    Parameters
    ----------
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates

    Returns
    -------
    tuple_slc_new : list of tuples of int
        Updated List with tuples including remaining slices.
    numCont_new : list of int
        Updated List of contour numbers in the remaining contour groups.

    """
    numCont_o = heartLayer['info']['numCont_group']
    numCont_o = [int(num) for num in numCont_o[1:-1].split(',')]
    slcCont_o = heartLayer['info']['slcCont_group']
    slcCont_o = [int(num) for num in slcCont_o[1:-1].split(',')]

    q_slcCont_new = ask4input('Select the individual slices in which you would like to select the contours \n\t\t(e.g. to select contours of slices 15,20,21,22 type 15,20-22): ',str)
    slcCont_new = []
    comma_split = q_slcCont_new.split(',')
    for string in comma_split:
        if string == '':
            pass
        else:
            slcCont_new.append(int(string))

    numCont_new = []
    for slc in slcCont_new:
        try:
            index_start = slcCont_o.index(slc)
            numCont_new.append(numCont_o[index_start])
        except:
            slcCont_add = slcCont_o.copy()
            slcCont_add.append(slc)
            slcCont_add.sort()
            index_start = slcCont_add.index(slc)
            numCont_new.append(numCont_o[index_start-1])

    tuple_slc_new = [(slc, slc+1) for slc in slcCont_new]

    return tuple_slc_new, numCont_new

#%% func - new_indiv_tuples
def new_indiv_tuples(heartLayer):
    """
    Function that returns a new group of tuples if user wants to fill the dictionary of individual slices

    Parameters
    ----------
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates

    Returns
    -------
    tuple_slc_new : list of tuples of int
        Updated List with tuples including remaining slices.
    numCont_new : list of int
        Updated List of contour numbers in the remaining contour groups.

    """
    numCont_o = heartLayer['info']['numCont_group']
    numCont_o = [int(num) for num in numCont_o[1:-1].split(',')]
    slcCont_o = heartLayer['info']['slcCont_group']
    slcCont_o = [int(num) for num in slcCont_o[1:-1].split(',')]

    q_slcCont_new = ask4input('Select the individual slices in which you would like to select the contours \n\t\t(e.g. to select contours of slices 15,20,21,22 type 15,20-22): ',str)
    slcCont_new = []
    comma_split = q_slcCont_new.split(',')
    for string in comma_split:
        if '-' in string:
                minus_split = string.split('-')
                #print(minus_split)
                for n in list(range(int(minus_split[0]),int(minus_split[1])+1,1)):
                    #print(n)
                    slcCont_new.append(n)
        elif string == '':
            pass
        else:
            slcCont_new.append(int(string))

    numCont_new = []
    for slc in slcCont_new:
        try:
            index_start = slcCont_o.index(slc)
            numCont_new.append(numCont_o[index_start])
        except:
            slcCont_add = slcCont_o.copy()
            slcCont_add.append(slc)
            slcCont_add.sort()
            index_start = slcCont_add.index(slc)
            numCont_new.append(numCont_o[index_start-1])

    tuple_slc_new = [(slc, slc+1) for slc in slcCont_new]

    return tuple_slc_new, numCont_new

#%% func - getClicks
def getClicks (clicks, myIm, scale, text):
    """
    Function to get clicks of prompted image slice.

    Parameters
    ----------
    clicks : list of tuples with coordinates (clicks)
        Initialised list of clicks coordinates.
    myIm : numpy array
        Imaged being processed
    scale : float
        Number that scales the image to be prompted.
    text : str
        Title of prompted window with image indicating the process being performed (cleaning/closing/etc...).

    Returns
    -------
    clicks : list of tuples with coordinates (clicks)
        Final list of clicks coordinates.

    """

    print("- Getting clicks... Press ENTER when done")

    window_width = int(myIm.shape[1] * scale)
    window_height = int(myIm.shape[0] * scale)

    def on_mouse(event, x, y, flags, params):
        if event == cv2.EVENT_LBUTTONDOWN:
            #print ('\r    Seed: ' + str(y) + ', ' + str(x), myIm[y,x])
            clicks.append((y,x))

    cv2.namedWindow(text, cv2.WINDOW_NORMAL)#,cv2.WINDOW_NORMAL)
    cv2.resizeWindow(text, window_width, window_height)
    cv2.setMouseCallback(text, on_mouse, 0, )
    cv2.imshow(text, myIm)
    cv2.waitKey()
    cv2.destroyAllWindows()

    return clicks

#%% func - slcNum_def
def slcNum_def (slc):
    """
    Function to define the key for the slices within dictionary

    Parameters
    ----------
    slc : int
        Slice number.

    Returns
    -------
    sliceNum : str
        Key corresponding to slice number

    """

    # Creation of key for slice number
    if slc < 10:
        sliceNum = "slc00"+str(slc)
    elif slc < 100:
        sliceNum = "slc0"+str(slc)
    else:
        sliceNum = "slc"+str(slc)

    return sliceNum

#%% - CONTOURS CLEAN AND CLOSURE
#%% A. func - automCloseStackContours
def automCloseStackContours(myStack, ch, slices, new_stack, plotEvery, level):
    """
    Funtion that automatically closes the contours of the input slice between the range given by 'slices'

    Parameters
    ----------
    myStack : numpy array
        Stack being processed.
    ch : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    slices : tuple
        start_slice, end_slice
    new_stack : numpy array
        Initialised closed stack.
    plotEvery : int
        Number defining the interval in which slices will be plotted.
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    new_stack : numpy array
        Automatically closed stack.
    done_autom : boolean
        True if automatic closure of contours has already been performed, else False.

    """

    print('\n')
    bar = Bar('Automatically closing contours', max=slices[1]-slices[0], suffix = suffix, check_tty=False, hide_cursor=False)
    for index, slc in enumerate(range(slices[0], slices[1])):

        plotshow = False
        mult = index % plotEvery
        if mult == 0:
            plotshow = True
            print("\n------------- Channel "+str(ch)+" / Slice "+str(slc)+" -------------")

        myIm = myStack[slc][:][:]
        # Get the contours of slice
        contours, numCont = getContExpCont (myIm, minLenContour = 100, level = level)
        # Sort the contours by length (bigger to smaller)
        contours = sorted(contours, key = len, reverse=True)
        # Get properties of each contour
        props = getContProps(myIm, ch, slc, cont_sort = contours, num_contours = numCont, plotshow = False)
        # Select only the contours that have a max intensity greater than min_expInt
        filt_cont, filt_props = filterContours(contours, props, min_expInt = 15000, printData = plotshow)
        # Get properties of filtered contours
        # if plotshow:
        #     print('\n-> FILTERED CONTOURS')
        _ = getContProps(myIm, ch, slc, cont_sort = filt_cont, num_contours = len(filt_cont), plotshow = plotshow)
        # Find distance between all contours and save information of those whose are at a distance less than minDist
        data2Connect = distBtwAllCont (contours = filt_cont, minDist = 15, printData = plotshow)
        # Draw lines between closer contours
        myIm_closedCont = automCloseContours (myIm, data2Connect)
        # Show new closed contours
        new_contours, new_numCont = getContExpCont (myIm_closedCont, minLenContour = 100, level = level)
        new_contours = sorted(new_contours, key = len, reverse=True)
        # if plotshow:
        #     print('\n-> FINAL CONTOURS')
        _ = getContProps(myIm_closedCont, ch, slc, cont_sort = new_contours, num_contours = new_numCont, plotshow = plotshow)

        new_stack[slc][:][:] = myIm_closedCont
        bar.next()

    bar.finish()
    print('- FINISHED Automatic closure of contours')
    done_autom = True
    alert('wohoo',1)

    return new_stack, done_autom

#%% func - getContExpCont
def getContExpCont (myIm, minLenContour, level):
    """
    Function that gets and returns the contours of a particular slice (slcNum)

    Parameters
    ----------
    myIm : numpy array
        Imaged being processed
    minLenContour : int
        Minimum number of points making up contoura and used to filter contours found in each slice
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    arr_contour_slc : list of arrays
        list of numpy array with coordinates making up the contours of the input slice.
    num_contour : int
        Number of contours found in slice.

    """

    # Create an empty array to save all the contours of each slice individually
    arr_contour_slc = []

    # Find all the contours of the image
    contours = measure.find_contours(myIm, level, 'high', 'high')

    # Variable to save the number of contours found
    num_contour = 0

    # Go through all the contours
    for n, contour in enumerate(contours):
        # Get only the contours made up of more than the designated number of points
        if len(contour)>minLenContour:
            # Append contour to the array
            arr_contour_slc.append(contour)
            num_contour += 1

    return arr_contour_slc, num_contour

#%% func - getContProps
def getContProps (myIm, ch, slc, cont_sort, num_contours, plotshow):
    """
    Funtion that returns the properties of the contours given as input

    Parameters
    ----------
    myIm : numpy array
        Imaged being processed
    ch : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    slc : int
        Slice number being processed.
    cont_sort : list of arrays
        list of arrays with contours' coordinates sorted from longer to shorter.
    num_contours : int
        Number of contours found in slice.
    plotshow : boolean
        True to show plot, else False

    Returns
    -------
    props_all : list of numpy array
        List of arrays with values of properties associated to each of the contours found in slice
        -->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]

    """

    if plotshow and num_contours > 0:
        #Print the channel and slice being analysed
        #print("Channel "+str(ch)+" / Slice "+str(slc))

        # Define the figure properties (columns, rows, image size)
        cols = 5
        # Limit the number of contours to plot to 30
        if num_contours > 30:
            num_contours = 30
            cont_plot = cont_sort[0:30]
        else:
            cont_plot = cont_sort

        rows = num_contours // cols
        if num_contours%cols != 0:
            rows = rows + 1
        imSize = 3
        colorImSize = 4
        fig11 = plt.figure(figsize=(cols*imSize+colorImSize, rows*imSize), constrained_layout=True)

        # gridspec inside gridspec
        outer_grid = fig11.add_gridspec(1,2)
        # Grid where color image will be placed
        color_grid = outer_grid[0].subgridspec(1,1, wspace=0, hspace=0)
        ax = fig11.add_subplot(color_grid[0])
        ax.imshow(myIm, cmap=plt.cm.gray)

        # Go through all the contours
        for num, contour in enumerate(cont_plot):
            ax.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color=colors[num])
            txt = "Cont"+str(num)
            ax.text(0.95,(0.97-0.035*(num+1)), txt,
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes,
                        color=colors[num], fontsize=10, weight = 'semibold')

        ax.set(xlabel = "Channel "+str(ch)+" / Slice "+str(slc))

        # Grid where subplots of each contour will be placed
        all_grid = outer_grid[1].subgridspec(rows, cols,  wspace=0.3, hspace=0.3)

    # Array to save the metrics of all the contours
    props_all = [None]*len(cont_sort)
    # Iterate through sorted contours
    for index, contList in enumerate(cont_sort):

        #-->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]
        props = maskContour(myIm, contList)
        props_all[index] = props

        area = format(props[0], '.0f')
        max_int = format(props[2], '.0f')
        mean_int = format(props[3], '.2f')
        lgth = format(props[4],'.0f')
        # per = format(props[5], '.0f')
        # sol = format(props[6], '.2f')

        if plotshow and index < 30 :
            ax = fig11.add_subplot(all_grid[index])
            ax.imshow(myIm, cmap=plt.cm.gray)
            ax.plot(contList[:, 1], contList[:, 0], linewidth=1.5, color = colors[index])
            ax.set_title("Contour "+str(index), fontsize=10, weight = 'semibold', color = colors[index])
            ax.set(xlabel = "Area:"+str(area)+" / L:"+str(lgth), ylabel = "Max:"+str(max_int)+" / Mean:"+str(mean_int))
            ax.set_xticks([])
            ax.set_yticks([])

    if plotshow:
        plt.show()

    return props_all

#%% func - maskContour
def maskContour (myIm, contour):
    """
    Function that measures the properties of the contour received as input

    Parameters
    ----------
    myIm : numpy array
        Imaged being processed
    contour : numpy array
        Coordinates of contour to use to mask image

    Returns
    -------
    props_all : list of numpy array
        List of arrays with values of properties associated to each of the contours found in slice
        -->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]

    """

    # Create an empty image to store the masked array
    r_mask = np.zeros_like(myIm, dtype='bool')
    # Create a contour masked image by using the contour coordinates rounded to their nearest integer value
    r_mask[np.round(contour[:, 1]).astype('int'), np.round(contour[:, 0]).astype('int')] = 1
    # Fill in the holes created by the contour boundary
    r_mask = ndimage.binary_fill_holes(r_mask)

    # Change the mask type to integer
    r_mask_int = r_mask.astype(np.int)

    # label image regions
    label_r_mask = label(r_mask_int)
    label_r_mask = np.transpose(label_r_mask)
    # print(label_r_mask.shape, myIm.shape)
    
    props = regionprops(label_r_mask, intensity_image = myIm)
    
    area = props[0].area
    centroid = props[0].centroid
    max_int = props[0].max_intensity
    mean_int = props[0].mean_intensity
    #cx_area = props[0].convex_area
    #ecc = props[0].eccentricity
    lgth = len(contour)
    per = props[0].perimeter
    sol = props[0].solidity
    bbox = props[0].bbox

    props_all = np.array([area, centroid, max_int, mean_int, lgth, per, sol, bbox], dtype=object)

    return props_all

#%% func - filterContours
def filterContours (contours, props, min_expInt, printData):
    """
    Funtion that filters contours by minimum and mean intensities

    Parameters
    ----------
    contours : list of arrays
        list of numpy array with coordinates making up the contours of the slice being processed.
    props : list of numpy array
        List of arrays with values of properties associated to each of the contours found in slice
        -->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]
    min_expInt : int
        Minimum expected intensity value.
    printData : boolean
        True to print number of initial and final contours, else False

    Returns
    -------
    filt_cont : list of arrays
        Filtered list of numpy array with coordinates making up the contours of the slice being processed.
    filt_props : list of numpy array
        Filtered list of arrays with values of properties associated to each of the filtered contours found in slice
        -->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]

    """

    filt_cont = []
    filt_props = []
    for num, cont in enumerate(contours):
        if props[num][2] > min_expInt and props[num][3] > 5000:
            filt_cont.append(cont)
            filt_props.append(props[num])

    if printData:
        print('- Number of initial contours: ', len(contours))
        print('- Number of final contours: ', len(filt_cont))

    return filt_cont, filt_props

#%% func - distBtwAllCont
def distBtwAllCont (contours, minDist, printData):
    """
    Function used to find minimum distance between a list of input contours

    Parameters
    ----------
    contours : list of arrays
        list of numpy array with coordinates making up the contours of the slice being processed.
    minDist : int
        Minimum distance in pixels
    printData : boolean
        True to print the number of the contours to connect, else False

    Returns
    -------
    final_list : list of lists
        List with data to be used to join contours [i, j+1, min_dist, indA, ptA, indB, ptB].

    """

    mtx_dist = np.ones(shape = (len(contours), len(contours)))*10000
    np.fill_diagonal(mtx_dist, 0)

    ij_list = []
    final_list = []

    for i, contA in enumerate(contours):
        # print('i:',i)
        for j, contB in enumerate(contours[i+1:]):
            # print('j:',j+1)
            #Get minimum distance between contA and contB
            data_func = minDistBtwContoursAB (contA, contB)
            min_dist, indA, ptA, indB, ptB = data_func
            data_func = i, j+1, min_dist, indA, ptA, indB, ptB
            mtx_dist[i,i+j+1] = min_dist

            if min_dist < minDist:
                ij_tuple = (i,i+j+1)
                ij_list.append(ij_tuple)
                final_list.append(data_func)
                if printData:
                    print('- Contours to connect: ', ij_tuple, ' / Distance [px]: ', format(min_dist, '.2f'))
                    # print('Points coordinates: ', ptA, '/', ptB)

    if printData:
        print('')
    return final_list#, ij_list, np.triu(mtx_dist)

#%% func - minDistBtwContoursAB
def minDistBtwContoursAB (contourA, contourB):
    """
    Function that measues the actual distance in pixels between contourA and contourB

    Parameters
    ----------
    contourA : numpy array
        Array with coordinates of contour A.
    contourB : numpy array
        Array with coordinates of contour B.

    Returns
    -------
    min_dist : float
        Minimum distance found between contours A and B.
    index_ptContA : int
        Index of point within array of contour A that results in minimum distance
    ptContA : numpy array
        Coordinates of point A
    index_ptContB : int
        Index of point within array of contour B that results in minimum distance
    ptContB : numpy array
        Coordinates of point B

    """

    # Find euclidian distance beteen all points in contours
    dist = cdist(contourA,contourB)
    # Get indexes where distance is minimum (closest points)
    min_dist = np.amin(dist)
    index4min = np.where(dist == min_dist)

    index_ptContA = index4min[0][0]
    ptContA = contourA[index_ptContA]
    index_ptContB = index4min[1][0]
    ptContB = contourB[index_ptContB]

    return min_dist, index_ptContA, ptContA, index_ptContB, ptContB

#%% func - automCloseContours
def automCloseContours (myIm, data2Connect):
    """
    Function that connects contours using the data2Connect given as input - (return of distBtwAllCont)

    Parameters
    ----------
    myIm : numpy array
        Imaged being processed
    data2Connect  :list of lists
        List with data to be used to join contours [i, j+1, min_dist, indA, ptA, indB, ptB].

    Returns
    -------
    myIm : numpy array
        Resulting processed imaged with connected contours

    """

    for num, connection in enumerate(data2Connect):
        ptAx = int(connection[4][0])
        ptAy = int(connection[4][1])
        ptBx = int(connection[6][0])
        ptBy = int(connection[6][1])

        rr, cc, val = line_aa(ptAx, ptAy, ptBx, ptBy)
        myIm[rr, cc] = val * 50000

    return myIm

#%% B. func - manuallyCloseContours
def manuallyCloseContours (stack_closed, stack_o, stack_m, slices, n_rows, chStr, exit_code, level):
    """
    Function used to manually close the contours of the stack

    Parameters
    ----------
    stack_closed : numpy array
        Input closed stack.
    stack_o : numpy array
        Copy of closed stack.
    stack_m : numpy array
        Copy of original masked stack to process.
    slices : tuple
        start_slice, end_slice
    n_rows : int
        Number of rows in subplot.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    exit_code : boolean
        True if user wants to exit code, else False.
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    stack_closed : numpy array
        Resulting closed stack.
    last_slc :  int
        Last processed slice
    exit_code : boolean
        True if user wants to exit code, else False.

    """
    print('\n- Closing manually slices: ',slices[0], '-',slices[1]-1)

    exit_code = False
    slcs_per_im = n_rows*4

    slices_first = list(range(slices[0],slices[1]+1,slcs_per_im))
    # print('slices_first: ',slices_first)
    slices_last =  list(range(slices[0]+slcs_per_im,slices[1]+1,slcs_per_im))
    # print('slices_last: ',slices_last)

    if slices_last != slices[-1]:
        slices_last.append(slices[-1])

    # print('slices_first:', shape(slices_first))

    for i in range(len(slices_first)):#shape(slices_first)[0]):
        slc_tuple = (slices_first[i], slices_last[i])
        plotSlcsRange(stack_closed, slc_tuple, 'INITIAL', slcs_per_im, level)
        slc_list, slc_end = getSlices(slc_tuple, 'you would like to close')
        # print('AJA - slc_list:', slc_list)

        while len(slc_list) != 0:
            exit_code, slc_end, last_slc, stack_closed = manuallyCloseContoursTuple (slc_list, stack_closed, stack_o, 
                                                                                     stack_m, chStr, exit_code, level)

            if exit_code:
                alert('error',1)
                print("\n- Exit script - last slice ", last_slc)
                break

            plotSlcsRange(stack_closed, slc_tuple, 'CHECKING (after having closed)', slcs_per_im, level)
            q_done = str(input('> Are you done CLOSING the contours for this - tuple ('+ str(slc_tuple[0])+','+str(slc_tuple[1]-1)+')?: \n - [0]:no/[1/ ]:yes! >>>>> ')).lower()
            if q_done == '1' or q_done == '':
                # print('slc list: ', slc_list)
                break
            #     if i == len(slices_first)-1:
            #         q_not_done_all = str(input('> Do you want to close any other slices? [0]:no/[1/ ]:yes! >>>>> '))
            #         if q_not_done_all == '1' or q_not_done_all == '': 
            #             slc_list, slc_end = getSlices((0,stack_closed.shape[0]), 'you would like to close')
            #         else: 
            #             exit_code = True
            #             break
            #     else: 
            #         break
            elif q_done == '0':
                slc_list, slc_end = getSlices(slc_tuple, 'you would like to close')
            
        if exit_code:
            break

    if slc_end == slices[1]:
        last_slc = slices[1]
        # alert('wohoo',1)
        # print("- All Done - Contours have been manually closed!")

    return stack_closed, last_slc, exit_code

#%% func - manuallyCloseContoursTuple
def manuallyCloseContoursTuple (slc_list, stack_closed, stack_o, stack_m, chStr, exit_code, level):
    """
    Function used to manually close the contours of the slices in the slc_list given as input

    Parameters
    ----------
    slc_list : list of int
        List of slices to be processed/checked.
    stack_closed : numpy array
        Input closed stack.
    stack_o : numpy array
        Copy of closed stack.
    stack_m : numpy array
        Copy of original masked stack to process.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    exit_code : boolean
        True if user wants to exit code, else False.
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    exit_code : boolean
        True if user wants to exit code, else False.
    slc : int
        Slice number being processed.
    last_slc : int
        Last processed slice
    stack_closed : numpy array
        Resulting closed stack.

    """
    select_process = 'esc'
    if len(slc_list) > 0: 
        for slc in slc_list:
            if exit_code:
                print('- EXIT! - slc: '+ str(slc)+' - Process: '+ select_process)
                break
    
            print("- Processing contours - slice ", str(slc))
            # Get image
    
            myIm = stack_closed[slc][:][:]
            # Get contours for that slice
            _ = getContExpCont_plt (myIm, slc, chStr, minLenContour=100, figSize=10, level = level)
            myIm_closed = myIm
    
            # closeContours per slice
            k_size = 60
            myIm_closed, done_close  = closeContoursSlc(myIm_closed, slc, chStr, k_size, k_size, level)
            if done_close == 'esc': #done_close in ['esc', 'n', 'N']
                last_slc = slc
                exit_code = True
                break
    
            last_slc = slc
            add_step = ['', '1']
            while not exit_code and done_close in add_step:   #exit_code == False:
                list_close = ['3', '4', '5']
                crop_size = [(120,120), (50,50), (100,300)]
                print('\n - Additional processes for - Slc '+str(slc)+':')
                print('   -[1]:draw black/clean slice\t\t-[2]:draw white')
                print('   -[3]:close (120x120)\t\t\t\t-[4]:close (50x50)')
                print('   -[5]:close (300x100)\t\t\t\t-[6]:close (user input size)')
                print('   -[7]:reset slice (to original after automatic contour closure)')
                print('   -[8]:reset slice (to original image - just masked)')
                print('   -[9/ ]:save (slice done!)\t\t-[esc]:exit')
                select_process = str(input('> Select: ')).lower()
    
                if select_process == 'esc':
                    done = select_process
                    last_slc = slc;
                    exit_code = True
                    break
                # Draw contour black/Clean contours
                elif select_process == '1':
                    myIm_closed, done = close_draw(myIm_closed, slc, chStr, 'black', level)
                # Draw contour white
                elif select_process == '2':
                    myIm_closed, done = close_draw(myIm_closed, slc, chStr, 'white', level)
                # Close contours
                elif select_process in list_close:
                    index_selected = list_close.index(select_process)
                    myIm_closed, done = closeContoursSlc(myIm_closed, slc, chStr, crop_size[index_selected][0], crop_size[index_selected][1], level)
                elif select_process == '6':
                    input_size = ask4input('Enter the rectangle size to crop and close contours separated by a comma - x,y -: ', str)
                    x_size, y_size = input_size.split(',')
                    myIm_closed, done = closeContoursSlc(myIm_closed, slc, chStr, int(y_size), int(x_size), level)
                # Done
                elif select_process == '9' or select_process == '':
                    done = 'OK'
                    stack_closed[slc][:][:] = myIm_closed
                    last_slc = slc
                    break
                #Reset slc
                elif select_process == '7':
                    myIm_closed = stack_o[slc][:][:]
                    done = 'reset_o_OK'
                    _ = getContExpCont_plt (myIm_closed, slc, chStr, 250, 10, level, True)
                    
                elif select_process == '8':
                    myIm_closed = stack_m[slc][:][:]
                    done = 'reset_m_OK'
                    _ = getContExpCont_plt (myIm_closed, slc, chStr, 250, 10, level, True)
                    
                else:
                    print('Error: Wrong input. Select only from the given options')
    
                if done == 'esc':
                    last_slc = slc
                    exit_code = True
                    break
    else: 
        exit_code = False
        last_slc = 0
        slc = 0

    return exit_code, slc, last_slc, stack_closed

#%% func - closeContoursSlc
def closeContoursSlc(myIm_closed, slc, chStr, kw, kh, level):
    """
    Function that closes the contours of the input image using the cropNcloseCont function

    Parameters
    ----------
    myIm_closed : numpy array
        Image to close
    slc : int
        Slice number being processed.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    kw : int
        width of the rectangular region to crop around the clicked point.
    kh : int
        length of the rectangular region to crop around the clicked point.
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    myIm_closed : numpy array
        Closed image
    done_close : boolean
        True if slice has already been closed, else False.

    """
    options_break = ['1','2','esc','']
    #Get clicks of positions to close contours
    while True:
        clicks = getClicks([], myIm_closed, scale=1, text='CLOSING CONTOURS')
        #Close contours and get Image
        myIm_closed = cropNcloseCont (clicks, myIm_closed, kw, kh, level)
        #Plot image with closed contours
        _ = getContExpCont_plt (myIm_closed, slc, chStr, 250, 10, level)
        #Ask if done
        # print('\nAre you done CLOSING the contours for this -slice '+ str(slc)+'?: \n')
        # done_close = input('> Are you done CLOSING the contours for this -slice '+ str(slc)+'?: \n - [0]: no, [1/ ]: yes!, [esc]: exit >>>>> ')
        done_close = str(input('> Are you done CLOSING the contours for this -slice '+ str(slc)+'?: \n - [0]:no/[1/ ]:yes, but I would like to do something else to it/[2]:yes, continue with the next slice!/[esc]:exit >>>>> ')).lower()
        if done_close in options_break:
            break

    return myIm_closed, done_close

#%% func - cropNcloseCont
def cropNcloseCont (clicks, myIm, wh, ww, level, plot_show = True):
    """
    Function that crops a rectangular region of the input image centered by the click, gets the contours present in that
    region and connects the two closest contours with a white line

    Parameters
    ----------
    clicks : list of tuples with coordinates (clicks)
        Initialised list of clicks coordinates to be used to close image.
    myIm : numpy array
        Imaged being processed
    wh :  int
        length of the rectangular region to crop around the clicked point.
    ww : int
        width of the rectangular region to crop around the clicked point.
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours
    plot_show : boolean, optional
        True to show plot, else False. The default is True.

    Returns
    -------
    myIm : numpy array
        Resulting processed imaged

    """
    wh = wh//2
    ww = ww//2

    for n, click in enumerate(clicks):
        #print('click:',click)
        y0, x0 = click

        #Crop image in a square with center: click
        ymin = y0-wh
        xmin = x0-ww
        if ymin < 0:
            ymin = 0
        if xmin < 0:
            xmin = 0
        imCrop = myIm[ymin:y0+wh,xmin:x0+ww]

        if plot_show:
            fig, ax = plt.subplots(1,3, figsize=(10,3))
            ax[0].imshow(imCrop, cmap=plt.cm.gray)
            ax[0].set_title("Original")

        #Get contours of the cropped image
        contours_crop = measure.find_contours(imCrop,level, 'high', 'high')
        #Organize contours in terms of number of points to get the two biggest
        contours_crop = sorted(contours_crop, key = len, reverse=True)
        
        if len(contours_crop) > 1:
            #Find euclidian distance beteen all points in contours
            dist = cdist(contours_crop[0],contours_crop[1])
            #Get indexes where distance is minimum (closest points)
            index4min = np.where(dist == np.amin(dist))

            index_ptCont0 = index4min[0][0]
            ptCont0 = contours_crop[0][index_ptCont0]
            index_ptCont1 = index4min[1][0]
            ptCont1 = contours_crop[1][index_ptCont1]

            if plot_show:
                ax[1].imshow(imCrop, cmap=plt.cm.gray)
                ax[1].plot(ptCont0[1], ptCont0[0],'r.')
                ax[1].plot(ptCont1[1], ptCont1[0],'g.')
                for n, contour in enumerate(contours_crop):
                    ax[1].plot(contour[:, 1], contour[:, 0], linewidth=1)
                ax[1].set_title("Cont and Pts")

            rr, cc, val = line_aa(int(ptCont0[0]), int(ptCont0[1]), int(ptCont1[0]), int(ptCont1[1]))
            imCrop[rr, cc] = val * 50000

            if plot_show:
                ax[2].imshow(imCrop, cmap=plt.cm.gray)
                ax[2].plot(ptCont0[1], ptCont0[0],'r.')
                ax[2].plot(ptCont1[1], ptCont1[0],'g.')
                ax[2].set_title("Closed Cont")
                plt.show()

    return myIm

#%% func - close_draw
def close_draw (myIm_closed, slc, chStr, color_draw, level, plot_show = True):
    """
    Function that collects clicks positions given by the user and connects them using a white or black line given as
    input by 'color_draw' parameter

    Parameters
    ----------
    myIm_closed : numpy array
        Image to close drawing
    slc : int
        Slice number being processed.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    color_draw : str
        Color used when drawing line between clicks ('b'/'black': black, 'w'/'white': white)
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours
    plot_show : boolean, optional
        True to show plot, else False. The default is True.

    Returns
    -------
    myIm_closed : numpy array
        Closed/cleaned image
    done_draw : boolean
        True if drawing is done, else False.

    """

    while True:
        #Get clicks of positions to close contours
        clicks = getClicks([], myIm_closed, scale=1, text='DRAWING SLICE ('+color_draw+')')
        # Draw white/black line following the clicked pattern
        myIm_closed = drawLine(clicks, myIm_closed, color_draw)
        _ = getContExpCont_plt (myIm_closed, slc, chStr, 250, 10, level, plot_show)

        if plot_show:
            print("- Are you done drawing ("+color_draw+") the contours for this - slice "+ str(slc)+"?: ")
            done_draw = str(input('>  - [0]:no/[1/ ]:yes!/[esc]:exit >>>>> ')).lower()
            if done_draw == '1' or done_draw == 'esc' or done_draw == '':
                break
        else:
            done_draw = ''
            break

    return myIm_closed, done_draw

#%% func - cleanSlice - REMOVE?
def cleanSlice (myIm_closed, slc, chStr, level):
    """
    Function that collects clicks positions given by the user and connects them using a black line

    Parameters
    ----------
    myIm_closed : numpy array
        Image to clean
    slc : int
        Slice number being processed.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    myIm_closed : numpy array
        Cleaned image
    done_clean : boolean
        True if cleaning is done, else False.

    """

    _ = getContExpCont_plt (myIm_closed, slc, chStr, 250, 7, level)
    while True:
        #Get clicks of positions to close contours
        clicks = getClicks([], myIm_closed, scale=1, text='CLEANING SLICE')
        # Draw black line following the clicked pattern
        myIm_closed = drawLine(clicks, myIm_closed, "b")
        _ = getContExpCont_plt (myIm_closed, slc, chStr, 250, 7, level)
        #Ask if done
        print("- Are you done CLEANING the contours for this -slice "+ str(slc)+"?: ")
        done_clean = str(input('>  - [0]:no/[1/ ]:yes!/[esc]:exit >>>>> ')).lower()
        if done_clean == '1' or done_clean == 'esc' or done_clean == '':
            break

    return myIm_closed, done_clean

#%% func - drawLine
def drawLine (clicks, myIm, color_draw):
    """
    Function that draws white or black line connecting all the clicks received as input

    Parameters
    ----------
    clicks : list of tuples with coordinates (clicks)
        Initialised list of clicks coordinates to be used to draw.
    myIm : numpy array
        Imaged being processed
    color_draw : str
        Color used when drawing line between clicks ('b': black, 'w': white)

    Returns
    -------
    myIm : numpy array
        Resulting processed imaged

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



#%% C. func - closeInfOutfStack
def closeInfOutfStack (stack_closed, slices, chStr, exit_code, region, level):
    """
    Function used to close the inflow/outflow tract of the stack

    Parameters
    ----------
    stack_closed : numpy array
        Input closed stack.
    slices : tuple
        start_slice, end_slice
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    exit_code : boolean
        True if user wants to exit code, else False.
    region : str
        Text indicating the region being closed ('Inflow'/'Outflow').
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    stack_closed : numpy array
        Resulting closed stack.
    slices : tuple
        start_slice, end_slice
    exit_code : boolean
        True if user wants to exit code, else False.

    """

    stack_o = np.copy(stack_closed)
    # Starting and ending slices to close inflow or outflow tracts
    start_InO, end_InO = slices
    # Get info about point every X images
    list_slc = list(range(start_InO, end_InO, 12))

    bar = Bar('Closing '+region+' Contours', max=len(list_slc), suffix = suffix, check_tty=False, hide_cursor=False)
    for slc in list_slc:
        happy = False
        print(slc)
        if slc != list_slc[-1]:
            slc_end = slc+12
        else:
            slc_end = end_InO
        print(slc_end)
        while not happy: #happy == False:
             stack_closed, q_happy, last_slc, exit_code = closeInfOutfTuple(stack_closed, (slc,slc_end), chStr, exit_code, level)
             if not q_happy and exit_code == False:
                 #Reset slices
                 stack_closed[slc:slc_end][:][:] = stack_o[slc:slc_end][:][:]
             elif exit_code == True:
                 #Exit and save last slice
                 break
             elif q_happy:
                 happy = True
                 bar.next()
        if exit_code:
            break
    bar.finish()

    if exit_code:
        alert('error',1)
        print("- Exit script - last slice", last_slc)
        last_slc, end_InO = slices
    else:
        end_InO, end_InO = slices

    return stack_closed, slices, exit_code

#%% func - closeInfOutfTuple
def closeInfOutfTuple (stack_closed, slices, chStr, exit_code, level):
    """
    Function used to close the inflow/outflow tract of the slices within the range of 'slices' given as input

    Parameters
    ----------
    stack_closed : numpy array
        Input closed stack.
    slices : tuple
        start_slice, end_slice
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    exit_code : boolean
        True if user wants to exit code, else False.
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    stack_closed : numpy array
        Resulting closed stack.
    q_happy : boolean
        True if the user considers successful the closure of the inflow/outflow tract in the slice range, else False
    last_slc : int
        Last processed slice
    exit_code : boolean
        True if user wants to exit code, else False.

    """
    
    print('- Closing inflow/outflow tract for slices: '+str(slices)+'...')
    for slc in range(slices[0], slices[1], 1):
        # print(slc)
        # Get image
        myIm = stack_closed[slc][:][:]
        myIm_closed = myIm
        done = False
        # Run method selection for every XX slices in range
        if slc == slices[0]:
            while not done:# and exit_code: 
                _ = getContExpCont_plt (myIm, slc, chStr, 600,7, level)
                method = str(input("\n> What method do you want to use for closing the inflow/outflow tract?\n\t[1/ ]:ConvexHull\n\t[2]:Close (450x150)\n\t[3]:Close (500x200)\n\t[4]:Close (200x50)\n\t[5]:Close (user input size)\n\t[6]:Draw\n\t[esc]:exit >>>>> : ")).lower()
                if method == "1" or method == '':
                    # Use convexHull method to close inflow/outflow tract
                    alert('error',1)
                    print('- IMPORTANT NOTE: Remember to click a point outside the convex hull of the contours!')
                    myIm_closed, clicks_chull = close_convexHull(myIm_closed, slc, chStr, clicks_in = [], option = 'first', level = level)
                    done = True
                elif method == "2":
                    # Use rectangle method to close inflow/outflow tract
                    rect_xA = 150
                    rect_yA = 450
                    myIm_closed, clicks_close = closeInfOutSlc(myIm_closed, slc, chStr, rect_xA, rect_yA, clicks_in = [], option = 'first', level = level)
                    done = True
                elif method == "3":
                    # Use rectangle method to close inflow/outflow tract
                    rect_xC = 200
                    rect_yC = 500
                    myIm_closed, clicks_close = closeInfOutSlc(myIm_closed, slc, chStr, rect_xC, rect_yC, clicks_in = [], option = 'first', level = level)
                    done = True
                elif method == "4":
                    # Use rectangle method to close inflow/outflow tract
                    rect_xB = 50
                    rect_yB = 200
                    myIm_closed, clicks_close = closeInfOutSlc(myIm_closed, slc, chStr, rect_xB, rect_yB, clicks_in = [], option = 'first', level = level)
                    done = True
                elif method == "5":
                    # Use rectangle method to close inflow/outflow tract
                    input_size = ask4input('Enter the rectangle size to crop and close contours separated by a comma - x,y -: ', str)
                    rect_xD, rect_yD = input_size.split(',')
                    myIm_closed, clicks_close = closeInfOutSlc(myIm_closed, slc, chStr, int(rect_yD), int(rect_xD), clicks_in = [], option = 'first', level = level)
                    done = True
                elif method == "6":
                    myIm_closed, done = close_draw(myIm_closed, slc, chStr, 'w', level, plot_show = False)
                    done = True
                elif method == 'esc':
                    last_slc = slc
                    exit_code = True
                    q_happy = False
                    break
                else:
                    print('Error: Wrong input. Select only from the given options')

        # Use the previous selected method for the rest of the XX-1 slices
        else:
            if method == "1" or method == '':
                myIm_closed, _ = close_convexHull(myIm_closed, slc, chStr, clicks_in = clicks_chull, option = 'autom', level = level)
            elif method == "2":
                myIm_closed, _ = closeInfOutSlc(myIm_closed, slc, chStr, rect_xA, rect_yA, clicks_in = clicks_close, option = 'autom', level = level)
            elif method == "3":
                myIm_closed, _ = closeInfOutSlc(myIm_closed, slc, chStr, rect_xC, rect_yC, clicks_in = clicks_close, option = 'autom', level = level)
            elif method == "4":
                myIm_closed, _ = closeInfOutSlc(myIm_closed, slc, chStr, rect_xB, rect_yB, clicks_in = clicks_close, option = 'autom', level = level)
            elif method == "5":
                myIm_closed, _ = closeInfOutSlc(myIm_closed, slc, chStr, int(rect_yD), int(rect_xD), clicks_in = clicks_close, option = 'autom', level = level)
            elif method == "6":
                myIm_closed, done = close_draw(myIm_closed, slc, chStr, 'white', level, plot_show = False)

        if slc == slices[1]-1:
            #print('I am the last slice of the group')
            showGridContours(myStack = stack_closed, slices = (slices[0], slices[1]), n_rows = 3, level = level)
            q_happy = ask4input('Are you happy with the way the contours have been closed? \n\t [0]:no, let me select other method to close them\n\t [1]:yes, continue!: ',bool)
            last_slc = slices[0]

        if exit_code:
            break

        # Add the closed image to the stack
        stack_closed[slc][:][:] = myIm_closed

    return stack_closed, q_happy, last_slc, exit_code

#%% func - close_convexHull
def close_convexHull(myIm_closed, slc, chStr, clicks_in, option, level):
    """
    Function that closes the inflow/ouflow tract of the input slice using convex hull

    Parameters
    ----------
    myIm_closed : numpy array
        Image to close using convex hull
    slc : int
        Slice number being processed.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    clicks_in : list of tuples with coordinates (clicks)
        Initialised list of clicks coordinates to close the inflow/outflow with convex Hull.
    option : str
        Text that indicates if the slice being processed is the first one or not. If 'first' user is prompted to click
        a point outside the contours close to the region than needs to be closed, else uses the click given in the first slice
        of that group of slices being closed which enters as input
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    myIm_closed : numpy array
        Closed image
    clicks : list of tuple with coordinates of last (clicks)
        Last clicked coordinates.

    """

    if option == 'first':
        plot_show = True
    else:
        plot_show = False
    
    contours_or = getContExpCont_plt (myIm_closed, slc, chStr, 250, 7, level, plot_show)
    contours_or = sorted(contours_or, key = len, reverse=True)
    ind_contours = []
    for i, cont in enumerate(contours_or): 
        props = maskContour(myIm_closed, cont)
        if props[2] > 15000:
            ind_contours.append(i)

    # print('ind_contours:', ind_contours)
    im_height, im_width = myIm_closed.shape
    black_array = np.uint16(np.zeros((150,im_width), dtype=int))

    myIm_closed = np.vstack((black_array,myIm_closed,black_array))
    myIm_closed = np.uint16(myIm_closed)

    contours_new = getContExpCont_plt (myIm_closed, slc, chStr, 250, 7, level, plot_show)
    contours_new = sorted(contours_new, key = len, reverse=True)
    
    # print('New contours: ',type(contours_new), len(contours_new))#, contours_new.shape)
    contours = []
    for index in ind_contours:
        contours.append(contours_new[index])
    # print('F contours: ',type(contours), len(contours))
    # print('Get contours: ',type(contours), len(contours))
    xy_contours = xy_allContours(contours)
    # print(xy_contours.shape, type(xy_contours))

    clicks = []
    #Get click to use convex hull to close
    if option == 'first':
        print("\n- Closing Inflow/Outflow tract - slice ", str(slc))
        # Get click for point to create convex hull from
        while len(clicks) == 0:
            clicks = getClicks(clicks_in, myIm_closed, scale=0.6, text='CONVEX HULL')
        # print('clicks first: ',str(clicks), type(clicks))
        clicks = clicks[-1]
        
    elif option =='autom':
        clicks = clicks_in
        # print('clicks autom: ',str(clicks), type(clicks))
            
    clicks_correct = False
    while not clicks_correct: 
        # Last point is considered the seed
        if len(clicks) > 0:
            # print ('clicks', clicks)
            # seed = clicks[0]
            # print ('seed', seed)
            y0, x0 = clicks
            point2add = np.array([[y0],[x0]])
            # print('point:', point2add.shape)
            xy_contours_new = np.concatenate((xy_contours, np.transpose(point2add)))
            qg_num = 'QG'+str(len(xy_contours_new)-1)
            # print(qg_num)
            hull = ConvexHull(points=xy_contours_new,qhull_options=qg_num)
            merge = hull.simplices[hull.good]
            closing_pt1, closing_pt2 = selectHull (merge, xy_contours_new)
            # print('closing pts: ', closing_pt1, closing_pt2, type(closing_pt1), type(closing_pt2))
            if type(closing_pt1) == tuple and type(closing_pt2) == tuple: 
                clicks_correct = True
                # print('clicks correct')
                break
            else: 
                alert('error',2)
                print('- ALERT!: Make sure you are clicking on a point outside the convex hull of the contours!')
                clicks = getClicks(clicks_in, myIm_closed, scale=0.6, text='CONVEX HULL')
                clicks = clicks[-1]

    if plot_show:
        fig, ax = plt.subplots(1,3, figsize=(10,3))
        ax[0].imshow(myIm_closed, cmap=plt.cm.gray)
        ax[0].set(xlabel= "y", ylabel = "x")
        ax[0].set_title("Original")
        for visible_facet in hull.simplices[hull.good]:
            ax[0].plot(hull.points[visible_facet, 1],
                    hull.points[visible_facet, 0],
                    color='red', lw=2)

        ax[1].imshow(myIm_closed, cmap=plt.cm.gray)
        ax[1].plot(closing_pt1[1], closing_pt1[0],'ro')
        ax[1].plot(closing_pt2[1], closing_pt2[0],'go')
        ax[1].set_title("Points to close")

    rr, cc, val = line_aa(int(closing_pt1[0]), int(closing_pt1[1]),
                          int(closing_pt2[0]), int(closing_pt2[1]))
    myIm_closed[rr, cc] = val * 50000

    contours_closed = measure.find_contours(myIm_closed, level, 'high', 'high')

    if plot_show:
        ax[2].imshow(myIm_closed, cmap=plt.cm.gray)
        for n, contour in enumerate(contours_closed):
            if len(contour)>250:
                ax[2].plot(contour[:, 1], contour[:, 0], linewidth=1)
        ax[2].set_title("Closed")
        plt.show()

    _ = getContExpCont_plt (myIm_closed, slc, chStr, 250, 10, level, plot_show)

    myIm_closed = myIm_closed[150:150+im_height]

    return myIm_closed, clicks

#%% func - xy_allContours
def xy_allContours (contours):
    """
    Function to create a unique array including all the points that make up a list of contours

    Parameters
    ----------
    contours : list of arrays
        list of numpy array with coordinates making up the contours of the slice being processed.

    Returns
    -------
    coordsXY : numpy array
        unique numpy array with all XY coordinates that make up all contours.

    """

    coordsXY = []

    if len(contours) == 1:
        # print('len(contours) = 1')
        coordsXY = np.array(contours)[0]
    else:
        for num, cont in enumerate(contours):
            coords2add = cont
            if len(coordsXY) == 0:
                coordsXY = coords2add
            else:
                coordsXY = np.concatenate((coordsXY,coords2add))

    # print('coordsXY:', type(coordsXY), coordsXY.shape)

    return coordsXY

#%% func - selectHull
def selectHull (merge, xy_contours):
    """
    Function that selects the longest good hull

    Parameters
    ----------
    merge :
        hull.simplices[hull.good].
    xy_contours : list of arrays
        list of numpy array with coordinates making up the contours of the slice being processed.

    Returns
    -------
    closing_pt1 : tuple with coordinates
        XY coordinates of 1st point to close the contour.
    closing_pt2 : tuple with coordinates
        XY coordinates of 2nd point to close the contour.

    """

    eu_dist = []
    pts1 = []
    pts2 = []
    if len(merge) > 0:
        for num, pair in enumerate(merge):
            #print(num,pair)
            x1 = xy_contours[pair[0]][0]
            y1 = xy_contours[pair[0]][1]
            xy1 = (x1,y1)
    
            x2 = xy_contours[pair[1]][0]
            y2 = xy_contours[pair[1]][1]
            xy2 = (x2,y2)
    
            eu_dist.append(distance.euclidean(xy1,xy2))
            pts1.append(xy1)
            pts2.append(xy2)
    
            index4max = np.where(eu_dist == np.amax(eu_dist))
            closing_pt1 = pts1[index4max[0][0]]
            closing_pt2 = pts2[index4max[0][0]]
    else: 
        closing_pt1 = []
        closing_pt2 = []
    

    return closing_pt1, closing_pt2

#%% func - closeInfOutSlc
def closeInfOutSlc(myIm_closed, slc, chStr, kw, kh, clicks_in, option, level):
    """
    Function that closes the inflow/outflow tract of the input image using the cropNcloseCont function

    Parameters
    ----------
    myIm_closed : numpy array
        Image to close using cropNcloseCont
    slc : int
        Slice number being processed.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    kw : int
        width of the rectangular region to crop around the clicked point.
    kh : int
        length of the rectangular region to crop around the clicked point.
    clicks_in : list of tuples with coordinates (clicks)
        Initialised list of clicks coordinates to close the inflow/outflow with cropNcloseCont
    option : str
        Text that indicates if the slice being processed is the first one or not. If 'first' user is prompted to click
        a point in between the contours to close, else uses the clicks given in the first slice
        of that group of slices being closed which enters as input
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    myIm_closed : numpy array
        Closed image
    final_click : list of tuple with coordinates of clicks
        Clicked coordinates.

    """

    #Get clicks of positions to close contours
    if option == 'first':
        plot_show = True
    else:
        plot_show = False

    while True:
        if option == 'first':
            print("\n- Closing Inflow/Outflow tract - slice ", str(slc))
            clicks = getClicks([], myIm_closed, scale=1, text='CLOSING INF/OUTFLOW')
        else:
            clicks = clicks_in
        #Close contours and get Image
        myIm_closed = cropNcloseCont (clicks, myIm_closed, kw, kh, level)
        #Plot image with closed contours
        _ = getContExpCont_plt (myIm_closed, slc, chStr, 250, 10, level, plot_show)
        break

    final_click = clicks

    return myIm_closed, final_click

#%% - DICTIONARIES
#%% func - dictCreation
def dictCreation (filename, file, minLenContour, chStr, slcCont, numCont,
                      xy_Scaling_um, z_Scaling_um, heartLayer, update, curr_v):
    """
    Creation of the heartLayer dictionary with basic info about the file being processed

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    file : str
        Loaded file name (.tif/.npy)
    minLenContour : int
        Minimum number of points making up contoura and used to filter contours found in each slice
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    slcCont : list of int
        List of slice numbers where contour groups start.
    numCont : list of int
        List of contour numbers in each of the contour groups.
    xy_Scaling_um : float
        x,y scaling values of the images taken. This information is taken from the metadata of the original file.
    z_Scaling_um : float
        z scaling values of the images taken. This information is taken from the metadata of the original file.
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates
    update : boolean
        True if the dictionary has already been created and is going to be updated, false otherwise.
    curr_v : int
        Number indicating the version of the dictionary.

    Returns
    -------
    heartLayer : dictionary
        Updated dictionary with information about the heart layer being processed including selected contours number and coordinates
    vNum: int
        Number indicating the version of the dictionary.
    update : boolean
        True, the dictionary has already been created and if this function is called again, it will only be updated.

    """

    #Create heartLayer dict
    if update == False:
        #Create heartLayer dict
        heartLayer = dict()
        #Create Info dict inside heartLayer
        infoDict = heartLayer["info"]= dict()

        #Save initial info
        infoDict["AAA_LS/F_ref"] = filename
        infoDict["AAA_chAnalysed"] = chStr
        infoDict["xy_Scaling_um"] = xy_Scaling_um
        infoDict["z_Scaling_um"] = z_Scaling_um
        infoDict["numCont_group"] = str(numCont)
        infoDict["slcCont_group"] = str(slcCont)
        vNum = "0"
        msg = "- heartLayer dictionary created! (v_"+vNum+")"

    else:
        vNum = str(curr_v+1)
        msg = "- heartLayer dictionary updated! (v_"+vNum+")"

    #Define titles in case a new version of the dictionary is to be saved
    proc_file_tt = "Proc_File_v"+vNum
    minNumPts_tt = "minNunPts_v"+vNum
    dateNtime_tt = "Date&Time_v"+vNum

    if update == False:
        #Save info
        infoDict[proc_file_tt] = file
        infoDict[minNumPts_tt] = minLenContour
        infoDict[dateNtime_tt] = str(datetime.now())

    else:
        #Input data into dictionary
        heartLayer["info"][proc_file_tt] = file
        heartLayer["info"][minNumPts_tt] = minLenContour
        heartLayer["info"][dateNtime_tt] = str(datetime.now())

    alert('wohoo',1)
    print(msg)
    update = True

    return heartLayer, int(vNum), update

#%% func - dictFill
def dictFill (tuple_slc, numContours, stack, chStr, heartLayer, minLenContour, level):
    """
    Function to fill the dictionary with info about the contours

    Parameters
    ----------
    tuple_slc : list of tuples of int
        List with tuples including all slices.
    numContours : list of int
        List of contour numbers in each of the contour groups.
    stack : numpy array
        Stack being processed.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates
    minLenContour : int
        Minimum number of points making up contoura and used to filter contours found in each slice
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    heartLayer : dictionary
        Updated dictionary with information about the heart layer being processed including selected contours number and coordinates
    tuple_slc : list of tuples of int
        List with tuples including all slices or a selected group of slices.
    numContours : list of int
        Updated List of contour numbers in each of the contour groups in which the heartLayer dict hasn't been filled yet.
    exit_txt : boolean
        True if user wants to exit code, else False.

    """

    exit_txt = False

    tic = perf_counter()
    chNum = int(chStr[-1])
    #print("Analysis of channel: ", chNum)

    for num, pair in enumerate(tuple_slc):
        #print(num, pair[0], pair[1])
        numContours_pair = numContours[num]
        #print("The number of contours expected in slcs", pair, "is: ", numContours_pair)

        #If there are no contours in the slc
        if numContours_pair == 0:
            for slc in range(pair[0],pair[1],1):
                print("- Save empty dict! - slc:", slc)
                myIm = stack[slc][:][:]
                # Creation of a dictionary for each slice within the heartLayer dictionary
                sliceNum = slcNum_def (slc)
                slcDict = heartLayer[sliceNum] = dict()

                #Save empty data in dictionary
                slcDict["FirstSlc"] = "Empty"
                slcDict["image_o"] = myIm
                slcDict["intContours"] = np.nan
                slcDict["extContours"] = np.nan
                slcDict["allContours"] = np.nan
                slcDict["index_intContours"] = ""
                slcDict["index_extContours"] = ""
                slcDict["index_allContours"] = ""
                slcDict["imIntFilledCont"] = np.zeros([myIm.shape[0],myIm.shape[1]])
                slcDict["imExtFilledCont"] = np.zeros([myIm.shape[0],myIm.shape[1]])
                slcDict["imAllFilledCont"] = np.zeros([myIm.shape[0],myIm.shape[1]])
                slcDict["pixIntFilledXY"] = np.nan
                slcDict["pixExtFilledXY"] = np.nan
                slcDict["pixAllFilledXY"] = np.nan

        else:
            # print('\n')
            bar = Bar('- Selecting contours', max=pair[1]-pair[0], suffix = suffix, check_tty=False, hide_cursor=False)
            for slc in range(pair[0],pair[1],1):
                myIm = stack[slc][:][:]
                #print("slc:", slc)

                # Creation of a dictionary for each slice within heartLayer dictionary
                sliceNum = slcNum_def (slc)
                slcDict = heartLayer[sliceNum] = dict()

                #Get the contours of slice
                # --getContExpCont_slc (myIm, sliceNum, chNum, minLenContour, plotshow, imageEvery)
                contours, numCont = getContExpCont(myIm, minLenContour, level)
                #Sort the contours by length (bigger to smaller)
                contours = sorted(contours, key = len, reverse=True)

                if slc == pair[0]:
                    slcDict["FirstSlc"] = "Yes"
                    # print("\n- Pair ", pair[0], ",", pair[1]-1," / slc", slc)
                    print("\n >> Pair ("+ str(pair[0])+","+str(pair[1]-1)+")/slc"+str(slc)+" - Number of expected contours: "+str(numContours_pair))
                    #Run function to select contours and get its metrics
                    # --selectContours_slc (ch, slc, myIm, cont_sort, num_contours, intNext, ask2save)
                    selectCont, txtSelect, propSelect, exit_txt = selectContours (myIm, slc, chNum, contours, numCont, numContours_pair, True)
                    if exit_txt:
                        print('\r- Exiting filling pair!')
                        break
                    else:
                        int_contNprops, ext_contNprops = shape_contNprops(selectCont, propSelect, False)
                        index_intContours = txtSelect[0]
                        index_extContours = txtSelect[1]
                        list_index = index_intContours+ " " +index_extContours

                else:
                    slcDict["FirstSlc"] = "No"
                    #print("Automatically selecting the contours of the next slices within the group - slc", slc, "/",pair[1]-1)
                    #print("Get properties...")
                    props_all = getProperties (myIm, contours)
                    #print("Select contours automatically...")
                    list_index = automSelectContours (propSelect, props_all, numContours_pair)
                    #print("list_index: ", list_index, " - slc:", slc)
                    contNprops = extractContours(list_index, contours, props_all, False, False)
                    #print("Group contours as int/ext -- slc:", slc)
                    index_intContours, index_extContours = classifyCont (contNprops[1])
                    #print("Get corresponding contours")
                    int_contNprops = extractContours(index_intContours, contNprops[0], contNprops[1], False, False)
                    ext_contNprops = extractContours(index_extContours, contNprops[0], contNprops[1], False, False)


                #Get data from selected contours for slice slc
                # All contours
                allContours = [y for x in [int_contNprops[0], ext_contNprops[0]] for y in x]

                #print("Get the images with filled contour")
                imIntFilledCont, pixIntFilledXY = fillContours(myIm, int_contNprops[0])
                imExtFilledCont, pixExtFilledXY = fillContours(myIm, ext_contNprops[0])
                imAllFilledCont, pixAllFilledXY = fillContours(myIm, allContours)
                plotFilledCont(myIm, allContours, imIntFilledCont, imExtFilledCont, imAllFilledCont, True, slc)
                propSelect.clear()
                propSelect.append(int_contNprops[1])
                propSelect.append(ext_contNprops[1])

                #Save all data in dictionary
                slcDict["image_o"] = myIm
                slcDict["intContours"] = int_contNprops[0]
                slcDict["extContours"] = ext_contNprops[0]
                slcDict["allContours"] = allContours
                slcDict["index_intContours"] = " ".join(str(item) for item in index_intContours)
                slcDict["index_extContours"] = " ".join(str(item) for item in index_extContours)
                slcDict["index_allContours"] = " ".join(str(item) for item in list_index)
                slcDict["imIntFilledCont"] = imIntFilledCont
                slcDict["imExtFilledCont"] = imExtFilledCont
                slcDict["imAllFilledCont"] = imAllFilledCont
                slcDict["pixIntFilledXY"] = pixIntFilledXY
                slcDict["pixExtFilledXY"] = pixExtFilledXY
                slcDict["pixAllFilledXY"] = pixAllFilledXY
                bar.next()

            bar.finish()

        if exit_txt:
            print('\r- Exiting filling tuples!')
            tuple_slc = tuple_slc[num:]
            numContours = numContours[num:]
            break

    # if exit_txt:
    #     print('Exiting filling heartLayer!')
    #     break

    if not exit_txt:
        toc = perf_counter()
        time = toc-tic
        tuple_slc = []
        alert('wohoo',1)
        print("- All Done - heartLayer dictionary was filled!")
        print("- Time taken to fill dictionary = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")

    return heartLayer, tuple_slc, numContours, exit_txt

#%% func - askModifyDict
def askModifyDict(filename, stack, stack_o, stack_m, heartLayer, channel, processDict, directories, level):
    """
    Function that asks the user if the dictionary created wants to be modified

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    stack : numpy array
        Stack being processed.
    stack_o : numpy array
        Copy of closed stack.
    stack_m : numpy array
        Copy of original masked stack to process.
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates
    channel : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    processDict : dictionary
        Dictionary with information about the processed carried out when closing contours
    directories : list of paths
        List of paths to the folders where the files are saved.
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours
    
        
    Returns
    -------
    heartLayer : dictionary
        Updated dictionary with information about the heart layer being processed including selected contours number and coordinates
    stack : numpy array
        Stack being processed.
    processDict : dictionary
        Updated dictionary with information about the processed carried out when closing contours

    """

    #% Check slices
    plotSelectedContours(imageEvery = 1, stack = stack, heartLayer = heartLayer, plot_all = False)
    
    #% Correct dictionary if needed
    q_modifyDict = ask4input("Do you want to modify the created dictionary (change any of the selected contours)? [0]:no/[1]:yes: ",bool)
    if q_modifyDict:
        while q_modifyDict:
            q_firstClose = ask4input("Do you need to close any contours before you modify the dictionary? [0]:no/[1]:yes: ",bool)
            if q_firstClose:
                stack, processDict, _ = main_manuallyCloseContours(filename, channel, directories, stack, 
                                                                         stack_o, stack_m, processDict, 7, level, True)
            heartLayer = modifyDict(stack = stack, chStr = channel, heartLayer = heartLayer,
                                                  minLenContour = 250, level = level)
            q_modifyDict = ask4input("Do you want to modify any other part of the dictionary? [0]:no/[1]:yes: ",bool)
            if not q_modifyDict:
                q_plotAll = ask4input('Do you want to plot the final selected contour masks for '+channel+'? [0]:no/[1]:yes: ',bool)
                if q_plotAll:
                    imageEvery = ask4input('Plot selected contours every X number of slices. X= ',int)
                    print('\n- Plotting selected contours ... (Plots tab)')
                    plotSelectedContours(imageEvery = imageEvery, stack = stack, heartLayer = heartLayer, plot_all = True)
                break
            else:
                plotSelectedContours(imageEvery = 1, stack = stack, heartLayer = heartLayer)

    return heartLayer, stack, processDict

#%% func - modifyDict
def modifyDict(stack, chStr, heartLayer, minLenContour, level):
    """
    Function to modify the contour info saved in the dictionary

    Parameters
    ----------
    stack : numpy array
        Stack being processed.
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates
    minLenContour : int
        Minimum number of points making up contoura and used to filter contours found in each slice
    level : float 
        Value along which to find contours in the array. Parameter for skimage.measure.find_contours

    Returns
    -------
    heartLayer : dictionary
        Updated dictionary with information about the heart layer being processed including selected contours number and coordinates

    """

    alert('frog',1); print("\n- Running script to modify dict...!")

    exit_txt = False
    not_new = True
    while not_new:
        q_newTuples = ask4input('Do you want to select the contours for a new group of slices OR for individual slices? \n\t[0]:no, I am done selecting contours\n\t[1]:yes, group of slices\n\t[2]:yes, individual slices: ', int)
        if q_newTuples == 1:
            tuple_slc_up, numCont_upf = new_tuples_from_originals(heartLayer)
            not_new = False
        elif q_newTuples == 2:
            tuple_slc_up, numCont_upf = new_indiv_tuples(heartLayer)
            not_new = False
        elif q_newTuples == 0:
            tuple_slc_up = []
            numCont_upf = []
            not_new = False
        else:
            print('ERROR: -'+q_newTuples+'- The number entered needs to be a [0],[1] or [2]!')

    print(tuple_slc_up)
    print(numCont_upf)

    while len(tuple_slc_up) != 0:
        heartLayer, tuple_slc_up, numCont_upf, exit_txt = dictFill (tuple_slc = tuple_slc_up, numContours = numCont_upf,
                                            stack = stack, chStr = chStr, heartLayer = heartLayer, minLenContour = 250, level = level)
        if exit_txt:
            print('EXIT - Modifying dictionary!!');
            break

    if not exit_txt:
        alert('wohoo',1)
        print("- All Done - heartLayer dictionary was modified!")
    else:
        print('-EXITING!')

    return heartLayer

#%% func - smallDict2Save
def smallDict2Save (heartLayer):
    """
    Function to save a small version of the dictionary created to recover contours selected in future

    Parameters
    ----------
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates

    Returns
    -------
    heartLayer2S : dictionary
        heartLayer2S dictionary filled with an small version.

    """

    #Define the keys to save of each slice
    keys2save = ["FirstSlc", "intContours", "extContours", "allContours",
                     "index_intContours", "index_extContours", "index_allContours"]

    heartLayer2S = dict()
    for slc, keySlc in enumerate(heartLayer.keys()):
        infoDict2S = heartLayer2S["info"]= dict()
        for data, keyInfo in enumerate(heartLayer["info"].keys()):
            infoDict2S[keyInfo] = heartLayer["info"][keyInfo]

        slcDict2S = heartLayer2S[keySlc] = dict()
        if keySlc[0:3] == "slc":
            for num, key in enumerate(keys2save):
                slcDict2S[key] = heartLayer[keySlc][key]

    alert('wohoo',1)
    print("- heartLayer dictionary TO SAVE has been created!!!")

    return heartLayer2S

#%% - CONTOURS SELECTION (INT AND EXT)
#%% func - selectContours
def selectContours (myIm, slcNum, chNum, cont_sort, num_contours, numContours_pair, ask2save):
    """
    Function to select the contours to add (all/int&ext)

    Parameters
    ----------
    myIm : numpy array
        Imaged being processed
    slcNum : int
        Slice number being processed.
    chNum : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    cont_sort : list of arrays
        list of arrays with contours' coordinates sorted from longer to shorter.
    num_contours : int
        Number of contours found in slice.
    numContours_pair : int
        Number of contours found in the contour group to which the processed slice belongs.
    ask2save : boolean
        True if you want the user to double check the selected contours before continuing, else False.

    Returns
    -------
    None.

    """
    #Print the channel and slice being analysed
    #print("Channel "+str(chNum)+" / Slice "+str(slcNum))

    #Define the figure properties (columns, rows, image size)
    cols = 3
    rows = num_contours // cols
    if num_contours%cols != 0:
        rows = rows + 1
    imSize = 3
    colorImSize = 4
    fig11 = plt.figure(figsize=(cols*imSize+colorImSize, rows*imSize), constrained_layout=True)

    # gridspec inside gridspec
    outer_grid = fig11.add_gridspec(1,2)
    # Grid where color image will be placed
    color_grid = outer_grid[0].subgridspec(1,1, wspace=0, hspace=0)
    ax = fig11.add_subplot(color_grid[0])
    ax.imshow(myIm, cmap=plt.cm.gray)

    # Go through all the contours
    for num, contour in enumerate(cont_sort):
        ax.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color=colors[num])
        txt = "Cont"+str(num)
        ax.text(0.95,(0.97-0.035*(num+1)), txt,
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=ax.transAxes,
                    color=colors[num], fontsize=10, weight = 'semibold')
#    ax.set_xticks([])
#    ax.set_yticks([])
    ax.set(xlabel = "Channel "+str(chNum)+" / Slice "+str(slcNum) + " / Contours expected: "+ str(numContours_pair))
    ax.set_title("Channel "+str(chNum)+" / Slice "+str(slcNum) + " \nContours expected: "+ str(numContours_pair), fontsize = 14, fontweight='bold')

    # Grid where subplots of each contour will be placed
    all_grid = outer_grid[1].subgridspec(rows, cols,  wspace=0.3, hspace=0.3)

    # Array to save the metrics of all the contours
    props_all = [None]*len(cont_sort)
    #Iterate through sorted contours
    for index, contList in enumerate(cont_sort):

        #-->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]
        props = maskContour(myIm, contList)
        props_all[index] = props

        max_int = format(props[2], '.2f')
        mean_int = format(props[3], '.2f')
        # sol = format(props[6], '.2f')
        per = format(props[5], '.0f')
        area = format(props[0], '.0f')
        #cx_area = format(props[2], '.0f')
        # lgth = format(props[4],'.0f')

        ax = fig11.add_subplot(all_grid[index])
        ax.imshow(myIm, cmap=plt.cm.gray)
        ax.plot(contList[:, 1], contList[:, 0], linewidth=1.5, color = colors[index])
        ax.set_title("Contour "+str(index), fontsize=10, weight = 'semibold', color = colors[index])
        ax.set(xlabel = "Area:"+str(area)+" / Per:"+str(per), ylabel = "MaxInt:"+str(max_int)+" / MeanInt:"+str(mean_int))
        ax.set_xticks([])
        ax.set_yticks([])

    plt.show()

    #Create list of contours to export
    cont2Exp = [None]*2
    #Save the number of the contours being exported
    txt2Exp = [None]*2
    #Save the properties of the contours being exported
    props2Exp = [None]*2

    exit_txt = False

    if ask2save == True:
       alert('bubble', 1)
       sureBoth = False

       while not sureBoth:
           # Internal
           #print('- \n[Ch'+str(chNum)+'/Slc'+str(slcNum)+']')
           print('- Enter the contour numbers for slice '+str(slcNum)+':')
           selecContInt = str(input("\t> INTERNAL contours (separated by SPACES)/'esc': ")).lower()
           if selecContInt == 'esc':
               exit_txt = True
               break
           # External
           selecContExt = str(input("\t> EXTERNAL contours (separated by SPACES)/'esc': ")).lower()
           if selecContExt == 'esc':
               exit_txt = True
               break
           q_sureBoth = str(input("> Are you sure?  - INTERNAL: "+ selecContInt + " / EXTERNAL: "+selecContExt+" - >>> [y]:yes/[n]:no/[esc]:exit: ")).lower()
           if q_sureBoth == "y":
               sureBoth = True
               print('\n')
               break
           elif q_sureBoth == 'esc':
               exit_txt = True
               break

       if exit_txt:
           print("- Exiting!")
           #break
       else:
           txt2Exp[0] = selecContInt
           intCont2add = selecContInt.split()
           intCont2add = [int(j) for j in intCont2add]

           #Empty lists to save internal contours and its metrics
           intContours = [None] * len(intCont2add)
           intProps = [None] * len(intCont2add)
           for index, cont2addInt in enumerate(intCont2add):
               intContours[index] = cont_sort[cont2addInt]
               intProps[index] = props_all[cont2addInt]

           cont2Exp[0] = intContours
           props2Exp[0] = intProps

           txt2Exp[1] = selecContExt
           extCont2add = selecContExt.split()
           extCont2add = [int(k) for k in extCont2add]

           #Empty lists to save external contours and its metrics
           extContours = [None] * len(extCont2add)
           extProps = [None] * len(extCont2add)
           for index, cont2addExt in enumerate(extCont2add):
               extContours[index] = cont_sort[cont2addExt]
               extProps[index] = props_all[cont2addExt]

           cont2Exp[1] = extContours
           props2Exp[1] = extProps

    # if ask2save == True:
    #     alert('bubble', 1)
    #     sureInt = False
    #     sureExt = False
    #     # Internal
    #     while True:
    #         selecContInt = str(input("> Enter INTERNAL contours to add [Ch"+str(chNum)+"/Slc"+str(slcNum)+"] (separated by SPACES)/'esc': ")).lower()
    #         if selecContInt == 'esc':
    #             exit_txt = True
    #             break
    #         sureInt = str(input("> Are you sure these are the INTERNAL contour numbers: -"+ selecContInt + "?- [y]:yes/[n]:no/[esc]:exit: ")).lower()
    #         if sureInt == "y":
    #             break
    #         elif sureInt == 'esc':
    #             exit_txt = True
    #             break

    #     if exit_txt:
    #         print("- Exiting!")
    #         #break
    #     else:
    #         txt2Exp[0] = selecContInt
    #         intCont2add = selecContInt.split()
    #         intCont2add = [int(j) for j in intCont2add]

    #         #Empty lists to save internal contours and its metrics
    #         intContours = [None] * len(intCont2add)
    #         intProps = [None] * len(intCont2add)
    #         for index, cont2addInt in enumerate(intCont2add):
    #             intContours[index] = cont_sort[cont2addInt]
    #             intProps[index] = props_all[cont2addInt]

    #         cont2Exp[0] = intContours
    #         props2Exp[0] = intProps


    #     #alert('bubble', 1)
    #     # External
    #     while True:
    #         selecContExt = str(input("> Enter EXTERNAL contours to add [Ch"+str(chNum)+"/Slc"+str(slcNum)+"] (separated by a SPACES)/'esc': ")).lower()
    #         if selecContExt == 'esc':
    #             exit_txt = True
    #             break
    #         sureExt = str(input("> Are you sure these are the EXTERNAL contour numbers: -"+ selecContInt + "?- [y]:yes/[n]:no/[esc]:exit: ")).lower()
    #         if sureExt == "y":
    #             print('\n')
    #             break
    #         elif sureExt == 'esc':
    #             exit_txt = True
    #             break

    #     if exit_txt:
    #         print("- Exiting!")
    #         #break
    #     else:
    #         txt2Exp[1] = selecContExt
    #         extCont2add = selecContExt.split()
    #         extCont2add = [int(k) for k in extCont2add]

    #         #Empty lists to save external contours and its metrics
    #         extContours = [None] * len(extCont2add)
    #         extProps = [None] * len(extCont2add)
    #         for index, cont2addExt in enumerate(extCont2add):
    #             extContours[index] = cont_sort[cont2addExt]
    #             extProps[index] = props_all[cont2addExt]

    #         cont2Exp[1] = extContours
    #         props2Exp[1] = extProps

    return cont2Exp, txt2Exp, props2Exp, exit_txt

#%% func - shape_contNprops
def shape_contNprops (selectCont, propSelect, empty):
    """
    Function to save selected contours in same way as extractContours function

    Parameters
    ----------
    selectCont : list of arrays
        list of numpy array with coordinates making up the selected contours
    propSelect : list of numpy array
        List of arrays with values of properties associated to each of the selected contours found in slice
        -->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]
    empty : boolean
        rue if you want to return empty structures, else False.

    Returns
    -------
    int_contNprops : list of numpy array
        List of arrays with values of properties associated to the INTERNAL contours found in slice
    ext_contNprops : list of numpy array
        List of arrays with values of properties associated to the EXTERNAL contours found in slice

    """


    int_contNprops = []
    ext_contNprops = []

    if empty:
        empty_list = [None]
        for num in range(2):
            int_contNprops.append(empty_list)
            ext_contNprops.append(empty_list)
    else:

        int_contNprops.append(selectCont[0])
        int_contNprops.append(propSelect[0])

        ext_contNprops.append(selectCont[1])
        ext_contNprops.append(propSelect[1])

    return int_contNprops, ext_contNprops

#%% func - getProperties
def getProperties (myIm, contours):
    """
    Function to get properties of all the identified contours

    Parameters
    ----------
    myIm : numpy array
        Imaged being processed
    contours : list of arrays
        list of numpy array with coordinates making up the contours of the slice being processed.

    Returns
    -------
    propsAllCont : list of numpy array
        List of arrays with values of properties associated to each of the identified contour in slice
        -->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]

    """
    # Array to save the metrics of all the contours
    propsAllCont = [None]*len(contours)
    #Iterate through sorted contours
    for index, contList in enumerate(contours):

        #[0. area, 1. centroid, 2. cx_area, 3. ecc, 4. lgth, 5. per, 6. sol, 7.bbox])
        props = maskContour(myIm, contList)
        propsAllCont[index] = props

    return propsAllCont

#%% func - automSelectContours
def automSelectContours (propContSelected, propsAllCont, numContours_pair):
    """
    Function to automatically select contours by area, centroid position, mean intensity and perimeter

    Parameters
    ----------
    propContSelected : list of numpy array
        List of arrays with values of properties associated to the selected contours found in previous slice
        -->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]
    propsAllCont : list of numpy array
        List of arrays with values of properties associated to all the contours found in current slice
        -->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]
    numContours_pair : int
        Number of contours found in the contour group to which the processed slice belongs.

    Returns
    -------
    list_index : list of int
        List of int with the indexes of the contours automatically selected.

    """
    lenPropsAllCont = len(propsAllCont)
    #print('lenPropsAllCont:',lenPropsAllCont)
    #0.area, 1.centroid, 2.meanInt, 3.perimeter
    scale_imp = [0.20,0.70,0,0.10]
    numC = 0

    index = []
    #Iterate through the selected contours of the previous slice and find a match for each
    for j in range(len(propContSelected)):
        #Enter if there were actually contours selected
        if len(propContSelected[j]) != 0 and numC < numContours_pair:
            for numContSel, propContSel in enumerate(propContSelected[j]):
                #print('numContours added: ', numC)
                # Get properties of each of the selected processDict
                area_sel = propContSel[0]
                cent_sel = propContSel[1]
                meanInt_sel = propContSel[3]
                per_sel = propContSel[3]

                # Create empty array to save distances
                bigNum = 10**20
                dif_area = np.ones(lenPropsAllCont)*bigNum
                dif_cent = np.ones(lenPropsAllCont)*bigNum
                dif_meanInt = np.ones(lenPropsAllCont)*bigNum
                dif_per = np.ones(lenPropsAllCont)*bigNum

                max_area = 0; max_cent = 0; #max_maxInt = 0;
                max_meanInt  = 0; max_per = 0

                tot = np.zeros(lenPropsAllCont)
                for numCont, propCont in enumerate(propsAllCont):
                    if not numCont in index:
                        # Get properties of each of the contours in next slide
                        area_cont = propCont[0]
                        cent_cont = propCont[1]
                        # maxInt_cont = propCont[2]
                        meanInt_cont = propCont[3]
                        per_cont = propCont[5]

                        # Save difference in properties
                        dif_area[numCont] = abs(area_sel-area_cont)
                        if dif_area[numCont] > max_area:
                            max_area = dif_area[numCont]
                        dif_cent[numCont] = abs(distance.euclidean(cent_sel, cent_cont))
                        if dif_cent[numCont] > max_cent:
                            max_cent = dif_cent[numCont]
                        dif_meanInt[numCont] = abs(meanInt_sel-meanInt_cont)
                        if dif_meanInt[numCont] > max_meanInt:
                            max_meanInt = dif_meanInt[numCont]
                        dif_per[numCont] = abs(per_sel-per_cont)
                        if dif_per[numCont] > max_per:
                            max_per = dif_per[numCont]

                try: 
                    dif_area = dif_area/max_area
                except: 
                    dif_area = 0
                try: 
                    dif_cent = dif_cent/max_cent
                except: 
                    dif_cent = 0
                try:
                    dif_meanInt = dif_meanInt/max_meanInt
                except: 
                    dif_meanInt = 0
                try: 
                    dif_per = dif_per/max_per
                except: 
                    dif_per = 0

                # print(dif_area)
                # print(dif_cent)
                # print(dif_maxInt)
                # print(dif_meanInt)
                # print(dif_per)

                for num in range(lenPropsAllCont):
                    if not num in index:
                        #0.area, 1.centroid, 2.meanInt, 3.perimeter
                        tot[num] = dif_area[num]*scale_imp[0]+dif_cent[num]*scale_imp[1]+dif_meanInt[num]*scale_imp[2]+dif_per[num]*scale_imp[3]
                    else:
                        tot[num] = 10**30
                # print('tot:', tot)
                index.append(np.where(tot == min(tot))[0][0])

    # print('index:', index)
    list_index = index[0:numContours_pair]
    # print("\t- Selected contours:", list_index)

    return list_index

# def automSelectContours (propContSelected, propsAllCont, numContours_pair):

#     """ automSelectContours_slc
#     Function to select all/int&ext contours per slice

#     Parameters:
#         - propContSelected: list - list of numpy arrays with properties of each of the selected contours
#         - propsAllCont: list - list of numpy arrays with properties associated to each contour
#                      [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]

#     Returns:
#         - list_index: list - list of int with the indexes of the contours automatically selected

#      Version: April 15, 2020

#     """

#     maxNumCont = len(propContSelected[0])+len(propContSelected[1])
#     list_index = [None]*numContours_pair#maxNumCont

#     #print('list_index:',list_index)
#     #print(len(list_index))
#     scale_imp = [0.20,0.70,0.05,0.05]

#     numC = 0
#     #Iterate through the selected contours of the previous slice and find a match for each
#     for j in range(len(propContSelected)):
#         #print('Number of contours to find match: ', len(propContSelected[j]))
#         #Enter if there were actually contours selected
#         if len(propContSelected[j]) != 0 and numC < maxNumCont:

#             for numContSel, propContSel in enumerate(propContSelected[j]):
#                 #print('numContours added: ', numC)
#                 # Get properties of each of the selected processDict
#                 area_sel = propContSel[0]
#                 cent_sel = propContSel[1]
#                 lgth_sel = propContSel[4]
#                 per_sel = propContSel[5]

#                 # Create empty array to save distances
#                 dif_area = np.ones(len(propsAllCont))*10**20
#                 dif_cent = np.ones(len(propsAllCont))*10**20
#                 dif_lgth = np.ones(len(propsAllCont))*10**20
#                 dif_per = np.ones(len(propsAllCont))*10**20

#                 for numCont, propCont in enumerate(propsAllCont):
#                     if not numCont in list_index:
#                         # Get properties of each of the contours in next slide
#                         area_cont = propCont[0]
#                         cent_cont = propCont[1]
#                         lgth_cont = propCont[4]
#                         per_cont = propCont[5]

#                         # Save difference in properties
#                         dif_area[numCont] = abs(area_sel-area_cont)
#                         dif_cent[numCont] = abs(distance.euclidean(cent_sel, cent_cont))
#                         dif_lgth[numCont] = abs(lgth_sel-lgth_cont)
#                         dif_per[numCont] = abs(per_sel-per_cont)

#                 dif_area = dif_area/max(dif_area)
#                 dif_cent = dif_cent/max(dif_cent)
#                 dif_lgth = dif_lgth/max(dif_lgth)
#                 dif_per = dif_per/max(dif_per)

#                 #print("dif_area:", dif_area)
#                 #print("dif_cent:", dif_cent)
#                 #print("dif_lgth:", dif_lgth)
#                 #print("dif_per:", dif_per)

#                 tot = np.zeros(len(propsAllCont))
#                 for num in range(len(propsAllCont)):
#                     if not num in list_index:
#                         tot[num] = dif_area[num]*scale_imp[0]+dif_cent[num]*scale_imp[1]+dif_lgth[num]*scale_imp[2]+dif_per[num]*scale_imp[3]
#                     else:
#                         tot[num] = 10**30
#                 #print("tot:", tot)

#                 index = np.where(tot == min(tot))[0][0]
#                 #print("index:", index)

#                 if numC <= maxNumCont-1:
#                     list_index[numC] = index
#                     print("\t- Selected contours:", list_index)
#                     numC += 1

#     return list_index

#%% func - extractContours
def extractContours (list_index, contours, propsAllCont, justContours, listArrayType):
    """
    Function to get contours and properties of list_index

    Parameters
    ----------
    list_index : list of int
        List of int with the indexes of the contours automatically selected.
    contours : list of arrays
        list of numpy array with coordinates making up the contours of the slice being processed.
    propsAllCont : list of numpy array
        List of arrays with values of properties associated to all the contours found in current slice
        -->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]
    justContours : boolean
        If True, then the function returns just contours. If False, the function returns contours and its properties
    listArrayType : boolean
        If True, the list_index comes as [(array([3], dtype=int64),),...]. If False, the the list_index comes as [1,0]

    Returns
    -------
    val2return : List
        list containing:
        - if jutsContours == True   -- val2return contains in [0] a list of the extracted contours
                                    -- val2return contains in [1] a list of the properties of the extracted contours
        - if jutsContours == False  -- val2return contains in [0] a list of the extracted contours
                                    -- val2return contains in [1] an empty list .

    """
    contExtracted = []
    propsExtracted = []

    if listArrayType:
        for num in range(len(list_index)):
            index = list_index[num][0][0]
            #print("index: ", index)
            contExtracted.append(contours[index])
            if not justContours:
                propsExtracted.append(propsAllCont[index])

    else:
        for index in list_index:
            contExtracted.append(contours[index])
            if not justContours:
                propsExtracted.append(propsAllCont[index])

    val2return = [contExtracted, propsExtracted]

    return val2return

#%% func - classifyCont
def classifyCont (propsCont2class):
    """
    Function to classify contours as external or internal

    Parameters
    ----------
    propsCont2class : list of numpy arrays
        List of arrays with values of properties associated to the selected contours that need to be classified as
        internal/external
        -->>#2 [0. area, 1. centroid, 2. max_int, 3. mean_int, 4. lgth, 5. per, 6. sol, 7. bbox]

    Returns
    -------
    int_contours : list of arrays
        resulting list of numpy array with coordinates making up the INTERNAL contours of the slice being processed.
    ext_contours : list of arrays
        resulting list of numpy array with coordinates making up the EXTERNAL contours of the slice being processed.

    """
    # Create an empty array to save all the contours of each slice individually
    bbox_region = []
    cent_region = []
    area_region = []

    ext_contours = []
    int_contours = []
    # Go through all the properties and extract them into arrays
    for num, props in enumerate(propsCont2class):
        row_c, col_c = props[1]
        cent_region.append([row_c,col_c])
        [min_row, min_col, max_row, max_col] = props[-1]
        bbox_reg = [min_row, max_row,min_col, max_col]
        bbox_region.append(bbox_reg)
        area_region.append(props[0])

    #print("bbox:",bbox_region)
    #print("cent:",cent_region)
    #print("area:",area_region)

    #Sort the area vector from biggest to smallest
    max_area = sorted(area_region, reverse=True)
    # Get the organised indices
    org_indexes = []
    for val in max_area:
        index_val = np.where(area_region == val)
        org_indexes.append(index_val[0][0])
    #print("organized_indexes: ",org_indexes)

    #Value to cut org_index array
    it = 0
    for index_big in org_indexes:
        if index_big not in int_contours:
            #print("INDEX_BIG:",index_big)
            # Get bounding box of the first region
            [min_row_b, max_row_b, min_col_b, max_col_b] = bbox_region[index_big]
            #Get the rest of the organised indexes to iterate through those
            org_indexes_sm = org_indexes[it+1:]
            #print(">> org_indexes_sm:", org_indexes_sm)
            for index_sm in org_indexes_sm:
                if index_big not in int_contours:
                    #print("INDEX_SM:",index_sm)
                    # Get bounding box of the small region
                    [min_row_s, max_row_s, min_col_s, max_col_s] = bbox_region[index_sm]

                    #Compare the bbox of the big index with the bbox of small index
                    if min_row_s > min_row_b  and max_row_s < max_row_b:
                        #print("     Row analysis: cont", index_big, "is outside cont", index_sm)
                        if min_col_s > min_col_b  and max_col_s < max_col_b:
                            #print("     Col analysis: cont", index_big, "is outside cont", index_sm)
                            int_contours.append(index_sm)
                            #print(" >>Internal:", index_sm, ", External:", index_big)

            #Get the centroid of the small index
            # [cent_row_s, cent_col_s] = cent_region[index_sm]
            # #Compare the bbox of the big index with the centroid position of the small index
            # if cent_row_s > min_row_b  and cent_row_s < max_row_b:
            #     #print("     Row analysis: cont", index_big, "is outside cont", index_sm)
            #     if cent_col_s > min_col_b  and cent_col_s < max_col_b:
            #         #print("     Col analysis: cont", index_big, "is outside cont", index_sm)
            #         int_contours.append(index_sm)
            #         #print(" >>Internal:", index_sm, ", External:", index_big)
            #else:
                #print("     Contour", index_big, "is NOT related to contour ", index_sm)

        it += 1

    for indexx in org_indexes:
        if not indexx in int_contours:
            ext_contours.append(indexx)

    #print("int_contours:", int_contours)
    #print("ext_contours:", ext_contours)

    return int_contours, ext_contours

#%% func - fillContours
def fillContours (myIm, contours):
    """
    Function to mask the image with the selected contours, fill them and get coordinates

    Parameters
    ----------
    myIm : numpy array
        Imaged being processed
    contours : list of arrays
        list of numpy array with coordinates making up the contours of the slice being processed.

    Returns
    -------
    imFilledCont : numpy array of booleans [0,1]
        Array with filled contours of myIm
    coordsXY : numpy array
        numpy array with X and Y coordinates that make up the filled contour.

    """

    if len(contours) == 0: #all(val is None for val in contours):
        imFilledCont = np.zeros([myIm.shape[0],myIm.shape[1]])
        coordsXY = np.nan

    else:
        for n, cont in enumerate(contours):
            # Create an empty image to store the masked array
            r_mask = np.zeros_like(myIm, dtype='bool')
            # Create a contour masked image by using the contour coordinates rounded to their nearest integer value
            r_mask[np.round(cont[:, 1]).astype('int'), np.round(cont[:, 0]).astype('int')] = 1
            # Fill in the holes created by the contour boundary
            r_mask = ndimage.binary_fill_holes(r_mask)
            r_mask = np.transpose(r_mask)

            if n == 0:
                resulting_mask = r_mask
            if n > 0:
                resulting_mask = np.logical_xor(resulting_mask,r_mask)

        imFilledCont = resulting_mask.astype(int)
        coordsXY = np.where(resulting_mask)
        coordsXY = np.transpose(np.asarray(coordsXY))

    return imFilledCont, coordsXY

#%% - STACK CREATION
#%% func - s3_create
def s3_create (stack_shape, im_contours, heartLayer):
    """
    Function to create a 3D array of a particular channel and group of contours

    Parameters
    ----------
    stack_shape : numpy array
        Stack shape [x,y,z]
    im_contours : str
        String indicating the imFilledContours to use.
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates

    Returns
    -------
    s3 : numpy array
        3D array with filled contours of each slice.

    """

    x_dim = stack_shape[0]
    y_dim = stack_shape[1]
    z_dim = stack_shape[2]

    s3 = np.empty((x_dim,y_dim,z_dim+2))

    for pos, keySlc in enumerate(heartLayer.keys()):
        if keySlc[0:3] == "slc":
            slcNum = int(keySlc[3:6])
            im_FilledCont = heartLayer[keySlc][im_contours]
            s3[:,:,slcNum+1] = im_FilledCont

    s3 = s3.astype('uint8')
    # print("s3_create, DONE - [uint8]")

    return s3

#%% func - ch_clean
def ch_clean (mask_s3, toClean_s3, option):
    """
    Function to clean channel using the other as a mask

    Parameters
    ----------
    mask_s3 : numpy array
        Stack of images to use as mask.
    toClean_s3 : numpy array
        Stack of images to clean.
    option : str
        "clean"/"cardiacjelly" that changes the title to use on top of the subplots.

    Returns
    -------
    toRemove_s3 : numpy array
         Resulting stack of images with regions to be removed.
    cleaned_s3 : numpy array
         Resulting stack of cleaned images.

    """

    if option == "cardiacjelly":
        print('- Extracting cardiac jelly')
    elif option == "clean":
        print('- Cleaning endocardium')

    toRemove_s3 = np.empty((toClean_s3.shape[0],toClean_s3.shape[1],toClean_s3.shape[2]))
    cleaned_s3 = np.empty((toClean_s3.shape[0],toClean_s3.shape[1],toClean_s3.shape[2]))

    bar = Bar('- Cleaning endocardial slices', max=toClean_s3.shape[2], suffix = suffix, check_tty=False, hide_cursor=False)
    for slc in range(toClean_s3.shape[2]):
        mask_slc = mask_s3[:,:,slc]
        toClean_slc = toClean_s3[:,:,slc]

        toRemove_slc = np.logical_and(toClean_slc, mask_slc)
        cleaned_slc = np.logical_xor(toClean_slc, toRemove_slc)

        toRemove_s3[:,:,slc] = toRemove_slc
        cleaned_s3[:,:,slc] = cleaned_slc
        bar.next()

    bar.finish()
    toRemove_s3 = toRemove_s3.astype('uint8')
    cleaned_s3 = cleaned_s3.astype('uint8')

    return toRemove_s3, cleaned_s3

#%% func - clean_wInvIntMyoc
def clean_wInvIntMyoc(s3_ch0_ext, s1_ch1_2clean):
    """
    Function to clean endocardial channels using the inverted internal myocardium as a mask

    Parameters
    ----------
    s3_ch0_ext : numpy array
        Stack of filled external contours of myocardial channel to invert and use as mask.
    s1_ch1_2clean : numpy array
        Stack of images to clean.

    Returns
    -------
    toRemove_s3 : numpy array
         Resulting stack of images with regions to be removed.
    s1_ch1_clean_s3 : numpy array
         Resulting stack of cleaned images.
    ch0_ext_inv_s3 : numpy array
         Resulting stack of inverted s3_ch0_ext images.

    """

    ch0_ext_inv_s3 = np.empty((s3_ch0_ext.shape[0],s3_ch0_ext.shape[1],s3_ch0_ext.shape[2]), dtype = 'uint8')
    s1_ch1_clean_s3 = np.empty((s1_ch1_2clean.shape[0],s1_ch1_2clean.shape[1],s1_ch1_2clean.shape[2]), dtype = 'uint8')
    toRemove_s3 = np.empty((s1_ch1_2clean.shape[0],s1_ch1_2clean.shape[1],s1_ch1_2clean.shape[2]), dtype = 'uint8')

    for slc in range(s1_ch1_2clean.shape[2]):
        ch0_ext_slc = s3_ch0_ext[:,:,slc]
        ch1_orig_slc = s1_ch1_2clean[:,:,slc]

        # Invert ch0_ext
        ch0_inv_slc = np.where((ch0_ext_slc==0)|(ch0_ext_slc==1), ch0_ext_slc^1, ch0_ext_slc)
        # inverted_ch0_ext AND ch1_2clean
        toRemove_slc = np.logical_and(ch1_orig_slc, ch0_inv_slc)
        # Keep only the clean bit
        cleaned_slc = np.logical_xor(ch1_orig_slc, toRemove_slc)

        ch0_ext_inv_s3[:,:,slc] = ch0_inv_slc
        s1_ch1_clean_s3[:,:,slc] = cleaned_slc
        toRemove_s3[:,:,slc] = toRemove_slc


    ch0_ext_inv_s3 = ch0_ext_inv_s3.astype('uint8')
    s1_ch1_clean_s3 = s1_ch1_clean_s3.astype('uint8')
    toRemove_s3 = toRemove_s3.astype('uint8')

    return toRemove_s3, s1_ch1_clean_s3, ch0_ext_inv_s3

#%% - SAVING
#%% func - save_s3s_fromDict
def save_s3s_fromDict (filename, chStr, stack_shape, heartLayer, dir_txtNnpy, save):
    """
    Function that saves the created s3s

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    stack_shape : numpy array
        Stack shape [x,y,z]
    heartLayer : dictionary
        Dictionary with information about the heart layer being processed including selected contours number and coordinates
    dir_txtNnpy : path
        Path to the folder where the np arrays are saved.
    save : boolean
        True if you want to save the final s3s, else False.

    Returns
    -------
    s3_int : numpy array of booleans
        Array with filled internal contours. The size of this array corresponds to the size of the original stack.
    s3_ext : numpy array of booleans
        Array with filled external contours. The size of this array corresponds to the size of the original stack.
    s3_all : numpy array of booleans
        Array with filled heart layer. The size of this array corresponds to the size of the original stack.

    """
    # Internal
    s3_title_int = filename+"_s3_"+chStr+"_int"
    s3_int =  s3_create(stack_shape, "imIntFilledCont", heartLayer)
    # External
    s3_title_ext = filename+"_s3_"+chStr+"_ext"
    s3_ext =  s3_create(stack_shape, "imExtFilledCont", heartLayer)
    # All
    s3_title_all = filename+"_s3_"+chStr+"_all"
    s3_all =  s3_create(stack_shape, "imAllFilledCont", heartLayer)
    print("- s3_created for all - [uint8]!")

    if save:
        # Export numpy file to create mesh
        # Internal
        s3_dir_int = os.path.join(dir_txtNnpy,s3_title_int)
        np.save(s3_dir_int, s3_int)
        print("- Saved s3 internal!")

        # External
        s3_dir_ext = os.path.join(dir_txtNnpy,s3_title_ext)
        np.save(s3_dir_ext, s3_ext)
        print("- Saved s3 external!")

        # All
        s3_dir_all = os.path.join(dir_txtNnpy,s3_title_all)
        np.save(s3_dir_all, s3_all)
        print("- Saved s3 all!")
        alert("countdown",1)
    alert('wohoo',1)

    return s3_int, s3_ext, s3_all

#%% func - saveStackAsNPY
def saveStackAsNPY(myStack, filename, chStr, stage, dir2save):
    """
    Function that saves the input myStack as a numpy array in the dir2save given as input

    Parameters
    ----------
    myStack : numpy array
        Stack being processed.
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    chStr : str
        Text indicating channel being processed ('ch0': myocardium/ 'ch1': endocardium).
    stage : str
        String indicating the process stage of the stack.
    dir2save : path
        Path to the folder where the npy arrays are saved.

    Returns
    -------
    None.

    """

    title = filename+"_St_"+chStr+"_"+stage
    dir2save = os.path.join(dir2save,title)
    np.save(dir2save, myStack)
    print("- Stack was saved as npy array! - "+ title)
    alert("countdown",1)

#%% func - save_s3s
def save_s3s(filename, s3_all, s3_int, s3_ext, dir_txtNnpy, layer):
    """
    Function that saves all the s3s generated as a numpy array in the dir_txtNnpy given as input

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    s3_all : numpy array of booleans
        Array with filled heart layer. The size of this array corresponds to the size of the original stack.
    s3_int : numpy array of booleans
        Array with filled internal contours. The size of this array corresponds to the size of the original stack.
    s3_ext : numpy array of booleans
        Array with filled external contours. The size of this array corresponds to the size of the original stack.
    dir_txtNnpy : path
        Path to the folder where the np arrays are saved.
    layer : str
        Text indicating the layer bwing saved to include as part of filename.

    Returns
    -------
    None.

    """

    print('- Saving '+ layer + ' s3s (all/int/ext) ')
    if not isinstance(s3_all, str):
        #All
        s3_title = filename+"_s3_"+layer
        s3_dir = os.path.join(dir_txtNnpy,s3_title)
        np.save(s3_dir, s3_all)

    if not isinstance(s3_int, str):
        # Internal
        s3_title_int = filename+"_s3_"+layer+"_int"
        s3_dir_int = os.path.join(dir_txtNnpy,s3_title_int)
        np.save(s3_dir_int, s3_int)

    if layer != 'cj' and not isinstance(s3_ext, str):
        # External
        s3_title_ext = filename+"_s3_"+layer+"_ext"
        s3_dir_ext = os.path.join(dir_txtNnpy,s3_title_ext)
        np.save(s3_dir_ext, s3_ext)

    alert('wohoo',1)
    print('- All s3s have been saved for '+layer+'!')

#%% func - save_s3
def save_s3(filename, s3, dir_txtNnpy, layer):
    """
    Function that saves all the s3s generated as a numpy array in the dir_txtNnpy given as input

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    s3 : numpy array of booleans
        Array with heart layer to save. The size of this array corresponds to the size of the original stack.
    dir_txtNnpy : path
        Path to the folder where the np arrays are saved.
    layer : str
        Text indicating the layer being saved to include as part of filename (chX_int or chX_ext or chX).

    Returns
    -------
    None.

    """

    s3_title = filename+"_s3_"+layer
    s3_dir = os.path.join(dir_txtNnpy,s3_title)
    np.save(s3_dir, s3)

#%% - ALERT WHEN IMPORTED
print ("IMPORTED: morphoHeart_funcContours")
alert('jump',1)
