# -*- coding: utf-8 -*-
"""
morphoHeart - A. CONTOURS

@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
from time import perf_counter

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
initial_runD = True

#%% Start A_Contours
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcContours as fcCont
    from morphoHeart_modules.morphoHeart_funcMeshes import createLayerMesh

    # Creating global variable
    global heartLayer

    #%% Get directories and file
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse)
    # Get file to process and directories
    folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse, 6. dir_LS_Folder selected
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')
    # Import the metadata to know pixel size and distance between slices
    xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
    res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]

    #%% i. IMPORT
    # Get channel to process
    channel = fcCont.selectChannel()
    # Import stack
    stack, stack_o, file, stack_shape  =fcCont.main_importImages(filename, channel, directories, n_rows = 6)

    #%% ii. CLOSE CONTOURS
    q_ABC = fcBasics.ask4input('Do you want to automatically close contours, manually close contours or \n\tclose inflow/outflow tracts of this stack? [0]:no/[1]:yes: ',bool)
    if q_ABC:
        ticABC = perf_counter()
        # >> Automatically close contours
        stack_closed, processDict, done_autom = fcCont.main_automCloseCont(filename, channel, directories, stack, plotEvery = 100, n_rows = 6)
        # >> Manually close remaining contours
        stack_closed, processDict, done_manual = fcCont.main_manuallyCloseContours(filename, channel, directories, stack_closed, stack_o, processDict, 6)
        if done_autom and done_manual:
            # >> Close inflow and outflow tracts
            stack_closed, processDict, done_infOutf = fcCont.main_closeInfAndOutfTract(filename, channel, directories, stack_closed, processDict, 6)
            if done_infOutf:
                q_ABC_done = True
                # >> Save stack before continuing to selectContours
                fcCont.saveStackAsNPY(stack_closed, filename, channel, 'closedCJ', directories[1])
                first = processDict[channel]['G-Slc_tissueLayerFirst']
                last = processDict[channel]['G-Slc_tissueLayerLast']
                fcCont.showGridContours(myStack = stack_closed, slices = (first,last), n_rows = 6)

                tocABC = perf_counter()
                fcBasics.printTime(ticABC, tocABC, 'close contours')
    else:
        stack_closed = stack

    #%% iii. SELECT CONTOURS
    if q_ABC_done:
        exit_txt = False
        q_D = fcBasics.ask4input('Do you want to select the layer contours or modify the selected contours for '+channel+'? [0]:no/[1]:yes: ',bool)
        if q_D:
            if initial_runD:
                ticD = perf_counter()
                # Define number of contours per slice
                slcCont_o, numCont_o = fcCont.getSlicesContNum(stack_closed)
                # Show tuples according to slcContours
                tuple_slc, numCont, slcCont = fcCont.tuple_pairs(numCont_o, slcCont_o, False, 20)
                # Creation of heartLayer dictionary with basic info
                heartLayer, version_N, update_dict = fcCont.dictCreation(filename = filename, file = file, minLenContour = 250,
                                                        chStr = channel, slcCont = slcCont_o, numCont = numCont_o, xy_Scaling_um = xy_Scaling_um,
                                                        z_Scaling_um = z_Scaling_um, heartLayer = dict(), update = False, curr_v= 0)
                initial_runD = False
            else:
                tuple_slc, numCont = fcCont.update_tuple_pair(tuple_slc, numCont, heartLayer)
            # Fill dictionary of selected contours
            while len(tuple_slc) != 0:
                heartLayer, tuple_slc, numCont, exit_txt = fcCont.dictFill (tuple_slc = tuple_slc, numContours = numCont,
                                                    stack = stack_closed, chStr = channel, heartLayer = heartLayer, minLenContour = 250)
                if exit_txt:
                    print('EXIT!!');
                    break
            if len(tuple_slc) == 0:
                heartLayer = fcCont.askModifyDict(filename, stack, heartLayer, channel)
                q_fullDict = fcBasics.ask4input('Is the dictionary full for '+channel+'? [0]:no/[1]:yes: ',bool)
                if q_fullDict:
                    #% Create dictionary to save and save it
                    heartLayer2S = fcCont.smallDict2Save(heartLayer)
                    fcBasics.saveDict(filename, heartLayer2S , "heartlayer2S_"+channel, directories[0])
                    q_save_s3s = fcBasics.ask4input("Do you want to save the stacks created? [0]:no/[1]:yes: ",bool)
                    _, _, s3_all = fcCont.save_s3s_fromDict(filename = filename, chStr = channel, stack_shape = stack_shape,
                                                     heartLayer = heartLayer, dir_txtNnpy = directories[1], save = q_save_s3s)
                    #Find surfaces in 3D for channel analysed
                    mesh = createLayerMesh(filename = filename, s3 = s3_all, resolution = res, layer = 'Channel '+str(channel),
                                                name = 'Channel '+str(channel), colour = 'cornflowerblue', alpha = 1, plotshow=True)
                tocD = perf_counter()
                fcBasics.printTime(ticD, tocD, 'select contours')
    else:
        print('- You need to have closed all the contours and inflow/outflow tracts of the stack to continue with the selection of contours process')

#%% Other functions to run individually in case
others = False
if others:
    #%% Save the stack you have processed so far
    fcCont.saveStackAsNPY(stack_closed, filename, channel, 'closedCJ', directories[1])

    #%% Manually close the contours or clean any slices you might have missed
    stack_closed, processDict = fcCont.main_manuallyCloseContours(filename, channel, directories, stack_closed, stack_o, processDict, 6)

    #%% Plot the contours from slices 'first' to 'last' in a grid
    first = 0
    last = 281 # any number within the number of slices in the stack
    fcCont.showGridContours(myStack = stack_closed, slices = (first,last), n_rows = 6)

    #%% Plot the selected contours (masks) of all the slices to double check
    fcCont.plotSelectedContours(imageEvery = 1, stack = stack_closed, heartLayer = heartLayer)
