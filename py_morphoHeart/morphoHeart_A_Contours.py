# -*- coding: utf-8 -*-
"""
morphoHeart - A. CONTOURS
Welcome to the first code of morphoHeart!
The objective of this script is to close and select the contours of the myocardial and endocardial channels of the heart
so that in future scripts extraction of the cardiac jelly volume and morphological measurements of the heart and each 
of its tissue layers can be carried out! At the end of this script you will end with three numpy arrays containing 
all the information of the internal, external and tissue layer contours of either the myocardium or the endocardium. 
This scipt thus need to be run twice, one time for each tissue layer. You can find the final arrays saved in a folder 
called txt_npy within the main folder of the heart being processed.

Happy contour closing and selecting!

@author: Juliana Sanchez-Posada
Version: 13th April, 2021
"""

#%% Importing python packages
import os
from time import perf_counter
from vedo import Plotter

init = False
# Verify working dir
def setWorkingDir (root_path, init):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    # init = True
    print("Current working directory: {0}".format(os.getcwd()))

    return root_path, init

root_path, init = setWorkingDir(os.getcwd(),init)

#%% Start A_Contours
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcContours as fcCont
    from morphoHeart_modules.morphoHeart_funcMeshes import createLayerMesh

    # Creating global variable
    global heartLayer, stack_closed, processDict

    #%% Set-up variables to start running code 
    initial_runABC = True; initial_runD = True
    q_ABC_done = False; q_fullDict = False
    n_rows = fcBasics.ask4input('Enter the number of rows you want to see when plotting slices with contours in grid \n\t (Note: if running in analysis computer the recommended is max 5):', int)
    
    #%% SELECT FILE AND GET METADATA 
    #   This section allows the user to select file to process and get its main directories and metadata
    #   ================================================================================================================
    
    # Get main directories 
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse)
    # Get file to process and directories
    folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse, 6. dir_LS_Folder selected
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')
    # Import the metadata to know pixel size and distance between slices
    xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
    res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]

    #%% SELECT CHANNEL AND LOAD FILES 
    #   This section allows the user to select the channel/heart layer to process and load corresponding files [tif/npy]
    #   [ch0: myocardium/ch1: endocardium]
    #   ================================================================================================================
    
    # Get channel to process
    channel = fcCont.selectChannel()
    # Import stack
    stack_closed, stack_o, stack_m, file, stack_shape, processDict = fcCont.main_importImages(filename, channel, directories, n_rows = n_rows)

    #%% ii. CLOSE CONTOURS
    #   This section allows the user to run either/all of the three different functions to close the channel's contours.
    #   1. Automatic closure of the contours
    #   2. Manual closure / cleaning of the contours
    #   3. Inflow/Outflow tract closure
    #   Throughout the process a dictionary will be saved keeping track of the processes the user has completed for each 
    #   of the file's channel and at the end will save the final closed stack as an numpy array.
    #   NOTE: Intermediate steps of this process can be saved.
    #   ================================================================================================================

    q_ABC = fcBasics.ask4input('Do you want to run any of these processes: \n\t- Automatically close contours \n\t- Manually close contours or \n\t- Close inflow/outflow tracts of this stack? \n\t >[0]:no/[1]:yes: ',bool)
    if q_ABC:
        ticABC = perf_counter()
        # >> Automatically close contours
        stack_closed, processDict, done_autom = fcCont.main_automCloseCont(filename, channel, directories, stack_closed, 
                                                                           plotEvery = 100, n_rows = n_rows, processDict = processDict)
        # >> Manually close remaining contours
        stack_closed, processDict, done_manual = fcCont.main_manuallyCloseContours(filename, channel, directories, 
                                                               stack_closed, stack_o, stack_m, processDict, n_rows = n_rows)
        if done_autom and done_manual:
            # >> Close inflow and outflow tracts
            stack_closed, processDict, done_infOutf = fcCont.main_closeInfAndOutfTract(filename, channel, 
                                                           directories, stack_closed, processDict, n_rows = n_rows)
            # >> Manually close additional contours (if needed)
            stack_closed, processDict, done_manual = fcCont.main_manuallyCloseContours(filename, channel, 
                                       directories, stack_closed, stack_o, stack_m, processDict, n_rows = n_rows, checking = True)
            # if done_infOutf:
            q_ABC_done = fcBasics.ask4input('Checking: Are you done closing the contours and inflow/outflow tracts? [0]:no/[1]:yes!:', bool)
            if q_ABC_done:
                # >> Save stack before continuing to selectContours
                fcCont.saveStackAsNPY(stack_closed, filename, channel, 'closedCJ', directories[1])
                first = processDict[channel]['G-Slc_tissueLayerFirst']
                last = processDict[channel]['G-Slc_tissueLayerLast']
                fcCont.showGridContours(myStack = stack_closed, slices = (first,last), n_rows = n_rows)
                tocABC = perf_counter()
                fcBasics.printTime(ticABC, tocABC, 'close contours')
    else:
        q_ABC_done = fcBasics.ask4input('Checking: Are you done closing the contours and inflow/outflow tracts? [0]:no/[1]:yes!:', bool)

    #%% iii. SELECT CONTOURS
    #   This section allows the user to select the internal and external contours of the heart layer being processed. 
    #   Before starting the user should scan through the plotted closed slices and write down a list with:
    #   1. all the (first) slices in which a group of consecutive slices have the same number of contours 
    #   2. the corresponding number of contours in each group
    #   At the beggininng of the process the user will be asked to enter this information. Then the contour groups will 
    #   be subdivided into groups containing each a maximum of 40 slices. The user will be asked to select the internal 
    #   and external contours for the first slice of each group and the next 39 (or less) slices the contours will be
    #   automatically selected. 
    #   When the contours for all slices have been selected, the user will be able to modify the selected contours in 
    #   case some of them were incorrectly automatically selected.
    #   At the end of this step, numpy arrays with the masks of the internal, external and layer contours will be saved
    #   as well as a dictionary containing information about the selected contours per slice.
    #   Finally a surface reconstruction of the heart tissue layer will be created and plotted.
    #   NOTE: Intermediate step of this process cannot be saved, so please make sure to run it 
    #   completely before closing. 
    #   ================================================================================================================
    
    if q_ABC_done:
        exit_txt = False
        q_D = fcBasics.ask4input('Do you want to select the layer contours or modify the selected contours for '+channel+'? [0]:no/[1]:yes: ',bool)
        if q_D:
            while not q_fullDict:
                if initial_runD:
                    ticD = perf_counter()
                    # Define number of contours per slice
                    slcCont_o, numCont_o = fcCont.getSlicesContNum(stack_closed)
                    # Show tuples according to slcContours
                    tuple_slc, numCont, slcCont = fcCont.tuple_pairs(numCont_o, slcCont_o, False, 40)
                    # Creation of heartLayer dictionary with basic info
                    heartLayer, version_N, update_dict = fcCont.dictCreation(filename = filename, file = file, minLenContour = 250,
                                                            chStr = channel, slcCont = slcCont_o, numCont = numCont_o, 
                                                            xy_Scaling_um = xy_Scaling_um,z_Scaling_um = z_Scaling_um, 
                                                            heartLayer = dict(), update = False, curr_v= 0)
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
                    heartLayer, stack, processDict = fcCont.askModifyDict(filename, stack_closed, heartLayer, channel, processDict, directories)
                    q_fullDict = fcBasics.ask4input('Is the dictionary full for '+channel+'? [0]:no/[1]:yes: ',bool)
                    if q_fullDict:
                        #% Create dictionary to save and save it
                        heartLayer2S = fcCont.smallDict2Save(heartLayer)
                        fcBasics.saveDict(filename, heartLayer2S , "heartlayer2S_"+channel, directories[0])
                        _, _, s3_all = fcCont.save_s3s_fromDict(filename = filename, chStr = channel, stack_shape = stack_shape,
                                                         heartLayer = heartLayer, dir_txtNnpy = directories[1], save = True)
                        #Find surfaces in 3D for channel analysed
                        mesh = createLayerMesh(filename = filename, s3 = s3_all, resolution = res, layer = 'Channel '+str(channel),
                                                    name = 'Channel '+str(channel), colour = 'cornflowerblue', alpha = 1, plotshow=True)
                    tocD = perf_counter()
                    fcBasics.printTime(ticD, tocD, 'select contours')
    else:
        print('- You need to have closed all the contours and inflow/outflow tracts of the stack to continue with the selection-of-contours process')

#%% OTHER FUNCTIONS
#   This section allows the user to re-run individual functions in case they are needed. Ignore otherwise. 
#   ====================================================================================================================
others = False
if others:
    #%% Save PNGs of the closed stack as individual PNGs in a folder inside Im_LSXX_FXX folder
    name = fcBasics.ask4input('Enter reference version (integer number) of the slices you want to save as PNGs: ', int)
    fcCont.savePltContours(dir_ims2Analyse = directories[5], filename = filename,
                                      myStack = stack_closed, chStr = channel,
                                      slices = (0,len(stack_closed)), contVersion = '0')
    
    #%% Save the stack you have processed so far
    fcCont.saveStackAsNPY(stack_closed, filename, channel, 'closedCJ', directories[1])

    #%% Manually close the contours or clean any slices you might have missed
    stack_closed, processDict, _ = fcCont.main_manuallyCloseContours(filename, channel, directories, stack_closed, 
                                                                     stack_o, stack_m, processDict, 7, True)

    #%% Plot the contours from slices 'first' to 'last' in a grid
    # Enter in 'first' and 'last' any slice number within the number of slices in the stack
    first = 110
    last = 130
    fcCont.showGridContours(myStack = stack_m, slices = (first,last+1), n_rows = n_rows)
    
    #%% Plot all the contours of the stack in a grid
    first = 0
    last = stack_closed.shape[0]
    fcCont.showGridContours(myStack = stack_closed, slices = (first,last), n_rows = n_rows)

    #%% Plot the selected contours (masks) of a selected group of slices to double check
    fcCont.plotSelectedContours(imageEvery = 1, stack = stack_closed, heartLayer = heartLayer)

    #%% Plot again resulting mesh
    vp = Plotter(N=1, axes =13)
    vp.show(mesh, at=0,interactive = True)

#%% Init
init = True