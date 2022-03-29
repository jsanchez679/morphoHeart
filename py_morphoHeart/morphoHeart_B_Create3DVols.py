# -*- coding: utf-8 -*-
"""
morphoHeart - B. CREATE 3D VOLUMES AND MESHES TO EXTRACT CENTRELINE
Welcome to the second code of morphoHeart!
If you are running this code you must have already closed and selected the contours of both the myocardium and 
endocardium of one of your hearts and are looking forward to see how the cardiac jelly of this heart looks in 3D!!!
The objective of this script is to create the volume reconstructions of both the myocardium and endocardium,  
extract the cardiac jelly and create its 3D volume and save a mesh (or meshes) from which to extract the heart 
centreline in the next script.
To do this we need to initially clean the endocardium using the internal contours of the myocardium. Next,
cut the inflow and outflow tracts of the endocardium (particularly remove the aortic arches and the common
cardinal vein), and using the same method cut also the myocardium so that we end up analysing only the heart tissues. 
Finally we smooth, clean and save the mesh (or meshes) from which we are going to extract the centreline. 
At the end of this script you will end with all the created meshes saved in the 'meshes' folder and the centreline mesh(es) 
saved in the 'centreline' folder.

Important note: 
Intermediate steps of this process can be saved, but I suggest you to run it complete as it doesn't take too long.

Happy running!

@author: Juliana Sanchez-Posada
Version: 15th November, 2021
"""

#%% Importing python packages
import os
from time import perf_counter
# from vedo import *
from vedo import embedWindow, settings, Plotter, Text2D
embedWindow(False)
settings.legendSize = .3
# settings.legendPos = 1
settings.legendFont="VTK"

# Verify working dir
def setWorkingDir (root_path, init = False):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    init = not bool(int(input('Do you want to execute the script all at once or run it by cells? \n\t[0]: all at once (recommended if you are already familiar with the script)\n\t[1]: by cells (recommended if you are NOT yet familiar with the script OR \n\t\t you just need to create the mesh(es) to extract centreline). >>>: ')))
    print("Current working directory: {0}".format(os.getcwd()))
    if not init: 
        print('\nIMPORTANT NOTES:\n- Remember to start running from cell %Start B_Create3DVols.\n- NEVER run as an individual cell the cell called %Importing python packages')
    return root_path, init

root_path, init = setWorkingDir(os.getcwd())
save = True

#%% Start B_Create3DVols (*** 1)
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcContours as fcCont
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
    tic = perf_counter()

    #%% SELECT FILE AND GET METADATA (*** 2)
    #   This section allows the user to select file to process, get its main directories and metadata,
    #   and define some properties
    #   ================================================================================================================

    # Get main directories
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse)
    # Get file to process and directories
    folder, df_file, file_num, blind = fcBasics.selectFile(df_dataset); filename = folder[0:-3]; dORv = filename[9:10]
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse, 6. dir_LS_Folder selected
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')
    # Import the metadata to know pixel size and distance between slices
    xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
    res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]
    # Initialise variables
    plotshow = False
    
    [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes] = fcBasics.import_dicts('mH_B1', filename, directories)
    # try: 
    #     # Import dictionaries
    #     dicts = fcBasics.loadDicts(filename = filename, dicts_name = ['dict_obj'], directories = [directories[0]])
    #     [dict_planes, dict_pts, dict_kspl, dict_colour] = fcMeshes.splitDicts(dicts[0])
    # except: 
    #     dict_planes = dict(); dict_pts = dict(); dict_kspl = dict(); dict_colour = dict()
    txt = Text2D(filename, c="k", font= 'VTK'); initial_smooth = True
    if dORv == 'D' or 'CJ' in filename:
        azimuth = -90
    else: 
        azimuth = 0
        
    #%% SELECT NUMBER OF LAYERS AND CHAMBERS TO ANALYSE
    tissLay2p, dict_process = fcBasics.selectLayer2process('mH_B1',filename, directories)
    
    # chambers2p = fcBasics.selectChambers2process()
    
    #%% CLEAN ENDOCARDIUM
    #   This section will import the external contour mask of the endocardium and the endocardial mask and clean them 
    #   using the user selected mask (1: the myocardial tissue mask, or 2: the inverted internal myocardial mask, 
    #   which is the recommended option). Use the 1st option if you want to keep pieces of the endocardium that might be 
    #   present outside of the myocardium, but keep in mind that these pieces will be disconnected from the big 
    #   endocardial tissue layer.
    #   At the end of this step, numpy arrays of the cleaned masks will be saved and a 3D interactive plot with the
    #   myocardium and endocardium (original and clean), will pop-up.
    #   ================================================================================================================

    if dict_process['Processes']['B_C3DV']['CleanEndocardium']:
        # Load stacks
        s3s, stackShape = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1],
                                            end_name = ['ch0_int','ch0_all','ch1_ext','ch1_all'])
        s3_ch0_int, s3_ch0, s3_ch1_ext, s3_ch1 = s3s
    
        #% Clean endocardium
        # Plot s3s - channel side by side (if plotshow = True)
        fcCont.plt_s3(start_slc = 0, end_slc = s3_ch0.shape[2]-1, im_every = 10,
                          s3_int = s3_ch0, s3_ext = s3_ch1, plotshow = plotshow, option = "both ch")
        # Create Ext Myocardial mesh - Ch0
        myoc_o = fcMeshes. createLayerMesh(filename = filename, s3 = s3_ch0, resolution = res, layer = 'Myoc',
                                           name = 'Myocardium', colour = 'darkcyan', alpha = 0.1, plotshow = plotshow)
        # Create Ext Endocardial mesh - Ch1
        endo_o = fcMeshes.createExtLayerMesh(filename = filename, s3_ext = s3_ch1_ext, resolution = res, layer = 'Endo',
                                             info = 'Original', plotshow = plotshow)
        
        m_selected = False
        while not m_selected: 
            method = fcBasics.ask4input('Select the mask you would like to use to clean the endocardial layer: \n\t[1]: Just the myocardial tissue layer \n\t[2]: (Recommended) The inverted internal myocardium (more profound cleaning). >>>: ', int)
            # Clean Endocardium - Option 1 (Using just myocardial layer, so that the two tissue layers don't colocalise in space)
            if method == 1:
                m_selected = True
                s3_endo_rem, s3_ch1_cl = fcCont.ch_clean(s3_ch0, s3_ch1, option = "clean")
                # Plot s3s of clean endocardium - (if plotshow = True)
                fcCont.ch_clean_plt(mask_s3 = s3_ch0, toClean_s3 = s3_ch1, toRemove_s3 = s3_endo_rem, 
                                          cleaned_s3 = s3_ch1_cl, plotshow = plotshow, im_every = 15, option = "clean")
                # Clean Ext Endocardium
                _, s3_ch1_ext_cl = fcCont.ch_clean(mask_s3 = s3_ch0, toClean_s3 = s3_ch1_ext, option = "clean")
                extractLargest = fcBasics.ask4input('Do you want to just keep the largest surface of the cleaned endocardium? \n\t[0]: no, keep all the pieces\n\t[1]: yes, just keep the biggest surface : ', bool)
            # Clean Endocardium - Option 2 (Using inverted internal myocardium).
            elif method == 2:
                m_selected = True
                s3_endo_rem, s3_ch1_cl, s3_invIntMyoc = fcCont.clean_wInvIntMyoc(s3_ch0_int, s3_ch1)
                # Plot s3s of clean endocardium - (if plotshow = True)
                fcCont.ch_clean_plt(mask_s3 = s3_invIntMyoc, toClean_s3 = s3_ch1, toRemove_s3 = s3_endo_rem,
                                          cleaned_s3 = s3_ch1_cl, plotshow = plotshow, im_every = 15, option = "clean")
                # Clean Ext Endocardium
                _, s3_ch1_ext_cl, _ = fcCont.clean_wInvIntMyoc(s3_ch0_int, s3_ch1_ext)
                extractLargest = True
                
        fcCont.save_s3(filename = filename, s3 = s3_ch1_cl, dir_txtNnpy = directories[1], layer = 'ch1_cut')
        fcCont.save_s3(filename = filename, s3 = s3_ch1_ext_cl, dir_txtNnpy = directories[1], layer = 'ch1_cut_ext')
        # Re-Create Ext Endocardial mesh - Ch1
        endo_ext = fcMeshes.createExtLayerMesh(filename = filename, s3_ext = s3_ch1_ext_cl, resolution = res, 
                                               layer = 'Endo', info = 'Cleaned', plotshow = plotshow, 
                                               extractLargest = extractLargest)
    
        # Plot Ext Myocardium and Clean Ext Endocardium
        myoc_o.alpha(0.1); endo_o.color('mediumorchid').legend('Orig.Ext.Endo'); settings.legendSize = .2
        vp = Plotter(N=4, axes=10)
        vp.show(endo_o, txt, at=0, zoom=1)
        vp.show(myoc_o, endo_o, at=1, zoom=1)
        vp.show(myoc_o, endo_ext, at=2, zoom=1)
        vp.show(endo_ext, at=3, zoom=1, azimuth = azimuth, interactive=True)

        del s3_ch0_int, s3_ch0, s3_ch1_ext, s3_ch1, s3_endo_rem, s3_ch1_cl, s3_ch1_ext_cl, 
        
    else: 
        if 'Myocardium' in tissLay2p:
            # Create Ext Myocardial mesh - Ch0
            myoc_o = fcMeshes. createLayerMesh(filename = filename, s3 = s3_ch0, resolution = res, layer = 'Myoc',
                                           name = 'Myocardium', colour = 'darkcyan', alpha = 0.1, plotshow = plotshow)
        if 'Endocardium' in tissLay2p:
            # Re-Create Ext Endocardial mesh - Ch1
            endo_ext = fcMeshes.createExtLayerMesh(filename = filename, s3_ext = s3_ch1_ext_cl, resolution = res, 
                                               layer = 'Endo', info = 'Cleaned', plotshow = plotshow,
                                               extractLargest = extractLargest)
    
    #%% CUT INFLOW AND OUTFLOW TRACTS OF BOTH TISSUE LAYERS
    #   This section will allow the user to define planes to cut both the inflow and outflow regions of the myocardial
    #   and endocardial tissue layers. This process is needed as it will remove the aortic arches and create sharp cuts
    #   at both ends of the heart layers to obtain a clean reconstruction of the cardiac jelly. At the end of the process
    #   volume reconstructions of the new cut meshes will be generated.
    #   Initially, a 3D interactive plot with both heart layers will pop-up to allow the user to see if a cut is needed
    #   in the inflow tract and if so, select which of the two (if not the two) tissue layers wants to be cut.
    #   Next, a new 3D interactive plot will pop-up to allow the user rotate and translate a pre-defined grey plane to
    #   the position where the cut wants to be made. When the user is happy with the selected plane, he/she should close
    #   the window.
    #   A new pop-up window will appear showing the selected plane and both tissue layers. When the window is closed the
    #   user will be asked if he/she is happy with the defined plane. If so, the user will need to go through the same
    #   process for the outflow tract. Once the user had given all the input, the selected tissue layers will be cut,
    #   and a pop-up window showing the resulting meshes will be shown.
    #   NOTE: If the result of the cuts performed want to be seen agan, un-comment the commented plot.
    #   ================================================================================================================

    if dict_process['Processes']['B_C3DV']['CutIFandOFTs']:
        #% Cut inflow and outflow regions of endocardium and myocardium and save all final s3s
        meshes_cut, dict_planes = fcMeshes.selectCutS3sOptMxLoad(filename = filename,
                                            m_endo = endo_ext, m_myoc = myoc_o,
                                            dict_planes = dict_planes, resolution = res,
                                            dir_txtNnpy = directories[1], save = save)
    
        myoc_cut, myoc_cut_int, myoc_cut_ext, endo_cut, endo_cut_int, endo_cut_ext = meshes_cut
    
        # settings.legendSize = .3
        # vp = Plotter(N=6, axes=10)
        # vp.show(myoc_cut, txt, at=0, zoom=1)
        # vp.show(myoc_cut_int, at=1, zoom=1)
        # vp.show(myoc_cut_ext, at=2,  zoom=1)
        # vp.show(endo_cut, at=3,  zoom=1)
        # vp.show(endo_cut_int, at=4)
        # vp.show(endo_cut_ext, at=5, zoom=1, azimuth = azimuth, interactive=True)
    
        if save:
            dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [myoc_cut, myoc_cut_int, myoc_cut_ext, endo_cut, endo_cut_int, endo_cut_ext],
                                names = ['myoc', 'myoc_int', 'myoc_ext', 'endo', 'endo_int', 'endo_ext'], dict_colour = dict_colour,
                                dir_stl = directories[2], extension = 'vtk')

    #%% EXTRACT CARDIAC JELLY FROM CONTOURS
    #   This section will extract the contours of the cardiac jelly 'subtracting' the filled external contours mask of
    #   the endocardium from the filled internal contours mask of the myocardium, and return volume reconstructions of
    #   the cardiac jelly. At the end of this process, a 3D interactive plot will pop-up, with the 3D reconstructions of
    #   the external and internal cardiac jelly as well as the cardiac jelly alone. Finally, a 3D interactive plot with
    #   all the heart tissue layers will appear.
    #   ================================================================================================================

    if dict_process['Processes']['B_C3DV']['ExtractCJ']:
        # Extract cardiac jelly from contours
        try: 
            [s3_ch0_int_cut, s3_ch1_ext_cut], _ = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1],
                                                                    end_name = ['ch0_cut_int','ch1_cut_ext' ])
        except:
            [s3_ch0_int_cut, s3_ch1_ext_cut], _ = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1],
                                                                    end_name = ['ch0_int','ch1_cut_ext' ])
            
        # Plot s3s - channel side by side (if plotshow = True)
        fcCont.plt_s3(start_slc = 0, end_slc = s3_ch0_int_cut.shape[2]-1, im_every = 10,
                          s3_int = s3_ch0_int_cut, s3_ext = s3_ch1_ext_cut, plotshow = plotshow, option = "cardiacjelly")
        # Get cardiac jelly per slice
        s3_cjIn, s3_cj = fcCont.ch_clean(mask_s3 = s3_ch1_ext_cut, toClean_s3 = s3_ch0_int_cut, option = "cardiacjelly")
        fcCont.ch_clean_plt(mask_s3 = s3_ch0_int_cut, toClean_s3 = s3_ch1_ext_cut, toRemove_s3 = s3_cjIn,
                                  cleaned_s3 = s3_cj, plotshow = plotshow, im_every = 5, option = "cardiacjelly")
        # Create meshes - CJ, CJ Internal surface and CJ External surface
        cj_all, cj_in, cj_out = fcMeshes.createAll3LayerMeshes(filename = filename, s3_all = s3_cj, s3_in = s3_cjIn,
                                            s3_out = s3_ch0_int_cut, resolution = res, layer = 'CJ')
        # Save Heart layers (cj) and s3s
        if save:
            fcCont.save_s3s(filename = filename, s3_all = s3_cj, s3_int = s3_cjIn, s3_ext = s3_ch0_int_cut,
                                    dir_txtNnpy = directories[1], layer = 'cj')
            dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [cj_all, cj_in, cj_out],
                                names = ['cj', 'cj_in', 'cj_out'], dict_colour = dict_colour,
                                dir_stl = directories[2], extension = 'vtk')
    
        # Plot all layers
        settings.legendSize = .3
        vp = Plotter(N=6, axes=13)
        vp.show(myoc_cut, txt, at=0, zoom=1)
        vp.show(endo_cut, at=1, zoom=1)
        vp.show(myoc_cut_int.color('turquoise'), at=2,  zoom=1)
        vp.show(cj_all.clone().alpha(0.05), at=3,  zoom=1)
        vp.show(myoc_cut.clone().alpha(0.1), endo_cut, cj_all, at=4)
        vp.show(endo_cut_ext.color('orchid'), at=5, zoom=1, azimuth = azimuth, interactive=True)
    
        del s3_ch0_int_cut, s3_ch1_ext_cut, s3_cjIn, s3_cj

    #%% MEASURE SURFACE AREA
    #   Now that all the volumetric reconstructions of all the heart tissue layers have been generated, surface area
    #   measurements will be taken and added to the measurements dataframe.
    #   ================================================================================================================
    
    if dict_process['Processes']['B_C3DV']['MeasureSurfAreaWholeTissues']:
        # Save SurfArea and Int/Ext Volume of meshes to dataframe and save df
        df_res = fcMeshes.addSurfArea2df(df_res = df_file,  file_num = file_num,
                                meshes = [myoc_cut, myoc_cut_int, myoc_cut_ext, endo_cut, endo_cut_int, endo_cut_ext, cj_all, cj_in, cj_out])
    
        df_res = fcMeshes.addIntExtVol2df(df_res = df_res, file_num = file_num, meshes = [myoc_cut_int, myoc_cut_ext, endo_cut_int, endo_cut_ext])
        if save:
            fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
    
        # Save
        if save:
            dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour],
                                             names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour'],
                                             dir2save = directories[0])
        
    #%% CREATE, CUT AND EXPORT MESHES TO OBTAIN CENTRELINE (*** 3)
    #   This section allows the user to create one or two meshes from which the centreline will be extracted. As the
    #   package used to create the centrelines (vmtk) needs a smooth mesh with blunt cuts in the inflow and outflow
    #   tracts, this section will smooth the meshes and ask the user to define planes to cut both the inflow and outflow
    #   tracts (following a similar process to the one previously used to cut the meshes).
    #   This section will save the resulting meshes (smoothed and cut).
    #   NOTE: If when running this process the computer freezes and the kernel needs to be restarted, restart the kernel, 
    #   make sure the rootpath is set correctly by right-cliking on this script's file name (tab on top) and selecting 
    #   'Set console working directory'. Then run the sections of the code that have a (***) at the end of the title in 
    #   the order given by the numbers (first code to run is at the end of the script). 
    #   This process will allow you just to create the meshes for the centreline without having to run again all the code.
    #   ================================================================================================================
    
    if dict_process['Processes']['B_C3DV']['CreateMeshes4CL']:
        if 'myoc_cut_int' not in locals():
            # Import dictionaries
            dicts = fcBasics.loadDicts(filename = filename, dicts_name = ['dict_obj'], directories = [directories[0]])
            dict_obj = fcMeshes.splitDicts(dicts[0])
            if len(dict_obj) == 4:
                [dict_planes, dict_pts, dict_kspl, dict_colour] = dict_obj
            else:
                [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes] = dict_obj
        
            # Import meshes
            initial_smooth = True
            [myoc_cut_int, endo_cut_ext] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc_int', 'myoc_ext'],
                                                                      extension = 'vtk', dir_stl = directories[2],
                                                                      alpha = [1,1], dict_colour = dict_colour)
    
        # Create meshes for centreline
        happyWithCuts = False
        while not happyWithCuts:
            if initial_smooth:
                meshes4cl, names_exp = fcMeshes.createMeshes4CL(filename = filename, meshes = [myoc_cut_int, endo_cut_ext],
                                                                plotshow = True); initial_smooth = False
                meshes4cl_o = meshes4cl.copy()
            else: 
                meshes4cl = []; meshes4cl_names = []
                q_selectMeshes = fcBasics.ask4input('Select the meshes you would like to export to extract centreline \n\t\t[0]: Int.Myocardium/[1]: Ext.Endocardium/[2]: both: ', int)
                if q_selectMeshes in [0,2]:
                    meshes4cl.append(meshes4cl_o[0])
                    meshes4cl_names.append(names_exp[0])
                if q_selectMeshes in [1,2]:
                    meshes4cl.append(meshes4cl_o[1])
                    meshes4cl_names.append(names_exp[1])
        
            # Cut inflow and outflow tract of meshes for centreline
            meshes4clf, dicts = fcMeshes.cutMeshes4CL(filename = filename, meshes = meshes4cl,
                                                      cuts = ['inflow', 'outflow'], cut_direction = [True, False],
                                                      dicts = [dict_planes, dict_pts, dict_kspl], plotshow = True)
            happyWithCuts = fcBasics.ask4input('Are you happy with the cuts made? \n\t[0]: no, I would like to recut the meshes\n\t[1]: yes, continue with next steps! >>>:', bool)
        
        # Save
        if save:
            dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = meshes4clf, names = names_exp,
                                          dict_colour = dict_colour, dir_stl = directories[3], extension='stl')
    
        # Plot and save all meshes
        settings.legendSize = .2
        if 'myoc_cut' in locals():
            vp = Plotter(N=4+len(meshes4clf), axes=13)
            vp.show(myoc_cut, txt, at=0, zoom=1)
            vp.show(endo_cut, at=1, zoom=1)
            vp.show(meshes4clf[0], at=2,  zoom=1)
            vp.show(cj_all, at=3,  zoom=1)
            if len(meshes4clf) == 1:
                vp.show(myoc_cut, endo_cut, cj_all, at=4, zoom=1, azimuth = azimuth, interactive=True)
            else:
                vp.show(myoc_cut, endo_cut, cj_all, at=4)
                vp.show(meshes4clf[1], at=5, zoom=1, azimuth = azimuth, interactive=True)
        else: 
            vp = Plotter(N=len(meshes4clf), axes=13)
            if len(meshes4clf) == 1:
                vp.show(meshes4clf[0], at=0, zoom=1, azimuth = azimuth, interactive=True)
            else:
                vp.show(meshes4clf[0], at=0, zoom=1)
                vp.show(meshes4clf[1], at=1, zoom=1, azimuth = azimuth, interactive=True)

    #%% SAVE ALL AND PRINT INSTRUCTIONS (*** 4)
    #   This section allows the user to save the dictionaries created with all the planes, splines and spheres created
    #   and print the instructions the user needs to follow to extract the centreline(s).
    #   ================================================================================================================

    # Merge and save all dictionaries
    if save:
        dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour],
                                         names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour'],
                                         dir2save = directories[0])

    # Instructions for VMTK
    _, _ = fcBasics.code4vmtkCL(filename = filename, mesh_name = names_exp,
                           dir_cl = directories[3], printshow = True)

    toc = perf_counter()
    fcBasics.printTime(tic, toc, 'Create 3D Volumes')

#%% Init 
init = True