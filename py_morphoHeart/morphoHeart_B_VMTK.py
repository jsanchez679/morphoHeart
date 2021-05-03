# -*- coding: utf-8 -*-
"""
morphoHeart - B. CALCULATE CENTRELINES USING VMTK AND EXPORT RESULTS
Welcome to the third code of morphoHeart!
If you are running this code you must have already created one or two meshes to extract the centreline from and have
already used MeshLab to reconstruct them and cut both its inflow and outflow tracts to make sure they have blunt ends.
If so, the objective of this script is to extract centreline(s) from this (these) mesh(es) and save it(them) as (a) 
dictionary(ies) with numpy arrays, so that we can use them in future scripts. 

Happy centreline extraction!

@author: Juliana Sanchez-Posada
Version: 13th April, 2021
"""
#%% Importing python packages
import os
import vmtk
from vmtk import pypes, vmtkscripts

# Verify working dir
def setWorkingDir (root_path, init = False):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    init = not bool(int(input('Do you want to execute the script all at once or run it by cells? \n\t[0]: all at once (recommended if you are already familiar with the script)\n\t[1]: by cells (recommended if you are NOT yet familiar with the script). >>>: ')))
    print("Current working directory: {0}".format(os.getcwd()))
    if not init: 
        print('\nIMPORTANT NOTES:\n- Remember to start running from cell %Start B_VMTK.\n- NEVER run as an individual cell the cell called %Importing python packages')
    return root_path, init

root_path, init = setWorkingDir(os.getcwd())

#%% Start B_VMTK
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics

    #%% SELECT FILE AND GET METADATA 
    #   This section allows the user to select file to process, get its main directories and metadata, 
    #   and define some properties
    #   ================================================================================================================
    
    # Get main directories
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse)
    # Get file to process and directories
    folder, df_file, file_num, blind = fcBasics.selectFile(df_dataset); filename = folder[0:-3]
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')

    #%% GET MESHLAB MESHES AND EXTRACT CENTRELINES
    #   This section will find the mesh(es) that was(were) processed in Meshlab, and use it(them) to get the heart
    #   centreline. When running each mesh the code will show the user a window with a 3D interactive plot of the 
    #   selected mesh and a numbered cap in each blunt cut (inflow and outflow tract). (Make sure both blunt cuts
    #   are numbered and appear as it they had been closed/capped. If not, it is recommended to cancel the execution of 
    #   the code -clicking on the red square on the console-, go back to Meshlab, cut a little bit more that heart 
    #   section, resave the mesh and re-run this code.)
    #   When pressing 'q', the user will be asked to introduce the cap number (profile IDs) corresponding to each, 
    #   the inflow and outflow tracts of the tube (in this case the heart), after which the centrelines will start to be
    #   calculated. When this process is finished, the pop-up window will load again, and when pressing 'q' the user 
    #   will be able to see the resulting centrelines. At the end of the code, the information related to the resulting 
    #   centrelines will be saved as dictionaries and will be imported in the next script.
    #   ================================================================================================================
    
    mesh_name = []
    while len(mesh_name) == 0: 
        vmtktxt, mesh_name = fcBasics.code4vmtkCL(filename = filename, mesh_name = ['myoc_int','endo_ext'],
                                                    dir_cl = directories[3], printshow = False)
        if len(mesh_name) == 0:
            _ = fcBasics.ask4input('No meshes are recognised in the path: "'+ str('//'.join(os.path.normpath(directories[3]).split(os.sep)[-4:])) + '". \n\t Make sure you have named your cleaned meshes correctly after running the processing in Meshlab and press -Enter- when ready.', str)
    
    
    for n, mesh in enumerate(mesh_name):
        #Create centrelines
        myArguments = vmtktxt[n]
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

        fcBasics.saveDict(filename = filename, dict2save = centrelines_dict, name = mesh+'_npcl', dir2save = directories[3])

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

#%%
init = True