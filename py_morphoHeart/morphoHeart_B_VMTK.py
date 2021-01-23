# -*- coding: utf-8 -*-
"""
morphoHeart - B. VTP CENTRELINE FORMAT TO NUMPY ARRAY
Script to extract the information from the centreline and save it as a
dictionary with numpy arrays

Only to be run on activated environment vmtk_env (Python 3.6)

@author: Juliana Sanchez-Posada
"""
#%% Importing python packages
import os
import platform
import numpy as np
import json
from vmtk import vmtkscripts

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

#%% Start B_VMTK
if init:
    # Class to save dictionary
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

    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    # alert,ask4input,getMainDirectories,exportDatasetCSV,selectFile,createDirectories2Save

    #%% Get main directories (check which ones are actually used)
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse)
    # Get file to process and directories
    folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')

    #%% Get vtp files and export them as np.arrays
    cl_dir = fcBasics.code4vmtkCL(filename = filename, mesh_name = ['myoc_int','endo_ext'],
                           dir_cl = directories[3], printshow = False)

    mesh_name = ['myoc_int','endo_ext']
    for n, mesh in enumerate(mesh_name):
        # Read the vtp file as a dictionary
        centerlineReader = vmtkscripts.vmtkSurfaceReader()
        centerlineReader.InputFileName = cl_dir[n]
        centerlineReader.Execute()
        clNumpyAdaptor = vmtkscripts.vmtkCenterlinesToNumpy()
        clNumpyAdaptor.Centerlines = centerlineReader.Surface
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

        # Saving a small version of heartLayer dictionary
        jsonCentr_name = filename+"_"+mesh+"_npcl.json"
        json2save_dir = os.path.join(directories[3],jsonCentr_name)

        with open(json2save_dir, "w") as write_file:
            json.dump(centrelines_dict, write_file, cls=NumpyArrayEncoder)

        fcBasics.alert('countdown',1)
        print("\n>>>>>>> numpyCenterlines saved correctly for "+mesh+"!\n\n")

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
