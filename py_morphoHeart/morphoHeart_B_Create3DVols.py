# -*- coding: utf-8 -*-
"""
morphoHeart - B. CREATE 3D VOLUMES AND MESHES TO EXTRACT CENTRELINE
@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
import platform

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

from vtkplotter import *
from vtkplotter import embedWindow
embedWindow(False)
save = True

#%% Importing morphoHeart packages
import morphoHeart_funcBasics as fcBasics
import morphoHeart_funcContours as fcCont
import morphoHeart_funcMeshes as fcMeshes

#%% Get main directories (check which ones are actually used)
_, _, dir_lsOngoing, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
df_dataset = fcBasics.exportDatasetCSV(dir_lsOngoing, dir_data2Analyse)
# Get file to process and directories 
folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]
# directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse
dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')
# Import the metadata to know pixel size and distance between slices
xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]
# Initialise variables
plotshow = False
dict_planes = dict(); dict_pts = dict(); dict_kspl = dict(); dict_colour = dict()

#%% Load stacks
s3s, stackShape = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1], 
                                    end_name = ['ch0_int','ch0_ext','ch0_all','ch1_int','ch1_ext','ch1_all' ])
s3_ch0_int, s3_ch0_ext, s3_ch0, s3_ch1_int, s3_ch1_ext, s3_ch1 = s3s

#%% Clean endocardium
# Plot s3s - channel side by side
fcCont.plt_s3(start_slc = 0, end_slc = s3_ch0.shape[2]-1, im_every = 10, 
                  s3_int = s3_ch0, s3_ext = s3_ch1, plotshow = plotshow, option = "both ch")
# Create Ext Myocardial mesh - Ch0
myoc_ext = fcMeshes.createExtLayerMesh(filename = filename, s3_ext = s3_ch0_ext, resolution = res, layer = 'Myoc', info = 'Original', plotshow = plotshow)
# Create Ext Endocardial mesh - Ch1
endo_o = fcMeshes.createExtLayerMesh(filename = filename, s3_ext = s3_ch1_ext, resolution = res, layer = 'Endo', info = 'Original', plotshow = plotshow)
# Clean Endocardium
s3_endo_rem, s3_ch1_cl = fcCont.ch_clean(s3_ch0, s3_ch1, option = '')
fcCont.ch_clean_plt(mask_s3 = s3_ch0, toClean_s3 = s3_ch1, toRemove_s3 = s3_endo_rem, 
                          cleaned_s3 = s3_ch1_cl, plotshow = plotshow, im_every = 15, option = "clean")
# Clean Ext Endocardium
s3_ch1_rem, s3_ch1_ext_cl = fcCont.ch_clean(mask_s3 = s3_ch0, toClean_s3 = s3_ch1_ext, option = "clean")
# Re-Create Ext Endocardial mesh - Ch1
endo_ext = fcMeshes.createExtLayerMesh(filename = filename, s3_ext = s3_ch1_ext_cl, resolution = res, layer = 'Endo', info = 'Cleaned', plotshow = plotshow)

# Plot Ext Myocardium and Clean Ext Endocardium
myoc_ext.alpha(0.01); endo_o.color('mediumorchid').legend('Orig.Ext.Endo')
vp = Plotter(N=4, axes=4)
vp.show(endo_o, at=0, zoom=1)
vp.show(myoc_ext, endo_o, at=1, zoom=1)
vp.show(endo_ext, at=2, zoom=1)
vp.show(myoc_ext, endo_ext, at=3, zoom=1, interactive=True)

#%% Cut inflow and outflow regions of endocardium and myocardium and save all final s3s
s3s_cut, meshes_cut, dict_planes = fcMeshes.selectCutS3sOptMx(filename = filename, 
                                    s3s2cut = [s3_ch0_int, s3_ch0_ext, s3_ch0, s3_ch1_int, s3_ch1_ext_cl, s3_ch1_cl], 
                                    m_endo = endo_ext, m_myoc = myoc_ext, 
                                    dict_planes = dict_planes, resolution = res, 
                                    dir_txtNnpy = directories[1], save = save)

myoc_cut, myoc_cut_int, myoc_cut_ext, endo_cut, endo_cut_int, endo_cut_ext = meshes_cut

# vp = Plotter(N=6, axes=7)
# vp.show(myoc_cut, at=0, zoom=1)
# vp.show(myoc_cut_int, at=1, zoom=1)
# vp.show(myoc_cut_ext, at=2,  zoom=1)
# vp.show(endo_cut, at=3,  zoom=1)
# vp.show(endo_cut_int, at=4)
# vp.show(endo_cut_ext, at=5, zoom=1, interactive=True)

s3_ch0_int_cut, _, _, _, s3_ch1_ext_cut, _ = s3s_cut
if save:
    dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [myoc_cut, myoc_cut_int, myoc_cut_ext, endo_cut, endo_cut_int, endo_cut_ext],
                        names = ['myoc', 'myoc_int', 'myoc_ext', 'endo', 'endo_int', 'endo_ext'], dict_colour = dict_colour,
                        dir_stl = directories[2], extension = 'vtk') 

#%% Extract cardiac jelly from contours
# Plot s3s - channel side by side
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
vp = Plotter(N=6, axes=7)
vp.show(myoc_cut, at=0, zoom=1)
vp.show(endo_cut, at=1, zoom=1)
vp.show(myoc_cut_int.color('turquoise'), at=2,  zoom=1)
vp.show(cj_all.clone().alpha(0.05), at=3,  zoom=1)
vp.show(myoc_cut.clone().alpha(0.1), endo_cut, cj_all, at=4)
vp.show(endo_cut_ext.color('orchid'), at=5, zoom=1, interactive=True)

#%% Save all meshes and plane dict
# Save SurfArea and Int/Ext Volume of meshes to dataframe and save df
df_res = fcMeshes.addSurfArea2df(df_res = df_file,  file_num = file_num,
                        meshes = [myoc_cut, myoc_cut_int, myoc_cut_ext, endo_cut, endo_cut_int, endo_cut_ext, cj_all, cj_in, cj_out])

df_res = fcMeshes.addIntExtVol2df(df_res = df_res, file_num = file_num, meshes = [myoc_cut_int, myoc_cut_ext, endo_cut_int, endo_cut_ext])
if save:
    fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)

#%% Create, cut and export meshes for centreline
# Create meshes for centreline
meshes4cl = fcMeshes.createMeshes4CL(filename = filename, meshes = [myoc_cut_int, endo_cut_ext], 
                            names = ['Int.Myoc(Cut)', 'Ext.Endo(Cut)'], mesh_colors = ['springgreen', 'violet'], 
                            plotshow = True)

# Cut inflow and outflow tract of meshes for centreline
meshes4clf, dicts = fcMeshes.cutMeshes4CL(filename = filename, meshes = meshes4cl, names = ['Int.Myoc(Cut)', 'Ext.Endo(Cut)'],
                                 cuts = ['inflow', 'outflow'], cut_direction = [True, False], 
                                 mark_colors = ['deepskyblue', 'tomato'], mesh_colors = ['springgreen', 'violet'], 
                                 dicts = [dict_planes, dict_pts, dict_kspl], plotshow = True)
myoc_int_CL, endo_ext_CL = meshes4clf

#%% Save meshes 4 centreline
if save:
    dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [myoc_int_CL, endo_ext_CL], names =['myoc_int_cut4cl','endo_ext_cut4cl'], 
                                  dict_colour = dict_colour, dir_stl = directories[3], extension='stl')

#%% Plot and save all meshes 
vp = Plotter(N=6, axes=7)
vp.show(myoc_cut, at=0, zoom=1)
vp.show(endo_cut, at=1, zoom=1)
vp.show(myoc_int_CL, at=2,  zoom=1)
vp.show(cj_all, at=3,  zoom=1)
vp.show(myoc_cut, endo_cut, cj_all, at=4)
vp.show(endo_ext_CL, at=5, zoom=1, interactive=True)

#%% Merge and save all dictionaries
dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour], 
                                     names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour'], dir2save = directories[0])

#%% Instructions for VMTK
fcBasics.code4vmtkCL(filename = filename, mesh_name = ['myoc_int','endo_ext'],
                       dir_cl = directories[3], printshow = True)