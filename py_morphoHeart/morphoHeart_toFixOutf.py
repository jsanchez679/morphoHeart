# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:56:04 2020

@author: mdp18js
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

from vtkplotter import *
from vtkplotter import embedWindow
embedWindow(False)
first = True

c="k"
font= 'CallingCode'
save = True
azimuth = 0

#%% Importing morphoHeart packages
import morphoHeart_funcBasics as fcBasics 
import morphoHeart_funcContours as fcCont
import morphoHeart_funcMeshes as fcMeshes

#%% Get main directories (check which ones are actually used)
_, _, dir_lsOngoing, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
df_dataset = fcBasics.exportDatasetCSV(dir_lsOngoing, dir_data2Analyse)
# Get file to process and directories 
folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]; dORv = filename[9:10]
#stage = df_file.loc[file_num,'Stage']
# directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse
dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')

# Import the metadata to know pixel size and distance between slices
xy_Scaling_um, z_Scaling_um = fcBasics.metadataExt(filename,dir_data2Analyse)
res = [xy_Scaling_um,xy_Scaling_um,z_Scaling_um]
# Import df_results
df_res = fcBasics.loadDF(filename = filename, file = 'ResultsDF', dir_results = dir_results)
file_num = df_res[df_res['Folder']==folder].index.values[0]
if dORv == 'D':
    azimuth = -90
# Initialise variables
dict_shapes = dict()

#%% Load data
# Import dictionaries
[dict_obj, myoc_int_npcl, endo_ext_npcl] = fcBasics.loadDicts(filename = filename, dicts_name = ['dict_obj','myoc_int_npcl','endo_ext_npcl'],
                                                                directories = [directories[0], directories[3], directories[3]])
[dict_planes, dict_pts, dict_kspl, dict_colour, _] = fcMeshes.splitDicts(dict_obj)
# Import meshes
[m_myoc, m_endo, m_cj, m_cjOut, m_cjIn] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc','endo','cj','cj_out','cj_in'],
                                                              extension = 'vtk', dir_stl = directories[2],
                                                              alpha = [0.5,0.5,0.5,1,1], dict_colour = dict_colour)
[m_myocInt, m_myocExt, m_endoInt, m_endoExt] = fcMeshes.openMeshes(filename = filename, meshes_names = ['myoc_int','myoc_ext','endo_int', 'endo_ext'],
                                                              extension = 'vtk', dir_stl = directories[2],
                                                              alpha = [1,1,1,1], dict_colour = dict_colour)

# Plot meshes
text = str(filename); txt = Text2D(text, c=c, font=font)
vp = Plotter(N=6, axes=7)
vp.show(m_myoc, txt, at=0, zoom=1)
vp.show(m_endo, at=1, zoom=1)
vp.show(m_cj, at=2, zoom=1)
vp.show(m_cjIn, at=3, zoom=1)
vp.show(m_cjOut, at=4, zoom=1)
vp.show(m_myoc, m_endo, m_cj, at=5, azimuth = azimuth, interactive=True)

#%% #TO DELETEEE!!!
# Myocardium
s3s2cut, stackShape = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1],
                                        end_name = ['ch0_cut', 'ch0_int', 'ch0_ext'])

s3_ch0_cut2, s3_ch0_int_cut2, s3_ch0_ext_cut2 = s3s2cut
# pl_outflow = dict_planes['pl2CutIm_outflow']
# pl_outfNormal = np.asarray(pl_outflow['pl_normal'])
# pl_outfCentre = np.asarray(pl_outflow['pl_centre'])

pl_inflow = dict_planes['pl2CutIm_inflow']
pl_infNormal = np.asarray(pl_inflow['pl_normal'])
pl_infCentre = np.asarray(pl_inflow['pl_centre'])

resolution = res
# d_cuts, normal_cuts, d_red, d_green = fcMeshes.getCuttingPlanesInfo([pl_infNormal, pl_outfNormal], [pl_infCentre, pl_outfCentre])
d_cuts, normal_cuts, d_red, d_green = fcMeshes.getCuttingPlanesInfo([pl_infNormal], [pl_infCentre])
# s3_ch0_cut, s3_ch0_int_cut, s3_ch0_ext_cut  = fcMeshes.cutInfAndOutfOpt(s3_ch0_cut2, s3_ch0_int_cut2, s3_ch0_ext_cut2, 
#                                                                         d_cuts, normal_cuts, d_red, d_green, resolution, '(Myoc)')
# Cut myocardial s3_all, s3_int, s3_ext
s3_ch0_cut, s3_ch0_int_cut, s3_ch0_ext_cut = fcMeshes.cutInfOrOutfOpt (s3_ch0_cut2, s3_ch0_int_cut2, s3_ch0_ext_cut2, 
                                                                      d_cut = d_cuts[0], normal_cut = normal_cuts[0], 
                                                                      d_red = d_red[0], d_green = d_green[0], resolution = resolution, 
                                                                      option = 'inflow', mesh_name = '(Myoc)')

#Create new mesh myocardial s3_all
myoc_cut = fcMeshes.getCutMesh(filename = filename, s3_cut = s3_ch0_cut, resolution = resolution, 
                                        mesh_original = m_myoc, layer = 'Myoc', plotshow = True)
myoc_cut.color('darkcyan').alpha(0.5)
m_myoc = myoc_cut
#Create new mesh myocardial s3_int
myoc_cut_int = fcMeshes.getCutMesh(filename = filename, s3_cut = s3_ch0_int_cut, resolution = resolution, 
                                mesh_original = '', layer = 'Int.Myoc', plotshow = False)
#Create new mesh myocardial s3_ext
myoc_cut_ext = fcMeshes.getCutMesh(filename = filename, s3_cut = s3_ch0_ext_cut, resolution = resolution, 
                                mesh_original = '', layer = 'Ext.Myoc', plotshow = False)

vp = Plotter(N=2, axes = 7)
vp.show(myoc_cut_ext, at = 0)
vp.show(myoc_cut_int, at = 1, interactive = True)

if save:
    fcMeshes.save_s3s(filename = filename, s3_all = s3_ch0_cut, s3_int = s3_ch0_int_cut, s3_ext = s3_ch0_ext_cut, 
            dir_txtNnpy = directories[1], layer = 'ch0_cut')
    
    dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [myoc_cut, myoc_cut_int, myoc_cut_ext],
                        names = ['myoc', 'myoc_int', 'myoc_ext'], dict_colour = dict_colour,
                        dir_stl = directories[2], extension = 'vtk') 
    
#%%
# Endocardium
s3s2cut_endo, stackShape = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1],
                                        end_name = ['ch1_cut', 'ch1_int', 'ch1_ext'])
s3_ch1_cut2, s3_ch1_int_cut2, s3_ch1_ext_cut2 = s3s2cut_endo

# s3_ch1_cut, s3_ch1_int_cut, s3_ch1_ext_cut  = fcMeshes.cutInfAndOutfOpt(s3_ch1_cut2, s3_ch1_int_cut2, s3_ch1_ext_cut2, 
#                                                                         d_cuts, normal_cuts, d_red, d_green, resolution, '(Endo)')
s3_ch1_cut, s3_ch1_int_cut, s3_ch1_ext_cut = fcMeshes.cutInfOrOutfOpt (s3_ch1_cut2, s3_ch1_int_cut2, s3_ch1_ext_cut2, 
                                                                      d_cut = d_cuts[0], normal_cut = normal_cuts[0], 
                                                                      d_red = d_red[0], d_green = d_green[0], resolution = resolution, 
                                                                      option = 'inflow', mesh_name = '(Endo)')

#Create new mesh endocardial s3_all
endo_cut = fcMeshes.getCutMesh(filename = filename, s3_cut = s3_ch1_cut, resolution = resolution, 
                                        mesh_original = m_endo, layer = 'Endo', plotshow = False)
endo_cut.color('darkmagenta').alpha(0.5)
m_endo = endo_cut
#Create new mesh endocardial s3_int
endo_cut_int = fcMeshes.getCutMesh(filename = filename, s3_cut = s3_ch1_int_cut, resolution = resolution, 
                                mesh_original = '', layer = 'Int.Endo', plotshow = False)
#Create new mesh endocardial s3_ext
endo_cut_ext = fcMeshes.getCutMesh(filename = filename, s3_cut = s3_ch1_ext_cut, resolution = resolution, 
                                mesh_original = '', layer = 'Ext.Endo', plotshow = False)

vp = Plotter(N=2, axes = 7)
vp.show(endo_cut_int, at = 0)
vp.show(endo_cut_ext, at = 1, interactive = True)

if save:
    fcMeshes.save_s3s(filename = filename, s3_all = s3_ch1_cut, s3_int = s3_ch1_int_cut, s3_ext = s3_ch1_ext_cut, 
            dir_txtNnpy = directories[1], layer = 'ch1_cut')
    
    dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [endo_cut, endo_cut_int, endo_cut_ext],
                        names = ['endo', 'endo_int', 'endo_ext'], dict_colour = dict_colour,
                        dir_stl = directories[2], extension = 'vtk') 

#%% Cardiac jelly
# [s3_ch0_int_cut, s3_ch1_ext_cut], stackShape = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1],
#                                         end_name = ['ch0_cut_int', 'ch1_cut_ext'])
# # Get cardiac jelly per slice
# s3_cjIn, s3_cj = fcCont.ch_clean(mask_s3 = s3_ch1_ext_cut, toClean_s3 = s3_ch0_int_cut, option = "cardiacjelly")
# fcCont.ch_clean_plt(mask_s3 = s3_ch0_int_cut, toClean_s3 = s3_ch1_ext_cut, toRemove_s3 = s3_cjIn, 
#                           cleaned_s3 = s3_cj, plotshow = True, im_every = 50, option = "cardiacjelly")
# # Create meshes - CJ, CJ Internal surface and CJ External surface
# cj_all, cj_in, cj_out = fcMeshes.createAll3LayerMeshes(filename = filename, s3_all = s3_cj, s3_in = s3_cjIn,
#                                     s3_out = s3_ch0_int_cut, resolution = res, layer = 'CJ')
# # Save Heart layers (cj) and s3s
# if save:
#     fcCont.save_s3s(filename = filename, s3_all = s3_cj, s3_int = s3_cjIn, s3_ext = s3_ch0_int_cut, 
#                             dir_txtNnpy = directories[1], layer = 'cj')
#     dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [cj_all, cj_in, cj_out],
#                         names = ['cj', 'cj_in', 'cj_out'], dict_colour = dict_colour,
#                         dir_stl = directories[2], extension = 'vtk') 
    
#%%
s3s2cut, stackShape = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1],
                                        end_name = ['ch0_cut_ext', 'ch1_cut_int'])

names2cut = ['Ext.Myoc', 'Int.Endo']
# Create empty lists to save atrium and ventricles
atr_meshes = []
atr_color = ['darkblue', 'indigo']
vent_meshes = []
vent_color = ['indigo', 'darkblue']

pl_Chamber = dict_planes['pl2CutIm_Chamber']
pl_imCh_normal = np.asarray(pl_Chamber['pl_normal'])
pl_imCh_centre = np.asarray(pl_Chamber['pl_centre'])

from itertools import count

for n, s3, name in zip(count(), s3s2cut, names2cut):
    # Mask s3s vent and atrium
    print('- Cutting s3 (', name,')')
    s3_vent, s3_atr = fcMeshes.maskChamberS3s(s3_mask = s3, pl_normal = pl_imCh_normal, pl_centre = pl_imCh_centre, resolution = resolution)
    # Create chamber meshes 
    atr = fcMeshes.createLayerMesh(filename = filename, s3 = s3_atr, resolution = resolution, layer = name, name = name+'_Atr',
                                    colour = atr_color[n], alpha = 0.01, plotshow = False)
    atr_meshes.append(atr)
    vent = fcMeshes.createLayerMesh(filename = filename, s3 = s3_vent, resolution = resolution, layer = name, name = name+'_Vent',
                                    colour = vent_color[n], alpha = 0.01, plotshow = False)
    vent_meshes.append(vent)

print('- All layers have been cut!')
m_atrExtMyo, m_atrIntEnd = atr_meshes
m_ventExtMyo, m_ventIntEnd = vent_meshes

vp = Plotter(N=2, axes = 10)
vp.show(atr_meshes[0], vent_meshes[0], at = 0, zoom = 2.5)
vp.show(atr_meshes[1], vent_meshes[1], at = 1, zoom = 2.5, azimuth = azimuth, interactive=True)
    
# Add chamber volume info for each layer and save df
meshes = [m_atrExtMyo, m_ventExtMyo, m_atrIntEnd, m_ventIntEnd]
    
names = ['Atr.ExtMyoc', 'Vent.ExtMyoc', 
             'Atr.IntEndo', 'Vent.IntEndo']
for n, name, mesh in zip(count(), names, meshes):
    volume = mesh.volume()
    df_res.loc[file_num,'Vol_'+name] = volume
    
fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = dir_results)
fcBasics.saveFilledDF(filename = filename, df_res = df_res, dir2save = os.path.join(dir_data2Analyse, 'R_All'))