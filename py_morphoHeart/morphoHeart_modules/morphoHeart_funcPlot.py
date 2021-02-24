# -*- coding: utf-8 -*-
"""
morphoHeart_funcPlot
Functions included:
    -

Version: Nov, 2020
@author: Juliana Sanchez-Posada

"""

#%% Importing python packages
import os
import numpy as np
from itertools import count
import math
from progress.bar import Bar
suffix = '%(index)d/%(max)d - %(elapsed)ds'

from vedo import *
from vedo import embedWindow#, Plotter, Text2D, settings, load
embedWindow(False)

c="k"
font= 'CallingCode'

#%% Importing morphoHeart packages
from .morphoHeart_funcBasics import alert, ask4input, loadDF, loadDicts
from .morphoHeart_funcMeshes import splitDicts

#%% func - createDictPaths
def createDictPaths(names, df_files, dir_data2Analyse):
    """


    Parameters
    ----------
    names : TYPE
        DESCRIPTION.
    df_files : TYPE
        DESCRIPTION.
    dir_data2Analyse : TYPE
        DESCRIPTION.

    Returns
    -------
    dict_paths : TYPE
        DESCRIPTION.

    """

    names_meshes = ['myoc','myoc_int','myoc_ext','myoc_atr','myoc_vent',
                    'endo','endo_int','endo_ext','endo_atr','endo_vent',
                    'cj','cj_in','cj_out','cj_atr','cj_vent']
    names_thickness = ['myoc_extBall','myoc_intBall','myoc_thickness','endo_thickness', 'cj_thickness']

    dict_paths = dict()
    for file in df_files['Folder'].tolist():
        filename = file[2:]
        dict_paths[file] = file_dict = dict()
        file_dict['path_results'] =  os.path.join(dir_data2Analyse, 'R_'+filename,'Results_'+filename)
        dir_meshes = os.path.join(dir_data2Analyse, 'R_'+filename,'Results_'+filename, 'meshes')
        dir_txtNnpy = os.path.join(dir_data2Analyse, 'R_'+filename,'Results_'+filename, 'txt_npy')
        dir_dict = os.path.join(dir_data2Analyse, 'R_'+filename,'Results_'+filename, 'dictionaries')
        dir_cl = os.path.join(dir_data2Analyse, 'R_'+filename,'Results_'+filename, 'centreline')
        file_dict['path_dict'] =  os.path.join(dir_dict)
        file_dict['path_stl'] =  os.path.join(dir_meshes)
        for name in names:
            if name in names_meshes:
                file_dict['path_'+name] = os.path.join(dir_meshes,filename+'_'+name+'.vtk')

            elif name in names_thickness:
                file_dict['path_'+name] = os.path.join(dir_meshes,filename+'_'+name+'.vtk')
                file_dict['thickness_'+name] = os.path.join(dir_txtNnpy,filename+'_'+name+'.npy')

            else:
                file_dict['path_myocInt_CL'] = os.path.join(dir_cl,filename+'_myoc_int_npcl.json')
                file_dict['path_endoExt_CL'] = os.path.join(dir_cl,filename+'_endo_ext_npcl.json')

    return dict_paths

    #Later to create all
#     [dict_obj, myoc_int_npcl, endo_ext_npcl] = fcBasics.loadDicts(filename = filename, dicts_name = ['dict_obj','myoc_int_npcl','endo_ext_npcl'],
#                                                                 directories = [directories[0], directories[3], directories[3]])
# [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes] = fcMeshes.splitDicts(dict_obj)


#%% func - plotTitle
def plotTitle(file, df_files, additional_title = ''):
    """


    Parameters
    ----------
    file : TYPE
        DESCRIPTION.
    df_files : TYPE
        DESCRIPTION.
    additional_title : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    txt : TYPE
        DESCRIPTION.

    """

    file_num = df_files[df_files['Folder']==file].index.values[0]
    text_file = file[2:]
    text_strain = "\n >> Strain: " + df_files.loc[file_num,'Strain']
    text_stage = "\n >> Stage: " + df_files.loc[file_num,'Stage']+" hpf"
    text_genotype = "\n >> Genotype " + df_files.loc[file_num,'Gene_A']+': '+df_files.loc[file_num,'Genotype_A']

    if df_files.loc[file_num,'Gene_B'] != '-':
        text_genotype = text_genotype + ' / Genotype '+ df_files.loc[file_num,'Gene_B']+': '+df_files.loc[file_num,'Genotype_B']

    text = text_file+text_strain+text_stage+text_genotype

    if additional_title != '':
        text = additional_title + "\n "

    txt = Text2D(text, c=c, font= font)

    return txt

#%% func - selectObjects
def selectObjects ():
    """


    Returns
    -------
    objs2loop : TYPE
        DESCRIPTION.

    """

    objs = ['myoc','myoc_int','myoc_ext','myoc_atr','myoc_vent',
              'endo','endo_int','endo_ext','endo_atr','endo_vent',
              'cj','cj_in','cj_out','cj_atr','cj_vent',
              'myoc_extBall','myoc_intBall','myoc_thickness','endo_thickness', 'cj_thickness']
    objsN = []
    print('\nObjects:')
    for c, obj in zip(count(), objs):
        if obj != '':
            txt = str(c)+'. '+ obj
        else:
            txt = obj
        objsN.append(txt)
    list_columns(objsN)

    input_obj = input('Select the objects you would like to plot: ')

    if input_obj == 'All':
        obj_num = list(range(0,len(objs),1))

    else:
        obj_num = []
        comma_split = input_obj.split(',')

        for string in comma_split:
            if '-' in string:
                minus_split = string.split('-')
                #print(minus_split)
                for n in list(range(int(minus_split[0]),int(minus_split[1])+1,1)):
                    #print(n)
                    obj_num.append(n)
            else:
                obj_num.append(int(string))

    objs2loop = []
    for i, num in enumerate(obj_num):
        objs2loop.append(objs[num])

    return objs2loop


#%% func - list_columns
def list_columns(obj, cols=4, columnwise=True, gap=4):
    """
    Print the given list in evenly-spaced columns.

    Parameters
    ----------
    obj : list
        The list to be printed.
    cols : int
        The number of columns in which the list should be printed.
    columnwise : bool, default=True
        If True, the items in the list will be printed column-wise.
        If False the items in the list will be printed row-wise.
    gap : int
        The number of spaces that should separate the longest column
        item/s from the next column. This is the effective spacing
        between columns based on the maximum len() of the list items.
    """
    # #https://stackoverflow.com/questions/1524126/how-to-print-a-list-more-nicely
    # for a,b,c,d,e in zip(objsN[::5],objsN[1::5],objsN[2::5],objsN[3::5],objsN[4::5]):
    #     print ('{:<25}{:<25}{:<25}{:<25}{:<}'.format(a,b,c,d,e))

    sobj = [str(item) for item in obj]
    if cols > len(sobj): cols = len(sobj)
    max_len = max([len(item) for item in sobj])
    if columnwise: cols = int(math.ceil(float(len(sobj)) / float(cols)))
    plist = [sobj[i: i+cols] for i in range(0, len(sobj), cols)]
    if columnwise:
        if not len(plist[-1]) == cols:
            plist[-1].extend(['']*(len(sobj) - len(plist[-1])))
        plist = zip(*plist)
    printer = '\n'.join([''.join([c.ljust(max_len + gap) for c in p]) for p in plist])

    print(printer)

#%% func - loadMultMeshes # CONTINUE!!!
def loadMultMeshes(names, dict_paths):
    """


    Parameters
    ----------
    names : TYPE
        DESCRIPTION.
    dict_paths : TYPE
        DESCRIPTION.

    Returns
    -------
    list_meshes : TYPE
        DESCRIPTION.

    """

    names_meshes = ['myoc','myoc_int','myoc_ext','myoc_atr','myoc_vent',
                    'endo','endo_int','endo_ext','endo_atr','endo_vent',
                    'cj','cj_in','cj_out','cj_atr','cj_vent']

    names_thickness = ['myoc_extBall','myoc_intBall','myoc_thickness','endo_thickness', 'cj_thickness']

    # - Names and legends
    names_all = ['myoc','myoc_ext','myoc_int','myoc_atr','myoc_vent',
                  'endo','endo_ext','endo_int','endo_atr','endo_vent',
                  'cj','cj_out','cj_in','cj_atr','cj_vent',
                  'myoc_intBall', 'myoc_extBall', 'myoc_thickness','endo_thickness','cj_thickness']

    legend_all = ['Myocardium','Ext.Myoc', 'Int.Myoc','Atrium(Myoc)','Ventricle(Myoc)',
                  'Endocardium', 'Ext.Endo', 'Int.Endo','Atrium(Endo)','Ventricle(Endo)',
                  'CardiacJelly','Ext.CJ','Int.CJ','Atrium(CJ)','Ventricle(CJ)',
                  'Int.Myoc Ball.','Ext.Myoc Ball.','Myoc.Thickness','Endo.Thickness','CJ.Thickness']

    bar_names_all = ['Int.Myoc\nBalloning\n[um]','Ext.Myoc\nBalloning\n[um]','Myoc.Thickness\n[um]','Endo.Thickness\n[um]','CJ.Thickness\n[um]']
    alpha_all = [1,1,0.1,1,1]

    list_meshes = []
    alphas = []; mins = []; maxs = []
    m = 0

    bar = Bar('- Loading meshes...', max = len(dict_paths.keys())*len(names), suffix = suffix, check_tty=False, hide_cursor=False)
    for n, file in enumerate(dict_paths.keys()):
        # print('File: ',file)
        # Image taken from side (D) or Ventral
        dORv = file[11:12]
        if dORv == 'V':
            rotY = 0
        elif dORv == 'D':
            rotY = 90

        df_res = loadDF(filename = file[2:], file = 'ResultsDF', dir_results = dict_paths['R_'+file[2:]]['path_results'])
        file_num = df_res[df_res['Folder']==file[2:]+'_2A'].index.values[0]
        rotX = df_res.loc[file_num,'ang_Heart']
        # myoc_title = file[2:]+"_myoc.vtk"
        # myoc = load(os.path.join(dict_paths[file]['path_stl'], myoc_title)).rotateY(rotY).rotateX(rotX)
        # xc,yc,zc = myoc.centerOfMass().tolist()
        # print('xyzc:', xc, yc, zc)

        # print('ang_Heart:', rotX)
        [dict_obj] = loadDicts(filename = file[2:], dicts_name = ['dict_obj'], directories = [dict_paths['R_'+file[2:]]['path_dict']], print_txt = False)
        [_, _, _, dict_colour, _] = splitDicts(dict_obj)

        if n == 1:
            alphas = alphas*len(dict_paths.keys())
            mins = mins*len(dict_paths.keys())
            maxs = maxs*len(dict_paths.keys())

        for name in names:
            index = names_all.index(name)
            path = dict_paths[file]['path_'+name]
            mesh_out = load(path)
            # print('name: ',name, index)

            if name in names_meshes:
                if n == 0:
                    q_alpha = ask4input('Alpha value for -'+ name + '-: ', float)
                    alphas.append(q_alpha)
                    # print ('\n')
                    mins.append(0)
                    maxs.append(20)

                mesh_colour = dict_colour[name]['colour']
                mesh_out.alpha(alphas[m]).legend(legend_all[index]).wireframe().color(mesh_colour).rotateY(rotY).rotateX(rotX)
                # x0,y0,z0 = mesh_out.centerOfMass().tolist()
                # mesh_out = mesh_out.x(x0-xc).y(y0-yc).z(z0-zc)
                # print(x0, y0, z0)
                #print(mesh_out.centreOfMass())

            elif name in names_thickness:
                if n == 0:
                    q_min = ask4input('Minimum value for -'+ name + '-: ',float)
                    mins.append(q_min)
                    q_max = ask4input('Maximum value for -'+ name + '-: ',float)
                    maxs.append(q_max)
                    # print ('\n')
                    alphas.append(1)

                # Load colour array
                path_thickness = dict_paths[file]['thickness_'+name]
                sp_colour = np.load(path_thickness)
                mesh_out.pointColors(sp_colour, cmap="jet", vmin = mins[m], vmax =maxs[m])
                mesh_out.addScalarBar(title=bar_names_all[index-15])
                mesh_out.alpha(alpha_all[index-15]).legend(legend_all[index]).wireframe().rotateY(rotY).rotateX(rotX)
                # x0,y0,z0 = mesh_out.centerOfMass().tolist()
                # mesh_out = mesh_out.x(x0-xc).y(y0-yc).z(z0-zc)
                # print(x0, y0, z0)
                #print(mesh_out.centreOfMass())
                mesh_out.mapper().SetScalarRange(mins[m], maxs[m])

            list_meshes.append(mesh_out)
            bar.next()
            m += 1
                # print(list_meshes)
    bar.finish()

    return list_meshes

#%% func - yieldMultMeshes
def yieldMultMeshes(names, dict_paths):
    """


    Parameters
    ----------
    names : TYPE
        DESCRIPTION.
    dict_paths : TYPE
        DESCRIPTION.

    Yields
    ------
    list_meshes : TYPE
        DESCRIPTION.

    """
    #https://www.programiz.com/python-programming/generator

    names_meshes = ['myoc','myoc_int','myoc_ext','myoc_atr','myoc_vent',
                    'endo','endo_int','endo_ext','endo_atr','endo_vent',
                    'cj','cj_in','cj_out','cj_atr','cj_vent']
    names_thickness = ['myoc_extBall','myoc_intBall','myoc_thickness','endo_thickness', 'cj_thickness']

    # - Names and legends
    legend_all = ['Myocardium','Ext.Myoc', 'Int.Myoc','Atrium(Myoc)','Ventricle(Myoc)',
                  'Endocardium', 'Ext.Endo', 'Int.Endo','Atrium(Endo)','Ventricle(Endo)',
                  'CardiacJelly','Ext.CJ','Int.CJ','Atrium(CJ)','Ventricle(CJ)']
    names_all = ['myoc','myoc_ext','myoc_int','myoc_atr','myoc_vent',
                  'endo','endo_ext','endo_int','endo_atr','endo_vent',
                  'cj','cj_out','cj_in','cj_atr','cj_vent']

    for n, file in enumerate(dict_paths.keys()):
        list_meshes = []
        #Check best way to get this df either going into each folder or getting the big collated df?
        df_res = loadDF(filename = file[2:], file = 'ResultsDF', dir_results = dict_paths['R_'+file[2:]]['path_results'])
        [dict_obj] = loadDicts(filename = file[2:], dicts_name = ['dict_obj'], directories = [dict_paths['R_'+file[2:]]['path_dict']])
        [_, _, _, dict_colour, _] = splitDicts(dict_obj)

        for name in names:
            if name in names_meshes:
                index = names_all.index(name)
                path = dict_paths[file]['path_'+name]
                mesh_out = load(path)
                mesh_colour = dict_colour[name]['colour']
                mesh_out.alpha(1).legend(legend_all[index]).wireframe().color(mesh_colour)

                # vp = Plotter(N=1, axes=10)
                # vp.show(mesh_out, at=0, interactive=True)

                list_meshes.append(mesh_out)

    # return list_meshes
        print(file)
        print(list_meshes)

        yield list_meshes

#%% func - getVarsANDLabels
def getVarsANDLabels (variables, labels):
    """
    

    Parameters
    ----------
    variables : TYPE
        DESCRIPTION.
    labels : TYPE
        DESCRIPTION.

    Returns
    -------
    vars2loop : TYPE
        DESCRIPTION.
    labels2loop : TYPE
        DESCRIPTION.

    """

    print('\nVariables:')
    for c, value in enumerate(variables, 1):
        print(c-1, value)
    input_var = input('Select the variables you would like to process: ')

    if input_var == 'All':
        var_num = list(range(0,len(variables),1))

    else:
        var_num = []
        comma_split = input_var.split(',')

        for string in comma_split:
            if '-' in string:
                minus_split = string.split('-')
                #print(minus_split)
                for n in list(range(int(minus_split[0]),int(minus_split[1])+1,1)):
                    #print(n)
                    var_num.append(n)
            else:
                var_num.append(int(string))

    vars2loop = []
    labels2loop = []
    for i, num in enumerate(var_num):
        vars2loop.append(variables[num])
        labels2loop.append(labels[num])

    return vars2loop, labels2loop


#%% - ALERT WHEN IMPORTED
print ("IMPORTED: morphoHeart_funcPlot")
alert('jump',1)
