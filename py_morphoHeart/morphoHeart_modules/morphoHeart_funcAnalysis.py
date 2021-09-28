# -*- coding: utf-8 -*-
"""
morphoHeart_funcAnalysis

Version: Nov, 2020
@author: Juliana Sanchez-Posada

"""
#%% Importing python packages
import os
import numpy as np
from sklearn.neighbors import KernelDensity
# from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
# import random

import pandas as pd
import seaborn as sns
from itertools import count, combinations
from progress.bar import Bar
import math
import locale
locale.setlocale(locale.LC_ALL, 'en_US')  

from scipy import stats 
import scikit_posthocs as sp
from statannot import add_stat_annotation

suffix = '%(index)d/%(max)d - %(elapsed)ds'
font = {'family':'sans-serif', 'color':'gray', 'weight':'light', 'size':10, 'horizontalalignment':'center'}

#%% Importing morphoHeart packages
from .morphoHeart_funcBasics import alert, ask4input, getInputNumbers, loadDF

# import matplotlib.font_manager
# matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
#%% - GENERAL LOAD/SAVE
#%% func - q_savePlot
def q_savePlot():
    save = ask4input('Do you want to save the figures? [0]:no/[1]:yes: ', bool)
    if save: 
        info = ask4input('Add info for plot name (eg. hapln1a241s_wtmt, wt): ', str)
        svgpng = ['png','svg']
        svgORpng = ask4input('Do you want to save the created figures as [0]:png, [1]:svg, [2]:both? : ', int)
        if svgORpng == 2:
            svgORpng = svgpng
        else:
            svgORpng = [svgpng[svgORpng]]
    else: 
        info = ''
        svgORpng = ''
    
    return save, info, svgORpng
        
#%% func - getGenotypeAll
def getGenotypeAll(df_input):
    df_input['GenotA'] = df_input['Gene_A']+':'+df_input['Genotype_A']
    df_input['GenotB'] = df_input['Gene_B']+':'+df_input['Genotype_B']
    df_input['GenotypeAll'] = df_input['GenotA']+'/'+df_input['GenotB']
    df_input['Ref'] = df_input['LS_Session']+'_'+df_input['Fish_ref']
    
    # Create new column with complete genotype
    genotypeAll = [] 
    for i, genotA, genotB in zip(count(), df_input["GenotA"], df_input["GenotB"]): 
        if genotB == "-:-": 
            genotypeAll.append(genotA) 
        else: 
            genotypeAll.append(genotA+'/'+genotB) 
    df_input["GenotypeAll"] = genotypeAll 
    
    return df_input

#%% func - getMainStrain
def getMainStrain(df_meas):
    
    strain_o = []
    for i, strain in zip(count(), df_meas['Strain']):
        if '(' in strain:
            strain_o.append(strain.split('(')[0] + 'InX')
        else: 
            strain_o.append(strain)

    df_meas['Strain_o'] = strain_o
    
    return df_meas

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

#%% func - getVariables
def getVariables (list_vars, name):
    """


    Parameters
    ----------
    list_vars : TYPE
        DESCRIPTION.
    name : TYPE
        DESCRIPTION.

    Returns
    -------
    variables2loop : TYPE
        DESCRIPTION.

    """

    print('\n- Variables:')
    for c, value in enumerate(list_vars, 1):
        print(c-1, value)
    input_var = ask4input('Select the '+ name +' you would like to process: ', str)

    if input_var == 'all':
        var_num = list(range(0,len(list_vars),1))

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

    variables2loop = []
    for i, num in enumerate(var_num):
        variables2loop.append(list_vars[num])

    return variables2loop
    
#%% func - selectVariables
def selectVariables (vars_dict):

    variables = list(vars_dict.keys())
    num_variables = [str(num)+'. '+var for num, var in enumerate(variables)]
    # print('\n- Variables:')
    # for c, value in enumerate(variables, 1):
    #     print(c-1, value)
    list_columns(num_variables, cols = 4)
    input_var = ask4input('Select the variables you would like to process: ', str)

    if input_var == 'all':
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
        labels2loop.append(vars_dict[variables[num]])

    return vars2loop, labels2loop

#%% func - selectVariables_auto
def selectVariables_auto (vars_dict, group):

    pl_groups = plot_groups()
    variables = list(vars_dict.keys())
    
    var_num = []
    for gr in group:
        var = pl_groups[gr]['vars']
        for v in var:
            var_num.append(variables.index(v))

    vars2loop = []
    labels2loop = []
    for i, num in enumerate(var_num):
        vars2loop.append(variables[num])
        labels2loop.append(vars_dict[variables[num]])
        
    return vars2loop, labels2loop

#%% - SORT/FILTER DATAFRAMES
#%% func - sortDFCols
def sortDFCols(df_meas):
    
    first_cols = ['Folder', 'LS_Session', 'Fish_ref', 'Ref', 'Strain', 'Strain_o', 'Stage', 'Manip',
       'Gene_A', 'Genotype_A','GenotA', 'Gene_B', 'Genotype_B', 'GenotB', 'GenotypeAll']

    other_cols = []
    for col in df_meas.columns:
        if col == 'spAnalysis':
            first_cols.append(col)
        elif col not in first_cols:
            other_cols.append(col)
    
    other_cols = sorted(other_cols,  key=lambda s: s.lower())
    
    df_meas_sorted = df_meas.reindex(first_cols+other_cols, axis =1)
    
    return df_meas_sorted

#%% func - filterDF2Plot
def filterDF2Plot (df_input, df2plot, df_type = 'meas'):
    """
    Define which df to plot

    Parameters
    ----------
    df_input : TYPE
        DESCRIPTION.

    Returns
    -------
    df_filt : TYPE
        DESCRIPTION.
    genots : TYPE
        DESCRIPTION.
    strains : TYPE
        DESCRIPTION.
    stages : TYPE
        DESCRIPTION.

    """
        
    if isinstance(df2plot, pd.DataFrame):
        print('\n>> FILTERED DATAFRAME: ')
        genots, strains, strains_o, stages = printDFINfo(df2plot, q_return = True, df_type = df_type)
        q_useSameDF = ask4input('Do you want to use the same filtered dataframe? \n\t-[0]: no, I want to  filter the dataset again.\n\t-[1]: yes, use the same one! >>: ',bool)
        if q_useSameDF:
            df_filt = df2plot
        else: #re-filter dataframe
            print('\n>> ORIGINAL DATAFRAME: ')
            df_filt, genots, strains, strains_o, stages = filterDF(df_input, True, df_type)
            
    else: # df2plot == []
        filtDF = ask4input('Do you want to filter dataframe to contain just one group of data? \n\t(e.g. only wild-types, only one strain or stage) [0]:no/[1]:yes: ', bool)
        print('\n>> ORIGINAL DATAFRAME: ')
        df_filt, genots, strains, strains_o, stages = filterDF(df_input, filtDF, df_type)
    
    return df_filt, genots, strains, strains_o, stages

#%% func - filterDF
def filterDF(df2filter, filtDF, df_type):
        # filtDF = ask4input('Do you want to filter dataframe to contain just one group of data? \n\t(e.g. only wild-types, only one strain or stage) [0]:no/[1]:yes: ', bool)
        genots, strains, strains_o, stages = printDFINfo(df2filter, q_return = True, df_type = df_type)
        
        while filtDF:
            filter_type = ask4input('Select the category by which you want to filter the data \n\t[0]: Genotype (wild-type, homozygous mutants)\n\t[1]: Strain (including F#)\n\t[2]: Strain_o (NOT including F#)\n\t[3]: Stage\t>> :', int)
            
            if filter_type == 0:
                str_print = 'Select genotypes to include '
                for i, gen in enumerate(genots):
                    str_print = str_print+'\n\t['+str(i)+']: '+ gen
                str_print = str_print+'\n\t[all]: All >> : '
                gen_filt = ask4input(str_print, str)
                filt_genot = getInputNumbers(gen_filt, genots)
                sel_genots = [genots[j] for j in filt_genot]
                if df_type != 'kde':
                    df2filter = df2filter[df2filter['GenotypeAll'].isin(sel_genots)]
                else: 
                    df2filter = df2filter[df2filter['Genotype'].isin(sel_genots)]
                    
            if filter_type == 1:
                str_print = 'Select strains to include '
                for i, strain in enumerate(strains):
                    str_print = str_print+'\n\t['+str(i)+']: '+ strain
                str_print = str_print+'\n\t[all]: All >> : '
                strain_filt = ask4input(str_print, str)
                filt_strain = getInputNumbers(strain_filt, strains)
                sel_strains = [strains[j] for j in filt_strain]
                df2filter = df2filter[df2filter['Strain'].isin(sel_strains)]
            
            if filter_type == 2:
                str_print = 'Select main strains to include '
                for i, strain_o in enumerate(strains_o):
                    str_print = str_print+'\n\t['+str(i)+']: '+ strain_o
                str_print = str_print+'\n\t[all]: All >> : '
                strain_o_filt = ask4input(str_print, str)
                filt_strain_o = getInputNumbers(strain_o_filt, strains_o)
                sel_strains_o = [strains_o[j] for j in filt_strain_o]
                df2filter = df2filter[df2filter['Strain_o'].isin(sel_strains_o)]
                
            if filter_type == 3:
                str_print = 'Select stages to include '
                for i, stg in enumerate(stages):
                    str_print = str_print+'\n\t['+str(i)+']: '+ stg
                str_print = str_print+'\n\t[all]: All >> : '
                stage_filt = ask4input(str_print, str)
                filt_stage = getInputNumbers(stage_filt, stages)
                sel_stages = [stages[j] for j in filt_stage]
                df2filter = df2filter[df2filter['Stage'].isin(sel_stages)]
            
            if df_type != 'kde':
                print('\n- Filtered dataframe: ')
                print(df2filter[['Ref','Strain','Strain_o', 'Stage','GenotypeAll']])
            else: 
                print('\n- Sample Filtered dataframe: ')
                print(df2filter[['Strain','Strain_o', 'Stage','Genotype']].sample(10))
                
            if df_type != 'kde':
                genots = sorted(df2filter.GenotypeAll.unique(), reverse = True)
            else: 
                genots = sorted(df2filter.Genotype.unique(), reverse = True)
                
            strains = sorted(df2filter.Strain.unique())
            strains_o = sorted(df2filter.Strain_o.unique())
            stages = sorted(df2filter.Stage.unique())
            
            filtDF = ask4input('Do you want to filter the dataframe using another category [0]:no/[1]:yes? : ', bool)
    
        return df2filter, genots, strains, strains_o, stages
    
#%% func - filterR_Autom 
def filterR_Autom (df_input, filters, group, col_out = ''):
        
    df_pivot = df_input.groupby(filters)
    df_out = df_pivot.get_group(group)
    if col_out != '':
        df_out = df_out[col_out]
    
    return df_out

#%% func - printDFINfo
def printDFINfo(df2plot, q_return = False, df_type = 'meas'):
    print('\n>> Information found in the imported dataframe:')
    if df_type != 'kde':
        genot = sorted(df2plot.GenotypeAll.unique(), reverse = True)
    else: 
        genot = sorted(df2plot.Genotype.unique(), reverse = True)
    strains = sorted(df2plot.Strain.unique())
    strains_o = sorted(df2plot.Strain_o.unique())
    stages = sorted(df2plot.Stage.unique())
    
    print(' - Genotypes: ', genot)
    print(' - Strains: ', strains)
    print(' - Strains_o: ', strains_o)
    print(' - Stages: ', stages)
    
    input()
    if q_return:
        return genot, strains, strains_o, stages
    
#%% - CALCULATE/DEFINE VARIABLES, LEGENDS AND TITLES
#%% func - getVarRatios
def getVarRatios(df_meas):
    
    df_meas['Looping_Ratio_Myoc'] = df_meas['Length_CL_Int.Myoc(Cut)']/ df_meas['linLine_Int.Myoc(Cut)']
    df_meas['Looping_Ratio_Endo'] = df_meas['Length_CL_Ext.Endo(Cut)']/ df_meas['linLine_Ext.Endo(Cut)']
    df_meas['Vol_Tissue'] = df_meas['Vol_Ext.Myoc'] - df_meas['Vol_Int.Endo']
    df_meas['Vol_Atr.Tissue'] = df_meas['Vol_Atr.ExtMyoc'] - df_meas['Vol_Atr.IntEndo']
    df_meas['Vol_Vent.Tissue'] = df_meas['Vol_Vent.ExtMyoc'] - df_meas['Vol_Vent.IntEndo']
    
    
    df_meas['Ratio_VolMyoc2VolExtMyoc'] = df_meas['Vol_Myoc']/ df_meas['Vol_Ext.Myoc']
    df_meas['Ratio_VolEndo2VolExtMyoc'] = df_meas['Vol_Myoc']/ df_meas['Vol_Ext.Myoc']
    df_meas['Ratio_VolCJ2VolExtMyoc'] = df_meas['Vol_CJ']/ df_meas['Vol_Ext.Myoc']
    df_meas['Ratio_VolLumen2VolExtMyoc'] = df_meas['Vol_Int.Endo']/ df_meas['Vol_Ext.Myoc']
    
    df_meas['Ratio_VolMyoc2VolTissue'] = df_meas['Vol_Myoc']/ df_meas['Vol_Tissue']
    df_meas['Ratio_VolEndo2VolTissue'] = df_meas['Vol_Myoc']/ df_meas['Vol_Tissue']
    df_meas['Ratio_VolCJ2VolTissue'] = df_meas['Vol_CJ']/ df_meas['Vol_Tissue']
    df_meas['Ratio_VolLumen2VolTissue'] = df_meas['Vol_Int.Endo']/ df_meas['Vol_Tissue']
    
    df_meas['Ratio_VolLeftCJ2VolRightCJ'] = df_meas['Vol_CJ.Left']/ df_meas['Vol_CJ.Right']
    df_meas['Ratio_VolAtrCJ2VolVentCJ'] = df_meas['Vol_Atr.CJ']/ df_meas['Vol_Vent.CJ']
    
    
    
    return df_meas

#%% func - def_variables
def def_variables(plot_type):

    if plot_type == 'strip_plots':
        
        vars_dict = {   # Angles
                       "ang_Heart" : "Angle Heart wrt Sample (\N{DEGREE SIGN})",
                       "ang_AtrS" : "Sagittal Atrial Angle (\N{DEGREE SIGN})",
                       "ang_VentS" : "Sagittal Ventricular Angle (\N{DEGREE SIGN})",
                       "ang_BtwChambersS" : "Sagittal Angle between Chambers (\N{DEGREE SIGN})",
                       
                       "ang_AtrV" : "Ventral Atrial Angle (\N{DEGREE SIGN})",
                       "ang_VentV" : "Ventral Ventricular Angle (\N{DEGREE SIGN})",
                       "ang_BtwChambersV" : "Ventral Angle between Chambers (\N{DEGREE SIGN})",
                       
                       # Ellipse measurements
                       "EllipAtr_Asphericity" : "Atrial Ellipse Asphericity",
                       "EllipAtr_Depth" : "Atrial depth [um]",
                       "EllipAtr_Length" : "Atrial length [um]",
                       "EllipAtr_Width" : "Atrial width [um]",
                       "EllipVent_Asphericity" : "Ventricular Ellipse Asphericity",
                       "EllipVent_Depth" :"Ventricular depth [um]",
                       "EllipVent_Length" : "Ventricular length [um]",
                       "EllipVent_Width" : "Ventricular width [um]",
                       
                       # Surface Areas
                       "SurfArea_Myoc" : "Surface Area\nMyocardium [um$^2$]",
                       "SurfArea_Int.Myoc" : "Surface Area\nInt.Myocardium [um$^2$]",
                       "SurfArea_Ext.Myoc" : "Surface Area\nExt.Myocardium [um$^2$]",
                       "SurfArea_Endo" : "Surface Area\nEndocardium [um$^2$]",
                       "SurfArea_Int.Endo" : "Surface Area\nInt.Endocardium [um$^2$]",
                       "SurfArea_Ext.Endo" : "Surface Area\nExt. Endocardium [um$^2$]",
                       "SurfArea_CJ" : "Surface Area\nCardiac Jelly [um$^2$]",
                       "SurfArea_Int.CJ" : "Surface Area\nInt.Cardiac Jelly [um$^2$]",
                       "SurfArea_Ext.CJ" : "Surface Area\nExt.Cardiac Jelly [um$^2$]",
                       "SurfArea_Atr.ExtCJ" : "Atrial Surface Area\nExternal Cardiac Jelly [um$^2$]",
                       "SurfArea_Atr.ExtMyoc" : "Atrial Heart\nSurface Area [um$^2$]",
                       "SurfArea_Atr.IntEndo" : "Atrial Lumen\nSurface Area [um$^2$]",
                       "SurfArea_Vent.ExtCJ" : "Ventricular Surface Area\nExternal Cardiac Jelly [um$^2$]",
                       "SurfArea_Vent.ExtMyoc" : "Ventricular Heart\nSurface Area [um$^2$]",
                       "SurfArea_Vent.IntEndo" : "Ventricular Lumen\nSurface Area [um$^2$]",
                        
                       # Looping ratio
                       "linLine_Int.Myoc(Cut)" : "Linear Heart Length\n(Int.Myoc) [um]",
                       "linLine_Ext.Endo(Cut)" : "Linear Heart Length\n(Ext.Endo) [um]",
                       "Length_CL_Int.Myoc(Cut)" : "Looped Heart Length\n(Int.Myoc) [um]",
                       "Length_CL_Ext.Endo(Cut)" : "Looped Heart Length\n(Ext.Endo) [um]",
                       'Looping_Ratio_Myoc' : 'Looping Ratio (Int.Myoc)',
                       'Looping_Ratio_Endo' : 'Looping Ratio (Ext.Endo)',
                        
                       # Volumes
                       "Vol_Int.Myoc" : "Volume Int.Myocardium [um$^3$]",
                       "Vol_Ext.Myoc" : "Heart Volume [um$^3$]",
                       "Vol_Int.Endo" : "Heart Lumen Volume [um$^3$]",
                       "Vol_Ext.Endo" : "Volume Ext.Endocardium [um$^3$]",
                       "Vol_Myoc" : "Volume Myocardium [um$^3$]",
                       "Vol_Atr.Myoc" : "Atrial Volume\nMyocardium [um$^3$]",
                       "Vol_Vent.Myoc" : "Ventricular Volume\nMyocardium [um$^3$]",
                       "Vol_Endo" : "Volume Endocardium [um$^3$]",
                       "Vol_Atr.Endo" : "Atrial Volume\nEndocardium [um$^3$]",
                       "Vol_Vent.Endo" : "Ventricular Volume\nEndocardium [um$^3$]",
                       "Vol_CJ" : "Volume Cardiac Jelly [um$^3$]",
                       "Vol_Atr.CJ" : "Atrial Volume\nCardiac Jelly [um$^3$]",
                       "Vol_Vent.CJ" : "Ventricular Volume\nCardiac Jelly [um$^3$]",
                       'Vol_Atr.ExtMyoc' : 'Atrial Volume [um$^3$]',
                       'Vol_Vent.ExtMyoc' : 'Ventricular Volume [um$^3$]',
                       'Vol_Atr.IntEndo' : 'Atrial Lumen Volume [um$^3$]',
                       'Vol_Vent.IntEndo' : 'Ventricular Lumen Volume [um$^3$]',
                       
                       'Vol_CJ.Left' : "Left Volume\nCardiac Jelly [um$^3$]",
                       'Vol_CJ.Right' : "Right Volume\nCardiac Jelly [um$^3$]",
                       
                       'Vol_Tissue' : 'Tissue Volume [um$^3$]',
                       'Vol_Atr.Tissue' : 'Atrial Tissue Volume [um$^3$]', 
                       'Vol_Vent.Tissue' : 'Ventricular Tissue Volume [um$^3$]',
                       
                       # Volume Ratios 
                       'Ratio_VolMyoc2VolExtMyoc' : 'Myocardial Volume / \nHeart Volume', 
                       'Ratio_VolEndo2VolExtMyoc' : 'Endocardial Volume / \nHeart Volume',
                       'Ratio_VolCJ2VolExtMyoc' : 'Cardiac Jelly Volume / \nHeart Volume', 
                       'Ratio_VolLumen2VolExtMyoc' : 'Lumen Volume / \nHeart Volume',
                       
                       'Ratio_VolMyoc2VolTissue' : 'Myocardial Volume / \nTissue Volume', 
                       'Ratio_VolEndo2VolTissue' : 'Endocardial Volume / \nTissue Volume',
                       'Ratio_VolCJ2VolTissue' : 'Cardiac Jelly Volume / \nTissue Volume', 
                       'Ratio_VolLumen2VolTissue' : 'Lumen Volume / \nTissue Volume', 
                       
                       'Ratio_VolLeftCJ2VolRightCJ' : 'Left Vol.Cardiac Jelly / \nRight Vol.Cardiac Jelly', 
                       'Ratio_VolAtrCJ2VolVentCJ' : 'Atrial Vol.Cardiac Jelly / \nVentricular Vol.Cardiac Jelly'}
        
    
    elif plot_type == 'bar_plots':
        
        vars_dict = {"Vol_Int.Myoc" : "Volume Int.Myocardium [um$^3$]",
                       "Vol_Ext.Myoc" : "Heart\nVolume ",
                       "Vol_Int.Endo" : "Heart\nLumen",
                       "Vol_Ext.Endo" : "Volume Ext.Endocardium [um$^3$]",
                       "Vol_Myoc" : "Myocardium",
                       "Vol_Atr.Myoc" : "Atrial\nMyocardium",
                       "Vol_Vent.Myoc" : "Ventricular\nMyocardium",
                       "Vol_Endo" : "Endocardium",
                       "Vol_Atr.Endo" : "Atrial\nEndocardium",
                       "Vol_Vent.Endo" : "Ventricular\nEndocardium",
                       "Vol_CJ" : "Cardiac\nJelly",
                       "Vol_Atr.CJ" : "Atrial\nCardiac Jelly",
                       "Vol_Vent.CJ" : "Ventricular\nCardiac Jelly",
                       'Vol_Atr.ExtMyoc' : 'Atrium',
                       'Vol_Vent.ExtMyoc' : 'Ventricle',
                       'Vol_Atr.IntEndo' : 'Atrial\nLumen',
                       'Vol_Vent.IntEndo' : 'Ventricular\nLumen',
                       "Vol_Myoc_Change" : "Myocardium",
                       "Vol_Endo_Change" : "Endocardium",
                       "Vol_CJ_Change" : "Cardiac Jelly",
                       "Vol_Atr.Myoc_Change" : "Atrial\nMyocardium",
                       "Vol_Atr.CJ_Change" : "Atrial\nCardiac Jelly",
                       "Vol_Atr.Endo_Change" : "Atrial\nEndocardium",
                       "Vol_Vent.Myoc_Change" : "Ventricular\nMyocardium",
                       "Vol_Vent.CJ_Change" : "Ventricular\nCardiac Jelly",
                       "Vol_Vent.Endo_Change" : "Ventricular\nEndocardium",
                       "Vol_Int.Endo_Change" : "Heart\nLumen",
                       "Vol_Atr.IntEndo_Change" : "Atrial\nLumen",
                       "Vol_Vent.IntEndo_Change" : "Ventricular\nLumen"}

    return vars_dict

#%% func - plot_groups
def plot_groups():
    pl_groups = {'Surface Area': 
                     {'title': 'Surface Areas', 
                      'vars' : ['SurfArea_Ext.Myoc','SurfArea_Int.Endo'],
                      'n_cols': 2, 'yticks_lab':'1e3 - d.', 'ylim' : ''},
                 'Surface Area CJ': 
                     {'title': 'Cardiac Jelly Thickness Analysis', 
                      'vars' : ['SurfArea_Ext.CJ'],
                      'n_cols': 1, 'yticks_lab':'th,', 'ylim' : ''},
                 'Heart Size': 
                     {'title': 'Heart Size', 
                      'vars' : ['Vol_Ext.Myoc', 'Vol_Atr.ExtMyoc', 'Vol_Vent.ExtMyoc'],
                      'n_cols': 3, 'yticks_lab':'1e6 - d.', 'ylim' : ''},
                 'Lumen Size' : 
                     {'title': 'Lumen Size', 
                      'vars' : ['Vol_Int.Endo', 'Vol_Atr.IntEndo', 'Vol_Vent.IntEndo'],
                      'n_cols': 3, 'yticks_lab':'1e6 - d.', 'ylim' : ''}, 
                 'Heart Looping': 
                     {'title': 'Heart Looping', 
                      'vars' : ['linLine_Int.Myoc(Cut)', 'Length_CL_Int.Myoc(Cut)', 'Looping_Ratio_Myoc'],
                      'n_cols': 3, 'yticks_lab':'d.', 'ylim' : ''}, 
                 'Tissue Myocardium':  
                     {'title': 'Tissue Layer Volumes - Myocardium', 
                      'vars' : ['Vol_Myoc', 'Vol_Atr.Myoc', 'Vol_Vent.Myoc'],
                      'n_cols': 3, 'yticks_lab':'1e3 - d.', 'ylim' : ''},
                 'Tissue Endocardium': 
                     {'title': 'Tissue Layer Volumes - Endocardium', 
                      'vars' : ['Vol_Endo', 'Vol_Atr.Endo', 'Vol_Vent.Endo'],
                      'n_cols': 3, 'yticks_lab':'1e3 - d.', 'ylim' : ''},
                 'Tissue Cardiac Jelly (AtrVent)' : 
                     {'title': 'Tissue Layer Volumes - Cardiac Jelly (Atr&Vent)', 
                      'vars' : ['Vol_CJ', 'Vol_Atr.CJ', 'Vol_Vent.CJ'],
                      'n_cols':3, 'yticks_lab':'1e3 - d.', 'ylim' : ''},
                 'Tissue Cardiac Jelly (LR)' : 
                     {'title': 'Tissue Layer Volumes - Cardiac Jelly (L&R)', 
                      'vars' : ['Vol_CJ','Vol_CJ.Left', 'Vol_CJ.Right'],
                      'n_cols':3, 'yticks_lab':'1e3 - d.', 'ylim' : ''},
                 'Ratios Cardiac Jelly' : 
                     {'title': 'Ratios Tissue Layer Volumes - Cardiac Jelly', 
                      'vars' : ['Ratio_VolAtrCJ2VolVentCJ', 'Ratio_VolLeftCJ2VolRightCJ'],
                      'n_cols': 2, 'yticks_lab':'d.', 'ylim' : ((2,12), (0,35))},
                 'AtrEllipse' : 
                     {'title': 'Atrial Ellipse', 
                      'vars' : ['EllipAtr_Depth','EllipAtr_Length', 'EllipAtr_Width'],
                      'n_cols': 3, 'yticks_lab':'d. - 0', 'ylim' : ''},
                 'VentEllipse' : 
                     {'title': 'Ventricular Ellipse', 
                      'vars' : ['EllipVent_Depth','EllipVent_Length', 'EllipVent_Width'],
                      'n_cols': 3, 'yticks_lab':'d. - 0', 'ylim' : ''},
                 'ChamberAsphericity' : 
                     {'title': 'Chamber Asphericity', 
                      'vars' : ['EllipAtr_Asphericity','EllipVent_Asphericity'],
                      'n_cols': 2, 'yticks_lab':'d.', 'ylim' : ''},
                 'Sagittal Angles' :
                     {'title': 'Sagittal Angles', 
                      'vars' : ['ang_AtrS', 'ang_VentS', 'ang_BtwChambersS'],
                      'n_cols': 3, 'yticks_lab':'d. - 0', 'ylim' : ''},
                 'Ventral Angles' :
                     {'title': 'Ventral Angles', 
                      'vars' : ['ang_AtrV', 'ang_VentV', 'ang_BtwChambersV'],
                      'n_cols': 3, 'yticks_lab':'d. - 0', 'ylim' : ''},
                 'Volume Percentages' : 
                     {'title': 'Heart Composition', 
                      'vars' : ['Ratio_VolMyoc2VolExtMyoc', 'Ratio_VolEndo2VolExtMyoc', 'Ratio_VolCJ2VolExtMyoc', 'Ratio_VolLumen2VolExtMyoc'],
                      'n_cols': 4, 'yticks_lab':'d.', 'ylim' : ''},
                 'Tissue Volume Percentages' : 
                     {'title': 'Tissue Volume Composition', 
                      'vars' : ['Ratio_VolMyoc2VolTissue', 'Ratio_VolEndo2VolTissue', 'Ratio_VolCJ2VolTissue', 'Ratio_VolLumen2VolTissue'],
                      'n_cols': 4, 'yticks_lab':'d.', 'ylim' : ''}}
    
    return pl_groups

#%% func - def_legends
def def_legends(df_input, df_type = 'meas'):
    """
    

    Parameters
    ----------
    genots : TYPE
        DESCRIPTION.
    strains : TYPE
        DESCRIPTION.
    stages : TYPE
        DESCRIPTION.

    Returns
    -------
    out_genots : TYPE
        DESCRIPTION.
    out_strains : TYPE
        DESCRIPTION.
    out_stages : TYPE
        DESCRIPTION.

    """
    
    # GENOTYPES
    if df_type != 'kde':
        try: 
            genots = sorted(df_input.GenotypeAll.unique(), reverse = True)
            bool_genots = True
        except: 
            bool_genots = False
    else: 
        try: 
            genots = sorted(df_input.Genotype.unique(), reverse = True)
            bool_genots = True
        except: 
            bool_genots = False
    
    out_genots = []
    if bool_genots: 
        all_genots = ['hapln1a:wt', 'hapln1a:ht', 'hapln1a:mt', 
                      'hapln1a:wt/spaw:wt', 'hapln1a:wt/spaw:ht', 'hapln1a:wt/spaw:mt', 
                      'hapln1a:ht/spaw:wt', 'hapln1a:ht/spaw:ht', 'hapln1a:ht/spaw:mt',
                      'hapln1a:mt/spaw:wt', 'hapln1a:mt/spaw:ht', 'hapln1a:mt/spaw:mt',
                      'galff:+/uas:+', 'galff:+/uas:-', 'galff:-/uas:+', 'galff:-/uas:-',
                      'vcana:wt', 'vcana:ht', 'vcana:mt', 'wt:wt']
        
        leg_genots = ['$hapln1a^{+/+}$', '$hapln1a^{+/-}$', '$hapln1a^{-/-}$', 
                      '$hapln1a^{+/+}/spaw^{+/+}$', '$hapln1a^{+/+}/spaw^{+/-}$', '$hapln1a^{+/+}/spaw^{-/-}$', 
                      '$hapln1a^{+/-}/spaw^{+/+}$', '$hapln1a^{+/-}/spaw^{+/-}$', '$hapln1a^{+/-}/spaw^{-/-}$',
                      '$hapln1a^{-/-}/spaw^{+/+}$', '$hapln1a^{-/-}/spaw^{+/-}$', '$hapln1a^{-/-}/spaw^{-/-}$',
                      'galff:+/uas:+', 'galff:+/uas:-', 'galff:-/uas:+', 'hapln1a_OE:wt',
                      '$vcana^{+/+}$', '$vcana^{+/-}$', '$vcana^{-/-}$', 'wt']
        
        for gen in genots:
            out_genots.append(leg_genots[all_genots.index(gen)])
    
    # STRAINS
    try: 
        strains = sorted(df_input.Strain.unique())
        bool_strains = True
    except: 
        bool_strains = False
    
    out_strains = []
    if bool_strains: 
        all_strains = ['hapln1a prom187/+ (F2s) InX','hapln1a prom187/+ (F3s) InX', 
                       'hapln1a prom241/+ (F2s) InX', 'hapln1a prom241/+ (F3s) InX',
                       'spaw+/-; hapln1a prom241/+ InX', 'vcana prom365/+ (F2s) InX',
                       'myl7lckGFP_kdrlRasCherry',
                       'myl7:galFF; UAS:TFP x UAS:hapln1a, cryaa:CFP']
        
        leg_strains = ["$hapln1a^{\Delta 187} (F2s)$", "$hapln1a^{\Delta 187} (F3s)$", 
                       "$hapln1a^{\Delta 241} (F2s)$", "$hapln1a^{\Delta 241} (F3s)$",
                       "$spaw^{+/-}; hapln1a^{\Delta 241}$", "$vcana^{\Delta 365}$",
                       'myl7:lck-GFP; kdrl:rasCherry',
                       'myl7:galFF; UAS:TFP x UAS:hapln1a, cryaa:CFP']
    
        for strain in strains:
            out_strains.append(leg_strains[all_strains.index(strain)])
        
    # MAIN STRAINS
    try: 
        strains_o = sorted(df_input.Strain_o.unique())
        bool_strains_o = True
    except: 
        bool_strains_o = False
    
    out_strains_o = []
    if bool_strains_o: 
        all_strains_o = ['hapln1a prom187/+ InX', 'hapln1a prom241/+ InX',
                       'spaw+/-; hapln1a prom241/+ InX', 'vcana prom365/+ (F2s) InX',
                       'myl7lckGFP_kdrlRasCherry',
                       'myl7:galFF; UAS:TFP x UAS:hapln1a, cryaa:CFP']
        
        leg_strains_o = ["$hapln1a^{\Delta 187}$","$hapln1a^{\Delta 241}$", 
                       "$spaw^{+/-}; hapln1a^{\Delta 241}$", "$vcana^{\Delta 365}$",
                       'myl7:lck-GFP; kdrl:rasCherry',
                       'myl7:galFF; UAS:TFP x UAS:hapln1a, cryaa:CFP']
        
        for strain_o in strains_o:
            out_strains_o.append(leg_strains_o[all_strains_o.index(strain_o)])
    
    # STAGES
    try: 
        stages = sorted(df_input.Stage.unique())
        bool_stages = True
    except: 
        bool_stages = False
    
    out_stages = []
    if bool_stages:
        if df_type != 'changes':
            all_stages = ['32-34', '48-50', '58-60', '72-74']
            leg_stages = ['32-34hpf', '48-50hpf', '58-60hpf', '72-74hpf']
        else: 
            all_stages = ['32->48', '48->74']
            leg_stages = ["32$\Rightarrow$48hpf", "48$\Rightarrow$74hpf"]
        
        for stg in stages:
            try: 
                out_stages.append(leg_stages[all_stages.index(stg)])
            except :
                out_stages.append('')
    
    # TIME-POINTS
    try: 
        time_points = sorted(df_input.time_point.unique())
        bool_tp = True
    except: 
        bool_tp = False
    
    out_tp = []
    if bool_tp:
        # if df_type != 'changes':
        #     all_stages = ['32-34', '48-50', '58-60', '72-74']
        #     leg_stages = ['32-34hpf', '48-50hpf', '58-60hpf', '72-74hpf']
        # else: 
        #     all_stages = ['32->48', '48->74']
        #     leg_stages = ["32$\Rightarrow$48hpf", "48$\Rightarrow$74hpf"]
        
        for tp in time_points:
            # try: 
            out_tp.append(tp)#leg_stages[all_stages.index(stg)])
            # except :
            #     out_tp.append('')
    
    
    x_labels = {'GenotypeAll': 'Genotype', 'Genotype': 'Genotype',
                    'Strain' : 'Strain', 'Strain_o' : 'Strain_o',
                    'Stage': 'Stage [hpf]', 'time_point': 'Heart phase contraction'}
                
    dict_legends = {'GenotypeAll': out_genots, 'Genotype': out_genots,
                    'Strain' : out_strains, 'Strain_o' : out_strains_o,
                    'Stage': out_stages, 'time_point' : out_tp, 'xlabels': x_labels}

    return dict_legends 

#%% func - def_var_names 
def def_var_names(plot_type, labels):
    vars_dict =  def_variables(plot_type)
    names = []
    for lab in labels: 
        names.append(vars_dict[lab])
    
    return names

#%% func - get_txt_title
def get_txt_title(vars_opt, df2plot):
    txt_title = '\n'
    for var in vars_opt:
        if len(df2plot[var].unique()) <= 1:
            txt_title = txt_title +'['+ var + ':'+df2plot[var].unique()[0]+'] '
    
    return txt_title

#%% - MODIFY DFs TO PLOT
#%% func - get_dfPctChange
def get_dfPctChange(df2plot, vars2plot, group_vars):
    
    def get_pct_change(df_all_stages, df_pct_change):
        if (df_all_stages.index == ['32-34','48-50','72-74']).all():
            # print(df_all_stages)
            # print(df_all_stages.diff())
            df_change = df_all_stages.pct_change()*100
            df_change = df_change.dropna()
            df_change = df_change.rename(columns = lambda s: s+'_Change', index={'48-50':'32->48', '72-74':'48->74'}) # index={'48-50':'32->50', '72-74':'48->74'})
            # print(df_change)
            if bool_genot:
                df_change ['GenotypeAll'] = [gen for nn in range(len(stages)-1)]
            if bool_strains_o:
                df_change['Strain_o'] = [strain for nn in range(len(stages)-1)]
            df_change = df_change.reset_index().set_index(new_index)
            # print(df_change)
            df_pct_change.append(df_change)
            
    bool_strains_o = False
    bool_genot = False
    
    df4plot = df2plot.groupby(group_vars)[vars2plot].mean()
    new_index = []
    if 'Strain_o' in group_vars: 
        strains_o = df2plot.Strain_o.unique()
        bool_strains_o = True
        new_index.append('Strain_o')
    if 'GenotypeAll' in group_vars:
        genots = df2plot.GenotypeAll.unique()
        bool_genot = True
        new_index.append('GenotypeAll')
    new_index.append('Stage')
    
    stages = sorted(df2plot.Stage.unique())
    
    df_pct_change = []
    if bool_genot:
        for gen in genots:
            if bool_strains_o:
                for strain in strains_o:
                    # print('\n',gen, strain)
                    df_all_stages = df4plot.loc[gen,strain]
                    get_pct_change(df_all_stages, df_pct_change)
            else: 
                # print('\n',gen)
                df_all_stages = df4plot.loc[gen]
                get_pct_change(df_all_stages, df_pct_change)
    else: 
        if bool_strains_o:
            for strain in strains_o:
                # print('\n',strain)
                df_all_stages = df4plot.loc[strain]
                get_pct_change(df_all_stages, df_pct_change)
    
        
    # for gen in genots:
    #     for strain in strains_o:
    #         print('\n',gen, strain)
    #         df_all_stages = df4plot.loc[gen,strain]
    #         if (df_all_stages.index == ['32-34','48-50','72-74']).all():
    #             print(df_all_stages)
    #             print(df_all_stages.diff())
    #             df_change = df_all_stages.pct_change()*100
    #             df_change = df_change.dropna()
    #             df_change = df_change.rename(columns = lambda s: s+'_Change', index={'48-50':'32->48', '72-74':'48->74'}) # index={'48-50':'32->50', '72-74':'48->74'})
    #             print(df_change)
    #             df_change ['GenotypeAll'] = [gen for nn in range(len(stages)-1)]
    #             df_change['Strain_o'] = [strain for nn in range(len(stages)-1)]
    #             df_change = df_change.reset_index().set_index(['GenotypeAll', 'Strain_o', 'Stage'])
    #             print(df_change)
    #             df_pct_change.append(df_change)
    df4plot_pct_change = pd.concat(df_pct_change)
    
    return df_pct_change, df4plot_pct_change

#%% func - meanHM
def meanHM(df_dataset_hm, filters, groups, chamber, variable, opt_norm,  dir2load_df, dir2save_hmf, dir_data2Analyse, save):
    
    # Get minimum and max values for each group of heatmaps 
    # (_o: originals, _W: normalised using whole heart, _C: normalised per chamber)
    normalise, perChamber, norm_type = opt_norm
    min_vals_o = []; max_vals_o = []; min_vals_N = []; max_vals_N = []
    for group in groups:
        folders = filterR_Autom (df_dataset_hm, filters, group, col_out = 'Folder')
        
        print('\n >> Variable: ', variable,' - Group: ', group,' - Chamber: ', chamber)
        print(folders)
        if normalise: 
            [dfs_o, dfs_hmN], num, _  = getHeatmaps2Unify(folders, chamber, variable, dir2load_df, 
                                                       dir_data2Analyse, normalise = normalise, 
                                                       opt_norm = norm_type, perChamber = perChamber)
            df_hmW = concatHeatmaps(dfs_hmN, operation = 'mean')
            min_vals_N.append(df_hmW.min().min())
            max_vals_N.append(df_hmW.max().max())
            
        else: 
            [dfs_o], num  = getHeatmaps2Unify(folders, chamber, variable, dir2load_df, 
                                                       dir_data2Analyse, normalise = normalise, 
                                                       opt_norm = norm_type, perChamber = perChamber)

        df_of = concatHeatmaps(dfs_o, operation = 'mean')
        min_vals_o.append(df_of.min().min())
        max_vals_o.append(df_of.max().max())
        
    min_o = math.floor(min(min_vals_o))
    max_o = math.ceil(max(max_vals_o))
    if normalise: 
        max_N = math.ceil(max(max_vals_N))

    for group in groups:
        folders = filterR_Autom (df_dataset_hm, filters, group, col_out = 'Folder')
        print('\n >> Variable: ', variable,' - Group: ', group,' - Chamber: ', chamber)
        print('\t - Maximum value_o:', max_o)
        if normalise: 
            print('\t - Maximum value_N:', max_N)
            [dfs_o, dfs_hmN], num, norm_vals  = getHeatmaps2Unify(folders, chamber, variable, dir2load_df, 
                                                       dir_data2Analyse, normalise = normalise, 
                                                       opt_norm = norm_type, perChamber = perChamber)
            df_hmW = concatHeatmaps(dfs_hmN, operation = 'mean')
            norm_vals = ['{:.2f}'.format(val) for val in norm_vals]
            print('\t - Normalisation: ', norm_type, ' - Norm_val:', norm_vals)
        else: 
            [dfs_o], num  = getHeatmaps2Unify(folders, chamber, variable, dir2load_df, 
                                                       dir_data2Analyse, normalise = normalise, 
                                                       opt_norm = norm_type, perChamber = perChamber)
        df_of = concatHeatmaps(dfs_o, operation = 'mean')

        gen_info = list(group[1])
        if '/' in gen_info:
            ind_sep = gen_info.index('/') 
            gen_info = list(gen_info); gen_info[ind_sep] = '_'; 
        ind_gen = list(filter(lambda x: gen_info[x] == ':', range(len(gen_info)))) 
        for ind in ind_gen:
            gen_info[ind+1] = gen_info[ind+1].upper()
        if len(ind_gen) == 2:
            ind_gen[1] -= 1
        [gen_info.pop(ind) for ind in ind_gen]
        gen_info  = ''.join(gen_info)
        unifyHeatmap(df_of, chamber, genotype=group[1], gen_info = gen_info, stage=group[0], thickness= variable, 
                  vmin=min_o, vmax=max_o, n_val = num, normalise = 'o_all', dir2save = dir2save_hmf, save = save, cmap = 'turbo')
        if normalise:
            if perChamber:
                txt_title = 'normCh-'+norm_type
            else: 
                txt_title = 'normWh-'+norm_type
            unifyHeatmap(df_hmW, chamber, genotype=group[1], gen_info = gen_info, stage=group[0], thickness= variable, 
                      vmin=0, vmax=max_N, n_val = num, normalise = txt_title,  dir2save = dir2save_hmf, 
                      save = save, cmap = 'inferno')

#%% func - concatHeatmaps
def concatHeatmaps(df2concat, operation = 'mean'):
    
    if operation == 'std':
        df_concat = pd.concat(df2concat).groupby(level=0).std()
    elif operation == 'mean':
        df_concat = pd.concat(df2concat).groupby(level=0).mean()
    elif operation == 'sem':
        df_concat = pd.concat(df2concat).groupby(level=0).sem()

    return df_concat

#%% - STATISTICAL ANALYSIS
#%% func - normalityTests
def normalityTests (alpha, factor, fact_levels, variable, data, test):
    
    pval_norm = np.ones(((len(fact_levels)),3))
    
    all_data_normal = True
    txt_normtest = test +' Normality tests ('+r'$\alpha$:'+str(alpha)+'), \n- p-val: - '
    data_groups = []
    
    for num, level in enumerate(fact_levels):
        data2testNorm = data.groupby(factor)[variable].get_group(level)
        print(' > Factor level:',level)
        # print(data2testNorm)

        # D’Agostino and Pearson’s - stats.normaltest()
        if test == 'D''Agostino & Pearson':
            k1, p1 = stats.normaltest(data2testNorm)
            pval_norm[num,0] = p1
            p = p1
        # Shapiro-Wilk test - stats.shapiro()
        elif test == 'Shapiro-Wilk':
            k2, p2 = stats.shapiro(data2testNorm)
            pval_norm[num,1] = p2
            p = p2
        # Kolmogorov-Smirnov test - stats.kstest(rvs, cdf, args=(), N=20, alternative='two-sided', mode='approx')
        elif test == 'Kolmogorov-Smirnov':
            k3, p3 = stats.kstest(data2testNorm, 'norm')
            pval_norm[num,2] = p3
            p = p3
        
        if p >= alpha:  # null hypothesis: factor comes from a normal distribution
            print("\t-",level, ": Data is normally distributed,: p-value={:.5f}".format(p))
        else:
            print("\t-",level, ": Data is NOT normally distributed: p-value={:.5f}".format(p))
            all_data_normal = False
        
        normtxt = level+':'+str(format(p,'.3f'))
        if level != fact_levels[-1]:
            txt_normtest = txt_normtest + normtxt +' / '
        else: 
            txt_normtest = txt_normtest + normtxt + ' - '  
        
        data_groups.append(data2testNorm)
    print(" >> Data normally distributed?", all_data_normal)
    
    return pval_norm, txt_normtest, all_data_normal, data_groups
    
# D’Agostino and Pearson’s - stats.normaltest()
# Shapiro-Wilk test - stats.shapiro()
# Kolmogorov-Smirnov test - stats.kstest(rvs, cdf, args=(), N=20, alternative='two-sided', mode='approx')

#%% func - selectStatisticalTest
def selectStatisticalTest (fact_levels, all_data_normal, variable):
    #Define variables for which the analysis should always be non-parametric
    vars4nonParam = ['Looping_Ratio', 'Index_Value']
    
    # Define statistical test to use 
    if variable not in vars4nonParam:
        if len(fact_levels) > 2:
            if all_data_normal == False: 
                test2use = 'Kruskal'
                test_txt = "Non-Parametric test \n Kruskal-Wallis"
                multcomp_txt = "Dunn's (Bonferroni Corr.)"
            else: 
                test2use = 'One-way-ANOVA'
                test_txt = 'Parametric test -\n One-way-ANOVA'
                multcomp_txt = "Tukey's Multiple Comparison"
        
        elif len(fact_levels) == 2:
            if all_data_normal == False: 
                test2use = 'Mann-Whitney'
                test_txt = 'Non-Parametric test -\n Mann-Whitney-Wilcoxon \ntest two-sided'#', Bonferroni Correction'
                multcomp_txt = ''
            else: 
                test2use = 't-test_ind'
                test_txt = 'Parametric test -\n Independend t-test'
                multcomp_txt = ''
    else:
        print(" - Non-parametric variable!")
        if len(fact_levels) > 2:
            test2use = 'Kruskal'
            test_txt = "Non-Parametric test \n Kruskal-Wallis"
            multcomp_txt = "Dunn's (Bonferroni Corr.)"
        elif len(fact_levels) == 2:
            test2use = 'Mann-Whitney'
            test_txt = 'Non-Parametric test \n Mann-Whitney-Wilcoxon \ntest two-sided'#', Bonferroni Correction'
            multcomp_txt = ''
        
    print(" >> Test selected to use: ", test2use)
    txt_testSelected = 'Test: '+test_txt
    
    return test2use, multcomp_txt, txt_testSelected

#%% func - runStatisticalTests
def runStatisticalTests(data, filters, norm_test, box_pairs_all, box_pairs_f, vars_group, alpha):
    
    list_of_lists = True
    var_filt = []
    opt_filt = []
    for filt in filters: 
        for filt_val in sorted(data[filt].unique()):
            var_filt.append(filt)
            opt_filt.append(filt_val)
    
    dict_stats = dict()
    for var in vars_group: 
        print("\n>>>> Running normality tests for variable:", var)#, "/ allele: ", alleleName, '\n')
        dict_spStatisticalRes = dict()
        pvals_mc = []
        txt_selected = []
        txt_norm = []
        for n, pair_f, pair_all in zip(count(), box_pairs_f, box_pairs_all): 
            
            # print('pair_all:', pair_all)
            inter = list(set(pair_all[0]).intersection(*pair_all[1:]))[0]
            filt_var = var_filt[opt_filt.index(inter)]
            print(' >> Groups being analysed:',pair_all)#, '\n - inter:',inter, '- filt_var:', filt_var)
            data_out = filterR_Autom (data, [filt_var], inter)
            factor = list(set(filters) -set([filt_var]))[0]
            fact_levels = sorted(data[factor].unique())
            # print(' - factor:',factor, '- fact_levels:', fact_levels)
            
            # print('pair_f:',pair_f) 
            indiv_pairs = get_indiv_pairs(pair_f, inter)
            # print('indiv_pairs:', indiv_pairs)
            
            dict_spStatisticalRes_in = dict()
            print('\n > Group: ', inter)
            
            pvalues_norm, txt_normtest, all_data_normal, data_groups = normalityTests (alpha, factor, fact_levels, var, data_out, norm_test)
            test2use, multcomp_txt, txt_testSelected = selectStatisticalTest (fact_levels, all_data_normal, var)
            
            dict_spStatisticalRes_in['pvalues_norm'] = pvalues_norm
            dict_spStatisticalRes_in['txt_normtest'] = txt_normtest
            dict_spStatisticalRes_in['all_data_normal'] = all_data_normal
            dict_spStatisticalRes_in['test2use'] = test2use
            dict_spStatisticalRes_in['multcomp_txt'] = multcomp_txt
            
            if test2use == 'One-way-ANOVA':
                pval_spTest, pval_mcTest, pval_multComp = stTest_OWANOVA (factor, var, indiv_pairs, data_out, data_groups)
                dict_spStatisticalRes_in['pval_test'] = pval_spTest
                dict_spStatisticalRes_in['pval_multComp_o'] = pval_mcTest
                dict_spStatisticalRes_in['pval_multComp'] = pval_multComp
                dict_spStatisticalRes_in['stat_test'] = False
                dict_spStatisticalRes_in['test'] = []
                txt_testSelected = txt_testSelected+" (ANOVA's p-val: "+str(format(pval_spTest,'.3f'))+'), \n'+multcomp_txt
                
            elif test2use == 'Kruskal':
                pval_spTest, pval_mcTest, pval_multComp = stTest_Kruskal (factor, var, indiv_pairs, data_out, data_groups)
                dict_spStatisticalRes_in['pval_test'] = pval_spTest
                dict_spStatisticalRes_in['pval_multComp_o'] = pval_mcTest
                dict_spStatisticalRes_in['pval_multComp'] = pval_multComp
                dict_spStatisticalRes_in['stat_test'] = False
                dict_spStatisticalRes_in['test'] = []
                txt_testSelected = txt_testSelected+" (Dunn's 'p-val: "+str(format(pval_spTest,'.3f'))+'), \n'+multcomp_txt
            
            # Two variables
            elif test2use == 'Mann-Whitney': 
                pval_spTest, pval_multComp = stTest_Mann_Whitney(data_groups)
                dict_spStatisticalRes_in['pval_test'] = pval_spTest
                dict_spStatisticalRes_in['pval_multComp_o'] = pval_spTest
                dict_spStatisticalRes_in['pval_multComp'] = pval_multComp
                dict_spStatisticalRes_in['stat_test'] = False
                dict_spStatisticalRes_in['test'] = []
                txt_testSelected = txt_testSelected+' (U-statistic: '+str(format(pval_spTest,'.3f'))+'),'#+multcomp_txt
                list_of_lists = False
            
            elif test2use == 't-test_ind':
                pval_spTest, pval_multComp = stTest_t_test_indVar(data_groups)
                dict_spStatisticalRes_in['pval_test'] = pval_spTest
                dict_spStatisticalRes_in['pval_multComp_o'] = pval_spTest
                dict_spStatisticalRes_in['pval_multComp'] = pval_multComp
                dict_spStatisticalRes_in['stat_test'] = False
                dict_spStatisticalRes_in['test'] = []
                txt_testSelected = txt_testSelected+' (t-statistic: '+str(format(pval_spTest,'.3f'))+'),'#+multcomp_txt
                list_of_lists = False
                
            else: 
                pval_multComp = []
                dict_spStatisticalRes_in['stat_test'] = True
                dict_spStatisticalRes_in['test'] = test2use
                dict_spStatisticalRes_in['pval_test'] = None
                dict_spStatisticalRes_in['pval_multComp'] = pval_multComp
                
            dict_spStatisticalRes_in['txt_testSelected'] = txt_testSelected
            dict_spStatisticalRes[inter] = dict_spStatisticalRes_in
            # print(dict_spStatisticalRes_in)
            # input()
            pvals_mc.append(pval_multComp)
            txt_selected.append(txt_testSelected)
            txt_norm.append(txt_normtest)
            print('\n')
        
        if list_of_lists:
            pval_multComp_all = [item for sublist in pvals_mc for item in sublist]
        else: 
            pval_multComp_all = pvals_mc
            
        box_pairs = [item for sublist in box_pairs_f for item in sublist]
        # print('box_pairs:', box_pairs)
        # print('pval_multComp_all:', pval_multComp_all)
        
        dict_spStatisticalRes['box_pairs'] = box_pairs
        dict_spStatisticalRes['pval_multComp_all'] = pval_multComp_all
        dict_spStatisticalRes['txt_testSelected_all'] = txt_selected
        dict_spStatisticalRes['txt_normtest_all'] = txt_norm
        
        # input()
        
        dict_stats[var] = dict_spStatisticalRes
        
    return dict_stats

#%% func - stTest_OWANOVA
def stTest_OWANOVA (factor, var, indiv_pairs, data_out, data_groups):
    
    if len(data_groups) == 3:
        statOWA, pval_OWA = stats.f_oneway(data_groups[0], data_groups[1], data_groups[2])
        
    elif len(data_groups) == 4:
        statOWA, pval_OWA = stats.f_oneway(data_groups[0], data_groups[1], data_groups[2], data_groups[3])
    
    pval_Tukey_o = sp.posthoc_tukey(data_out, var, factor)
    print ("p-values Tukey's test:\n", pval_Tukey_o)
    
    pval_Tukey = []
    for pair in indiv_pairs:
        # print(pair[0], pair[1])
        pval = pval_Tukey_o.loc[pair[0],pair[1]]
        pval_Tukey.append(list(pval.iloc[0])[0])
    
    return pval_OWA, pval_Tukey_o, pval_Tukey

#%% func - stTest_Kruskal
def stTest_Kruskal (factor, var, indiv_pairs, data_out, data_groups): #(x: factor, y:variable, vals, data
    
    if len(data_groups) == 3:
        statKruskal, pval_Kruskal = stats.kruskal(data_groups[0], data_groups[1], data_groups[2])
        
    elif len(data_groups) == 4:
        statKruskal, pval_Kruskal = stats.kruskal(data_groups[0], data_groups[1], data_groups[2], data_groups[3])
        
    pval_Dunn_o = sp.posthoc_dunn(data_out, var, factor, 'bonferroni')
    print (" - p-values Dunn's test:\n", pval_Dunn_o)
    
    pval_Dunn = []
    for pair in indiv_pairs:
        # print(pair[0], pair[1])
        pval = pval_Dunn_o.loc[pair[0],pair[1]]
        pval_Dunn.append(list(pval.iloc[0])[0])
    
    return pval_Kruskal, pval_Dunn_o, pval_Dunn

#%% func - stTest_Mann_Whitney
def stTest_Mann_Whitney(data):
    
    u_stat, pval = stats.mannwhitneyu(data[0], data[1], alternative='two-sided')
    print (" -Mann-Whitney U-statistic: ", format(pval,'.3f'))
    
    return u_stat, pval

#%% fun - stTest_t_test_indVar
def stTest_t_test_indVar(data):
    
    stat, pval = stats.ttest_ind(a=data[0], b=data[1])
    print (" - t-statistic t-test independent samples: ", format(pval,'.3f'))
    
    return stat, pval
    
#%% func - def_box_pairs
def def_box_pairs(data, x_val, hue_val, btw_x = True, btw_hue = False):
    
    if x_val =='GenotypeAll':
        reverse_x = True
    else: 
        reverse_x = False
    x_values = sorted(data[x_val].unique(), reverse=reverse_x)
    
    if hue_val =='GenotypeAll':
        reverse_hue = True
    else: 
        reverse_hue = False
    hue_values = sorted(data[hue_val].unique(), reverse=reverse_hue)
    
    box_pairs_all = []
    if btw_x:
        box_pairs = []
        for x_val in x_values: 
            box_units = []
            for hue in hue_values: 
                box_units.append((x_val, hue))
            box_pairs.append(tuple(box_units))
        box_pairs_all.append(box_pairs)
    
    if btw_hue: 
        box_pairs = []
        for hue in hue_values: 
            box_units = []
            for x_val in x_values: 
                box_units.append((x_val, hue))
            box_pairs.append(tuple(box_units))
        box_pairs_all.append(box_pairs)
    # print(box_pairs_all)
    
    box_pairs_all_f = []
    for box in box_pairs_all[0]:
        box_pairs_f = []
        # print('box:',box)
        if len(box) > 2:
            boxes = list(combinations(box, 2))
            for b in boxes:
                box_pairs_f.append(b)
                # print('b:',b)
        else: #elif len(box) == 2: 
            box_pairs_f.append(box)
        box_pairs_all_f.append(box_pairs_f)
            
    return box_pairs_all[0], box_pairs_all_f

#%% func - get_indiv_pairs 
def get_indiv_pairs(pair_f, inter):
    
    indiv_pairs = []
    for pair in pair_f:
        # print(pair)
        in_pair = []
        for p in pair:
            p = list(p)
            # print(p)
            p.remove(inter)
            in_pair.append(p)
            # print('p_r:',p)
        indiv_pairs.append(in_pair)
        
    return indiv_pairs
    
#%% func - txtMultComp
def txtMultComp (box_pairs, pval, x_values, x_labels, test_selected, test_norm, alpha):
    txt_multcomp = '\n - - - - - - - - - - \n'
    num = int(len(box_pairs)/len(x_values))
    chunks = [box_pairs[x:x+num] for x in range(0, len(box_pairs), num)]
    
    n = 0
    for j, chunk, x_val, x_lab, test, norm in zip(count(), chunks, x_values, x_labels, test_selected, test_norm):
        txt_multcomp = txt_multcomp + x_lab + ': '+ norm +'\n'+ test + '\n - p-val: '
        for i, pair in enumerate(chunk):
            a_pair = list(pair[0]); b_pair = list(pair[1])
            a_pair.remove(x_val);b_pair.remove(x_val)
            pval_txt = str(format(pval[n],'.3f'))
            if pval[n] <= alpha:
                pval_txt = pval_txt +'*'
            pairtxt =  str(a_pair)+ ' vs. '+ str(b_pair) +': '+pval_txt
            n += 1
            if i < len(chunk)-1:
                txt_multcomp = txt_multcomp + pairtxt + ' \n '
            else: 
                txt_multcomp = txt_multcomp +  pairtxt + ' - '
        
        # print(j, len(x_values))
        # if j < len(x_values)-1:
        txt_multcomp = txt_multcomp+'\n - - - - - - - - - - \n'
            
    # print(txt_multcomp)
    
    # for i, pair in enumerate(box_pairs):
    #     pairtxt = str(pair[0])+ ' vs. '+ str(pair[1]) +': '+str(format(pval[i],'.3f'))
    #     if i < len(box_pairs)-1:
    #         txt_multcomp = txt_multcomp + pairtxt + ' \n '
    #     else: 
    #         txt_multcomp = txt_multcomp +  pairtxt + ' - '
    # print(txt_multcomp)
        
    return txt_multcomp

#%% - CREATE PLOTS
styles = ['o', '^', 's', 'v', 'D', '<', 'p', '>'] 
# styles = ['o', 'o', 's', 's', 'D', '<', 'p', '>'] # https://matplotlib.org/stable/api/markers_api.html
palettes = ['deeppink','royalblue','mediumturquoise', 'darkmagenta','darkorange','limegreen', 'gold']
# def plotInGroups(plot_type, input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
#                      h_plot, w_plot, save, dir2save, info, dpi = 300, sharey = False, h_add = 5, w_add = 1, ext = 'png'):

#%% func - plotInGroups
def plotInGroups (df2plot, vars2plot, x_var, hue_var, shape_var, title, labels2plot, ips, dir2save,
                  n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th,', info ='', 
                  save = True, dpi = 300, ext = 'png'):
    
    print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    sns.set_context('poster') # notebook, talk, poster, paper
    # Get legends
    dict_legends = def_legends(df2plot)
    
    # Set up the matplotlib figure
    num_vars = len(vars2plot)
    n_rows = math.ceil(num_vars/n_cols)
    # print('n_rows:' , n_rows)
    
    index_right_col = list(range(n_cols,(n_cols+1)*n_rows,4))
    index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
    index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))

    for index in index_no_plot:
        vars2plot.insert(index, '')
        labels2plot.insert(index, '')
    
    # Set up the matplotlib figure
    h_plot, w_plot = ips
    if num_vars == 1:
        h_add = 0; w_add = 0
    size_col = (n_cols+1)*h_plot+h_add
    size_row = n_rows*w_plot+w_add
    
    # Genotypes and Strains being plotted 
    values = []
    for var in [x_var, hue_var, shape_var]:
        if var =='GenotypeAll':
            reverse = True
        else: 
            reverse = False
        values.append(sorted(df2plot[var].unique(), reverse=reverse))
    x_values, hue_values, shape_values = values
    
    # - number of x_var
    n_x = len(x_values)
    
    #  Create figure  - plt.clf()
    gridkw = dict(width_ratios=[1]*n_cols+[0.2])
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
    sns.set_style("ticks")
    sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                    "ytick.labelsize":8, "lines.linewidth": 2.5,
                                                    "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
    
    # Define legends for shape and hue
    hue_legend = dict_legends[hue_var]
    legend_elem_hue = []
    for aa, hue_val in enumerate(hue_legend):
        legend_elem_hue.append(Line2D([0], [0], marker='h', color='w', label=hue_val,
                                markerfacecolor=palettes[aa], markersize=20))
    space = [Line2D([0], [0], marker='o', color='w', label='',
                                markerfacecolor='w', markersize=20)]
    shape_legend = dict_legends[shape_var]
    legend_elem_shape = []
    for n_str, shape_val, mark in zip(count(), shape_legend, styles):
        legend_elem_shape.append(Line2D([0], [0], marker=mark, color='w', label=shape_val,
                                markerfacecolor='k', markersize=15))
        
    handle_new = legend_elem_hue+space+legend_elem_shape
    legend_new = hue_legend+['']+shape_legend
    
    marker_size = 12
    dodge = True
    for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
        # Create list to contain info of yticks and ist highest value
        y_vals_all = []
        max_y_vals = []
    
        if n == 0: 
            for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
                print('\t- '+svar+': ', value)
                
        if n in index_no_plot:
            if n == n_cols:
                ax.set_axis_off()
                ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
            else: 
                ax.remove()
        else: 
            m = sns.swarmplot(data=df2plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
                              palette = palettes, dodge = dodge, size = marker_size)
            box = ax.get_position()
            if yticks_lab == '1e6 - d.':
                ylabel = ylabel +' x 10$^6$'
            elif yticks_lab == '1e3 - d.':
                ylabel = ylabel +' x 10$^3$'
            ax.set(xlabel=dict_legends['xlabels'][x_var], ylabel=ylabel);
            # ax.set_xticks(ax.get_xticks())
            ax.set_xticklabels(dict_legends[x_var], rotation=0)
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            ax.get_legend().remove()
            sns.despine()
            
            if n == 0:
                handles, labels = m.get_legend_handles_labels()

            y_vals = ax.get_yticks()
            y_vals_all.append(y_vals)
            max_y_vals.append(y_vals[-1])
            if n_x > 1: 
                vline_pos = list(range(1,n_x,1))
                for pos in vline_pos:
                    ax.axvline(pos - 0.5, ymax = 0.95, color='dimgrey', ls='-.', linewidth=0.8)
            
            # Define axes based on higherst bar
            max_y_index = max_y_vals.index(max(max_y_vals))
            ax.set_yticks(y_vals_all[max_y_index])
            if yticks_lab == '1e6 - d.':
                ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
            elif yticks_lab == '1e3 - d.':
                ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'th,':
                ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd. - 0':
                ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd. - 1':
                ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd.':
                ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])

    fig.suptitle(title+'\n', fontsize = 30, y=1)
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'pl_meas', 'R_')
            if info != '':
                fig_title = dir2savef+info+"_"+title+"."+extf
            else: 
                fig_title = dir2savef+title+"."+extf

            plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)

#%% func - plotInGroupsShape
def plotInGroupsShape (df2plot, vars2plot, x_var, hue_var, shape_var, title, labels2plot, ips, dir2save,
                  n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th,', ylim = '', info ='', 
                  save = True, dpi = 300, ext = 'png'):
    
    print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    sns.set_context('poster') # notebook, talk, poster, paper
    # Get legends
    dict_legends = def_legends(df2plot)
    
    # Set up the matplotlib figure
    num_vars = len(vars2plot)
    n_rows = math.ceil(num_vars/n_cols)
    # print('n_rows:' , n_rows)
    
    index_right_col = list(range(n_cols,(n_cols+1)*n_rows,4))
    index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
    index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))

    for index in index_no_plot:
        vars2plot.insert(index, '')
        labels2plot.insert(index, '')
    
    # Set up the matplotlib figure
    h_plot, w_plot = ips
    
    # Genotypes and Strains being plotted 
    values = []
    for var in [x_var, hue_var, shape_var]:
        if var =='GenotypeAll':
            reverse = True
        else: 
            reverse = False
        values.append(sorted(df2plot[var].unique(), reverse=reverse))
    x_values, hue_values, shape_values = values
    
    # - number of x_var
    n_x = len(x_values)
    if n_x == 1:
        h_plot = 3.5
    if num_vars == 1:
        h_add = 0; w_add = 0
    size_col = (n_cols+1)*h_plot+h_add
    size_row = n_rows*w_plot+w_add
    
    #  Create figure  - plt.clf()
    gridkw = dict(width_ratios=[1]*n_cols+[0.2])
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
    sns.set_style("ticks")
    sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                    "ytick.labelsize":8, "lines.linewidth": 2.5,
                                                    "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
    
    # Define legends for shape and hue
    hue_legend = dict_legends[hue_var]
    legend_elem_hue = []
    for aa, hue_val in enumerate(hue_legend):
        legend_elem_hue.append(Line2D([0], [0], marker='h', color='w', label=hue_val,
                                markerfacecolor=palettes[aa], markersize=20))
    space = [Line2D([0], [0], marker='o', color='w', label='',
                                markerfacecolor='w', markersize=20)]
    shape_legend = dict_legends[shape_var]
    legend_elem_shape = []
    for n_str, shape_val, mark in zip(count(), shape_legend, styles):
        legend_elem_shape.append(Line2D([0], [0], marker=mark, color='w', label=shape_val,
                                markerfacecolor='k', markersize=15))
        
    handle_new = legend_elem_hue+space+legend_elem_shape
    legend_new = hue_legend+['']+shape_legend
    
    marker_size = 10
    dodge = True
    jitter = 0.3
    
    plot_no = 0
    for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
        # Create list to contain info of yticks and ist highest value
        y_vals_all = []
        max_y_vals = []
    
        if n == 0: 
            for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
                print('\t- '+svar+': ', value)
                
        if n in index_no_plot:
            if n == n_cols:
                ax.set_axis_off()
                ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
            else: 
                ax.remove()
        else: 
            for j, val, style in zip(count(), shape_values, styles):
                df_plot = df2plot[df2plot[shape_var] == val]
                # print(x_var, var, hue_var, hue_values, x_values, style, palettes)
                m = sns.stripplot(data=df_plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
                                  marker = style, palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
                box = ax.get_position()
                if yticks_lab == '1e6 - d.':
                    ylabel = ylabel +' x 10$^6$'
                elif yticks_lab == '1e3 - d.':
                    ylabel = ylabel +' x 10$^3$'
                ax.set(xlabel=dict_legends['xlabels'][x_var], ylabel=ylabel);
                # ax.set_xticks(ax.get_xticks())
                ax.set_xticklabels(dict_legends[x_var], rotation=0)
                ax.set_position([box.x0, box.y0, box.width*1, box.height])
                ax.get_legend().remove()
                sns.despine()
                
                if ylim != '':
                    # print('ylim:', ylim[plot_no][0],'-', ylim[plot_no][1])
                    ax.set_ylim(ylim[plot_no][0], ylim[plot_no][1])
                
                if n == 0:
                    handles, labels = m.get_legend_handles_labels()

                y_vals = ax.get_yticks()
                y_vals_all.append(y_vals)
                max_y_vals.append(y_vals[-1])
                if n_x > 1: 
                    vline_pos = list(range(1,n_x,1))
                    for pos in vline_pos:
                        ax.axvline(pos - 0.5, ymax = 0.95, color='dimgrey', ls='-.', linewidth=0.8)
            
            # Define axes based on higherst bar
            max_y_index = max_y_vals.index(max(max_y_vals))
            ax.set_yticks(y_vals_all[max_y_index])
            if yticks_lab == '1e6 - d.':
                ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
            elif yticks_lab == '1e3 - d.':
                ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'th,':
                ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd. - 0':
                ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd. - 1':
                ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd.':
                ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])
                
        plot_no +=1

    fig.suptitle(title+'\n', fontsize = 30, y=1)
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'pl_meas', 'R_')
            if info != '':
                fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
            else: 
                fig_title = dir2savef+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf

            plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
    
#%% func - plotInGroupsStats
def plotInGroupsStats(df2plot, vars2plot, x_var, hue_var, shape_var, title, labels2plot, ips, dir2save, stats_set,
                      n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th,', ylim = '', info ='', 
                      save = True, dpi = 300, ext = 'png'):
    
    if stats_set[0]: 
        tests_res = []
        box_pairs_all = []
        stats = True
        dict_stats = stats_set[1]
        alpha = stats_set[2]
        txt_multcomp_all = []
        
    print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    sns.set_context('poster') # notebook, talk, poster, paper
    # Get legends
    dict_legends = def_legends(df2plot)
    
    # Set up the matplotlib figure
    num_vars = len(vars2plot)
    n_rows = math.ceil(num_vars/n_cols)+1
    # print('n_rows:' , n_rows)
    # print('n_cols:' , n_cols)
    
    index_right_col = list(range(n_cols,(n_cols+1)*n_rows,n_cols+1))
    index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
    index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))
    index_stats = sorted(list(set(index_no_graph).difference(set(index_right_col))))
    # print(index_stats)

    for index in index_no_plot:
        vars2plot.insert(index, '')
        labels2plot.insert(index, '')
    
    # Set up the matplotlib figure
    h_plot, w_plot = ips
    
    # Genotypes and Strains being plotted 
    values = []
    for var in [x_var, hue_var, shape_var]:
        if var =='GenotypeAll':
            reverse = True
        else: 
            reverse = False
        values.append(sorted(df2plot[var].unique(), reverse=reverse))
    x_values, hue_values, shape_values = values
    
    # - number of x_var
    n_x = len(x_values)
    if n_x == 1:
        h_plot = 6
        hspace = 0.2; wspace = 1
    else: 
        hspace = 0.5; wspace = 0.5
        
    if num_vars == 1:
        h_add = 0; w_add = 0
    size_col = (n_cols+1)*h_plot+h_add
    size_row = n_rows*w_plot+w_add
    
    #  Create figure  - plt.clf()
    gridkw = dict(width_ratios=[1]*n_cols+[0.2], height_ratios = [1,1.2])
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
    fig.subplots_adjust(hspace=hspace, wspace=wspace)
    
    sns.set_style("ticks")
    sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                    "ytick.labelsize":8, "lines.linewidth": 2.5,
                                                    "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
    
    # Define legends for shape and hue
    hue_legend = dict_legends[hue_var]
    legend_elem_hue = []
    for aa, hue_val in enumerate(hue_legend):
        legend_elem_hue.append(Line2D([0], [0], marker='h', color='w', label=hue_val,
                                markerfacecolor=palettes[aa], markersize=20))
    space = [Line2D([0], [0], marker='o', color='w', label='',
                                markerfacecolor='w', markersize=20)]
    shape_legend = dict_legends[shape_var]
    legend_elem_shape = []
    for n_str, shape_val, mark in zip(count(), shape_legend, styles):
        legend_elem_shape.append(Line2D([0], [0], marker=mark, color='w', label=shape_val,
                                markerfacecolor='k', markersize=15))
        
    handle_new = legend_elem_hue+space+legend_elem_shape
    legend_new = hue_legend+['']+shape_legend
    
    marker_size = 8
    dodge = True
    plot_no = 0
    for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
        
        # Create list to contain info of yticks and ist highest value
        y_vals_all = []
        max_y_vals = []
    
        if n == 0: 
            for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
                print('\t- '+svar+': ', value)
        
        if n in index_no_plot:
            if n == n_cols:
                ax.set_axis_off()
                ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
            else: 
                if n in index_stats:
                    text2add = txt_multcomp_all[n-n_cols-1]
                    ax.text(0.5, 0.5, text2add, size=20, ha='center', va='center')
                    # ax.remove()
                    ax.set_axis_off()
                    # ax.get_xaxis().set_visible(False)
                    # ax.get_yaxis().set_visible(False)
                else: 
                    ax.remove()
                
        else: 
            print(' > ',var,'\n')
            m = sns.stripplot(data=df2plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
                              palette = palettes, dodge = dodge, size = marker_size, jitter = 0.2)
            # m = sns.swarmplot(data=df2plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
            #                   palette = palettes, dodge = dodge, size = marker_size)
            if stats: 
                box_pairs = dict_stats[var]['box_pairs']
                stat_test = False
                test = None
                p_val = dict_stats[var]['pval_multComp_all']
                txt_testSelected_all = dict_stats[var]['txt_testSelected_all']
                txt_normtest_all = dict_stats[var]['txt_normtest_all']
                
                m1, test_results = add_stat_annotation(ax, data=df2plot, x=x_var, y=var, hue=hue_var, 
                                                       hue_order = hue_values, order = x_values,
                                                        box_pairs=box_pairs, perform_stat_test = stat_test, test = test,
                                                        pvalues=p_val, comparisons_correction=None, #'bonferroni',
                                                        line_offset_to_box=0.4, line_offset=0.1,
                                                        line_height=0.015, text_offset=5,
                                                        text_format='star', loc='inside', fontsize='small', verbose=0);
                tests_res.append(p_val)
                box_pairs_all.append(box_pairs)
                txt_multcomp = txtMultComp (box_pairs, p_val, x_values, dict_legends[x_var], txt_testSelected_all, txt_normtest_all, alpha)
                txt_multcomp_all.append(txt_multcomp)
                # print('txt_multcomp:', txt_multcomp)
                
            box = ax.get_position()
            if yticks_lab == '1e6 - d.':
                ylabel = ylabel +' x 10$^6$'
            elif yticks_lab == '1e3 - d.':
                ylabel = ylabel +' x 10$^3$'
            ax.set(xlabel=dict_legends['xlabels'][x_var], ylabel=ylabel);
            # ax.set_xticks(ax.get_xticks())
            ax.set_xticklabels(dict_legends[x_var], rotation=0)
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            ax.get_legend().remove()
            sns.despine()
            if ylim != '':
                ax.set_ylim(ylim[plot_no][0], ylim[plot_no][1])
            
            if n == 0:
                handles, labels = m.get_legend_handles_labels()

            y_vals = ax.get_yticks()
            y_vals_all.append(y_vals)
            max_y_vals.append(y_vals[-1])
            if n_x > 1: 
                vline_pos = list(range(1,n_x,1))
                for pos in vline_pos:
                    ax.axvline(pos - 0.5, ymax = 0.95, color='dimgrey', ls='-.', linewidth=0.8)
            
            # Define axes based on higherst bar
            max_y_index = max_y_vals.index(max(max_y_vals))
            ax.set_yticks(y_vals_all[max_y_index])
            if yticks_lab == '1e6 - d.':
                ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
            elif yticks_lab == '1e3 - d.':
                ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'th,':
                ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd. - 0':
                ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd. - 1':
                ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd.':
                ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])
            
            plot_no +=1

    fig.suptitle(title+'\n', fontsize = 30, y=0.95)
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'pl_stats', 'R_')
            if info != '':
                fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
            else: 
                fig_title = dir2savef+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf

            plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
    
    if stats: 
        return tests_res, box_pairs_all, x_values, dict_legends[x_var]
    
#%% func - plotInGroupsStatsBU
def plotInGroupsStatsBU(df2plot, vars2plot, x_var, hue_var, shape_var, title, labels2plot, ips, dir2save, stats_set,
                      n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th,', info ='', 
                      save = True, dpi = 300, ext = 'png'):
    
    if stats_set[0]: 
        tests_res = []
        stats = True
        dict_stats = stats_set[1]
        
    print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    sns.set_context('poster') # notebook, talk, poster, paper
    # Get legends
    dict_legends = def_legends(df2plot)
    
    # Set up the matplotlib figure
    num_vars = len(vars2plot)
    n_rows = math.ceil(num_vars/n_cols)+1
    print('n_rows:' , n_rows)
    
    index_right_col = list(range(n_cols,(n_cols+1)*n_rows,4))
    index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
    index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))
    print(index_no_plot)

    for index in index_no_plot:
        vars2plot.insert(index, '')
        labels2plot.insert(index, '')
    
    # Set up the matplotlib figure
    h_plot, w_plot = ips
    if num_vars == 1:
        h_add = 0; w_add = 0
    size_col = (n_cols+1)*h_plot+h_add
    size_row = n_rows*w_plot+w_add
    
    # Genotypes and Strains being plotted 
    values = []
    for var in [x_var, hue_var, shape_var]:
        if var =='GenotypeAll':
            reverse = True
        else: 
            reverse = False
        values.append(sorted(df2plot[var].unique(), reverse=reverse))
    x_values, hue_values, shape_values = values
    
    # - number of x_var
    n_x = len(x_values)
    
    #  Create figure  - plt.clf()
    gridkw = dict(width_ratios=[1]*n_cols+[0.2], height_ratios = [1,1])
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
    sns.set_style("ticks")
    sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                    "ytick.labelsize":8, "lines.linewidth": 2.5,
                                                    "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
    
    # Define legends for shape and hue
    hue_legend = dict_legends[hue_var]
    legend_elem_hue = []
    for aa, hue_val in enumerate(hue_legend):
        legend_elem_hue.append(Line2D([0], [0], marker='h', color='w', label=hue_val,
                                markerfacecolor=palettes[aa], markersize=20))
    space = [Line2D([0], [0], marker='o', color='w', label='',
                                markerfacecolor='w', markersize=20)]
    shape_legend = dict_legends[shape_var]
    legend_elem_shape = []
    for n_str, shape_val, mark in zip(count(), shape_legend, styles):
        legend_elem_shape.append(Line2D([0], [0], marker=mark, color='w', label=shape_val,
                                markerfacecolor='k', markersize=15))
        
    handle_new = legend_elem_hue+space+legend_elem_shape
    legend_new = hue_legend+['']+shape_legend
    
    marker_size = 10
    dodge = True
    for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
        
        # Create list to contain info of yticks and ist highest value
        y_vals_all = []
        max_y_vals = []
    
        if n == 0: 
            for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
                print('\t- '+svar+': ', value)
        
        if n in index_no_plot:
            if n == n_cols:
                ax.set_axis_off()
                ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
            else: 
                ax.remove()
        else: 
            print(' > ',var,'\n')
            m = sns.swarmplot(data=df2plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
                              palette = palettes, dodge = dodge, size = marker_size)
            if stats: 
                box_pairs = dict_stats[var]['box_pairs']
                stat_test = False
                test = None
                p_val = dict_stats[var]['pval_multComp_all']
                
                m1, test_results = add_stat_annotation(ax, data=df2plot, x=x_var, y=var, hue=hue_var, 
                                                       hue_order = hue_values, order = x_values,
                                                        box_pairs=box_pairs, perform_stat_test = stat_test, test = test,
                                                        pvalues=p_val, comparisons_correction=None, #'bonferroni',
                                                        line_offset_to_box=0.4, line_offset=0.1,
                                                        line_height=0.015, text_offset=5,
                                                        text_format='star', loc='inside', verbose=0);
                tests_res.append(test_results)
                # txt_multcomp = txtMultComp(box_pairs, test_results)
                # print('txt_multcomp:', txt_multcomp)
                
            box = ax.get_position()
            if yticks_lab == '1e6 - d.':
                ylabel = ylabel +' x 10$^6$'
            elif yticks_lab == '1e3 - d.':
                ylabel = ylabel +' x 10$^3$'
            ax.set(xlabel=dict_legends['xlabels'][x_var], ylabel=ylabel);
            # ax.set_xticks(ax.get_xticks())
            ax.set_xticklabels(dict_legends[x_var], rotation=0)
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            ax.get_legend().remove()
            sns.despine()
            
            if n == 0:
                handles, labels = m.get_legend_handles_labels()

            y_vals = ax.get_yticks()
            y_vals_all.append(y_vals)
            max_y_vals.append(y_vals[-1])
            if n_x > 1: 
                vline_pos = list(range(1,n_x,1))
                for pos in vline_pos:
                    ax.axvline(pos - 0.5, ymax = 0.95, color='dimgrey', ls='-.', linewidth=0.8)
            
            # Define axes based on higherst bar
            max_y_index = max_y_vals.index(max(max_y_vals))
            ax.set_yticks(y_vals_all[max_y_index])
            if yticks_lab == '1e6 - d.':
                ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
            elif yticks_lab == '1e3 - d.':
                ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'th,':
                ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd. - 0':
                ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd. - 1':
                ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd.':
                ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])

    fig.suptitle(title+'\n', fontsize = 30, y=1)
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'pl_meas', 'R_')
            if info != '':
                fig_title = dir2savef+info+"_"+title+"."+extf
            else: 
                fig_title = dir2savef+title+"."+extf

            plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
    
    if stats: 
        return tests_res

#%% func - relPlotInGroups
def relPlotInGroups (df2plot, vars2plot, x_var, hue_var, shape_var, title, labels2plot, ips, dir2save,
                     n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th.', info ='', 
                     save = True, dpi = 300, ext = 'png'):
    
    print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    sns.set_context('poster') # notebook, talk, poster, paper
    # Get legends
    dict_legends = def_legends(df2plot)
    
    # Set up the matplotlib figure
    num_vars = len(vars2plot)
    n_rows = math.ceil(num_vars/n_cols)
    # print('n_rows:' , n_rows)
    
    index_right_col = list(range(n_cols,(n_cols+1)*n_rows,4))
    index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
    index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))

    for index in index_no_plot:
        vars2plot.insert(index, '')
        labels2plot.insert(index, '')
    
    # Set up the matplotlib figure
    h_plot, w_plot = ips
    if num_vars == 1:
        h_add = 0; w_add = 0
    size_col = (n_cols+1)*h_plot+h_add
    size_row = n_rows*w_plot+w_add
    
    # Genotypes and Strains being plotted 
    values = []
    for var in [x_var, hue_var, shape_var]:
        if var =='GenotypeAll':
            reverse = True
        else: 
            reverse = False
        values.append(sorted(df2plot[var].unique(), reverse=reverse))
    x_values, hue_values, shape_values = values
    
    # - number of x_var
    # n_x = len(x_values)
    
    #  Create figure  - plt.clf()
    gridkw = dict(width_ratios=[1]*n_cols+[0.2])
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
    sns.set_style("ticks")
    sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                    "ytick.labelsize":8, "lines.linewidth": 2.5,
                                                    "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
    
    # Define legends for shape and hue
    hue_legend = dict_legends[hue_var]
    legend_elem_hue = []
    for aa, hue_val in enumerate(hue_legend):
        legend_elem_hue.append(Line2D([0], [0], marker='h', color='w', label=hue_val,
                                markerfacecolor=palettes[aa], markersize=20))
    space = [Line2D([0], [0], marker='o', color='w', label='',
                                markerfacecolor='w', markersize=20)]
    shape_legend = dict_legends[shape_var]
    legend_elem_shape = []
    for n_str, shape_val, mark in zip(count(), shape_legend, styles):
        legend_elem_shape.append(Line2D([0], [0], marker=mark, color='w', label=shape_val,
                                markerfacecolor='k', markersize=15))
        
    handle_new = legend_elem_hue+space+legend_elem_shape
    legend_new = hue_legend+['']+shape_legend
    
    # marker_size = 12
    # dodge = True
    # jitter = 0.3
    for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
        # Create list to contain info of yticks and ist highest value
        y_vals_all = []
        max_y_vals = []
    
        if n == 0: 
            for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
                print('\t- '+svar+': ', value)

        if n in index_no_plot:
            if n == n_cols:
                ax.set_axis_off()
                ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
            else: 
                ax.remove()
        else: 
            m = sns.lineplot(x=x_var, y=var,estimator=None, data= df2plot, 
                             palette = 'Set2', ax = ax)
            box = ax.get_position()
            if yticks_lab == '1e6 - d.':
                ylabel = ylabel +' x 10$^6$'
            elif yticks_lab == '1e3 - d.':
                ylabel = ylabel +' x 10$^3$'
            # print(dict_legends['xlabels'][x_var], ylabel)
            
            ax.set(xlabel=dict_legends['xlabels'][x_var], ylabel=ylabel);
            ax.set_xticks(ax.get_xticks())
            ax.set_xticklabels(dict_legends[x_var], rotation=0)
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            sns.despine()
            
            if n == 0:
                handles, labels = m.get_legend_handles_labels()

            y_vals = ax.get_yticks()
            y_vals_all.append(y_vals)
            max_y_vals.append(y_vals[-1])
            
            # Define axes based on higherst bar
            max_y_index = max_y_vals.index(max(max_y_vals))
            ax.set_yticks(y_vals_all[max_y_index])
            if yticks_lab == '1e6 - d.':
                ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
            elif yticks_lab == '1e3 - d.':
                ax.set_yticklabels(['{:.2f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'th,':
                ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
            elif yticks_lab == 'd.':
                ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])

    fig.suptitle(title+'\n', fontsize = 30, y=1)
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'pl_meas', 'R_')
            if info != '':
                fig_title = dir2savef+info+"_"+title+"."+extf
            else: 
                fig_title = dir2savef+title+"."+extf

            plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)

#%% func - plotPerVariableLabels
def plotPerVariableLabels(script, input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
                     h_plot, w_plot, save, dir2save, info, dpi = 300, h_add = 5, w_add = 1):
    
    styles = ['o', '^', 's', 'v', 'D', '<', 'p', '>'] # https://matplotlib.org/stable/api/markers_api.html
    variables, ylabels = def_variables(script)
    
    for i, input_var, title in zip(count(), input_vars, titles):
        vars2plot, labels2plot = selectVariables_auto(variables, ylabels, input_var)
        
        # Set up the matplotlib figure
        stages = sorted(df2plot.Stage.unique())
        plots_per_col = len(stages)
        plots_per_row = 1
        stages.append('')
        
        # Set up the matplotlib figure
        size_col = (plots_per_col+1)*h_plot-12
        size_row = plots_per_row*w_plot+w_add
        
        # Genotypes and Strains being plotted 
        genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
        strains = sorted(df2plot.Strain.unique())

        index_no_plot = list(range(3,(plots_per_col+1)*plots_per_row,4))

        if i == 0: 
            print('- Genotypes: ', genots)
            print('- Strains: ', strains)
            print('- Stages: ', stages)
        
        palettes = ['mediumturquoise', 'darkmagenta']
        
        for nn, var, ylabel in zip(count(), vars2plot, labels2plot):
            #  Create figure  - plt.clf()
            gridkw = dict(width_ratios=[1,1,1,0.2])
            fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
            fig.subplots_adjust(hspace=1.5, wspace=0.05)
            
            sns.set_style("ticks")
            sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                            "ytick.labelsize":'small', "lines.linewidth": 2.5,
                                                            "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
            # Define legends for strain and genotype
            legend_elem_gen = []
            for gen, gen_lab in enumerate(gen_legend):
                legend_elem_gen.append(Line2D([0], [0], marker='h', color='w', label=gen_lab,
                                        markerfacecolor=palettes[gen], markersize=10))
            handle_new = legend_elem_gen
            legend_new = gen_legend
            
            marker_size = 10
            dodge = False
            jitter = 0.2
            for n, ax, stg in zip(count(), axes.flatten(), stages):
                # print(n)
                if n in index_no_plot:
                    if n == 3:
                        ax.set_axis_off()
                        ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(0, 1), frameon = False)
                    else: 
                        ax.remove()
                else: 
                    df_plot = df2plot[df2plot['Stage'] == stg]
                    m = sns.stripplot(x="Strain", y=var, hue="GenotypeAll", hue_order = genots, data=df_plot, ax = ax, order=strains,
                                  marker = styles[0], palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
                    box = ax.get_position()
                    m.set_xticklabels(strain_legend, rotation=50)
                    ax.set(xlabel="\nStage: "+stg+'hpf')
                    ymin, ymax = ax.get_ylim()
                    up = (ymax-ymin)//20
                    # print(up)
                    df_plot = df_plot.sort_values(by=['Stage','Strain', var])
                    fish_refs = df_plot.Ref.unique()
                    ha = ['right', 'left', 'center']*50
                    for j, ref in enumerate(fish_refs):
                        strain_pos =strains.index(df_plot['Strain'].values[j])
                        gen_pos = genots.index(df_plot['GenotypeAll'].values[j])
                        # print(ref, ha[j])
                        # rand_x = random.uniform(0.95, 1.05)
                        # rand_y = random.uniform(0.95, 1.05)
                        # ax.text(x=strain_pos*rand_x, y=df_plot[var].values[j]*rand_y, s=ref, horizontalalignment=ha[j], size=10, color=palettes[gen_pos])
                        ax.text(x=strain_pos, y=df_plot[var].values[j]+up, s=ref, horizontalalignment=ha[j], size=10, color=palettes[gen_pos])
    
                    if n == 0:
                        ax.set(ylabel='\n'+ylabel+'\n')
                    else: 
                        ax.set(ylabel = '')
                        ax.tick_params(axis = 'y', labelcolor='w', width=0.1)
                        ax.spines['left'].set_color('gray')
                        ax.spines['left'].set_linestyle((0,(4,4)))#"dashed")
                    ax.set_position([box.x0, box.y0, box.width*1, box.height])
                    ax.get_legend().remove()
                    
                    # ax.legend(handles0[:],labels0[:],loc='center left', bbox_to_anchor=(1, 0.5))
                    sns.despine()
                    
                    if n == 0:
                        handles, labels = m.get_legend_handles_labels()
                        # print(handles)
                        print(var, '- Genotypes: ',labels)
    
            fig.suptitle(title, fontsize = 30, y=1)
            dir2savef = os.path.join(dir2save, 'pl_labels', 'R_')
            if info != '':
                fig_title = dir2savef+"Lab_"+info+"_"+var+".png"
            else: 
                fig_title = dir2savef+"Lab_"+var+".png"
            
            if save: 
                plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)

#%% func - barPlots
def barPlots(df2plot, vars2plot, group_vars, x_var, col_var, colours, title, txt_title, ylabel, 
             dir2save, yticks_lab, info='', stack100 = True, sub_bar_lab = True, save = True, ext = ['png']): 
    
    mpl.rcParams.update(mpl.rcParamsDefault)
    print('\n> '+title)
    print('\t - x_var: '+x_var+ ' - col_var:'+col_var)
    # Create list to contain info of bars and ist highest value
    y_vals_all = []
    max_y_vals = []
    
    # Filter dataframe 
    df4plot = df2plot.groupby(group_vars)[vars2plot].mean()
    df4plot_count = df2plot.groupby(group_vars)[vars2plot].count()
    dict_legend = def_legends(df2plot.reset_index())
    
    # Count and define figure settings accordingly
    # - number of columns define legend position 
    # Define legends for columns 
    reverse = False
    ascending = False
    if col_var == 'GenotypeAll':
        reverse = True
        ascending = True
        col_legend = dict_legend[col_var]#gen_legend
    elif col_var == 'Stage':
        ascending = False
        stage_legend = dict_legend[col_var]
        col_legend = ['Stage: ' + s for s in stage_legend]#'Stage:'+stage_legend
    else: 
        reverse = False
        ascending = False
        print('HELP A')
    
    col_values = sorted(df2plot[col_var].unique(), reverse = reverse)
    n_col = len(col_values)
    print(n_col)
    for nn in range(n_col): col_values.append('')
    
    # - number of variables/ stacked bars define number of columns in legend
    if n_col >= 3:
        leg_pos = 4
        loc='center'
        bbox_to_anchor=(0.5, 0.5)
        ncols_leg = len(vars2plot)
    else:
        leg_pos = n_col
        loc='center'
        if len(vars2plot) == 3:
            ncols_leg = 3
        else: 
            ncols_leg = 2
        if n_col == 1:
            bbox_to_anchor=(0.5, 0.5)
        else: 
            bbox_to_anchor=(1, 0.5)
            
    # - number of x_var
    n_x = len(df2plot[x_var].unique())
    
    
    # If bar plots are stacked to 100 transform data to percentages
    if stack100: 
        df4plot = df4plot.apply(lambda zz: zz*100/sum(zz), axis=1)

    # Find minimum bar height to later define offset  
    min_val = df4plot.min().min()
    
    #Create figure
    gridkw = dict(width_ratios=[1]*n_col, height_ratios=[1,0.3])
    fig_size = (math.ceil(2.5*(n_x/2)*(n_col)),7)
    fig, axes = plt.subplots(ncols=n_col,nrows=2, figsize=fig_size, sharey=True, gridspec_kw=gridkw)
    #print('Figure size: ', fig_size)
    for n, ax, col_val in zip(count(), axes.flatten(), col_values):
        # First row will contain the graphs
        if n < n_col:
            df_col = df4plot.iloc[df4plot.index.get_level_values(col_var) == col_val]
            df_col = df_col.droplevel(col_var, axis="index").sort_index(key=lambda xx: xx.str.lower(), ascending=ascending)
            df_col_count = df4plot_count.iloc[df4plot.index.get_level_values(col_var) == col_val]
            df_col_count = df_col_count.mean(axis=1).droplevel(col_var, axis="index").sort_index(key=lambda yy: yy.str.lower(), ascending=ascending)
            
            if df_col.index.nlevels > 1:
                joined_titles = df_col.index.map(('\n'.join))
                df2plot_titles = [jt+'\nn='+str(int(num)) for i, jt, num in zip(count(), joined_titles, df_col_count)]
            else: 
                df_col_copy = df_col.copy()
                dict_legend = def_legends(df_col_copy.reset_index())
                df2plot_titles = [leg_val+'\nn='+str(int(num)) for i, leg_val, num in zip(count(), dict_legend[x_var], df_col_count)]
               
            print('\n - '+col_val +'\n'); print(df_col)#, df2plot_titles)
            # Initialize the bottom at zero for the first set of bars.
            bottom = np.zeros(len(df_col))
            
            # Plot each layer of the bar, adding each bar to the "bottom" so
            # the next bar starts higher.
            for i, colm in enumerate(df_col.columns):
                ax.bar(df2plot_titles, df_col[colm], bottom=bottom, label=colm, color=colours[i], width = 0.7)
                bottom += np.array(df_col[colm])
                ax.set_xticks(df2plot_titles)
                ax.set_xticklabels(df2plot_titles, rotation=30)
            y_vals = ax.get_yticks()
            y_vals_all.append(y_vals)
            max_y_vals.append(y_vals[-1])
            
            if n == 0:
                handles, labels = ax.get_legend_handles_labels()
                if not stack100 and yticks_lab == '1e3 - d.':
                    ylabel = ylabel +' x 10$^3$'
                ax.set(ylabel=ylabel)
                
            # Sum up the rows of our data to get the total value of each bar.
            totals = df_col.sum(axis=1)
            if not stack100:
                # Set an offset that is used to bump the label up a bit above the bar.
                y_offset = round(0.2*min_val)
                # Add labels to each bar.
                for i, total in enumerate(totals):
                    if total < 100:
                          ax.text(i, total + y_offset, 
                                  ['{:.1f}'.format(w) for w in total], 
                                  ha='center', weight='bold', size =10)
                    else: 
                          # ax.text(i, total + y_offset ,locale.format("%d", total, grouping=True), ha='center', weight='bold', size =9)
                          if yticks_lab == '1e3 - d.':
                              ax.text(i, total + y_offset ,
                                      '{:.1f}'.format(total/1e3),
                                      ha='center', weight='bold', size =9)
                          elif yticks_lab == 'th,':
                              ax.text(i, total + y_offset ,
                                      locale.format("%d", total, grouping=True),
                                      ha='center', weight='bold', size =9)
                              
            # Let's put the annotations inside the bars themselves by using a
            # negative offset.
            if sub_bar_lab: 
                if stack100:
                    y_offset = -0.5
                else: 
                    y_offset = round(-0.4*min_val)
                # For each patch (basically each rectangle within the bar), add a label.
                for bar in ax.patches:
                    if bar.get_height() < 100:
                        bar_height = "{0:.1f}".format(bar.get_height())
                    else: 
                        if yticks_lab == '1e3 - d.':
                            bar_height = '{:.1f}'.format(bar.get_height()/1e3)
                        elif yticks_lab == 'th,':
                            bar_height = locale.format("%d",  bar.get_height(), grouping=True)
                    ax.text(
                          # Put the text in the middle of each bar. get_x returns the start
                          # so we add half the width to get to the middle.
                          bar.get_x() + bar.get_width() / 2,
                          # Vertically, add the height of the bar to the start of the bar,
                          # along with the offset.
                          bar.get_height()//2 + bar.get_y()+ y_offset,
                          # This is actual value we'll show. original: round(bar.get_height())
                          bar_height,
                          # Center the labels and style them a bit.
                          ha='center', color='w', weight='bold', size=9)
                  
            ax.set_title(col_legend[n])
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
        
        # Second row will contain the legend
        if n >= n_col:
            ax.set_axis_off()
            ax.set_title('')
            if n== leg_pos:
                legends = def_var_names('bar_plots', labels)
                ax.legend(handles, legends, loc=loc, bbox_to_anchor = bbox_to_anchor, 
                          frameon = True, fontsize=10, ncol = ncols_leg)
    
    # Define axes based on higherst bar
    if not stack100:
        max_y_index = max_y_vals.index(max(max_y_vals))
        ax.set_yticks(y_vals_all[max_y_index])
        # ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
        if yticks_lab == '1e3 - d.':
            ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
        elif yticks_lab == 'th,':
            ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
        # ax.ticklabel_format(axis='y', style='plain')

    fig.suptitle(title+txt_title+'\n', fontsize = 11, y=1.01)
    
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'pl_bars', 'R_')
            title2save = title.replace(' per Stage and Genotype','')
            title2save = title2save.replace(' per Genotype and Stage','')
            title2save = title2save.replace(' ','_')
            if stack100: 
                title2save = 'Perc'+title2save
            if info != '':
                fig_title = dir2savef+info+"_"+title2save+"."+extf
            else: 
                fig_title = dir2savef+title2save+"."+extf
            plt.savefig(fig_title, dpi=300, bbox_inches='tight', transparent=True)

#%% func - pctChange_barPlots
def pctChange_barPlots(df2plot, vars2plot, group_vars, x_var, col_var, colours, title, txt_title, ylabel, 
                       dir2save, info='', sub_bar_lab = True, save = True, ext = ['png']): 
    
    mpl.rcParams.update(mpl.rcParamsDefault)
    print('\n> '+title.replace('\n',' '))
    print('\t - x_var: '+x_var+ ' - col_var:'+col_var)
    # Create list to contain info of bars and ist highest value
    y_vals_all = []
    max_y_vals = []
    
    # Filter dataframe 
    df4plot = df2plot.groupby(group_vars)[vars2plot].mean()
    dict_legend = def_legends(df2plot.reset_index(), df_type = 'changes')
    # print(dict_legend)

    # Count and define figure settings accordingly
    # - number of x_var
    x_values = sorted(df2plot.index.unique(level = x_var))
    n_x = len(x_values)
        
    # Define legends
    col_legend = dict_legend[col_var]
    reverse = False
    if col_var == 'GenotypeAll':
        reverse = True
        
    # - number of columns define legend position 
    col_values = sorted(df2plot.index.unique(level = col_var), reverse = reverse)
    n_col = len(df2plot.index.unique(level = col_var))
    for nn in range(n_col): col_values.append('')
    
    # - number of variables/ stacked bars define number of columns in legend
    if n_col >= 2:
        leg_pos = 2
        loc='center'
        bbox_to_anchor=(1, 0.5)
        ncols_leg = len(vars2plot)
    else: 
        leg_pos = n_col
        loc='center'
        bbox_to_anchor=(0.5, 0.5)
        if len(vars2plot) == 3:
            ncols_leg = 3
        else: 
            ncols_leg = 2
        # if n_col == 1:
        #     bbox_to_anchor=(0.5, 0.5)
        # else: 
        #     bbox_to_anchor=(1, 0.5)
    
    #Create figure
    gridkw = dict(width_ratios=[1]*n_col, height_ratios=[1,0.3])
    fig_size = (math.ceil(2.5*(n_x/2)*(n_col)),7)
    fig, axes = plt.subplots(ncols=n_col,nrows=2, figsize=fig_size, sharey=True, gridspec_kw=gridkw)
    #print('Figure size: ', fig_size)
    for n, ax, col_val in zip(count(), axes.flatten(), col_values):
        # First row will contain the graphs
        if n < n_col:
            df_col = df4plot.iloc[df4plot.index.get_level_values(col_var) == col_val]
            df_col = df_col.droplevel(col_var, axis="index").sort_index(key=lambda xx: xx.str.lower(), ascending=True)
            if df_col.index.nlevels > 1:
                df2plot_titles = df_col.index.map(('\n'.join))
            else: 
                df_col_copy = df_col.copy()
                dict_legend = def_legends(df_col_copy.reset_index(), df_type = 'changes')
                df2plot_titles = dict_legend[x_var]
                
            print('\n - '+col_val +'\n'); print(df_col)#, df2plot_titles)
            # Initialize the bottom at zero for the first set of bars.
            bottom_positive = np.zeros(len(df_col)) 
            bottom_negative = np.zeros(len(df_col))
            bottom = np.zeros(len(df_col))
            
            # Plot each layer of the bar, adding each bar to the "bottom" so
            # the next bar starts higher.
            for i, colm in enumerate(df_col.columns):
                #print('ACAAA'); print(df_col[colm])
                for v, val in enumerate(df_col[colm]):
                    #print(val)
                    if val > 0: 
                        bottom[v] = bottom_positive[v]
                        bottom_positive[v] += val
                    else: 
                        bottom[v] = bottom_negative[v]
                        bottom_negative[v] += val
                #print('bottom', bottom)
                ax.bar(df2plot_titles, df_col[colm], bottom=bottom, label=colm, color=colours[i], width = 0.7)
                bottom += np.array(df_col[colm])
                ax.set_xticks(df2plot_titles)
                ax.set_xticklabels(df2plot_titles, rotation=30)
            y_vals = ax.get_yticks()
            y_vals_all.append(y_vals)
            max_y_vals.append(y_vals[-1])
            
            if n == 0:
                handles, labels = ax.get_legend_handles_labels()
                ax.set(ylabel=ylabel)

            # Let's put the annotations inside the bars themselves by using a
            # negative offset.
            if sub_bar_lab:
                y_offset = -0.5
                # For each patch (basically each rectangle within the bar), add a label.
                for bar in ax.patches:
                    if bar.get_height() < 100:
                        bar_height = "{0:.1f}".format(bar.get_height())
                    else: 
                        bar_height = locale.format("%d",  bar.get_height(), grouping=True)
                    ax.text(
                          # Put the text in the middle of each bar. get_x returns the start
                          # so we add half the width to get to the middle.
                          bar.get_x() + bar.get_width() / 2,
                          # Vertically, add the height of the bar to the start of the bar,
                          # along with the offset.
                          bar.get_height()//2 + bar.get_y()+ y_offset,
                          # This is actual value we'll show. original: round(bar.get_height())
                          bar_height,
                          # Center the labels and style them a bit.
                          ha='center', color='w', weight='bold', size=10)
            
            ax.set_title(col_legend[n])
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.axhline(0, color='dimgrey', linewidth=0.8)
        
        # Second row will contain the legend
        if n >= n_col:
            ax.set_axis_off()
            ax.set_title('')
            if n== leg_pos:
                legends = def_var_names('bar_plots', labels)
                ax.legend(handles, legends, loc=loc, bbox_to_anchor = bbox_to_anchor, frameon = True, fontsize=10, ncol = ncols_leg)
    
    # Define axes based on higherst bar
    max_y_index = max_y_vals.index(max(max_y_vals))
    ax.set_yticks(y_vals_all[max_y_index])
    ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
    # ax.ticklabel_format(axis='y', style='plain')

    fig.suptitle(title+txt_title+'\n', fontsize = 11, y=1.01)
    
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'pl_chgbars', 'R_')
            title2save = title.replace(' per Stage and Genotype','')
            title2save = title2save.replace(' per Genotype and Stage','')
            title2save = title2save.replace('\n','_')
            title2save = title2save.replace(' ','_')
            if info != '':
                fig_title = dir2savef+info+"_"+title2save+"."+extf
            else: 
                fig_title = dir2savef+title2save+"."+extf
            # print(fig_title)
            plt.savefig(fig_title, dpi=300, bbox_inches='tight', transparent=True)

#%% func - plotKDEs
def plotKDEs(classif, classif_lab, df_PDF, save, dir2save, info, ext, dpi = 300):
    
    stages = sorted(df_PDF.Stage.unique())
    n_stages = len(stages)
    genots = sorted(df_PDF.Genotype.unique(), reverse = True)
    dict_legends = def_legends(df_PDF)
    stages.append('')

    pal = [['tomato','gold'],['deepskyblue', 'darkblue'],['greenyellow','darkgreen']]
    err_pal = [['lightcoral','goldenrod','lightcoral','goldenrod'],['cyan', 'steelblue','cyan', 'steelblue'],['yellowgreen','limegreen','yellowgreen','limegreen']]
    
    width_ratios=[1 for n in range(n_stages)]+[0.1]
    gridkw = dict(width_ratios=width_ratios)
    for j, cl, cl_lab in zip(count(), classif, classif_lab): 
        xgrid_max = df_PDF['x_grid'].max()
        # print('xgrid_max:', xgrid_max)
        fig, axes = plt.subplots(nrows=1, ncols=n_stages+1, figsize=(13*n_stages, 10), gridspec_kw=gridkw, sharex = False, sharey = True)
        fig.subplots_adjust(wspace=0.1)
        for n, ax, stg in zip(count(), axes.flatten(), stages):
            if n < n_stages:
                df_plot = df_PDF[df_PDF['Stage'] == stg]
                m = sns.lineplot(data=df_plot, x="x_grid", y=cl, hue="Genotype", hue_order = genots, 
                                 estimator='mean', ci=95, ax = ax, palette = pal[j], err_kws = dict(alpha = 0.3))
                fill_under_lines(color = err_pal[j], ax = ax, alpha = 0.2)
                ax.set_xlim(0, xgrid_max)
                ax.set_ylim(0, 1)        
                ax.get_legend().remove()
                ax.set(xlabel='Cardiac jelly thickness [um]', ylabel='\nThickness distribution\n'+ cl_lab + '\n')
                ax.title.set_text("Stage: "+stg+'hpf\n')
                sns.despine()
                if n == 0:
                    handles, labels = m.get_legend_handles_labels()
                    print(cl_lab,'- Genotypes: ',labels)
                # if n == 2: 
                #     ax.legend(handles, gen_legend, loc='lower right', frameon = False)
            else: 
                ax.set_axis_off()
                ax.legend(handles, dict_legends['GenotypeAll'], loc='upper left', frameon = False)
        
        dir2savef = os.path.join(dir2save, 'pl_kde', 'R_')
        if save: 
            for extf in ext:  
                plt.savefig(dir2savef+info+"kdeAll_"+cl+"."+extf, dpi=300, bbox_inches='tight', transparent=True)

#%% func - kde_sklearn
def kde_sklearn(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scikit-learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)

#%% func - kdeThPlots
def kdeThPlots(filename, df_file, file_num, variable, thData, dir2save, save = True):
    """
    Function that creates Kernel Density Estimation Plots for each region of the heart (e.g, atrium, ventricle, 
    left, right, dorsal and ventral and then

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    df_file : TYPE
        DESCRIPTION.
    file_num : TYPE
        DESCRIPTION.
    variable : TYPE
        DESCRIPTION.
    thData : TYPE
        DESCRIPTION.
    dir2save : TYPE
        DESCRIPTION.
    save : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    df_pdfs : TYPE
        DESCRIPTION.
        
    Some links for reference:
        #https://stackabuse.com/kernel-density-estimation-in-python-using-scikit-learn/
        #https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
        #https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html

    """
    if variable == 'cj_thickness':
        title = filename + ' - Cardiac Jelly Thickness [um]'
        xlabel = 'Cardiac Jelly Thickness [um]'
        step = 0.05
    elif variable == 'myoc_intBall' :
        title = filename + ' - Myoc.Int Ballooning [um]'
        xlabel = 'Myoc.Int Ballooning [um]'
        step = 0.25
    
    regions_div = [['AtrVent','atrium', 'ventricle'], ['LeftRight','left','right'],['DorsVent','dorsal','ventral']]
    color = [['','tomato','gold'],['','deepskyblue', 'darkblue'],['','greenyellow', 'darkgreen']]
    num_linspace = int((max(round(thData[variable]))/step)+2)
    x_grid = np.linspace(0, max(round(thData[variable]))+step, num_linspace)

    genotype = df_file.loc[file_num,'Gene_A']+':'+df_file.loc[file_num,'Genotype_A']
    if df_file.loc[file_num,'Gene_B'] != '-':
        genotype = str(genotype+'/'+df_file.loc[file_num,'Gene_B']+':'+df_file.loc[file_num,'Genotype_B'])
        
    df_pdfs = pd.DataFrame(columns = ['Filename','Strain','Stage','Genotype', 'x_grid','AtrVent','atrium', 'ventricle', 'LeftRight','left','right','DorsVent','dorsal','ventral'])
    df_pdfs['Filename'] = [filename for j in range(len(x_grid))]
    df_pdfs['Strain'] = [df_file.loc[file_num,'Strain'] for j in range(len(x_grid))]
    df_pdfs['Stage'] = [df_file.loc[file_num,'Stage'] for j in range(len(x_grid))]
    df_pdfs['Genotype'] = [genotype for j in range(len(x_grid))]
    df_pdfs['x_grid'] = x_grid
    plot_dir = os.path.join(dir2save, filename+"_")
    
    print('\n- Creating density plots... this process takes a while, about 8-12 min/plot out of 3 plots, be patient :)' )
    bar = Bar('- Creating density plots', max=3, suffix = suffix, check_tty=False, hide_cursor=False)
    for i, reg, col in zip(count(), regions_div, color):
        df_one = thData[thData[reg[0]] == reg[1]]
        th_one = df_one[variable]
        bw_one = len(th_one)**(-1./(1+4))*np.std(th_one)
        pdf_one = kde_sklearn(th_one, x_grid, bandwidth=bw_one)
        
        df_two = thData[thData[reg[0]] == reg[2]]
        th_two = df_two[variable]
        bw_two = len(th_two)**(-1./(1+4))*np.std(th_two)
        pdf_two = kde_sklearn(th_two, x_grid, bandwidth=bw_two)
        
        print('\n - bandwidths:', format(bw_one,'.2f'), '-',format(bw_two, '.2f'))
        pdf_comb = np.add(pdf_one, pdf_two)
        
        pct_one = pdf_one/pdf_comb
        # pct_two = pdf_two/pdf_comb
    
        ones = np.ones((len(x_grid),))
        
        fig, ax = plt.subplots(figsize=(8,5))
        ax.fill_between(x_grid, pct_one, alpha=0.5, label= reg[1], color = col[1])
        ax.fill_between(x_grid, pct_one, ones, alpha=0.5, label= reg[2], color= col[2])
        # ax.plot(x_grid, pct_one, linewidth=1, alpha=0.5, label= reg[1])
        # ax.plot(x_grid, pct_two, linewidth=1, alpha=0.5, label= reg[2])
        ax.set_xlim(0, x_grid[-1])
        ax.set_ylim(0, 1)        
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.set_title(title, fontsize = 10)
        ax.set_xlabel(xlabel, fontsize=10)
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        if save: 
            plt.savefig(plot_dir+"kde"+reg[0]+".png", dpi=300, bbox_inches='tight', transparent=True)
        plt.show()

        # fig, ax = plt.subplots()
        # ax.plot(x_grid, pdf_one, linewidth=1, alpha=0.5, label= 'pdf_'+reg[1])
        # ax.plot(x_grid, pdf_two, linewidth=1, alpha=0.5, label= 'pdf_'+reg[2])
        # ax.legend(loc='upper right')
        
        df_pdfs[reg[0]] = pct_one
        df_pdfs[reg[1]] = pdf_one
        df_pdfs[reg[2]] = pdf_two
        
        bar.next()
        
    bar.finish()
    alert('wohoo',1)
    
    return df_pdfs

#%% func - fill_under_lines
def fill_under_lines(color, ax=None, alpha=.2, **kwargs):
        if ax is None:
            ax = plt.gca()
        # print('ja:',len(ax.lines))
        for i, line in enumerate(ax.lines):
            x, y = line.get_xydata().T
            # print(color[i])
            ax.fill_between(x, 0, y, color=color[i], alpha=alpha, **kwargs)
            
#%% func - plotKDEIndiv
def plotKDEIndiv(classif, classif_lab, df_PDF, save, dir2save, info, ext, dpi = 300):
    
    stages_o = sorted(df_PDF.Stage.unique())
    n_stages = len(stages_o)
    stages = []
    for n in range(n_stages):
        stages.append(stages_o[n])
        stages.append('')
        
    # strains = sorted(df_PDF.Strain.unique())
    # genots = sorted(df_PDF.Genotype.unique(), reverse = True)
    # # gen_legend, strain_legend, stage_legend = 
    # dict_legends = def_legends(df_PDF)
    gridkw = dict(height_ratios=[1 for n in range(n_stages)], width_ratios=[1,0.3])
    
    for j, cl, cl_lab in zip(count(), classif, classif_lab): 
        xgrid_max = df_PDF['x_grid'].max()
        # print('xgrid_max:', xgrid_max)
        fig, axes = plt.subplots(nrows=n_stages, ncols=2, figsize=(12, 10*n_stages), gridspec_kw=gridkw, sharex = False, sharey = True)
        fig.subplots_adjust(hspace=0.5)
        for n, ax, stg in zip(count(), axes.flatten(), stages):
            # print(n)
            if n in [0,2,4,6]:
                df_plot = df_PDF[df_PDF['Stage'] == stg]
                m = sns.lineplot(data=df_plot, x="x_grid", y=cl, hue='Filename', style = "Genotype",# hue_order = genots, 
                                 estimator='mean', ci=95, ax = ax, #palette = pal[j], 
                                 err_kws = dict(alpha = 0.3))
                # fill_under_lines(color = err_pal[j], ax = ax, alpha = 0.2)
                ax.set_xlim(0, xgrid_max)
                ax.set_ylim(0, 1)        
                ax.get_legend().remove()
                ax.set(xlabel='Cardiac jelly thickness [um]', ylabel='\nThickness distribution\n'+ cl_lab + '\n')
                ax.title.set_text("Stage: "+stg+'hpf\n')
                sns.despine()
                handles, labels = m.get_legend_handles_labels()
            else: 
                # print('aja')
                ax.set_axis_off()
                ax.legend(handles, labels, loc='upper left', fontsize = 'xx-small', frameon = False)
                
        dir2savef = os.path.join(dir2save, 'pl_kde', 'R_')
        if save: 
            for extf in ext:  
                plt.savefig(dir2savef+info+"kdeIndivAll_"+cl+"."+extf, dpi=300, bbox_inches='tight', transparent=True)

#%% func - getHeatmaps2Unify
def getHeatmaps2Unify(folders, chamber, thickness, dir2load_df, dir_data2Analyse, normalise = False, 
                      opt_norm = 'opt_div', perChamber = False):
    
    # List for original dataframes
    dfs_o = []
    # List for normalised dataframes
    if normalise:
        dfs_hmf = []
        norm_vals = []
    
    num = 0
    # print('thickness:', thickness)
    for n, file in enumerate(folders):
        # print(file)
        name = 'unloop'+chamber+'_'+thickness
        hmf_file = 'hmf_'+name
        # print('hmf_file:',hmf_file)
        try: 
            df_load = loadDF(file[2:], hmf_file, dir2load_df)
            dfs_o.append(df_load)
            if normalise: 
                dir_results = os.path.join(dir_data2Analyse, 'R_all', 'df_all', 'df_meas')
                df_res = loadDF(file[2:], 'ResultsDF', dir_results)
                file_num = df_res[df_res['Folder']==file[2:]+'_2A'].index.values[0]
                # print(file[2:], file_num)
                # print(df_res.head(5))
                if perChamber:
                    dict_norm = {'CjTh': {'unloopAtr_CjTh': df_res.loc[file_num,'Vol_Atr.CJ'], 'unloopVent_CjTh': df_res.loc[file_num,'Vol_Vent.CJ']},
                                  'MyocTh': {'unloopAtr_MyocTh': df_res.loc[file_num,'Vol_Atr.Myoc'], 'unloopVent_MyocTh': df_res.loc[file_num,'Vol_Vent.Myoc']},
                                  'EndoTh': {'unloopAtr_EndoTh': df_res.loc[file_num,'Vol_Atr.Endo'], 'unloopVent_EndoTh': df_res.loc[file_num,'Vol_Vent.Endo']},
                                  'myocIntBall': {'unloopAtr_myocIntBall': df_res.loc[file_num,'Vol_Atr.ExtMyoc'], 'unloopVent_myocIntBall': df_res.loc[file_num,'Vol_Vent.ExtMyoc']}}
                else: 
                    dict_norm = {'CjTh': {'unloopAtr_CjTh': df_res.loc[file_num,'Vol_CJ'], 'unloopVent_CjTh': df_res.loc[file_num,'Vol_CJ']},
                                  'MyocTh': {'unloopAtr_MyocTh': df_res.loc[file_num,'Vol_Myoc'], 'unloopVent_MyocTh': df_res.loc[file_num,'Vol_Myoc']},
                                  'EndoTh': {'unloopAtr_EndoTh': df_res.loc[file_num,'Vol_Endo'], 'unloopVent_EndoTh': df_res.loc[file_num,'Vol_Endo']},
                                  'myocIntBall': {'unloopAtr_myocIntBall': df_res.loc[file_num,'Vol_Ext.Myoc'], 'unloopVent_myocIntBall': df_res.loc[file_num,'Vol_Ext.Myoc']}}
                # print(dict_norm)
                if opt_norm == 'opt_div':
                    norm_val = 1000000/dict_norm[thickness][name]
                elif opt_norm == 'opt_mult':
                    norm_val = dict_norm[thickness][name]/1000000
                norm_vals.append(norm_val)
                df_load = df_load*norm_val
                dfs_hmf.append(df_load)
            num += 1
        
        except: 
            print('-No hmf heatmap dataframe found for '+thickness+' within '+file+' results folder!')
            continue
        
    if normalise: 
        return [dfs_o, dfs_hmf], num, norm_vals
    else: 
        return [dfs_o], num
    
#%% func - unifyHeatmap
def unifyHeatmap(df, chamber, stage, genotype, gen_info, thickness, vmin, vmax, n_val, normalise, dir2save, save, cmap = 'turbo'):
    
    sns.set_context('notebook', font_scale=1.25)
    stage = stage+'hpf'
    if thickness == 'CjTh':
        title = 'Cardiac jelly thickness [um] - ('+chamber+', '+stage+', '+genotype+' - n='+str(n_val)+')_'+normalise+'\n'
    elif thickness == 'myocIntBall': 
        title = 'Myocardium ballooning [um] - ('+chamber+', '+stage+', '+genotype+' - n='+str(n_val)+')_'+normalise+'\n'
    elif thickness == 'MyocTh':
        title = 'Myocardial thickness [um] - ('+chamber+', '+stage+', '+genotype+' - n='+str(n_val)+')_'+normalise+'\n'
    elif thickness == 'EndoTh':
        title = 'Endocardial thickness [um] - ('+chamber+', '+stage+', '+genotype+' - n='+str(n_val)+')_'+normalise+'\n'
    
    # Make figure 
    df_max = format(df.max().max(), '.2f')
    fig, ax = plt.subplots(figsize=(16, 10))
    if 'std' in gen_info:
        title = '(STD - max:'+df_max+'um) '+ title 
    elif 'sem' in gen_info:
        title = '(SEM - max:'+df_max+'um) '+ title 
    ax = sns.heatmap(df, cmap=cmap, vmin = vmin, vmax = vmax)
    
    max_val = df.index.max()
    if max_val > 1.1: 
        y_labels = [1,2]
        y_text = '[Inflow tract >> Valve]'
    else: # Ventricle
        y_labels = [0,1]
        y_text = '[Valve >> Outflow tract]'
        
    x_pos = ax.get_xticks()
    # x_lab = ax.get_xticklabels()
    x_pos_new = np.linspace(x_pos[0], x_pos[-1], 19)
    x_lab_new = np.arange(-180,200,20)
    ax.set_xticks(x_pos_new) 
    ax.set_xticklabels(x_lab_new, rotation=30)
    
    y_pos = ax.get_yticks()
    # y_lab = ax.get_yticklabels()
    y_pos_new = np.linspace(y_pos[0], y_pos[-1], 11)
    y_lab_new = np.linspace(y_labels[0],y_labels[1],11)
    y_lab_new = [format(y,'.2f') for y in y_lab_new]
    ax.set_yticks(y_pos_new) 
    ax.set_yticklabels(y_lab_new, rotation=0)
    
    plt.ylabel('Centreline position '+y_text+'\n', fontsize=15)
    plt.xlabel('Angle (\N{DEGREE SIGN}) [Dorsal >> Right >> Ventral >> Left >> Dorsal]', fontsize=15)
    plt.title(title, fontsize = 15)
    
    dir4heatmap = os.path.join(dir2save,'pl_hmnorm', 'hmfAll_'+gen_info+'_'+thickness+'_'+chamber+'_'+stage+'_'+normalise+'.png')
    # print(dir4heatmap)
    if save: 
        plt.savefig(dir4heatmap, dpi=300, bbox_inches='tight', transparent=True)

#%% Palette stuff
    # Save a palette to a variable:
    # palette = sns.color_palette("husl", 8*n_gen*n_strain)
    # sns.palplot(palette)
    # for gen in range(n_gen):
    #     palettes.append(palette[gen*len(palette)//n_gen+2])
    
#%% - ALERT WHEN IMPORTED
print ("IMPORTED: morphoHeart_funcAnalysis")
alert('jump',1)

#%% - OTHERS TO ORGANIZE
#%% Others -  Plots of all groups of variables x3 (stages)

    # for i, input_var, title in zip(count(), input_vars, titles):
    #     vars2plot, labels2plot = fcAn.selectVariables_auto(variables, ylabels, input_var)
        
    #     # Set up the matplotlib figure
    #     num_vars = len(vars2plot)
    #     plots_per_col = len(stages)
    #     plots_per_row = num_vars
    #     stages = sorted(df2plot.Stage.unique())
    #     n_stages = len(stages)
        
    #     vars2plotX3 = []
    #     labels2plotX3 = []
    #     for ii, var, lab in zip(count(), vars2plot, labels2plot): 
    #         for n in range(n_stages):
    #             vars2plotX3.append(var)
    #             labels2plotX3.append(lab)
        
    #     index_first = list(range(0,(plots_per_col+1)*plots_per_row,4))
    #     stagesX3 = stages*num_vars
    #     index_no_plot = list(range(3,(plots_per_col+1)*plots_per_row,4))
    #     for index in index_no_plot:
    #         vars2plotX3.insert(index, '')
    #         labels2plotX3.insert(index, '')
    #         stagesX3.insert(index,'')
            
    #     # Set up the matplotlib figure
    #     size_col = (plots_per_col+1)*h_plot+h_add-5
    #     size_row = plots_per_row*w_plot+w_add+20
        
    #     # Number of genotypes: 
    #     genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
    #     n_gen = len(genots)
    #     strains = sorted(df2plot.Strain.unique())
    #     n_strain = len(strains)

    #     if i == 0: 
    #         print('Genotypes: ', genots)
    #         print('Strains: ', strains)
    #         print('Stages: ', stages)
        
    #     palettes = ['mediumturquoise', 'darkmagenta']
        
    #     #  Create figure  - plt.clf()
    #     gridkw = dict(width_ratios=[1,1,1,0.2])
    #     fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
    #     fig.subplots_adjust(hspace=1.5, wspace=0)
        
    #     sns.set_style("ticks")
    #     sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
    #                                                     "ytick.labelsize":'small', "lines.linewidth": 2.5,
    #                                                     "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
    #     # Define legends for strain and genotype
    #     legend_elem_gen = []
    #     for gen, gen_lab in enumerate(gen_legend):
    #         legend_elem_gen.append(Line2D([0], [0], marker='h', color='w', label=gen_lab,
    #                                 markerfacecolor=palettes[gen], markersize=20))
    #     handle_new = legend_elem_gen
    #     legend_new = gen_legend
        
    #     marker_size = 12
    #     dodge = False
    #     jitter = 0.2
    #     for n, ax, stg, var, ylabel in zip(count(), axes.flatten(), stagesX3, vars2plotX3, labels2plotX3):
    #         if n in index_no_plot:
    #             if n == 3:
    #                 ax.set_axis_off()
    #                 ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(0, 1), frameon = False)
    #             else: 
    #                 ax.remove()
    #         else: 
    #             df_plot = df2plot[df2plot['Stage'] == stg]
    #             m = sns.stripplot(x="Strain", y=var, hue="GenotypeAll", hue_order = genots, data=df_plot, ax = ax, order=strains,
    #                           marker = style, palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
    #             box = ax.get_position()
                
    #             m.set_xticklabels(strain_legend, rotation=45)
    #             # m.set_axis_labels(strains)
    #             if n in index_first:
    #                 ax.set(ylabel=ylabel)
    #             else: 
    #                 ax.set(ylabel = '')
    #                 ax.tick_params(axis = 'y', labelcolor='w', width=0.1)
    #                 ax.spines['left'].set_color('gray')
    #                 ax.spines['left'].set_linestyle("dashed")
    #             ax.set(xlabel="\nStage: "+stg+'hpf')
    #             ax.set_position([box.x0, box.y0, box.width*1, box.height])
    #             ax.get_legend().remove()
                
    #             # ax.legend(handles0[:],labels0[:],loc='center left', bbox_to_anchor=(1, 0.5))
    #             sns.despine()
                
    #             if n == 0:
    #                 handles, labels = m.get_legend_handles_labels()
    #                 print(handles)
    #                 print(labels)

    #     fig.suptitle(title, fontsize = 30, y=0.9)
        
        # if save: 
        #     plt.savefig(dir_pl_meas+title+"StgStrainSplit.png", dpi=300, bbox_inches='tight', transparent=True)
        
# func - compareKDE
# import seaborn as sns
# sns.set_theme(style="darkgrid")

# # Load an example dataset with long-form data
# fmri = sns.load_dataset("fmri")

# # Plot the responses for different events and regions
# sns.lineplot(x="timepoint", y="signal",
#              hue="region", style="event",
#              data=fmri)

# flights = sns.load_dataset("flights")
# flights.head()

# flights_wide = flights.pivot("year", "month", "passengers")
# flights_wide.head()
# sns.lineplot(data=flights_wide)

# sns.lineplot(data=flights, x="year", y="passengers")
# sns.lineplot(data=flights, x="year", y="passengers", hue="month")

#%% TO REMOVE
#%%func - plotStatisticsOLD
# def plotStatisticsOLD(data, x, y, x_label, y_label, box_pairs, figsize, genot2filter,
#                     colors_bgd, colors_pts, order2plot,
#                     alleleName, txt_normtest, multcomp_txt, txt_testSelected, txt_note,
#                     test2use, pval_test, pval_multComp,
#                     plots_dir, saveFig, ext='png'):
    
#     if test2use == 'One-way-ANOVA' or test2use == 'Kruskal':
#         stat_test = False
#         p_val = pval_multComp
#         test = None
#         txt_testSelected = txt_testSelected+' (p-value: '+str(format(pval_test,'.3f'))+'), '+multcomp_txt
#     else: 
#         stat_test = True
#         p_val = None
#         test = test2use
    
#     data_filtered = data[data[y].notnull()]
#     exp_ref_order = sorted(data_filtered.Exp_Ref.unique().tolist())
    
#     n_wt = str(len(data_filtered[data_filtered['Genotype']=='wt']))
#     n_ht = str(len(data_filtered[data_filtered['Genotype']=='ht']))
#     n_mt = str(len(data_filtered[data_filtered['Genotype']=='mt']))
    
#     #Define xtickslabel and figure sizes
#     if genot2filter == 'A':
#         xticklabel = ["$hapln1a^{+/+}$\nn ="+n_wt, "$hapln1a^{+/-}$\nn ="+n_ht, "$hapln1a^{-/-}$\nn ="+n_mt]
#         figsize = (6, 8.5)
#     else:
#         xticklabel = ["$hapln1a^{+/+}$\nn ="+n_wt, "$hapln1a^{-/-}$\nn ="+n_mt]
#         figsize = (5, 8.5)
    
#     gridkw = dict(height_ratios=[4, 1])
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, gridspec_kw=gridkw)
#     # Axis 1
#     g = sns.violinplot(x=x, y=y, palette = colors_bgd, inner='quartile', scale='count',
#                         width=0.5, linewidth=0.5, legend_out=True, data=data, 
#                         order=order2plot, ax=ax1)
#     sns.swarmplot(x=x, y=y, hue="Exp_Ref", hue_order= exp_ref_order, palette = colors_pts, size=5, 
#                   order=order2plot, data=data, ax=ax1).set(xlabel=x_label, ylabel=y_label)
    
#     ax1, test_results = add_stat_annotation(g, data=data, x=x, y=y,
#                                             perform_stat_test=stat_test, 
#                                             test=test,
#                                             pvalues=p_val,
#                                             comparisons_correction=None,#'bonferroni'
#                                             order=order2plot, box_pairs=box_pairs,
#                                             line_offset_to_box=0.1, line_offset=0.01,
#                                             line_height=0.015, 
#                                             text_format='star', verbose=0)
#     ax1.legend(bbox_to_anchor=(1, 1), loc=0)
    
#     txt_multcomp = txtMultComp(box_pairs, test_results)
        
#     g.set_xticklabels(xticklabel)
    
#     # Axis 2
#     ax2.set_axis_off()        
#     t = (txt_normtest+'\n'+txt_testSelected+'\n'+txt_multcomp+'\n'+txt_note)
#     ax2.text(0.5, 0.25, t, fontdict=font, wrap=True)
#     plt.suptitle(alleleName+" - "+y_label, y=0.92)
    
#     if saveFig:
#         plt.savefig(plots_dir+y+"_st."+ext, dpi=300, bbox_inches='tight')
    
#     return test_results




#%% func - plotInGroups2
# def plotInGroups2 (plot_type, input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
#                      h_plot, w_plot, save, dir2save, info, dpi = 300, sharey = False, h_add = 5, w_add = 1, ext = 'png'):
    
#     vars_dict = def_variables(plot_type)
    
#     for i, input_var, title in zip(count(), input_vars, titles):
#         vars2plot, labels2plot = selectVariables_auto(vars_dict, input_var)
#         sns.set_context('poster') # notebook, talk, poster, paper
#         # Set up the matplotlib figure
#         num_vars = len(vars2plot)
#         plots_per_col = 3
#         plots_per_row = math.ceil(num_vars/plots_per_col)
        
#         index_no_plot = list(range(3,(plots_per_col+1)*plots_per_row,4))
#         for index in index_no_plot:
#             vars2plot.insert(index, '')
#             labels2plot.insert(index, '')
            
#         # Set up the matplotlib figure
#         size_col = (plots_per_col+1)*h_plot+h_add
#         size_row = plots_per_row*w_plot+w_add
        
#         # Genotypes and Strains being plotted 
#         genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
#         strains = sorted(df2plot.Strain.unique())
#         stages = sorted(df2plot.Stage.unique())
        
#         if i == 0: 
#             print('- Genotypes: ', genots)
#             print('- Strains: ', strains)
#             print('- Stages: ', stages)
        
#         palettes = ['deeppink','royalblue','mediumturquoise', 'darkmagenta','darkorange','limegreen', 'gold']
       
#         #  Create figure  - plt.clf()
#         gridkw = dict(width_ratios=[1,1,1,0.2])
#         fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
#         fig.subplots_adjust(hspace=0.5, wspace=0.5)
        
#         sns.set_style("ticks")
#         sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
#                                                         "ytick.labelsize":8, "lines.linewidth": 2.5,
#                                                         "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
#         # Define legends for strain and genotype
#         legend_elem_gen = []
#         for gen, gen_lab in enumerate(gen_legend):
#             legend_elem_gen.append(Line2D([0], [0], marker='h', color='w', label=gen_lab,
#                                     markerfacecolor=palettes[gen], markersize=20))
#         space = [Line2D([0], [0], marker='o', color='w', label='',
#                                     markerfacecolor='w', markersize=20)]
#         legend_elem_strain = []
#         for n_str, str_lab, mark in zip(count(), strain_legend, styles):
#             legend_elem_strain.append(Line2D([0], [0], marker=mark, color='w', label=str_lab,
#                                     markerfacecolor='k', markersize=15))
            
#         handle_new = legend_elem_gen+space+legend_elem_strain
#         legend_new = gen_legend+['']+strain_legend
        
#         marker_size = 12
#         dodge = True
#         jitter = 0.3
#         for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
#             if n in index_no_plot:
#                 if n == 3:
#                     ax.set_axis_off()
#                     ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
#                 else: 
#                     ax.remove()
#             else: 
#                 for j, strain, style in zip(count(), strains, styles):
#                     df_plot = df2plot[df2plot['Strain'] == strain]
#                     m = sns.stripplot(x="Stage", y=var, hue="GenotypeAll", hue_order = genots, data=df_plot, ax = ax, order=['32-34','48-50','72-74'],
#                                   marker = style, palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
#                     box = ax.get_position()
#                     ax.set(xlabel="Stage [hpf]", ylabel=ylabel);
#                     ax.set_position([box.x0, box.y0, box.width*1, box.height])
#                     ax.get_legend().remove()
                    
#                     # ax.legend(handles0[:],labels0[:],loc='center left', bbox_to_anchor=(1, 0.5))
#                     sns.despine()
                    
#                     if n == 0:
#                         handles, labels = m.get_legend_handles_labels()
#                         #print(handles)
#                         print(labels)

#         fig.suptitle(title, fontsize = 30, y=1)
#         if save: 
#             for extf in ext: 
#                 dir2savef = os.path.join(dir2save, 'meas_all', 'R_')
#                 if info != '':
#                     fig_title = dir2savef+info+"_"+title+"."+extf
#                 else: 
#                     fig_title = dir2savef+title+"."+extf

#                 plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
            
#%% func - plotRelPlotInGroups2
# def plotRelPlotInGroups2(plot_type, vars2plot, labels2plot, titles, df2plot, gen_legend, strain_legend , stage_legend,
#                      h_plot, w_plot, save, dir2save, info, dpi = 300, sharey = False, 
#                      h_add = 5, w_add = 1, ext = 'png'):
    
#     # vars_dict = def_variables(plot_type)
    
#     for i, input_var, title in zip(count(), vars2plot, titles):
        
#         # vars2plot, labels2plot = selectVariables_auto(vars_dict, input_var)
#         sns.set_context('poster') # notebook, talk, poster, paper
#         # Set up the matplotlib figure
#         num_vars = len(vars2plot)
#         plots_per_col = 3
#         plots_per_row = math.ceil(num_vars/plots_per_col)
        
#         index_no_plot = list(range(3,(plots_per_col+1)*plots_per_row,4))
#         for index in index_no_plot:
#             vars2plot.insert(index, '')
#             labels2plot.insert(index, '')
            
#         # Set up the matplotlib figure
#         size_col = (plots_per_col+1)*h_plot+h_add
#         size_row = plots_per_row*w_plot+w_add
        
#         # Genotypes and Strains being plotted 
#         genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
#         strains = sorted(df2plot.Strain.unique())
#         stages = sorted(df2plot.Stage.unique())
        
#         if i == 0: 
#             print('- Genotypes: ', genots)
#             print('- Strains: ', strains)
#             print('- Stages: ', stages)
        
#         palettes = ['deeppink','royalblue','mediumturquoise', 'darkmagenta','darkorange','limegreen', 'gold']
       
#         #  Create figure  - plt.clf()
#         gridkw = dict(width_ratios=[1,1,1,0.2])
#         fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
#         fig.subplots_adjust(hspace=0.5, wspace=0.5)
        
#         sns.set_style("ticks")
#         sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
#                                                         "ytick.labelsize":8, "lines.linewidth": 2.5,
#                                                         "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
#         # Define legends for strain and genotype
#         legend_elem_gen = []
#         for gen, gen_lab in enumerate(gen_legend):
#             legend_elem_gen.append(Line2D([0], [0], marker='h', color='w', label=gen_lab,
#                                     markerfacecolor=palettes[gen], markersize=20))
#         space = [Line2D([0], [0], marker='o', color='w', label='',
#                                     markerfacecolor='w', markersize=20)]
#         legend_elem_strain = []
#         for n_str, str_lab, mark in zip(count(), strain_legend, styles):
#             legend_elem_strain.append(Line2D([0], [0], marker=mark, color='w', label=str_lab,
#                                     markerfacecolor='k', markersize=15))
            
#         handle_new = legend_elem_gen+space+legend_elem_strain
#         legend_new = gen_legend+['']+strain_legend
        
#         for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
#             # if n == 0:
#             #     min_val = 0.8*df2plot[var].min()
#             #     max_val = 1.2*df2plot[var].max()
                
#             if n in index_no_plot:
#                 if n == 3:
#                     ax.set_axis_off()
#                     ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
#                 else: 
#                     ax.remove()
#             else: 
#                 m = sns.lineplot(x="time_point", y=var,
#                     estimator=None, data= df2plot, palette = 'Set2', 
#                     ax = ax)
#                 box = ax.get_position()
#                 # ax.set_ylim([min_val, max_val])
#                 ax.set(xlabel="Heart phase contraction", ylabel=ylabel);
#                 ax.set_position([box.x0, box.y0, box.width*1, box.height])
#                 try:
#                     ax.label_outer()
#                 except:
#                     pass
#                 sns.despine()
        
#                 # if n == 0:
#                 #     handles, labels = m.get_legend_handles_labels()
#                 #     #print(handles)
#                 #     print(labels)

#         fig.suptitle(title, fontsize = 30, y=1)
#         dir2savef = os.path.join(dir2save, 'meas_all', 'R_')
#         if info != '':
#             fig_title = dir2savef+info+"_"+title+"."+ext
#         else: 
#             fig_title = dir2savef+title+"."+ext
        
#         if save: 
#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
            
#%% func - plotPerVariable
# def plotPerVariable(script, input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
#                      h_plot, w_plot, save, dir2save, info, dpi = 300, h_add = 5, w_add = 1, ext = 'png'):
    
#     styles = ['o', '^', 's', 'v', 'D', '<', 'p', '>'] # https://matplotlib.org/stable/api/markers_api.html
#     variables, ylabels = def_variables(script)
    
#     for i, input_var, title in zip(count(), input_vars, titles):
#         vars2plot, labels2plot = selectVariables_auto(variables, ylabels, input_var)
        
#         # Set up the matplotlib figure
#         stages = sorted(df2plot.Stage.unique())
#         plots_per_col = len(stages)
#         plots_per_row = 1
#         stages.append('')
        
#         # Set up the matplotlib figure
#         size_col = (plots_per_col+1)*h_plot-12
#         size_row = plots_per_row*w_plot+w_add
        
#         # Genotypes and Strains being plotted 
#         genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
#         strains = sorted(df2plot.Strain.unique())

#         index_no_plot = list(range(3,(plots_per_col+1)*plots_per_row,4))

#         if i == 0: 
#             print('- Genotypes: ', genots)
#             print('- Strains: ', strains)
#             print('- Stages: ', stages)
        
#         palettes = ['deeppink','royalblue','mediumturquoise', 'darkmagenta','darkorange','limegreen', 'gold']
        
#         for nn, var, ylabel in zip(count(), vars2plot, labels2plot):
#             #  Create figure  - plt.clf()
#             gridkw = dict(width_ratios=[1,1,1,0.2])
#             fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
#             fig.subplots_adjust(hspace=1.5, wspace=0.2)
            
#             sns.set_style("ticks")
#             sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
#                                                             "ytick.labelsize":'small', "lines.linewidth": 2.5,
#                                                             "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
#             # Define legends for strain and genotype
#             legend_elem_gen = []
#             for gen, gen_lab in enumerate(gen_legend):
#                 legend_elem_gen.append(Line2D([0], [0], marker='h', color='w', label=gen_lab,
#                                         markerfacecolor=palettes[gen], markersize=10))
#             handle_new = legend_elem_gen
#             legend_new = gen_legend
            
#             marker_size = 10
#             dodge = False
#             jitter = 0.3
#             for n, ax, stg in zip(count(), axes.flatten(), stages):
#                 # print(n)
#                 if n in index_no_plot:
#                     if n == 3:
#                         ax.set_axis_off()
#                         ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(0, 1), frameon = False)
#                     else: 
#                         ax.remove()
#                 else: 
#                     df_plot = df2plot[df2plot['Stage'] == stg]
#                     m = sns.stripplot(x="Strain", y=var, hue="GenotypeAll", hue_order = genots, data=df_plot, ax = ax, order=strains,
#                                   marker = styles[0], palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
#                     box = ax.get_position()
#                     m.set_xticklabels(strain_legend, rotation=50)
#                     # m.set_axis_labels(strains)
#                     ax.set(xlabel="\nStage: "+stg+'hpf')
#                     if n == 0:
#                         ax.set(ylabel='\n'+ylabel+'\n')
#                     else: 
#                         ax.set(ylabel = '')
#                         ax.tick_params(axis = 'y', labelcolor='w', width=0.1)
#                         ax.spines['left'].set_color('gray')
#                         ax.spines['left'].set_linestyle((0,(4,4)))#"dashed")
#                     ax.set_position([box.x0, box.y0, box.width*1, box.height])
#                     ax.get_legend().remove()
                    
#                     # ax.legend(handles0[:],labels0[:],loc='center left', bbox_to_anchor=(1, 0.5))
#                     sns.despine()
                    
#                     if n == 0:
#                         handles, labels = m.get_legend_handles_labels()
#                         # print(handles)
#                         print(var, '- Genotypes: ',labels)
    
#             fig.suptitle(title, fontsize = 30, y=1)
            
#             if save: 
#                 for extf in ext: 
#                     dir2savef = os.path.join(dir2save, 'meas_Ind', 'R_')
#                     if info != '':
#                         fig_title = dir2savef+info+"_Ind_"+var+"."+extf
#                     else: 
#                         fig_title = dir2savef+"Ind_"+var+"."+ext
    
#                     plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)

#%% func - getHeatmaps2Unify_BU
# def getHeatmaps2Unify_BU(folders, chamber, thickness, dir2save, operation = 'mean'):
    
#     dfs_hmf = []
#     num = 0
#     for n, file in enumerate(folders):
#         print(file)
#         hmf_file = 'hmf_unloop'+chamber+thickness
#         # hmf_file = 'hmN_unloop'+chamber+thickness
#         try: 
#             dfs_hmf.append(loadDF(file[2:], hmf_file, dir2save))
#             num += 1
#         except: 
#             print('-No hmf heatmap dataframe found for '+thickness+' within '+file+' results folder!')
#             continue
        
#     if operation == 'std':
#         df_hmf = pd.concat(dfs_hmf).groupby(level=0).std()
#     elif operation == 'mean':
#         df_hmf = pd.concat(dfs_hmf).groupby(level=0).mean()
#     elif operation == 'sem':
#         df_hmf = pd.concat(dfs_hmf).groupby(level=0).sem()
    
    
#     return df_hmf, num

#%% func - barPlotsBU
# def barPlotsBU(df2plot, vars2plot, colours, title, ylabel, stack100 = True, totals_lab = True, sub_bar_lab = True): 
    
#     y_vals_all = []
#     max_y_vals = []
#     strain_o = True #ask4input('Strain_o? >>:', bool)
    
#     vars_opt = ['Stage','GenotypeAll','Manip']
#     if strain_o: 
#         vars_opt.append('Strain_o')
#     else: 
#         vars_opt.append('Strain')
        
#     group_vars = []
#     txt_title = '\n'
#     for var in vars_opt:
#         if len(df2plot[var].unique()) > 1:
#             group_vars.append(var)
#         else: 
#             txt_title = txt_title +'['+ var + ':'+df2plot[var].unique()[0]+'] '
            
#     df4plot = df2plot.groupby(group_vars)[vars2plot].mean()
#     stages = sorted(df2plot['Stage'].unique())
#     stages.append('')
#     n_genot = len(df2plot['GenotypeAll'].unique())
    
#     if stack100: 
#         df4plot = df4plot.apply(lambda x: x*100/sum(x), axis=1)
#     min_val = df4plot.min().min()
    
#     #Create figure
#     gridkw = dict(width_ratios=[1,1,1,0.2])
#     fig, axes = plt.subplots(ncols=len(stages),nrows=1, figsize=(3*(n_genot/2)*(len(stages)-0.8),6), sharey=True, gridspec_kw=gridkw)
#     # print(3*(n_genot/2)*(len(stages)-0.8))
#     for n, ax, stg in zip(count(), axes.flatten(), stages):
#         if n < len(stages):
#             df_stage = df4plot.iloc[df4plot.index.get_level_values('Stage') == stg]
#             df_stage = df_stage.droplevel('Stage', axis="index")
#             df_stage = df_stage.sort_index(key=lambda x: x.str.lower(), ascending=False)
#             if df_stage.index.nlevels > 1:
#                 df2plot_titles = df_stage.index.map(('\n'.join))
#             else: 
#                 df_stage_copy = df_stage.copy()
#                 gen_legend, _, strains_o_legend, _ = def_legends(df_stage_copy.reset_index())
#                 df2plot_titles = gen_legend
#                 #df2plot_titles = df_stage.index
#             # Initialize the bottom at zero for the first set of bars.
#             bottom = np.zeros(len(df_stage))
            
#             # Plot each layer of the bar, adding each bar to the "bottom" so
#             # the next bar starts higher.
#             for i, col in enumerate(df_stage.columns):
#               ax.bar(df2plot_titles, df_stage[col], bottom=bottom, label=col, color=colours[i])
#               bottom += np.array(df_stage[col])
#               ax.set_xticks(df2plot_titles)
#               ax.set_xticklabels(df2plot_titles, rotation=30)
#             y_vals = ax.get_yticks()
#             # print(y_vals)
#             y_vals_all.append(y_vals)
#             max_y_vals.append(y_vals[-1])
#             # print(max_y_vals)
            
#             if n == 0:
#                 handles, labels = ax.get_legend_handles_labels()
#                 ax.set(ylabel=ylabel)
#             # Sum up the rows of our data to get the total value of each bar.
#             totals = df_stage.sum(axis=1)
#             if not stack100 or totals_lab:
#                 # Set an offset that is used to bump the label up a bit above the bar.
#                 if stack100:
#                     y_offset = 2
#                 else: 
#                     y_offset = round(0.2*min_val)
#                 # Add labels to each bar.
#                 for i, total in enumerate(totals):
#                     if total < 100:
#                           ax.text(i, total + y_offset, total, ha='center',
#                                   weight='bold', size =10)
#                     else: 
#                           ax.text(i, total + y_offset ,locale.format("%d", total, grouping=True), ha='center',
#                                   weight='bold', size =10)
              
#             # Let's put the annotations inside the bars themselves by using a
#             # negative offset.
#             if sub_bar_lab: 
#                 # if stack100:
#                 #     y_offset = -5
#                 # else: 
#                 y_offset = round(-0.4*min_val)
#                 print(y_offset, min_val)
#                 # For each patch (basically each rectangle within the bar), add a label.
#                 for bar in ax.patches:
#                     if bar.get_height() < 100:
#                         bar_height = "{0:.1f}".format(bar.get_height())
#                     else: 
#                         bar_height = locale.format("%d",  bar.get_height(), grouping=True)
#                     ax.text(
#                           # Put the text in the middle of each bar. get_x returns the start
#                           # so we add half the width to get to the middle.
#                           bar.get_x() + bar.get_width() / 2,
#                           # Vertically, add the height of the bar to the start of the bar,
#                           # along with the offset.
#                           bar.get_height()//2 + bar.get_y()+ y_offset,
#                           # This is actual value we'll show. original: round(bar.get_height())
#                           bar_height,
#                           # Center the labels and style them a bit.
#                           ha='center', color='w', weight='bold', size=10)
                  
#             ax.set_title('Stage: '+stg+'hpf')
#             ax.spines['right'].set_visible(False)
#             ax.spines['top'].set_visible(False)
            
#         if n == len(stages)-1:
#             ax.set_axis_off()
#             ax.set_title('')
#             legends = def_var_names('bar_plots', labels)
#             ax.legend(handles, legends, loc='upper left', bbox_to_anchor=(-0.5, 1), frameon = False, fontsize=10)
#             print(labels)
            
#     if not stack100:
#         # print(y_vals_all)
#         # print(max_y_vals)
#         max_y_index = max_y_vals.index(max(max_y_vals))
#         # print(max_y_index)
#         ax.set_yticks(y_vals_all[max_y_index])
#         ax.set_yticklabels([locale.format("%d", x, grouping=True) for x in y_vals_all[max_y_index]])
#         # ax.ticklabel_format(axis='y', style='plain')

#     fig.suptitle(title+txt_title+'\n', fontsize = 12, y=1.01)
    
#%% func - barPlotsBU_lastVersion
# def barPlotsBU_lastVersion(df2plot, vars2plot, group_vars, colours, title, txt_title, ylabel, dir2save, info='', stack100 = True, 
#                          sub_bar_lab = True, save = True, ext = ['png']): 
    
#     print('\n> '+title)
#     # Create list to contain info of bars and ist highest value
#     y_vals_all = []
#     max_y_vals = []
    
#     # Filter dataframe 
#     df4plot = df2plot.groupby(group_vars)[vars2plot].mean()
#     df4plot_count = df2plot.groupby(group_vars)[vars2plot].count()
#     # Count and define figure settings accordingly
#     # - number of stages define legend position 
#     stages = sorted(df2plot['Stage'].unique())
#     n_stages_o = len(stages)
#     for nn in range(n_stages_o): stages.append('')
#     if n_stages_o >= 3:
#         leg_pos = 4
#         loc='center'
#     else: 
#         leg_pos = n_stages_o
#         loc='center left'
        
#     # - number of genotypes
#     n_genot = len(df2plot['GenotypeAll'].unique())
    
#     # - number of variables/ stacked bars define number of columns in legend
#     if len(vars2plot) == 3:
#         ncols_leg = 3
#     else: 
#         ncols_leg = 2
    
#     # If bar plots are stacked to 100 transform data to percentages
#     if stack100: 
#         df4plot = df4plot.apply(lambda x: x*100/sum(x), axis=1)

#     # Find minimum bar height to later define offset  
#     min_val = df4plot.min().min()
    
#     #Create figure
#     gridkw = dict(width_ratios=[1,1,1], height_ratios=[1,0.3])
#     fig_size = (math.ceil(2.5*(n_genot/2)*(n_stages_o)),7)
#     fig, axes = plt.subplots(ncols=n_stages_o,nrows=2, figsize=fig_size, sharey=True, gridspec_kw=gridkw)
#     #print('Figure size: ', fig_size)
#     for n, ax, stg in zip(count(), axes.flatten(), stages):
#         # First row will contain the graphs
#         if n < n_stages_o:
#             df_stage = df4plot.iloc[df4plot.index.get_level_values('Stage') == stg]
#             df_stage = df_stage.droplevel('Stage', axis="index").sort_index(key=lambda x: x.str.lower(), ascending=False)
#             df_stage_count = df4plot_count.iloc[df4plot.index.get_level_values('Stage') == stg]
#             df_stage_count = df_stage_count.mean(axis=1).droplevel('Stage', axis="index").sort_index(key=lambda x: x.str.lower(), ascending=False)
#             if df_stage.index.nlevels > 1:
#                 joined_titles = df_stage.index.map(('\n'.join))
#                 df2plot_titles = [jt+'\nn='+str(int(num)) for i, jt, num in zip(count(), joined_titles, df_stage_count)]
#             else: 
#                 df_stage_copy = df_stage.copy()
#                 gen_legend, _, strains_o_legend, _ = def_legends(df_stage_copy.reset_index())
#                 df2plot_titles = [gen+'\nn='+str(int(num)) for i, gen, num in zip(count(), gen_legend, df_stage_count)]
#             print(df_stage, df2plot_titles)
#             # Initialize the bottom at zero for the first set of bars.
#             bottom = np.zeros(len(df_stage))
            
#             # Plot each layer of the bar, adding each bar to the "bottom" so
#             # the next bar starts higher.
#             for i, col in enumerate(df_stage.columns):
#               ax.bar(df2plot_titles, df_stage[col], bottom=bottom, label=col, color=colours[i])
#               bottom += np.array(df_stage[col])
#               ax.set_xticks(df2plot_titles)
#               ax.set_xticklabels(df2plot_titles, rotation=30)
#             y_vals = ax.get_yticks()
#             y_vals_all.append(y_vals)
#             max_y_vals.append(y_vals[-1])
            
#             if n == 0:
#                 handles, labels = ax.get_legend_handles_labels()
#                 ax.set(ylabel=ylabel)
#             # Sum up the rows of our data to get the total value of each bar.
#             totals = df_stage.sum(axis=1)
#             if not stack100:
#                 # Set an offset that is used to bump the label up a bit above the bar.
#                 if stack100:
#                     y_offset = 2
#                 else: 
#                     y_offset = round(0.2*min_val)
#                 # Add labels to each bar.
#                 for i, total in enumerate(totals):
#                     if total < 100:
#                           ax.text(i, total + y_offset, total, ha='center', weight='bold', size =10)
#                     else: 
#                           ax.text(i, total + y_offset ,locale.format("%d", total, grouping=True), ha='center',
#                                   weight='bold', size =10)
              
#             # Let's put the annotations inside the bars themselves by using a
#             # negative offset.
#             if sub_bar_lab: 
#                 if stack100:
#                     y_offset = -0.5
#                 else: 
#                     y_offset = round(-0.4*min_val)
#                 # For each patch (basically each rectangle within the bar), add a label.
#                 for bar in ax.patches:
#                     if bar.get_height() < 100:
#                         bar_height = "{0:.1f}".format(bar.get_height())
#                     else: 
#                         bar_height = locale.format("%d",  bar.get_height(), grouping=True)
#                     ax.text(
#                           # Put the text in the middle of each bar. get_x returns the start
#                           # so we add half the width to get to the middle.
#                           bar.get_x() + bar.get_width() / 2,
#                           # Vertically, add the height of the bar to the start of the bar,
#                           # along with the offset.
#                           bar.get_height()//2 + bar.get_y()+ y_offset,
#                           # This is actual value we'll show. original: round(bar.get_height())
#                           bar_height,
#                           # Center the labels and style them a bit.
#                           ha='center', color='w', weight='bold', size=10)
                  
#             ax.set_title('Stage: '+stg+'hpf')
#             ax.spines['right'].set_visible(False)
#             ax.spines['top'].set_visible(False)
        
#         # Second row will contain the legend
#         if n >= n_stages_o:
#             ax.set_axis_off()
#             ax.set_title('')
#             if n== leg_pos:
#                 legends = def_var_names('bar_plots', labels)
#                 ax.legend(handles, legends, loc=loc, frameon = True, fontsize=10, ncol = ncols_leg)
    
#     # Define axes based on higherst bar
#     if not stack100:
#         max_y_index = max_y_vals.index(max(max_y_vals))
#         ax.set_yticks(y_vals_all[max_y_index])
#         ax.set_yticklabels([locale.format("%d", x, grouping=True) for x in y_vals_all[max_y_index]])
#         # ax.ticklabel_format(axis='y', style='plain')

#     fig.suptitle(title+txt_title+'\n', fontsize = 12, y=1.01)
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'meas_all', 'R_')
#             title2save = title.replace(' per Stage and Genotype','')
#             title2save = title2save.replace(' ','_')
#             if stack100: 
#                 title2save = 'Perc'+title2save
#             if info != '':
#                 fig_title = dir2savef+title2save+"_"+info+"."+extf
#             else: 
#                 fig_title = dir2savef+title2save+"."+extf
#             plt.savefig(fig_title, dpi=300, bbox_inches='tight', transparent=True)
            
# #-----
#         # if len(df2plot.GenotypeAll.unique()) > 1: 
#         #     fcAn.barPlotsBU_lastVersion(df2plot = df2plot, vars2plot = vars2plot, group_vars = group_vars, colours = colours, 
#         #               title = title_sp +'Tissue Composition per Stage and Genotype', txt_title = txt_title,
#         #               ylabel = title_sp +'Tissue Composition (%)', dir2save = dir_pl_meas,
#         #               info = info, stack100 = True, sub_bar_lab = True, save = save)
#         #     fcAn.barPlots(df2plot = df2plot, vars2plot = vars2plot, group_vars = group_vars, colours = colours, 
#         #               title = title_sp +'Tissue Composition per Stage and Genotype', txt_title = txt_title,
#         #               ylabel = title_sp +'Tissue Composition [um$^3$]', dir2save = dir_pl_meas,
#         #               info = info, stack100 = False, sub_bar_lab = True, save = save)
#         #     fcAn.barPlots(df2plot = df2plot, vars2plot = vars2plot + [add_var],group_vars = group_vars, colours = colours+['tomato'], 
#         #               title = title_sp +'Composition per Stage and Genotype', txt_title = txt_title,
#         #               ylabel = title_sp +'Composition (%)', dir2save = dir_pl_meas,
#         #               info = info, stack100 = True, sub_bar_lab = True, save = save)
#         #     fcAn.barPlots(df2plot = df2plot, vars2plot = vars2plot + [add_var],group_vars = group_vars, colours = colours+['tomato'], 
#         #               title = title_sp +'Composition per Stage and Genotype', txt_title = txt_title,
#         #               ylabel = title_sp +'Composition [um$^3$]', dir2save = dir_pl_meas,
#         #               info = info, stack100 = False, sub_bar_lab = True, save = save)
            
#%% func - pctChange_barPlotsBU
# def pctChange_barPlotsBU(df2plot, vars2plot, group_vars, colours, title, txt_title, ylabel, dir2save, info='',
#                          sub_bar_lab = True, save = True, ext = ['png']): 
    
#     print('\n> '+title.replace('\n',' '))
#     # Create list to contain info of bars and ist highest value
#     y_vals_all = []
#     max_y_vals = []
    
#     # Filter dataframe 
#     df4plot = df2plot.groupby(group_vars)[vars2plot].mean()

#     # Count and define figure settings accordingly
#     stages = sorted(df2plot.index.unique(level = 'Stage'))
#     n_stages_o = len(stages)
        
#     # - number of genotypes
#     # - number of genotypes define legend position 
#     genot = sorted(df2plot.index.unique(level = 'GenotypeAll'), reverse = True)
#     gen_legend, _, _, _ = def_legends(df2plot.reset_index(), df_type = 'changes')
#     n_genot = len(df2plot.index.unique(level = 'GenotypeAll'))
#     for nn in range(n_genot): genot.append('')
#     if n_genot >= 2:
#         leg_pos = 2
#         loc='center'
#         bbox_to_anchor=(1, 0.5)
#     else: 
#         leg_pos = n_genot
#         loc='center'
#         bbox_to_anchor=(0.5, 0.5)
        
#     # - number of variables/ stacked bars define number of columns in legend
#     if len(vars2plot) == 3:
#         ncols_leg = 3
#     else: 
#         ncols_leg = 2
    
#     #Create figure
#     gridkw = dict(width_ratios=[1]*n_genot, height_ratios=[1,0.3])
#     fig_size = (math.ceil(2.5*(n_genot/2)*(n_stages_o)),7)
#     fig, axes = plt.subplots(ncols=n_genot,nrows=2, figsize=fig_size, sharey=True, gridspec_kw=gridkw)
#     #print('Figure size: ', fig_size)
#     for n, ax, gen in zip(count(), axes.flatten(), genot):
#         # First row will contain the graphs
#         if n < n_genot:
#             df_genot = df4plot.iloc[df4plot.index.get_level_values('GenotypeAll') == gen]
#             df_genot = df_genot.droplevel('GenotypeAll', axis="index").sort_index(key=lambda x: x.str.lower(), ascending=True)
#             if df_genot.index.nlevels > 1:
#                 df2plot_titles = df_genot.index.map(('\n'.join))
#             else: 
#                 df_genot_copy = df_genot.copy()
#                 _, _, _, stage_legend = def_legends(df_genot_copy.reset_index(), df_type = 'changes')
#                 df2plot_titles =  stage_legend
#             print('\n - '+gen +'\n'); print(df_genot)#, df2plot_titles)
#             # Initialize the bottom at zero for the first set of bars.
#             bottom_positive = np.zeros(len(df_genot)) 
#             bottom_negative = np.zeros(len(df_genot))
#             bottom = np.zeros(len(df_genot))
            
#             # Plot each layer of the bar, adding each bar to the "bottom" so
#             # the next bar starts higher.
#             for i, col in enumerate(df_genot.columns):
#                 #print('ACAAA'); print(df_genot[col])
#                 for v, val in enumerate(df_genot[col]):
#                     #print(val)
#                     if val > 0: 
#                         bottom[v] = bottom_positive[v]
#                         bottom_positive[v] += val
#                     else: 
#                         bottom[v] = bottom_negative[v]
#                         bottom_negative[v] += val
#                 #print('bottom', bottom)
#                 ax.bar(df2plot_titles, df_genot[col], bottom=bottom, label=col, color=colours[i])
#                 bottom += np.array(df_genot[col])
#                 ax.set_xticks(df2plot_titles)
#                 ax.set_xticklabels(df2plot_titles, rotation=30)
#             y_vals = ax.get_yticks()
#             y_vals_all.append(y_vals)
#             max_y_vals.append(y_vals[-1])
            
#             if n == 0:
#                 handles, labels = ax.get_legend_handles_labels()
#                 ax.set(ylabel=ylabel)
            
#             # # Sum up the rows of our data to get the total value of each bar.
#             # totals = df_genot.sum(axis=1)

#             # # Set an offset that is used to bump the label up a bit above the bar.
#             # # Add labels to each bar.
#             # for i, total in enumerate(totals):
#             #     if total > 0:
#             #         ax.text(i, total + 2,"{0:.1f}".format(total), ha='center', weight='bold', size =10)
#             #     else: 
#             #         ax.text(i, total - 7,"{0:.1f}".format(total), ha='center', weight='bold', size =10)

#             # Let's put the annotations inside the bars themselves by using a
#             # negative offset.
#             if sub_bar_lab:
#                 y_offset = -0.5
#                 # For each patch (basically each rectangle within the bar), add a label.
#                 for bar in ax.patches:
#                     if bar.get_height() < 100:
#                         bar_height = "{0:.1f}".format(bar.get_height())
#                     else: 
#                         bar_height = locale.format("%d",  bar.get_height(), grouping=True)
#                     ax.text(
#                           # Put the text in the middle of each bar. get_x returns the start
#                           # so we add half the width to get to the middle.
#                           bar.get_x() + bar.get_width() / 2,
#                           # Vertically, add the height of the bar to the start of the bar,
#                           # along with the offset.
#                           bar.get_height()//2 + bar.get_y()+ y_offset,
#                           # This is actual value we'll show. original: round(bar.get_height())
#                           bar_height,
#                           # Center the labels and style them a bit.
#                           ha='center', color='w', weight='bold', size=10)
            
#             ax.set_title(gen_legend[n])
#             ax.spines['right'].set_visible(False)
#             ax.spines['top'].set_visible(False)
#             ax.axhline(0, color='dimgrey', linewidth=0.8)
#             # ax.spines['bottom'].set_visible(False)
        
#         # Second row will contain the legend
#         if n >= n_genot:
#             ax.set_axis_off()
#             ax.set_title('')
#             if n== leg_pos:
#                 legends = def_var_names('bar_plots', labels)
#                 ax.legend(handles, legends, loc=loc, bbox_to_anchor = bbox_to_anchor, frameon = True, fontsize=10, ncol = ncols_leg)
    
#     # Define axes based on higherst bar
#     max_y_index = max_y_vals.index(max(max_y_vals))
#     ax.set_yticks(y_vals_all[max_y_index])
#     ax.set_yticklabels([locale.format("%d", x, grouping=True) for x in y_vals_all[max_y_index]])
#     # ax.ticklabel_format(axis='y', style='plain')

#     fig.suptitle(title+txt_title+'\n', fontsize = 11, y=1.01)
    
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'meas_all', 'R_')
#             title2save = title.replace(' per Genotype and Stage','')
#             title2save = title2save.replace('\n','_')
#             title2save = title2save.replace(' ','_')
#             if info != '':
#                 fig_title = dir2savef+title2save+"_"+info+"."+extf
#             else: 
#                 fig_title = dir2savef+title2save+"."+extf
#             # print(fig_title)
#             plt.savefig(fig_title, dpi=300, bbox_inches='tight', transparent=True)
            
            
    # # Options splitting stages rather than genotypes ----------------
    # print('> '+title)
    # # Create list to contain info of bars and ist highest value
    # y_vals_all = []
    # max_y_vals = []
    
    # # Filter dataframe 
    # df4plot = df2plot.groupby(group_vars)[vars2plot].mean()

    # # Count and define figure settings accordingly
    # # - number of stages define legend position 
    # stages = sorted(df2plot.index.unique(level = 'Stage'))
    # n_stages_o = len(stages)
    # for nn in range(n_stages_o): stages.append('')
    # if n_stages_o >= 2:
    #     leg_pos = 3
    #     loc='center'
    # else: 
    #     leg_pos = n_stages_o
    #     loc='center left'
        
    # # - number of genotypes
    # n_genot = len(df2plot.index.unique(level = 'GenotypeAll'))
    
    # # - number of variables/ stacked bars define number of columns in legend
    # if len(vars2plot) == 3:
    #     ncols_leg = 1
    # else: 
    #     ncols_leg = 2

    # # Find minimum bar height to later define offset  
    # min_val = df4plot.abs().min().min()
    
    # #Create figure
    # gridkw = dict(width_ratios=[1,1], height_ratios=[1,0.3])
    # fig_size = (math.ceil(2.5*(n_genot/2)*(n_stages_o)),7)

    # fig, axes = plt.subplots(ncols=n_stages_o,nrows=2, figsize=fig_size, sharey=True, gridspec_kw=gridkw)
    # #print('Figure size: ', fig_size)
    # for n, ax, stg in zip(count(), axes.flatten(), stages):
    #     # First row will contain the graphs
    #     if n < n_stages_o:
    #         df_stage = df4plot.iloc[df4plot.index.get_level_values('Stage') == stg]
    #         df_stage = df_stage.droplevel('Stage', axis="index").sort_index(key=lambda x: x.str.lower(), ascending=False)
    #         if df_stage.index.nlevels > 1:
    #             df2plot_titles = df_stage.index.map(('\n'.join))
    #         else: 
    #             df_stage_copy = df_stage.copy()
    #             gen_legend, _, strains_o_legend, _ = def_legends(df_stage_copy.reset_index())
    #             df2plot_titles =  gen_legend

    #         # Initialize the bottom at zero for the first set of bars.
    #         bottom = np.zeros(len(df_stage))
            
    #         # Plot each layer of the bar, adding each bar to the "bottom" so
    #         # the next bar starts higher.
    #         for i, col in enumerate(df_stage.columns):
    #           ax.bar(df2plot_titles, df_stage[col], bottom=bottom, label=col, color=colours[i])
    #           bottom += np.array(df_stage[col])
    #           ax.set_xticks(df2plot_titles)
    #           ax.set_xticklabels(df2plot_titles, rotation=30)
    #         y_vals = ax.get_yticks()
    #         y_vals_all.append(y_vals)
    #         max_y_vals.append(y_vals[-1])
            
    #         if n == 0:
    #             handles, labels = ax.get_legend_handles_labels()
    #             ax.set(ylabel=ylabel)
    #         # Sum up the rows of our data to get the total value of each bar.
    #         totals = df_stage.sum(axis=1)

    #         # Set an offset that is used to bump the label up a bit above the bar.
    #         # Add labels to each bar.
    #         for i, total in enumerate(totals):
    #             if total > 0:
    #                 ax.text(i, total + 2,"{0:.1f}".format(total), ha='center', weight='bold', size =10)
    #             else: 
    #                 ax.text(i, total - 7,"{0:.1f}".format(total), ha='center', weight='bold', size =10)

    #         # Let's put the annotations inside the bars themselves by using a
    #         # negative offset.
    #         if sub_bar_lab:
    #             y_offset = -0.5
    #             # For each patch (basically each rectangle within the bar), add a label.
    #             for bar in ax.patches:
    #                 if bar.get_height() < 100:
    #                     bar_height = "{0:.1f}".format(bar.get_height())
    #                 else: 
    #                     bar_height = locale.format("%d",  bar.get_height(), grouping=True)
    #                 ax.text(
    #                       # Put the text in the middle of each bar. get_x returns the start
    #                       # so we add half the width to get to the middle.
    #                       bar.get_x() + bar.get_width() / 2,
    #                       # Vertically, add the height of the bar to the start of the bar,
    #                       # along with the offset.
    #                       bar.get_height()//2 + bar.get_y()+ y_offset,
    #                       # This is actual value we'll show. original: round(bar.get_height())
    #                       bar_height,
    #                       # Center the labels and style them a bit.
    #                       ha='center', color='w', weight='bold', size=10)
                  
    #         ax.set_title('Stage: '+stg+'hpf')
    #         ax.spines['right'].set_visible(False)
    #         ax.spines['top'].set_visible(False)
        
    #     # Second row will contain the legend
    #     if n >= n_stages_o:
    #         ax.set_axis_off()
    #         ax.set_title('')
    #         if n== leg_pos:
    #             legends = def_var_names('bar_plots', labels)
    #             ax.legend(handles, legends, loc=loc, frameon = True, fontsize=10, ncol = ncols_leg)
    
    # # Define axes based on higherst bar
    # if not stack100:
    #     max_y_index = max_y_vals.index(max(max_y_vals))
    #     ax.set_yticks(y_vals_all[max_y_index])
    #     ax.set_yticklabels([locale.format("%d", x, grouping=True) for x in y_vals_all[max_y_index]])
    #     # ax.ticklabel_format(axis='y', style='plain')

    # fig.suptitle(title+txt_title+'\n', fontsize = 12, y=1.01)

#%% func - barPlots_oneGenot
def barPlots_oneGenot(df2plot, vars2plot, group_vars, colours, title, txt_title, ylabel, strain, dir2save, info='', stack100 = True, 
                        sub_bar_lab = True, save = True, ext = ['png']): 
    
    print('\n>>> '+title)
    # Create list to contain info of bars and ist highest value
    y_vals_all = []
    max_y_vals = []
    
    # Filter dataframe 
    df4plot = df2plot.groupby(group_vars)[vars2plot].mean()
    df4plot_count = df2plot.groupby(group_vars)[vars2plot].count()
    
    # Count and define figure settings accordingly
    # - number of stages define legend position 
    stages = sorted(df2plot['Stage'].unique())
    n_stages_o = len(stages)
    for nn in range(n_stages_o): stages.append('')
    if n_stages_o >= 3:
        leg_pos = 4
        loc='center'
    else: 
        leg_pos = n_stages_o
        loc='center left'
    
    # - number of strains define subplots width
    if strain != '': 
        n_strain = len(df2plot[strain].unique())
    else:
        n_strain = 1.5
        
    # - number of variables/ stacked bars define number of columns in legend
    if len(vars2plot) == 3:
        ncols_leg = 3
    else: 
        ncols_leg = 2
        
    # If bar plots are stacked to 100 transform data to percentages
    if stack100: 
        df4plot = df4plot.apply(lambda x: x*100/sum(x), axis=1)
    
    # Find minimum bar height to later define offset  
    min_val = df4plot.min().min()

    # Create figure
    gridkw = dict(width_ratios=[1,1,1], height_ratios=[1,0.3])
    fig_size = (math.ceil(2.5*(n_strain/2)*(n_stages_o)),8)
    fig, axes = plt.subplots(ncols=n_stages_o,nrows=2, figsize=fig_size, sharey=True, gridspec_kw=gridkw)
    #print('Figure size: ', fig_size)
    for n, ax, stg in zip(count(), axes.flatten(), stages):
        # First row will contain the graphs
        if n < n_stages_o:
            df_stage = df4plot.iloc[df4plot.index.get_level_values('Stage') == stg]
            df_stage_count = df4plot_count.iloc[df4plot.index.get_level_values('Stage') == stg]
            if df_stage.index.nlevels > 1:
                df_stage = df_stage.droplevel('Stage', axis="index").sort_index(key=lambda x: x.str.lower(), ascending=False)
                df_stage_count = df_stage_count.mean(axis=1).droplevel('Stage', axis="index").sort_index(key=lambda x: x.str.lower(), ascending=False)
                if df_stage.index.nlevels > 1:
                    joined_titles = df_stage.index.map(('\n'.join))
                    df2plot_titles = [jt+'\nn='+str(int(num)) for i, jt, num in zip(count(), joined_titles, df_stage_count)]
                else: 
                    df2plot_titles = df_stage.index
            else: 
                df_stage_copy = df_stage.copy()
                _, _, strains_o_legend, _ = def_legends(df_stage_copy.reset_index())
                df_stage_count = df_stage_count.mean(axis=1)
                if strains_o_legend == []:
                    df2plot_titles = [strain+'\nn='+str(int(num)) for i, strain, num in zip(count(), [info], df_stage_count)] 
                else: 
                    df2plot_titles = [strain+'\nn='+str(int(num)) for i, strain, num in zip(count(), strains_o_legend, df_stage_count)]
                    
            # Initialize the bottom at zero for the first set of bars.
            bottom = np.zeros(len(df_stage))
            
            # Plot each layer of the bar, adding each bar to the "bottom" so
            # the next bar starts higher.
            for i, col in enumerate(df_stage.columns):
              ax.bar(df2plot_titles, df_stage[col], bottom=bottom, label=col, color=colours[i])
              bottom += np.array(df_stage[col])
              ax.set_xticks(df2plot_titles)
              ax.set_xticklabels(df2plot_titles, rotation=30)
            y_vals = ax.get_yticks()
            y_vals_all.append(y_vals)
            max_y_vals.append(y_vals[-1])
            
            if n == 0:
                handles, labels = ax.get_legend_handles_labels()
                ax.set(ylabel=ylabel)
            # Sum up the rows of our data to get the total value of each bar.
            totals = df_stage.sum(axis=1)
            if not stack100:
                # Set an offset that is used to bump the label up a bit above the bar.
                if stack100:
                    y_offset = 2
                else: 
                    y_offset = round(0.2*min_val)
                # Add labels to each bar.
                for i, total in enumerate(totals):
                    if total < 100:
                          ax.text(i, total + y_offset, total, ha='center', weight='bold', size =10)
                    else: 
                          ax.text(i, total + y_offset ,locale.format("%d", total, grouping=True), ha='center',
                                  weight='bold', size =10)
              
            # Let's put the annotations inside the bars themselves by using a
            # negative offset.
            if sub_bar_lab: 
                if stack100:
                    y_offset = -0.5
                else: 
                    y_offset = round(-0.4*min_val)
                # print(y_offset, min_val)
                # For each patch (basically each rectangle within the bar), add a label.
                for bar in ax.patches:
                    if bar.get_height() < 100:
                        bar_height = "{0:.1f}".format(bar.get_height())
                    else: 
                        bar_height = locale.format("%d",  bar.get_height(), grouping=True)
                    ax.text(
                          # Put the text in the middle of each bar. get_x returns the start
                          # so we add half the width to get to the middle.
                          bar.get_x() + bar.get_width() / 2,
                          # Vertically, add the height of the bar to the start of the bar,
                          # along with the offset.
                          bar.get_height()//2 + bar.get_y()+ y_offset,
                          # This is actual value we'll show. original: round(bar.get_height())
                          bar_height,
                          # Center the labels and style them a bit.
                          ha='center', color='w', weight='bold', size=10)
                  
            ax.set_title('Stage: '+stg+'hpf')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
        
        # Second row will contain the legend
        if n >= n_stages_o:
            ax.set_axis_off()
            ax.set_title('')
            if n== leg_pos:
                legends = def_var_names('bar_plots', labels)
                ax.legend(handles, legends, loc=loc, frameon = True, fontsize=10, ncol = ncols_leg)
    
    # Define axes based on higherst bar
    if not stack100:
        max_y_index = max_y_vals.index(max(max_y_vals))
        ax.set_yticks(y_vals_all[max_y_index])
        ax.set_yticklabels([locale.format("%d", x, grouping=True) for x in y_vals_all[max_y_index]])
        # ax.ticklabel_format(axis='y', style='plain')

    fig.suptitle(title+txt_title+'\n', fontsize = 12, y=1.01)
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'meas_all', 'R_')
            title2save = title.replace(' per Stage','')
            title2save = title2save.replace(' ','_')
            if stack100: 
                title2save = 'Perc'+title2save
            if info != '':
                fig_title = dir2savef+title2save+"_"+info+"."+extf
            else: 
                fig_title = dir2savef+title2save+"."+extf
            # print(fig_title)
            plt.savefig(fig_title, dpi=300, bbox_inches='tight', transparent=True)
            
#%% func - strainORstrain_o
def strainORstrain_o(df2plot):
    
    vars_opt = ['Stage','GenotypeAll']#,'Manip']
    strain = ''
    include_strain = ask4input('Include strain division? >>:', bool)
    if include_strain:
        strain_o = ask4input('Strain_o? >>:', bool)
        if strain_o: 
            vars_opt.append('Strain_o')
            strain = 'Strain_o'
        else: 
            vars_opt.append('Strain')
            strain = 'Strain'
    
    group_vars = []
    txt_title = '\n'
    for var in vars_opt:
        if len(df2plot[var].unique()) > 1:
            group_vars.append(var)
        else: 
            txt_title = txt_title +'['+ var + ':'+df2plot[var].unique()[0]+'] '
    
    return group_vars, txt_title, strain


#%% END