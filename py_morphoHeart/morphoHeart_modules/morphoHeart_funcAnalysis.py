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
import math
import locale
locale.setlocale(locale.LC_ALL, 'en_US')  

from scipy import stats 
import scikit_posthocs as sp
from statannot import add_stat_annotation

from progress.bar import Bar
suffix = '%(index)d/%(max)d - %(elapsed)ds'

#%% Importing morphoHeart packages
from .morphoHeart_funcBasics import alert, ask4input, getInputNumbers, loadDF, saveDF

#%% Interesting links
# Colors for graphs
# https://chartio.com/learn/charts/how-to-choose-colors-data-visualization/
# https://projects.susielu.com/viz-palette?colors=[%22#dc133b%22,%22#ffa500%22,%22#6495ed%22,%22#ade64f%22,%22#c36bea%22,%22#214d4e%22,%22#36cbd3%22,%22#a23e27%22,%22#839292%22,%22#8b008b%22,%22#3ff44c%22,%22#000080%22]&backgroundColor=%22white%22&fontColor=%22black%22&mode=%22normal%22
# http://vrl.cs.brown.edu/color
# https://www.rapidtables.com/web/color/RGB_Color.html#color-table
#

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
    
    # Add hapln1a multiple mutations to gene and genotype
    genotA = []
    for aa, geneA, genotypeA, strain in zip(count(), df_input['Gene_A'], df_input['Genotype_A'], df_input['Strain']):
        if 'prom187' in strain:
            genotA.append(geneA+'187:'+genotypeA)
        elif 'prom241' in strain:
            genotA.append(geneA+'241:'+genotypeA)
        elif '3bp' in strain:
            genotA.append(geneA+'3bp:'+genotypeA)
        elif 'STOP' in strain:
            genotA.append(geneA+'STOP:'+genotypeA)
        elif 'prom365' in strain:
            genotA.append(geneA+'365:'+genotypeA)
        elif 'myl7:lifeActGFP/+; fli1a:AcTagRFP/fli1a:AcTagRFP' in strain:
            genotA.append(geneA+'_tc:'+genotypeA)
        else: 
            genotA.append(geneA+':'+genotypeA)
    df_input['GenotA'] = genotA
    
    # Merge genotypes
    df_input['GenotB'] = df_input['Gene_B']+':'+df_input['Genotype_B']
    df_input['GenotypeAll'] = df_input['GenotA']+'/'+df_input['GenotB']
    
    # Create new column with complete genotype and another to merge all wild-types indistinctive of strain
    genotypeAll = [] 
    genotypeF = []
    for i, genotA, genotB in zip(count(), df_input["GenotA"], df_input["GenotB"]): 
        if ':wt' in genotA:
            if ':wt' in genotB:
                genotypeAll.append(genotA+'/'+genotB) 
                genotypeF.append('wt:wt')
            elif '-:-' in genotB:
                genotypeAll.append(genotA) 
                genotypeF.append('wt:wt')
            else:
                genotypeAll.append(genotA+'/'+genotB) 
                genotypeF.append(genotA+'/'+genotB) 
                
        elif '-:-' in genotB:
            genotypeAll.append(genotA) 
            genotypeF.append(genotA) 
            
        else: 
            genotypeAll.append(genotA+'/'+genotB) 
            genotypeF.append(genotA+'/'+genotB) 
        # if genotB == "-:-": 
        #     genotypeAll.append(genotA) 
        # else: 
        #     genotypeAll.append(genotA+'/'+genotB) 
        
    df_input["GenotypeAll"] = genotypeAll 
    df_input["GenotypeF"] = genotypeF
    
    # Add ref for Session and Fish/Embryo
    df_input['Ref'] = df_input['LS_Session']+'_'+df_input['Fish_ref']
    
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

#%% func - spAnalysis
def spAnalysis(df_meas):
    
    spAn = []
    for i, txt in zip(count(), df_meas['spAnalysis']):
        if 'spaw' in txt:
            if 'lf' in txt: 
                spAn.append('Left looper')
            if 'rt' in txt: 
                spAn.append('Right looper')
            if 'ct' in txt: 
                spAn.append('No looper')
        else: #if 'normal' in txt:
            spAn.append('normal')
            
    df_meas['spAnalysis'] = spAn
    
    return df_meas

#%% func - modifySibsGenot
def modifySibsGenot(df2plot):
    askQ = False
    genotF = df2plot['GenotypeF']
    for gen in genotF:
        if gen == 'hapln1a241:wt/spaw:ht':
            askQ = True
    if askQ:
        modif = ask4input('Are you sure you want to modify sibs (hapln1a241:wt,spaw:ht) genotype? \n\t[0]: no, [1]: yes! >:', bool)
        if modif: 
            newGenotF = []
            for genot in genotF:
                if genot == 'hapln1a241:wt/spaw:ht':
                    newGenotF.append('wt:wt')
                else: 
                    newGenotF.append(genot)
            df2plot['GenotypeF'] = newGenotF
    
    return df2plot

#%% func - modifyOEControlGenot
def modifyOEControlGenot(df2plot):
    askQ = False
    genotF = df2plot['GenotypeF']
    for gen in genotF:
        if gen == 'galff:-/uas:+' or gen == 'galff:+/uas:-':
            askQ = True
    if askQ:
        modif = ask4input('Are you sure you want to modify controls (GAL4:- or UAS:-) genotype? \n\t[0]: no, [1]: yes! >:', bool)
        if modif: 
            # print('AJAAA')
            newGenotF = []
            for genot in genotF:
                if '-' in genot:
                    # print(genot)
                    newGenotF.append('wt:ctrl')
                else: 
                    # print(genot, 'Aja')
                    newGenotF.append(genot)
            df2plot['GenotypeF'] = newGenotF
    
    return df2plot
    
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
def selectVariables_auto (vars_dict, group, plot_type):

    if plot_type == 'group':
        pl_groups = plot_groups()
    elif plot_type == 'indiv':
        pl_groups = plot_indiv()
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
       'Gene_A', 'Genotype_A','GenotA', 'Gene_B', 'Genotype_B', 'GenotB', 'GenotypeAll', 'GenotypeF']

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
    df_meas['Ratio_VolEndo2VolExtMyoc'] = df_meas['Vol_Endo']/ df_meas['Vol_Ext.Myoc']
    df_meas['Ratio_VolCJ2VolExtMyoc'] = df_meas['Vol_CJ']/ df_meas['Vol_Ext.Myoc']
    df_meas['Ratio_VolLumen2VolExtMyoc'] = df_meas['Vol_Int.Endo']/ df_meas['Vol_Ext.Myoc']
    
    df_meas['Ratio_VolMyoc2VolTissue'] = df_meas['Vol_Myoc']/ df_meas['Vol_Tissue']
    df_meas['Ratio_VolEndo2VolTissue'] = df_meas['Vol_Endo']/ df_meas['Vol_Tissue']
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
                       "EllipAtr_Depth" : "Atrial depth [$\mu$m]",
                       "EllipAtr_Length" : "Atrial length [$\mu$m]",
                       "EllipAtr_Width" : "Atrial width [$\mu$m]",
                       "EllipVent_Asphericity" : "Ventricular Ellipse Asphericity",
                       "EllipVent_Depth" :"Ventricular depth [$\mu$m]",
                       "EllipVent_Length" : "Ventricular length [$\mu$m]",
                       "EllipVent_Width" : "Ventricular width [$\mu$m]",
                       
                       # Surface Areas
                       "SurfArea_Myoc" : "Surface Area\nMyocardium [$\mu$m$^2$]",
                       "SurfArea_Int.Myoc" : "Surface Area\nInt.Myocardium [$\mu$m$^2$]",
                       "SurfArea_Ext.Myoc" : "Surface Area\nExt.Myocardium [$\mu$m$^2$]",
                       "SurfArea_Endo" : "Surface Area\nEndocardium [$\mu$m$^2$]",
                       "SurfArea_Int.Endo" : "Surface Area\nInt.Endocardium [$\mu$m$^2$]",
                       "SurfArea_Ext.Endo" : "Surface Area\nExt. Endocardium [$\mu$m$^2$]",
                       "SurfArea_CJ" : "Surface Area\nCardiac Jelly [$\mu$m$^2$]",
                       "SurfArea_Int.CJ" : "Surface Area\nInt.Cardiac Jelly [$\mu$m$^2$]",
                       "SurfArea_Ext.CJ" : "Surface Area\nExt.Cardiac Jelly [$\mu$m$^2$]",
                       "SurfArea_Atr.ExtCJ" : "Atrial Surface Area\nExternal Cardiac Jelly [$\mu$m$^2$]",
                       "SurfArea_Atr.ExtMyoc" : "Atrial Heart\nSurface Area [$\mu$m$^2$]",
                       "SurfArea_Atr.IntEndo" : "Atrial Lumen\nSurface Area [$\mu$m$^2$]",
                       "SurfArea_Vent.ExtCJ" : "Ventricular Surface Area\nExternal Cardiac Jelly [$\mu$m$^2$]",
                       "SurfArea_Vent.ExtMyoc" : "Ventricular Heart\nSurface Area [$\mu$m$^2$]",
                       "SurfArea_Vent.IntEndo" : "Ventricular Lumen\nSurface Area [$\mu$m$^2$]",
                        
                       # Looping ratio
                       "linLine_Int.Myoc(Cut)" : "Linear Heart Length [$\mu$m]",
                       "linLine_Ext.Endo(Cut)" : "Linear Heart Length\n(Ext.Endo) [$\mu$m]",
                       "Length_CL_Int.Myoc(Cut)" : "Looped Heart Length [$\mu$m]",
                       "Length_CL_Ext.Endo(Cut)" : "Looped Heart Length\n(Ext.Endo) [$\mu$m]",
                       'Looping_Ratio_Myoc' : 'Looping Ratio',
                       'Looping_Ratio_Endo' : 'Looping Ratio (Ext.Endo)',
                        
                       # Volumes
                       "Vol_Int.Myoc" : "Volume Int.Myocardium [$\mu$m$^3$]",
                       "Vol_Ext.Myoc" : "Heart Volume [$\mu$m$^3$]",
                       "Vol_Int.Endo" : "Heart Lumen Volume [$\mu$m$^3$]",
                       "Vol_Ext.Endo" : "Volume Ext.Endocardium [$\mu$m$^3$]",
                       "Vol_Myoc" : "Volume Myocardium [$\mu$m$^3$]",
                       "Vol_Atr.Myoc" : "Atrial Volume\nMyocardium [$\mu$m$^3$]",
                       "Vol_Vent.Myoc" : "Ventricular Volume\nMyocardium [$\mu$m$^3$]",
                       "Vol_Endo" : "Volume Endocardium [$\mu$m$^3$]",
                       "Vol_Atr.Endo" : "Atrial Volume\nEndocardium [$\mu$m$^3$]",
                       "Vol_Vent.Endo" : "Ventricular Volume\nEndocardium [$\mu$m$^3$]",
                       "Vol_CJ" : "Volume Cardiac Jelly [$\mu$m$^3$]",
                       "Vol_Atr.CJ" : "Atrial Volume\nCardiac Jelly [$\mu$m$^3$]",
                       "Vol_Vent.CJ" : "Ventricular Volume\nCardiac Jelly [$\mu$m$^3$]",
                       'Vol_Atr.ExtMyoc' : 'Atrial Volume [$\mu$m$^3$]',
                       'Vol_Vent.ExtMyoc' : 'Ventricular Volume [$\mu$m$^3$]',
                       'Vol_Atr.IntEndo' : 'Atrial Lumen Volume [$\mu$m$^3$]',
                       'Vol_Vent.IntEndo' : 'Ventricular Lumen Volume [$\mu$m$^3$]',
                       
                       'Vol_CJ.Left' : "Left Volume\nCardiac Jelly [$\mu$m$^3$]",
                       'Vol_CJ.Right' : "Right Volume\nCardiac Jelly [$\mu$m$^3$]",
                       
                       'Vol_Tissue' : 'Tissue Volume [$\mu$m$^3$]',
                       'Vol_Atr.Tissue' : 'Atrial Tissue Volume [$\mu$m$^3$]', 
                       'Vol_Vent.Tissue' : 'Ventricular Tissue Volume [$\mu$m$^3$]',
                       
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
        
        vars_dict = {"Vol_Int.Myoc" : "Volume Int.Myocardium [$\mu$m$^3$]",
                       "Vol_Ext.Myoc" : "Heart\nVolume ",
                       "Vol_Int.Endo" : "Heart\nLumen",
                       "Vol_Ext.Endo" : "Volume Ext.Endocardium [$\mu$m$^3$]",
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
                      'n_cols': 4, 'yticks_lab':'d.', 'ylim' : ''}
                     }
    
    return pl_groups

#%% func - plot_indiv
def plot_indiv():
    pl_indiv = {
                'Vol_Ext.Myoc': 
                    {'graph_no': '01', 'title': 'Heart Size', 
                     'vars' : ['Vol_Ext.Myoc'],
                     'n_cols': 1, 'yticks_lab':'1e6 - d.', 
                     'ylim' : (0.75e6,3e6), 'yset' : 'round'},
                'Vol_Atr.ExtMyoc': 
                    {'graph_no': '02', 'title': 'Atrial Size', 
                     'vars' : ['Vol_Atr.ExtMyoc'],
                     'n_cols': 1, 'yticks_lab':'1e6 - d.', 
                     'ylim' : (0.2e6,2e6), 'yset' : 'round'},
                'Vol_Vent.ExtMyoc': 
                    {'graph_no': '03', 'title': 'Ventricular Size', 
                     'vars' : ['Vol_Vent.ExtMyoc'],
                     'n_cols': 1, 'yticks_lab':'1e6 - d.', 
                     'ylim' : (0.2e6,2e6), 'yset' : 'round'},
                    
                'Vol_Int.Endo': 
                    {'graph_no': '04', 'title': 'Lumen Size', 
                     'vars' : ['Vol_Int.Endo'],
                     'n_cols': 1, 'yticks_lab':'1e6 - d.', 
                     'ylim' :  (0,1.6e6), 'yset' : 'round'},
                'Vol_Atr.IntEndo': 
                    {'graph_no': '05', 'title': 'Atrial Lumen Size', 
                     'vars' : ['Vol_Atr.IntEndo'],
                     'n_cols': 1, 'yticks_lab':'1e6 - d.', 
                     'ylim' :  (0,1.2e6), 'yset' : 'round'},
                'Vol_Vent.IntEndo': 
                    {'graph_no': '06', 'title': 'Ventricular Lumen Size', 
                     'vars' : ['Vol_Vent.IntEndo'],
                     'n_cols': 1, 'yticks_lab':'1e6 - d.', 
                     'ylim' : (0,0.6e6), 'yset' : 'round'},#(0,1.2e6)
                     # 'ylim' : (0,1.2e6), 'yset' : 'round'},
                    
                'linLine_Int.Myoc(Cut)': 
                    {'graph_no': '07', 'title': 'Linear Heart Length', 
                     'vars' : ['linLine_Int.Myoc(Cut)'],
                     'n_cols': 1, 'yticks_lab':'d. - 0', 
                     'ylim' : (100,280), 'yset' : 'round'},
                'Length_CL_Int.Myoc(Cut)': 
                    {'graph_no': '08', 'title': 'Looped Heart Length', 
                     'vars' : ['Length_CL_Int.Myoc(Cut)'],
                     'n_cols': 1, 'yticks_lab':'d. - 0', 
                     'ylim' : (180,360), 'yset' : 'round'},
                'Looping_Ratio_Myoc': 
                    {'graph_no': '09', 'title': 'Looping Ratio', 
                     'vars' : ['Looping_Ratio_Myoc'],
                     'n_cols': 1, 'yticks_lab':'d.', 
                     'ylim' : (1,2.6), 'yset' : 'dec'},
                    
                'Vol_Myoc': 
                    {'graph_no': '10', 'title': 'Myocardial Volume', 
                     'vars' : ['Vol_Myoc'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.',
                     'ylim' : (250e3,750e3), 'yset' : 'round'},
                'Vol_Atr.Myoc': 
                    {'graph_no': '11', 'title': 'Atrial Myocardium Volume', 
                     'vars' : ['Vol_Atr.Myoc'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' : (100e3,500e3), 'yset' : 'round'},
                'Vol_Vent.Myoc': 
                    {'graph_no': '12', 'title': 'Ventricular Myocardium Volume', 
                     'vars' : ['Vol_Vent.Myoc'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' : (100e3,500e3), 'yset' : 'round'},
                
                'Vol_Endo': 
                    {'graph_no': '13', 'title': 'Endocardial Volume', 
                     'vars' : ['Vol_Endo'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' : (100e3,450e3), 'yset' : 'round'},
                'Vol_Atr.Endo': 
                    {'graph_no': '14', 'title': 'Atrial Endocardium Volume', 
                     'vars' : ['Vol_Atr.Endo'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' : (50e3,260e3), 'yset' : 'round'},
                'Vol_Vent.Endo': 
                    {'graph_no': '15', 'title': 'Ventricular Endocardium Volume', 
                     'vars' : ['Vol_Vent.Endo'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' : (50e3,260e3), 'yset' : 'round'},
                
                'Vol_CJ': 
                    {'graph_no': '16', 'title': 'Cardiac Jelly Volume', 
                     'vars' : ['Vol_CJ'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' : (0,700e3), 'yset' : 'round'},
                'Vol_Atr.CJ': 
                    {'graph_no': '17', 'title': 'Atrial Cardiac Jelly Volume', 
                     'vars' : ['Vol_Atr.CJ'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' : (0,600e3), 'yset' : 'round'},
                'Vol_Vent.CJ': 
                    {'graph_no': '18', 'title': 'Ventricular Cardiac Jelly Volume', 
                     'vars' : ['Vol_Vent.CJ'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' :(0,140e3), 'yset' : 'round'}, # (0,600e3)
                     # 'ylim' :(0,600e3), 'yset' : 'round'},
                'Vol_CJ.Left': 
                    {'graph_no': '19', 'title': 'Left Cardiac Jelly Volume', 
                     'vars' : ['Vol_CJ.Left'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' : (0,700e3), 'yset' : 'round'},
                'Vol_CJ.Right': 
                    {'graph_no': '20', 'title': 'Right Cardiac Jelly Volume', 
                     'vars' : ['Vol_CJ.Right'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' : (0,700e3), 'yset' : 'round'},
                'Vol_CJ.LeftvsRight': 
                    {'graph_no': '21', 'title': 'Cardiac Jelly Volume leftOverRight', 
                     'vars' : ['Ratio_VolLeftCJ2VolRightCJ'],
                     'n_cols': 1, 'yticks_lab':'d.', 
                     'ylim' : '', 'yset' : 'dec'},
                'Vol_CJ.Right2': 
                    {'graph_no': '22', 'title': 'Right Cardiac Jelly Volume', 
                     'vars' : ['Vol_CJ.Right'],
                     'n_cols': 1, 'yticks_lab':'1e3 - d.', 
                     'ylim' : (0,400e3), 'yset' : 'round'},
                
                'EllipAtr_Depth': 
                    {'graph_no': '25', 'title': 'Atrial Depth', 
                     'vars' : ['EllipAtr_Depth'],
                     'n_cols': 1, 'yticks_lab':'d. - 0', 
                     'ylim' : (60,180), 'yset' : 'round'},
                'EllipAtr_Length': 
                    {'graph_no': '26', 'title': 'Atrial Length', 
                     'vars' : ['EllipAtr_Length'],
                     'n_cols': 1, 'yticks_lab':'d. - 0', 
                     'ylim' : (80,240), 'yset' : 'round'},
                'EllipAtr_Width': 
                    {'graph_no': '27', 'title': 'Atrial Width', 
                     'vars' : ['EllipAtr_Width'],
                     'n_cols': 1, 'yticks_lab':'d. - 0', 
                     'ylim' : (70,160), 'yset' : 'round'},
                'EllipAtr_Asphericity': 
                    {'graph_no': '28', 'title': 'Atrial Asphericity', 
                      'vars' : ['EllipAtr_Asphericity'],
                      'n_cols': 1, 'yticks_lab':'d.', 
                      'ylim' : (0,0.4), 'yset' : 'dec'},
                    
                'EllipVent_Depth': 
                    {'graph_no': '29', 'title': 'Ventricular Depth', 
                     'vars' : ['EllipVent_Depth'],
                     'n_cols': 1, 'yticks_lab':'d. - 0', 
                     'ylim' : (60,180), 'yset' : 'round'},
                'EllipVent_Length': 
                    {'graph_no': '30', 'title': 'Ventricular Length', 
                     'vars' : ['EllipVent_Length'],
                     'n_cols': 1, 'yticks_lab':'d. - 0', 
                     'ylim' : (80,240), 'yset' : 'round'},
                'EllipVent_Width': 
                    {'graph_no': '31', 'title': 'Ventricular Width', 
                     'vars' : ['EllipVent_Width'],
                     'n_cols': 1, 'yticks_lab':'d. - 0', 
                     'ylim' : (70,160), 'yset' : 'round'},
                'EllipVent_Asphericity': 
                    {'graph_no': '32', 'title': 'Ventricular Asphericity', 
                      'vars' : ['EllipVent_Asphericity'],
                      'n_cols': 1, 'yticks_lab':'d.', 
                      'ylim' : (0,0.4), 'yset' : 'dec'},
                    
                'ang_AtrS': 
                    {'graph_no': '40', 'title': 'Atrium Orientation (Saggital Plane)', 
                      'vars' : ['ang_AtrS'],
                      'n_cols': 1, 'yticks_lab':'d. - 0',
                      'ylim' : (0,60), 'yset' : 'round'},
                'ang_VentS': 
                    {'graph_no': '41', 'title': 'Ventricle Orientation (Saggital Plane)', 
                      'vars' : ['ang_VentS'],
                      'n_cols': 1, 'yticks_lab':'d. - 0', 
                      'ylim' : (120,190), 'yset' : 'round'},
                'ang_BtwChambersS': 
                    {'graph_no': '42', 'title': 'Angle Between Chambers  (Saggital Plane)', 
                      'vars' : ['ang_BtwChambersS'],
                      'n_cols': 1, 'yticks_lab':'d. - 0', 
                      'ylim' : (80,200), 'yset' : 'round'},
                
                'ang_AtrV': 
                    {'graph_no': '43', 'title': 'Atrium Orientation (Ventral Plane)', 
                      'vars' : ['ang_AtrV'],
                      'n_cols': 1, 'yticks_lab':'d. - 0', 
                      'ylim' : (-5,45), 'yset' : 'round'},
                'ang_VentV': 
                    {'graph_no': '44', 'title': 'Ventricle Orientation (Ventral Plane)', 
                      'vars' : ['ang_VentV'],
                      'n_cols': 1, 'yticks_lab':'d. - 0', 
                      'ylim' : (110,190), 'yset' : 'round'},
                'ang_BtwChambersV': 
                    {'graph_no': '45', 'title': 'Angle Between Chambers  (Ventral Plane)', 
                      'vars' : ['ang_BtwChambersV'],
                      'n_cols': 1, 'yticks_lab':'d. - 0', 
                      'ylim' : (100,190), 'yset' : 'round'},
                    
                'Ratio_VolMyoc2VolExtMyoc': 
                    {'graph_no': '50', 'title': 'Heart Composition_Myoc', 
                      'vars' : ['Ratio_VolMyoc2VolExtMyoc'],
                      'n_cols': 1, 'yticks_lab':'d.', 
                      'ylim' : (0,1), 'yset' : 'dec'},
                'Ratio_VolEndo2VolExtMyoc': 
                    {'graph_no': '51', 'title': 'Heart Composition_Endo', 
                      'vars' : ['Ratio_VolEndo2VolExtMyoc'],
                      'n_cols': 1, 'yticks_lab':'d.', 
                      'ylim' : (0,1), 'yset' : 'dec'},
                'Ratio_VolCJ2VolExtMyoc': 
                    {'graph_no': '52', 'title': 'Heart Composition_CJ', 
                      'vars' : ['Ratio_VolCJ2VolExtMyoc'],
                      'n_cols': 1, 'yticks_lab':'d.', 
                      'ylim' : (0,1), 'yset' : 'dec'},
                'Ratio_VolLumen2VolExtMyoc': 
                    {'graph_no': '53', 'title': 'Heart Composition_Lumen', 
                      'vars' : ['Ratio_VolLumen2VolExtMyoc'],
                      'n_cols': 1, 'yticks_lab':'d.', 
                      'ylim' : (0,1), 'yset' : 'dec'},
                    
                'Ratio_VolMyoc2VolTissue': 
                    {'graph_no': '55', 'title': 'Tissue Composition_Myoc', 
                      'vars' : ['Ratio_VolMyoc2VolTissue'],
                      'n_cols': 1, 'yticks_lab':'d.', 
                      'ylim' : (0,1), 'yset' : 'dec'},
                'Ratio_VolEndo2VolTissue': 
                    {'graph_no': '56', 'title': 'Tissue Composition_Endo', 
                      'vars' : ['Ratio_VolEndo2VolTissue'],
                      'n_cols': 1, 'yticks_lab':'d.', 
                      'ylim' : (0,1), 'yset' : 'dec'},
                'Ratio_VolCJ2VolTissue': 
                    {'graph_no': '57', 'title': 'Tissue Composition_CJ', 
                      'vars' : ['Ratio_VolCJ2VolTissue'],
                      'n_cols': 1, 'yticks_lab':'d.', 
                      'ylim' : (0,1), 'yset' : 'dec'},
                'Ratio_VolLumen2VolTissue': 
                    {'graph_no': '58', 'title': 'Tissue Composition_Lumen', 
                      'vars' : ['Ratio_VolLumen2VolTissue'],
                      'n_cols': 1, 'yticks_lab':'d.', 
                      'ylim' : (0,1), 'yset' : 'dec'},
                    
        
                # '': 
                #     {'graph_no': 01, 'title': '', 
                #       'vars' : [''],
                #       'n_cols': 1, 'yticks_lab':'1e3 - d.', 'ylim' : ''},
                # '': 
                #     {'graph_no': 01, 'title': '', 
                #       'vars' : [''],
                #       'n_cols': 1, 'yticks_lab':'1e3 - d.', 'ylim' : ''},
                # '': 
                #     {'graph_no': 01, 'title': '', 
                #       'vars' : [''],
                #       'n_cols': 1, 'yticks_lab':'1e3 - d.', 'ylim' : ''},
                

                    
                 }
    return pl_indiv
        
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
    if df_type == 'meas':
        try: 
            genots = sorted(df_input.GenotypeAll.unique(), reverse = True)
            genotsF = sorted(df_input.GenotypeF.unique(), reverse = True)
            bool_genots = True
            bool_genotsF = True
        except: 
            bool_genots = False
            bool_genotsF = False
            
    elif df_type == 'changes':
        try: 
            genots = sorted(df_input.GenotypeAll.unique(), reverse = True)
            bool_genots = True
        except: 
            bool_genots = False
        try: 
            genotsF = sorted(df_input.GenotypeF.unique(), reverse = True)
            bool_genotsF = True
        except: 
            bool_genotsF = False
            
    elif df_type == 'changesF':
        try: 
            genots = sorted(df_input.GenotypeF.unique(), reverse = True)
            bool_genotsF = True
            bool_genots = False
        except: 
            bool_genotsF = False
            bool_genots = False
    else:  # 'kde'
        try: 
            genots = sorted(df_input.Genotype.unique(), reverse = True)
            genotsF = sorted(df_input.GenotypeF.unique(), reverse = True)
            bool_genots = True
            bool_genotsF = False
        except: 
            bool_genots = False
            bool_genotsF = False
    
    # leg_genots = {'hapln1a187:wt': {'legend': '$wild$-$type_{187KO}$', 'color': '#000080'}, 
    #                   'hapln1a241:wt': {'legend': '$wild$-$type_{241KO}$', 'color': '#000080'}, 
    #                   'hapln1a3bp:wt': {'legend': '$wild$-$type_{3KO}$', 'color': '#000080'},  
    #                   'hapln1aSTOP:wt': {'legend': '$wild$-$type_{STKO}$', 'color': '#000080'}, 
                    
    #                   'hapln1a187:ht': {'legend': '$hapln1a^{+/-}_{\Delta187}$', 'color': '#839292'}, 
    #                   'hapln1a241:ht': {'legend': '$hapln1a^{+/-}_{\Delta241}$', 'color': '#839292'}, 
                      
    #                   'hapln1a187:mt': {'legend': '$hapln1a^{-/-}_{\Delta187}$', 'color': '#8B008B'}, 
    #                   'hapln1a241:mt': {'legend': '$hapln1a^{-/-}_{\Delta241}$', 'color': '#dc133b'}, 
    #                   'hapln1a3bp:mt': {'legend': '$hapln1a^{-/-}_{\Delta3bp}$', 'color': '#c36bea'}, 
    #                   'hapln1aSTOP:mt': {'legend': '$hapln1a^{-/-}_{STOP}$', 'color': '#72e5ef'}, 
                      
    #                   'hapln1a241:wt/spaw:wt': {'legend': '$wild$-$type_{sh}$', 'color': '#000080'},
    #                   'hapln1a241:wt/spaw:ht': {'legend': '$hapln1a^{+/+}_{\Delta241};spaw^{+/-}$', 'color': '#000080'},
    #                   'hapln1a241:wt/spaw:mt': {'legend': '$spaw^{-/-}$', 'color': '#36cbd3'}, 
    #                   'hapln1a241:mt/spaw:wt': {'legend': '$hapln1a^{-/-}_{\Delta241sh}$', 'color': '#dc133b'}, 
    #                   'hapln1a241:mt/spaw:mt': {'legend': '$hapln1a^{-/-}_{\Delta241};spaw^{-/-}$', 'color': '#3ff44c'},
                      
    #                   'galff:+/uas:+': {'legend': '$galff^+;uas^+$', 'color': '#ade64f'}, 
    #                   'galff:+/uas:-': {'legend': '$galff^+;uas^-$', 'color': '#6495ed'}, 
    #                   'galff:-/uas:+': {'legend': '$galff^-;uas^+$', 'color': '#214d4e'}, 
    #                   'galff:-/uas:-': {'legend': '$wild$-$type_{OE}$', 'color': '#000080'},
                      
    #                   'vcana365:wt': {'legend': '$wild$-$type_{365KO}$', 'color': '#000080'}, 
    #                   'vcana365:ht': {'legend': '$vcana^{+/-}_{\Delta365}$', 'color': '#a23e27'}, 
    #                   'vcana365:mt': {'legend': '$vcana^{-/-}_{\Delta365}$', 'color': '#ffa500'}, 
    #                   'wt:wt': {'legend': '$wild$-$type$','color': '#000080'},
    #                   'wt_tc:wt': {'legend': '$wild$-$type_{tc}$','color': '#000080'},
    #                   'wt:ctrl': {'legend': '$wild$-$type_{OE}$','color': '#000080'}}
    
    leg_genots = {'hapln1a187:wt': {'legend': '$sibling$', 'color': '#000080'}, 
                      'hapln1a241:wt': {'legend': '$sibling$', 'color': '#000080'}, 
                      'hapln1a3bp:wt': {'legend': '$sibling', 'color': '#000080'},  
                      'hapln1aSTOP:wt': {'legend': '$sibling$', 'color': '#000080'}, 
                    
                      'hapln1a187:ht': {'legend': '$hapln1a^{+/-}_{\Delta187}$', 'color': '#839292'}, 
                      'hapln1a241:ht': {'legend': '$hapln1a^{+/-}_{\Delta241}$', 'color': '#839292'}, 
                      
                      'hapln1a187:mt': {'legend': '$hapln1a^{\Delta187/\Delta187}$', 'color': '#8B008B'}, 
                      'hapln1a241:mt': {'legend': '$hapln1a^{\Delta241/\Delta241}$', 'color': '#dc133b'}, 
                      'hapln1a3bp:mt': {'legend': '$hapln1a^{\Delta3/\Delta3}$', 'color': '#c36bea'}, 
                      'hapln1aSTOP:mt': {'legend': '$hapln1a^{STOP/STOP}$', 'color': '#72e5ef'}, 
                      
                      'hapln1a241:wt/spaw:wt': {'legend': '$sibling$', 'color': '#000080'},
                      'hapln1a241:wt/spaw:ht': {'legend': '$sibling$', 'color': '#000080'},
                      'hapln1a241:wt/spaw:mt': {'legend': '$spaw^{-/-}$', 'color': '#36cbd3'}, 
                      'hapln1a241:mt/spaw:wt': {'legend': '$hapln1a^{\Delta241/\Delta241}$', 'color': '#dc133b'}, 
                      'hapln1a241:mt/spaw:mt': {'legend': '$hapln1a^{\Delta241/\Delta241};spaw^{-/-}$', 'color': '#3ff44c'},
                      
                      'galff:+/uas:+': {'legend': '$Gal4^+$; $UAS^+$', 'color': '#ade64f'}, 
                      'galff:+/uas:-': {'legend': '$Gal4^+$; $UAS^-$', 'color': '#6495ed'}, 
                      'galff:-/uas:+': {'legend': '$Gal4^-$; $UAS^+$', 'color': '#214d4e'}, 
                      'galff:-/uas:-': {'legend': '$sibling$', 'color': '#000080'},
                      
                      'vcana365:wt': {'legend': '$sibling$', 'color': '#000080'}, 
                      'vcana365:ht': {'legend': '$vcana^{\Delta365/+}$', 'color': '#a23e27'}, 
                      'vcana365:mt': {'legend': '$vcana^{\Delta365/\Delta365}$', 'color': '#ffa500'}, 
                      'wt:wt': {'legend': '$wild$-$type$','color': '#000080'},
                      'wt_tc:wt': {'legend': '$wild$-$type_{tc}$','color': '#000080'},
                      'wt:ctrl': {'legend': '$wild$-$type_{OE}$','color': '#000080'}}
    
    out_genots = dict()
    out_genotsF = dict()
    if bool_genots: 
        for gen in genots:
            out_genots[gen] = leg_genots[gen]
    if bool_genotsF: 
        for genf in genotsF:
            out_genotsF[genf] = leg_genots[genf]
    
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
                       'myl7:galFF; UAS:TFP x UAS:hapln1a, cryaa:CFP', 
                       'myl7:lifeActGFP/+; fli1a:AcTagRFP/fli1a:AcTagRFP']
        
        leg_strains = ["$hapln1a^{\Delta 187} (F2s)$", "$hapln1a^{\Delta 187} (F3s)$", 
                       "$hapln1a^{\Delta 241} (F2s)$", "$hapln1a^{\Delta 241} (F3s)$",
                       "$spaw^{+/-}; hapln1a^{\Delta 241}$", "$vcana^{\Delta 365}$",
                       'myl7:lck-GFP; kdrl:rasCherry',
                       "$hapln1a^{OE}$", "wt (G&R)"]
    
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
                       'spaw+/-; hapln1a prom241/+ InX', 'vcana prom365/+ InX',
                       'myl7lckGFP_kdrlRasCherry',
                       'myl7:galFF; UAS:TFP x UAS:hapln1a, cryaa:CFP',
                       'myl7:lifeActGFP/+; fli1a:AcTagRFP/fli1a:AcTagRFP']
        
        leg_strains_o = ["$hapln1a^{\Delta 187}$","$hapln1a^{\Delta 241}$", 
                       "$spaw^{+/-}; hapln1a^{\Delta 241}$", "$vcana^{\Delta 365}$",
                       'myl7:lck-GFP; kdrl:rasCherry',
                       "$hapln1a^{OE}$", "wt (G&R)"]
        
        for strain_o in strains_o:
            out_strains_o.append(leg_strains_o[all_strains_o.index(strain_o)])
    
    # STAGES
    try: 
        stages = sorted(df_input.Stage.unique())
        bool_stages = True
    except: 
        bool_stages = False
    
    out_stages = dict()
    if bool_stages:
        leg_stages = {'32-34': {'legend': '32-34hpf', 'color': '#72e5ef'}, 
                      '48-50': {'legend': '48-50hpf', 'color': '#af437c'},
                      '58-60': {'legend': '58-60hpf', 'color': '#94d86f'},
                      '72-74': {'legend': '72-74hpf', 'color': '#5f4ac2'},
                      '32->48': {'legend': '32$\Rightarrow$48hpf', 'color': '#FF6347'},
                      '48->74': {'legend': '48$\Rightarrow$74hpf', 'color': '#FFFF00'}}
        
        for stg in stages:
            out_stages[stg] = leg_stages[stg]
            
    # SPECIFIC ANALYSIS
    try: 
        spAn = sorted(df_input.spAnalysis.unique())
        bool_spAn = True
    except: 
        bool_spAn = False
    
    out_spAn = dict()
    if bool_spAn:
        leg_spAn = {'normal': {'legend':'WT', 'color': '#72e5ef'}, 
                      'No looper': {'legend': 'No looper', 'color': '#af437c'},
                      'Right looper': {'legend': 'Right Looper', 'color': '#94d86f'},
                      'Left looper': {'legend': 'Left Looper', 'color': '#5f4ac2'}}

        for spA in spAn:
            out_spAn[spA] = leg_spAn[spA]
    
    # TIME-POINTS
    try: 
        time_points = sorted(df_input.time_point.unique())
        bool_tp = True
    except: 
        bool_tp = False
    
    out_tp = []
    if bool_tp:
        for tp in time_points:
            out_tp.append(tp)#leg_stages[all_stages.index(stg)])
            
    # FISH Reference
    try: 
        fish_ref = sorted(df_input.Fish_ref.unique())
        bool_fr = True
    except: 
        bool_fr = False
    
    out_fr = []
    if bool_fr:
        for fr in fish_ref:
            out_fr.append(fr)
    
    x_labels = {'GenotypeAll': 'Genotype', 'GenotypeF': 'Genotype_f',
                    'Strain' : 'Strain', 'Strain_o' : 'Strain_o',
                    'Stage': 'Stage [hpf]', 'time_point': 'Heart phase contraction',
                    'Fish_ref' : 'Fish Reference', 'spAnalysis': 'Specific Analysis'}
                
    dict_legends = {'GenotypeAll': out_genots, 'GenotypeF': out_genotsF,
                    'Strain' : out_strains, 'Strain_o' : out_strains_o,
                    'Stage': out_stages, 'time_point' : out_tp, 'Fish_ref': out_fr,
                    'spAnalysis': out_spAn, 'xlabels': x_labels}

    return dict_legends 

#%% func - genot_order
def genot_order(values, mix_wts = False):
    #['hapln1a:wt', 'hapln1a:ht', 'hapln1a:mt', 
        #               'hapln1a:wt/spaw:wt', 'hapln1a:wt/spaw:mt', 'hapln1a:mt/spaw:wt', 'hapln1a:mt/spaw:mt',
        #               'galff:+/uas:+', 'galff:+/uas:-', 'galff:-/uas:+', 'galff:-/uas:-',
        #               'vcana:wt', 'vcana:ht', 'vcana:mt', 'wt:wt']
        
    # print('-> Genot Values: ', values )
    if not mix_wts:
        poss_combinations = [['hapln1a241:wt','hapln1a241:mt'],
                             ['hapln1a241:wt','hapln1a241:ht','hapln1a241:mt'],
                             ['hapln1a187:wt','hapln1a187:mt'],
                             ['hapln1a187:wt','hapln1a187:ht','hapln1a187:mt'],
                             ['hapln1a241:wt','hapln1a241:mt','hapln1a187:wt','hapln1a187:mt'],
                             ['hapln1a241:wt','hapln1a241:mt','hapln1a187:mt'],
                             
                             ['hapln1a241:wt/spaw:wt','hapln1a241:wt/spaw:mt'],
                             ['hapln1a241:wt/spaw:wt','hapln1a241:wt/spaw:mt','hapln1a241:mt/spaw:wt','hapln1a241:mt/spaw:mt'],
                             ['hapln1a241:wt/spaw:wt','hapln1a241:mt','hapln1a241:wt/spaw:mt','hapln1a241:mt/spaw:mt'],
                             ['hapln1a241:wt','hapln1a241:mt','hapln1a241:wt/spaw:mt','hapln1a241:mt/spaw:mt'],
                             ['hapln1a241:wt','hapln1a241:mt','hapln1a241:wt/spaw:mt'],
                             
                             ['galff:-/uas:-','galff:+/uas:-', 'galff:-/uas:+','galff:+/uas:+'],
                             ['galff:-/uas:-','galff:+/uas:-', 'galff:-/uas:+','galff:+/uas:+','hapln1a241:mt'],
                             ['galff:+/uas:-', 'galff:-/uas:+','galff:+/uas:+'],
                             ['galff:+/uas:-', 'galff:-/uas:+','galff:+/uas:+','hapln1a241:mt'],
                             ['vcana365:wt', 'vcana365:ht', 'vcana365:mt'],
                             ['vcana365:wt', 'vcana365:mt'],
                             ['wt:wt','wt_tc:wt','hapln1a241:wt','hapln1a241:ht','hapln1a241:mt', 
                              'hapln1a187:wt','hapln1a187:ht','hapln1a187:mt',
                              'hapln1a3bp:wt','hapln1a3bp:ht','hapln1a3bp:mt',
                              'hapln1aSTOP:wt','hapln1aSTOP:ht','hapln1aSTOP:mt',
                              'hapln1a241:wt/spaw:wt', 'hapln1a241:wt/spaw:ht', 'hapln1a241:mt/spaw:wt','hapln1a241:wt/spaw:mt', 'hapln1a241:mt/spaw:mt',
                              'galff:-/uas:-','galff:+/uas:-', 'galff:-/uas:+', 'wt:ctrl','galff:+/uas:+',
                              'vcana365:wt', 'vcana365:ht','vcana365:mt']]
        
        for poss in poss_combinations:
            if set(values) == set(poss):
                res_order = poss
                # print('-> Resulting Order for Values: ', res_order)
            else: 
                res_order = []
                for val in poss_combinations[-1]:
                    if val in values:
                        res_order.append(val)
                # print('ERROR: No match with possible comination! Update combination list!')
            
    return res_order
    
#%% func - props_ordered
def props_ordered(df2plot, x_var, hue_var, shape_var):
    values = []
    for var in [x_var, hue_var, shape_var]:
        # print(var)
        if var =='GenotypeAll':
            vals2app = genot_order(df2plot[var].unique(), mix_wts = False)
            # print(vals2app)
            values.append(vals2app)
        elif var =='GenotypeF':
            vals2app = genot_order(df2plot[var].unique(), mix_wts = False)
            # print(vals2app)
            values.append(vals2app)
        else: 
            reverse = False
            values.append(sorted(df2plot[var].unique(), reverse=reverse))
    x_values, hue_values, shape_values = values
    
    return x_values, hue_values, shape_values

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
            if bool_genotAll:
                df_change ['GenotypeAll'] = [gen for nn in range(len(stages)-1)]
            if bool_genotF:
                df_change ['GenotypeF'] = [gen for nn in range(len(stages)-1)]
            if bool_strains_o:
                df_change['Strain_o'] = [strain for nn in range(len(stages)-1)]
            df_change = df_change.reset_index().set_index(new_index)
            # print(df_change)
            df_pct_change.append(df_change)
            
    bool_strains_o = False
    bool_genotAll = False
    bool_genotF = False
    
    df4plot = df2plot.groupby(group_vars)[vars2plot].mean()
    new_index = []
    if 'Strain_o' in group_vars: 
        strains_o = df2plot.Strain_o.unique()
        bool_strains_o = True
        new_index.append('Strain_o')
    if 'GenotypeAll' in group_vars:
        genots = df2plot.GenotypeAll.unique()
        bool_genotAll = True
        new_index.append('GenotypeAll')
    if 'GenotypeF' in group_vars:
        genots = df2plot.GenotypeF.unique()
        bool_genotF = True
        new_index.append('GenotypeF')
    new_index.append('Stage')
    
    stages = sorted(df2plot.Stage.unique())
    
    df_pct_change = []
    if bool_genotAll or bool_genotF:
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



#%% - STATISTICAL ANALYSIS
#%% func - normalityTests
def normalityTests (alpha, factor, fact_levels, variable, data, test):
    
    pval_norm = np.ones(((len(fact_levels)),3))
    
    all_data_normal = True
    txt_normtest = test +' Normality tests \n('+'alpha:'+str(alpha)+'), \n- p-val: - '
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
        txt_normtest = txt_normtest + normtxt 
        if level != fact_levels[-1]:
            txt_normtest = txt_normtest +', \n '
        else: 
            txt_normtest = txt_normtest + ' - '  
        
        data_groups.append(data2testNorm)
    print(" >> Data normally distributed?", all_data_normal)
    # print(txt_normtest)
    
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
            pvalues_norm, txt_normtest, all_data_normal, data_groups = normalityTests(alpha, factor, fact_levels, var, data_out, norm_test)
            # print(var)
            test2use, multcomp_txt, txt_testSelected = selectStatisticalTest(fact_levels, all_data_normal, var)
            
            dict_spStatisticalRes_in['pvalues_norm'] = pvalues_norm
            dict_spStatisticalRes_in['txt_normtest'] = txt_normtest
            dict_spStatisticalRes_in['all_data_normal'] = all_data_normal
            dict_spStatisticalRes_in['test2use'] = test2use
            dict_spStatisticalRes_in['multcomp_txt'] = multcomp_txt
            
            if test2use == 'One-way-ANOVA':
                pval_spTest, pval_mcTest, pval_multComp = stTest_OWANOVA(factor, var, indiv_pairs, data_out, data_groups)
                dict_spStatisticalRes_in['pval_test'] = pval_spTest
                dict_spStatisticalRes_in['pval_multComp_o'] = pval_mcTest
                dict_spStatisticalRes_in['pval_multComp'] = pval_multComp
                dict_spStatisticalRes_in['stat_test'] = False
                dict_spStatisticalRes_in['test'] = []
                txt_testSelected = txt_testSelected+" (ANOVA's p-val: "+str(format(pval_spTest,'.3f'))+'), \n'+multcomp_txt
                
            elif test2use == 'Kruskal':
                pval_spTest, pval_mcTest, pval_multComp = stTest_Kruskal(factor, var, indiv_pairs, data_out, data_groups)
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
    elif len(data_groups) == 5:
        statOWA, pval_OWA = stats.f_oneway(data_groups[0], data_groups[1], data_groups[2], data_groups[3], 
                                           data_groups[4])
    elif len(data_groups) == 6:
        statOWA, pval_OWA = stats.f_oneway(data_groups[0], data_groups[1], data_groups[2], data_groups[3], 
                                           data_groups[4], data_groups[5])
    elif len(data_groups) == 7:
        statOWA, pval_OWA = stats.f_oneway(data_groups[0], data_groups[1], data_groups[2], data_groups[3], 
                                           data_groups[4], data_groups[5], data_groups[6])
    
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
    elif len(data_groups) == 5:
        statKruskal, pval_Kruskal = stats.kruskal(data_groups[0], data_groups[1], data_groups[2], data_groups[3],
                                                  data_groups[4])
    elif len(data_groups) == 6:
        statKruskal, pval_Kruskal = stats.kruskal(data_groups[0], data_groups[1], data_groups[2], data_groups[3], 
                                                  data_groups[4], data_groups[5])
    elif len(data_groups) == 7:
        statKruskal, pval_Kruskal = stats.kruskal(data_groups[0], data_groups[1], data_groups[2], data_groups[3],
                                                  data_groups[4], data_groups[5], data_groups[6])
        
        
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
def txtMultComp (box_pairs, pval, x_values, x_labels, test_selected, test_norm, alpha, indiv = False):
    txt_multcomp = '\n -- STATISTICS -- \n'
    
    if indiv:
        chunks = [box_pairs]
        x_values = [x_values]
        x_labels = [x_labels]
    else:
        num = int(len(box_pairs)/len(x_values))
        chunks = [box_pairs[x:x+num] for x in range(0, len(box_pairs), num)]
    
    n = 0
    for j, chunk, x_val, x_lab, test, norm in zip(count(), chunks, x_values, x_labels, test_selected, test_norm):
        txt_multcomp = txt_multcomp + x_lab + ': '+ norm +'\n'+ test + '\n - p-val: '
        for i, pair in enumerate(chunk):
            # print(pair, n, i)
            a_pair = list(pair[0]); b_pair = list(pair[1])
            # print(x_val,'-',a_pair, '-',b_pair)
            a_pair.remove(x_val);b_pair.remove(x_val)
            pval_txt = str(format(pval[n],'.5f'))
            if pval[n] > alpha:
                pval_txt = pval_txt +' (ns)'
            elif pval[n] <= 1e-4:
                pval_txt = pval_txt +' (****)'
            elif pval[n] <= 1e-3:
                pval_txt = pval_txt +' (***)'
            elif pval[n] <= 1e-2:
                pval_txt = pval_txt +' (**)'
            elif pval[n] <= 0.05:
                pval_txt = pval_txt +' (*)'
            pairtxt =  str(a_pair)+ ' vs. '+ str(b_pair) +': '+pval_txt
            n += 1
            txt_multcomp = txt_multcomp + pairtxt
            if i < len(chunk)-1:
                txt_multcomp = txt_multcomp + ' \n '
            else: 
                txt_multcomp = txt_multcomp + ' - '

        txt_multcomp = txt_multcomp+'\n - - - - - - - - - - \n'
     
    return txt_multcomp

#%% - CREATE PLOTS
class Thesis_contxt():
    def __init__(self,context):
        self.context = context
        
    def get_context(self):
        return self.context
    
    def set_context(self, new_context):
        self.context = new_context

#%% Plot properties!
#%% func - setSNSContext
def setSNSContext(thesis_contxt):
    fontname = 'Myriad Pro'
    
    if thesis_contxt: # for thesis
        contxt = 'poster'
        font_scale = 1
        rc_dict = {'font.size': 8, 
                   'fontname' : 'Myriad Pro',
                   # 'lines.linewidth': 3.0,
                   'axes.linewidth': 0.7, 
                   'axes.labelsize': 8.0, #size for the axis labels (names)
                   'axes.titlesize': 8, #'axes.titleweight': 200,
                   'xtick.bottom' : True, 'ytick.left' : True,
                   'xtick.major.size': 3, 'ytick.major.size': 3,
                   'xtick.major.width': 0.25, 'ytick.major.width': 0.25,
                   # 'xtick.minor.width': 0.25, 'ytick.minor.width': 0.25,
                   'xtick.labelsize': 6, 'ytick.labelsize': 8, # size of the labels of each tick
                   'legend.fontsize': 8, 'legend.title_fontsize': 8,
                   'figure.dpi': 100, 'savefig.dpi': 300, 'svg.fonttype': None}
                   # 'figure.figsize': (3,2)}
        
        boxprops = dict(linestyle = '-', linewidth = 0.25, alpha = 0.4)
        flierprops = dict(marker='*', markersize=2.5, markerfacecolor='#FFD700', markeredgecolor = '#000000',markeredgewidth = 0.15)#linestyle = 'none')
        whiskerprops = dict(linewidth = 0.35, alpha = 0.6, color='#000000')#'708090')
        capprops = dict(linewidth = 0.35, color='#708090')
        medianprops = dict(linewidth=0.2, linestyle = '--', color = '#708090') #'white')#
        meanprops = dict(linewidth=0.1, linestyle = '-.', color = '#696969')
        
    else: 
        contxt = 'poster'
        font_scale = 1
        rc_dict = {'font.size': 20, 
                    'fontname' : 'Myriad Pro',
                    # 'lines.linewidth': 3.0,
                    'axes.linewidth': 1.5, 
                    'axes.labelsize': 16.0, #size for the axis labels (names)
                    'axes.titlesize': 10, #'axes.titleweight': 200,
                    'xtick.bottom' : True, 'ytick.left' : True,
                    'xtick.major.size': 10, 'ytick.major.size': 8,
                    'xtick.major.width': 1, 'ytick.major.width': 1,
                    'xtick.labelsize': 14, 'ytick.labelsize': 14, # size of the labels of each tick
                    'legend.fontsize': 16, 'legend.title_fontsize': 16,
                    'figure.dpi': 100, 'savefig.dpi': 300, 'svg.fonttype': None}
        
        boxprops = dict(linestyle = '-', linewidth = 1.2, alpha = 0.4)
        flierprops = dict(marker='*', markersize=8, markerfacecolor='#FFD700', markeredgecolor = '#000000',markeredgewidth = 0.8)#linestyle = 'none')
        whiskerprops = dict(linewidth = 1, color='#A9A9A9')
        capprops = dict(color='#A9A9A9')
        medianprops = dict(linewidth=1, linestyle = '-', color = 'white')#'#696969')
        meanprops = dict(linewidth=1.5, linestyle = '-.', color = '#696969')
        
    all_box_props = [boxprops, flierprops, whiskerprops, capprops, medianprops, meanprops]
    
    return contxt, font_scale, rc_dict, all_box_props, fontname

# from matplotlib import rcParams
# rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Tahoma', 'Verdana', 'DejaVu Sans',
#                                'Lucida Grande']
# To print context parameters
# print(sns.plotting_context())
# https://matplotlib.org/devdocs/tutorials/introductory/customizing.html

# Set context for all figures
contxt = 'poster'
font_scale = 1
rc_dict = {'font.size': 20, 
            'fontname' : 'Myriad Pro',
            # 'lines.linewidth': 3.0,
            'axes.linewidth': 1.5, 
            'axes.labelsize': 16.0, #size for the axis labels (names)
            'axes.titlesize': 10, #'axes.titleweight': 200,
            'xtick.bottom' : True, 'ytick.left' : True,
            'xtick.major.size': 10, 'ytick.major.size': 8,
            'xtick.major.width': 1, 'ytick.major.width': 1,
            'xtick.labelsize': 14, 'ytick.labelsize': 14, # size of the labels of each tick
            'legend.fontsize': 16, 'legend.title_fontsize': 16,
            'figure.dpi': 100, 'savefig.dpi': 300, 'svg.fonttype': None}

boxprops = dict(linestyle = '-', linewidth = 1.2, alpha = 0.3)
flierprops = dict(marker='*', markersize=8, markerfacecolor='#FFD700', markeredgecolor = '#000000',markeredgewidth = 0.8)#linestyle = 'none')
whiskerprops = dict(linewidth = 1, color='#A9A9A9')
capprops = dict(color='#A9A9A9')
medianprops = dict(linewidth=1, linestyle = '-', color = 'white')#'#696969')
meanprops = dict(linewidth=1.5, linestyle = '-.', color = '#696969')


# contxt = 'poster'
# font_scale = 1
# rc_dict = {'font.size': 8, 
#            'fontname' : 'Myriad Pro',
#            # 'lines.linewidth': 3.0,
#            'axes.linewidth': 0.7, 
#            'axes.labelsize': 8.0, #size for the axis labels (names)
#            'axes.titlesize': 8, #'axes.titleweight': 200,
#            'xtick.bottom' : True, 'ytick.left' : True,
#            'xtick.major.size': 3, 'ytick.major.size': 3,
#            'xtick.major.width': 0.25, 'ytick.major.width': 0.25,
#            'xtick.labelsize': 7, 'ytick.labelsize': 8, # size of the labels of each tick
#            'legend.fontsize': 8, 'legend.title_fontsize': 8,
#            'figure.dpi': 100, 'savefig.dpi': 300, 'svg.fonttype': None}
#            # 'figure.figsize': (3,2)}

# boxprops = dict(linestyle = '-', linewidth = 0.15, alpha = 0.4)
# flierprops = dict(marker='*', markersize=2.5, markerfacecolor='#FFD700', markeredgecolor = '#000000',markeredgewidth = 0.15)#linestyle = 'none')
# whiskerprops = dict(linewidth = 0.25, color='#708090')
# capprops = dict(linewidth = 0.25, color='#708090')
# medianprops = dict(linewidth=0.2, linestyle = '--', color = '#708090') #'white')#
# meanprops = dict(linewidth=0.1, linestyle = '-.', color = '#696969')

# Setting font!
# https://jonathansoma.com/lede/data-studio/matplotlib/changing-fonts-in-matplotlib/
# https://jonathansoma.com/lede/data-studio/matplotlib/list-all-fonts-available-in-matplotlib-plus-samples/
# https://bastibe.de/2016-05-30-matplotlib-font-cache.html
# import matplotlib.font_manager
# matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
# matplotlib.font_manager._rebuild()

fontname = 'Myriad Pro' #'Comic Sans MS'

# Palette stuff
# Save a palette to a variable:
# palette = sns.color_palette("husl", 8*n_gen*n_strain)
# sns.palplot(palette)
# for gen in range(n_gen):
#     palettes.append(palette[gen*len(palette)//n_gen+2])
    
eight_palette = ["#dc133b", "#ffa500", "#6495ed", "#ade64f", "#8b6fed", "#859947", "#994440", "#36cbd3"]
# twelve_palete = ["#dc133b", "#ffa500", "#6495ed", "#ade64f", "#c36bea", "#214d4e", "#36cbd3", "#a23e27", "#839292", "#8b008b", "#3ff44c", "#000080"]

# twelve pallete with metadata
# [{"color":"#dc133b", "name":"dark red","textColor":"white"},
# {"color":"#ffa500","name":"orange","textColor":"black"},
# {"color":"#6495ed","name":"blue","textColor":"black"},
# {"color":"#ade64f","name":"yellow green","textColor":"black"},
# {"color":"#c36bea","name":"pink","textColor":"black"},
# {"color":"#214d4e","name":"teal", "textColor":"white"},
# {"color":"#36cbd3", "name":"turquoise", "textColor":"black"},
# {"color":"#a23e27", "name":"brown", "textColor":"white"},
# {"color":"#839292", "name":"grey", "textColor":"black"},
# {"color":"#8b008b", "name":"magenta", "textColor":"white"},
# {"color":"#3ff44c", "name":"green","textColor":"black"},
# {"color":"#000080","name":"indigo","textColor":"white"}]

styles = ['o', '^', 's', 'v', 'D', '<', 'p', '>']*20 # https://matplotlib.org/stable/api/markers_api.html
palettes = ['deeppink','royalblue','mediumturquoise', 'darkmagenta','darkorange','limegreen', 'gold']*20

# Other links
# https://towardsdatascience.com/matplotlib-vs-ggplot2-c86dd35a9378
# https://www.r-graph-gallery.com/boxplot.html
# https://medium.com/analytics-vidhya/r-style-visualizations-in-python-560c6bbfb14a
# https://matplotlib.org/2.0.2/examples/pylab_examples/spine_placement_demo.html
# https://datascienceplus.com/seaborn-categorical-plots-in-python/
# https://towardsdatascience.com/scattered-boxplots-graphing-experimental-results-with-matplotlib-seaborn-and-pandas-81f9fa8a1801
# https://github.com/cfcooney/medium_posts/blob/master/scattered_boxplots.ipynb
# https://htmlcolorcodes.com/

# Padding 
# https://dfrieds.com/data-visualizations/style-plots-python-matplotlib.html
# https://matplotlib.org/stable/tutorials/introductory/customizing.html

#%% ->No statistics
import matplotlib.ticker as ticker
#%% func - plotIndivperX
def plotIndivperX (# General Plot Settings
                    df2plot, vars2plot, x_var, hue_var, shape_var, 
                    # Size, title, labels and legends
                    title, labels2plot, dict_legends, ips, 
                    # Other plot settings
                    suptitle = True, right_legend = False, 
                    ctx_style = 0, uni_color = False, 
                    yticks_lab = 'th,', ylim = [], yset = '', box_plot = True, show_fliers = False,
                    # Saving settings
                    save = True, dpi = 300, ext = ['png'], info = '', dir2save = ''):

    contxt, font_scale, rc_dict, all_box_props, fontname  = setSNSContext(ctx_style)
    boxprops, flierprops, whiskerprops, capprops, medianprops, meanprops = all_box_props
    
    print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    
    # Genotypes and Strains being plotted 
    values = props_ordered(df2plot, x_var, hue_var, shape_var)
    x_values, hue_values, shape_values = values
    # print(x_values, hue_values, shape_values)
    
    # Set up the matplotlib figure 
    n_rows = 1
    n_cols = len(hue_values)
    addVars = False
    for index in range(n_cols):
        if index == 0 and len(vars2plot) == 1:
            addVars = True
            var2copy = vars2plot[0]
            label2copy = labels2plot[0]
        
        if addVars and index != 0:
            vars2plot.insert(index, var2copy)
            labels2plot.insert(index, label2copy)
    
    # As a right legend wants to be added, we need to add '' to all labels 
    if right_legend: 
        index_no_plot = [n_cols]
        for index in index_no_plot:
            vars2plot.insert(index, '')
            labels2plot.insert(index, '')
        width_ratios = [1]*n_cols+[0.2]
        n_cols = n_cols+1
    else: 
        width_ratios = [1]*n_cols
        index_no_plot = [1000]
        
    # Set up the matplotlib figure
    h_plot, w_plot = ips
    
    
    if len(x_values) >= 6 and len(hue_values) >= 2:
        h_add = 1.4
        wspace = 0.05
        box_width = 0.4
        bbox_lf = -0.5
        w_add = 0.2
        print('A', len(x_values), len(hue_values))
    
    elif len(x_values) >= 3 and len(hue_values) >= 2:
        h_add = 1
        wspace = 0.05
        box_width = 0.4
        bbox_lf = -0.5
        w_add = 0.2
        print('A', len(x_values), len(hue_values))
    
    elif len(x_values) >= 3 and len(hue_values) < 2:
        h_add = 0.4
        wspace = 0.05
        box_width = 0.4
        bbox_lf = -0.5
        w_add = 0.2
        print('C', len(x_values), len(hue_values))
        
    else:
        h_add = 0.6
        h_plot = 0.2#0.8
        wspace = 0.15#0.15
        box_width = 0.4
        bbox_lf = -1.8
        w_add = 0.2#1
        print('B', len(x_values), len(hue_values))

    size_col = (n_cols)*(h_plot*len(x_values))+h_add
    size_row = n_rows*w_plot+w_add
    print('-size: ', size_col, size_row)
    
    # Define legends for x and hue
    v_legend = []; v_color = []
    for zz, v_var, v_values in zip(count(), [x_var, hue_var], [x_values, hue_values]):
        if isinstance(dict_legends[v_var], dict):
            legend = []
            color = []
            for v_dict in v_values: 
                legend.append(dict_legends[v_var][v_dict]['legend'])
                color.append(dict_legends[v_var][v_dict]['color'])
        else: 
            legend = dict_legends[v_var]
            color = ''
        v_legend.append(legend)
        v_color.append(color)
    x_legend, hue_legend = v_legend
    x_color, _ = v_color
    print(x_color, x_values)
    
    if uni_color:
        for pp, hue_val in enumerate(dict_legends[hue_var]):
            print(hue_val)
            color = dict_legends[hue_var][hue_val]['color']
        print(pp+1, color)
        x_color = [color for n in range(len(x_values))]
        print(x_color)
        
    if right_legend:
        # Define legends for x
        legend_elem = []
        for aa, xval, xcol in zip(count(), x_legend, x_color):
            legend_elem.append(Line2D([0], [0], marker='o', color='w', label=xval,
                                    markerfacecolor=xcol, markersize=12))
        handle_new = legend_elem
        
        for index in index_no_plot:
            hue_legend.insert(index, '')
            hue_values.insert(index, '')
            
    marker_size = 2.5; dodge = True; jitter = 0.3
    
    for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
        print('\t- '+svar+': ', value)
        
    ##  CREATE FIGURE 
    gridkw = dict(width_ratios=width_ratios)
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
    fig.subplots_adjust(hspace=1, wspace=wspace)
    
    sns.set_style("ticks")
    sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    
    # print(n_cols, n_rows, axes, right_legend, vars2plot, labels2plot, hue_values)
    if n_cols == 1 and n_rows == 1 and not right_legend:
        axes_fl = [axes]
    else: 
        axes_fl = axes.flatten()
    
    # Create list to contain info of yticks and its highest value
    # y_vals_all = []; 
    max_y_vals = []; min_y_vals = []
    for n, ax, var, ylabel, hue_value in zip(count(), axes_fl, vars2plot, labels2plot, hue_values):
        if n in index_no_plot and right_legend:
                if n == n_cols-1:
                    ax.set_axis_off()
                    ax.legend(handle_new, x_legend, loc='upper left', bbox_to_anchor=(bbox_lf, 1), frameon = False)
                else: 
                    ax.remove()
        else: 
            df_xfilt = df2plot[df2plot[hue_var] == hue_value]
            m = sns.stripplot(data=df_xfilt, x=x_var, y=var, ax = ax, order=x_values,
                              marker = 'o', palette = x_color, size = marker_size, 
                              linewidth=0.3, jitter = jitter)
            if box_plot: 
                m = sns.boxplot(data=df_xfilt, x=x_var, y=var, ax = ax, order=x_values,
                                   dodge = dodge, width= box_width, showfliers = show_fliers, palette = x_color,
                                   boxprops = boxprops, whiskerprops = whiskerprops, capprops = capprops, 
                                   flierprops = flierprops, medianprops = medianprops)#, 
                                   # meanline = True, showmeans = False, meanprops = meanprops,)
    
            box = ax.get_position()
            # ax.tick_params(which='major', width=0.25, length=1)
            # ax.tick_params(which='minor', width=0.25, length=0.5)
            # ax.yaxis.set_major_locator(plt.MaxNLocator(6))
            # ax.yaxis.set_minor_locator(plt.MaxNLocator(0))
            # ax.locator_params(axis='y', nbins=4)
            
            if n == 0:
                if yticks_lab == '1e6 - d.':
                    ylabel = ylabel +' x 10$^6$'
                elif yticks_lab == '1e3 - d.':
                    ylabel = ylabel +' x 10$^3$'
                ax.set_xlabel(hue_legend[n], fontname = fontname)
                ax.set_ylabel(ylabel, fontname = fontname)
            else:
                ax.set_xlabel(hue_legend[n], fontname = fontname)
                ax.set_ylabel('', fontname = fontname)
    
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            xticks = ax.get_xticks()
            ax.set_xticks(xticks)
            ax.set_xticklabels(x_legend, rotation=45, horizontalalignment='right', fontname = fontname)
            sns.despine()    
            
            if isinstance(ylim, tuple):
                print('ylim:', ylim[0],'-', ylim[1])
                ax.set_ylim(ylim[0], ylim[1])
            
            if n == 0:
                handles, labels = m.get_legend_handles_labels()
            else:
                ax.spines['left'].set_linestyle('-.')
                ax.spines['left'].set_color('#696969')
                ax.spines['left'].set_linewidth(0.6)
                ax.tick_params(left = False)
    
            y_vals = ax.get_yticks()
            max_y_vals.append(y_vals[-1])
            min_y_vals.append(y_vals[0])
    
    # Find min and maximum value for each axis 
    # print(min_y_vals, max_y_vals)
    min_y = min(min_y_vals); max_y = max(max_y_vals)

    if yset == 'round':
        rounding_set = -1
        y_step = round((max_y - min_y)/8,rounding_set)
        if y_step == 0: 
            rounding_set = 0
            alert('frog',1)
        y_step = round((max_y - min_y)/8,rounding_set)
        y_vals_set = list(range(int(min_y), int(max_y), int(y_step)))
        y_vals_set.append(y_vals_set[-1]+int(y_step))
        
    elif yset == 'dec':
        rounding_set = 2
        y_step = round((max_y - min_y)/8,rounding_set)
        # print(y_step)
        num_y = (max_y - min_y)//y_step
        # print(num_y)
        y_vals_set = np.linspace(min_y, max_y,int(num_y), endpoint = True)
    # print(y_vals_set)
    
    y_valsF = y_vals_set#y_vals_all[max_y_index]
    
    for n, ax in zip(count(), axes_fl):
        ax.set_yticks(y_valsF)
        # Define axes based on higherst bar
        if yticks_lab == '1e6 - d.':
            ax.set_yticks(ax.get_yticks())  #
            ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_valsF])
        elif yticks_lab == '1e3 - d.':
            ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_valsF])
        elif yticks_lab == 'th,':
            ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_valsF])
        elif yticks_lab == 'd. - 0':
            ax.set_yticklabels(['{:.0f}'.format(w) for w in y_valsF])
        elif yticks_lab == 'd. - 1':
            ax.set_yticklabels(['{:.1f}'.format(w) for w in y_valsF])
        elif yticks_lab == 'd.':
            ax.set_yticklabels(['{:.2f}'.format(w) for w in y_valsF])
        
        for tick in ax.get_xticklabels():
            tick.set_fontname(fontname)
        for tick in ax.get_yticklabels():
            tick.set_fontname(fontname)
        try:
            ax.yaxis.get_offset_text().set_fontname(fontname)
            # print(r' -> YAY')
        except:
            print(r' -> Nop')

    if suptitle:
        fig.suptitle(title+'\n', y=1, fontname = fontname)
    # fig.tight_layout()
    
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'R_')
            # fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
            fig_title = dir2savef+info+"_"+title+"."+extf

            plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)

#%% func - plotInGroups
# def plotInGroups(# General Plot Settings
#                  df2plot, vars2plot, x_var, hue_var, shape_var, 
#                  # Size, title, labels and legends
#                  title, labels2plot, dict_legends, ips, 
#                  # Other plot settings
#                  n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th,', ylim = '', 
#                  # Saving settings
#                  save = True, dpi = 300, ext = 'png', info ='',dir2save = ''):
    
#     print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    
#     # Set up the matplotlib figure
#     num_vars = len(vars2plot)
#     n_rows = math.ceil(num_vars/n_cols)
#     # print('n_rows:' , n_rows)
    
#     index_right_col = list(range(n_cols,(n_cols+1)*n_rows,4))
#     index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
#     index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))

#     for index in index_no_plot:
#         vars2plot.insert(index, '')
#         labels2plot.insert(index, '')
    
#     # Set up the matplotlib figure
#     h_plot, w_plot = ips
    
#     # Genotypes and Strains being plotted 
#     values = props_ordered(df2plot, x_var, hue_var, shape_var)
#     x_values, hue_values, shape_values = values
    
#     # - number of x_var
#     n_x = len(x_values)
#     if n_x == 1:
#         h_plot = 3.5
        
#     if num_vars == 1:
#         h_add = 0; w_add = 0
#     size_col = (n_cols+1)*h_plot+h_add
#     size_row = n_rows*w_plot+w_add
    
#     # Define legends for x and hue
#     v_legend = []; v_color = []
#     for zz, v_var, v_values in zip(count(), [x_var, hue_var], [x_values, hue_values]):
#         if isinstance(dict_legends[v_var], dict):
#             legend = []
#             color = []
#             for v_dict in v_values: 
#                 legend.append(dict_legends[v_var][v_dict]['legend'])
#                 color.append(dict_legends[v_var][v_dict]['color'])
#         else: 
#             legend = dict_legends[v_var]
#             color = ''
#         v_legend.append(legend)
#         v_color.append(color)
#     x_legend, hue_legend = v_legend
#     x_color, hue_color = v_color

#     legend_elem_hue = []
#     for aa, hue_val, hue_col in zip(count(), hue_legend, hue_color):
#         legend_elem_hue.append(Line2D([0], [0], marker='o', color='w', label=hue_val,
#                                 markerfacecolor=hue_col, markersize=20))
        
#     handle_new = legend_elem_hue
#     legend_new = hue_legend
    
#     marker_size = 10; dodge = True; jitter = 0.3
#     plot_no = 0
    
#     #  CREATE FIGURE
#     gridkw = dict(width_ratios=[1]*n_cols+[0.2])
#     fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
#     fig.subplots_adjust(hspace=0.5, wspace=0.2)
    
#     sns.set_style("ticks")
#     sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    
#     for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
        
#         # Create list to contain info of yticks and ist highest value
#         y_vals_all = []
#         max_y_vals = []
        
#         if n == 0: 
#             for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
#                 print('\t- '+svar+': ', value)
                
#         if n in index_no_plot:
#             if n == n_cols:
#                 ax.set_axis_off()
#                 ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-0.5, 1), frameon = False)
#             else: 
#                 ax.remove()
#         else: 
#             print(' > ',var,'\n')
#             m = sns.boxplot(data=df2plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
#                             dodge = dodge, width= 0.8, showfliers = True, palette = hue_color,
#                             boxprops = boxprops, whiskerprops = whiskerprops, capprops = capprops, 
#                             flierprops = flierprops, medianprops = medianprops, 
#                             meanline = True, meanprops = meanprops, showmeans = True)

#             m = sns.stripplot(data=df2plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
#                               marker = 'o', palette = hue_color, dodge = dodge, size = marker_size, 
#                               linewidth=1, jitter = jitter)

#             box = ax.get_position()
#             if yticks_lab == '1e6 - d.':
#                 ylabel = ylabel +' x 10$^6$'
#             elif yticks_lab == '1e3 - d.':
#                 ylabel = ylabel +' x 10$^3$'
                
#             ax.set_xlabel(hue_legend, fontname = fontname)
#             ax.set_ylabel(ylabel, fontname = fontname)
#             ax.set_xticklabels(x_legend, rotation=45, horizontalalignment='right', fontname = fontname)
#             ax.set_position([box.x0, box.y0, box.width*1, box.height])
#             ax.get_legend().remove()
#             sns.despine()
            
#             if ylim != '':
#                 # print('ylim:', ylim[plot_no][0],'-', ylim[plot_no][1])
#                 ax.set_ylim(ylim[plot_no][0], ylim[plot_no][1])
            
#             if n == 0:
#                 handles, labels = m.get_legend_handles_labels()

#             y_vals = ax.get_yticks()
#             y_vals_all.append(y_vals)
#             max_y_vals.append(y_vals[-1])
#             if n_x > 1: 
#                 vline_pos = list(range(1,n_x,1))
#                 for pos in vline_pos:
#                     ax.axvline(pos - 0.5, ymax = 0.95, color='dimgrey', ls='-.', linewidth=0.8)
            
#             # Define axes based on higherst bar
#             max_y_index = max_y_vals.index(max(max_y_vals))
#             ax.set_yticks(y_vals_all[max_y_index])
#             if yticks_lab == '1e6 - d.':
#                 ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == '1e3 - d.':
#                 ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'th,':
#                 ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 0':
#                 ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 1':
#                 ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd.':
#                 ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])
                
#             for tick in ax.get_xticklabels():
#                 tick.set_fontname(fontname)
#             for tick in ax.get_yticklabels():
#                 tick.set_fontname(fontname)
#             try:
#                 ax.yaxis.get_offset_text().set_fontname(fontname)
#             except:
#                 print(r' -> Nop')
                
#         plot_no +=1

#     fig.suptitle(title+'\n', y=1, fontname = fontname)
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'R_')
#             if info != '':
#                 fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
#             else: 
#                 fig_title = dir2savef+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf

#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)

#%% ->With statistics        
#%% func - plotIndivStatsperX
# def plotIndivStatsperX (# General Plot Settings
#                     df2plot, vars2plot, x_var, hue_var, shape_var, 
#                     # Size, title, labels and legends
#                     title, labels2plot, dict_legends, ips, 
#                     # Statistic Settings
#                     stats_set, statsTxt = True, 
#                     # Other plot settings
#                     suptitle = True, right_legend = True,
#                     yticks_lab = 'th,', ylim = '', yset = '', box_plot = True, show_fliers = False,
#                     # Saving settings
#                     save = True, dpi = 300, ext = 'png', info = '', dir2save = ''):

#     if stats_set[0]: 
#         tests_res = []
#         box_pairs_all = []
#         statist = True
#         dicts_stats = stats_set[1]
#         alpha = stats_set[2]
#         txt_multcomp_all = []
        
#     print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    
#     # Genotypes and Strains being plotted 
#     values = props_ordered(df2plot, x_var, hue_var, shape_var)
#     x_values, hue_values, shape_values = values
#     # print(x_values, hue_values, shape_values)
    
#     # Set up the matplotlib figure 
#     if statsTxt:
#         n_rows = 2
#         height_ratios = [1,1]
#     else: 
#         n_rows = 1
#         height_ratios = [1]
#     n_cols = len(hue_values)

#     addVars = False
#     for index in range(n_cols):
#         if index == 0 and len(vars2plot) == 1:
#             addVars = True
#             var2copy = vars2plot[0]
#             label2copy = labels2plot[0]
        
#         if addVars and index != 0:
#             vars2plot.insert(index, var2copy)
#             labels2plot.insert(index, label2copy)
    
#     # As a right legend wants to be added, we need to add '' to all labels 
#     index_no_plot = []
#     if right_legend: 
#         index_no_plot = [n_cols]
#         for index in index_no_plot:
#             vars2plot.insert(index, '')
#             labels2plot.insert(index, '')
#         width_ratios = [1]*n_cols+[0.2]
#         n_cols = n_cols+1
#     else: 
#         width_ratios = [1]*n_cols
    
#     # print('index_no_plot:', index_no_plot)
#     index_no_plot2 = []
#     if statsTxt:
#         index_no_plot2 = [n_cols+aaa for aaa in range(n_cols)]
#         for index in index_no_plot2:
#             vars2plot.insert(index, '')
#             labels2plot.insert(index, '')
    
#     # print(vars2plot, labels2plot)
#     # print('index_no_plot2:', index_no_plot2)
    
#     index_plot = list(range(n_cols*n_rows))
#     index_plot = list(set(index_plot) - set(index_no_plot) - set(index_no_plot2))
#     # print('index_plot:', index_plot)
    
#     # Set up the matplotlib figure
#     h_plot, w_plot = ips
#     h_add = 1; w_add = 1

#     if len(x_values) > 3:
#         h_add = 0
#         wspace = 0.05
#         box_width = 0.6
#         bbox_lf = -0.5
#     else:
#         h_plot = 2#0.8
#         wspace = 0.15
#         box_width = 0.6
#         bbox_lf = -1.8
    
#     size_col = (n_cols)*(h_plot*len(x_values))+h_add
#     size_row = n_rows*w_plot+w_add
    
#     # Define legends for x and hue
#     v_legend = []; v_color = []
#     for zz, v_var, v_values in zip(count(), [x_var, hue_var], [x_values, hue_values]):
#         if isinstance(dict_legends[v_var], dict):
#             legend = []
#             color = []
#             for v_dict in v_values: 
#                 legend.append(dict_legends[v_var][v_dict]['legend'])
#                 color.append(dict_legends[v_var][v_dict]['color'])
#         else: 
#             legend = dict_legends[v_var]
#             color = ''
#         v_legend.append(legend)
#         v_color.append(color)
#     x_legend, hue_legend = v_legend
#     x_color, _ = v_color
        
#     if right_legend:
#         # Define legends for x
#         legend_elem = []
#         for aa, xval, xcol in zip(count(), x_legend, x_color):
#             legend_elem.append(Line2D([0], [0], marker='o', color='w', label=xval,
#                                     markerfacecolor=xcol, markersize=12))
#         handle_new = legend_elem
        
#     for index in index_no_plot+index_no_plot2:
#         hue_legend.insert(index, '')
#         hue_values.insert(index, '')
#         dicts_stats.insert(index,'')
            
#     marker_size = 8; dodge = True; jitter = 0.3
    
#     for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
#         print('\t- '+svar+': ', value)
        
#     ##  CREATE FIGURE 
#     gridkw = dict(width_ratios=width_ratios, height_ratios=height_ratios)
#     fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
#     fig.subplots_adjust(hspace=0.25, wspace=wspace)
    
#     sns.set_style("ticks")
#     sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    
#     # print(n_cols, n_rows, axes, right_legend, vars2plot, labels2plot, hue_values)
#     if n_cols == 1 and n_rows == 1 and not right_legend and not statsTxt:
#         axes_fl = [axes]
#     else: 
#         axes_fl = axes.flatten()
    
#     if statsTxt:
#         index_no_plot_all = index_no_plot+[index_no_plot2[-1]]
#     else:
#         index_no_plot_all = index_no_plot
#     # Create list to contain info of yticks and its highest value
#     # y_vals_all = []; 
#     max_y_vals = []; min_y_vals = []
#     for n, ax, var, ylabel, hue_value, dict_stats in zip(count(), axes_fl, vars2plot, labels2plot, hue_values, dicts_stats):
#         # print(n, len(axes_fl))
#         if n in index_no_plot_all:
#                 if n == n_cols-1 and right_legend:
#                     ax.set_axis_off()
#                     ax.legend(handle_new, x_legend, loc='upper left', bbox_to_anchor=(bbox_lf, 1), frameon = False)
#                 else: 
#                     ax.remove()
#         elif n in index_no_plot2 and statsTxt:
#             # print(n, 'aja')# elif n == 1 and statsTxt and not right_legend: 
#             ax.text(0.5,10,txt_multcomp_all[n-n_cols], size=11, ha='center', va='baseline', fontname = fontname)
#             ax.set_axis_off()
            
#         else: 
#             df_xfilt = df2plot[df2plot[hue_var] == hue_value]
#             m = sns.stripplot(data=df_xfilt, x=x_var, y=var, ax = ax, order=x_values,
#                               marker = 'o', palette = x_color, size = marker_size, 
#                               linewidth=0.5, jitter = jitter)
#             if box_plot: 
#                 m = sns.boxplot(data=df_xfilt, x=x_var, y=var, ax = ax, order=x_values,
#                                    dodge = dodge, width= box_width, showfliers = show_fliers, palette = x_color,
#                                    boxprops = boxprops, whiskerprops = whiskerprops, capprops = capprops, 
#                                    flierprops = flierprops, medianprops = medianprops, 
#                                    meanline = True, meanprops = meanprops, showmeans = True)
    
#             box = ax.get_position()
            
#             if statist: 
#                 box_pairs = dict_stats[var]['box_pairs']
#                 stat_test = False
#                 test = None
#                 p_val = dict_stats[var]['pval_multComp_all']
#                 txt_testSelected_all = dict_stats[var]['txt_testSelected_all']
#                 txt_normtest_all = dict_stats[var]['txt_normtest_all']
                
#                 # print(box_pairs, '\n ',p_val, '\n ',txt_testSelected_all, '\n ',txt_normtest_all)
#                 # Added fontname as an imput to add_stat_annotation function in statannot
#                 # m1, test_results = add_stat_annotation(ax, data=df_xfilt, x=x_var, y=var, hue=hue_var, hue_order = hue_values,
#                 #                                        order = x_values,
#                 #                                        box_pairs = box_pairs, perform_stat_test = stat_test, test = test,
#                 #                                        pvalues = p_val, comparisons_correction=None, #'bonferroni',
#                 #                                        line_offset_to_box=0.4, line_offset=0.2,
#                 #                                        line_height=0.015, text_offset=5,
#                 #                                        text_format='star', loc='inside', fontsize='x-small', verbose=0,
#                 #                                        fontname = fontname)
#                 tests_res.append(p_val)
#                 box_pairs_all.append(box_pairs)
#                 txt_multcomp = txtMultComp(box_pairs, p_val, hue_values[n], hue_legend[n], txt_testSelected_all, 
#                                             txt_normtest_all, alpha, indiv = True)
#                 txt_multcomp_all.append(txt_multcomp)
#                 # print(n, txt_multcomp)
                
#             if n == 0:
#                 if yticks_lab == '1e6 - d.':
#                     ylabel = ylabel +' x 10$^6$'
#                 elif yticks_lab == '1e3 - d.':
#                     ylabel = ylabel +' x 10$^3$'
#                 ax.set_xlabel(hue_legend[n], fontname = fontname)
#                 ax.set_ylabel(ylabel, fontname = fontname)
#             else:
#                 ax.set_xlabel(hue_legend[n], fontname = fontname)
#                 ax.set_ylabel('', fontname = fontname)
    
#             ax.set_position([box.x0, box.y0, box.width*1, box.height])
#             xticks = ax.get_xticks()
#             ax.set_xticks(xticks)
#             ax.set_xticklabels(x_legend, rotation=45, horizontalalignment='right', fontname = fontname)
#             sns.despine()
            
#             if isinstance(ylim, tuple):
#                 addY = (ylim[1] - ylim[0])/10 
#                 print('ylim:', ylim[0],'-', ylim[1]+ addY)
#                 ax.set_ylim(ylim[0], ylim[1]+addY)
            
#             if n == 0:
#                 handles, labels = m.get_legend_handles_labels()
#             else:
#                 ax.spines['left'].set_linestyle('-.')
#                 ax.spines['left'].set_color('#696969')
#                 ax.spines['left'].set_linewidth(0.8)
#                 ax.tick_params(left = False)
    
#             y_vals = ax.get_yticks()
#             # print(y_vals)
#             max_y_vals.append(y_vals[-1])
#             min_y_vals.append(y_vals[0])
    
#     # # Find min and maximum value for each axis 
#     # print(min_y_vals, max_y_vals)
#     min_y = min(min_y_vals); max_y = max(max_y_vals)

#     if yset == 'round':
#         rounding_set = -1
#         y_step = round((max_y - min_y)/8,rounding_set)
#         if y_step == 0:
#             y_step = round((max_y - min_y)/8,0)
#         # print(y_step)
#         y_vals_set = list(range(int(min_y), int(max_y), int(y_step)))
#         y_vals_set.append(y_vals_set[-1]+int(y_step))
        
#     elif yset == 'dec':
#         rounding_set = 2
#         y_step = round((max_y - min_y)/8,rounding_set)
#         if y_step == 0:
#             y_step = round((max_y - min_y)/8,3)
#         # print(y_step)
#         num_y = (max_y - min_y)//y_step
#         # print(num_y)
#         y_vals_set = np.linspace(min_y, max_y,int(num_y), endpoint = True)
#     # print(y_vals_set)
    
#     y_valsF = y_vals_set#y_vals_all[max_y_index]
    
#     for n, ax in zip(count(), axes_fl[index_plot]):
#         ax.set_yticks(y_valsF)
#         # Define axes based on higherst bar
#         if yticks_lab == '1e6 - d.':
#             ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_valsF])
#         elif yticks_lab == '1e3 - d.':
#             ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_valsF])
#         elif yticks_lab == 'th,':
#             ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_valsF])
#         elif yticks_lab == 'd. - 0':
#             ax.set_yticklabels(['{:.0f}'.format(w) for w in y_valsF])
#         elif yticks_lab == 'd. - 1':
#             ax.set_yticklabels(['{:.1f}'.format(w) for w in y_valsF])
#         elif yticks_lab == 'd.':
#             ax.set_yticklabels(['{:.2f}'.format(w) for w in y_valsF])
        
#         for tick in ax.get_xticklabels():
#             tick.set_fontname(fontname)
#         for tick in ax.get_yticklabels():
#             tick.set_fontname(fontname)
#         try:
#             ax.yaxis.get_offset_text().set_fontname(fontname)
#             # print(r' -> YAY')
#         except:
#             print(r' -> Nop')
    
#     # if statsTxt and n_cols > 1:
#     #     gs = axes[1,0].get_gridspec()
#     #     for axis in axes[1,:]:
#     #         axis.remove()
#     #     axbig = fig.add_subplot(gs[1,:])
    
#     #     axbig.annotate(txt_multcomp, (0.5, 0.225),
#     #                xycoords='figure fraction', va='center', fontname = fontname, fontsize = 14, 
#     #                horizontalalignment = 'center')
#     #     axbig.axis('off')
    
#     if suptitle:
#         fig.suptitle(title+'\n', y=1, fontname = fontname)
#     # fig.tight_layout()
    
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'R_')
#             if info != '':
#                 fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
#             else: 
#                 fig_title = dir2savef+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf

#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
            
#     if statist: 
#         return tests_res, box_pairs_all

#%% func - plotIndivStatTxtsperX
# def plotIndivStatsTxtperX (# General Plot Settings
#                     df2plot, vars2plot, x_var, hue_var, shape_var, 
#                     # Size, title, labels and legends
#                     title, labels2plot, dict_legends, ips, 
#                     # Statistic Settings
#                     stats_set, statsTxt = True, 
#                     # Other plot settings
#                     suptitle = True, right_legend = True,
#                     yticks_lab = 'th,', ylim = '', yset = '', box_plot = True, show_fliers = False,
#                     # Saving settings
#                     save = True, dpi = 300, ext = 'png', info = '', dir2save = ''):

#     if stats_set[0]: 
#         tests_res = []
#         box_pairs_all = []
#         statist = True
#         dicts_stats = stats_set[1]
#         alpha = stats_set[2]
#         txt_multcomp_all = []
        
#     print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    
#     # Genotypes and Strains being plotted 
#     values = props_ordered(df2plot, x_var, hue_var, shape_var)
#     x_values, hue_values, shape_values = values
#     # print(x_values, hue_values, shape_values)
    
#     # Set up the matplotlib figure 
#     if statsTxt:
#         n_rows = 2
#         height_ratios = [1,1]
#     else: 
#         n_rows = 1
#         height_ratios = [1]
#     n_cols = len(hue_values)

#     addVars = False
#     for index in range(n_cols):
#         if index == 0 and len(vars2plot) == 1:
#             addVars = True
#             var2copy = vars2plot[0]
#             label2copy = labels2plot[0]
        
#         if addVars and index != 0:
#             vars2plot.insert(index, var2copy)
#             labels2plot.insert(index, label2copy)
    
#     # As a right legend wants to be added, we need to add '' to all labels 
#     index_no_plot = []
#     if right_legend: 
#         index_no_plot = [n_cols]
#         for index in index_no_plot:
#             vars2plot.insert(index, '')
#             labels2plot.insert(index, '')
#         width_ratios = [1]*n_cols+[0.2]
#         n_cols = n_cols+1
#     else: 
#         width_ratios = [1]*n_cols
    
#     # print('index_no_plot:', index_no_plot)
#     index_no_plot2 = []
#     if statsTxt:
#         index_no_plot2 = [n_cols+aaa for aaa in range(n_cols)]
#         for index in index_no_plot2:
#             vars2plot.insert(index, '')
#             labels2plot.insert(index, '')
    
#     # print(vars2plot, labels2plot)
#     # print('index_no_plot2:', index_no_plot2)
    
#     index_plot = list(range(n_cols*n_rows))
#     index_plot = list(set(index_plot) - set(index_no_plot) - set(index_no_plot2))
#     # print('index_plot:', index_plot)
    
#     # Set up the matplotlib figure
#     h_plot, w_plot = ips
#     h_add = 1; w_add = 1

#     if len(x_values) > 3:
#         h_add = 0
#         wspace = 0.05
#         box_width = 0.6
#         bbox_lf = -0.5
#     else:
#         h_plot = 2#0.8
#         wspace = 0.15
#         box_width = 0.6
#         bbox_lf = -1.8
    
#     size_col = (n_cols)*(h_plot*len(x_values))+h_add
#     size_row = n_rows*w_plot+w_add
    
#     # Define legends for x and hue
#     v_legend = []; v_color = []
#     for zz, v_var, v_values in zip(count(), [x_var, hue_var], [x_values, hue_values]):
#         if isinstance(dict_legends[v_var], dict):
#             legend = []
#             color = []
#             for v_dict in v_values: 
#                 legend.append(dict_legends[v_var][v_dict]['legend'])
#                 color.append(dict_legends[v_var][v_dict]['color'])
#         else: 
#             legend = dict_legends[v_var]
#             color = ''
#         v_legend.append(legend)
#         v_color.append(color)
#     x_legend, hue_legend = v_legend
#     x_color, _ = v_color
        
#     if right_legend:
#         # Define legends for x
#         legend_elem = []
#         for aa, xval, xcol in zip(count(), x_legend, x_color):
#             legend_elem.append(Line2D([0], [0], marker='o', color='w', label=xval,
#                                     markerfacecolor=xcol, markersize=12))
#         handle_new = legend_elem
        
#     for index in index_no_plot+index_no_plot2:
#         hue_legend.insert(index, '')
#         hue_values.insert(index, '')
#         dicts_stats.insert(index,'')
            
#     marker_size = 8; dodge = True; jitter = 0.3
    
#     for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
#         print('\t- '+svar+': ', value)
        
#     ##  CREATE FIGURE 
#     gridkw = dict(width_ratios=width_ratios, height_ratios=height_ratios)
#     fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(size_col, size_row), sharex=False, sharey=False, gridspec_kw=gridkw)
#     fig.subplots_adjust(hspace=0.25, wspace=wspace)
    
#     sns.set_style("ticks")
#     sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    
#     # print(n_cols, n_rows, axes, right_legend, vars2plot, labels2plot, hue_values)
#     if n_cols == 1 and n_rows == 1 and not right_legend and not statsTxt:
#         axes_fl = [axes]
#     else: 
#         axes_fl = axes.flatten()
    
#     if statsTxt:
#         index_no_plot_all = index_no_plot+[index_no_plot2[-1]]
#     else:
#         index_no_plot_all = index_no_plot
#     # Create list to contain info of yticks and its highest value
#     # y_vals_all = []; 
#     max_y_vals = []; min_y_vals = []
#     for n, ax, var, ylabel, hue_value, dict_stats in zip(count(), axes_fl, vars2plot, labels2plot, hue_values, dicts_stats):
#         # print(n, len(axes_fl))
#         if n in index_no_plot_all:
#                 if n == n_cols-1 and right_legend:
#                     ax.set_axis_off()
#                     ax.legend(handle_new, x_legend, loc='upper left', bbox_to_anchor=(bbox_lf, 1), frameon = False)
#                 else: 
#                     ax.remove()
#         elif n in index_no_plot2 and statsTxt:
#             # print(n, 'aja')# elif n == 1 and statsTxt and not right_legend: 
#             ax.set_ylim(0,1)
#             ax.set_xlim(0,1)
#             ax.text(0.5,0.5,txt_multcomp_all[n-n_cols], size=11, ha='center', va='baseline', fontname = fontname)
#             ax.set_axis_off()
            
#         else: 
#             df_xfilt = df2plot[df2plot[hue_var] == hue_value]
#             m = sns.stripplot(data=df_xfilt, x=x_var, y=var, ax = ax, order=x_values,
#                               marker = 'o', palette = x_color, size = marker_size, 
#                               linewidth=0.5, jitter = jitter)
#             if box_plot: 
#                 m = sns.boxplot(data=df_xfilt, x=x_var, y=var, ax = ax, order=x_values,
#                                    dodge = dodge, width= box_width, showfliers = show_fliers, palette = x_color,
#                                    boxprops = boxprops, whiskerprops = whiskerprops, capprops = capprops, 
#                                    flierprops = flierprops, medianprops = medianprops, 
#                                    meanline = True, meanprops = meanprops, showmeans = True)
    
#             box = ax.get_position()
            
#             if statist: 
#                 box_pairs = dict_stats[var]['box_pairs']
#                 stat_test = False
#                 test = None
#                 p_val = dict_stats[var]['pval_multComp_all']
#                 txt_testSelected_all = dict_stats[var]['txt_testSelected_all']
#                 txt_normtest_all = dict_stats[var]['txt_normtest_all']
                
#                 # print(box_pairs, '\n ',p_val, '\n ',txt_testSelected_all, '\n ',txt_normtest_all)
#                 # Added fontname as an imput to add_stat_annotation function in statannot
#                 # m1, test_results = add_stat_annotation(ax, data=df_xfilt, x=x_var, y=var, hue=hue_var, hue_order = hue_values,
#                 #                                        order = x_values,
#                 #                                        box_pairs = box_pairs, perform_stat_test = stat_test, test = test,
#                 #                                        pvalues = p_val, comparisons_correction=None, #'bonferroni',
#                 #                                        line_offset_to_box=0.4, line_offset=0.2,
#                 #                                        line_height=0.015, text_offset=5,
#                 #                                        text_format='star', loc='inside', fontsize='x-small', verbose=0,
#                 #                                        fontname = fontname)
#                 tests_res.append(p_val)
#                 box_pairs_all.append(box_pairs)
#                 txt_multcomp = txtMultComp(box_pairs, p_val, hue_values[n], hue_legend[n], txt_testSelected_all, 
#                                             txt_normtest_all, alpha, indiv = True)
#                 txt_multcomp_all.append(txt_multcomp)
#                 # print(n, txt_multcomp)
                
#             if n == 0:
#                 if yticks_lab == '1e6 - d.':
#                     ylabel = ylabel +' x 10$^6$'
#                 elif yticks_lab == '1e3 - d.':
#                     ylabel = ylabel +' x 10$^3$'
#                 ax.set_xlabel(hue_legend[n], fontname = fontname)
#                 ax.set_ylabel(ylabel, fontname = fontname)
#             else:
#                 ax.set_xlabel(hue_legend[n], fontname = fontname)
#                 ax.set_ylabel('', fontname = fontname)
    
#             ax.set_position([box.x0, box.y0, box.width*1, box.height])
#             xticks = ax.get_xticks()
#             ax.set_xticks(xticks)
#             ax.set_xticklabels(x_legend, rotation=45, horizontalalignment='right', fontname = fontname)
#             sns.despine()
            
#             if isinstance(ylim, tuple):
#                 addY = (ylim[1] - ylim[0])/10 
#                 print('ylim:', ylim[0],'-', ylim[1]+ addY)
#                 ax.set_ylim(ylim[0], ylim[1]+addY)
            
#             if n == 0:
#                 handles, labels = m.get_legend_handles_labels()
#             else:
#                 ax.spines['left'].set_linestyle('-.')
#                 ax.spines['left'].set_color('#696969')
#                 ax.spines['left'].set_linewidth(0.8)
#                 ax.tick_params(left = False)
    
#             y_vals = ax.get_yticks()
#             # print(y_vals)
#             max_y_vals.append(y_vals[-1])
#             min_y_vals.append(y_vals[0])
    
#     # # Find min and maximum value for each axis 
#     # print(min_y_vals, max_y_vals)
#     min_y = min(min_y_vals); max_y = max(max_y_vals)

#     if yset == 'round':
#         rounding_set = -1
#         y_step = round((max_y - min_y)/8,rounding_set)
#         if y_step == 0:
#             y_step = round((max_y - min_y)/8,0)
#         # print(y_step)
#         y_vals_set = list(range(int(min_y), int(max_y), int(y_step)))
#         y_vals_set.append(y_vals_set[-1]+int(y_step))
        
#     elif yset == 'dec':
#         rounding_set = 2
#         y_step = round((max_y - min_y)/8,rounding_set)
#         if y_step == 0:
#             y_step = round((max_y - min_y)/8,3)
#         # print(y_step)
#         num_y = (max_y - min_y)//y_step
#         # print(num_y)
#         y_vals_set = np.linspace(min_y, max_y,int(num_y), endpoint = True)
#     # print(y_vals_set)
    
#     y_valsF = y_vals_set#y_vals_all[max_y_index]
    
#     for n, ax in zip(count(), axes_fl[index_plot]):
#         ax.set_yticks(y_valsF)
#         # Define axes based on higherst bar
#         if yticks_lab == '1e6 - d.':
#             ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_valsF])
#         elif yticks_lab == '1e3 - d.':
#             ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_valsF])
#         elif yticks_lab == 'th,':
#             ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_valsF])
#         elif yticks_lab == 'd. - 0':
#             ax.set_yticklabels(['{:.0f}'.format(w) for w in y_valsF])
#         elif yticks_lab == 'd. - 1':
#             ax.set_yticklabels(['{:.1f}'.format(w) for w in y_valsF])
#         elif yticks_lab == 'd.':
#             ax.set_yticklabels(['{:.2f}'.format(w) for w in y_valsF])
        
#         for tick in ax.get_xticklabels():
#             tick.set_fontname(fontname)
#         for tick in ax.get_yticklabels():
#             tick.set_fontname(fontname)
#         try:
#             ax.yaxis.get_offset_text().set_fontname(fontname)
#             # print(r' -> YAY')
#         except:
#             print(r' -> Nop')
    
#     # if statsTxt and n_cols > 1:
#     #     gs = axes[1,0].get_gridspec()
#     #     for axis in axes[1,:]:
#     #         axis.remove()
#     #     axbig = fig.add_subplot(gs[1,:])
    
#     #     axbig.annotate(txt_multcomp, (0.5, 0.225),
#     #                xycoords='figure fraction', va='center', fontname = fontname, fontsize = 14, 
#     #                horizontalalignment = 'center')
#     #     axbig.axis('off')
    
#     if suptitle:
#         fig.suptitle(title+'\n', y=1, fontname = fontname)
#     # fig.tight_layout()
    
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'R_')
#             if info != '':
#                 fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
#             else: 
#                 fig_title = dir2savef+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf

#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
            
#     if statist: 
#         return tests_res, box_pairs_all

#%% func - plotIndivStatTxtAllX
def plotIndivStatsTxtAllX (# General Plot Settings
                    df2plot, vars2plot, x_var, hue_var, shape_var, 
                    # Size, title, labels and legends
                    title, dict_legends, 
                    # Statistic Settings
                    stats_set, 
                    # Other plot settings
                    suptitle = True, ctx_style = 0,
                    # Saving settings
                    save = True, dpi = 300, ext = 'png', info = '', dir2save = ''):

    contxt, font_scale, rc_dict, all_box_props, fontname  = setSNSContext(ctx_style)
    boxprops, flierprops, whiskerprops, capprops, medianprops, meanprops = all_box_props
    
    if stats_set[0]: 
        tests_res = []
        box_pairs_all = []
        statist = True
        dicts_stats = stats_set[1]
        alpha = stats_set[2]
        txt_multcomp_all = []
        
    print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    
    # Genotypes and Strains being plotted 
    values = props_ordered(df2plot, x_var, hue_var, shape_var)
    x_values, hue_values, shape_values = values
    
    n_cols = len(hue_values)

    addVars = False
    for index in range(n_cols):
        if index == 0 and len(vars2plot) == 1:
            addVars = True
            var2copy = vars2plot[0]
        
        if addVars and index != 0:
            vars2plot.insert(index, var2copy)
    
    # Define legends for x and hue
    v_legend = []
    for zz, v_var, v_values in zip(count(), [x_var, hue_var], [x_values, hue_values]):
        if isinstance(dict_legends[v_var], dict):
            legend = []
            for v_dict in v_values: 
                legend.append(dict_legends[v_var][v_dict]['legend'])
        else: 
            legend = dict_legends[v_var]
        v_legend.append(legend)
    x_legend, hue_legend = v_legend
    
    for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
        print('\t- '+svar+': ', value)
    
    sns.set_style("ticks")
    sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    
    for n, var, dict_stats in zip(count(), vars2plot, dicts_stats):
        box_pairs = dict_stats[var]['box_pairs']
        # stat_test = False
        # test = None
        p_val = dict_stats[var]['pval_multComp_all']
        txt_testSelected_all = dict_stats[var]['txt_testSelected_all']
        txt_normtest_all = dict_stats[var]['txt_normtest_all']

        tests_res.append(p_val)
        box_pairs_all.append(box_pairs)
        txt_multcomp = txtMultComp(box_pairs, p_val, hue_values[n], hue_legend[n], txt_testSelected_all, 
                                    txt_normtest_all, alpha, indiv = True)
        txt_multcomp_all.append(txt_multcomp)
    
    txt_multcomp_f = '\n'.join(txt_multcomp_all)
    txt_multcomp_f = title+'\n-.-.-.-.-.-.-.-.-.-.-.-\n'+txt_multcomp_f
    print(txt_multcomp_f)
    
    fig = plt.figure()
    ax = fig.add_subplot()
    fig.subplots_adjust(top=0.85)
    # fig.suptitle(title+'\n', y=1.2, fontname = fontname, fontsize=14)
    ax.text(0.5,0.5,txt_multcomp_f, size=11, ha='center', va='center', fontname = fontname)
    ax.set_axis_off()
    # plt.text(0.5, 0.5, txt_multcomp, fontname = fontname, fontsize = 14, color='green')
    # plt.show()

    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'R_')
            # fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")_stats."+extf
            fig_title = dir2savef+info+"_"+title+"_stats."+extf
            plt.savefig(fig_title, dpi=100, bbox_inches='tight', transparent=True)
            
    plt.show()
            
    if statist: 
        return tests_res, box_pairs_all

#%% func - plotInGroupsStats
# def plotInGroupsStats(# General Plot Settings
#                       df2plot, vars2plot, x_var, hue_var, shape_var, 
#                       # Size, title, labels and legends
#                       title, labels2plot, dict_legends, ips, 
#                       # Statistic Settings
#                       stats_set,
#                       # Other plot settings
#                       n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th,', ylim = '', 
#                       # Saving settings
#                       save = True, dpi = 300, ext = 'png', info ='', dir2save=''):
    
#     if stats_set[0]: 
#         tests_res = []
#         box_pairs_all = []
#         statist = True
#         dict_stats = stats_set[1]
#         alpha = stats_set[2]
#         txt_multcomp_all = []
        
#     print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    
#     # Set up the matplotlib figure
#     num_vars = len(vars2plot)
#     n_rows = math.ceil(num_vars/n_cols)+1
    
#     index_right_col = list(range(n_cols,(n_cols+1)*n_rows,n_cols+1))
#     index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
#     index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))
#     index_stats = sorted(list(set(index_no_graph).difference(set(index_right_col))))

#     for index in index_no_plot:
#         vars2plot.insert(index, '')
#         labels2plot.insert(index, '')
    
#     # Set up the matplotlib figure
#     h_plot, w_plot = ips
    
#     # Genotypes and Strains being plotted 
#     values = props_ordered(df2plot, x_var, hue_var, shape_var)
#     x_values, hue_values, shape_values = values
    
#     # - number of x_var
#     n_x = len(x_values)
#     if n_x == 1:
#         h_plot = 6
#         hspace = 0.2; wspace = 1
#     else: 
#         hspace = 0.5; wspace = 0.5
        
#     if num_vars == 1:
#         h_add = 0; w_add = 0
#     size_col = (n_cols+1)*h_plot+h_add
#     size_row = n_rows*w_plot+w_add
    
#     # Define legends for x and hue
#     v_legend = []; v_color = []
#     for zz, v_var, v_values in zip(count(), [x_var, hue_var], [x_values, hue_values]):
#         if isinstance(dict_legends[v_var], dict):
#             legend = []
#             color = []
#             for v_dict in v_values: 
#                 legend.append(dict_legends[v_var][v_dict]['legend'])
#                 color.append(dict_legends[v_var][v_dict]['color'])
#         else: 
#             legend = dict_legends[v_var]
#             color = ''
#         v_legend.append(legend)
#         v_color.append(color)
#     x_legend, hue_legend = v_legend
#     x_color, hue_color = v_color

#     legend_elem_hue = []
#     for aa, hue_val, hue_col in zip(count(), hue_legend, hue_color):
#         legend_elem_hue.append(Line2D([0], [0], marker='o', color='w', label=hue_val,
#                                 markerfacecolor=hue_col, markersize=20))
        
#     handle_new = legend_elem_hue
#     legend_new = hue_legend

#     marker_size = 8; dodge = True; plot_no = 0
    
#     #  CREATE FIGURE
#     gridkw = dict(width_ratios=[1]*n_cols+[0.2], height_ratios = [1,1.2])
#     fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
#     fig.subplots_adjust(hspace=hspace, wspace=wspace)
    
#     sns.set_style("ticks")
#     sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    
#     for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
#         # Create list to contain info of yticks and ist highest value
#         y_vals_all = []
#         max_y_vals = []
    
#         if n == 0: 
#             for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
#                 print('\t- '+svar+': ', value)
        
#         if n in index_no_plot:
#             if n == n_cols:
#                 ax.set_axis_off()
#                 ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
#             else: 
#                 if n in index_stats:
#                     text2add = '- - - - - - - - - - \n STATISTICS: \n'+txt_multcomp_all[n-n_cols-1]
#                     ax.text(0.5, 0.5, text2add, size=16, ha='center', va='center', fontname = fontname)
#                     ax.set_axis_off()
#                 else: 
#                     ax.remove()
                
#         else: 
#             print(' > ',var,'\n')
#             m = sns.boxplot(data=df2plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
#                             dodge = dodge, width =0.8, showfliers = True, palette = hue_color,
#                             boxprops = boxprops, whiskerprops = whiskerprops, capprops = capprops, 
#                             flierprops = flierprops, medianprops = medianprops, 
#                             meanline = True, meanprops = meanprops, showmeans = True)
                                  
#             m = sns.stripplot(data=df2plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
#                               palette = hue_color, dodge = dodge, size = marker_size, jitter = 0.3)
            
#             if statist: 
#                 box_pairs = dict_stats[var]['box_pairs']
#                 stat_test = False
#                 test = None
#                 p_val = dict_stats[var]['pval_multComp_all']
#                 txt_testSelected_all = dict_stats[var]['txt_testSelected_all']
#                 txt_normtest_all = dict_stats[var]['txt_normtest_all']
                
#                 m1, test_results = add_stat_annotation(ax, data=df2plot, x=x_var, y=var, hue=hue_var, 
#                                                        hue_order = hue_values, order = x_values,
#                                                         box_pairs=box_pairs, perform_stat_test = stat_test, test = test,
#                                                         pvalues=p_val, comparisons_correction=None, #'bonferroni',
#                                                         line_offset_to_box=0.4, line_offset=0.1,
#                                                         line_height=0.015, text_offset=5,
#                                                         text_format='star', loc='inside', fontsize='small', verbose=0,
#                                                         fontname = fontname);
#                 tests_res.append(p_val)
#                 box_pairs_all.append(box_pairs)
#                 txt_multcomp = txtMultComp (box_pairs, p_val, x_values, x_legend, txt_testSelected_all, txt_normtest_all, alpha)
#                 txt_multcomp_all.append(txt_multcomp)
#                 print('txt_multcomp:', txt_multcomp)
                
#             box = ax.get_position()
#             if yticks_lab == '1e6 - d.':
#                 ylabel = ylabel +' x 10$^6$'
#             elif yticks_lab == '1e3 - d.':
#                 ylabel = ylabel +' x 10$^3$'
            
#             ax.set_xlabel(dict_legends['xlabels'][x_var], fontname = fontname)
#             ax.set_ylabel(ylabel, fontname = fontname)
#             ax.set_xticklabels(dict_legends[x_var], rotation=0, fontname = fontname)
#             ax.set_position([box.x0, box.y0, box.width*1, box.height])
#             ax.get_legend().remove()
#             sns.despine()
#             if ylim != '':
#                 ax.set_ylim(ylim[plot_no][0], ylim[plot_no][1])
            
#             if n == 0:
#                 handles, labels = m.get_legend_handles_labels()

#             y_vals = ax.get_yticks()
#             y_vals_all.append(y_vals)
#             max_y_vals.append(y_vals[-1])
#             if n_x > 1: 
#                 vline_pos = list(range(1,n_x,1))
#                 for pos in vline_pos:
#                     ax.axvline(pos - 0.5, ymax = 0.95, color='dimgrey', ls='-.', linewidth=0.8)
            
#             # Define axes based on higherst bar
#             max_y_index = max_y_vals.index(max(max_y_vals))
#             ax.set_yticks(y_vals_all[max_y_index])
#             if yticks_lab == '1e6 - d.':
#                 ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == '1e3 - d.':
#                 ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'th,':
#                 ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 0':
#                 ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 1':
#                 ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd.':
#                 ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])
            
#             for tick in ax.get_xticklabels():
#                 tick.set_fontname(fontname)
#             for tick in ax.get_yticklabels():
#                 tick.set_fontname(fontname)
#             try:
#                 ax.yaxis.get_offset_text().set_fontname(fontname)
#             except:
#                 print(r' -> Nop')
                
#             plot_no +=1

#     fig.suptitle(title+'\n', fontsize = 30, y=0.95, fontname = fontname)
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'pl_stats', 'R_')
#             if info != '':
#                 fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
#             else: 
#                 fig_title = dir2savef+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf

#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
    
#     if statist: 
#         return tests_res, box_pairs_all

#%% -> Time course - no statistics
#%% func - plotIndivTimeCourse
def plotIndivTimeCourse (# General Plot Settings
                    df2plot, vars2plot, x_var, hue_var,
                    # Size, title, labels and legends
                    title, labels2plot, dict_legends, ips, 
                    # Other plot settings
                    suptitle = True, right_legend = True, ctx_style = 0,
                    yticks_lab = 'th,', ylim = '', 
                    # Saving settings
                    save = True, dpi = 300, ext = 'png', info = '', dir2save = ''):

    contxt, font_scale, rc_dict, all_box_props, fontname  = setSNSContext(ctx_style)
    boxprops, flierprops, whiskerprops, capprops, medianprops, meanprops = all_box_props
    
    
    print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var)
    
    # Genotypes and Strains being plotted 
    values = props_ordered(df2plot, x_var, hue_var, 'Strain')
    x_values, hue_values, _ = values
    
    # Set up the matplotlib figure 
    n_rows = 1
    n_cols = 1
    
    # As a right legend wants to be added, we need to add '' to all labels 
    if right_legend: 
        index_no_plot = [n_cols]
        for index in index_no_plot:
            vars2plot.insert(index, '')
            labels2plot.insert(index, '')
        width_ratios = [1]*n_cols+[0.2]
        n_cols = n_cols+1
    else: 
        width_ratios = [1]*n_cols
        index_no_plot = [1000]
        
    # Set up the matplotlib figure
    h_plot, w_plot = ips
    h_add = 0; w_add = 0
    print(len(x_values), len(hue_values))
    if len(x_values) > 3:
        h_add = 0
        wspace = 0.05
        bbox_lf = -0.5
    else:
        # h_plot = 0.8
        wspace = 0.05
        bbox_lf = -0.4
        
    size_col = (n_cols)*(h_plot*len(x_values))+h_add
    size_row = n_rows*w_plot+w_add
    print('fig_size:', size_col, size_row)
    
    # Define legends for x and hue
    v_legend = []
    for zz, v_var, v_values in zip(count(), [x_var, hue_var], [x_values, hue_values]):
        if isinstance(dict_legends[v_var], dict):
            legend = []
            for v_dict in v_values: 
                legend.append(dict_legends[v_var][v_dict]['legend'])
        else: 
            legend = dict_legends[v_var]
        v_legend.append(legend)
    x_legend, hue_legend = v_legend
    # print('x_legend, hue_legend')
    # print(x_legend, hue_legend)
        
    if right_legend:
        # Define legends for x
        legend_elem = []
        palette_elem = []
        palette_black = []
        for aa, xval in zip(count(), hue_legend):
            legend_elem.append(Line2D([0], [0], color=eight_palette[aa], label=xval,lw = 1))
            palette_elem.append(eight_palette[aa])
            palette_black.append('#A9A9A9')
        handle_new = legend_elem
            
    marker_size = 3
    
    for k, svar, value in zip(count(), [x_var, hue_var], values):
        print('\t- '+svar+': ', value)
        
    ##  CREATE FIGURE 
    gridkw = dict(width_ratios=width_ratios)
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
    fig.subplots_adjust(hspace=1, wspace=wspace)
    
    sns.set_style("ticks")
    sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    
    # Create list to contain info of yticks and its highest value
    for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
        if n in index_no_plot and right_legend:
            if n == n_cols-1:
                ax.set_axis_off()
                ax.legend(handle_new, hue_legend, loc='upper left', bbox_to_anchor=(bbox_lf, 1), frameon = False)
            else: 
                ax.remove()
        else: 
            m = sns.lineplot(data=df2plot, x=x_var, y=var, hue=hue_var, ax = ax, 
                             markers=True, palette = palette_elem)
            m = sns.stripplot(data=df2plot, x=x_var, y=var, ax = ax, order=x_values,
                              marker = 'o',size = marker_size, palette = palette_black, 
                              linewidth=0.5, jitter = False)
            box = ax.get_position()
            ax.get_legend().remove()
            
            for num in range(len(hue_legend)):
                ax.lines[num].set_linestyle("--")
                ax.lines[num].set_linewidth(1)
    
            if n == 0:
                if yticks_lab == '1e6 - d.':
                    ylabel = ylabel +' x 10$^6$'
                elif yticks_lab == '1e3 - d.':
                    ylabel = ylabel +' x 10$^3$'
                xlabel = dict_legends['xlabels'][x_var]
                ax.set_xlabel(xlabel, fontname = fontname)
                ax.set_ylabel(ylabel, fontname = fontname)
            else:
                ax.set_xlabel('', fontname = fontname)
                ax.set_ylabel('', fontname = fontname)
    
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            xticks = ax.get_xticks()
            ax.set_xticks(xticks)
            ax.set_xticklabels(x_legend, rotation=45, horizontalalignment='center', fontname = fontname)
            sns.despine()
            
            if ylim != '':
                print('ylim:', ylim[0],'-', ylim[1])
                ax.set_ylim(ylim[0], ylim[1])
            
            if n == 0:
                handles, labels = m.get_legend_handles_labels()
            else:
                ax.spines['left'].set_linestyle('-.')
                ax.spines['left'].set_color('#696969')
                ax.spines['left'].set_linewidth(0.8)
                ax.tick_params(left = False)
    
            y_vals = ax.get_yticks()
            ax.set_yticks(y_vals)
    
            # Define axes based on higherst bar
            if yticks_lab == '1e6 - d.':
                ax.set_yticks(ax.get_yticks())  #
                ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals])
            elif yticks_lab == '1e3 - d.':
                ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals])
            elif yticks_lab == 'th,':
                ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals])
            elif yticks_lab == 'd. - 0':
                ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals])
            elif yticks_lab == 'd. - 1':
                ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals])
            elif yticks_lab == 'd.':
                ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals])
            
            for tick in ax.get_xticklabels():
                tick.set_fontname(fontname)
            for tick in ax.get_yticklabels():
                tick.set_fontname(fontname)
            try:
                ax.yaxis.get_offset_text().set_fontname(fontname)
                # print(r' -> YAY')
            except:
                print(r' -> Nop')

    if suptitle:
        fig.suptitle(title+'\n', y=1, fontname = fontname)
    # fig.tight_layout()
    
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'R_')
            fig_title = dir2savef+info+"_"+title+"."+extf


            plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)

#%% func - plotMultTimeCourse
def plotMultTimeCourse (# General Plot Settings
                    df2plot, vars2plot, x_var, hue_var,
                    # Size, title, labels and legends
                    title, labels2plot, dict_legends, ips, 
                    # Other plot settings
                    suptitle = True, right_legend = True, ctx_style = 0,
                    yticks_lab = 'th,', ylim = '', 
                    # Saving settings
                    save = True, dpi = 300, ext = 'png', info = '', dir2save = ''):

    contxt, font_scale, rc_dict, all_box_props, fontname  = setSNSContext(ctx_style)
    boxprops, flierprops, whiskerprops, capprops, medianprops, meanprops = all_box_props
    
    print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var)
    
    # Genotypes and Strains being plotted 
    values = props_ordered(df2plot, x_var, hue_var, 'Strain')
    x_values, hue_values, _ = values
    
    # Set up the matplotlib figure 
    n_rows = 1
    n_cols = 1
    
    # As a right legend wants to be added, we need to add '' to all labels 
    if right_legend: 
        index_no_plot = [n_cols]
        for index in index_no_plot:
            vars2plot.insert(index, '')
            labels2plot.insert(index, '')
        width_ratios = [1]*n_cols+[0.2]
        n_cols = n_cols+1
    else: 
        width_ratios = [1]*n_cols
        index_no_plot = [1000]
        
    # Set up the matplotlib figure
    h_plot, w_plot = ips
    h_add = 0; w_add = 0
    print(len(x_values), len(hue_values))
    if len(x_values) > 3:
        h_add = 0
        wspace = 0.05
        bbox_lf = -0.5
    else:
        # h_plot = 0.8
        wspace = 0.05
        bbox_lf = -0.4
        
    size_col = (n_cols)*(h_plot*len(x_values))+h_add
    size_row = n_rows*w_plot+w_add
    print('fig_size:', size_col, size_row)
    print(x_values, hue_values)
    # Define legends for x and hue
    v_legend = []; v_color = []
    for zz, v_var, v_values in zip(count(), [x_var, hue_var], [x_values, hue_values]):
        if isinstance(dict_legends[v_var], dict):
            legend = []
            color = []
            for v_dict in v_values: 
                legend.append(dict_legends[v_var][v_dict]['legend'])
                color.append(dict_legends[v_var][v_dict]['color'])
        else: 
            legend = dict_legends[v_var]
            color = ''
        v_legend.append(legend)
        v_color.append(color)
    x_legend, hue_legend = v_legend
    x_color, hue_color = v_color
    # print('x_legend, hue_legend')
    # print(x_legend, hue_legend)
    # print(x_color, hue_color)
        
    if right_legend:
        # Define legends for x
        legend_elem = []
        for aa, xval, xcol in zip(count(), hue_legend, hue_color):
            legend_elem.append(Line2D([0], [0], color=xcol, label=xval,lw = 1))
        handle_new = legend_elem
    
    for k, svar, value in zip(count(), [x_var, hue_var], values):
        print('\t- '+svar+': ', value)
        
    ##  CREATE FIGURE 
    gridkw = dict(width_ratios=width_ratios)
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
    fig.subplots_adjust(hspace=1, wspace=wspace)
    
    sns.set_style("ticks")
    sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    
    # Create list to contain info of yticks and its highest value
    for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
        print(var, ylabel)
        if n in index_no_plot and right_legend:
            if n == n_cols-1:
                ax.set_axis_off()
                ax.legend(handle_new, hue_legend, loc='upper left', bbox_to_anchor=(bbox_lf, 1), frameon = False, fontsize=7)
            else: 
                ax.remove()
        else: 
            # m = sns.lineplot(data=df2plot, x=x_var, y=var, hue=hue_var, ax = ax, 
            #                  markers=True, ci = None, palette = palette_elem)
            m = sns.pointplot(data=df2plot, x=x_var, y=var, hue=hue_var, ax = ax, order=x_values,
                              hue_order = hue_values,
                              scale = 0.15, 
                              marker = 'o', palette = hue_color,#palette_elem, 
                              dodge=True, linestyles='--', linewidth = 0.1, errwidth=0.4, capsize=.12)
            box = ax.get_position()
            # plt.setp(ax.collections, sizes=[10])
            ax.get_legend().remove()
            # lw = m.lines[0].get_linewidth() # lw of first line
            plt.setp(m.lines,linewidth=0.4) 
            # for num in range(len(hue_legend)):
            #     ax.lines[num].set_linestyle("--")
            #     ax.lines[num].set_linewidth(1)
    
            if n == 0:
                if yticks_lab == '1e6 - d.':
                    ylabel = ylabel +' x 10$^6$'
                elif yticks_lab == '1e3 - d.':
                    ylabel = ylabel +' x 10$^3$'
                xlabel = dict_legends['xlabels'][x_var]
                ax.set_xlabel(xlabel, fontname = fontname)
                ax.set_ylabel(ylabel, fontname = fontname)
            else:
                ax.set_xlabel('', fontname = fontname)
                ax.set_ylabel('', fontname = fontname)
    
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            xticks = ax.get_xticks()
            ax.set_xticks(xticks)
            ax.set_xticklabels(x_legend, rotation=45, horizontalalignment='center', fontname = fontname)
            sns.despine()
            
            if ylim != '':
                print('ylim:', ylim[0],'-', ylim[1])
                ax.set_ylim(ylim[0], ylim[1])
            
            if n == 0:
                handles, labels = m.get_legend_handles_labels()
            else:
                ax.spines['left'].set_linestyle('-.')
                ax.spines['left'].set_color('#696969')
                ax.spines['left'].set_linewidth(0.8)
                ax.tick_params(left = False)
    
            y_vals = ax.get_yticks()
            ax.set_yticks(y_vals)
    
            # Define axes based on higherst bar
            if yticks_lab == '1e6 - d.':
                ax.set_yticks(ax.get_yticks())  #
                ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals])
            elif yticks_lab == '1e3 - d.':
                ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals])
            elif yticks_lab == 'th,':
                ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals])
            elif yticks_lab == 'd. - 0':
                ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals])
            elif yticks_lab == 'd. - 1':
                ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals])
            elif yticks_lab == 'd.':
                ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals])
            
            for tick in ax.get_xticklabels():
                tick.set_fontname(fontname)
            for tick in ax.get_yticklabels():
                tick.set_fontname(fontname)
            try:
                ax.yaxis.get_offset_text().set_fontname(fontname)
                # print(r' -> YAY')
            except:
                print(r' -> Nop')

    if suptitle:
        fig.suptitle(title+'\n', y=1, fontname = fontname)
    # fig.tight_layout()
    
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'R_')
            fig_title = dir2savef+info+"_"+title+"."+extf


            plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
            
#%% -> RelPlots
#%% func - relPlotInGroups (test if working!)
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
            
#%% func - barPlots (Working!)
def barPlots(# General Plot Settings
             df2plot, vars2plot, group_vars, x_var, col_var,
             # General Plot Settings
             title, txt_title, ylabel, colours, dict_legends,
             # Other plot settings
             yticks_lab, ylim, stack100 = True, sub_bar_lab = True, 
             ctx_style = 0, bot_legend = True, 
             # Saving settings
             save = True, ext = ['png'], info='', dir2save = ''): 
    
    
    contxt, font_scale, rc_dict, all_box_props, fontname  = setSNSContext(ctx_style)
    boxprops, flierprops, whiskerprops, capprops, medianprops, meanprops = all_box_props
    
    mpl.rcParams.update(mpl.rcParamsDefault)
    print('\n> '+title)
    print('\t - x_var: '+x_var+ ' - col_var:'+col_var)
    # Create list to contain info of bars and ist highest value
    y_vals_all = []
    max_y_vals = []
    
    # Filter dataframe 
    df4plot = df2plot.groupby(group_vars)[vars2plot].mean()
    df4plot_count = df2plot.groupby(group_vars)[vars2plot].count()
    # dict_legend = def_legends(df2plot.reset_index())
    
    # Genotypes and Strains being plotted 
    values = props_ordered(df2plot, x_var, col_var, 'Strain_o')
    x_values, col_values, shape_values = values
    
    # Count and define figure settings accordingly
    # - number of columns define legend position 
    # Define legends for columns 
    v_legend = []
    for zz, v_var, v_values in zip(count(), [x_var, col_var], [x_values, col_values]):
        if isinstance(dict_legends[v_var], dict):
            legend = []
            for v_dict in v_values: 
                legend.append(dict_legends[v_var][v_dict]['legend'])
        else: 
            legend = dict_legends[v_var]
        v_legend.append(legend)
    x_legend, col_legend = v_legend
    
    # reverse = False
    ascending = False
    if col_var == 'GenotypeAll':
        ascending = True
    elif col_var == 'Stage':
        ascending = False
        col_legend = ['Stage: ' + s for s in col_legend]
    else: 
        ascending = False
        print('HELP A')
        
    
    n_col = len(col_values)
    print(n_col)
    for nn in range(n_col): col_values.append('')
    
    # - number of variables/ stacked bars define number of columns in legend
    if n_col >= 3:
        leg_pos = 4
        loc='lower center'
        bbox_to_anchor=(0.5, 0.5)
        ncols_leg = len(vars2plot)
    else:
        leg_pos = n_col
        loc='lower center'
        if len(vars2plot) == 3:
            ncols_leg = 3
        else: 
            ncols_leg = 2
        if n_col == 1:
            bbox_to_anchor=(0.5, 0.5)
            print('A')
        else: 
            if len(colours) <= 3:
                bbox_to_anchor=(1, -1.3)
            else:
                bbox_to_anchor=(1, -1.7)
                print('B')
            
    # - number of x_var
    n_x = len(df2plot[x_var].unique())
    print(n_x)

    # If bar plots are stacked to 100 transform data to percentages
    if stack100: 
        df4plot = df4plot.apply(lambda zz: zz*100/sum(zz), axis=1)

    # Find minimum bar height to later define offset  
    min_val = df4plot.min().min()
    
    #Create figure
    gridkw = dict(width_ratios=[1]*n_col, height_ratios=[1,0.3])
    if bot_legend:
        fig_size = (math.ceil(1*(n_x/2)*(n_col)),2.4)
    else: 
        fig_size = (math.ceil(1*(n_x/2)*(n_col)),2.4) # for h1aOE subtract -0.5 from x 
    print(fig_size)
    fig, axes = plt.subplots(ncols=n_col,nrows=2, figsize=fig_size, sharey=True, gridspec_kw=gridkw)
    
    sns.set_style("ticks")
    sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    # print(rc_dict)
    
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
                ax.bar(df2plot_titles, df_col[colm], bottom=bottom, label=colm, color=colours[i], width = 0.75)
                bottom += np.array(df_col[colm])
                ax.set_xticks(df2plot_titles)
                ax.set_xticklabels(df2plot_titles, rotation=30, horizontalalignment='center', 
                                   fontname = fontname, size =8)
            
            if x_var == 'Stage':
                ax.set_xlabel('Stage [hpf]', fontname = fontname, size =8)
            y_vals = ax.get_yticks()
            y_vals_all.append(y_vals)
            max_y_vals.append(y_vals[-1])
            
            if n == 0:
                handles, labels = ax.get_legend_handles_labels()
                if not stack100 and yticks_lab == '1e3 - d.':
                    ylabel = ylabel +' x 10$^3$'
                ax.set_ylabel(ylabel, fontname = fontname, size =8)
                
            # Sum up the rows of our data to get the total value of each bar.
            totals = df_col.sum(axis=1)
            if not stack100:
                ax.set_ylim(ylim[0], ylim[1])
                # Set an offset that is used to bump the label up a bit above the bar.
                y_offset = round(0.2*min_val)
                # Add labels to each bar.
                for i, total in enumerate(totals):
                    if total < 100:
                          ax.text(i, total + y_offset, 
                                  ['{:.1f}'.format(w) for w in total], 
                                  ha='center', weight='bold', size =6, fontname = fontname)
                    else: 
                          # ax.text(i, total + y_offset ,locale.format("%d", total, grouping=True), ha='center', weight='bold', size =9)
                          if yticks_lab == '1e3 - d.':
                              ax.text(i, total + y_offset ,
                                      '{:.1f}'.format(total/1e3),
                                      ha='center', weight='bold', size =6, fontname = fontname)
                          elif yticks_lab == 'th,':
                              ax.text(i, total + y_offset ,
                                      locale.format("%d", total, grouping=True),
                                      ha='center', weight='normal', size =6, fontname = fontname)
                              
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
                          ha='center', color='w', weight='bold', size=6,fontname = fontname)
                  
            ax.set_title(col_legend[n])
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
        
        # Second row will contain the legend
        if n >= n_col:
            ax.set_axis_off()
            ax.set_title('')
            if n== leg_pos and bot_legend:
                legends = def_var_names('bar_plots', labels)
                ax.legend(handles, legends, loc=loc, bbox_to_anchor = bbox_to_anchor, 
                          frameon = True, fontsize=7, ncol = ncols_leg)
    
    # Define axes based on higherst bar
    if not stack100:
        for n, ax in zip(count(), axes.flatten()):
            max_y_index = max_y_vals.index(max(max_y_vals))
            ax.set_yticks(y_vals_all[max_y_index])
            # ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
            if yticks_lab == '1e3 - d.':
                ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]], 
                                   fontname = fontname, size =8);# print('aa')
            elif yticks_lab == 'th,':
                ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]],
                                   fontname = fontname, size =8); #print('bb')
            # ax.ticklabel_format(axis='y', style='plain')

            # ax.tick_params(axis='y', which='major', labelsize=8)
            for tick in ax.get_xticklabels():
                tick.set_fontname(fontname)
            for tick in ax.get_yticklabels():
                tick.set_fontname(fontname)
                # print('IN')
        
    else: 
        for n, ax in zip(count(), axes.flatten()):
            ax.tick_params(axis='y', which='major', labelsize=8)
            for tick in ax.get_xticklabels():
                tick.set_fontname(fontname)
            for tick in ax.get_yticklabels():
                tick.set_fontname(fontname)
                # print('IN')
        
    fig.suptitle(title+txt_title+'\n', fontsize = 8, y=1.1, fontname = fontname)
    
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'pl_bars', 'R_')
            title2save = title.replace(' per Stage and Genotype','')
            title2save = title2save.replace(' per Genotype and Stage','')
            title2save = title2save.replace(' ','_')
            if stack100: 
                title2save = 'Perc'+title2save

            fig_title = dir2savef+info+"_"+title2save+"."+extf
          
            plt.savefig(fig_title, dpi=300, bbox_inches='tight', transparent=True)

#%% func - pctChange_barPlots (test if working!)
def pctChange_barPlots(# General Plot Settings
                       df2plot, vars2plot, group_vars, x_var, col_var, 
                       # General Plot Settings
                       title, txt_title, ylabel, colours, dict_legends,
                       # Other plot settings
                       sub_bar_lab = True,
                       ctx_style = 0, bot_legend = True, 
                       # Saving settings
                       save = True, ext = ['png'], info='', dir2save = ''):
    
    contxt, font_scale, rc_dict, all_box_props, fontname  = setSNSContext(ctx_style)
    boxprops, flierprops, whiskerprops, capprops, medianprops, meanprops = all_box_props
    
    mpl.rcParams.update(mpl.rcParamsDefault)
    print('\n> '+title.replace('\n',' '))
    print('\t - x_var: '+x_var+ ' - col_var:'+col_var)
    # Create list to contain info of bars and ist highest value
    y_vals_all = []
    max_y_vals = []
    
    # Filter dataframe 
    df4plot = df2plot.groupby(group_vars)[vars2plot].mean()
    
    # Count and define figure settings accordingly
    # - number of x_var
    x_values = sorted(df2plot.index.unique(level = x_var))
    n_x = len(x_values)
    
    reverse = False
    if col_var == 'GenotypeAll':
        reverse = True
    # - number of columns define legend position 
    col_values = sorted(df2plot.index.unique(level = col_var), reverse = reverse)
    n_col = len(df2plot.index.unique(level = col_var))
    
    # Define legends for columns 
    v_legend = []
    for zz, v_var, v_values in zip(count(), [x_var, col_var], [x_values, col_values]):
        if isinstance(dict_legends[v_var], dict):
            legend = []
            for v_dict in v_values: 
                # print(v_var, v_dict)
                legend.append(dict_legends[v_var][v_dict]['legend'])
        else: 
            legend = dict_legends[v_var]
        v_legend.append(legend)
    x_legend, col_legend = v_legend

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
    
    #Create figure
    gridkw = dict(width_ratios=[1]*n_col, height_ratios=[1,0.3])
    fig_size = (math.ceil(1*(n_x/2)*(n_col)),2.4)
    fig, axes = plt.subplots(ncols=n_col,nrows=2, figsize=fig_size, sharey=True, gridspec_kw=gridkw)
    
    sns.set_style("ticks")
    sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    
    for n, ax, col_val in zip(count(), axes.flatten(), col_values):
        # First row will contain the graphs
        if n < n_col:
            df_col = df4plot.iloc[df4plot.index.get_level_values(col_var) == col_val]
            df_col = df_col.droplevel(col_var, axis="index").sort_index(key=lambda xx: xx.str.lower(), ascending=True)
            if df_col.index.nlevels > 1:
                df2plot_titles = df_col.index.map(('\n'.join))
            else: 
                # df_col_copy = df_col.copy()
                df2plot_titles = x_legend #dict_legend[x_var]
                
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
                ax.bar(df2plot_titles, df_col[colm], bottom=bottom, label=colm, color=colours[i], width = 0.75)
                bottom += np.array(df_col[colm])
                ax.set_xticks(df2plot_titles)
                ax.set_xticklabels(df2plot_titles, rotation=30, 
                                   fontname = fontname, size =8)
            y_vals = ax.get_yticks()
            y_vals_all.append(y_vals)
            max_y_vals.append(y_vals[-1])
            
            if n == 0:
                handles, labels = ax.get_legend_handles_labels()
                ax.set_ylabel(ylabel, fontname = fontname, size =8)

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
                          ha='center', color='w', weight='bold', size=6,  fontname = fontname)
            
            ax.set_title(col_legend[n],  fontname = fontname)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.axhline(0, color='dimgrey', linewidth=0.8)
        
        # Second row will contain the legend
        if n >= n_col:
            ax.set_axis_off()
            ax.set_title('')
            if n== leg_pos and bot_legend:
                legends = def_var_names('bar_plots', labels)
                ax.legend(handles, legends, loc=loc, bbox_to_anchor = bbox_to_anchor, 
                          frameon = True, fontsize=7, ncol = ncols_leg)
                # ax.legend.set_name(fontname)
        
        ax.tick_params(axis='y', which='major', labelsize=8)
        for tick in ax.get_xticklabels():
            tick.set_fontname(fontname)
        for tick in ax.get_yticklabels():
            tick.set_fontname(fontname)
        try:
            ax.yaxis.get_offset_text().set_fontname(fontname)
            # print(r' -> YAY')
        except:
            print(r' -> Nop')
        
    # # Define axes based on higherst bar
    # max_y_index = max_y_vals.index(max(max_y_vals))
    # ax.set_yticks(y_vals_all[max_y_index])
    # ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]], fontname = fontname)

    fig.suptitle(title+txt_title+'\n', fontsize = 8, y=1.2,  fontname = fontname)
    
    if save: 
        for extf in ext: 
            dir2savef = os.path.join(dir2save, 'pl_chgbars', 'R_')
            title2save = title.replace(' per Stage and Genotype','')
            title2save = title2save.replace(' per Genotype and Stage','')
            title2save = title2save.replace('\n','_')
            title2save = title2save.replace(' ','_')
            fig_title = dir2savef+info+"_"+title2save+"."+extf

            # print(fig_title)
            plt.savefig(fig_title, dpi=300, bbox_inches='tight', transparent=True)

#%% func - barPlots_oneGenot - remove?
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
            
#%% func - strainORstrain_o - remove?
# def strainORstrain_o(df2plot):
    
#     vars_opt = ['Stage','GenotypeAll']#,'Manip']
#     strain = ''
#     include_strain = ask4input('Include strain division? >>:', bool)
#     if include_strain:
#         strain_o = ask4input('Strain_o? >>:', bool)
#         if strain_o: 
#             vars_opt.append('Strain_o')
#             strain = 'Strain_o'
#         else: 
#             vars_opt.append('Strain')
#             strain = 'Strain'
    
#     group_vars = []
#     txt_title = '\n'
#     for var in vars_opt:
#         if len(df2plot[var].unique()) > 1:
#             group_vars.append(var)
#         else: 
#             txt_title = txt_title +'['+ var + ':'+df2plot[var].unique()[0]+'] '
    
#     return group_vars, txt_title, strain

#%% HEATMAPS!
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
            print('-No hmf heatmap dataframe found for '+thickness+'-'+file+' inside the results folder!')
            continue
        
    if normalise: 
        return [dfs_o, dfs_hmf], num, norm_vals
    else: 
        return [dfs_o], num
    
#%% func - unifyHeatmap
def unifyHeatmap(df, chamber, stage, strain, genotype, gen_info, thickness, vmin, vmax, n_val, normalise, dir2save, save, info, cmap = 'turbo'):
    
    sns.set_context('notebook', font_scale=1.25)#, rc = rc_dict)
    stage = stage+'hpf'
    if thickness == 'CjTh':
        title = 'Cardiac jelly thickness [$\mu$m] - ('+chamber+', '+stage+', '+strain+', '+genotype+' - n='+str(n_val)+')_'+normalise+'\n'
    elif thickness == 'myocIntBall': 
        title = 'Myocardium ballooning [$\mu$m] - ('+chamber+', '+stage+', '+strain+', '+genotype+' - n='+str(n_val)+')_'+normalise+'\n'
    elif thickness == 'MyocTh':
        title = 'Myocardial thickness [$\mu$m] - ('+chamber+', '+stage+', '+strain+', '+genotype+' - n='+str(n_val)+')_'+normalise+'\n'
    elif thickness == 'EndoTh':
        title = 'Endocardial thickness [$\mu$m] - ('+chamber+', '+stage+', '+strain+', '+genotype+' - n='+str(n_val)+')_'+normalise+'\n'
    
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
    ax.set_xticklabels(x_lab_new, rotation=30, fontname = fontname)
    
    y_pos = ax.get_yticks()
    # y_lab = ax.get_yticklabels()
    y_pos_new = np.linspace(y_pos[0], y_pos[-1], 11)
    y_lab_new = np.linspace(y_labels[0],y_labels[1],11)
    y_lab_new = [format(y,'.2f') for y in y_lab_new]
    ax.set_yticks(y_pos_new) 
    ax.set_yticklabels(y_lab_new, rotation=0, fontname = fontname)
    
    plt.ylabel('Centreline position '+y_text+'\n', fontsize=15, fontname = fontname)
    plt.xlabel('Angle (\N{DEGREE SIGN}) [Dorsal >> Right >> Ventral >> Left >> Dorsal]', fontsize=15, fontname = fontname)
    plt.title(title, fontsize = 15, fontname = fontname)
    
    dir4heatmap = os.path.join(dir2save,'pl_hmnorm', 'hmfAll_'+info+'-'+gen_info+'_'+thickness+'_'+chamber+'_'+stage+'_'+normalise+'.png')
    # print(dir4heatmap)
    if save: 
        plt.savefig(dir4heatmap, dpi=300, bbox_inches='tight', transparent=True)

#%% func - meanHM
def meanHM(df_dataset_hm, filters, groups, chamber, variable, opt_norm,  dir2load_df, dir2save_hmf, dir_data2Analyse, save, info):
    
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
    max_o = 0.7*math.ceil(max(max_vals_o))
    if normalise: 
        max_N = math.ceil(max(max_vals_N))
        
    max_o = ask4input('Define max value for heatmaps: ', int)

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
        unifyHeatmap(df_of, chamber, genotype=group[1], gen_info = gen_info, stage=group[0], strain =group[2], thickness= variable, 
                  vmin=min_o, vmax=max_o, n_val = num, normalise = 'o_all', dir2save = dir2save_hmf, save = save, info = info, cmap = 'turbo')
        if normalise:
            if perChamber:
                txt_title = 'normCh-'+norm_type
            else: 
                txt_title = 'normWh-'+norm_type
            unifyHeatmap(df_hmW, chamber, genotype=group[1], gen_info = gen_info, stage=group[0], strain =group[2], thickness= variable, 
                      vmin=0, vmax=max_N, n_val = num, normalise = txt_title,  dir2save = dir2save_hmf, 
                      save = save, info = info, cmap = 'inferno')

#%% func - concatHeatmaps
def concatHeatmaps(df2concat, operation = 'mean'):
    
    if operation == 'std':
        df_concat = pd.concat(df2concat).groupby(level=0).std()
    elif operation == 'mean':
        df_concat = pd.concat(df2concat).groupby(level=0).mean()
    elif operation == 'sem':
        df_concat = pd.concat(df2concat).groupby(level=0).sem()

    return df_concat


#%% func - defineHMRegions
def defineHMRegions(angle_div, position_div, chamber): 
    
    # Define angle region boundaries
    angles_f = [-180, -180+angle_div//2]
    while angles_f[-1]+angle_div <180:
        angles_f.append(angles_f[-1]+angle_div)
    angles_f.append(180)
    
    ang_cols = [str(ang)+'->'+str(angles_f[a+2]) for a, ang in enumerate(angles_f[1:-2])]
    ang_cols.append(str(angles_f[-2])+'->'+str(angles_f[1]))
    angles_tup = [(ang,angles_f[a+2]) for a, ang in enumerate(angles_f[1:-2])]
    angles_tup.append((angles_f[-2],angles_f[1]))
    
    # Define position region boundaries 
    if 'Atr' in chamber: 
        pos_f = np.linspace(1,2,num = int(1/position_div)+1, endpoint = True)
    else: # 'Vent'
        pos_f = np.linspace(0,1,num = int(1/position_div)+1, endpoint = True)
    
    pos_rows = [str(pos)+'->'+str(pos_f[p+1]) for p, pos in enumerate(pos_f[0:-1])]
    pos_tup = [(pos,pos_f[p+1]) for p, pos in enumerate(pos_f[0:-1])]
    
    return ang_cols, angles_tup, pos_rows, pos_tup

#%% func - defineIndex4Regions
def defineIndex4Regions(hmf, angles_tup, pos_tup): 
    
    hmf_rows = list(hmf.index)
    hmf_cols = list(hmf.columns.values)
    hmf_colsf = [float(val) for val in hmf_cols]
    # Get indexes to crop heatmap
    aa_tup = []
    for a, angle in zip(count(), angles_tup):
        # print('\n Cols', angle)
        found_ci = False; found_cii = False
        for cc, col in enumerate(hmf_colsf):
            if col >= angle[0] and not found_ci:
                ci = cc; found_ci = True
            if col >= angle[1] and not found_cii:
                cii = cc; found_cii = True
        aa_tup.append((ci,cii))
        # print('->', ci, '/', cii)
        
    pp_tup = []
    for p, pos in zip(count(), pos_tup):
        # print('\n Rows', pos)
        found_ri = False; found_rii = False
        for rr, row in enumerate(hmf_rows):
            if row >= pos[0] and not found_ri:
                ri = rr; found_ri = True
            if row >= pos[1] and not found_rii:
                rii = rr; found_rii = True
        if not found_rii:
            rii = -1
        pp_tup.append((ri,rii))
        # print('->', ri, '/', rii)
        
    return aa_tup, pp_tup

#%% func - regionaliseDF - Function based on filterUnloopedDF (fcMeshes)
def regionaliseDF(hmf, chamber, angle_div, position_div, operations =['mean','max','min','median']):
    """
    #https://stackoverflow.com/questions/38940946/average-of-multiple-dataframes-with-the-same-columns-and-indices
    Function based on filterUnloopedDF (fcMeshes)

    """
    # Define angle and position region boundaries
    ang_cols, angles_tup, pos_rows, pos_tup = defineHMRegions(angle_div, position_div, chamber)
    
    hmf_mean = pd.DataFrame(columns=['z_plane']+ang_cols)
    hmf_mean['z_plane'] = pos_rows
    hmf_max = pd.DataFrame(columns=['z_plane']+ang_cols)
    hmf_max['z_plane'] = pos_rows
    hmf_min = pd.DataFrame(columns=['z_plane']+ang_cols)
    hmf_min['z_plane'] = pos_rows
    hmf_median = pd.DataFrame(columns=['z_plane']+ang_cols)
    hmf_median['z_plane'] = pos_rows
    
    aa_tup, pp_tup = defineIndex4Regions(hmf, angles_tup, pos_tup)
    
    region_labels = []
    n_region = 0
    
    for aa, atup, ang_lab in zip(count(), aa_tup, ang_cols):
        # print(ang_lab)
        mean_col = []; max_col = []; min_col = []; median_col = []
        for pp, ptup, pos_lab in zip(count(), pp_tup, pos_rows):
            # print('\n aa:',aa,atup, ang_lab, '-pp:',pp,ptup, pos_lab)
            if aa != len(aa_tup)-1:
                # print(ang_lab)
                hmf_crop = hmf.iloc[ptup[0]:ptup[1], atup[0]:atup[1]]
            else: 
                hmf_crop_right = hmf.iloc[ptup[0]:ptup[1], atup[0]:]
                hmf_crop_left = hmf.iloc[ptup[0]:ptup[1], 0:atup[1]]
                hmf_crop = pd.concat([hmf_crop_right, hmf_crop_left], axis=1)
                
            vmean = hmf_crop.mean().mean(); mean_col.append(vmean)
            vmax = hmf_crop.max().max(); max_col.append(vmax); 
            vmin = hmf_crop.min().min(); min_col.append(vmin); 
            np_crop = hmf_crop.to_numpy().flatten()
            vmedian = np.nanmedian(np_crop); median_col.append(vmedian); 
            # vstdev = np.nanstd(np_crop)
            # print('mean:', vmean, '- max:',vmax, '- min:',vmin, '- median:', vmedian, '- stdev:', vstdev)

            # _ = plt.hist(np_crop[~np.isnan(np_crop)], bins='auto')
            # plt.title('Histogram '+str(atup)+'/'+str(ptup))
            # plt.show()
            
            if n_region < 10:
                region_labels.append((pos_lab,ang_lab,'R-0'+str(n_region)))
            else: 
                region_labels.append((pos_lab,ang_lab,'R-'+str(n_region)))
            n_region +=1

        hmf_mean[ang_lab] = mean_col
        hmf_max[ang_lab] = max_col
        hmf_min[ang_lab] = min_col
        hmf_median[ang_lab] = median_col
        
    hmf_all = []
    for op in operations: 
        if op == 'mean': 
            hmf_mean = hmf_mean.set_index('z_plane')
            hmf_mean.columns.names = ['angle_range']
            hmf_all.append(hmf_mean)
        if op == 'max': 
            hmf_max = hmf_max.set_index('z_plane')
            hmf_max.columns.names = ['angle_range']
            hmf_all.append(hmf_max)
        if op == 'min': 
            hmf_min= hmf_min.set_index('z_plane')
            hmf_min.columns.names = ['angle_range']
            hmf_all.append(hmf_min)
        if op == 'mean': 
            hmf_median = hmf_median.set_index('z_plane')
            hmf_median.columns.names = ['angle_range']
            hmf_all.append(hmf_median)

    df_reglab = pd.DataFrame(region_labels, columns=['z_plane','angle_range','hm_region'])
    
    return hmf_all, df_reglab

#%% func - stackRegHM()
def stackRegHM(angle_div, position_div,
               df_dataset_hm, df_meas, var2an, 
               dir2load_df, save, dir2save):
    
    chambers = ['Atr', 'Vent']
    variables = ['CjTh', 'myocIntBall','MyocTh', 'EndoTh']
    coVars2add = pd.DataFrame([['Vol_Atr.CJ','Vol_Atr.ExtMyoc','Vol_Atr.Myoc', 'Vol_Atr.Endo'], 
                               ['Vol_Vent.CJ','Vol_Vent.ExtMyoc','Vol_Vent.Myoc', 'Vol_Vent.Endo']],
                                    index=chambers, columns=variables)
    df_m = df_meas.set_index('Folder')
    
    ini = False
    print('\n>>> Dividing '+ var2an+' heatmaps into regions...')
    folders = df_dataset_hm['Folder']#filterR_Autom (df_dataset_hm, filters, group, col_out = 'Folder')
    df_data = df_dataset_hm.set_index('Folder')
    for chamber in chambers: 
        print('\n >> Variable: ', var2an,' - Chamber: ', chamber, '\n')
        hmf_file = 'hmf_'+'unloop'+chamber+'_'+var2an
        bar = Bar('- Processing ', max=len(folders), suffix = suffix, check_tty=False, hide_cursor=False)
        for n, file in enumerate(folders):
            try: 
                hmf = loadDF(file[2:], hmf_file, dir2load_df)
                # print(r'Filename: ',file, ' - heatmap shape: ', hmf.shape)
                hmf_all, df_reglab = regionaliseDF(hmf, chamber, angle_div, position_div, operations =['mean'])#,'max','min','median'])
                hmf_mean = hmf_all[0] # hmf_mean, hmf_max, hmf_min, hmf_median = hmf_all
                df_reglab = df_reglab.set_index(['z_plane', 'angle_range'])
                reglab = list(df_reglab['hm_region'])
                # input()
                st_mean = hmf_mean.stack()
                st_mean = st_mean.to_frame()
                st_mean = st_mean.rename_axis(['z_plane', 'angle_range'])
                st_mean = st_mean.rename(columns={0: 'value'})
                st_mean = st_mean.join(df_reglab, on=['z_plane', 'angle_range'])
                st_mean_reg = st_mean.set_index('hm_region')
                st_mean_regT = st_mean_reg.T
                
                st_mean_regT['Chamber'] = chamber
                st_mean_regT['var2An'] = var2an
                st_mean_regT['Folder'] = file
                st_mean_regT['Strain_o'] = df_data.loc[file,'Strain_o']
                st_mean_regT['Stage'] = df_data.loc[file,'Stage']
                st_mean_regT['Manip'] = df_data.loc[file,'Manip']
                st_mean_regT['GenotypeAll'] = df_data.loc[file,'GenotypeAll']
                st_mean_regT['GenotypeF'] = df_data.loc[file,'GenotypeF']
                coVar2add = coVars2add.loc[chamber, var2an]
                st_mean_regT['coVar2add'] = coVar2add
                st_mean_regT['coVar_value'] = df_m.loc[file[2:]+'_2A',coVar2add]
                
                # print(st_mean_regT)
                if not ini: 
                    st_meanF = st_mean_regT.copy()
                    ini = True
                else: 
                    st_meanF = pd.concat([st_meanF, st_mean_regT], ignore_index=True)
            except: 
                print('\n-No hmf heatmap dataframe found for '+var2an+'-'+file+' inside the results folder!')
                
            bar.next()
        bar.finish()
        
    first_cols = ['Folder', 'Strain_o', 'Stage', 'Manip','GenotypeAll', 'GenotypeF',
                  'Chamber','var2An','coVar2add', 'coVar_value']
    
    st_meanF = st_meanF.reindex(first_cols+reglab, axis =1)
    
    if save: 
        saveDF(filename = 'hm_reg', df2save = st_meanF, df_name = var2an, dir2save = dir2save)
        
    alert('frog', 1)

    return st_meanF, reglab

#%% KDE PLOTS
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
                ax.set(xlabel='Cardiac jelly thickness [$\mu$m]', ylabel='\nThickness distribution\n'+ cl_lab + '\n')
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
        title = filename + ' - Cardiac Jelly Thickness [$\mu$m]'
        xlabel = 'Cardiac Jelly Thickness [$\mu$m]'
        step = 0.05
    elif variable == 'myoc_intBall' :
        title = filename + ' - Myoc.Int Ballooning [$\mu$m]'
        xlabel = 'Myoc.Int Ballooning [$\mu$m]'
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
                ax.set(xlabel='Cardiac jelly thickness [$\mu$m]', ylabel='\nThickness distribution\n'+ cl_lab + '\n')
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


#%% ------
#%% OLD - But Working!
# #%% func - plotStatsperXBU (working, text of tests not included)
# def plotStatsperXBU (# General Plot Settings
#                     df2plot, vars2plot, x_var, hue_var, shape_var, 
#                     # Size, title, labels and legends
#                     title, labels2plot, dict_legends, ips, 
#                     # Statistic Settings
#                     stats_set,
#                     # Other settings
#                     suptitle = True, right_legend = False,
#                     yticks_lab = 'th,', ylim = '', box_plot = True, show_fliers = False,
#                     # Saving settings
#                     save = True, dpi = 300, ext = 'png', info = '', dir2save = ''):

#     if stats_set[0]: 
#         tests_res = []
#         box_pairs_all = []
#         statist = True
#         dicts_stats = stats_set[1]
#         alpha = stats_set[2]
#         txt_multcomp_all = []
        
#     print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    
#     # Genotypes and Strains being plotted 
#     values = props_ordered(df2plot, x_var, hue_var, shape_var)
#     x_values, hue_values, shape_values = values
    
#     # Set up the matplotlib figure 
#     n_rows = 1
#     n_cols = len(hue_values)

#     for index in range(n_cols):
#         if index == 0 and len(vars2plot) == 1:
#             addVars = True
#             var2copy = vars2plot[0]
#             label2copy = labels2plot[0]
        
#         if addVars and index != 0:
#             vars2plot.insert(index, var2copy)
#             labels2plot.insert(index, label2copy)
    
#     # As a right legend wants to be added, we need to add '' to all labels 
#     if right_legend: 
#         index_no_plot = [n_cols]
#         for index in index_no_plot:
#             vars2plot.insert(index, '')
#             labels2plot.insert(index, '')
#         width_ratios = [1]*n_cols+[0.2]
#         n_cols = n_cols+1
#     else: 
#         width_ratios = [1]*n_cols
#         index_no_plot = [1000]
        
#     # Set up the matplotlib figure
#     h_plot, w_plot = ips
#     h_add = 1; w_add = 1

#     if len(x_values) > 3:
#         h_add = 0
#         wspace = 0.05
#         box_width = 0.8
#         bbox_lf = -0.5
#     else:
#         h_plot = 0.8
#         wspace = 0.15
#         box_width = 0.6
#         bbox_lf = -1.8
    
#     size_col = (n_cols)*(h_plot*len(x_values))+h_add
#     size_row = n_rows*w_plot+w_add
    
#     # Define legends for x and hue
#     v_legend = []; v_color = []
#     for zz, v_var, v_values in zip(count(), [x_var, hue_var], [x_values, hue_values]):
#         if isinstance(dict_legends[v_var], dict):
#             legend = []
#             color = []
#             for v_dict in v_values: 
#                 legend.append(dict_legends[v_var][v_dict]['legend'])
#                 color.append(dict_legends[v_var][v_dict]['color'])
#         else: 
#             legend = dict_legends[v_var]
#             color = ''
#         v_legend.append(legend)
#         v_color.append(color)
#     x_legend, hue_legend = v_legend
#     x_color, _ = v_color
        
#     if right_legend:
#         # Define legends for x
#         legend_elem = []
#         for aa, xval, xcol in zip(count(), x_legend, x_color):
#             legend_elem.append(Line2D([0], [0], marker='o', color='w', label=xval,
#                                     markerfacecolor=xcol, markersize=12))
#         handle_new = legend_elem
        
#         for index in index_no_plot:
#             hue_legend.insert(index, '')
#             hue_values.insert(index, '')
            
#     marker_size = 8; dodge = True; jitter = 0.3
    
#     for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
#         print('\t- '+svar+': ', value)
        
#     ##  CREATE FIGURE 
#     gridkw = dict(width_ratios=width_ratios)
#     fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
#     fig.subplots_adjust(hspace=1, wspace=wspace)
    
#     sns.set_style("ticks")
#     sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
    
#     # Create list to contain info of yticks and its highest value
#     y_vals_all = []; max_y_vals = []
#     tests_res_all = []; box_pairs_allf = []
#     for n, ax, var, ylabel, hue_value, dict_stats in zip(count(), axes.flatten(), vars2plot, labels2plot, hue_values, dicts_stats):
        
#         if n in index_no_plot and right_legend:
#                 if n == n_cols-1:
#                     ax.set_axis_off()
#                     ax.legend(handle_new, x_legend, loc='upper left', bbox_to_anchor=(bbox_lf, 1), frameon = False)
#                 else: 
#                     ax.remove()
#         else: 
#             df_xfilt = df2plot[df2plot[hue_var] == hue_value]
#             m = sns.stripplot(data=df_xfilt, x=x_var, y=var, ax = ax, order=x_values,
#                               marker = 'o', palette = x_color, size = marker_size, 
#                               linewidth=0.5, jitter = jitter)
#             if box_plot: 
#                 m = sns.boxplot(data=df_xfilt, x=x_var, y=var, ax = ax, order=x_values,
#                                    dodge = dodge, width= box_width, showfliers = show_fliers, palette = x_color,
#                                    boxprops = boxprops, whiskerprops = whiskerprops, capprops = capprops, 
#                                    flierprops = flierprops, medianprops = medianprops, 
#                                    meanline = True, meanprops = meanprops, showmeans = True)
    
#             box = ax.get_position()
            
#             if statist: 
#                 box_pairs = dict_stats[var]['box_pairs']
#                 stat_test = False
#                 test = None
#                 p_val = dict_stats[var]['pval_multComp_all']
#                 txt_testSelected_all = dict_stats[var]['txt_testSelected_all']
#                 txt_normtest_all = dict_stats[var]['txt_normtest_all']
                
#                 # Added fontname as an imput to add_stat_annotation function in statannot
#                 m1, test_results = add_stat_annotation(ax, data=df2plot, x=x_var, y=var, hue=hue_var, 
#                                                        hue_order = hue_values, order = x_values,
#                                                        box_pairs=box_pairs, perform_stat_test = stat_test, test = test,
#                                                        pvalues=p_val, comparisons_correction=None, #'bonferroni',
#                                                        line_offset_to_box=0.4, line_offset=0.1,
#                                                        line_height=0.015, text_offset=5,
#                                                        text_format='star', loc='inside', fontsize='small', verbose=0,
#                                                        fontname = fontname)
#                 tests_res.append(p_val)
#                 box_pairs_all.append(box_pairs)
#                 txt_multcomp = txtMultComp(box_pairs, p_val, hue_values[n], hue_legend[n], txt_testSelected_all, 
#                                             txt_normtest_all, alpha, indiv = True)
#                 txt_multcomp_all.append(txt_multcomp)
                
#             if n == 0:
#                 if yticks_lab == '1e6 - d.':
#                     ylabel = ylabel +' x 10$^6$'
#                 elif yticks_lab == '1e3 - d.':
#                     ylabel = ylabel +' x 10$^3$'
#                 ax.set_xlabel(hue_legend[n], fontname = fontname)
#                 ax.set_ylabel(ylabel, fontname = fontname)
#             else:
#                 ax.set_xlabel(hue_legend[n], fontname = fontname)
#                 ax.set_ylabel('', fontname = fontname)
    
#             ax.set_position([box.x0, box.y0, box.width*1, box.height])
#             xticks = ax.get_xticks()
#             ax.set_xticks(xticks)
#             ax.set_xticklabels(x_legend, rotation=45, horizontalalignment='right', fontname = fontname)
#             sns.despine()
            
#             if ylim != '':
#                 print('ylim:', ylim[n][0],'-', ylim[n][1])
#                 ax.set_ylim(ylim[n][0], ylim[n][1])
            
#             if n == 0:
#                 handles, labels = m.get_legend_handles_labels()
#             else:
#                 ax.spines['left'].set_linestyle('-.')
#                 ax.spines['left'].set_color('#696969')
#                 ax.spines['left'].set_linewidth(0.8)
#                 ax.tick_params(left = False)
    
#             y_vals = ax.get_yticks()
#             y_vals_all.append(y_vals)
#             max_y_vals.append(y_vals[-1])
                
#             # Define axes based on higherst bar
#             if yticks_lab == '1e6 - d.':
#                 ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals])
#             elif yticks_lab == '1e3 - d.':
#                 ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals])
#             elif yticks_lab == 'th,':
#                 ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals])
#             elif yticks_lab == 'd. - 0':
#                 ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals])
#             elif yticks_lab == 'd. - 1':
#                 ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals])
#             elif yticks_lab == 'd.':
#                 ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals])
            
#             for tick in ax.get_xticklabels():
#                 tick.set_fontname(fontname)
#             for tick in ax.get_yticklabels():
#                 tick.set_fontname(fontname)
#             try:
#                 ax.yaxis.get_offset_text().set_fontname(fontname)
#                 # print(r' -> YAY')
#             except:
#                 print(r' -> Nop')

#     if suptitle:
#         fig.suptitle(title+'\n', y=1, fontname = fontname)
#     # fig.tight_layout()
    
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'pl_meas', 'R_')
#             if info != '':
#                 fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
#             else: 
#                 fig_title = dir2savef+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf

#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
            
#     if statist: 
#         return tests_res, box_pairs_all
            
# #%% func - plotNoStatsperXBU (working fine, no legend included)
# def plotNoStatsperXBU (# General Plot Settings
#                     df2plot, vars2plot, x_var, hue_var, shape_var, 
#                     # Size, title, labels and legends
#                     title, labels2plot, dict_legends, ips, suptitle = True,
#                     # Other settings
#                     yticks_lab = 'th,', ylim = '', box_plot = True, show_fliers = False,
#                     # Saving settings
#                     save = True, dpi = 300, ext = 'png', info = '', dir2save = ''):

    
#     print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
    
#     # Genotypes and Strains being plotted 
#     values = props_ordered(df2plot, x_var, hue_var, shape_var)
#     x_values, hue_values, shape_values = values
    
#     # Set up the matplotlib figure 
#     n_rows = 1
#     n_cols = len(hue_values)

#     for index in range(n_cols):
#         if index == 0 and len(vars2plot) == 1:
#             addVars = True
#             var2copy = vars2plot[0]
#             label2copy = labels2plot[0]
        
#         if addVars and index != 0:
#             vars2plot.insert(index, var2copy)
#             labels2plot.insert(index, label2copy)
    
#     # Set up the matplotlib figure
#     h_plot, w_plot = ips
#     h_add = 1; w_add = 1

#     if len(x_values) > 3:
#         h_add = 0
#         wspace = 0.05
#         box_width = 0.8
#     else:
#         h_plot = 0.8
#         wspace = 0.15
#         box_width = 0.6
    
#     size_col = (n_cols)*(h_plot*len(x_values))+h_add
#     size_row = n_rows*w_plot+w_add
    
#     # Define legends for x and hue
#     v_legend = []; v_color = []
#     for zz, v_var, v_values in zip(count(), [x_var, hue_var], [x_values, hue_values]):
#         if isinstance(dict_legends[v_var], dict):
#             legend = []
#             color = []
#             for v_dict in v_values: 
#                 print(v_var, v_dict)#dict_legends[x_var].keys():
#                 legend.append(dict_legends[v_var][v_dict]['legend'])
#                 color.append(dict_legends[v_var][v_dict]['color'])
#         else: 
#             legend = dict_legends[v_var]
#             color = ''
            
#         v_legend.append(legend)
#         v_color.append(color)
    
#     x_legend, hue_legend = v_legend
#     x_color, hue_color = v_color
        
#     # hue_legend = dict_legends[hue_var]
    
#     marker_size = 8; dodge = True; jitter = 0.3
    
#     for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
#         print('\t- '+svar+': ', value)
        
#     ##  CREATE FIGURE 
#     gridkw = dict(width_ratios=[1]*n_cols)
#     fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
#     fig.subplots_adjust(hspace=1, wspace=wspace)
    
#     sns.set_style("ticks")
#     sns.set_context(contxt, font_scale = font_scale, rc = rc_dict)
#     # print(sns.plotting_context())
    
#     y_vals_all = []; max_y_vals = []
#     for n, ax, var, ylabel, hue_value in zip(count(), axes.flatten(), vars2plot, labels2plot, hue_values):
#         # Create list to contain info of yticks and its highest value
#         print(' > Plotting ',hue_value,'\n')
#         df_xfilt = df2plot[df2plot[hue_var] == hue_value]
#         m = sns.stripplot(data=df_xfilt, x=x_var, y=var, ax = ax, order=x_values,
#                           marker = 'o', palette = x_color, size = marker_size, 
#                           linewidth=0.5, jitter = jitter)
#         if box_plot: 
#             m = sns.boxplot(data=df_xfilt, x=x_var, y=var, ax = ax, order=x_values,
#                                dodge = dodge, width= box_width, showfliers = show_fliers, palette = x_color,
#                                boxprops = boxprops, whiskerprops = whiskerprops, capprops = capprops, 
#                                flierprops = flierprops, medianprops = medianprops, 
#                                meanline = True, meanprops = meanprops, showmeans = True)

#         box = ax.get_position()

#         if n == 0:
#             if yticks_lab == '1e6 - d.':
#                 ylabel = ylabel +' x 10$^6$'
#             elif yticks_lab == '1e3 - d.':
#                 ylabel = ylabel +' x 10$^3$'
#             ax.set_xlabel(hue_legend[n], fontname = fontname)
#             ax.set_ylabel(ylabel, fontname = fontname)
#         else:
#             ax.set_xlabel(hue_legend[n], fontname = fontname)
#             ax.set_ylabel('', fontname = fontname)

#         ax.set_position([box.x0, box.y0, box.width*1, box.height])
#         xticks = ax.get_xticks()
#         ax.set_xticks(xticks)
#         ax.set_xticklabels(x_legend, rotation=45, horizontalalignment='right', fontname = fontname)#, fontsize = 16)
#         sns.despine()
        
#         if ylim != '':
#             print('ylim:', ylim[n][0],'-', ylim[n][1])
#             ax.set_ylim(ylim[n][0], ylim[n][1])
        
#         if n == 0:
#             handles, labels = m.get_legend_handles_labels()
#         else:
#             ax.spines['left'].set_linestyle('-.')
#             ax.spines['left'].set_color('#696969')
#             ax.spines['left'].set_linewidth(0.8)
#             ax.tick_params(left = False)

#         y_vals = ax.get_yticks()
#         y_vals_all.append(y_vals)
#         max_y_vals.append(y_vals[-1])
            
#         # Define axes based on higherst bar
#         if yticks_lab == '1e6 - d.':
#             ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals])
#         elif yticks_lab == '1e3 - d.':
#             ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals])
#         elif yticks_lab == 'th,':
#             ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals])
#         elif yticks_lab == 'd. - 0':
#             ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals])
#         elif yticks_lab == 'd. - 1':
#             ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals])
#         elif yticks_lab == 'd.':
#             ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals])
        
#         for tick in ax.get_xticklabels():
#             tick.set_fontname(fontname)
#         for tick in ax.get_yticklabels():
#             tick.set_fontname(fontname)
#         try:
#             ax.yaxis.get_offset_text().set_fontname(fontname)
#             print(r' -> YAY')
#         except:
#             print(r' -> Nop')

#     if suptitle:
#         fig.suptitle(title+'\n', y=1, fontname = fontname)
#     # fig.tight_layout()
    
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'pl_meas', 'R_')
#             if info != '':
#                 fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
#             else: 
#                 fig_title = dir2savef+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf

#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)

#%% - ALERT WHEN IMPORTED
print ("IMPORTED: morphoHeart_funcAnalysis")
alert('jump',1)

#%% - OTHERS TO ORGANIZE?
#%% Plots variables x3 (stages)

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

#%% END
#%% TO REMOVE
#%% func - plotInGroups2
# def plotInGroups2 (df2plot, vars2plot, x_var, hue_var, shape_var, title, labels2plot, ips, dir2save,
#                   n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th,', info ='', 
#                   save = True, dpi = 300, ext = 'png'):
    
#     print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
#     sns.set_context('poster') # notebook, talk, poster, paper
#     # Get legends
#     dict_legends = def_legends(df2plot)
    
#     # Set up the matplotlib figure
#     num_vars = len(vars2plot)
#     n_rows = math.ceil(num_vars/n_cols)
#     # print('n_rows:' , n_rows)
    
#     index_right_col = list(range(n_cols,(n_cols+1)*n_rows,4))
#     index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
#     index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))

#     for index in index_no_plot:
#         vars2plot.insert(index, '')
#         labels2plot.insert(index, '')
    
#     # Set up the matplotlib figure
#     h_plot, w_plot = ips
#     if num_vars == 1:
#         h_add = 0; w_add = 0
#     size_col = (n_cols+1)*h_plot+h_add
#     size_row = n_rows*w_plot+w_add
    
#     # Genotypes and Strains being plotted 
#     values = []
#     for var in [x_var, hue_var, shape_var]:
#         if var =='GenotypeAll':
#             reverse = True
#         else: 
#             reverse = False
#         values.append(sorted(df2plot[var].unique(), reverse=reverse))
#     x_values, hue_values, shape_values = values
    
#     # - number of x_var
#     n_x = len(x_values)
    
#     #  Create figure  - plt.clf()
#     gridkw = dict(width_ratios=[1]*n_cols+[0.2])
#     fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
#     fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
#     sns.set_style("ticks")
#     sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left": True,
#                                                     "ytick.labelsize": 8, "lines.linewidth": 2.5,
#                                                     "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize": "large"})
    
#     # Define legends for shape and hue
#     hue_legend = dict_legends[hue_var]
#     legend_elem_hue = []
#     for aa, hue_val in enumerate(hue_legend):
#         legend_elem_hue.append(Line2D([0], [0], marker='h', color='w', label=hue_val,
#                                 markerfacecolor=palettes[aa], markersize=20))
#     space = [Line2D([0], [0], marker='o', color='w', label='',
#                                 markerfacecolor='w', markersize=20)]
#     shape_legend = dict_legends[shape_var]
#     legend_elem_shape = []
#     for n_str, shape_val, mark in zip(count(), shape_legend, styles):
#         legend_elem_shape.append(Line2D([0], [0], marker=mark, color='w', label=shape_val,
#                                 markerfacecolor='k', markersize=15))
        
#     handle_new = legend_elem_hue+space+legend_elem_shape
#     legend_new = hue_legend+['']+shape_legend
    
#     marker_size = 12
#     dodge = True
#     for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
#         # Create list to contain info of yticks and ist highest value
#         y_vals_all = []
#         max_y_vals = []
    
#         if n == 0: 
#             for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
#                 print('\t- '+svar+': ', value)
                
#         if n in index_no_plot:
#             if n == n_cols:
#                 ax.set_axis_off()
#                 ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
#             else: 
#                 ax.remove()
#         else: 
#             m = sns.swarmplot(data=df2plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
#                               palette = palettes, dodge = dodge, size = marker_size)
#             box = ax.get_position()
#             if yticks_lab == '1e6 - d.':
#                 ylabel = ylabel +' x 10$^6$'
#             elif yticks_lab == '1e3 - d.':
#                 ylabel = ylabel +' x 10$^3$'
#             ax.set(xlabel=dict_legends['xlabels'][x_var], ylabel=ylabel);
#             # ax.set_xticks(ax.get_xticks())
#             ax.set_xticklabels(dict_legends[x_var], rotation=0)
#             ax.set_position([box.x0, box.y0, box.width*1, box.height])
#             ax.get_legend().remove()
#             sns.despine()
            
#             if n == 0:
#                 handles, labels = m.get_legend_handles_labels()

#             y_vals = ax.get_yticks()
#             y_vals_all.append(y_vals)
#             max_y_vals.append(y_vals[-1])
#             if n_x > 1: 
#                 vline_pos = list(range(1,n_x,1))
#                 for pos in vline_pos:
#                     ax.axvline(pos - 0.5, ymax = 0.95, color='dimgrey', ls='-.', linewidth=0.8)
            
#             # Define axes based on higherst bar
#             max_y_index = max_y_vals.index(max(max_y_vals))
#             ax.set_yticks(y_vals_all[max_y_index])
#             if yticks_lab == '1e6 - d.':
#                 ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == '1e3 - d.':
#                 ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'th,':
#                 ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 0':
#                 ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 1':
#                 ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd.':
#                 ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])

#     fig.suptitle(title+'\n', fontsize = 30, y=1)
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'pl_meas', 'R_')
#             if info != '':
#                 fig_title = dir2savef+info+"_"+title+"."+extf
#             else: 
#                 fig_title = dir2savef+title+"."+extf

#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)


#%% func - plotInGroupsShape
# https://towardsdatascience.com/matplotlib-vs-ggplot2-c86dd35a9378
# https://www.r-graph-gallery.com/boxplot.html
# https://medium.com/analytics-vidhya/r-style-visualizations-in-python-560c6bbfb14a
# https://matplotlib.org/2.0.2/examples/pylab_examples/spine_placement_demo.html
# https://datascienceplus.com/seaborn-categorical-plots-in-python/
# https://towardsdatascience.com/scattered-boxplots-graphing-experimental-results-with-matplotlib-seaborn-and-pandas-81f9fa8a1801
# https://github.com/cfcooney/medium_posts/blob/master/scattered_boxplots.ipynb
# https://htmlcolorcodes.com/

# def plotInGroupsShape (df2plot, vars2plot, x_var, hue_var, shape_var, title, labels2plot, ips, dir2save,
#                   n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th,', ylim = '', info ='', 
#                   save = True, dpi = 300, ext = 'png'):
    
#     print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
#     sns.set_context('poster') # notebook, talk, poster, paper
#     # Get legends
#     dict_legends = def_legends(df2plot)
    
#     # Set up the matplotlib figure
#     num_vars = len(vars2plot)
#     n_rows = math.ceil(num_vars/n_cols)
#     # print('n_rows:' , n_rows)
    
#     index_right_col = list(range(n_cols,(n_cols+1)*n_rows,4))
#     index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
#     index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))

#     for index in index_no_plot:
#         vars2plot.insert(index, '')
#         labels2plot.insert(index, '')
    
#     # Set up the matplotlib figure
#     h_plot, w_plot = ips
    
#     # Genotypes and Strains being plotted 
#     values = []
#     for var in [x_var, hue_var, shape_var]:
#         if var =='GenotypeAll':
#             reverse = True
#         else: 
#             reverse = False
#         values.append(sorted(df2plot[var].unique(), reverse=reverse))
#     x_values, hue_values, shape_values = values
    
#     # - number of x_var
#     n_x = len(x_values)
#     if n_x == 1:
#         h_plot = 3.5
#     if num_vars == 1:
#         h_add = 0; w_add = 0
#     size_col = (n_cols+1)*h_plot+h_add
#     size_row = n_rows*w_plot+w_add
    
#     #  Create figure  - plt.clf()
#     gridkw = dict(width_ratios=[1]*n_cols+[0.2])
#     fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
#     fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
#     sns.set_style("ticks")
#     sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
#                                                     "ytick.labelsize":8, "lines.linewidth": 2.5,
#                                                     "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
    
#     # Define legends for shape and hue
#     hue_legend = dict_legends[hue_var]
#     legend_elem_hue = []
#     for aa, hue_val in enumerate(hue_legend):
#         legend_elem_hue.append(Line2D([0], [0], marker='h', color='w', label=hue_val,
#                                 markerfacecolor=palettes[aa], markersize=20))
#     space = [Line2D([0], [0], marker='o', color='w', label='',
#                                 markerfacecolor='w', markersize=20)]
#     shape_legend = dict_legends[shape_var]
#     legend_elem_shape = []
#     for n_str, shape_val, mark in zip(count(), shape_legend, styles):
#         legend_elem_shape.append(Line2D([0], [0], marker=mark, color='w', label=shape_val,
#                                 markerfacecolor='k', markersize=15))
        
#     handle_new = legend_elem_hue+space+legend_elem_shape
#     legend_new = hue_legend+['']+shape_legend
    
#     marker_size = 10
#     dodge = True
#     jitter = True#0.3
    
#     plot_no = 0
#     for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
#         # Create list to contain info of yticks and its highest value
#         y_vals_all = []
#         max_y_vals = []
#         if n == 0: 
#             for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
#                 print('\t- '+svar+': ', value)
                
#         if n in index_no_plot:
#             if n == n_cols:
#                 ax.set_axis_off()
#                 ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
#             else: 
#                 ax.remove()
#         else: 
#             for j, val, style in zip(count(), shape_values, styles):
#                 df_plot = df2plot[df2plot[shape_var] == val]
#                 # print(x_var, var, hue_var, hue_values, x_values, style, palettes)
                
#                 m = sns.boxplot(data=df_plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
#                                   dodge = dodge, width=1, boxprops = boxprops, whiskerprops = whiskerprops, capprops = capprops,
#                                   flierprops = flierprops, medianprops = medianprops, showmeans = False)
#                 m = sns.stripplot(data=df_plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
#                                   marker = style, palette = palettes, jitter=jitter, dodge = dodge, size = marker_size, 
#                                   linewidth=1)
                
#                 for pp,boxs in enumerate(m.artists):
#                 #     box.set_edgecolor('black')
#                     boxs.set_facecolor('white')

#                 box = ax.get_position()
#                 if j == 0:
#                     if yticks_lab == '1e6 - d.':
#                         ylabel = ylabel +' x 10$^6$'
#                     elif yticks_lab == '1e3 - d.':
#                         ylabel = ylabel +' x 10$^3$'
#                 ax.set(xlabel=dict_legends['xlabels'][x_var], ylabel=ylabel);
#                 # ax.set_xticks(ax.get_xticks())
#                 ax.set_xticklabels(dict_legends[x_var], rotation=0)
#                 ax.set_position([box.x0, box.y0, box.width*1, box.height])
#                 ax.get_legend().remove()
#                 sns.despine()
                
#                 if ylim != '':
#                     # print('ylim:', ylim[plot_no][0],'-', ylim[plot_no][1])
#                     ax.set_ylim(ylim[plot_no][0], ylim[plot_no][1])
                
#                 if n == 0:
#                     handles, labels = m.get_legend_handles_labels()

#                 y_vals = ax.get_yticks()
#                 y_vals_all.append(y_vals)
#                 max_y_vals.append(y_vals[-1])
#                 if n_x > 1: 
#                     vline_pos = list(range(1,n_x,1))
#                     for pos in vline_pos:
#                         ax.axvline(pos - 0.5, ymax = 0.95, color='dimgrey', ls='-.', linewidth=0.8)
            
#             # Define axes based on higherst bar
#             max_y_index = max_y_vals.index(max(max_y_vals))
#             ax.set_yticks(y_vals_all[max_y_index])
#             if yticks_lab == '1e6 - d.':
#                 ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == '1e3 - d.':
#                 ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'th,':
#                 ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 0':
#                 ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 1':
#                 ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd.':
#                 ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])
                
#         plot_no +=1

#     fig.suptitle(title+'\n', fontsize = 30, y=1)
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'pl_meas', 'R_')
#             if info != '':
#                 fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
#             else: 
#                 fig_title = dir2savef+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf

#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)

#%% func - plotInGroupsShapeBU
# def plotInGroupsShapeBU (df2plot, vars2plot, x_var, hue_var, shape_var, title, labels2plot, ips, dir2save,
#                   n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th,', ylim = '', info ='', 
#                   save = True, dpi = 300, ext = 'png'):
    
#     print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
#     sns.set_context('poster') # notebook, talk, poster, paper
#     # Get legends
#     dict_legends = def_legends(df2plot)
    
#     # Set up the matplotlib figure
#     num_vars = len(vars2plot)
#     n_rows = math.ceil(num_vars/n_cols)
#     # print('n_rows:' , n_rows)
    
#     index_right_col = list(range(n_cols,(n_cols+1)*n_rows,4))
#     index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
#     index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))

#     for index in index_no_plot:
#         vars2plot.insert(index, '')
#         labels2plot.insert(index, '')
    
#     # Set up the matplotlib figure
#     h_plot, w_plot = ips
    
#     # Genotypes and Strains being plotted 
#     values = []
#     for var in [x_var, hue_var, shape_var]:
#         if var =='GenotypeAll':
#             reverse = True
#         else: 
#             reverse = False
#         values.append(sorted(df2plot[var].unique(), reverse=reverse))
#     x_values, hue_values, shape_values = values
    
#     # - number of x_var
#     n_x = len(x_values)
#     if n_x == 1:
#         h_plot = 3.5
#     if num_vars == 1:
#         h_add = 0; w_add = 0
#     size_col = (n_cols+1)*h_plot+h_add
#     size_row = n_rows*w_plot+w_add
    
#     #  Create figure  - plt.clf()
#     gridkw = dict(width_ratios=[1]*n_cols+[0.2])
#     fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
#     fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
#     sns.set_style("ticks")
#     sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
#                                                     "ytick.labelsize":8, "lines.linewidth": 2.5,
#                                                     "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
    
#     # Define legends for shape and hue
#     hue_legend = dict_legends[hue_var]
#     legend_elem_hue = []
#     for aa, hue_val in enumerate(hue_legend):
#         legend_elem_hue.append(Line2D([0], [0], marker='h', color='w', label=hue_val,
#                                 markerfacecolor=palettes[aa], markersize=20))
#     space = [Line2D([0], [0], marker='o', color='w', label='',
#                                 markerfacecolor='w', markersize=20)]
#     shape_legend = dict_legends[shape_var]
#     legend_elem_shape = []
#     for n_str, shape_val, mark in zip(count(), shape_legend, styles):
#         legend_elem_shape.append(Line2D([0], [0], marker=mark, color='w', label=shape_val,
#                                 markerfacecolor='k', markersize=15))
        
#     handle_new = legend_elem_hue+space+legend_elem_shape
#     legend_new = hue_legend+['']+shape_legend
    
#     marker_size = 10
#     dodge = True
#     jitter = 0.3
    
#     plot_no = 0
#     for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
#         # Create list to contain info of yticks and ist highest value
#         y_vals_all = []
#         max_y_vals = []
    
#         if n == 0: 
#             for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
#                 print('\t- '+svar+': ', value)
                
#         if n in index_no_plot:
#             if n == n_cols:
#                 ax.set_axis_off()
#                 ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
#             else: 
#                 ax.remove()
#         else: 
#             for j, val, style in zip(count(), shape_values, styles):
#                 df_plot = df2plot[df2plot[shape_var] == val]
#                 # print(x_var, var, hue_var, hue_values, x_values, style, palettes)
#                 m = sns.stripplot(data=df_plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
#                                   marker = style, palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
#                 box = ax.get_position()
#                 if yticks_lab == '1e6 - d.':
#                     ylabel = ylabel +' x 10$^6$'
#                 elif yticks_lab == '1e3 - d.':
#                     ylabel = ylabel +' x 10$^3$'
#                 ax.set(xlabel=dict_legends['xlabels'][x_var], ylabel=ylabel);
#                 # ax.set_xticks(ax.get_xticks())
#                 ax.set_xticklabels(dict_legends[x_var], rotation=0)
#                 ax.set_position([box.x0, box.y0, box.width*1, box.height])
#                 ax.get_legend().remove()
#                 sns.despine()
                
#                 if ylim != '':
#                     # print('ylim:', ylim[plot_no][0],'-', ylim[plot_no][1])
#                     ax.set_ylim(ylim[plot_no][0], ylim[plot_no][1])
                
#                 if n == 0:
#                     handles, labels = m.get_legend_handles_labels()

#                 y_vals = ax.get_yticks()
#                 y_vals_all.append(y_vals)
#                 max_y_vals.append(y_vals[-1])
#                 if n_x > 1: 
#                     vline_pos = list(range(1,n_x,1))
#                     for pos in vline_pos:
#                         ax.axvline(pos - 0.5, ymax = 0.95, color='dimgrey', ls='-.', linewidth=0.8)
            
#             # Define axes based on higherst bar
#             max_y_index = max_y_vals.index(max(max_y_vals))
#             ax.set_yticks(y_vals_all[max_y_index])
#             if yticks_lab == '1e6 - d.':
#                 ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == '1e3 - d.':
#                 ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'th,':
#                 ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 0':
#                 ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 1':
#                 ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd.':
#                 ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])
                
#         plot_no +=1

#     fig.suptitle(title+'\n', fontsize = 30, y=1)
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'pl_meas', 'R_')
#             if info != '':
#                 fig_title = dir2savef+info+"_"+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf
#             else: 
#                 fig_title = dir2savef+title+"_(x_var_"+x_var+"-hue_var_"+hue_var+")."+extf

#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
            

#%% func - plotInGroupsStatsBU
# def plotInGroupsStatsBU(df2plot, vars2plot, x_var, hue_var, shape_var, title, labels2plot, ips, dir2save, stats_set,
#                       n_cols = 3, h_add = 5, w_add = 1, sharey = False, yticks_lab = 'th,', info ='', 
#                       save = True, dpi = 300, ext = 'png'):
    
#     if stats_set[0]: 
#         tests_res = []
#         statist = True
#         dict_stats = stats_set[1]
        
#     print('\n>> '+ title+' - x_var: '+x_var+ ' - hue_var: '+hue_var+ ' - shape_var: '+shape_var)
#     sns.set_context('poster') # notebook, talk, poster, paper
#     # Get legends
#     dict_legends = def_legends(df2plot)
    
#     # Set up the matplotlib figure
#     num_vars = len(vars2plot)
#     n_rows = math.ceil(num_vars/n_cols)+1
#     print('n_rows:' , n_rows)
    
#     index_right_col = list(range(n_cols,(n_cols+1)*n_rows,4))
#     index_no_graph = list(range(num_vars, (n_cols+1)*n_rows))
#     index_no_plot = sorted(list(set(index_right_col).union(set(index_no_graph))))
#     print(index_no_plot)

#     for index in index_no_plot:
#         vars2plot.insert(index, '')
#         labels2plot.insert(index, '')
    
#     # Set up the matplotlib figure
#     h_plot, w_plot = ips
#     if num_vars == 1:
#         h_add = 0; w_add = 0
#     size_col = (n_cols+1)*h_plot+h_add
#     size_row = n_rows*w_plot+w_add
    
#     # Genotypes and Strains being plotted 
#     values = []
#     for var in [x_var, hue_var, shape_var]:
#         if var =='GenotypeAll':
#             reverse = True
#         else: 
#             reverse = False
#         values.append(sorted(df2plot[var].unique(), reverse=reverse))
#     x_values, hue_values, shape_values = values
    
#     # - number of x_var
#     n_x = len(x_values)
    
#     #  Create figure  - plt.clf()
#     gridkw = dict(width_ratios=[1]*n_cols+[0.2], height_ratios = [1,1])
#     fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
#     fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
#     sns.set_style("ticks")
#     sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
#                                                     "ytick.labelsize":8, "lines.linewidth": 2.5,
#                                                     "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
    
#     # Define legends for shape and hue
#     hue_legend = dict_legends[hue_var]
#     legend_elem_hue = []
#     for aa, hue_val in enumerate(hue_legend):
#         legend_elem_hue.append(Line2D([0], [0], marker='h', color='w', label=hue_val,
#                                 markerfacecolor=palettes[aa], markersize=20))
#     space = [Line2D([0], [0], marker='o', color='w', label='',
#                                 markerfacecolor='w', markersize=20)]
#     shape_legend = dict_legends[shape_var]
#     legend_elem_shape = []
#     for n_str, shape_val, mark in zip(count(), shape_legend, styles):
#         legend_elem_shape.append(Line2D([0], [0], marker=mark, color='w', label=shape_val,
#                                 markerfacecolor='k', markersize=15))
        
#     handle_new = legend_elem_hue+space+legend_elem_shape
#     legend_new = hue_legend+['']+shape_legend
    
#     marker_size = 10
#     dodge = True
#     for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
        
#         # Create list to contain info of yticks and ist highest value
#         y_vals_all = []
#         max_y_vals = []
    
#         if n == 0: 
#             for k, svar, value in zip(count(), [x_var, hue_var, shape_var], values):
#                 print('\t- '+svar+': ', value)
        
#         if n in index_no_plot:
#             if n == n_cols:
#                 ax.set_axis_off()
#                 ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
#             else: 
#                 ax.remove()
#         else: 
#             print(' > ',var,'\n')
#             m = sns.swarmplot(data=df2plot, x=x_var, y=var, hue = hue_var, hue_order = hue_values, ax = ax, order=x_values,
#                               palette = palettes, dodge = dodge, size = marker_size)
#             if statist: 
#                 box_pairs = dict_stats[var]['box_pairs']
#                 stat_test = False
#                 test = None
#                 p_val = dict_stats[var]['pval_multComp_all']
                
#                 m1, test_results = add_stat_annotation(ax, data=df2plot, x=x_var, y=var, hue=hue_var, 
#                                                        hue_order = hue_values, order = x_values,
#                                                         box_pairs=box_pairs, perform_stat_test = stat_test, test = test,
#                                                         pvalues=p_val, comparisons_correction=None, #'bonferroni',
#                                                         line_offset_to_box=0.4, line_offset=0.1,
#                                                         line_height=0.015, text_offset=5,
#                                                         text_format='star', loc='inside', verbose=0);
#                 tests_res.append(test_results)
#                 # txt_multcomp = txtMultComp(box_pairs, test_results)
#                 # print('txt_multcomp:', txt_multcomp)
                
#             box = ax.get_position()
#             if yticks_lab == '1e6 - d.':
#                 ylabel = ylabel +' x 10$^6$'
#             elif yticks_lab == '1e3 - d.':
#                 ylabel = ylabel +' x 10$^3$'
#             ax.set(xlabel=dict_legends['xlabels'][x_var], ylabel=ylabel);
#             # ax.set_xticks(ax.get_xticks())
#             ax.set_xticklabels(dict_legends[x_var], rotation=0)
#             ax.set_position([box.x0, box.y0, box.width*1, box.height])
#             ax.get_legend().remove()
#             sns.despine()
            
#             if n == 0:
#                 handles, labels = m.get_legend_handles_labels()

#             y_vals = ax.get_yticks()
#             y_vals_all.append(y_vals)
#             max_y_vals.append(y_vals[-1])
#             if n_x > 1: 
#                 vline_pos = list(range(1,n_x,1))
#                 for pos in vline_pos:
#                     ax.axvline(pos - 0.5, ymax = 0.95, color='dimgrey', ls='-.', linewidth=0.8)
            
#             # Define axes based on higherst bar
#             max_y_index = max_y_vals.index(max(max_y_vals))
#             ax.set_yticks(y_vals_all[max_y_index])
#             if yticks_lab == '1e6 - d.':
#                 ax.set_yticklabels(['{:.2f}'.format(w/1e6) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == '1e3 - d.':
#                 ax.set_yticklabels(['{:.0f}'.format(w/1e3) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'th,':
#                 ax.set_yticklabels([locale.format("%d", w, grouping=True) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 0':
#                 ax.set_yticklabels(['{:.0f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd. - 1':
#                 ax.set_yticklabels(['{:.1f}'.format(w) for w in y_vals_all[max_y_index]])
#             elif yticks_lab == 'd.':
#                 ax.set_yticklabels(['{:.2f}'.format(w) for w in y_vals_all[max_y_index]])

#     fig.suptitle(title+'\n', fontsize = 30, y=1)
#     if save: 
#         for extf in ext: 
#             dir2savef = os.path.join(dir2save, 'pl_meas', 'R_')
#             if info != '':
#                 fig_title = dir2savef+info+"_"+title+"."+extf
#             else: 
#                 fig_title = dir2savef+title+"."+extf

#             plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
    
#     if statist: 
#         return tests_res

#%% func - plotPerVariableLabels
# def plotPerVariableLabels(script, input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
#                      h_plot, w_plot, save, dir2save, info, dpi = 300, h_add = 5, w_add = 1):
    
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
        
#         palettes = ['mediumturquoise', 'darkmagenta']
        
#         for nn, var, ylabel in zip(count(), vars2plot, labels2plot):
#             #  Create figure  - plt.clf()
#             gridkw = dict(width_ratios=[1,1,1,0.2])
#             fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
#             fig.subplots_adjust(hspace=1.5, wspace=0.05)
            
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
#             jitter = 0.2
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
#                     ax.set(xlabel="\nStage: "+stg+'hpf')
#                     ymin, ymax = ax.get_ylim()
#                     up = (ymax-ymin)//20
#                     # print(up)
#                     df_plot = df_plot.sort_values(by=['Stage','Strain', var])
#                     fish_refs = df_plot.Ref.unique()
#                     ha = ['right', 'left', 'center']*50
#                     for j, ref in enumerate(fish_refs):
#                         strain_pos =strains.index(df_plot['Strain'].values[j])
#                         gen_pos = genots.index(df_plot['GenotypeAll'].values[j])
#                         # print(ref, ha[j])
#                         # rand_x = random.uniform(0.95, 1.05)
#                         # rand_y = random.uniform(0.95, 1.05)
#                         # ax.text(x=strain_pos*rand_x, y=df_plot[var].values[j]*rand_y, s=ref, horizontalalignment=ha[j], size=10, color=palettes[gen_pos])
#                         ax.text(x=strain_pos, y=df_plot[var].values[j]+up, s=ref, horizontalalignment=ha[j], size=10, color=palettes[gen_pos])
    
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
#             dir2savef = os.path.join(dir2save, 'pl_labels', 'R_')
#             if info != '':
#                 fig_title = dir2savef+"Lab_"+info+"_"+var+".png"
#             else: 
#                 fig_title = dir2savef+"Lab_"+var+".png"
            
#             if save: 
#                 plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
                
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


#%% func - plotInGroups3
# def plotInGroups3 (plot_type, input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
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




