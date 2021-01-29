# -*- coding: utf-8 -*-
"""
morphoHeart_funcAnalysis

Version: Nov, 2020
@author: Juliana Sanchez-Posada

"""
#%% Importing python packages


#%% Importing morphoHeart packages
from .morphoHeart_funcBasics import alert

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

    print('\nVariables:')
    for c, value in enumerate(list_vars, 1):
        print(c-1, value)
    input_var = input('Select the '+ name +' you would like to process: ')

    if input_var == 'All':
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

#%% func - def_variables
def def_variables(module):
    if module == 'morphoHeart_D_AnalyseData':
        variables = ["SurfArea_Myoc","SurfArea_Int.Myoc","SurfArea_Ext.Myoc",
                     "SurfArea_Endo","SurfArea_Int.Endo","SurfArea_Ext.Endo",
                     "SurfArea_CJ","SurfArea_Int.CJ","SurfArea_Ext.CJ",
                     "Vol_Int.Myoc","Vol_Ext.Myoc",
                     "Vol_Int.Endo","Vol_Ext.Endo",
                     "linLine_Int.Myoc(Cut)","linLine_Ext.Endo(Cut)",
                     "Length_CL_Int.Myoc(Cut)","Length_CL_Ext.Endo(Cut)",
                     "Vol_Myoc","Vol_Atr.Myoc","Vol_Vent.Myoc",
                     "Vol_Endo","Vol_Atr.Endo","Vol_Vent.Endo",
                     "Vol_CJ","Vol_Atr.CJ","Vol_Vent.CJ",
                     "ang_Heart","ang_Atr","ang_Vent","ang_BtwChambers",
                     'Looping_Ratio_Myoc','Looping_Ratio_Endo',
                     'Vol_Atr.ExtMyoc','Vol_Vent.ExtMyoc','Vol_Atr.IntEndo','Vol_Vent.IntEndo']

        ylabels  = ["Surface Area Myocardium [um$^2$]","Surface Area Int.Myocardium [um$^2$]","Surface Area Ext.Myocardium [um$^2$]",
                     "Surface Area Endocardium [um$^2$]","Surface Area Int.Endocardium [um$^2$]","Surface Area Ext. Endocardium [um$^2$]",
                     "Surface Area CJ [um$^2$]","Surface Area Int.CJ [um$^2$]","Surface Area Ext.CJ [um$^2$]",
                     "Volume Int.Myocardium [um$^3$]","Heart Volume [um$^3$]",
                     "Heart Lumen Volume [um$^3$]","Volume Ext.Endocardium [um$^3$]",
                     "Linear Heart Length (Int.Myoc) [um]","Linear Heart Length (Ext.Endo) [um]",
                     "Looped Heart Length (Int.Myoc) [um]","Linear Heart Length (Ext.Endo) [um]",
                     "Volume Myocardium [um$^3$]","Atrial Volume Myocardium [um$^3$]","Ventricular Volume Myocardium [um$^3$]",
                     "Volume Endocardium [um$^3$]","Atrial Volume Endocardium [um$^3$]","Ventricular Volume Endocardium [um$^3$]",
                     "Volume CJ [um$^3$]","Atrial Volume CJ [um$^3$]","Ventricular Volume CJ [um$^3$]",
                     "Angle Heart wrt Sample (\N{DEGREE SIGN})","Atrial Angle (\N{DEGREE SIGN})","Ventricular Angle (\N{DEGREE SIGN})","Angle between Chambers (\N{DEGREE SIGN})",
                     'Looping Ratio (Int.Myoc)','Looping Ratio (Ext.Endo)',
                     'Atrial Volume [um$^3$]','Ventricular Volume [um$^3$]','Atrial Lumen Volume [um$^3$]','Ventricular Lumen Volume [um$^3$]']
    
    return variables, ylabels
    