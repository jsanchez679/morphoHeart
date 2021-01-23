# -*- coding: utf-8 -*-
"""
morphoHeart_funcAnalysis

Version: Nov, 2020    
@author: Juliana Sanchez-Posada

"""
#%% Importing python packages


#%% Importing morphoHeart packages
from morphoHeart_funcBasics import alert

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
    for c, value in enumerate(variables, 1):
        print(c-1, value)
    input_var = input('Select the '+ name +' you would like to process: ')
    
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
   
    variables2loop = []
    for i, num in enumerate(var_num):
        variables2loop.append(variables[num])

    return variables2loop



