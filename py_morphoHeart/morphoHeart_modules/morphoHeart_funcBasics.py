# -*- coding: utf-8 -*-
"""
morphoHeart_funcBasics

Version: Jan 22, 2021
@author: Juliana Sanchez-Posada

"""

#%% Importing python packages
import os
import pandas as pd
from playsound import playsound
from pathlib import Path
import numpy as np
from progress.bar import Bar
from itertools import count
import json
import shutil
import psutil

suffix = '%(index)d/%(max)d - %(elapsed)ds'

#https://thispointer.com/python-pandas-how-to-display-full-dataframe-i-e-print-all-rows-columns-without-truncation/
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

#%% class - NumpyArrayEncoder
# Definition of class to save dictionary
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

#%% func - alert
def alert(sound, times):
    """
    Function to play selected sound as an alert

    Parameters
    ----------
    sound : str
        'frog'/'clown'/'countdown'/'wohoo'/'whistle'/'jump'/'beep', etc.
    times : int
        number of times to play the sound.

    Returns
    -------
    None.

    """

    if sound == 'frog':
        sound2play ='Sounds/Frog-sound-ribbit.mp3'
    elif sound == 'clown':
        sound2play = 'Sounds/Mistake-sound-effect.mp3'
    elif sound == 'countdown':
        sound2play ='Sounds/Countdown-start.mp3'
    elif sound == 'wohoo':
        sound2play = 'Sounds/Cartoon-woohoo.mp3'
    elif sound == 'whistle':
        sound2play = 'Sounds/Swanee-whistle-down.mp3'
    elif sound == 'jump':
        sound2play = 'Sounds/Jumping-sound-effect.mp3'
    elif sound == 'beep':
        sound2play = 'Sounds/Error-beep-sound-effect.mp3'
    elif sound == 'bubble':
        sound2play = 'Sounds/Bubble-sound-effect.mp3'
    elif sound == 'dig_camera':
        sound2play = 'Sounds/Digital-camera-click-sound.mp3'
    elif sound == 'dig_pocket':
        sound2play = 'Sounds/Digital-pocket-camera-click-sound.mp3'
    elif sound == 'light':
        sound2play = 'Sounds/Light-speed-sound-effect.mp3'
    else:
        sound2play = 'Sounds/Error-sound.mp3'

    for i in range(times):
        playsound(sound2play)

#%% func - getRAMUse
def getRAMUse():
    """
    Funtion that prints RAM memory being used

    Returns
    -------
    None.

    """
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2.**30  # memory use in GB...I think
    print('RAM Memory Use:', memoryUse)

#%% func - ask4input
def ask4input(text, type_response, keep = False):
    """
    Function that ask for user input and transforms it into the expected type

    Parameters
    ----------
    text : str
        Text asking the question and giving possible options.
    type_response : data types int/str/float/boolean
        data type of the expected response
    keep : Boolean
        If True, leaves the string as it is (with upper and lower cases), else, changes the whole string to lower case.

    Returns
    -------
    response : int/str/float/bool
        returns an object with the corresponding type_response given as input.

    """
    alert('beep',1)
    exit_now = False
    while exit_now == False:
        response = input('> '+text+' ')
        if type_response == int:
            try:
                response = int(response)
                exit_now = True
            except:
                print('ERROR: -'+response+'- The number entered needs to be an integer!')
                pass
        elif type_response == float:
            try:
                response = float(response)
                exit_now = True
            except:
                print('ERROR: -'+response+'- The number entered needs to be a float!')
                pass
        elif type_response == str:
            if keep == False:
                response = response.lower()
            exit_now = True
        elif type_response == bool:
            try:
                if int(response) in [0,1]:
                    response = bool(int(response))
                    exit_now = True
            except:
                print('ERROR: -'+response+'- The number entered needs to be a [0]:no or [1]:yes!')
                pass

    return response

#%% func - getMainDirectories
def getMainDirectories(root_path):
    """
    Function to get the main directories of py_morphoHeart

    Parameters
    ----------
    root_path : path
        Path to py_morphoHeart folder.

    Returns
    -------
    dir_pyLSAnalysis : path
        Path to py_LSAnalysis folder.
    dir_lsAnalysis : path
        Path to lsAnalysis folder.
    dir_lsOngoing : path
        Path to lsOngoing folder.
    dir_data2Analyse : path
        Path to data2Analyse folder (LS_MorphoHeart).

    """
    # dir_pyLSAnalysis = Path(root_path).parent
    # dir_lsAnalysis = Path(dir_pyLSAnalysis).parent
    # dir_lsOngoing = Path(dir_lsAnalysis).parent
    # dir_data2Analyse = Path(dir_lsAnalysis, "LS_MorphoHeart")

    dir_morphoHeart = Path(root_path).parent
    dir_parent = Path(dir_morphoHeart).parent
    # dir_lsOngoing = Path(dir_parent).parent
    dir_im_morphoHeart = Path(dir_parent, "im_morphoHeart")


    # return dir_pyLSAnalysis, dir_lsAnalysis, dir_lsOngoing, dir_data2Analyse
    return dir_morphoHeart, dir_parent, dir_im_morphoHeart

#%% func - exportDatasetCSV
def exportDatasetCSV(dir_data2Analyse, end_name = '2A'):
    """
    Function that finds all the folders within the dir_data2Analyse folder whose name ends with 2A (a.k.a: to analyse)
    and creates a dataframe including the embryo information (strain, stage, manipulation, genotype and lightsheet ref)
    The resulting dataframe is saved as a '.csv' file within the dir_data2Analyse folder.

    Parameters
    ----------
    dir_lsOngoing : path
        Path to lsOngoing folder where the genotyping record is saved
    dir_data2Analyse : path
        Path to data2Analyse folder (LS_MorphoHeart).
    end_name : str
        String telling which type of folder to look for.

    Returns
    -------
    df_dataset : dataframe
        Resulting dataframe containing all the embryos information of the folders '_2A'.

    """

    # - Location of genotyping record
    genot_file = 'LS_Genotyping Record.xlsx'
    genot_dir = os.path.join(dir_data2Analyse, genot_file)
    # - Get general info of all LS Sessions
    info_cols = ['LS No','Strain','Stage','Manipulation']
    df_info = pd.read_excel(genot_dir, sheet_name = 'RAW_All', header = 0,
                            usecols=info_cols, engine='openpyxl')
    df_info = df_info.set_index('LS No')
    # - Get genotyping data of all LS Sessions
    df_genot = pd.read_excel(genot_dir, sheet_name = 'All_Genotype', header = 1, engine='openpyxl')
    df_genot = df_genot.set_index('Pos')

    if end_name == 'R':
        glob_txt = '*R_L*'
    else:
        glob_txt = '*_2A*'

    df_dataset = pd.DataFrame()
    for folder_2A in dir_data2Analyse.glob(glob_txt):
        folder = os.path.basename(folder_2A)
        # print(folder)
        # folder = str(folder)
        if end_name == '2A':
            ls_Session = folder.split('_')[0]
            fish_ref = folder.split('_')[1]
        elif end_name == 'R':
            ls_Session = folder.split('_')[1]
            fish_ref = folder.split('_')[2]
        # print(ls_Session, fish_ref)

        stage = df_info.loc[ls_Session,'Stage']
        strain = df_info.loc[ls_Session,'Strain']
        manip = df_info.loc[ls_Session,'Manipulation']

        # - Get genotype(s) of fish
        genotype_A = df_genot.loc[fish_ref,ls_Session]
        gene_A = df_genot.loc['Gene',ls_Session]

        ls_Session_B = ls_Session + '.1'
        if ls_Session_B in df_genot.columns:
            genotype_B = df_genot.loc[fish_ref,ls_Session_B]
            gene_B = df_genot.loc['Gene',ls_Session_B]
        else:
            genotype_B = '-'
            gene_B = '-'

        new_row = {'Folder': folder, 'LS_Session': ls_Session, 'Fish_ref': fish_ref,
                    'Strain': strain, 'Stage': stage, 'Manip': manip,
                    'Gene_A': gene_A, 'Genotype_A': genotype_A,
                    'Gene_B': gene_B, 'Genotype_B': genotype_B}

        df_dataset = df_dataset.append(new_row, ignore_index=True)
        df_dataset = df_dataset[['Folder','LS_Session', 'Fish_ref','Strain','Stage', 'Manip',
                                  'Gene_A', 'Genotype_A','Gene_B', 'Genotype_B']]

        dir_dfDataset = os.path.join(dir_data2Analyse,'df_dataset'+end_name+'.csv')

    df_dataset.to_csv(dir_dfDataset)

    return df_dataset

#%% func - loadDF
def loadDF(filename, file, dir_results):
    """
    Function that loads a dataframe

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    file : str
        Particular name given to the dataframe that wants to be loaded (e.g. ResultsDF).
    dir_results : path
        Path to the folder where the dataframe is saved.

    Returns
    -------
    df_results : dataframe
        Loaded dataframe.

    """

    name_csv = filename+'_'+file+'.csv'
    dir_filledDF = os.path.join(dir_results, name_csv)
    df_results = pd.read_csv(dir_filledDF, index_col=0)

    return df_results

#%% func - loadNPY
def loadNPY(filename, names, dir_txtNnpy, print_txt = True):
    """
    Function that loads a (group of) npy array(s)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    names : list of str
        List of particular names given to the npy arrays that want to be loaded (e.g. ['stackShape', 's3_ch0']).
    dir_txtNnpy :  path
        Path to the npy arrays folder.
    print_txt : Bool
        True if text is to be printed, else False. The default is True.

    Returns
    -------
    npys : list of arrays
        List of Loaded array(s).

    """
    # print('- Loading np arrays...')
    npys = []
    if print_txt:
        bar = Bar('- Loading np arrays', max=len(names), suffix = suffix, check_tty=False, hide_cursor=False)
    for i, name in enumerate(names):
        npy_title = filename+'_'+name+'.npy'
        npy_dir = os.path.join(dir_txtNnpy, npy_title)
        npy = np.load(npy_dir)
        npys.append(np.asarray(npy))
        if print_txt:
            bar.next()
            
    if print_txt:
        bar.finish()
        alert('jump',1)

    return npys

#%% func - loadDicts
def loadDicts(filename, dicts_name, directories, print_txt = True):
    """
    Function that loads a (group of) dictionary(ies)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    dicts_name :  list of str
        List of particular names given to the dictionaries that want to be loaded (e.g. ['dict_obj']).
    directories : list of paths
        List of paths to the folders where the dictionaries are saved.
    print_txt : Bool
        True if text is to be printed, else False. The default is True.

    Returns
    -------
    dicts_out : list of dictionaries
        List of Loaded dictionaries.

    """
    if print_txt:
        print('- Importing dictionaries')
    dicts_out = []
    for i, name, sp_dir in zip(count(), dicts_name, directories):

        jsonDict_name = filename+"_"+name+".json"
        json2open_dir = os.path.join(sp_dir,jsonDict_name)
        #print("Started Reading JSON file")
        with open(json2open_dir, "r") as read_file:
            if print_txt:
                print("\t>> "+name+": Converting JSON encoded data into Numpy Array")
            decodedArray = json.load(read_file)

        dicts_out.append(decodedArray)

    alert('jump',1)

    return dicts_out

#%% func - selectFile
def selectFile (df_dataset):
    """
    Function that allows the user to select the file to be processed, printing the names of all the folders in the
    input dataframe and prompting the user to input the file number to process.

    Parameters
    ----------
    df_dataset : dataframe
        Dataframe containing all the embryos information of the folders '_2A'.

    Returns
    -------
    filename : str
        Reference name given to the images of the embryo selected to be processed (LSXX_FXX_X_XX_XXXX).
    df_file : dataframe
        Dataframe just with the information of the filename selected by the user.
    input_num :  int
        Index number selected from the input dataframe.

    """

    df_datasetTemp = df_dataset.copy()
    df_datasetTemp['Folder_2A'] = df_dataset['Folder'].str.slice(0,10)
    df_datasetTemp['Gene A'] = df_dataset['Gene_A']+': '+df_dataset['Genotype_A']
    df_datasetTemp['Gene B'] = df_dataset['Gene_B']+': '+df_dataset['Genotype_B']

    print('\n--- HEARTS TO BE PROCESSED ---')
    df_datasetTemp = df_datasetTemp.sort_values(by=['Stage','Strain','Gene_A','Genotype_A','Gene_B','Genotype_B'],
                                                ascending = (True, True, True, False, True, False))
    print(df_datasetTemp[['Folder_2A','Strain','Stage','Gene A','Gene B']])

    list_folder = df_dataset['Folder'].tolist()

    # print('\nFiles to process:')
    # for c, value in enumerate(list_folder, 1):
    #     print('-',c-1,':\t', value)

    alert('beep',1)
    input_num = int(input('> Select the file you want to process: '))
    filename = list_folder[input_num]

    print('\n\t\t\tFILE SELECTED ')
    df_file = df_dataset.loc[df_dataset['Folder'] == filename]
    print(df_file.T,'\n')

    return filename, df_file, input_num

#%% func - selectHearts
def selectHearts (df_dataset):
    """
    Function that allows the user to select the hearts to plot, printing a modified version of the
    input dataframe and prompting the user to input the file number(s)

    Parameters
    ----------
    df_dataset : dataframe
        Dataframe containing all the embryos information of the folders 'R_'.

    Returns
    -------
    df_file : dataframe
        Dataframe just with the information of the filename(s) selected by the user

    """


    df_datasetTemp = df_dataset.copy()
    df_datasetTemp['Ref'] = df_dataset['Folder'].str.slice(2,10)
    df_datasetTemp['Gene A'] = df_dataset['Gene_A']+': '+df_dataset['Genotype_A']
    df_datasetTemp['Gene B'] = df_dataset['Gene_B']+': '+df_dataset['Genotype_B']

    # q_filter = ask4input('Do you want to filter hearts by genotype and/or stage? \n\t[0]:yes, genotype / [1]:yes, stage / [2]:yes, both / [3]:no, manually select hearts: ', int)
    # if q_filter in [0,2]:
    #     genes = df_dataset[["Gene_A", "Gene_B"]].values.ravel()
    #     genes = pd.unique(genes)
    #     genes = genes.tolist()
    #     genes.remove('-')

    #     print('\nGenes:')
    #     for c, value in enumerate(genes, 1):
    #         print('-',c-1,':\t', value)
    #     input_var = input('Select the genes by which you would like to filter: ')
    #     genes_num = getInputNumbers(input_var, genes)
    #     genes_selected = [genes[num] for i, num in enumerate(genes_num)]

    #     numOfGenes = len(genes_selected)
    #     if numOfGenes == 1:
    #         df_filtered = df_dataset[(df_dataset['Gene_A'] == genes_selected[0]) & (df_dataset['Gene_B'] == '-')]
    #     elif numOfGenes == 2:
    #         df_filtered = df_dataset[(df_dataset['Gene_A'] == genes_selected[0]) & (df_dataset['Gene_B'] ==  genes_selected[1])]

    #     q_plotWtHtMt = ask4input("- Which genotypes would you like to include? \n\t-(e.g. write 'wt,ht,mt' if you want them all, or 'wt' in case you only want wild-types) : ",str)
    #     ### CONTINUE...
    #     #END list_folder

    # elif q_filter == 3:
    print('-------- PROCESSED HEARTS --------')
    df_datasetTemp = df_datasetTemp.sort_values(by=['Stage','Strain','Gene_A','Genotype_A','Gene_B','Genotype_B'],
                                                ascending = (True, True, True, False, True, False))
    print(df_datasetTemp[['Ref','Strain','Stage','Gene A','Gene B']])
    list_folder = df_dataset['Folder'].tolist()
    alert('beep',1)
    input_num = ask4input('Enter the numbers of the hearts you want to plot: ', str)
    hearts_num = getInputNumbers(input_num, list_folder)

    print('\n-------- SELECTED HEARTS --------')
    df_file = df_dataset.loc[hearts_num, :]

    print(df_file.T,'\n')

    return df_file

#%% func - getInputNumbers
def getInputNumbers(input_var, variables):
    """
    Function that gets the variables corresponding to the user input.

    Parameters
    ----------
    input_var : str
        String with input given by the user.
    variables : list of str
        List with variables to select.

    Returns
    -------
    var_num : list of str
        List of selected variables.

    """
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

    return var_num

#%% func - createDirectories2Save
def createDirectories2Save (filename, dir_data2Analyse, end_name):
    """
    Function to create directories to save data within original folder (v. Sept 03, 2020)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    dir_data2Analyse : path
        Path to data2Analyse folder (LS_MorphoHeart).
    end_name : str
        String with '2A' or 'R', according to the type of folder in which the other directories are to be created.

    Returns
    -------
    dir_results : path
        Path to the Results directory.
    directories : list of paths
        List of paths: 0. dir_dict, 1. dir_txtNnpy, 2. dir_stl, 3. dir_cl, 4. dir_imsNvideos, 5. dir_ims2Analyse, 6. dir_folder2A.

    """

    res_folder = "Results_"+filename

    if end_name == '2A':
        folder2A = filename+'_'+end_name
    elif end_name == 'R':
        folder2A = end_name+'_'+filename

    dir_folder2A = os.path.join(dir_data2Analyse, folder2A)
    dir_results = os.path.join(dir_data2Analyse, folder2A, res_folder)
    if os.path.isdir(dir_results) == False:
        os.mkdir(dir_results)
        print(res_folder, " was created as a directory!")
    else:
        print(res_folder, " is an existing directory!\n")

    dirsResults = ["dictionaries", "txt_npy", "meshes", "centreline", "imgs_videos"]
    directories = []

    for num, direc in enumerate(dirsResults):
        dir2create = os.path.join(dir_data2Analyse, folder2A, res_folder, direc)
        directories.append(dir2create)
        if os.path.isdir(dir2create) == False:
            os.mkdir(dir2create)
            print("\t-",direc, " was created as a directory!")
        # else:
        #     print("\t-",direc, " is an existing directory!")

    dir_ims2Analyse = os.path.join(dir_data2Analyse,folder2A,'Im_'+filename)
    directories.append(dir_ims2Analyse)
    directories.append(dir_folder2A)

    dir_dict, dir_txtNnpy, dir_stl, dir_cl, dir_imsNvideos, dir_ims2Analyse, dir_folder2A = directories

    return dir_results, directories

#%% func - metadataExt
def metadataExt (filename, dir_data2Analyse):
    """
    Function that gets the metadata from the filename given as input and exports the x, y and z spacing as a npy array. (v. April 17, 2020)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    dir_data2Analyse : path
        Path to data2Analyse folder (LS_MorphoHeart).

    Returns
    -------
    xy_Scaling_um : float
        x and y scaling taken from the metadata of the file.
    z_Scaling_um : float
        z scaling taken from the metadata of the file.

    """

    Metadata = "METADATA_RAW_"+filename[0:4]+".txt"
    dir_Metadata = os.path.join(dir_data2Analyse,'Z_METADATA',str(Metadata))
    usecols_MD = ["ImageTitle","FileName", "ScalingX", "ScalingY", "ScalingZ", "DimensionZ"]
    df_Metadata = pd.read_csv(dir_Metadata, delimiter="\t", header = 1, index_col = 0, usecols = usecols_MD)
    index_stack =  filename+".czi"

    # xyz scaling from Metadata
    xy_Scaling = df_Metadata.at[index_stack, "ScalingX"]
    z_Scaling = df_Metadata.at[index_stack, "ScalingZ"]
    # xyz scaling in micrometers (um)
    xy_Scaling_um = xy_Scaling*(10**6)
    z_Scaling_um = z_Scaling*(10**6)

    return xy_Scaling_um, z_Scaling_um

#%% func - code4vmtkCL (Changes needed - list/for)
def code4vmtkCL(filename, mesh_name, dir_cl, printshow = True):
    """
    Function that gets directories information and prints a series of instructions to process the meshes to obtain centreline and
    run the vmtk code.

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    mesh_name : list of str
        List of mesh names.
    dir_cl : path
        Path to the centreline directory.
    printshow : Bool
        True if text is to be printed, else False.

    Returns
    -------
    vmtktxts : list of str
        Codes to run VMTK using pypes to extract centreline of the internal myoc and/or external endoc

    """
    
    mesh_titles = []
    vmtktxts = []
    names = []
    meshML_dirs = []
    
    vmtktxtsf = []
    namesf = []

    for i in range(len(mesh_name)):
        if printshow:
            name = mesh_name[i][0:8]
        else: 
            name = mesh_name[i]
        names.append(name)
        mesh_title = filename+"_"+name+"_cut4cl.stl"
        mesh_titles.append(mesh_title)
        meshML_title = filename+"_"+name+"_cut4clML.stl"
        meshML_dir = '"'+os.path.join(dir_cl, meshML_title)+'"'
        meshML_dirs.append(os.path.join(dir_cl, meshML_title))
        cl_title = filename+"_"+name+"_cl.vtp"
        cl_dir = os.path.join(dir_cl, cl_title)
        cl_dir_str = '"'+cl_dir+'"'
        vmtktxt = "vmtksurfacereader -ifile "+ meshML_dir +" --pipe vmtksurfacesmoothing -passband 0.1 -iterations 30 --pipe vmtkcenterlines -seedselector openprofiles -ofile"+ cl_dir_str+ " --pipe vmtkrenderer --pipe vmtksurfaceviewer -opacity 0.25 --pipe vmtksurfaceviewer -i @vmtkcenterlines.o -array MaximumInscribedSphereRadius"
        vmtktxts.append(vmtktxt)

    if printshow:
        print("\nYou are done in python for a little while... \n\t\t\tto get the centreline of each of the selected meshes follow the next steps:")
        print(" >>> 1. Open the file(s):  -", mesh_titles," - in Meshlab")
        print(" >>> 2. Run Filters > Remeshing, Simplification.. > Screened Poisson Surf Reco (check Pre-clean)")
        print(" >>> 3. Cut inflow/outflow tract and export the resulting surface adding 'ML' at the end of the filename \n\t\t\t(e.g _cut4clML.stl) in the same folder")
        print(" >>> 4. Back in Spyder, open script morphoHeart_B_VMTK.py and run it...")
        # print(str(vmtktxtA), '\n\n\n', str(vmtktxtB))
        
        vmtktxtsf = vmtktxts
        namesf = names
    
    else:
        for j, ml_dir in enumerate(meshML_dirs):
            # print(cldir)
            if os.path.exists(ml_dir):
                vmtktxtsf.append(vmtktxts[j])
                # print('yay')
                namesf.append(names[j])

    alert("wohoo",1)
        
    return vmtktxtsf, namesf

#%% func - saveFilledDF
def saveFilledDF(filename, df_res, dir2save):
    """
    Function that exports the filled dataframe containing all the measured data obtained in the processing of the images
    (e.g. surface area, volumes, orientation angles, etc) as a '.csv' file in the directory given as input

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    df_res : dataframe
        Dataframe containing all the measured data.
    dir2save : path
        Path to the folder where the dataframe is to  be saved.

    Returns
    -------
    None.

    """

    print('\n- File dataframe for '+filename+' has been saved!')
    # print(df_res.T.head(1:2),'\n')

    name_csv = filename+'_ResultsDF.csv'
    dir_filledDF = os.path.join(dir2save, name_csv)
    df_res.to_csv(dir_filledDF)

    print('-', name_csv, ' has been saved!')
    alert('wohoo',1)

#%% func - saveDF
def saveDF(filename, df2save, df_name, dir2save):
    """
    Function that exports a dataframe given as input as a '.csv' file in the directory given as input as well.

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    df2save : dataframe
        Dataframe to save.
    df_name : str
        Specific name given to the dataframe to be saved.
    dir2save : path
        Path to the folder where the dataframe is to  be saved.

    Returns
    -------
    None.

    """

    name_csv = filename+'_'+df_name+'.csv'
    dir_DF = os.path.join(dir2save, name_csv)
    df2save.to_csv(dir_DF)

    print('-', df_name, ' has been saved!')
    alert('wohoo',1)

#%% func - saveDict
def saveDict(filename, dict2save, name, dir2save, print_txt = True):
    """
    Functions that saves dict given as input

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    dict2save : dict
        Dictionary to be saved.
    name :  str
        Specific name given to the dictionary to be saved.
    dir2save : path
        Path to the folder where the dictionary is to be saved.
    print_txt : bool, optional
        True if confirmation of action is needed, else False. The default is True.

    Returns
    -------
    None.

    """

    jsonDict_name = filename+"_"+name+".json"
    json2save_dir = os.path.join(dir2save,jsonDict_name)

    with open(json2save_dir, "w") as write_file:
        json.dump(dict2save, write_file, cls=NumpyArrayEncoder)
    if print_txt:
        print("\t>> Dictionary saved correctly!\n\t> File: "+jsonDict_name);
        alert("countdown",1)

#%% func - printTime
def printTime(tic, toc, txt):
    """
    Function that prints a text with the time taken to perform a process

    Parameters
    ----------
    tic : float
        Time for starting stopwatch.
    toc : float
        Time for ending stopwatch.
    txt : str
        Process being performed

    Returns
    -------
    None.

    """
    time = toc-tic
    print("\t>> Total time taken to "+txt+" = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")


#%% func - changeDirName
def changeDirName(filename, dir_o):
    """
    Function that changes the name of the folder_2A to R_folder, once the processing of that particular file is complete
    (points have been classified into regions and all measurements have been made)

    Parameters
    ----------
    filename : str
        Reference name given to the images of the embryo being processed (LSXX_FXX_X_XX_XXXX).
    dir_o : path
        Path of the original folder.

    Returns
    -------
    None.

    """
    change_name = ask4input("Are you done processing -"+filename+"- and want to change the folder's name? [0]:no/[1]:yes: ", bool)
    if change_name:
        dir_data2Analyse = Path(dir_o).parent
        new_name = 'R_'+filename
        dir_new = os.path.join(dir_data2Analyse, new_name)
        shutil.move(dir_o, dir_new)

        print('- Processed folder has changed name - ', new_name)

#%% Code to create documentation
# from pathlib import Path
# # pdoc --html <path to code> --output-dir <path to documentation>
# (morphoHeart) D:\Documents JSP\Dropbox\Dropbox_Juliana\PhD_Thesis\Data_ongoing\LS_ongoing\A_LS_Analysis\morphoHeart\py_morphoHeart>
#       pdoc --html --output-dir docs_morphoHeart morphoHeart_modules

#%% Code to create reqs
# (morphoHeart) D:\Documents JSP\Dropbox\Dropbox_Juliana\PhD_Thesis\Data_ongoing\LS_ongoing\A_LS_Analysis\morphoHeart>
#       pipreqs .

#%% - ALERT WHEN IMPORTED
print ("IMPORTED: morphoHeart_funcBasics")
alert('jump',1)
