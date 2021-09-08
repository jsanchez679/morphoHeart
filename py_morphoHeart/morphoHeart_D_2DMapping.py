# -*- coding: utf-8 -*-
"""
morphoHeart - D. TRANSFORM ALL THICKNESS DATA INTO 2D
Welcome to the fifth code of morphoHeart!
If you are running this code you must have already divided the heart tissue layers into atrium and ventricle and created
heatmap visualisation of the layers' thickness and heart ballooning. The objective of this code is to transform all the 
thickness and balooning data  you got in the last one into a 2D representation, so that then you can compare it easier 
to other hearts with different genotypes or manipulations. In the end, you can also create distribution plots to compare
the cardiac jelly distribution between different regions of the heart. All the plots created in this script will be saved 
in the 'imgs_videos' folder of each heart, and the numerical data in a csv file within the Results folder. 

When you are done with this code celebrate as this is the last script you need to run with each heart!!!! Woohoo!!
Happy 2D mapping!

@author: Juliana Sanchez-Posada
Version: 15th April, 2021
"""
#%% Importing python packages
import os
import numpy as np
from time import perf_counter
from vedo import Plotter, Cylinder, settings, Text2D
from vedo import embedWindow
embedWindow(False)
settings.legendSize = .3
settings.legendFont="VTK"

# Verify working dir
def setWorkingDir (root_path, init = False):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    init = not bool(int(input('Do you want to execute the script all at once or run it by cells? \n\t[0]: all at once (recommended if you are already familiar with the script)\n\t[1]: by cells (recommended if you are NOT yet familiar with the script). >>>: ')))
    print("- Current working directory: {0}".format(os.getcwd()))
    if not init: 
        print('\nIMPORTANT NOTE:\n- Remember to start running from cell %Start D_2DMapping.\n- NEVER run as an individual cell the cell called %Importing python packages')
    return root_path, init

root_path, init = setWorkingDir(os.getcwd())

c="k"; font= 'VTK';
save = True; plot = True; plotshow = False
azimuth = 0

#%% Start D_2DMapping
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes
    from morphoHeart_modules import morphoHeart_funcAnalysis as fcAn
    tic = perf_counter()

    #%% SELECT FILE AND GET METADATA
    #   This section allows the user to select file to process, get its main directories and metadata,
    #   and define some properties
    #   ================================================================================================================

    # Get main directories
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse)
    # Get file to process and directories
    folder, df_file, file_num, blind = fcBasics.selectFile(df_dataset); filename = folder[0:-3]; dORv = filename[9:10]
    #stage = df_file.loc[file_num,'Stage']
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse, 6. dir_LS_Folder selected
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
    txt = Text2D(filename, c="k", font= 'CallingCode')
    df_res = fcBasics.spAnalysis(df_res, file_num)

    #%% LOAD MESHES, CENTRELINES AND OBJECTS
    #   This section will load all the objects and dataframes needed to unloop the selected heart and create 
    #   distribution plots. It will first ask you to select the thicknesses/ballooning values you would like to
    #   unloop, and then the chamber(s) you would like to unloop. Finally it will load the datasets needed according to 
    #   your selection.
    #   ================================================================================================================

    # Get existing cl_dictionaries
    _, mesh_name = fcBasics.code4vmtkCL(filename = filename, mesh_name = ['myoc_int','endo_ext'],
                                                dir_cl = directories[3], printshow = False)
    # Import dictionaries
    dicts = fcBasics.loadDicts(filename = filename, dicts_name = ['dict_obj']+[txt+'_npcl' for txt in mesh_name],
                                                                    directories = [directories[0]]+ [directories[3]]*len([txt+'_npcl' for txt in mesh_name]))
    dict_obj = fcMeshes.splitDicts(dicts[0])
    if len(dict_obj) == 4:
        [dict_planes, dict_pts, dict_kspl, dict_colour] = dict_obj
    else:
        [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes] = dict_obj
    
    # Select analysis to run
    tissue_analysis, tissue_opt, m_myoc, dict_unloop = fcMeshes.load_tissues2unloop(filename, directories, dir_results, dict_colour)
    
    # Create centreline
    kspl_CL, linLines, ksplSph_o, dict_kspl = fcMeshes.createCLs(dict_cl = dicts[1:], dict_pts = dict_pts, 
                                                                 dict_kspl = dict_kspl, dict_planes = dict_planes, 
                                                                 colors = ['deepskyblue', 'tomato'], myoc = m_myoc)
    
    # Get centreline ribbon
    cl_ribbonV, kspl_ext, _, _, _ = fcMeshes.createCLRibbon(filename = filename, file_num = file_num, 
                                                   df_res = df_res, kspl_CL2use = kspl_CL[0], linLine = linLines[0],
                                                   mesh = m_myoc, dict_kspl = dict_kspl, dict_shapes = dict_shapes, 
                                                   dict_planes = dict_planes, clRib_type = 'extV', plotshow = True)
    
    # Get ch_cylinder
    cyl_chamber = dict_shapes['cyl2CutChambers_final']
    disk = Cylinder(pos = cyl_chamber['cyl_centre'],r = cyl_chamber['radius_max'], height = 2*0.225, 
                    axis = cyl_chamber['cyl_axis'], c = 'purple', cap = True, res = 300)

    #%% UNLOOP HEART
    #   Now we can unloop the heart to create 2D heatmaps of any of the thickness or ballooning measurements already 
    #   taken from this heart. To do this, initially the centreline of the heart is going to be divided into an atrial
    #   and a ventricular centrelines, each containing about 600 points. Each of these centrelines will be then
    #   used to create 150 cross sectional planes of the chamber thickness mesh, all perpendicular to the centreline 
    #   and equally spaced. Thickness information about the points of the mesh that are cut by each of the cross-
    #   sectional planes will be saved as well as information regarding its position (e.g. angular position with 
    #   respect to the ventral center of the heart, as well as cross-sectional plane in which the points were found). 
    #
    #   If plotshow is True in function unloopChambers, as the heart is cross-sectioned, you will be able to see 
    #   a selection of cross-sectional planes and points that are being cut by that plane. You needs to close those 
    #   3D plots for the code to continue processing. 
    #   Once all the information has been acquired for both chambers, a 3D interactive plot will appear showing a 
    #   selection of the planes used to cross-section each chamber. When closed, the code will transform all the 
    #   acquired data and create 2D heatmaps for each chamber and thickness measurement. When all the heatmap
    #   plots have been created, a new 3D interactive plot will appear where you will be able to compare the thickness 
    #   measurements on the 2D heatmap and on the color-coded mesh. 
    #   Note: Make sure you select the 'Plots' Pane in Spyder to be able to see the resulting heatmaps and compare 
    #   them to the color coded meshes. 
    #   The heatmap creation can take longer than 5min/each, and four need to be created. 
    #   Again, be patient, it is worth it! :)
    #   ================================================================================================================
    
    saveHM = True; savePlot = True; plotshow = False; print_txt = False;
        
    # tissue_analysis, tissue_opt, m_myoc, dict_unloop = fcMeshes.load_tissues2unloop(filename, directories, dir_results, dict_colour)
   
    for tissue in tissue_analysis:
        print('>>> Unlooping ', tissue)
        df_AtrVent = np.asarray(dict_unloop[tissue]['df_AtrVent'])
        if plotshow: 
            fcMeshes.plotPtClassif(filename = filename, mesh = dict_unloop[tissue]['mesh'].alpha(1), 
                                    pts_whole = dict_unloop[tissue]['pts_whole'], 
                                    pts_class = dict_unloop[tissue]['pts_class'])

        # Unlooping the tissue chambers
        return_list, kspl_vSurf, _ = fcMeshes.unloopChambers(filename = filename, 
                                        mesh = dict_unloop[tissue]['mesh'].alpha(0.05), 
                                        kspl_CL = kspl_CL[0], kspl_ext = kspl_ext,
                                        no_planes = 150, pl_CLRibbon =  dict_planes['pl_Parallel2LinLine'], 
                                        param = dict_unloop[tissue]['param'], 
                                        param_name = dict_unloop[tissue]['param_name'],
                                        df_AtrVent = dict_unloop[tissue]['df_AtrVent'], 
                                        selected_chambers = dict_unloop[tissue]['selected_chambers'],
                                        dict_kspl = dict_kspl, dict_shapes = dict_shapes, 
                                        dict_pts = dict_pts, dict_planes = dict_planes, 
                                        dir_results = dir_results, 
                                        save_names = dict_unloop[tissue]['save_names'],
                                        plotshow = plotshow, tol=0.05, print_txt=print_txt)
        spheres_zeroDeg, arr_vectZeroDeg = return_list 
        
        print('- Creating heatmaps for each chamber... Note: this can take a while.')
        heatmaps_th, scale_th = fcMeshes.heatmapUnlooped(filename = filename, 
                                                         val2unloop = dict_unloop[tissue]['param_name'], 
                                                         dir_results = dir_results, dirImgs = directories[4], 
                                                         save_names= dict_unloop[tissue]['save_names'],
                                                         hm_names = dict_unloop[tissue]['hm_names'],
                                                         saveHM = saveHM, savePlot = savePlot, cmap = 'turbo')
        del heatmaps_th
    
        _ = fcMeshes.filterUnloopedDF(filename = filename, val2unloop = dict_unloop[tissue]['param_name'], 
                                              dir_results = dir_results, dir_data2Analyse = dir_data2Analyse,
                                              save_names= dict_unloop[tissue]['save_names'],
                                              hm_names = dict_unloop[tissue]['hm_names'],
                                              saveHM = saveHM, cmap='turbo')
        
        if 'kspl_vSurf' not in locals() or 'return_list' not in locals():
            kspl_vSurf = []; return_list = []
        m_Th  = dict_unloop[tissue]['mesh']
        m_Th.pointColors(dict_unloop[tissue]['param'], cmap="turbo", vmin=scale_th[0][0], vmax=scale_th[0][1]).addScalarBar()
        m_Th.mapper().SetScalarRange(scale_th[0][0],scale_th[0][1])
        settings.legendSize = .2
        txt = Text2D(filename+"\n\n >> Plot to validate heatmaps", c="k", font= font)
        vp = Plotter(N=2, axes = 4)
        vp.show(m_Th.alpha(0.01), cl_ribbonV, kspl_vSurf, spheres_zeroDeg, arr_vectZeroDeg, txt, at = 0)
        vp.show(m_Th.clone().alpha(1), at = 1, interactive = True)
    
    #%% GET THICKNESS REGION PLOTS (PROBALITY DENSITY ESTIMATION PLOTS)
    #   Finally, we can create distribution plots that will allow us to compare the cardiac thickness distribution between 
    #   regions of the heart (left vs right, dorsal vs ventral, atrium vs ventricle). 
    #   NOTE: This code can be modified to also get the distribution plots of the thickness of the other tissue layers. 
    #   Just let me know if you want me to include it! :)
    #   ================================================================================================================
    cjPDF = True
    if cjPDF: 
        try: 
            thData = fcBasics.loadDF(filename = filename, file = 'df_cjThNmyocIntBall', dir_results = os.path.join(dir_results, 'csv_all'))
        except: 
            thData = fcBasics.loadDF(filename = filename, file = 'df_cjTh', dir_results = os.path.join(dir_results, 'csv_all'))
                
        dir_kdeIm = fcBasics.new_dir(directories[4], 'kde')
        dir_cjPDFs = fcBasics.new_dir(dir_results, 'csv_all')
        file_num = df_file[df_file['Folder']==folder].index.values[0]
        df_cjPDFs = fcAn.kdeThPlots(filename = filename, df_file = df_file, file_num = file_num, variable = 'cj_thickness', 
                                    thData = thData, dir2save = dir_kdeIm, save = True)
        
        if save: 
            fcBasics.saveDF(filename = filename, df2save = df_cjPDFs, df_name = 'cjPDFs', dir2save = dir_cjPDFs)
            fcBasics.saveDF(filename = filename, df2save = df_cjPDFs, df_name = 'cjPDFs', dir2save = os.path.join(dir_data2Analyse, 'R_All','df_all','df_cjPDFs'))
        
        if save: 
            # Append all dicts to one object dict
            dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
                                                 names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'],
                                                 dir2save = directories[0])
                
    toc = perf_counter()
    fcBasics.printTime(tic, toc, 'Transform thickness data into 2D')

#%% GET THICKNESS REGION PLOTS (PROBALITY DENSITY ESTIMATION PLOTS)
    #   Finally, we can create distribution plots that will allow us to compare the cardiac thickness distribution between 
    #   regions of the heart (left vs right, dorsal vs ventral, atrium vs ventricle). 
    #   NOTE: This code can be modified to also get the distribution plots of the thickness of the other tissue layers. 
    #   Just let me know if you want me to include it! :)
    #   ================================================================================================================
    # cjPDF = False
    # if cjPDF:
    #     file_num = df_file[df_file['Folder']==folder].index.values[0]
    #     df_cjPDFs = fcAn.kdeThPlots(filename = filename, df_file = df_file, file_num = file_num, variable = 'cj_thickness', 
    #                                 thData = df_cjThNmyocIntBall, dir2save = directories[4], save = True)
        
    #     if save: 
    #         fcBasics.saveDF(filename = filename, df2save = df_cjPDFs, df_name = 'cjPDFs', dir2save = dir_results)
    #         fcBasics.saveDF(filename = filename, df2save = df_cjPDFs, df_name = 'cjPDFs', dir2save = os.path.join(dir_data2Analyse, 'R_All','df_cjPDFs'))
        
    # if save: 
    #     # Append all dicts to one object dict
    #     dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
    #                                              names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])
            
    # toc = perf_counter()
    # fcBasics.printTime(tic, toc, 'Transform thickness data into 2D')
    
    #%% func - kdeThPlots
# def kdeThPlots(filename, df_file, file_num, variable, thData, dir2save, save = True):
#     """
#     Function that creates Kernel Density Estimation Plots for each region of the heart (e.g, atrium, ventricle, 
#     left, right, dorsal and ventral and then

#     Parameters
#     ----------
#     filename : TYPE
#         DESCRIPTION.
#     df_file : TYPE
#         DESCRIPTION.
#     file_num : TYPE
#         DESCRIPTION.
#     variable : TYPE
#         DESCRIPTION.
#     thData : TYPE
#         DESCRIPTION.
#     dir2save : TYPE
#         DESCRIPTION.
#     save : TYPE, optional
#         DESCRIPTION. The default is True.

#     Returns
#     -------
#     df_pdfs : TYPE
#         DESCRIPTION.
        
#     Some links for reference:
#         #https://stackabuse.com/kernel-density-estimation-in-python-using-scikit-learn/
#         #https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
#         #https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html

#     """
#     if variable == 'cj_thickness':
#         title = filename + ' - Cardiac Jelly Thickness [um]'
#         xlabel = 'Cardiac Jelly Thickness [um]'
#         step = 0.05
#     elif variable == 'myoc_intBall' :
#         title = filename + ' - Myoc.Int Ballooning [um]'
#         xlabel = 'Myoc.Int Ballooning [um]'
#         step = 0.25
    
#     regions_div = [['AtrVent','atrium', 'ventricle'], ['LeftRight','left','right'],['DorsVent','dorsal','ventral']]
#     color = [['','tomato','gold'],['','deepskyblue', 'darkblue'],['','greenyellow', 'darkgreen']]
#     num_linspace = int((max(round(thData[variable]))/step)+2)
#     x_grid = np.linspace(0, max(round(thData[variable]))+step, num_linspace)

#     genotype = df_file.loc[file_num,'Gene_A']+':'+df_file.loc[file_num,'Genotype_A']
#     if df_file.loc[file_num,'Gene_B'] != '-':
#         genotype = str(genotype+'/'+df_file.loc[file_num,'Gene_B']+':'+df_file.loc[file_num,'Genotype_B'])
        
#     df_pdfs = pd.DataFrame(columns = ['Filename','Strain','Stage','Genotype', 'x_grid','AtrVent','atrium', 'ventricle', 'LeftRight','left','right','DorsVent','dorsal','ventral'])
#     df_pdfs['Filename'] = [filename for j in range(len(x_grid))]
#     df_pdfs['Strain'] = [df_file.loc[file_num,'Strain'] for j in range(len(x_grid))]
#     df_pdfs['Stage'] = [df_file.loc[file_num,'Stage'] for j in range(len(x_grid))]
#     df_pdfs['Genotype'] = [genotype for j in range(len(x_grid))]
#     df_pdfs['x_grid'] = x_grid
#     plot_dir = os.path.join(dir2save, filename+"_")
    
#     print('\n- Creating density plots... this process takes a while, about 8-12 min/plot out of 3 plots, be patient :)' )
#     bar = Bar('- Creating density plots', max=3, suffix = suffix, check_tty=False, hide_cursor=False)
#     for i, reg, col in zip(count(), regions_div, color):
#         df_one = thData[thData[reg[0]] == reg[1]]
#         th_one = df_one[variable]
#         bw_one = len(th_one)**(-1./(1+4))*np.std(th_one)
#         pdf_one = kde_sklearn(th_one, x_grid, bandwidth=bw_one)
        
#         df_two = thData[thData[reg[0]] == reg[2]]
#         th_two = df_two[variable]
#         bw_two = len(th_two)**(-1./(1+4))*np.std(th_two)
#         pdf_two = kde_sklearn(th_two, x_grid, bandwidth=bw_two)
        
#         print('\n - bandwidths:', format(bw_one,'.2f'), '-',format(bw_two, '.2f'))
#         pdf_comb = np.add(pdf_one, pdf_two)
        
#         pct_one = pdf_one/pdf_comb
#         # pct_two = pdf_two/pdf_comb
    
#         ones = np.ones((len(x_grid),))
        
#         fig, ax = plt.subplots(figsize=(8,5))
#         ax.fill_between(x_grid, pct_one, alpha=0.5, label= reg[1], color = col[1])
#         ax.fill_between(x_grid, pct_one, ones, alpha=0.5, label= reg[2], color= col[2])
#         # ax.plot(x_grid, pct_one, linewidth=1, alpha=0.5, label= reg[1])
#         # ax.plot(x_grid, pct_two, linewidth=1, alpha=0.5, label= reg[2])
#         ax.set_xlim(0, x_grid[-1])
#         ax.set_ylim(0, 1)        
#         box = ax.get_position()
#         ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#         ax.set_title(title, fontsize = 10)
#         ax.set_xlabel(xlabel, fontsize=10)
#         # Put a legend to the right of the current axis
#         ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#         if save: 
#             plt.savefig(plot_dir+"kde"+reg[0]+".png", dpi=300, bbox_inches='tight', transparent=True)
#         plt.show()

#         # fig, ax = plt.subplots()
#         # ax.plot(x_grid, pdf_one, linewidth=1, alpha=0.5, label= 'pdf_'+reg[1])
#         # ax.plot(x_grid, pdf_two, linewidth=1, alpha=0.5, label= 'pdf_'+reg[2])
#         # ax.legend(loc='upper right')
        
#         df_pdfs[reg[0]] = pct_one
#         df_pdfs[reg[1]] = pdf_one
#         df_pdfs[reg[2]] = pdf_two
        
#         bar.next()
        
#     bar.finish()
#     alert('wohoo',1)
    
#     return df_pdfs

#%% Init
init = True

    #%%
    # # Get centreline ribbon
    # cl_ribbon, kspl_ext, dict_kspl, dict_shapes, dict_planes = fcMeshes.createCLRibbon(filename = filename, kspl_CL2use = kspl_CL[0], linLine = linLines[0],
    #                                                              mesh = m_cjTh, dict_kspl = dict_kspl, dict_shapes = dict_shapes, dict_planes = dict_planes)
    
    # [m_cjThLnR] = fcMeshes.divideMeshesLnR(filename = filename, meshes = [m_cjTh], cl_ribbon = cl_ribbon)
    
    #%% TESTS
    # kspl_vSurf = fcMeshes.getExtCLonSurf(filename,  m_cjTh.alpha(0.05), kspl_ext, dict_planes['pl_Parallel2LinLine'])
    
    # kspl_CLnew = fcMeshes.getExtCLHighRes(filename, m_cjTh.alpha(0.05), kspl_ext, kspl_CL[0], dict_planes)
    
    # vent, atr, dict_pts, sph_cut = fcMeshes.ksplChamberCut(m_cjTh.alpha(0.05), kspl_CLnew, dict_shapes, dict_pts)