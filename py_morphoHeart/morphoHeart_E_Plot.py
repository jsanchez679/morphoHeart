# -*- coding: utf-8 -*-
"""
morphoHeart - F. PLOT

@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import count
from vedo import *
from vedo import embedWindow#, settings, Plotter, Text2D
embedWindow(False)
settings.legendSize = .3
# settings.legendPos = 1
settings.legendFont="VTK"

init = False
# Verify working dir
def setWorkingDir (root_path, init):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    # init = True
    print("Current working directory: {0}".format(os.getcwd()))

    return root_path, init

root_path, init = setWorkingDir(os.getcwd(),init)

c="k"; font= 'VTK' # 'CallingCode'
azimuth = 0

#%% Start E_Plot
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    # from morphoHeart_modules import morphoHeart_funcContours as fcCont
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes

    #%% Get directories and file
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = '')

    #%% Plot things for just one heart
    # Get file to process and directories
    folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[2:]; dORv = filename[9:10]
    #stage = df_file.loc[file_num,'Stage']
    # directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse
    dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = 'R')

    # Import df_results
    df_res = fcBasics.loadDF(filename = filename, file = 'ResultsDF', dir_results = dir_results)
    file_num = df_res[df_res['Folder']== filename+'_2A'].index.values[0]
    if dORv == 'D':
        azimuth = -90
    elevation = df_res.loc[file_num,'ang_Heart']

    #%% Load data
    # Import dictionaries
    [dict_obj, myoc_int_npcl] = fcBasics.loadDicts(filename = filename, dicts_name = ['dict_obj','myoc_int_npcl'],
                                                                    directories = [directories[0], directories[3]])
    [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes] = fcMeshes.splitDicts(dict_obj)
    # Import meshes
    # m_names = ['myoc','myoc_atr','myoc_vent','endo','endo_atr','endo_vent','cj','cj_atr','cj_vent','cj_in']
    # m_all = fcMeshes.openMeshes(filename = filename, meshes_names = m_names, extension = 'vtk',
    #                             dir_stl = directories[2], alpha = [1]*len(m_names), dict_colour = dict_colour)
    # m_myoc, m_atrMyoc, m_ventMyoc, m_endo, m_atrEndo, m_ventEndo, m_cj, m_atrCJ, m_ventCJ, m_cjIn = m_all
    m_names = ['myoc','endo', 'cj']
    m_all = fcMeshes.openMeshes(filename = filename, meshes_names = m_names, extension = 'vtk',
                                dir_stl = directories[2], alpha = [1]*len(m_names), dict_colour = dict_colour)
    m_myoc, m_endo, m_cj = m_all
    # mTh_names = ['cj_thickness','myoc_thickness','endo_thickness','myoc_intBall']
    # [m_thAll, colour_thAll] = fcMeshes.openThicknessMeshes(filename = filename, meshes_names = mTh_names, extension = 'vtk',
    #                               dir_stl = directories[2], dir_txtNnpy = directories[1])
    # m_cjTh, m_myocTh, m_endoTh, m_myocIntBall = m_thAll
    mTh_names = ['cj_thickness']
    [m_thAll, colour_thAll] = fcMeshes.openThicknessMeshes(filename = filename, meshes_names = mTh_names, extension = 'vtk',
                                  dir_stl = directories[2], dir_txtNnpy = directories[1])
    m_cjTh = m_thAll[0]
    
    kSplinesCuts = fcMeshes.createKSpls(dict_kspl, kspl_list = ['ksplCut4CL_outflow-Int.Myoc(Cut)', 'ksplCut4CL_inflow-Int.Myoc(Cut)'])
    sphCuts = fcMeshes.createSpheres(dict_pts, pts_list = ['sph_Cut4CL_inflow-Int.Myoc(Cut)', 'sph_Cut4CL_outflow-Int.Myoc(Cut)'])
    kspl_CL, linLines, sph_CL, sph_CL_colour, dict_shapes, dict_kspl = fcMeshes.createCLs(dict_cl = [myoc_int_npcl], dict_pts = dict_pts,
                                                                                          dict_kspl = dict_kspl, dict_shapes = dict_shapes,
                                                                                          colors = ['tomato'])

    #%%
    # fcMeshes.saveMultVideos(filename = filename, info = ['myoc_endo','cj', 'cj_thickness'],
    #                         meshes4video = [[m_myoc.alpha(0.01), m_endo.alpha(0.01), kspl_CL[0]],[m_cj.alpha(0.01)],[m_cjTh]],
    #                         rotAngle= df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow = False, alpha_cube = 0)
    
    mTh_names = ['cj_thickness']
    [m_thAll, colour_thAll] = fcMeshes.openThicknessMeshes(filename = filename, meshes_names = mTh_names, extension = 'vtk',
                                  dir_stl = directories[2], dir_txtNnpy = directories[1])
    m_cjTh = m_thAll[0]
    
    dir2save = os.path.join(dir_data2Analyse, 'R_all', 'Videos', 'cube')
    [cj_thickness] = fcBasics.loadNPY(filename = filename, names = ['cj_thickness'], dir_txtNnpy = directories[1])
    m_cjTh.pointColors(cj_thickness, cmap="jet", vmin=0, vmax=25)
    m_cjTh.addScalarBar()
    m_cjTh.alpha(1)
    m_cjTh.mapper().SetScalarRange(0,25)

    vp = Plotter(N=1, axes=13)
    vp.show(m_cjTh, at=0, azimuth = azimuth, elevation = elevation, interactive=True)

    fcMeshes.saveVideo (filename = filename, info = 'cj_thickness0to25', meshes4video = [m_cjTh],
                        rotAngle  = df_res.loc[file_num,'ang_Heart'], dir2save = dir2save, alpha_cube = 0, plotshow=True)
    
    # fcMeshes.saveMultVideos(filename = filename, info = ['cj', 'cj_thickness'],
    #                         meshes4video = [[m_cj.alpha(0.01)],[m_cjTh]],
    #                         rotAngle= df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow = True, alpha_cube = 0)

#%%
    m_names = ['myoc_int']
    m_all = fcMeshes.openMeshes(filename = filename, meshes_names = m_names, extension = 'vtk',
                                dir_stl = directories[2], alpha = [1]*len(m_names), dict_colour = dict_colour)
    [m_myoc_int] = m_all

    vp = Plotter(N=1, axes=13)
    vp.show(m_myoc_int, at=0, interactive=True)

    fcMeshes.saveMesh(filename, m_myoc_int, 'myoc_int', directories[2], 'stl')
    #%%
    # plot_dir = os.path.join(directories[4], filename+"_")
    vp = Plotter(N=9, axes=13)
    text = filename + "\n\n >> Heart layers and chambers"; txt = Text2D(text, c=c, font=font)
    vp.show(m_myoc, txt, at=0)
    vp.show(m_atrMyoc, at=3)
    vp.show(m_endo, at=1)
    vp.show(m_atrEndo, at=4)
    vp.show(m_cj, at=2)
    vp.show(m_atrCJ, at=5)
    vp.show(m_ventMyoc, at=6)
    vp.show(m_ventEndo, at=7)
    vp.show(m_ventCJ, at=8,elevation = elevation, azimuth = azimuth, interactive=True)

    vp = Plotter(N=3, axes=13)
    text = filename + "\n\n >> Heart layers and chambers"; txt = Text2D(text, c=c, font=font)
    vp.show(m_myoc, txt, at=0)
    vp.show(m_endo, at=1)
    vp.show(m_cj, at=2, elevation = elevation, azimuth = azimuth, interactive=True).screenshot(filename='screenshot.png')

    #%%
    kSplinesCuts = fcMeshes.createKSpls(dict_kspl, kspl_list = ['ksplCut4CL_inflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Int.Myoc(Cut)', 'ksplCut4CL_inflow-Int.Myoc(Cut)'])
    sphCuts = fcMeshes.createSpheres(dict_pts, pts_list = ['sph_Cut4CL_inflow-Int.Myoc(Cut)', 'sph_Cut4CL_inflow-Ext.Endo(Cut)', 'sph_Cut4CL_outflow-Int.Myoc(Cut)', 'sph_Cut4CL_outflow-Ext.Endo(Cut)'])
    kspl_CL, linLines, sph_CL, sph_CL_colour, dict_shapes, dict_kspl = fcMeshes.createCLs(dict_cl = [myoc_int_npcl,endo_ext_npcl], dict_pts = dict_pts,
                                                                                          dict_kspl = dict_kspl, dict_shapes = dict_shapes,
                                                                                          colors = ['deepskyblue', 'tomato'])
    vp = Plotter(N=1, axes=10)
    vp.show(txt, kSplinesCuts, sphCuts, kspl_CL, linLines, m_myoc.alpha(0.01), m_endo.alpha(0.01), at=0, azimuth = azimuth, interactive=True)


    plane_Ch = Plane(pos=dict_planes['pl2CutMesh_Chamber']['pl_centre'],normal=dict_planes['pl2CutMesh_Chamber']['pl_normal'], sx=500)
    # Cut cl with plane
    ksplCL_cut = kspl_CL[0].clone().cutWithMesh(plane_Ch)
    # Find point of centreline closer to last point of kspline cut
    ksplCL_cutPt, num_pt = fcMeshes.findClosestPt(ksplCL_cut.points()[0], kspl_CL[0].points())

    #%%
    # df_res_0 = df_res
    sph_orient, lines_orient, dict_pts, dict_kspl, df_res = fcMeshes.getChambersOrientation(filename = filename, file_num = file_num, num_pt = num_pt, kspl_CL2use = kspl_CL[0],
                                                                                    myoc_meshes = [m_myoc, m_atrMyoc, m_ventMyoc], linLine = linLines[0],
                                                                                    dict_pts = dict_pts, dict_kspl = dict_kspl, df_res = df_res)
    print('ang_heart_original:',df_res_0.loc[file_num, 'ang_Heart'])
    print('ang_heart_new:',df_res.loc[file_num, 'ang_Heart'])

    #%%
    # plot_dir = os.path.join(directories[4], filename+"_")
    vp = Plotter(N=9, axes=13)
    text = filename + "\n\n >> Heart layers and chambers"; txt = Text2D(text, c=c, font=font)
    vp.show(m_myoc, txt, at=0)
    vp.show(m_atrMyoc, at=3)
    vp.show(m_endo, at=1)
    vp.show(m_atrEndo, at=4)
    vp.show(m_cj, at=2)
    vp.show(m_atrCJ, at=5)
    vp.show(m_ventMyoc, at=6)
    vp.show(m_ventEndo, at=7)
    vp.show(m_ventCJ, at=8, elevation = elevation, interactive=True)


    m_atrMyocC = m_ventMyoc
    m_atrMyocC.legend(m_atrMyoc._legend)
    m_ventMyocC = m_atrMyoc
    m_ventMyocC.legend('Ventricle(Myoc)')

    m_atrEndoC = m_ventEndo
    m_atrEndoC.legend(m_atrEndo._legend)
    m_ventEndoC = m_atrEndo
    m_ventEndoC.legend('Ventricle(Endo)')

    m_atrCJC = m_ventCJ
    m_atrCJC.legend(m_atrCJ._legend)
    m_ventCJC = m_atrCJ
    m_ventCJC.legend('Ventricle(CJ)')

    vp = Plotter(N=9, axes=13)
    text = filename + "\n\n >> Heart layers and chambers"; txt = Text2D(text, c=c, font=font)
    vp.show(m_myoc, txt, at=0)
    vp.show(m_atrMyocC, at=3)
    vp.show(m_endo, at=1)
    vp.show(m_atrEndoC, at=4)
    vp.show(m_cj, at=2)
    vp.show(m_atrCJC, at=5)
    vp.show(m_ventMyocC, at=6)
    vp.show(m_ventEndoC, at=7)
    vp.show(m_ventCJC, at=8, elevation = elevation, interactive=True)

    # screenshot(filename=plot_dir+'Heart_Layers_Chambers.png')

    # vp = Plotter(N=5, axes=13)
    # text = filename + "\n\n >> Thickness"; txt = Text2D(text, c=c, font=font)
    # vp.show(m_cjTh, m_cjIn.color('white'), txt, at=0)
    # vp.show(m_myocTh, at=1)
    # vp.show(m_endoTh, at=2)
    # vp.show(m_myocIntBall, at=3)
    # vp.show(m_myocExtBall, at = 4, elevation = elevation, interactive=True)
    # vp.screenshot(filename=plot_dir+'ThicknessPlot.png')

    dict_colour = fcMeshes.saveMeshes(filename = filename, meshes = [m_atrMyocC, m_atrEndoC, m_atrCJC, m_ventMyocC, m_ventEndoC, m_ventCJC],
                            names = ['myoc_atr', 'endo_atr', 'cj_atr', 'myoc_vent', 'endo_vent', 'cj_vent'],
                            dict_colour = dict_colour, dir_stl = directories[2], extension = 'vtk')

    #%% #CHECK names of elements in dict
    # text = filename + "\n\n >> Heart layers and chambers"; txt = Text2D(text, c=c, font=font)
    # kspl_all = fcMeshes.createKSpls(dict_kspl, kspl_list =['CL_Ext.Endo(Cut)', 'CL_Int.Myoc(Cut)','linLine_Ext.Endo(Cut)', 'linLine_Int.Myoc(Cut)'])#, 'CL_ext.Endo(Cut)'
    # cl_endo, cl_myoc, lin_endo, lin_myoc = kspl_all
    # clExt = fcMeshes.createKSpls(dict_kspl, kspl_list =['kspl_CLExtD','kspl_CLExtV'])

    # sph_all = fcMeshes.createSpheres(dict_pts, pts_list = ['sph_AtrCentre-Orient', 'sph_VentCentre-Orient', 'sph_InflowCentre-Orient',
    #                                    'sph_OutflowCentre-Orient', 'sph_AVCCentre-Orient', 'sph_ChamberCut'])
    # [sph_AtrC, sph_VenC, sph_infC, sph_outfC, sph_valve, sph_chCut] = sph_all
    # cl_ribbon = Ribbon(clExt[0], clExt[1], alpha=0.2, res=(500, 150)).wireframe(True).legend("rib_ExtCL(D-V)").color('purple')

    # vp = Plotter(N=1, axes=13)
    # vp.show(txt, cl_ribbon, kspl_all, sph_all, m_myoc.alpha(0.01), m_endo.alpha(0.01), at=0, azimuth = azimuth, elevation = elevation, interactive=True)


    #%% cj thickness fixed range
    dir2save = os.path.join(dir_data2Analyse, 'R_all', 'Videos', 'new')
    [cj_thickness] = fcBasics.loadNPY(filename = filename, names = ['cj_thickness'], dir_txtNnpy = directories[1])
    m_cjTh.pointColors(cj_thickness, cmap="jet", vmin=0, vmax=20)
    m_cjTh.addScalarBar()
    m_cjTh.alpha(1)
    m_cjIn.color("white").alpha(1).wireframe()
    m_cjTh.mapper().SetScalarRange(0,20)

    vp = Plotter(N=1, axes=13)
    vp.show(m_cjTh, at=0, azimuth = azimuth, elevation = elevation, interactive=True)

    fcMeshes.saveVideo (filename = filename, info = 'cj_thickness0to25', meshes4video = [m_cjTh],
                        rotAngle  = df_res.loc[file_num,'ang_Heart'], dir2save = dir2save, plotshow=True)
    
    # fcMeshes.saveVideo (filename = filename, info = 'heartAll', meshes4video = [m_myoc.alpha(0.01), m_endo.alpha(0.01)],
    #                     rotAngle  = df_res.loc[file_num,'ang_Heart'], dir2save = dir2save, plotshow=True)

#%%

    vp = Plotter(N=1, axes=13)
    vp.show(m_endo.alpha(0.5), m_cj.alpha(1), at=0, azimuth = azimuth, elevation = elevation, interactive=True)

    fcMeshes.saveVideo (filename = filename, info = 'heartEndoCJ', meshes4video = [m_endo.alpha(0.5), m_cj.alpha(1)],
                        rotAngle  = df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow=True)

    #%%
    kspl_CL, linLines, _, _, _, _ = fcMeshes.createCLs(dict_cl = [myoc_int_npcl,endo_ext_npcl], dict_pts = dict_pts,
                                                dict_kspl = dict_kspl, dict_shapes = dict_shapes,
                                                colors = ['deepskyblue', 'tomato'])
    sph_ballonning = fcMeshes.sphInSpline(kspl_CL = kspl_CL[0], name = 'sphs_ball', every = 0.6)
    [myoc_intBall] = fcBasics.loadNPY(filename = filename, names = ['myoc_intBall'], dir_txtNnpy = directories[1])
    m_myocIntBall.pointColors(myoc_intBall, cmap="jet", vmin=0, vmax=75)
    m_myocIntBall.addScalarBar()
    m_myocIntBall.alpha(1)
    # m_myocIntBall.color("white").alpha(1).wireframe()
    m_myocIntBall.mapper().SetScalarRange(0,75)

    # vp = Plotter(N=1, axes=13)
    # vp.show(sph_ballonning, m_myocIntBall, at=0, azimuth = azimuth, elevation = elevation, interactive=True)

    fcMeshes.saveVideo (filename = filename, info = 'myoc_intBall0to75', meshes4video = [m_myocIntBall],
                        rotAngle  = df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow=True)


    # fcMeshes.saveVideo(filename = filename, info = 'cj', meshes4video = [m_cj.alpha(0.01)],
    #                         rotAngle = df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow = True)

    # fcMeshes.saveVideo(filename = filename, info = 'myoc', meshes4video = [m_myoc],
    #                         rotAngle = df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow = False)
    # fcMeshes.saveVideo(filename = filename, info = 'endo', meshes4video = [m_endo],
    #                         rotAngle = df_res.loc[file_num,'ang_Heart'], dir2save = directories[4], plotshow = False)

    #%%
    # meshes4video = [kspl_all[0], kspl_all[1],kspl_all[2],kspl_all[3],
    #                 sph_all[0],sph_all[1],sph_all[2],sph_all[3],sph_all[4],sph_all[5],
    #                 m_myoc.alpha(0.01), m_endo.alpha(0.01)]
    # fcMeshes.saveVideo(filename = filename, info = 'heartAll', meshes4video = meshes4video,
    #                       rotAngle = df_res.loc[file_num,'ang_Heart'], dir2save = directories[4])

    
    


    #%%
    # fig, ax = plt.subplots(figsize=(8, 5))
    # # palette = sns.color_palette('jet')
    # g = sns.scatterplot(ax=ax, x='z_plane', y='theta', hue='cj_thickness', marker='o', data=df_filt, legend = False, cmap = 'jet')
    # # g.legend(bbox_to_anchor=(1, 1), ncol=1)
    # # g.set(xlim = (50000,250000))
    # ylabels = ['{:,.1f}'.format(y) for y in g.get_yticks()]
    # xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()/100]
    # g.set_xticklabels(xlabels)
    # g.set_yticklabels(ylabels)
    # # g.set(xlim = (0,1))
    # g.set(ylim = (-180,180))

    # #%%
    # #generate some psuedo data
    # df = pd.DataFrame({'num':[50000, 75000, 100000, 125000], 'Rent/Sqft':np.random.randn(4), 'Region':list('abcd')})

    # import matplotlib.pyplot as plt
    # import matplotlib.ticker as ticker
    # import seaborn as sns
    # import pandas as pd
    # sns.set(style="white")
    # fig, ax = plt.subplots(figsize=(8, 5))
    # palette = sns.color_palette("bright", 4)
    # g = sns.scatterplot(ax=ax, x="num", y="Rent/Sqft", hue="Region", marker='o', data=df, s=100, palette= palette)
    # g.legend(bbox_to_anchor=(1, 1), ncol=1)
    # g.set(xlim = (50000,250000))
    # xlabels = ['{:,.2f}'.format(x) + 'K' for x in g.get_xticks()/1000]
    # g.set_xticklabels(xlabels)

    # #%%
    # from matplotlib import pyplot as PLT
    # from matplotlib import cm as CM
    # from matplotlib import mlab as ML
    # import numpy as NP

    # n = 1e5
    # x = y = NP.linspace(-5, 5, 100)
    # X, Y = NP.meshgrid(x, y)
    # Z1 = ML.bivariate_normal(X, Y, 2, 2, 0, 0)
    # Z2 = ML.bivariate_normal(X, Y, 4, 1, 1, 1)
    # ZD = Z2 - Z1
    # x = X.ravel()
    # y = Y.ravel()
    # z = ZD.ravel()
    # gridsize=30
    # PLT.subplot(111)

    # # if 'bins=None', then color of each hexagon corresponds directly to its count
    # # 'C' is optional--it maps values to x-y coordinates; if 'C' is None (default) then
    # # the result is a pure 2D histogram

    # PLT.hexbin(x, y, C=z, gridsize=gridsize, cmap=CM.jet, bins=None)
    # PLT.axis([x.min(), x.max(), y.min(), y.max()])

    # cb = PLT.colorbar()
    # cb.set_label('mean value')
    # PLT.show()

    # #%%
    # import matplotlib.pyplot as plt
    # plt.style.use('seaborn-white')
    # import numpy as np

    # X = df_filt['z_plane']
    # Y = df_filt['theta']
    # Z = df_filt['cj_thickness']
    # plt.contourf(X, Y, Z, 20, cmap='RdGy')
    # plt.colorbar();

    # #%%
    # df = sns.load_dataset('iris')
    # ax = sns.kdeplot(x, y)
    # ax = sns.kdeplot(x, y, shade=True)
    # ax = sns.kdeplot(x, y, cbar=True)

    # iris = sns.load_dataset("iris")
    # >>> setosa = iris.loc[iris.species == "setosa"]
    # >>> virginica = iris.loc[iris.species == "virginica"]
    # >>> ax = sns.kdeplot(setosa.sepal_width, setosa.sepal_length,
    # ...                  cmap="Reds", shade=True, shade_lowest=False)
    # >>> ax = sns.kdeplot(virginica.sepal_width, virginica.sepal_length,
    # ...                  cmap="Blues", shade=True, shade_lowest=False)


    # grid_kws = {"height_ratios": (.9, .05), "hspace": .5}
    # f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
    # ax = sns.heatmap(heatmap, ax=ax,
    #                   cbar_ax=cbar_ax,
    #                   cbar_kws={"orientation": "horizontal"})
    # # f.figure.suptitle(filename, fontsize = 12)
    # plt.xlabel('Centreline position [Atrium >> Ventricle]', fontsize=10)
    # plt.savefig(plot_dir+x+"Unlooped.png", dpi=300, bbox_inches='tight')

    #%%
    # if dORv == 'V':
    #     azimuth = 0
    #     linLineX = lin_myoc.clone().projectOnPlane('x').c(lin_myoc.color()).x(0)
    # elif dORv == 'D': #Acaaa
    #     azimuth = -90
    #     linLineX = lin_myoc.clone().projectOnPlane('z').c(lin_myoc.color()).z(0)
    # ptsPl_linLine = Points([linLineX.points()[0], lin_myoc.points()[0], lin_myoc.points()[1]])
    # pl_linLine = fitPlane(ptsPl_linLine.points()).scale(4).c('mediumaquamarine').alpha(1).legend('pl_Parallel2LinLine')
    # pl_linLine_normal = pl_linLine.normal
    # pl_linLine_centre = pl_linLine.center
    # dict_planes = fcMeshes.addPlane2Dict (plane = pl_linLine, pl_centre = pl_linLine_centre,
    #                                         pl_normal = pl_linLine_normal, info = '', dict_planes = dict_planes)

    # dict_obj = fcMeshes.fillNsaveObjDict(filename = filename, dicts = [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes],
    #                                           names = ['dict_planes', 'dict_pts', 'dict_kspl', 'dict_colour', 'dict_shapes'], dir2save = directories[0])


    #%%
    # df_atr = df_classCJTh.loc[df_classCJTh.AtrVent == 'atrium', 'cj_thickness']
    # df_vent = df_classCJTh.loc[df_classCJTh.AtrVent == 'ventricle', 'cj_thickness']
    # df_left = df_classCJTh.loc[df_classCJTh.LeftRight == 'left', 'cj_thickness']
    # df_right = df_classCJTh.loc[df_classCJTh.LeftRight == 'right', 'cj_thickness']
    # df_dorsal = df_classCJTh.loc[df_classCJTh.DorsVent == 'dorsal', 'cj_thickness']
    # df_ventral = df_classCJTh.loc[df_classCJTh.DorsVent == 'ventral', 'cj_thickness']

    # # Plot
    # kwargs = dict(hist_kws={'alpha':.6}, kde_kws={'linewidth':2})

    # plt.figure(figsize=(10,7), dpi= 80)
    # sns.distplot(df_atr, color="dodgerblue", label="atrium", **kwargs)
    # sns.distplot(df_vent, color="springgreen", label="ventricle", **kwargs)
    # plt.xlim(-5,25)
    # plt.ylim(0,0.6)
    # # plt.suptitle(title, y=0.95, size=12)
    # plt.legend();

    # plt.figure(figsize=(10,7), dpi= 80)
    # sns.distplot(df_left, color="crimson", label="left", **kwargs)
    # sns.distplot(df_right, color="tomato", label="right", **kwargs)
    # plt.xlim(-5,25)
    # plt.ylim(0,0.6)
    # # plt.suptitle(title, y=0.95, size=12)
    # plt.legend();

    # plt.figure(figsize=(10,7), dpi= 80)
    # sns.distplot(df_dorsal, color="darkviolet", label="dorsal", **kwargs)
    # sns.distplot(df_ventral, color="deeppink", label="ventral", **kwargs)
    # plt.xlim(-5,25)
    # plt.ylim(0,0.6)
    # # plt.suptitle(title, y=0.95, size=12)
    # plt.legend();

    #%% Others
    # grid = sns.FacetGrid(df_classCJTh, row="AtrVent", col="LeftRight", hue = "DorsVent", margin_titles=True)
    # grid.map(plt.hist, "cj_thickness", bins=np.linspace(0.1, 20, 50));

    # dfs_ALD = df_classCJTh[(df_classCJTh['AtrVent'] == 'atrium') & (df_classCJTh['LeftRight'] == 'left')
    #                       & (df_classCJTh['DorsVent'] == 'dorsal')]['cj_thickness']
    # dfs_ALV = df_classCJTh[(df_classCJTh['AtrVent'] == 'atrium') & (df_classCJTh['LeftRight'] == 'left')
    #                       & (df_classCJTh['DorsVent'] == 'ventral')]['cj_thickness']
    # dfs_ARD = df_classCJTh[(df_classCJTh['AtrVent'] == 'atrium') & (df_classCJTh['LeftRight'] == 'right')
    #                       & (df_classCJTh['DorsVent'] == 'dorsal')]['cj_thickness']
    # dfs_ARV = df_classCJTh[(df_classCJTh['AtrVent'] == 'atrium') & (df_classCJTh['LeftRight'] == 'right')
    #                       & (df_classCJTh['DorsVent'] == 'ventral')]['cj_thickness']

    # dfs_VLD = df_classCJTh[(df_classCJTh['AtrVent'] == 'ventricle') & (df_classCJTh['LeftRight'] == 'left')
    #                       & (df_classCJTh['DorsVent'] == 'dorsal')]['cj_thickness']
    # dfs_VLV = df_classCJTh[(df_classCJTh['AtrVent'] == 'ventricle') & (df_classCJTh['LeftRight'] == 'left')
    #                       & (df_classCJTh['DorsVent'] == 'ventral')]['cj_thickness']
    # dfs_VRD =  df_classCJTh[(df_classCJTh['AtrVent'] == 'ventricle') & (df_classCJTh['LeftRight'] == 'right')
    #                       & (df_classCJTh['DorsVent'] == 'dorsal')]['cj_thickness']
    # dfs_VRV = df_classCJTh[(df_classCJTh['AtrVent'] == 'ventricle') & (df_classCJTh['LeftRight'] == 'right')
    #                       & (df_classCJTh['DorsVent'] == 'ventral')]['cj_thickness']

    # df_classALLs = [dfs_ALD, dfs_ALV, dfs_ARD, dfs_ARV, dfs_VLD, dfs_VLV, dfs_VRD, dfs_VRV]

    # for col in df_classALLs:
    #     sns.kdeplot(col, shade=True)

    # sns.kdeplot(
    #     data=tips, x="total_bill", hue="time",
    #     cumulative=True, common_norm=False, common_grid=True,
    # )

init = True