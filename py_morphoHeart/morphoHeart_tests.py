# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 11:22:39 2021

@author: mdp18js
"""
#%% Importing python packages
import os
from vedo import *
#from vedo import Plotter, Cube, settings, Text2D
from vedo import embedWindow
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

c="k"; font= 'VTK';
save = True; plot = True
azimuth = 0

#%% Start C_CutAndMeasure
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcContours as fcCont
    from morphoHeart_modules import morphoHeart_funcMeshes as fcMeshes

    #%% SELECT FILE AND GET METADATA
    #   This section allows the user to select file to process, get its main directories and metadata,
    #   and define some properties
    #   ================================================================================================================

    # Get main directories
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse)
    # Get file to process and directories
    folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]; dORv = filename[9:10]
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
    
#%% Load data
    # Import dictionaries
    [dict_obj, myoc_int_npcl] = fcBasics.loadDicts(filename = filename, dicts_name = ['dict_obj','myoc_int_npcl'],
                                                                    directories = [directories[0], directories[3]])
    [dict_planes, dict_pts, dict_kspl, dict_colour, dict_shapes] = fcMeshes.splitDicts(dict_obj)
    # Import meshes
   
    mTh_names = ['myoc_intBall']
    [m_thAll, colour_thAll] = fcMeshes.openThicknessMeshes(filename = filename, meshes_names = mTh_names, extension = 'vtk',
                                  dir_stl = directories[2], dir_txtNnpy = directories[1])
    m_myocIntBall = m_thAll[0]


    #%% Create ksplines, points, lines and centrelines
    #   This section will create spline(s) of the centreline(s) using the information from the loaded dictionary(ies) and
    #   return first a 3D interactive plot with the myocardium, endocardium and centreline(s) and second a new plot
    #   with two different visualisations of the maximum inscribed spheres with which the centreline(s) were calculated.
    #   Note: The resulting centreline(s) will comprise 300 points (important for next step)
    #   ================================================================================================================

    kSplinesCuts = fcMeshes.createKSpls(dict_kspl, kspl_list = ['ksplCut4CL_inflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Ext.Endo(Cut)', 'ksplCut4CL_outflow-Int.Myoc(Cut)', 'ksplCut4CL_inflow-Int.Myoc(Cut)'])
    sphCuts = fcMeshes.createSpheres(dict_pts, pts_list = ['sph_Cut4CL_inflow-Int.Myoc(Cut)', 'sph_Cut4CL_inflow-Ext.Endo(Cut)', 'sph_Cut4CL_outflow-Int.Myoc(Cut)', 'sph_Cut4CL_outflow-Ext.Endo(Cut)'])
    kspl_CL, linLines, sph_CL, sph_CL_colour, dict_shapes, dict_kspl = fcMeshes.createCLs(dict_cl = [myoc_int_npcl], dict_pts = dict_pts,
                                                                                          dict_kspl = dict_kspl, dict_shapes = dict_shapes,
                                                                                          colors = ['deepskyblue', 'tomato'])
    
    
    #%%
    vp = Plotter(N=1, axes = 13)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0], at = 0, interactive=True)
    
     # Divide heart layers into chambers and save data
    pl_Chambers, dict_planes, dict_pts, num_pt = fcMeshes.getInfo2CutChambers(filename = filename, kspl_CL = kspl_CL[0],
                                                                              mesh2cut = m_myocIntBall, dict_planes = dict_planes, dict_pts = dict_pts)

    pl_Chambers_normal = dict_planes['pl2CutMesh_Chamber']['pl_normal']
    pl_Chambers_centre = dict_planes['pl2CutMesh_Chamber']['pl_centre']
    
    pts2cut, data2cut = fcMeshes.getPointsAtPlane(points = m_myocIntBall.points(), pl_normal = pl_Chambers_normal,
                                        pl_centre = pl_Chambers_centre, tol = 1, addData = myoc_intBall)
    
    ordpts, angpts = fcMeshes.order_pts(points = pts2cut)
    cl_point = kspl_CL[0].points()[num_pt]
    sph_c = Sphere(pos = cl_point, c = 'darkorange', r = 2)
    sph_atr = Sphere(pos =  kspl_CL[0].points()[100], c = 'lime', r = 2)

    # Create spline around cut
    kspl = KSpline(ordpts, continuity=0, tension=0, bias=0, closed=True)
    kspl.color('magenta').legend('ksplCutChambers').lw(2)

    vp = Plotter(N=1, axes = 13)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0],kspl, sph_c, sph_atr, at = 0, interactive=True)
    
    import numpy as np
    colour_title = filename+"_myoc_intBall.npy"
    colour_dir = os.path.join(directories[1], colour_title)
    myoc_intBall = np.load(colour_dir)
        
    dist = []
    for pt in ordpts:
        dist.append(fcMeshes.findDist(cl_point, pt))
    
    r_circle_max = min(dist)*2
    r_circle_min = min(dist)*0.8
    normal_unit = fcMeshes.unit_vector(pl_Chambers_normal)*10
    
    step_rad = int(((r_circle_max-r_circle_min)/0.225)+1)
    for j, rad in enumerate(np.linspace(r_circle_min, r_circle_max, 10)):
        for i,h in enumerate(np.linspace(0.225/2,0.225*2,9)):
            cyl = Cylinder(pos = cl_point,r = rad, height = h, axis = normal_unit, c = 'lime', cap = True, res = 2000)#.wireframe(True)
            if i == 0 and j == 0:
                cyl_pts = cyl.points()
            else: 
                cyl_pts = np.concatenate((cyl_pts, cyl.points()))
    
    vp = Plotter(N=1, axes = 1)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0], cyl, sph_c, at = 0, interactive=True)
    
    # cyl_pts = cyl.points()
    cyl_points_rot = np.zeros_like(cyl_pts)
    for i, pt in enumerate(cyl_pts):
        cyl_points_rot[i] = (np.dot(fcMeshes.rotation_matrix(axis = [0,0,1], theta = np.radians(90)),pt))
    
    # sphs_cyl = Spheres(cyl_points_rot, r=0.1, c='red')
    # sph_origin = Sphere([0,0,0], r=3, c='black')
    
    # vp = Plotter(N=1, axes = 1)
    # vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0],kspl, cyl, sph_c, sphs_cyl, sph_origin, m_myocIntBall.clone().rotateZ(90), at = 0, interactive=True)
    
    cyl_pix = np.transpose(np.asarray([cyl_points_rot[:,i]//res[i] for i in range(len(res))]))
    cyl_pix = cyl_pix.astype(int)
    cyl_pix = np.unique(cyl_pix, axis =0)
    
    [s3_mask], _ = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1], end_name = ['ch0_all'], print_txt = False)
    
    xdim, ydim, zdim = s3_mask.shape
    s3_cyl = np.zeros_like(s3_mask)
    s3_cyl[cyl_pix[:,0],cyl_pix[:,1],cyl_pix[:,2]] = 1
    
    # s3_mask[cyl_pix[:,0],cyl_pix[:,1],cyl_pix[:,2]] = 0
    
    # fcCont.plt_s3(start_slc = 0, end_slc = s3_mask.shape[2]-1, im_every = 10,
    #                   s3_int = s3_mask, s3_ext = s3_cyl, plotshow = True, option = "both ch")
    
    import random
    
    # bar = Bar('Masking cut', max=zdim, suffix = suffix, check_tty=False, hide_cursor=False)
    for slc in range(zdim):
        im_cyl =s3_cyl[:,:,slc]
        pos_pts = np.where(im_cyl == 1)
        im = s3_mask[:,:,slc]
    
        clicks = [(pos_pts[0][i], pos_pts[1][i]) for i in range(pos_pts[0].shape[0])]
        clicks_random = random.sample(clicks, len(clicks))
        myIm = fcCont.drawLine(clicks_random, im, '0')
        s3_mask[:,:,slc] = myIm
    
    # fcCont.plt_s3(start_slc = 0, end_slc = s3_mask.shape[2]-1, im_every = 1,
                      # s3_int = s3_mask, s3_ext = s3_cyl, plotshow = True, option = "both ch")
    fcCont.plt_s3(start_slc = 155, end_slc = 322, im_every = 20,
                      s3_int = s3_mask, s3_ext = s3_cyl, plotshow = True, option = "both ch")
    
    test_mesh = fcMeshes.createLayerMesh(filename, s3_mask, res, layer='test', name='test', colour = 'cornflowerblue', alpha= 0.1, plotshow = True)
    
    # fcMeshes.saveMesh(filename, mesh= test_mesh, mesh_name= 'test_split', dir_stl = directories[2], extension = 'vtk')
    
    # vp = Plotter(N=1, axes = 1)
    # vp.show(test_mesh, at = 0, interactive=True)
    
    splitem = test_mesh.splitByConnectivity(maxdepth=2)
    
    vp = Plotter(N=1, axes = 1)
    vp.show(splitem, at = 0, interactive=True)
    
    splitem[0].isInside(point = kspl_CL[0].points()[100])
    splitem[1].isInside(point = kspl_CL[0].points()[100])
    
    splitem[1].legend('Atr')
    
    # pos_pts_new = np.where(myIm == 1)
    
    # crop = myIm[300:550, 300:400]
    
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(1,2, figsize = (8,8))
    # ax[0].imshow(crop)
    # ax[1].imshow(myIm)
    # plt.show()
    
    #%%
    # from vedo import *

    car = load(datadir+"porsche.ply").alpha(0.2)
    
    line = [(-9.,0.,0.), (0.,1.,0.), (9.,0.,0.)]
    tube = Tube(line).triangulate().c("violet",0.2)
    
    contour = car.intersectWith(tube).lineWidth(4).c('black')
    
    show(car, tube, contour, __doc__, axes=7)
    #%%
    cyls = []
    for j, h in enumerate(np.linspace(0.225/2,0.225*5,9)):
        cyl = Cylinder(pos = cl_point,r = r_circle_max, height = h, axis = normal_unit, c = 'lime', cap = True, res = 2000).wireframe(True)
        cyls.append(cyl)
        if i == 0 and j == 0:
            cyl_pts = cyl.points()
        else: 
            cyl_pts = np.concatenate((cyl_pts, cyl.points()))
    
    vp = Plotter(N=1, axes = 1)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0], cyls, sph_c, at = 0, interactive=True)
    
    cyl_points_rot = np.zeros_like(cyl_pts)
    for i, pt in enumerate(cyl_pts):
        cyl_points_rot[i] = (np.dot(fcMeshes.rotation_matrix(axis = [0,0,1], theta = np.radians(90)),pt))
    
    sph_origin = Sphere([0,0,0], r=3, c='black')
    sphs_cyl = Spheres(cyl_points_rot, r=0.1, c='red')
    
    vp = Plotter(N=1, axes = 1)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0],kspl, cyls,sph_c, sphs_cyl, sph_origin, m_myocIntBall.clone().rotateZ(90), at = 0, interactive=True)
    
    
    
    
    cyls = []
    for i, rad in enumerate(np.linspace(r_circle_min,r_circle_max,100)):
        for j, h in enumerate(np.linspace(0.225,0.225*5,9)):
            cyl = Cylinder(pos = cl_point,r = rad, height = h, axis = normal_unit, c = 'lime', cap = True, res = 100).wireframe(True)
            cyls.append(cyl)
            if i == 0 and j == 0:
                cyl_pts = cyl.points()
            else: 
                cyl_pts = np.concatenate((cyl_pts, cyl.points()))
        
    sph_c = Sphere(pos = cl_point, c = 'darkorange', r = 2)
    
    vp = Plotter(N=1, axes = 1)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0],kspl, cyls,sph_c, at = 0, interactive=True)
    
    cyl_points_rot = np.zeros_like(cyl_pts)
    for i, pt in enumerate(cyl_pts):
        cyl_points_rot[i] = (np.dot(fcMeshes.rotation_matrix(axis = [0,0,1], theta = np.radians(90)),pt))
    
    sph_origin = Sphere([0,0,0], r=3, c='black')
    sphs_cyl = Spheres(cyl_points_rot, r=0.1, c='red')
    
    vp = Plotter(N=1, axes = 1)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0],kspl, cyls,sph_c, sphs_cyl, sph_origin, m_myocIntBall.clone().rotateZ(90), at = 0, interactive=True)
    
    cyl_pix = np.transpose(np.asarray([cyl_points_rot[:,i]//res[i] for i in range(len(res))]))
    cyl_pix = cyl_pix.astype(int)
    cyl_pix = np.unique(cyl_pix, axis =0)
    
    [s3_mask], _ = fcCont.loadStacks(filename = filename, dir_txtNnpy = directories[1], end_name = ['ch0_all'], print_txt = False)
    
    xdim, ydim, zdim = s3_mask.shape
    s3_cyl = np.zeros_like(s3_mask)
    s3_cyl[cyl_pix[:,0],cyl_pix[:,1],cyl_pix[:,2]] = 1
    
    s3_mask[cyl_pix[:,0],cyl_pix[:,1],cyl_pix[:,2]] = 0
    
    fcCont.plt_s3(start_slc = 0, end_slc = s3_mask.shape[2]-1, im_every = 10,
                      s3_int = s3_mask, s3_ext = s3_cyl, plotshow = True, option = "both ch")
    
    #%%%
    #Create a ring in 0,0,0
    pts_cyl = np.zeros((2000,3))
    n = 0
    x = np.asarray([r_circle_max*cos(np.radians(theta)) for theta in np.linspace(0,360,360)])
    y = np.asarray([r_circle_max*sin(np.radians(theta)) for theta in np.linspace(0,360,360)])
    z = np.ones((360,))*3
    
    xyz_cyl = np.stack((x,y,z), axis = 1)
    
    xyz_sph = Spheres(xyz_cyl, r = 1, c = 'cyan')
    vp = Plotter(N=1, axes = 1)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0], cyls,sph_c, sph_origin, xyz_sph, at = 0, interactive=True)
    
    orig_normal = [0,0,1]
    rotX = np.dot([pl_Chambers_normal[2], pl_Chambers_normal[1]],[orig_normal[2],orig_normal[1]])
    rotY = np.dot([pl_Chambers_normal[2], pl_Chambers_normal[0]],[orig_normal[2],orig_normal[0]])
    rotZ = np.dot([pl_Chambers_normal[1], pl_Chambers_normal[0]],[orig_normal[1],orig_normal[0]])
    
    pts_cyl_rot = np.zeros_like(xyz_cyl)
    for i, pt in enumerate(xyz_cyl):
        pts_cyl_rot[i] = (np.dot(fcMeshes.rotation_matrix(axis = [1,0,0], theta = np.pi-rotX),pt))
    
    for i, pt in enumerate(pts_cyl_rot):
        pts_cyl_rot[i] = (np.dot(fcMeshes.rotation_matrix(axis = [0,1,0], theta = np.pi-rotY),pt))

    for i, pt in enumerate(pts_cyl_rot):
        pts_cyl_rot[i] = (np.dot(fcMeshes.rotation_matrix(axis = [0,0,1], theta = np.pi-rotZ),pt))
        
        
    xyz_sph_rot= Spheres(pts_cyl_rot, r = 2, c = 'gold')
    vp = Plotter(N=1, axes = 1)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0],sph_c, sph_origin, xyz_sph, xyz_sph_rot, at = 0, interactive=True)
    
    #%%%
    #find two points in the plane and use equation defined by https://math.stackexchange.com/questions/73237/parametric-equation-of-a-circle-in-3d-space/73242
    # to paint circle 
    
    import numpy as np
    pt_a = ordpts[np.where(angpts == min(angpts))[0]][0]
    # np.where((dists >= r) & (dists <= r + dr))
    pt_b = ordpts[0]#[np.where((angpts > -5) & (angpts < 5))[0][-1]]
    pt_c = cl_point
    
    sph_a = Sphere(pos = pt_a, c = 'darkblue', r = 2)
    sph_b = Sphere(pos = pt_b, c = 'red', r = 2)
    sph_c = Sphere(pos = pt_c, c = 'darkorange', r = 2)
    
    
    vp = Plotter(N=1, axes = 13)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0],kspl, sph_a, sph_b, sph_c, at = 0, interactive=True)
    
    ptu_a = fcMeshes.unit_vector(pt_a)
    ptu_b = fcMeshes.unit_vector(pt_a)
    ptu_c = fcMeshes.unit_vector(pt_c)
    
    sphu_a = Sphere(pos = ptu_a, c = 'darkblue', r = 0.5)
    sphu_b = Sphere(pos = ptu_b, c = 'red', r = 0.5)
    sphu_c = Sphere(pos = ptu_c, c = 'darkorange', r =0.5)
    
    vp = Plotter(N=1, axes = 13)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0],kspl, sph_a, sph_b, sph_c, sphu_a, sphu_b, sphu_c, at = 0, interactive=True)
    
    ring_pts = np.empty([len(range(0,360,5)), 3])
    for i, theta in enumerate( range(0,360,5)):
        theta_rad = np.radians(theta)
        ring_pts[i,0] = r_circle*np.cos(theta_rad) #pt_c[0] + r_circle*np.cos(theta_rad)*ptu_a[0] + r_circle*np.sin(theta_rad)*ptu_b[0]
        ring_pts[i,1] =  r_circle*np.sin(theta_rad) #pt_c[1] + r_circle*np.cos(theta_rad)*ptu_a[1] + r_circle*np.sin(theta_rad)*ptu_b[1]
        ring_pts[i,2] = 0 #pt_c[2] + r_circle*np.cos(theta_rad)*ptu_a[2] + r_circle*np.sin(theta_rad)*ptu_b[2]
        
    kspl_ring = KSpline(ring_pts,continuity=0, tension=0, bias=0, closed=True)
    kspl_ring.color('darkorange').legend('ksplRing').lw(2)
        
    vp = Plotter(N=1, axes = 1)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0],kspl, sph_a, sph_b, sph_c, kspl_ring, at = 0, interactive=True)    
    
    orig_normal = [0,0,1]
    rotX = np.dot([pl_Chambers_normal[1], pl_Chambers_normal[2]],[orig_normal[1],orig_normal[2]])
    rotY = np.dot([pl_Chambers_normal[0], pl_Chambers_normal[2]],[orig_normal[0],orig_normal[2]])
    rotZ = np.dot([pl_Chambers_normal[0], pl_Chambers_normal[1]],[orig_normal[0],orig_normal[1]])
    
    
    rot_ring = np.empty([len(range(0,360,5)),3])
    for i, pt in enumerate(ring_pts):
        pt_rotX = (np.dot(fcMeshes.rotation_matrix(axis = [1,0,0], theta = rotX), pt)) #+pt_c[0] 
        pt_rotY = (np.dot(fcMeshes.rotation_matrix(axis = [0,1,0], theta = rotY), pt_rotX))#+pt_c[1]
        rot_ring[i,:] = (np.dot(fcMeshes.rotation_matrix(axis = [0,0,1], theta = rotZ), pt_rotY))#+pt_c[2]
        
    kspl_ring_rot = KSpline(rot_ring,continuity=0, tension=0, bias=0, closed=True)
    kspl_ring_rot.color('red').legend('ksplRing_rot').lw(2)
    
    vp = Plotter(N=1, axes = 1)
    vp.show(m_myocIntBall.alpha(0.01), kspl_CL[0],kspl, sph_a, sph_b, sph_c, kspl_ring, kspl_ring_rot, at = 0, interactive=True) 
        
    
    
    
    s3_div =  fcMeshes.maskCutInChamberS3s (s3_mask = s3_mask, pl_normal = pl_Chambers_normal, pl_centre = pl_Chambers_centre, 
                                            r_min = r_circle, cl_pt = cl_point, resolution = res, tol = 2)
     
    mesh = fcMeshes.createLayerMesh(filename, s3_mask, res, layer='test', name = 'test', colour = 'cyan', alpha = 1, plotshow = True)
    
init = True
