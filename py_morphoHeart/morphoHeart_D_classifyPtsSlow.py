# -*- coding: utf-8 -*-
"""
morphoHeart -  D. CLASSIFY POINTS SLOW Version (ATRIUM/VENTRICLE, DORSAL/VENTRAL, LEFT/RIGHT) 
@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os; import platform

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

suffix = '%(index)d/%(max)d - %(elapsed)ds'
save = True

#%% Importing morphoHeart packages
import morphoHeart_funcBasics as fcBasics 
import morphoHeart_funcMeshes as fcMeshes

#%% Get main directories (check which ones are actually used)
_, _, dir_lsOngoing, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
df_dataset = fcBasics.exportDatasetCSV(dir_lsOngoing, dir_data2Analyse)
# Get file to process and directories 
folder, df_file, file_num = fcBasics.selectFile(df_dataset); filename = folder[0:-3]
# directories = 0.dir_dict, 1.dir_txtNnpy, 2.dir_stl, 3.dir_cl, 4.dir_imsNvideos, 5.dir_ims2Analyse, 6.dir_folder2A
dir_results, directories = fcBasics.createDirectories2Save (filename, dir_data2Analyse, end_name = '2A')

#%% Load data for cjTh and myocIntBall
# Dictionary and decode it
[dict_cjThNmyocIntBall] = fcBasics.loadDicts(filename = filename, dicts_name = ['cjTh_myocIntBall2class'], directories = [directories[0]])
[AnV, DnV_Atr, DnV_Vent, pts_left, pts_whole, meas_param] = fcMeshes.decodeDict (dict2classify = dict_cjThNmyocIntBall, 
                                                                                 info = ['cj_thickness','myoc_intBall'])
# npy arrays
cj_thickness, myoc_intBall = meas_param

#%% Classify points
ptsClass_cjThNmyocIntBall = fcMeshes.classifyPts(AnV = AnV, DnV_Atr = DnV_Atr, DnV_Vent = DnV_Vent,
                                     pts_left =pts_left, pts_whole = pts_whole, info = 'cjTh_myocIntBall')

df_cjThNmyocIntBall = fcMeshes.ptsClass2df(pts_classFinal = ptsClass_cjThNmyocIntBall, meas_param = [cj_thickness,myoc_intBall], 
                                           names = ['cj_thickness','myoc_intBall'])

if save: 
    fcBasics.saveDF(filename = filename, df2save = df_cjThNmyocIntBall, df_name = 'df_cjThNmyocIntBall', dir2save = dir_results)

#%% Change folder name, meaning you finished processing the current heart
fcBasics.changeDirName(filename, directories[6])


#What was previously in morphoHeart_C_CutAndMeasure to then classify the points
if save: 
    # Classify Planes Sides
    AnV, DnV_Atr, DnV_Vent = fcMeshes.classifyPlanesSides(filename = filename, mesh = m_myoc, dict_planes = dict_planes)
    
    # CJ thickness - Int Myoc Ballooning ('Same external mesh)
    dict_cjThNmyocIntBall = fcMeshes.dict2classifyPts(AnV = AnV, DnV_Atr = DnV_Atr, DnV_Vent = DnV_Vent,
                                                   pts_left = m_cjThLnR.points(), pts_whole = m_cjTh.points(), 
                                                   meas_param = [cj_thickness,myoc_intBall], names = ['cj_thickness','myoc_intBall'],
                                                   info = 'cjTh_myocIntBall')
    fcMeshes.saveDict(filename = filename, dict2save = dict_cjThNmyocIntBall, name = 'cjTh_myocIntBall2class', dir2save = directories[0])
    # Ext Myoc Ballooning
    dict_extMyocBall2classify = fcMeshes.dict2classifyPts(AnV = AnV, DnV_Atr = DnV_Atr, DnV_Vent = DnV_Vent,
                                                   pts_left = m_myocExtLnR.points(), pts_whole = m_myocExtBall.points(), 
                                                   meas_param = [myoc_extBall], names = ['myoc_extBall'], info = 'myoc_extBall')
    fcMeshes.saveDict(filename = filename, dict2save = dict_extMyocBall2classify, name = 'myocExtBall2class', dir2save = directories[0])


#%% Functions

#%% func - classifyPlanesSides
def classifyPlanesSides(filename, mesh, dict_planes):
    
    # Classify marking points at either side of planes 
    # - AnV
    d_AnV, normal_AnV, centre_AnV, classSorted_AnV = getPlaneData2Classify(filename = filename, mesh = mesh, dict_planes = dict_planes, name = 'pl2CutMesh_Chamber', info = 'Atrium-Ventricle')
    AnV = [d_AnV, normal_AnV, classSorted_AnV]
    # - DnV_Atr
    d_DnV_Atr, normal_DnV_Atr, centre_DnV_Atr, classSorted_DnV_Atr = getPlaneData2Classify(filename = filename, mesh = mesh, dict_planes = dict_planes, name = 'pl_AtrCoronal', info = 'Dorsal-Ventral')
    DnV_Atr = [d_DnV_Atr, normal_DnV_Atr, classSorted_DnV_Atr]
    # - DnV_Vent
    d_DnV_Vent, normal_DnV_Vent, centre_DnV_Vent, classSorted_DnV_Vent = getPlaneData2Classify(filename = filename, mesh = mesh, dict_planes = dict_planes, name = 'pl_VentCoronal', info = 'Dorsal-Ventral')
    DnV_Vent = [d_DnV_Vent, normal_DnV_Vent, classSorted_DnV_Vent]
    
    return AnV, DnV_Atr, DnV_Vent

#%% func - dict2classifyPts
def dict2classifyPts(AnV, DnV_Atr, DnV_Vent, pts_left, pts_whole, meas_param, names, info):
    print('- Creating dict for '+info)
     # - AnV
    d_AnV, normal_AnV, classSorted_AnV = AnV
    # - DnV_Atr
    d_DnV_Atr, normal_DnV_Atr, classSorted_DnV_Atr = DnV_Atr
    # - DnV_Vent
    d_DnV_Vent, normal_DnV_Vent, classSorted_DnV_Vent = DnV_Vent
    
    dict2ClassifyPts = dict()
    
    # Classify marking points at either side of planes 
    # - AnV
    dictAnV = dict2ClassifyPts['AnV'] = dict()
    dictAnV['d'] = d_AnV
    dictAnV['normal'] = normal_AnV
    dictAnV['classSorted'] = classSorted_AnV
    
    # - DnV_Atr
    dictDnV_Atr = dict2ClassifyPts['DnV_Atr'] = dict()
    dictDnV_Atr['d'] = d_DnV_Atr
    dictDnV_Atr['normal'] = normal_DnV_Atr
    dictDnV_Atr['classSorted'] = classSorted_DnV_Atr
    
    # - DnV_Vent
    dictDnV_Vent = dict2ClassifyPts['DnV_Vent'] = dict()
    dictDnV_Vent['d'] = d_DnV_Vent
    dictDnV_Vent['normal'] = normal_DnV_Vent
    dictDnV_Vent['classSorted'] = classSorted_DnV_Vent

    # Save points in dict2ClassifyPts
    dict2ClassifyPts['pts_Left'] = pts_left
    dict2ClassifyPts['pts_Whole'] = pts_whole
    
    for i, param, name in zip(count(), meas_param, names):
        dict2ClassifyPts['param_'+name] = param
    
    return dict2ClassifyPts

#%% func - getPlaneData2Classify
def getPlaneData2Classify(filename, mesh, dict_planes, name, info):
    
    print('- Classifying sides of mesh as '+ info +'...')
    
    pl_normal = np.asarray(dict_planes[name]['pl_normal'])
    pl_centre = np.asarray(dict_planes[name]['pl_centre'])
    
    # Make normal of plane unitary 
    normal_unit = unit_vector(pl_normal)
    normalx10 = np.asarray([k*100 for k in normal_unit])
    
    # Create points at either side of the plane
    pt_plus = pl_centre+normalx10
    sph_plus = Sphere(pos=pt_plus, r=5, c='red')#.legend("pt(+)")
    pt_minus = pl_centre-normalx10
    sph_minus = Sphere(pos=pt_minus, r=5, c='green')#.legend("pt(-)")
    
    sph_centre = Sphere(pos=pl_centre, r=5, c='navy')#.legend("pt(o)")
    
    plane = Plane(pos = pl_centre, normal = pl_normal, sx = 300).legend(info).color('plum').alpha(1)
    
    info_split = info.split('-')
    text = filename+'\n\n >> Have a look at the position of the red and green spheres with respect to the \n     myocardium and classify them as "'+ info_split[0] + '" or "' + info_split[1] + '" \n     after closing the window'
    txt = Text2D(text, c="k", font= 'CallingCode')
    
    vp = Plotter(N=1, axes=4)
    vp.show(mesh, plane, sph_plus, sph_minus, sph_centre, txt, at=0, interactive=1)
    
    if info == "Atrium-Ventricle":
        classifier = ask4input("Select option -[0]:red=atrium/green=ventricle or [1]:red=ventricle/green=atrium-: ",int)
        if classifier == 0:
            class_pt_plus = "atrium"
            class_pt_minus = "ventricle"
        elif classifier == 1:
            class_pt_plus = "ventricle"
            class_pt_minus = "atrium"
            
    elif info == "Dorsal-Ventral":
        classifier = ask4input("Select option -[0]:red=dorsal/green=ventral or [1]:red=ventral/green=dorsal-: ",int)
        if classifier == 0:
            class_pt_plus = "dorsal"
            class_pt_minus = "ventral"
        elif classifier == 1:
            class_pt_plus = "ventral"
            class_pt_minus = "dorsal"
    
    # Get plane values
    d = normal_unit.dot(pl_centre)
    a, b, c = normal_unit
    #print("Normal_unit", normal_unit)
    #print("Plane d", d)
    
    x_plus,y_plus,z_plus = pt_plus
    dotProd_pt_plus = dot((a,b,c), (x_plus,y_plus,z_plus))
    #print("dotProd_pt_plus:", dotProd_pt_plus)
    x_minus,y_minus,z_minus = pt_minus
    dotProd_pt_minus = dot((a,b,c), (x_minus,y_minus,z_minus))
    #print("dotProd_pt_minus:", dotProd_pt_minus)
    
    dotProds = [dotProd_pt_minus, d, dotProd_pt_plus]
    classes = [class_pt_minus, 'center', class_pt_plus]
    #print (classes)
    classes_sorted = [a for _, a in sorted(zip(dotProds,classes))]
    #print(classes_sorted)
    
    return [d, dotProds], normal_unit, pl_centre, classes_sorted

#%% func - classifyPts
def classifyPts(AnV, DnV_Atr, DnV_Vent, pts_left, pts_whole, info):
    
    print('- Classifying points for ', info)
    tic = perf_counter()
    
    # Create empty lists
    pts_classLnR = []
    pts_classAtnVent = []
    pts_classDnV = []
    
    # Classify marking points at either side of planes 
    # - AnV
    d_AnV, normal_AnV, classSorted_AnV = AnV
    # - DnV_Atr
    d_DnV_Atr, normal_DnV_Atr, classSorted_DnV_Atr = DnV_Atr
    # - DnV_Vent
    d_DnV_Vent, normal_DnV_Vent, classSorted_DnV_Vent = DnV_Vent
    
    # Classify mesh points at either side of extended centreline (left/right)
    
    # - Get x, y, and z arrays of pts that make up the left side of the mesh
    x_pts_left = pts_left[:,0]
    y_pts_left = pts_left[:,1]
    z_pts_left = pts_left[:,2]
    
    bar = Bar('Classifying pts' , max = 100, suffix = suffix, check_tty=False, hide_cursor=False)
    num_pts = round(len(pts_whole)/100)
    
    for index, pt_whole in enumerate(pts_whole): 
        #print('index:', index)
        # look first for each individual coordinate
        x = pt_whole[0]
        where_x = np.where(x_pts_left == x)[0]
        y = pt_whole[1]
        where_y = np.where(y_pts_left == y)[0]
        z = pt_whole[2]
        where_z = np.where(z_pts_left == z)[0]
        # Unify information from the three coordinates 
        where_final = intersecOfSets(where_x, where_y, where_z)
       
        # Classify point as left/right -> If there is an intersection, point belongs to the left side, if not it belongs to the right
        if len(where_final) == 0:
            pts_classLnR.append('right')
        else: 
            pts_classLnR.append('left')
        
        # Classify point as atrium/ventricle
        pt_classAnV = classifyPtAtSidesofPlane(pt = [x,y,z], d = d_AnV, normal = normal_AnV, 
                                               classes_sorted = classSorted_AnV)
        pts_classAtnVent.append(pt_classAnV)
        
        # Classify point as dorsal/ventral using atrial or ventricular coronal plane
        if pt_classAnV == 'atrium':
            pt_classDnV = classifyPtAtSidesofPlane(pt = [x,y,z], d = d_DnV_Atr, normal = normal_DnV_Atr, 
                                                   classes_sorted = classSorted_DnV_Atr)
        elif pt_classAnV == 'ventricle':
            pt_classDnV = classifyPtAtSidesofPlane(pt = [x,y,z], d = d_DnV_Vent, normal = normal_DnV_Vent, 
                                                   classes_sorted = classSorted_DnV_Vent)
        pts_classDnV.append(pt_classDnV)
        
        if index % num_pts == 0:
            bar.next()
    
    bar.finish()
    toc = perf_counter()
    time = toc-tic
    
    alert('whistle',3)
    print("- All Done - points have been classified!")
    print("- Time taken to classify = ",format(time,'.2f'), "s/", format(time/60,'.2f'), "m/", format(time/3600,'.2f'), "h")     

    pts_classFinal = [pts_classLnR, pts_classAtnVent, pts_classDnV]
    
    return pts_classFinal

#%% func - ptsClass2df 
def ptsClass2df(pts_classFinal, meas_param, names):
    
    pts_classLnR, pts_classAtnVent, pts_classDnV = pts_classFinal
    
    df_ptsClass = {'AtrVent':pts_classAtnVent,'DorsVent': pts_classDnV, 'LeftRight': pts_classLnR}
    df_ptsClass = pd.DataFrame(df_ptsClass, columns = ['AtrVent','DorsVent','LeftRight'])
    
    order_list = ['AtrVent', 'DorsVent','LeftRight']
    for i, param, name in zip(count(),meas_param, names):
        df_ptsClass[name] = param 
        order_list.append(name)
    
    df_ptsClass = df_ptsClass[order_list]
    
    print(df_ptsClass.sample(10))
    
    return df_ptsClass