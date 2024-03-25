'''
morphoHeart_API

@author: Juliana Sanchez-Posada
'''
#%% Imports - ########################################################

#%% morphoHeart Imports - ##################################################
from .modules import mH_funcBasics as fcB
from .modules import mH_funcContours as fcC
from .modules import mH_funcMeshes as fcM
from .modules.mH_classes_new import *
from .gui.config import mH_config
from .gui.gui_classes import *

# Load logo
path_mHImages = mH_config.path_mHImages
path_logo = path_mHImages / 'logo-07.jpg'
logo = vedo.Picture(str(path_logo))

#%% API functions
def set_process(controller, running_process):
    proc_set = ['wf_info', 'active_process']
    controller.organ.update_settings(proc_set, running_process, 'mH')
    setattr(controller.main_win, 'running_process', running_process)

def mask_channel(controller, ch_name): 
    # Masking channel
    set_process(controller, 'masking_'+ch_name)
    # Workflow process
    workflow = controller.organ.workflow['morphoHeart']
    process = ['ImProc', ch_name, 'A-MaskChannel','Status']
    #Initial message
    controller.main_win.win_msg('Masking Channel '+str(ch_name[-1])+'!')
    #Enable and make visible close cont buttons
    enable_close_functions(controller=controller, process = 'mask', ch_name=ch_name)
    #Get channel
    im_ch = controller.organ.obj_imChannels[ch_name]
    controller.main_win.win_msg('Masking Channel '+str(ch_name[-1]))
    #Mask channel
    im_ch.maskIm()
    #Update Status in GUI and in CH Progress 
    status_btn = getattr(controller.main_win, 'mask_'+ch_name+'_status')
    if getattr(controller.main_win, 'mask_with_npy_'+ch_name).isChecked():
        if getattr(controller.main_win, 'set_npy_mask_'+ch_name).isChecked(): 
            wfn = controller.organ.workflow['morphoHeart']['ImProc'][ch_name]['A-MaskNPY']['Status']
            if wfn == 'DONE':
                controller.main_win.update_status(workflow, process, status_btn)
                controller.main_win.update_ch_progress()
            else: 
                controller.main_win.update_status(None, 'Initialised', status_btn, override=True)
        else: 
            controller.main_win.update_status(None, 'Initialised', status_btn, override=True)
    else: 
        controller.main_win.update_status(workflow, process, status_btn)
        controller.main_win.update_ch_progress()
    #Toggle button
    getattr(controller.main_win, 'mask_'+ch_name+'_play').setChecked(True)
    #Enable reset
    getattr(controller.main_win, 'reset_mask_'+ch_name).setEnabled(True)
    #Win msg 
    controller.main_win.win_msg('Channel '+str(ch_name[-1])+' has been masked!')
    #Plot masked stack
    controller.main_win.plot_all_slices(ch = ch_name)

def mask_npy_channel(controller, ch_name): 
    # Masking channel
    set_process(controller, 'masking_'+ch_name)
    # Workflow process
    workflow = controller.organ.workflow['morphoHeart']
    process = ['ImProc', ch_name, 'A-MaskNPY','Status']
    #Initial message
    controller.main_win.win_msg('Masking Channel '+str(ch_name[-1])+' with S3 from another Channel!')
    #Enable and make visible close cont buttons
    enable_close_functions(controller=controller, process = 'mask', ch_name=ch_name)
    #Get channel
    im_ch = controller.organ.obj_imChannels[ch_name]
    fcC.mask_with_npy(organ = controller.organ, 
                      im_ch=im_ch, gui_mask_npy = controller.main_win.gui_mask_npy[ch_name])
    #Update status
    status_btn = getattr(controller.main_win, 'mask_'+ch_name+'_status')
    if controller.organ.mH_settings['setup']['mask_ch'][ch_name]: 
        wfn = controller.organ.workflow['morphoHeart']['ImProc'][ch_name]['A-MaskChannel']['Status']
        if wfn != 'DONE': 
            controller.main_win.update_status(None, 'Initialised', status_btn, override=True)
            controller.main_win.update_ch_progress()
        else: 
            controller.main_win.update_status(workflow, process,status_btn)
    else: 
        controller.main_win.update_status(workflow, process, status_btn)
        cS_mask = getattr(controller.main_win, 'cS_'+ch_name+'_mask')
        controller.main_win.update_status(workflow, process, cS_mask)
    
    #Toggle button
    getattr(controller.main_win, 'npy_mask_'+ch_name+'_play').setChecked(True)
    #Enable reset
    getattr(controller.main_win, 'reset_npy_mask_'+ch_name).setEnabled(True)
    #Win msg 
    controller.main_win.win_msg('Channel '+str(ch_name[-1])+' has been masked with another Channel S3!')
    #Plot masked stack
    start = controller.main_win.gui_mask_npy[ch_name]['start_slc']
    end = controller.main_win.gui_mask_npy[ch_name]['end_slc']+1
    controller.main_win.plot_all_slices(ch = ch_name, slice_range=(start, end))

def autom_close_contours(controller, ch_name): 
    #Automatically close contours
    set_process(controller, 'autom_'+ch_name)
    #Uncheck DONE
    getattr(controller.main_win, 'autom_close_'+ch_name+'_done').setChecked(False)
    #Initial message    
    controller.main_win.win_msg('Automatic Closure of Contours has started for Channel '+str(ch_name[-1])+'!')
    #Enable and make visible close cont buttons
    enable_close_functions(controller=controller, process = 'autom', ch_name=ch_name)
    #Get channel
    im_ch = controller.organ.obj_imChannels[ch_name]
    # Automatically Close Contours
    im_ch.closeContours_auto(gui_param = controller.main_win.gui_autom_close_contours[ch_name], 
                              gui_plot = controller.main_win.plot_contours_settings[ch_name], 
                              win = controller.main_win)
    #Toggle button
    getattr(controller.main_win, 'autom_close_'+ch_name+'_play').setChecked(True)
    #Enable reset
    getattr(controller.main_win, 'reset_autom_close_'+ch_name).setEnabled(True)
    #Update organ workflow
    controller.main_win.user_done(process='autom_close', ch_name=ch_name)
    #Plot autom closed stack
    start = controller.main_win.gui_autom_close_contours[ch_name]['start_slc']
    end = controller.main_win.gui_autom_close_contours[ch_name]['end_slc']+1
    controller.main_win.plot_all_slices(ch = ch_name, slice_range=(start, end))

#Manually close contours
def manual_close_contours(controller, ch_name):
    #Close contours manually
    set_process(controller, 'manual_'+ch_name)
    #Uncheck DONE
    getattr(controller.main_win, 'manual_close_'+ch_name+'_done').setChecked(False)
    #Uncheck Saved
    getattr(controller.main_win, 'save_manually_closed_'+ch_name).setChecked(False)
    #Enable and make visible close cont buttons
    enable_close_functions(controller=controller, process = 'manual', ch_name=ch_name)
    #Get channel and images and save them as attribute
    controller.main_win.set_im_proc(ch_name)
    controller.main_win.slc_py = None
    controller.main_win.slc_list = None
    # Get slices
    slc_first_py = controller.main_win.gui_manual_close_contours[ch_name]['start_slc']
    slc_last_py = controller.main_win.gui_manual_close_contours[ch_name]['end_slc']
    if slc_last_py+1 > controller.main_win.im_proc[ch_name].shape[0]:
        pass
    else:
        slc_last_py = slc_last_py+1
        print('Last slice added 1')
    # Get Plot settings
    n_rows = int(controller.main_win.plot_contours_settings[ch_name]['n_rows'])
    n_cols = int(controller.main_win.plot_contours_settings[ch_name]['n_cols'])
    slcs_per_im = n_rows*n_cols
    #Get manual slices tuples
    manual_slices = fcC.set_tuples(slc_first=slc_first_py, slc_last=slc_last_py, slcs_per_im=slcs_per_im)
    controller.main_win.manual_slices = manual_slices
    #Message
    controller.main_win.win_msg('Manually closing contours for Channel '+str(ch_name[-1])+'. No. slices to close: '+str(slc_last_py-slc_first_py)+')')
    #Toggle button
    getattr(controller.main_win, 'manual_close_'+ch_name+'_play').setChecked(True)
    #Update organ workflow
    controller.main_win.user_done(process='manual_close', ch_name=ch_name)
    #Plot initial tuple
    controller.main_win.slc_tuple = controller.main_win.manual_slices[0]
    plot_tuple_to_manually_close(controller, ch_name, controller.main_win.slc_tuple, controller.main_win.im_proc[ch_name], initial=True)

def plot_tuple_to_manually_close(controller, ch_name, slc_tuple, im_proc, initial=True):
    #Get settings
    level = controller.main_win.gui_manual_close_contours[ch_name]['level']
    min_contour_len = controller.main_win.gui_manual_close_contours[ch_name]['min_contour_len']

    # > Fill widget data
    getattr(controller.main_win, 'tuple_analysed_'+ch_name).setText(str(slc_tuple[0]+1)+' - '+str(slc_tuple[1]))
    stack_cut = copy.deepcopy(im_proc[:][:][slc_tuple[0]:slc_tuple[1]+1])
    params = {'ch_name': ch_name, 'stack': stack_cut, 
                'slices_plot': slc_tuple, 'text': 'Contours', 
                'level': level, 'min_contour_length': min_contour_len}
    controller.main_win.add_thumbnail(function ='plot_slc_range', params = params, 
                        name = 'Cont Slcs '+str(slc_tuple[0]+1)+'-'+str(slc_tuple[1]-1+1)+'')
    controller.main_win.plot_slc_range(params)

    #Reset sly_py
    controller.main_win.slc_py = None
    controller.main_win.closing_slc.setText('')
    if initial: 
        #Enable buttons to close
        enable_close_functions(controller=controller, process='manual', ch_name=ch_name, widgets=False)
        #Enable widget
        getattr(controller.main_win, 'close_tuples_'+ch_name+'_widget').setEnabled(True)
        controller.main_win.win_msg('!Channel '+ch_name[-1]+': Select the slices you would like to close from tuple '+str(slc_tuple[0]+1)+'-'+str(slc_tuple[1])+'. [e.g: to close slices 5, 9, 10, and 11 type: "5,9-11"]')
        lineEdit = getattr(controller.main_win, 'slices_to_close_'+ch_name)
        lineEdit.clear()
        lineEdit.setFocus()
    else: 
        getattr(controller.main_win, 'prev_slice_'+ch_name).setEnabled(False)
        getattr(controller.main_win, 'next_slice_'+ch_name).setEnabled(False)
        controller.main_win.win_msg('!Channel '+ch_name[-1]+': Are you done closing the contours in this tuple '+str(slc_tuple[0]+1)+'-'+str(slc_tuple[1])+'? If so, continue closing contours in the next tuple. If not, type down the slices you would like to close.')

def close_slcs_tuple(controller, ch_name):
    main_win = controller.main_win
    main_win.slc_list = fcC.get_slices(lineEdit = getattr(main_win, 'slices_to_close_'+ch_name), 
                                                  slc_tuple=main_win.slc_tuple, 
                                                  win=main_win)
    getattr(main_win, 'prev_slice_'+ch_name).setEnabled(True)
    getattr(main_win, 'next_slice_'+ch_name).setEnabled(True)

    if len(main_win.slc_list) > 0: 
        print('slc_list:', main_win.slc_list)
        #Start with first slice
        main_win.win_msg('!Channel '+ch_name[-1]+': Manually closing contours of Slice '+str(main_win.slc_list[0]+1))
        # Get slice to close
        main_win.slc_py = main_win.slc_list[0]
        print('Initial slc_py:', main_win.slc_py)
        slc_user = main_win.slc_py+1
        # Get image from stack
        main_win.myIm = main_win.im_proc[ch_name][:][:][main_win.slc_py]
        #Get settings
        level = main_win.gui_manual_close_contours[ch_name]['level']
        min_contour_len = main_win.gui_manual_close_contours[ch_name]['min_contour_len']

        # Get contours for that slice
        params = {'myIm': copy.deepcopy(main_win.myIm), 'slc_user': slc_user, 'ch': ch_name, 
                'level': level, 'min_contour_length': min_contour_len}
        main_win.add_thumbnail(function ='fcC.plot_contours_slc', params = params, 
                                name = 'Cont Slc '+str(slc_user))
        main_win.plot_contours_slc(params)
        main_win.closing_slc.setText(str(slc_user))
        main_win.close_draw_btns_widget.setEnabled(True)
    else: 
        return

def next_prev_slice_in_tuple(next:bool, controller, ch_name:str):
    main_win = controller.main_win
    #Find the position of slc_py in slc_list
    index = main_win.slc_list.index(main_win.slc_py)
    if next: 
        new_index = index+1
        if new_index>= len(main_win.slc_list):
            plot_tuple_to_manually_close(controller, ch_name, controller.main_win.slc_tuple, 
                                         controller.main_win.im_proc[ch_name], initial=False)
            return
        else: 
            pass
    else: #previous
        new_index = index-1

    #Set new slc_py
    main_win.slc_py = main_win.slc_list[new_index]
    print('New slc_py:', main_win.slc_py)
    slc_user = main_win.slc_py+1
    # Get image from stack
    main_win.myIm = main_win.im_proc[ch_name][:][:][main_win.slc_py]
    main_win.win_msg('!Channel '+ch_name[-1]+': Manually closing contours of Slice '+str(slc_user))
    #Get settings
    level = main_win.gui_manual_close_contours[ch_name]['level']
    min_contour_len = main_win.gui_manual_close_contours[ch_name]['min_contour_len']

    # Get contours for that slice
    params = {'myIm': copy.deepcopy(main_win.myIm), 'slc_user': slc_user, 'ch': ch_name, 
            'level': level, 'min_contour_length': min_contour_len}
    main_win.add_thumbnail(function ='fcC.plot_contours_slc', params = params, 
                            name = 'Cont Slc '+str(slc_user))
    main_win.plot_contours_slc(params)
    main_win.closing_slc.setText(str(slc_user))

def next_prev_tuple_to_manually_close(next:bool, controller, ch_name):
    main_win = controller.main_win
    save_after_tuple = getattr(main_win, 'save_after_tuple_'+ch_name).isChecked()
    #Find the position of slc_tuple in manual_slices
    index = main_win.manual_slices.index(main_win.slc_tuple)
    if next: 
        new_index = index+1
        if new_index+1 > len(main_win.manual_slices):
            main_win.win_msg('!You have reached the end of the slice tuples. The closed channel will be saved and all the slices within the input range will be plotted for you to check!')
            main_win.save_closed_channel(ch=ch_name)
            slc_first_py = controller.main_win.gui_manual_close_contours[ch_name]['start_slc']
            slc_last_py = controller.main_win.gui_manual_close_contours[ch_name]['end_slc']
            main_win.plot_all_slices(ch = ch_name, slice_range = (slc_first_py, slc_last_py+1))
            return
        else:
            pass
    else: 
        new_index = index-1
    
    if save_after_tuple:
        main_win.win_msg('!Saving channel before continuing with next tuple.')
        main_win.save_closed_channel(ch=ch_name)
        #Check Save ch
        getattr(main_win, 'save_manually_closed_'+ch_name).setChecked(True)
    #Plot next tuple
    main_win.slc_tuple = main_win.manual_slices[new_index]
    print('New slc_tuple:', main_win.slc_tuple)
    plot_tuple_to_manually_close(controller, ch_name, main_win.slc_tuple, controller.main_win.im_proc[ch_name], initial=True)

#Select contours
def select_contours(controller, ch_name): 
    #Select contours
    set_process(controller, 'select_'+ch_name)
    #Uncheck DONE
    getattr(controller.main_win, 'select_contours_'+ch_name+'_done').setChecked(False)
    #Uncheck Saved S3
    getattr(controller.main_win, 'save_select_contours_'+ch_name).setChecked(False)
    #Get channel and save is as attribute
    controller.main_win.set_im_proc(ch_name, close=False)
    #Check if the channels3 already exist, and if not create empty channelS3 
    if len(controller.main_win.im_ch[ch_name].contStack)<3:
        print('No mask stacks have been created or saved for this channel') 
        controller.main_win.im_ch[ch_name].create_chS3s(layerDict={}, win=controller.main_win)
    #Load s3s and save them as attributes
    controller.main_win.im_ch[ch_name].load_chS3s(cont_types = ['int', 'tiss', 'ext'])
    if ch_name not in controller.main_win.s3_cont.keys(): 
        controller.main_win.s3_cont[ch_name] = {}
    for cont in ['int', 'tiss', 'ext']: 
        s3 = getattr(controller.main_win.im_ch[ch_name], 's3_'+cont).s3()
        controller.main_win.s3_cont[ch_name][cont] = s3
        # setattr(controller.main_win.s3_cont, 's3_'+cont, s3)
    #Enable and make visible close cont buttons
    enable_close_functions(controller=controller, process = 'selecting', ch_name=ch_name)
    #Create attribute to save selected contours 
    controller.main_win.dict_s3s = {}
    #Toggle button
    getattr(controller.main_win, 'select_contours_'+ch_name+'_play').setChecked(True)
    #########
    if getattr(controller.main_win, 'select_state_'+ch_name) != 'Finished': 
        #Get tuples_out
        gui_select_cont = controller.main_win.gui_select_contours[ch_name]
        controller.main_win.tuples_out = fcC.tuple_pairs(gui_select_cont)
        #Initialise progress bar
        controller.main_win.prog_bar_range(0,len(controller.main_win.tuples_out))
        #Message
        controller.main_win.win_msg('Selecting contours for Channel '+str(ch_name[-1])+'.')
        #Update organ workflow
        controller.main_win.user_done(process='select_contours', ch_name=ch_name)
        #Plot initial tuple
        at_least_one = False
        if getattr(controller.main_win, 'select_state_'+ch_name) == None: 
            set_state(controller, ch_name, 'Initialised')
            set_index_active(controller, ch_name, 0)
        else: 
            #Find the index of the tuple active and compare it the index_active saved
            print(getattr(controller.main_win, 'index_active_'+ch_name))
            controller.main_win.prog_bar_update(getattr(controller.main_win, 'index_active_'+ch_name))

        while not at_least_one: 
            try: 
                initial_index = controller.main_win.tuples_out[getattr(controller.main_win, 'index_active_'+ch_name)]
            except: 
                set_index_active(controller, ch_name, 0); print('Index_active_'+ch_name+' was reset to 0')
                initial_index = controller.main_win.tuples_out[getattr(controller.main_win, 'index_active_'+ch_name)]

            if initial_index['int_cont']+initial_index['ext_cont'] > 0:
                plot_props_to_select(controller, ch_name, getattr(controller.main_win, 'index_active_'+ch_name))
                at_least_one = True
                break
            else: 
                set_index_active(controller, ch_name, getattr(controller.main_win, 'index_active_'+ch_name)+1)
    else: 
        #Message
        controller.main_win.win_msg('!The contours have already been selected for Channel '+str(ch_name[-1])+'. Use the "Check / Modify Selection" module to check and modify the selected contours if needed.')

def set_state(controller, ch_name, select_state):
    proc_set = ['wf_info', 'select_contours', ch_name, 'select_state']
    controller.organ.update_settings(proc_set, select_state, 'mH')
    setattr(controller.main_win, 'select_state_'+ch_name, select_state)

def set_index_active(controller, ch_name, index_active):
    proc_set = ['wf_info', 'select_contours', ch_name, 'index_active']
    controller.organ.update_settings(proc_set, index_active, 'mH')
    setattr(controller.main_win, 'index_active_'+ch_name, index_active)

def plot_props_to_select(controller, ch_name, index_active):
    main_win = controller.main_win
    # > Fill widget data
    tuple_active = main_win.tuples_out[index_active] 

    first, last = tuple_active['tuple_pair']
    getattr(main_win, 'select_tuple_'+ch_name).setText(str(first+1)+' - '+str(last+1-1))
    getattr(main_win, 'select_fisrt_slc_'+ch_name).setText(str(first+1))
    getattr(main_win, 'int_of_'+ch_name).setText('/ '+str(tuple_active['int_cont']))
    getattr(main_win, 'ext_of_'+ch_name).setText('/ '+str(tuple_active['ext_cont']))
    main_win.slc_py = first

    print('Checking updates:', controller.organ.mH_settings['wf_info']['select_contours'][ch_name])
    #Get slice contours and plot contours 
    get_slc_cont_plot_props(controller=controller, ch_name=ch_name, tuple_active=tuple_active)
    
    #Enable widget
    getattr(main_win, 'select_contours_'+ch_name+'_widget').setEnabled(True)
    main_win.win_msg('!Channel '+ch_name[-1]+': Enter the internal and external contours for slice '+str(first+1)+". If no internal/external contours, leave the space empty.")
    lineEdit = getattr(main_win, 'int_cont_'+ch_name)
    lineEdit.clear()
    lineEdit.setFocus()

def get_slc_cont_plot_props(controller, ch_name, tuple_active):

    #Get settings
    level = controller.main_win.gui_select_contours[ch_name]['level']
    min_contour_len = controller.main_win.gui_select_contours[ch_name]['min_contour_len']

    #Get slice
    im_proc = controller.main_win.im_proc[ch_name]
    controller.main_win.myIm = im_proc[:][:][controller.main_win.slc_py]
    slc_user = controller.main_win.slc_py+1
    # Show new closed contours
    new_contours = get_contours(controller.main_win.myIm, min_contour_length = min_contour_len, level = level)
    # Sort the contours by length (bigger to smaller)
    controller.main_win.new_contours = sorted(new_contours, key = len, reverse=True)
    params_props = {'myIm': copy.deepcopy(controller.main_win.myIm), 'ch': ch_name, 'slc':slc_user, 
                    'cont_sort': controller.main_win.new_contours, 'tuple_active': tuple_active, 
                    'win':controller.main_win}
    fcC.plot_props(params_props)
    controller.main_win.add_thumbnail(function='fcC.plot_props', params = params_props, 
                                        name='Conts. Slc'+str(slc_user))

def select_slcs_tuple(controller, ch_name):
    main_win = controller.main_win
    #Uncheck SAVE s3
    getattr(main_win, 'save_select_contours_'+ch_name).setChecked(False)
    #Get contour numbers
    tuples_out_slc = controller.main_win.tuples_out[getattr(main_win, 'index_active_'+ch_name)]
    num_contours = fcC.get_contour_num(lineEdit_int = getattr(main_win, 'int_cont_'+ch_name),
                                        lineEdit_ext = getattr(main_win, 'ext_cont_'+ch_name),
                                        tuples_out_slc = tuples_out_slc, 
                                        num_contours = len(main_win.new_contours), 
                                        win=main_win)
    
    if num_contours != None: 
        #Diable Play
        getattr(main_win, 'slc_tuple_select_'+ch_name+'_play').setEnabled(False)
        #Get tuples 
        tuple_range = main_win.tuples_out[getattr(main_win, 'index_active_'+ch_name)]['tuple_pair']
        getattr(main_win, 'progress_select_'+ch_name).setRange(tuple_range[0], tuple_range[1])
        main_win.win_msg('!Selecting contours for Slices '+str(tuple_range[0]+1)+'-'+str(tuple_range[1])+'...')
        
        input_numbers = [len(num_contours[key])> 0 for key in num_contours.keys()]
        #Use this if to check if the right number of contours where entered
        if any(flag == True for flag in input_numbers): 
            #Get the s3s with the selected contours
            selected_cont, all_cont, s3s = s3_with_contours(controller = controller, 
                                                            ch_name = ch_name,
                                                            num_contours=num_contours, 
                                                            contours=main_win.new_contours, 
                                                            slc=main_win.slc_py+1)
            #Plot that image with filled contours
            # params_filled = {'myIm': copy.deepcopy(main_win.myIm), 'slc':main_win.slc_py+1, 
            #                 'ch': ch_name, 's3s': s3s, 'win': main_win, 'all_cont': all_cont}
            # fcC.plot_filled_contours(params_filled)
            # main_win.add_thumbnail(function='fcC.plot_filled_contours', params = params_filled, 
            #                         name='FilledCont. Slc'+str(main_win.slc_py+1))

            # print('selected_cont:', selected_cont)
            #Get number of contours
            no_int = main_win.tuples_out[getattr(main_win, 'index_active_'+ch_name)]['int_cont']
            no_ext = main_win.tuples_out[getattr(main_win, 'index_active_'+ch_name)]['ext_cont']
            level = controller.main_win.gui_select_contours[ch_name]['level']
            min_contour_len = controller.main_win.gui_select_contours[ch_name]['min_contour_len']

            dict_plot = []
            dict_plot.append({'myIm': copy.deepcopy(main_win.myIm), 'slc':main_win.slc_py+1, 
                            'ch': ch_name, 's3s': s3s, 'all_cont': all_cont})
        
            slc_o = tuple_range[0]+1
            getattr(main_win, 'progress_select_'+ch_name).setValue(slc_o)
            if len(range(slc_o, tuple_range[1], 1)) > 0: 
                for slc in range(slc_o, tuple_range[1], 1):
                    print('slc:', slc, 'slc_user: ', slc+1)
                    main_win.myIm = main_win.im_proc[ch_name][:][:][slc]
                    #Get contours, sort them and get their properties
                    contours = get_contours(myIm=main_win.myIm, min_contour_length = min_contour_len, level = level)
                    contours = sorted(contours, key = len, reverse=True)
                    props_all = fcC.get_cont_props(myIm = main_win.myIm, cont_sort=contours)
                    ################
                    # #Plot contours (TO DELETE)
                    # params_props = {'myIm': copy.deepcopy(myIm), 'ch': ch_name, 'slc':slc+1, 
                    #                     'cont_sort': contours, 'win':controller.main_win}
                    # fcC.plot_props(params_props)
                    # controller.main_win.add_thumbnail(function='fcC.plot_props', params = params_props, 
                    #                                     name='Conts. Slc'+str(slc+1))
                    ################

                    #Automatically select the contours
                    selected_out = fcC.autom_select_contours(props_first = selected_cont,
                                                            props_myIm = props_all)
                    
                    selected_outf = fcC.confirm_selection(selected_out=selected_out,
                                                        props_myIm = props_all,
                                                        num_conts = {'ext': no_ext, 'int': no_int}, 
                                                        slc = slc)
                    selected_cont = None
                    #Get the s3s with the automatically selected contours
                    selected_cont, all_out, s3s_out = s3_with_contours(controller = controller,
                                                                       ch_name = ch_name,
                                                                        num_contours=selected_outf,
                                                                        contours=contours,  
                                                                        slc = slc+1)
                    #Plot that image with filled contours
                    # params_filled_out = {'myIm': copy.deepcopy(myIm), 'slc':slc+1, 
                    #             'ch': ch_name, 's3s': s3s_out, 'win': main_win, 'all_cont': all_out}
                    # fcC.plot_filled_contours(params_filled_out)
                    # main_win.add_thumbnail(function='fcC.plot_filled_contours', params = params_filled_out, 
                    #                         name='FilledCont. Slc'+str(slc+1))

                    #Make a plot for 12 slc to minimise slcs shown 
                    if len(dict_plot)== 12: 
                        dict_plot = []
                        slc_o = slc+1

                    params_slc = {'myIm': copy.deepcopy(main_win.myIm), 'slc':slc+1, 
                                'ch': ch_name, 's3s': s3s_out, 'all_cont': all_out}
                    dict_plot.append(params_slc)

                    if len(dict_plot)== 12 or slc == tuple_range[1]-1: 
                        slc_f = slc
                        params_group = {'win': main_win, 'dict_plot': dict_plot}
                        fcC.plot_group_filled_contours(params = params_group)
                        main_win.add_thumbnail(function='fcC.plot_group_filled_contours', params = params_group, 
                                            name='Filled Slcs'+str(slc_o)+'-'+str(slc_f+1))
                        
                    getattr(main_win, 'progress_select_'+ch_name).setValue(slc)
                
            else: 
                params_group = {'win': main_win, 'dict_plot': dict_plot}
                fcC.plot_group_filled_contours(params = params_group)
                main_win.add_thumbnail(function='fcC.plot_group_filled_contours', params = params_group, 
                                    name='Filled Slcs'+str(slc_o))

            alert('bubble')
            main_win.win_msg('!Contours have been selected for Slices '+str(tuple_range[0]+1)+'-'+str(tuple_range[1])+'. To continue selecting contours go to the next slice group.')
            getattr(main_win, 'progress_select_'+ch_name).setValue(tuple_range[1])
    
    getattr(main_win, 'next_group_'+ch_name).setEnabled(True)
    main_win.prog_bar_update(getattr(main_win, 'index_active_'+ch_name)+1)
    
def next_tuple_select(controller, ch_name): 
    main_win = controller.main_win
    getattr(main_win, 'next_group_'+ch_name).setEnabled(False)
    getattr(main_win, 'slc_tuple_select_'+ch_name+'_play').setEnabled(True)
    save_after_tuple = getattr(main_win, 'save_after_group_'+ch_name).isChecked()
    #Get the next tuple
    at_least_one = False
    getattr(main_win, 'int_cont_'+ch_name).clear()
    getattr(main_win, 'ext_cont_'+ch_name).clear()
    while not at_least_one: 
        new_index = getattr(main_win, 'index_active_'+ch_name) + 1
        if new_index not in main_win.tuples_out.keys(): 
            main_win.win_msg('!You have reached the end of the stack.')
            getattr(main_win, 'tick_modify_select_'+ch_name).setEnabled(True)
            getattr(main_win, 'progress_select_'+ch_name).setValue(100)
            main_win.prog_bar_update(100)
            set_state(controller, ch_name, 'Finished')
            return
        else: 
            if save_after_tuple: 
                main_win.win_msg('!Saving channel before continuing with next group.')
                main_win.save_closed_channel(ch=main_win.im_ch.channel_no, s3s=True)
                #Check SAVE s3
                getattr(main_win, 'save_select_contours_'+ch_name).setChecked(True)

            set_index_active(controller, ch_name, new_index)
            next_index = main_win.tuples_out[getattr(main_win, 'index_active_'+ch_name)]
            if next_index['int_cont']+next_index['ext_cont'] > 0:
                #Plot next tuple
                plot_props_to_select(controller, ch_name, getattr(main_win, 'index_active_'+ch_name))
                at_least_one = True
                getattr(main_win, 'progress_select_'+ch_name).setValue(len(main_win.tuples_out))
                break
            else: 
                set_index_active(controller, ch_name, getattr(main_win, 'index_active_'+ch_name)+1)

def s3_with_contours(controller, ch_name, num_contours, contours, slc):
    
    main_win = controller.main_win
    selected_cont = {}; all_cont = {'contours':[]}; s3s = {}
    for ctype in ['internal', 'external']: 
        selected_cont[ctype] = {'contours': [], 'props': []}
        for cont in num_contours[ctype]: 
            sp_cont = contours[cont]
            selected_cont[ctype]['contours'].append(sp_cont)
            all_cont['contours'].append(sp_cont)
            # Get properties of the contour
            sp_props = fcC.get_cont_props(main_win.myIm, [sp_cont])
            selected_cont[ctype]['props'].append(sp_props[0])
        
        slc_s3 = main_win.s3_cont[ch_name][ctype[:3]][:,:,slc] #getattr(main_win, 's3_'+ctype[:3])[:,:,slc]
        if len(num_contours[ctype])>0: 
            slc_s3 = fcC.fill_contours(selected_cont[ctype], slc_s3)
        else: 
            slc_s3 = np.zeros_like(slc_s3, dtype='bool').astype(int)
        main_win.s3_cont[ch_name][ctype[:3]][:,:,slc] = slc_s3 
        #getattr(main_win, 's3_'+ctype[:3])[:,:,slc] = slc_s3
        s3s[ctype[:3]] = slc_s3
    
    slc_tiss = main_win.s3_cont[ch_name]['tiss'][:,:,slc] #getattr(main_win, 's3_tiss')[:,:,slc]
    if len(all_cont['contours'])>0: 
        slc_tiss = fcC.fill_contours(all_cont, slc_tiss)
    else: 
        slc_tiss = np.zeros_like(slc_tiss, dtype='bool').astype(int)
    main_win.s3_cont[ch_name]['tiss'][:,:,slc] = slc_tiss
    # getattr(main_win, 's3_tiss')[:,:,slc] = slc_tiss
    s3s['tiss'] = slc_tiss
    main_win.dict_s3s[slc-1] = selected_cont

    return selected_cont, all_cont, s3s

def modify_selected_contours(controller, ch_name): 
    main_win = controller.main_win
    lineEdit = getattr(controller.main_win, 'man_int_cont_'+ch_name)
    lineEdit.clear()
    lineEdit.setFocus()
    main_win.slc_to_select = fcC.get_slices(lineEdit = getattr(main_win, 'select_manually_slcs_'+ch_name), 
                                            slc_tuple=(1,int(getattr(main_win, 'total_stack_slices_'+ch_name).text())), 
                                            win=main_win)
    #Get tuples_out
    gui_select_cont = main_win.gui_select_contours[ch_name]
    main_win.tuples_out = fcC.tuple_pairs(gui_select_cont)
    
    #Initialise progress bar
    main_win.prog_bar_range(0,len(main_win.slc_to_select))
    main_win.win_msg('Selecting contours for Channel '+str(ch_name[-1])+'.')
    main_win.index_slc_active = 0

    getattr(main_win, 'slc_manually_select_'+ch_name+'_play').setEnabled(True)
    getattr(main_win, 'next_slc_select_'+ch_name).setEnabled(False)
    plot_props_to_slc_select(controller, ch_name, main_win.index_slc_active)
    
def plot_props_to_slc_select(controller, ch_name, index_active): 
    
    #Slc active
    slc_active = controller.main_win.slc_to_select[index_active]
    controller.main_win.slc_py = slc_active
    #Find the tuple in which this slc fits 
    controller.main_win.tuple_active = fcC.find_slc_within_tuples(slc_active, controller.main_win.tuples_out)
    #Set the values in gui
    getattr(controller.main_win, 'select_slc_'+ch_name).setText(str(slc_active+1))
    getattr(controller.main_win, 'man_int_of_'+ch_name).setText('/ '+str(controller.main_win.tuple_active['int_cont']))
    getattr(controller.main_win, 'man_ext_of_'+ch_name).setText('/ '+str(controller.main_win.tuple_active['ext_cont']))
    #Get slice contours and plot contours 
    get_slc_cont_plot_props(controller=controller, ch_name=ch_name, tuple_active=controller.main_win.tuple_active)
    controller.main_win.win_msg('!Channel '+ch_name[-1]+': Enter the internal and external contours for slice '+str(slc_active+1)+". If no internal/external contours, leave the space empty.")
    lineEdit = getattr(controller.main_win, 'man_int_cont_'+ch_name)
    lineEdit.clear()
    lineEdit.setFocus()

def select_slc(controller, ch_name):
    main_win = controller.main_win
    #Get contour numbers
    num_contours = fcC.get_contour_num(lineEdit_int = getattr(main_win, 'man_int_cont_'+ch_name),
                                        lineEdit_ext = getattr(main_win, 'man_ext_cont_'+ch_name),
                                        tuples_out_slc = main_win.tuple_active, 
                                        num_contours = len(main_win.new_contours), 
                                        win=main_win, ignore=True)
    
    if num_contours != None: 
        #Diable Play
        getattr(main_win, 'slc_manually_select_'+ch_name+'_play').setEnabled(False)
        main_win.win_msg('!Selecting contours for Slice '+str(controller.main_win.slc_py+1)+'...')
        
        #Get the s3s with the selected contours
        selected_cont, all_cont, s3s = s3_with_contours(controller = controller, 
                                                        ch_name = ch_name,
                                                        num_contours=num_contours, 
                                                        contours=main_win.new_contours, 
                                                        slc=main_win.slc_py+1)
        # Plot that image with filled contours
        params_filled = {'myIm': copy.deepcopy(main_win.myIm), 'slc':main_win.slc_py+1, 
                        'ch': ch_name, 's3s': s3s, 'win': main_win, 'all_cont': all_cont}
        fcC.plot_filled_contours(params_filled)
        main_win.add_thumbnail(function='fcC.plot_filled_contours', params = params_filled, 
                                name='FilledCont. Slc'+str(main_win.slc_py+1))
        main_win.prog_bar_update(controller.main_win.index_slc_active+1)
        getattr(main_win, 'next_slc_select_'+ch_name).setEnabled(True)
        
def next_slc_select(controller, ch_name):
    main_win = controller.main_win
    getattr(main_win, 'next_slc_select_'+ch_name).setEnabled(False)
    getattr(main_win, 'slc_manually_select_'+ch_name+'_play').setEnabled(True)
    getattr(main_win, 'man_int_cont_'+ch_name).clear()
    getattr(main_win, 'man_ext_cont_'+ch_name).clear()
    #Get the next slice
    at_least_one = False
    while not at_least_one: 
        new_index = main_win.index_slc_active + 1
        if new_index >= len(main_win.slc_to_select): 
            main_win.win_msg('!You have reached the end of the slices.')
            getattr(main_win, 'progress_select_'+ch_name).setValue(100)
            main_win.prog_bar_update(100)
            return
        else: 
            #Plot next slice
            main_win.index_slc_active = new_index
            plot_props_to_slc_select(controller, ch_name, main_win.index_slc_active)
            at_least_one = True
            main_win.prog_bar_update(new_index)
            break

def enable_close_functions(controller, process, ch_name, widgets=True): 

    scrollArea = getattr(controller.main_win, 'scrollArea_'+ch_name)
    ignore=False
    if process == 'mask': 
        widget_central = getattr(controller.main_win, 'mask_'+ch_name+'_widget')
    elif process == 'autom': 
        widget_central = getattr(controller.main_win, 'autom_close_'+ch_name+'_widget')
    elif process == 'manual': 
        widget_central = getattr(controller.main_win, 'manual_close_'+ch_name+'_widget')
        if widgets: 
            #Show and enable
            controller.main_win.close_draw_btns_widget.setVisible(True)
            controller.main_win.close_draw_btns_widget.setEnabled(False)
            getattr(controller.main_win, 'next_slice_'+ch_name).setShortcut("Ctrl+M")
            getattr(controller.main_win, 'prev_slice_'+ch_name).setShortcut("Ctrl+N")
            #Open Buttons Section
            controller.main_win.functions_btns_open.setChecked(False)
            controller.main_win.open_section(name = 'functions_btns')
            #Enable buttons in manually close subsection
            getattr(controller.main_win, 'save_after_tuple_'+ch_name).setEnabled(True)
            getattr(controller.main_win, 'save_manually_closed_'+ch_name).setEnabled(True)
            getattr(controller.main_win, 'manual_close_'+ch_name+'_done').setEnabled(True)
            #UnCheck DONE
            getattr(controller.main_win, 'manual_close_'+ch_name+'_done').setChecked(False)
        else:
            #Enable buttons to close
            controller.main_win.close_draw_btns_widget.setEnabled(True)
    elif process == 'selecting':
        widget_central = getattr(controller.main_win, 'select_contours_all_'+ch_name+'_widget')
        controller.main_win.close_draw_btns_widget.setVisible(False)
        getattr(controller.main_win, 'close_tuples_'+ch_name+'_widget').setEnabled(False)
        if widgets: 
            #Enable
            getattr(controller.main_win, 'select_contours_'+ch_name+'_widget').setEnabled(True)
            getattr(controller.main_win, 'select_plot_slc_'+ch_name).setEnabled(True)
            getattr(controller.main_win, 'select_plot_all_'+ch_name).setEnabled(True)
            getattr(controller.main_win, 'save_select_contours_'+ch_name).setEnabled(True)
            getattr(controller.main_win, 'select_contours_'+ch_name+'_done').setEnabled(True)
            #UnCheck DONE
            getattr(controller.main_win, 'select_contours_'+ch_name+'_done').setChecked(False)
            getattr(controller.main_win, 'plot_final_s3s_'+ch_name).setEnabled(False)
            getattr(controller.main_win, 'plot_final_s3s_'+ch_name).setVisible(False)
            if getattr(controller.main_win, 'select_state_'+ch_name) != 'Finished': 
                getattr(controller.main_win, 'next_group_'+ch_name).setShortcut("Ctrl+Shift+.")
                #Close Buttons Section
                controller.main_win.functions_btns_open.setChecked(True)
                controller.main_win.open_section(name = 'functions_btns')
                #Enable buttons in selecting subsection
                getattr(controller.main_win, 'save_after_group_'+ch_name).setEnabled(True)
                controller.main_win.closing_slc.setText('')
                #Diable 
                getattr(controller.main_win, 'next_group_'+ch_name).setEnabled(False)
            else: 
                # activate modify option
                getattr(controller.main_win, 'manually_select_widget_'+ch_name).setEnabled(True)
                getattr(controller.main_win, 'tick_modify_select_'+ch_name).setEnabled(True)
                getattr(controller.main_win, 'tick_modify_select_'+ch_name).setChecked(True)
                getattr(controller.main_win, 'slc_manually_select_'+ch_name+'_play').setEnabled(False)
                controller.main_win.enable_modify(ch_name)
        scrollArea.verticalScrollBar().setValue(scrollArea.verticalScrollBar().maximum())
        ignore=True
    else: 
        widget_central = getattr(controller.main_win,'plot_slices_'+ch_name+'_widget')
    if not ignore:
        scrollArea.ensureWidgetVisible(widget_central)
    
    #Shortcuts not working...
    print('')

#ANALYSIS TAB
def run_keeplargest(controller):
    workflow = controller.organ.workflow
    set_process(controller, 'keeplargest')
    #Check channels have already been created: 
    if len(controller.organ.obj_imChannels.keys()) == controller.organ.mH_settings['setup']['no_chs']:
        fcM.s32Meshes(organ = controller.organ, gui_keep_largest=controller.main_win.gui_keep_largest, 
                    win = controller.main_win, rotateZ_90=controller.organ.mH_settings['setup']['rotateZ_90'])

        #Update Status in GUI
        controller.main_win.update_status(None, 'DONE', controller.main_win.keeplargest_status, override = True)
        #Enable button for plot all
        getattr(controller.main_win, 'keeplargest_plot').setEnabled(True)
        #Toggle button
        getattr(controller.main_win, 'keeplargest_play').setChecked(True)
        #Update progress in main_win
        controller.main_win.update_workflow_progress()
        #Add meshes to plot_user
        controller.main_win.fill_comboBox_all_meshes()
        
    else: 
        title = 'Channels not closed / Contours not selected!'
        msg = 'You are not done closing/selecting the contours of the input channels! \nPlease go back to  -mH: Segment Channels-  Tab and continue processing the channels before running the processes in this tab'
        prompt = Prompt_ok(title, msg, parent=controller.welcome_win)
        prompt.exec()
        print('output:', prompt.output)
        return

    print('\nEND Keeplargest')
    print('organ.mH_settings:', controller.organ.mH_settings)
    print('organ.workflow:', workflow)

def run_cleanup(controller):
    set_process(controller, 'cleanup')
    if controller.main_win.keeplargest_play.isChecked(): 
        workflow = controller.organ.workflow
        if controller.main_win.gui_clean['plot2d']:
            plot_settings = (True, controller.main_win.gui_clean['n_slices'])
        else: 
            plot_settings = (False, None) 

        fcM.clean_ch(organ = controller.organ, 
                    gui_clean = controller.main_win.gui_clean, 
                    win=controller.main_win, plot_settings=plot_settings)

        #Enable button for plot all
        getattr(controller.main_win, 'clean_plot').setEnabled(True)
        #Toggle button
        getattr(controller.main_win, 'cleanup_play').setChecked(True)
        #Update progress in main_win
        controller.main_win.update_workflow_progress()
        #Add meshes to plot_user
        controller.main_win.fill_comboBox_all_meshes()

    else: 
        controller.main_win.win_msg('*To clean-up the tissues make sure you have at least run the  -Keep Largest-  section.')
            
def run_trimming(controller):
    set_process(controller, 'trimming')
    if controller.main_win.keeplargest_play.isChecked(): 
        workflow = controller.organ.workflow
        meshes, no_cut, cuts_out = get_trimming_planes(organ = controller.organ,
                                                        gui_trim = controller.main_win.gui_trim,
                                                        win = controller.main_win)
        
        fcM.trim_top_bottom_S3s(organ = controller.organ, meshes = meshes, 
                                no_cut = no_cut, cuts_out = cuts_out,
                                win = controller.main_win)
        #Enable button for plot all
        getattr(controller.main_win, 'trimming_plot').setEnabled(True)
        #Toggle button
        getattr(controller.main_win, 'trimming_play').setChecked(True)
        #Update progress in main_win
        controller.main_win.update_workflow_progress()
        #Add meshes to plot_user
        controller.main_win.fill_comboBox_all_meshes()
    else: 
        controller.main_win.win_msg('*To trim the tissues make sure you have at least run the  -Keep Largest-  section.')

def get_trimming_planes(organ, gui_trim, win): 
    filename = organ.user_organName
    #Get meshes to cut
    meshes = []
    no_cut = []
    settings = {'color': {}, 'name':{}}
    aa = 0
    for ch in organ.obj_imChannels.keys():
        for cont in ['tiss', 'ext', 'int']:
            if gui_trim['top']['chs'][ch][cont] or gui_trim['bottom']['chs'][ch][cont]:
                meshes.append(organ.obj_meshes[ch+'_'+cont])
                settings['color'][aa] = organ.obj_meshes[ch+'_'+cont].color
                settings['name'][aa] = organ.obj_meshes[ch+'_'+cont].legend
                aa+=1
                break
            else: 
                no_cut.append(ch+'_'+cont)
    print('settings:', settings)

    # User user input to select which meshes need to be cut
    cuts_names = {'top': {'heart_def': 'outflow tract','other': 'top'},
                'bottom': {'heart_def': 'inflow tract','other': 'bottom'}}
    cuts_out = copy.deepcopy(gui_trim)

    cut_top = []; cut_bott = []; #cut_chs = {}
    cuts_flat = flatdict.FlatDict(gui_trim)
    print('A:',cuts_flat)
    for key in cuts_flat.keys():
        if 'top' in key and 'object' not in key: 
            cut_top.append(cuts_flat[key])
        if 'bot' in key and 'object' not in key: 
            cut_bott.append(cuts_flat[key])
                                
    # print('cut_chs:', cut_chs)
    print('cut_top:', cut_top)
    print('cut_bott:', cut_bott)
            
    if mH_config.heart_default:
        name_dict =  'heart_def'     
    else: 
        name_dict = 'other'

    #Define plane to cut bottom
    if any(cut_bott):
        happy = False
        #Define plane to cut bottom
        while not happy: 
            plane_bott, pl_dict_bott = fcM.get_plane(filename=filename, 
                                                txt = 'cut '+cuts_names['bottom'][name_dict],
                                                meshes = meshes, settings=settings)#, win=win)  
            title = 'Happy with the defined plane?' 
            msg = 'Are you happy with the defined plane to cut '+cuts_names['bottom'][name_dict]+'?'
            items = {0: {'opt':'no, I would like to define a new plane.'}, 1: {'opt':'yes, continue!'}}
            prompt = Prompt_ok_cancel_radio(title, msg, items, parent=win)
            prompt.exec()
            print('output:', prompt.output, '\n')  
            if prompt.output[0] == 1: 
                happy = True

        cuts_out['bottom']['plane_info_mesh'] = pl_dict_bott
        # Reorient plane to images (s3)
        plane_bottIm, pl_dict_bottIm = fcM.rotate_plane2im(pl_dict_bott['pl_centre'], 
                                                            pl_dict_bott['pl_normal'])
        cuts_out['bottom']['plane_info_image'] = pl_dict_bottIm
        
    #Define plane to cut top
    if any(cut_top):
        happy = False
        #Define plane to cut top
        while not happy: 
            plane_top, pl_dict_top = fcM.get_plane(filename=filename, 
                                                txt = 'cut '+cuts_names['top'][name_dict],
                                                meshes = meshes, settings=settings)#, win=win)

            title = 'Happy with the defined plane?' 
            msg = 'Are you happy with the defined plane to cut '+cuts_names['top'][name_dict]+'?'
            items = {0: {'opt':'no, I would like to define a new plane.'}, 1: {'opt':'yes, continue!'}}
            prompt = Prompt_ok_cancel_radio(title, msg, items, parent=win)
            prompt.exec()
            print('output:', prompt.output, '\n')  
            if prompt.output[0] == 1: 
                happy = True
        
        cuts_out['top']['plane_info_mesh'] = pl_dict_top
        # Reorient plane to images (s3)
        plane_topIm, pl_dict_topIm = fcM.rotate_plane2im(pl_dict_top['pl_centre'], 
                                                        pl_dict_top['pl_normal'])
        cuts_out['top']['plane_info_image'] = pl_dict_topIm
        
    print('cuts_out:', cuts_out)
    return meshes, no_cut, cuts_out

def run_axis_orientation(controller, only_roi=False):
    set_process(controller, 'orientation')
    if controller.main_win.keeplargest_play.isChecked(): 
        #CHECK HERE WHICH ONE HAS BEEN RUN AND IF STACK WAS RUN, 
        # PROGRAM SO THAT THE STACK DATA DOESN'T GET REMOVED 
        workflow = controller.organ.workflow['morphoHeart']
        if not only_roi: 
            fcM.get_stack_orientation(organ = controller.organ,  
                                    gui_orientation = controller.main_win.gui_orientation, 
                                    win = controller.main_win)
        
        try: 
            print('tryA')
            gui_orientation = controller.main_win.gui_orientation
        except: 
            print('exceptA')
            gui_orientation = controller.organ.mH_settings['wf_info']['orientation']

        on_hold = fcM.get_roi_orientation(organ = controller.organ,
                                            gui_orientation = gui_orientation,  
                                            win = controller.main_win)
        
        #Update Status in GUI
        process = ['MeshesProc', 'A-Set_Orientation', 'Status']
        controller.main_win.update_status(workflow, process, controller.main_win.orient_status)

        #Toggle button
        print('organ.on_hold:', on_hold)
        select_btn = getattr(controller.main_win, 'orientation_play')
        if not on_hold:
            controller.organ.on_hold = on_hold
            select_btn.setChecked(True)
        else: 
            controller.organ.on_hold = on_hold
            select_btn.setChecked(False)

        #Update progress in main_win
        controller.main_win.update_workflow_progress()
    else: 
        controller.main_win.win_msg('*To set the stack and roi axes make sure you have at least run the  -Keep Largest-  section.')

def run_chNS(controller):
    set_process(controller, 'chNS')
    if controller.main_win.keeplargest_play.isChecked(): 
        if controller.main_win.gui_chNS['plot2d']:
            plot_settings = (True, controller.main_win.gui_chNS['n_slices'])
        else: 
            plot_settings = (False, None) 

        fcM.extract_chNS(organ = controller.organ, 
                        rotateZ_90 = controller.organ.mH_settings['setup']['rotateZ_90'],
                        win = controller.main_win, 
                        plot_settings = plot_settings)

        #Toggle button
        getattr(controller.main_win, 'chNS_play').setChecked(True)
        getattr(controller.main_win, 'summary_whole_plot_chNS').setEnabled(True)

        #Update progress in main_win
        controller.main_win.update_workflow_progress()

        #Add meshes to plot_user
        controller.main_win.fill_comboBox_all_meshes()

    else: 
        controller.main_win.win_msg('*To extract the channel from the negative space make sure you have at least run the  -Keep Largest-  section.')

# > CENTRELINE
def run_centreline_clean(controller):
    set_process(controller, 'centreline')
    if controller.main_win.keeplargest_play.isChecked(): 
        workflow = controller.organ.workflow['morphoHeart']

        m4clf = fcM.proc_meshes4cl(controller.organ, 
                                    win=controller.main_win) 
                                                                    
        if not hasattr(controller.organ, 'obj_temp'):
            controller.organ.obj_temp = {}
        controller.organ.obj_temp['centreline'] = {'SimplifyMesh': m4clf}
        print('obj_temp: ', controller.organ.obj_temp)

        #Enable button for ML
        getattr(controller.main_win, 'centreline_ML_play').setEnabled(True)

        #Toggle button
        getattr(controller.main_win, 'centreline_clean_play').setChecked(True)

        #Update Status in GUI
        process = process =  ['MeshesProc','C-Centreline','Status']
        controller.main_win.update_status(workflow, process, controller.main_win.centreline_status)

        prompt_meshLab(controller)

        #Update progress in main_win
        controller.main_win.update_workflow_progress()

    else: 
        controller.main_win.win_msg('*To extract the centreline from tissues make sure you have at least run the  -Keep Largest-  section.')

def run_centreline_ML(controller): 
    set_process(controller, 'centreline')
    workflow = controller.organ.workflow['morphoHeart']
    proc_simp = ['MeshesProc', 'C-Centreline', 'SimplifyMesh', 'Status']
    if get_by_path(workflow, proc_simp) == 'DONE':
        #Check all the MeshLab meshes have been created
        cl_names = list(controller.organ.mH_settings['measure']['CL'].keys())
        all_saved = []
        for nn, cl in enumerate(cl_names): 
            ch, cont, _ = cl.split('_')
            #Get name and path
            directory = controller.organ.dir_res(dir ='centreline')
            name_ML = controller.organ.mH_settings['wf_info']['centreline']['dirs'][ch][cont]['dir_meshLabMesh']
            mesh_dir = directory / name_ML
            if mesh_dir.is_file():
                controller.main_win.win_msg('MeshLab cut mesh ' +str(name_ML)+' was found!')
                all_saved.append(True)

                # Update organ workflow
                status_sq = getattr(controller.main_win, 'meshLab_status'+str(nn+1))
                controller.main_win.update_status(None, 'DONE', status_sq, override=True)

                #Enable button for vmtk
                getattr(controller.main_win, 'centreline_vmtk_play').setEnabled(True)

                #Toggle button
                getattr(controller.main_win, 'centreline_ML_play').setChecked(True)

                #Update progress in main_win
                controller.main_win.update_workflow_progress()

            else:
                error = '*'+str(name_ML)+' has not been created! Clean this mesh in MeshLab to proceed.'
                controller.main_win.win_msg(error)
                msg_add = [str(name_ML)+' was not found in the centreline folder. Make sure you have named your cleaned meshes correctly after running the processing in Meshlab.', 
                            'To clean up the meshes with MeshLab follow the next steps:']
                prompt_meshLab(controller, msg_add=msg_add)
                controller.main_win.centreline_ML_play.setChecked(False)
                return
            
            controller.main_win.win_msg('All MeshLab cut Meshes were successfully found!')

def run_centreline_vmtk(controller): 
    set_process(controller, 'centreline')
    if controller.main_win.centreline_ML_play.isChecked(): 
        workflow = controller.organ.workflow['morphoHeart']
        try: 
            print('try')
            voronoi = controller.main_win.gui_centreline['vmtk_CL']['voronoi']
        except:
            print('except') 
            voronoi = controller.organ.mH_settings['wf_info']['centreline']['vmtk_CL']['voronoi']
        done = fcM.extract_cl(organ=controller.organ, 
                                    win=controller.main_win, 
                                    voronoi=voronoi)
        if not done: 
            title = 'Something went wrong in VMTK...'
            msg = 'Something went wrong when trying to extract the centreline using VMTK. Please go back and check all the meshes created in MeshLab and re-run this process to extract the centrelines.' 
            prompt = Prompt_ok(title, msg, parent=controller.main_win)
            prompt.exec()
            print('output:',prompt.output, '\n')
            return
        
        #Enable button for select
        getattr(controller.main_win, 'centreline_select').setEnabled(True)

        #Toggle button
        getattr(controller.main_win, 'centreline_vmtk_play').setChecked(True)

        #Update progress in main_win
        controller.main_win.update_workflow_progress()
    else: 
        controller.main_win.win_msg('*To extract the centrelines using VMTK make sure you have run the  -MeshLab-  cleanup first.')

def run_centreline_select(controller):
    set_process(controller, 'centreline')
    if controller.main_win.centreline_vmtk_play.isChecked(): 

        nPoints = controller.main_win.gui_centreline['buildCL']['nPoints']
        same_plane = controller.main_win.gui_centreline['SimplifyMesh']['same_planes']

        workflow = controller.organ.workflow['morphoHeart']
        process = ['MeshesProc','C-Centreline','buildCL','Status']
        cl_names = list(controller.organ.mH_settings['measure']['CL'].keys())
        nn = 0
        for name in cl_names: 
            ch, cont, _ = name.split('_')
            print('name:', name)
            proc_wft = ['MeshesProc', 'C-Centreline', 'buildCL', ch, cont, 'Status']

            dict_clOpt = fcM.create_CLs(organ=controller.organ, 
                                        name=name,
                                        nPoints = nPoints)
            
            title = 'Select best centreline for tissue-contour' 
            msg = 'Select the preferred centreline for processing this tissue-contour ('+ch+'-'+cont+')'
            items = {1: {'opt': 'Option 1'}, 2: {'opt': 'Option 2'}, 3: {'opt': 'I would like to manually define the centreline'}}
            prompt = Prompt_ok_cancel_radio(title, msg, items, parent=controller.main_win)
            prompt.exec()
            print('output:', prompt.output, '\n')  
            opt_selected = prompt.output[0]
            cl_selected = items[opt_selected]['opt']

            proc_set = ['wf_info','centreline','buildCL','connect_cl',ch+'_'+cont]
            if opt_selected == 3:
                #Here, add option to align centreline
                happy = False
                while not happy: 
                    cl_final, dict_clOpt = fcM.define_CL(organ=controller.organ, 
                                                        name=name,
                                                        nPoints = nPoints)
                    
                    title = 'Check centreline...'
                    msg = 'Are you happy with the created centreline?! \n If so, select  -OK-, else select  -Cancel- and redefine it.'
                    prompt2 = Prompt_ok_cancel(title, msg, parent=controller.main_win)
                    prompt2.exec()
                    print('output:', prompt2.output)
                    happy = prompt2.output
                    update = 'User Modif.'
                    #Add text to Opt Centreline
                    getattr(controller.main_win, 'opt_cl'+str(nn+1)).setText(update)

            else: 
                update = cl_selected +'-'+dict_clOpt[items[opt_selected]['opt']]['description']
                cl_final = dict_clOpt[cl_selected]['kspl']
                #Add text to Opt Centreline
                getattr(controller.main_win, 'opt_cl'+str(nn+1)).setText(prompt.output[1])
            
            controller.organ.update_settings(process = proc_set, update = update, mH='mH')
            
            #Add centreline to organ
            controller.organ.add_object(cl_final, proc='Centreline', class_name=ch+'_'+cont, name='KSpline')
            controller.organ.obj_meshes[ch+'_'+cont].set_centreline()
            
            # Update organ workflow
            controller.organ.update_mHworkflow(process = proc_wft, update = 'DONE')
            status_sq = getattr(controller.main_win, 'opt_cl_status'+str(nn+1))
            controller.main_win.update_status(workflow, proc_wft, status_sq)

            #Enable button for plot cl
            getattr(controller.main_win, 'cl_plot'+str(nn+1)).setEnabled(True)

            #Update progress in main_win
            controller.main_win.update_workflow_progress()
            nn+=1
        
        # Update organ workflow
        controller.organ.update_mHworkflow(process = process, update = 'DONE')
        controller.organ.update_mHworkflow(process = ['MeshesProc','C-Centreline','Status'], update = 'DONE')
        controller.organ.check_status(process='MeshesProc')

        #Update Status in GUI
        controller.main_win.update_status(workflow, ['MeshesProc','C-Centreline','Status'], 
                                        controller.main_win.centreline_status)

        #Toggle button
        getattr(controller.main_win, 'centreline_select').setChecked(True)

        #Check if user selected measuring centreline length
        cl_measurements = controller.organ.mH_settings['setup']['params']['2']['measure']
        cl_measure = [cl_measurements[key] for key in cl_measurements.keys()]
        if any(cl_measure): 
            fcM.measure_centreline(organ=controller.organ, nPoints=nPoints)
            controller.main_win.fill_results()

        #If orientation is on_hold
        if controller.organ.on_hold: 
            run_axis_orientation(controller=controller, only_roi=True)

    else: 
        controller.main_win.win_msg('*To select the centreline to use make sure you have created the centrelines using VMTK first.')
    
def prompt_meshLab(controller, msg_add=None):

    # Prompt with instructions
    title = 'Instructions for cleaning up meshes with MeshLab!'
    msg = []
    if msg_add != None: 
        msg = msg+msg_add
    else: 
        msg.append('You are done in morphoHeart for a little while. To get the centreline of each of the selected meshes follow the next steps:')
    msg.append("   > 1. Open the .stl file(s) in Meshlab")
    msg.append("   > 2. Run Filters > Remeshing, Simplification.. > Screened Poisson Surf Reco (check Pre-clean)")
    msg.append("   > 3. Cut inflow and outflow tract as close as the cuts from the original mesh you opened in Meshlab.") 
    msg.append("   > 4. Export the resulting surface adding 'ML' at the end of the filename (e.g _cut4clML.stl) in the same folder")
    msg.append("   > 5. Come back, press OK and continue processing!...")
    prompt = Prompt_ok_cancel_Big(title, msg, parent=controller.main_win)
    prompt.exec()
    print('output:',prompt.output, '\n')

# > HEATMAPS
def run_heatmaps3D(controller, btn):
    set_process(controller, 'heatmaps')
    if controller.main_win.keeplargest_play.isChecked(): 
        workflow = controller.organ.workflow['morphoHeart']
        thck_values = {'i2e': {'short': 'th_i2e',
                                'method': 'D-Thickness_int>ext', 
                                'param': 'thickness int>ext', 
                                'n_type': 'int>ext'},
                    'e2i': {'short': 'th_e2i', 
                                'method': 'D-Thickness_ext>int', 
                                'param': 'thickness ext>int',
                                'n_type': 'ext>int'}}

        hm3d_list = list(controller.main_win.hm_btns.keys())
        if btn != None: 
            for key in hm3d_list: 
                num_key = controller.main_win.hm_btns[key]['num']
                if int(num_key) == int(btn): 
                    hm3get = key
                    break
            hm3d_set = [hm3get] #[segm_list[int(num)-1]]
        else: 
            hm3d_set = hm3d_list
        print('hm3d_set:',hm3d_set)

        for hmitem in hm3d_set:
            short, ch_info = hmitem.split('[') #short = th_i2e, th_e2i, ball
            ch_info = ch_info[:-1]
            if 'th' in short: 
                _, th_val = short.split('_')
                ch, cont = ch_info.split('-')
                method = thck_values[th_val]['method']
                mesh_tiss = controller.organ.obj_meshes[ch+'_tiss'].legend
                print('\n>> Extracting thickness information for '+mesh_tiss+'... \nNOTE: it takes about 5min to process each mesh... just be patient :) ')
                controller.main_win.win_msg('Extracting thickness information for '+mesh_tiss+'... NOTE: it takes about 5min to process each mesh... just be patient :)')
                setup = controller.main_win.gui_thickness_ballooning[hmitem]
                fcM.get_thickness(organ = controller.organ, name = (ch, cont), 
                                    thck_dict = thck_values[th_val], 
                                    setup = setup)

            else: # if 'ball' in short
                ch_cont, cl_info = ch_info.split('(')
                ch, cont = ch_cont.split('-')

                cl_info = cl_info[:-1].split('.')[1]
                cl_ch, cl_cont = cl_info.split('-')
                mesh2ball = controller.organ.obj_meshes[ch+'_'+cont].legend
                controller.main_win.win_msg('Extracting centreline>tissue information for '+mesh2ball+'... NOTE: it takes about 10-15 to process each mesh... just be patient :)')
                setup = controller.main_win.gui_thickness_ballooning[hmitem]

                fcM.extract_ballooning(organ = controller.organ, name = (ch, cont),
                                    name_cl = (cl_ch, cl_cont), setup = setup)

            #Enable buttons to plot heatmaps
            controller.main_win.hm_btns[hmitem]['plot'].setEnabled(True)
            hm2d_btn = controller.main_win.hm_btns[hmitem]['play2d']
            num = controller.main_win.hm_btns[hmitem]['num']
            d3d2_btn = getattr(controller.main_win, 'd3d2_'+str(num))
            #Hereee!
            if d3d2_btn.isChecked() and controller.main_win.thickness2D_set.isChecked(): 
                hm2d_btn.setEnabled(True)

            # controller.main_win.prog_bar_update(nn)
            #Update progress in main_win
            controller.main_win.update_workflow_progress()

        # Update organ workflow
        all_all_done = []
        processes = ['D-Thickness_int>ext', 'D-Thickness_ext>int', 'D-Ballooning']
        for proc in processes:
            all_done = []
            if len(workflow['MeshesProc'][proc].keys())>0: 
                for ch in workflow['MeshesProc'][proc].keys():
                    if ch != 'Status':
                        for cont in workflow['MeshesProc'][proc][ch].keys():
                            print(proc, ch, cont, workflow['MeshesProc'][proc][ch][cont]['Status'])
                            all_done.append(workflow['MeshesProc'][proc][ch][cont]['Status'])

                if all(flag == 'DONE' for flag in all_done): 
                    proc_wft = ['MeshesProc', proc, 'Status']
                    controller.organ.update_mHworkflow(process = proc_wft, update = 'DONE')
                    all_all_done.append('DONE')
                elif any(flag == 'DONE' for flag in all_done):
                    proc_wft = ['MeshesProc', proc, 'Status']
                    controller.organ.update_mHworkflow(process = proc_wft, update = 'Initialised')
                    all_all_done.append('Initialised')
                else: 
                    pass
            else: 
                all_all_done.append(True)

        # Update mH_settings
        proc_set = ['wf_info']
        update = controller.main_win.gui_thickness_ballooning
        controller.organ.update_settings(proc_set, update, 'mH', add='heatmaps')

        #Update Status in GUI
        if all(flag == 'DONE' for flag in all_all_done): 
            process = ['MeshesProc', processes[0], 'Status']
            toggle = True
        elif any(flag == 'DONE' for flag in all_all_done):
            for proc in processes: 
                test_proc = ['MeshesProc', proc, 'Status']
                if get_by_path(workflow, test_proc) == 'Initialised' or get_by_path(workflow, test_proc) == 'NI': 
                    process = test_proc
                    break
            toggle = False
        else: 
            process = ['MeshesProc', processes[0], 'Status']
            toggle = False
        controller.main_win.update_status(workflow, process, controller.main_win.heatmaps_status)
                
        #Toggle button
        if toggle: 
            getattr(controller.main_win, 'heatmaps3D_play').setChecked(True)
        else: 
            getattr(controller.main_win, 'heatmaps3D_play').setChecked(False)

        print('\nEND Heatmaps')
        print('organ.mH_settings:', controller.organ.mH_settings)
        print('organ.workflow:', workflow)

        #Add meshes to plot_user
        controller.main_win.fill_comboBox_all_meshes()

    else: 
        controller.main_win.win_msg('*To extract the thickness meshes make sure you have at least run the  -Keep Largest-  section.')

def run_heatmaps2D(controller, btn):
    set_process(controller, 'heatmaps')
    if controller.main_win.keeplargest_play.isChecked():
        if controller.main_win.thickness2D_set.isChecked():  
            workflow = controller.organ.workflow['morphoHeart']
            hm2d_list = list(controller.main_win.hm_btns.keys())
            if btn != None: 
                for key in hm2d_list: 
                    num_key = controller.main_win.hm_btns[key]['num']
                    if int(num_key) == int(btn): 
                        hm2get = key
                        break
                hm2d_set = [hm2get] #[segm_list[int(num)-1]]
            else: 
                hm2d_set = hm2d_list
            print('hm2d_set:',hm2d_set)

            gui_heatmaps2d = controller.main_win.gui_thickness_ballooning['heatmaps2D']

            for hmitem in hm2d_set: 
                short, ch_info = hmitem.split('[') #short = th_i2e, th_e2i, ball
                ch_info = ch_info[:-1]
                if 'th' in short: 
                    ch, contf = ch_info.split('-')
                    if 'i2e' in short: 
                        n_type = 'intTOext'
                        cont = 'ext'
                    else:# 'e2i' in short: 
                        n_type = 'extTOint'
                        cont = 'int'
                    array_name = 'thck('+n_type+')'
                else: 
                    # hm2d_set: ['ball[ch1-int(CL.ch1-int)]']
                    ch_cont, cl_info = ch_info.split('(CL.')
                    ch, contf = ch_cont.split('-')
                    from_cl, from_cl_type = cl_info.split('-')
                    array_name = 'ballCL('+from_cl+'_'+from_cl_type
                    cont = contf

                print('array_name:', array_name)
                #Get whole mesh and heatmap3d data
                whole_mesh = controller.organ.obj_meshes[ch+'_'+contf]
                array_mesh = controller.organ.obj_meshes[ch+'_'+cont]
                
                #Get data from heatmap 3D - LS52_F02_ch1_tiss_CMthck(intTOext)
                print(whole_mesh.dirs['arrays'])
                
                title = str(whole_mesh.dirs['arrays'][array_name])+'.npy'
                dir_npy = whole_mesh.parent_organ.dir_res(dir='csv_all') / title
                npy_array = np.load(dir_npy)
                print(short+'> whole_mesh: '+ch+'_'+contf+' -- array_mesh: '+ch+'_'+cont+' - loaded_npy:', dir_npy)

                #Get the extended and cut centrelines
                # Get the centreline ribbon
                cl_ribbon, kspl_ext = controller.main_win.plot_cl_ext1D(plotshow=False)

                # KSpline on surface
                ext_plane = controller.main_win.gui_thickness_ballooning['heatmaps2D']['direction']['plane_normal']
                kspl_vSurf = fcM.get_extCL_on_surf(mesh = array_mesh.mesh, kspl_ext = kspl_ext, direction = ext_plane)

                # Create new extended and cut kspline with higher resolution
                kspl_CLnew = fcM.get_extCL_highRes(organ = controller.organ, mesh = array_mesh.mesh, 
                                                    kspl_ext = kspl_ext) 
                
                #If the user wants to use the segments to aid heatmap unlooping....
                if gui_heatmaps2d['use_segms']: 
                    cut2use = gui_heatmaps2d['segms'].split(': ')[0]
                    segm_setup = controller.organ.mH_settings['setup']['segm'][cut2use]
                    segm_cuts_info = controller.organ.mH_settings['wf_info']['segments']['setup'][cut2use]['cut_info']
                    print('segm_setup:',segm_setup)

                    #See if the meshes have already been cut
                    obj_segm = {}
                    for n_segm in range(1,segm_setup['no_segments']+1,1):
                        segm_name = cut2use+'_'+ch+'_'+cont+'_segm'+str(n_segm)
                        if segm_name in controller.organ.submeshes.keys():
                            key = 'segm'+str(n_segm)+':'+segm_setup['name_segments']['segm'+str(n_segm)]
                            print('segm_name:', segm_name)
                            obj_segm[key] = controller.organ.obj_subm[segm_name]
                        else: 
                            controller.main_win.win_msg('*The segments for this tissue ('+ch+'_'+cont+") have not been obtained yet! Please create this tissue's segments to be able to run this process.") 
                            controller.main_win.hm_btns[hmitem]['play2d'].setChecked(False)
                            print('*The segments for this tissue have not been obtained yet! ('+controller.main_win.gui_thickness_ballooning['heatmaps2D']['segms']+')')
                            return
                        
                    print('len(obj_segm):', len(obj_segm), obj_segm)

                    data = {hmitem: npy_array}
                    #Use the whole tissue to create 2D heatmap
                    if len(obj_segm) != segm_setup['no_segments']: 
                        controller.main_win.win_msg('*The number of segments saved for this tissue ('+ch+'_'+cont+") is different than the number of segments expected. Please re-run segment creation for this tissue to be able to run this process.") 
                        print('*The number of segments saved for this tissue ('+ch+'_'+cont+") is different than the number of segments expected. Please re-run segment creation for this tissue to be able to run this process.")
                        return

                    else: 
                        # Classify points
                        print('Classify points into segments')
                        df_classPts, class_name = fcM.classify_heart_pts(array_mesh.mesh, obj_segm, data)

                        # Create kspline for each segment
                        print('controller.main_win.ordered_kspl:',controller.main_win.ordered_kspl)
                        ordered_kspl = copy.deepcopy(controller.organ.mH_settings['wf_info']['heatmaps']['heatmaps2D']['div'])
                        ordered_kspl = fcM.kspl_chamber_cut(organ = controller.organ, 
                                                            mesh = array_mesh.mesh, 
                                                            kspl_CLnew = kspl_CLnew, 
                                                            segm_cuts_info=segm_cuts_info, 
                                                            cut=cut2use, ordered_segm=ordered_kspl)
                        print('ordered_kspl (segm):',ordered_kspl)
                        order_segm_upside = list(ordered_kspl.keys())
                        order_segm_upside.reverse()
                
                #If the user wants to unloop the whole tissue without segments aid....
                else:
                    data = {hmitem: npy_array} 
                    #Use the whole tissue to create 2D heatmap
                    print('Use the whole tissue to create 2D heatmap')
                    df_classPts, class_name = fcM.classify_heart_pts(array_mesh.mesh, obj_segm=[], data=data)
                    # Create kspline for whole tissue
                    ordered_kspl = copy.deepcopy(controller.organ.mH_settings['wf_info']['heatmaps']['heatmaps2D']['div'])
                    print(ordered_kspl)
                    kspl_CLnewf = vedo.KSpline(kspl_CLnew.points()[::-1], res = 600).lw(5).color('deeppink').legend('HighResCL')
                    ordered_kspl['div1']['kspl'] = kspl_CLnewf
                    ordered_kspl['div1']['num_pts_range'] = (0, kspl_CLnewf.NPoints())
                    print('ordered_kspl (whole):',ordered_kspl)
                    order_segm_upside = ordered_kspl

                ##################################################################################################
                # Now Unloop the chambers (or not if whole) and get heatmap
                dirs = {}
                for div in order_segm_upside: 
                    print('\n\n- Unlooping the heart chambers for '+ordered_kspl[div]['name']+'...')
                    # print('div:', div, 'ordered_kspl[div]:', ordered_kspl[div])
                    # print(array_name, short, hmitem)
                    df_unloopedf, title_df = fcM.unloop_chamber(organ = controller.organ,
                                                                mesh = array_mesh.mesh, 
                                                                kspl_CLnew = kspl_CLnew,
                                                                kspl_vSurf = kspl_vSurf,
                                                                df_classPts = df_classPts,
                                                                labels = (hmitem, class_name),
                                                                gui_heatmaps2d = gui_heatmaps2d, 
                                                                kspl_data=ordered_kspl[div])
                    dirs[div] = Path(title_df)
                    
                    fcM.heatmap_unlooped(organ = controller.organ, kspl_data = ordered_kspl[div], 
                                                df_unloopedf = df_unloopedf, hmitem= hmitem, ch = ch,
                                                gui_thball = controller.main_win.gui_thickness_ballooning)
                    
                    #Update dirs to check if things have been made
                    controller.main_win.gui_thickness_ballooning[hmitem]['hm2d_dirs'] = dirs
                                
                #Enable buttons to plot heatmaps
                plot_btn = controller.main_win.hm_btns[hmitem]['plot2d']
                plot_btn.setEnabled(True)
                #Enable button to filter heatmaps
                controller.main_win.filter_2dhm.setEnabled(True)

            print('\nEND Heatmaps')
            print('organ.mH_settings:', controller.organ.mH_settings)
            print('organ.workflow:', workflow)

        else: 
            controller.main_win.win_msg('*Set the 2D Heatmap settings first to be able to run this process')
    else: 
        controller.main_win.win_msg('*To extract the 2D Heatmaps make sure you have run the  -Keep Largest-  section, extracted the Centreline and created the 3D Heatmap.')

# > SEGMENTS
def run_segments(controller, btn): 
    set_process(controller, 'segments')
    workflow = controller.organ.workflow['morphoHeart']
    segm_list = list(controller.main_win.segm_btns.keys())
    if btn != None: 
        cut, num = btn.split('_')
        for key in segm_list: 
            cut_key = key.split(':')[0]
            num_key = controller.main_win.segm_btns[key]['num']
            if cut_key == cut and int(num_key) == int(num): 
                segm2cut = key
                break
        segm_set = [segm2cut] #[segm_list[int(num)-1]]
    else: 
        segm_set = segm_list
    print('segm_set:',segm_set)

    #Setup everything to cut
    if controller.main_win.gui_segm['use_centreline']: 
        cl_name = controller.main_win.gui_segm['centreline'].split('(')[1][:-1]
        cl_ch, cl_cont = cl_name.split('_')
        if workflow['MeshesProc']['C-Centreline']['buildCL'][cl_ch][cl_cont]['Status'] != 'DONE':
            controller.main_win.win_msg('*Finish obtaining the centreline ('+cl_name+') to continue running this section!')
            controller.main_win.segm_btns[segm_set[0]]['play'].setChecked(False)
            return 
        else: 
            #Get centreline
            nPoints = controller.organ.mH_settings['wf_info']['centreline']['buildCL']['nPoints']
            cl = controller.organ.obj_meshes[cl_name].get_centreline(nPoints = nPoints)
            spheres_spl = fcM.sphs_in_spline(kspl = cl, colour = True)
            cl_spheres = {'centreline': cl, 
                        'spheres': spheres_spl, 
                        'nPoints' : nPoints}
    else: 
        cl_spheres = None

    # print('cl_spheres: ', cl_spheres)
    print('organ.obj_temp:', controller.organ.obj_temp)

    #Loop through all the tissues that are going to be segmented
    for segm in segm_set: 
        print('Cutting segm:', segm)
        #Find cut
        cut, ch_cont = segm.split(':')
        ch, cont = ch_cont.split('_')
        #Extract info for cut
        colors_all = controller.organ.mH_settings['setup']['segm'][cut]['colors']
        palette = [colors_all[key] for key in colors_all.keys()]
        segm_names = controller.organ.mH_settings['setup']['segm'][cut]['name_segments']
        
        #Find method to cut
        method = controller.organ.mH_settings['wf_info']['segments']['setup'][cut]['ch_info'][ch][cont]
        mesh2cut = controller.organ.obj_meshes[ch+'_'+cont]
        print('Cutting into segments:', mesh2cut.name, '- method: ', method)

        #Get usernames string
        user_names = '('+', '.join([segm_names[val] for val in segm_names])+')'
        print('\n- Dividing '+mesh2cut.legend+' into segments '+user_names)
        controller.main_win.win_msg('Dividing '+mesh2cut.legend+' into segments '+user_names)

        if method == 'ext-ext' or method == 'indep-ext': #TEST!!!
            # -> Get the discs that are going to be used to cut
            get_segm_discs(controller.organ, 
                            cut = cut, ch=ch, cont=cont, 
                            cl_spheres=cl_spheres, win=controller.main_win)
            # -> Create masks of discs
            # Get s3 dimensions
            mesh2cut.imChannel.load_chS3s([cont])
            s3_shape = getattr(mesh2cut.imChannel, 's3_'+cont).shape_s3
            fcM.create_disc_mask(controller.organ, cut = cut, s3_shape = s3_shape, h_min = 0.1125)
            ext_subsgm, meshes_segm = fcM.segm_ext_ext(controller.organ, mesh2cut, cut, 
                                                      segm_names, palette, win=controller.main_win)
            if ext_subsgm == None and meshes_segm == None: 
                btn2uncheck = controller.main_win.segm_btns[segm]['play']
                controller.main_win.win_msg('!Dividing '+mesh2cut.legend+' resulted in less segments than expected. Please re-run this process and make sure the defined disc(s) is(are) splitting the mesh into the expected number of segments.', btn2uncheck)
                title = 'Something went wrong when cutting the '+mesh2cut.legend+' into segments...'
                msg = 'Dividing '+mesh2cut.legend+' resulted in less segments than expected. Please re-run this process and make sure the defined disc(s) is(are) splitting the mesh into the expected number of segments.' 
                prompt = Prompt_ok(title, msg, parent=controller.main_win)
                prompt.exec()
                print('output:',prompt.output, '\n')
                return
            else: 
                #Add submeshes of ext_ext as attribute to organ
                controller.organ.ext_subsgm = ext_subsgm

                #Enable Plot Buttons
                controller.main_win.segm_btns[segm]['plot'].setEnabled(True)
                print('wf:', controller.organ.workflow['morphoHeart']['MeshesProc'])

                #Enable play buttons of meshes with other methods
                for sgmt in segm_list:
                    if sgmt != segm: 
                        play_btn = controller.main_win.segm_btns[sgmt]['play']
                        play_btn.setEnabled(True)
            
        elif method == 'cut_with_ext-ext' or method == 'cut_with_other_ext-ext':
            #Loading external subsegments 
            try: 
                ext_subsgm = controller.organ.ext_subsgm
                print('try ext_subsgm')
            except: 
                ext_subsgm = controller.organ.get_ext_subsgm(cut)
                print('except ext_subsgm')
            print('ext_subsgm: ',ext_subsgm)

            # -> Get segments using ext segments
            meshes_segm = fcM.get_segments(controller.organ, mesh2cut, cut, 
                                            segm_names, palette, ext_subsgm,  win=controller.main_win)
            
            #Enable Plot Buttons
            btn = controller.main_win.segm_btns[segm]['plot']
            btn.setEnabled(True)
            print('wf:', controller.organ.workflow['morphoHeart']['MeshesProc'])

        else: #method == 'cut_with_ext-indep'
            print('No functions for this method:', method)
            alert('error_beep')
            meshes_segm = None
            controller.main_win.win_msg('*Process to cut indep other than ext not yet coded: '+ method)
            controller.main_win.segm_btns[segm]['play'].setChecked(False)
            return 

        #Save meshes temporarily within segment buttons
        controller.main_win.segm_btns[segm]['meshes'] = meshes_segm
        print(controller.main_win.segm_btns)
        print('meshes_segm: ',meshes_segm)
        print(controller.main_win.segm_btns[segm])

        #Update progress in main_win
        controller.main_win.update_workflow_progress()
        #Fill-up results table
        controller.main_win.fill_results()
        #Check button 
        controller.main_win.segm_btns[segm]['play'].setChecked(True)
        #Check segm-sections
        if hasattr(controller.main_win, 'segm_sect_btns'):
            btn_final = 'NA'
            for btn_ss in controller.main_win.segm_sect_btns:
                if 's'+cut.title() in btn_ss and ch_cont in btn_ss:
                    btn_final = btn_ss
                    break
            if btn_final != 'NA': 
                reg_cut = btn_final.split(':')[0].split('_o_')[1]
                reg_btn = reg_cut+':'+ch_cont
                if controller.main_win.sect_btns[reg_btn]['plot'].isEnabled(): 
                    controller.main_win.segm_sect_btns[btn_final]['play'].setEnabled(True)
        
        if controller.organ.analysis['morphoCell']: 
            if controller.organ.mC_settings['setup']['segm_mC'] != False: 
                if controller.organ.mC_settings['setup']['segm_mC'][cut]['use_mH_settings']:
                    controller.main_win.update_status(None, 'DONE', getattr(controller.main_win, 'mH_segm_'+cut.lower()+'_done'), override=True)
        
    # Update organ workflow and GUI Status
    flat_semg_wf = flatdict.FlatDict(copy.deepcopy(workflow['MeshesProc']['E-Segments']))
    all_done = []
    for key in flat_semg_wf.keys(): 
        key_split = key.split(':')
        if len(key_split) > 1: 
            all_done.append(flat_semg_wf[key])

    proc_wft = ['MeshesProc', 'E-Segments', 'Status']
    if all(flag == 'DONE' for flag in all_done): 
        controller.organ.update_mHworkflow(process = proc_wft, update = 'DONE')
    elif any(flag == 'DONE' for flag in all_done): 
        controller.organ.update_mHworkflow(process = proc_wft, update = 'Initialised')
    else: 
        pass
    controller.main_win.update_status(workflow, proc_wft, controller.main_win.segments_status)
    
    #Add meshes to plot_user
    controller.main_win.fill_comboBox_all_meshes()

    print('organ.obj_temp:', controller.organ.obj_temp)

def get_segm_discs(organ, cut, ch, cont, cl_spheres, win): 

    if cl_spheres != None: 
        segm_using_cl = True
        cl = cl_spheres['centreline']
        spheres_spl = cl_spheres['spheres']
        nPoints = cl_spheres['nPoints']
    else:
        segm_using_cl = False
        cl = []

    # Get user segment names
    dict_names = organ.mH_settings['setup']['segm'][cut]['name_segments']
    user_names = '('+', '.join([dict_names[val] for val in dict_names])+')'

    # Number of discs expected to be created
    no_cuts_4segments = organ.mH_settings['setup']['segm'][cut]['no_cuts_4segments']
    #Get mesh
    mesh2cut = organ.obj_meshes[ch+'_'+cont]
    res = mesh2cut.resolution

    for n in range(no_cuts_4segments):
        happyWithDisc = False

        win.win_msg('Creating Disc No.'+str(n)+' for cutting tissues into segments!')
        while not happyWithDisc: 
            if segm_using_cl:
                msg = '\n> Define the centreline point number to use to initialise Disc No.'+str(n)+' to divide heart into segments '+user_names+' \n[NOTE: Spheres appear in centreline every 10 points, but you can select intermediate points, eg. 142].'
                txt = [(0, organ.user_organName + msg)]
                obj = [(mesh2cut.mesh.alpha(0.05), cl, spheres_spl)]
                plot_grid(obj=obj, txt=txt, axes=5, sc_side=max(organ.get_maj_bounds()))
                
                num_pt = None
                while num_pt == None: 
                    msg = 'Enter the centreline point number you want to use to initialise Disc No.'+str(n)+' to divide the tissue into segments '+user_names+':'
                    title = 'Centreline point number to initialise Disc'
                    prompt = Prompt_user_input(msg = msg, title = title, info = (0, nPoints), parent = win)
                    prompt.exec()
                    num_pt = prompt.output
                    print('prompt.output:', prompt.output)
                
                if prompt.output != None:
                    del prompt
                    cl_orient = None
                    while cl_orient == None: 
                        items = {0: {'opt': 'Use the centreline orientation at the selected point to define disc orientation', 'lineEdit': False, 'regEx': 'all'}, 
                                1: {'opt': 'Initialise the disc in the selected position but in plane normal to the y-z plane', 'lineEdit': False}}
                        title = 'Initilise disk orientation!'
                        msg = 'Select the orientation in which you would like to initialise the disc: '
                        prompt = Prompt_ok_cancel_radio(title, msg, items, parent = win)
                        prompt.exec()
                        cl_orient = prompt.output
                        print('output:', prompt.output, '\n')

                    if cl_orient[0] == 0: 
                        pl_normal, pl_centre = fcM.get_plane_normal2pt(pt_num = num_pt, points = cl.points())
                    else:# prompt.output[0] == 1: 
                        pl_centre = cl.points()[num_pt]
                        pl_normal = [1,0,0]
                    del prompt

            else: 
                pl_centre = mesh2cut.mesh.center_of_mass()
                pl_normal = [1,0,0]
            
            height = 2*0.225; disc_color = 'purple'; disc_res = 300
            # Modify (rotate and move cylinder/disc)
            radius = organ.mH_settings['wf_info']['segments']['radius'][cut]
            cyl_test, sph_test, rotX, rotY, rotZ = fcM.modify_disc(filename = organ.user_organName,
                                                                txt = 'cut tissues into segments '+user_names+' - Disc No.'+str(n+1)+'/'+str(no_cuts_4segments), 
                                                                mesh = mesh2cut.mesh,
                                                                option = [True,True,True,True,True,True],
                                                                def_pl = {'pl_normal': pl_normal, 'pl_centre': pl_centre},
                                                                radius = radius, height = height, 
                                                                color = disc_color, res = disc_res, 
                                                                zoom=0.8)

            # Get new normal of rotated disc
            pl_normal_corrected = fcM.new_normal_3DRot(normal = pl_normal, rotX = rotX, rotY = rotY, rotZ = rotZ)
            normal_unit = fcM.unit_vector(pl_normal_corrected)*10
            # Get central point of newly defined disc
            pl_centre_new = sph_test.pos()
            
            # Newly defined centreline point
            sph_cut = vedo.Sphere(pos = pl_centre_new, r=4, c='gold').legend('sph_ChamberCut')
    
            # Build new disc to confirm
            cyl_final = vedo.Cylinder(pos = pl_centre_new, r = radius, height = height, axis = normal_unit, c = disc_color, cap = True, res = disc_res)
            cyl_final.legend('Disc')

            msg = '\n> Check the position and the radius of Disc No.'+str(n)+' to cut the tissue into segments '+user_names+'.\nMake sure it is cutting the tissue effectively separating it into individual segments.\nClose the window when done.'
            txt = [(0, organ.user_organName + msg)]
            obj = [(mesh2cut.mesh, cyl_final, cl)]
            plot_grid(obj=obj, txt=txt, axes=5, sc_side=max(organ.get_maj_bounds()))

            happy = None
            while happy == None: 
                disc_radius = radius
                items = {0: {'opt': 'no, I would like to define a new position for the disc'}, 
                         1: {'opt': 'yes, but I would like to redefine the disc radius [um]', 'lineEdit': True, 'regEx': "int3d"}, 
                         2: {'opt': 'yes, I am happy with both, disc position and radius'}}
                title = 'Happy with the defined Disc No.'+str(n)+'?'
                msg = 'Are you happy with the position of the disc [radius: '+str(disc_radius)+'um] to cut tissue into segments  '+user_names+'?'
                prompt = Prompt_ok_cancel_radio(title, msg, items, parent = win)
                prompt.exec()
                happy = prompt.output
                del prompt
            if happy[0] == 1: 
                happy_rad = False
                disc_radius = int(happy[2])
                while not happy_rad: 
                    cyl_final = vedo.Cylinder(pos = pl_centre_new, r = disc_radius, height = height, axis = normal_unit, c = disc_color, cap = True, res = disc_res)
                    msg = '\n> New radius:'+str(disc_radius)+'um. \nCheck the radius of Disc No.'+str(n)+' to cut the tissue into segments '+user_names+'. \nMake sure it is cutting the tissue effectively separating it into individual segments.\nClose the window when done.'
                    txt = [(0, organ.user_organName + msg)]
                    obj = [(mesh2cut.mesh.alpha(1), cl, cyl_final, sph_cut)]
                    plot_grid(obj=obj, txt=txt, axes=5, sc_side=max(organ.get_maj_bounds()))

                    items = {0: {'opt': 'no, I would like to change its value', 'lineEdit': True, 'regEx': "int3d"}, 
                             1: {'opt': 'yes, it cuts the tissue without disrupting too much the segments!'}}
                    title = 'Happy with the defined Disc No.'+str(n)+'?'
                    msg = 'Is the selected radius ['+str(disc_radius)+'um] for Disc No.'+str(n)+' sufficient to cut the tissue into segments '+user_names+'?'
                    prompt = Prompt_ok_cancel_radio(title, msg, items, parent = win)
                    prompt.exec()
                    output = prompt.output
                    print('output:', output)
                    del prompt
                    if output[0] == 0: 
                        disc_radius = int(output[2])
                    else: 
                        happy_rad = True
                happyWithDisc = True
            elif happy[0] == 2:
                happyWithDisc = True
                
        # Save disc info
        cyl_name = 'Disc No.'+str(n)

        cyl_final.legend(cyl_name)
        cyl_dict = {'radius': disc_radius,
                    'normal_unit': normal_unit,
                    'pl_centre': pl_centre_new,
                    'height': height, 
                    'color': disc_color, 
                    'res': res}
        
        # Create new key in mH_settings to save disc info
        proc = ['wf_info', 'segments','setup', cut, 'cut_info']
        if 'cut_info' not in organ.mH_settings['wf_info']['segments']['setup'][cut].keys():
            organ.update_settings(proc, update = {}, mH = 'mH')
        organ.update_settings(proc+[cyl_name], update = cyl_dict, mH = 'mH')
        print('wf_info:', organ.mH_settings['wf_info']['segments']['setup'])

# > Angles Segments
def run_angles(controller, btn): 

    #Get df_res 
    df_res = fcM.df_reset_index(df = controller.organ.mH_settings['df_res'], mult_index=['Parameter', 'Tissue-Contour', 'User (Tissue-Contour)'])
    
    #Get direction of cut
    dir, cut = btn.split('_')
    axes = controller.main_win.gui_angles[cut]['axes'].split(',')
    ax_selected = axes[int(dir[-1])-1].strip()
    mtype = controller.main_win.gui_angles[cut]['axis_lab']
    axis = getattr(controller.main_win, 'angle_'+dir+'_'+cut).text()

    #Get cube
    cubes = getattr(controller.organ, mtype+'_cube')
    orient_cube = cubes['cube'].clone().alpha(0.1)
    cube_pts = orient_cube.points()

    #Get plane 
    pl_normal = controller.main_win.gui_orientation[mtype]['planar_views'][ax_selected]['pl_normal']
    id_cell = controller.main_win.gui_orientation[mtype]['planar_views'][ax_selected]['idcell']
    pl_centre = orient_cube.cell_centers()[id_cell]
    # sph = vedo.Sphere(pos = pl_centre, r=1, c='tomato')
    plane = vedo.Plane(pos=pl_centre, normal = pl_normal, s=(300,300)).alpha(0.5)
    if 'ref_vector' not in controller.main_win.gui_orientation[mtype]['planar_views'][ax_selected].keys(): 
        controller.main_win.win_msg('*Please re-run the "Organ/ROI Axis Labels" section to update Reference Axes and come back to run this section!')
        return False
    else:
        ref_vector_cube = np.array(controller.main_win.gui_orientation[mtype]['planar_views'][ax_selected]['ref_vector'])

    #Get the corner closer to the 0,0,0
    unique_corner_pts = np.unique(cube_pts, axis=0)
    min_dist = 10000000
    for n, pt in enumerate(unique_corner_pts):
        dist = np.linalg.norm(pt)
        if dist < min_dist: 
            min_dist = dist
            nf = n

    zero_corner = unique_corner_pts[nf]

    #Get external mesh
    mesh_name = controller.main_win.gui_angles[cut]['mesh_to_use']
    #Get the segments and iterate through them
    divs = controller.main_win.gui_angles[cut]['div']
    # div_settings = {}
    div_settings = {'proj_plane': {'pl_centre': pl_centre, 
                                    'pl_normal': pl_normal}}
    
    orig_line = {}
    for div in divs: 
        segm_no = divs[div]['segm']
        subm_name = cut.title()+'_'+mesh_name+'_'+segm_no
        subsegm = controller.organ.obj_subm[subm_name]

        #Create the line 
        p0 = divs[div]['segment_pts']['start']
        p1 = divs[div]['segment_pts']['end']
        orig_line[div] = {'p0': p1, 'p1': p0}

        sph_p0 = vedo.Sphere(pos = p0, r=2, c='black')
        sph_p1 = vedo.Sphere(pos = p1, r=2, c='gold')

        sph_proj_p0 = sph_p0.clone().project_on_plane(plane = plane).color('black').alpha(1)
        proj_p0 = sph_proj_p0.center_of_mass()
        sph_proj_p1 = sph_p1.clone().project_on_plane(plane = plane).color('gold').alpha(1)
        proj_p1 = sph_proj_p1.center_of_mass()

        #Measure angles
        vector_zero = np.array(([0.0,0.0,0.0], proj_p0-proj_p1))
        ref_vect = np.array(([0.0,0.0,0.0], ref_vector_cube))
        angle = fcM.find_angle_btw_pts(ref_vect, vector_zero)
        angle2 = fcM.find_angle_btw_pts_atan(vref = ref_vector_cube, 
                                             vect = proj_p0-proj_p1, 
                                             normal = pl_normal)
        print('angle-acos: ', angle, ' - angle-atan2:', angle2)

        div_settings[div] = {'proj_line': {'p0': proj_p1, 'p1': proj_p0}, 
                                'proj_vect_corner_box': {'p0': zero_corner, 'p1': zero_corner+proj_p0-proj_p1},
                                'ref_vector_corner_box': {'p0': zero_corner, 'p1': zero_corner+(ref_vector_cube*50)},
                                'vector_zero': vector_zero,
                                'vector_ref': ref_vect, 
                                'angle': angle2} 

        #Add angle to the df
        tiss, cont, segm = subsegm.sub_legend.split('_')
        cont_dict = {'ext': 'external', 'int': 'internal', 'tiss': 'tissue'}
        index = ('Angles: Segment-'+axis.title(), cut.title()+'_'+mesh_name+'_'+segm_no+'_'+dir, cut.title()+': '+ tiss+'-'+cont_dict[cont]+' ('+segm+'-'+axis+')')
        df_res = fcB.df_add_value(df = df_res, index = index, value=angle2)

        #Get things to plot
        if mH_config.dev: 
            arrow_normal = vedo.Arrow(start_pt=pl_centre, end_pt=pl_centre+(np.array(pl_normal)*50), c='fuchsia', s=0.1)
            sub_mesh = subsegm.get_segm_mesh()
            arrow_or = vedo.Arrow(start_pt=p1, end_pt=p0, c='darkred', s=0.1)
            clear_cube = cubes['clear'].clone().alpha(0.1)
            proj_arrow_or = vedo.Arrow(start_pt=proj_p1, end_pt=proj_p0, s=0.1, c='darkpurple').alpha(1)
            moved_proj_arrow_or = vedo.Arrow(start_pt = zero_corner, end_pt=zero_corner+proj_p0-proj_p1, s=0.1, c='darkpurple').alpha(1)
            col = controller.main_win.gui_orientation[mtype]['planar_views'][ax_selected]['color'][0:3]
            ref_vector_corner = vedo.Arrow(start_pt = zero_corner, end_pt=zero_corner+(ref_vector_cube*50), s=0.1, c=col).alpha(1)

            plt = vedo.Plotter(N=1, axes=1)
            plt.add_icon(logo, pos=(0.1,1), size=0.25)
            plt.show(sub_mesh, arrow_or, clear_cube, plane, arrow_normal, proj_arrow_or, moved_proj_arrow_or, ref_vector_corner, 
                     sph_p0, sph_p1, sph_proj_p0, sph_proj_p1, at=0, interactive=True)
        
    print('div_settings:', div_settings)
    #Resave df_res 
    controller.organ.mH_settings['df_res'] = df_res
    #Fill-up results table
    controller.main_win.fill_results()
    #Update original line
    controller.main_win.gui_angles[cut]['original_line'] = orig_line
    #Update gui_angles 
    if 'axis_settings' in controller.main_win.gui_angles[cut].keys(): 
        controller.main_win.gui_angles[cut]['axis_settings'][axis] = div_settings
    else: 
        controller.main_win.gui_angles[cut]['axis_settings'] = {axis : div_settings}

    #Update organ mH_settings
    proc_set = ['wf_info']
    update = controller.main_win.gui_angles
    controller.organ.update_settings(proc_set, update, 'mH', add='segm_angles')

    #Toggle button
    getattr(controller.main_win, 'play_angle_'+btn).setChecked(True)
    #Enable plot btn 
    getattr(controller.main_win, 'plot_angle_'+btn).setEnabled(True)
   
    return True

def run_ellipsoids(controller):

    #Get df_res 
    df_res = fcM.df_reset_index(df = controller.organ.mH_settings['df_res'], mult_index=['Parameter', 'Tissue-Contour', 'User (Tissue-Contour)'])
        
    #Get Cube, Planes to project things
    vect_ref = controller.main_win.gui_ellips['extended_dir']['plane_normal']
    ref_vect = fcM.unit_vector(vect_ref)
    mtype_cube = controller.main_win.gui_ellips['axis']

    #Get planar views
    planar_views = controller.organ.mH_settings['wf_info']['orientation'][mtype_cube]['planar_views']

    segm_btns = controller.main_win.segm_btns
    segm_ellips = []
    cont_dict = {'ext': 'external', 'int': 'internal', 'tiss': 'tissue'}
    for btn in list(segm_btns.keys()): 
        if segm_btns[btn]['play'].isChecked() and '_ext' in btn: 
            #Run Ellipsoid for this segm
            controller.main_win.win_msg('Creating ellipsoids for segments comprising '+btn)
            #Get segments
            cut, ch_cont = btn.split(':')
            ch, cont = ch_cont.split('_')
            divs = controller.main_win.gui_angles[cut.lower()]['div']
            
            for subm in controller.organ.mH_settings['setup']['segm'][cut]['name_segments'].keys():
                submesh_name = cut.title()+'_'+ch+'_'+cont+'_'+subm
                submesh = controller.organ.obj_subm[submesh_name]
                mesh_segm = submesh.get_segm_mesh()
                tiss, cont, segm_name = submesh.sub_legend.split('_')
                print(tiss, cont, segm_name)
                segm_name = controller.organ.mH_settings['setup']['segm'][cut]['name_segments'][subm]
                dims= fcM.measure_ellipse(organ = controller.organ, 
                                          mesh_segm = mesh_segm, segm_no = subm, 
                                            planar_views = planar_views, divs = divs, 
                                            cut = cut.lower(), ref_vect = ref_vect, 
                                            gui_angles = controller.main_win.gui_angles, 
                                            name = submesh.sub_legend,
                                            plot = False)
                
                print('Ellipsoid '+tiss+'-'+cont+'-'+segm_name+':',dims)
                segm_ellips.append(cut.title()+': '+ch_cont+'_'+subm+' ('+segm_name+')')
                for val in dims.keys(): 
                    ind_param = 'Ellipsoid: '+val
                    ind_tiss = cut.title()+'_'+ch_cont+'_'+subm
                    ind_tiss_user = cut.title()+': '+ tiss +'-'+cont_dict[cont]+' ('+segm_name+')'
                    index = (ind_param, ind_tiss, ind_tiss_user)
                    print(index)
                    df_res = fcB.df_add_value(df = df_res, index = index, value=dims[val])
    
    #Fill comboBox with measured ellipsoids
    controller.main_win.view_ellipsoid.clear()
    controller.main_win.view_ellipsoid.addItems(['----']+segm_ellips)
    controller.main_win.view_ellipsoid.setEnabled(True)
    #Resave df_res 
    controller.organ.mH_settings['df_res'] = df_res
    #Fill-up results table
    controller.main_win.fill_results()
    #Update organ mH_settings
    controller.main_win.gui_ellips['segm_ellips'] = segm_ellips
    proc_set = ['wf_info']
    update = controller.main_win.gui_ellips
    controller.organ.update_settings(proc_set, update, 'mH', add='gui_ellips')

    #Toggle button
    getattr(controller.main_win, 'ellipsoids_play').setChecked(True)
    #Enable plot btn 
    getattr(controller.main_win, 'plot_ellipsoid').setEnabled(True)
                
# > SECTIONS
def run_sections(controller, btn): 
    if controller.main_win.sections_set.isChecked(): 
        set_process(controller, 'regions')
        workflow = controller.organ.workflow['morphoHeart']
        sect_list = list(controller.main_win.sect_btns.keys())
        if btn != None: 
            cut, num = btn.split('_')
            for key in sect_list: 
                cut_key = key.split(':')[0]
                num_key = controller.main_win.sect_btns[key]['num']
                if cut_key == cut and int(num_key) == int(num): 
                    sect2cut = key
                    break
            sect_set = [sect2cut] 
        else: 
            sect_set = sect_list
        print('sect_set:',sect_set)

        #Setup everything to cut
        clRib_type = 'ext2sides'
        sect_settings = controller.organ.mH_settings['wf_info']['sections']

        #Loop through all the tissues that are going to be sectioned
        for sect in sect_set: 
            #Find cut
            cut, ch_cont = sect.split(':')
            ch, cont = ch_cont.split('_')
            #Extract info for cut
            colors_all = controller.organ.mH_settings['setup']['sect'][cut]['colors']
            palette = [colors_all[key] for key in colors_all.keys()]
            sect_names = controller.organ.mH_settings['setup']['sect'][cut]['name_sections']
            
            #Find mesh to cut
            mesh2cut = controller.organ.obj_meshes[ch+'_'+cont]

            #Check if the mask has already been saved
            if 'mask_name' in sect_settings[cut.title()].keys(): 
                #Add mask_name as attribute to organ in case it is not
                mask_name = sect_settings[cut.title()]['mask_name']
                if not hasattr(controller.organ, 'mask_sect_'+cut.lower()): 
                    setattr(controller.organ, 'mask_sect_'+cut.lower(), mask_name)
            else: 
                #Find centreline and settings to cut
                cl_name = controller.main_win.gui_sect[cut]['centreline'].split('(')[1][:-1]
                nPoints = controller.main_win.gui_sect[cut]['nPoints']
                mesh_cl = controller.organ.obj_meshes[cl_name]
                nRes = controller.main_win.gui_sect[cut]['nRes']
                filling_method = controller.main_win.gui_sect[cut]['filling_method']

                #Create mask_cube and save
                ext_plane = getattr(controller.main_win, 'extend_dir_'+cut.lower())['plane_normal']
                # -> Create ribbon
                if controller.main_win.reg_same_centreline.isChecked(): 
                    #Find if the othe KSpline extended has already been created
                    if cut == 'Cut1':
                        cuto = 'Cut2'
                    else: 
                        cuto = 'Cut1'
                    if 'ext_nPoints' in controller.main_win.gui_sect[cuto].keys(): 
                        ext_points = controller.main_win.gui_sect[cuto]['ext_pts']
                        use_prev = True
                    else: 
                        ext_points = None
                        use_prev = False
                else: 
                    ext_points = None
                    use_prev = False
                    
                #test this for a case in which the same centreline is used for both cuts
                #  (is it asking for happy the second time around?)
                happy = False
                while not happy: 
                    cl_ribbon, kspl_ext = mesh_cl.get_clRibbon(nPoints=nPoints, nRes=nRes, 
                                                                pl_normal=ext_plane, clRib_type=clRib_type, 
                                                                use_prev=use_prev, ext_points = ext_points)
                    
                    title = 'Check extended centreline...'
                    msg = 'Are you happy with the extended centreline/ribbon created?! \n If so, select  -OK-, else select  -Cancel- and redefine extended centreline.'
                    prompt = Prompt_ok_cancel(title, msg, parent=controller.main_win)
                    prompt.exec()
                    print('output:', prompt.output)
                    happy = prompt.output
                
                controller.main_win.gui_sect[cut]['ext_pts'] = kspl_ext.points()
                # proc_kspl_ext = ['wf_info', 'sections', cut, 'ext_pts']
                # controller.organ.update_settings(proc_kspl_ext, kspl_ext.points(), 'mH')

                # obj = [(mesh2cut.mesh, cl_ribbon)]
                # txt = [(0, controller.organ.user_organName+'- Extended Centreline Ribbon to cut organ into sections')]
                # plot_grid(obj=obj, txt=txt, axes=8, sc_side=max(controller.organ.get_maj_bounds()))
                
                # -> Create high resolution ribbon
                controller.main_win.win_msg('Creating high resolution centreline ribbon for '+cut.title())
                controller.organ.obj_imChannels[ch].load_chS3s([cont])
                s3_shape = getattr(controller.organ.obj_imChannels[ch], 's3_'+cont).shape_s3
                s3_filledCube, test_rib = fcM.get_stack_clRibbon(organ = controller.organ,
                                                                    s3_shape = s3_shape, 
                                                                    mesh_cl = mesh_cl, 
                                                                    cl_ribbon = cl_ribbon, 
                                                                    win=controller.main_win)
                
                # obj = [(test_rib, mesh_cl.mesh)]
                # txt = [(0, controller.organ.user_organName)]
                # plot_grid(obj=obj, txt=txt, axes=5, sc_side=max(controller.organ.get_maj_bounds()))

                # -> Create cube of ribbon and mask one side
                print('Creating cube sections for masking ('+cut.title()+')')
                controller.main_win.win_msg('Creating cube sections for masking ('+cut.title()+')')
                mask_cube_split, s3_filledCubes = fcM.get_cube_clRibbon(organ = controller.organ,
                                                                        s3_shape = s3_shape, 
                                                                        cut = cut,  
                                                                        s3_filledCube = s3_filledCube,
                                                                        res = mesh_cl.resolution,  
                                                                        pl_normal = ext_plane, 
                                                                        filling_method = filling_method)
                
                # obj = [(mask_cube_split[0], mask_cube_split[1], test_rib, mesh_cl.mesh)]
                # txt = [(0, controller.organ.user_organName)]
                # plot_grid(obj=obj, txt=txt, axes=5, sc_side=max(controller.organ.get_maj_bounds()))
                
                #Select the side of the ribbon that corresponds to section 1
                selected_side = fcM.select_ribMask(controller.organ, cut, mask_cube_split, mesh_cl.mesh, kspl_ext)
                
                if selected_side['name'] == 'NS': 
                    name_sect1 = controller.organ.mH_settings['setup']['sect'][cut.title()]['name_sections']['sect1']
                    title = name_sect1.title()+' section not selected!' 
                    msg = 'Select the mesh that corresponds to Section No.1 ('+name_sect1.upper()+'):'
                    items = {0: {'opt':'dark blue mesh (A)'}, 1: {'opt':'light blue mesh (B)'}}
                    prompt = Prompt_ok_cancel_radio(title, msg, items, parent=controller.main_win)
                    prompt.exec()
                    print('output:', prompt.output, '\n')  
                    if prompt.output[0] == 0: 
                        selected_side = {'name': 'Filled CLRibbon SideA'}
                    else: #prompt.output[0] == 1: 
                        selected_side = {'name': 'Filled CLRibbon SideB'}

                else: 
                    pass
                print('final selected_mesh:', selected_side)

                fcM.save_ribMask_side(organ = controller.organ, 
                                        cut = cut, 
                                        selected_side=selected_side, 
                                        s3_filledCubes = s3_filledCubes)
                    
                #Enable plot button for centreline extension
                cl_ext_btn = getattr(controller.main_win, 'cl_ext_'+cut.lower()).setEnabled(True)
            
            #Cut input tissue into sections
            print(type(mesh2cut))
            meshes_sect = fcM.get_sections(controller.organ, mesh2cut, cut, 
                                            sect_names, palette, win=controller.main_win)
            
            #Save meshes temporarily within section buttons
            controller.main_win.sect_btns[sect]['meshes'] = meshes_sect
            print(controller.main_win.sect_btns)
            print('meshes_sect: ',meshes_sect)
            print(controller.main_win.sect_btns[sect])

            #Enable Plot Buttons
            btn = controller.main_win.sect_btns[sect]['plot']
            btn.setEnabled(True)
            print('wf:', controller.organ.workflow['morphoHeart']['MeshesProc'])

            #Update progress in main_win
            controller.main_win.update_workflow_progress()
            #Fill-up results table
            controller.main_win.fill_results()
            #Check button 
            controller.main_win.sect_btns[sect]['play'].setChecked(True)
            #Check segm-sections
            if hasattr(controller.main_win, 'segm_sect_btns'):
                btn_final = 'NA'
                for btn_ss in controller.main_win.segm_sect_btns:
                    if '_o_'+cut.title() in btn_ss and ch_cont in btn_ss:
                        btn_final = btn_ss
                        break
                if btn_final != 'NA': 
                    seg_cut = btn_final.split('_o_')[0][1:]
                    seg_btn = seg_cut+':'+ch_cont
                    if controller.main_win.segm_btns[seg_btn]['plot'].isEnabled(): 
                        controller.main_win.segm_sect_btns[btn_final]['play'].setEnabled(True)
        
        # Update organ workflow and GUI Status
        flat_sect_wf = flatdict.FlatDict(copy.deepcopy(workflow['MeshesProc']['E-Sections']))
        all_done = []
        for key in flat_sect_wf.keys(): 
            key_split = key.split(':')
            if len(key_split) > 1: 
                all_done.append(flat_sect_wf[key])

        proc_wft = ['MeshesProc', 'E-Sections', 'Status']
        if all(flag == 'DONE' for flag in all_done): 
            controller.organ.update_mHworkflow(process = proc_wft, update = 'DONE')
        elif any(flag == 'DONE' for flag in all_done): 
            controller.organ.update_mHworkflow(process = proc_wft, update = 'Initialised')
        else: 
            pass
        controller.main_win.update_status(workflow, proc_wft, controller.main_win.sections_status)

        #Add meshes to plot_user
        controller.main_win.fill_comboBox_all_meshes()
    else: 
        controller.main_win.win_msg('*You need to re-set the settings of the regions to be able to run this process.')

# > SEGM-SECT
def run_segm_sect(controller, btn): 
    set_process(controller, 'segments-regions')
    workflow = controller.organ.workflow['morphoHeart']
    segm_sect_list = list(controller.main_win.segm_sect_btns.keys())
    if btn != None: 
        scuto, rcut_num = btn.split('_o_')
        scut = scuto[1:]
        rcut, num = rcut_num.split('_')
        for key in segm_sect_list: 
            ch_cont = key.split(':')[1]
            seg_cut = key.split('_o_')[0][1:]
            reg_cut = key.split(':')[0].split('_o_')[1]
            num_key = controller.main_win.segm_sect_btns[key]['num']
            if seg_cut == scut and reg_cut == rcut and int(num_key) == int(num): 
                segmsect2cut = key
                break
        segm_sect_set = [segmsect2cut] 
    else: 
        segm_sect_set = segm_sect_list
    print('segm_sect_set:',segm_sect_set)

    #Loop through all the tissues that are going to be sectioned
    for segm_sect in segm_sect_set: 
        ch_cont = segm_sect.split(':')[1]
        ch, cont = ch_cont.split('_')
        seg_cut = segm_sect.split('_o_')[0][1:]
        reg_cut = segm_sect.split(':')[0].split('_o_')[1]
        print(ch, cont, seg_cut, reg_cut)

        #Extract info for cut
        colors_all = controller.organ.mH_settings['setup']['segm-sect']['s'+seg_cut][reg_cut]['colors']
        segm_names = controller.organ.mH_settings['setup']['segm'][seg_cut]['name_segments']
        sect_names = controller.organ.mH_settings['setup']['sect'][reg_cut]['name_sections']
        print(colors_all, '\n', segm_names, '\n', sect_names)

        #Get mask to use
        sect_settings = controller.organ.mH_settings['wf_info']['sections']
        if 'mask_name' in sect_settings[reg_cut.title()].keys(): 
            #Add mask_name as attribute to organ in case it is not
            mask_name = sect_settings[reg_cut.title()]['mask_name']
            if not hasattr(controller.organ, 'mask_sect_'+reg_cut.lower()): 
                setattr(controller.organ, 'mask_sect_'+reg_cut.lower(), mask_name)
        else: 
            print('error no mask? ')

        # Cut input tissue into sections
        meshes_segm_sect = fcM.get_segm_sects(controller.organ, 
                                                ch_cont = ch_cont, 
                                                cuts = (seg_cut, reg_cut), 
                                                names = (segm_names, sect_names), 
                                                palette=colors_all,
                                                win=controller.main_win)
        
        #Save meshes temporarily within section buttons
        controller.main_win.segm_sect_btns[segm_sect]['meshes'] = meshes_segm_sect
        print('meshes_segm_sect: ',meshes_segm_sect)
        print(controller.main_win.segm_sect_btns[segm_sect])

        #Enable Plot Buttons
        btn = controller.main_win.segm_sect_btns[segm_sect]['plot']
        btn.setEnabled(True)
        print('wf:', controller.organ.workflow['morphoHeart']['MeshesProc'])

    #Update progress in main_win
    controller.main_win.update_workflow_progress()

    #Fill-up results table
    controller.main_win.fill_results()

    # Update organ workflow and GUI Status
    flat_sect_wf = flatdict.FlatDict(copy.deepcopy(workflow['MeshesProc']['E-Segments_Sections']))
    all_done = []
    for key in flat_sect_wf.keys(): 
        key_split = key.split(':')
        if len(key_split) > 1: 
            all_done.append(flat_sect_wf[key])

    proc_wft = ['MeshesProc', 'E-Segments_Sections', 'Status']
    if all(flag == 'DONE' for flag in all_done): 
        controller.organ.update_mHworkflow(process = proc_wft, update = 'DONE')
    elif any(flag == 'DONE' for flag in all_done): 
        controller.organ.update_mHworkflow(process = proc_wft, update = 'Initialised')
    else: 
        pass
    controller.main_win.update_status(workflow, proc_wft, controller.main_win.segm_sect_status)
    
    #Add meshes to plot_user
    controller.main_win.fill_comboBox_all_meshes()

#> MEASURE
def run_measure(controller): 
    set_process(controller, 'measure')
    if controller.main_win.keeplargest_play.isChecked(): 
        organ = controller.organ
        measurements = organ.mH_settings['measure']
        all_done = []; whole_done = []; other_done = []
        df_res = fcB.df_reset_index(df=organ.mH_settings['df_res'], 
                                     mult_index= ['Parameter', 'Tissue-Contour'])
        print(controller.main_win.index_param)
        #{'Looping direction', 'Centreline: Linear Length', 'Ellipsoid: Depth', 
        # 'Ellipsoid: Length', 'Ellipsoid: Asphericity', 'Volume: Segm-Reg', 'Volume: Segment', 
        # 'Surface Area', 'Centreline: Looped Length', 'Embryo Notes', 'Aortic length', 
        # 'Ellipsoid: Width', 'Volume: Region', 'Volume', 'Angles: Segment'}

        for param in measurements.keys():
            for item in measurements[param]:
                if 'whole' in item: 
                    ch, cont, _ = item.split('_')
                    try: 
                        mesh2meas = organ.obj_meshes[ch+'_'+cont]
                        mesh = mesh2meas.mesh
                        if param == 'Vol': 
                            vol = mesh.volume()
                            df_res = fcB.df_add_value(df=df_res, index=('Volume', ch+'_'+cont+'_whole'), value=vol)
                            # organ.mH_settings['measure']['Vol'][mesh2meas.name+'_whole'] = vol
                        if param == 'SA': 
                            area = mesh.area()
                            df_res = fcB.df_add_value(df=df_res, index=('Surface Area', ch+'_'+cont+'_whole'), value=area)
                            # organ.mH_settings['measure']['SA'][mesh2meas.name+'_whole'] = area
                        all_done.append(True); whole_done.append(True)
                    except: 
                        print('Ch-Cont: '+ch+'_'+cont+' not found!')
                        all_done.append(False); whole_done.append(False)
                else: 
                    print('Variable to measure:', item, param)
                    all_done.append(False)
                    pass

        # Update organ workflow and GUI Status
        # Whole 
        if all(whole_done): 
            controller.main_win.update_status(None, 'DONE', controller.main_win.measure_whole_status, override=True)
        elif any(whole_done): 
            controller.main_win.update_status(None, 'Initialised', controller.main_win.measure_whole_status, override=True)
        else: 
            controller.main_win.update_status(None, 'NI', controller.main_win.measure_whole_status, override=True)

        # All 
        if all(all_done): 
            controller.main_win.update_status(None, 'DONE', controller.main_win.measure_status, override=True)
        elif any(all_done): 
            controller.main_win.update_status(None, 'Initialised', controller.main_win.measure_status, override=True)
        else: 
            controller.main_win.update_status(None, 'NI', controller.main_win.measure_status, override=True)
        
        #Fill-up results table
        df_res = fcB.df_reset_index(df=df_res, mult_index= ['Parameter', 'Tissue-Contour', 'User (Tissue-Contour)'])
        organ.mH_settings['df_res'] = df_res
        controller.main_win.fill_results()
    
    else: 
        controller.main_win.win_msg('*To measure whole tissue make sure you have at least run the  -Keep Largest-  section.')

#MORPHOCELL TAB
def run_isosurface(controller, btn): 

    #set_process(controller, 'isosurface_mC')
    workflow = controller.organ.workflow['morphoCell']
    res = fcM.create_iso_volume(organ = controller.organ, ch = btn)
    if not res: 
        controller.main_win.win_msg('*Somthing went wrong when loading image file for '+btn.title()+'. Please check to continue.')
        getattr(controller.main_win, btn+'_play').setChecked(False)
        return 
    
    proc = ['A-SetExtraChs', btn, 'Status']
    controller.main_win.organ.update_mCworkflow(process = proc, update = 'DONE')

    getattr(controller.main_win, btn+'_play').setChecked(True)
    getattr(controller.main_win, 'iso_plot_'+btn).setEnabled(True)
    all_done = []
    for ch in ['chB', 'chC', 'chD']: 
        btn_o = getattr(controller.main_win, ch+'_play')
        all_done.append(btn_o.isVisible() and btn_o.isChecked())
    
    proc_wft = ['A-SetExtraChs', 'Status']
    if all(all_done):
        controller.organ.update_mCworkflow(process = proc_wft, update = 'DONE')
    elif any(all_done): 
        controller.organ.update_mCworkflow(process = proc_wft, update = 'DONE')
    else: 
        pass
    controller.main_win.update_status(workflow, proc_wft, controller.main_win.isosurface_cells_status)    

def run_remove_cells(controller): 

    #set_process(controller, 'remove_cells_mC')
    workflow = controller.organ.workflow['morphoCell']
    # Check which channels want to be seen
    controller.organ.cellsMC['chA'].remove_cells()
    
    getattr(controller.main_win, 'remove_cells_play').setChecked(True)
    proc_wft = ['A-CleanCells', 'Status']
    controller.organ.update_mCworkflow(process = proc_wft, update = 'DONE')
    controller.main_win.update_status(workflow, proc_wft, controller.main_win.remove_cells_status)    

def run_segments_mC(controller, btn): 

    #set_process(controller, 'segments_mC')
    workflow = controller.organ.workflow['morphoCell']
    segm_cell_list = list(controller.main_win.cell_segm_btns.keys())
    if btn != None: 
        segm_cell_set = [btn+':chA']
    else: 
        segm_cell_set = segm_cell_list
    print('segm_cell_set:',segm_cell_set)

    #Loop through all the cuts that can be done
    for segm in segm_cell_set: 
        print('Cutting segm_cell:', segm)
        cut, ch = segm.split(':')
        #Extract info for cut
        segm_names = controller.organ.mC_settings['setup']['segm_mC'][cut]['name_segments']
        #Get usernames string
        user_names = '('+', '.join([segm_names[val] for val in segm_names])+')'
        print('\n- Dividing Cells into segments '+user_names)
        controller.main_win.win_msg('Dividing Cells into segments '+user_names)
        
        cells_position = controller.organ.cellsMC['chA'].df_cells()
        column_name = 'Segment-'+cut
        run_all = False; run_reclass = False; cells_out = None
        if len(controller.organ.mC_settings['wf_info']['segments_cells'][cut]['cut_info']) > 0: 
            if column_name in list(cells_position.columns): 
                #Ask if the used wants to re-run process: 
                items = {0: {'opt': 'I would like to redefine the plane and continue from there'}, 
                            1: {'opt': 'I would like to only re-classify some cells that are still misclassified'}}
                title = 'Cells have already been classified for '+cut+'...'
                msg = 'The cells from this organ have already been classified for '+cut+': '+user_names[1:-1]+'. Select from which process you would like to re-run cell classification:'
                prompt = Prompt_ok_cancel_radio(title, msg, items, parent = controller.main_win)
                prompt.exec()
                print(prompt.output)
                if prompt.output[0] == 0: 
                    run_all = True
                else: 
                    run_all = False
                    run_reclass = True
            else: 
                run_all = True
        else: 
            run_all = True
        
        if run_all: 
            get_segm_planes(organ = controller.organ, cut = cut, win = controller.main_win)
            cells_out, segm_class = fcM.get_cells_within_planes(controller = controller, 
                                                                organ = controller.organ,
                                                                cells_position = cells_position, 
                                                                cut=cut)
        
        if run_all or run_reclass: 
            if cells_out == None: 
                segm_class = cells_position['Segment-'+cut]
                color_class = controller.organ.cellsMC['chA'].get_colour_class(cut, mtype='segm')
                cells_out = controller.organ.cellsMC['chA'].colour_cells(sphs_pos = cells_position, 
                                                                        color_class = color_class)
            
            cells_out, segm_classf = fcM.modify_cell_class(organ = controller.organ, 
                                                            cells = cells_out, cut=cut, 
                                                            cells_class = segm_class)
            
            controller.organ.cellsMC['chA'].cells = cells_out
            cells_position = controller.organ.cellsMC['chA'].assign_class(cells_position, segm_classf, col_name = 'Segment-'+cut)
            controller.organ.cellsMC['chA'].save_cells(cells_position)

            fcM.count_cells(controller.organ, cells_position, cut, mtype = 'segm', col_name = 'Segment-'+cut)

            #Update workflow
            controller.main_win.cell_segm_btns[segm]['play'].setChecked(True)
            proc_wf = ['B-Segments', cut,'Status']
            controller.organ.update_mCworkflow(process = proc_wf, update = 'DONE')

            #Enable Plot Buttons
            controller.main_win.cell_segm_btns[segm]['plot'].setEnabled(True)
    
    all_done = []
    for cut in workflow['B-Segments'].keys(): 
        if cut != 'Status':
            all_done.append(workflow['B-Segments'][cut]['Status'])
    
    proc_wft = ['B-Segments','Status']
    if all(flag == 'DONE' for flag in all_done): 
        controller.organ.update_mCworkflow(process = proc_wft, update = 'DONE')
    elif any(flag == 'DONE' for flag in all_done): 
        controller.organ.update_mCworkflow(process = proc_wft, update = 'Initialised')
    else:
        controller.organ.update_mCworkflow(process = proc_wft, update = 'NI')

    controller.main_win.update_status(workflow, proc_wft, controller.main_win.cell_segments_status)
            
def get_segm_planes(organ, cut, win): 

    # Get user segment names
    dict_names = organ.mC_settings['setup']['segm_mC'][cut]['name_segments']
    user_names = '('+', '.join([dict_names[val] for val in dict_names])+')'

    # Number of discs expected to be created
    no_cuts_4segments = organ.mC_settings['setup']['segm_mC'][cut]['no_cuts_4segments']
    #Get iso_vols
    iso_vols, vol_settings = organ.organ_vol_iso()
    
    for n in range(no_cuts_4segments):
        happyWithPlane = False
        win.win_msg('Creating Plane No.'+str(n)+' for cutting tissues into segments!')

        while not happyWithPlane:
            try: 
                pl_centre = iso_vols[0].center_of_mass()
            except: # this except doesn't work - find the boundaries of the spheres in get_plane?
                cells_position = organ.cellsMC['chA'].df_cells()
                filt_segm = cells_position[["Position X", "Position Y", "Position Z"]]
                pl_centre = list(filt_segm.mean())

            pl_normal = (1,0,0)
            def_pl = {'pl_normal': pl_normal, 'pl_centre': pl_centre}
            _, pl_dict = fcM.get_plane(filename = organ.user_organName, 
                                            txt = 'cut cells into segments '+user_names+' - Plane No.'+str(n+1)+'/'+str(no_cuts_4segments), 
                                            meshes = iso_vols, 
                                            settings = vol_settings, 
                                            def_pl = def_pl)
            
            items = {0: {'opt': 'no, I would like to define a new position for the plane'}, 
                        1: {'opt': 'yes, I am happy with the defined plane'}}
            title = 'Happy with the defined Plane No.'+str(n)+'?'
            msg = 'Are you happy with the position of the plane to cut cells into segments  '+user_names+'?'
            prompt = Prompt_ok_cancel_radio(title, msg, items, parent = win)
            prompt.exec()
            happyWithPlane = prompt.output
            
        plane_name = 'Disc No.'+str(n)
        
        # Create new key in mC_settings to save disc info
        proc = ['wf_info', 'segments_cells', cut, 'cut_info']
        if 'cut_info' not in organ.mC_settings['wf_info']['segments_cells'][cut].keys():
            organ.update_settings(proc, update = {}, mH = 'mC')
        organ.update_settings(proc+[plane_name], update = pl_dict, mH = 'mC')
        print('wf_info:', organ.mC_settings['wf_info']['segments_cells'][cut])

def run_IND_segm(controller, plot=True): 

    wf_info = controller.organ.mC_settings['wf_info']['ind_segments_cells']
    measure = controller.organ.mC_settings['measure']['mC_segm']

    if 'cluster_group' not in wf_info.keys(): 
        wf_info['cluster_group'] = {}

    for cut in measure: 
        print('cut:', cut)
        for segm in measure[cut]:
            ss, user_name = segm.split(':')
            n_closest_cells = wf_info[cut][ss]
            segm_count = measure[cut][segm]['total']
            
            df_clusters_ALL, sphs_distColour_ALL, sil_distColour_ALL = fcM.extract_segm_IND(controller.organ, cut, segm, 
                                                                                            n_closest_cells, segm_count, plot)

        
            #Add prompt that asks which option to select
            if len(df_clusters_ALL) > 1:
                items = {}
                for opt in range(len(df_clusters_ALL)): 
                    items[opt] = {'opt': 'Option '+str(opt+1)}
                title = 'Clustering Options for '+cut+'-'+segm.title()
                msg = 'Select preferred clustering option for cells in '+cut+'-'+segm.title()
                prompt = Prompt_ok_cancel_radio(title, msg, items, parent = controller.main_win)
                prompt.exec()

                sel_option = int(prompt.output[0])
            else: 
                sel_option = 0
            
            #Add selected option to wf_info
            if cut not in wf_info['cluster_group']: 
                wf_info['cluster_group'][cut] = {ss: sel_option}
            else: 
                wf_info['cluster_group'][cut][ss] = sel_option
            
            #Save results in df?
            df_clusterf = df_clusters_ALL[sel_option]
            controller.organ.mC_settings['measure']['mC_segm'][cut][segm]['IND'] = df_clusterf

            if not hasattr(controller.organ, 'ind'): 
                controller.organ.ind = {'mC_segm':{cut:{ss: [sphs_distColour_ALL[sel_option],sil_distColour_ALL[sel_option]]}}}
            else: 
                if 'mC_segm' not in controller.organ.ind:
                    controller.organ.ind['mC_segm'] = {cut:{ss: [sphs_distColour_ALL[sel_option],sil_distColour_ALL[sel_option]]}}
                else: 
                    if cut not in controller.organ.ind['mC_segm']:
                        controller.organ.ind['mC_segm'][cut] = {ss: [sphs_distColour_ALL[sel_option],sil_distColour_ALL[sel_option]]}
                    else: 
                        controller.organ.ind['mC_segm'][cut][ss] = [sphs_distColour_ALL[sel_option],sil_distColour_ALL[sel_option]]
            
            #Enable plot button
            getattr(controller.main_win, cut.lower()+'_IND_'+ss+'_plot').setEnabled(True)

def run_zones(controller, zone): 
    
    # set_process(controller, 'segments')
    workflow = controller.organ.workflow['morphoCell']
    zone_cell_list = list(controller.main_win.cell_zone_btns.keys())
    if zone != None: 
        zone_cell_set = [zone+':chA']
    else: 
        zone_cell_set = zone_cell_list
    print('zone_cell_set:',zone_cell_set)
    #Loop through all the cuts that can be done
    for zone in zone_cell_set: 
        print('Adding cells to '+zone)
        zz, ch = zone.split(':')
        #Get setings for that zone
        zone_set = controller.organ.mC_settings['wf_info']['zone_cells'][zz]
        if zone_set['use_segm']: 
            segm_cut_info = zone_set['ind_segm']
            cut, names = segm_cut_info.split(':')
            #Check if the cut has been run
            if getattr(controller.main_win, cut.lower()+'_segm_cell_play').isChecked():
                #Get all names
                dict_names = controller.organ.mC_settings['setup']['segm_mC'][cut]['name_segments']
                segm_names = {}
                for key, item in dict_names.items():
                    segm_names[key] = key+':'+item
                #See if the segms selected are more than one
                all_names = names.split(', ')
                #Create an empty dictionary where values of all segm will be saved
                all_data_to_run = {}

                dict_sp = controller.organ.mC_settings['measure']['mC_segm'][cut]
                #Option 1: Use all segments for cut
                if len(all_names) > 1:
                    #Using all segments - check if they have been obtained 
                    for segm in segm_names: 
                        print(segm)
                        if 'IND' in dict_sp[segm_names[segm]]:
                            all_data_to_run[segm] = {}
                        else: 
                            controller.main_win.win_msg('*Extract the IND values for '+cut+':'+segm_names[segm].title()+' to run this process.', 
                                                        getattr(controller.main_win, zz.lower()+'_cell_play'))
                            return None
                        
                #Option 2: Only one segment
                else: 
                    #Get the segm that will be used
                    for ss in dict_names: 
                        if all_names[0] == dict_names[ss]: 
                            segm = ss
                            break
                    if 'IND' in dict_sp[segm_names[segm]]:
                        all_data_to_run[segm] = {}
                    else: 
                        controller.main_win.win_msg('*Extract the IND values for '+cut+':'+segm_names[segm].title()+' to run this process.', 
                                                    getattr(controller.main_win, zz.lower()+'_cell_play'))
                        return None

                for ss in all_data_to_run: 
                    #Run a for for all the segm in all_data_to_run
                    n_closest_cells = controller.organ.mC_settings['wf_info']['ind_segments_cells'][cut][ss]
                    segm_count = controller.organ.mC_settings['measure']['mC_segm'][cut][segm_names[ss]]['total']
                    index_closest_cells, min_dist_to_cells, avg_min_dist = extract_segm_IND(controller.organ, cut, 
                                                                                            segm_names[ss], n_closest_cells, 
                                                                                            segm_count, plot=False, 
                                                                                            process = 'zones')
                    sphs_distColour = controller.organ.ind['mC_segm'][cut][ss][0]
                    all_data_to_run[ss]['index_closest_cells'] = index_closest_cells
                    all_data_to_run[ss]['min_dist_to_cells'] = min_dist_to_cells
                    all_data_to_run[ss]['avg_min_dist'] = avg_min_dist
                    all_data_to_run[ss]['sphs_distColour'] = sphs_distColour

                happy = False
                while not happy: 
                    df_zones, sphs_zones, sil_zones = fcM.select_cell_for_zones(controller.organ, cut = cut, zone = zz, 
                                                        segm_names = segm_names, all_data = all_data_to_run)
                    title = 'Happy?'
                    msg = 'Are you happy with the selected cells/zone? If so, select  -OK-,  else select  -Cancel-  and select again the cells for all zones' 
                    prompt = Prompt_ok_cancel(title, msg, parent=controller.main_win)
                    prompt.exec()
                    print('output:',prompt.output, '\n')
                    if prompt.output: 
                        happy = True

                #Save results in df?
                if 'mC_zone' not in controller.organ.mC_settings['measure']:
                    controller.organ.mC_settings['measure']['mC_zone'] = {zz: df_zones}
                else: 
                    if zz not in controller.organ.mC_settings['measure']['mC_zone']:
                        controller.organ.mC_settings['measure']['mC_zone'][zz] = df_zones

                if not hasattr(controller.organ, 'ind'): 
                    controller.organ.ind = {'mC_zone':{zz: [sphs_zones, sil_zones]}}
                else: 
                    if 'mC_zone' not in controller.organ.ind:
                        controller.organ.ind['mC_zone'] = {zz: [sphs_zones, sil_zones]}
                    else: 
                        controller.organ.ind['mC_zone'][zz] = [sphs_zones, sil_zones]

                #Update workflow
                controller.main_win.cell_zone_btns[zz+':chA']['play'].setChecked(True)
                proc_wf = ['B-Zones', zz,'Status']
                controller.organ.update_mCworkflow(process = proc_wf, update = 'DONE')

                #Enable plot button
                controller.main_win.cell_zone_btns[zz+':chA']['plot'].setEnabled(True)

            else: 
                controller.main_win.win_msg('*Please run the  -Segments-  module to be able to use its settings.')

    all_done = []
    for zzn in workflow['B-Zones'].keys(): 
        if zzn != 'Status':
            all_done.append(workflow['B-Zones'][zzn]['Status'])
    
    proc_wft = ['B-Zones','Status']
    if all(flag == 'DONE' for flag in all_done): 
        controller.organ.update_mCworkflow(process = proc_wft, update = 'DONE')
    elif any(flag == 'DONE' for flag in all_done): 
        controller.organ.update_mCworkflow(process = proc_wft, update = 'Initialised')
    else:
        controller.organ.update_mCworkflow(process = proc_wft, update = 'NI')

    controller.main_win.update_status(workflow, proc_wft, controller.main_win.cell_zones_status)
            