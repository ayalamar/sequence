# block structure
def initexp(cfg):

    #LOGICAL_CONFIG
    cfg['exp_name'] = '2'
    cfg['num_blocks'] = 10
    cfg ['block_names'] = ['Reach with cursor', 'Reach without cursor', 'Reach with cursor', 'Reach with cursor',
                           'Reach without cursor', 'Reach with cursor', 'Reach without cursor','Reach with cursor', 'Reach without cursor',  'Reach without cursor']
    cfg['trial_numbers'] = [60,12,6,180,12,180,12,12,12,12]
    cfg['rotations'] = [0, 0, 0,30, 0,30, 0, 30, 0,0,30, 0,0]
    cfg['garage_distance'] = ['far', 'far', 'far','far', 'far', 'far', 'far', 'far', 'far', 'far']
    cfg['garage_left_right'] = ['left', 'left', 'left','left', 'left','left', 'left', 'left', 'left','left']
    cfg['target_angles'] = [75, 105, 90]

    #DISPLAY_CONFIG
    cfg['target_distance'] = 12
    cfg['target_radius'] = 0.0375 * 2 / 3  #relative to the height of screen
    cfg['close_garage_distance'] = 0.1 #relative to the width of screen
    cfg['far_garage_distance'] = 9
    cfg['visible_home_range'] = 2
    cfg['acquire_target_limit'] = cfg['target_radius']  #relative to the height of screen
    cfg['speed_limit'] = 9
    cfg['acquire_target_time'] = 50
    cfg['acquire_garage_limit'] = 2
    cfg['garage_size'] = 3

    return cfg

