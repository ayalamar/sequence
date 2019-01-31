import collections
import os
from utility import counterbalance, realign, find_diff, speed, writedata, calibrate, steps
from psychopy import event, core
from gui import targetstim, garagestim, instructions

def initblocks(cfg):
    blocks = []

    for block in range(cfg['starting_block'], cfg['num_blocks']):
        blockinfo = collections.defaultdict(list)
        blockinfo['filename'] = cfg['exp_name'] + os.path.sep + 'p%03d' % int(cfg['P_number']) + str(block)
        blockinfo['block_name'] = cfg['block_names'][block]
        blockinfo['num_trials'] = cfg['trial_numbers'][block]
        blockinfo['target_angles'] = counterbalance(cfg['target_angles'], cfg['trial_numbers'][block])
        blockinfo['target_stim'] = []
        if cfg['garage_distance'][block] is 'both':
            blockinfo['garage_distance'] = counterbalance(['far', 'close'], cfg['trial_numbers'][block])
        else:
            blockinfo['garage_distance'] = [cfg['garage_distance'][block]] * cfg['trial_numbers'][block]

        if cfg['garage_left_right'][block] is 'both':
            blockinfo['garage_left_right'] = counterbalance(['left', 'right'], cfg['trial_numbers'][block])
        else:
            blockinfo['garage_left_right'] = [cfg['garage_left_right'][block]] * cfg['trial_numbers'][block]

        for trial_number, garage_location in enumerate(blockinfo['garage_left_right']):
            blockinfo['target_stim'].append(targetstim(cfg, blockinfo['target_angles'][trial_number]))
            blockinfo['garage_stim'].append(garagestim(cfg, blockinfo['target_stim'][trial_number].pos,
                                                            blockinfo['garage_distance'][trial_number],
                                                            garage_location))

            if garage_location is 'left':
                blockinfo['rotation'].append(cfg['rotations'][block] * -1)

            elif garage_location is 'right':
                blockinfo['rotation'].append(cfg['rotations'][block] * 1)

        blockinfo['number'] = block

        blocks.append(blockinfo)

    return blocks

def runblock(block, cfg):
    print(block['number'])
    if block['number'] != 3:
       instructions(cfg, block['block_name'])
    for trial_number in range(0, block['num_trials']):
        try:
            runtrial(cfg, block, trial_number)
        except:
            exit()


def runtrial(cfg, block, trial_number):

    #Global Variables
    trial_data = collections.defaultdict(list)
    trial_data['trial_number'] = trial_number
    trial_data['hand_pos'] = [cfg['mouse'].getPos()]
    trial_data['hand_pos'][0] = calibrate(trial_data['hand_pos'][0])
    trial_data['step'] = 'start'
    trial_data['cursor_target_diff'] = 10000
    trial_data['cursor_home_diff'] = 0
    trial_data['garage_home_diff'] = find_diff(cfg['homepos'], block['garage_stim'][trial_number][2].end)
    trial_data['core_start'] = core.getTime()
    trial_data['time'] = [core.getTime() - trial_data['core_start']]
    trial_data['clamp'] = 'off'
    trial_data['cursor_pos'] = [realign(cfg, block['rotation'][trial_number], trial_data)]
    trial_data['left_target'] = 0

    #Visual settings
    trial_data['target_angle'] = block['target_angles'][trial_number]
    trial_data['target_pos'] = block['target_stim'][trial_number].pos
    trial_data['garage_pos'] = block['garage_stim'][trial_number][2].end
    trial_data['clamp_movement'] = steps(trial_data['target_pos'][0], trial_data['garage_pos'][0],100)


    cfg['mouse'].setVisible(False)
    cfg['home_stim'].draw()
    cfg['cursorstim'].draw()
   # cfg['arrow'].draw()

    block['target_stim'][trial_number].draw()
    for stim in block['garage_stim'][trial_number]:
        stim.draw()

    cfg['win'].flip()
    trial_data['num_step'] =[]

    while trial_data['step'] is not 'end':
        response = event.getKeys()
        if 'escape' in response:
            cfg['win'].close()
            exit()

        cfg['cursorstim'].pos = trial_data['cursor_pos'][-1]
        cfg['home_stim'].draw()
        cfg['cursorstim'].draw()


        if trial_data['step'] is not 'start':
            block['target_stim'][trial_number].draw()

        for stim in block['garage_stim'][trial_number]:
            stim.draw()

        cfg['win'].flip()
        trial_data['hand_pos'].append(calibrate(cfg['mouse'].getPos()))
        trial_data['time'].append(core.getTime() - trial_data['core_start'])
        trial_data['cursor_pos'].append(realign(cfg, block['rotation'][trial_number], trial_data))
        trial_data['cursor_home_diff'] = find_diff(cfg['homepos'], trial_data['cursor_pos'][-1])
        #if trial_data['left_target'] == 1:
        #    trial_data['cursor_target_diff'] = find_diff(trial_data['moved_goalposts'],trial_data['cursor_pos'][-1])
        #else:
        trial_data['cursor_target_diff'] = find_diff(trial_data['target_pos'],trial_data['cursor_pos'][-1])
        trial_data['cursor_garage_diff'] = find_diff(block['garage_stim'][trial_number][2].end,
                                                     trial_data['cursor_pos'][-1])
        trial_data['target_garage_diff'] = find_diff(block['target_stim'][trial_number].pos,
                                                     block['garage_stim'][trial_number][2].end)



        if trial_data['step'] is 'start':
            trial_data['num_step'].append(1)
            if trial_data['trial_number'] != 0:
                trial_data['step'] = 'leave_home'
            else:
                if cfg['mouse'].getPressed()[0] == 1:
                    if trial_data['cursor_home_diff'] < cfg['acquire_target_limit'] * cfg['height']:
                        trial_data['step'] = 'leave_home'

        elif trial_data['step'] is 'leave_home':
            trial_data['num_step'].append(1)
            if trial_data['cursor_home_diff'] > cfg['visible_home_range']:
                    trial_data['step'] = 'reach_target'
                    trial_data['reach_target_timer_start'] = trial_data['time'][-1]
                    cfg['home_stim'].opacity = 0

        elif trial_data['step'] is 'reach_target':
            trial_data['num_step'].append(2)
            block['target_stim'][trial_number].pos = trial_data['target_pos']
            #block['target_stim'][trial_number].opacity = 1
            if trial_data['time'][-1] >= trial_data['reach_target_timer_start'] + 1.5 and trial_data['left_target']  ==0:
                cfg['textstim'].draw()

            if trial_data['left_target'] == 1:
                cfg['slowstim'].draw()

            if block['block_name'] == 'Reach without cursor':
                cfg['cursorstim'].opacity = 0

            if acquire_target(cfg, block, trial_data):
                trial_data['step'] = 'waiting_on_target'
                trial_data['cursor_astarget_position'] = trial_data['cursor_pos'][-1]


        elif trial_data['step'] is 'waiting_on_target':
            trial_data['num_step'].append(3)
            if acquire_target(cfg, block, trial_data):
                trial_data['step'] = 'reach_garage'
                block['garage_stim'][trial_number] = color_garage(block['garage_stim'][trial_number], 'green')
                trial_data['hold_on_targetstart'] = trial_data['time'][-1]
                clock_position = 1

        elif trial_data['step'] is 'reach_garage':
                trial_data['num_step'].append(4)

                block['target_stim'][trial_number].opacity = 0

                if acquire_target(cfg, block, trial_data):
                    clock_position = int(1.75*(trial_data['time'][-1] - trial_data['hold_on_targetstart']) / 0.010)
                    new_pos = [trial_data['clamp_movement'][clock_position], trial_data['target_pos'][1]]
                    block['target_stim'][trial_number].setPos(new_pos)

                    if clock_position >= 90:
                        trial_data['step'] = 'reach_home'
                        cfg['cursorstim'].opacity = 0
                        cfg['home_stim'].opacity = 1
                else:
                    trial_data['step'] = 'reach_target'
                    cfg['sound'].play()
                    block['target_stim'][trial_number].pos = trial_data['target_pos']
                    block['target_stim'][trial_number].opacity = 1
                    trial_data['left_target'] = 1

        elif trial_data['step'] is 'reach_home':
            trial_data['num_step'].append(5)
            trial_data['clamp'] = 'off'
            if block['rotation'][trial_number] != 0 \
                    or (block is 'Reach without cursor' and \
                    trial_data['cursor_home_diff'] < 2/3 * trial_data['garage_home_diff']):
                        pass
                #cfg = set_arrow(cfg, trial_data['cursor_pos'][-1])
                #cfg['arrow'].draw()

            if trial_data['cursor_home_diff'] < trial_data['garage_home_diff']/5:
                cfg['cursorstim'].opacity = 1
            else:
                cfg['cursorstim'].opacity = 0



            if trial_data['cursor_home_diff'] < cfg['acquire_target_limit']*cfg['height']:
                trial_data['step'] = 'end'


    writedata(cfg, block, trial_data)




def acquire_target(cfg, block, trial_data):
    if block['block_name'] == 'Reach with cursor':
        if trial_data['step'] is 'reach_target' or trial_data['step'] is 'waiting_on_target' :
            return trial_data['cursor_target_diff'] < cfg['acquire_target_limit'] * cfg['height']*2 and\
                   speed(trial_data['cursor_pos'],trial_data['time'], cfg['speed_limit'],block['block_name'])

        if trial_data['step'] is 'reach_garage':
            return find_diff(trial_data['cursor_astarget_position'], trial_data['cursor_pos'][-1]) < cfg['target_radius'] * cfg['height']*2
    else:
        if trial_data['step'] is 'reach_target' or trial_data['step'] is 'waiting_on_target' :
            return trial_data['cursor_home_diff'] > 6 and \
               speed(trial_data['cursor_pos'], trial_data['time'], cfg['speed_limit'],block['block_name'])

        if trial_data['step'] is 'reach_garage':
            return find_diff(trial_data['cursor_astarget_position'], trial_data['cursor_pos'][-1]) < cfg['target_radius'] * cfg['height']*2



def acquire_garage(cfg, block, trial_data):
    if block['block_name'] == 'Reach with cursor':
        return trial_data['target_garage_diff'] < cfg['acquire_garage_limit']


def color_garage(garage, color):
    for stim in garage:
        stim.setLineColor(color=color, colorSpace='rgb255')
        stim.opacity = 1

    return garage


def set_arrow(cfg, cursor_pos):
    arrow_vector = [0, 0]
    arrow_vector[0] = cursor_pos[0] - cfg['homepos'][0]
    arrow_vector[1] = cursor_pos[1] - cfg['homepos'][1]
    arrow_size = find_diff(cursor_pos, cfg['homepos'])
    unit_arrow = arrow_vector / (arrow_size / cfg['home_stim'].radius)
    unit_arrow[1] = unit_arrow[1] + cfg['homepos'][1]
    cfg['arrow'].end = unit_arrow
    cfg['arrow'].opacity = 0
    return cfg
