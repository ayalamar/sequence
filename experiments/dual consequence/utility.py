from math import sqrt, cos, sin, pi
import csv
from numpy import matrix, matmul, linspace, poly1d
from random import shuffle


def counterbalance(values, size):
    lst = []
    shuffle(values)
    for i in range(0, size / len(values)):
        shuffle(values)
        for item in values:
            lst.append(item)
    return lst


def find_diff(a, b):
    return sqrt((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2)


def realign(cfg, rotation, trial_data):
    if trial_data['clamp'] is 'off':
        return rotator(cfg, rotation, trial_data['hand_pos'][-1])
    elif trial_data['clamp'] is 'horizontal_clamp':
        hmove = trial_data['hand_pos'][-1][0] - trial_data['hand_pos'][-2][0]
        return [trial_data['cursor_pos'][-1][0] + hmove, trial_data['cursor_pos'][-1][1]]


def writedata(cfg, block, trial_data):
    outfilename = 'Data/' + cfg['data_folder'] + '/P' + str(cfg['P_number']).zfill(3) + '-' + str(block['number']) + '-' \
                  + str(trial_data['trial_number']) + '.txt'
    trial_number = trial_data['trial_number']
    with open(outfilename, 'wb') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(
            ['time', 'trial', 'cursorx', 'cursory', 'hand_x', 'hand_Y', 'rotation', 'step', 'targetangle',
             'targetposx', 'targetposy', 'garagelocation', 'garageposx',
             'garageposy', 'homex', 'homey', 'group', 'trial_type'])
        for index, (time, cursor, hand) in enumerate(
                zip(trial_data['time'], trial_data['cursor_pos'], trial_data['hand_pos'])):
            if block['garage_left_right'][trial_number] == 'left':
                output = -1
            else:
                output = 1
            if block['block_name'] == 'Reach with cursor':
                t_output = 1
            else:
                t_output = -1
            csv_writer.writerow(
                [time, trial_number, cursor[0], cursor[1], hand[0], hand[1], block['rotation'][trial_number],
                 trial_data['num_step'][index - 1], trial_data['target_angle'],
                 block['target_stim'][trial_number].pos[0],
                 block['target_stim'][trial_number].pos[1], output,
                 trial_data['garage_pos'][0], trial_data['garage_pos'][1],
                 cfg['home_stim'].pos[0], cfg['home_stim'].pos[1],
                 cfg['exp_name'], t_output])
        csvfile.close()


def rotator(cfg, angle, position):
    homex = cfg['homepos'][0]
    homey = cfg['homepos'][1]
    [x, y] = (position[0] - homex, position[1] - homey)
    angle = pi * (angle / 180.0)
    temp = matrix([[1, 2, 3]])
    temp = matmul(matrix([[x, y, 1]]), matrix([[cos(angle), -sin(angle), 0], [sin(angle), cos(angle), 0], [0, 0, 1]]))
    xx = temp.item(0)
    yy = temp.item(1) + homey
    return [xx, yy]


def speed(pos, time, limit, block_name):
    current_time = time[-1]
    first_time = next(x[0] for x in enumerate(time) if x[1] >= current_time - 0.2)
    if block_name == 'Reach with cursor':
        if find_diff(pos[first_time], pos[-1]) <= 0:
            return 1
        else:
            return 0
    else:
        if find_diff(pos[first_time], pos[-1]) <= 2:
            return 1
        else:
            return 0


def calibrate(data):
    x = data[0]
    y = data[1]

    # Aspect Ratio
    x = x / 1.6

    # Scale UP
    x = x * 1.15
    y = y * 1.15

    data = [x, y]
    return data


def add(x, y): return x + y


def steps(start, end, size):
    return poly1d(linspace(end, start, size))
