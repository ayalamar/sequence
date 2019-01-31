from psychopy import visual, data, event, monitors, tools, gui, sound
import math

#Colors
red = (255,0,0)
green = (0,255,0)
blue = (0,0,255)
darkBlue = (0,0,128)
white = (255,255,255)
black = (0,0,0)
pink = (255,200,200)
lightgrey = (83, 79, 75)
darkgrey = (222, 212, 201)

def getparticipant(cfg):
    expName = cfg['exp_name']
    expInfo = {'participant initials:': '', 'participant number': '2', 'starting block': '0'}
    dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
    cfg['date'] = data.getDateStr()
    cfg['P_number'] = int(expInfo['participant number'])
    cfg['starting_block'] = int(expInfo['starting block'])
    cfg['expName'] = expName
    return cfg


def window(cfg):
    mon = monitors.Monitor('Ali', width=47)
    mon.setSizePix([1680, 1050])
    #win = visual.Window(size=[500 , 500], color=(0, 0, 0), units='cm', fullscr=False, monitor=mon)
    win = visual.Window(size=[1680 , 1050], color=(0, 0, 0), units='cm', fullscr=True, monitor=mon)
    cfg['win'] = win
    cmsize = tools.monitorunittools.pix2cm(win.size, mon)
    cfg['width'] = cmsize[0]
    cfg['height'] = cmsize[1]
    cfg['winpos'] = win.pos
    return cfg


def homestim(cfg):
    homepos = [0, -.25*cfg['height']]
    homestim = visual.Circle(win=cfg['win'], pos=homepos, radius=0.0375 * 2/3 * cfg['height'])
    homestim.setFillColor(color=red, colorSpace='rgb255')
    homestim.setLineColor(color=red, colorSpace='rgb255')
    cfg['homepos'] = homepos
    cfg['init_cursorpos'] = [0, cfg['homepos'][1]-10]
    cfg['home_stim'] = homestim
    return cfg


def cursorstim(cfg):
    homepos = [0, -.25*cfg['height']]
    cfg['cursorstim'] = visual.Circle(win=cfg['win'], pos=homepos, radius=0.0375 * 2/3 * cfg['height'])
    cfg['cursorstim'].setFillColor(color=(93, 220,245), colorSpace='rgb255')
    cfg['cursorstim'].setLineColor(color=(93, 220, 245), colorSpace='rgb255')
    cfg['cursorpos'] = [0, homepos[1]-10]
    cfg['cursorrecord'] = homepos
    return cfg


def soundstim(cfg):
    sound.init(rate=44100, stereo=True, buffer=128)
    cfg['sound'] = sound.Sound('ding.wav', secs = 1)
    return cfg


def mousestim(cfg):
    cfg['mouse'] = event.Mouse(newPos=cfg['homepos'], win=cfg['win'])
    return cfg


def arrowstim(cfg, endpoint):
    cfg['arrow'] = visual.Line(win=cfg['win'], start=cfg['homepos'], end=endpoint, opacity=0)
    cfg['arrow'].opacity = 0
    return cfg

def targetstim(cfg, angle):
    targetpos = [cfg['target_distance'] * math.cos(angle * (math.pi / 180)),
                 (cfg['target_distance'] * math.sin(angle * (math.pi / 180))) + cfg['homepos'][1]]

    targetstim = visual.Circle(win=cfg['win'], pos=targetpos, radius=cfg['target_radius'] * cfg['height'])
    targetstim.setFillColor(color=(255, 255, 0), colorSpace='rgb255')
    targetstim.setLineColor(color=(255, 255, 0), colorSpace='rgb255')
    return targetstim


def garagestim(cfg, target_position, garage_distance, garage_location):
    ##leftgarage


    if garage_distance is 'far':
        garage_distance = cfg['far_garage_distance']
    elif garage_distance is 'near':
        garage_distance = 0
    
    if garage_location is 'left':
        garage_distance = -1 * garage_distance
        garage_left_corner = -1

    else:
        garage_left_corner = 1

    garagetopright = [target_position[0] + garage_distance, target_position[1] + 1]
    garagetopleft = [garagetopright[0] - 2*garage_left_corner, target_position[1] + 1]
    garagebottomright = [target_position[0] + garage_distance, target_position[1] - 1]
    garagebottomleft = [garagebottomright[0] - 2*garage_left_corner, target_position[1] - 1]

    garage_top = visual.Line(win=cfg['win'], lineWidth=10000, start=(garagetopleft[0], garagetopleft[1]),
                             end=(garagetopright[0], garagetopright[1]))
    garage_bottom = visual.Line(win=cfg['win'], lineWidth=10000, start=(garagebottomleft[0], garagebottomleft[1]),
                                end=(garagebottomright[0], garagebottomright[1]))
    garage_left = visual.Line(win=cfg['win'], lineWidth=10000, start=(garagetopleft[0], garagetopleft[1]),
                              end=(garagebottomleft[0], garagebottomleft[1]))
    garage_right = visual.Line(win=cfg['win'], lineWidth=1000, start=(garagetopright[0], garagetopright[1]),
                               end=(garagebottomright[0], garagebottomright[1]))



    garage_left = garage_right

    return [garage_top, garage_bottom, garage_left]



def textstim(cfg):
     cfg['textstim'] = visual.TextStim(cfg['win'], text="TOO SLOW!", pos=(0,0), bold=True)
     cfg['slowstim'] = visual.TextStim(cfg['win'], text="STAY STILL!", pos=(0,0), bold=True)
     return cfg

def instructions(cfg, blockname):
    txt = visual.TextStim(cfg['win'], blockname)
    txt.setColor(white, 'rgb255')
    txt.draw()
    cfg['win'].flip()
    event.waitKeys(keyList=['space'])
