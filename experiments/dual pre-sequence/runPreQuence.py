import pandas as pd
import scipy as sp
import time
import os
from random import shuffle, seed
from psychopy import gui, visual, event, core, sound
from ctypes import *

def runExp():
	# set up config object
	cfg = {}  
	
	# ask for participant info
	cfg = getParticipant(cfg)
	
	# set up psychopy environment:
	cfg = openEnvironment(cfg)
	
	cfg = makeStimuli(cfg)
	
	# set random number seed (so participants have different trial orders)
	seed(cfg['ppno']*999)
	
	
	# generate block & trial order:
	cfg = makeBlocks(cfg)
	
	# loop through the blocks:
	for blockno in range(cfg['starting_block'], len(cfg['blocks'])):
		
		cfg['blockno'] = blockno
		
		# if the instruction is not "empty": ''
		#    show the instruction
		
		if len(cfg['blocks'][blockno]['instruction']) > 0:	
			cfg = showInstruc(cfg)
			
		
		# loop through trials:
		ntrials = cfg['blocks'][blockno]['trials'].shape[0]
		if (cfg['mode']=='test'):
			ntrials = 3
		
		for trialno in range(ntrials):
			cfg['trialno'] = trialno
			cfg = runTrial(cfg)
			
		# concatenate all trial data into one task data file?
		taskdata = 0
		
		for trialno in range(ntrials):
			
			filename = "p%03d-%d-%d.csv" % (cfg['ppno'], cfg['blockno'],trialno)
			filepath = os.path.join('data',filename)
			trialdata = pd.DataFrame.from_csv(path = filepath, index_col = None )
			
			if isinstance(taskdata, pd.DataFrame):
				taskdata = pd.concat([taskdata, trialdata])
			else:
				taskdata = trialdata
				
		filename = "p%03d-%d.csv" % (cfg['ppno'], cfg['blockno'])
		filepath = os.path.join('data',filename)
		taskdata.to_csv(filepath, index = False )
				


	cfg = closeEnvironment(cfg)
	
	
def runTrial(cfg):
	
	blockno = cfg['blockno']
	trialno = cfg['trialno']
	
	trialDef = cfg['blocks'][blockno]['trials'].iloc[trialno]
	cursorvisible = trialDef.CursorVis
	prehomeposition = trialDef.Prepos
	target = trialDef.Target
	rotation = trialDef.Rotation
	
	#cfg['home'].draw()

	theta = target * sp.pi/180
	targetx = sp.cos(theta) * cfg['targetdistance']
	targety = sp.sin(theta) * cfg['targetdistance']
	targetx = targetx + cfg['home'].pos[0]
	targety = targety + cfg['home'].pos[1]
	
	cfg['target'].pos = [targetx,targety]
	#cfg['target'].draw()
	
	# rotate, position & then draw prehome
	
	# prehomeposition:
	# -1: left, ori=0
	#  1: right, ori=180
	
	cfg['preHome'].pos = [trialDef.Prepos*cfg['prehomedistance'],cfg['home'].pos[1]] 
	prehomepos = trialDef.Prepos*cfg['prehomedistance'],cfg['home'].pos[1]
	cfg['preHome'].ori = 180*((trialDef.Prepos + 1)/2)
	# cfg['preHome'].draw()
	
	cfg['arrow'].ori = 180*((trialDef.Prepos + 1)/2)
	
	
	#cfg['cursor'].pos = cfg['mouse'].getPos()
	#cfg['cursor'].draw()
	
	#cfg['win'].flip()
	
	# steps:
	# -2: reach home
	# -1: reach pre-home
	#  0: hold home
    #  1: leave home
	#  2: reach target
	
	#  5: reach home
	
	# set up variables for trial control:
	step = -2
	
	# create vectors to collect data samples:
	cursorx_px = []
	cursory_px = []
	mousex_px  = []
	mousey_px  = []
	time_ms    = []
	step_vec   = []
	
	while (step < 6):
		
		# where would the cursor be?
		mousepos = cfg['mouse'].Pos()
		
		# get cursor position
		cursorpos = mousepos[0:2]
		
		# rotate cursor position
		#if trialDef.Rotation != 0:
		#	cursorpos = rotatepos(point = cursorpos, ref = cfg['home'].pos, rotation = trialDef.Rotation)
		
		cfg['cursor'].pos = cursorpos

		
		if (step == -2):
			# participant needs to go to the home position
			# usually already true, except on the first trial of a block...
			
			# show home position
			cfg['home'].draw()
			
			if cfg['mode'] == 'test':
				cfg['cursor'].draw()
			
			# if distance between cursor and home position is lower than XYZ, move on to step -1 
			# a = cursorx - homex ; b = cursory - homey ; c^2 = a^2 + b^2; c = distance between cursor and home
			if sp.sqrt((cursorpos[0] - cfg['home'].pos[0])**2 + (cursorpos[1] - cfg['home'].pos[1])**2) < cfg['radius']:
				step = -1
				
			
		if (step == -1):
			# participants needs to go to the point to start the lead-in movement
			
			# show home
			cfg['home'].draw()
			# show arrow
			cfg['arrow'].draw()
			
			if cfg['mode'] == 'test':
				cfg['cursor'].draw()
				cfg['preHome'].draw()
			
			# check if participant stayed "clamped", if not, restart trial (step = -2), and BEEP
			if mousepos[1] > (cfg['home'].pos[1] + cfg['radius']):
				step = -2
				cfg['sound'].play()
			
			# if distance between cursor and pre-home is lower than XYZ, move on to step 0
			if sp.sqrt((mousepos[0] - prehomepos[0])**2 + (mousepos[1] - prehomepos[1])**2) < cfg['radius']:
				step = 0
			
			
		if (step == 0):
			# participants need to make lead-in movement to home
			# show home
			cfg['home'].draw()
			# show cursor 
			cfg['cursor'].draw()
			
			# stay clamped during lead-in movement
			if mousepos[1] > (cfg['home'].pos[1] + cfg['radius']):
				step = -2
				cfg['sound'].play()
			
			# if distance between cursor and home is lower than XYZ for 100 ms, move on to step 1
			
			# determine which samples to use:
			sample_idx = sp.array([s >= (time_ms[-1]-100) for s in time_ms]).nonzero()[0]
			curX = sp.array(cursorx_px)[sample_idx] - cfg['home'].pos[0]
			curY = sp.array(cursory_px)[sample_idx] - cfg['home'].pos[1]
			
			maxdist = max((curX**2 + curY**2)**0.5)

			if maxdist < cfg['radius']:
				step = 1		
	
		if (step == 1):
			# participants see target and leave the home position
			# rotate cursor position
			if trialDef.Rotation != 0:
				cursorpos = rotatepos(point = cursorpos, ref = cfg['home'].pos, rotation = trialDef.Rotation)
			cfg['cursor'].pos = cursorpos
			
			# show target
			cfg['target'].draw()
			# show cursor
			cfg['cursor'].draw()
			
			# if distance between cursor and home is higher than XYZ, move on to step 2
			if sp.sqrt((cursorpos[0] - cfg['home'].pos[0])**2 + (cursorpos[1] - cfg['home'].pos[1])**2) > cfg['radius']:
				step = 2
				
		if (step == 2):
			# participant needs to reach to target (or stop moving for no cursor blocks)
			# rotate cursor position
			if trialDef.Rotation != 0:
				cursorpos = rotatepos(point = cursorpos, ref = cfg['home'].pos, rotation = trialDef.Rotation)
			cfg['cursor'].pos = cursorpos

			# show target
			cfg['target'].draw()

			#   if distance between cursor and target is lower than XYZ, move on to step 5
			if trialDef.CursorVis:
				cfg['cursor'].draw()
				if sp.sqrt((cursorpos[0] - cfg['target'].pos[0])**2 + (cursorpos[1] - cfg['target'].pos[1])**2) < cfg['radius']:
					step = 3
					
			else:
				if sp.sqrt((cursorpos[0] - cfg['home'].pos[0])**2 + (cursorpos[1] - cfg['home'].pos[1])**2) > (cfg['targetdistance']/2.):
					sample_idx = sp.array([s >= (time_ms[-1]-100) for s in time_ms]).nonzero()[0]
					curX = sp.array(cursorx_px)[sample_idx]
					curY = sp.array(cursory_px)[sample_idx]
					
					distances = ((sp.diff(curX)**2 + sp.diff(curY)**2)**0.5) / cfg['pixpercm']
					timediffs = sp.diff(sp.array(time_ms)[sample_idx])/1000
					
					#print(distances)
					#print(timediffs)
					#print(distances / timediffs)
					
					avg_speed = sp.mean(distances / timediffs)
					
					#print('speed: %0.2f -- criterion: %0.2f' %(avg_speed,1*cfg['pixpercm']))
					
					# speed below 1 cm/s ?
					if avg_speed < 1:
						step = 3
				#else:
					#print('distance: %0.2f -- criterion: %0.2f' %(sp.sqrt((cursorpos[0] - cfg['home'].pos[0])**2 + (cursorpos[1] - cfg['home'].pos[1])**2),cfg['targetdistance']/2.))
				
				if cfg['mode'] == 'test':
					cfg['cursor'].draw()
			
		if (step == 3):
			step = 4
		if (step == 4):
			step = 5
			
		if (step == 5):
			
			# show home
			cfg['home'].draw()
			
			if cfg['mode'] == 'test':
				cfg['cursor'].draw()
			
			# if distance to home < XYZ: show cursor
			
			
			# if distance between cursor and home is lower than XYZ, set step to 6
			if sp.sqrt((cursorpos[0] - cfg['home'].pos[0])**2 + (cursorpos[1] - cfg['home'].pos[1])**2) > cfg['radius']:
				step = 6
			
		cfg['win'].flip()
		# add data to mouse X / Y vectors, cursor X / Y vectors, time vector
		cursorx_px.append(cursorpos[0]) 
		cursory_px.append(cursorpos[1])
		mousex_px.append(mousepos[0]) 
		mousey_px.append(mousepos[1])
		time_ms.append(mousepos[2]*1000)
		step_vec.append(step)
		

		
	# #####################################
	# convert all pixel data to centimeters!
	# #####################################
	
	d = {'step': step_vec, 
		'time_ms': time_ms, 
		'mousex_px': mousex_px,
		'mousey_px': mousey_px,
		'cursorx_px': cursorx_px,
		'cursory_px':cursory_px,
		'homex_px': cfg['home'].pos[0],
		'homey_px': cfg['home'].pos[1],
		'rotation': trialDef.Rotation,
		'participant': cfg['ppno'],
		'targetangle_deg': trialDef.Target,
		'targetx_px': cfg['target'].pos[0],
		'targety_px': cfg['target'].pos[1],
		'trial': cfg['trialno'],
		'prehome': trialDef.Prepos,
		'prehomex_px': prehomepos[0],
		'prehomey_px':prehomepos[1]}
	
	trialdata = pd.DataFrame(d)
	
	filename = "p%03d-%d-%d.csv" % (cfg['ppno'], cfg['blockno'],cfg['trialno'])
	filepath = os.path.join('data',filename)
	
	trialdata.to_csv(filepath, index=False)
	
	# buf = "A = %d\n , B = %s\n" % (a, b)
	# p002-0-0.csv
	return(cfg)



def getParticipant(cfg):
	expInfo = {'1) participant number':' ', '2) starting block': '0', '3) mode': 'real'}
	dlg = gui.DlgFromDict(dictionary = expInfo, title = 'welcome to experiment')
	cfg['ppno'] = int(expInfo['1) participant number'])
	cfg['starting_block'] = int(expInfo['2) starting block'])
	cfg['mode'] = expInfo['3) mode']
	return cfg


def showInstruc(cfg):
	# print(cfg['blocks'][cfg['blockno']]['instruction'])
	# show this on screen!
	
	instruction = visual.TextStim(cfg['win'], text = cfg['blocks'][cfg['blockno']]['instruction'], height = 16)
	waitingForSpace = True
	
	while waitingForSpace: 
		pressed = event.getKeys(keyList = ['space'])
		instruction.draw()
		cfg['win'].flip()
		if len(pressed) > 0:
			cfg['win'].flip()
			waitingForSpace = False

	
	return(cfg)


def openEnvironment(cfg):
	
	# 47.4 * 29.6 
	
	# 1680 / 47.4 ~ 34.443 pix per cm
	# 1050 / 29.6 ~ 35.4729 pix per cm
	
	# let's say 35 pix / cm?
	
	cfg['pixpercm'] = 35
	
	winSize = [640,480]
	#winSize = [1680, 1050]
	
	cfg['win'] = visual.Window(size = winSize, color =(0,0,0), units ='pix', fullscr=False)
	#cfg['win'] = visual.Window(size = [1680, 1050], color =(0,0,0), units ='pix', fullscr=True)
	
	cfg['winSize'] = winSize
	
	cfg['psyMouse'] = event.Mouse(visible = False, newPos = None, win = cfg['win'])
	
	cfg = addMouse(cfg)
	
	return(cfg)

def addMouse(cfg):
	
	# the X coordinate is scaled, so that on the 16:9 widescreen monitor, and square tablet, reaches are still proportional
	# factor: 1.6
	# both X and Y are then scaled up, so that the movement in centimeters is equal on the screen and tablet
	# factor: 1.05
	
	try:
		class myMouse:
			Xlib = CDLL("libX11.so.6")
			display = Xlib.XOpenDisplay(None)
			if display == 0: sys.exit(2) # no display or can't be accessed...
			w = Xlib.XRootWindow(display, c_int(cfg['monitorIndex']-1))
			(root_id, child_id) = (c_uint32(), c_uint32())
			(root_x, root_y, win_x, win_y) = (c_int(), c_int(), c_int(), c_int())
			mask = c_uint()
			      
    	  	def Pos(self):
    	  		#print('X11 mouse')
    	  		ret = self.Xlib.XQueryPointer(self.display, c_uint32(self.w), byref(self.root_id), byref(self.child_id), byref(self.root_x), byref(self.root_y), byref(self.win_x), byref(self.win_y), byref(self.mask))
    	  		if ret == 0: sys.exit(1)
    	  		return [(self.root_x.value - (cfg['width']/2))/1.6*1.05, -1 * (self.root_y.value - (cfg['height']/2))*1.05, time.time()] # c_int can't be used by regular Python to do math, but the values of c_ints are ints - also, we return the current time
      	  
	except:
		# Xlib is not available (not on Linux?)
 	 	# use Psychopy:
 	 	class myMouse:
  	  	  
 		    def Pos(self):
 		    	#print('PsychoPy mouse')
 		     	[X,Y] = cfg['psyMouse'].getPos()
 		     	return [X/1.6*1.05,Y*1.05,time.time()]
      	  

    cfg['mouse'] = myMouse()
  	
	return(cfg)


def closeEnvironment(cfg):
	
	cfg['psyMouse'].setVisible(True)
	cfg['win'].close() 
	return(cfg)
	
	
def makeStimuli(cfg):
	
	radius = 0.0375*2/3*cfg['winSize'][1]
	cfg['radius'] = radius
	
	# add home position
	cfg['home'] = visual.Circle(win = cfg['win'], pos = [0,-0.25*cfg['winSize'][1]], radius = radius)
	cfg['home'].setFillColor(color=(255,0,0), colorSpace = 'rgb255')
	cfg['home'].setLineColor(color=(255,0,0), colorSpace = 'rgb255')
	
	
	# add target
	cfg['target'] = visual.Circle(win = cfg['win'], pos = [0,-0.25*cfg['winSize'][1]], radius = radius)
	cfg['target'].setFillColor(color=(255,255,0), colorSpace = 'rgb255')
	cfg['target'].setLineColor(color=(255,255,0), colorSpace = 'rgb255')
	
	if (cfg['winSize'][1] == 1050):
		cfg['targetdistance'] = cfg['pixpercm'] * 12
	else:
		cfg['targetdistance'] = cfg['pixpercm'] * 5
	
	# add cursor
	cfg['cursor'] = visual.Circle(win = cfg['win'], pos = [0,-0.25*cfg['winSize'][1]], radius = radius)
	cfg['cursor'].setFillColor(color=(93,220,245), colorSpace = 'rgb255')
	cfg['cursor'].setLineColor(color=(93,220,245), colorSpace = 'rgb255')
	
	# add pre-home target
	cfg['preHome'] = visual.ShapeStim( win=cfg['win'], vertices= ((.5,.5),(-.5,.5),(-.5,-.5),(.5,-.5)) ,closeShape = False, size = 2*cfg['pixpercm'])
	
	# add arrow on home position 
	cfg['arrow'] = visual.ShapeStim( win=cfg['win'], vertices= ((-.5,0),(.5,.5),(.25,0),(.5,-.5)) ,closeShape = True, size = radius*1.5, pos = cfg['home'].pos)
	cfg['arrow'].setFillColor(color=(0,0,0), colorSpace = 'rgb')
	
	if (cfg['winSize'][1] == 1050):
		cfg['prehomedistance'] = cfg['pixpercm'] * 7.5
	else:
		cfg['prehomedistance'] = cfg['pixpercm'] * 3.5
		
	# add BEEP
	sound.init(rate=44100, stereo=True, buffer=128)
	cfg['sound'] = sound.Sound('ding.wav', secs = 1)
		
	return(cfg)

def makeBlocks(cfg):
	
	# empty list with blocks to do:
	blocks = []
	
	# set properties of blocks:
	cursorVisibilityList = [1,0,1,1,0,1,0,1,0,0,1,0,0]
	rotated = [0,0,0,1,0,1,0,1,0,0,1,0,0]
	preposind = [0,0,0,0,0,0,0,0,0,0,0,0,0]
	
	trialSets = [10,4,1,60,2,60,2,4,2,2,4,2,2]
	trialTargets = sp.array([75,90,105,75,90,105])
	trialRotations = sp.array([30,30,30,-30,-30,-30])
	trialPrepos = [sp.array([1,1,1,-1,-1,-1]), sp.repeat(1,6), sp.repeat(-1,6)]
	
	
	# counter-balance include(1)/exclude(0) strategy:
	orderIndex = cfg['ppno'] % 4 # either 0 or 1
	instrorder  = [[0,1],[1,0],[0,1],[1,0]][orderIndex]
	preposorder = [[1,2],[2,1],[2,1],[1,2]][orderIndex]
	
	instruct = ['reach with cursor',
				'reach without cursor',
				'reach with cursor',
	            '',
	            'reach without cursor',
				'reach with cursor',
				'reach without cursor',
				'reach with cursor',
		        ['reach without cursor (include strategy)', 'reach without cursor (exclude strategy)'][instrorder[0]],
		        ['reach without cursor (include strategy)', 'reach without cursor (exclude strategy)'][instrorder[1]],
		        'reach with cursor',
		        ['reach without cursor (include strategy)', 'reach without cursor (exclude strategy)'][instrorder[1]],
		        ['reach without cursor (include strategy)', 'reach without cursor (exclude strategy)'][instrorder[0]]
		        
	]
	
	preposind[8] = preposorder[0]
	preposind[9] = preposorder[1]
	preposind[11] = preposorder[1]
	preposind[12] = preposorder[0]
	
	for blockno in range(13): 
		
		blockDef = {}
		blockDef['instruction'] = instruct[blockno]
		
		target = []
		cursorVisibility = []
		rotation = []
		prepos = []
		
		for trialSet in range(trialSets[blockno]):
			# there are 6 trial types
			# that have to appear in random order in each trial set
			# each target appears once before repeating
			
			# first and second set of three trials each have all three targets:
			first = []
			second = []
			
			
			for indices in range(3):
				trialindices = [[0,3],[1,4],[2,5]][indices]
				
				# but with the lead-in randomized
				shuffle(trialindices)
				
				# the first set of three gets one lead-in
				first.append(trialindices[0])
				# the second gets the other
				second.append(trialindices[1])
				
			# within each three trials we shuffle the order: 
			shuffle(first)
			shuffle(second)
			
			# and all 6 determine 6 trials in this trialSet
			trialInd = first + second
			
			# use the random trial order to append trial descriptors to the block lists:
			target = target + list(trialTargets[trialInd])
			rotation = rotation + list(trialRotations[trialInd] * rotated[blockno])
			prepos = prepos + list(trialPrepos[preposind[blockno]][trialInd])
		
			cursorVisibility = cursorVisibility + list(sp.repeat(cursorVisibilityList[blockno], 6))
			
		# convert the block lists to a dataframe:
		trialDefs = pd.DataFrame({'Target':target,'Rotation':rotation,'Prepos':prepos,'CursorVis':cursorVisibility})
		blockDef['trials'] = trialDefs
		# and add it to the list of block definitions:
		blocks.append(blockDef)
		
	cfg['blocks'] = blocks
	return(cfg)


def rotatepos(point, ref=[0,0], rotation = 0):
	
	t = rotation*sp.pi/180
	
	R = sp.array([[sp.cos(t), -1*sp.sin(t)], [sp.sin(t), sp.cos(t)]])
	
	point = sp.array(point) - sp.array(ref)
	
	return(list(sp.dot(R,point) + sp.array(ref)))
	
	
	


runExp()
