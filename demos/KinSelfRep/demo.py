from lxml import etree
import subprocess as sub
import numpy as np


RECORD_HISTORY = True
DRAW_LINKS = False
DRAW_WALLS = False

SWARM_SIZE = 4  # if this is greater than 4 then add more coordinates to (lines 103 and 104)

WORLD_SIZE = 107
BODY_SIZE = (8, 8, 7)
# if body size changes, or if the stiffness/density of body material changes, 
# then the cilia force of the material will need to be recalibrated

MIN_BOT_SIZE = 64
EVAL_PERIOD = 4
SETTLE_TIME = 0.5

SPACE_BETWEEN_DEBRIS = 2
DEBRIS_MAT = 2
REPLENISH_DEBRIS_EVERY = EVAL_PERIOD + SETTLE_TIME

ATTACH_WATCH_DISTANCE = 1.25
ATTACH_BOUNDING_RADIUS = 0.8

SEED = 1
np.random.seed(SEED)

bx, by, bz = BODY_SIZE
wx, wy, wz = (WORLD_SIZE, WORLD_SIZE, bz)

spacing_between_bodies = int(4.25*bx)-1  
# if spacing is too big compared to worldsize then you'll get an error when trying to broadcast the body onto the world

# controller
BASE_CILIA_FORCE = np.ones((wx, wy, wz, 3))  * -1  # pointing downward
BASE_CILIA_FORCE[:, :, :, :2] = 2 * np.random.rand(wx, wy, wz, 2) - 1  # initial forces

# morphology
BODY_PLAN = np.ones(BODY_SIZE, dtype=np.int)

# create data folder if it doesn't already exist
sub.call("mkdir data{}".format(SEED), shell=True)
sub.call("cp base.vxa data{}/base.vxa".format(SEED), shell=True)

# clear old .vxd robot files from the data directory
sub.call("rm data{}/*.vxd".format(SEED), shell=True)

# remove old sim output.xml if we are saving new stats
if not RECORD_HISTORY:
    sub.call("rm output{}.xml".format(SEED), shell=True)

# start vxd file
root = etree.Element("VXD")

vxa_min_bot_size = etree.SubElement(root, "MinimumBotSize")
vxa_min_bot_size.set('replace', 'VXA.Simulator.MinimumBotSize')
vxa_min_bot_size.text = str(MIN_BOT_SIZE)

vxa_debris_spacing = etree.SubElement(root, "SpaceBetweenDebris")
vxa_debris_spacing.set('replace', 'VXA.Simulator.SpaceBetweenDebris')
vxa_debris_spacing.text = str(SPACE_BETWEEN_DEBRIS)

vxa_replenish_debris_every = etree.SubElement(root, "ReplenishDebrisEvery")
vxa_replenish_debris_every.set('replace', 'VXA.Simulator.ReplenishDebrisEvery')
vxa_replenish_debris_every.text = str(REPLENISH_DEBRIS_EVERY)

vxa_eval_period = etree.SubElement(root, "ReinitializeInitialPositionAfterThisManySeconds")
vxa_eval_period.set('replace', 'VXA.Simulator.ReinitializeInitialPositionAfterThisManySeconds')
vxa_eval_period.text = str(EVAL_PERIOD)

vxa_settle_time = etree.SubElement(root, "SettleTimeBeforeNextRoundOfReplication")
vxa_settle_time.set('replace', 'VXA.Simulator.SettleTimeBeforeNextRoundOfReplication')
vxa_settle_time.text = str(SETTLE_TIME) 

vxa_debris_mat = etree.SubElement(root, "DebrisMat")
vxa_debris_mat.set('replace', 'VXA.Simulator.DebrisMat')
vxa_debris_mat.text = str(DEBRIS_MAT)

vxa_world_size = etree.SubElement(root, "WorldSize")
vxa_world_size.set('replace', 'VXA.Simulator.WorldSize')
vxa_world_size.text = str(wx)

attach_detach = etree.SubElement(root, "AttachDetach")
attach_detach.set('replace', 'VXA.Simulator.AttachDetach')
etree.SubElement(attach_detach, "watchDistance").text = str(ATTACH_WATCH_DISTANCE)
etree.SubElement(attach_detach, "boundingRadius").text = str(ATTACH_BOUNDING_RADIUS)

# set seed for browain cilia motion
vxa_seed = etree.SubElement(root, "RandomSeed")
vxa_seed.set('replace', 'VXA.Simulator.RandomSeed')
vxa_seed.text = str(SEED)


if RECORD_HISTORY:
    # sub.call("rm a{0}_gen{1}.hist".format(seed, pop.gen), shell=True)
    history = etree.SubElement(root, "RecordHistory")
    history.set('replace', 'VXA.Simulator.RecordHistory')
    etree.SubElement(history, "RecordStepSize").text = '100'
    etree.SubElement(history, "RecordVoxel").text = '1'
    etree.SubElement(history, "RecordLink").text = str(int(DRAW_LINKS))
    etree.SubElement(history, "RecordFixedVoxels").text = str(int(DRAW_WALLS))
    etree.SubElement(history, "RecordCoMTraceOfEachVoxelGroupfOfThisMaterial").text = '0'  # draw CoM trace=
    

structure = etree.SubElement(root, "Structure")
structure.set('replace', 'VXA.VXC.Structure')
structure.set('Compression', 'ASCII_READABLE')
etree.SubElement(structure, "X_Voxels").text = str(wx)
etree.SubElement(structure, "Y_Voxels").text = str(wy)
etree.SubElement(structure, "Z_Voxels").text = str(wz)


world = np.zeros((wx, wy, wz), dtype=np.int8)

bodies = [BODY_PLAN] * SWARM_SIZE

s = spacing_between_bodies
a = [s, 2*s, s, 2*s]
b = [s, 2*s, 2*s, s]

for n, (ai, bi) in enumerate(zip(a,b)):
    try:  
        world[ai:ai+bx, bi:bi+by, :] = bodies[n]
    except IndexError:
        pass

world = np.swapaxes(world, 0,2)
# world = world.reshape([wz,-1])
world = world.reshape(wz, wx*wy)

for i in range(wx):
    for j in range(wy):
        if (i == 0) or (j == 0) or (i == wx-1) or (j == wy-1):
            world[:, i*wx+j] = 4  # wall

for i in range(2, wx, SPACE_BETWEEN_DEBRIS+1): 
    for j in range(2, wy, SPACE_BETWEEN_DEBRIS+1):
        for k in range(1): # DEBRIS_HEIGHT = 1
            try:
                if ((world[k, i*wx+j] == 0) 
                and (world[k-1, i*wx+j] in [0, DEBRIS_MAT]) and (world[k+1, i*wx+j] in [0, DEBRIS_MAT])
                and (world[k, (i+1)*wx+j] == 0) and (world[k, (i-1)*wx+j] == 0) 
                and (world[k, i*wx+j+1] == 0) and (world[k, i*wx+j-1] == 0) ):

                    world[k, i*wx+j] = DEBRIS_MAT  # pellet

            except IndexError:
                pass

data = etree.SubElement(structure, "Data")
for i in range(world.shape[0]):
    layer = etree.SubElement(data, "Layer")
    str_layer = "".join([str(c) for c in world[i]])
    layer.text = etree.CDATA(str_layer)

# cilia motion
base_cilia_force = np.swapaxes(BASE_CILIA_FORCE, 0,2)
base_cilia_force = base_cilia_force.reshape(wz, 3*wx*wy)

data = etree.SubElement(structure, "BaseCiliaForce")
for i in range(base_cilia_force.shape[0]):
    layer = etree.SubElement(data, "Layer")
    str_layer = "".join([str(c) + ", " for c in base_cilia_force[i]])
    layer.text = etree.CDATA(str_layer)

# save the vxd to data folder
with open('data'+str(SEED)+'/bot_0.vxd', 'wb') as vxd:
    vxd.write(etree.tostring(root))

# ok let's finally evaluate all the robots in the data directory

if RECORD_HISTORY:
    sub.call("./voxcraft-sim -i data{0} > a.hist".format(SEED), shell=True)
else:
    sub.call("./voxcraft-sim -i data{0} -o output{0}.xml".format(SEED), shell=True)
