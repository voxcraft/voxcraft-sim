from lxml import etree
import subprocess as sub
import numpy as np


RECORD_HISTORY = True
DRAW_LINKS = False
DRAW_WALLS = False

SEED = 2
np.random.seed(SEED)

WORLD_SIZE = 68  # 107  # needs to be increased for SWARM_SIZE >25
BODY_SIZE = (6, 6, 5)  # (8, 8, 7)
# if body size changes, or if the stiffness/density of body material changes, 
# then the cilia force of the material will need to be recalibrated
bx, by, bz = BODY_SIZE
wx, wy, wz = (WORLD_SIZE, WORLD_SIZE, bz)

SWARM_SIZE = 2**2  # this should be a perfect square: 4, 9, 16, 25, etc.
ROWS = int(np.sqrt(SWARM_SIZE))
spacing_between_bodies = (wx - bx*ROWS) / (ROWS+1)

MIN_BOT_SIZE = 64 #16
EVAL_PERIOD = 4  # 4
SETTLE_TIME = 0.5  # 0.5

DEBRIS_HEIGHT = 2
HIGH_DEBRIS_CONCENTRATION = False

SPACE_BETWEEN_DEBRIS = 2 
DEBRIS_MAT = 2
REPLENISH_DEBRIS_EVERY = EVAL_PERIOD + SETTLE_TIME

RANDMONIZE_CILIA_EVERY = 1

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

# delete old hist file
sub.call("rm a.hist", shell=True)

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

vxa_debris_height = etree.SubElement(root, "DebrisHeight")
vxa_debris_height.set('replace', 'VXA.Simulator.DebrisHeight')
vxa_debris_height.text = str(DEBRIS_HEIGHT)

vxa_high_debris_concentration = etree.SubElement(root, "HighDebrisConcentration")
vxa_high_debris_concentration.set('replace', 'VXA.Simulator.HighDebrisConcentration')
vxa_high_debris_concentration.text = str(int(HIGH_DEBRIS_CONCENTRATION))

vxa_world_size = etree.SubElement(root, "WorldSize")
vxa_world_size.set('replace', 'VXA.Simulator.WorldSize')
vxa_world_size.text = str(wx)

vxa_randomize_cilia_every = etree.SubElement(root, "RandomizeCiliaEvery")
vxa_randomize_cilia_every.set('replace', 'VXA.Simulator.RandomizeCiliaEvery')
vxa_randomize_cilia_every.text = str(RANDMONIZE_CILIA_EVERY)

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
a = []
b = []
for r in range(ROWS):
    for c in range(SWARM_SIZE/ROWS):
        a += [(r+1)*s+r*bx+int(wx%s>0)]
        b += [(c+1)*s+c*bx+int(wx%s>0)]

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


def empty(i, j, k):
    if ((world[k, i*wx+j] == 0) 
    and (world[k-1, i*wx+j] in [0, DEBRIS_MAT]) and (world[k+1, i*wx+j] in [0, DEBRIS_MAT])
    and (world[k, (i+1)*wx+j] == 0) and (world[k, (i-1)*wx+j] == 0) 
    and (world[k, i*wx+j+1] == 0) and (world[k, i*wx+j-1] == 0) ):
        return True
    else:
        return False


for i in range(2, wx, SPACE_BETWEEN_DEBRIS+1): 
    for j in range(2, wy, SPACE_BETWEEN_DEBRIS+1):
        # for k in range(DEBRIS_HEIGHT):
        k = DEBRIS_HEIGHT-1
        try:
            if empty(i, j, k):
                world[k, i*wx+j] = DEBRIS_MAT  # pellet

            if empty(i+1, j+1, k/2) and HIGH_DEBRIS_CONCENTRATION:
                world[k/2, (i+1)*wx+(j+1)] = DEBRIS_MAT  # pellet

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
