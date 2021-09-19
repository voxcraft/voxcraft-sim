from lxml import etree
import subprocess as sub
import numpy as np

SEED = 0
np.random.seed(SEED)

RECORD_HISTORY = True

WORLD_SIZE = 80
WORLD_HEIGHT = 5
BODY_SIZES = [(7, 7, 5),]*12 # (6, 6, 5)  # (8, 8, 7)
# if body size changes, or if the stiffness/density of body material changes, 
# then the cilia force of the material will need to be recalibrated
wx, wy, wz = (WORLD_SIZE, WORLD_SIZE, WORLD_HEIGHT)

EVAL_PERIOD = 1  # 4
SETTLE_TIME = 0  # 0.5

RANDMONIZE_CILIA_EVERY = 0.25 # 5

# controller
BASE_CILIA_FORCE = np.ones((wx, wy, wz, 3))  * -1  # pointing downward
BASE_CILIA_FORCE[:, :, :, :2] = 2 * np.random.rand(wx, wy, wz, 2) - 1  # initial forces

# light source corner
lx = 0
ly = 0
lz = 2
l_size = 2
LIGHT_BULB = np.ones((l_size,)*3, dtype=np.int)*3  # materials: cilia, no cilia, lightbulb

# data
world = np.zeros((wx, wy, wz), dtype=np.int8)

world[lx:lx+l_size, ly:ly+l_size, lz:lz+l_size] = LIGHT_BULB 

for (bx, by, bz) in BODY_SIZES:
    body = np.ones((bx, by, bz), dtype=np.int8)

    sphere = np.zeros((by+2,)*3, dtype=np.int8) 
    radius = by//2+1
    r2 = np.arange(-radius, radius+1)**2
    dist2 = r2[:, None, None] + r2[:, None] + r2
    sphere[dist2 <= radius**2] = 1

    # max_size = 0
    for layer in range(bz):
        if layer > bz//2:
            pad = (bz-1) - (by-bz)//2
        else:
            pad = (by-bz)//2
        body[:, :, layer] *= sphere[1:bx+1, 1:by+1, layer+pad]
        # max_size += np.sum(sphere[1:bx+1, 1:by+1, layer+pad])

    while True:  # shift down until in contact with surface plane
        if np.sum(body[:, :, 0]) == 0:
            body[:, :, :-1] = body[:, :, 1:]
            body[:, :, -1] = np.zeros_like(body[:, :, -1])
        else:
            break

    attepts = 0
    while True:
        attepts += 1
        corners = np.random.randint(l_size+1, wx-bx, 2)
        if np.sum(world[corners[0]-1:corners[0]+bx+1, corners[1]-1:corners[1]+by+1, :bz]) == 0:
            world[corners[0]:corners[0]+bx, corners[1]:corners[1]+by, :bz] = body
            break
        if attepts > 500:
            break


world = np.swapaxes(world, 0,2)
# world = world.reshape([wz,-1])
world = world.reshape(wz, wx*wy)

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

vxa_world_size = etree.SubElement(root, "LightPosX")
vxa_world_size.set('replace', 'VXA.Simulator.LightPosX')
vxa_world_size.text = str(lx+l_size//2-0.5)
vxa_world_size = etree.SubElement(root, "LightPosY")
vxa_world_size.set('replace', 'VXA.Simulator.LightPosX')
vxa_world_size.text = str(ly+l_size//2-0.5)
vxa_world_size = etree.SubElement(root, "LightPosZ")
vxa_world_size.set('replace', 'VXA.Simulator.LightPosX')
vxa_world_size.text = str(lz+l_size//2-0.5)
print("light pos: " + str(lx+l_size//2-0.5) + ", " + str(ly+l_size//2-0.5) + ", " + str(lz+l_size//2-0.5))

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
    etree.SubElement(history, "RecordLink").text = '1'
    etree.SubElement(history, "RecordFixedVoxels").text = '1'
    etree.SubElement(history, "RecordCoMTraceOfEachVoxelGroupfOfThisMaterial").text = '0'  # draw CoM trace
    

structure = etree.SubElement(root, "Structure")
structure.set('replace', 'VXA.VXC.Structure')
structure.set('Compression', 'ASCII_READABLE')
etree.SubElement(structure, "X_Voxels").text = str(wx)
etree.SubElement(structure, "Y_Voxels").text = str(wy)
etree.SubElement(structure, "Z_Voxels").text = str(wz)

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
