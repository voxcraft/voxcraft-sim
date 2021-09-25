from lxml import etree
import subprocess as sub
import numpy as np
import sys
from cilia_utils import restricted_cilia
from shape_utils import make_one_shape_only

# inputs
# 1: seed

SEED = 0
np.random.seed(int(sys.argv[1]))

N_CUTS = 4
N_PATCHES = 4

RECORD_HISTORY = True

DEBUG = False  # straight cilia vectors, instead of random angles

bx, by, bz = (9, 9, 7)
wx, wy, wz = (50, 9, 7)

USING_LIGHT_SOURCE = False
CILIA_DELAY_IN_LIGHT = 0
CILIA_DECAY_IN_DARK = 0
CILIA_FACTOR_IN_LIGHT = 1

# data
world = np.zeros((wx, wy, wz), dtype=np.int8)

# controller
BASE_CILIA_FORCE = np.zeros((wx, wy, wz, 3))

if USING_LIGHT_SOURCE:
    l_size = 2
    lx = wx//2-1 # light source min x
    ly = wy//2-1 # min y
    lz = wz-l_size # min z
    LIGHT_BULB = np.ones((l_size,)*3, dtype=np.int)*3  # materials: cilia, no cilia, lightbulb
    world[lx:lx+l_size, ly:ly+l_size, lz:lz+l_size] = LIGHT_BULB 

# morpholoy
body = np.ones((bx, by, bz), dtype=np.int8)

sphere = np.zeros((by+2,)*3, dtype=np.int8) 
radius = by//2+1
r2 = np.arange(-radius, radius+1)**2
dist2 = r2[:, None, None] + r2[:, None] + r2
sphere[dist2 <= radius**2] = 1

# remove the min and max layers and as many middle layers as necessary
for layer in range(bz):
    if layer > bz//2:
        pad = 1+by-bz
    else:
        pad = 1
    body[:, :, layer] *= sphere[1:bx+1, 1:by+1, layer+pad]

# carve away random chunks
for cut in range(N_CUTS):
    noise = np.random.rand(bx,by,bz)
    chunk = make_one_shape_only(noise)
    chunk = chunk.astype(np.int8)
    chunk *= body
    if np.sum(body) > np.sum(chunk):
        body -= chunk

# material distribution
for patch in range(N_PATCHES):
    noise = np.random.rand(bx,by,bz)
    chunk = make_one_shape_only(noise)
    chunk = chunk.astype(np.int8)
    chunk *= body
    body[chunk > 1] = 2

# shift down until in contact with surface plane
while True:
    if np.sum(body[:, :, 0]) == 0:
        body[:, :, :-1] = body[:, :, 1:]
        body[:, :, -1] = np.zeros_like(body[:, :, -1])
    else:
        break

for x in range(0, wx, bx+1):
    world[x:x+bx, :by, :bz] = body  # same body
    BASE_CILIA_FORCE[x:x+bx, :by, :bz] = restricted_cilia(body, DEBUG)  # different cilia

world = np.swapaxes(world, 0,2)
# world = world.reshape([wz,-1])
world = world.reshape(wz, wx*wy)

# get new voxcraft build
sub.call("cp ../../build/voxcraft-sim .", shell=True)
sub.call("cp ../../build/vx3_node_worker .", shell=True)

# create data folder if it doesn't already exist
sub.call("mkdir data{}".format(SEED), shell=True)
sub.call("cp base.vxa data{}/base.vxa".format(SEED), shell=True)

# clear old .vxd robot files from the data directory
sub.call("rm data{}/*.vxd".format(SEED), shell=True)

# delete old hist file
sub.call("rm a.hist", shell=True)

# delete old workspace
sub.call("rm -r workspace", shell=True)

# remove old sim output.xml if we are saving new stats
if not RECORD_HISTORY:
    sub.call("rm output{}.xml".format(SEED), shell=True)

# start vxd file
root = etree.Element("VXD")

vxa_using_light = etree.SubElement(root, "UsingLightSource")
vxa_using_light.set('replace', 'VXA.Simulator.UsingLightSource')
vxa_using_light.text = str(int(USING_LIGHT_SOURCE))

if USING_LIGHT_SOURCE:
    vxa_cilia_delay = etree.SubElement(root, "CiliaDelayInLight")
    vxa_cilia_delay.set('replace', 'VXA.Simulator.CiliaDelayInLight')
    vxa_cilia_delay.text = str(CILIA_DELAY_IN_LIGHT)

    vxa_cilia_decay = etree.SubElement(root, "CiliaDecayInDark")
    vxa_cilia_decay.set('replace', 'VXA.Simulator.CiliaDecayInDark')
    vxa_cilia_decay.text = str(CILIA_DECAY_IN_DARK)

    vxa_cilia_factor = etree.SubElement(root, "CiliaFactorInLight")
    vxa_cilia_factor.set('replace', 'VXA.Simulator.CiliaFactorInLight')
    vxa_cilia_factor.text = str(CILIA_FACTOR_IN_LIGHT)

    vxa_light_pos_x = etree.SubElement(root, "LightPosX")
    vxa_light_pos_x.set('replace', 'VXA.Simulator.LightPosX')
    vxa_light_pos_x.text = str(lx+l_size/2-0.5)

    vxa_light_pos_y = etree.SubElement(root, "LightPosY")
    vxa_light_pos_y.set('replace', 'VXA.Simulator.LightPosY')
    vxa_light_pos_y.text = str(ly+l_size/2-0.5)

    vxa_light_pos_z = etree.SubElement(root, "LightPosZ")
    vxa_light_pos_z.set('replace', 'VXA.Simulator.LightPosZ')
    vxa_light_pos_z.text = str(lz+l_size/2-0.5)

    print("light pos: " + str(lx+l_size/2-0.5) + ", " + str(ly+l_size/2-0.5) + ", " + str(lz+l_size/2-0.5) )


if RECORD_HISTORY:
    # sub.call("rm a{0}_gen{1}.hist".format(seed, pop.gen), shell=True)
    history = etree.SubElement(root, "RecordHistory")
    history.set('replace', 'VXA.Simulator.RecordHistory')
    etree.SubElement(history, "RecordStepSize").text = '100'
    etree.SubElement(history, "RecordVoxel").text = '1'
    etree.SubElement(history, "RecordLink").text = '0'
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