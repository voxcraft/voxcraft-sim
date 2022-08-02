from lxml import etree
import subprocess as sub
import numpy as np
import sys
from cilia_utils import restricted_cilia

# python demo.py 0 50 6 12 3 3 2000 10 0
# python demo.py 0 100 15 12 11 9 1000 10 0

# uses base_demo.vxa (line 103)

# inputs
# 1: DEBUG cilia: 0=random push angle, 1=exactly perpendicular
# 2: world size
# 3: world height
# 4: num of bots
# 5: bot length and width
# 6: bot height
# 7: cilia time delay in time steps when in light until full "charge" (100% of cilia factor is applied)
# 8: cilia decay % of 1 to remove from charge each time step when in dark
# 9: cilia factor when fully charged 

SEED = 1
np.random.seed(SEED)

RECORD_HISTORY = True

DEBUG = int(sys.argv[1])  # straight cilia vectors, instead of random angles

WORLD_SIZE = int(sys.argv[2])
WORLD_HEIGHT = int(sys.argv[3])
BODY_SIZES =  int(sys.argv[4]) * [(int(sys.argv[5]), int(sys.argv[5]), int(sys.argv[6])),]  # + [(2,)*3]*200  # + [(11, 11, 9),]*2  + [(5, 5, 4),]*10
# if body size changes, or if the stiffness/density of body material changes, 
# then the cilia force of the material will need to be recalibrated
wx, wy, wz = (WORLD_SIZE, WORLD_SIZE, WORLD_HEIGHT)

CILIA_DELAY_IN_LIGHT = int(sys.argv[7])
CILIA_DECAY_IN_DARK = float(sys.argv[8])/100
CILIA_FACTOR_IN_LIGHT = int(sys.argv[9])

# controller
BASE_CILIA_FORCE = np.zeros((wx, wy, wz, 3))
# BASE_CILIA_FORCE = np.ones((wx, wy, wz, 3))  * -1  # pointing downward
# BASE_CILIA_FORCE[:, :, :, :2] = 2 * np.random.rand(wx, wy, wz, 2) - 1  # unrestricted forces

# light source corner
l_size = 2
lx = wx//2-1
ly = wy//2-1
lz = wz-l_size
LIGHT_BULB = np.ones((l_size,)*3, dtype=np.int)*3  # materials: cilia, no cilia, lightbulb

# data
world = np.zeros((wx, wy, wz), dtype=np.int8)

world[lx:lx+l_size, ly:ly+l_size, lz:lz+l_size] = LIGHT_BULB 

for (bx, by, bz) in BODY_SIZES:
    body = np.ones((bx, by, bz), dtype=np.int8)

    if bx > 5:

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
            BASE_CILIA_FORCE[corners[0]:corners[0]+bx, corners[1]:corners[1]+by, :bz] = restricted_cilia(body, DEBUG)
            break
        if attepts > 500:
            break


world = np.swapaxes(world, 0,2)
# world = world.reshape([wz,-1])
world = world.reshape(wz, wx*wy)

# get new voxcraft build
sub.call("cp ../../build/voxcraft-sim .", shell=True)
sub.call("cp ../../build/vx3_node_worker .", shell=True)

# create data folder if it doesn't already exist
sub.call("mkdir data{}".format(SEED), shell=True)
sub.call("cp base_demo.vxa data{}/base.vxa".format(SEED), shell=True)

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
