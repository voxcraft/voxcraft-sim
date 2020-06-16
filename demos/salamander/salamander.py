import os
from lxml import etree

import numpy as np
np.random.seed(1)
import binvox_rw

print("Start")

os.system("sh rebuild.sh")

with open('./salamander_100.binvox', 'rb') as f:
    model = binvox_rw.read_as_3d_array(f)

body = model.data
r = np.random.random(size=body.shape)
body[r<0.8] = 0

X_Voxels, Y_Voxels, Z_Voxels = body.shape

phaseoffset = np.random.random(size=body.shape)
cilia = np.zeros(shape=[3, X_Voxels, Y_Voxels, Z_Voxels], dtype=float)
cilia_shift = np.zeros_like(cilia)


world = np.zeros(body.shape, dtype=int)
world_phaseoffset = np.zeros_like(world, dtype=float)
world_cilia = np.zeros(shape=[3, world.shape[0], world.shape[1], world.shape[2]], dtype=float)
world_cilia_shift = np.zeros_like(world_cilia)

def put_into(world, data, offset = [0,0,0], prefix=0):
    x = offset[0]
    x1 = offset[0]+data.shape[0+prefix]
    y = offset[1]
    y1 = offset[1]+data.shape[1+prefix]
    z = offset[2]
    z1 = offset[2]+data.shape[2+prefix]
    if prefix==0:
        world[x:x1,y:y1,z:z1] = data
    else:
        world[:,x:x1,y:y1,z:z1] = data

put_into(world, body)
put_into(world_phaseoffset, phaseoffset)
put_into(world_cilia, cilia, prefix=1)
put_into(world_cilia_shift, cilia_shift, prefix=1)

def generate_structure(VXD, body, phaseoffset, cilia, cilia_shift):
    X_Voxels = body.shape[0]
    Y_Voxels = body.shape[1]
    Z_Voxels = body.shape[2]

    # TODO: tricky, too tricky, next time, please don't use 2-d to store 3-d/4-d data! confusing!
    # Please use good 1d structure to store n-d data!
    body = np.swapaxes(body, 0,1)
    phaseoffset = np.swapaxes(phaseoffset, 0,1)
    cilia = np.swapaxes(cilia, 1,2)
    cilia_shift = np.swapaxes(cilia_shift, 1,2)
    
    body_flatten = body.reshape(X_Voxels*Y_Voxels, Z_Voxels)
    phaseoffset_flatten = phaseoffset.reshape(X_Voxels*Y_Voxels, Z_Voxels)
    cilia_flatten = cilia.reshape(3, X_Voxels*Y_Voxels, Z_Voxels)
    cilia_shift_flatten = cilia_shift.reshape(3, X_Voxels*Y_Voxels, Z_Voxels)
    
    Structure = etree.SubElement(VXD, "Structure")
    Structure.set("replace", "VXA.VXC.Structure")
    Structure.set("Compression", "ASCII_READABLE")
    etree.SubElement(Structure, "X_Voxels").text = f"{X_Voxels}"
    etree.SubElement(Structure, "Y_Voxels").text = f"{Y_Voxels}"
    etree.SubElement(Structure, "Z_Voxels").text = f"{Z_Voxels}"
    
    Data = etree.SubElement(Structure, "Data")
    for i in range(Z_Voxels):
        string = "".join([f"{c}" for c in body_flatten[:,i]])
        etree.SubElement(Data, "Layer").text = etree.CDATA(string)

    PhaseOffset = etree.SubElement(Structure, "PhaseOffset")
    for i in range(Z_Voxels):
        string = ",".join([f"{c}" for c in phaseoffset_flatten[:,i].reshape(-1)])
        etree.SubElement(PhaseOffset, "Layer").text = etree.CDATA(string)

    # BaseCiliaForce = etree.SubElement(Structure, "BaseCiliaForce")
    # for i in range(Z_Voxels):
    #     string = ", ".join([f"{c}" for c in np.swapaxes(cilia_flatten[:,:,i], 0,1).reshape(-1)])
    #     etree.SubElement(BaseCiliaForce, "Layer").text = etree.CDATA(string)
    
    # ShiftCiliaForce = etree.SubElement(Structure, "ShiftCiliaForce")
    # for i in range(Z_Voxels):
    #     string = ", ".join([f"{c}" for c in np.swapaxes(cilia_shift_flatten[:,:,i], 0,1).reshape(-1)])
    #     etree.SubElement(ShiftCiliaForce, "Layer").text = etree.CDATA(string)
    
def generate_vxd(body, phaseoffset, cilia, cilia_shift):
    VXD = etree.Element("VXD")
    generate_structure(VXD, body, phaseoffset, cilia, cilia_shift)
    file_content = etree.tostring(VXD, pretty_print=True).decode("utf-8")
    return file_content

print("Generating")
file_content = generate_vxd(world, world_phaseoffset, world_cilia, world_cilia_shift)
print("Writing to file")
with open("data/robot.vxd", "w") as f:
    print(file_content, file=f)
# print("==[data/robot.vxd]==")
# print(file_content)

print("Executing sim...")
os.system("./voxcraft-sim -i data > a.history")

c = os.stat('a.history').st_size
print(f"a.history file size: {c}")
if c<1000:
    print("==[a.history]==")
    with open("a.history", "r") as f:
        content = f.read()
        print(content)
else:
    print("Executing viz...")
    os.system("./voxcraft-viz a.history")
