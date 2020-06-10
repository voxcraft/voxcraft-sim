import os
from lxml import etree

import numpy as np
np.random.seed(1)

print("Start")

os.system("sh rebuild.sh")

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

def generate_structure(VXD, body,  cilia):
    X_Voxels = body.shape[0]
    Y_Voxels = body.shape[1]
    Z_Voxels = body.shape[2]

    body = np.swapaxes(body, 0,1)
    cilia = np.swapaxes(cilia, 1,2)
    
    body_flatten = body.reshape(X_Voxels*Y_Voxels, Z_Voxels)
    cilia_flatten = cilia.reshape(3, X_Voxels*Y_Voxels, Z_Voxels)
    
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

    BaseCiliaForce = etree.SubElement(Structure, "BaseCiliaForce")
    for i in range(Z_Voxels):
        string = ", ".join([f"{c}" for c in np.swapaxes(cilia_flatten[:,:,i], 0,1).reshape(-1)])
        etree.SubElement(BaseCiliaForce, "Layer").text = etree.CDATA(string)
    
def generate_vxd(body, cilia):
    VXD = etree.Element("VXD")
    generate_structure(VXD, body, cilia)
    file_content = etree.tostring(VXD, pretty_print=True).decode("utf-8")
    return file_content

X_Voxels, Y_Voxels, Z_Voxels = 2,3,1
body = np.ones(shape=[X_Voxels,Y_Voxels,Z_Voxels], dtype=int)
# body[1,1,0] = 0
cilia = np.zeros(shape=[3, X_Voxels, Y_Voxels, Z_Voxels], dtype=float)
cilia[0,0,1,0] = 1

cilia1 = np.zeros(shape=[3, X_Voxels, Y_Voxels, Z_Voxels], dtype=float)
cilia1[0,1,1,0] = -1

world = np.zeros([10,10,1], dtype=int)
world_cilia = np.zeros(shape=[3, world.shape[0], world.shape[1], world.shape[2]], dtype=float)

put_into(world, body)
put_into(world, body[::-1,:,:], offset=[6,0,0])
put_into(world_cilia, cilia, prefix=1)
put_into(world_cilia, cilia1, offset=[6,0,0], prefix=1)


print("Generating")
file_content = generate_vxd(world, world_cilia)
print("Writing to file")
with open("robot.vxd", "w") as f:
    print(file_content, file=f)
# print("==[data/robot.vxd]==")
# print(file_content)

print("Executing sim...")
os.system("./voxcraft-sim -i . > a.history")

c = os.stat('a.history').st_size
print(f"a.history file size: {c}")
if c<1000:
    print("==[a.history]==")
    with open("a.history", "r") as f:
        content = f.read()
        print(content)
else:
    print("Executing viz...")
    os.system("voxcraft-viz a.history")
