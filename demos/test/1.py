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

world = np.zeros([30,30,200], dtype=int)
world_cilia = np.zeros(shape=[3, world.shape[0], world.shape[1], world.shape[2]], dtype=float)
Senario="trivial"
# Senario="shoot"
if Senario=="trivial":
    body = np.ones(shape=[1,1,1], dtype=int)
    # body[0,0,1] = 0
    put_into(world, body)
elif Senario=="two links together":
    # test ok
    x,y,z = 2,1,4
    body_tiny = np.ones(shape=[x,y,z], dtype=int)
    body_tiny[0,0,1] = body_tiny[1,0,2] =0
    tiny_cilia = np.zeros(shape=[3,x,y,z])
    # tiny_cilia = np.random.random(size=[3,x,y,z]) * 0.2
    put_into(world, body_tiny)
    put_into(world_cilia, tiny_cilia, prefix=1)
elif Senario=="fall giggling":
    # 
    x,y,z = 2,2,5
    body_tiny = np.zeros(shape=[x,y,z], dtype=int)
    body_tiny[0,0,0] = body_tiny[0,0,1] = body_tiny[1,0,1] = body_tiny[1,1,4] =1
    # tiny_cilia = np.random.random(size=[3,x,y,z]) * 0.2
    put_into(world, body_tiny)
elif Senario=="regular_fall":
    # group calculating error
    x,y,z = 5,5,5
    body_tiny = np.ones(shape=[x,y,z], dtype=int)
    for i,x in enumerate(np.nditer(body_tiny, op_flags=['readwrite'])):
        x[...] = (i%2==0)
    # body_tiny[:,:,3:] = 0
    tiny_cilia = np.zeros(shape=[3,x,y,z])
    # tiny_cilia = np.random.random(size=[3,x,y,z]) * 0.2
    put_into(world, body_tiny)
    put_into(world_cilia, tiny_cilia, prefix=1)
elif Senario=="voxels fall":
    # group calculating error
    x,y,z = 20,20,200
    body_tiny = np.ones(shape=[x,y,z], dtype=int)
    body_tiny[np.random.random(size=[x,y,z])<0.9] = 0
    # body_tiny[:,:,3:] = 0
    tiny_cilia = np.zeros(shape=[3,x,y,z])
    # tiny_cilia = np.random.random(size=[3,x,y,z]) * 0.2
    put_into(world, body_tiny)
    put_into(world_cilia, tiny_cilia, prefix=1)
elif Senario=="shoot":
    # 
    body_tiny = np.ones(shape=[5,1,1], dtype=int)
    v = np.ones(shape=[1,1,1], dtype=int)
    c = np.zeros(shape=[3,1,1,1])
    c[1,0,0,0] = -0.4
    c[0,0,0,0] = -0.2

    # tiny_cilia = np.random.random(size=[3,2,2,3]) * 0.5
    put_into(world, body_tiny)
    put_into(world, v, offset=[2,2,0])
    put_into(world_cilia, c, offset=[2,2,0], prefix=1)

elif Senario=="tiny blob":
    # test ok
    body_tiny = np.ones(shape=[2,2,3], dtype=int)
    # body_tiny[np.random.random(size=[2,2,3])<0.3] = 0
    body_tiny[0,:,0]=0
    tiny_cilia = np.zeros(shape=[3,2,2,3])
    # tiny_cilia = np.random.random(size=[3,2,2,3]) * 0.5
    put_into(world, body_tiny)
    put_into(world_cilia, tiny_cilia, prefix=1)

elif Senario=="normal attach":
    # test ok
    X_Voxels, Y_Voxels, Z_Voxels = 2,3,1
    body = np.ones(shape=[X_Voxels,Y_Voxels,Z_Voxels], dtype=int)
    body[1,1,0] = 0
    cilia = np.zeros(shape=[3, X_Voxels, Y_Voxels, Z_Voxels], dtype=float)
    cilia[0,0,1,0] = 0.6

    put_into(world, body)
    put_into(world, body[::-1,:,:], offset=[4,0,0])
    put_into(world_cilia, cilia, prefix=1)
    put_into(world_cilia, -cilia[:,::-1,:,:], offset=[4,0,0], prefix=1)

elif Senario=="verticle attach":
    # test ok
    X_Voxels, Y_Voxels, Z_Voxels = 2,3,1
    body = np.ones(shape=[X_Voxels,Y_Voxels,Z_Voxels], dtype=int)
    body[1,1,0] = 0
    cilia = np.zeros(shape=[3, X_Voxels, Y_Voxels, Z_Voxels], dtype=float)
    cilia[1,0,1,0] = 0.6

    put_into(world, body)
    put_into(world, body[::-1,:,:], offset=[4,0,0])
    put_into(world_cilia, cilia, prefix=1)
    put_into(world_cilia, -cilia[:,::-1,:,:], offset=[4,0,0], prefix=1)
    world = np.swapaxes(world, 0,1)
    world_cilia = np.swapaxes(world_cilia,1,2)

elif Senario=="two Y_NEG":
    # still bad
    X_Voxels, Y_Voxels, Z_Voxels = 4,1,1
    body = np.ones(shape=[X_Voxels,Y_Voxels,Z_Voxels], dtype=int)
    body1 = np.ones(shape=[1,1,1], dtype=int)
    cilia = np.zeros(shape=[3, X_Voxels, Y_Voxels, Z_Voxels], dtype=float)
    cilia[1,0,0,0] = 0.18

    put_into(world, body)
    put_into(world, body1, offset=[5,0,0])
    put_into(world_cilia, cilia, prefix=1)

elif Senario=="Orthogonal attach":
    # still bad
    X_Voxels, Y_Voxels, Z_Voxels = 2,4,1
    body = np.ones(shape=[X_Voxels,Y_Voxels,Z_Voxels], dtype=int)
    # body[1,1,0] = 0
    cilia = np.zeros(shape=[3, X_Voxels, Y_Voxels, Z_Voxels], dtype=float)
    cilia[0,0,0,0] = 0.6

    cilia1 = np.zeros(shape=[3, X_Voxels, Y_Voxels, Z_Voxels], dtype=float)
    cilia1[0,1,-1,0] = 0

    put_into(world, body)
    put_into(world, body[::-1,:,:], offset=[4,0,0])
    put_into(world_cilia, cilia, prefix=1)
    put_into(world_cilia, cilia1, offset=[4,0,0], prefix=1)

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
