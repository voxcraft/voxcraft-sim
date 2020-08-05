# For easy modification of VXA/VXD files
# Author: Sida Liu 2020

#
# Usage:
#
# import robot
# r = robot.Robot("in_filename")
# b = r.get_body()
# b is a numpy array, and we can modify it
# r.set_body(b)
# r.write("out_filename")
#

import numpy as np
from lxml import etree


class Robot:
    tree = None

    def __init__(self, filename):
        self.read(filename)

    def read(self, filename):
        with open(filename, "r") as f:
            self.tree = etree.parse(f)

    def get_body(self):
        X_Voxels = int(self.tree.find("*/X_Voxels").text)
        Y_Voxels = int(self.tree.find("*/Y_Voxels").text)
        Z_Voxels = int(self.tree.find("*/Z_Voxels").text)
        print(X_Voxels)
        data = self.tree.find("*/Data")
        layers = data.xpath("Layer")

        assert(len(layers) == Z_Voxels)

        ret = []
        print(len(layers))
        for layer in layers:
            text = layer.text
            assert(X_Voxels*Y_Voxels == len(text))
            arr = []
            for c in text:
                arr.append(int(c))
            arr = np.array(arr).reshape(X_Voxels, Y_Voxels)
            ret.append(arr)
        ret = np.array(ret).swapaxes(0,1).swapaxes(1,2)
        return ret

    def set_body(self, body):
        X_Voxels, Y_Voxels, Z_Voxels = body.shape
        body_flatten = body.reshape(-1, Z_Voxels)
        self.tree.find("*/X_Voxels").text = str(X_Voxels)
        self.tree.find("*/Y_Voxels").text = str(Y_Voxels)
        self.tree.find("*/Z_Voxels").text = str(Z_Voxels)
        Data = self.tree.find("*/Data")
        for child in list(Data):
            Data.remove(child)
        for i in range(Z_Voxels):
            string = "".join([f"{c}" for c in body_flatten[:,i]])
            etree.SubElement(Data, "Layer").text = etree.CDATA(string)
    
    def write(self, filename):
        self.tree.write(filename)