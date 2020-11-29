from lxml import etree


def tensor_to_cdata(X):
    """Converts a tensor representation of voxels to a cdata representation.

    A tensor representation of voxels in R^4 is defined by:
        f(x, y, z) = m
    where m is the material type in which 0 represents empty space.
    The corresponding cdata representation transforms the z axis
    to the first dimension and stacks the data column-major.

    Example (four voxel pillar with gap):
        X = [[[0, 0, 1],[0, 0, 1],[0, 0, 0],[0, 0, 1]]]
        C = [[1],[1],[0],[1]]

    """
   
    return X.reshape(X.shape[0], -1)


def add_cdata_to_xml(C, x_size, y_size, z_size, file_path):
    """Writes the extent and CData to XML file.

    """

    VXD = etree.Element("VXD")
    Structure = etree.SubElement(VXD, "Structure")
    # Does this need to happen?
    Structure.set("replace", "VXA.VXC.Structure")
    Structure.set("Compression", "ASCII_READABLE")
    etree.SubElement(Structure, "X_Voxels").text = f"{x_size}"
    etree.SubElement(Structure, "Y_Voxels").text = f"{y_size}"
    etree.SubElement(Structure, "Z_Voxels").text = f"{z_size}"
    
    Data = etree.SubElement(Structure, "Data")
    for i in range(z_size):
        string = "".join([f"{int(c)}" for c in C[i, :]])
        etree.SubElement(Data, "Layer").text = etree.CDATA(string)

    file_content = etree.tostring(VXD, pretty_print=True).decode("utf-8")
    with open(file_path, "w") as f:
        print(file_content, file=f)

    # PhaseOffset = etree.SubElement(Structure, "PhaseOffset")
    # for i in range(Z_Voxels):
    #     string = ",".join([f"{c}" for c in phaseoffset_flatten[:,i].reshape(-1)])
    #     etree.SubElement(PhaseOffset, "Layer").text = etree.CDATA(string)
    # 
