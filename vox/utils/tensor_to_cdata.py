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
