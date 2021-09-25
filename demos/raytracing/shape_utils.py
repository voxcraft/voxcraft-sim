import numpy as np
from scipy.ndimage.measurements import label


def make_circle(diameter):
    circle = np.zeros((diameter+2,)*2, dtype=np.int8) 
    radius = diameter//2+1
    r2 = np.arange(-radius, radius+1)**2
    dist2 = r2[:, None] + r2
    circle[dist2 <= radius**2] = 1
    return circle[1:-1, 1:-1, 1:-1]


def make_sphere(diameter):
    sphere = np.zeros((diameter+2,)*3, dtype=np.int8) 
    radius = diameter//2+1
    r2 = np.arange(-radius, radius+1)**2
    dist2 = r2[:, None, None] + r2[:, None] + r2
    sphere[dist2 <= radius**2] = 1
    return sphere[1:-1, 1:-1, 1:-1]


def make_one_shape_only(output_state):
    """Find the largest continuous arrangement of True elements after applying boolean mask.
    Avoids multiple disconnected softbots in simulation counted as a single individual.
    Parameters
    ----------
    output_state : numpy.ndarray
        Network output
    Returns
    -------
    part_of_ind : bool
        True if component of individual
    """
    if np.sum(output_state) == 0:
        return output_state

    # find coordinates
    array = output_state > 0
    labeled, ncomponents = label(array)

    largest_count = 0
    largest_label = 0
    for n in range(ncomponents+1):
        this_count = np.sum(labeled == n)
        vox_count = np.sum(array[labeled == n])
        if (this_count > largest_count) and (vox_count > 0):
            largest_label = n
            largest_count = this_count

    return labeled == largest_label


