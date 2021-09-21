# Thanks to Caitlin Grasso for the following code:
import numpy as np

def front(x,y,z,d=1):
    return x,y-d,z 

def back(x,y,z,d=1):
    return x,y+d,z

def left(x,y,z,d=1):
    return x-d,y,z

def right(x,y,z,d=1):
    return x+d,y,z

def get_empty_neighbor_positons(body, pos):
    empty_neigh = []
    for direction in [front, back, left, right]:
        neigh_pos = direction(pos[0], pos[1], pos[2])

        # Checking if the neighboring voxels is in array bounds
        if neigh_pos[0]>=0 and neigh_pos[0]<body.shape[0] and neigh_pos[1]>=0 and neigh_pos[1]<body.shape[1]:

            if body[neigh_pos] == 0: # in bounds and empty
                empty_neigh.append(neigh_pos)

        else: # out of array bounds so by default there is an empty neighbor
            empty_neigh.append(neigh_pos)

    return empty_neigh

def restricted_cilia(body, DEBUG=False):
    # Assumes voxels do not have cilia forces in the z-direction

    cilia = np.zeros((body.shape[0], body.shape[1], body.shape[2], 3))

    rad45 = (45/180)*np.pi
    
    # iterate through ciliated cells
    for x in range(body.shape[0]):
        for y in range(body.shape[1]):
            for z in range(body.shape[2]):
                
                if body[x,y,z]:
                    curr_pos = (x,y,z)
                    # Get neighboring empty voxel locations
                    empty_neigh = get_empty_neighbor_positons(body, curr_pos)

                    # Compute vectors to directions of the empty neighbors
                    vectors = []
                    for empty_neigh_pos in empty_neigh:

                        # voxels can push only
                        x_comp = curr_pos[0]-empty_neigh_pos[0]
                        y_comp = curr_pos[1]-empty_neigh_pos[1]
                        z_comp = curr_pos[2]-empty_neigh_pos[2] # should always be 0
                        assert z_comp==0
                        
                        # by default all of the vectors are unit vectors because the distance between voxels is 1
                        vectors.append([x_comp,y_comp,z_comp])

                    # DEBUG keeps thrusters in the center
                    if DEBUG:
                        if len(vectors)>0:
                            cilia[x,y,z,:] = vectors[0]
                        continue

                    # Compute range of angles the cilia force vector can lie in
                    # +/- 45 degrees of the vector to each empty neighboring voxel
                    if len(vectors)>0:
                        bounds = []
                        for vector in vectors:

                            angle_in_radians = np.arctan2(vector[1],vector[0])
                            lb = angle_in_radians-rad45 # lower bound
                            ub = angle_in_radians+rad45 # upper bound
                            
                            # convert angle to positive degrees around the unit circle if the angles are negative
                            # print(lb,ub)
                            if lb < 0 and ub < 0: 
                                lb = 2*np.pi + lb
                                ub = 2*np.pi + ub
                                bounds.append([lb,ub])

                            elif lb<0 and ub > 0:
                                lb = 2*np.pi + lb
                                bounds.append([lb,np.pi*2])
                                bounds.append([0,ub])

                            elif lb>0 and ub<0:
                                ub = 2*np.pi + ub
                                bounds.append([ub,np.pi*2])
                                bounds.append([0,lb])

                            else:
                                bounds.append([lb,ub])

                        # print(bounds)

                        # first choose a random range to choose from if more than one 
                        # important if the ranges are not touching
                        # i.e. missing voxels to the left and right but not up and down
                        range_index = np.random.randint(len(bounds))
                        angle_range = bounds[range_index]

                        # Choose a random angle in the range (on the unit circle)
                        cilia_force_angle = (angle_range[1] - angle_range[0]) * np.random.random() + angle_range[0]
                        
                        # Compute the x and y components of the unit vector given the chosen angle
                        cilia_x_comp = np.sin(cilia_force_angle)
                        cilia_y_comp = np.cos(cilia_force_angle)
                        cilia_z_comp = 0
                        cilia[x,y,z,:] = [cilia_x_comp, cilia_y_comp, cilia_z_comp]              

    return cilia