from l_system import ProbabilisticLSystem, PytorchGrowthFunction
from vox.utils.tensor_to_cdata import tensor_to_cdata, add_cdata_to_xml


l_system = ProbabilisticLSystem(growth_iterations=7, 
                                search_radius=2,
                                max_voxels=2)

# The output size for the growth function is the
# number of possible configurations that can be grown.
input_size = len(l_system.materials)
output_size = len(l_system.configuration_map)

growth_function = PytorchGrowthFunction(input_size, output_size)

# Grow the creature.
l_system.expand(growth_function)

# Prepare creature for simulation.
X = l_system.to_tensor()
C = tensor_to_cdata(X)
add_cdata_to_xml(C, X.shape[0], X.shape[1], X.shape[2], "data/robot.vxd")
