from l_system import ProbabilisticLSystem, PytorchGrowthFunction
from vox.utils.tensor_to_cdata import tensor_to_cdata, add_cdata_to_xml
import argparse


parser = argparse.ArgumentParser(
    description="Generate a VXD file for the voxcraft simulation."
)
parser.add_argument("--growth-iterations", type=int)
parser.add_argument("--search-radius", type=int)
parser.add_argument("--max-voxels", type=int)
parser.add_argument("--record-history", 
                    action="store_true")

args = parser.parse_args()

l_system = ProbabilisticLSystem(
    growth_iterations=args.growth_iterations,
    search_radius=args.search_radius,
    max_voxels=args.max_voxels,
)

# The output size for the growth function is the
# number of possible configurations that can be grown.
input_size = len(l_system.materials)
output_size = len(l_system.configuration_map)

growth_function = PytorchGrowthFunction(input_size, output_size, sample=True)

# Grow the creature.
l_system.expand(growth_function)

# Prepare creature for simulation.
X = l_system.to_tensor()
C = tensor_to_cdata(X)
add_cdata_to_xml(C, X.shape[0], X.shape[1], X.shape[2], "data/robot.vxd",
                 record_history=args.record_history)
