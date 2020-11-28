from collections import deque
import numpy as np
from vox.utils.tensor_to_cdata import tensor_to_cdata


class Voxel:
    """A voxel of material type, build level and connections.

    IE a graph representation of a voxel.

    """

    def __init__(self, material):
        self.material = material
        self.negative_x = None
        self.positive_x = None
        self.positive_z = None
        self.negative_z = None
        self.positive_y = None
        self.negative_y = None
        self.level = None


class ProbabilisticLSystem:
    """Generate patterns given patterns.

    Use the local context to decide what pattern to generate next.
    IE the the configuration of voxels added depend on the proportion
    of the voxel types.

    """

    def __init__(
        self,
        materials=(0, 1),
        directions=(
            "negative_x",
            "positive_x",
            "negative_y",
            "positive_y",
            "negative_z",
            "positive_z",
        ),
        growth_iterations=3,
        max_voxels=5,
        search_radius=1,
    ):
        self.axiom = Voxel(1, 0)
        self.growth_iterations = growth_iterations
        self.materials = materials
        self.directions = directions
        self.max_voxels = max_voxels
        self.initialize_configurations()
        self.search_radius = search_radius

    def expand(self, growth_function):
        """Expand the axiom and grow the body.

        Expand out the axiom given the condtional probability of
        new growth given nearby previous growth.

        """

        def attach_voxels(configuration, current_voxel):
            """Attach a configuration of voxels to the current voxel.

            Attach a configuration of voxels (IE a
            combination of voxels of a given material and placements)
            to to the current voxel.

            """

            voxels = []
            for c in configuration:
                direction = c[0]
                voxel = c[1]
                if direction == "negative_x":
                    current_voxel.negative_x = voxel
                    voxel.positive_x = current_voxel
                elif direction == "positive_x":
                    current_voxel.positive_x = voxel
                    voxel.negative_x = current_voxel
                elif direction == "negative_y":
                    current_voxel.negative_y = voxel
                    voxel.positive_y = current_voxel
                elif direction == "positive_y":
                    current_voxel.positive_y = voxel
                    voxel.negative_y = current_voxel
                elif direction == "negative_z":
                    current_voxel.negative_z = voxel
                    voxel.positive_z = current_voxel
                elif direction == "positive_z":
                    current_voxel.positive_z = voxel
                    voxel.negative_z = current_voxel
                else:
                    raise ValueError("Invalid direction.")
                voxel.level = current_voxel.level + 1
                voxels.append(voxel)
            return voxels

        self.axiom.level = 0
        body = deque([self.axiom])
        while len(body) > 0:
            voxel = body.pop()
            if voxel.level > self.growth_iterations:
                break
            X = self.get_function_input(voxel)
            configuration_number = growth_function(X)
            configuration = self.configuration[configuration_number]
            voxels = attach_voxels(configuration)
            body.append(voxels)

    def get_function_input(self, voxel):
        """Get the material proportions of nearby voxels"""

        initial_level = voxel.level
        total_voxels = 0
        material_proportions = {}
        for m in self.materials:
            material_proportions[m] = 0
        search_voxels = deque([voxel])
        while len(search_voxels) > 0:
            voxel = search_voxels.pop()
            if np.abs(initial_level - voxel.level) > self.search_radius:
                break
            material_proportions[voxel.material] += 1
            if voxel.negative_x:
                search_voxels.appendleft(voxel.negative_x)
            elif voxel.positive_x:
                search_voxels.appendleft(voxel.positive_x)
            elif voxel.negative_y:
                search_voxels.appendleft(voxel.negative_y)
            elif voxel.positive_y:
                search_voxels.appendleft(voxel.positive_y)
            elif voxel.negative_z:
                search_voxels.appendleft(voxel.negative_z)
            elif voxel.positive_z:
                search_voxels.appendleft(voxel.positive_z)
            total_voxels += 1
        for m in material_proportions:
            material_proportions[m] /= total_voxels
        return material_proportions.values()

    def initialize_configurations(self):
        """Map every possible configuration.

        Create a map from i = 0 to n of every possible way in which
        voxels with materials could be placed on three-dimensional
        surfaces.

        """

        self.configuration_map = {}
        self.configuration_map[0] = None
        i = 0
        for m in self.materials:
            for d in self.directions:
                for v in range(1, self.max_voxels + 1):
                    voxels = []
                    for _ in v:
                        voxels.appendleft((d, Voxel(m)))
                    self.configuration_map[i] = voxels
                    i += 1

    def to_tensor(self):
        """Convert the graph representation of the body to a tensor.

        Fill a three-dimensional tensor with the material types of
        each voxel IE:
            X[i][j][k] = m

        """

        extent = (2 * self.growth_iterations) + 1
        X = np.zeros((extent, extent, extent))
        middle = int(np.floor(extent / 2))
        x, y, z = middle, middle, middle
        to_process = deque([(x, y, z, self.axiom)])
        while len(to_process) > 0:
            x, y, z, voxel = to_process.pop()
            X[x, y, z] = voxel.material
            if voxel.negative_x:
                x -= 1
                to_process.appendleft((x, y, z, voxel.negative_x))
            elif voxel.positive_x:
                x += 1
                to_process.appendleft((x, y, z, voxel.positive_x))
            elif voxel.negative_y:
                y -= 1
                to_process.appendleft((x, y, z, voxel.negative_y))
            elif voxel.positive_y:
                y += 1
                to_process.appendleft((x, y, z, voxel.positive_y))
            elif voxel.negative_z:
                z -= 1
                to_process.appendleft((x, y, z, voxel.negative_z))
            elif voxel.positive_z:
                z += 1
                to_process.appendleft((x, y, z, voxel.positive_z))
        return X
