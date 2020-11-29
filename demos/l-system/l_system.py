from collections import deque
from mlp import MLP
import torch
import numpy as np
import itertools


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
        self.searched = False
        self.processed = False

    def __str__(self):
        return f"Voxel with material: {self.material} and level: {self.level}"


class PytorchGrowthFunction():

    def __init__(self, input_size, output_size):
        self.growth_function = MLP(input_size, output_size)
        
    def predict(self, X):
        """Predict the configuration index based on nearby voxels.

        Attributes:
            X: List of material proportions.
            
        Returns:
            index of configuration as an int.

        """

        X = torch.tensor([X])
        probabilities = self.growth_function(X).squeeze()
        max_value = 0
        index = 0
        for i in range(probabilities.shape[0]):
            if probabilities[i] > max_value:
                max_value = probabilities[i]
                index = i
        return int(index)


class ProbabilisticLSystem:
    """Generate patterns given patterns.

    Use the local context to decide what pattern to generate next.
    IE the the configuration of voxels added depend on the proportion
    of the voxel types.

    """

    def __init__(
        self,
        materials=(0, 1, 2),
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
        self.axiom = Voxel(1)
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
                material = c[0]
                direction = c[1]
                voxel = Voxel(material)
                if direction == "negative_x":
                    current_voxel.negative_x = voxel
                    voxel.positive_x = current_voxel
                if direction == "positive_x":
                    current_voxel.positive_x = voxel
                    voxel.negative_x = current_voxel
                if direction == "negative_y":
                    current_voxel.negative_y = voxel
                    voxel.positive_y = current_voxel
                if direction == "positive_y":
                    current_voxel.positive_y = voxel
                    voxel.negative_y = current_voxel
                if direction == "negative_z":
                    current_voxel.negative_z = voxel
                    voxel.positive_z = current_voxel
                if direction == "positive_z":
                    current_voxel.positive_z = voxel
                    voxel.negative_z = current_voxel
                voxel.level = current_voxel.level + 1
                voxels.append(voxel)
            print("Added voxels:")
            for v in voxels:
                print(v)
            return voxels

        self.axiom.level = 0
        body = deque([self.axiom])
        while len(body) > 0:
            voxel = body.pop()
            if voxel.level == self.growth_iterations:
                break
            X = self.get_function_input(voxel)
            configuration_index = growth_function.predict(X)
            configuration = self.configuration_map[configuration_index]
            voxels = attach_voxels(configuration, voxel)
            for v in voxels:
                body.appendleft(v)

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
            voxel.searched = True

            if np.abs(initial_level - voxel.level) > self.search_radius:
                break

            material_proportions[voxel.material] += 1
            if voxel.negative_x and not voxel.negative_x.searched:
                search_voxels.appendleft(voxel.negative_x)
            if voxel.positive_x and not voxel.positive_x.searched:
                search_voxels.appendleft(voxel.positive_x)
            if voxel.negative_y and not voxel.negative_y.searched:
                search_voxels.appendleft(voxel.negative_y)
            if voxel.positive_y and not voxel.positive_y.searched:
                search_voxels.appendleft(voxel.positive_y)
            if voxel.negative_z and not voxel.negative_z.searched:
                search_voxels.appendleft(voxel.negative_z)
            if voxel.positive_z and not voxel.positive_z.searched:
                search_voxels.appendleft(voxel.positive_z)
            total_voxels += 1

        for m in material_proportions:
            material_proportions[m] /= total_voxels
        return list(material_proportions.values())

    def initialize_configurations(self):
        """Map every possible configuration.

        Create a map from i = 0 to n of every possible way in which
        voxels with materials could be placed on three-dimensional
        surfaces.

        """
        # Get all possible voxels.
        possible_voxels = []
        for m in self.materials[1:]:
            for d in self.directions:
                possible_voxels.append((m, d))

        self.configuration_map = {}
        self.configuration_map[0] = None
        i = 0
        for num_voxels in range(1, self.max_voxels + 1):
            for subset in itertools.combinations(possible_voxels, num_voxels):
                self.configuration_map[i] = subset
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
            voxel.processed = True
            print(f"Added voxel at {x}, {y}, {z} with material {voxel.material}")
            X[x, y, z] = voxel.material

            if voxel.negative_x and not voxel.negative_x.processed:
                x -= 1
                to_process.appendleft((x, y, z, voxel.negative_x))
            if voxel.positive_x and not voxel.positive_x.processed:
                x += 1
                to_process.appendleft((x, y, z, voxel.positive_x))
            if voxel.negative_y and not voxel.negative_y.processed:
                y -= 1
                to_process.appendleft((x, y, z, voxel.negative_y))
            if voxel.positive_y and not voxel.positive_y.processed:
                y += 1
                to_process.appendleft((x, y, z, voxel.positive_y))
            if voxel.negative_z and not voxel.negative_z.processed:
                z -= 1
                to_process.appendleft((x, y, z, voxel.negative_z))
            if voxel.positive_z and not voxel.positive_z.processed:
                z += 1
                to_process.appendleft((x, y, z, voxel.positive_z))
        return X
